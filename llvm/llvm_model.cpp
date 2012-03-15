#include "llvm/llvm_model.hpp"
#include "interface/log2_dot.hpp"
#include "interface/plugin.hpp"
#include "interface/distribution.hpp"

#include <iostream>

#include "llvm/Module.h"
#include "llvm/DerivedTypes.h"
#include "llvm/LLVMContext.h"
#include "llvm/ExecutionEngine/ExecutionEngine.h"

#include "llvm/Support/IRBuilder.h"

using namespace std;
using namespace theta;
using namespace llvm;
using namespace theta::utils;

namespace {
    // conversion utilities fom theta to llvm data structures: the first three are for
    // converting data / histograms
    size_t get_total_nbins(const ObsIds & oids, const VarIdManager & vm){
        size_t total_nbins = 0;
        for(ObsIds::const_iterator it=oids.begin(); it!=oids.end(); ++it){
            total_nbins += vm.get_nbins(*it);
            if(total_nbins % 2) ++total_nbins;
        }
        return total_nbins;
    }
    
    // make the concatenated histogram, based on the observables in oids and binning info in vm.
    // assumes that h_concat has been allocated correctly with get_total_nbins(!!)
    void get_concatenated_from_data(Histogram1D & h_concat, const Data & dat, const ObsIds & oids){
        size_t offset = 0;
        for(ObsIds::const_iterator it=oids.begin(); it!=oids.end(); ++it){
            size_t nbins = dat[*it].get_nbins();
            memcpy(h_concat.get_data() + offset, dat[*it].get_data(), nbins * sizeof(double));
            offset += nbins;
            if(offset % 2) ++offset;
        }
    }
    
    void get_data_from_concatenated(Data & d, const Histogram1D & h_concat, const ObsIds & oids, const VarIdManager & vm){
        theta_assert(get_total_nbins(oids, vm) == h_concat.size());
        size_t offset = 0;
        for(ObsIds::const_iterator it=oids.begin(); it!=oids.end(); ++it){
            size_t nbins = vm.get_nbins(*it);
            const std::pair<double, double> & range = vm.get_range(*it);
            d[*it] = Histogram1D(nbins, range.first, range.second);
            memcpy(d[*it].get_data(), h_concat.get_data() + offset, nbins * sizeof(double));
            if(nbins % 2){
                theta_assert(*(h_concat.get_data() + offset + nbins) == 0);
            }
            offset += nbins;
            if(offset % 2) ++offset;
        }
    }
    
    // convert ParValues into double * ( use &par_values[0])
    void get_par_values(vector<double> & par_values, const ParValues & values, const ParIds & par_ids){
        par_values.resize(par_ids.size());
        size_t i=0;
        for(ParIds::const_iterator it=par_ids.begin(); it!=par_ids.end(); ++it, ++i){
            par_values[i] = values.get_unchecked(*it);
        }
    }
    
    
    // par_values must have correct size!
    void get_par_values(vector<double> & par_values, const ParValues & values, const vector<ParId> & par_ids){
        //par_values.resize(par_ids.size());
        size_t i=0;
        for(vector<ParId>::const_iterator it=par_ids.begin(); it!=par_ids.end(); ++it, ++i){
            par_values[i] = values.get_unchecked(*it);
        }
    }
}


class LLVMPluginBuilder: public PluginBuilder<theta::Function>{
public:
    virtual std::auto_ptr<theta::Function> build(const Configuration & cfg, const std::string & type){
        string new_type = type;
        if(type=="exp_function" || type == "multiply") new_type = "llvm_" + type;
        return PluginManager<theta::Function>::build(cfg, new_type);
    }
};

struct llvm_pb_sentinel{
   llvm_pb_sentinel(){
       std::auto_ptr<PluginBuilder<theta::Function> > pb(new LLVMPluginBuilder());
       PluginManager<theta::Function>::set_plugin_builder(pb);
   }
   ~llvm_pb_sentinel(){
       PluginManager<theta::Function>::reset_plugin_builder();
   }
};


llvm_model::~llvm_model(){}

void llvm_model::set_prediction(const ObsId & obs_id, boost::ptr_vector<theta::Function> & coeffs_, boost::ptr_vector<HistogramFunction> & histos_){
    observables.insert(obs_id);
    const size_t n = coeffs_.size();
    if(n!=coeffs_.size()) throw invalid_argument("Model::setPrediction: number of histograms and coefficients do not match");
    if(histos[obs_id].size()>0 || coeffs[obs_id].size()>0)
        throw invalid_argument("Model::setPrediction: prediction already set for this observable");
    coeffs[obs_id].transfer(coeffs[obs_id].end(), coeffs_.begin(), coeffs_.end(), coeffs_);
    histos[obs_id].transfer(histos[obs_id].end(), histos_.begin(), histos_.end(), histos_);
    for(boost::ptr_vector<theta::Function>::const_iterator it=coeffs[obs_id].begin(); it!=coeffs[obs_id].end(); ++it){
        ParIds pids = (*it).get_parameters();
        parameters.insert(pids.begin(), pids.end());
    }
    size_t nbins = 0;
    double xmin = NAN, xmax = NAN;
    bool first = true;
    for(boost::ptr_vector<HistogramFunction>::const_iterator it=histos[obs_id].begin(); it!=histos[obs_id].end(); ++it){
        if(first){
            it->get_histogram_dimensions(nbins, xmin, xmax);
            first = false;
        }
        else{
            size_t nbins_tmp = 0;
            double xmin_tmp = NAN, xmax_tmp = NAN;
            it->get_histogram_dimensions(nbins_tmp, xmin_tmp, xmax_tmp);
            if(nbins!=nbins_tmp || xmin!=xmin_tmp || xmax!=xmax_tmp){
                throw invalid_argument("llvm_model::set_prediction: histogram dimensions mismatch");
            }
        }
        const ParIds & pids = (*it).get_parameters();
        parameters.insert(pids.begin(), pids.end());
    }
}


template<typename HT>
void llvm_model::get_prediction_impl(DataT<HT> & result, const ParValues & parameters) const{
    histos_type::const_iterator h_it = histos.begin();
    coeffs_type::const_iterator c_it = coeffs.begin();
    for(; h_it != histos.end(); ++h_it, ++c_it){
        const ObsId & oid = h_it->first;
        histos_type::const_mapped_reference hfs = *(h_it->second);
        coeffs_type::const_mapped_reference h_coeffs = *(c_it->second);
        theta_assert(hfs.size() > 0 && hfs.size() == h_coeffs.size());
        // overwrite result[oid] with first term:
        hfs[0].apply_functor(copy_to<HT>(result[oid]), parameters);
        result[oid] *= h_coeffs[0](parameters);
        // add the rest:
        for (size_t i = 1; i < hfs.size(); i++) {
            hfs[i].apply_functor(add_with_coeff_to<HT>(result[oid], h_coeffs[i](parameters)), parameters);
        }
    }
}

void llvm_model::get_prediction(Data & result, const ParValues & parameters) const {
    if(llvm_always){
        if(model_get_prediction==0) generate_llvm();
        Histogram1D pred_concat(get_total_nbins(observables, vm));
        vector<double> parameter_values;
        get_par_values(parameter_values, parameters, this->parameters);
        model_get_prediction(&parameter_values[0], pred_concat.get_data());
        get_data_from_concatenated(result, pred_concat, observables, vm);
        return;
    }
    get_prediction_impl<Histogram1D>(result, parameters);
}

void llvm_model::get_prediction(DataWithUncertainties & result, const theta::ParValues & parameters) const{
    get_prediction_impl<Histogram1DWithUncertainties>(result, parameters);
}

// the generated function is
// void <prefix>_get_prediction(const double * par_values, double * data);
// where data must be set to zero prior to calling this function.
void llvm_model::generate_llvm() const {
    module.reset(new llvm_module(parameters));
    LLVMContext & context = module->module->getContext();
    const Type * double_t = Type::getDoubleTy(context);
    const Type * i32_t = Type::getInt32Ty(context);
    const Type * void_t = Type::getVoidTy(context);
    std::vector<const Type*> arg_types(2, double_t->getPointerTo());
    FunctionType * FT = FunctionType::get(void_t, arg_types, false);
    llvm::Function * F = llvm::Function::Create(FT, llvm::Function::ExternalLinkage, "model_evaluate", module->module);
    llvm::Function::arg_iterator iter = F->arg_begin();
    // the function parameters:
    Value * par_values = iter++;
    Value * data = iter;
    IRBuilder<> Builder(context);
    BasicBlock * BB = BasicBlock::Create(context, "entry", F);
    Builder.SetInsertPoint(BB);
    // generate code for all Functions and HistorgamFunctions:
    histos_type::const_iterator h_it = histos.begin();
    coeffs_type::const_iterator c_it = coeffs.begin();
    size_t iobs = 0;
    size_t obs_offset = 0;
    for(; h_it != histos.end(); ++h_it, ++c_it, ++iobs){
        theta_assert(h_it->first == c_it->first);
        histos_type::const_mapped_reference h_functions = *(h_it->second);
        coeffs_type::const_mapped_reference h_coeffs = *(c_it->second);
        theta_assert(h_functions.size() == h_coeffs.size());
        // data_iobs = data + obs_offset:
        Constant * index = ConstantInt::get(i32_t, obs_offset);
        GetElementPtrInst* data_iobs = GetElementPtrInst::Create(data, index, "", BB);
        for(size_t i=0; i<h_functions.size(); ++i){
            stringstream ss_prefix;
            ss_prefix << "c" << iobs << "_p" << i;
            llvm::Function * coeff_function = create_llvm_function(&h_coeffs[i], *module, ss_prefix.str());
            set_private_inline(coeff_function);
            Value * coeff = Builder.CreateCall(coeff_function, par_values);
            llvm::Function * histo_function = create_llvm_histogram_function(&h_functions[i], *module, ss_prefix.str());
            set_private_inline(histo_function);
            Builder.CreateCall3(histo_function, coeff, par_values, data_iobs);
        }
        obs_offset += vm.get_nbins(h_it->first);
        if(obs_offset % 2) ++obs_offset;
    }
    Builder.CreateRetVoid();
    //module->module->dump();
    module->optimize();
    //cout << "\n\nafter optimization:\n\n";
    //module->module->dump();
    model_get_prediction = reinterpret_cast<t_model_get_prediction>(module->getFunctionPointer(F));
}

std::auto_ptr<NLLikelihood> llvm_model::get_nllikelihood(const Data & data) const{
    if(not(data.get_observables()==observables)){
        throw invalid_argument("Model::createNLLikelihood: observables of model and data mismatch!");
    }
    if(model_get_prediction == 0){
        generate_llvm();
    }
    return std::auto_ptr<NLLikelihood>(new llvm_model_nll(*this, data, model_get_prediction, nbins_total));
}


llvm_model::llvm_model(const Configuration & ctx): vm(*ctx.pm->get<VarIdManager>()),
 llvm_always(false), model_get_prediction(0){
    SettingWrapper s = ctx.setting;
    if(s.exists("bb_uncertainties")){
        bool bb_uncertainties =  s["bb_uncertainties"];
        if(bb_uncertainties) throw ConfigurationException("llvm_model does not implement bb_uncertainties!");
    }
    ObsIds observables = vm.get_all_observables();
    llvm_pb_sentinel b;
    //go through observables to find the template definition for each of them:
    for (ObsIds::const_iterator obsit = observables.begin(); obsit != observables.end(); obsit++) {
        string obs_name = vm.get_name(*obsit);
        if(not s.exists(obs_name)) continue;
        SettingWrapper obs_setting = s[obs_name];
        boost::ptr_vector<HistogramFunction> histos;
        boost::ptr_vector<theta::Function> coeffs;
        for (size_t i = 0; i < obs_setting.size(); i++) {
            Configuration context(ctx, obs_setting[i]["histogram"]);
            auto_ptr<HistogramFunction> hf = PluginManager<HistogramFunction>::build(context);
            Configuration cfg(ctx, obs_setting[i]["coefficient-function"]);
            auto_ptr<theta::Function> coeff_function = PluginManager<theta::Function>::build(cfg);
            coeffs.push_back(coeff_function);
            theta_assert(coeff_function.get()==0);
            histos.push_back(hf);
        }
        set_prediction(*obsit, coeffs, histos);
    }
    parameter_distribution = PluginManager<Distribution>::build(Configuration(ctx, s["parameter-distribution"]));
    if(not (parameter_distribution->get_parameters() == get_parameters())){
        stringstream ss;
        ss << "'parameter-distribution' has to define the same set of parameters the model depends on. However";
        ParIds dist_pars = parameter_distribution->get_parameters();
        const ParIds & m_pars = parameters;
        ParIds all_pars = parameters;
        all_pars.insert(dist_pars.begin(), dist_pars.end());
        for(ParIds::const_iterator p_it=all_pars.begin(); p_it!=all_pars.end(); ++p_it){
            if(m_pars.contains(*p_it) && dist_pars.contains(*p_it)) continue;
            if(m_pars.contains(*p_it)){
               ss << ", the model depends on '"<< vm.get_name(*p_it) << "' which the parameter distribution does not include";
            }
            else ss << ", the parameter distribution depends on '" << vm.get_name(*p_it) << "' which the model does not depend on";
        }
        throw ConfigurationException(ss.str());
    }
    if(ctx.setting.exists("llvm_always")){
        llvm_always = ctx.setting["llvm_always"];
    }
    if(ctx.setting.exists("rvobs-distribution")){
        rvobservable_distribution = PluginManager<Distribution>::build(Configuration(ctx, ctx.setting["rvobs-distribution"]));
        rvobservables = rvobservable_distribution->get_parameters();
    }
    // type checking:
    for(ParIds::const_iterator it=parameters.begin(); it!=parameters.end(); ++it){
        if(vm.get_type(*it) != "par"){
            throw ConfigurationException("Type error: parameter '" + vm.get_name(*it) + "' is used as model parameter, but was not declared as such.");
        }
    }
    for(ParIds::const_iterator it=rvobservables.begin(); it!=rvobservables.end(); ++it){
        if(vm.get_type(*it) != "rvobs"){
            throw ConfigurationException("Type error: parameter '" + vm.get_name(*it) + "' is used as real-valued observable, but was not declared as such.");
        }
    }
    nbins_total = get_total_nbins(observables, vm);
}


/* llvm_model_nll */
void llvm_model_nll::fill_par_ids_vec(){
    par_ids_vec.clear();
    for(ParIds::const_iterator it = par_ids.begin(); it!=par_ids.end(); ++it){
        par_ids_vec.push_back(*it);
    }
    parameter_values.resize(par_ids_vec.size());
}


llvm_model_nll::llvm_model_nll(const llvm_model & m, const Data & dat, t_model_get_prediction model_get_prediction_, size_t nbins_total): model(m),
  data_concatenated(nbins_total, 0, 1), model_get_prediction(model_get_prediction_){
    theta::Function::par_ids = model.get_parameters();
    fill_par_ids_vec();
    get_concatenated_from_data(data_concatenated, dat, model.get_observables());
    pred_concatenated = Histogram1D(data_concatenated.size(), pred_concatenated.get_xmin(), pred_concatenated.get_xmax());
    rvobs_values = dat.get_rvobs_values();
}

void llvm_model_nll::set_additional_term(const boost::shared_ptr<theta::Function> & term){
    additional_term = term;
    par_ids = model.get_parameters();
    if(additional_term.get()){
         const ParIds & pids = additional_term->get_parameters();
         par_ids.insert(pids.begin(), pids.end());
    }
    fill_par_ids_vec();
}

void llvm_model_nll::set_override_distribution(const boost::shared_ptr<Distribution> & d){
    override_distribution = d;
}

double llvm_model_nll::operator()(const ParValues & values) const{
    //0. convert values to vector:
    get_par_values(parameter_values, values, par_ids_vec);
    double result = 0.0;
    //1. the model prior first, because if we are out of bounds, we should not evaluate
    //   the likelihood of the templates ...
    if(override_distribution){
        result += override_distribution->eval_nl(values);
    }
    else{
        result += model.get_parameter_distribution().eval_nl(values);
    }
    //2. get the prediction of the model:
    pred_concatenated.set_all_values(0.0);
    model_get_prediction(&parameter_values[0], pred_concatenated.get_data());
    //3. the template likelihood
    result += template_nllikelihood(data_concatenated.get_data(), pred_concatenated.get_data(), data_concatenated.size());
    //4. the likelihood part for the real-valued observables, if set:
    const Distribution * rvobs_dist = model.get_rvobservable_distribution();
    if(rvobs_dist){
        ParValues all_values(values);
        all_values.set(rvobs_values);
        result += rvobs_dist->eval_nl(all_values);
    }
    //5. The additional likelihood terms, if set:
    if(additional_term){
       result += (*additional_term)(values);
    }
    return result;
}


REGISTER_PLUGIN(llvm_model)

