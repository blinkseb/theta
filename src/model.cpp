#include "interface/model.hpp"
#include "interface/log2_dot.hpp"
#include <fstream>

using namespace std;
using namespace theta;
using namespace theta::utils;

REGISTER_PLUGIN_BASETYPE(Model);

Model::Model(const Model & model, const boost::shared_ptr<VarIdManager> & vm_):
                   vm(vm_), parameters(model.parameters), observables(model.observables){
}

ParIds Model::getParameters() const{
    return parameters;
}

ObsIds Model::getObservables() const{
    return observables;
}


/* default_model */
void default_model::set_prediction(const ObsId & obs_id, boost::ptr_vector<Function> & coeffs_, boost::ptr_vector<HistogramFunction> & histos_){
    observables.insert(obs_id);
    const size_t n = coeffs_.size();
    if(n!=coeffs_.size()) throw InvalidArgumentException("Model::setPrediction: number of histograms and coefficients do not match");
    if(histos[obs_id].size()>0)
        throw InvalidArgumentException("Model::setPrediction: prediction already set for this observable");
    histos[obs_id].transfer(histos[obs_id].end(), histos_.begin(), histos_.end(), histos_);
    if(coeffs[obs_id].size()>0)
        throw InvalidArgumentException("Model::setPrediction: prediction already set for this observable (coeffs)");
    coeffs[obs_id].transfer(coeffs[obs_id].end(), coeffs_.begin(), coeffs_.end(), coeffs_);
    for(boost::ptr_vector<Function>::const_iterator it=coeffs[obs_id].begin(); it!=coeffs[obs_id].end(); ++it){
        ParIds pids = (*it).getParameters();
        parameters.insert(pids.begin(), pids.end());
    }
    Histogram1D h;
    for(boost::ptr_vector<HistogramFunction>::const_iterator it=histos[obs_id].begin(); it!=histos[obs_id].end(); ++it){
        if(h.get_nbins()==0){
            h = it->get_histogram_dimensions();
        }
        else{
            h.check_compatibility(it->get_histogram_dimensions());
        }
        ParIds pids = (*it).getParameters();
        parameters.insert(pids.begin(), pids.end());
    }
}

void default_model::get_prediction(Data & result, const ParValues & parameters) const {
    if(!(parameters.getParameters()==this->parameters)){
        throw InvalidArgumentException("deault_model::get_prediction: parameters was incomplete");
    }
    histos_type::const_iterator h_it = histos.begin();
    coeffs_type::const_iterator c_it = coeffs.begin();
    for(; h_it != histos.end(); ++h_it, ++c_it){
        const ObsId & oid = h_it->first;
        theta_assert(c_it->first == oid);
        histos_type::const_mapped_reference h_producers = *(h_it->second);
        coeffs_type::const_mapped_reference h_coeffs = *(c_it->second);
        bool result_init = false;
        for (size_t i = 0; i < h_producers.size(); i++) {
            if(result_init){
                result[oid].add_with_coeff(h_coeffs[i](parameters), h_producers[i](parameters));
            }
            else{
                result[oid] = h_producers[i](parameters);
                result[oid] *= h_coeffs[i](parameters);
                result_init = true;
            }
        }
        //return empty prediction as correct zero template:
        if(not result_init){
            const pair<double, double> & o_range = vm->get_range(oid);
            result[oid].reset_n(vm->get_nbins(oid));
            result[oid].reset_range(o_range.first, o_range.second);
        }
    }
}

void default_model::codegen(std::ostream & out, const std::string & prefix, const PropertyMap & pm) const{
    // generate code for all Functions and HistorgamFunctions:
    histos_type::const_iterator h_it = histos.begin();
    coeffs_type::const_iterator c_it = coeffs.begin();
    size_t iobs = 0;
    for(; h_it != histos.end(); ++h_it, ++c_it, ++iobs){
        histos_type::const_mapped_reference h_functions = *(h_it->second);
        coeffs_type::const_mapped_reference h_coeffs = *(c_it->second);
        for(size_t i=0; i<h_functions.size(); ++i){
            stringstream ss_prefix;
            ss_prefix << prefix << "_hf_" << iobs << "_" << i;
            h_coeffs[i].codegen(out, ss_prefix.str(), pm);
            h_functions[i].codegen(out, ss_prefix.str(), pm);
        }
    }
    // now generate the evaluate function itself:
    h_it = histos.begin();
    c_it = coeffs.begin();
    out << "void " << prefix << "_get_prediction(const double * par_values, double * pred_data){" << endl;
    iobs = 0;
    for(; h_it != histos.end(); ++h_it, ++c_it, ++iobs){
        histos_type::const_mapped_reference h_functions = *(h_it->second);
        for(size_t i=0; i<h_functions.size(); ++i){
            stringstream ss_prefix;
            ss_prefix << prefix << "_hf_" << iobs << "_" << i;
            out << "    " << ss_prefix.str() + "_hf_add_with_coeff(" + ss_prefix.str() + "_evaluate(par_values), par_values, pred_data + obs_offset[" << iobs << "]);" << endl;
        }
    }
    out << "}" << endl;
}

void default_model::get_prediction_randomized(Random & rnd, Data & result, const ParValues & parameters) const{
    for(ObsIds::const_iterator obsit=observables.begin(); obsit!=observables.end(); ++obsit){
        histos_type::const_iterator it = histos.find(*obsit);
        theta_assert(it!=histos.end());
        histos_type::const_mapped_reference h_producers = *(it->second);
        coeffs_type::const_iterator it2 = coeffs.find(*obsit);
        coeffs_type::const_mapped_reference h_coeffs = *(it2->second);
        bool result_init = false;
        for (size_t i = 0; i < h_producers.size(); i++) {
            if(result_init){
                result[*obsit].add_with_coeff(h_coeffs[i](parameters), h_producers[i].getRandomFluctuation(rnd, parameters));
            }
            else{
                result[*obsit] = h_producers[i].getRandomFluctuation(rnd, parameters);
                result[*obsit] *= h_coeffs[i](parameters);
                result_init = true;
            }
        }
        if(not result_init){
            const pair<double, double> & o_range = vm->get_range(*obsit);
            result[*obsit].reset_n(vm->get_nbins(*obsit));
            result[*obsit].reset_range(o_range.first, o_range.second);
        }
    }
}

std::auto_ptr<NLLikelihood> default_model::getNLLikelihood(const Data & data) const{
    if(not(data.getObservables()==observables)){
        throw FatalException("Model::createNLLikelihood: observables of model and data mismatch!");
    }
    return std::auto_ptr<NLLikelihood>(new default_model_nll(*this, data, observables));
    //return getCompiledNLLikelihood(data);
}

std::auto_ptr<NLLikelihood> default_model::getCompiledNLLikelihood(const Data & data) const{
    if(not(data.getObservables()==observables)){
        throw FatalException("Model::createNLLikelihood: observables of model and data mismatch!");
    }
    if(model_get_prediction == 0){
        // generate code:
        ofstream out("compiled_model.cpp");
        PropertyMap pm;
        pm.set("default", vm);
        codegen::header(out, pm, parameters, observables);
        codegen(out, "model", pm);
        codegen::footer(out);
        out.close();
        // compile and load:
        void * handle = codegen::compile_load_so("compiled_model.cpp");
        model_get_prediction = reinterpret_cast<codegen::t_model_get_prediction>(dlsym(handle, "model_get_prediction"));
        if(model_get_prediction==0){
            throw Exception("model prediction symbol not found");
        }
    }
    return std::auto_ptr<NLLikelihood>(new compiled_model_nll(*this, data, model_get_prediction));
}

default_model::default_model(const Configuration & ctx): Model(ctx.pm->get<VarIdManager>()), model_get_prediction(0){
    SettingWrapper s = ctx.setting;
    ObsIds observables = vm->getAllObsIds();
    //go through observables to find the template definition for each of them:
    for (ObsIds::const_iterator obsit = observables.begin(); obsit != observables.end(); obsit++) {
        string obs_name = vm->getName(*obsit);
        if(not s.exists(obs_name)) continue;
        SettingWrapper obs_setting = s[obs_name];
        boost::ptr_vector<HistogramFunction> histos;
        boost::ptr_vector<Function> coeffs;
        for (size_t i = 0; i < obs_setting.size(); i++) {
            Configuration context(ctx, obs_setting[i]["histogram"]);
            auto_ptr<HistogramFunction> hf = PluginManager<HistogramFunction>::instance().build(context);
            Configuration cfg(ctx, obs_setting[i]["coefficient-function"]);
            auto_ptr<Function> coeff_function = PluginManager<Function>::instance().build(cfg);
            coeffs.push_back(coeff_function);
            theta_assert(coeff_function.get()==0);
            histos.push_back(hf);
        }
        set_prediction(*obsit, coeffs, histos);
    }
    parameter_distribution = PluginManager<Distribution>::instance().build(Configuration(ctx, s["parameter-distribution"]));
    if(not (parameter_distribution->getParameters() == getParameters())){
        stringstream ss;
        ss << "'parameter-distribution' has to define the same set of parameters the model depends on. However";
        ParIds dist_pars = parameter_distribution->getParameters();
        ParIds m_pars = getParameters();
        ParIds all_pars = m_pars;
        all_pars.insert(dist_pars.begin(), dist_pars.end());
        for(ParIds::const_iterator p_it=all_pars.begin(); p_it!=all_pars.end(); ++p_it){
            if(m_pars.contains(*p_it) && dist_pars.contains(*p_it)) continue;
            if(m_pars.contains(*p_it)){
               ss << ", the model depends on '"<< vm->getName(*p_it) << "' which the parameter distribution does not include";
            }
            else ss << ", the parameter distribution depends on '" << vm->getName(*p_it) << "' which the model does not depend on";
        }
        throw ConfigurationException(ss.str());
    }
}

default_model::default_model(const default_model & model, const theta::PropertyMap & pm): Model(model, pm.get<VarIdManager>()), model_get_prediction(0){
    for(histos_type::const_iterator it=model.histos.begin(); it!=model.histos.end(); ++it){
        boost::ptr_vector<HistogramFunction> & vec = histos[it->first];
        for(size_t i=0; i<it->second->size(); ++i){
            vec.push_back(it->second->at(i).clone());
        }
    }
    for(coeffs_type::const_iterator it=model.coeffs.begin(); it!=model.coeffs.end(); ++it){
        boost::ptr_vector<Function> & vec = coeffs[it->first];
        for(size_t i=0; i<it->second->size(); ++i){
            vec.push_back(it->second->at(i).clone());
        }
    }
    parameter_distribution = model.parameter_distribution->clone();
}

std::auto_ptr<Model> default_model::clone(const theta::PropertyMap & pm) const{
    return std::auto_ptr<Model>(new default_model(*this, pm));
}


/* default_model_nll */
default_model_nll::default_model_nll(const default_model & m, const Data & dat, const ObsIds & obs): model(m),
        data(dat), obs_ids(obs){
    Function::par_ids = model.getParameters();
    for(ParIds::const_iterator it=par_ids.begin(); it!=par_ids.end(); ++it){
        ranges[*it] = model.get_parameter_distribution().support(*it);
    }
}

void default_model_nll::set_additional_term(const boost::shared_ptr<Function> & term){
    additional_term = term;
}

void default_model_nll::set_override_distribution(const boost::shared_ptr<Distribution> & d){
    override_distribution = d;
    if(d){
        for(ParIds::const_iterator it=par_ids.begin(); it!=par_ids.end(); ++it){
            ranges[*it] = d->support(*it);
        }
    }else{
        for(ParIds::const_iterator it=par_ids.begin(); it!=par_ids.end(); ++it){
            ranges[*it] = model.get_parameter_distribution().support(*it);
        }
    }
}


double default_model_nll::operator()(const ParValues & values) const{
    double result = 0.0;
    //1. the model prior first, because if we are out of bounds, we should not evaluate
    //   the likelihood of the templates ...
    if(override_distribution){
        result += override_distribution->evalNL(values);
    }
    else{
        result += model.get_parameter_distribution().evalNL(values);
    }
    //2. get the prediction of the model:
    model.get_prediction(predictions, values);
    //3. the template likelihood    
    for(ObsIds::const_iterator obsit=obs_ids.begin(); obsit!=obs_ids.end(); obsit++){
        const double * pred_data = predictions[*obsit].getData();
        const double * data_data = data[*obsit].getData();
        result += template_nllikelihood(data_data, pred_data, data[*obsit].get_nbins());
    }
    //3. The additional likelihood terms, if set:
    if(additional_term){
       result += (*additional_term)(values);
    }
    return result;
}

/* compiled_model_nll */
compiled_model_nll::compiled_model_nll(const default_model & m, const Data & dat, codegen::t_model_get_prediction model_get_prediction_): model(m),
  model_get_prediction(model_get_prediction_){
    Function::par_ids = model.getParameters();
    size_t total_nbins = 0;
    ObsIds observables = dat.getObservables();
    for(ObsIds::const_iterator it=observables.begin(); it!=observables.end(); ++it){
        total_nbins += model.vm->get_nbins(*it);
        if(total_nbins % 2) ++total_nbins;
    }
    data_concatenated.reset_n(total_nbins);
    size_t offset = 0;
    for(ObsIds::const_iterator it=observables.begin(); it!=observables.end(); ++it){
        size_t nbins = model.vm->get_nbins(*it);
        memcpy(data_concatenated.getData() + offset, dat[*it].getData(), nbins * sizeof(double));
        offset += nbins;
        if(offset % 2) ++offset;
    }
    pred_concatenated.reset_n(total_nbins);
    parameter_values.resize(par_ids.size());
}

void compiled_model_nll::set_additional_term(const boost::shared_ptr<Function> & term){
    additional_term = term;
}

void compiled_model_nll::set_override_distribution(const boost::shared_ptr<Distribution> & d){
    override_distribution = d;
}

double compiled_model_nll::operator()(const ParValues & values) const{
    //0. convert values to vector:
    size_t i=0;
    for(ParIds::const_iterator it=par_ids.begin(); it!=par_ids.end(); ++it, ++i){
        parameter_values[i] = values.get_unchecked(*it);
    }
    double result = 0.0;
    //1. the model prior first, because if we are out of bounds, we should not evaluate
    //   the likelihood of the templates ...
    if(override_distribution){
        result += override_distribution->evalNL(values);
    }
    else{
        result += model.get_parameter_distribution().evalNL(values);
    }
    //2. get the prediction of the model:
    pred_concatenated.set_all_values(0.0);
    model_get_prediction(&parameter_values[0], pred_concatenated.getData());
    //3. the template likelihood
    result += template_nllikelihood(data_concatenated.getData(), pred_concatenated.getData(), data_concatenated.size());
    //3. The additional likelihood terms, if set:
    if(additional_term){
       result += (*additional_term)(values);
    }
    return result;
}


REGISTER_PLUGIN_DEFAULT(default_model)

