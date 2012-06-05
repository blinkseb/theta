#include "interface/model.hpp"
#include "interface/histogram-function.hpp"
#include "interface/distribution.hpp"
#include "interface/plugin.tcc"
#include "interface/log2_dot.hpp"

#include <limits>

using namespace std;
using namespace theta;

REGISTER_PLUGIN_BASETYPE(Model);

const ParIds & Model::get_parameters() const{
    return parameters;
}

const ParIds & Model::get_rvobservables() const{
    return rvobservables;
}

const ObsIds & Model::get_observables() const{
    return observables;
}


/* default_model */
void default_model::set_prediction(const ObsId & obs_id, boost::ptr_vector<Function> & coeffs_, boost::ptr_vector<HistogramFunction> & histos_){
    observables.insert(obs_id);
    const size_t n = coeffs_.size();
    if(n!=coeffs_.size()) throw invalid_argument("number of histograms and coefficients do not match");
    if(histos[obs_id].size()>0 || coeffs[obs_id].size()>0){
        throw invalid_argument("prediction already set for this observable");
    }
    if(n==0) throw invalid_argument("empty prediction set");
    coeffs[obs_id].transfer(coeffs[obs_id].end(), coeffs_.begin(), coeffs_.end(), coeffs_);
    histos[obs_id].transfer(histos[obs_id].end(), histos_.begin(), histos_.end(), histos_);
    for(boost::ptr_vector<Function>::const_iterator it=coeffs[obs_id].begin(); it!=coeffs[obs_id].end(); ++it){
        parameters.insert_all(it->get_parameters());
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
                throw invalid_argument("default_model::set_prediction: histogram dimensions mismatch");
            }
        }
        parameters.insert_all(it->get_parameters());
    }
}

template<typename HT>
void default_model::get_prediction_impl(DataT<HT> & result, const ParValues & parameters) const{
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

void default_model::get_prediction(DataWithUncertainties & result, const ParValues & parameters) const {
    get_prediction_impl<Histogram1DWithUncertainties>(result, parameters);
}

void default_model::get_prediction(Data & result, const ParValues & parameters) const {
    get_prediction_impl<Histogram1D>(result, parameters);
}

std::auto_ptr<NLLikelihood> default_model::get_nllikelihood(const Data & data) const{
    if(not(data.get_observables() == observables)){
        throw invalid_argument("default_model::get_nllikelihood: observables of model and data mismatch!");
    }
    if(not(data.get_rvobs_values().contains_all(rvobservables))){
        throw invalid_argument("default_model::get_nllikelihood: real-values observables of model and data mismatch!");
    }
    if(bb_uncertainties){
        return std::auto_ptr<NLLikelihood>(new default_model_bbadd_nll(*this, data, observables));
    }
    return std::auto_ptr<NLLikelihood>(new default_model_nll(*this, data, observables));
}

default_model::default_model(const Configuration & ctx): bb_uncertainties(false) {
    Setting s = ctx.setting;
    boost::shared_ptr<VarIdManager> vm = ctx.pm->get<VarIdManager>();
    if(s.exists("bb_uncertainties")){
        bb_uncertainties =  s["bb_uncertainties"];
    }
    //go through observables to find the template definition for each of them:
    ObsIds observables = vm->get_all_observables();
    for (ObsIds::const_iterator obsit = observables.begin(); obsit != observables.end(); obsit++) {
        string obs_name = vm->get_name(*obsit);
        if(not s.exists(obs_name)) continue;
        Setting obs_setting = s[obs_name];
        if(obs_setting.size()==0) throw ConfigurationException("observable '" + vm->get_name(*obsit) + "' is empty");
        boost::ptr_vector<HistogramFunction> histos;
        boost::ptr_vector<Function> coeffs;
        for (size_t i = 0; i < obs_setting.size(); i++) {
            auto_ptr<HistogramFunction> hf = PluginManager<HistogramFunction>::build(Configuration(ctx, obs_setting[i]["histogram"]));
            auto_ptr<Function> coeff_function = PluginManager<Function>::build(Configuration(ctx, obs_setting[i]["coefficient-function"]));
            coeffs.push_back(coeff_function);
            histos.push_back(hf);
        }
        set_prediction(*obsit, coeffs, histos);
    }
    if(ctx.setting.exists("rvobs-distribution")){
        rvobservable_distribution = PluginManager<Distribution>::build(Configuration(ctx, ctx.setting["rvobs-distribution"]));
        rvobservables = rvobservable_distribution->get_parameters();
        // add parameters:
        parameters.insert_all(rvobservable_distribution->get_distribution_parameters());
    }
    // type checking for rvobs ParIds vs. parameter ParIds:
    for(ParIds::const_iterator it=parameters.begin(); it!=parameters.end(); ++it){
        if(vm->get_type(*it) != "par"){
            throw ConfigurationException("Type error: parameter '" + vm->get_name(*it) + "' is used as model parameter, but was not declared as such.");
        }
    }
    for(ParIds::const_iterator it=rvobservables.begin(); it!=rvobservables.end(); ++it){
        if(vm->get_type(*it) != "rvobs"){
            throw ConfigurationException("Type error: parameter '" + vm->get_name(*it) + "' is used as real-valued observable, but was not declared as such.");
        }
    }
    
    // parameter distribution:
    if(parameters.size() == 0){
        parameter_distribution.reset(new EmptyDistribution());
    }
    parameter_distribution = PluginManager<Distribution>::build(Configuration(ctx, s["parameter-distribution"]));
    if(not (parameter_distribution->get_parameters() == parameters)){
        stringstream ss;
        ss << "'parameter-distribution' has to define the same set of parameters the model depends on. However";
        ParIds dist_pars = parameter_distribution->get_parameters();
        ParIds all_pars = parameters;
        all_pars.insert_all(dist_pars);
        for(ParIds::const_iterator p_it=all_pars.begin(); p_it!=all_pars.end(); ++p_it){
            if(parameters.contains(*p_it) && dist_pars.contains(*p_it)) continue;
            if(parameters.contains(*p_it)){
               ss << ", the model depends on '"<< vm->get_name(*p_it) << "' which the parameter distribution does not include";
            }
            else ss << ", the parameter distribution depends on '" << vm->get_name(*p_it) << "' which the model does not depend on";
        }
        throw ConfigurationException(ss.str());
    }
}

default_model::~default_model(){
}

/* default_model_nll */
default_model_nll::default_model_nll(const default_model & m, const Data & dat, const ObsIds & obs): model(m),
        data(dat), obs_ids(obs){
    par_ids = model.get_parameters();
}

void default_model_nll::set_additional_term(const boost::shared_ptr<Function> & term){
    additional_term = term;
    par_ids = model.get_parameters();
    if(additional_term.get()){
         par_ids.insert_all(additional_term->get_parameters());
    }
}

void default_model_nll::set_override_distribution(const boost::shared_ptr<Distribution> & d){
    override_distribution = d;
}


double default_model_nll::operator()(const ParValues & values) const{
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
    model.get_prediction(predictions, values);
    //3. the template likelihood    
    for(ObsIds::const_iterator obsit=obs_ids.begin(); obsit!=obs_ids.end(); obsit++){
        const double * pred_data = predictions[*obsit].get_data();
        const double * data_data = data[*obsit].get_data();
        result += template_nllikelihood(data_data, pred_data, data[*obsit].get_nbins());
    }
    //4. the likelihood part for the real-valued observables, if set:
    const Distribution * rvobs_dist = model.get_rvobservable_distribution();
    if(rvobs_dist){
        ParValues all_values(values);
        all_values.set(data.get_rvobs_values());
        result += rvobs_dist->eval_nl(all_values);
    }
    //5. The additional likelihood terms, if set:
    if(additional_term){
       result += (*additional_term)(values);
    }
    return result;
}

// bbadd
default_model_bbadd_nll::default_model_bbadd_nll(const default_model & m, const Data & dat, const ObsIds & obs): default_model_nll(m, dat, obs){
}

double default_model_bbadd_nll::operator()(const ParValues & values) const{
    double result = 0.0;
    //1. the model prior first, because if we are out of bounds, we should not evaluate
    //   the likelihood of the templates ...
    if(override_distribution){
        result += override_distribution->eval_nl(values);
    }
    else{
        result += model.get_parameter_distribution().eval_nl(values);
    }
    //2. get the prediction of the model, with uncertainties:
    model.get_prediction(predictions_wu, values);
    //3. the template likelihood. This is the only thing different w.r.t. the "non-bb" version ...
    for(ObsIds::const_iterator obsit=obs_ids.begin(); obsit!=obs_ids.end(); obsit++){
        const Histogram1DWithUncertainties & pred_obs = predictions_wu[*obsit];
        const Histogram1D & data_obs = data[*obsit];
        const size_t nbins = data_obs.get_nbins();
        theta_assert(nbins == pred_obs.get_nbins());
        for(size_t ibin=0; ibin < nbins; ++ibin){
            const double p = pred_obs.get_value(ibin);
            const double d = data_obs.get(ibin);
            const double p_unc2 = pred_obs.get_uncertainty2(ibin);
            double beta = 0.0;
            if(p_unc2 > 0.0){
                double dummy;
                theta::utils::roots_quad(dummy, beta, p + p_unc2, p_unc2 * (p - d));
                result += 0.5 * beta * beta / p_unc2;
            }
            const double new_pred = beta + p;
            // As special case, new_pred == 0.0 can happen (for p == p_unc2 and data == 0). In this case,
            // the log-term can be skipped as it has vanishing contribution to the nll.
            result += new_pred;
            if(d > 0.0){
                if(new_pred <= 0.0){
                    return numeric_limits<double>::infinity();
                }
                result -= d * utils::log(new_pred);
            }
        }
    }
    //4. the likelihood part for the real-valued observables, if set:
    const Distribution * rvobs_dist = model.get_rvobservable_distribution();
    if(rvobs_dist){
        ParValues all_values(values);
        all_values.set(data.get_rvobs_values());
        result += rvobs_dist->eval_nl(all_values);
    }
    //5. The additional likelihood terms, if set:
    if(additional_term){
        result += (*additional_term)(values);
    }
    return result;
}


REGISTER_PLUGIN_DEFAULT(default_model)
