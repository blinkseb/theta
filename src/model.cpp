#include "interface/model.hpp"
#include "interface/histogram-function.hpp"
#include "interface/distribution.hpp"
#include "interface/plugin.tcc"
#include "interface/log2_dot.hpp"

#include <limits>

using namespace std;
using namespace theta;

REGISTER_PLUGIN_BASETYPE(Model);

const ParIds & Model::getParameters() const{
    return parameters;
}

const ParIds & Model::getRVObservables() const{
    return rvobservables;
}

const ObsIds & Model::getObservables() const{
    return observables;
}


/* default_model */
void default_model::set_prediction(const ObsId & obs_id, boost::ptr_vector<Function> & coeffs_, boost::ptr_vector<HistogramFunction> & histos_){
    observables.insert(obs_id);
    const size_t n = coeffs_.size();
    if(n!=coeffs_.size()) throw invalid_argument("Model::setPrediction: number of histograms and coefficients do not match");
    if(histos[obs_id].size()>0 || coeffs[obs_id].size()>0)
        throw invalid_argument("Model::setPrediction: prediction already set for this observable");
    coeffs[obs_id].transfer(coeffs[obs_id].end(), coeffs_.begin(), coeffs_.end(), coeffs_);
    histos[obs_id].transfer(histos[obs_id].end(), histos_.begin(), histos_.end(), histos_);
    for(boost::ptr_vector<Function>::const_iterator it=coeffs[obs_id].begin(); it!=coeffs[obs_id].end(); ++it){
        ParIds pids = (*it).getParameters();
        parameters.insert(pids.begin(), pids.end());
    }
    Histogram1DWithUncertainties h;
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

void default_model::get_prediction(DataWithUncertainties & result, const ParValues & parameters) const {
    /*if(!(parameters.contains_all(this->parameters))){
        throw invalid_argument("default_model::get_prediction: not all parameters set!");
    }*/
    histos_type::const_iterator h_it = histos.begin();
    coeffs_type::const_iterator c_it = coeffs.begin();
    for(; h_it != histos.end(); ++h_it, ++c_it){
        const ObsId & oid = h_it->first;
        //theta_assert(c_it->first == oid);
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
        //theta_assert(result_init);
    }
}

std::auto_ptr<NLLikelihood> default_model::getNLLikelihood(const Data & data) const{
    if(not(data.getObservables() == observables)){
        throw invalid_argument("Model::getNLLikelihood: observables of model and data mismatch!");
    }
    if(not(data.getRVObsValues().getParameters() == rvobservables)){
        throw invalid_argument("Model::getNLLikelihood: real-values observables of model and data mismatch!");
    }
    if(bb_uncertainties){
        return std::auto_ptr<NLLikelihood>(new default_model_bbadd_nll(*this, data, observables));
    }
    else{
        return std::auto_ptr<NLLikelihood>(new default_model_nll(*this, data, observables));
    }
    throw logic_error("logic_error in default_model::getNLLikelihood");
}

default_model::default_model(const Configuration & ctx): bb_uncertainties(false) {
    SettingWrapper s = ctx.setting;
    boost::shared_ptr<VarIdManager> vm = ctx.pm->get<VarIdManager>();
    if(s.exists("bb_uncertainties")){
        bb_uncertainties =  s["bb_uncertainties"];
    }
    //go through observables to find the template definition for each of them:
    ObsIds observables = vm->getAllObsIds();
    for (ObsIds::const_iterator obsit = observables.begin(); obsit != observables.end(); obsit++) {
        string obs_name = vm->getName(*obsit);
        if(not s.exists(obs_name)) continue;
        SettingWrapper obs_setting = s[obs_name];
        boost::ptr_vector<HistogramFunction> histos;
        boost::ptr_vector<Function> coeffs;
        for (size_t i = 0; i < obs_setting.size(); i++) {
            Configuration context(ctx, obs_setting[i]["histogram"]);
            auto_ptr<HistogramFunction> hf = PluginManager<HistogramFunction>::build(context);
            Configuration cfg(ctx, obs_setting[i]["coefficient-function"]);
            auto_ptr<Function> coeff_function = PluginManager<Function>::build(cfg);
            coeffs.push_back(coeff_function);
            theta_assert(coeff_function.get()==0);
            histos.push_back(hf);
        }
        set_prediction(*obsit, coeffs, histos);
    }
    if(ctx.setting.exists("rvobs-distribution")){
        rvobservable_distribution = PluginManager<Distribution>::build(Configuration(ctx, ctx.setting["rvobs-distribution"]));
        rvobservables = rvobservable_distribution->getParameters();
        // add parameters:
        const ParIds & dist_pars = rvobservable_distribution->getDistributionParameters();
        parameters.insert(dist_pars.begin(), dist_pars.end());
    }
    // type checking for rvobs ParIds vs. parameter ParIds:
    for(ParIds::const_iterator it=parameters.begin(); it!=parameters.end(); ++it){
        if(vm->getType(*it) != "par"){
            throw ConfigurationException("Type error: parameter '" + vm->getName(*it) + "' is used as model parameter, but was not declared as such.");
        }
    }
    for(ParIds::const_iterator it=rvobservables.begin(); it!=rvobservables.end(); ++it){
        if(vm->getType(*it) != "rvobs"){
            throw ConfigurationException("Type error: parameter '" + vm->getName(*it) + "' is used as real-valued observable, but was not declared as such.");
        }
    }
    
    // parameter distribution:
    if(parameters.size() == 0){
        parameter_distribution.reset(new EmptyDistribution());
    }
    parameter_distribution = PluginManager<Distribution>::build(Configuration(ctx, s["parameter-distribution"]));
    if(not (parameter_distribution->getParameters() == parameters)){
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

default_model::~default_model(){
}

/* default_model_nll */
default_model_nll::default_model_nll(const default_model & m, const Data & dat, const ObsIds & obs): model(m),
        data(dat), obs_ids(obs){
    par_ids = model.getParameters();
}

void default_model_nll::set_additional_term(const boost::shared_ptr<Function> & term){
    additional_term = term;
    par_ids = model.getParameters();
    if(additional_term.get()){
         const ParIds & pids = additional_term->getParameters();
         par_ids.insert(pids.begin(), pids.end());
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
        result += override_distribution->evalNL(values);
    }
    else{
        result += model.get_parameter_distribution().evalNL(values);
    }
    //2. get the prediction of the model:
    model.get_prediction(predictions, values);
    //3. the template likelihood    
    for(ObsIds::const_iterator obsit=obs_ids.begin(); obsit!=obs_ids.end(); obsit++){
        const double * pred_data = predictions[*obsit].get_values().getData();
        const double * data_data = data[*obsit].getData();
        result += template_nllikelihood(data_data, pred_data, data[*obsit].get_nbins());
    }
    //4. the likelihood part for the real-valued observables, if set:
    const Distribution * rvobs_dist = model.get_rvobservable_distribution();
    if(rvobs_dist){
        ParValues all_values(values);
        all_values.set(data.getRVObsValues());
        result += rvobs_dist->evalNL(all_values);
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
        result += override_distribution->evalNL(values);
    }
    else{
        result += model.get_parameter_distribution().evalNL(values);
    }
    //2. get the prediction of the model:
    model.get_prediction(predictions, values);
    //3. the template likelihood. This is the only thing different w.r.t. the "non-bb" version ...
    for(ObsIds::const_iterator obsit=obs_ids.begin(); obsit!=obs_ids.end(); obsit++){
        const Histogram1DWithUncertainties & pred_obs = predictions[*obsit];
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
        all_values.set(data.getRVObsValues());
        result += rvobs_dist->evalNL(all_values);
    }
    //5. The additional likelihood terms, if set:
    if(additional_term){
        result += (*additional_term)(values);
    }
    return result;
}


REGISTER_PLUGIN_DEFAULT(default_model)
