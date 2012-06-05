#include "plugins/deltanll_hypotest.hpp"
#include "plugins/asimov_likelihood_widths.hpp"
#include "interface/plugin.hpp"
#include "interface/minimizer.hpp"
#include "interface/histogram.hpp"
#include "interface/distribution.hpp"
#include "interface/variables-utils.hpp"
#include "interface/model.hpp"

#include <sstream>

using namespace theta;

namespace{

bool in_support(const ParValues & values, const std::map<theta::ParId, std::pair<double, double> > & support){
    for(std::map<theta::ParId, std::pair<double, double> >::const_iterator it=support.begin(); it!=support.end(); ++it){
        double val = values.get(it->first);
        if(val < it->second.first || val > it->second.second) return false;
    }
    return true;
}

}

void deltanll_hypotest::produce(const theta::Data & data, const theta::Model & model){
    if(not init){
        ParIds model_pars = model.get_parameters();
        if(not (s_plus_b->get_parameters() == model_pars) or not (b_only->get_parameters() == model_pars)){
            throw std::invalid_argument("parameters in s+b / b only distributions do not coincide with model parameters");
        }
        s_plus_b_width.set(asimov_likelihood_widths(model, s_plus_b, additional_nll_term));
        b_only_width.set(asimov_likelihood_widths(model, b_only, additional_nll_term));
        init = true;
    }
    std::auto_ptr<NLLikelihood> nll = get_nllikelihood(data, model);

    double nll_sb, nll_b;

    // For more robustness, try to fit the more restrictive model first, then the other one, starting at the same parameter values, if possible.
    // This is hard to tell in general, but usually, the b-only model is more restrictive.
    // In case of restrict_poi however, s+b is more restrictive, so switch fitting in that case.
    // It's

    if(restrict_poi){
        theta_assert(!std::isnan(poi_value));
        //s+b first:
        nll->set_override_distribution(s_plus_b);
        s_plus_b_support[*restrict_poi].first =  s_plus_b_support[*restrict_poi].second = poi_value;
        s_plus_b_mode.set(*restrict_poi, poi_value);
        MinimizationResult minres = minimizer->minimize(*nll, s_plus_b_mode, s_plus_b_width, s_plus_b_support);
        nll_sb = minres.fval;

        //now b-only:
        nll->set_override_distribution(b_only);
        b_only_support[*restrict_poi].second = poi_value;
        if(in_support(minres.values, b_only_support)){
            minres = minimizer->minimize(*nll, minres.values, b_only_width, b_only_support);
        }
        else{
            minres = minimizer->minimize(*nll, b_only_mode, b_only_width, b_only_support);
        }
        nll_b = minres.fval;
        products_sink->set_product(c_poi, poi_value);
    }
    else{
        // b-only first:
        nll->set_override_distribution(b_only);
        MinimizationResult minres = minimizer->minimize(*nll, b_only_mode, b_only_width, b_only_support);
        nll_b = minres.fval;
        // now s+b:
        nll->set_override_distribution(s_plus_b);
        if(in_support(minres.values, s_plus_b_support)){
            // start at the b-only minimum. This is more robust than starting at the original point again.
            minres = minimizer->minimize(*nll, minres.values, s_plus_b_width, s_plus_b_support);
        }
        else{
            minres = minimizer->minimize(*nll, s_plus_b_mode, s_plus_b_width, s_plus_b_support);
        }
        nll_sb = minres.fval;
    }
    
    products_sink->set_product(c_nll_sb, nll_sb);
    products_sink->set_product(c_nll_b, nll_b);
    products_sink->set_product(c_nll_diff, nll_b - nll_sb);
}

void deltanll_hypotest::set_parameter_values(const theta::ParValues & values){
    if(restrict_poi){
        poi_value = values.get(*restrict_poi);
    }
}


deltanll_hypotest::deltanll_hypotest(const theta::Configuration & cfg):
        ParameterDependentProducer(cfg), init(false) {
    Setting s = cfg.setting;
    minimizer = theta::PluginManager<Minimizer>::build(theta::Configuration(cfg, s["minimizer"]));
    s_plus_b = theta::PluginManager<Distribution>::build(theta::Configuration(cfg, s["signal-plus-background-distribution"]));
    b_only = theta::PluginManager<Distribution>::build(theta::Configuration(cfg, s["background-only-distribution"]));
    if(s.exists("restrict_poi")){
        boost::shared_ptr<VarIdManager> vm = cfg.pm->get<VarIdManager>();
        restrict_poi = vm->get_par_id(s["restrict_poi"]);
        par_ids.insert(*restrict_poi);
        default_poi_value = NAN;
        if(s.exists("default_poi_value")) default_poi_value = s["default_poi_value"];
        poi_value = default_poi_value;
    }
    fill_mode_support(s_plus_b_mode, s_plus_b_support, *s_plus_b);
    fill_mode_support(b_only_mode, b_only_support, *b_only);
    if(not (s_plus_b->get_parameters()==b_only->get_parameters())){
        throw ConfigurationException("parameters of the distributions 'signal-plus-background' and 'background-only' do not match");
    }
    c_nll_b = products_sink->declare_product(*this, "nll_b", theta::typeDouble);
    c_nll_sb = products_sink->declare_product(*this, "nll_sb", theta::typeDouble);
    c_nll_diff = products_sink->declare_product(*this, "nll_diff", theta::typeDouble);
    if(restrict_poi){
       c_poi = products_sink->declare_product(*this, "poi", theta::typeDouble);
    }
}

REGISTER_PLUGIN(deltanll_hypotest)

