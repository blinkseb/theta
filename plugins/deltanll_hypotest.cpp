#include "plugins/deltanll_hypotest.hpp"
#include "plugins/asimov_likelihood_widths.hpp"
#include "interface/plugin.hpp"
#include "interface/minimizer.hpp"
#include "interface/histogram.hpp"
#include "interface/distribution.hpp"
#include "interface/variables-utils.hpp"

#include "interface/redirect_stdio.hpp" //TODO: remove

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
        ParIds model_pars = model.getParameters();
        if(not (s_plus_b->getParameters() == model_pars) or not (b_only->getParameters() == model_pars)){
            throw FatalException(Exception("parameters in s+b / b only distributions do not coincide with model parameters"));
        }
        s_plus_b_width.set(asimov_likelihood_widths(model, s_plus_b));
        b_only_width.set(asimov_likelihood_widths(model, b_only));
        init = true;
    }
    std::auto_ptr<NLLikelihood> nll = get_nllikelihood(data, model);
    // fit b-only first:
    nll->set_override_distribution(b_only);
    MinimizationResult minres = minimizer->minimize(*nll, b_only_mode, b_only_width, b_only_support);
    double nll_b = minres.fval;

    nll->set_override_distribution(s_plus_b);    
    if(in_support(minres.values, s_plus_b_support)){
        // start at the b-only minimum. This is more robust than starting at the original point again.
        minres = minimizer->minimize(*nll, minres.values, s_plus_b_width, s_plus_b_support);
    }
    else{
        minres = minimizer->minimize(*nll, s_plus_b_mode, s_plus_b_width, s_plus_b_support);
    }
    
    double nll_sb = minres.fval;
    
    products_sink->set_product(c_nll_sb, nll_sb);
    products_sink->set_product(c_nll_b, nll_b);
    products_sink->set_product(c_nll_diff, nll_b - nll_sb);
}

std::auto_ptr<theta::Producer> deltanll_hypotest::clone(const theta::PropertyMap & pm) const{
    return std::auto_ptr<Producer>(new deltanll_hypotest(*this, pm));
}

deltanll_hypotest::deltanll_hypotest(const deltanll_hypotest & rhs, const theta::PropertyMap & pm): Producer(rhs, pm), init(rhs.init),
  s_plus_b_mode(rhs.s_plus_b_mode), b_only_mode(rhs.b_only_mode), s_plus_b_width(rhs.s_plus_b_width), b_only_width(rhs.b_only_width),
  s_plus_b_support(rhs.s_plus_b_support), b_only_support(rhs.b_only_support){
    s_plus_b = rhs.s_plus_b->clone();
    b_only = rhs.b_only->clone();
    minimizer = rhs.minimizer->clone(pm);
    
    c_nll_b = products_sink->declare_product(*this, "nll_b", theta::typeDouble);
    c_nll_sb = products_sink->declare_product(*this, "nll_sb", theta::typeDouble);
    c_nll_diff = products_sink->declare_product(*this, "nll_diff", theta::typeDouble);
}


deltanll_hypotest::deltanll_hypotest(const theta::Configuration & cfg):
        Producer(cfg), init(false){
    SettingWrapper s = cfg.setting;
    minimizer = theta::PluginManager<Minimizer>::instance().build(theta::Configuration(cfg, s["minimizer"]));
    s_plus_b = theta::PluginManager<Distribution>::instance().build(theta::Configuration(cfg, s["signal-plus-background-distribution"]));
    b_only = theta::PluginManager<Distribution>::instance().build(theta::Configuration(cfg, s["background-only-distribution"]));
    DistributionUtils::fillModeSupport(s_plus_b_mode, s_plus_b_support, *s_plus_b);
    DistributionUtils::fillModeSupport(b_only_mode, b_only_support, *b_only);
    if(not (b_only_mode.getParameters()==s_plus_b_mode.getParameters())){
        throw ConfigurationException("parameters of the distributions 'signal-plus-background' and 'background-only' do not match");
    }
    c_nll_b = products_sink->declare_product(*this, "nll_b", theta::typeDouble);
    c_nll_sb = products_sink->declare_product(*this, "nll_sb", theta::typeDouble);
    c_nll_diff = products_sink->declare_product(*this, "nll_diff", theta::typeDouble);
}

REGISTER_PLUGIN(deltanll_hypotest)

