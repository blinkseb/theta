#include "plugins/nll_scan.hpp"
#include "plugins/asimov_likelihood_widths.hpp"
#include "plugins/reduced_nll.hpp"
#include "interface/plugin.hpp"
#include "interface/histogram.hpp"
#include "interface/distribution.hpp"

#include <sstream>

using namespace theta;
using namespace std;
using namespace libconfig;

void nll_scan::produce(const Data & data, const Model & model) {
    std::auto_ptr<NLLikelihood> nll = get_nllikelihood(data, model);
    if(not start_step_ranges_init){
        const Distribution & d = nll->get_parameter_distribution();
        fill_mode_support(m_start, m_ranges, d);
        m_step.set(asimov_likelihood_widths(model, override_parameter_distribution, additional_nll_term));
        start_step_ranges_init = true;
    }
    MinimizationResult minres = minimizer->minimize(*nll, m_start, m_step, m_ranges);
    products_sink->set_product(c_maxl, minres.values.get(pid));
    ReducedNLL nll_r(*nll, pid, minres.values, re_minimize ? minimizer.get() : 0, m_start, m_step, m_ranges);
    nll_r.set_offset_nll(minres.fval);
    
    theta::Histogram1D result(n_steps, start, start + n_steps * step);
    for(unsigned int i=0; i<n_steps; ++i){
        double x = start + i * step;
        result.set(i, nll_r(x));
    }
    products_sink->set_product(c_nll, result);
}


nll_scan::nll_scan(const theta::Configuration & cfg): Producer(cfg), pid(cfg.pm->get<VarIdManager>()->get_par_id(cfg.setting["parameter"])),
   re_minimize(true), start_step_ranges_init(false){
    SettingWrapper s = cfg.setting;
    minimizer = PluginManager<Minimizer>::build(theta::Configuration(cfg, s["minimizer"]));
    if(s.exists("re-minimize")){
        re_minimize = s["re-minimize"];
    }
    start = s["parameter-values"]["start"];
    stop = s["parameter-values"]["stop"];
    n_steps = s["parameter-values"]["n-steps"];
    if(n_steps<2){
        throw ConfigurationException("nll_scan: n-steps must be >= 2");
    }
    if(start >= stop){
        throw ConfigurationException("nll_scan: start < stop must hold");
    }
    step = (stop - start) / n_steps;
    c_nll = products_sink->declare_product(*this, "nll", theta::typeHisto);
    c_maxl = products_sink->declare_product(*this, "maxl", theta::typeDouble);
}

REGISTER_PLUGIN(nll_scan)
