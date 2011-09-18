#include "plugins/minimizer_chain.hpp"
#include "interface/plugin.hpp"

using namespace theta;
using namespace theta::plugin;
using namespace std;


MinimizationResult minimizer_chain::minimize(const Function & f, const ParValues & start,
                const ParValues & step, const map<ParId, pair<double, double> > & ranges){
    MinimizationResult res;
    bool success = false;
    for(size_t i=0; i < minimizers.size(); ++i){
        try{
            res = minimizers[i].minimize(f, start, step, ranges);
            success = true;
        }
        catch(MinimizationException & ex){
            // is this was the last attempt: re-throw, otherwise silently ignore.
            if(i+1==minimizers.size()) throw;
        }
        if(success) break;
    }
    if(last_minimizer.get()){
        ParValues step2 = step;
        // set step2 to the errors from the minimization, if available:
        const ParIds & pids = res.errors_plus.getParameters();
        for(ParIds::const_iterator it=pids.begin(); it!=pids.end(); ++it){
            double width = res.errors_plus.get(*it);
            if(width > 0)
               step2.set(*it, width);
        }
        res = last_minimizer->minimize(f, res.values, step2, ranges);
    }
    return res;
}

std::auto_ptr<theta::Minimizer> minimizer_chain::clone(const PropertyMap & pm) const{
    return std::auto_ptr<theta::Minimizer>(new minimizer_chain(*this, pm));
}

minimizer_chain::minimizer_chain(const minimizer_chain & rhs, const PropertyMap & pm){
    minimizers.reserve(rhs.minimizers.size());
    for(size_t i=0; i<rhs.minimizers.size(); ++i){
        minimizers.push_back(rhs.minimizers[i].clone(pm));
    }
    if(rhs.last_minimizer.get())
       last_minimizer = rhs.last_minimizer->clone(pm);
}

minimizer_chain::minimizer_chain(const theta::plugin::Configuration & cfg){
    SettingWrapper s = cfg.setting;
    const size_t n = s["minimizers"].size();
    minimizers.reserve(n);
    size_t n_minimizers = 0;
    for(size_t i=0; i<n; ++i){
        minimizers.push_back(PluginManager<Minimizer>::instance().build(plugin::Configuration(cfg, s["minimizers"][i])));
        ++n_minimizers;
    }
    if(s.exists("last_minimizer")){
        last_minimizer = PluginManager<Minimizer>::instance().build(plugin::Configuration(cfg, s["last_minimizer"]));
        ++n_minimizers;
    }
    if(n_minimizers==0) throw ConfigurationException("no minimizers specified; required is at least one (counting last_minimizer)");
}

REGISTER_PLUGIN(minimizer_chain)

