#include "plugins/multiply.hpp"

using namespace libconfig;
using namespace theta;
using namespace theta::plugin;

multiply::multiply(const Configuration & cfg): literal_factor(1.0){
    size_t n = cfg.setting["factors"].size();
    for(size_t i=0; i<n; ++i){
        Setting::Type t = cfg.setting["factors"][i].getType();
        if(t==Setting::TypeFloat){
            literal_factor *= static_cast<double>(cfg.setting["factors"][i]);
        }
        else if(t==Setting::TypeString){
           ParId pid = cfg.vm->getParId(cfg.setting["factors"][i]);
           v_pids.push_back(pid);
           par_ids.insert(pid);
        }
        else if(t==Setting::TypeGroup){
            std::auto_ptr<Function> f = PluginManager<Function>::instance().build(Configuration(cfg, cfg.setting["factors"][i]));
            const ParIds & f_p = f->getParameters();
            par_ids.insert(f_p.begin(), f_p.end());
            functions.push_back(f);
        }
        else{
           std::stringstream ss;
           ss << "Invalid config setting type in setting 'factors' at index " << i;
           throw ConfigurationException(ss.str());
        }
    }
}

double multiply::operator()(const ParValues & v) const{
    double result = literal_factor;
    for(size_t i=0; i<v_pids.size(); ++i){
        result *= v.get(v_pids[i]);
    }
    for(size_t i=0; i<functions.size(); ++i){
        result *= functions[i](v);
    }
    return result;
}


REGISTER_PLUGIN(multiply)

