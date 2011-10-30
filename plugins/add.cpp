#include "plugins/add.hpp"
#include "interface/plugin.hpp"
#include "interface/redirect_stdio.hpp"

using namespace theta;
using namespace std;

add::add(const Configuration & cfg): literal_addend(0.0){
    boost::shared_ptr<VarIdManager> vm = cfg.pm->get<VarIdManager>();
    size_t n = cfg.setting["addends"].size();
    for(size_t i=0; i<n; ++i){
        libconfig::Setting::Type t = cfg.setting["addends"][i].getType();
        if(t==libconfig::Setting::TypeFloat){
            literal_addend += static_cast<double>(cfg.setting["addends"][i]);
        }
        else if(t==libconfig::Setting::TypeString){
           ParId pid = vm->getParId(cfg.setting["addends"][i]);
           v_pids.push_back(pid);
           par_ids.insert(pid);
        }
        else if(t==libconfig::Setting::TypeGroup){
            std::auto_ptr<Function> f = PluginManager<Function>::build(Configuration(cfg, cfg.setting["addends"][i]));
            const ParIds & f_p = f->getParameters();
            par_ids.insert(f_p.begin(), f_p.end());
            functions.push_back(f);
        }
        else{
           std::stringstream ss;
           ss << "Invalid config setting type in setting 'addends' at index " << i;
           throw ConfigurationException(ss.str());
        }
    }
}

double add::operator()(const ParValues & v) const{
    double result = literal_addend;
    for(size_t i=0; i<v_pids.size(); ++i){
        result += v.get_unchecked(v_pids[i]);
    }
    for(size_t i=0; i<functions.size(); ++i){
        result += functions[i](v);
    }
    return result;
}


REGISTER_PLUGIN(add)

