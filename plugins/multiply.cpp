#include "plugins/multiply.hpp"
#include "interface/plugin.hpp"
#include "interface/redirect_stdio.hpp"

using namespace theta;
using namespace std;

multiply::multiply(const Configuration & cfg): literal_factor(1.0){
    boost::shared_ptr<VarIdManager> vm = cfg.pm->get<VarIdManager>();
    string type = cfg.setting["type"];
    if(type == "mult"){
        theta::out << "Warning: function plugin with type='mult' is obsolete. Use type='multiply' instead and adapt configuration accordingly (see documentation; in particular use 'factors' setting instead of 'parameters')." << endl;
        //compatibility mode: search for "parameters", instead of "factors"
        size_t n = cfg.setting["parameters"].size();
        for(size_t i=0; i<n; ++i){
            ParId pid = vm->get_par_id(cfg.setting["parameters"][i]);
            v_pids.push_back(pid);
            par_ids.insert(pid);
        }
        return;
    }
    size_t n = cfg.setting["factors"].size();
    for(size_t i=0; i<n; ++i){
        Setting::Type t = cfg.setting["factors"][i].get_type();
        if(t==Setting::TypeFloat){
            literal_factor *= static_cast<double>(cfg.setting["factors"][i]);
        }
        else if(t==Setting::TypeString){
           ParId pid = vm->get_par_id(cfg.setting["factors"][i]);
           v_pids.push_back(pid);
           par_ids.insert(pid);
        }
        else if(t==Setting::TypeGroup){
            std::auto_ptr<Function> f = PluginManager<Function>::build(Configuration(cfg, cfg.setting["factors"][i]));
            par_ids.insert_all(f->get_parameters());
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
    std::vector<theta::ParId>::const_iterator it_end = v_pids.end();
    for(std::vector<theta::ParId>::const_iterator it = v_pids.begin(); it!= it_end; ++it){
        result *= v.get_unchecked(*it);
    }
    boost::ptr_vector<theta::Function>::const_iterator fit_end = functions.end();
    for(boost::ptr_vector<theta::Function>::const_iterator fit = functions.begin(); fit != fit_end; ++fit){
        result *= (*fit)(v);
    }
    return result;
}


REGISTER_PLUGIN(multiply)
REGISTER_PLUGIN_NAME(multiply,mult)
