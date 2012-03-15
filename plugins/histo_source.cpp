#include "plugins/histo_source.hpp"
#include "interface/histogram-function.hpp"

using namespace theta;
using namespace std;

histo_source::histo_source(const Configuration & cfg): DataSource(cfg){
    size_t n = cfg.setting.size();
    boost::shared_ptr<VarIdManager> vm = cfg.pm->get<VarIdManager>();
    for(size_t i=0; i<n; ++i){
        SettingWrapper s = cfg.setting[i];
        if(not s.exists("type")) continue;
        string obs_name = s.get_name();
        ObsId obs_id = vm->get_obs_id(obs_name);
        std::auto_ptr<HistogramFunction> hf = PluginManager<HistogramFunction>::build(Configuration(cfg, s));
        if(hf->get_parameters().size() > 0){
            throw ConfigurationException("histo_source: given histogram depends on parameters, which is not allowed");
        }
        hf->apply_functor(copy_to<Histogram1D>(data[obs_id]), ParValues());
    }
    if(cfg.setting.exists("rvobs-values")){
        ParValues rvobs_values;
        for(size_t i=0; i<cfg.setting["rvobs-values"].size(); ++i){
            string name = cfg.setting["rvobs-values"][i].get_name();
            ParId id = vm->get_par_id(name);
            // type checking:
            if(vm->get_type(id)!="rvobs"){
                throw ConfigurationException("type error: parameter '" + name + "' used as real-valued observable, but was not declared as such.");
            }
            rvobs_values.set(id, cfg.setting["rvobs-values"][i]);
        }
        data.set_rvobs_values(rvobs_values);
    }
}

void histo_source::fill(Data & dat){
    dat = data;
}

REGISTER_PLUGIN(histo_source)
