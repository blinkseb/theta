#include "plugins/histo_source.hpp"
#include "interface/histogram-function.hpp"

using namespace theta;
using namespace theta::plugin;
using namespace std;

histo_source::histo_source(const Configuration & cfg): DataSource(cfg){
    size_t n = cfg.setting.size();
    for(size_t i=0; i<n; ++i){
        SettingWrapper s = cfg.setting[i];
        if(not s.exists("type")) continue;
        string obs_name = s.getName();
        ObsId obs_id = cfg.vm->getObsId(obs_name);
        std::auto_ptr<HistogramFunction> hf = PluginManager<HistogramFunction>::instance().build(Configuration(cfg, s));
        if(hf->getParameters().size() > 0){
            throw ConfigurationException("histo_source: given histogram depends on parameters, which is not allowed");
        }
        data[obs_id] = (*hf)(ParValues());
    }
}

void histo_source::fill(Data & dat){
    dat = data;
}

REGISTER_PLUGIN(histo_source)
