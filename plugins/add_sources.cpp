#include "plugins/add_sources.hpp"

using namespace theta;
using namespace theta::plugin;

add_sources::add_sources(const theta::plugin::Configuration & cfg): DataSource(cfg){
    SettingWrapper s = cfg.setting["sources"];
    size_t n = s.size();
    if(n==0) ConfigurationException("have to specify at least one source in 'sources'");
    for(size_t i=0; i<n; ++i){
        sources.push_back(PluginManager<DataSource>::instance().build(Configuration(cfg, s[i])));
    }
}

void add_sources::fill(theta::Data & dat){
    dat.reset();
    sources[0].fill(dat);
    ObsIds obs_result = dat.getObservables();
    ObsIds obs_all = obs_result;
    for(size_t i=1; i<sources.size(); ++i){
        theta::Data data_tmp;
        sources[i].fill(data_tmp);
        ObsIds obs_tmp = data_tmp.getObservables();
        obs_all.insert(obs_tmp.begin(), obs_tmp.end());
        for(ObsIds::const_iterator it=obs_all.begin(); it!=obs_all.end(); ++it){
            // add data_tmp[*it]  to  dat[*it]. However, some might be empty, so check:
            if(!obs_tmp.contains(*it)) continue;
            if(obs_result.contains(*it)){
               dat[*it] += data_tmp[*it];
            }
            else{
               dat[*it] = data_tmp[*it];
            }
        }
    }
}

REGISTER_PLUGIN(add_sources)
