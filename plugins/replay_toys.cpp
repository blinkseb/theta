#include "plugins/replay_toys.hpp"
#include "interface/data.hpp"
#include "interface/histogram.hpp"

using namespace theta;
using namespace std;

void replay_toys::fill(theta::Data & dat){
    if(!res->has_data()){
        throw logic_error("replay_toys: no more input data available");
    }
    for(size_t i=0; i<observables.size(); ++i){
        dat[observables[i]] = res->get_histogram(i);
    }
    ++(*res);
}

replay_toys::replay_toys(const theta::Configuration & cfg): DataSource(cfg){
    boost::shared_ptr<VarIdManager> vm = cfg.pm->get<VarIdManager>();
    if(cfg.setting.exists("observables")){
        size_t n = cfg.setting["observables"].size();
        observables.reserve(n);
        for(size_t i=0; i<n; i++){
            observables.push_back(vm->get_obs_id(cfg.setting["observables"][i]));
        }
    }
    else{
        ObsIds oids = vm->get_all_observables();
        observables.reserve(oids.size());
        for(ObsIds::const_iterator oit=oids.begin(); oit!=oids.end(); ++oit){
            observables.push_back(*oit);
        }
    }
    input_database = PluginManager<DatabaseInput>::build(Configuration(cfg, cfg.setting["input_database"]));
    
    std::vector<std::string> colnames; // column names in the same order as "observables"
    std::string pdw_name("pdw");
    if(cfg.setting.exists("pdw_name")){
        pdw_name = static_cast<string>(cfg.setting["pdw_name"]);
        for(size_t i=0; i<observables.size(); ++i){
            colnames.push_back(pdw_name + "_data_" + vm->get_name(observables[i]));
        }
    }else{
        std::vector<std::pair<std::string, data_type> > columns = input_database->get_all_columns("products");
        for(size_t i=0; i<observables.size(); ++i){
            string colend = "_data_" + vm->get_name(observables[i]);
            size_t n_found = 0;
            for(size_t icol=0; icol<columns.size(); ++icol){
                size_t p = columns[icol].first.find(colend);
                if(p!=string::npos && p==columns[icol].first.size() - colend.size()){
                    ++n_found;
                    if(n_found>1){
                        throw ConfigurationException("'pdw_name' not given, so guessed columns to use from which to take data. However, there are multiple columns ending with'" + colend+"'");
                    }
                    colnames.push_back(columns[icol].first);
                }
            }
        }
    }
    theta_assert(colnames.size() == observables.size());
    res = input_database->query("products", colnames);
    if(!res->has_data()){
        throw ConfigurationException("no input data available");
    }
}

REGISTER_PLUGIN(replay_toys)
