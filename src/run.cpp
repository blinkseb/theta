#include "interface/run.hpp"
#include "interface/histogram.hpp"
#include "interface/phys.hpp"
#include "interface/model.hpp"
#include "interface/redirect_stdio.hpp"

#include <iomanip>


using namespace theta;
using namespace std;


void Run::run(){
    //log the start of the run:
    //use eventid = 0 to indicate a "run-scoped" entry
    logtable->append(runid, 0, LogTable::info, "run start");
   
    Data data;
    //main event loop:
    for (int eventid = 1; eventid <= n_event; eventid++) {
        if(stop_execution)break;
        try{
            data_source->fill(data);
        }
        catch(DataSource::DataUnavailable &){
            break;
        }
        catch(theta::Exception & ex){
           ex.message += " (in Run::run while throwing toy data)";
           throw;
        }
        logtable->append(runid, eventid, LogTable::info, "start");
        bool error = false;
        for (size_t j = 0; j < producers.size(); j++) {
            try {
                producers[j].produce(data, *model);
            } catch (Exception & ex) {
                error = true;
                std::stringstream ss;
                ss << "Producer '" << producers[j].getName() << "' failed: " << ex.message << ".";
                logtable->append(runid, eventid, LogTable::error, ss.str());
                break;
            }
            catch(FatalException & f){
                stringstream ss;
                ss << "Producer '" << producers[j].getName() << "': " << f.message;
                f.message = ss.str();
                throw;
            }
        }
        //only add a row if no error ocurred to prevent NULL values and similar things ...
        if(!error){
            products_table->add_row(runid, eventid);
        }
        logtable->append(runid, eventid, LogTable::info, "end");
        if(progress_listener) progress_listener->progress(eventid, n_event);
    }
    
    logtable->append(runid, 0, LogTable::info, "run end");
    if(log_report){
        const int* n_messages = logtable->get_n_messages();
        LogTable::e_severity s = logtable->get_loglevel();
        theta::cout << endl << endl << "Log report:" << endl;
        theta::cout << "  errors:   " << setw(6) << n_messages[0] << endl;
        if(s > 0)
            theta::cout << "  warnings: " << setw(6) << n_messages[1] << endl;
        if(s > 1)
            theta::cout << "  infos:    " << setw(6) << n_messages[2] << endl;
        if(s > 2)
            theta::cout << "  debug:    " << setw(6) << n_messages[3] << endl;
    }
}

Run::Run(const plugin::Configuration & cfg){
    SettingWrapper s = cfg.setting;
    
    //0. set default values for members:
    vm = cfg.vm;
    log_report = true;
    runid = 1;
    n_event = s["n-events"];
    
    //1. setup database and tables:
    db = plugin::PluginManager<Database>::instance().build(plugin::Configuration(cfg, s["output_database"]));

    std::auto_ptr<Table> logtable_underlying = db->create_table("log");
    logtable.reset(new LogTable(logtable_underlying));
    
    std::auto_ptr<Table> rndinfo_table_underlying = db->create_table("rndinfo");
    rndinfo_table.reset(new RndInfoTable(rndinfo_table_underlying));
    cfg.pm->set("default", rndinfo_table);
    
    std::auto_ptr<Table> products_table_underlying = db->create_table("products");
    products_table.reset(new ProductsTable(products_table_underlying));
    cfg.pm->set<ProductsSink>("default", products_table);
    
    boost::shared_ptr<int> ptr_runid(new int(runid));
    cfg.pm->set("runid", ptr_runid);
        
    //2. model and data_source
    model = plugin::PluginManager<Model>::instance().build(plugin::Configuration(cfg, s["model"]));
    data_source = plugin::PluginManager<DataSource>::instance().build(plugin::Configuration(cfg, s["data_source"]));
    
    //3. logging stuff
    LogTable::e_severity level = LogTable::warning;
    if(s.exists("log-level")){
        std::string loglevel = s["log-level"];
        if(loglevel=="error") level = LogTable::error;
        else if(loglevel=="warning") level = LogTable::warning;
        else if(loglevel=="info")level = LogTable::info;
        else if(loglevel=="debug")level = LogTable::debug;
        else{
            std::stringstream ss;
            ss << "log level given in " << s["log-level"].getPath() << " unknown (given '" << loglevel << 
                  "'; only allowed values are 'error', 'warning', 'info' and 'debug')";
            throw ConfigurationException(ss.str());
        }
    }
    logtable->set_loglevel(level);
    if(s.exists("log-report")){
        log_report = s["log-report"];
    }
    
    //4. producers:
    size_t n_p = s["producers"].size();
    if (n_p == 0)
        throw ConfigurationException("no producers specified!");
    for (size_t i = 0; i < n_p; i++) {
         producers.push_back(plugin::PluginManager<Producer>::instance().build(plugin::Configuration(cfg, s["producers"][i])));
    }
}

REGISTER_PLUGIN_DEFAULT(Run)

