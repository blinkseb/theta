#include "test/utils.hpp"

#include "interface/plugin.hpp"

#include <string>
#include <iostream>
#include <boost/test/unit_test.hpp>

using namespace theta;
using namespace theta::plugin;

using namespace std;


ConfigCreator::ConfigCreator(const std::string & cfg_string, const boost::shared_ptr<theta::VarIdManager> & vm):
      dummy(setup_config(cfg_string)), rec(new SettingUsageRecorder()), cfg(vm, SettingWrapper(config.getRoot(), config.getRoot(), rec)){
}

int ConfigCreator::setup_config(const std::string & cfg_string){
    try{
        config.readString(cfg_string);
    }
    catch(libconfig::ParseException & ex){
        cerr << "ConfigCreator: parse Exception: " << ex.getError() << " on line " << ex.getLine() << endl;
    }
    return 0;
}


void load_core_plugins(){
    static bool loaded(false);
    if(loaded) return;
    BOOST_TEST_CHECKPOINT("loading core plugin");
    try{
        PluginLoader::load("lib/core-plugins.so");
    }
    catch(Exception & ex){
      std::cout << ex.message << std::endl;
      throw;
    }
    BOOST_TEST_CHECKPOINT("loaded core plugin");
    loaded = true;
}

bool load_root_plugins(){
    static bool loaded(false);
    if(loaded) return true;
    BOOST_TEST_CHECKPOINT("loading root plugin");
    try{
        PluginLoader::load("lib/root.so");
    }
    catch(Exception & ex){
        return false;
    }
    BOOST_TEST_CHECKPOINT("loaded root plugin");
    loaded = true;
    return true;
}

