#include "test/utils.hpp"

#include "interface/plugin.hpp"
#include "interface/variables.hpp"

#include <string>
#include <iostream>
#include <boost/test/unit_test.hpp>

using namespace theta;


using namespace std;


ConfigCreator::ConfigCreator(const std::string & cfg_string, const boost::shared_ptr<theta::VarIdManager> & vm):
      rec(new SettingUsageRecorder()), cfg(setup_cfg(cfg_string)){
    cfg.pm->set("default", vm);
}

Configuration ConfigCreator::setup_cfg(const string & cfg_string){
    try{
        Setting root = LibconfigSetting::parse(cfg_string, rec);
        return Configuration(root);
    }
    catch(Exception & ex){
        std::cerr << "ConfigCreator: " << ex.what() << endl;
        throw;
    }
}

void load_core_plugins(){
    static bool loaded(false);
    if(loaded) return;
    BOOST_TEST_CHECKPOINT("loading core plugin");
    try{
        PluginLoader::load("lib/core-plugins.so");
    }
    catch(exception & ex){
      std::cout << "std::exception in load_core_plugins: " << ex.what() << std::endl;
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
    catch(exception & ex){
        return false;
    }
    BOOST_TEST_CHECKPOINT("loaded root plugin");
    loaded = true;
    return true;
}

bool load_llvm_plugins(){
    static bool loaded(false);
    if(loaded) return true;
    BOOST_TEST_CHECKPOINT("loading llvm plugins");
    try{
        PluginLoader::load("lib/llvm-plugins.so");
    }
    catch(exception & ex){
        //std::cout << "error loadin llvm plugins: " << ex.message << endl;
        return false;
    }
    BOOST_TEST_CHECKPOINT("loaded llvm plugin");
    loaded = true;
    return true;
}


