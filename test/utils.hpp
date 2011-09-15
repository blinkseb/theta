#ifndef TEST_UTILS_HPP
#define TEST_UTILS_HPP

#include "interface/plugin.hpp"
#include "interface/cfg-utils.hpp"
#include "libconfig/libconfig.h++"

void load_core_plugins();

//returns false if loading root plugins was not successful
bool load_root_plugins();

class ConfigCreator{
public:
    
    ConfigCreator(const std::string & cfg_string, const boost::shared_ptr<theta::VarIdManager> & vm);
    
    const theta::plugin::Configuration & get(){
        return cfg;
    }
    
private:
    int setup_config(const std::string &);
    
    libconfig::Config config;
    int dummy;
    
    boost::shared_ptr<theta::SettingUsageRecorder> rec;
    theta::plugin::Configuration cfg;
};

#endif
