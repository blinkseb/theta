#include "interface/cfg-utils.hpp"
#include "interface/exception.hpp"
#include <sstream>
#include <limits>

using namespace std;
using namespace libconfig;
using namespace theta;

void SettingUsageRecorder::markAsUsed(const libconfig::Setting & s){
    used_paths.insert(s.getPath());
}

void SettingUsageRecorder::get_unused(std::vector<std::string> & unused, const libconfig::Setting & aggregate_setting) const{
    int n = aggregate_setting.getLength();
    for(int i=0; i<n; ++i){
        std::string path = aggregate_setting[i].getPath();
        bool a_unused = false;
        if(used_paths.find(path) == used_paths.end()){
            unused.push_back(path);
            a_unused = true;
        }
        //don't descend if already aggregate was reported as unused ...
        if(aggregate_setting[i].isAggregate() && !a_unused){
            get_unused(unused, aggregate_setting[i]);
        }
    }
}

double SettingWrapper::get_double_or_inf() const {
    rec->markAsUsed(setting);
    if(setting.getType()==Setting::TypeFloat) return setting;
    string infstring = setting;
    if(infstring == "inf" || infstring == "+inf") return numeric_limits<double>::infinity();
    if(infstring == "-inf") return -numeric_limits<double>::infinity();
    stringstream error;
    error << "error reading double (or \"inf\") from configuration path " << getPath();
    throw InvalidArgumentException(error.str());
}

const Setting & SettingWrapper::resolve_link(const Setting & setting, const Setting & root, const boost::shared_ptr<SettingUsageRecorder> & rec){
    try{
        std::string next_path;
        //hard-code maximum redirection level of 10:
        for(int i=0; i <= 10; ++i){
            const libconfig::Setting * p_s = &setting;
            
            if(i>0){
                p_s = &root;
                do{
                    size_t dotpos = next_path.find('.');
                    if(dotpos==string::npos) dotpos = next_path.size();
                    p_s = &((*p_s)[next_path.substr(0, dotpos)]);
                    next_path.erase(0, dotpos+1);
                }while(next_path.size());
            }

            const libconfig::Setting & s = *p_s;
            
            if(s.getType() != libconfig::Setting::TypeString) return s;
            std::string link = s;
            if(link.size()==0 || link[0]!='@'){
                return s;
            }
            link.erase(0, 1);
            next_path = link;
            //mark any intermediate link as used:
            rec->markAsUsed(s);
        }
    }
    catch(Exception & ex){
        std::stringstream ss;
        ss << "Exception while trying to resolve link at " << setting.getPath() << ": " << ex.message;
        ex.message = ss.str();
        throw;
    }
    std::stringstream ss;
    ss << "While trying to resolve link at " << setting.getPath() << ": link level is too deep";
    throw ConfigurationException(ss.str());
}

SettingWrapper::SettingWrapper(const libconfig::Setting & s, const libconfig::Setting & root,
                             const boost::shared_ptr<SettingUsageRecorder> & recorder):
           rootsetting(root), rec(recorder), setting(resolve_link(s, rootsetting, rec)), setting_name("<noname>"){
    const char * name = s.getName();
    if(name) setting_name = name;
}

