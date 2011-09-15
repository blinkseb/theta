#include "interface/variables.hpp"
#include "interface/variables-utils.hpp"
#include "interface/plugin.hpp"
#include "interface/exception.hpp"

#include <sstream>
#include <cmath>

using namespace theta;
using namespace theta::utils;
using namespace std;
using namespace libconfig;

ParId VarIdManager::createParId(const std::string & name) {
    if (parNameExists(name)) {
        stringstream ss;
        ss << "VarIdManager::createParId: parameter '"<< name <<"' defined twice";
        throw InvalidArgumentException(ss.str());
    }
    ParId result(next_pid_id);
    ++next_pid_id;
    pid_to_name[result.id] = name;
    name_to_pid[name] = result.id;
    return result;
}

ObsId VarIdManager::createObsId(const std::string & name, size_t nbins, double min, double max) {
    if (obsNameExists(name)) {
        stringstream ss;
        ss << "VarIdManager::createObsId: observable '" << name << "' defined twice";
        throw InvalidArgumentException(ss.str());
    }
    if (min >= max) {
        stringstream ss;
        ss << "Observable " << name << " has min >= max, i.e., empty range";
        throw InvalidArgumentException(ss.str());
    }
    if(nbins==0){
        stringstream ss;
        ss << "Observable '" << name << "' has no bins";
        throw InvalidArgumentException(ss.str());
    }
    
    ObsId result(next_oid_id);
    ++next_oid_id;
    oid_to_name[result.id] = name;
    name_to_oid[name] = result.id;
    oid_to_range[result.id] = make_pair(min, max);
    oid_to_nbins[result.id] = nbins;
    return result;
}

bool VarIdManager::parNameExists(const std::string & name) const {
    return name_to_pid.find(name) != name_to_pid.end();
}

bool VarIdManager::obsNameExists(const std::string & name) const {
    return name_to_oid.find(name) != name_to_oid.end();
}

std::string VarIdManager::getName(const ParId & id) const {
    std::map<size_t, std::string>::const_iterator it = pid_to_name.find(id.id);
    if (it == pid_to_name.end()) {
        throw NotFoundException("VarIdManager::getVarname: did not find given VarId.");
    }
    return it->second;
}

std::string VarIdManager::getName(const ObsId & id) const {
    std::map<size_t, std::string>::const_iterator it = oid_to_name.find(id.id);
    if (it == oid_to_name.end()) {
        throw NotFoundException("VarIdManager::getVarname: did not find given VarId.");
    }
    return it->second;
}

ParId VarIdManager::getParId(const std::string & name) const {
    std::map<std::string, size_t>::const_iterator it = name_to_pid.find(name);
    if (it == name_to_pid.end()) {
        stringstream ss;
        ss << __FUNCTION__ << ": did not find variable '" << name << "'";
        throw NotFoundException(ss.str());
    }
    return ParId(it->second);
}

ObsId VarIdManager::getObsId(const std::string & name) const {
    std::map<std::string, size_t>::const_iterator it = name_to_oid.find(name);
    if (it == name_to_oid.end()) {
        stringstream ss;
        ss << __FUNCTION__ << ": did not find variable '" << name << "'";
        throw NotFoundException(ss.str());
    }
    return ObsId(it->second);
}

size_t VarIdManager::get_nbins(const ObsId & id) const{
    std::map<size_t, size_t>::const_iterator it = oid_to_nbins.find(id.id);
    if (it == oid_to_nbins.end()) {
        throw NotFoundException("VarIdManager::getNbins: did not find given variable id.");
    }
    return it->second;
}

const pair<double, double> & VarIdManager::get_range(const ObsId & id) const{
    std::map<size_t, pair<double, double> >::const_iterator it = oid_to_range.find(id.id);
    if (it == oid_to_range.end()) {
        throw NotFoundException("VarIdManager::getRange: did not find given variable id.");
    }
    return it->second;
}

ObsIds VarIdManager::getAllObsIds() const{
    std::map<size_t, pair<double, double> >::const_iterator it = oid_to_range.begin();
    ObsIds result;
    for(; it!= oid_to_range.end(); ++it){
       result.insert(ObsId(it->first));
    }
    return result;
}

ParIds VarIdManager::getAllParIds() const{
    std::map<size_t, string>::const_iterator it = pid_to_name.begin();
    ParIds result;
    for(; it!= pid_to_name.end(); ++it){
       result.insert(ParId(it->first));
    }
    return result;
}

/* ParValues */
ParIds ParValues::getParameters() const {
    ParIds result;
    for (size_t i=0; i<values.size(); i++) {
        if(not isnan(values[i])){
            result.insert(ParId(i));
        }
    }
    return result;
}

void ParValues::fail_get(const ParId & pid) const{
    std::stringstream ss;
    ss << "ParValues::get: given VarId " << pid.id << " not found";
    throw NotFoundException(ss.str());
}


void VarIdManagerUtils::apply_settings(plugin::Configuration & ctx){
    SettingWrapper s = ctx.setting;
    if(s.exists("observables")){
        size_t nobs = s["observables"].size();
        for (size_t i = 0; i < nobs; ++i) {
            string obs_name = s["observables"][i].getName();
            double min = s["observables"][i]["range"][0].get_double_or_inf();
            double max = s["observables"][i]["range"][1].get_double_or_inf();
            unsigned int nbins = s["observables"][i]["nbins"];
            ctx.vm->createObsId(obs_name, static_cast<size_t>(nbins), min, max);
        }
    }
    if(s.exists("parameters")){
        //get parameters:
        size_t npar = s["parameters"].size();
        for (size_t i = 0; i < npar; i++) {
            string par_name = s["parameters"][i];
            ctx.vm->createParId(par_name);
        }
    }
}

