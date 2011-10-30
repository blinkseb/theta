#include "plugins/add_squared.hpp"
#include "interface/plugin.hpp"
#include "interface/redirect_stdio.hpp"

using namespace theta;
using namespace std;

add_squared::add_squared(const Configuration & cfg){
    boost::shared_ptr<VarIdManager> vm = cfg.pm->get<VarIdManager>();
    size_t n = cfg.setting["parameters"].size();
    for(size_t i=0; i<n; ++i){
        ParId pid = vm->getParId(cfg.setting["parameters"][i]);
        v_pids.push_back(pid);
        par_ids.insert(pid);
        weights.push_back(cfg.setting["weights"][i]);
        if(weights.back() < 0.0) throw ConfigurationException("negative weight is not allowed");
    }
    theta_assert(weights.size() == v_pids.size());
}

double add_squared::operator()(const ParValues & v) const{
    double result = 0.0;
    for(size_t i=0; i<v_pids.size(); ++i){
        result += weights[i] * pow(v.get(v_pids[i]), 2);
    }
    return sqrt(result);
}


REGISTER_PLUGIN(add_squared)

