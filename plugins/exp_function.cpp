#include "plugins/exp_function.hpp"
#include "interface/plugin.hpp"
#include "interface/utils.hpp"

using namespace std;
using namespace theta;

exp_function::exp_function(const theta::Configuration & cfg){
    ParValues val_lambdas_plus, val_lambdas_minus;
    boost::shared_ptr<VarIdManager> vm = cfg.pm->get<VarIdManager>();
    if(cfg.setting.exists("parameters")){
        for(size_t i=0; i<cfg.setting["parameters"].size(); ++i){
            ParId pid = vm->get_par_id(cfg.setting["parameters"][i]);
            par_ids.insert(pid);
            val_lambdas_plus.set(pid, cfg.setting["lambdas_plus"][i]);
            val_lambdas_minus.set(pid, cfg.setting["lambdas_minus"][i]);
        }
    }
    else{
        ParId pid = vm->get_par_id(cfg.setting["parameter"]);
        par_ids.insert(pid);
        if(cfg.setting.exists("lambda_minus")){
            val_lambdas_minus.set(pid, cfg.setting["lambda_minus"]);
            val_lambdas_plus.set(pid, cfg.setting["lambda_plus"]);
        }
        else{
            double lambda = cfg.setting["lambda"];
            val_lambdas_minus.set(pid, lambda);
            val_lambdas_plus.set(pid, lambda);
        }
    }
    //convert to the "vector" version:
    lambdas_plusminus.reserve(2 * par_ids.size());
    n = par_ids.size();
    v_pids.reserve(par_ids.size());
    for(ParIds::const_iterator it=par_ids.begin(); it!=par_ids.end(); ++it){
        v_pids.push_back(*it);
        lambdas_plusminus.push_back(val_lambdas_plus.get(*it));
        lambdas_plusminus.push_back(val_lambdas_minus.get(*it));
    }
    
}

double exp_function::operator()(const theta::ParValues & values) const{
    double exponent_total = 0.0;
    for(size_t i=0; i<n; ++i){
        double val = values.get_unchecked(v_pids[i]);
        int offset = val > 0 ? 0: 1;
        exponent_total += lambdas_plusminus[2*i + offset] * val;
    }
    return theta::utils::exp(exponent_total);
}

REGISTER_PLUGIN(exp_function)
