#include "plugins/exp_function.hpp"
#include "interface/plugin.hpp"

using namespace std;

exp_function::exp_function(const theta::plugin::Configuration & cfg): pid(cfg.vm->getParId(cfg.setting["parameter"])){
    par_ids.insert(pid);
    if(cfg.setting.exists("lambda_minus")){
        lambda_minus = cfg.setting["lambda_minus"];
        lambda_plus = cfg.setting["lambda_plus"];
    }
    else{
        lambda_plus = lambda_minus = cfg.setting["lambda"];
    }
}

double exp_function::operator()(const theta::ParValues & values) const{
    double val = values.get(pid);
    double lambda = val < 0 ? lambda_minus: lambda_plus;
    return exp(lambda * val);
}

REGISTER_PLUGIN(exp_function)

