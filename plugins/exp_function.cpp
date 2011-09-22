#include "plugins/exp_function.hpp"
#include "interface/plugin.hpp"

using namespace std;
using namespace theta;

exp_function::exp_function(const theta::Configuration & cfg){
    ParValues val_lambdas_plus, val_lambdas_minus;
    boost::shared_ptr<VarIdManager> vm = cfg.pm->get<VarIdManager>();
    if(cfg.setting.exists("parameters")){
        for(size_t i=0; i<cfg.setting["parameters"].size(); ++i){
            ParId pid = vm->getParId(cfg.setting["parameters"][i]);
            par_ids.insert(pid);
            val_lambdas_plus.set(pid, cfg.setting["lambdas_plus"][i]);
            val_lambdas_minus.set(pid, cfg.setting["lambdas_minus"][i]);
        }
    }
    else{
        ParId pid = vm->getParId(cfg.setting["parameter"]);
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
    lambdas_plus.reserve(par_ids.size());
    lambdas_minus.reserve(par_ids.size());
    v_pids.reserve(par_ids.size());
    for(ParIds::const_iterator it=par_ids.begin(); it!=par_ids.end(); ++it){
        v_pids.push_back(*it);
        lambdas_plus.push_back(val_lambdas_plus.get(*it));
        lambdas_minus.push_back(val_lambdas_minus.get(*it));
    }
}

double exp_function::operator()(const theta::ParValues & values) const{
    double exponent_total = 0.0;
    const size_t n = v_pids.size();
    for(size_t i=0; i<n; ++i){
        double val = values.get_unchecked(v_pids[i]);
        double lambda = val < 0 ? lambdas_minus[i] : lambdas_plus[i];
        exponent_total += lambda * val;
    }
    return exp(exponent_total);
}

void exp_function::codegen(std::ostream & out, const std::string & prefix, const theta::PropertyMap & pm) const{
    boost::shared_ptr<VarIdManager> vm = pm.get<VarIdManager>();
    out << "double " << prefix << "_evaluate(const double * par_values){" << endl
        << "    double exponent_total = 0.0, val;" << endl;
        for(size_t i=0; i<v_pids.size(); ++i){
            out << "    val = par_values[pindex_" << vm->getName(v_pids[i]) << "];" << endl;
            if(lambdas_minus[i]==lambdas_plus[i]){
                out << "    exponent_total += " << codegen::dtos(lambdas_minus[i]) << " * val;" << endl;
            }
            else{
                out << "    exponent_total += (val < 0 ? " << codegen::dtos(lambdas_minus[i]) << ":" << codegen::dtos(lambdas_plus[i]) << ") * val;" << endl;
            }
        }
    out << "    return utils::exp(exponent_total);" << endl;
    out << "}" << endl << endl;
}

std::auto_ptr<theta::Function> exp_function::clone() const{
    return std::auto_ptr<theta::Function>(new exp_function(*this));
}

REGISTER_PLUGIN(exp_function)
