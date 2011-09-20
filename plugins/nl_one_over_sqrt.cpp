#include "plugins/nl_one_over_sqrt.hpp"
#include "interface/exception.hpp"
#include "interface/plugin.hpp"

using namespace std;
using namespace theta;

nl_one_over_sqrt::nl_one_over_sqrt(const theta::Configuration & cfg): pid(cfg.pm->get<VarIdManager>()->getParId(cfg.setting["parameter"])){
    par_ids.insert(pid);
}

double nl_one_over_sqrt::operator()(const theta::ParValues & values) const{
    double val = values.get(pid);
    if(val < 0.0) throw theta::MathException("nl_one_over_sqrt: negative argument");
    return 0.5 * val;
}

std::auto_ptr<theta::Function> nl_one_over_sqrt::clone() const{
    return std::auto_ptr<theta::Function>(new nl_one_over_sqrt(*this));
}

REGISTER_PLUGIN(nl_one_over_sqrt)
