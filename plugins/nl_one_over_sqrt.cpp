#include "plugins/nl_one_over_sqrt.hpp"
#include "interface/exception.hpp"
#include "interface/plugin.hpp"

using namespace std;

nl_one_over_sqrt::nl_one_over_sqrt(const theta::plugin::Configuration & cfg): pid(cfg.vm->getParId(cfg.setting["parameter"])){
    par_ids.insert(pid);
}

double nl_one_over_sqrt::operator()(const theta::ParValues & values) const{
    double val = values.get(pid);
    if(val < 0.0) throw theta::MathException("nl_one_over_sqrt: negative argument");
    return 0.5 * val;
}

REGISTER_PLUGIN(nl_one_over_sqrt)
