#include "interface/minimizer.hpp"
#include "interface/plugin.tcc"

REGISTER_PLUGIN_BASETYPE(theta::Minimizer);

void theta::MinimizationResult::operator=(const MinimizationResult& rhs){
    fval = rhs.fval;
    values.set(rhs.values);
    errors_plus.set(rhs.errors_plus);
    errors_minus.set(rhs.errors_minus);
    covariance = rhs.covariance;
}

theta::Minimizer::~Minimizer(){}
