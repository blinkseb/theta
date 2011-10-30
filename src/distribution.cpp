#include "interface/distribution.hpp"
#include "interface/exception.hpp"
#include "interface/plugin.tcc"

REGISTER_PLUGIN_BASETYPE(theta::Distribution);

void theta::DistributionUtils::fillModeSupport(theta::ParValues & mode,
                std::map<theta::ParId, std::pair<double, double> > & support, const theta::Distribution & d){
    ParIds pids = d.getParameters();
    d.mode(mode);
    theta_assert(mode.getParameters()==pids);
    for(ParIds::const_iterator p_it=pids.begin(); p_it!=pids.end(); ++p_it){
        support[*p_it] = d.support(*p_it);
    }
}

theta::Distribution::~Distribution(){}


theta::EmptyDistribution::~EmptyDistribution(){}


const std::pair<double, double> & theta::EmptyDistribution::support(const ParId & p) const{
    throw std::invalid_argument("EmptyDistribution::support not implemented");
}
