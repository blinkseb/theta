#include "interface/distribution.hpp"

REGISTER_PLUGIN_BASETYPE(theta::Distribution);

void theta::DistributionUtils::fillModeSupport(theta::ParValues & mode,
                std::map<theta::ParId, std::pair<double, double> > & support, const theta::Distribution & d){
    ParIds pids = d.getParameters();
    d.mode(mode);
    assert(mode.getParameters()==pids);
    for(ParIds::const_iterator p_it=pids.begin(); p_it!=pids.end(); ++p_it){
        support[*p_it] = d.support(*p_it);
    }
}
