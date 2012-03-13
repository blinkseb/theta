#include "interface/data.hpp"

using namespace theta;

Data theta::strip_uncertainties(const DataWithUncertainties & data){
    Data result;
    ObsIds observables = data.getObservables();
    for(ObsIds::const_iterator it=observables.begin(); it!=observables.end(); ++it){
        result[*it] = data[*it].get_values_histogram();
    }
    return result;
}

