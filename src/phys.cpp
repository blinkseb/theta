#include "interface/phys.hpp"

using namespace theta;

REGISTER_PLUGIN_BASETYPE(Function);
REGISTER_PLUGIN_BASETYPE(DataSource);

/* DATA */
ObsIds Data::getObservables() const{
    ObsIds result;
    std::vector<Histogram>::const_iterator it = data.begin();
    size_t i=0;
    for(;it!=data.end(); ++it, ++i){
        if(it->get_nbins()!=0) result.insert(ObsId(i));
    }
    return result;
}

void Data::fail_get(const ObsId & oid) const{
    throw NotFoundException("Data::operator[]() const: no data found for given ObsId");
}
