#include "interface/phys.hpp"

using namespace theta;

REGISTER_PLUGIN_BASETYPE(Function);
REGISTER_PLUGIN_BASETYPE(DataSource);


ObsIds Data::getObservables() const{
    ObsIds result;
    std::vector<Histogram1D>::const_iterator it = data.begin();
    size_t i=0;
    for(;it!=data.end(); ++it, ++i){
        if(it->get_nbins()!=0) result.insert(ObsId(i));
    }
    return result;
}

void Data::fail_get(const ObsId & oid) const{
    throw NotFoundException("Data::operator[]() const: no data found for given ObsId");
}



void Function::codegen(std::ostream & out, const std::string & prefix, const PropertyMap & pm) const{
    out  << "double " << prefix << "_evaluate(const double * par_values){" << std::endl
         << "    return reinterpret_cast<Function*>(" << reinterpret_cast<size_t>(this) << ")->operator()(ParValues(par_values, all_par_ids));" << std::endl
         << "}" << std::endl;
}



