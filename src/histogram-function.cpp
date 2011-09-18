#include "interface/histogram-function.hpp"
#include "interface/random.hpp"
#include "interface/cfg-utils.hpp"

REGISTER_PLUGIN_BASETYPE(theta::HistogramFunction);

using namespace theta;

const Histogram1D &  ConstantHistogramFunctionError::getRandomFluctuation(Random & rnd, const ParValues & values) const{
    const size_t nbins = h.get_nbins();
    for(size_t i=0; i<nbins; ++i){
        double c = h.get(i);
        double err_i = err.get(i);
        if(err_i==0.0){
            fluc.set(i, c);
        }
        else{
            double factor = -1.0;
            //factor is gaussian around 1.0 with the current error, truncated at zero:
            while(factor < 0.0){
                factor = 1.0 + rnd.gauss(err_i);
            }
            fluc.set(i, factor * c);
        }
    }
    return fluc;
}
