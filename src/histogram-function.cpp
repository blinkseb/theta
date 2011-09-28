#include "interface/histogram-function.hpp"
#include "interface/random.hpp"
#include "interface/cfg-utils.hpp"
#include "interface/plugin.tcc"

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


const Histogram1D & ConstantHistogramFunction::operator()(const ParValues & values) const{
    return h;
}

Histogram1D ConstantHistogramFunction::get_histogram_dimensions() const{
    return h;
}

void ConstantHistogramFunction::set_histo(const Histogram1D & h_){
    h = h_;
}

ConstantHistogramFunction::ConstantHistogramFunction(){}


ConstantHistogramFunctionError::ConstantHistogramFunctionError(){}

const Histogram1D & ConstantHistogramFunctionError::operator()(const ParValues & values) const{
    return h;
}

Histogram1D ConstantHistogramFunctionError::get_histogram_dimensions() const{
    return h;
}


void ConstantHistogramFunctionError::set_histos(const Histogram1D & histo, const Histogram1D & error){
    histo.check_compatibility(error); //throws if not compatible
    h = histo;
    err = error;
    fluc = h;
    //check that errors are positive:
    for(size_t i=0; i<h.get_nbins(); ++i){
        if(error.get(i)<0.0) throw InvalidArgumentException("ConstantHistogramFunctionError: error histogram contains negative entries");
    }
}

