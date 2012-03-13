#include "interface/histogram-function.hpp"
#include "interface/random.hpp"
#include "interface/cfg-utils.hpp"
#include "interface/plugin.tcc"

REGISTER_PLUGIN_BASETYPE(theta::HistogramFunction);

using namespace theta;

const Histogram1DWithUncertainties & ConstantHistogramFunction::operator()(const ParValues & values) const{
    return h;
}

Histogram1DWithUncertainties ConstantHistogramFunction::get_histogram_dimensions() const{
    return h;
}

void ConstantHistogramFunction::set_histo(const Histogram1DWithUncertainties & h_){
    h = h_;
}

ConstantHistogramFunction::ConstantHistogramFunction(){}

