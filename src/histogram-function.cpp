#include "interface/histogram-function.hpp"
#include "interface/cfg-utils.hpp"
#include "interface/exception.hpp"
#include "interface/plugin.tcc"

REGISTER_PLUGIN_BASETYPE(theta::HistogramFunction);

using namespace theta;

void ConstantHistogramFunction::apply_functor(const functor<Histogram1DWithUncertainties> & f, const ParValues & values) const{
    f(h_wu);
}

void ConstantHistogramFunction::apply_functor(const functor<Histogram1D> & f, const ParValues & values) const{
    f(h);
}
        
void ConstantHistogramFunction::get_histogram_dimensions(size_t & nbins, double & xmin, double & xmax) const{
    nbins = h.get_nbins();
    xmin = h.get_xmin();
    xmax = h.get_xmax();
}

void ConstantHistogramFunction::set_histo(const Histogram1DWithUncertainties & h_){
    h = h_.get_values_histogram();
    h_wu = h_;
}

ConstantHistogramFunction::ConstantHistogramFunction(){}

Histogram1DWithUncertainties theta::get_constant_histogram(const Configuration & cfg){
    std::auto_ptr<HistogramFunction> hf = PluginManager<HistogramFunction>::build(cfg);
    if(hf->get_parameters().size()!=0){
        throw std::invalid_argument("Histogram defined in path '" + cfg.setting.get_path() + "' is not constant");
    }
    Histogram1DWithUncertainties res;
    hf->apply_functor(copy_to<Histogram1DWithUncertainties>(res), ParValues());
    return res;
}
