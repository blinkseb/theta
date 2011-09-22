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


void HistogramFunction::codegen(std::ostream & out, const std::string & prefix, const PropertyMap & pm) const{
    out  << "void " << prefix << "_hf_add_with_coeff(double coeff, const double * par_values, double * histodata){" << std::endl
         << "    Histogram1D h = reinterpret_cast<HistogramFunction*>(" << reinterpret_cast<size_t>(this) << ")->operator()(ParValues(par_values, all_par_ids));" << std::endl
         << "    add_fast_with_coeff(histodata, h.getData(), coeff, h.size());" << std::endl
         << "}" << std::endl;
}



void ConstantHistogramFunction::codegen(std::ostream & out, const std::string & prefix, const PropertyMap & pm) const{
    using std::endl;
    const size_t nbins = h.get_nbins();
    out << endl;
    out << "/* start with ConstantHistogramFunction '" << prefix << "' */" << endl;
    out << "DoubleVector " << prefix << "_hf_h_init(){" << endl
        << "    const double data[" << nbins << "] = {";
        for(size_t i=0; i<nbins; ++i){
            if(i > 0) out << ", ";
            out << codegen::dtos(h.get(i));
        }
        out << "};" << endl;
        out << "    DoubleVector result(" << nbins << ");" << endl;
        out << "    memcpy(result.getData(), data, " << nbins * sizeof(double) << ");" << endl;
        out << "    return result;" << endl;
        out << "}" << endl << endl;
    out  << "const DoubleVector " << prefix << "_hf_h(" << prefix << "_hf_h_init());" << endl << endl;
    out  << "void " << prefix << "_hf_add_with_coeff(double coeff, const double * par_values, double * histodata){" << std::endl
         << "    add_fast_with_coeff(histodata, " << prefix << "_hf_h.getData(), coeff, " << nbins << ");" << std::endl
         << "}" << std::endl;
    out << "/* end of with ConstantHistogramFunction '" << prefix << "' */" << endl;
}


