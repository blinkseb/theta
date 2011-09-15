#include "interface/histogram.hpp"
#include "interface/random.hpp"
#include "interface/exception.hpp"
#include "interface/utils.hpp"

#include <cmath>
#include <sstream>
#include <limits>
#include <new>

using namespace theta;

namespace{
   int n_allocs = 0;
   int n_frees = 0;

   double * allocate_histodata(size_t nbins){
      ++n_allocs;
      double * result = 0;
      const size_t nbins_orig = nbins;
      //always allocate an even number of bins:
      if(nbins_orig % 2) ++nbins;
      //for the add_fast routine, which might use SSE optimizations, we need this alignment. And
      // while we at it, we should make sure double is as expected:
      BOOST_STATIC_ASSERT(sizeof(double)==8);
      int err = posix_memalign(reinterpret_cast<void**>(&result), 16, sizeof(double) * (nbins + 2));
      if(err!=0){
        throw std::bad_alloc();
      }
      //set the extra allocated double to zero to make sure no time-consuming garbage is there ...
      if(nbins_orig % 2) result[nbins + 1] = 0.0;
      return result;
   }
   
   void free_histodata(double * histodata){
      ++n_frees;
      free(histodata);
   }
}

void get_allocs_frees(int & n_alls, int & n_frs){
    n_alls = n_allocs;
    n_frs = n_frees;
}


Histogram::Histogram(size_t b, double x_min, double x_max) : histodata(0), nbins(b), xmin(x_min), xmax(x_max) {
    if(xmin >= xmax) throw InvalidArgumentException("Histogram: xmin >= xmax not allowed");
    histodata = allocate_histodata(nbins);
    reset();
}

Histogram::Histogram(const Histogram& rhs) {
    initFromHisto(rhs);
}

void Histogram::operator=(const Histogram & rhs) {
    if (&rhs == this) return;
    if (nbins != rhs.nbins || xmin != rhs.xmin || xmax != rhs.xmax) {
        free_histodata(histodata);
        initFromHisto(rhs);
    } else {
        memcpy(histodata, rhs.histodata, sizeof (double) *(nbins + 2));
    }
}

Histogram::~Histogram() {
    free_histodata(histodata);
}

void Histogram::reset(size_t b, double x_min, double x_max) {
    //only re-allocate if there where changes.
    if (b > 0 && (b != nbins || x_min != xmin || x_max != xmax)) {
        nbins = b;
        xmin = x_min;
        xmax = x_max;
        if(xmin >= xmax) throw InvalidArgumentException("Histogram: xmin >= xmax not allowed");
        free_histodata(histodata);
        histodata = allocate_histodata(nbins);
    }
    memset(histodata, 0, sizeof (double) *(nbins + 2));
}

void Histogram::reset_to_1(){
   for(size_t i=0; i<=nbins+1; i++){
      histodata[i] = 1.0;
   }
}

void Histogram::multiply_with_ratio_exponented(const Histogram & nominator, const Histogram & denominator, double exponent){
   check_compatibility(nominator);
   check_compatibility(denominator);
   const double * n_data = nominator.histodata;
   const double * d_data = denominator.histodata;
   double s = 0.0;
   for(size_t i=0; i<=nbins+1; i++){
      if(d_data[i]>0.0)
         histodata[i] *= pow(n_data[i] / d_data[i], exponent);
      s += histodata[i];
   }
}

void Histogram::initFromHisto(const Histogram & h) {
    nbins = h.nbins;
    xmin = h.xmin;
    xmax = h.xmax;
    histodata = allocate_histodata(nbins);
    memcpy(histodata, h.histodata, sizeof (double) *(nbins + 2));
}

void Histogram::fill(double xvalue, double weight) {
    int bin = static_cast<int> ((xvalue - xmin) * nbins / (xmax - xmin) + 1);
    if (bin < 0)
        bin = 0;
    if (static_cast<size_t> (bin) > nbins + 1)
        bin = nbins + 1;
    histodata[bin] += weight;
}

void Histogram::fail_check_compatibility(const Histogram & h) const {
    if (nbins != h.nbins || xmin!=h.xmin || xmax != h.xmax){
        std::stringstream s;
        s <<  "Histogram::check_compatibility: Histograms are not compatible (nbins, xmin, xmax) are: "
              " (" << nbins << ", " << xmin << ", " << xmax << ") and "
              " (" << h.nbins << ", " << h.xmin << ", " << h.xmax << ")";
        throw InvalidArgumentException(s.str());
    }
}

void Histogram::operator*=(const Histogram & h) {
    check_compatibility(h);
    const double * hdata = h.histodata;
    for (size_t i = 0; i <= nbins + 1; ++i) {
        histodata[i] *= hdata[i];
    }
}

