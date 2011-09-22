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

   double * allocate_doubles(size_t n){
      ++n_allocs;
      double * result = 0;
      const size_t n_orig = n;
      //always allocate an even number of bins:
      if(n_orig % 2) ++n;
      //for the add_fast routine, which might use SSE optimizations, we need this alignment. And
      // while we at it, we should make sure double is as expected:
      BOOST_STATIC_ASSERT(sizeof(double)==8);
      int err = posix_memalign(reinterpret_cast<void**>(&result), 16, sizeof(double) * n);
      if(err!=0){
        throw std::bad_alloc();
      }
      //set the extra allocated double to zero to make sure no time-consuming garbage is there ...
      if(n_orig % 2) result[n-1] = 0.0;
      return result;
   }
   
   // note: data can be NULL
   void free_doubles(double * data){
      ++n_frees;
      free(data);
   }
}

void get_allocs_frees(int & n_alls, int & n_frs){
    n_alls = n_allocs;
    n_frs = n_frees;
}

DoubleVector::DoubleVector(size_t n): data(0), n_data(n){
    if(n_data > 0){
       data = allocate_doubles(n_data);
       set_all_values(0.0);
    }
}

DoubleVector::~DoubleVector(){
    free_doubles(data);
}

DoubleVector::DoubleVector(const DoubleVector & rhs): data(0), n_data(rhs.n_data){
   if(n_data > 0){
       data = allocate_doubles(n_data);
       memcpy(data, rhs.data, sizeof(double) * n_data);
   }
}

void DoubleVector::operator=(const DoubleVector & rhs){
    if(&rhs == this) return;
    if(n_data != rhs.n_data){
        free_doubles(data);
        if(rhs.n_data > 0)
            data = allocate_doubles(rhs.n_data);
        else
            data = 0;
        n_data = rhs.n_data;
    }
    if(n_data > 0)
       memcpy(data, rhs.data, sizeof(double) * n_data);
}

void DoubleVector::reset_n(size_t new_n){
    if(n_data!=new_n){
        free_doubles(data);
        if(new_n > 0)
            data = allocate_doubles(new_n);
        else
            data = 0;
        n_data = new_n;
    }
    set_all_values(0.0);
}

Histogram1D::Histogram1D(size_t b, double x_min, double x_max) : DoubleVector(b), xmin(x_min), xmax(x_max) {
    if(xmin >= xmax) throw InvalidArgumentException("Histogram: xmin >= xmax not allowed");
    set_all_values(0.0);
}

Histogram1D::Histogram1D(const Histogram1D& rhs): DoubleVector(rhs), xmin(rhs.xmin), xmax(rhs.xmax){
}

void Histogram1D::operator=(const Histogram1D & rhs) {
    DoubleVector::operator=(rhs);
    xmin = rhs.xmin;
    xmax = rhs.xmax;
}

Histogram1D::~Histogram1D() {
}

void Histogram1D::reset_range(double x_min, double x_max){
    if(size() > 0 && x_min >= x_max) throw InvalidArgumentException("Histogram: xmin >= xmax not allowed");
    xmin = x_min;
    xmax = x_max;
}

void Histogram1D::multiply_with_ratio_exponented(const Histogram1D & nominator, const Histogram1D & denominator, double exponent){
   check_compatibility(nominator);
   check_compatibility(denominator);
   const double * n_data = nominator.getData();
   const double * d_data = denominator.getData();
   double * data = getData();
   const size_t n = size();
   for(size_t i=0; i<n; i++){
      if(d_data[i]>0.0)
         data[i] *= pow(n_data[i] / d_data[i], exponent);
   }
}

void Histogram1D::fill(double xvalue, double weight) {
    int bin = static_cast<int> ((xvalue - xmin) * get_nbins() / (xmax - xmin));
    if (bin < 0 || (bin==0 && xvalue < xmin)) return;
    if (static_cast<size_t> (bin) >= size()) return;
    getData()[bin] += weight;
}

void Histogram1D::fail_check_compatibility(const Histogram1D & h) const {
    if (get_nbins() != h.get_nbins() || xmin!=h.xmin || xmax != h.xmax){
        std::stringstream s;
        s <<  "Histogram::check_compatibility: Histograms are not compatible (nbins, xmin, xmax) are: "
              " (" << get_nbins() << ", " << xmin << ", " << xmax << ") and "
              " (" << h.get_nbins() << ", " << h.xmin << ", " << h.xmax << ")";
        throw InvalidArgumentException(s.str());
    }
}

void Histogram1D::operator*=(const Histogram1D & h) {
    check_compatibility(h);
    const double * hdata = h.getData();
    double * data = getData();
    const size_t n = size();
    for (size_t i = 0; i < n; ++i) {
        data[i] *= hdata[i];
    }
}
