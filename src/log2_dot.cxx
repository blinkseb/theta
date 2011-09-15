// fallback implementation for log2_dot.s which is implemented specifically for 64-bit architectures

#include "interface/log2_dot.hpp"
#include <cstddef>
#include <limits>
#include <math.h>

#include "interface/utils.hpp"

double log2_dot(const double * x, const double * y, unsigned int n){
    double result = 0.0;
    for(size_t i=0; i<n; ++i){
        result += y[i] * log2(x[i]);
    }
    return result;
}

double template_nllikelihood(const double * data, const double * pred, unsigned int n){
   double result = 0.0;
   for(unsigned int i=0; i<n; ++i){
        result += pred[i];
        if(pred[i] > 0.0){
             if(data[i] > 0.0){
                 result -= data[i] * theta::utils::log(pred[i]);
             }
         }else if(data[i] > 0.0){
             return std::numeric_limits<double>::infinity();
         }
    }
    return result;
}

