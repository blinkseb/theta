#ifndef UTILS_HPP
#define UTILS_HPP

#include <cmath>
#include <string>

#ifdef USE_CRLIBM
#include "crlibm/crlibm.h"
#endif

#ifdef USE_AMDLIBM
#include "amdlibm/include/amdlibm.h"
#endif

#ifdef __SSE2__
#include <emmintrin.h>
#endif



namespace theta { namespace utils{

extern std::string theta_dir;
void fill_theta_dir(char** argv);

/// Replaces the string "$THETA_DIR" by the theta directory; to be used by plugins resolving filenames
std::string replace_theta_dir(const std::string & path);

double phi_inverse(double p);

/** \brief The lngamma function
 *
 * Forwards to the boost implementation which is thread save (note that
 * C99 implementations need not be therad save).
 */
double lngamma(double x);

/** \brief add 2 vectors, possibly with sse optimization
 *
 * Calculates x+=y. Requires x and y to be aligned at a 16-byte address. It is assumed
 * that an even number of doubles has been allocated for x, even if n is odd.
 */
inline void add_fast(double * x, const double * y, const size_t n){
#ifndef __SSE2__
   for(size_t i=0; i<n; ++i){
      x[i] += y[i];
   }
#else
  for(size_t i=0; i<n; i+=2){
     __m128d XMM0 = _mm_load_pd(x + i);
     __m128d XMM1 = _mm_load_pd(y + i);
     XMM0 = _mm_add_pd(XMM0, XMM1);
     _mm_store_pd(x+i  , XMM0);
  }
#endif
}

/** \brief multiply a vector with a constant, possibly with sse optimization
 *
 * Calculates x*=cy. Requires x to be aligned at a 16-byte address. It is assumed
 * that an even number of doubles has been allocated for x, even if n is odd.
 *
 */
inline void mul_fast(double * x, double c, const size_t n){
#ifndef __SSE2__
   for(size_t i=0; i<n; ++i){
      x[i] *= c;
   }
#else
  __m128d XMM2 = _mm_set1_pd(c);
  for(size_t i=0; i<n; i+=2){
     __m128d XMM0 = _mm_load_pd(x + i);
     XMM0 = _mm_mul_pd(XMM0, XMM2);
     _mm_store_pd(x+i, XMM0);
  }
#endif
}


/** \brief add 2 vectors, possibly with sse optimization
 *
 * Calculates x+=c * y. Requires x and y to be aligned at a 16-byte address. It is assumed
 * that an even number of doubles has been allocated for x, even if n is odd.
 */
inline void add_fast_with_coeff(double * x, const double * y, double c, const size_t n){
#ifndef __SSE2__
   for(size_t i=0; i<n; ++i){
      x[i] += c * y[i];
   }
#else
  __m128d XMM2 = _mm_set1_pd(c);
  for(size_t i=0; i<n; i+=2){
     __m128d XMM0 = _mm_load_pd(x + i);
     __m128d XMM1 = _mm_load_pd(y + i);
     XMM1 = _mm_mul_pd(XMM1, XMM2);
     XMM0 = _mm_add_pd(XMM0, XMM1);
     _mm_store_pd(x+i, XMM0);
  }
#endif
}

/** \brief possible redirections of log
 *
 * Tests have shown that the common log function is slow. In order
 * to allow easy switching for testing, all code in %theta should
 * use this log function.
 */
inline double log(double x){
#ifdef USE_CRLIBM
    return log_rn(x);
#else
#ifdef USE_AMDLIBM
    return amd_log(x);
#else
    return ::log(x);
#endif
#endif
}

inline double exp(double x){
#ifdef USE_AMDLIBM
    return amd_exp(x);
#else
    return ::exp(x);
#endif
}


}}

#endif
