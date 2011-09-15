#ifndef LOG2_DOT_HPP
#define LOG2_DOT_HPP

#include <cmath>

extern "C"{
    /** \brief Calculate inner product of a vector with log of the vector
     *
     * result = sum_i=0^{n-1}  y[i] * log2(x[i])
     */
    double log2_dot(const double * x, const double * y, unsigned int n);
    
    /** \brief calculate the negative-log-likelihood for the data and prediction
     *
     * calculates
     * sum_{i=0}^{n-1}  pred[i] - data[i] * log(pred[i])
     *
     * If, for some i, data[i] and pred[i] are both zero, this is i is skipped in
     * the sum (i.e., summand is treated as 0).
     *
     * If, for some i, data[i] > 0.0 and pred[i] <= 0.0, +infinity is returned.
     */
    double template_nllikelihood(const double * data, const double * pred, unsigned int n);
}

#endif
