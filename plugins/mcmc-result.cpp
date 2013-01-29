#include "plugins/mcmc-result.hpp"
#include "interface/exception.hpp"

#include <cmath>
#include <sstream>

using namespace theta;

Result::Result(size_t n) :
    npar(n), count(0), count_different_points(0), means(n, 0.0), count_covariance(n, n) {
}

void Result::fill(const double * p, double nll, size_t weight){
    if(!std::isfinite(nll)) return;
    //factor for the covariance ...
    double factor = count * 1.0 / (count + 1);
    for(size_t i=1; i<weight; ++i){
        factor += (count*count*1.0) / ((count+i) * (count+i+1));
    }
    for(size_t i=0; i<npar; ++i){
        double diff_i = p[i] - means[i];
        for(size_t j=i; j<npar; ++j){
            double diff_j = p[j] - means[j];
            count_covariance(i,j) += factor * diff_i * diff_j;
        }
        //the old means[i] is now no longer needed, as j>=i ...
        means[i] += weight * 1.0 / (count + weight) * diff_i;
    }
    count += weight;
    count_different_points++;
    fill2(p, nll, weight);
}

size_t Result::getnpar() const {
    return npar;
}

size_t Result::getCount() const {
    return count;
}

size_t Result::getCountDifferent() const {
    return count_different_points;
}

std::vector<double> Result::getMeans() const {
    return means;
}

Matrix Result::getCov() const {
    Matrix result(npar, npar);
    for (size_t i = 0; i < npar; i++) {
        for (size_t j = i; j < npar; j++) {
            result(j, i) = result(i, j) = count_covariance(i,j) / count;
        }
    }
    return result;
}

