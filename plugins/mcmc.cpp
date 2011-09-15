#include "plugins/mcmc.hpp"
#include "plugins/mcmc-result.hpp"
#include "plugins/asimov_likelihood_widths.hpp"
#include "interface/utils.hpp"
#include "interface/distribution.hpp"
#include "interface/model.hpp"
#include "interface/redirect_stdio.hpp"

#include <algorithm>
#include <cassert>

#include <limits>

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <sstream>

using namespace std;

namespace theta{

void get_cholesky(const Matrix & cov, Matrix & result, int expect_reduced){
    size_t npar_reduced = cov.getRows();
    size_t npar = npar_reduced;
    result.reset(npar, npar);
    //get the overall "scale" of the matrix by looking at the largest diagonal element:
    double m = fabs(cov(0,0));
    for(size_t i=0; i<cov.getRows(); i++){
        m = max(m, fabs(cov(i,i)));
    }
    bool par_has_zero_cov[npar];
    for (size_t i = 0; i < npar; i++) {
       if((par_has_zero_cov[i] = utils::close_to(cov(i,i), 0, m)))
            npar_reduced--;
    }
    if(npar_reduced==0){
        throw InvalidArgumentException("get_cholesky: number of reduced dimensions is zero (all parameters fixed?)");
    }
    if(expect_reduced > 0 && static_cast<size_t>(expect_reduced)!=npar_reduced){
        throw InvalidArgumentException("get_cholesky: number of reduced dimensions not as expected");
    }
    //calculate cholesky decomposition    cov = L L^T    of the reduced covariance matrix:
    Matrix cov_c(npar_reduced, npar_reduced);
    size_t row = 0;
    for (size_t i = 0; i < npar; i++) {
        if(par_has_zero_cov[i])continue;
        size_t col = 0;
        for (size_t j = 0; j < npar; j++) {
            if(par_has_zero_cov[j])continue;
            cov_c(row, col) = cov(i,j);
            col++;
        }
        row++;
    }
    cov_c.cholesky_decomposition();
    //store the result in result, correctly filling the zero vectors ...
    row = 0;
    for(size_t i=0; i<npar; i++){
        if(par_has_zero_cov[i])continue;
        size_t col = 0;
        for (size_t j = 0; j < npar; j++) {
            if(par_has_zero_cov[j])continue;
            result(i,j) = cov_c(row, col);
            col++;
        }
        row++;
    }
}

bool jump_rates_converged(const vector<double> & jump_rates){
    if(jump_rates.size() < 2) return false;
    const size_t n = jump_rates.size();
    double last_rate = jump_rates[n-1];
    double nlast_rate = jump_rates[n-2];
    if(last_rate < 0.05) return false;
    if(last_rate > 0.5) return false;
    //if jump rate looks reasonable and did not change too much in the last iterations: break:
    if(fabs((nlast_rate - last_rate) / max(last_rate, nlast_rate)) < 0.10) return true;
    //TODO: otherwise, it could still have converged well enough, but it is difficult to have a very good
    // criterion here ...
    return false;
}


Matrix get_sqrt_cov2(Random & rnd, const Model & model, std::vector<double> & startvalues,
                    const boost::shared_ptr<theta::Distribution> & override_parameter_distribution,
                    const boost::shared_ptr<VarIdManager> & vm){
    const size_t max_passes = 50;
    const size_t iterations = 8000;
    const size_t n = model.getParameters().size();
    Matrix sqrt_cov(n, n);
    Matrix cov(n, n);
    startvalues.resize(n);
    
    ParValues widths = asimov_likelihood_widths(model, override_parameter_distribution);
    
    ParIds par_ids = model.getParameters();
    size_t k=0;
    int n_fixed_parameters = 0;
    
    const Distribution & dist = override_parameter_distribution.get()? *override_parameter_distribution : model.get_parameter_distribution();
    ParValues pv_start;
    dist.mode(pv_start);
    Data asimov_data;
    model.get_prediction(asimov_data, pv_start);
    ParIds fixed_pars;
    for(ParIds::const_iterator it = par_ids.begin(); it!=par_ids.end(); ++it, ++k){
        double width = widths.get(*it) * 2.38 / sqrt(n);
        startvalues[k] = pv_start.get(*it);
        if(width==0.0){
            ++n_fixed_parameters;
            fixed_pars.insert(*it);
        }
        cov(k, k) = width*width;
    }
    assert(k==n);
    get_cholesky(cov, sqrt_cov, static_cast<int>(n) - n_fixed_parameters);
    
    std::auto_ptr<NLLikelihood> p_nll = model.getNLLikelihood(asimov_data);
    p_nll->set_override_distribution(override_parameter_distribution);
    NLLikelihood & nll = *p_nll;
    vector<double> jump_rates;
    jump_rates.reserve(max_passes);
    Result res(n);
    for (size_t i = 0; i < max_passes; i++) {
        res.reset();
        metropolisHastings(nll, res, rnd, startvalues, sqrt_cov, iterations, iterations/10);
        startvalues = res.getMeans();
        cov = res.getCov();
        get_cholesky(cov, sqrt_cov, static_cast<int>(n) - n_fixed_parameters);
        jump_rates.push_back(static_cast<double>(res.getCountDifferent()) / res.getCount());
        if(jump_rates_converged(jump_rates)) break;
    }
    if(jump_rates.size()==max_passes){
        theta::cout << "WARNING in get_sqrt_cov: covariance estimate did not really converge; jump rates were: ";
        for(size_t i=0; i<max_passes; ++i){
            theta::cout << jump_rates[i] << "; ";
        }
    }
    return sqrt_cov;
}

}
