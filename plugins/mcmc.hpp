#ifndef PLUGINS_MCMC_HPP
#define PLUGINS_MCMC_HPP

#include "interface/matrix.hpp"
#include "interface/exception.hpp"
#include "interface/random.hpp"
#include "interface/phys.hpp"

#include <vector>
#include <cmath>

#include <boost/scoped_array.hpp>

namespace theta {
/** Run the metropolis-Hastings Markov-Chain Monte-Carlo algorithm.
 *
 * @param nllikelihood is a function object which must implement double operator()(const double*) const which returns
 *  the negative logarithm of the likelihood (rather: the posterior) to integrate.
 * @param res is a "container" where the Markov Chain will be saved. It must implement the methods:
 * <ul>
 *   <li><tt>size_t getnpar()</tt> the number of parameters this result object was constructed for (required here
 *      for consistency checks only)</li>
 *   <li><tt>void fill(const double * x, double nll, size_t n)</tt>: this function will be called to indicate that the
 *       chain has <tt>n</tt> points at parameter values <tt>x</tt> and negative logarithm value of <tt>nll</tt>. What
 *       the result object does with this information is not of our concern here. Typically, it will either save the complete
 *       chain or fill the x values into a histogram</li>
 * </ul>
 * @param rand is the random number generator to use to compute the next candidate point of the chain.
 *   It has to implement the <tt>double gauss()</tt> and <tt>double uniform()</tt>
 *   methods which shall return a random number distributed according to a gaussian around 0 with width 1 and a uniform distribution
 *   in the intervall [0, 1], respectively (whether or not the endpoints are actually included in the uniform case does not matter here).
 * @param startvalues is the point where to start the Markov Chain.
 * @param sqrt_cov is (apart from a overall factor) the matrix used to transform the uniform-Gaussian
 *   jumping kernel before adding it to the current point. It should be set to the Cholesky decomposition of the
 *   covariance matrix of the likelihood (or an approximation thereof), hence the name "square root of covariance".
 * @param iterations is the number of iterations for the Markov Chain which will be reported to the \c res object.
 * @param burn_in is the number of Markov-Chain iterations run at the beginning which are <em>not</em> reported to the \c res
 *   object.
 */
template<class resulttype>
void metropolisHastings(const Function & nllikelihood, resulttype &res, Random & rand,
        const std::vector<double> & startvalues, const Matrix & sqrt_cov, size_t iterations, size_t burn_in, bool ignore_inf_nll = false) {    
    const size_t npar = startvalues.size();
    if(npar != sqrt_cov.get_n_rows() || npar!=sqrt_cov.get_n_cols() || npar!=nllikelihood.get_parameters().size() || npar!=res.getnpar())
        throw std::invalid_argument("metropolisHastings: dimension/size of arguments mismatch");
    size_t npar_reduced = npar;
    for(size_t i=0; i<npar; i++){
        if(sqrt_cov(i,i)==0) --npar_reduced;
    }
    theta_assert(npar_reduced > 0);
    double factor = 2.38 / sqrt(npar_reduced);

    //keep the lower triangle   L    of the sqrt_cov to multiply that with deltaX:
    boost::scoped_array<double> lm(new double[npar * (npar + 1) / 2]);
    size_t z = 0;
    for (size_t i = 0; i < npar; i++) {
        for (size_t j = 0; j <= i; j++) {
            lm[z] = sqrt_cov(i,j) * factor;
            theta_assert(std::isfinite(lm[z]));
            z++;
        }
    }

    //0. initialization:
    boost::scoped_array<double> x(new double[npar]);
    boost::scoped_array<double> x_new(new double[npar]);
    boost::scoped_array<double> dx(new double[npar]);
    
    //set the starting point:
    std::copy(startvalues.begin(), startvalues.end(), &x[0]);
    double nll = nllikelihood(x.get());
    if(!std::isfinite(nll) and not ignore_inf_nll) throw Exception("first nll value was inf");
    // if this exception is thrown, it means that the likelihood function at the start values was inf or NAN.
    // One common reason for this is that the model produces a zero prediction in some bin where there is >0 data which should
    // be avoided by the model (e.g., by filling in some small number in all bins with content zero).

    const size_t iter = burn_in + iterations;
    size_t weight = 1;
    //note: splitting up this for-loop
    // into two loops (one for burn-in, one for
    // recording) seems to save some time for the saved if(it>=burn_in)
    // but it is actually slower ... (tested with gcc 4.3.3, -O3).
    for (size_t it = 1; it < iter; it++) {
        for (size_t i = 0; i < npar; i++) {
            dx[i] = rand.gauss();
        }
        //multiply lm with dx and add to x to get x_new:
        z = 0;
        for (size_t i = 0; i < npar; ++i) {
            x_new[i] = x[i];
            for (size_t j = 0; j <= i; ++j) {
                x_new[i] += lm[z] * dx[j];
                z++;
            }
        }
        double nll_new = nllikelihood(x_new.get());
        if ((nll_new <= nll) || (rand.uniform() < exp(nll - nll_new))) {
            if(it > burn_in){
                res.fill(x.get(), nll, weight);
                weight = 1;
            }
            x.swap(x_new);
            nll = nll_new;
        } else if(it > burn_in){
            ++weight;
        }
    }
    res.fill(x.get(), nll, weight);
}

/** \brief estimate the square root (cholesky decomposition) of the covariance matrix of the likelihood function
 *
 * The method will start a Markov chain at the given startvalues with the \c iterations iterations.
 * The found covariance matrix is used in a next pass where the Markov chain jumping rules
 * are adjusted to the covariance found in the previous step. This procedure is repeated until the
 * estimate is stable.
 *
 * \param rnd is the random number generator to use
 * \param model is the Model to use to get the asimov data from
 * \param[out] startvalues will contain the suggested startvalues. The contents when calling this function will be ignored.
 */
Matrix get_sqrt_cov2(Random & rnd, const Model & model, std::vector<double> & startvalues,
                    const boost::shared_ptr<theta::Distribution> & override_parameter_distribution,
                    const boost::shared_ptr<theta::Function> & additional_nll_term);

/** \brief Calculate the cholesky decomposition, but allow zero eigenvalues.
 *
 *
 * Writes the cholesky decomposition into \c result. This function treats the case of fixing a parameter through
 * setting its covariance diagonal entry to zero correctly in only decomposing the non-zero part of the matrix
 * and keeping the zero entries where they were.
 *
 * If \c expect_reduced is non-negative, it will be checked whether it matches the determined number
 * of reduced (=non-fixed) dimensions. If it does not, an Exception will be thrown. If \c expect_reduced
 * is negative, no check will be done.
 */
void get_cholesky(const Matrix & cov, Matrix & result, int expect_reduced = -1);

}

#endif
