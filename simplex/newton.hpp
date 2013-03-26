#ifndef NEWTON_HPP
#define NEWTON_HPP

#include "interface/minimizer.hpp"
#include "simplex/newton-utils.hpp"
#include <vector>
#include <map>

namespace newton_internal{
    class Linesearch;
}


/** \brief Quasi-Newton Minimizer
 * 
 * Minimizer implementing the Quasi-Newton method of calculating updated estimates for the
 * covariance matrix.
 * 
 * Configured via a setting group like
 * \begincode
 * minimizer = {
 *    type = "newton_minimizer";
 *    par_eps = 1e-4; // optional; default is 1e-3 (1e-6) for improve_cov = false (true)
 *    maxit = 100000; // optional; default is 10k
 *    improve_cov = true; // optional; default is false
 * };
 * \endcode
 * 
 * Any quasi-Newton method performs minimization in several steps:
 * 1. At the current point x0, using the gradient g and current estimate of the Hessian, propose a search vector v (direction and magnitude). If the
 *    function is actually quadratic with the current Hessian estimate and the gradient is accurate, x0 + v is the actual minimum.
 * 2. using the direction from 1., attempt minimizing the function along this direction [note that this search does not have to be very accurate.];
 *    this leads to an improved estimate of the minimum
 * 3. At the new minimum estimate, calculate the gradient and use it to update the covariance estimate; start over at 1. at the new point.
 * 
 * The ingredients to specify the algorithm are:
 *  * The linesearch routine (which is a 1D minimization). Currently, the only supported type is "brent", which uses Brent's algorithm to find the minimum
 *    along the line with an accuracy better than eps * (line search step).
 *  * The algorithm used to update the estimate of the Hessian matrix. Here, SR1 is used; see http://en.wikipedia.org/wiki/SR1_formula
 *  * The stopping criterion. Here iteration is stopped either after maxit iterations (in this case, a failure is reported), or if the step from the old to the new point
 *    is "small" in the sense that it is smaller than par_eps / step in each parameter, where step is the initial step size for this parameter.
 * 
 * Given that, most parameters have a straight-forward meaning:
 * 
 * \c par_eps controls the stopping criterion: once the step size in all parameters is smaller than this value (compared to this parameter's step size), the iteration stops.
 * 
 * \c maxit is the maximum number of iterations of the algorithm. One iteration executes all steps 1. to 3. as explained above.
 * 
 * \c improve_cov : if true, an additional step is performed after the actual minimization to improve the estimate for the covariance matrix by exploring
 *  function points and gradients around the found minimum. The tested points are at around 100 times the accuracy, so you should set
 *  \c par_eps to a small value when using this option.
 * 
 * Note that the algorithm typically makes n function evaluations at each point to calculate the gradient, where n is the number of parameters of f.
 * For the line search algorithm, another m ~ O(10) function evaluations are required. maxit controls the number of complete iterations of the
 * algorithm, leading to a maximum number of function evaluations of about maxit * (n + m).
 */
class newton_minimizer: public theta::Minimizer{
public:
    struct options{
        int maxit;
        double par_eps;
        bool debug;
        bool improve_cov;
        options(): maxit(10000), par_eps(1e-3), debug(false), improve_cov(false){}
    };

    newton_minimizer(const options & opts_);
    newton_minimizer(const theta::Configuration & cfg);
    virtual theta::MinimizationResult minimize(const theta::Function & f, const theta::ParValues & start,
                                               const theta::ParValues & step, const theta::Ranges & ranges);
    virtual theta::MinimizationResult minimize2(const theta::Function & f, const theta::FunctionInfo & info, const theta::ParValues & fixed_parameters);
    virtual boost::shared_ptr<theta::FunctionInfo> create_nll_function_info(const theta::Model & m, const boost::shared_ptr<theta::Distribution> & override_parameter_distribution,
                                                                            const theta::ParValues & fixed_parameters = theta::ParValues());
private:
    std::auto_ptr<newton_internal::Linesearch> ls;
    options opts;
};

#endif
