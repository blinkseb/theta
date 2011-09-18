#ifndef PLUGINS_MCMC_MINIMIZER_HPP
#define PLUGINS_MCMC_MINIMIZER_HPP

#include "interface/decls.hpp"
#include "interface/variables.hpp"
#include "interface/database.hpp"
#include "interface/minimizer.hpp"
#include "interface/random-utils.hpp"
#include "interface/matrix.hpp"

#include <string>

/** \brief Find the minimum using Markov-Chain Monte-Carlo
 *
 * Configuration is done via a setting group like
 * \code
 * minimizer = {
 *   type = "mcmc_minimizer";
 *   name = "min0";
 *   iterations = 10000;
 *   burn-in = 100; //optional. default is iterations / 10
 *   stepsize_factor = 0.1; //optional, default is 1.0
 *   after_minimizer = {...}; //optional.
 * };
 *
 * \endcode
 *
 * \c type is always "mcmc_minimizer" to select this Minimizer.
 * 
 * \c name is a name required for saving the random seed in the RndInfoTable.
 * 
 * \c iterations is the number of MCMC iterations.
 *
 * \c burn_in is the number of MCMC iterations to do at the beginning and throw away.
 * 
 * \c width_factor control the MCMC jump kernel size, see below.
 * 
 * \c after_minimizer is a Minimizer specification for a minimizer to run after the MCMC chain, see below.
 * 
 * The Markov-Chain Monte-Carlo method will be run first. As jump kernel, a diagonal
 * covariance matrix is used, based on the \c step parameter of the \c minimize() method. These
 * are scaled by \c stepsize_factor.
 * at the minimum are saved in MinimizationResult::values. The covariance
 * matrix for the whole chain (excluding burn-in) is used to set MinimizationResult::errors_plus
 * and errors_minus (to the same value).
 * 
 * Due to the randomness of the method (and the non-adaptiveness regarding the step size), the
 * result is not expected to be very accurate. Therefore, it is recommended to run another, more
 * precise Minimization method. This can be done in various ways, for example using \c last_minimizer in
 * \link minimizer_chain minimizer_chain \endlink or by setting \c after_minimizer. If set, this minimizer
 * will be run using the result of the MCMC method as starting point.
 */
class mcmc_minimizer: public theta::Minimizer, public theta::RandomConsumer{
public:
    /// \brief Constructor used by the plugin system to build an instance from settings in a configuration file
    mcmc_minimizer(const theta::plugin::Configuration & ctx);
    
    virtual theta::MinimizationResult minimize(const theta::Function & f, const theta::ParValues & start,
                const theta::ParValues & step, const std::map<theta::ParId, std::pair<double, double> > & ranges);
    
    virtual std::auto_ptr<theta::Minimizer> clone(const PropertyMap & pm) const;
    
private:
    mcmc_minimizer(const mcmc_minimizer & rhs, const PropertyMap & pm);
    
    std::auto_ptr<theta::Minimizer> after_minimizer;
    std::string name;
    unsigned int iterations;
    unsigned int burn_in;
    double stepsize_factor;
};

#endif
