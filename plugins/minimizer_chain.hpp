#ifndef PLUGINS_MINIMIZER_CHAIN_HPP
#define PLUGINS_MINIMIZER_CHAIN_HPP

#include "interface/decls.hpp"
#include "interface/minimizer.hpp"

#include <boost/ptr_container/ptr_vector.hpp>

#include <memory>

/** \brief Run a sequence of minimizers, until the first succeeds
 *
 * Configuration is done via a setting group like
 * \code
 * minimizer = {
 *   type = "minimizer_chain";
 *   minimizers = ("@minimizer1", "@minimizer2", "@minimizer3");
 *   last_minimizer = "@minimizer1"; // optional.
 * };
 * 
 * minimizer1 = { ... }; // some other minimizers ...
 * minimizer2 = { ... };
 * minimizer3 = { ... };
 *
 * \endcode
 *
 * \c type is always "minimizer_chain" to select this theta::Minimizer.
 * 
 * \c minimizers is the list of Minimizers to try, in the given order
 * 
 * \c last_minimizer is a Minimizer to run as last step
 *
 * Each of the minimizers given will be used to find a minimum, in the given order.
 * The minimum of the first successful minimization attempt will be returned. If
 * all minimizers fail, a MinimizationException according to the last minimizer in the chain
 * will be thrown.
 * 
 * If \c last_minimizer is given, it will be always run last, after another Minimizer succeeded.
 * as start value for \c last_minimizer, the parameter values at the minimum found so far is used; as step
 * sizes either the errors or (if these are not available), the original step sizes.
 * 
 * The main use case for \c last_minimizer is to provide error estimates: the minimizers in the chain
 * will probably not all provide an error estimate. Specifying a Minimizer here which provides errors
 * as \c last_minimizer ensures that all (successful) minimizations have an error estimate for the parameters.
 */
class minimizer_chain: public theta::Minimizer {
public:
    /// \brief Constructor used by the plugin system to build an instance from settings in a configuration file
    minimizer_chain(const theta::Configuration & cfg);
    
    virtual theta::MinimizationResult minimize(const theta::Function & f, const theta::ParValues & start,
                const theta::ParValues & step, const std::map<theta::ParId, std::pair<double, double> > & ranges);
    
private:
    std::auto_ptr<theta::Minimizer> last_minimizer;
    boost::ptr_vector<theta::Minimizer> minimizers;
};

#endif
