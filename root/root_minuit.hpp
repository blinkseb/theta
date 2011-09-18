#ifndef THETA_ROOT_MINUIT_HPP
#define THETA_ROOT_MINUIT_HPP

#include "interface/exception.hpp"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/IFunction.h"

#include "interface/plugin.hpp"
#include "interface/phys.hpp"
#include "interface/minimizer.hpp"

/** \brief Minimizer using the MINUIT minimizer from root
 *
 * Configuration with a setting like:
 * \code
 * {
 *  type = "root_minuit";
 *
 *  printlevel = 1; // optional. Default is 0
 *  method = "simplex"; //optional. Default is "migrad"
 *  tolerance = 0.001; //optional. Default as in ROOT::Minuit2
 *  max_iterations = 10000; // optional. Default as in ROOT::Minuit2
 *  max_function_calls = 100000; //optional. Default as in ROOT::Minuit2
 *  n_retries = 10; // optional; the default is 2
 * }
 * \endcode
 *
 * \c printlevel is the verbosity level of the minimizer. The default of 0 does not print anything.
 *  Increase this value in case you are debugging a problem and suspect that it has to do with the minimization.
 *  The value is passed to ROOT::Minuit2::Minuit2Minimizer::SetPrintLevel().
 *
 * \c method must be either "simplex" or "migrad". Refer to the MINUIT documentation on details of these methods.
 *
 * \c tolerance is the Tolerance as should be documented in ROOT::Minuit2::Minuit2Minimizer::SetTolerance.
 *  Default is the one used by ROOT::Minuit2::Minuit2Minimizer.
 *
 * \c max_iterationc and \c max_function_calls are the according settings of ROOT::Minuit2::Minuit2Minimizer.
 *
 * \c n_retries is the maximum number of retries in case the first minimization attempt fails. In some cases,
 *    the success rate for the minimization can be increased by re-running the minimization starting at the parameter values
 *    of the current (failed) attempt.
 *
 * Please note that this plugin relies on the Minuit2 implementation of ROOT which is poorly documented. Minuit2
 * is a C++ proxy to the fortran MINUIT for which you can find more documentation.
 */
class root_minuit: public theta::Minimizer{
public:
    /** \brief Constructor used by the plugin system to build an instance from a configuration file.
     */
    root_minuit(const theta::plugin::Configuration & cfg);

    /** \brief Implement the Minimizer::minimize routine.
     *
     * See documentation of Minimizer::minimize for an introduction
     *
     * The minimizer forwards the task to a ROOT::Minuit2::Minuit2Minimizer using its
     * functions also for setting limits on the parameters. If minimization fails, it
     * is attempted up to three times through repeating calls of ROOT::Minuit2::Minuit2Minimizer::Minimize.
     * If still an error is reported, a \ref theta::MinimizationException is thrown which contains the status
     * code returned by ROOT::Minuit2::Minuit2Minimizer::Status().
     */
    virtual theta::MinimizationResult minimize(const theta::Function & f, const theta::ParValues & start,
            const theta::ParValues & step, const std::map<theta::ParId, std::pair<double, double> > & ranges);
    
    virtual std::auto_ptr<theta::Minimizer> clone(const PropertyMap & pm) const;
private:
    
    ROOT::Minuit2::EMinimizerType type;
    double tolerance;
    int printlevel;
    int max_iterations, max_function_calls, n_retries;
};

#endif
