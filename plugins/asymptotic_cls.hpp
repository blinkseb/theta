#ifndef PLUGINS_ASYMPTOTIC_CLS
#define PLUGINS_ASYMPTOTIC_CLS

#include "interface/decls.hpp"
#include "interface/main.hpp"
#include "interface/variables.hpp"


/** \brief Calculate Asymptotic CLs limits
 * 
 * Calculate expected and observed CLs limits according to the formulas given in
 * http://arxiv.org/abs/1007.1727 with the test statistic for upper limits q_mu-tilde.
 * It uses the likelihood ratio of the asimov dataset for the expected test statistic distribution (called
 * "sigma_A" in the paper).
 *
 * \code
 * main = {
 *   type = "asymptotic_cls";
 *   parameter = "beta_signal";
 *   cl = 0.9; // optional, default is 0.95
 *   model = "@model";
 *   minimizer = { ... };
 * 
 *   //options for observed limits:
 *   data = {...}; // optional
 * 
 *   // options for expected limits:
 *   quantiles_expected = (0.025, 0.16, 0.5, 0.84, 0.975); // optional; this is the default
 *   parameter_value_expected = 1.0; // optional, default is 0.0
 * };
 * \endcode
 * 
 * \c parameter is the parameter name to calculate the asymptotic CLs limits for.
 * 
 * \c cl is the confidence level of the limit to calculate
 * 
 * \c model is the statistical model
 * 
 * \c minimizer is a Minimizer specification for the minimizer to use for the profile likelihood ratio test statistic calculation
 * 
 * \c data is a DataSource specification for the observed data. If missing, no "observed" limit will be calculated.
 * 
 * \c quantiles_expected are the quantiles of the expected limit distribution to calculate. The default setting given above corresponds to computing
 *   the (typical) median, central 1sigma and central 2sigma intervals. All values must be &gt; 0 and &lt; 1. Can be set to an empty list
 *   to suppress calculation of expected limits.
 * 
 * \c parameter_value_expected is the value of the \c parameter to be used for the expected limit calculation. The default of 0.0
 *   means to compute the expected limit without any signal.
 * 
 * The products table will contain the columns "q" and "limit", both of \c typeDouble. "q" is set to the quantile for the expected limit, using 0.0 for observed data.
 * The "limit" is set to the calculated asymptotic limit.
 */
class asymptotic_cls: public theta::Main{
public:
    explicit asymptotic_cls(const theta::Configuration & cfg);
    virtual void run();
    
private:
    theta::ParId parameter;
    double cl;
    std::auto_ptr<theta::Model> model;
    std::auto_ptr<theta::Minimizer> minimizer;
    
    std::auto_ptr<theta::DataSource> data;
    
    std::vector<double> quantiles_expected;
    double parameter_value_expected;
};


#endif
