#include "interface/phys.hpp"

/** \brief Add (weighted) squares of parameters and take square root
 *
 * Example configuration:
 * \code
 * f = {
 *   type = "add_squared";
 *   parameters = ("p1", "p2");
 *   weights = (1.0, 1.2);
 * };
 * \endcode
 *
 * \c parameters is a list of parameters to add.
 * 
 * \c weights is a list of floating point values.
 * 
 * If the parameters are p_i and the weights w_i, the function value is
 * sqrt(w_1 * p_1 ^2 + w_2 * p_2^2 + ...)
 *
 * This Function can be useful to split a single nuisance parameter with width 1.0 arouns 0.0 into "subcomponents": Replace the
 * original parameters by a n of new parameters with the same Gaussian prior.
 * Everywhere where you used the original parameter, you would now use this function
 * which depends on all the n new parameters with weights such that the sum of the squared weights is 1.0.
 */
class add_squared: public theta::Function{
public:
    add_squared(const theta::Configuration & cfg);
    virtual double operator()(const theta::ParValues & v) const;
private:
    std::vector<theta::ParId> v_pids;
    std::vector<double> weights;
};

