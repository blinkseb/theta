#ifndef PLUGINS_ASIMOV_LIKELIHOOD_WIDTHS_HPP
#define PLUGINS_ASIMOV_LIKELIHOOD_WIDTHS_HPP

#include "interface/decls.hpp"
#include <boost/shared_ptr.hpp>

/** \brief Calculate parameter widths from the asimov data for the given model
 *
 * The model's distribution will be used to determine the most probable values / modes of all parameters.
 * The prediction using these parameter values (without Poisson smearing = Asimov data) is used to construct
 * a nllikelihood function. From this likelihood function, approximate 1sigma parameter widths are extracted by scanning
 * the function in all directions until the function value reaches (nll value at minimum)+1/2. No minimisation is performed
 * in this process.
 * If the "+1/2 point" cannot be reached within the parameter boundaries given by the Distribution support,
 * the width for that parameter is set to the total width of the parameter's support. This can be 0, e.g., for
 * delta functions. It can also be infinity, e.g., for parameters with a flat distribution which do not really affect the model prediction.
 * Callers of this method should make sure to treat such cases properly. The zero case is often valid (=parameter is fixed), whereas
 * the large value / infinity case usually indicates a problem bacause it means the likelihood function (and thus the model prediction)
 * hardly depends on this parameter.
 *
 * If override_parameter_distribution is given, it will be used instead of the model's Distribution for both the Asimov
 * data construction as well as for the nllikelihood definition.
 *
 * This function is meant to be used to choose initial step sizes for minimization and MCMC integration.
 */
theta::ParValues asimov_likelihood_widths(const theta::Model & model, const boost::shared_ptr<theta::Distribution> & override_parameter_distribution,
     const boost::shared_ptr<theta::Function> & additional_nll_term);

#endif
