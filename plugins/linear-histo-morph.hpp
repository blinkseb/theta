#ifndef LINEAR_HISTO_MORPH_HPP
#define LINEAR_HISTO_MORPH_HPP

#include "interface/histogram-function.hpp"
#include "interface/plugin.hpp"


/** \brief A HistogramFunction which interpolates a "zero" Histogram and several "distorted" Histograms as generic method to treat systematic uncertainties.
 *
 * The configuration block is something like (mass is some observable previously defined;
 * s, delta1, delta2 are parameters):
 * \code
 *   histogram = {
 *      type = "linear_histo_morph";
 *      parameters = ("delta1", "delta2");
 *      nominal-histogram = { / * fixed-histogram-specification * / };
 *      nominal-expectation = 120.0;
 *      delta1-kappa-plus-histogram = { / * fixed-histogram-specification * / };
 *      delta1-kappa-minus-histogram = { / * fixed-histogram-specification * / };
 *      delta1-minus-relexp = 0.17;  // relative change in expectation
 *      delta1-plus-relexp = -0.13;
 *      delta2-kappa-plus-histogram = { / * fixed-histogram-specification * / };
 *      delta2-kappa-minus-histogram = { / * fixed-histogram-specification * / };
 *   };
 * \endcode
 * Here, <tt>fixed-histogram-specification</tt> is a Histogram Setting block that returns a Histogram
 * which does not depend on any parameters. That is, any Histogram of type="fixed_*".
 *
 * No random fluctuations are done: the method \c getRandomFluctuation is not overriden from \link theta::HistogramFunction HistogramFunction \endlink
 * and will therefore return the same values as the usual, non-random, evaluation operator, \c operator().
 *
 * A histogram \c h0 is interpolated by parameters p_i specified in the \c parameters linearly:
 *  -# if all p_i = 0, the original Histogram h0 is returned
 *  -# for all other values, the histogram contents are interpolated binwise linearly: to get the bin value of bin k, first
 *    the sum over all kappa[i][k] * fabs(delta[i]) is calculated where i runs over all uncertainties (i.e., all parameters
 *    given in the \c parameters setting). Then, this is multiplied binwise with the nominal histogram.
 *  -# The normalization of the result is calculated as product of (1 + abs(parameters[i]) * relexp[i]) where i runs over all parameters.
 *
 * It is valid to give only one or even no kappa histogram. In this case, only the rate uncertainty is applied. If no
 * explicit rate uncertainty is given, 0.0 is assumed.
 */
class linear_histo_morph : public theta::HistogramFunction {
public:
    
    /** \brief Constructor used by the plugin system to build an instance from settings in a configuration file
     */
    linear_histo_morph(const theta::plugin::Configuration & ctx);
        
    /** Returns the interpolated Histogram as documented in the class documentation.
     * throws a NotFoundException if a parameter is missing.
     */
    virtual const theta::Histogram & operator()(const theta::ParValues & values) const;

    /// Return a Histogram of the same dimenions as the one returned by operator()
    virtual theta::Histogram get_histogram_dimensions() const{
        return h;
    }

private:
    /** \brief Build a (constant) Histogram from a Setting block.
    *
    * Will throw an InvalidArgumentException if the Histogram is not constant.
    */
    static theta::Histogram getConstantHistogram(const theta::plugin::Configuration & ctx, theta::SettingWrapper s);
    
    theta::Histogram h0;
    double h0exp;
    
    std::vector<theta::Histogram> kappa_plus;
    std::vector<theta::Histogram> kappa_minus;
    std::vector<double> plus_relexp;
    std::vector<double> minus_relexp;
    
    //the interpolation parameters used to interpolate between hplus and hminus.
    std::vector<theta::ParId> parameters;
    
    //the Histogram returned by operator(). Defined as mutable to allow operator() to be const.
    mutable theta::Histogram h;
};

#endif
