#ifndef DISTRIBUTION_HPP
#define DISTRIBUTION_HPP

#include "interface/decls.hpp"
#include "interface/variables.hpp"
#include "interface/plugin.hpp"

namespace theta{

    /** \brief A probability distribution of parameters in one or more dimensions
     *
     * The distribution class provides methods for generating random numbers according to the distribution
     * it is representing. Further, it provides the negative logarithm of the probability density, including
     * derivatives at a given point.
     *
     * The intended use of this class is twofold:
     * <ol>
     *   <li>As priors in the pseudo data generation. Here, only the sample method is important.</li>
     *   <li>As additional terms in the likelihood function. Here, the other methods
     *       (mainly mpv, width, evalNL, support) are used.</li>
     * </ol>
     *
     * This de-coupling allows some seemingly inconsistent definitions of a Distribution which, for example,
     * a Distribution always sampling the same value but has a non-trivial density. Another example would
     * be an improper density in one and (for example, a flat density on infinite range).
     */
    class Distribution{
    public:
        /// Define us as the base_type for derived classes; required for the plugin system
        typedef Distribution base_type;
        
        /** Sample values from this distribution using \c rnd as random number generator
         *  and respecting limits set for the parameters in \c vm.
         *
         * The ParValues instance \c result will be used to set the random values.
         * This will set values for all variables this distribution
         * is defined on (i.e. all returned by getVariables()). Other parameter
         * values will not be touched.
         *
         * \param[out] result Fill the sampled values here.
         * \param rnd Proxy to the random number generator to use for sampling.
         */
        virtual void sample(ParValues & result, Random & rnd) const = 0;

        /** \brief Provides the mode (most probable values)
         *
         * All parameters returned by getParameters() are set to their most
         * probable value.
         *
         * Derived classes must ensure that if calling evalNL with all parameter values
         * set to their mode as returned by this function has a non-zero probability
         * (i.e., a non-infinite evalNL result).
         *
         * This function is mainly used to select valid start values for
         * algorithms like minimizations or markov chains.
         */
        virtual void mode(ParValues & result) const = 0;

        /** \brief The negative logarithm of the probability density.
         * 
         * The density is not guaranteed to be normalized, i.e., the negative
         * logarithm can be shifted by an arbitrary (but constant) value.
         *
         * \c values must contain (at least) the variable ids in getParameters().
         * Otherwise, a NotFoundException is thrown.
         * 
         * \param values The point in parameter space the density is to be evaluated.
         * \return The density at \c values.
         */
        virtual double evalNL(const ParValues & values) const = 0;
        
       /** \brief The negative logarithm of the probability, and derivatives thereof.
        * 
        * Returns the negative logarithm of the probability density at \c values, just as \c evalNL.
        * Additionally, the partial derivatives are filled into \c derivatives.
        * 
        * \param values The point in parameter space to use to evaluate the density and its derivatives.
        * \param[out] derivatives The container that will be filled with the partial derivatives.
        * \return The density at \c values. This is the same value as \c evalNL(values) would return.
        */
        virtual double evalNL_withDerivatives(const ParValues & values, ParValues & derivatives) const = 0;

        /** \brief Get the support of a parameter
         *
         * The support is the set of values on which the density is non-vanishing. If the support is not
         * an interval, this method should return the smallest interval containing the support.
         *
         * This is mainly used to set constraints for that parameter in a minimization procedure.
         *
         * If \c p is not in getParameters(), the behaviour is undefined (i.e., derived one-dimensional classes
         *  need not check whether p is the correct ParId).
         */
        virtual const std::pair<double, double> & support(const ParId & p) const = 0;

        /** \brief Get the parameters this Distribution depends on and provides values for
         */
        const ParIds & getParameters() const{
            return par_ids;
        }
        
        /// declare destructor as virtual, as polymorphic access will happen
        virtual ~Distribution(){};
    protected:
        ParIds par_ids;
    };
    
    /// \brief namespace for free functions closely related to the \link Distribution Distribution\endlink class
    namespace DistributionUtils{
        
        /** \brief Fill mode and and support from a Distribution instance
         *
         * This is a utility routine calling the Distribution::mode and
         * and Distribution::support routines for all parameters of the Distribution and filling
         * the result into the parameters \c mode, \c width and \c support
         */
        void fillModeSupport(theta::ParValues & mode, std::map<theta::ParId, std::pair<double, double> > & support, const theta::Distribution & d);
    }
    
}

#endif
