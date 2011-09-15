#ifndef MINIMIZER_HPP
#define MINIMIZER_HPP

#include "interface/variables.hpp"
#include "interface/matrix.hpp"
#include "interface/phys.hpp"

#include <map>

namespace theta{

    /** \brief The result of a minimization process, returned by \c Minimizer::minimize().
     */
    struct MinimizationResult{
        
        /** \brief The function value at the minimum.
         */
        double fval;

        /** \brief The parameter values at the function minimum.
         *
         * Contains exactly the parameters the function to minimize depends
         * on.
         */
        ParValues values;

        /** \brief The errors at the function minimum, in positive direction.
         *
         * Contains all parameters the function to minimize depends
         * on. How this is calculated depends on the minimizer used. In cases where
         * the minimizer provides symmetrical errors, the entries are equal to \c errors_minus.
         *
         * Set to -1 if not provided by the minimizer.
         */
        ParValues errors_plus;

        /** \brief The errors at the function minimum, in negative direction.
         *
         * Contains all parameters the function to minimize depends
         * on. How this is calculated depends on the minimizer used. In cases where
         * the minimizer provides symmetrical errors, the entries are equal to \c errors_plus.
         *
         * Note that while these are the errors in negative direction, the
         * entries are always positive in case it contains valid errors.
         *
         * Set to -1 if not provided by the minimizer.
         */
        ParValues errors_minus;

        /** \brief Contains the error matrix at the minimum.
         *
         * It is quadratic and has values.size() rows.
         * The index convention is such that 0 corresponds to the first ParId in
         * the sorted (=as iterated) list of parameter ids contained in \c values.
         *
         * It is the negative unity matrix of the correct size in case the
         * minimization does not provide errors.
         */
        Matrix covariance;
        
        /// Define explicitely as ParValues::operator= is private
        void operator=(const MinimizationResult& rhs){
            fval = rhs.fval;
            values.set(rhs.values);
            errors_plus.set(rhs.errors_plus);
            errors_minus.set(rhs.errors_minus);
            covariance = rhs.covariance;
        }
    };


    /** \brief Abstract interface to different minimizers.
     *
     * The possible settings are documented at derived classes.
     */
    class Minimizer{
    public:
        
        /// Define us as the base_type for derived classes; required for the plugin system
        typedef Minimizer base_type;

        /// declare destructor virtual as we expect polymorphic access to derived classes
        virtual ~Minimizer(){}

        /** \brief Attempt to minimize the function.
         *
         * The function f is attempted to be minimized, with specified start values, step sizes and
         * ranges.
         *
         * If a serious error occurs during minimization and the minimization fails,
         * a MinimizationException is thrown. The reasons for such a failure are manifold
         * and also depend on the particular minimization algorithm and should be documented
         * in derived classes.
         * 
         * If either step is 0.0 or the range contains only one value, the parameter should be
         * considered as fixed.
         */
        virtual MinimizationResult minimize(const theta::Function & f, const theta::ParValues & start,
                const theta::ParValues & step, const std::map<theta::ParId, std::pair<double, double> > & ranges) = 0;
                
        const VarIdManager & get_vm() const{
            return *vm;
        }

    protected:
        /// Pointer to the relevant VarIdManager instance. Used to control parameter limits
        const boost::shared_ptr<theta::VarIdManager> vm;
        
        /// The configured tolerance. Meaning depends on the derived class
        double tolerance;

        /// Construct Minimizer from a Configuration instance, setting the VarIdManager vm
        Minimizer(const theta::plugin::Configuration & cfg): vm(cfg.vm){}
    };
    
}

#endif
