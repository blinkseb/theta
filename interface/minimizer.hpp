#ifndef MINIMIZER_HPP
#define MINIMIZER_HPP

#include "interface/decls.hpp"
#include "interface/distribution.hpp"
#include "interface/variables.hpp"
#include "interface/matrix.hpp"

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
         * Empty if not provided by the minimizer.
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
         * Empty if not provided by the minimizer.
         */
        ParValues errors_minus;

        /** \brief Contains the error matrix at the minimum.
         *
         * It is quadratic and has values.size() rows.
         * The index convention is such that 0 corresponds to the first ParId in
         * the sorted (=as iterated) list of parameter ids contained in \c values.
         *
         * It is set to a Matrix of size 0 in case the
         * minimization does not provide errors.
         */
        Matrix covariance;
        
        /// Define explicitely as ParValues::operator= is private
        void operator=(const MinimizationResult& rhs);
    };
    
    /** \brief A minimization problem: function, parameter bounds, and start/step values
     * 
     * Successful construction guarantees some basic properties: start, step, and range values are available for
     * all parameters of the function; all step sizes are &gt;= 0.0; the start value is always within the range.
     * 
     * The matrix is an estimate for the covariance matrix, using the function parameters to convert ParIds to indices
     * consistently.
     */
    // this could also be done in parameters of Minimizer::minimize, but this is more general ...
    class MinimizationProblem{
    public:
        MinimizationProblem(const theta::Function & f, const theta::ParValues & start, const theta::ParValues & step, const std::map<theta::ParId, std::pair<double, double> > & ranges);
        MinimizationProblem(const theta::Function & f, const theta::ParValues & start, const theta::Matrix & matrix, const std::map<theta::ParId, std::pair<double, double> > & ranges);
        
        /// The function to minimize
        const theta::Function & get_f() const{
            return f;
        }
        
        /// The start values for the minimization.
        theta::ParValues get_start() const{
            return start;
        }
        
        /// The parameter ranges. Guaranteed to contain all parameters of f.
        const Ranges & get_ranges() const{
            return ranges;
        }
        
        /// the step values contain approximately the standard deviation of the parameters to estimate
        theta::ParValues get_steps() const{
            return step;
        }
        
        /// the step matrix holds approximately the covariance matrix (if known).
        Matrix get_step_matrix() const{
            return matrix;
        }
        
    private:
        const Function & f;
        ParValues start;
        ParValues step;
        Ranges ranges;
        Matrix matrix;
        
        void check_consistency() const;
    };


    /** \brief Abstract interface to different minimizers.
     *
     * The possible settings are documented at derived classes.
     */
    class Minimizer{
    public:
        
        /// Define this class as the base_type for derived classes; required for the plugin system
        typedef Minimizer base_type;

        /// declare destructor virtual as we expect polymorphic access to derived classes
        virtual ~Minimizer();

        /** \brief Find the minimum of a function.
         *
         * The minimum of f is found, starting at the specified start values, step sizes and
         * ranges. The step size should be an estimate of the parameter uncertainty.
         *
         * If a serious error occurs during minimization and the minimization fails,
         * a MinimizationException is thrown. The reasons for such a failure 
         * depend on the particular minimization algorithm and should be documented
         * in derived classes.
         * 
         * If for a parameter either step is 0.0 or the range contains only one value, this parameter is considered
         * as constant and it is not varied during minimization.
         */
        virtual MinimizationResult minimize(const theta::Function & f, const theta::ParValues & start,
                const theta::ParValues & step, const Ranges & ranges) = 0;
        
        /** \brief Alternative minimize method with more information
         *
         * Same as the minimize method, but can make use of more information about the minimization problem
         * via the MinimizationProblem instance.
         * 
         * Derived classes should implement this method. The default implementation is a simple forwarding to
         * the minimize method above. In derived classes, it makes sense to implement the main
         * function in this method and make minimize a forward to this one, using the appropriate constructor
         * of MinimizationProblem.
         */
        virtual MinimizationResult minimize2(const MinimizationProblem & mp);
        
    };
    
}

#endif
