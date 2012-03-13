#ifndef HISTOGRAM_FUNCTION_HPP
#define HISTOGRAM_FUNCTION_HPP

#include "interface/decls.hpp"
#include "interface/histogram.hpp"
#include "interface/histogram-with-uncertainties.hpp"
#include "interface/variables.hpp"

namespace theta {

    /** \brief A Histogram-valued function which depends on zero or more parameters.
     *
     * This class is used extensively for model building: a physical model is given by specifying
     * the expected observation in one or more observables and this expectation in turn is specified
     * by histograms which depend on the model parameters. As this can be seen as a histogram-valuesd function,
     * the class is called \c HistogramFunction.
     */
    class HistogramFunction{
    public:
        
        /// Define us as the base_type for derived classes; required for the plugin system
        typedef HistogramFunction base_type;
        
        /** \brief Returns the Histogram for the given parameter values.
         *
         * The returned reference is only guaranteed to be valid as long as this HistogramFunction object.
         */
        virtual const Histogram1DWithUncertainties & operator()(const ParValues & values) const = 0;

        /** \brief Returns the parameters which this HistogramFunction depends on.
         */
        const ParIds & getParameters() const{
            return par_ids;
        }

        /** \brief Get a Histogram of the dimensions (nbins, xmin, xmax) also returned by the evaluation operator
         *
         * The content of the returned Histogram does not matter.
         *
         * This function is used as part of the consistency checks to make sure that the Histogram dimensions match; to save
         * time, it is not usually not used during usual likelihood evaluation, etc.
         */
        virtual Histogram1DWithUncertainties get_histogram_dimensions() const = 0;

        /// Declare the destructor virtual as there will be polymorphic access to derived classes
        virtual ~HistogramFunction(){}
        
    protected:
        /// To be filled by derived classes:
        ParIds par_ids;
    };
    

    /** \brief A simple HistogramFunction which always returns the same Histogram, independent of any parameters.
     *
     * It does not implement any kind of error, i.e., getRandomFluctuation() returns always the same Histogram.
     */
    class ConstantHistogramFunction: public HistogramFunction{
    public:

        /** \brief Returns the Histogram \c h set via set_histo
         */
        virtual const Histogram1DWithUncertainties & operator()(const ParValues & values) const;
        
        /// Return a Histogram of the same dimensions as the one returned by operator()
        virtual Histogram1DWithUncertainties get_histogram_dimensions() const;

    protected:
        /** \brief Set the constant Histogram to return
         *
         * This method is meant for derived classes which can use it to set the constant Histogram to
         * be returned by operator()
         */
        void set_histo(const Histogram1DWithUncertainties & h);
        
        /** \brief Default constructor to be used by derived classes
         */
        ConstantHistogramFunction();
     private:
        Histogram1DWithUncertainties h;
    };

}

#endif
