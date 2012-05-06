#ifndef HISTOGRAM_FUNCTION_HPP
#define HISTOGRAM_FUNCTION_HPP

#include "interface/decls.hpp"
#include "interface/histogram.hpp"
#include "interface/histogram-with-uncertainties.hpp"
#include "interface/variables.hpp"

namespace theta {
    
    template<typename T>
    class functor{
    public:
        virtual void operator()(const T&) const = 0;
        virtual ~functor(){}
    };
    
    template<typename T>
    class copy_to: public functor<T>{
    private:
        T & into;
    public:
        copy_to(T & into_): into(into_){}
        virtual void operator()(const T& t) const {
            into = t;
        }
    };
    
    template<typename T>
    class add_with_coeff_to: public functor<T>{
    private:
        T & h0;
        double coeff;
    public:
        add_with_coeff_to(T & h0_, double coeff_): h0(h0_), coeff(coeff_){}
        virtual void operator()(const T& t) const {
            theta_assert(h0.get_nbins() == t.get_nbins());
            h0.add_with_coeff(coeff, t);
        }
    };


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
        
        //@{
        /** \brief Apply a functor on the resulting Histogram1D / Histogram1DWithUncertainties
         *
         * This construction is an efficient generalization of the more straight-forward approach of providing
         * evaluation operators which directly return Histogram1D or Histogram1DWithUncertainties: the former
         * implementation would require copying the result to the return value, while this implementation can avoid this
         * copy completely, if the functor does not perform such a copy.
         *
         * To just "get the result histogram", do:
         * \code
         * const ParValues & values;
         * const HistogramFunction & hf;
         * ...
         * Histogram1D h;
         * hf.apply_functor(copy_to<Histogram1D>(h), values);
         * \endcode
         *
         * There are two version: with and without bin-by-bin uncertainties.
         */
        virtual void apply_functor(const functor<Histogram1DWithUncertainties> & f, const ParValues & values) const = 0;
        virtual void apply_functor(const functor<Histogram1D> & f, const ParValues & values) const = 0;
        //@}

        /** \brief Returns the parameters which this HistogramFunction depends on.
         */
        const ParIds & get_parameters() const{
            return par_ids;
        }

        /** \brief Get the dimensions of the Histogram (nbins, xmin, xmax) filled by the evaluation operators
         *
         * This function is used as part of the consistency checks to make sure that the Histogram dimensions match; to save
         * time, it is not usually not used during usual likelihood evaluation, etc.
         */
        virtual void get_histogram_dimensions(size_t & nbins, double & xmin, double & xmax) const = 0;

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

        virtual void apply_functor(const functor<Histogram1DWithUncertainties> & f, const ParValues & values) const;
        virtual void apply_functor(const functor<Histogram1D> & f, const ParValues & values) const;
        virtual void get_histogram_dimensions(size_t & nbins, double & xmin, double & xmax) const;

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
        Histogram1DWithUncertainties h_wu;
        Histogram1D h;
    };
 
    /** \brief Build a HistogramFunction according to the given configuration and return the result
     *
     * This assumes that the HistogramFunction specified by cfg does not depend on any parameters.
     * If this is not the case, an invalid_argument exception will be thrown.
     */
    Histogram1DWithUncertainties get_constant_histogram(const Configuration & cfg);
    
    
}

#endif
