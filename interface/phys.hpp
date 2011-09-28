#ifndef PHYS_HPP
#define PHYS_HPP

#include "interface/decls.hpp"
#include "interface/variables.hpp"
#include "interface/histogram.hpp"
#include "interface/producer.hpp"

#include <vector>


namespace theta {
    
    /** \brief A real-valued function which depends on one or more parameters
     *
     * This is the base class for function plugins.
     */
    class Function{
    public:
        /// Define this as the base_type for derived classes; required for the plugin system
        typedef Function base_type;
        
        /** \brief Evaluate the function at the given parameter values.
         *
         * @return The function value at \c v.
         */
        virtual double operator()(const ParValues & v) const = 0;

        /** \brief Evaluate the function, using the parameter values given as array of doubles.
         *
         * This does the same as operator()(const ParValues&), however, it takes a pointer
         * to an array of doubles instead.
         *
         * This function is provided to make it easier to provide an interface to other libraries, which
         * do not need to know about the theta-style of variable handling (i.e.,
         * \c ParId, \c ParValue, \c VarIdManager classes).
         *
         * The translation from this "external" array format into the internal (i.e.,
         * \c ParValues) format is done by simultaneously iterating over the
         * ParIds as returned by getParameters() and the given array \c x.
         * (As iteration over a \c ParIds object always yields the same order of parameters,
         * this mapping is well-defined.)
         */
        double operator()(const double * x) const{
            size_t i=0;
            for(ParIds::const_iterator it=par_ids.begin(); it!=par_ids.end(); ++it, ++i){
                theta_assert(!std::isnan(x[i]));
                pv.set(*it, x[i]);
            }
            return operator()(pv);
        }


        /** \brief Returns the parameters this function depends on
         */
        const ParIds & getParameters() const{
            return par_ids;
        }
        
        /// Declare destructor virtual as polymorphic access to derived classes will happen.
        virtual ~Function(){}
        
    protected:
        /** \brief The parameters this function depends on
         *
         * Has to be set correctly by derived classes
         */
        ParIds par_ids;
        
    private:
        mutable ParValues pv; //saving this class-wide and not in operator()(const double*) saves quiet some time ...
    };
    
    
    /** \brief Contains data for one or more observables
     *  
     * A data object can be constructed:
     * -# "by hand": use the default constructor and set data Histograms for a number of observables
     * -# by using the DataFactory, which parses a configuration setting and typically reads data from root files
     * -# by sampling from a \c Model using the \c samplePseudoData method.
     *
     * After construction, the Data object is typically used to get a likelihood function from a Model.
     *
     * \sa Model::getLikelihood(Data) DataFactory Model::samplePseudoData
     */
    class Data {
    private:
        void fail_get(const ObsId & oid) const;
    public:
        /** \brief Returns all observable ids for which there is data
         * 
         * Returns those observable ids for which there was previously
         * a Histogram1D saved via the operator[].
         */
        ObsIds getObservables() const;

        /** \brief Access the histogram with an observable
         *
         * The const version is usually only used to read a previously set
         * Histogram. If no Histogram is saved for the supplied observable id,
         * a NotFoundException will be thrown from the const version.
         */
        ///@{
        Histogram1D & operator[](const ObsId & id){
            if(id.id >= data.size()) data.resize(id.id + 1);
            return data[id.id];
        }
        const Histogram1D & operator[](const ObsId & id) const{
            if(id.id >= data.size() || data[id.id].get_nbins()==0) fail_get(id);
            return data[id.id];
        }
        //@}
        
        
        /// \brief reset all current Histograms, i.e., set to zero entry
        void reset(){
            std::vector<Histogram1D>::iterator it = data.begin();
            for(; it!=data.end(); ++it){
               it->set_all_values(0.0);
            }
        }

    private:
        std::vector<Histogram1D> data;
    };
    
    /** \brief A data-providing class, can be used as base class in the plugin system
     *
     * DataSource classes are used as part of a run, which, for each pseudo
     * experiment, calls the DataSource::fill function to get the pseudo data.
     */
    class DataSource: public ProductsSource{
    public:
        
        /// Define this as the base_type for derived classes; required for the plugin system
        typedef DataSource base_type;
        
        /** \brief Exception class to be thrown by fill
         */
        class DataUnavailable{};
        
        /** \brief Fill the provided Data instance with data
         *
         * This method is called exactly once per event. It sets
         * the histograms of the provided observables and lets everything else unaffected.
         *
         * It can throw a DataUnavailable exception if no more data is available and the
         * run should end.
         */
        virtual void fill(Data & dat) = 0;
        
        /// Declare destructor virtual as polymorphic access to derived classes will happen.
        virtual ~DataSource(){}
        
    protected:
        /// proxy to ProductsTableWriter constructor for derived classes
        DataSource(const theta::Configuration & cfg): ProductsSource(cfg){}
    };
    

}

#endif
