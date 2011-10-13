#ifndef PRODUCER_HPP
#define PRODUCER_HPP

#include "interface/decls.hpp"
#include "interface/variables.hpp"
#include "interface/data_type.hpp"

#include <string>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>

namespace theta {

/** \brief Container for the by-event products to be used by ProductsSources
 *
 * Abstract base class which takes products from a ProductsSource.
 */
class ProductsSink: private boost::noncopyable{
public:
    virtual Column declare_product(const ProductsSource & source, const std::string & product_name, const data_type & type) = 0;
    virtual void set_product(const Column & c, double d) = 0;
    virtual void set_product(const Column & c, int i) = 0;
    virtual void set_product(const Column & c, const std::string & s) = 0;
    virtual void set_product(const Column & c, const Histogram1D & h) = 0;
    virtual ~ProductsSink(){}
};


/** \brief Base class for all classes writing products to a ProductsSink
 *
 * products_sink is filled via Configuration cfg.pm.
 *
 * Each ProductsSource has a name in order to identify
 * the products in case of multiple same producers. This name
 * is set explicitely in the configuration file via the 'name' setting.
 */
class ProductsSource{
public:
    /// Get the name as configured via the configuration file
    const std::string & getName() const;
    
protected:
    /// To be used by derived classes, to fill name and products_sink
    ProductsSource(const Configuration & cfg);
    ProductsSource(const std::string & name_, const boost::shared_ptr<ProductsSink> & sink);
    
    std::string name;
    boost::shared_ptr<ProductsSink> products_sink;
};

/** \brief The abstract base class for all statistical methods or other objects creating per-event data
 *
 * It is called "producer" as it produces results, given Data and
 * a Model. Every Producer belongs to exactly one Run, which calls its produce method (typically
 * repeatedly on some pseudo data).
 *
 * Common to all producers are the settings "override-parameter-distribution" and "additional-nll-term".
 * The former is a distribution for all model parameters to be used for the likelihood function instead
 * of the ones from the model. The latter is a Function to add to the negative log-likelihood such as external
 * constraints or priors which cannot be expressed as parameter distributions.
 */
class Producer: public ProductsSource{
public:
    /// Define us as the base_type for derived classes; required for the plugin system
    typedef Producer base_type;
    
    /** Declare the destructor as virtual, as we expect polymorphic
     *  access to derived classes.
     */
    virtual ~Producer();
    
    /** \brief Run a statistical algorithm on the data and model and write out the results
     *
     * The result should be written to \c products_sink by calling ProductsSink::set_product on
     * product columns previously defined via ProductsSink::decalre_product in the constructor.
     *
     * Derived classes may assume that all calls are done with the same \c model.
     *
     * In case of an error, the method should through a class derived from \link theta::Exception \endlink
     * or \link theta::FatalException \endlink.
     */
    virtual void produce(const Data & data, const Model & model) = 0;
    
protected:
    /** \brief Construct from a Configuration instance
     *
     * Parses the settings "add-nllikelihood-functions".
     */
    Producer(const Configuration & cfg);
    
    /** \brief Get the likelihood for the provided Data and Model, including the setting of \c override-parameter-distribution, if applicable
     * 
     * Derived classes should always use this method and never construct the NLLIkelihood
     * directly from a Model instance to ensure consistent treatment of the additional likelihood terms.
     */
    std::auto_ptr<NLLikelihood> get_nllikelihood(const Data & data, const Model & model);
    
    boost::shared_ptr<theta::Distribution> override_parameter_distribution;
    boost::shared_ptr<theta::Function> additional_nll_term;
};

/** \brief Base class for Producers whose result depend on parameter values
 *
 * If a producer depends on the value of certain model parameters, it should inherit from this class
 * instead of inheriting from Producer directly.
 *
 * Derived classes should fill par_ids to indicate which parameters the class depends on (which can also be
 * an empty set).
 */
class ParameterDependentProducer: public Producer{
public:
    virtual ~ParameterDependentProducer();
    virtual void setParameterValues(const ParValues & values) = 0;
    ParIds getParameters() const {
        return par_ids;
    }
protected:
    theta::ParIds par_ids;

    // forward to Producer(cfg)
    ParameterDependentProducer(const Configuration & cfg): Producer(cfg){}
};


}

#endif

