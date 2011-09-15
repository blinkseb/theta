#ifndef PLUGIN_CORE_HPP
#define PLUGIN_CORE_HPP

#include "interface/plugin.hpp"
#include "interface/phys.hpp"

/** \brief A data source using a list of constant Histograms
 *
 * It will always yield the same data: which is given as constant histogram functions (i.e.,
 * HistogramFunction instances which do not depend on any parameter).
 *
 * Configured via a setting group like
 * \code
 * data_source = {
 *   type = "histo_source";
 *   obs1 = { // assuming "obs1" is an observable
 *        type = "root_histogram"; filename = "DATA.root"; histoname = "h_obs1";
 *    }; // or some other constant histogram specification
 *   obs2 = { ... }; // some constant histogram specification to use for observable obs2
 * };
 * \endcode
 *
 * \c type must be "histo_source" to select this DataSource type
 *
 * For all observables, a setting specifying the (constant) HistogramFunction must be specified,
 * with the name of that observable.
 *
 * This DataSource does not save any values to the EventTable.
 */
class histo_source: public theta::DataSource{
public:

    /// Construct from a Configuration; required by the plugin system
    histo_source(const theta::plugin::Configuration & cfg);

    /** \brief Fills the provided Data instance with data from the model
     *
     * Will only throw the DataUnavailable Exception if the parameter distribution instance
     * (i.e., either the model's instance or the override-parameter-distribution instance) throws an Exception.
     */
    virtual void fill(theta::Data & dat);
    
private:
    theta::Data data;
};

#endif

