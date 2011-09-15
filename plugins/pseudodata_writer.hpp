#ifndef PLUGIN_PSEUDODATA_WRITER_HPP
#define PLUGIN_PSEUDODATA_WRITER_HPP

#include "interface/decls.hpp"
#include "interface/database.hpp"
#include "interface/producer.hpp"
#include "interface/variables.hpp"

#include <vector>
#include <string>

/** \brief Producer to save the pseudo data used in the pseudo experiments
 *
 * This is meant mainly for cross checks, visualization and debugging. The producer should
 * not be used in large-scale production as it will save a lot of data.
 *
 * Configuration is done with a settings block like:
 * \code
 * pd = {
 *  type = "pseudodata_writer";
 *  name = "pd";
 *  observables = ("o0", "o1");
 *  write-data = true;
 * }
 * \endcode
 *
 * \c type has always to be "pseudodata_writer" in order to use this producer
 *
 * \c name is a name chosen by the user used to construct unique column names in the result table (this name and two underscores are
 *   prepended to the column names explained below).
 *
 * \c observables is a list of observable names for which to write out the pseudodata information
 *
 * \c write-data is a boolean which specifies whether or not to actually write the data: if \c true, the
 *    data will be saved as a Histogram for each observable. If set to \c false, only
 *    the normalization of the data will be saved.
 *
 * For each observable, the result table will contain a column of the name "n_events_<observable name>" and -- if \c write-data
 * is true -- a column "data_<observable name>" with the histogram.
 */
class pseudodata_writer: public theta::Producer {
public:

    /// \brief Constructor used by the plugin system to build an instance from settings in a configuration file
    pseudodata_writer(const theta::plugin::Configuration & cfg);
    virtual void produce(const theta::Data & data, const theta::Model & model);
    
private:
    boost::shared_ptr<theta::VarIdManager> vm;
    std::vector<theta::ObsId> observables;
    boost::ptr_vector<theta::Column> n_events_columns;
    boost::ptr_vector<theta::Column> data_columns;
    bool write_data;
};

#endif
