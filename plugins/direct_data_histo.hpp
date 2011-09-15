#ifndef DIRECT_DATA_HISTO_HPP
#define DIRECT_DATA_HISTO_HPP

#include "interface/histogram-function.hpp"
#include "interface/plugin.hpp"


/** \brief A constant Histogram by specifying all bin contents directly in the configuration file
 *
 * Configuration:
 * \code
 *   histogram = {
 *      type = "direct_data_histo";
 *      range = [0.0, 10.0];
 *      nbins = 9;
 *      data = (1.0, 2.0, 3.0, 2.0, 10.0, 2.0, 0.01, 0.01, 2.0);
 *   };
 * \endcode
 *
 * \c range specifies the range used to construct the histogram
 *
 * \c nbins is the number of bins within the range (i.e., not counting overflow and underflow bin)
 *
 * \c data specifies the bin content. This list/array must contain exactly nbins entries; note that underflow and overflow bins
 *  are set to 0.
 */
class direct_data_histo : public theta::ConstantHistogramFunction {
public:
    
    /** \brief Constructor used by the plugin system to build an instance from settings in a configuration file
     */
    direct_data_histo(const theta::plugin::Configuration & ctx);
};

#endif
