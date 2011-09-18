#ifndef RANDOM_UTILS_HPP
#define RANDOM_UTILS_HPP

#include "interface/plugin.hpp"
#include "interface/random.hpp"
#include <string>

namespace theta{

/// \brief Base class for plugins using a random number generator.
class RandomConsumer{
protected:
   /** \brief Constructor to be used by derived classes
    *
    * Will save the random seed in the RndInfoTable of the cfg.pm, if this is set.
    */
   RandomConsumer(const theta::plugin::Configuration & cfg, const std::string & name);
   
   /// random seed used
   int seed;
   
   /// random number generator instance to be used by derived classes
   std::auto_ptr<Random> rnd_gen;
   
   // pseudo copy-constructor for clones
   RandomConsumer(const RandomConsumer & rhs, const PropertyMap & pm, const std::string & name);
};


/** \brief Replace data by a Poisson value
 *
 * Replaces each data value by a Poisson with mean equal to the original value.
 */
void randomize_poisson(DoubleVector & d, Random & rnd);


}


#endif

