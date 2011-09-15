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
};


/** \brief Dice a poisson for each bin of a given Histogram
 *
 * Replaces the bin contents of the given Histogram \c h with a random variable drawn from a Poisson distribution.
 * As mean, the original bin content is used.
 */
void randomize_poisson(Histogram & h, Random & rnd);


}


#endif

