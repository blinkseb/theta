#ifndef PLUGINS_MCMC_QUANTILES_HPP
#define PLUGINS_MCMC_QUANTILES_HPP

#include "interface/decls.hpp"
#include "interface/variables.hpp"
#include "interface/database.hpp"
#include "interface/producer.hpp"
#include "interface/random-utils.hpp"
#include "interface/matrix.hpp"

#include <string>

/** \brief Write out the Markov Chain in text files
 *
 * The result can be used to construct many informations from the chain directly, such as
 * multidimensional posteriors, etc.
 *
 * It will create one text file per run.
 *
 * Configuration is done via a setting group like
 * \code
 * chain = {
 *   type = "mcmc_chain_txt";
 *   outfile_prefix = "chain";
 *   iterations = 10000;
 *   burn-in = 100; //optional. default is iterations / 10
 *   re-init = 1; //optional. Default is 0
 * };
 *
 * \endcode
 *
 * \c type is always "mcmc_chain_txt" to select this producer.
 *
 * \c outfile_prefix is the filename prefix for the produces .txt files. The files will be constructed by the \c outfile_prefix given here, and an incrementing counter,
 *   separated by an underscore ("_"). The filename ends with ".txt".
 *
 * \c iterations is the number of MCMC iterations. See additional comments about runtime and suggested robustness tests
 *     in the documentation of \link mcmc_posterior_ratio mcmc_posterior_ratio \endlink.
 *
 * \c burn-in is the number of MCMC iterations to do at the beginning and throw away. See additional comments in the
 *     documentation of \link mcmc_posterior_ratio mcmc_posterior_ratio \endlink
 *
 * \c re-init is an optional integer which controls re-initialisation of the jumping kernel width. The deault of 0 never re-initialises. For a value N > 0, 
 * re-initialisation is done every N toys.
 *
 * <b>Important:</b> This plugin does not write anything to the output database.
 */
class mcmc_chain_txt: public theta::Producer, public theta::RandomConsumer{
public:
    mcmc_chain_txt(const theta::Configuration & ctx);
    virtual void produce(const theta::Data & data, const theta::Model & model);
    
private:
    
    static std::string construct_name();
    
    //whether sqrt_cov* and startvalues* have been initialized:
    bool init;
    

    int re_init, itoy;
    
    boost::shared_ptr<theta::VarIdManager> vm;
    std::vector<std::string> parameter_names;
    std::string outfile_prefix;
    
    //MCMC parameters:
    unsigned int iterations;
    unsigned int burn_in;
    theta::Matrix sqrt_cov;
    std::vector<double> startvalues;
};

#endif
