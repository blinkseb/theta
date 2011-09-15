#include "plugins/mcmc_quantiles.hpp"
#include "plugins/mcmc.hpp"
#include "interface/plugin.hpp"
#include "interface/model.hpp"
#include "interface/histogram.hpp"
#include "interface/distribution.hpp"

#include <sstream>
#include <iomanip>

using namespace theta;
using namespace std;
using namespace libconfig;

//the result class for the metropolisHastings routine.
class MCMCPosteriorQuantilesResult{
    public:
        //ipar_ is the parameter of interest
        MCMCPosteriorQuantilesResult(size_t npar_, size_t ipar_, size_t n_iterations_): npar(npar_), ipar(ipar_), n_iterations(n_iterations_){
            par_values.reserve(n_iterations);
        }
        
        size_t getnpar() const{
            return npar;
        }
        
        //just save the parameter value we are interested in
        void fill(const double * x, double, size_t n_){
            for(size_t i=0; i<n_; ++i){
               par_values.push_back(x[ipar]);
            }
        }
        
        //return the quantile q
        double get_quantile(double q){
            if(par_values.size()!=n_iterations){
                throw InvalidArgumentException("MCMCPosteriorQuantilesResult: called get_quantile before chain has finished!");
            }
            int index = static_cast<int>(q * n_iterations);
            if(index >= static_cast<int>(n_iterations)) index = n_iterations-1;
            std::nth_element(par_values.begin(), par_values.begin() + index, par_values.end());
            return par_values[index];
        }
        
    private:
        size_t npar;
        size_t ipar;
        size_t n_iterations;
        vector<double> par_values;
};

void mcmc_quantiles::produce(const Data & data, const Model & model) {
    if(!init){
        try{
            sqrt_cov = get_sqrt_cov2(*rnd_gen, model, startvalues, override_parameter_distribution, vm);
            //find the number of the parameter of interest:
            ParIds model_pars = model.getParameters();
            ipar=0;
            for(ParIds::const_iterator it=model_pars.begin(); it!=model_pars.end(); ++it, ++ipar){
                if(*it == par_id) break;
            }
            //now ipar has the correct value ...
            init = true;
        }
        catch(Exception & ex){
            ex.message = "initialization failed: " + ex.message;
            throw FatalException(ex);
        }
    }
    
    std::auto_ptr<NLLikelihood> nll = get_nllikelihood(data, model);
    MCMCPosteriorQuantilesResult result(nll->getnpar(), ipar, iterations);
    metropolisHastings(*nll, result, *rnd_gen, startvalues, sqrt_cov, iterations, burn_in);
    
    for(size_t i=0; i<quantiles.size(); ++i){
        products_sink->set_product(columns[i], result.get_quantile(quantiles[i]));
    }
}

mcmc_quantiles::mcmc_quantiles(const theta::plugin::Configuration & cfg): Producer(cfg), RandomConsumer(cfg, getName()),
   init(false), par_id(cfg.vm->getParId(cfg.setting["parameter"])){
    vm = cfg.vm;
    SettingWrapper s = cfg.setting;
    string parameter = s["parameter"];
    size_t n = s["quantiles"].size();
    if(n==0){
        throw ConfigurationException("mcmc_quantiles: list of requested quantiles is empty");
    }
    quantiles.reserve(n);
    for(size_t i=0; i<n; ++i){
        quantiles.push_back(s["quantiles"][i]);
        if(quantiles[i]<=0.0 || quantiles[i]>=1.0){
            throw ConfigurationException("mcmc_quantiles: quantiles out of range (must be strictly between 0 and 1)");
        }
    }
    iterations = s["iterations"];
    if(s.exists("burn-in")){
        burn_in = s["burn-in"];
    }
    else{
        burn_in = iterations / 10;
    }
    for(size_t i=0; i<quantiles.size(); ++i){
        stringstream ss;
        ss << "quant" << setw(5) << setfill('0') << static_cast<int>(quantiles[i] * 10000 + 0.5);
        columns.push_back(products_sink->declare_product(*this, ss.str(), theta::typeDouble));
    }
}

REGISTER_PLUGIN(mcmc_quantiles)

