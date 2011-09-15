#include "plugins/mcmc_posterior_ratio.hpp"
#include "plugins/mcmc.hpp"
#include "interface/plugin.hpp"
#include "interface/model.hpp"
#include "interface/histogram.hpp"
#include "interface/distribution.hpp"

using namespace theta;
using namespace std;
using namespace libconfig;

//the result class for the metropolisHastings routine.
class MCMCPosteriorRatioResult{
    public:
        MCMCPosteriorRatioResult(size_t npar_): npar(npar_), min_nll_value(numeric_limits<double>::infinity()), n_total(0){}
        
        size_t getnpar() const{
            return npar;
        }
        
        void fill(const double *, double nll, size_t n_){
            nll_values.push_back(nll);
            if(nll < min_nll_value) min_nll_value = nll;
            n.push_back(n_);
            n_total += n_;
        }
        
        //return the negative logarithm of the average posterior
        double get_nl_average_posterior(){
            double posterior_sum = 0.0;
            //instead of calculating
            // - log (   1/N  * sum_{i=1}^N exp (-nll_i)    )
            // calculate
            // - log (  exp(-min_nll) * 1/N * sum_{i=1}^N  exp(-nll_i + min_nll)   ) = min_nll - log (  1/N * sum_{i=1}^N  exp (-nll_i + min_nll)   )
            // which is the same , but numerically much better
            for(size_t i=0; i<n.size(); ++i){
                posterior_sum += n[i] * exp(min_nll_value - nll_values[i]);
            }
            return min_nll_value - log(posterior_sum / n_total);
        }
        
        double get_jump_rate(){
            return (1.0 * n.size()) / n_total;
        }
    private:
        size_t npar;
        vector<double> nll_values;
        double min_nll_value;
        vector<size_t> n;
        size_t n_total;
};

void mcmc_posterior_ratio::produce(const theta::Data & data, const theta::Model & model) {
    if(!init){
        try{
            sqrt_cov_sb = get_sqrt_cov2(*rnd_gen, model, startvalues_sb, s_plus_b, vm);
            sqrt_cov_b = get_sqrt_cov2(*rnd_gen, model, startvalues_b, b_only, vm);
            init = true;
        }catch(Exception & ex){
            ex.message = "initialization failed: " + ex.message;
            throw FatalException(ex);
        }
    }
    
    std::auto_ptr<NLLikelihood> nll = get_nllikelihood(data, model);
    
    //a. calculate s plus b:
    MCMCPosteriorRatioResult res_sb(nll->getnpar());
    metropolisHastings(*nll, res_sb, *rnd_gen, startvalues_sb, sqrt_cov_sb, iterations, burn_in);
    double nl_posterior_sb = res_sb.get_nl_average_posterior();

    //b. calculate b only:
    MCMCPosteriorRatioResult res_b(nll->getnpar());
    metropolisHastings(*nll, res_b, *rnd_gen, startvalues_b, sqrt_cov_b, iterations, burn_in);
    double nl_posterior_b = res_b.get_nl_average_posterior();

    if(std::isnan(nl_posterior_sb) || std::isnan(nl_posterior_b)){
        throw Exception("average posterior was NAN");
    }
    products_sink->set_product(*c_nl_posterior_sb, nl_posterior_sb);
    products_sink->set_product(*c_nl_posterior_b, nl_posterior_b);
}

mcmc_posterior_ratio::mcmc_posterior_ratio(const theta::plugin::Configuration & cfg): Producer(cfg), RandomConsumer(cfg, getName()), init(false){
    SettingWrapper s = cfg.setting;
    vm = cfg.vm;
    
    s_plus_b = theta::plugin::PluginManager<Distribution>::instance().build(theta::plugin::Configuration(cfg, s["signal-plus-background-distribution"]));
    b_only = theta::plugin::PluginManager<Distribution>::instance().build(theta::plugin::Configuration(cfg, s["background-only-distribution"]));    
    
    iterations = s["iterations"];
    if(s.exists("burn-in")){
        burn_in = s["burn-in"];
    }
    else{
        burn_in = iterations / 10;
    }
    c_nl_posterior_sb = products_sink->declare_product(*this, "nl_posterior_sb", theta::typeDouble);
    c_nl_posterior_b =  products_sink->declare_product(*this, "nl_posterior_b",  theta::typeDouble);
}

REGISTER_PLUGIN(mcmc_posterior_ratio)
