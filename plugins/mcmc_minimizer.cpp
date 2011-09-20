#include "plugins/mcmc_minimizer.hpp"
#include "plugins/mcmc.hpp"
#include "plugins/mcmc-result.hpp"
#include "interface/plugin.hpp"
#include "interface/model.hpp"
#include "interface/histogram.hpp"
#include "interface/distribution.hpp"

using namespace theta;
using namespace std;

//the result class for the metropolisHastings routine.
class MCMCMinResult: public Result{
    public:
        //ipar_ is the parameter of interest
        MCMCMinResult(size_t npar, const ParIds & pids_): Result(npar), pids(pids_), min_nll(numeric_limits<double>::infinity()){
        }
                
        virtual void fill2(const double * x, double nll, size_t n_){
            if(nll < min_nll){
                min_nll = nll;
                values_at_min.set(ParValues(x, pids));
            }
        }
        
        const ParValues & values_at_minimum() const{
            return values_at_min;
        }
        
        double get_min_nll() const{
            return min_nll;
        }
        
    private:
        ParIds pids;
        double min_nll;
        ParValues values_at_min;
};


MinimizationResult mcmc_minimizer::minimize(const Function & f, const ParValues & start,
                const ParValues & step, const map<ParId, pair<double, double> > & ranges){
    const ParIds & pids = f.getParameters();
    size_t n = pids.size();
    Matrix sqrt_cov(n,n);
    ParIds::const_iterator it=pids.begin();
    vector<double> v_start(n);
    for(size_t i=0; i<n; ++it, ++i){
        sqrt_cov(i,i) = step.get(*it) * stepsize_factor;
        v_start[i] = start.get(*it);
    }
    MCMCMinResult result(n, pids);
    metropolisHastings(f, result, *rnd_gen, v_start, sqrt_cov, iterations, burn_in);
    //now step 2: call the after_minimizer:
    const ParValues & start2 = result.values_at_minimum();
    ParValues step2;
    Matrix cov = result.getCov();
    it=pids.begin();
    for(size_t i=0; i<n; ++it, ++i){
        step2.set(*it, sqrt(cov(i, i)));
    }
    if(after_minimizer.get()){
        return after_minimizer->minimize(f, start2, step2, ranges);
    }
    else{
        MinimizationResult res;
        res.fval = result.get_min_nll();
        res.values.set(start2);
        res.errors_minus.set(step2);
        res.errors_plus.set(step2);
        res.covariance = cov;
        return res;
    }
}

std::auto_ptr<theta::Minimizer> mcmc_minimizer::clone(const theta::PropertyMap & pm) const{
    //throw InvalidArgumentException("mcmc_minimizer does not support clone");
    return std::auto_ptr<theta::Minimizer>(new mcmc_minimizer(*this, pm));
}

mcmc_minimizer::mcmc_minimizer(const mcmc_minimizer & rhs, const theta::PropertyMap & pm): RandomConsumer(rhs, pm, rhs.name), name(rhs.name),
        iterations(rhs.iterations), burn_in(rhs.burn_in), stepsize_factor(rhs.stepsize_factor){
    after_minimizer = rhs.after_minimizer->clone(pm);
}

mcmc_minimizer::mcmc_minimizer(const theta::Configuration & cfg): RandomConsumer(cfg, cfg.setting["name"]), name(cfg.setting["name"]),
  stepsize_factor(1.0){
    SettingWrapper s = cfg.setting;
    iterations = s["iterations"];
    if(s.exists("burn-in")){
        burn_in = s["burn-in"];
    }
    else{
        burn_in = iterations / 10;
    }
    if(s.exists("after_minimizer")){
        after_minimizer = PluginManager<Minimizer>::instance().build(Configuration(cfg, s["after_minimizer"]));
    }
    if(s.exists("stepsize_factor")){
        stepsize_factor = s["stepsize_factor"];
    }
}

REGISTER_PLUGIN(mcmc_minimizer)

