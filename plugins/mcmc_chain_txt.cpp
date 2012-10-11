#include "plugins/mcmc_chain_txt.hpp"
#include "plugins/mcmc.hpp"
#include "interface/plugin.hpp"
#include "interface/model.hpp"
#include "interface/histogram.hpp"
#include "interface/distribution.hpp"

#include <sstream>
#include <fstream>
#include <iomanip>

using namespace theta;
using namespace std;

//the result class for the metropolisHastings routine.
class MCMCChainResult{
    public:
        MCMCChainResult(const vector<string> & parameter_names_, const string & filename): parameter_names(parameter_names_){
            outfile.open(filename.c_str(), ios_base::out | ios_base::trunc);
            if(!outfile.good()){
                throw invalid_argument("mcmc_chain_txt: could not open file '" + filename + "' for writing.");
            }
            // write the "header":
            outfile << "# weight";
            for(size_t i=0; i<parameter_names.size(); ++i){
                outfile << " " << parameter_names[i];
            }
            outfile << endl;
        }
        
        size_t getnpar() const{
            return parameter_names.size();
        }
        
        void fill(const double * x, double, size_t n_){
            outfile << n_;
            for(size_t i=0; i<parameter_names.size(); ++i){
                outfile << " " << x[i];
            }
            outfile << endl;
        }
        
        ~MCMCChainResult(){
            outfile.close();
        }
        
    private:
        vector<string> parameter_names;
        ofstream outfile;
};

void mcmc_chain_txt::produce(const Data & data, const Model & model) {
    std::auto_ptr<NLLikelihood> nll = get_nllikelihood(data, model);
    
    if(!init || (re_init > 0 && itoy % re_init == 0)){
        try{
            sqrt_cov = get_sqrt_cov2(*rnd_gen, model, startvalues, override_parameter_distribution, additional_nll_term);
            parameter_names.clear();
            const ParIds & pars = nll->get_parameters();
            parameter_names.reserve(pars.size());
            for(ParIds::const_iterator pid = pars.begin(); pid!=pars.end(); ++pid){
                parameter_names.push_back(vm->get_name(*pid));
            }
            init = true;
        }
        catch(Exception & ex){
            ex.message = "initialization failed: " + ex.message;
            throw invalid_argument(ex.message);        
        }
    }
    ++itoy;
    
    stringstream ss;
    ss << outfile_prefix << "_" << itoy << ".txt";
    MCMCChainResult result(parameter_names, ss.str());
    metropolisHastings(*nll, result, *rnd_gen, startvalues, sqrt_cov, iterations, burn_in);
}

string mcmc_chain_txt::construct_name(){
    static int i = 0;
    stringstream ss;
    ss << "mcmc_chain_txt" << i;
    ++i;
    return ss.str();
}


mcmc_chain_txt::mcmc_chain_txt(const theta::Configuration & cfg): Producer(cfg, construct_name()), RandomConsumer(cfg, get_name()),
   init(false), itoy(0), vm(cfg.pm->get<VarIdManager>()){
    Setting s = cfg.setting;
    outfile_prefix = static_cast<string>(s["outfile_prefix"]);
    iterations = s["iterations"];
    if(s.exists("burn-in")){
        burn_in = s["burn-in"];
    }
    else{
        burn_in = iterations / 10;
    }
    if(s.exists("re-init")){
        re_init = s["re-init"];
    }
}

REGISTER_PLUGIN(mcmc_chain_txt)
