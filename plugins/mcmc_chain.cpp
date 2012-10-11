#include "plugins/mcmc_chain.hpp"
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
class MCMCResultTextfile{
    public:
        MCMCResultTextfile(const vector<string> & parameter_names_, const string & filename): parameter_names(parameter_names_){
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
        
        ~MCMCResultTextfile(){
            outfile.close();
        }
        
    private:
        vector<string> parameter_names;
        ofstream outfile;
};


class MCMCResultDatabase{
    public:
        MCMCResultDatabase(const vector<string> & parameter_names, Table & table_): npar(parameter_names.size()), table(table_){
            c_weight = table.add_column("weight", typeInt);
            for(size_t i=0; i<parameter_names.size(); ++i){
                c_parameters.push_back(table.add_column(parameter_names[i], typeDouble));
            }
        }
        
        size_t getnpar() const{
            return npar;
        }
        
        void fill(const double * x, double, size_t n_){
            Row r;
            r.set_column(c_weight, (int)n_);
            for(size_t i=0; i<c_parameters.size(); ++i){
                r.set_column(c_parameters[i], x[i]);
            }
            table.add_row(r);
        }
        
    private:
        size_t npar;
        Column c_weight;
        vector<Column> c_parameters;
        Table & table;
};

void mcmc_chain::produce(const Data & data, const Model & model) {
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
    
    if(database.get()){
        stringstream ss;
        ss << "chain_" << itoy;
        std::auto_ptr<Table> table = database->create_table(ss.str());
        MCMCResultDatabase result(parameter_names, *table);
        metropolisHastings(*nll, result, *rnd_gen, startvalues, sqrt_cov, iterations, burn_in);
    }
    else{
        stringstream ss;
        ss << outfile_prefix << "_" << itoy << ".txt";
        MCMCResultTextfile result(parameter_names, ss.str());
        metropolisHastings(*nll, result, *rnd_gen, startvalues, sqrt_cov, iterations, burn_in);
    }
}

string mcmc_chain::construct_name(){
    static int i = 0;
    stringstream ss;
    ss << "mcmc_chain_txt" << i;
    ++i;
    return ss.str();
}


mcmc_chain::mcmc_chain(const theta::Configuration & cfg): Producer(cfg, construct_name()), RandomConsumer(cfg, get_name()),
   init(false), itoy(0), vm(cfg.pm->get<VarIdManager>()){
    Setting s = cfg.setting;
    if(s.exists("outfile_prefix")){
        if(s.exists("output_database")){
            throw ConfigurationException("You have to specify excatly one of the settings 'outfile_prefix', 'output_database', but found both");
        }
        outfile_prefix = static_cast<string>(s["outfile_prefix"]);
    }
    else{
        database = PluginManager<theta::Database>::build(Configuration(cfg, s["output_database"]));
    }
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

REGISTER_PLUGIN(mcmc_chain)
