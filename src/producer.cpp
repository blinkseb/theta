#include "interface/producer.hpp"
#include "interface/model.hpp"
#include "interface/plugin.tcc"
#include "interface/phys.hpp"
#include "interface/distribution.hpp"
#include "interface/variables-utils.hpp"

#include <sstream>

using namespace theta;

REGISTER_PLUGIN_BASETYPE(Producer);

namespace{
bool nameOk(const std::string & name){
    if(name=="") return false;
    for(size_t i=0; i<name.size(); ++i){
        char c = name[i];
        if((c >= 'A') && (c <= 'Z')) continue;
        if((c >= 'a') && (c <= 'z')) continue;
        if((c >= '0') && (c <= '9')) continue;
        if(c=='_') continue;
        return false;
    }
    return true;
}

}


const std::string & ProductsSource::getName() const{
   return name;
}

ProductsSource::ProductsSource(const std::string & name_, const boost::shared_ptr<ProductsSink> & sink): name(name_), products_sink(sink){}

ProductsSource::ProductsSource(const Configuration & cfg): name(cfg.setting["name"]){
    products_sink = cfg.pm->get<ProductsSink>();
    if(not nameOk(name)){
       throw ConfigurationException("name '" + name + "' is not a valid product name");
   }
}


Producer::Producer(const Configuration & cfg): ProductsSource(cfg){
    if(cfg.setting.exists("override-parameter-distribution")){
        override_parameter_distribution = PluginManager<Distribution>::build(Configuration(cfg, cfg.setting["override-parameter-distribution"]));
    }
    if(cfg.setting.exists("additional-nll-term")){
        additional_nll_term = PluginManager<Function>::build(Configuration(cfg, cfg.setting["additional-nll-term"]));
    }
}

Producer::~Producer(){}

std::auto_ptr<NLLikelihood> Producer::get_nllikelihood(const Data & data, const Model & model){
    std::auto_ptr<NLLikelihood> nll = model.getNLLikelihood(data);
    if(override_parameter_distribution){
        ParIds pars = model.getParameters();
        if(additional_nll_term.get()){
            const ParIds & additional_pars  = additional_nll_term->getParameters();
            pars.insert(additional_pars.begin(), additional_pars.end());
        }
        if(!(override_parameter_distribution->getParameters()==pars)){
            std::stringstream ss;
            ss << "Producer " + getName() + ": override-parameter-distribution must define exactly the parameter models and those of"
                  " 'additional-nll-term' (if present).";
            throw std::invalid_argument(ss.str());
        }
        nll->set_override_distribution(override_parameter_distribution);
    }
    nll->set_additional_term(additional_nll_term);
    return nll;
}

ParameterDependentProducer::~ParameterDependentProducer(){}

