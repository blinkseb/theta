#include "interface/producer.hpp"
#include "interface/model.hpp"
#include "interface/phys.hpp"

using namespace theta;
using namespace theta::plugin;

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

ProductsSource::ProductsSource(const plugin::Configuration & cfg): name(cfg.setting["name"]){
    products_sink = cfg.pm->get<ProductsSink>();
    if(not nameOk(name)){
       throw InvalidArgumentException("name '" + name + "' is not a valid product name. ");
   }
}


Producer::Producer(const Configuration & cfg): ProductsSource(cfg){
    if(cfg.setting.exists("override-parameter-distribution")){
        override_parameter_distribution = PluginManager<Distribution>::instance().build(Configuration(cfg, cfg.setting["override-parameter-distribution"]));
    }
    if(cfg.setting.exists("additional-nll-term")){
        additional_nll_term = PluginManager<Function>::instance().build(Configuration(cfg, cfg.setting["additional-nll-term"]));
    }
}

std::auto_ptr<NLLikelihood> Producer::get_nllikelihood(const Data & data, const Model & model){
    std::auto_ptr<NLLikelihood> nll = model.getNLLikelihood(data);
    if(override_parameter_distribution){
        if(!(override_parameter_distribution->getParameters()==model.getParameters())){
            throw FatalException("producer " + getName() + ": override parameter distribution does not define the model parameters");
        }
        nll->set_override_distribution(override_parameter_distribution);
    }
    nll->set_additional_term(additional_nll_term);
    return nll;
}


