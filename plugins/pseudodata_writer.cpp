#include "plugins/pseudodata_writer.hpp"
#include "interface/plugin.hpp"
#include "interface/phys.hpp"
#include "interface/histogram.hpp"
#include <sstream>

using namespace theta;
using namespace std;
using namespace libconfig;

void pseudodata_writer::produce(const Data & data, const Model & model) {
    for(size_t i=0; i<observables.size(); ++i){
        const Histogram1D & h = data[observables[i]];
        double n_event = h.get_sum();
        products_sink->set_product(n_events_columns[i], n_event);
        if(write_data){
            products_sink->set_product(data_columns[i], h);
        }
    }
}

pseudodata_writer::pseudodata_writer(const theta::Configuration & cfg): Producer(cfg){
    size_t n = cfg.setting["observables"].size();
    boost::shared_ptr<VarIdManager> vm = cfg.pm->get<VarIdManager>();
    observables.reserve(n);
    for(size_t i=0; i<n; i++){
        observables.push_back(vm->getObsId(cfg.setting["observables"][i]));
    }
    write_data = cfg.setting["write-data"];
    declare_products(vm);
}

void pseudodata_writer::declare_products(const boost::shared_ptr<VarIdManager> & vm){
   for(size_t i=0; i<observables.size(); ++i){
        n_events_columns.push_back(products_sink->declare_product(*this, "n_events_" + vm->getName(observables[i]), theta::typeDouble));
        if(write_data)
            data_columns.push_back(products_sink->declare_product(*this, "data_" + vm->getName(observables[i]), theta::typeHisto));
    }
}

std::auto_ptr<theta::Producer> pseudodata_writer::clone(const theta::PropertyMap & pm) const{
    return std::auto_ptr<theta::Producer>(new pseudodata_writer(*this, pm));
}

pseudodata_writer::pseudodata_writer(const pseudodata_writer & rhs, const theta::PropertyMap & pm): Producer(rhs, pm),
  observables(rhs.observables), write_data(rhs.write_data){
    boost::shared_ptr<VarIdManager> vm = pm.get<VarIdManager>();
    declare_products(vm);
}

REGISTER_PLUGIN(pseudodata_writer)
