#include "plugins/direct_data_histo.hpp"

direct_data_histo::direct_data_histo(const theta::plugin::Configuration & cfg){
   unsigned int nbins = cfg.setting["nbins"];
   double xmin = cfg.setting["range"][0];
   double xmax = cfg.setting["range"][1];
   theta::Histogram h(nbins, xmin, xmax);
   if(cfg.setting["data"].size() != nbins){
      throw theta::ConfigurationException("The length of " + cfg.setting["data"].getPath() + " and the nbins setting are inconsistent.");
   }
   for(unsigned int i=0; i<nbins; ++i){
      h.set(i+1, cfg.setting["data"][i]);
   }
   set_histo(h);
}

REGISTER_PLUGIN(direct_data_histo)
