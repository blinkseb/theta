#include "interface/random-utils.hpp"
#include "interface/random.hpp"
#include "interface/histogram.hpp"
#include "interface/database.hpp"
#include "interface/plugin.hpp"

#include <boost/date_time/local_time/local_time.hpp>
#include <unistd.h>

using namespace theta;

RandomConsumer::RandomConsumer(const theta::Configuration & cfg, const std::string & name): seed(-1) {
   std::auto_ptr<RandomSource> rnd_source;
   std::string source_type = "taus";
   if(cfg.setting.exists("rnd_gen")){
       SettingWrapper s = cfg.setting["rnd_gen"];
       if(s.exists("source_type")){
           source_type = static_cast<std::string>(s["source_type"]);
       }
       if(s.exists("seed")){
          seed = s["seed"];
       }
   }
   if(source_type=="taus"){
      rnd_source.reset(new RandomSourceTaus());
   }
   else if(source_type == "mt"){
      rnd_source.reset(new RandomSourceMersenneTwister());
   }
   else{
      throw ConfigurationException("unknown source_type given for rnd_gen (valid values are 'taus' and 'mt')");
   }
   if(seed == -1){
       using namespace boost::posix_time;
       using namespace boost::gregorian;
       ptime t(microsec_clock::universal_time());
       time_duration td = t - ptime(date(1970, 1, 1));
       seed = td.total_microseconds();
       // to avoid clashes with other RandomConsumers initialized in the same microsecond / clock resolution
       // interval: use also the RandomConsumer's name which should be unique within one theta configuration.
       for(size_t i=0; i<name.size(); ++i){
           seed = seed * 33 + (int)name[i];
       }
       // to avoid clashes in case of batch system usage with jobs starting in the same clock resolution
       // interval with the same configuration (=same name), also use the hostname for the seed:
       char hname[HOST_NAME_MAX + 1];
       gethostname(hname, HOST_NAME_MAX + 1);
       // In case the hostname does not fit into hname, the name is truncated but no error is returned.
       // This should not happen, as we use HOST_NAME_MAX. On the other hand, we do not check for
       // any errors potentially returned by gethostname ...
       hname[HOST_NAME_MAX] = '\0';
       int c;
       size_t i=0;
       while((c = hname[i++])){
           seed = seed * 33 + (int)hname[i];
       }
       
   }
   rnd_gen.reset(new Random(rnd_source));
   int runid = *(cfg.pm->get<int>("runid"));
   rnd_gen->set_seed(seed + runid - 1);
   cfg.pm->get<RndInfoTable>()->append(runid, name, seed);
}

RandomConsumer::~RandomConsumer(){}

void theta::randomize_poisson(DoubleVector & d, Random & rnd){
    const size_t n = d.size();
    double * data = d.getData();
    for(size_t i=0; i<n; ++i){
        double mu = data[i];
        if(mu > 0.){
            data[i] = rnd.poisson(mu);
        }
    }
}
