#include "plugins/mle.hpp"
#include "plugins/asimov_likelihood_widths.hpp"
#include "interface/plugin.hpp"
#include "interface/minimizer.hpp"
#include "interface/histogram.hpp"
#include "interface/distribution.hpp"

#include <sstream>

using namespace theta;
using namespace std;
using namespace libconfig;
using namespace theta::plugin;

namespace{

int get_index(const ParId & pid, const ParIds & pids){
   int result = 0;
   for(ParIds::const_iterator it=pids.begin(); it!=pids.end(); ++it, ++result){
      if((*it)==pid) return result;
   }
   throw NotFoundException("get_index: no such parameter");
}

}

void mle::produce(const theta::Data & data, const theta::Model & model) {
    std::auto_ptr<NLLikelihood> nll = get_nllikelihood(data, model);
    if(not start_step_ranges_init){
        const Distribution & d = nll->get_parameter_distribution();
        DistributionUtils::fillModeSupport(start, ranges, d);
        step.set(asimov_likelihood_widths(model, override_parameter_distribution));
        start_step_ranges_init = true;
    }
    MinimizationResult minres = minimizer->minimize(*nll, start, step, ranges);
    products_sink->set_product(*c_nll, minres.fval);
    for(size_t i=0; i<save_ids.size(); ++i){
        products_sink->set_product(parameter_columns[i], minres.values.get(save_ids[i]));
        products_sink->set_product(error_columns[i], 0.5 * (minres.errors_plus.get(save_ids[i]) + minres.errors_minus.get(save_ids[i])) );
    }
    if(write_covariance){
       const size_t N = save_ids.size();
       Histogram h(N*N, 0, N*N);
       ParIds pars = nll->getParameters();
       for(size_t i=0; i<N; ++i){
           int index_i = get_index(save_ids[i], pars);
           for(size_t j=0; j<N; ++j){
               int index_j = get_index(save_ids[j], pars);
               h.set(i*N + j + 1, minres.covariance(index_i,index_j));
           }
       }
       products_sink->set_product(*c_covariance, h);
    }
    if(write_ks_ts){
        ObsIds obs = data.getObservables();
        Data pred;
        model.get_prediction(pred, minres.values);
        double ks_ts = 0.0;
        for(ObsIds::const_iterator it=obs.begin(); it!=obs.end(); ++it){
            const Histogram & data_o = data[*it];
            const Histogram & pred_o = pred[*it];
            data_o.check_compatibility(pred_o);
            double sum_d=0, sum_p=0;
            for(size_t i=1; i<=data_o.get_nbins(); ++i){
                sum_d += data_o.get(i);
                sum_p += pred_o.get(i);
                ks_ts = max(ks_ts, fabs(sum_d - sum_p));
            }
        }
        products_sink->set_product(*c_ks_ts, ks_ts);
    }
    if(write_bh_ts){
        ObsIds obs = data.getObservables();
        Data pred;
        model.get_prediction(pred, minres.values);
        double bh_ts = 0.0;
        const Histogram & data_o = data[*bh_ts_obsid];
        const Histogram & pred_o = pred[*bh_ts_obsid];
        data_o.check_compatibility(pred_o);
        for(size_t i=1; i<=data_o.get_nbins(); ++i){
            double bump = 0.0;
            for(size_t j=i; j<=data_o.get_nbins(); ++j){
                bump += data_o.get(j) - pred_o.get(j);
                bh_ts = max(bh_ts, bump);
            }
        }
        products_sink->set_product(*c_bh_ts, bh_ts);
    }
}

mle::mle(const theta::plugin::Configuration & cfg): Producer(cfg), start_step_ranges_init(false), write_covariance(false), write_ks_ts(false), write_bh_ts(false){
    SettingWrapper s = cfg.setting;
    minimizer = PluginManager<Minimizer>::instance().build(Configuration(cfg, s["minimizer"]));
    size_t n_parameters = s["parameters"].size();
    for (size_t i = 0; i < n_parameters; i++) {
        string par_name = s["parameters"][i];
        save_ids.push_back(cfg.vm->getParId(par_name));
        parameter_names.push_back(par_name);
    }
    if(s.exists("write_covariance")){
       write_covariance = s["write_covariance"];
    }
    if(s.exists("write_ks_ts")){
       write_ks_ts = s["write_ks_ts"];
    }
    if(s.exists("bh")){
        bh_ts_obsid.reset(new ObsId(cfg.vm->getObsId(s["bh"])));
        write_bh_ts = true;
    }
    c_nll = products_sink->declare_product(*this, "nll", theta::typeDouble);
    for(size_t i=0; i<save_ids.size(); ++i){
        parameter_columns.push_back(products_sink->declare_product(*this, parameter_names[i], theta::typeDouble));
        error_columns.push_back(products_sink->declare_product(*this, parameter_names[i] + "_error", theta::typeDouble));
    }
    if(write_covariance){
       c_covariance = products_sink->declare_product(*this, "covariance", theta::typeHisto);
    }
    if(write_ks_ts){
       c_ks_ts = products_sink->declare_product(*this, "ks_ts", theta::typeDouble);
    }
    if(write_bh_ts){
       c_bh_ts = products_sink->declare_product(*this, "bh_ts", theta::typeDouble);
    }
}

REGISTER_PLUGIN(mle)

