#include "plugins/cubiclinear_histomorph.hpp"
#include "interface/random.hpp"

using namespace std;
using namespace theta;

const Histogram1D & cubiclinear_histomorph::operator()(const ParValues & values) const {
    h = h0;
    const size_t n_sys = hplus_diff.size();
    for (size_t isys = 0; isys < n_sys; isys++) {
        const double delta = values.get(vid[isys]) * parameter_factors[isys];
        if(delta==0.0) continue;
        //linear extrpolation beyond 1 sigma:
        if(fabs(delta) > 1){
            const Histogram1D & t_sys = delta > 0 ? hplus_diff[isys] : hminus_diff[isys];
            h.add_with_coeff(fabs(delta), t_sys);
        }
        else{
            //cubic interpolation:
            diff_total = diff[isys];
            diff_total *= 0.5 * delta;
            diff_total.add_with_coeff(delta * delta - 0.5 * pow(fabs(delta), 3), sum[isys]);
            h += diff_total;
        }
    }
    for(size_t i=0; i<h.get_nbins(); ++i){
       h.set(i, max(h.get(i), 0.0));
    }
    if(normalize_to_nominal){
       h *= h0_sum / h.get_sum();
    }
    return h;
}


const Histogram1D & cubiclinear_histomorph::getRandomFluctuation(Random & rnd, const ParValues & values) const{
    if(h0err.get_nbins()==0 && decorr_delta_width==0.0) return operator()(values);
    h = h0;
    // MC stat. uncertainty (for nominal only) if h0err filled:
    if(h0err.size()){
        for(size_t i=0; i<h.get_nbins(); ++i){
            double err_i = h0err.get(i);
            h.set(i, rnd.gauss(err_i) + h.get(i));
        }
    }
    const size_t n_sys = hplus_diff.size();
    double par_factor = decorr_delta_width >=1.0 ? 0.0 : sqrt(1.0 - decorr_delta_width * decorr_delta_width);
    for (size_t isys = 0; isys < n_sys; isys++) {
        const double delta0 = values.get(vid[isys]) * par_factor; 
        for(size_t ibin=0; ibin < h.get_nbins(); ++ibin){
            double delta = (delta0 + rnd.gauss(decorr_delta_width)) * parameter_factors[isys];
            //see above in operator(); same here, just for each bin with a different delta ...
            if(fabs(delta) > 1){
                const Histogram1D & t_sys = delta > 0 ? hplus_diff[isys] : hminus_diff[isys];
                h.set(ibin, h.get(ibin) + fabs(delta) * t_sys.get(ibin));
            }
            else{
                double d = diff[isys].get(ibin);
                d *= 0.5 * delta;
                d += (delta * delta - 0.5 * pow(fabs(delta), 3)) * sum[isys].get(ibin);
                h.set(ibin, h.get(ibin) + d);
            }
        }
    }
    //same as above: truncate at zero and renormalize:
    for(size_t i=0; i<h.get_nbins(); ++i){
       h.set(i, max(h.get(i), 0.0));
    }
    if(normalize_to_nominal){
       h *= h0_sum / h.get_sum();
    }
    return h;
}

cubiclinear_histomorph::cubiclinear_histomorph(const Configuration & ctx): normalize_to_nominal(false), decorr_delta_width(0.0) {
    boost::shared_ptr<VarIdManager> vm = ctx.pm->get<VarIdManager>();
    //build nominal histogram:
    h0 = getConstantHistogram(ctx, ctx.setting["nominal-histogram"]);
    if(ctx.setting.exists("normalize_to_nominal")){
        normalize_to_nominal = ctx.setting["normalize_to_nominal"];
    }
    //h0err contains the absolute uncertainties of each bin
    if(ctx.setting.exists("nominal-uncertainty-histogram")){
        h0err = getConstantHistogram(ctx, ctx.setting["nominal-uncertainty-histogram"]);
        h0err.check_compatibility(h0);
    }
    //IMPORTANT: decorr_delta_width assumes that parameters are distributed with a width of 1.0:
    // the parameter value is scaled down by sqrt(1 - decorr_width**2) and a random number with width decorr_width is added
    // to recover the fluctuation in this bin. Then, as LAST step, the parameter_factor is applied.
    if(ctx.setting.exists("decorr_delta_width")){
        decorr_delta_width = ctx.setting["decorr_delta_width"];
        if(decorr_delta_width < 0.0 || decorr_delta_width > 1.0) throw ConfigurationException("decorr_delta_width must be between 0.0 and 1.0");
    }
    SettingWrapper psetting = ctx.setting["parameters"];
    size_t n = psetting.size();
    parameter_factors.resize(n, 1.0);
    bool have_parameter_factors = ctx.setting.exists("parameter_factors");
    for(size_t i=0; i<n; i++){
        string par_name = psetting[i];
        ParId pid = vm->getParId(par_name);
        par_ids.insert(pid);
        vid.push_back(pid);
        string setting_name;
        //plus:
        setting_name = par_name + "-plus-histogram";
        hplus_diff.push_back(getConstantHistogram(ctx, ctx.setting[setting_name]));
        hplus_diff.back().check_compatibility(h0);
        hplus_diff.back().add_with_coeff(-1.0, h0);
        //minus:
        setting_name = par_name + "-minus-histogram";
        hminus_diff.push_back(getConstantHistogram(ctx, ctx.setting[setting_name]));
        hminus_diff.back().check_compatibility(h0);
        hminus_diff.back().add_with_coeff(-1.0, h0);
        
        sum.push_back(hplus_diff.back());
        sum.back() += hminus_diff.back();
        diff.push_back(hplus_diff.back());
        diff.back().add_with_coeff(-1, hminus_diff.back());
        
        if(have_parameter_factors){
            parameter_factors[i] = ctx.setting["parameter_factors"][i];
        }
    }
    h = h0;
    h0_sum = h0.get_sum();
}

Histogram1D cubiclinear_histomorph::getConstantHistogram(const Configuration & cfg, SettingWrapper s){
    std::auto_ptr<HistogramFunction> hf = PluginManager<HistogramFunction>::build(Configuration(cfg, s));
    if(hf->getParameters().size()!=0){
        stringstream ss;
        ss << "Histogram defined in path " << s.getPath() << " is not constant (but has to be).";
        throw invalid_argument(ss.str());
    }
    return (*hf)(ParValues());
}

REGISTER_PLUGIN(cubiclinear_histomorph)
