#include "plugins/cubiclinear_histomorph.hpp"
#include "interface/random.hpp"

using namespace std;
using namespace theta;

const Histogram1DWithUncertainties & cubiclinear_histomorph::operator()(const ParValues & values) const {
    h_wu = h0;
    const size_t n_sys = hplus_diff.size();
    for (size_t isys = 0; isys < n_sys; isys++) {
        const double delta = values.get(vid[isys]) * parameter_factors[isys];
        if(delta==0.0) continue;
        //linear extrpolation beyond 1 sigma:
        if(fabs(delta) > 1){
            const Histogram1D & t_sys = delta > 0 ? hplus_diff[isys] : hminus_diff[isys];
            h_wu.add_with_coeff(fabs(delta), t_sys);
        }
        else{
            //cubic interpolation:
            diff_total = diff[isys];
            diff_total *= 0.5 * delta;
            diff_total.add_with_coeff(delta * delta - 0.5 * pow(fabs(delta), 3), sum[isys]);
            h_wu += diff_total;
        }
    }
    double h_sum = 0.0;
    for(size_t i=0; i<h_wu.get_nbins(); ++i){
        double val = h_wu.get_value(i);
        if(val < 0.0){
            h_wu.set(i, 0.0, h_wu.get_uncertainty(i));
        }
        else{
            h_sum += val;
        }
    }
    if(normalize_to_nominal && h_sum > 0.0){
       h_wu *= h0_sum / h_sum;
    }
    return h_wu;
}

cubiclinear_histomorph::cubiclinear_histomorph(const Configuration & ctx): normalize_to_nominal(false) {
    boost::shared_ptr<VarIdManager> vm = ctx.pm->get<VarIdManager>();
    //build nominal histogram:
    h0 = getConstantHistogram(ctx, ctx.setting["nominal-histogram"]);
    // to check compatibility and for the construction of jhplus and hminus, use this temporary Histogram1D with no uncertainties:
    Histogram1D h0_tmp = h0.get_values_histogram();
    if(ctx.setting.exists("normalize_to_nominal")){
        normalize_to_nominal = ctx.setting["normalize_to_nominal"];
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
        hplus_diff.push_back(getConstantHistogram(ctx, ctx.setting[setting_name]).get_values_histogram());
        hplus_diff.back().check_compatibility(h0_tmp);
        hplus_diff.back().add_with_coeff(-1.0, h0_tmp);
        //minus:
        setting_name = par_name + "-minus-histogram";
        hminus_diff.push_back(getConstantHistogram(ctx, ctx.setting[setting_name]).get_values_histogram());
        hminus_diff.back().check_compatibility(h0_tmp);
        hminus_diff.back().add_with_coeff(-1.0, h0_tmp);
        
        sum.push_back(hplus_diff.back());
        sum.back() += hminus_diff.back();
        diff.push_back(hplus_diff.back());
        diff.back().add_with_coeff(-1, hminus_diff.back());
        
        if(have_parameter_factors){
            parameter_factors[i] = ctx.setting["parameter_factors"][i];
        }
    }
    h0_sum = 0;
    for(size_t i=0; i<h0.get_nbins(); ++i){
        h0_sum += h0.get_value(i);
    }
}

Histogram1DWithUncertainties cubiclinear_histomorph::getConstantHistogram(const Configuration & cfg, SettingWrapper s){
    std::auto_ptr<HistogramFunction> hf = PluginManager<HistogramFunction>::build(Configuration(cfg, s));
    if(hf->getParameters().size()!=0){
        stringstream ss;
        ss << "Histogram defined in path " << s.getPath() << " is not constant (but has to be).";
        throw invalid_argument(ss.str());
    }
    return (*hf)(ParValues());
}

REGISTER_PLUGIN(cubiclinear_histomorph)
