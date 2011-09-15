#include "plugins/cubiclinear_histomorph.hpp"

using namespace std;
using namespace theta;
using namespace theta::plugin;

const Histogram & cubiclinear_histomorph::operator()(const ParValues & values) const {
    h = h0;
    const size_t n_sys = hplus_diff.size();
    for (size_t isys = 0; isys < n_sys; isys++) {
        const double delta = values.get(vid[isys]);
        if(delta==0.0) continue;
        //linear extrpolation beyond 1 sigma:
        if(fabs(delta) > 1){
            const Histogram & t_sys = delta > 0 ? hplus_diff[isys] : hminus_diff[isys];
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
    for(size_t i=1; i<=h.get_nbins(); ++i){
       h.set(i, max(h.get(i), 0.0));
    }
    h.set(0,0);
    h.set(h.get_nbins() + 1,0);
    return h;
}

cubiclinear_histomorph::cubiclinear_histomorph(const Configuration & ctx){
    SettingWrapper psetting = ctx.setting["parameters"];
    //build nominal histogram:
    h0 = getConstantHistogram(ctx, ctx.setting["nominal-histogram"]);
    size_t n = psetting.size();
    for(size_t i=0; i<n; i++){
        string par_name = psetting[i];
        ParId pid = ctx.vm->getParId(par_name);
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
    }
    assert(hplus_diff.size()==hminus_diff.size());
    assert(vid.size()==hminus_diff.size());
    assert(vid.size()==n);
    h = h0;
}

Histogram cubiclinear_histomorph::getConstantHistogram(const Configuration & cfg, SettingWrapper s){
    std::auto_ptr<HistogramFunction> hf = PluginManager<HistogramFunction>::instance().build(Configuration(cfg, s));
    if(hf->getParameters().size()!=0){
        stringstream ss;
        ss << "Histogram defined in path " << s.getPath() << " is not constant (but has to be).";
        throw InvalidArgumentException(ss.str());
    }
    return (*hf)(ParValues());
}

REGISTER_PLUGIN(cubiclinear_histomorph)
