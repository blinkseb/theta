//#include "core.hpp"
#include "plugins/interpolating-histogram.hpp"

using namespace std;
using namespace theta;
using namespace libconfig;

const Histogram1DWithUncertainties & interpolating_histo::operator()(const ParValues & values) const {
    h.set_all_values(1.0);
    const size_t n_sys = hplus.size();
    for (size_t isys = 0; isys < n_sys; isys++) {
        const double delta = values.get(vid[isys]);
        const Histogram1D & t_sys = delta > 0 ? hplus[isys] : hminus[isys];
        if (t_sys.get_nbins() == 0)continue;
        h.multiply_with_ratio_exponented(t_sys, h0, fabs(delta));
    }
    h *= h0;
    h_wu.set(h);
    return h_wu;
}

interpolating_histo::interpolating_histo(const Configuration & ctx){
    SettingWrapper psetting = ctx.setting["parameters"];
    boost::shared_ptr<VarIdManager> vm = ctx.pm->get<VarIdManager>();
    //build nominal histogram:
    h0 = getConstantHistogram(ctx, ctx.setting["nominal-histogram"]);
    size_t n = psetting.size();
    //note: allow n==0 to allow the user to disable systematics.
    // In case of unintentional type error (parameters="delta1,delta2";), user will get a warning about
    // the unused delta*-{plus,minus}-histogram blocks anyway ...
    for(size_t i=0; i<n; i++){
        string par_name = psetting[i];
        ParId pid = vm->getParId(par_name);
        par_ids.insert(pid);
        vid.push_back(pid);
        stringstream setting_name;
        //plus:
        setting_name << par_name << "-plus-histogram";
        hplus.push_back(getConstantHistogram(ctx, ctx.setting[setting_name.str()] ));
        //minus:
        setting_name.str("");
        setting_name << par_name << "-minus-histogram";
        hminus.push_back(getConstantHistogram(ctx, ctx.setting[setting_name.str()] ));
    }
    h = h0;
    
    const size_t nsys = hplus.size();
    std::set<ParId> pid_set;
    for(size_t i=0; i<nsys; i++){
        pid_set.insert(vid[i]);
        h0.check_compatibility(hplus[i]);
        h0.check_compatibility(hminus[i]);
    }
    //to make calculation of derivatives easier, we do not allow the same parameter twice.
    if(pid_set.size()!=nsys){
        throw invalid_argument("interpolating_histo::interpolating_histo: having one parameter parametrizing two interpolations is not supported.");
    }
    h = h0;
    h_wu.set(h0);
}

Histogram1D interpolating_histo::getConstantHistogram(const Configuration & cfg, SettingWrapper s){
    std::auto_ptr<HistogramFunction> hf = PluginManager<HistogramFunction>::build(Configuration(cfg, s));
    if(hf->getParameters().size()!=0){
        stringstream ss;
        ss << "Histogram defined in path " << s.getPath() << " is not constant (but has to be).";
        throw invalid_argument(ss.str());
    }
    return (*hf)(ParValues()).get_values_histogram();//copies the internal reference, so this is ok.
}

REGISTER_PLUGIN(interpolating_histo)
