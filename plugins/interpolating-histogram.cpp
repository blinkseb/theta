//#include "core.hpp"
#include "plugins/interpolating-histogram.hpp"

using namespace std;
using namespace theta;
using namespace theta::plugin;
using namespace libconfig;

const Histogram & interpolating_histo::operator()(const ParValues & values) const {
    h.reset_to_1();
    const size_t n_sys = hplus.size();
    for (size_t isys = 0; isys < n_sys; isys++) {
        const double delta = values.get(vid[isys]);
        const Histogram & t_sys = delta > 0 ? hplus[isys] : hminus[isys];
        if (t_sys.get_nbins() == 0)continue;
        h.multiply_with_ratio_exponented(t_sys, h0, fabs(delta));
    }
    h *= h0;
    return h;
}

const Histogram & interpolating_histo::gradient(const ParValues & values, const ParId & pid) const{
    const size_t n_sys = hplus.size();
    bool found = false;
    for (size_t isys = 0; isys < n_sys; isys++) {
        if(vid[isys]!=pid)continue;
        found = true;
        const double delta = values.get(vid[isys]);
        const Histogram & t_sys = delta > 0 ? hplus[isys] : hminus[isys];
        if (t_sys.get_nbins() == 0)continue;
        h.multiply_with_ratio_exponented(t_sys, h0, fabs(delta));
        const size_t nbins = t_sys.get_nbins();
        for(size_t i=1; i<=nbins; ++i){
            if(h0.get(i)>0.0)
                h.set(i, h.get(i) * theta::utils::log(t_sys.get(i) / h0.get(i)));
        }
        break;
    }
    if(not found){
        h.reset();
    }
    else{
        h *= h0;
    }
    return h;
}

interpolating_histo::interpolating_histo(const Configuration & ctx){
    SettingWrapper psetting = ctx.setting["parameters"];
    //build nominal histogram:
    h0 = getConstantHistogram(ctx, ctx.setting["nominal-histogram"]);
    size_t n = psetting.size();
    //note: allow n==0 to allow the user to disable systematics.
    // In case of unintentional type error (parameters="delta1,delta2";), user will get a warning about
    // the unused delta*-{plus,minus}-histogram blocks anyway ...
    for(size_t i=0; i<n; i++){
        string par_name = psetting[i];
        ParId pid = ctx.vm->getParId(par_name);
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
    assert(hplus.size()==hminus.size());
    assert(vid.size()==hminus.size());
    assert(vid.size()==n);
    h = h0;
    
    const size_t nsys = hplus.size();
    std::set<ParId> pid_set;
    for(size_t i=0; i<nsys; i++){
        pid_set.insert(vid[i]);
        h0.check_compatibility(hplus[i]);
        h0.check_compatibility(hminus[i]);
        //set overflow and underflow to zero. Those are not included in the likelihood and in order
        // to get the normalization right, we have to subtract those:
        hplus[i].set(0,0);
        hplus[i].set(hplus[i].get_nbins()+1,0);
        hminus[i].set(0,0);
        hminus[i].set(hminus[i].get_nbins()+1,0);
    }
    //to make calculation of derivatives easier, we do not allow the same parameter twice.
    if(pid_set.size()!=nsys){
        throw InvalidArgumentException("interpolating_histo::interpolating_histo: having one parameter parametrizing two interpolations is not supported.");
    }
    h0.set(0,0);
    h0.set(h0.get_nbins()+1,0);
}

Histogram interpolating_histo::getConstantHistogram(const Configuration & cfg, SettingWrapper s){
    std::auto_ptr<HistogramFunction> hf = PluginManager<HistogramFunction>::instance().build(Configuration(cfg, s));
    if(hf->getParameters().size()!=0){
        stringstream ss;
        ss << "Histogram defined in path " << s.getPath() << " is not constant (but has to be).";
        throw InvalidArgumentException(ss.str());
    }
    return (*hf)(ParValues());//copies the internal reference, so this is ok.
}

REGISTER_PLUGIN(interpolating_histo)
