#include "plugins/linear-histo-morph.hpp"
#include "interface/utils.hpp"

using namespace std;
using namespace theta;
using namespace theta::plugin;
using namespace libconfig;

const Histogram & linear_histo_morph::operator()(const ParValues & values) const {
    h.reset_to_1();
    const size_t n_sys = kappa_plus.size();
    //1. interpolate linearly in each bin; also calculate normalization
    double scale_unc = 1;
    for (size_t isys = 0; isys < n_sys; ++isys) {
        const double delta = values.get(parameters[isys]);
        if(delta==0.0) continue;
        const Histogram & kappa_sys = delta > 0 ? kappa_plus[isys] : kappa_minus[isys];
        if(kappa_sys.get_nbins() > 0)
           h.add_with_coeff(fabs(delta), kappa_sys);
        double relexp = delta > 0 ? plus_relexp[isys] : minus_relexp[isys];
        // multiply with 1 + |\delta| relexp, but cutoff at 0:
        scale_unc *= max(0.0, 1.0 + fabs(delta) * relexp);
    }
    h *= h0;
    //2. lower bin cutoff
    for(size_t i=1; i <= h.get_nbins(); ++i){
       if(h.get(i) < 0.0){
         h.set(i, 0.0);
         //throw UnphysicalPredictionException();
       }
    }
    //3.a. rescale to nominal:
    h *= h0exp / h.get_sum_of_bincontents();
    //3.b. apply scale uncertainty
    h *= scale_unc;
    return h;
}

linear_histo_morph::linear_histo_morph(const Configuration & ctx){
    SettingWrapper psetting = ctx.setting["parameters"];
    size_t n = psetting.size();
    for(size_t i=0; i<n; i++){
        string par_name = psetting[i];
        ParId pid = ctx.vm->getParId(par_name);
        par_ids.insert(pid);
        parameters.push_back(pid);
        if(ctx.setting.exists(par_name + "-kappa-plus-histogram"))
           kappa_plus.push_back(getConstantHistogram(ctx, ctx.setting[par_name + "-kappa-plus-histogram"] ));
        else
           kappa_plus.push_back(Histogram());
        if(ctx.setting.exists(par_name + "-kappa-minus-histogram"))
            kappa_minus.push_back(getConstantHistogram(ctx, ctx.setting[par_name + "-kappa-minus-histogram"] ));
        else
           kappa_minus.push_back(Histogram());
        if(ctx.setting.exists(par_name + "-plus-relexp"))
            plus_relexp.push_back(ctx.setting[par_name + "-plus-relexp"]);
        else
            plus_relexp.push_back(0.0);
        if(ctx.setting.exists(par_name + "-minus-relexp"))
            minus_relexp.push_back(ctx.setting[par_name + "-minus-relexp"]);
        else
            minus_relexp.push_back(0.0);
    }
    assert(parameters.size()==n);
    assert(n==kappa_minus.size());
    assert(n==kappa_plus.size());
    assert(n==plus_relexp.size());
    assert(n==minus_relexp.size());
    
    const size_t nsys = kappa_plus.size();
    std::set<ParId> pid_set;
    h0 = getConstantHistogram(ctx, ctx.setting["nominal-histogram"]);
    h = h0;
    for(size_t i=0; i < nsys; i++){
        pid_set.insert(parameters[i]);
        if(kappa_plus[i].get_nbins() > 0)
           h0.check_compatibility(kappa_plus[i]);
        if(kappa_minus[i].get_nbins() > 0)
           h0.check_compatibility(kappa_minus[i]);
        kappa_plus[i].set(0,0);
        kappa_plus[i].set(kappa_plus[i].get_nbins()+1,0);
        kappa_minus[i].set(0,0);
        kappa_minus[i].set(kappa_minus[i].get_nbins()+1,0);
    }
    if(pid_set.size()!=nsys){
        throw InvalidArgumentException("linear_histo_morph: duplicate parameter in parameter list.");
    }
    h0.set(0,0);
    h0.set(h0.get_nbins()+1,0);
    h0exp = ctx.setting["nominal-expectation"];
    h0 *= h0exp / h0.get_sum_of_bincontents();
}

Histogram linear_histo_morph::getConstantHistogram(const Configuration & cfg, SettingWrapper s){
    std::auto_ptr<HistogramFunction> hf = PluginManager<HistogramFunction>::instance().build(Configuration(cfg, s));
    if(hf->getParameters().size()!=0){
        stringstream ss;
        ss << "Histogram defined in path " << s.getPath() << " is not constant (but has to be).";
        throw InvalidArgumentException(ss.str());
    }
    return (*hf)(ParValues());
}

REGISTER_PLUGIN(linear_histo_morph)
