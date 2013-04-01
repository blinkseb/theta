#include "plugins/cubiclinear_histomorph.hpp"
#include "interface/random.hpp"

using namespace std;
using namespace theta;

namespace{
    template<typename HT, typename HT2>
    void add_with_coeff2(HT & t, double c1, const HT2 & v1, double c2, const HT2 & v2){
        t.add_with_coeff(c1, v1);
        t.add_with_coeff(c2, v2);
    }

    template<>
    void add_with_coeff2<Histogram1D, Histogram1D>(Histogram1D & t, double c1, const Histogram1D & v1, double c2, const Histogram1D & v2){
        t.add_with_coeff2(c1, v1, c2, v2);
    }
}

template<typename HT>
void cubiclinear_histomorph::add_morph_terms(HT & t, const ParValues & values) const{
    const size_t n_sys = hplus_diff.size();
    for (size_t isys = 0; isys < n_sys; isys++) {
        const double delta = values.get_unchecked(vid[isys]) * parameter_factors[isys];
        if(delta==0.0) continue;
        //linear extrapolation beyond 1 sigma:
        if(fabs(delta) > 1){
            const Histogram1D & t_sys = delta > 0 ? hplus_diff[isys] : hminus_diff[isys];
            t.add_with_coeff(fabs(delta), t_sys);
        }
        else{
            //cubic interpolation:
            const double d2 = delta * delta;
            const double d3 = d2 * fabs(delta);
            add_with_coeff2(t, 0.5*delta, diff[isys], d2 - 0.5 * d3, sum[isys]);
        }
    }
    double h_sum = 0.0;
    const size_t n = t.get_nbins();
#ifndef __SSE2__
    for(size_t i=0; i < n; ++i){
        double val = t.get(i);
        if(val < 0.0){
            t.set(i, 0.0);
        }
        else{
            h_sum += val;
        }
    }
#else
    double * hdata = t.get_data();
    const __m128d zero = _mm_setzero_pd();
    __m128d sum = zero;
    for(size_t i=0; i < n; i+=2){
        __m128d data = _mm_load_pd(hdata + i);
        __m128d truncated = _mm_max_pd(data, zero);
        _mm_store_pd(hdata + i, truncated);
        sum = _mm_add_pd(sum, truncated);
    }
    double sa[2];
    _mm_storeu_pd(sa, sum);
    h_sum = sa[0] + sa[1];
#endif
    if(normalize_to_nominal && h_sum > 0.0){
       t *= h0_sum / h_sum;
    }
}

void cubiclinear_histomorph::add_with_coeff_to(Histogram1DWithUncertainties & hres, double coeff, const ParValues & values) const{
    h_wu = h0_wu;
    add_morph_terms(h_wu, values);
    hres.add_with_coeff(coeff, h_wu);
}

void cubiclinear_histomorph::add_with_coeff_to(Histogram1D & hres, double coeff, const ParValues & values) const{
    h = h0;
    add_morph_terms(h, values);
    hres.add_with_coeff(coeff, h);
}

void cubiclinear_histomorph::get_histogram_dimensions(size_t & nbins, double & xmin, double & xmax) const{
    nbins = h0.get_nbins();
    xmin = h0.get_xmin();
    xmax = h0.get_xmax();
}

cubiclinear_histomorph::cubiclinear_histomorph(const Configuration & ctx): normalize_to_nominal(false) {
    boost::shared_ptr<VarIdManager> vm = ctx.pm->get<VarIdManager>();
    //build nominal histogram:
    h0_wu = get_constant_histogram(Configuration(ctx, ctx.setting["nominal-histogram"]));
    h0 = h0_wu.get_values_histogram();
    if(ctx.setting.exists("normalize_to_nominal")){
        normalize_to_nominal = ctx.setting["normalize_to_nominal"];
    }
    Setting psetting = ctx.setting["parameters"];
    size_t n = psetting.size();
    parameter_factors.resize(n, 1.0);
    bool have_parameter_factors = ctx.setting.exists("parameter_factors");
    for(size_t i=0; i<n; i++){
        string par_name = psetting[i];
        ParId pid = vm->get_par_id(par_name);
        par_ids.insert(pid);
        vid.push_back(pid);
        string setting_name;
        //plus:
        setting_name = par_name + "-plus-histogram";
        hplus_diff.push_back(get_constant_histogram(Configuration(ctx, ctx.setting[setting_name])).get_values_histogram());
        hplus_diff.back().check_compatibility(h0);
        hplus_diff.back().add_with_coeff(-1.0, h0);
        //minus:
        setting_name = par_name + "-minus-histogram";
        hminus_diff.push_back(get_constant_histogram(Configuration(ctx, ctx.setting[setting_name])).get_values_histogram());
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
    h0_sum = h0.get_sum();
    h = h0;
    h_wu = h0_wu;
}

REGISTER_PLUGIN(cubiclinear_histomorph)
