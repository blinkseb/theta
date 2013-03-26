#include "root/root_minuit1.hpp"
#include "interface/plugin.hpp"
#include "interface/phys.hpp"

#include "interface/exception.hpp"
#include "TMinuit.h"

using namespace theta;
using namespace std;

namespace{
    
// for TMinuit, function evaluation is provided by deriving from TMinuit and overloading Eval.
class MyTMinuit: public TMinuit{
private:
    const theta::Function & f;
    size_t ndim;
    vector<double> min_;
    double f_at_min;
    
public:
    
    explicit MyTMinuit(const Function & f_): TMinuit(f_.get_parameters().size()), f(f_), ndim(f.get_parameters().size()){
        min_.resize(ndim);
    }
    
    const vector<double> & get_min() const{
        return min_;
    }
    
    double get_f_at_min() const{
        return f_at_min;
    }
    
    // see http://root.cern.ch/root/html/TMinuit.html#TMinuit:Eval
    virtual Int_t Eval(Int_t npar, Double_t * grad, Double_t & fval, Double_t * par, Int_t flag){
        if(flag==3){
            // this means we are done: save the parameters:
            std::copy(par, par + ndim, min_.begin());
            f_at_min = f(par);
        }
        else if(flag==1){
            theta_assert(npar >= 0 and static_cast<size_t>(npar) == ndim);
            for(size_t i=0; i<ndim; ++i){
                if(std::isnan(par[i])){
                    throw MinimizationException("minuit called likelihood function with NAN argument!");
                }
            }
            fval = f(par);
            if(std::isinf(fval)){
                throw MinimizationException("function to minimize was infinity during minimization");
            }
        }
        else if(flag==2){ // we should calculate gradient
            return 1; // is this how this is supposed to work?. No.
        }
        return 0;
    }
};

}


MinimizationResult root_minuit1::minimize(const theta::Function & f, const theta::ParValues & start,
        const theta::ParValues & steps, const Ranges & ranges){
    MyTMinuit min(f);
    //1. setup parameters, limits and initial step sizes
    ParIds parameters = f.get_parameters();
    int ivar=0;
    for(ParIds::const_iterator it=parameters.begin(); it!=parameters.end(); ++it, ++ivar){
        pair<double, double> range = ranges.get(*it);
        double def = start.get(*it);
        double step = steps.get(*it);
        stringstream ss;
        ss << "par" << ivar;
        std::string sname = ss.str();
        const char * name = sname.c_str();
        //use not the ranges directly, but a somewhat more narrow range (one permille of the respective border)
        // in order to avoid that the numerical evaluation of the numerical derivative at the boundaries pass these
        // boundaries ...
        if(step == 0.0){
            min.DefineParameter(ivar, name, def, 0.0, def, def);
            min.FixParameter(ivar);
        }
        else if(std::isinf(range.first)){
            if(std::isinf(range.second)){
                min.DefineParameter(ivar, name, def, step, def - infinity * step, def + infinity * step);
            }
            else{
                min.DefineParameter(ivar, name, def, step, def - infinity * step, range.second - fabs(range.second) * 0.001);
            }
        }
        else{
            if(std::isinf(range.second)){
                min.DefineParameter(ivar, name, def, step, range.first + fabs(range.first) * 0.001, def + infinity * step);
            }
            else{ // both ends are finite
                theta_assert(range.first < range.second);
                min.DefineParameter(ivar, name, def, step, range.first + fabs(range.first) * 0.001, range.second - fabs(range.second) * 0.001);
            }
        }
    }
    min.SetErrorDef(0.5);
    //min.SetPrintLevel(0);
    min.mngrad(); // let minuit see that we don't calculate the gradient. TODO: really needed?
    int res = min.Migrad();
    if(res!=0){
        for(unsigned int i=1; i<=n_retries; i++){
            res = min.Migrad();
            if(res==0) break;
        }
    }
    if(res!=0){
        stringstream s;
        s << "MINUIT returned " << res;
        throw MinimizationException(s.str());
    }
    if(hesse){
        min.mnhess();
    }
    
    // done with minimization, convert result:
    const vector<double> & xmin = min.get_min();
    MinimizationResult result;
    result.fval = min.get_f_at_min();
    ivar = 0;
    for(ParIds::const_iterator it=parameters.begin(); it!=parameters.end(); ++it, ++ivar){
        result.values.set(*it, xmin[ivar]);
        double eplus = 0.0, eminus = 0.0, eparab = 0.0, dummy;
        min.mnerrs(ivar, eplus, eminus, eparab, dummy); // should be set to 0.0 if not available
        if(eplus > 0.0 and eminus > 0.0){
            result.errors_plus.set(*it, eplus);
            result.errors_minus.set(*it, eminus);
        }
        else if(eparab > 0.0){
            result.errors_plus.set(*it, eparab);
            result.errors_minus.set(*it, eparab);
        }
        else{
            result.errors_plus.set(*it, -1.0);
            result.errors_minus.set(*it, -1.0);
        }
    }
    const size_t n = parameters.size();
    vector<double> emat(n*n);
    min.mnemat(&emat[0], n);
    result.covariance.reset(n, n);
    for(size_t i=0; i<parameters.size(); ++i){
        for(size_t j=0; j<parameters.size(); ++j){
            result.covariance(i,j) = emat[i*n+j];
        }
    }
    return result;
}

root_minuit1::root_minuit1(const Configuration & cfg): tolerance(-1), infinity(1e5), n_retries(2), hesse(false) {
    if(cfg.setting.exists("n_retries")){
        n_retries = cfg.setting["n_retries"];
    }
    if(cfg.setting.exists("infinity")){
        infinity = cfg.setting["infinity"];
    }
    if(cfg.setting.exists("hesse")){
        hesse = cfg.setting["hesse"];
    }
    if(cfg.setting.exists("tolerance")){
        tolerance = cfg.setting["tolerance"];
        if(tolerance <= 0){
            throw ConfigurationException("tolerance <= 0.0 not allowed");
        }
    }
}

REGISTER_PLUGIN(root_minuit1)

