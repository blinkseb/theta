#include "plugins/asimov_likelihood_widths.hpp"
#include "plugins/secant.hpp"
#include "interface/distribution.hpp"
#include "interface/model.hpp"
#include "interface/exception.hpp"

#include <sstream>


using namespace theta;
using namespace std;

namespace{

//function object depending on one double
// which evaluates the nll using parameter values at some given point
// except one parameter. Also subtracts value at the mode such that
// the value there is 0.
struct nll_mode_pid{
   nll_mode_pid(const ParValues & mode_, const ParId & pid_, const NLLikelihood & nll_, double subtract_value_): values(mode_), pid(pid_),
      subtract_value(subtract_value_), nll(nll_){
   }
   
   double operator()(double p) const{
       values.set(pid, p);
       return nll(values) - subtract_value;
   }
private:
   mutable ParValues values;
   const ParId pid;
   double subtract_value;
   const NLLikelihood & nll;
};

}

theta::ParValues asimov_likelihood_widths(const theta::Model & model, const boost::shared_ptr<Distribution> & override_parameter_distribution, const boost::shared_ptr<theta::Function> & additional_nll_term){
    const Distribution & dist = override_parameter_distribution.get()? *override_parameter_distribution: model.get_parameter_distribution();
    ParIds parameters = model.get_parameters();
    ParValues mode;
    dist.mode(mode);
    Data asimov_data;
    model.get_prediction(asimov_data, mode);
    
    const Distribution * rvobs_dist = model.get_rvobservable_distribution();
    if(rvobs_dist){
        rvobs_dist->mode(mode);
        asimov_data.set_rvobs_values(ParValues(mode, model.get_rvobservables()));
    }
    std::auto_ptr<NLLikelihood> nll = model.get_nllikelihood(asimov_data);
    //0 value has same semantics for NLLikelihood:
    nll->set_override_distribution(override_parameter_distribution);
    nll->set_additional_term(additional_nll_term);
    if(additional_nll_term.get()){
        parameters.insert_all(additional_nll_term->get_parameters());
    }
    double nll_at_min = (*nll)(mode);
    ParValues result;
    int k=0;
    for(ParIds::const_iterator it=parameters.begin(); it!=parameters.end(); ++it, ++k){
        ParId pid = *it;
        const double pid_mode = mode.get(pid);
        std::pair<double, double> support = dist.support(pid);
        theta_assert(support.first <= pid_mode && pid_mode <= support.second);
        if(support.first == support.second){
            result.set(pid, 0.0);
            continue;
        }
        nll_mode_pid f(mode, pid, *nll, nll_at_min + 0.5);
        //if one end is finite, try to use it. Save whether the interval end is considered
        // "fl0", i.e. the interval end itself is finite but the function value there is invalid (< 0).
        bool low_is_fl0 = false, high_is_fl0 = false;
        if(std::isfinite(support.second)){
            double f2 = f(support.second);
            if(f2==0.0){
                result.set(pid, fabs(pid_mode - support.second));
                continue;
            }
            if(!std::isfinite(f2) || f2 < 0){
               low_is_fl0 = true;
            }
            else{
               result.set(pid, fabs(pid_mode - secant(pid_mode, support.second, 0.0, -0.5, f2, 0.05, f)));
               continue;
            }
        }
        if(std::isfinite(support.first)){
            double f2 = f(support.first);
            if(f2==0.0){
                result.set(pid, fabs(pid_mode - support.first));
                continue;
            }
            if(!std::isfinite(f2) || f2 < 0){
               high_is_fl0 = true;
            }
            else{
               result.set(pid, fabs(pid_mode - secant(support.first, pid_mode, 0.0, f2, -0.5, 0.05, f)));
               continue;
            }
        }
        //the support was either infinite or the values at the borders were not sufficiently high.
        // Treat second case first:
        if(low_is_fl0 && high_is_fl0){
            result.set(pid, support.second - support.first);
            continue;
        }
        //Now, one of the interval ends has to be infinite, otherwise we would not be here.
        //Scan in that direction:
        theta_assert(std::isinf(support.first) || std::isinf(support.second));
        bool found = false;
        for(double sign = -1.0; sign <= 1.001; sign+=2.0){
            if(!std::isinf(support.first) && sign < 0) continue;
            if(!std::isinf(support.second) && sign > 0) continue;
            // as step size, try the parameter value, if it is not zero:
            double step = fabs(pid_mode);
            if(step==0) step = 1.0;
            for(int i=0; i<1000; ++i){
                double fval = f(pid_mode + sign * step);
                if(isinf(fval)){
                    step /= 1.5;
                    continue;
                }
                step *= 2.0;
                if(fval > 0){
                    double xlow, xhigh, flow, fhigh;
                    xlow = pid_mode; flow = -0.5;
                    xhigh = pid_mode + sign * step; fhigh = fval;
                    if(sign < 0){
                        std::swap(xlow, xhigh);
                        std::swap(flow, fhigh);
                    }
                    theta_assert(xlow <= xhigh);
                    result.set(pid, fabs(pid_mode - secant(xlow, xhigh, 0.0, flow, fhigh, 0.05, f)));
                    found = true;
                    break;
                }
            }
            if(found) break;
        }
        if(found) continue;
        stringstream ss;
        ss << "asimov_likelihood_widths: could not find width for parameter " << pid;
        throw Exception(ss.str());
    }
    return result;
}

