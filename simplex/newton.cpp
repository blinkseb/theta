#include "simplex/newton.hpp"
#include "interface/asimov-utils.hpp"
#include "interface/phys.hpp"
#include "interface/plugin.hpp"
#include "interface/variables.hpp"
#include "interface/matrix.hpp"
#include "interface/cfg-utils.hpp"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>

using namespace theta;
using namespace std;
using namespace newton_internal;


namespace newton_internal {

class Linesearch{
public:
    // do a line search for a minimum starting at x0, with an estimate for the minimum at x0 + step.
    // The result is filled in x.
    // x does not have to be on the line defined by x0 and step. (this allows an easier
    // handling of constraints). The result point x must be in the range of f.
    // All memory belongs to the caller.
    virtual void do_linesearch(const RangedFunction & f, const std::vector<double> & x0, const std::vector<double> & step, std::vector<double> & x) const = 0;
    virtual ~Linesearch(){}
};


}


namespace{


// Function info for the newton Minimizer. So far, we only add epsilon_f.
class NewtonFunctionInfo: public FunctionInfo{
private:
    double epsilon_f;
    
public:
    NewtonFunctionInfo(const ParValues & start, const ParValues & step, const Ranges & ranges, const ParValues & fixed_parameters, double epsilon_f_):
        FunctionInfo(start, step, ranges, fixed_parameters), epsilon_f(epsilon_f_){
        theta_assert(epsilon_f > 0.0);
    }
    
    NewtonFunctionInfo(const FunctionInfo & info, double epsilon_f_): FunctionInfo(info), epsilon_f(epsilon_f_){
        theta_assert(epsilon_f > 0.0);
    }
    
    double get_epsilon_f() const{
        return epsilon_f;
    }
};


class IntervalLinesearch: public Linesearch {
private:
    double eps;
    bool debug;
    
    struct eval_f{
        const RangedFunction & f;
        const vector<double> & x0, step;
        vector<double> x;
        
        eval_f(const RangedFunction & f_, const vector<double> & x0_, const vector<double> & step_): f(f_), x0(x0_), step(step_){
            x.resize(x0.size());
        }
        
        double operator()(double c){
            x = x0;
            add_with_coeff(x, step, c);
            f.trunc_to_range(x);
            return f(x);
        }
    };
    
    struct ab_small{
        const RangedFunction & f;
        const vector<double> & x0, step;
        double eps;
        vector<double> x, y;
        
        ab_small(const RangedFunction & f_, const vector<double> & x0_, const vector<double> & step_, double eps_): f(f_), x0(x0_), step(step_), eps(eps_){
            x.resize(x0.size());
            y.resize(x0.size());
        }
        
        bool operator()(double a, double b){
            x = x0;
            add_with_coeff(x, step, a);
            f.trunc_to_range(x);
            y = x0;
            add_with_coeff(y, step, b);
            f.trunc_to_range(y);
            add_with_coeff(x, y, -1.0);
            return pnorm(x, step) < eps;
        }
    };
    
    
public:
    explicit IntervalLinesearch(double eps_, bool debug_ = false): eps(eps_), debug(debug_){}
    
    void do_linesearch(const RangedFunction & f, const vector<double> & x0, const vector<double> & step, vector<double> & x_new) const{
        const double f0 = f(x0);
        double fbest = f0;
        double cbest = 0.0; // best coefficient c found so far for   x = x0 + c * step
        // 1. find start of coefficients (a,b) in which the minimum is contained.
        // We assume that step is already a good start, so test nearby points:
        // TODO: delta could be tuned (?!)
        const double delta0 = 0.1;
        const size_t n = f.ndim();
        double a = 1 - delta0;
        
        // 1.a. make sure the lower point a is within the function range:
        vector<double> x(x0);
        add_with_coeff(x, step, a);
        // search for a value of a which is not outside the function range (that usually does not happen very often):
        bool found_a = false;
        for(int i=0; i < int(1/delta0 + 2); ++i){
            if(f.trunc_to_range(x) < n){
                found_a = true;
                break;
            }
            a -= delta0;
            x = x0;
            add_with_coeff(x, step, a);
        }
        theta_assert(found_a);
        
        // 1.b. find a, cbest such that a < cbest   and f(cbest) <= f(a).
        // the 1d function, depending on the coefficient of the step vector:
        eval_f f1d(f, x0, step);
        double fa = f1d(a);
        double delta = delta0;
        // find  a  such that fa > fbest and a < cbest
        while(fa < fbest || a >= cbest){
            if(fa < fbest){
                fbest = fa;
                cbest = a;
            }
            a -= delta;
            fa = f1d(a);
            delta *= 2;
        }
        theta_assert(fa >= fbest and a < cbest);
        // note that in rare cases, a can become negative, but that's Ok.
        
        // 1.c. find b > cbest with f(xb) > f(xbest)
        delta = max(cbest - a, delta0);
        double b = cbest + delta;
        double fb = f1d(b);
        while(fb < fbest){
            fbest = fb;
            cbest = b;
            b += delta;
            fb = f1d(b);
            delta *= 2;
        }
        // 2. using a < c < b with fa >= fc <= fb, search for the minimum of f.
        ab_small small_enough(f, x0, step, eps);
        double cmin = find_argmin(boost::ref(f1d), a, b, cbest, boost::ref(small_enough), fa, fb, fbest, 10000, debug);
        x_new = x0;
        add_with_coeff(x_new, step, cmin);
        f.trunc_to_range(x_new);
    }
};


} // anon. namespace

   
std::ostream & operator<<(std::ostream & out, const std::vector<double> & x){
    const size_t n = x.size();
    for(size_t i=0; i<n; ++i){
        out << x[i] << " ";
    }
    return out;
}

// PROBLEM: operator<< is defined globally


newton_minimizer::newton_minimizer(const newton_minimizer:: options & opts_): opts(opts_){
    ls.reset(new IntervalLinesearch(opts.par_eps, opts.debug));
}

newton_minimizer::newton_minimizer(const Configuration & cfg) {
    Setting s = cfg.setting;
    if(s.exists("maxit")){
        opts.maxit = s["maxit"];
    }
    if(s.exists("improve_cov")){
        opts.improve_cov = s["improve_cov"];
        // change default par_eps; might be overridden by user settings below.
        opts.par_eps = 1e-6;
    }
    if(s.exists("par_eps")){
        opts.par_eps = s["par_eps"];
    }
    if(s.exists("debug")){
        opts.debug = s["debug"];
    }
    ls.reset(new IntervalLinesearch(opts.par_eps, opts.debug));
}

namespace{
    
bool update_ihessian(Matrix & ih, const vector<double> & x, const vector<double> & xnew, const vector<double> & g, const vector<double> & gnew, bool debug){
    const size_t n = x.size();
    vector<double> dx(xnew);
    add_with_coeff(dx, x, -1.0);
    // calculate z = (Delta x - inverse_hessian * (next_grad - grad))
    vector<double> z(n);
    for(size_t i=0; i<n; ++i){
        z[i] = dx[i];
        for(size_t j=0; j<n; ++j){
            z[i] -=  ih(i, j) * (gnew[j] - g[j]);
        }
    }
    // calculate the denominator, z * (next_grad - grad):
    double denom = 0.0;
    double norm_graddiff = 0.0;
    for(size_t i=0; i<n; ++i){
        double graddiff = gnew[i] - g[i];
        norm_graddiff += graddiff * graddiff;
        denom += z[i] * graddiff;
    }
    double norm_z = norm(z);
    norm_graddiff = sqrt(norm_graddiff);
    if(debug)cout << "denom = " << denom << "; |z| = " << norm_z << "; |g_k+1 - g_k| = " << norm_graddiff << endl;
    if(fabs(denom) > 1e-7 * norm_graddiff * norm_z){
        // calculate the update to inverse_hessian: add z * z^T / denom:
        for(size_t i=0; i<n; ++i){
            for(size_t j=0; j< n; ++j){
                ih(i,j) += z[i] * z[j] / denom;
            }
        }
        return true;
    }
    else{
        if(debug) cout << "denom is too small for update; continuing with next iteration directly " << endl;
        return false;
    }
}
    
}


MinimizationResult newton_minimizer::minimize2(const theta::Function & f_, const theta::FunctionInfo & info_, const theta::ParValues & fixed_parameters){
    const NewtonFunctionInfo & info = dynamic_cast<const NewtonFunctionInfo &>(info_);
    RangedThetaFunction f(f_, fixed_parameters, info.get_step(), info.get_ranges());
    f.set_epsilon_f(info.get_epsilon_f());
    const ParIds & all_pids = f_.get_parameters();
    const ParIds & fixed_pids = f.get_fixed_parameters();
    theta_assert(fixed_pids == info.get_fixed_parameters());
    ParIds nonfixed_pids;
    for(ParIds::const_iterator it=all_pids.begin(); it!=all_pids.end(); ++it){
        if(fixed_pids.contains(*it)) continue;
        nonfixed_pids.insert(*it);
    }
    const size_t n = f.ndim();
    theta_assert(n == nonfixed_pids.size());
    // fill the matrix with step**2 at the diagonals:
    // TODO: for the future, could do that in FunctionInfo already and use non-diagonal elements as well ...
    ParValues pv_step = info.get_step();
    Matrix inverse_hessian(n, n);
    vector<double> step(n);
    {
        size_t i=0;
        for(ParIds::const_iterator it=nonfixed_pids.begin(); i<n; ++i, ++it){
            double st = pv_step.get(*it);
            inverse_hessian(i,i) = pow(st, 2);
            step[i] = st;
        }
    }
    vector<double> x(n), grad(n);
    vector<double> direction(n), next_x(n), next_grad(n);
    ParValues pv_start = info.get_start();
    pv_start.set(fixed_parameters);
    pv_start.fill(&x[0], nonfixed_pids);
    double fx = f.eval_with_derivative(x, grad);
    int it = 0;
    for(; it < opts.maxit; ++it){
        if(opts.debug){
            // dump current status:
            cout << endl << "Starting iteration " << it << ": " << endl << "x = " << x << " --> " << fx << endl;
            cout << "  g = " << grad << endl;
            cout << endl << "h = " << endl;
            for(size_t k=0; k<n; ++k){
                for(size_t l=0; l<n; ++l){
                    printf(" %8.4g", inverse_hessian(k, l));
                }
                cout << endl;
            }
        }
        
        // the estimated minimum is   -inverse_hessian * grad
        mul(direction, inverse_hessian, grad);
        mul(direction, -1.0);
        if(opts.debug){
            cout << "Iteration " << it << ": search direction = " << direction << endl;
        }
        if(pnorm(direction, step) < opts.par_eps){
            if(opts.debug){
                cout << "Iteration " << it << ": estimated minimum close enough; stopping iteration." << endl;
            }
            break;
        }
        if(opts.debug)cout << "Iteration " << it << ": doing linesearch now" << endl;
        ls->do_linesearch(f, x, direction, next_x);
        if(opts.debug){
            cout << "Iteration " << it << ": linesearch proposes:" << endl << "x =" << next_x << endl;
        }
        // calculate dx = next_x - x
        vector<double> dx(next_x);
        add_with_coeff(dx, x, -1.0);
        if(pnorm(dx, step) < opts.par_eps){
            if(opts.debug){
                 cout << "Iteration " << it << ": actual dx was smaller than par_eps, stopping iteration." << endl;
            }
            break;
        }
        // do the update to inverse_hessian:
        double next_fx = f.eval_with_derivative(next_x, next_grad);
        update_ihessian(inverse_hessian, x, next_x, grad, next_grad, opts.debug);
        // prepare for next iteration:
        swap(x, next_x);
        swap(grad, next_grad);
        fx = next_fx;
    }
    if(it==opts.maxit){
        throw MinimizationException("maximum number of iterations reached");
    }
    if(opts.improve_cov){
        // calculate better estimate for covariance here by going to 0.1 sigma ("posterior") in the
        // direction of the parameter and calculating the update.
        theta_assert(f.trunc_to_range(x)==0);
        for(size_t i=0; i<n; ++i){
            vector<double> next_x(x);
            const double step = 0.1 * sqrt(inverse_hessian(i,i));
            next_x[i] = x[i] + step;
            if(f.trunc_to_range(next_x) > 0){
                next_x[i] = x[i] - step;
                if(f.trunc_to_range(next_x) > 0){
                    throw MinimizationException("for calculating covariance: range too small (<0.1 sigma)");
                }
            }
            f.eval_with_derivative(next_x, next_grad);
            update_ihessian(inverse_hessian, x, next_x, grad, next_grad, opts.debug);
        }
    }
    MinimizationResult res;
    res.fval = fx;
    res.values.set(pv_start); // to get the fixed parameters right
    res.values.set(ParValues(&x[0], nonfixed_pids));
    // write the covariance we have so far; make sure to use 0 for fixed parameters:
    res.covariance = Matrix(all_pids.size(), all_pids.size());
    size_t i_nf = 0, j_nf = 0; // count non-fixed parameters
    size_t i=0; // index to all parameters
    for(ParIds::const_iterator pit=all_pids.begin(); pit!=all_pids.end(); ++pit, ++i){
        size_t j=0;
        j_nf = 0;
        for(ParIds::const_iterator pit2=all_pids.begin(); pit2!=all_pids.end(); ++pit2, ++j){
            if(nonfixed_pids.contains(*pit) && nonfixed_pids.contains(*pit2)){
                res.covariance(i,j) = inverse_hessian(i_nf, j_nf);
                if(i==j){
                    double error = sqrt(res.covariance(i,j));
                    res.errors_plus.set(*pit, error);
                    res.errors_minus.set(*pit, error);
                }
            }
            else{
                res.covariance(i,j) = 0.0;
                if(i==j){
                    res.errors_plus.set(*pit, 0.0);
                    res.errors_minus.set(*pit, 0.0);
                }
            }
            if(nonfixed_pids.contains(*pit2)) ++j_nf;
        }
        theta_assert(j_nf == nonfixed_pids.size());
        theta_assert(j == all_pids.size());
        if(nonfixed_pids.contains(*pit)) ++i_nf;
    }
    theta_assert(i == all_pids.size());
    theta_assert(i_nf == nonfixed_pids.size());
    return res;
}


MinimizationResult newton_minimizer::minimize(const theta::Function & f_, const theta::ParValues & start,
                                              const theta::ParValues & step, const Ranges & ranges){
    ParValues fixed_parameters;
    RangedThetaFunction f(f_, fixed_parameters, step, ranges);
    // get epsilon_f:
    vector<double> x0(f.ndim());
    start.fill(&x0[0], f.get_nonfixed_parameters());
    double epsilon_f = f_accuracy(f, x0, 0);
    if(opts.debug) cout << "epsilon_f = " << epsilon_f << endl;
    NewtonFunctionInfo info(start, step, ranges, fixed_parameters, epsilon_f);
    return minimize2(f_, info, fixed_parameters);
}

boost::shared_ptr<theta::FunctionInfo> newton_minimizer::create_nll_function_info(const theta::Model & m,
                                                                                  const boost::shared_ptr<theta::Distribution> & override_parameter_distribution,
                                                                                  const theta::ParValues & fixed_parameters){
    boost::shared_ptr<FunctionInfo> info = Minimizer::create_nll_function_info(m, override_parameter_distribution, fixed_parameters);
    // get epsilon_f for the asimov likelihood:
    asimov_data_nll nll(m, override_parameter_distribution, fixed_parameters);
    RangedThetaFunction f(nll, fixed_parameters, info->get_step(), info->get_ranges());
    vector<double> x0(f.ndim());
    ParValues start = info->get_start();
    start.fill(&x0[0], f.get_nonfixed_parameters());
    double epsilon_f = f_accuracy(f, x0, 0);
    boost::shared_ptr<NewtonFunctionInfo> result(new NewtonFunctionInfo(*info, epsilon_f));
    return result;
}


REGISTER_PLUGIN(newton_minimizer)

