#include "simplex/newton-utils.hpp"
#include "interface/phys.hpp"
#include "interface/distribution.hpp"
#include <limits>
#include <iostream>

using namespace std;
using namespace theta;
using namespace newton_internal;

double RangedFunction::eval_with_derivative(const vector<double> & x0, vector<double> & grad) const{
    const size_t n = ndim();
    theta_assert(x0.size()==n);
    vector<double> x(x0);
    grad.resize(n);
    double f0 = operator()(x);
    for(size_t i=0; i<n; ++i){
        double x0i = x0[i];
        const double sigma = step[i];
        theta_assert(sigma > 0);
        // we choose h ~ sqrt(epsilon_f / f'') where f'' is the second derivative.
        // Here, f''  ~  1 / step**2
        double h_ = sqrt(epsilon_f) * sigma;
        double x0prime = x0i + h_;
        if(x0prime > range_max[i] || x0prime < range_min[i]){
            x0prime = x0i - h_;
            if(x0prime > range_max[i] || x0prime < range_min[i]){
                throw invalid_argument("parameter range is to small");
            }
        }
        theta_assert(!std::isnan(x0prime));
        x[i] = x0prime;
        double fprime = operator()(x);
        x[i] = x0[i];
        theta_assert(!std::isnan(fprime));
        volatile double h = x0prime - x0i;
        double g = (fprime - f0) / h;
        // estimate the error: roundoff error + truncation error:
        //double g_err = epsilon_f / h + h / pow(sigma, 2);
        grad[i] = g;
    }
    return f0;
}

size_t RangedFunction::trunc_to_range(vector<double> & x) const{
    const size_t n = ndim();
    x.resize(n);
    size_t result = 0;
    for(size_t i=0; i<n; ++i){
        if(x[i] > range_max[i]){
            x[i] = range_max[i];
            ++result;
        }
        else if(x[i] < range_min[i]){
            x[i] = range_min[i];
            ++result;
        }
    }
    return result;
}

RangedFunction::~RangedFunction(){}

RangedThetaFunction::RangedThetaFunction(const Function & f_, const ParValues & pv_fixed_values, const ParValues & pv_step, const Ranges & ranges): f(f_){
    epsilon_f = numeric_limits<double>::epsilon();
    const ParIds & pids = f.get_parameters();
    values.resize(pids.size());
    step.reserve(pids.size());
    nonfixed_indices.reserve(pids.size());
    range_min.reserve(pids.size());
    range_max.reserve(pids.size());
    size_t i=0;
    for(ParIds::const_iterator pit=pids.begin(); pit!=pids.end(); ++pit, ++i){
        if(pv_fixed_values.contains(*pit)){
            values[i] = pv_fixed_values.get_unchecked(*pit);
            fixed_parameters.insert(*pit);
        }
        else{
            const pair<double, double> & range = ranges.get(*pit);
            if(pv_step.get(*pit)==0.0){
                theta_assert(range.first == range.second);
                values[i] = range.first;
                fixed_parameters.insert(*pit);
            }
            else{
                nonfixed_indices.push_back(i);
                nonfixed_parameters.insert(*pit);
                step.push_back(pv_step.get(*pit));
                theta_assert(step.back() > 0.0);
                range_min.push_back(range.first);
                range_max.push_back(range.second);
            }
        }
    }
}
    
void RangedThetaFunction::set_epsilon_f(double newval){
    theta_assert(newval > 0);
    epsilon_f = newval;
}

const ParIds & RangedThetaFunction::get_fixed_parameters() const{
    return fixed_parameters;
}

const ParIds & RangedThetaFunction::get_nonfixed_parameters() const{
    return nonfixed_parameters;
}
    
size_t RangedThetaFunction::ndim() const{
    return nonfixed_indices.size();
}

// vals are the ndim() non-fixed parameters in conversion order defined by pids_ as passed to the constructor.
double RangedThetaFunction::operator()(const vector<double> & vals) const {
    theta_assert(vals.size() == nonfixed_indices.size());
    for(size_t i=0; i<nonfixed_indices.size(); ++i){
        values[nonfixed_indices[i]] = vals[i];
    }
    return f(&values[0]);
}


double newton_internal::find_argmin(const boost::function<double (double)> & f, double a, double b, double c,
                                    const boost::function<bool (double, double)> & ab_small_enough, double fa, double fb, double fc, unsigned int maxit, bool debug){
    theta_assert(a < c && c < b);
    theta_assert(fa >= fc and fb >= fc);
    unsigned int it = 0;
    size_t f_eval = 0;
    for(; it < maxit; ++it){
        if(ab_small_enough(a,b)){
            if(debug) cout << "interval linesearch done. f_eval = " << f_eval << endl;
            return c;
        }
        if(debug) cout << "interval linesearch iteration " << it << "; at: " << a << " " << c << " " << b << endl
                       << "  f = " << fa << " >= " << fc << " <= " << fb << endl
                       << "  f_eval = " << f_eval << endl;
        theta_assert(a <= c and c <= b);
        theta_assert(fa >= fc and fb >= fc);
        // for now, use interval bisection. In the future, could use Brent.
        // bisect either a--c or c--b, using the larger one, to make most gain:
        if(c - a > b - c){
            double c_trial = a + 0.5 * (c-a); // a < c_trial < c < b
            if(debug) cout << "  c_trial = " << c_trial << endl;
            double fct = f(c_trial);
            ++f_eval;
            if(debug) cout << "  f(c_trial) = " << fct << endl;
            if(fct <= fc){
                // we found an even smaller function value; the new search triplet would be   a < c_trial < c. However,
                // c can be very close to b in which case the interval hardly shrinks. To avoid that, test yet
                // another point at (c + c_trial) / 2 if c is close to b:
                if(b - c < 0.5 * (c - a)){ 
                    double c_trial2 = c + 0.5 * (c_trial - c); // a < c_trial < c_trial2 < c < b
                    if(debug) cout << "  c_trial2 = " << c_trial2 << endl;
                    double fct2 = f(c_trial2);
                    ++f_eval;
                    if(debug) cout << "  f(c_trial2) = " << fct2 << endl;
                    if(fct2 <= fct){
                        // new triplet is c_trial < c_trial2 < c
                        a = c_trial;
                        fa = fct;
                        b = c;
                        fb = fc;
                        c = c_trial2;
                        fc = fct2;
                    }
                    else{
                        // new triplet is a < c_trial < c_trial2
                        b = c_trial2;
                        fb = fct2;
                        c = c_trial;
                        fc = fct;
                    }
                }
                else{ // i.e. b and c are not so close: just use triplet a < c_trial < c:
                    b = c;
                    fb = fc;
                    c = c_trial;
                    fc = fct;
                }
            }
            else{
                // c_trial < c < b
                a = c_trial;
                fa = fct;
            }
        }
        else{ // test c_trial between c and b, so we have   a < c < c_trial < b
            double c_trial = c + 0.5 * (b-c);
            if(debug) cout << "  c_trial = " << c_trial << endl;
            double fct = f(c_trial);
            ++f_eval;
            if(debug) cout << "  f(c_trial) = " << fct << endl;
            if(fct <= fc){
                if(c - a < 0.5 * (b - c)){
                    // a < c < c_trial2 < c_trial < b
                    double c_trial2 = c + 0.5 * (c_trial - c);
                    if(debug) cout << "  c_trial2 = " << c_trial2 << endl;
                    double fct2 = f(c_trial2);
                    ++f_eval;
                    if(debug) cout << "  f(c_trial2) = " << fct2 << endl;
                    if(fct2 <= fct){
                        // new triplet is c < c_trial2 < c_trial
                        a = c;
                        fa = fc;
                        c = c_trial2;
                        fc = fct2;
                        b = c_trial;
                        fb = fct;
                    }
                    else{
                        // new triplet is c_trial2 < c_trial < b
                        a = c_trial2;
                        fa = fct2;
                        c = c_trial;
                        fc = fct;
                    }
                }
                else{
                    // new triplet: c < c_trial < b
                    a = c;
                    fa = fc;
                    c = c_trial;
                    fc = fct;
                }
            }
            else{
                // new triplet: a < c < c_trial
                b = c_trial;
                fb = fct;
            }
        }
    }
    
    throw Exception("too many iterations to find minimum");
}



double newton_internal::f_accuracy(const RangedFunction & f, const vector<double> & x0, size_t i, double f_scale){
    theta_assert(f_scale > 0.0);
    theta_assert(i < x0.size());
    theta_assert(f.ndim() == x0.size());
    vector<double> x(x0);
    size_t f_eval = 0;
    theta_assert(f.trunc_to_range(x) == 0);
    const double f0 = fabs(f(x)) + f_scale;
    ++f_eval;
    double h = numeric_limits<double>::epsilon() * fabs(f0);
    x[i] = x0[i] + h;
    if(f.trunc_to_range(x) > 0){
        h = -h;
        x[i] = x0[i] + h;
        if(f.trunc_to_range(x) > 0) throw invalid_argument("f_accuracy: parameter range to small");
    }
    double f1 = fabs(f(x)) + f_scale;
    ++f_eval;
    // two cases: either f0==f1, then we have to increase h, or f0!=f1, then we have to decrease h:
    bool increase_h = f0 == f1;    
    double f_eps = increase_h ? numeric_limits<double>::infinity() : fabs(f1 - f0);
    for(int j=0; j<1025; ++j){ // 2**+-1024 = double_max / double_min
        if(increase_h) h *= 2;
        else h /= 2;
        x[i] = x0[i] + h;
        if(f.trunc_to_range(x) > 0){
            throw invalid_argument("f_accuracy: parameter range to small");
        }
        double f1 = fabs(f(x)) + f_scale;
        ++f_eval;
        if(f1 != f0){
            f_eps = min(f_eps, fabs(f1 - f0));
            // in increasing mode, this is the first non-zero difference, and we are done:
            if(increase_h){
                return f_eps;
            }
        }
        else{
            // if in decresing mode, this is the first zero difference, and we are done:
            if(!increase_h){
                return f_eps;
            }
        }
    }
    throw logic_error("too many iterations to find f_eps (f constant?)");
}


