#include "simplex/newton.hpp"
#include "interface/phys.hpp"
#include "interface/variables.hpp"
#include "interface/matrix.hpp"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>

#include <map>

using namespace theta;
using namespace std;

namespace{
const double eps = numeric_limits<double>::epsilon();
const double s_eps = sqrt(numeric_limits<double>::epsilon());
const double inf = numeric_limits<double>::infinity();
}

class Ranges{
	map<ParId, pair<double, double> > ranges_;
public:
	Ranges(const map<ParId, pair<double, double> > & ranges): ranges_(ranges){
		map<ParId, pair<double, double> >::const_iterator it = ranges.begin();
		for(; it!=ranges.end(); ++it){
			theta_assert(it->second.first <= it->second.second);
		}
	}

	pair<double, double> get(const ParId & pid) const{
		map<ParId, pair<double, double> >::const_iterator it = ranges_.find(pid);
		if(it==ranges_.end()){
			throw invalid_argument("no range set for this parid");
		}
		return it->second;
	}

	bool fixed(const ParId & pid) const{
		map<ParId, pair<double, double> >::const_iterator it = ranges_.find(pid);
		if(it==ranges_.end()){
			throw invalid_argument("no range set for this parid");
		}
		return it->second.first == it->second.second;
	}

	void set(const ParId & pid, double low, double high){
		theta_assert(low <= high);
		ranges_[pid] = make_pair(low, high);
	}

	void set(const map<ParId, pair<double, double> > & ranges){
		map<ParId, pair<double, double> >::const_iterator it = ranges.begin();
		for(; it!=ranges.end(); ++it){
			set(it->first, it->second.first, it->second.second);
		}
	}

	// truncate the values to the ranges.
	void trunc(ParValues & values) const{
		for(map<ParId, pair<double, double> >::const_iterator it=ranges_.begin(); it!=ranges_.end(); ++it){
			if(!values.contains(it->first)) continue;
			double val = values.get(it->first);
			if(val < it->second.first) val = it->second.first;
			else if(val > it->second.second) val = it->second.second;
			values.set(it->first, val);
		}
	}
};

// evaluate function, respecting the parameter ranges
// values must be within the range. Does not update grad for fixed parameters.
double eval_with_grad(const Function & f, const Ranges & ranges, const ParValues & values, ParValues & grad){
	const ParIds & pids = f.get_parameters();
	const double f0 = f(values);
	for(ParIds::const_iterator it=pids.begin(); it!=pids.end(); ++it){
		const ParId & pid = *it;
		if(ranges.fixed(pid)) continue;
		pair<double, double> range = ranges.get(pid);
		double x0 = values.get(pid);
		double x0prime = x0 * (1 + s_eps);
		if(fabs(x0prime) < eps) x0prime = 2*x0;
		if(x0==0) x0 = eps;
		if(x0prime > range.second || x0prime < range.first){
			x0prime = x0 * (1 - s_eps);
			if(fabs(x0prime) < eps) x0prime = 2*x0;
			if(x0==0) x0 = -eps;
			if(x0prime > range.second || x0prime < range.first){
				throw invalid_argument("relative parameter range is to small");
			}
		}
		ParValues vals(values);
		theta_assert(!std::isnan(x0prime));
		vals.set(pid, x0prime);
		double fprime = f(vals);
		theta_assert(!std::isnan(fprime));
		volatile double h = x0prime - x0;
		double g = (fprime - f0) / h;
		if(std::isnan(g)){
			cout << fprime << " - " << f0 << " / " << h << " = " << g <<  endl;
		}
		theta_assert(!std::isnan(g));
		grad.set(pid, g);
	}
	return f0;
}


// returns the pair (f(x0), df(x0) / ddir). Make sure that direction is 0.0 for fixed parameters.
pair<double, double> eval_with_derivative(const Function & f, const Ranges & ranges, const ParValues & x0, const ParValues & direction){
	const ParIds & pids = f.get_parameters();
	const double f0 = f(x0);
	ParValues x(x0);
	double dir_norm = 0.0;
	double x0_norm = 0.0;
	for(ParIds::const_iterator it=pids.begin(); it!=pids.end(); ++it){
		double d = direction.get(*it);
		if(d != 0.0){
			theta_assert(not ranges.fixed(*it));
		}
		dir_norm +=  d*d;
		x0_norm += x0.get(*it) * x0.get(*it);
	}
	dir_norm = sqrt(dir_norm);
	x0_norm = sqrt(x0_norm);

	// if we are very close to zero, we should use some minimum value:
	if(x0_norm < s_eps) x0_norm = s_eps;

	// use x0 + s_eps * x0norm * direction / dir_norm. If along some axis, this would be outside
	// of the range, the step is taken in the other direction.
	double sign_switch = 1.0;
	for(ParIds::const_iterator it=pids.begin(); it!=pids.end(); ++it){
		const ParId & pid = *it;
		pair<double, double> range = ranges.get(pid);
		double xi = x0.get(pid);
		double dx = s_eps * x0_norm * direction.get(pid) / dir_norm;
		double xiprime = xi + dx;
		if(xiprime > range.second || xiprime < range.first){
			xiprime = xi - dx;
			sign_switch *= -1.0;
			if(xiprime > range.second || xiprime < range.first){
				cout << "very small range! " << endl;
				throw invalid_argument("range to small in the given direction");
			}
		}
		x.set(pid, xiprime);
	}
	double fprime = f(x);
	double der = sign_switch * (fprime - f0) / (s_eps * x0_norm);
	return make_pair(f0, der);
}


class Linesearch{
public:
	// do a line search for a minimum starting at x0 in the given direction, respecting the parameter
	// ranges given in r. Return the values of the "best" point in that direction.
	// x0 and x0 + dir are the proposed first search point; x0 is guaranteed to be within the range, while x0 + dir could be outside.
	//
	// The routine is not required to stay on the proposed line.
	virtual ParValues do_linesearch(const Function & f, const Ranges & r, const ParValues & x0, const ParValues & dir) const = 0;
	virtual ~Linesearch(){}
};


class NewtonLinesearch: public Linesearch{
private:
	static void add(ParValues & x0, const ParValues & dx, const ParIds & pids, double coeff = 1.0){
		for(ParIds::const_iterator it=pids.begin(); it!=pids.end(); ++it){
			x0.set(*it, x0.get(*it) + coeff * dx.get(*it));
		}
	}

	static double norm(const ParValues & n, const ParIds & pids){
		double result = 0.0;
		for(ParIds::const_iterator it=pids.begin(); it!=pids.end(); ++it){
			result += n.get(*it) * n.get(*it);
		}
		return sqrt(result);
	}

	static double normdiff(const ParValues & x0, const ParValues & x1, const ParIds & pids){
		double result = 0.0;
		for(ParIds::const_iterator it=pids.begin(); it!=pids.end(); ++it){
			double diff = x0.get(*it) - x1.get(*it);
			result += diff * diff;
		}
		return sqrt(result);
	}

	bool debug;
public:
	explicit NewtonLinesearch(bool debug_ = false): debug(debug_){

	}

	virtual ParValues do_linesearch(const Function & f, const Ranges & r, const ParValues & x0, const ParValues & dir) const{
		const ParIds & pars = f.get_parameters();
		pair<double, double> fdf0 = eval_with_derivative(f, r, x0, dir);
		if(debug)cout << "linesearch started with f = " << fdf0.first << "; g = " << fdf0.second << endl;
		ParValues x(x0);
		double coeff = 1.0;
		add(x, dir, f.get_parameters(), coeff);
		r.trunc(x);
		pair<double, double> fdf1 = eval_with_derivative(f, r, x, dir);
		size_t coeff_it = 0;
		const size_t max_coeff_it = 30;
		while(fdf1.first > fdf0.first && coeff_it < max_coeff_it){
			coeff *= 0.5;
			if(debug)cout << "Got larger function value after step, trying to decreasing step size to " << coeff << endl;
			x.set(x0);
			add(x, dir, f.get_parameters(), coeff);
			r.trunc(x);
			fdf1 = eval_with_derivative(f, r, x, dir);
			++coeff_it;
		}
		if(coeff_it == max_coeff_it){
			if(debug)cout << "WARNING: seems like function is always larger in the suggested direction, flipping direction!" << endl;
			coeff = -1.0;
			x.set(x0);
			add(x, dir, f.get_parameters(), coeff);
			r.trunc(x);
			fdf1 = eval_with_derivative(f, r, x, dir);
			coeff_it = 0;
			while(fdf1.first > fdf0.first && coeff_it < max_coeff_it){
				coeff *= 0.5;
				if(debug)cout << "Got larger function value after step, trying to decreasing step size to " << coeff << endl;
				x.set(x0);
				add(x, dir, f.get_parameters(), coeff);
				r.trunc(x);
				fdf1 = eval_with_derivative(f, r, x, dir);
				++coeff_it;
			}
			if(coeff_it == max_coeff_it){
				if(debug)cout << "WARNING: the other direction wasn't better either, seems we're at a good minimum already (?)" << endl;
				throw MinimizationException("function always larger");
			}
		}
		// now make one "newton step", x_n+1  =  x_n - f'(x_n) / f''(x_n), where f''(x_n) is approximated by the difference
		// of derivatives.
		if(debug)cout << "found valid function value at " << x.get(*pars.begin()) << ": " << fdf1.first << "; g = " << fdf1.second << endl;
		double fpp = (fdf1.second - fdf0.second) / normdiff(x, x0, pars);
		if(debug) cout << "fpp = " << fpp << endl;
		if(!std::isfinite(fpp)){ // only happens if points were too close ...
			if(debug)cout << "best point very close!" << endl;
			return x0;
		}
		if(fpp != 0.0){
			if(debug)cout << "going in that direction with factor " << -fdf0.second / fpp / norm(dir, pars) << endl;
			ParValues x_next(x0);
			add(x_next, dir, pars, -fdf0.second / fpp / norm(dir, pars));
			r.trunc(x_next);
			return x_next;
		}
		else{
			return x;
		}
	}
	virtual ~NewtonLinesearch(){}
};

newton_minimizer::newton_minimizer(const newton_minimizer:: options & opts_): opts(opts_){
	ls.reset(new NewtonLinesearch(opts.debug));
}


MinimizationResult newton_minimizer::minimize2(const theta::MinimizationProblem & mp){
	Matrix inverse_hessian = mp.get_step_matrix();
	// the fixed parameters have a step size of 0.0; to make the matrix inversion work, add a diagonal element there:
	const size_t n = inverse_hessian.get_n_cols();
	theta_assert(n == inverse_hessian.get_n_rows());
	/*vector<size_t> fixed_indices;
	for(size_t i=0; i<n; ++i){
		if(inverse_hessian(i,i) == 0.0){
			inverse_hessian(i,i) = 1.0;
			fixed_indices.push_back(i);
		}
	}*/
	//cout << "detected " << fixed_indices.size() << " fixed parameters, based on the step matrix" << endl;
	// TODO clear these rows and cols completely s.t. they stay fixed (or wrap function?!)
	Ranges r(mp.get_ranges());
	const theta::Function & f = mp.get_f();
	const ParIds & pids = f.get_parameters();
	ParValues x(mp.get_start());
	ParValues grad, zero;
	vector<ParId> vpids;
	for(ParIds::const_iterator it=pids.begin(); it!=pids.end(); ++it){
		grad.set(*it, 0.0);
		zero.set(*it, 0.0);
		vpids.push_back(*it);
	}
	theta_assert(vpids.size() == n);
	double fx = eval_with_grad(f, r, x, grad);
    int it = 0;
    double f_hist[opts.n_hist];
	for(; it < opts.maxit; ++it){
		if(opts.n_hist > 0){
			f_hist[it % opts.n_hist] = fx;
		}
		if(opts.debug){
			// dump current status:
			cout << endl << "At iteration " << it << ": " << endl << "x = ";
			for(size_t k=0; k<n; ++k){
				printf(" %8.4g", x.get(vpids[k]));
			}
			cout << " --> " << fx << endl << "g = ";
			for(size_t k=0; k<n; ++k){
				printf(" %8.4g", grad.get(vpids[k]));
			}
			cout << endl;
			cout << endl << "h = " << endl;
			for(size_t k=0; k<n; ++k){
				for(size_t l=0; l<n; ++l){
					printf(" %8.4g", inverse_hessian(k, l));
				}
				cout << endl;
			}
		}

		// the search direction is given by    -inverse_hessian * grad
		ParValues dir(zero);
		double dir_norm = 0.0;
		if(opts.debug) cout << "search dir = ";
		for(size_t i=0; i<n; ++i){
			double dir_i = 0.0;
			for(size_t j=0; j<n; ++j){
				dir_i -= inverse_hessian(i, j) * grad.get(vpids[j]);
			}
			dir_norm += dir_i * dir_i;
			dir.set(vpids[i], dir_i);
			if(opts.debug) printf(" %8.4g", dir_i);
		}
		if(opts.debug)cout << endl;
		dir_norm = sqrt(dir_norm);
		if(dir_norm < opts.par_eps){
			if(opts.debug)cout << "estimated minimum close enough: " << dir_norm << "; stopping iteration." << endl;
			break;
		}

		if(it >= opts.n_hist && opts.n_hist > 0){
			double min_f = inf, max_f = -inf;
			for(int k=0; k<opts.n_hist; ++k){
				min_f = min(min_f, f_hist[k]);
				max_f = max(max_f, f_hist[k]);
			}
			if(max_f - min_f < opts.f_eps){
				if(opts.debug)cout << "f_eps criterion fulfilled, stopping iteration." << endl;
				break;
			}
		}
		if(opts.debug)cout << "doing linesearch now" << endl;
		ParValues next_x = ls->do_linesearch(f, r, x, dir);
		if(opts.debug){
			cout << "linesearch proposes:" << endl << "x =";
			for(size_t k=0; k<n; ++k){
				printf(" %8.4g", next_x.get(vpids[k]));
			}
			cout << endl;
		}

		// do the update to inverse_hessian:
		ParValues new_grad(zero);
		double new_fx = eval_with_grad(f, r, next_x, new_grad);
		// calculate z = (Delta x - inverse_hessian * (new_grad - grad))
		ParValues z(zero);
		ParValues dx(zero);
		for(size_t i=0; i<n; ++i){
			double z_i = next_x.get(vpids[i]) - x.get(vpids[i]);
			dx.set(vpids[i], z_i);
			for(size_t j=0; j<n; ++j){
				z_i -= inverse_hessian(i, j) * (new_grad.get(vpids[j]) - grad.get(vpids[j]));
			}
			z.set(vpids[i], z_i);
		}
		// calculate the denominator, z * (new_grad - grad):
		double denom = 0.0;
		double norm_dx = 0.0;
		double norm_z = 0.0;
		double norm_graddiff = 0.0;
		for(size_t i=0; i<n; ++i){
			norm_dx += pow(dx.get(vpids[i]), 2);
			norm_z += pow(z.get(vpids[i]), 2);
			double graddiff = new_grad.get(vpids[i]) - grad.get(vpids[i]);
			norm_graddiff += graddiff * graddiff;
			denom += z.get(vpids[i]) * graddiff;
		}
		norm_dx = sqrt(norm_dx);
		if(norm_dx < opts.par_eps){
			if(opts.debug){
				cout << "actual dx was smaller than par_eps, stopping iteration." << endl;
				break;
			}
		}
		norm_z = sqrt(norm_z);
		norm_graddiff = sqrt(norm_graddiff);
		if(opts.debug)cout << "denom = " << denom << "; |z| = " << norm_z << "; |g_k+1 - g_k| = " << norm_graddiff << endl;
		// TODO: there should be rather a relative check on denom ...
		if(fabs(denom) > s_eps * sqrt(norm_graddiff * norm_z)){
			// calculate the update to inverse_hessian, i.e. add z * z^T / denom:
			for(size_t i=0; i<n; ++i){
				for(size_t j=0; j< n; ++j){
					inverse_hessian(i,j) += z.get(vpids[i]) * z.get(vpids[j]) / denom;
				}
			}
		}
		else{
			if(opts.debug) cout << "denom is too small for update; continuing with next iteration directly " << endl;
		}
		// prepare for next iteration:
		x.set(next_x);
		fx = new_fx;
		grad.set(new_grad);
	}
	if(it==opts.maxit){
		throw MinimizationException("maximum number of iterations reached");
	}
	MinimizationResult res;
	res.fval = fx;
	res.values = x;
	// TODO: covariance ... (but: is only very approximate).
	return res;
}


MinimizationResult newton_minimizer::minimize(const theta::Function & f, const theta::ParValues & start,
                                               const theta::ParValues & step, const std::map<theta::ParId, std::pair<double, double> > & ranges){
	const ParIds & pids = f.get_parameters();
	if(opts.debug){
		cout << "start / steps: ";
		for(ParIds::const_iterator it=pids.begin(); it!=pids.end(); ++it){
			cout << " " <<  start.get(*it) << "+-" << step.get(*it);
		}
		cout << endl;
	}

	MinimizationProblem mp(f, start, step, ranges);
	Matrix m = mp.get_step_matrix();
	if(opts.debug)cout << "maxtrix(0,0) = " << m(0,0) << endl;
	return minimize2(mp);
}
