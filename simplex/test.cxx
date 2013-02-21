#include <boost/test/unit_test.hpp>

#include "interface/matrix.hpp"
#include "interface/phys.hpp"
#include "interface/variables.hpp"
#include "interface/random.hpp"

#include "simplex/simplex.hpp"
#include "simplex/newton.hpp"

#include <vector>
#include <iostream>

// define some test functions:

using namespace theta;
using namespace std;

namespace{
	void dump(ostream & out, const ParValues & vals, const VarIdManager & vm){
		ParIds pids = vm.get_all_parameters();
		for(ParIds::const_iterator it=pids.begin(); it!=pids.end(); ++it){
			try{
				double value = vals.get(*it);
				std::string name = vm.get_name(*it);
				out << " " << name << " = " << value;
			}
			catch(invalid_argument &){
				//ignore
			}
		}
	}
}

class quadratic_function: public theta::Function{
	Matrix hesse;
	ParValues x0;
	mutable size_t n_eval;
public:
	quadratic_function(const ParIds & pids, const ParValues & x0_, const Matrix & hesse_): hesse(hesse_), x0(x0_), n_eval(0){
		theta_assert(hesse.get_n_cols() == hesse.get_n_rows());
		theta_assert(hesse.get_n_cols() == pids.size());
		theta_assert(x0.contains_all(pids));
		par_ids = pids;
	}

	size_t get_n_eval() const{
		return n_eval;
	}
	void reset_n_eval(){
		n_eval = 0;
	}

	using Function::operator();

	virtual double operator()(const ParValues & x) const{
		++n_eval;
		const size_t n = par_ids.size();
		std::vector<double> xx0(n); // x-x0
		size_t i=0;
		for(ParIds::const_iterator it=par_ids.begin(); it!=par_ids.end(); ++it, ++i){
			xx0[i] = x.get(*it) - x0.get(*it);
			theta_assert(std::isfinite(xx0[i]));
		}
		double result =0.;
		for(i=0; i<n; ++i){
			for(size_t j=0; j<n; ++j){
				result += 0.5 * hesse(i,j) * xx0[j] * xx0[i];
			}
		}
		theta_assert(!std::isnan(result));
		return result;
	}

};

namespace{
const double inf = std::numeric_limits<double>::infinity();
}


BOOST_AUTO_TEST_SUITE(simplex_tests)

BOOST_AUTO_TEST_CASE(test1d){
	VarIdManager vm;
	ParId pid = vm.create_par_id("p0");
	ParValues x0;
	x0.set(pid, 0.0);
	ParIds pids;
	pids.insert(pid);

	Matrix hesse(1,1);
	hesse(0,0) = 1.0;

	quadratic_function f(pids, x0, hesse);
	options opts;
	simplex_minimizer min(opts);

	ParValues start, step;
	start.set(pid, 1.0);
	step.set(pid, 1.0);
	std::map<ParId, pair<double, double> > ranges;
	ranges[pid] = make_pair(-inf, inf);
	MinimizationResult minres = min.minimize(f, start, step, ranges);
	BOOST_CHECK_LT(minres.fval, opts.f_eps);
	BOOST_CHECK_LT(f.get_n_eval(), 100);
	cout << "f(" << minres.values.get(pid) << ") = " << minres.fval << " (evals: " << f.get_n_eval() << ")"<< endl;

	// step size small:
	start.set(pid, 10.0);
	step.set(pid, 0.1);
	f.reset_n_eval();
	minres = min.minimize(f, start, step, ranges);
	BOOST_CHECK_LT(minres.fval, opts.f_eps);
	BOOST_CHECK_LT(f.get_n_eval(), 100);
	cout << "f(" << minres.values.get(pid) << ") = " << minres.fval << " (evals: " << f.get_n_eval() << ")"<< endl;

	// step size large:
	start.set(pid, 0.1);
	step.set(pid, 100.0);
	f.reset_n_eval();
	minres = min.minimize(f, start, step, ranges);
	BOOST_CHECK_LT(minres.fval, opts.f_eps);
	BOOST_CHECK_LT(f.get_n_eval(), 100);
	cout << "f(" << minres.values.get(pid) << ") = " << minres.fval << " (evals: " << f.get_n_eval() << ")"<< endl;
}


BOOST_AUTO_TEST_CASE(test10d){
	VarIdManager vm;
	std::vector<ParId> vpids;
	ParIds pids;
	const size_t ndim = 10;
	ParValues x0, start, step, x;
	Matrix hesse(ndim,ndim);
	map<ParId, pair<double, double> > ranges;
	for(size_t i=0; i<ndim; ++i){
		stringstream ss;
		ss << "p" << i;
		ParId pid = vm.create_par_id(ss.str());
		vpids.push_back(pid);
		x0.set(pid, 1.0);
		start.set(pid, 1.0);
		step.set(pid, 1.0);
		x.set(pid, 0.0);
		hesse(i,i) = 1.0;
		pids.insert(pid);
		ranges[pid] = make_pair(-inf, inf);
	}

	quadratic_function f(pids, x0, hesse);
	options opts;
	cout << "alpha = " << opts.alpha << "; beta = " << opts.beta << "; gamma = " << opts.gamma << "; delta = " << opts.delta << endl;
	simplex_minimizer min(opts);

	cout << "f(x0) = " << f(x0) << endl;
	cout << "f(0) = " << f(x) << endl;
	cout << "f(start) = " << f(start) << endl;

	// starting at the minimum:
	MinimizationResult minres = min.minimize(f, start, step, ranges);
	BOOST_CHECK_LT(minres.fval, 2*opts.f_eps);
	//BOOST_CHECK_LT(f.get_n_eval(), 100);
	cout << "minimum found at: ";
	dump(cout, minres.values, vm);
	cout << endl;
	cout << "f = " << minres.fval << " (evals: " << f.get_n_eval() << ")"<< endl;
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(newton)

BOOST_AUTO_TEST_CASE(test1d){
	BOOST_CHECKPOINT("test1d entry");
	VarIdManager vm;
	ParId pid = vm.create_par_id("p0");
	ParValues x0;
	x0.set(pid, 0.0);
	ParIds pids;
	pids.insert(pid);

	Matrix hesse(1,1);
	hesse(0,0) = 1.0;

	BOOST_CHECKPOINT("test1d");
	quadratic_function f(pids, x0, hesse);
	newton_minimizer::options opts;
	newton_minimizer min(opts);

	BOOST_CHECKPOINT("test1d");
	ParValues start, step;
	start.set(pid, 1.0);
	step.set(pid, 1.0);
	std::map<ParId, pair<double, double> > ranges;
	ranges[pid] = make_pair(-inf, inf);
	MinimizationResult minres = min.minimize(f, start, step, ranges);
	BOOST_CHECK(minres.fval < opts.f_eps || minres.values.get(pid) < opts.par_eps);
	BOOST_CHECK_LT(f.get_n_eval(), 100);
	//cout << "f(" << minres.values.get(pid) << ") = " << minres.fval << " (evals: " << f.get_n_eval() << ")"<< endl;

	BOOST_CHECKPOINT("test1d");
	// step size small:
	start.set(pid, 10.0);
	step.set(pid, 0.1);
	f.reset_n_eval();
	minres = min.minimize(f, start, step, ranges);
	BOOST_CHECK_LT(minres.fval, opts.f_eps);
	BOOST_CHECK_LT(f.get_n_eval(), 100);
	//cout << "f(" << minres.values.get(pid) << ") = " << minres.fval << " (evals: " << f.get_n_eval() << ")"<< endl;

	BOOST_CHECKPOINT("test1d");
	// step size large:
	start.set(pid, 0.1);
	step.set(pid, 100.0);
	f.reset_n_eval();
	minres = min.minimize(f, start, step, ranges);
	BOOST_CHECK_LT(minres.fval, opts.f_eps);
	BOOST_CHECK_LT(f.get_n_eval(), 100);
	//cout << "f(" << minres.values.get(pid) << ") = " << minres.fval << " (evals: " << f.get_n_eval() << ")"<< endl;
}

BOOST_AUTO_TEST_CASE(test10d){
	VarIdManager vm;
	std::vector<ParId> vpids;
	ParIds pids;
	const size_t ndim = 10;
	ParValues x0, start, step, x;
	Matrix hesse(ndim,ndim);
	map<ParId, pair<double, double> > ranges;
	for(size_t i=0; i<ndim; ++i){
		stringstream ss;
		ss << "p" << i;
		ParId pid = vm.create_par_id(ss.str());
		vpids.push_back(pid);
		x0.set(pid, 1.0);
		start.set(pid, 1.0);
		step.set(pid, 1.0);
		x.set(pid, 0.0);
		hesse(i,i) = 1.0;
		pids.insert(pid);
		ranges[pid] = make_pair(-inf, inf);
	}
	quadratic_function f(pids, x0, hesse);

	newton_minimizer::options opts;
	newton_minimizer min(opts);

	/*cout << "f(x0) = " << f(x0) << endl;
	cout << "f(0) = " << f(x) << endl;
	cout << "f(start) = " << f(start) << endl;*/

	// starting at the minimum:
	f.reset_n_eval();
	MinimizationResult minres = min.minimize(f, start, step, ranges);
	BOOST_CHECK_LT(minres.fval, opts.f_eps);
	BOOST_CHECK_LT(f.get_n_eval(), 500);
	/*cout << "minimum found at: ";
	dump(cout, minres.values, vm);
	cout << endl;
	cout << "f = " << minres.fval << " (evals: " << f.get_n_eval() << ")"<< endl;*/

	// starting somewhat off, with too large stepsize:
	for(size_t i=0; i<ndim; ++i){
		start.set(vpids[i], 1.0 + i);
		step.set(vpids[i], 10.0);
	}
	f.reset_n_eval();
	minres = min.minimize(f, start, step, ranges);
	BOOST_CHECK_LT(minres.fval, opts.f_eps);
	BOOST_CHECK_LT(f.get_n_eval(), 500);
	/*cout << "minimum found at: ";
	dump(cout, minres.values, vm);
	cout << endl;
	cout << "f = " << minres.fval << " (evals: " << f.get_n_eval() << ")"<< endl;*/
}

// 30 dimensions with different "scales", not aligned to axes
BOOST_AUTO_TEST_CASE(test30d_skewed){
	const size_t ndim = 10;
	Matrix hesse(ndim, ndim), covariance(ndim, ndim);
	// make the covariance matrix a 30*30 matrix by transforming the diagonal matrix of eigenvalues lambda with
	// a householder transform trafo: cov = trafo * lambda * trafo
	Matrix lambda(ndim, ndim), trafo(ndim, ndim);
	vector<double> v(ndim);
	auto_ptr<RandomSource> rnd_src(new RandomSourceTaus());
	Random rnd(rnd_src);
	double v_norm2 = 0.0;
	for(size_t i=0; i < ndim; ++i){
		lambda(i, i) = pow(1.0 + i, 2);
		v[i] = rnd.uniform() - 0.5;
		v_norm2 += v[i] * v[i];
	}
	// construct trafo = I - 2 v v^T:
	//printf("trafo:\n");
	for(size_t i=0; i<ndim; ++i){
		for(size_t j=0; j<ndim; ++j){
			trafo(i,j) = (i==j?1:0) - 2 * v[i] * v[j] / v_norm2;
			//printf(" %10.3g", trafo(i,j));
		}
		//printf("\n");
	}
	// covariance_ij = trafo_ik * lambda_kl * trafo_lj, but lambda is diagonal, so
	// covariance_ij = trafo_ik * lambda_kk * trafo_kj
	for(size_t i=0; i<ndim; ++i){
		for(size_t j=0; j<ndim; ++j){
			for(size_t k=0; k<ndim; ++k){
				covariance(i,j) += trafo(i,k) * lambda(k,k) * trafo(k,j);
			}
		}
	}

	// test:
	VarIdManager vm;
	std::vector<ParId> vpids;
	ParIds pids;
	ParValues x0, start, step;
	map<ParId, pair<double, double> > ranges;
	for(size_t i=0; i<ndim; ++i){
		stringstream ss;
		ss << "p" << i;
		ParId pid = vm.create_par_id(ss.str());
		vpids.push_back(pid);
		x0.set(pid, 0.0);
		start.set(pid, 1.0);
		step.set(pid, 1.0);
		pids.insert(pid);
		ranges[pid] = make_pair(-inf, inf);
	}
	hesse = covariance;
	hesse.invert_cholesky();

	/*printf("covariance:\n");
	for(size_t i=0; i<ndim; ++i){
		for(size_t j=0; j<ndim; ++j){
			printf(" %10.4g", covariance(i,j));
		}
		printf("\n");
	}
	printf("true hesse:\n");
	for(size_t i=0; i<ndim; ++i){
		for(size_t j=0; j<ndim; ++j){
			printf(" %10.4g", hesse(i,j));
		}
		printf("\n");
	}*/

	quadratic_function f(pids, x0, hesse);

	newton_minimizer::options opts;
	opts.f_eps = 0.0;
	newton_minimizer min(opts);

	f.reset_n_eval();
	MinimizationResult minres = min.minimize(f, start, step, ranges);
	BOOST_CHECK_LT(minres.fval, 1e-10);
	BOOST_CHECK_LT(f.get_n_eval(), 400);

	// large step sizes seem to be much les problematic than too small ones, test that:
	f.reset_n_eval();
	for(size_t i=0; i<ndim; ++i){
		step.set( vpids[i], 100.0);
	}
	minres = min.minimize(f, start, step, ranges);
	BOOST_CHECK_LT(minres.fval, 1e-10);
	BOOST_CHECK_LT(f.get_n_eval(), 1000);

	// test whether it's faster with correct covariance:
	f.reset_n_eval();
	minres = min.minimize2(MinimizationProblem(f, start, covariance, ranges));
	BOOST_CHECK_LT(minres.fval, 1e-10);
	BOOST_CHECK_LT(f.get_n_eval(), 100);
}

BOOST_AUTO_TEST_CASE(range1d){
	VarIdManager vm;
	ParId pid = vm.create_par_id("p0");
	ParValues x0;
	x0.set(pid, 0.0);
	ParIds pids;
	pids.insert(pid);

	Matrix hesse(1,1);
	hesse(0,0) = 1.0;

	quadratic_function f(pids, x0, hesse);
	newton_minimizer::options opts;
	newton_minimizer min(opts);

	// minimize a 1d function with a range at the minimum:
	ParValues start, step;
	start.set(pid, 1.0);
	step.set(pid, 1.0);
	std::map<ParId, pair<double, double> > ranges;
	ranges[pid] = make_pair(0, inf);
	MinimizationResult minres = min.minimize(f, start, step, ranges);
	BOOST_CHECK_LT(minres.fval, 1e-10);
	BOOST_CHECK_LT(f.get_n_eval(), 20);

	// set range at 0.2
	f.reset_n_eval();
	ranges[pid] = make_pair(0.2, inf);
	minres = min.minimize(f, start, step, ranges);
	//cout << minres.values.get(pid) << endl;
	BOOST_CHECK_LT(minres.fval, 0.2 * 0.2 * (1 + 1e-7));
	BOOST_CHECK_LT(f.get_n_eval(), 50);
}


BOOST_AUTO_TEST_CASE(range10d){
	VarIdManager vm;
	std::vector<ParId> vpids;
	ParIds pids;
	const size_t ndim = 10;
	ParValues x0, start, step, x;
	Matrix hesse(ndim,ndim);
	map<ParId, pair<double, double> > ranges;
	for(size_t i=0; i<ndim; ++i){
		stringstream ss;
		ss << "p" << i;
		ParId pid = vm.create_par_id(ss.str());
		vpids.push_back(pid);
		x0.set(pid, 1.0);
		start.set(pid, 1.0);
		step.set(pid, 1.0);
		x.set(pid, 0.0);
		hesse(i,i) = 1.0;
		pids.insert(pid);
		ranges[pid] = make_pair(i==0?1.3:-inf, inf);
	}
	quadratic_function f(pids, x0, hesse);

	newton_minimizer::options opts;
	newton_minimizer min(opts);

	/*cout << "f(x0) = " << f(x0) << endl;
	cout << "f(0) = " << f(x) << endl;
	cout << "f(start) = " << f(start) << endl;*/

	// starting close to the minimum:
	f.reset_n_eval();
	start.set(vpids[0], 2.0);
	MinimizationResult minres = min.minimize(f, start, step, ranges);
	BOOST_CHECK_LT(fabs(minres.values.get(vpids[0]) - 1.3), 1e-7);
	BOOST_CHECK_LT(minres.fval, 0.3 * 0.3 + opts.f_eps);
	BOOST_CHECK_LT(f.get_n_eval(), 500);

}

BOOST_AUTO_TEST_SUITE_END()
