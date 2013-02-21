#include "interface/minimizer.hpp"
#include "interface/phys.hpp"
#include "interface/plugin.tcc"

REGISTER_PLUGIN_BASETYPE(theta::Minimizer);

using namespace theta;
using namespace std;

void MinimizationResult::operator=(const MinimizationResult& rhs){
    fval = rhs.fval;
    values.set(rhs.values);
    errors_plus.set(rhs.errors_plus);
    errors_minus.set(rhs.errors_minus);
    covariance = rhs.covariance;
}

MinimizationProblem::MinimizationProblem(const theta::Function & f_, const theta::ParValues & start_,
                                         const theta::ParValues & step_, const std::map<theta::ParId, std::pair<double, double> > & ranges_): f(f_), start(start_), step(step_), ranges(ranges_){
     const ParIds & pars = f.get_parameters();
     matrix.reset(pars.size(), pars.size());
     size_t i=0;
     for(ParIds::const_iterator pit=pars.begin(); pit!=pars.end(); ++pit, ++i){
         matrix(i,i) = pow(step.get(*pit), 2);
     }
     check_consistency();
}

MinimizationProblem::MinimizationProblem(const theta::Function & f_, const theta::ParValues & start_,
                                         const theta::Matrix & matrix_, const std::map<theta::ParId, std::pair<double, double> > & ranges_): f(f_), start(start_), ranges(ranges_), matrix(matrix_){
    const ParIds & pars = f.get_parameters();
    theta_assert(pars.size() == matrix.get_n_rows() && pars.size() == matrix.get_n_cols());
    size_t i=0;
    for(ParIds::const_iterator pit=pars.begin(); pit!=pars.end(); ++pit, ++i){
        theta_assert(matrix(i,i) >= 0.0);
        step.set(*pit, sqrt(matrix(i,i)));
    }
    check_consistency();
}

// check the consistency of the parameters
void MinimizationProblem::check_consistency() const{
    // all sets of ParId are the same:
    const ParIds & f_pars = f.get_parameters();
    theta_assert(start.contains_all(f_pars));
    theta_assert(step.contains_all(f_pars));
    ParIds ranges_pids;
    theta_assert(f_pars.size() == matrix.get_n_rows());
    theta_assert(f_pars.size() == matrix.get_n_cols());
    // check that start is within ranges, and that steps are non-negative:
    for(ParIds::const_iterator it=f_pars.begin(); it!=f_pars.end(); ++it){
        double val = start.get(*it);
        theta_assert(std::isfinite(val));
        const pair<double, double> & range = ranges.get(*it);
        theta_assert(range.first <= range.second);
        theta_assert(val <= range.second && val >= range.first);
        double st = step.get(*it);
        theta_assert(st >= 0.0 && std::isfinite(st));
    }
    for(size_t i=0; i<matrix.get_n_rows(); ++i){
        double me = matrix(i,i);
        theta_assert(std::isfinite(me));
        theta_assert(me >= 0.0);
    }
}

MinimizationResult Minimizer::minimize2(const MinimizationProblem & mp){
    return minimize(mp.get_f(), mp.get_start(), mp.get_steps(), mp.get_ranges());
}

Minimizer::~Minimizer(){}
