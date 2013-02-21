#include <algorithm>
#include <vector>
#include <cassert>
#include <cmath>

#include "interface/minimizer.hpp"


struct options{
  double step_factor;
  // convergence criterion: the maximum difference of function values of the simplex are <= f_eps  
  double f_eps;
  // parameters for the update: reflection, contraction, expansion and shrinkage parameters:
  double alpha, beta, gamma, delta;
  
  // maximum number of iterations:
  size_t max_iter;
  
  bool debug;

  options(): step_factor(2), f_eps(1e-5), alpha(1.0), beta(0.5), gamma(2.0), delta(0.5), max_iter(50000), debug(false){}
};


class simplex_minimizer: public theta::Minimizer{
    options opts;
public:
    simplex_minimizer(const theta::Configuration & cfg);
    simplex_minimizer(const options & opts_);
    virtual theta::MinimizationResult minimize(const theta::Function & f, const theta::ParValues & start,
                                               const theta::ParValues & step, const std::map<theta::ParId, std::pair<double, double> > & ranges);
    
};

