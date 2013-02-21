#include "interface/minimizer.hpp"
#include <memory>

class Linesearch;

class newton_minimizer: public theta::Minimizer{
public:
    struct options{
    	int maxit;
    	double par_eps; // iteration stops if minimum is nearer than that, measured in parameter space
    	double f_eps;
    	int n_hist; // alternatively, iteration stops once it is found that f has not changed more than f_eps over the last n_hist iterations
    	bool debug;
    	options(): maxit(1000), par_eps(1e-8), f_eps(1e-3), n_hist(4), debug(false){}
    };


    //newton_minimizer(const options & opts_);
    newton_minimizer(const options & opts_);
    virtual theta::MinimizationResult minimize(const theta::Function & f, const theta::ParValues & start,
                                               const theta::ParValues & step, const std::map<theta::ParId, std::pair<double, double> > & ranges);
    virtual theta::MinimizationResult minimize2(const theta::MinimizationProblem & mp);

private:
    std::auto_ptr<Linesearch> ls;
    options opts;

};


