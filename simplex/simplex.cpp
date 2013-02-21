#include "simplex.hpp"
#include <iostream>

#include "interface/plugin.tcc"
#include "interface/phys.hpp"

using namespace std;
using namespace theta;

/*// derive from this class and implement the objective function:
class function{
public:
    virtual double operator()(const double * x) const = 0;
    virtual size_t get_ndim() const = 0;
    virtual ~function();
};


function::~function(){}*/


// vector of doubles with fast additions, etc.:
class vec{
    vector<double> data;
public:
    
    explicit vec(size_t ndim = 0): data(ndim){
    }
    
    void set(size_t size, const double * value){
        data.resize(size);
        if(value!=0){
            std::copy(value, value + size, data.begin());
        }
    }
    
    double & operator[](size_t i){
        return data[i];
    }
    
    double operator[](size_t i) const{
        return data[i];
    }
    
    void operator+=(const vec & rhs){
        assert(rhs.data.size() == data.size());
        for(size_t i=0; i<data.size(); ++i){
            data[i] += rhs.data[i];
        }
    }
    
    void operator*=(double coeff){
        for(size_t i=0; i<data.size(); ++i){
            data[i] *= coeff;
        }
    }
    
    void zeros(){
        for(size_t i=0; i<data.size(); ++i){
            data[i] = 0.0;
        }
    }
    
    void add_with_coeff(double coeff, const vec & rhs){
        assert(rhs.data.size() == data.size());
        for(size_t i=0; i<data.size(); ++i){
            data[i] += coeff * rhs.data[i];
        }
    }
    
    const double * get_data()const{
        return &(data[0]);
    }
    
    bool equals(const vec & rhs, double eps) const{
        if(data.size() != rhs.data.size()) return false;
        for(size_t i=0; i<data.size(); ++i){
            if(fabs(data[i] - rhs.data[i]) > eps) return false;
        }
        return true;
    }
    
};



class constrained_f{
    const theta::Function & f;
    const std::map<theta::ParId, std::pair<double, double> > & ranges;
    ParIds nonfixed_pids;
    ParValues fixed_values;
public:
    constrained_f(const theta::Function & f_, const std::map<theta::ParId, std::pair<double, double> > & ranges_): f(f_), ranges(ranges_){
        // what really counts are the parameters in the function; ranges might contain more pids
        const ParIds & pids = f.get_parameters();
        for(ParIds::const_iterator it=pids.begin(); it!=pids.end(); ++it){
            std::map<theta::ParId, std::pair<double, double> >::const_iterator itt=ranges.find(*it);
            if(itt == ranges.end()){
                throw invalid_argument("ranges does not contain all pids of the function");
            }
            if(itt->second.first == itt->second.second){
                fixed_values.set(*it, itt->second.second);
            }
            else{
                nonfixed_pids.insert(*it);
            }
        }
    }
    
    double operator()(const double * x){
        size_t i=0;
        for(ParIds::const_iterator it=nonfixed_pids.begin(); it!=nonfixed_pids.end(); ++i, ++it){
            const pair<double, double> & range = ranges.find(*it)->second;
            if(x[i] < range.first || x[i] > range.second){
                //cout << "outside of range: " << x[i] << "    range is[" << range.first << ", " << range.second << "]" << endl;
                return std::numeric_limits<double>::infinity();
            }
        }
        ParValues values(x, nonfixed_pids);
        values.set(fixed_values);
        return f(values);
    }
    
    size_t get_ndim()const{
        return nonfixed_pids.size();
    }
    
    const ParIds & get_parameters() const{
        return nonfixed_pids;
    }
    
};

simplex_minimizer::simplex_minimizer(const options & opts_): opts(opts_){

}

simplex_minimizer::simplex_minimizer(const theta::Configuration & cfg){
    if(cfg.setting.exists("alpha")){
        opts.alpha = cfg.setting["alpha"];
        if(opts.alpha <= 0){
            throw ConfigurationException("alpha > 0 required");
        }
    }
    if(cfg.setting.exists("beta")){
        opts.beta = cfg.setting["beta"];
        if(opts.beta <= 0 || opts.beta >= 1){
            throw ConfigurationException("0 < beta < 1 required");
        }
    }
    if(cfg.setting.exists("gamma")){
        opts.gamma = cfg.setting["gamma"];
        if(opts.gamma <= 1){
            throw ConfigurationException("gamma > 1 required");
        }
    }
    if(cfg.setting.exists("delta")){
        opts.delta = cfg.setting["delta"];
        if(opts.delta <= 0 || opts.delta >= 1){
            throw ConfigurationException("0 < delta < 1 required");
        }
    }
    if(cfg.setting.exists("max_iter")){
        opts.max_iter = static_cast<unsigned int>(cfg.setting["max_iter"]);
    }
    if(cfg.setting.exists("step_factor")){
        opts.step_factor = cfg.setting["step_factor"];
        if(opts.step_factor <= 0){
            throw ConfigurationException("step_factor <=0 not allowed");
        }
    }
    if(cfg.setting.exists("f_eps")){
        opts.f_eps = cfg.setting["f_eps"];
        if(opts.f_eps <= 0){
            throw ConfigurationException("f_eps <=0 not allowed");
        }
    }
}


theta::MinimizationResult simplex_minimizer::minimize(const theta::Function & theta_f, const theta::ParValues & start,
                                                      const theta::ParValues & step, const std::map<theta::ParId, std::pair<double, double> > & ranges){
    
    // check consistency of the input values:
    {
        const ParIds & pids = theta_f.get_parameters();
        ParIds ranges_pids;
        for(std::map<theta::ParId, std::pair<double, double> >::const_iterator rit=ranges.begin(); rit!=ranges.end(); ++rit){
            ranges_pids.insert(rit->first);
            if(rit->second.second < rit->second.first) throw invalid_argument("invalid range");
            if(start.get(rit->first) > rit->second.second || start.get(rit->first) < rit->second.first) throw invalid_argument("invalid start point (not in range)");
            double s = step.get(rit->first);
            if(s < 0.0) throw invalid_argument("negative step size given");
            if(s==0 && (rit->second.second != rit->second.first)) throw invalid_argument("zero step size given, but range is finite");
        }
        // note: we are a little more restrictive here than we have to be: it would be enough to ask that pids is a subset of ranges_pids instead of equality
        if(!start.contains_all(pids) || !step.contains_all(pids) || !(pids == ranges_pids)){
            throw invalid_argument("inconsistent parameters to simplex_minimizer");
        }
    }
    
    constrained_f f(theta_f, ranges);
    const size_t ndim = f.get_ndim();
    const size_t npoints = ndim + 1;

    vector<vec> current_simplex(npoints);
    // function values, in the same ordering as current_simplex:
    vector<double> fvals(npoints);
    const ParIds & pids = f.get_parameters();
    //cout << "setup start" << endl;
    for(size_t i=0; i<npoints; ++i){
        current_simplex[i].set(ndim, 0);
        size_t j=0;
        for(ParIds::const_iterator it=pids.begin(); j<ndim; ++j, ++it){
            current_simplex[i][j] = start.get(*it);
            if(i==j){
                double s = step.get(*it) * opts.step_factor;
                double newval = current_simplex[i][j] + s;
                const pair<double, double> & range = ranges.find(*it)->second;
                /*if(sign * s > 0){
                    if(newval > range.second){
                        newval = range.second;
                    }
                }
                else{
                    if(newval < range.first){
                        newval = range.first;
                    }
                }*/
                if(newval > range.second){
                    newval = range.second;
                    /*newval = current_simplex[i][j] - s;
                    if(newval < range.first){
                        if(range.second - current_simplex[i][j] > current_simplex[i][j] - range.first){
                            newval = range.second;
                        }
                        else{
                            newval = range.first;
                        }
                    }*/
                }
                //assert(current_simplex[i][j] != newval && newval >= range.first && newval <= range.second);
                current_simplex[i][j] = newval;
            }
            // for the last point, use     start - step*step_factor    in each direction:
            else if(i==npoints-1){
                double s = step.get(*it) * opts.step_factor;
                double newval = current_simplex[i][j] - s;
                const pair<double, double> & range = ranges.find(*it)->second;
                if(newval < range.first){
                    newval = range.first;
                }
                current_simplex[i][j] = newval;
            }
            //cout << i << "," << j << ":" << current_simplex[i][j] << endl;
        }
        assert(j==ndim);
        // set fval:
        //cout << "get fval " << i << endl;
        fvals[i] = f(current_simplex[i].get_data());
        //cout << "get fval end " << i << endl;
    }
    
    //cout << "setup end" << endl;
    // indices iw,isw,ib are of the worst, second worst and best vertex, resp.:
    size_t iw, isw, ib;
    
    // initial function evaluation and determination of iw, isw, ib. Note: probably not the most efficient method,
    // but it's save in case of equalities ...
    ib = iw = isw = npoints;
    for(size_t i=0; i<npoints; ++i){
        if(iw == npoints || fvals[i] > fvals[iw]) iw = i;
    }
    for(size_t i=0; i<npoints; ++i){
        if(i==iw)continue;
        if(ib == npoints || fvals[i] < fvals[ib]) ib = i;
        if(isw == npoints || fvals[i] > fvals[isw]) isw = i;
    }
    assert(ib!=iw && iw!=isw);
    // note: isw==ib is allowed.
    
    // compute centroid:
    vec c(ndim);
    for(size_t i=0; i<npoints; ++i){
        if(i==iw)continue;
        c.add_with_coeff(1.0 / (npoints - 1), current_simplex[i]);
    }
    
    // in each iteration, make sure to update
    // a) current_simplex, fvals
    // b) iw, isw, ib
    // c) c
    
    for(size_t it=0; it<opts.max_iter; ++it){
        // convergence test:
        if(fabs(fvals[ib] - fvals[iw]) < opts.f_eps){
            theta::MinimizationResult result;
            result.fval = fvals[ib];
            size_t j=0;
            for(ParIds::const_iterator it=pids.begin(); it!=pids.end(); ++j, ++it){
                result.values.set(*it, current_simplex[ib][j]);
            }
            return result;
        }
        
        bool debug = opts.debug;
        if(debug){
			cout << endl << "iteraton " << it << " iw = " << iw << "; isw = " << isw << "; ib = " << ib << endl;
			// consistency checks (can be removed once in production):
			vec c2(c);
			c2.zeros();
			for(size_t i=0; i<npoints; ++i){
				if(i==iw)continue;
				c2.add_with_coeff(1.0 / (npoints - 1), current_simplex[i]);
			}
			assert(c2.equals(c, 1e-7));

			assert(iw!=ib && iw!=isw);
			//cout << "iw: " << iw << "; ib: " << ib  << "; isw: " << isw << endl;
			for(size_t i=0; i<npoints; ++i){
				cout << "[" << i << "]: ";
				for(size_t j=0; j<ndim; ++j){
					cout << current_simplex[i][j] << ", ";
				}
				cout << "--> " << fvals[i] << endl;
				if(i==iw)continue;
				assert(fvals[isw] >= fvals[i]);
			}
			cout << "[c]: ";
			for(size_t j=0; j<ndim; ++j){
				cout << c[j] << ", ";
			}
			cout << endl;
			assert(fvals[isw] <= fvals[iw]);
			for(size_t i=0; i<npoints; ++i){
				assert(fvals[ib] <= fvals[i]);
				assert(fvals[iw] >= fvals[i]);
			}
        } // if debug

        //reflect worst point at the centroid:
        vec xr(c);
        xr *= 1 + opts.alpha;
        xr.add_with_coeff(-opts.alpha, current_simplex[iw]);
        const double fr = f(xr.get_data());
        // accept xr if it is better than the second worst, but not the best:
        if(fr < fvals[isw] && fr >= fvals[ib]){
            if(debug)cout << "case 1: accepting xr as it's better than second worst" << endl;
            // a)
            const size_t new_index = iw;
            current_simplex[iw] = xr;
            fvals[iw] = fr;
            // b)
            // new worst is what was second worst so far:
            
            iw = isw;
            // find new second worst; make sure to use an initial value isw != iw:
            isw = npoints;
            for(size_t i=0; i< npoints; ++i){
                if(i==iw) continue;
                if(isw==npoints || fvals[i] > fvals[isw]) isw = i;
            }
            // ib is unchanged.
            // c) the new centroid contains the new point, and subtract what is now the worst:
            c.add_with_coeff(-1.0 / (npoints - 1), current_simplex[iw]);
            c.add_with_coeff(1.0 / (npoints - 1), current_simplex[new_index]);
        }
        // if fr is better than the current best: expand
        else if(fr < fvals[ib]){
        	if(debug)cout << "case 2: xr is better than current best" << endl;
            // compute the expansion point:
            vec xe(c);
            xe *= 1 - opts.gamma;
            xe.add_with_coeff(opts.gamma, xr);
            const double fe = f(xe.get_data());
            // depending on fe < fr, accept either xe or xr and terminate iteration:
            if(fe < fr){
            	if(debug)cout << "case 2a: expanded xr (=xe) is even better" << endl;
                // accept xe. a)
                current_simplex[iw] = xe;
                fvals[iw] = fe;
            }
            else{
            	if(debug)cout << "case 2b: expanded xr (=xe) wasn't better, using xr" << endl;
                // accept xr. a)
                current_simplex[iw] = xr;
                fvals[iw] = fr;
            }
            ib = iw;
            iw = isw;
            //look for new isw:
            isw = npoints;
			for(size_t i=0; i< npoints; ++i){
				if(i==iw) continue;
				if(isw==npoints || fvals[i] > fvals[isw]) isw = i;
			}
            // c. remove new worst from simplex:
			c.add_with_coeff(-1.0 / (npoints - 1), current_simplex[iw]);
            c.add_with_coeff(1.0 / (npoints - 1), current_simplex[ib]);
        }
        else{
        	if(debug)cout << "case 3: xr is worse than second worst" << endl;
            // compute the contraction point xc:
            vec xc(c);
            xc *= 1 - opts.beta;
            xc.add_with_coeff(opts.beta, current_simplex[iw]);
            const double fc = f(xc.get_data());
            if(fc <= fvals[iw]){
            	if(debug)cout << "case 3a: contracted point was better than current worst" << endl;
                // accept xc and terminate iteration:
                // a)
                const size_t new_index = iw;
                current_simplex[iw] = xc;
                fvals[iw] = fc;
                // b)
                // maybe the new point is also better than the second worst:
                if(fc < fvals[isw]){
                    iw = isw;
                    // find new second worst
                    isw = npoints;
                    for(size_t i=0; i< npoints; ++i){
                        if(i==iw)continue;
                        if(isw==npoints || fvals[i] > fvals[isw]) isw = i;
                    }
                    if(fvals[new_index] < fvals[ib]) ib = new_index;
                    // c) the new centroid contains the new point, and subtract what is now the worst:
                    c.add_with_coeff(1.0 / (npoints - 1), current_simplex[new_index]);
                    c.add_with_coeff(-1.0 / (npoints - 1), current_simplex[iw]);
                }
            }
            else{
            	if(debug)cout << "case 3b: contracted point was worse than current worst; shrinking" << endl;
                //shrink all points toward the best:
                for(size_t i=0; i<npoints; ++i){
                    if(i==ib)continue;
                    // a)
                    current_simplex[i] *= opts.delta;
                    current_simplex[i].add_with_coeff(1 - opts.delta, current_simplex[ib]);
                    fvals[i] = f(current_simplex[i].get_data());
                }
                // b)
                ib = iw = isw = npoints;
                for(size_t i=0; i<npoints; ++i){
                    if(iw == npoints || fvals[i] > fvals[iw]) iw = i;
                }
                for(size_t i=0; i<npoints; ++i){
                    if(i==iw)continue;
                    if(ib == npoints || fvals[i] < fvals[ib]) ib = i;
                    if(isw == npoints || fvals[i] > fvals[isw]) isw = i;
                }
                assert(ib!=iw && iw!=isw);
                // c)
                c.zeros();
                for(size_t i=0; i<npoints; ++i){
                    if(i==iw)continue;
                    c.add_with_coeff(1.0 / (npoints - 1), current_simplex[i]);
                }
            }
        }
    }
    // error: too many iterations...
    throw MinimizationException("too many iterations");
}

REGISTER_PLUGIN(simplex_minimizer)
