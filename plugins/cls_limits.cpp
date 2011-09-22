#include "interface/cfg-utils.hpp"
#include "interface/plugin.hpp"
#include "interface/histogram.hpp"
#include "interface/variables-utils.hpp"
#include "interface/variables.hpp"
#include "interface/minimizer.hpp"
#include "interface/random.hpp"
#include "interface/database.hpp"
#include "interface/main.hpp"
#include "interface/redirect_stdio.hpp"
#include "interface/phys.hpp"
#include "interface/model.hpp"
#include "interface/matrix.hpp"
#include "interface/random-utils.hpp"
#include "plugins/asimov_likelihood_widths.hpp"

#include "libconfig/libconfig.h++"

#include <fstream>
#include <iomanip>

using namespace std;
using namespace theta;
using namespace libconfig;

// the negative log likelihood to fit the function
// y = target_cls * exp(lambda * (x - limit)), with fixed target_cls and parameters lambda, limit
// to a dataset with correlated uncertainties for y.
class exp_nll: public Function{
    ParId lambda, limit;
    double target_cls;
    vector<double> x_values, y_values;
    Matrix inverse_covariance;
public:
    exp_nll(const ParId & pid_lambda, const ParId & pid_limit, double target_cls_, const vector<double> & x_values_, const vector<double> & y_values_, const Matrix & inverse_covariance_):
      lambda(pid_lambda), limit(pid_limit), target_cls(target_cls_), x_values(x_values_), y_values(y_values_), inverse_covariance(inverse_covariance_){
        par_ids.insert(limit);
        par_ids.insert(lambda);
        theta_assert(inverse_covariance.getRows() == inverse_covariance.getCols());
        theta_assert(x_values.size()==inverse_covariance.getRows() && x_values.size()==y_values.size());
    }
    
    std::auto_ptr<Function> clone() const{
        return std::auto_ptr<Function>(new exp_nll(*this));
    }
    
    double operator()(const ParValues & values) const{
        double limit_v = values.get(limit);
        double lambda_v = values.get(lambda);
        double result = 0.0;
        const size_t N = x_values.size();
        for(size_t i=0; i < N; ++i){
            for(size_t j=0; j < N; ++j){
                result += (y_values[i] - target_cls * exp(lambda_v * (x_values[i] - limit_v))) * (y_values[j] - target_cls * exp(lambda_v * (x_values[j] -limit_v))) * inverse_covariance(i,j);
            }
        }
        return 0.5 * result;
    }
};

// forwards everything to multiple other sinks
class MultiplexingProductsSink: public ProductsSink{
private:
    std::vector<boost::shared_ptr<ProductsSink> > sinks;
    std::map<Column, std::vector<Column> > sink_columns; // for each Column we return, this holds the columns of the sinks
    int next_column_id;
public:
    MultiplexingProductsSink(const std::vector<boost::shared_ptr<ProductsSink> > & sinks_): sinks(sinks_), next_column_id(0){}

    virtual Column declare_product(const ProductsSource & source, const std::string & product_name, const data_type & type){
        Column result(next_column_id++);
        for(size_t i=0; i<sinks.size(); ++i){
            sink_columns[result].push_back(sinks[i]->declare_product(source, product_name, type));
        }
        return result;
    }
    
    virtual void set_product(const Column & c, double d){
        for(size_t i=0; i<sinks.size(); ++i){
            sinks[i]->set_product(sink_columns[c][i], d);
        }
    }
    virtual void set_product(const Column & c, int k){
        for(size_t i=0; i<sinks.size(); ++i){
            sinks[i]->set_product(sink_columns[c][i], k);
        }
    }
    virtual void set_product(const Column & c, const std::string & s){
        for(size_t i=0; i<sinks.size(); ++i){
            sinks[i]->set_product(sink_columns[c][i], s);
        }
    }
    virtual void set_product(const Column & c, const Histogram1D & h){
        for(size_t i=0; i<sinks.size(); ++i){
            sinks[i]->set_product(sink_columns[c][i], h);
        }
    }
};

class SaveDoubleColumn: public ProductsSink{
private:
    string column_name;
    double value;
public:
    virtual Column declare_product(const ProductsSource & source, const std::string & product_name, const data_type & type){
        string colname = source.getName() + "__" + product_name;
        return Column(column_name == colname?1:0);
    }
    
    virtual void set_product(const Column & c, double d){
        if(c.get_id()==1) value = d;
    }
    virtual void set_product(const Column & c, int i){
    }
    virtual void set_product(const Column & c, const std::string & s){
    }
    virtual void set_product(const Column & c, const Histogram1D & h){
    }
    
    double get_value() const{
        return value;
    }
    
    explicit SaveDoubleColumn(const string & name): column_name(name), value(NAN){}
};

// contains the (numerically derived) cls values and uncertainties as function of the truth value for a fixed ts value.
class cls_vs_truth_data {
private:
    vector<double> ptruth_values, pcls_values, pcls_uncertainties_n0, pcls_uncertainties_n;
public:
    cls_vs_truth_data(const vector<double> & truth_values_, const vector<double> & cls_values_, const vector<double> & cls_uncertainties_n0_,
       const vector<double> & cls_uncertainties_n_): ptruth_values(truth_values_),
       pcls_values(cls_values_), pcls_uncertainties_n0(cls_uncertainties_n0_), pcls_uncertainties_n(cls_uncertainties_n_){
    }
    const vector<double> & truth_values() const{
        return ptruth_values;
    }
    const vector<double> & cls_values() const{
        return pcls_values;
    }
    
    // absolute unsertainty on the CLs values due to finite number of s+b toys ("n")
    const vector<double> & cls_uncertainties_n() const{
        return pcls_uncertainties_n;
    }
    
    // absolute unsertainty on the CLs values due to finite number of b only toys ("n0")
    const vector<double> & cls_uncertainties_n0() const{
        return pcls_uncertainties_n0;
    }
    void write_txt(ostream & out) const{
        out << "# truth: cls +- errorn0 +- errorn" << endl;
        for(size_t i=0; i<ptruth_values.size(); ++i){
            out << ptruth_values[i] << "   " << pcls_values[i] << "   " << pcls_uncertainties_n0[i]  << "   " << pcls_uncertainties_n[i] << endl;
        }
    }
};


void binom_with_error(size_t nominator, size_t denominator, double & p, double & p_error){
    p = 1.0 * nominator / denominator;
    p_error = max(sqrt(p*(1-p) / denominator), 1.0 / denominator);
}

// container for toy outcomes, i.e., pairs of (truth, test statistic for s+b) and (truth, test statistic for b only) tuples,
// where many truth values are the same. Note that this construction allows for the test statistic to be truth-dependent.
class truth_ts_values{
private:
    map<double, multiset<double> > truth_to_ts_sb; // signal plus background test statistic values
    map<double, multiset<double> > truth_to_ts_b;  // background only test statistic values
public:
    void add_point_b(double truth, double ts){
        truth_to_ts_b[truth].insert(ts);
    }
    // bulk insertion (I guess this has much better performance ...)
    template <class InputIterator>
    void add_points_b(double truth, InputIterator first, InputIterator last){
        truth_to_ts_b[truth].insert(first, last);
    }
    
    
    void add_point_sb(double truth, double ts){
        truth_to_ts_sb[truth].insert(ts);
    }
    template <class InputIterator>
    void add_points_sb(double truth, InputIterator first, InputIterator last){
        truth_to_ts_sb[truth].insert(first, last);
    }
    
    double get_ts_quantile(double truth, double q) const{
        theta_assert(q>=0.0 && q < 1.0);
        map<double, multiset<double> >::const_iterator it = truth_to_ts_sb.find(truth);
        if(it==truth_to_ts_sb.end()) throw NotFoundException("truth_ts_value::get_ts_quantile");
        const multiset<double> & ts_values = it->second;
        multiset<double>::const_iterator itt = ts_values.begin();
        advance(itt, q*ts_values.size());
        return *itt;
    }
    
    bool contains_truth_value(double truth) const{
        return truth_to_ts_sb.find(truth) != truth_to_ts_sb.end() && truth_to_ts_b.find(truth) != truth_to_ts_b.end();
    }
    
    size_t get_n(double truth) const{
        map<double, multiset<double> >::const_iterator it = truth_to_ts_sb.find(truth);
        if(it==truth_to_ts_sb.end()) return 0;
        return it->second.size();
    }
    
    size_t get_n0(double truth) const{
        map<double, multiset<double> >::const_iterator it = truth_to_ts_b.find(truth);
        if(it==truth_to_ts_b.end()) return 0;
        return it->second.size();
    }
    
    struct cls_info{
        double clsb, clsb_uncertainty;
        double clb, clb_uncertainty;
        double cls, cls_uncertainty_n0, cls_uncertainty_n;
    };
    
    cls_info get_cls(double ts_value, double truth_value) const{
        map<double, multiset<double> >::const_iterator it_b = truth_to_ts_b.find(truth_value);
        map<double, multiset<double> >::const_iterator it_sb = truth_to_ts_sb.find(truth_value);
        if(it_b == truth_to_ts_b.end() || it_sb == truth_to_ts_sb.end()){
            throw InvalidArgumentException("truth_ts_values::get_cls_vs_truth: truth_to_ts_sb or truth_to_b do not contain given truth value");
        }
        cls_info result;
        multiset<double>::const_iterator it_ts = it_sb->second.lower_bound(ts_value);
        double p, p_error;
        size_t nominator = distance(it_ts, it_sb->second.end());
        binom_with_error(nominator, it_sb->second.size(), p, p_error);
        result.clsb = p;
        result.clsb_uncertainty = p_error;
        it_ts = it_b->second.lower_bound(ts_value);
        double p0, p0_error;
        nominator = distance(it_ts, it_b->second.end());
        binom_with_error(nominator, it_b->second.size(), p0, p0_error);
        result.clb = p0;
        result.clb_uncertainty = p0_error;
        if(p0==1.0){
            //regularise cls_value by inserting a lower limit instead of infinity:
            result.cls = (1 - p) / (1 - p0 + p0_error);
            result.cls_uncertainty_n = p_error / (1 - p0 + p0_error);
            result.cls_uncertainty_n0 = p0_error * result.cls / (1 - p0 + p0_error);
        }
        else{
            result.cls = (1 - p) / (1 - p0);
            result.cls_uncertainty_n = p_error / (1 - p0);
            result.cls_uncertainty_n0 = p0_error * result.cls / (1 - p0);
        }
        return result;
    }
    
    // It is an error if there is no data for truth=0 or the ts value is too low (=too background-like) for truth=0.
    //
    // For the function "cls versus truth", there are two error sources for the error on cls:
    // * the finite number of toys for truth=0  (n0), correlated between different truth values
    // * the finite number of toys for the current truth value (n) which is uncorrelated
    // the two covariance matrices (where the latter is diagonal ...) are calculated and stored separately. This is helpful for the decision on
    // where to perform more toys, at truth=0, or at the intersection ...
    cls_vs_truth_data get_cls_vs_truth(double ts_value) const{
         if(truth_to_ts_b.size() != truth_to_ts_sb.size()){
             throw InvalidArgumentException("truth_ts_values::get_cls_vs_truth: truth_to_ts_sb and truth_to_b have different size");
         }
         const size_t n_truth = truth_to_ts_sb.size();
         vector<double> truth_values(n_truth+1);
         vector<double> cls_values(n_truth+1);
         vector<double> cls_uncertainties_n(n_truth+1);
         vector<double> cls_uncertainties_n0(n_truth+1);
         truth_values[0] = 0;
         cls_values[0] = 1;
         cls_uncertainties_n[0] = 0;
         cls_uncertainties_n0[0] = 0;
         map<double, multiset<double> >::const_iterator it_b = truth_to_ts_b.begin();
         map<double, multiset<double> >::const_iterator it_sb = truth_to_ts_sb.begin();
         for(size_t i=1; it_b!=truth_to_ts_b.end(); ++it_b, ++it_sb, ++i){
            truth_values[i] = it_b->first;
            if(it_b->first != it_sb->first){
                throw InvalidArgumentException("truth_ts_values::get_cls_vs_truth: truth_to_ts_sb and truth_to_b contain different truth values");
            }
            cls_info res = get_cls(ts_value, it_b->first);
            cls_values[i] = res.cls;
            cls_uncertainties_n0[i] = res.cls_uncertainty_n0;
            cls_uncertainties_n[i] = res.cls_uncertainty_n;
         }
         return cls_vs_truth_data(truth_values, cls_values, cls_uncertainties_n0, cls_uncertainties_n);
    }
};


template<typename T>
const std::auto_ptr<std::ostream> & operator<<(const std::auto_ptr<std::ostream> & p_out, const T & t){
    if(p_out.get()){
        (*p_out) << t;
    }
    return p_out;
}

void flush(const std::auto_ptr<std::ostream> & p_out){
    if(p_out.get()){
        (*p_out) << std::flush;
    }
}



struct fitexp_result{
    double limit, limit_error;
};

struct fitexp_parameters{
    Minimizer & minimizer;
    ParId pid_limit, pid_lambda;
    fitexp_parameters(Minimizer & min_, const ParId & limit_, const ParId & lambda_): minimizer(min_), pid_limit(limit_), pid_lambda(lambda_){}
};

// fit from xmin to (including) xmax. xmin and xmax must be truth_values in data.
fitexp_result fitexp(const cls_vs_truth_data & data, double target_cls, fitexp_parameters & pars, double xmin, double xmax, std::auto_ptr<std::ostream> & debug_out){
    theta_assert(xmax > xmin);
    const vector<double> & truth_values = data.truth_values();
    const vector<double> & cls_values = data.cls_values();
    const vector<double> & cls_uncertainties_n0 = data.cls_uncertainties_n0();
    const vector<double> & cls_uncertainties_n = data.cls_uncertainties_n();
    
    size_t imin = find(truth_values.begin(), truth_values.end(), xmin) - truth_values.begin();
    size_t imax = find(truth_values.begin(), truth_values.end(), xmax) - truth_values.begin();
    if(imin==truth_values.size() || imax==truth_values.size()){
        throw InvalidArgumentException("fitexp: xmin, xmax not found in data");
    }
    if(imin >= imax){
        throw InvalidArgumentException("fitexp: xmin, xmax define empty range");
    }
    const size_t N_range = imax - imin + 1;
    vector<double> x_values(truth_values.begin() + imin, truth_values.begin() + imax + 1);
    vector<double> y_values(cls_values.begin() + imin, cls_values.begin() + imax + 1);
    debug_out << "fitexp: fitting to imin, imax = (" << imin << ", " << imax << ");  xmin, xmax = "  << x_values[0] << ", " << x_values[x_values.size()-1] << ")\n";
    Matrix cov_inv(N_range, N_range);
    // build inverse covariance matrix:
    for(size_t i=0; i<N_range; ++i){
        cov_inv(i, i) = 1.0 / (pow(cls_uncertainties_n0[imin+i], 2) + pow(cls_uncertainties_n[imin+i], 2));
    }
    // use first and last point in range to determine initial lambda value; take care not to use log(0), so always add the
    // uncertainty which is always > 0.
    double lambda0 = min((log(y_values[N_range-1] + sqrt(1.0 / cov_inv(N_range-1, N_range-1))) - log(y_values[0] + sqrt(1 / cov_inv(0,0)))) / (x_values[N_range-1] - x_values[0]), -1e-10);
    double limit0 = 0.5 * (x_values[0] + x_values[N_range-1]);
    debug_out << "fitexp: limit0, lambda0 = " << limit0 << ", " << lambda0 << "\n";
    exp_nll nll(pars.pid_lambda, pars.pid_limit, target_cls, x_values, y_values, cov_inv);
    ParValues start, step;
    start.set(pars.pid_lambda, lambda0);
    start.set(pars.pid_limit, limit0);
    step.set(pars.pid_lambda, fabs(lambda0) * 0.1);
    step.set(pars.pid_limit, fabs(limit0) * 0.1);
    const double inf = numeric_limits<double>::infinity();
    map<ParId, pair<double, double> > ranges;
    ranges[pars.pid_lambda] = make_pair(-inf, 0.0);
    ranges[pars.pid_limit] = make_pair(0.0, inf);
    fitexp_result result;
    result.limit = result.limit_error = NAN;
    try{
        MinimizationResult res = pars.minimizer.minimize(nll, start, step, ranges);
        result.limit = res.values.get(pars.pid_limit);
        result.limit_error = 0.5 * (res.errors_plus.get(pars.pid_limit) + res.errors_minus.get(pars.pid_limit));
    }
    catch(MinimizationException & ex){
        debug_out << "fitexp_range: MinimizationException: " << ex.message << "\n";
        //result will be NAN
    }
    return result;
}



/** \brief Calculate CLs limits
 *
 *
 * \code
 * main = {
 *   type = "cls_limits";
 *
 *   // same as in type = default (theta::Run):
 *   producers = ("@lr");
 *   model = "@model";
 *   output_database = {
 *      type = "sqlite_database";
 *      filename = "toys.db";
 *   };
 *
 *   // specific to cls_toys
 *   minimizer = { ... };
 *   ts_column = "lr__nll_diff"; // full column name, with optional '-'
 *   truth_parameter = "beta_signal";
 *
 *   // optional parameters
 *   cl = 0.90; //optional; default is 0.95
 *   ts_values = (1.0, 1.7); //optional
 *   ts_values_background_bands = true; //optional
 *   reltol_limit = 0.001; // optional; default is 0.05
 *   tol_cls = 0.005; //optional; default is 0.015
 *   limit_hint = (200.0, 240.0); //optional; default is finding out on its own ...
 *   reuse_toys = (
 *       {type = "sqlite_database_in"; filename = "toys0.db"; }
 *   ); //optional; default is the empty list (), i.e., not to reuse any toys
 *   debuglog = "debug.txt"; // optional, default is no debug output
 * };
 * \endcode
 * 
 * Toys are drawn from the \c model with different values for the signal parameter given in \c signal_parameter (all other
 * parameters are drawn randomly from the model's parameter distribution). For each toy, the \c ts_producer is run and the
 * value from the given \c ts_column uis used as test statistic for the construction of CLs limits.
 * 
 * Which values for the signal parameter are used and how many toys are draw for each parameter value is determined automatically and based
 * on \c ts_values, \c reltol_limit, and \c tol_cls: \c ts_values controls for which test statistic values the toys calculation should be done such that the CLs limit
 * has a relative 1sigma uncertainty of at most \c reltol_limit and the (absolute) tolerance for the CLs value is \c tol_cls.
 *
 * Valid values for \c ts_values are either (i) a list of values or (ii) a DataSource specification
 * in which case the ts producer is run once for that DataSource. If \c ts_values_background_bands is true,
 * the background band is calculated in addition, i.e., the test statistic values at the median, central 1sigma and 2sigma
 * for toys with signal_parameter=0.
 * The neither \c ts_values and \c ts_values_background_bands are given, only background bands
 * are calculated, i.e., \c ts_values_background_bands defaults to \c true in this case.
 * If \c ts_values is given, \c ts_values_background_bands defaults to \c false and has to be set to \c true explicitly.
 *
 * If \c limit_hint is given, it should be an interval where the limits are probably contained,
 * it does not have to be exact; teh interval does not have to contain the true value. If given, the convergence can be faster. 
 *
 * Both, the toys used for the construction as well as the resulting limits are written to the
 * \c output_database. For the toys, the tables are the same as for theta::Run. In addition, the table
 * "cls_limits" is created which contains the columns "index", "ts_value", "limit", and "limit_uncertainty".
 * If sorted by index, this table will contain the limits for the test statistic values specified via
 * \c ts_values first, and then (if applicable) the 5 limits defining the "expected" band, in increasing
 * order of the ts value (=lower edge of 2sigma band, lower edge of 1sigma band, median, upper 1sigma edge, upper 2sigma edge).
 * 
 * (TODO: this is not implemented yet!) The optional setting \c reuse_toys is a list of theta::DatabaseInput
 * specifications. If present, it is assumed that these contain toys
 * from previous runs with the same model and ts definition (in particular, with the same truth_parameter and ts_column name ...).
 * Reusing toys can speed up the calculation substiantially.
 *
 * \c debuglog is the filename of the debug log. Debug output is only written if a filename is given here.
 *
 * The indicated progress and errors refer to the number of toys produced. It is hard to tell how many toys are
 * necessary; for usual settings in the order of 5000-10000.
 */
class cls_limits: public Main{
public:
    cls_limits(const Configuration & cfg);
    virtual void run();
    
private:
    class data_filler: public theta::RandomConsumer, public theta::ProductsSource{
    private:
        boost::shared_ptr<Model> model;
        Column truth_column;
    public:
        data_filler(const Configuration & cfg, const boost::shared_ptr<Model> & model_): RandomConsumer(cfg, "source"), ProductsSource("source", cfg.pm->get<ProductsSink>()), model(model_){
            truth_column = products_sink->declare_product(*this, "truth", theta::typeDouble);
        }
        void fill(Data & dat, const ParId & truth_parameter, double truth_value){
            dat.reset();
            Random & rnd = *rnd_gen;
            ParValues values;
            model->get_parameter_distribution().sample(values, rnd);
            values.set(truth_parameter, truth_value);
            products_sink->set_product(truth_column, truth_value);
            model->get_prediction(dat, values);
            ObsIds observables = dat.getObservables();
            for (ObsIds::const_iterator it = observables.begin(); it != observables.end(); it++) {
                 randomize_poisson(dat[*it], rnd);
            }
        }
    };
    
    void run_single_truth(double truth, bool bkg_only, int n_event);
    
    // run toys at the given truth value (and at truth=0!) until either
    // * the uncertainty on the CLs value is below tol_cls
    //   or
    // * the CLs value is found to be more than 3 sigma away from the target CLs 1-cl.
    void run_single_truth_adaptive(double ts_value, double truth);

    boost::shared_ptr<VarIdManager> vm;
    boost::shared_ptr<Model> model;
    //note: the tables hold a shared_ptr to this database to ensure proper destruction order
    boost::shared_ptr<Database> db;

    std::auto_ptr<LogTable> logtable;
    boost::shared_ptr<RndInfoTable> rndinfo_table;

    //the producers to be run on the pseudo data:
    boost::ptr_vector<Producer> producers;
    boost::shared_ptr<MultiplexingProductsSink> mps;
    boost::shared_ptr<SaveDoubleColumn> sdc;
    double ts_factor; // either +1.0 or -1.0
    boost::shared_ptr<ProductsTable> products_table;
    
    std::auto_ptr<Table> cls_limits_table;
    Column cls_limits__index, cls_limits__ts_value, cls_limits__limit, cls_limits__limit_uncertainty;
    
    ParId truth_parameter;
    std::auto_ptr<data_filler> source;
    
    // for fitting the CLs versus truth curve, we need a minimizer and some parameters:
    std::auto_ptr<Minimizer> minimizer;
    ParId pid_limit, pid_lambda;
    
    truth_ts_values tts;

    int runid;
    
    bool ts_background_bands;
    std::auto_ptr<theta::DataSource> ts_data_source;
    vector<double> ts_values;
    
    pair<double, double> limit_hint;
    double reltol_limit, tol_cls;
    double cl;
    
    std::auto_ptr<std::ostream> debug_out;
    
    int n_toys, n_toy_errors;
};


void cls_limits::run_single_truth(double truth, bool bkg_only, int n_event){
    if(!isfinite(truth)) throw InvalidArgumentException("run_single_truth: truth not finite");
    Data data;
    vector<double> ts_values;
    ts_values.reserve(n_event);
    for (int eventid = 1; eventid <= n_event; eventid++) {
        if(stop_execution)break;
        source->fill(data, truth_parameter, bkg_only?0.0:truth);
        bool error = false;
        for (size_t j = 0; j < producers.size(); j++) {
            try {
                producers[j].produce(data, *model);
            } catch (Exception & ex) {
                error = true;
                std::stringstream ss;
                ss << "Producer '" << producers[j].getName() << "' failed: " << ex.message << ".";
                logtable->append(runid, eventid, LogTable::error, ss.str());
                ++n_toy_errors;
                break;
            }
            catch(FatalException & f){
                stringstream ss;
                ss << "Producer '" << producers[j].getName() << "': " << f.message;
                f.message = ss.str();
                throw;
            }
        }
        ++n_toys;
        if(!error){
            ts_values.push_back(ts_factor * sdc->get_value());
            products_table->add_row(runid, eventid);
        }
        if(progress_listener){
            progress_listener->progress(n_toys, -1, n_toy_errors);
        }
    }
    if(ts_values.size() * 1.0 / n_event < 0.8){
        throw FatalException("cls_limits: ts_producers fails in more than 20% of the cases");
    }
    if(bkg_only){
       tts.add_points_b(truth, ts_values.begin(), ts_values.end());
    }
    else{
       tts.add_points_sb(truth, ts_values.begin(), ts_values.end());
    }
}


void cls_limits::run_single_truth_adaptive(double ts_value, double truth){
    if(!tts.contains_truth_value(truth)){
        run_single_truth(truth, true, 100);
        run_single_truth(truth, false, 100);
    }
    for(int kk=0; kk<100; ++kk){
        truth_ts_values::cls_info res = tts.get_cls(ts_value, truth);
        double cls_uncertainty = sqrt(pow(res.cls_uncertainty_n0, 2) + pow(res.cls_uncertainty_n, 2));
        debug_out << "run_single_truth_adaptive (iteration " << kk << "):\n";
        debug_out << " clsb = " << res.clsb << " +- " << res.clsb_uncertainty << "\n";
        debug_out << " clb  = " << res.clb  << " +- " << res.clb_uncertainty << "\n";
        debug_out << " cls  = " << res.cls  << " +- " << cls_uncertainty << "\n";
        if(fabs(res.cls - (1 - cl)) / cls_uncertainty > 3) return;
        if(cls_uncertainty < tol_cls) return;
        double target_cls_uncertainty = max(tol_cls, fabs(res.cls - (1 - cl))/3) * 0.7;
        if(res.cls_uncertainty_n0 > target_cls_uncertainty){
            size_t n0 = tts.get_n0(truth);
            size_t n0_new = pow(res.cls_uncertainty_n0 / target_cls_uncertainty, 2) * n0;
            theta_assert(n0_new >= n0);
            size_t n_toys = max(50ul, min<size_t>(n0_new - n0, 500));
            debug_out << "run_single_truth_adaptive: making " << n_toys << " more toys for background only\n";
            flush(debug_out);
            run_single_truth(truth, true, n_toys);
        }
        if(res.cls_uncertainty_n > target_cls_uncertainty){
            size_t n = tts.get_n(truth);
            size_t n_new = pow(res.cls_uncertainty_n / target_cls_uncertainty, 2) * n;
            theta_assert(n_new >= n);
            size_t n_toys = max(50ul, min<size_t>(n_new - n, 500));
            debug_out << "run_single_truth_adaptive: making " << n_toys << " more toys for signal + background\n";
            flush(debug_out);
            run_single_truth(truth, false, n_toys);
        }
    }
    throw FatalException("run_single_truth_adaptive: could not reach target accuracy for CLs after 100 iterations");
}



cls_limits::cls_limits(const Configuration & cfg): vm(cfg.pm->get<VarIdManager>()), ts_factor(1.0), truth_parameter(vm->getParId(cfg.setting["truth_parameter"])),
  pid_limit(vm->createParId("__limit")), pid_lambda(vm->createParId("__lambda")), runid(1),
ts_background_bands(true), limit_hint(NAN, NAN), reltol_limit(0.05), tol_cls(0.015), cl(0.95), n_toys(0), n_toy_errors(0) {
    SettingWrapper s = cfg.setting;
    
    //1. setup database and tables:
    db = PluginManager<Database>::instance().build(Configuration(cfg, s["output_database"]));

    std::auto_ptr<Table> logtable_underlying = db->create_table("log");
    logtable.reset(new LogTable(logtable_underlying));
    
    std::auto_ptr<Table> rndinfo_table_underlying = db->create_table("rndinfo");
    rndinfo_table.reset(new RndInfoTable(rndinfo_table_underlying));
    cfg.pm->set("default", rndinfo_table);
    
    std::auto_ptr<Table> products_table_underlying = db->create_table("products");
    products_table.reset(new ProductsTable(products_table_underlying));
    string colname = s["ts_column"];
    while(colname.size() > 0 && (colname[0]==' ' or colname[0]=='-')){
        if(colname[0]=='-') ts_factor = -1.0;
        colname = colname.substr(1);
    }
    sdc.reset(new SaveDoubleColumn(colname));
    std::vector<boost::shared_ptr<ProductsSink> > sinks;
    sinks.push_back(products_table);
    sinks.push_back(sdc);
    mps.reset(new MultiplexingProductsSink(sinks));
    cfg.pm->set<ProductsSink>("default", mps);
    
    boost::shared_ptr<int> ptr_runid(new int(runid));
    cfg.pm->set("runid", ptr_runid);
    
    cls_limits_table = db->create_table("cls_limits");
    cls_limits__index = cls_limits_table->add_column("index", typeInt);
    cls_limits__ts_value = cls_limits_table->add_column("ts_value", typeDouble);
    cls_limits__limit = cls_limits_table->add_column("limit", typeDouble);
    cls_limits__limit_uncertainty = cls_limits_table->add_column("limit_uncertainty", typeDouble);
            
    //2. model, data_source and minimizer
    model = PluginManager<Model>::instance().build(Configuration(cfg, s["model"]));
    minimizer = PluginManager<Minimizer>::instance().build(Configuration(cfg, s["minimizer"]));
    source.reset(new data_filler(cfg, model));
    
    //3. logging stuff
    logtable->set_loglevel(LogTable::warning);
    
    //4. producers:
    size_t n_p = s["producers"].size();
    if (n_p == 0)
        throw ConfigurationException("list of producers is empty");
    for (size_t i = 0; i < n_p; i++) {
         producers.push_back(PluginManager<Producer>::instance().build(Configuration(cfg, s["producers"][i])));
    }
    
    //5. ts_values
    if(s.exists("ts_values")){
        ts_background_bands = false;
        libconfig::Setting::Type type = s["ts_values"].getType();
        if(type==libconfig::Setting::TypeArray or type==libconfig::Setting::TypeList){
            size_t n_ts = s["ts_values"].size();
            if (n_ts == 0){
                throw ConfigurationException("list of ts_values is empty");
            }
            for (size_t i = 0; i < n_ts; i++) {
                ts_values.push_back(s["ts_values"][i]);
            }
        }
        else if(type==libconfig::Setting::TypeGroup){
            ts_data_source = PluginManager<DataSource>::instance().build(Configuration(cfg, s["ts_values"]));
        }
        else{
            throw ConfigurationException("invalid type for 'ts_values': must be either a list/array or a setting group");
        }
    }
    if(s.exists("ts_values_background_bands")){
        ts_background_bands = s["ts_values_background_bands"];
    }
    
    //6. misc
    if(s.exists("reltol_limit")) reltol_limit = s["reltol_limit"];
    if(s.exists("limit_hint")){
        limit_hint.first = s["limit_hint"][0];
        limit_hint.second = s["limit_hint"][1];
    }
    if(s.exists("cl")){
        cl = s["cl"];
        if(cl >= 0.999 || cl <= 0) throw ConfigurationException("invalid value for cl. Valid range is (0, 0.999)");
    }
    if(s.exists("debuglog")){
        string fname = s["debuglog"];
        debug_out.reset(new ofstream(fname.c_str()));
    }
}


bool cls_is_significantly_larger(const cls_vs_truth_data & data, size_t i, double target_cls){
    if(i==0) return true;
    const vector<double> & cls_values = data.cls_values();
    const vector<double> & cls_uncertainties_n = data.cls_uncertainties_n();
    const vector<double> & cls_uncertainties_n0 = data.cls_uncertainties_n0();
    return cls_values[i] > target_cls && fabs(cls_values[i] - target_cls) / sqrt(pow(cls_uncertainties_n[i], 2) + pow(cls_uncertainties_n0[i], 2)) > 3;
}

bool cls_is_significantly_smaller(const cls_vs_truth_data & data, size_t i, double target_cls){
    const vector<double> & cls_values = data.cls_values();
    const vector<double> & cls_uncertainties_n = data.cls_uncertainties_n();
    const vector<double> & cls_uncertainties_n0 = data.cls_uncertainties_n0();
    return cls_values[i] < target_cls && fabs(cls_values[i] - target_cls) / sqrt(pow(cls_uncertainties_n[i], 2) + pow(cls_uncertainties_n0[i], 2)) > 3;
}

void cls_limits::run(){
    //0. determine ts_epsilon. This is important for degenerate test statistic (such as likelihood ratios where
    // the backgruond only prior fixes the signal and is free in the s+b prior), because mathematically degenerate
    // values will not be necessarily exactly degenerate numerically. A possible solution to that is
    // to add small values to all 'observed' test statistic values, i.e., make them more signal-like (although
    // only by a very small amount).
    run_single_truth(0.0, false, 200);
    const double ts_epsilon = fabs(tts.get_ts_quantile(0.0, 0.975) - tts.get_ts_quantile(0.0, 0.025)) * 1e-3;
    debug_out << "ts_epsilon is " << ts_epsilon << "\n";
    //1.a. determine the ts value for the data source, if set:
    if(ts_data_source.get()){
        Data data;
        ts_data_source->fill(data);
        for (size_t j = 0; j < producers.size(); j++) {
            try {
                producers[j].produce(data, *model);
            } catch (Exception & ex) {
                ex.message += " (while running ts_producer for data)";
                throw;
            }
        }
        ts_values.push_back(ts_factor * sdc->get_value() + ts_epsilon);
    }
    //1.b. determine the ts_background_bands:
    if(ts_background_bands){
        debug_out << "Requested expected limit. For this, running some toys at truth=0:\n";
        flush(debug_out);
        run_single_truth(0.0, false, 1000);
        run_single_truth(0.0, true, 50);
        ts_values.push_back(tts.get_ts_quantile(0.0, 0.025) + ts_epsilon);
        ts_values.push_back(tts.get_ts_quantile(0.0, 0.16) + ts_epsilon);
        ts_values.push_back(tts.get_ts_quantile(0.0, 0.5) + ts_epsilon);
        ts_values.push_back(tts.get_ts_quantile(0.0, 0.84) + ts_epsilon);
        ts_values.push_back(tts.get_ts_quantile(0.0, 0.975) + ts_epsilon);
    }
    //2. some other values in about the right region ...
    double signal_width = limit_hint.second - limit_hint.first;
    if(!isfinite(limit_hint.first) || !isfinite(limit_hint.second)){
        ParValues widths = asimov_likelihood_widths(*model);
        signal_width = widths.get(truth_parameter);
        debug_out << "signal_width = " << signal_width << "\n";
        if(signal_width <= 0.0){
            debug_out << "Warning: signal_width <=0.0. Setting signal_width to 1.0.\n";
            signal_width = 1.0;
        }
        flush(debug_out);
    }
    debug_out << "ts_values.size() == " << ts_values.size() << "\n";
    //3. calculate the CLs limit for all ts_values:
    fitexp_parameters pars(*minimizer, pid_limit, pid_lambda);
    const size_t N_maxit = 200;
    for(size_t k=0; k<ts_values.size(); ++k){
        debug_out << "starting ts value index " << k << ", " << ts_values[k] << "\n";
        flush(debug_out);
        //3.a. find a "seed", i.e., an interval of truth values significantly larger than target_cls and significantly smaller than target_cls:
        //3.a.i. make toys to get a truth value with significantly smaller cls value than target_cls:
        cls_vs_truth_data data = tts.get_cls_vs_truth(ts_values[k]);
        while(not cls_is_significantly_smaller(data, data.cls_values().size()-1, 1 - cl)){
             double next_value;
             if(isfinite(limit_hint.first) && isfinite(limit_hint.second) && data.cls_values().back() < limit_hint.second){
                 next_value = limit_hint.second;
             }
             else if(data.truth_values().back() > 0.0){
                 next_value = data.truth_values().back()*2.0;
             }
             else{
                 next_value = signal_width;
             }
             debug_out << "making toys at high truth values to find upper limit on limit; next is truth=" << next_value << "\n";
             flush(debug_out);
             run_single_truth_adaptive(ts_values[k], next_value);
             data = tts.get_cls_vs_truth(ts_values[k]);
        }
        //3.a.ii. find a seed in truth values (i_low to i_high, both inclusive):
        size_t i_low, i_high;
        i_high = 0;
        i_low = data.cls_values().size()-1;
        while(not cls_is_significantly_smaller(data, i_high, 1-cl)) ++i_high; //terminates because we created a significantly smaller point in 3.a.i.
        while(not cls_is_significantly_larger(data, i_low, 1-cl)) --i_low; //terminates because this is true for i_low==0 (truth==0.0)
        double truth_low = data.truth_values()[i_low];
        double truth_high = data.truth_values()[i_high];
        debug_out << k << " seed: interval [" << truth_low << ", " << truth_high << "]; i_low, ihigh = (" << i_low << ", " << i_high << ")\n";
        flush(debug_out);
        //3.b. iterate
        fitexp_result latest_res;
        for(size_t i=0; i<=N_maxit; ++i){
            debug_out << k << "." << i << "\n";
            if(debug_out.get()){
                //debug_out << "tts:\n";
                //tts.write_txt(*debug_out, 70, -2, 5);
                debug_out << "cls vs truth:\n";
                data.write_txt(*debug_out);
                flush(debug_out);
            }
            if(i==N_maxit){
                throw FatalException("CLs method did not converge: too many iterations necessary. See debuglog for details.");
            }
            latest_res = fitexp(data, 1 - cl, pars, truth_low, truth_high, debug_out);
            theta_assert(latest_res.limit >= truth_low && latest_res.limit <= truth_high);
            debug_out << k << "." << i << " result: limit = " << latest_res.limit << " +- " << latest_res.limit_error << "\n";
            if(latest_res.limit_error / latest_res.limit < reltol_limit){
                debug_out << k << "." << i << ": result is good enough; done with this ts value.\n";
                break;
            }
            else{
                debug_out << k << "." << i << ": result is NOT good enough yet: relative limit accuracy is "
                          << latest_res.limit_error / latest_res.limit << ", target accuracy is " << reltol_limit << ".\n";
            }
            flush(debug_out);
            // add a point at the estimated limit:
            run_single_truth_adaptive(ts_values[k], latest_res.limit);
            data = tts.get_cls_vs_truth(ts_values[k]);
            // make interval [truth_low, truth_high] smaller if border is far away (>2sigma) from limit AND the next point can serve as new interval border
            i_low = find(data.truth_values().begin(), data.truth_values().end(), truth_low) - data.truth_values().begin();
            i_high = find(data.truth_values().begin(), data.truth_values().end(), truth_high) - data.truth_values().begin();
            if(fabs(truth_low - latest_res.limit) / latest_res.limit_error > 2 && i_low+1 < i_high){
                 debug_out << "lower border " << truth_low << " far away from current limit --> test if can be removed ...\n";
                 //re-run next candidate point:
                 run_single_truth_adaptive(ts_values[k], data.truth_values()[i_low+1]);
                 if(cls_is_significantly_larger(data, i_low+1, 1 - cl)){
                     ++i_low;
                     truth_low = data.truth_values()[i_low];
                     debug_out << "--> yes, replaced by " << truth_low << ".\n";
                 }
                 else{
                     debug_out << "--> no, next point not valid.\n";
                 }
            }
            if(fabs(truth_high - latest_res.limit) / latest_res.limit_error > 2 && i_low < i_high-1){
                 debug_out << "upper border " << truth_high << " far away from current limit --> test if can be removed ...\n";
                 //re-run next candidate point:
                 run_single_truth_adaptive(ts_values[k], data.truth_values()[i_high-1]);
                 if(cls_is_significantly_smaller(data, i_high-1, 1 - cl)){
                     --i_high;
                     truth_high = data.truth_values()[i_high];
                     debug_out << "--> yes, replaced by " << truth_high << ".\n";
                 }
                 else{
                     debug_out << "--> no, next point not valid.\n";
                 }
            }
        }
        Row r;
        r.set_column(cls_limits__index, (int)k);
        r.set_column(cls_limits__ts_value, ts_values[k]);
        r.set_column(cls_limits__limit, latest_res.limit);
        r.set_column(cls_limits__limit_uncertainty, latest_res.limit_error);
        cls_limits_table->add_row(r);
    }
}

REGISTER_PLUGIN(cls_limits)

