#include "interface/cfg-utils.hpp"
#include "interface/plugin.hpp"
#include "interface/histogram.hpp"
#include "interface/variables-utils.hpp"
#include "interface/minimizer.hpp"
#include "interface/random.hpp"
#include "interface/database.hpp"
#include "interface/main.hpp"
#include "interface/redirect_stdio.hpp"
#include "interface/phys.hpp"
#include "interface/model.hpp"
#include "interface/distribution.hpp"
#include "interface/random-utils.hpp"
#include "plugins/asimov_likelihood_widths.hpp"

#include <fstream>
#include <iomanip>
#include <boost/optional.hpp>

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
    
    // absolute uncertainty on the CLs values due to finite number of s+b toys ("n")
    const vector<double> & cls_uncertainties_n() const{
        return pcls_uncertainties_n;
    }
    
    // absolute uncertainty on the CLs values due to finite number of b only toys ("n0")
    const vector<double> & cls_uncertainties_n0() const{
        return pcls_uncertainties_n0;
    }
    void write_txt(ostream & out) const{
        out << "# truth cls  cls_error_n0  cls_error_n" << endl;
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
    // bulk insertion (I guess this has much better performance ...)
    template <class InputIterator>
    void add_points_b(double truth, InputIterator first, InputIterator last){
        truth_to_ts_b[truth].insert(first, last);
    }
    template <class InputIterator>
    void add_points_sb(double truth, InputIterator first, InputIterator last){
        truth_to_ts_sb[truth].insert(first, last);
    }

    // get the test statistic quantile q (0 <= q < 1) for the given truth value for the b-only case.
    double get_ts_b_quantile(double truth, double q){
        theta_assert(q>=0.0 && q < 1.0);
        map<double, multiset<double> >::const_iterator it = truth_to_ts_b.find(truth);
        if(it==truth_to_ts_b.end()) throw NotFoundException("truth_ts_value::get_ts_b_quantile");
        const multiset<double> & ts_values = it->second;
        multiset<double>::const_iterator itt = ts_values.begin();
        advance(itt, q*ts_values.size());
        return *itt;
    }
    
    bool contains_truth_value(double truth) const{
        return truth_to_ts_sb.find(truth) != truth_to_ts_sb.end() && truth_to_ts_b.find(truth) != truth_to_ts_b.end();
    }

    // get all truth values we have toys for.
    set<double> truth_values() const{
        set<double> result;
        result.insert(0.0);
        for(map<double, multiset<double> >::const_iterator it=truth_to_ts_sb.begin(); it!=truth_to_ts_sb.end(); ++it){
            result.insert(it->first);
        }
        return result;
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
            throw InvalidArgumentException("truth_ts_values::get_cls: truth_to_ts_sb or truth_to_b do not contain given truth value");
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
    
    // truth_to_ts is a map which specifies -- for each truth value -- the test statistic value.
    // (This is important for cases in which the test statistic value depends on the truth parameter).
    // The caller has to make sure that the ts value in this map is set for all truth values, see truth_values, including
    // truth = 0.0 (!).
    //
    // For the function "cls versus truth", there are two error sources for the error on cls:
    // * the finite number of toys for truth=0  (n0), correlated between different truth values
    // * the finite number of toys for the current truth value (n) which is uncorrelated
    // the two covariance matrices (where the latter is diagonal ...) are calculated and stored separately. This is helpful for the decision on
    // where to perform more toys, at truth=0, or at the intersection ...
    cls_vs_truth_data get_cls_vs_truth(const map<double, double> & truth_to_ts) const{
         theta_assert(truth_to_ts_b.size() == truth_to_ts_sb.size());
         const size_t n_truth = truth_to_ts_sb.size();
         vector<double> truth_values(n_truth);
         vector<double> cls_values(n_truth);
         vector<double> cls_uncertainties_n(n_truth);
         vector<double> cls_uncertainties_n0(n_truth);
         truth_values[0] = 0;
         cls_values[0] = 1;
         cls_uncertainties_n[0] = 0.001;
         cls_uncertainties_n0[0] = 0.001;
         map<double, multiset<double> >::const_iterator it_b = truth_to_ts_b.begin();
         map<double, multiset<double> >::const_iterator it_sb = truth_to_ts_sb.begin();
         // skip truth=0:
         ++it_b;
         ++it_sb;
         for(size_t i=1; it_b!=truth_to_ts_b.end(); ++it_b, ++it_sb, ++i){
            truth_values[i] = it_b->first;
            if(it_b->first != it_sb->first){
                throw InvalidArgumentException("truth_ts_values::get_cls_vs_truth: truth_to_ts_sb and truth_to_b contain different truth values");
            }
            map<double, double>::const_iterator it_t_ts = truth_to_ts.find(truth_values[i]);
            theta_assert(it_t_ts != truth_to_ts.end());
            double ts = it_t_ts->second;
            cls_info res = get_cls(ts, it_b->first);
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
    step.set(pars.pid_limit, fabs(x_values.back() - x_values[0]) * 0.1);
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
 *   // (almost) same as in type = "default" (theta::Run):
 *   producer = "@lr";
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
 *   cl = 0.90; // optional; default is 0.95
 *   data_source = { some datasource specification }; //optional
 *   expected_bands = true; //optional
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
 * \c data_source and \c expected_bands control for which data / test statistic the limits are calculated. \c data_source
 * is a data source and used to caluclate the observed limit on actual data. If \c expected_bands is true, the central expected,
 * +-1sigma and +-2sigma bands of expected limits are calculated as well. \c expected_bands defaults to true if data_source is
 * given and to false otherwise.
 *
 * The calculation of CLs limits in general requires making a large number of toys at different truth values, to scan for the 
 * desired CLs value 1 - cl. How many toys are drawn and where is done automatically and is driven by \c reltol_limit:
 * this is the relative 1sigma uncertainty tolerance for the limit. \c tol_cls is an absolute accuracy for the CLs value. It is used
 * internally only and does not affect the accuracy of the result directly. However, it can affect robustness of the method. Set it
 * to a smaller value if you 
 *
 * If \c limit_hint is given, it should be an interval where the limits are probably contained.
 * It does not have to be exact, but wil be used for a starting point to make some toys which can make the
 * convergence considerably faster.
 *
 * Both, the toys used for the construction as well as the resulting limits are written to the
 * \c output_database. For the toys, the tables are the same as for theta::Run. In addition, the table
 * "cls_limits" is created which contains the columns "q", "limit", and "limit_uncertainty".
 * The column q is the quantile of the background-only distribution used as test statistic. It is
 * The values 0.025, 0.16, 0.5, 0.84, and 0.975 define the "expected" bands. the special value -1 is used
 * to indicate that data was used from the data_source.
 * 
 * (TODO: this is not implemented yet!) The optional setting \c reuse_toys is a list of theta::DatabaseInput
 * specifications. If present, it is assumed that these contain toys
 * from previous runs with the same model and ts definition (in particular, with the same truth_parameter and ts_column name ...).
 * Reusing toys can speed up the calculation substiantially.
 *
 * \c debuglog is the filename of the debug log. Debug output is only written if a filename is given here.
 *
 * The indicated progress and errors refer to the number of toys produced. It is hard to tell how many toys are
 * necessary; this also depends on the CLb value: a low CLb value requires much more toys to rech the desired accuracy.
 * (note that if calculating the expected limit band, low CLb values exist by construction).
 * In such a case (CLb around 0.03), about 20k to 50k toys are necessary. For a more typical case (CLb around 0.5),
 * usually 5-10k of toys are sufficient. This is also true in cases where the test statistic is degenerate (= has often the same value)
 * as in such a case, low CLb values do not occur.
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
    //void run_single_truth_adaptive(double ts_value, double truth);
    void run_single_truth_adaptive(double q, map<double, double> & truth_to_ts, double ts_epsilon, double truth);
    void update_truth_to_ts(map<double, double> & truth_to_ts, double q, double ts_epsilon);

    boost::shared_ptr<VarIdManager> vm;
    boost::shared_ptr<Model> model;
    //note: the tables hold a shared_ptr to this database to ensure proper destruction order
    boost::shared_ptr<Database> db;

    std::auto_ptr<LogTable> logtable;
    boost::shared_ptr<RndInfoTable> rndinfo_table;

    //the producer to be run on the pseudo data which provides the test statistic:
    std::auto_ptr<Producer> producer;
    ParameterDependentProducer * pp_producer; // points to producer, if type matched. 0 otherwise.
    boost::shared_ptr<MultiplexingProductsSink> mps;
    boost::shared_ptr<SaveDoubleColumn> sdc;
    boost::shared_ptr<ProductsTable> products_table;
    
    std::auto_ptr<Table> cls_limits_table;
    Column cls_limits__ts_value, cls_limits__limit, cls_limits__limit_uncertainty;
    
    ParId truth_parameter;
    std::auto_ptr<data_filler> source;
    
    // for fitting the CLs versus truth curve, we need a minimizer and some parameters:
    std::auto_ptr<Minimizer> minimizer;
    ParId pid_limit, pid_lambda;
    
    truth_ts_values tts;

    int runid;
    
    bool expected_bands;
    std::auto_ptr<theta::DataSource> data_source;
    
    pair<double, double> limit_hint;
    double reltol_limit, tol_cls;
    double cl;
    
    std::auto_ptr<std::ostream> debug_out;
    
    int n_toys, n_toy_errors;
};


void cls_limits::run_single_truth(double truth, bool bkg_only, int n_event){
    if(!isfinite(truth)) throw InvalidArgumentException("run_single_truth: truth not finite");
    //debug_out << "starting run_single_truth(truth = " << truth << ", bkg_only = " << bkg_only << ", n_event = " << n_event << ")\n";
    // if the producer is parameter dependent, pass the truth value to the producer:
    if(pp_producer){
        ParValues values;
        values.set(truth_parameter, truth);
        pp_producer->setParameterValues(values);
    }
    Data data;
    vector<double> ts_values;
    ts_values.reserve(n_event);
    for (int eventid = 1; eventid <= n_event; eventid++) {
        if(stop_execution) break;
        source->fill(data, truth_parameter, bkg_only?0.0:truth);
        bool error = false;
        try {
            producer->produce(data, *model);
        } catch (Exception & ex) {
            error = true;
            std::stringstream ss;
            ss << "Producer '" << producer->getName() << "' failed: " << ex.message << ".";
            logtable->append(runid, eventid, LogTable::error, ss.str());
            ++n_toy_errors;
            continue;
        }
	    catch(FatalException & f){
		    stringstream ss;
		    ss << "Producer '" << producer->getName() << "': " << f.message;
		    f.message = ss.str();
		    throw;
	    }
        ++n_toys;
        if(!error){
            ts_values.push_back(sdc->get_value());
            products_table->add_row(runid, eventid);
        }
        if(progress_listener){
            progress_listener->progress(n_toys, -1, n_toy_errors);
        }
    }
    if(bkg_only){
       tts.add_points_b(truth, ts_values.begin(), ts_values.end());
    }
    else{
       tts.add_points_sb(truth, ts_values.begin(), ts_values.end());
    }
    // check success rate and fail in case this is < 0.8. Do not check in case of stop_execution flag set:
    if(ts_values.size() * 1.0 / n_event < 0.8 && !stop_execution){
        throw FatalException("cls_limits: ts_producer fails in more than 20% of the cases");
    }
}


void cls_limits::run_single_truth_adaptive(double q, map<double, double> & truth_to_ts, double ts_epsilon, double truth){
    if(!tts.contains_truth_value(truth)){
        run_single_truth(truth, true, 100);
        run_single_truth(truth, false, 100);
        if(stop_execution) return;
    }
    for(int kk=0; kk<100; ++kk){
        if(stop_execution) return;
        update_truth_to_ts(truth_to_ts, q, ts_epsilon);
        double ts_value = truth_to_ts[truth];
        truth_ts_values::cls_info res = tts.get_cls(ts_value, truth);
        double cls_uncertainty = sqrt(pow(res.cls_uncertainty_n0, 2) + pow(res.cls_uncertainty_n, 2));
        debug_out << "run_single_truth_adaptive (iteration " << kk << "; truth = " << truth << "; ts = " << ts_value << "):\n";
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
            if(stop_execution) return;
        }
        if(res.cls_uncertainty_n > target_cls_uncertainty){
            size_t n = tts.get_n(truth);
            size_t n_new = pow(res.cls_uncertainty_n / target_cls_uncertainty, 2) * n;
            theta_assert(n_new >= n);
            size_t n_toys = max(50ul, min<size_t>(n_new - n, 500));
            debug_out << "run_single_truth_adaptive: making " << n_toys << " more toys for signal + background\n";
            flush(debug_out);
            run_single_truth(truth, false, n_toys);
            if(stop_execution) return;
        }
    }
    throw FatalException("run_single_truth_adaptive: could not reach target accuracy for CLs after 100 iterations");
}



cls_limits::cls_limits(const Configuration & cfg): vm(cfg.pm->get<VarIdManager>()), truth_parameter(vm->getParId(cfg.setting["truth_parameter"])),
  pid_limit(vm->createParId("__limit")), pid_lambda(vm->createParId("__lambda")), runid(1),
  expected_bands(true), limit_hint(NAN, NAN), reltol_limit(0.05), tol_cls(0.015), cl(0.95), n_toys(0), n_toy_errors(0) {
    SettingWrapper s = cfg.setting;
    
    //1. setup database and tables:
    db = PluginManager<Database>::build(Configuration(cfg, s["output_database"]));

    std::auto_ptr<Table> logtable_underlying = db->create_table("log");
    logtable.reset(new LogTable(logtable_underlying));
    
    std::auto_ptr<Table> rndinfo_table_underlying = db->create_table("rndinfo");
    rndinfo_table.reset(new RndInfoTable(rndinfo_table_underlying));
    cfg.pm->set("default", rndinfo_table);
    
    std::auto_ptr<Table> products_table_underlying = db->create_table("products");
    products_table.reset(new ProductsTable(products_table_underlying));
    string colname = s["ts_column"];
    sdc.reset(new SaveDoubleColumn(colname));
    std::vector<boost::shared_ptr<ProductsSink> > sinks;
    sinks.push_back(products_table);
    sinks.push_back(sdc);
    mps.reset(new MultiplexingProductsSink(sinks));
    cfg.pm->set<ProductsSink>("default", mps);
    
    boost::shared_ptr<int> ptr_runid(new int(runid));
    cfg.pm->set("runid", ptr_runid);
    
    cls_limits_table = db->create_table("cls_limits");
    cls_limits__ts_value = cls_limits_table->add_column("ts_value", typeDouble);
    cls_limits__limit = cls_limits_table->add_column("limit", typeDouble);
    cls_limits__limit_uncertainty = cls_limits_table->add_column("limit_uncertainty", typeDouble);
            
    //2. model, data_source and minimizer
    model = PluginManager<Model>::build(Configuration(cfg, s["model"]));
    minimizer = PluginManager<Minimizer>::build(Configuration(cfg, s["minimizer"]));
    source.reset(new data_filler(cfg, model));
    
    //3. logging stuff
    logtable->set_loglevel(LogTable::warning);
    
    //4. producer:
    producer = PluginManager<Producer>::build(Configuration(cfg, s["producer"]));
    pp_producer = dynamic_cast<ParameterDependentProducer*>(producer.get()); // either 0 or identical to producer ...
    if(pp_producer){
        ParIds pars = pp_producer->getParameters();
        if(pars.size() == 0){
           // not parameter dependent after all:
           pp_producer = 0;
        }
        else{
            if(pars.size() > 1 or not pars.contains(truth_parameter)){
                 throw ConfigurationException("cls_limits only works for Producers which depend on the truth parameter only (or have no parameter dependence at all)");
            }
        }
    }
    
    //5. data_source and expected_bands
    if(s.exists("data_source")){
        expected_bands = false;
        data_source = PluginManager<DataSource>::build(Configuration(cfg, s["data_source"]));
    }
    if(s.exists("expected_bands")){
        expected_bands = s["expected_bands"];
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
    if(s.exists("tol_cls")){
        tol_cls = s["tol_cls"];
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

// update the truth_to_ts map to contain ts values for all truth_values in tts.
//
// mode controls from where to take the ts value:
// < 0: from data
// 0 < q < 1: from the background-only distribution at that truth value
void cls_limits::update_truth_to_ts(map<double, double> & truth_to_ts, double q, double ts_epsilon){
    set<double> truth_values = tts.truth_values();
    Data data;
    bool data_filled = false;
    boost::optional<double> constant_ts; //filled only in case of constant, non-poi dependent ts, i.e., pp_producer == 0
    // in case of data and constant ts, check whether we have calculated the ts value already:
    if(pp_producer==0 && q < 0 && truth_to_ts.size() > 0){
        constant_ts = truth_to_ts.begin()->second;
    }
    for(set<double>::const_iterator it = truth_values.begin(); it!=truth_values.end(); ++it){
        if(q < 0){
            // in case we already know the (constant) ts, it's easy:
            if(constant_ts){
                truth_to_ts[*it] = *constant_ts;
                continue;
            }
            // otherwise, calculate ts value, be it constant or not:
            theta_assert(data_source.get());
            if(!data_filled) data_source->fill(data);
            if(pp_producer){
                ParValues vals;
                vals.set(truth_parameter, *it);
                pp_producer->setParameterValues(vals);
            }
            try {
                producer->produce(data, *model);
            } catch (Exception & ex) {
                ex.message += " (while running ts_producer for data)";
                throw;
            }
            double ts = sdc->get_value() + ts_epsilon;
            truth_to_ts[*it] = ts;
            // if ts value is constant, fill it to constant_ts, for the next iteration:
            if(pp_producer==0){
                constant_ts = ts;
            }
        }
        else{
            theta_assert(q > 0 && q < 1);
            double ts = tts.get_ts_b_quantile(*it, q);
            truth_to_ts[*it] = ts + ts_epsilon;
        }
    }
}


void cls_limits::run(){
    
    // 0. determine signal width
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


    run_single_truth(0.0, false, 50);
    run_single_truth(0.0, true, 50);

    //1. determine ts_epsilon. This is important for degenerate test statistic (such as likelihood ratios where
    // the backgruond only prior fixes the signal and is free in the s+b prior), because mathematically degenerate
    // values will not be necessarily exactly degenerate numerically. A possible solution to that is
    // to add small values to all 'observed' test statistic values, i.e., make them more signal-like (although
    // only by a very small amount).
    // truth0 is a guess for a 'high' truth value, which should be in about the right region of the limit.
    double truth0 = 2 * signal_width;
    if(isfinite(limit_hint.second)) truth0 = max(truth0, limit_hint.second);
    run_single_truth(truth0, false, 200);
    run_single_truth(truth0, true, 200);
    const double ts_epsilon = fabs(tts.get_ts_b_quantile(truth0, 0.975) - tts.get_ts_b_quantile(truth0, 0.025)) * 1e-3;
    debug_out << "ts_epsilon is " << ts_epsilon << "\n";


    fitexp_parameters pars(*minimizer, pid_limit, pid_lambda);
    const size_t N_maxit = 200;
    vector<double> qs;
    vector<double> all_limits;
    vector<double> all_limits_uncertainties;
    if(data_source.get()){
        qs.push_back(-1);
    }
    if(expected_bands){
        qs.push_back(0.025);
        qs.push_back(0.16);
        qs.push_back(0.5);
        qs.push_back(0.84);
        qs.push_back(0.975);
    }
    for(size_t iq=0; iq < qs.size(); ++iq){
        debug_out << "starting qs index " << iq << ", " << qs[iq] << "\n";
        flush(debug_out);
        map<double, double> truth_to_ts;
        // run an update, to get the ts values from previous iq iterations:
        update_truth_to_ts(truth_to_ts, qs[iq], ts_epsilon);
        //3.a. find a seed:
        //3.a.i. make some adaptive toys at the highest current truth value:
        cls_vs_truth_data data = tts.get_cls_vs_truth(truth_to_ts);
        run_single_truth_adaptive(qs[iq], truth_to_ts, ts_epsilon, data.truth_values().back());
        if(stop_execution) return;
        update_truth_to_ts(truth_to_ts, qs[iq], ts_epsilon);
        data = tts.get_cls_vs_truth(truth_to_ts);
        while(not cls_is_significantly_smaller(data, data.cls_values().size()-1, 1 - cl)){
             double next_value = data.truth_values().back()*2.0;
             debug_out << "making toys at high truth values to find upper limit on limit; next is truth=" << next_value << "\n";
             flush(debug_out);
             run_single_truth_adaptive(qs[iq], truth_to_ts, ts_epsilon, next_value);
             if(stop_execution) return;
             update_truth_to_ts(truth_to_ts, qs[iq], ts_epsilon);
             data = tts.get_cls_vs_truth(truth_to_ts);
        }
        size_t i_low, i_high;
        i_high = 0;
        i_low = data.cls_values().size()-1;
        while(not cls_is_significantly_smaller(data, i_high, 1-cl)) ++i_high; //terminates because we created a significantly smaller point in 3.a.i.
        while(not cls_is_significantly_larger(data, i_low, 1-cl)) --i_low; //terminates because this is true for i_low==0 (truth==0.0)
        double truth_low = data.truth_values()[i_low];
        double truth_high = data.truth_values()[i_high];
        debug_out << iq << " seed: interval [" << truth_low << ", " << truth_high << "]; i_low, ihigh = (" << i_low << ", " << i_high << ")\n";
        flush(debug_out);
        //3.b. iterate
        fitexp_result latest_res;
        for(size_t i=0; i<=N_maxit; ++i){
            debug_out << iq << "." << i << "\n";
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
            if(isnan(latest_res.limit)){
                debug_out << "exp fit did not work; fill in some random point in fitted interval, with large error\n";
                // u is a uniform random number between 0 and 1. Details don't play a role here, so just hard code a linear
                // congruent generator using i as seed:
                double u = ((i * 1103515245 + 12345) % 4294967296) * 1.0 / 4294967296;
                latest_res.limit = truth_low + (0.25 + 0.5 * u) * (truth_high - truth_low);
                latest_res.limit_error = 0.5 * (truth_high - truth_low);
            }
            debug_out << iq << "." << i << " result: limit = " << latest_res.limit << " +- " << latest_res.limit_error << "\n";
            theta_assert(latest_res.limit >= truth_low && latest_res.limit <= truth_high);
            if(latest_res.limit_error / latest_res.limit < reltol_limit){
                debug_out << iq << "." << i << ": result is good enough; done with this ts value.\n";
                break;
            }
            else{
                debug_out << iq << "." << i << ": result is NOT good enough yet: relative limit accuracy is "
                          << latest_res.limit_error / latest_res.limit << ", target accuracy is " << reltol_limit << ".\n";
            }
            flush(debug_out);
            // add a point at the estimated limit:
            run_single_truth_adaptive(qs[iq], truth_to_ts, ts_epsilon, latest_res.limit);
            if(stop_execution) return;
            update_truth_to_ts(truth_to_ts, qs[iq], ts_epsilon);
            data = tts.get_cls_vs_truth(truth_to_ts);
            // make interval [truth_low, truth_high] smaller if (border is far away (>2sigma) from limit or border is 0) AND the next point can serve as new interval border
            i_low = find(data.truth_values().begin(), data.truth_values().end(), truth_low) - data.truth_values().begin();
            i_high = find(data.truth_values().begin(), data.truth_values().end(), truth_high) - data.truth_values().begin();
            if((fabs(truth_low - latest_res.limit) / latest_res.limit_error > 2 or i_low==0) && i_low+1 < i_high){
                 debug_out << "lower border " << truth_low << " far away from current limit --> test if can be removed ...\n";
                 //re-run next candidate point:
                 run_single_truth_adaptive(qs[iq], truth_to_ts, ts_epsilon, data.truth_values()[i_low+1]);
                 if(stop_execution) return;
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
                 run_single_truth_adaptive(qs[iq], truth_to_ts, ts_epsilon, data.truth_values()[i_high-1]);
                 if(stop_execution) return;
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
        r.set_column(cls_limits__ts_value, qs[iq]);
        r.set_column(cls_limits__limit, latest_res.limit);
        r.set_column(cls_limits__limit_uncertainty, latest_res.limit_error);
        all_limits.push_back(latest_res.limit);
        all_limits_uncertainties.push_back(latest_res.limit_error);
        cls_limits_table->add_row(r);
    }
    debug_out << "limits:\n";
    for(size_t i=0; i<all_limits.size(); ++i){
        debug_out << "q = " << qs[i] << ": " << all_limits[i] << "  +-  " << all_limits_uncertainties[i] << "\n";
    }
}

REGISTER_PLUGIN(cls_limits)

