#ifndef PLUGINS_CLS_LIMITS
#define PLUGINS_CLS_LIMITS

#include "interface/decls.hpp"
#include "interface/main.hpp"
#include "interface/database.hpp"

#include <map>
#include <memory>
#include <boost/shared_ptr.hpp>

class MultiplexingProductsSink;
class SaveDoubleColumn;
class data_filler;
class truth_ts_values;

/*
 *   reuse_toys = (
 *       {type = "sqlite_database_in"; filename = "toys0.db"; }
 *   ); //optional; default is the empty list (), i.e., not to reuse any toys

 * (TODO: this is not implemented yet!) The optional setting \c reuse_toys is a list of theta::DatabaseInput
 * specifications. If present, it is assumed that these contain toys
 * from previous runs with the same model and ts definition (in particular, with the same truth_parameter and ts_column name ...).
 * Reusing toys can speed up the calculation substiantially.
*/



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
 *   ts_column = "lr__nll_diff";
 *   truth_parameter = "beta_signal";
 *
 *   // optional parameters
 *   cl = 0.90; // optional; default is 0.95
 *   data_source = { some datasource specification }; //optional
 *   expected_bands = true; //optional
 *
 *   reltol_limit = 0.001; // optional; default is 0.05
 *   tol_cls = 0.005; //optional; default is 0.015
 *   limit_hint = (200.0, 240.0); //optional; default is finding out on its own ...
 *
 *   debuglog = "debug.txt"; // optional, default is no debug output
 * };
 * \endcode
 * 
 * Toys are drawn from the \c model with different values for the signal parameter given in \c truth_parameter (all other
 * parameters are drawn randomly from the model's parameter distribution). For each toy, the \c producer is run and the
 * value from the given \c ts_column is used as test statistic for the construction of CLs limits.
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
 * The products produced by each toy is used immidiately for the construction of the CLs limits.
 * They are also written to the \c output_database. The tables are the same as for theta::Run. In addition, the table
 * "cls_limits" is created which contains the columns "q", "limit", and "limit_uncertainty".
 * The column "q" is the quantile of the background-only distribution used as test statistic. It is
 * The values 0.025, 0.16, 0.5, 0.84, and 0.975 define the "expected" bands. the special value -1 is used
 * to indicate that data was used from the data_source.
 * 
 * \c debuglog is the filename of the debug log. Debug output is only written if a filename is given here.
 *
 * The indicated progress and errors refer to the number of toys produced. It is hard to tell how many toys are
 * necessary; aside the more obvious things like \c reltol_limit, the number of required toys 
 * also depends on the CLb value: a low CLb value requires much more toys to reach the desired accuracy
 * (note that if calculating the expected limit band, low CLb values occur naturally).
 * In such a case, it is not uncommon that 20k to 50k toys are necessary. For a typical case of calculating the observed limit, i.e. if CLb is around 0.5,
 * usually 5-10k of toys are sufficient. This is also true in cases where the test statistic is degenerate as
 * in such a case, very low CLb values do not occur.
 */
class cls_limits: public theta::Main{
public:
    cls_limits(const theta::Configuration & cfg);
    virtual void run();
    
private:
    
    void run_single_truth(double truth, bool bkg_only, int n_event);
    
    // run toys at the given truth value (and at truth=0!) until either
    // * the uncertainty on the CLs value is below tol_cls
    //   or
    // * the CLs value is found to be more than 3 sigma away from the target CLs 1-cl.
    //void run_single_truth_adaptive(double ts_value, double truth);
    void run_single_truth_adaptive(double q, std::map<double, double> & truth_to_ts, double ts_epsilon, double truth);
    void update_truth_to_ts(std::map<double, double> & truth_to_ts, double q, double ts_epsilon);

    boost::shared_ptr<theta::VarIdManager> vm;
    boost::shared_ptr<theta::Model> model;
    //note: the tables hold a shared_ptr to this database to ensure proper destruction order
    boost::shared_ptr<theta::Database> db;

    std::auto_ptr<theta::LogTable> logtable;
    boost::shared_ptr<theta::RndInfoTable> rndinfo_table;

    //the producer to be run on the pseudo data which provides the test statistic:
    std::auto_ptr<theta::Producer> producer;
    theta::ParameterDependentProducer * pp_producer; // points to producer, if type matched. 0 otherwise.
    boost::shared_ptr<MultiplexingProductsSink> mps;
    boost::shared_ptr<SaveDoubleColumn> sdc;
    boost::shared_ptr<theta::ProductsTable> products_table;
    
    std::auto_ptr<theta::Table> cls_limits_table;
    theta::Column cls_limits__q, cls_limits__limit, cls_limits__limit_uncertainty;
    
    theta::ParId truth_parameter;
    std::auto_ptr<data_filler> source;
    
    // for fitting the CLs versus truth curve, we need a minimizer and some parameters:
    std::auto_ptr<theta::Minimizer> minimizer;
    theta::ParId pid_limit, pid_lambda;
    
    std::auto_ptr<truth_ts_values> tts;

    int runid;
    
    bool expected_bands;
    std::auto_ptr<theta::DataSource> data_source;
    
    std::pair<double, double> limit_hint;
    double reltol_limit, tol_cls;
    double cl;
    
    std::auto_ptr<std::ostream> debug_out;
    
    int n_toys, n_toy_errors;
};


#endif

