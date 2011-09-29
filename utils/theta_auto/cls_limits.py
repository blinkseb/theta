import config, utils, os.path, datetime, math, bisect

import scipy.special

from theta_interface import *
import Report
import bayesian

from utils import *

# returns a tuple (theta config dictionary, ts column name) for ts producer specification ts (either 'lr' or 'mle'). Refers
# to the global settings "signal_prior", "nuisance_prior" and "minimizer"
def ts_producer_dict(ts, signal_prior_bkg='flat'):
    signal_prior_bkg = signal_prior_dict(signal_prior_bkg)
    if ts == 'mle':
        ts_producer = {'type': 'mle', 'name': 'mle', 'minimizer': "@minimizer", 'parameter': 'beta_signal',
           'override-parameter-distribution': product_distribution("@signal_prior", "@nuisance_prior")}
        ts_colname = 'mle__beta_signal'
    elif ts == 'lr':
        ts_producer = {'type': 'deltanll_hypotest', 'name': 'lr', 'minimizer': "@minimizer",
        'signal-plus-background-distribution': product_distribution("@signal_prior", "@nuisance_prior"),
        'background-only-distribution': product_distribution(signal_prior_bkg, "@nuisance_prior")}
        ts_colname = 'lr__nll_diff'
    else: raise RuntimeError, 'unknown ts "%s"' % ts
    return ts_producer, ts_colname

## \brief Calculate the test statistic value for data
#
# \return dictionary (spid) -> ts value
def ts_data(model, ts = 'lr', signal_prior = 'flat', nuisance_prior = '', signal_prior_bkg='fix:0', signal_processes = None, **options):
    if signal_processes is None: signal_processes = [[sp] for sp in model.signal_processes]
    signal_prior = signal_prior_dict(signal_prior)
    signal_prior_bkg = signal_prior_dict(signal_prior_bkg)
    nuisance_prior = nuisance_prior_distribution(model, nuisance_prior)
    main = {'model': '@model', 'n-events': 1, 'producers': ('@ts_producer',), 'output_database': sqlite_database(), 'log-report': False}
    ts_producer, ts_colname = ts_producer_dict(ts, signal_prior_bkg)    
    toplevel_settings = {'signal_prior': signal_prior, 'main': main, 'ts_producer': ts_producer, "minimizer": minimizer(need_error = False)}
    toplevel_settings.update(get_common_toplevel_settings(**options))
    main['data_source'], toplevel_settings['model-distribution-signal'] = data_source_dict(model, 'data')
    for sp in signal_processes:
        model_parameters = model.get_parameters(sp)
        toplevel_settings['nuisance_prior'] = nuisance_prior.get_cfg(model_parameters)
        name = write_cfg(model, sp, 'ts_data', '', additional_settings = toplevel_settings, **options)
        cfg_names_to_run = []
        cfg_names_to_run.append(name)
    if 'run_theta' not in options or options['run_theta']:
        run_theta(cfg_names_to_run)
    else: return None
    
    cachedir = os.path.join(config.workdir, 'cache')
    result = {}
    for name in cfg_names_to_run:
        method, sp, dummy = name.split('-',2)
        sqlfile = os.path.join(cachedir, '%s.db' % name)
        data = sql(sqlfile, 'select "%s" from products' % ts_colname)
        result[sp] = data[0][0]
    return result
 
## \brief Calculate the test statistic for toys with certain values of beta_signal
#
# The result can be used later for cls_limits, discovery, etc.
# This function is usually not used directly but called internally, e.g., from \ref discovery.
#
# \return a dictionary (spid, beta signal value) -> (list of ts values)
def ts_toys(model, beta_signal_values, n, ts = 'lr', signal_prior = 'flat', nuisance_prior = '', signal_prior_bkg='fix:0', signal_processes = None, **options):
    if signal_processes is None: signal_processes = [[sp] for sp in model.signal_processes]
    signal_prior = signal_prior_dict(signal_prior)
    signal_prior_bkg = signal_prior_dict(signal_prior_bkg)
    nuisance_prior = nuisance_prior_distribution(model, nuisance_prior)
    main = {'model': '@model', 'n-events': n, 'producers': ('@ts_producer',), 'output_database': sqlite_database(), 'log-report': False}
    ts_producer, ts_colname = ts_producer_dict(ts, signal_prior_bkg)
    toplevel_settings = {'signal_prior': signal_prior, 'main': main, 'ts_producer': ts_producer, "minimizer": minimizer(need_error = False)}
    toplevel_settings.update(get_common_toplevel_settings(**options))
    main['data_source'], toplevel_settings['model-distribution-signal'] = data_source_dict(model, 'toys:0', **options)
    #print toplevel_settings
    cfg_names_to_run = []
    for beta_signal_value in beta_signal_values:
        for sp in signal_processes:
            model_parameters = model.get_parameters(sp)
            toplevel_settings['nuisance_prior'] = nuisance_prior.get_cfg(model_parameters)
            toplevel_settings['model-distribution-signal']['beta_signal'] = beta_signal_value
            name = write_cfg(model, sp, 'ts_toys', '%f' % beta_signal_value, additional_settings = toplevel_settings, **options)
            cfg_names_to_run.append(name)
    if 'run_theta' not in options or options['run_theta']:
        run_theta(cfg_names_to_run, **options)
    else: return None
    cachedir = os.path.join(config.workdir, 'cache')
    result = {}
    for name in cfg_names_to_run:
        method, sp, s_beta_signal, dummy = name.split('-',3)
        beta_signal = float(s_beta_signal)
        sqlfile = os.path.join(cachedir, '%s.db' % name)
        data = sql(sqlfile, 'select "%s" from products' % ts_colname)
        result[(sp,beta_signal)] = [row[0] for row in data]
    return result


# returns tuple (Z, Z_error), where Z_error is from the limited number of toys (the larger one in Z):
def p_to_Z(n, n0):
    p = n*1.0 / n0
    p_error = max(math.sqrt(p*(1-p) / n0), 1.0 / n0)
    p_plus = min(p + p_error, 1.0 - p_error)
    p_minus = max(p - p_error, p_error)
    Z = math.sqrt(2) * scipy.special.erfinv(1 - 2*p)
    Z_p_plus = math.sqrt(2) * scipy.special.erfinv(1 - 2*p_plus)
    Z_p_minus = math.sqrt(2) * scipy.special.erfinv(1 - 2*p_minus)
    return Z, max(abs(Z - Z_p_minus), abs(Z - Z_p_plus))

def count_above(the_list, threshold):
    n=0
    for e in the_list:
        if e >= threshold: n+=1
    return n

## \brief Determine p-value / "N sigma" from tail distribution of background-only test statistic
#
# The number of toys will be increased adaptively.
#
# \param Z_error_max The maximum error in the Z value
#
# \return A tuple (expected significance, observed significance). Each entry in the tuple is in turn a triple (central value, error_minus, error_plus)
# where the error for the expected significance are the expected 1sigma spread and for the observed one the error from limited toy MC.
def discovery(model, use_data = True, Z_error_max = 0.05, **options):
    options['signal_prior']='flat'
    options['signal_prior_bkg'] = 'fix:0'
    ts_expected = ts_toys(model, [1.0], 10000, **options)
    # make (spid) -> (median, -1sigma, +1sigma):
    expected = {}
    for spid, beta_signal in ts_expected:
        ts_sorted = sorted(ts_expected[(spid, beta_signal)])
        expected[spid] = (ts_sorted[int(0.5 * len(ts_sorted))], ts_sorted[int(0.16 * len(ts_sorted))], ts_sorted[int(0.84 * len(ts_sorted))])
    del ts_expected
    if use_data: observed = ts_data(model, **options)
    n = 10000
    # dict  spid -> [n, n0] for observed p-value
    observed_nn0 = {}
    # dict spid -> (median [n, n0], -1sigma [n, n0], +1sigma [n, n0])
    expected_nn0 = {}
    for spid in expected:
        expected_nn0[spid] = ([0,0], [0,0], [0,0])
        observed_nn0[spid] = [0,0]
    observed_significance = None
    for seed in range(1, 100):
        options['toydata_seed'] = seed
        ts_bkgonly = ts_toys(model, [0.0], n, **options)
        max_Z_error = 0.0
        for spid, beta_signal in ts_bkgonly:
            Z = [0,0,0]
            for i in range(3):
                expected_nn0[spid][i][1] += len(ts_bkgonly[spid, beta_signal])
                expected_nn0[spid][i][0] += count_above(ts_bkgonly[spid, beta_signal], expected[spid][i])
                Z[i], Z_error = p_to_Z(*expected_nn0[spid][i])
                max_Z_error = max(max_Z_error, Z_error)
            print "expected significance for process '%s' (with +-1sigma band): %f  +%f  -%f"  % (spid, Z[0], Z[2] - Z[0], Z[0]-Z[1])
            expected_significance = Z[0], Z[2] - Z[0], Z[0]-Z[1]
            if use_data:
                observed_nn0[spid][1] += len(ts_bkgonly[spid, beta_signal])
                observed_nn0[spid][0] += count_above(ts_bkgonly[spid, beta_signal], observed[spid])
                Z, Z_error = p_to_Z(*observed_nn0[spid])
                max_Z_error = max(max_Z_error, Z_error)
                observed_significance = Z, Z_error, Z_error
                print "observed significance for process '%s' (with stat. error from limited toys): %f +- %f" % (spid, Z, Z_error)
        #print "max Z error: ", max_Z_error
        if max_Z_error < Z_error_max: break
    return expected_significance, observed_significance


## \brief Calculate CLs limits.
#
# Calculate expected and/or observed CLs limits for a given model. The choice of test statistic
# is driven by the 'ts' option which can be either 'lr' for a likelihood ratio test statistic
# calculated with the deltanll_hypotest producer or 'mle' for the maximum likelihood estimate
# of beta_signal to be used as test statistic. In both cases, 'nuisance_prior' controls the nuisance
# prior to be used in the likelihood function definition and 'signal_prior' is the constraint for signal.
#
# In case of ts='lr', you can also set 'signal_prior_bkg' which is the constraint used for the evaluation of the
# profiled likelihood in the background only case.
#
# For ts='mle', signal_prior_bkg has no meaning. The nuisance prior for the s+b / b only likelihood is always the same (so far).
#
# what controls which limits are computed. Valid vaues are:
# * 'observed': compute observed limit on data only
# * 'expected': compute +-1sigma, +-2sigma limit bands for background only
# * 'all': both 'data' and 'expected'
#
# returns a tuple of two plotutil.plotdata instances. The first contains expected limit (including the band) and the second the 'observed' limit
# if 'what' is not 'all', one of the plotdata instances is replaced with None.
def cls_limits(model, what = 'all',  cl = 0.95, ts = 'lr', signal_prior = 'flat', nuisance_prior = '', signal_prior_bkg='fix:0', signal_processes = None, **options):
    if signal_processes is None: signal_processes = [[sp] for sp in model.signal_processes]
    # if the signal+background distribution is more restrictive than the background only, there is a sign flip ...
    #if signal_prior_bkg.startswith('flat') and signal_prior.startswith('fix'): ts_sign = '-'
    #else: ts_sign=''
    ts_sign = ''
    signal_prior = signal_prior_dict(signal_prior)
    if 'beta_signal' in signal_prior: signal_prior['beta_signal']['fix-sample-value'] = signal_prior['beta_signal']['range'][0]
    signal_prior_bkg = signal_prior_dict(signal_prior_bkg)
    nuisance_prior = nuisance_prior_distribution(model, nuisance_prior)
    main = {'type': 'cls_limits', 'model': '@model', 'producers': ('@ts_producer',),
        'output_database': sqlite_database(), 'truth_parameter': 'beta_signal', 'minimizer': minimizer(need_error = True), 'debuglog': 'debug.txt'}
    if what in ('expected', 'all'):
        main['ts_values_background_bands'] = True
    elif what != 'observed': raise RuntimeError, "unknown option what='%s'" % what
    if what in ('observed', 'all'):
        main['ts_values'], dummy = data_source_dict(model, 'data')
    if ts == 'mle':
        ts_producer = {'type': 'mle', 'name': 'mle', 'minimizer': minimizer(need_error = False), 'parameter': 'beta_signal',
           'override-parameter-distribution': product_distribution("@signal_prior", "@nuisance_prior")}
        main['ts_column'] = ts_sign + 'mle__beta_signal'
    elif ts == 'lr':
        ts_producer = {'type': 'deltanll_hypotest', 'name': 'lr', 'minimizer': minimizer(need_error = False),
        'signal-plus-background-distribution': product_distribution("@signal_prior", "@nuisance_prior"),
        'background-only-distribution': product_distribution(signal_prior_bkg, "@nuisance_prior")}
        main['ts_column'] = ts_sign + 'lr__nll_diff'
    else: raise RuntimeError, 'unknown ts "%s"' % ts
    toplevel_settings = {'signal_prior': signal_prior, 'main': main, 'ts_producer': ts_producer,
       'model-distribution-signal': {'type': 'flat_distribution', 'beta_signal': {'range': [0.0, float("inf")], 'fix-sample-value': 0.0}}}
    toplevel_settings.update(get_common_toplevel_settings(**options))
    cfg_names_to_run = []
    for sp in signal_processes:
        model_parameters = model.get_parameters(sp)
        toplevel_settings['nuisance_prior'] = nuisance_prior.get_cfg(model_parameters)
        name = write_cfg(model, sp, 'cls_limits', '', additional_settings = toplevel_settings, **options)
        cfg_names_to_run.append(name)
    if 'run_theta' not in options or options['run_theta']:
        run_theta(cfg_names_to_run, **options)
    else: return None
    cachedir = os.path.join(config.workdir, 'cache')
    
    
    #expected results maps (process name) -> (median, band1, band2)
    # where band1 and band2 are tuples for the central 68 and 95%, resp.
    expected_results = {}
    # observed result maps (process name) -> (limit, limit_uncertainty)
    observed_results = {}
    spids = set()    
    for name in cfg_names_to_run:
        method, sp, dummy = name.split('-',2)
        spids.add(sp)
        sqlfile = os.path.join(cachedir, '%s.db' % name)
        data = sql(sqlfile, 'select "limit", limit_uncertainty from cls_limits order by "index"')
        index = 0
        if what in ('all', 'observed'):
            observed_results[sp] = data[index][0], data[index][1]
            index += 1
        if what in ('all', 'expected'):
            expected_results[sp] = (data[index+2][0], (data[index+1][0], data[index+3][0]), (data[index][0], data[index+4][0]))
    x_to_sp = get_x_to_sp(spids, **options)
    
    pd_expected, pd_observed = None, None
    result_table = Report.table()
    result_table.add_column('process', 'signal process')
    if what in ('all', 'expected'):
        pd_expected = plotutil.plotdata()
        pd_expected.color = '#000000'
        pd_expected.as_function = True
        pd_expected.legend = 'expected limit'
        pd_expected.x = sorted(list(x_to_sp.keys()))
        pd_expected.y  = []
        pd_expected.bands = [([], [], '#00ff00'), ([], [], '#00aa00')]
        result_table.add_column('expected limit', 'expected limit (median)')
        result_table.add_column('expected limit (1sigma band)')
        result_table.add_column('expected limit (2sigma band)')
    if what in ('all', 'observed'):
        pd_observed = plotutil.plotdata()
        pd_observed.color = '#000000'
        pd_observed.as_function = True
        pd_observed.x = sorted(list(x_to_sp.keys()))
        pd_observed.y = []
        pd_observed.legend = 'observed limit'
        result_table.add_column('observed limit')
    
    for x in sorted(x_to_sp.keys()):
        sp = x_to_sp[x]
        result_table.set_column('process', sp)
        if what in ('all', 'expected'):
            expected, band1, band2 = expected_results[sp]
            pd_expected.y.append(expected)
            pd_expected.bands[1][0].append(band1[0])
            pd_expected.bands[1][1].append(band1[1])
            pd_expected.bands[0][0].append(band2[0])
            pd_expected.bands[0][1].append(band2[1])
            result_table.set_column('expected limit', '%.2f' % expected)
            result_table.set_column('expected limit (1sigma band)', '[%.2f, %.2f]' % band1)
            result_table.set_column('expected limit (2sigma band)', '[%.2f, %.2f]' % band2)
        if what in ('all', 'observed'):
            observed, observed_unc = observed_results[sp]
            pd_observed.y.append(observed)
            result_table.set_column('observed limit', '%.2f +- %.2f' % (observed, observed_unc))
        result_table.add_row()
    bayesian.report_limit_band_plot(pd_expected, pd_observed, 'CLs', 'cls')
    return pd_expected, pd_observed


