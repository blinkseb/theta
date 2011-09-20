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


# returns dictionary (spid0) -> ts value
def ts_data(model, ts = 'lr', signal_prior = 'flat', nuisance_prior = '', signal_prior_bkg='fix:0', signal_processes = None, **options):
    if signal_processes is None: signal_processes = [[sp] for sp in model.signal_processes]
    signal_prior = signal_prior_dict(signal_prior)
    signal_prior_bkg = signal_prior_dict(signal_prior_bkg)
    nuisance_prior = nuisance_prior_distribution(model, nuisance_prior)
    main = {'model': '@model', 'n-events': 1, 'producers': ('@ts_producer',), 'output_database': sqlite_database(), 'log-report': False}
    cfg_options = {'plugin_files': ('$THETA_DIR/lib/core-plugins.so','$THETA_DIR/lib/root.so')}
    ts_producer, ts_colname = ts_producer_dict(ts, signal_prior_bkg)    
    toplevel_settings = {'signal_prior': signal_prior, 'main': main, 'options': cfg_options, 'ts_producer': ts_producer, "minimizer": minimizer(need_error = False)}
    main['data_source'], toplevel_settings['model-distribution-signal'] = data_source_dict(model, 'data')
    for sp in signal_processes:
        model_parameters = model.get_parameters(sp)
        toplevel_settings['nuisance_prior'] = nuisance_prior.get_cfg(model_parameters)
        name = write_cfg(model, sp, 'ts_data', '', additional_settings = toplevel_settings)
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
    

# run toys at certain values of beta_signal to produce test statistic values; to be used later for cls_limits or discovery, etc.
# returns a dictionary (spid, beta signal value) -> (list of ts values)
def ts_toys(model, beta_signal_values, n, ts = 'lr', signal_prior = 'flat', nuisance_prior = '', signal_prior_bkg='fix:0', signal_processes = None, **options):
    if signal_processes is None: signal_processes = [[sp] for sp in model.signal_processes]
    signal_prior = signal_prior_dict(signal_prior)
    signal_prior_bkg = signal_prior_dict(signal_prior_bkg)
    nuisance_prior = nuisance_prior_distribution(model, nuisance_prior)
    main = {'model': '@model', 'n-events': n, 'producers': ('@ts_producer',), 'output_database': sqlite_database(), 'log-report': False}
    cfg_options = {'plugin_files': ('$THETA_DIR/lib/core-plugins.so','$THETA_DIR/lib/root.so')}
    ts_producer, ts_colname = ts_producer_dict(ts, signal_prior_bkg)
    toplevel_settings = {'signal_prior': signal_prior, 'main': main, 'options': cfg_options, 'ts_producer': ts_producer, "minimizer": minimizer(need_error = False)}
    main['data_source'], toplevel_settings['model-distribution-signal'] = data_source_dict(model, 'toys:0', **options)
    #print toplevel_settings
    cfg_names_to_run = []
    for beta_signal_value in beta_signal_values:
        for sp in signal_processes:
            model_parameters = model.get_parameters(sp)
            toplevel_settings['nuisance_prior'] = nuisance_prior.get_cfg(model_parameters)
            toplevel_settings['model-distribution-signal']['beta_signal'] = beta_signal_value
            name = write_cfg(model, sp, 'ts_toys', '%f' % beta_signal_value, additional_settings = toplevel_settings)
            cfg_names_to_run.append(name)
    if 'run_theta' not in options or options['run_theta']:
        run_theta(cfg_names_to_run)
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


# returns a tuple (expected significance, observed significance). Each entry in the tuple is in turn a triple (central value, error_minus, error_plus)
# where the error for the expected significance are the expected 1sigma spread and for the observed one the error from limited toy MC.
# If you think that latter error is too large, decrease Z_error_max.
def discovery(model, use_data = True, Z_error_max = 0.05, **options):
    options['signal_prior']='flat'
    options['signal_prior_bkg'] = 'fix:0'
    ts_expected = ts_toys(model, [1.0], 1000, **options)
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


# ts is either 'lr' for likelihood ratio or 'mle' for maximum likelihood estimate.
#
# for 'lr', there are different variants: (note: the names are taken from
# https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideHiggsAnalysisCombinedLimit )
# * 'LEP' will do a 'LEP-like' likelihood ratio in which no minimization is done at all; the likelihood ratio is evaluated using the
#    most probable values for the nuisance parameters and beta_signal=1, beta_signal=0.
#    this can be achieved by
#    nuisance_prior = 'shape:fix;rate:fix'
#    signal_prior = 'fix:1'
#    signal_prior_bkg = 'fix:0'
# * 'Tevatron' minimizes the likelihood ratio w.r.t. all nuisance parameters and use beta_signal=1, beta_signal=0
#    This can be achieved by using
#     nuisance_prior = ''
#     signal_prior = 'fix:1'
#     signal_prior_bkg = 'fix:0'
# * 'LHC' / 'ATLAS' method is similar to 'Tevatron' but will also let beta_signal free in the numerator of the likelihood ratio (and use beta_signal=0 for
#    the denominator). This is achieved by
#     nuisance_prior = ''
#     signal_prior = 'fix:1'
#     signal_prior_bkg = 'flat:[0,1]'
#
# For ts='mle', signal_prior_bkg has no meaning. The nuisance prior for the s+b / b only likelihood is always the same (so far).
#
# forwhat controls which limits are computed. Valid vaues are:
# * 'data': compute observed limit only
# * 'expected': compute +-1sigma, +-2sigma limit bands for background only
# * 'all': both 'data' and 'expected'
def cls_limits(model, forwhat = 'all',  cl = 0.95, ts = 'lr', signal_prior = 'flat', nuisance_prior = '', signal_prior_bkg='fix:0', signal_processes = None, **options):
    if signal_processes is None: signal_processes = [[sp] for sp in model.signal_processes]
    # if the signal+background distribution is more restrictive than the background only, there is a sign flip ...
    #if signal_prior_bkg.startswith('flat') and signal_prior.startswith('fix'): ts_sign = '-'
    #else: ts_sign=''
    ts_sign = ''
    signal_prior = signal_prior_dict(signal_prior)
    signal_prior_bkg = signal_prior_dict(signal_prior_bkg)
    nuisance_prior = nuisance_prior_distribution(model, nuisance_prior)
    minimizer = minimizer()
    main = {'type': 'cls_limits', 'model': '@model', 'producers': ('@ts_producer',),
        'output_database': sqlite_database(), 'truth_parameter': 'beta_signal', 'minimizer': minimizer, 'debuglog': 'debug.txt'}
    if forwhat in ('expected', 'all'):
        main['ts_values_background_bands'] = True
    elif forwhat != 'data': raise RuntimeError, "unknown option forwhat='%s'" % forwhat
    if forwhat in ('data', 'all'):
        main['ts_values'], dummy = data_source_dict(model, 'data')
    cfg_options = {'plugin_files': ('$THETA_DIR/lib/core-plugins.so','$THETA_DIR/lib/root.so')}
    if ts == 'mle':
        ts_producer = {'type': 'mle', 'name': 'mle', 'minimizer': minimizer, 'parameter': 'beta_signal',
           'override-parameter-distribution': product_distribution("@signal_prior", "@nuisance_prior")}
        main['ts_column'] = ts_sign + 'mle__beta_signal'
    elif ts == 'lr':
        ts_producer = {'type': 'deltanll_hypotest', 'name': 'lr', 'minimizer': minimizer,
        'signal-plus-background-distribution': product_distribution("@signal_prior", "@nuisance_prior"),
        'background-only-distribution': product_distribution(signal_prior_bkg, "@nuisance_prior")}
        main['ts_column'] = ts_sign + 'lr__nll_diff'
    else: raise RuntimeError, 'unknown ts "%s"' % ts
    toplevel_settings = {'signal_prior': signal_prior, 'main': main, 'options': cfg_options, 'ts_producer': ts_producer,
       'model-distribution-signal': {'type': 'flat_distribution', 'beta_signal': {'range': [0.0, float("inf")], 'fix-sample-value': 0.0}}}
    cfg_names_to_run = []
    for sp in signal_processes:
        model_parameters = model.get_parameters(sp)
        toplevel_settings['nuisance_prior'] = nuisance_prior.get_cfg(model_parameters)
        name = write_cfg(model, sp, 'cls_limits', '', additional_settings = toplevel_settings)
        cfg_names_to_run.append(name)
    if 'run_theta' not in options or options['run_theta']:
        run_theta(cfg_names_to_run, **options)
    else: return None
    
    cachedir = os.path.join(config.workdir, 'cache')
    
    result = {}
    #result_table = Report.table()
    #result_table.add_column('process', 'signal process')
    #result_table.add_column('limit', header)
    for name in cfg_names_to_run:
        method, sp, dummy = name.split('-',2)
        #result_table.set_column('process', sp)
        sqlfile = os.path.join(cachedir, '%s.db' % name)
        data = sql(sqlfile, 'select "limit", limit_uncertainty from cls_limits order by "index"')
        for row in data:
            print row
        #result_table.add_row()
    #config.report.new_section("CLs limits on ensemble '%s'" % input)
    #config.report.add_html(result_table.html())
    return result

