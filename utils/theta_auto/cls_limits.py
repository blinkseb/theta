# -*- coding: utf-8 -*-
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
    elif ts == 'lhclike':
        ts_producer = {'type': 'deltanll_hypotest', 'name': 'lr', 'minimizer': "@minimizer",
        'signal-plus-background-distribution': product_distribution(signal_prior_bkg, "@nuisance_prior"),
        'background-only-distribution': product_distribution(signal_prior_bkg, "@nuisance_prior"),
        'restrict_poi': 'beta_signal'}
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
# In case of ts='lr', you can also set 'signal_prior_bkg' which is the constraint used for the evaluation of the
# profiled likelihood in the background only case. Setting signal_prior_bkg to None chooses the default which is
#  * 'flat' in case ts=='lhclike'
#  * 'fix:0' otherwise
#
# The return value depends on 'ret': in case of 'data', returns a dictionary (spid, beta signal value) -> (list of ts values).
# If ret is 'cfg_names', returns a list of config file names.
# All other values of ret are ignored for now; in that case, None is returned
#
# bkg_only_toys is a special setting for the lhclike test statistic which depends on beta_signal. If set to true, the toys
# are always background only and the beta_signal_values are only used in the evaluation of the test statistic. Otherwise,
# the beta_signal_values are used both, for toy generation, and as parameter for evaluating the test statistic.
#
def ts_toys(model, beta_signal_values, n, ts = 'lr', signal_prior = 'flat', nuisance_prior = '', signal_prior_bkg = None, signal_processes = None, bkg_only_toys = False, ret = 'data', **options):
    if signal_processes is None: signal_processes = [[sp] for sp in model.signal_processes]
    if signal_prior_bkg is None:
        if ts == 'lhclike': signal_prior_bkg = 'flat'
        else: signal_prior_bkg = 'fix:0'
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
        if ts=='lhclike': ts_producer['default_poi_value'] = beta_signal_value
        for sp in signal_processes:
            model_parameters = model.get_parameters(sp)
            toplevel_settings['nuisance_prior'] = nuisance_prior.get_cfg(model_parameters)
            toplevel_settings['model-distribution-signal']['beta_signal'] = beta_signal_value
            if bkg_only_toys: toplevel_settings['model-distribution-signal']['beta_signal'] = 0.0
            name = write_cfg(model, sp, 'ts_toys', '%f' % beta_signal_value, additional_settings = toplevel_settings, **options)
            cfg_names_to_run.append(name)
    run_theta_ = options.get('run_theta', True)
    if run_theta_: run_theta(cfg_names_to_run, **options)
    else: return cfg_names_to_run
    if ret == 'data':
        cachedir = os.path.join(config.workdir, 'cache')
        result = {}
        for name in cfg_names_to_run:
            method, sp, s_beta_signal, dummy = name.split('-',3)
            beta_signal = float(s_beta_signal)
            sqlfile = os.path.join(cachedir, '%s.db' % name)
            data = sql(sqlfile, 'select "%s" from products' % ts_colname)
            result[(sp,beta_signal)] = [row[0] for row in data]
        return result
    elif ret=='cfg_names': return cfg_names_to_run


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


class truth_ts_values:
    def __init__(self):
        self.truth_to_ts_sb = {}
        self.truth_to_ts_b = {}

    def truth_values(self):
        result = set(self.truth_to_ts_sb.keys())
        if 0.0 in result: result.remove(0.0)
        return result

    def add_point_b(self, truth, ts):
        if truth not in self.truth_to_ts_b: self.truth_to_ts_b[truth] = []
        self.truth_to_ts_b[truth].append(ts)

    def add_point_sb(self, truth, ts):
        if truth not in self.truth_to_ts_sb: self.truth_to_ts_sb[truth] = []
        self.truth_to_ts_sb[truth].append(ts)


    def get_truth_to_ts(self, q):
        assert q > 0 and q < 1, str(q)
        result = {}
        for t in self.truth_values():
            b = sorted(self.truth_to_ts_b[t])
            ts = b[int(q*len(b))]
            result[t] = ts
        return result

    def get_cls_vs_truth(self, truth_to_ts):
        pd = plotutil.plotdata()
        pd.color = '#000000'
        pd.draw_line = False
        pd.x = sorted(list(self.truth_values()))
        pd.y = []
        pd.as_function = True
        pd.yerrors = []
        for t in pd.x:
            ts = truth_to_ts[t]
            n_b = len(self.truth_to_ts_b[t])
            clb = len([x for x in self.truth_to_ts_b[t] if x < ts]) * 1.0 / n_b
            clb_error = max(math.sqrt(clb * (1 - clb) / n_b), 1.0 / n_b)
            n_sb = len(self.truth_to_ts_sb[t])
            clsb = len([x for x in self.truth_to_ts_sb[t] if x < ts]) * 1.0 / n_sb
            clsb_error = max(math.sqrt(clsb * (1 - clsb) / n_sb), 1.0 / n_sb)
            if clb==0: clb += clb_error
            cls = clsb / clb
            pd.y.append(cls)
            cls_error = math.sqrt((clsb_error / clb)**2 + (clb_error * cls / clb)**2)
            pd.yerrors.append(cls_error)
        return pd
            

# make some plots of the test statistic for the report, using the given db filename
def debug_cls_plots(dbfile, ts_column):
    limits = sql(dbfile, 'select "index", "limit", limit_uncertainty from cls_limits')
    indices = [row[0] for row in limits]
    have_data = 0 in indices
    data = sql(dbfile, 'select runid, lr__poi, source__truth, "%s" from products' % ts_column)
    tts = truth_ts_values()
    for row in data:
        if row[0] == 0: continue
        if row[2] > 0:  tts.add_point_sb(row[2], row[3])
        else: tts.add_point_b(row[1], row[3])
    truth_values = sorted(list(tts.truth_values()))
    plotsdir = os.path.join(config.workdir, 'plots')
    config.report.new_section('debug_cls for file %s' % dbfile)
    result_table = Report.table()
    #result_table.add_column('q')
    result_table.add_column('limit')
    result_table.add_column('limit_uncertainty')
    for row in limits:
        #result_table.set_column('q', '%g' % row[0])
        result_table.set_column('limit', '%g' % row[1])
        result_table.set_column('limit_uncertainty', '%g' % row[2])
        result_table.add_row()
    config.report.add_html(result_table.html())
    pd = tts.get_cls_vs_truth(truth_to_ts_data)
    plotname = 'debug_cls-cls-vs-truth-data.png'
    plot([pd], 'beta_signal', 'CLs', os.path.join(plotsdir, plotname))
    config.report.add_html('<p>For data: <br/><img src="plots/%s" /></p>' % plotname)

"""
    # plot the test statistic value for data versus truth value (constant, unless lhclike test stat.)
    pds = []
    i = 0
    colors = ['#ff0000', '#00ff00',  '#0000ff', '#aaaaaa', '#000000', '#00ffff']
    for q in qs:
        if q == -1.0:
            truth_to_ts = {}
            tsdata = sorted([row for row in data if row[0]==0], cmp=lambda r1, r2: cmp(r1[1], r2[1]))
            for r in tsdata: truth_to_ts[r[1]] = r[3]
            truth_to_ts_data = truth_to_ts
        else:
            truth_to_ts = tts.get_truth_to_ts(q)
        pd_ts_vs_truth = plotutil.plotdata()
        pd_ts_vs_truth.x = sorted(truth_to_ts.keys())
        pd_ts_vs_truth.y = [truth_to_ts[t] for t in pd_ts_vs_truth.x]
        pd_ts_vs_truth.as_function = True
        pd_ts_vs_truth.color = colors[i]
        i+=1
        pd_ts_vs_truth.legend = 'q = %g' % q
        pds.append(pd_ts_vs_truth)

    plotname = 'debug_cls-ts-vs-beta_signal.png'
    plot(pds, 'beta_signal', 'ts', os.path.join(plotsdir, plotname))
    config.report.add_html('<p>test statistic vs beta_signal:<br/><img src="plots/%s" /></p>' % plotname)
    # plot the test statistic distributions for b-only and s+b; one plot per truth value
    i=0
    for t in truth_values: 
        l_sb = tts.truth_to_ts_sb[t]
        l_b = tts.truth_to_ts_b[t]
        ts_min, ts_max = min(min(l_b), min(l_sb)), max(max(l_b), max(l_sb))
        pd_sb = plotutil.plotdata()
        pd_sb.histogram(l_sb, ts_min, ts_max, 100)
        pd_sb.legend = 's+b'
        pd_sb.color = '#ff0000'
        pd_b = plotutil.plotdata()
        pd_b.histogram(l_b, ts_min, ts_max, 100)
        pd_b.legend = 'b'
        pd_b.color = '#00ff00'
        i+=1
        plotname = 'debug_cls-cls-histo%d.png' % i
        plot([pd_b, pd_sb], 'ts', 'N', os.path.join(plotsdir, plotname), logy=True, ymin=0.5)
        config.report.add_html('<p>For truth = %g; n_sb = %d; n_b = %d: <br/><img src="plots/%s" /></p>' % (t, len(l_sb), len(l_b),  plotname))
    i=0
    for q in qs:
        if q < 0: pd = tts.get_cls_vs_truth(truth_to_ts_data)
        else: pd = tts.get_cls_vs_truth(tts.get_truth_to_ts(q))
        i+=1
        plotname = 'debug_cls-cls-vs-truth%d.png' % i
        plot([pd], 'beta_signal', 'CLs', os.path.join(plotsdir, plotname))
        config.report.add_html('<p>For q = %g: <br/><img src="plots/%s" /></p>' % (q,  plotname))
"""


## \brief Make a grid for CLs limit calculation
#
#
# To make use of distributed grid calculation, three steps are necessary:
# 1. In your analysis.py, call cls_limits_grid with run_theta=False. This returns a dictionary containing all config filenames to run. As usual, you find them in the working directory.
# 2. Run (somewhere, distributed), all config filenames and put the output together with the original config file in the cache directory.
# 3. In analysis.py, call cls_limits, passing the dictionary from the (same!) call to cls_limits_grid as the reuse_toys parameter.
#
# Returns a dictionary (spid) -> (list of config names). This dictionary can be passed as the 'reuse_toys' parameter to cls_limits
#
# If options include 'run_theta = False', theta was not actually run. For each signal process
# group, n_configs configuration files are created.
#
def cls_limits_grid(model, beta_signal_range, n_beta_signal = 21, n_toys_sb = 2000, n_toys_b = 2000, n_configs = 1, ts = 'lr', signal_prior = 'flat',\
        nuisance_prior = '', signal_prior_bkg = None, signal_processes = None, **options):
    if signal_processes is None: signal_processes = [[sp] for sp in model.signal_processes]
    assert len(signal_processes) > 0
    if signal_prior_bkg is None:
        if ts == 'lhclike': signal_prior_bkg = 'flat'
        else: signal_prior_bkg = 'fix:0'
    signal_prior = signal_prior_dict(signal_prior)
    signal_prior_bkg = signal_prior_dict(signal_prior_bkg)
    nuisance_prior = nuisance_prior_distribution(model, nuisance_prior)
    main = {'type': 'cls_limits', 'model': '@model', 'producer': '@ts_producer', 'rnd_gen': {},
        'output_database': sqlite_database(), 'truth_parameter': 'beta_signal', 'mode': 'generate_grid',
        'truth_range': (float(beta_signal_range[0]), float(beta_signal_range[1])), 'n_truth': n_beta_signal,
        'n_sb_toys_per_truth': n_toys_sb, 'n_b_toys_per_truth': n_toys_b }
    ts_producer, tscol = ts_producer_dict(ts, signal_prior_bkg)
    toplevel_settings = {'signal_prior': signal_prior, 'main': main, 'ts_producer': ts_producer, 'minimizer': minimizer(need_error = False),
       'model-distribution-signal': {'type': 'flat_distribution', 'beta_signal': {'range': [0.0, float("inf")], 'fix-sample-value': 0.0}}}
    if ts == 'lhclike': del toplevel_settings['signal_prior']
    toplevel_settings.update(get_common_toplevel_settings(**options))
    cfg_names_to_run = []
    result = {}
    for sp in signal_processes:
        spid = ''.join(sp)
        if spid not in result: result[spid] = []
        model_parameters = model.get_parameters(sp)
        toplevel_settings['nuisance_prior'] = nuisance_prior.get_cfg(model_parameters)
        for i in range(n_configs):
            main['rnd_gen']['seed'] = i + 1
            name = write_cfg(model, sp, 'cls_limits_grid', '%d' % i, additional_settings = toplevel_settings, **options)
            cfg_names_to_run.append(name)
            result[spid].append(name)
    run_theta_ = options.get('run_theta', True)
    if run_theta_:
        run_theta(cfg_names_to_run, **options)
    return result

## \brief Calculate CLs limits.
#
# Calculate expected and/or observed CLs limits for a given model. The choice of test statistic
# is driven by the 'ts' option which can be either 'lr' for a likelihood ratio test statistic
# calculated with the deltanll_hypotest producer or 'mle' for the maximum likelihood estimate
# of beta_signal to be used as test statistic. In both cases, 'nuisance_prior' controls the nuisance
# prior to be used in the likelihood function definition and 'signal_prior' is the constraint for signal.
#
# In case of ts='lr', you can also set 'signal_prior_bkg' which is the constraint used for the evaluation of the
# profiled likelihood in the background only case. Setting signal_prior_bkg to None chooses the default which is
#  * 'flat' in case ts=='lhclike'
#  * 'fix:0' otherwise
#  
# For LHC-like test statistic, use ts = 'lhclike'. In this case, signal_prior has no effect: if fitting
# the s+b model, beta_signal is always fixed to r. 'signal_prior_bkg' is applied. You would usually set it to 'flat' in this case
# so that the b-only fit varies beta_signal on [0, r].
#
# For ts='mle', signal_prior_bkg has no meaning. The nuisance prior for the s+b / b only likelihood is always the same (so far).
#
# what controls which limits are computed. Valid vaues are:
# * 'observed': compute observed limit on data only
# * 'expected': compute +-1sigma, +-2sigma limit bands for background only
# * 'all': both 'data' and 'expected'
#
# reuse_toys should contain a dictionary (spid) -> (list of config filenames) as returned by cls_limits_grid. If set, toys from the grid db files
# in the cache directory will be used.
#
# If setting debug_cls to True, the report will contain all test statistic distributions and "CLs versus beta_signal"-plots
# used to determine the CLs limits. As the name suggest, this is useful mainly for debugging.
#
# returns a tuple of two plotutil.plotdata instances. The first contains expected limit (including the band) and the second the 'observed' limit
# if 'what' is not 'all', one of the plotdata instances is replaced with None.
def cls_limits(model, what = 'all',  cl = 0.95, ts = 'lr', signal_prior = 'flat', nuisance_prior = '', signal_prior_bkg = None,
   signal_processes = None, reuse_toys = {}, truth_max = None, debug_cls = False, **options):
    if signal_processes is None: signal_processes = [[sp] for sp in model.signal_processes]
    assert len(signal_processes) > 0
    if signal_prior_bkg is None:
        if ts == 'lhclike': signal_prior_bkg = 'flat'
        else: signal_prior_bkg = 'fix:0'
    signal_prior = signal_prior_dict(signal_prior)
    signal_prior_bkg = signal_prior_dict(signal_prior_bkg)
    nuisance_prior = nuisance_prior_distribution(model, nuisance_prior)
    main = {'type': 'cls_limits', 'model': '@model', 'producer': '@ts_producer', 'expected_bands' : 0,
        'output_database': sqlite_database(), 'truth_parameter': 'beta_signal', 'minimizer': minimizer(need_error = True),
        'tol_cls': 0.025, 'clb_cutoff': 0.001, 'debuglog': '@debuglog-name'}
    if truth_max is not None: main['truth_max'] = float(truth_max)
    if what in ('expected', 'all'):
        main['expected_bands'] = 2000
    elif what != 'observed': raise RuntimeError, "unknown option what='%s'" % what
    if what in ('observed', 'all'):
        main['data_source'], dummy = data_source_dict(model, 'data')
    ts_producer, main['ts_column'] = ts_producer_dict(ts, signal_prior_bkg)
    toplevel_settings = {'signal_prior': signal_prior, 'main': main, 'ts_producer': ts_producer, 'minimizer': minimizer(need_error = False), 'debuglog-name': 'XXX',
       'model-distribution-signal': {'type': 'flat_distribution', 'beta_signal': {'range': [0.0, float("inf")], 'fix-sample-value': 0.0}}}
    if ts == 'lhclike': del toplevel_settings['signal_prior']
    toplevel_settings.update(get_common_toplevel_settings(**options))
    cachedir = os.path.join(config.workdir, 'cache')
    cfg_names_to_run = []
    for sp in signal_processes:
        model_parameters = model.get_parameters(sp)
        spid = ''.join(sp)
        toplevel_settings['nuisance_prior'] = nuisance_prior.get_cfg(model_parameters)
        toplevel_settings['debuglog-name'] = 'debuglog' + spid + '.txt'
        spid = ''.join(sp)
        if spid in reuse_toys:
            reuse_names = reuse_toys[spid]
            filenames = map(lambda s: os.path.join(cachedir, s) + '.db', reuse_names)
            #TODO: check whether files exist and complain here already!
            main['reuse_toys'] = {'input_database': {'type': 'sqlite_database_in', 'filenames': filenames}}
        elif 'reuse_toys' in main: del main['reuse_toys']
        name = write_cfg(model, sp, 'cls_limits', '', additional_settings = toplevel_settings, **options)
        cfg_names_to_run.append(name)
    if 'run_theta' not in options or options['run_theta']:
        run_theta(cfg_names_to_run, **options)
    else: return None
    
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
        data = sql(sqlfile, 'select "index", "limit", limit_uncertainty from cls_limits order by "index"')
        if what in ('all', 'observed'):
            assert data[0][0] == 0, "cls limit calculation for data failed! See debug.txt."
            observed_results[sp] = data[0][1], data[0][2]
            del data[0] # so rest is expected ...
        if what in ('all', 'expected'):
            # sort by limit:
            limits = sorted([r[1] for r in data])
            limits_noinf = sorted([r[1] for r in data if r[1] != float("inf")])
            limits = limits_noinf
            n = len(limits)
            n_inf = len([x for x in limits if x==float("inf")])
            if n_inf * 1.0 / n >= 0.025: print "WARNING: too many results are infinity: %f%% (acceptable: 2.5%%)" % (n_inf * 100. / n)
            expected_results[sp] = (limits[n/2], (limits[int(0.16*n)], limits[int(0.84*n)]), (limits[int(0.025*n)], limits[int(0.975*n)]))                
        if debug_cls:
            debug_cls_plots(sqlfile, main['ts_column'])
    x_to_sp = get_x_to_sp(spids, **options)
    
    pd_expected, pd_observed = None, None
    if what in ('all', 'expected'):
        pd_expected = plotutil.plotdata()
        pd_expected.color = '#000000'
        pd_expected.as_function = True
        pd_expected.legend = 'expected limit'
        pd_expected.x = sorted(list(x_to_sp.keys()))
        pd_expected.y  = []
        pd_expected.bands = [([], [], '#00ff00'), ([], [], '#00aa00')]
    if what in ('all', 'observed'):
        pd_observed = plotutil.plotdata()
        pd_observed.color = '#000000'
        pd_observed.as_function = True
        pd_observed.x = sorted(list(x_to_sp.keys()))
        pd_observed.y = []
        pd_observed.yerrors = []
        pd_observed.legend = 'observed limit'
    
    for x in sorted(x_to_sp.keys()):
        sp = x_to_sp[x]
        if what in ('all', 'expected'):
            expected, band1, band2 = expected_results[sp]
            pd_expected.y.append(expected)
            pd_expected.bands[1][0].append(band1[0])
            pd_expected.bands[1][1].append(band1[1])
            pd_expected.bands[0][0].append(band2[0])
            pd_expected.bands[0][1].append(band2[1])
        if what in ('all', 'observed'):
            observed, observed_unc = observed_results[sp]
            pd_observed.y.append(observed)
            pd_observed.yerrors.append(observed_unc)
    bayesian.report_limit_band_plot(pd_expected, pd_observed, 'CLs', 'cls')
    return pd_expected, pd_observed


