# -*- coding: utf-8 -*-
import config, utils, os.path, datetime, math, bisect

import scipy.special

import Report
from utils import *


# container for toys made for the CLs construction
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
        
    # return a tuple (cls+b, clb, cls) of plotutil instances which are
    # the CLs+b, clb and cls curves as a function of the "truth" (beta_signal)
    def get_cl_vs_truth(self, truth_to_ts):
        pd_clsb, pd_clb, pd_cls = plotutil.plotdata(as_function = True, legend = 'CLs+b', color = '#0000ff'), plotutil.plotdata(as_function = True, legend = 'CLb'), plotutil.plotdata(as_function = True, color = '#ff0000', legend = 'CLs')
        pd_clsb.x = sorted(list(self.truth_values()))
        pd_clb.x, pd_cls.x = pd_clsb.x, pd_clsb.x
        for t in pd_clb.x:
            ts = truth_to_ts[t]
            n_b = len(self.truth_to_ts_b[t])
            clb = len([x for x in self.truth_to_ts_b[t] if x < ts]) * 1.0 / n_b
            clb_error = max(math.sqrt(clb * (1 - clb) / n_b), 1.0 / n_b)
            pd_clb.y.append(clb)
            n_sb = len(self.truth_to_ts_sb[t])
            clsb = len([x for x in self.truth_to_ts_sb[t] if x < ts]) * 1.0 / n_sb
            clsb_error = max(math.sqrt(clsb * (1 - clsb) / n_sb), 1.0 / n_sb)
            pd_clsb.y.append(clsb)
            if clb==0: clb += clb_error
            cls = clsb / clb
            pd_cls.y.append(cls)
            cls_error = math.sqrt((clsb_error / clb)**2 + (clb_error * cls / clb)**2)
        return pd_clsb, pd_clb, pd_cls
            

# make some plots of the test statistic for the report, using the given db filename
#
# returns a list of four plotutil.plotdata instances which contain curves as a function of the truth value (beta_signal):
# * cls+b for data
# * clb for data
# * cls for data
# * expected CLs curve, including 1sigma and 2sigma bands
def debug_cls_plots(dbfile, ts_column = 'lr__nll_diff'):
    limits = sql(dbfile, 'select "index", "limit", limit_uncertainty from cls_limits')
    indices = [row[0] for row in limits]
    have_data = 0 in indices
    data = sql(dbfile, 'select runid, eventid, lr__poi, source__truth, "%s" from products order by runid, eventid' % ts_column)
    tts = truth_ts_values()
    truth_to_ts_data = {}
    for row in data:
        if row[0] == 0 and row[1] == 0:
            truth_to_ts_data[row[2]] = row[4]
            continue
        if row[3] > 0:  tts.add_point_sb(row[3], row[4])
        else: tts.add_point_b(row[2], row[4])
    plotsdir = os.path.join(config.workdir, 'plots')
    config.report.new_section('debug_cls for file %s' % dbfile)
    pd_clsb, pd_clb, pd_cls = tts.get_cl_vs_truth(truth_to_ts_data)
    
    expected_pd_cls = plotutil.plotdata(as_function = True)
    expected_pd_cls.x = pd_cls.x[:]
    # use median as "expected":
    truth_to_ts = tts.get_truth_to_ts(0.5)
    pd_clsb, pd_clb, pd_cls = tts.get_cl_vs_truth(truth_to_ts)
    expected_pd_cls.y[:] = pd_cls.y
    #build the bands:
    band1s = [[], [], '#00ff00']
    band2s = [[], [], '#ffff00']
    truth_to_ts = tts.get_truth_to_ts(0.025)
    pd_clsb, pd_clb, pd_cls = tts.get_cl_vs_truth(truth_to_ts)
    band2s[0][:] = pd_cls.y[:]
    truth_to_ts = tts.get_truth_to_ts(0.975)
    pd_clsb, pd_clb, pd_cls = tts.get_cl_vs_truth(truth_to_ts)
    band2s[1][:] = pd_cls.y[:]
    truth_to_ts = tts.get_truth_to_ts(0.16)
    pd_clsb, pd_clb, pd_cls = tts.get_cl_vs_truth(truth_to_ts)
    band1s[0][:] = pd_cls.y[:]
    truth_to_ts = tts.get_truth_to_ts(0.84)
    pd_clsb, pd_clb, pd_cls = tts.get_cl_vs_truth(truth_to_ts)
    band1s[1][:] = pd_cls.y[:]        
    expected_pd_cls.bands = [band2s, band1s]    
    
    #draw all:
    pds = [pd_clsb, pd_clb, pd_cls, expected_pd_cls]
    plotname = 'debug_cl-vs-truth-data.png'
    plot(pds, 'beta_signal', 'p-value', os.path.join(plotsdir, plotname))
    config.report.add_html('<p>For data: <br/><img src="plots/%s" /></p>' % plotname)
    
    return pds
        

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
def cls_limits(model, what = 'all',  cl = 0.95, ts = 'lhclike', signal_prior = 'flat', nuisance_prior = '', signal_prior_bkg = None,
   signal_processes = None, reuse_toys = {}, truth_max = None, debug_cls = False, bootstrap_nuisancevalues = False, **options):
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
        'tol_cls': 0.025, 'clb_cutoff': 0.02, 'debuglog': '@debuglog-name', 'rnd_gen': {'seed': options.get('toydata_seed', 1)}}
    if truth_max is not None: main['truth_max'] = float(truth_max)
    if 'reltol_limit' in options: main['reltol_limit'] = float(options['reltol_limit'])
    if bootstrap_nuisancevalues: main['nuisancevalues-for-toys'] = 'datafit'
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
            limits = sorted([r[1] for r in data if r[1] != float("inf")])
            n = len(limits)
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
    report_limit_band_plot(pd_expected, pd_observed, 'CLs', 'cls')
    return pd_expected, pd_observed
