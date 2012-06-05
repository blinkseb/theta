# -*- coding: utf-8 -*-

# This file contains statistical method connected to frequentist inference, such as p-values via toys, etc.


# Returns a new model, leaves the "model" parameter unchanged.
#
# * for each nuisance parameter:
#   - add a real-valued observable "rvobs_<parameter name>"
#   - add a Gaussian distribution for the rvobs with mean = nuisance parameter, parameter = rvobs_...
# * make all nuisance parameter priors flat in the model distribution, s.t.
#   all constraints come from the real-valued observables, not from the priors.
# * add the default value according to model.distribution to the real-valued observables in data
def frequentize_model(model):
    result = copy.deepcopy(model)
    for p in model.distribution.get_parameters():
        prior_nuisance = model.distribution.get_distribution(p)
        # have to use the conjugate distribution here. gauss is self-conjugate, so no problem here in most cases:
        if prior_nuisance['typ'] != 'gauss': raise RuntimeError, "currently, only gaussian nuisance parameters are supported"
        rvobs = 'rvobs_%s' % p
        result.rvobs_distribution.set_distribution(rvobs, 'gauss', mean = p, width = prior_nuisance['width'], range = prior_nuisance['range'])
        result.distribution.set_distribution_parameters(p, width = inf)
        result.data_rvobsvalues[rvobs] = prior_nuisance['mean'] # usually 0.0
    return result


# FIXME: cleanup, from here!

"""
# returns a tuple (theta config dictionary, ts column name) for ts producer specification ts (either 'lr' or 'mle'). Refers
# to the global settings "signal_prior", "nuisance_prior" and "minimizer"
def ts_producer_dict(ts, signal_prior_bkg='flat'):
    assert ts in ('lr', 'lhclike')
    signal_prior_bkg = signal_prior_dict(signal_prior_bkg)
    if ts == 'lr':
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
def ts_data(model, ts, signal_prior = 'flat', nuisance_prior = '', signal_prior_bkg='fix:0', signal_processes = None, **options):
    if signal_processes is None: signal_processes = [[sp] for sp in model.signal_processes]
    signal_prior = signal_prior_dict(signal_prior)
    signal_prior_bkg = signal_prior_dict(signal_prior_bkg)
    nuisance_prior = nuisance_prior_distribution(model, nuisance_prior)
    main = {'model': '@model', 'n-events': 1, 'producers': ('@ts_producer',), 'output_database': sqlite_database(), 'log-report': False}
    ts_producer, ts_colname = ts_producer_dict(ts, signal_prior_bkg)    
    toplevel_settings = {'signal_prior': signal_prior, 'main': main, 'ts_producer': ts_producer, "minimizer": minimizer(need_error = False)}
    toplevel_settings.update(get_common_toplevel_settings(**options))
    main['data_source'], toplevel_settings['model-distribution-signal'] = data_source_dict(model, 'data')
    cfg_names_to_run = []
    for sp in signal_processes:
        model_parameters = model.get_parameters(sp)
        toplevel_settings['nuisance_prior'] = nuisance_prior.get_cfg(model_parameters)
        name = write_cfg(model, sp, 'ts_data', '', additional_settings = toplevel_settings, **options)
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
# bkg_only_toys is a special setting for the lhclike test statistic, as this test statistic in general depends on beta_signal. If bkg_only_toys
# is set to True, the toys are always background only and the beta_signal_values are only used in the evaluation of the test statistic. Otherwise,
# the beta_signal_values are used both, for toy generation, and as parameter for evaluating the test statistic.
#
def ts_toys(model, beta_signal_values, n, ts, signal_prior = 'flat', nuisance_prior = '', signal_prior_bkg = None,
   signal_processes = None, bkg_only_toys = False, ret = 'data', asimov = False, **options):
    if signal_processes is None: signal_processes = [[sp] for sp in model.signal_processes]
    if signal_prior_bkg is None:
        if ts == 'lhclike': signal_prior_bkg = 'flat'
        else: signal_prior_bkg = 'fix:0'
    signal_prior = signal_prior_dict(signal_prior)
    signal_prior_bkg = signal_prior_dict(signal_prior_bkg)
    nuisance_prior = nuisance_prior_distribution(model, nuisance_prior)
    main = {'model': '@model', 'n-events': n, 'producers': ('@ts_producer',), 'output_database': sqlite_database(), 'log-report': False}
    ts_producer, ts_colname = ts_producer_dict(ts, signal_prior_bkg)
    toplevel_settings = {'main': main, 'ts_producer': ts_producer, "minimizer": minimizer(need_error = False)}
    if ts != 'lhclike': toplevel_settings['signal_prior'] = signal_prior
    toplevel_settings.update(get_common_toplevel_settings(**options))
    inp_spec = 'toys:0'
    if asimov: inp_spec = 'toys-asimov:0'
    main['data_source'], toplevel_settings['model-distribution-signal'] = data_source_dict(model, inp_spec, **options)
    #print toplevel_settings
    cfg_names_to_run = []
    for beta_signal_value in beta_signal_values:
        if ts=='lhclike': ts_producer['default_poi_value'] = beta_signal_value
        for sp in signal_processes:
            model_parameters = model.get_parameters(sp)
            toplevel_settings['nuisance_prior'] = nuisance_prior.get_cfg(model_parameters)
            toplevel_settings['model-distribution-signal']['beta_signal'] = beta_signal_value
            if bkg_only_toys:
                toplevel_settings['model-distribution-signal']['beta_signal'] = 0.0
                id = 'b'
            else: id = 'sb'
            name = write_cfg(model, sp, 'ts_toys', '%f' % beta_signal_value, id = id, additional_settings = toplevel_settings, **options)
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
    
## \brief Determine p-value / "N sigma" from tail distribution of background-only test statistic
#
# The number of toys will be increased adaptively.
#
# \param Z_error_max The maximum error in the Z value
#
# \return A tuple (expected significance, observed significance). Each entry in the tuple is in turn a triple (central value, error_minus, error_plus)
# where the error for the expected significance are the expected 1sigma spread and for the observed one the error from limited toy MC.
def discovery(model, ts = 'lr', use_data = True, Z_error_max = 0.05, **options):
    options['signal_prior']='flat'
    options['signal_prior_bkg'] = 'fix:0'
    ts_expected = ts_toys(model, [1.0], 10000, ts = ts, **options)
    # make (spid) -> (median, -1sigma, +1sigma):
    expected = {}
    for spid, beta_signal in ts_expected:
        ts_sorted = sorted(ts_expected[(spid, beta_signal)])
        expected[spid] = (ts_sorted[int(0.5 * len(ts_sorted))], ts_sorted[int(0.16 * len(ts_sorted))], ts_sorted[int(0.84 * len(ts_sorted))])
    del ts_expected
    if use_data: observed = ts_data(model, ts = ts, **options)
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
        ts_bkgonly = ts_toys(model, [0.0], n, ts = ts, **options)
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

"""

