# -*- coding: utf-8 -*-

# This file contains statistical method connected to frequentist inference, such as p-values via toys, etc.

from theta_interface import *
from likelihood import mle
import bisect


## Returns a new model, leaves the "model" parameter unchanged.
#
# * for each nuisance parameter with a finite width:
#   - add a real-valued observable "rvobs_<parameter name>"
#   - add a Gaussian distribution for the rvobs with mean = nuisance parameter, parameter = rvobs_...
#   - make the nuisance parameter prior flat in the model distribution, s.t.
#     all constraints come from the real-valued observables, not from the priors.
#   - add the default value according to model.distribution to the real-valued observables in data
#
# Note that frequentize_model is idempotent.
def frequentize_model(model):
    result = model.copy()
    for p in model.distribution.get_parameters():
        prior_nuisance = model.distribution.get_distribution(p)
        # have to use the conjugate distribution here. gauss is self-conjugate, so no problem here in most cases:
        if prior_nuisance['typ'] != 'gauss': raise RuntimeError, "only gaussian nuisance parameters are supported"
        if prior_nuisance['width'] == inf: continue
        rvobs = 'rvobs_%s' % p
        result.rvobs_distribution.set_distribution(rvobs, 'gauss', mean = p, width = prior_nuisance['width'], range = prior_nuisance['range'])
        result.distribution.set_distribution_parameters(p, width = inf)
        result.data_rvobsvalues[rvobs] = prior_nuisance['mean'] # usually 0.0
    return result
    

def get_bootstrapped_model(model):
    """
    Return a new :class:`theta_auto.Model` instance in which the data values for the global
    observables are set to the best fit values for data.
    """
    model_freq = frequentize_model(model)
    pars = model_freq.get_parameters([])
    res = mle(model_freq, 'data', 1, signal_process_groups = {'': []})
    par_values = {}
    for p in pars:
        par_values[p] = res[''][p][0][0]
    for p in pars:
        model_freq.distribution.set_distribution_parameters(p, mean = par_values[p])
    return model_freq


    
def make_data(model, input, n, signal_process_groups = None, nuisance_prior_toys = None, options = None, seed = None):
    """
    Make toys data and save it to a file. This is useful to run different statistical models on the same data.
    
    Returns a dictionary with the signal process group id as key and the
    path to the .db file in which the toy data has been saved as value. This path
    can be used as ``input`` argument to many methods.
    """
    if signal_process_groups is None: signal_process_groups = model.signal_process_groups
    if options is None: options = Options()
    result = {}
    for spid, signal_processes in signal_process_groups.iteritems():
        r = Run(model, signal_processes, signal_prior = 'flat', input = input, n = n,
             producers = [PDWriter()], nuisance_prior_toys = nuisance_prior_toys, seed = seed)
        r.run_theta(options)
        result[spid] = r.get_db_fname()
    return result



def deltanll(model, input, n, signal_process_groups = None, nuisance_constraint = None, nuisance_prior_toys = None, signal_prior = 'flat', options = None,
  run_theta = True, seed = None, lhclike_beta_signal = None, pchi2 = False):
    """
    Calculate the delta-log-likelihood test statistic suitable for signal search for ``input``. The test statistic is
    
    .. math::
    
      q_0 = \\log \\frac{L_{s+b}}{L_b}
      
    where both nominator and denominator are maximized: for L_s+b, all parameters are free, including the signal strength parameter ``beta_signal``. For L_b,
    ``beta_signal`` is set to 0.0. In both cases, the likelihood contains terms according to ``nuisance_contraint``.
    
    For the parameters ``model``, ``input``, ``n``, ``signal_process_groups``, ``nuisance_constraint``, ``nuisance_prior_toys``, ``signal_prior`` and ``options`` refer to :ref:`common_parameters`.
    
    More parameters:
    
    * ``seed`` - this is the random seed to use for toy data generation. It is only relevant for ``input="toys..."``. The default value of ``None`` will use a seed which is
      different in each :program:`theta` run. While this is usually a good idea, it makes the result not exactly reproducible.
    * ``run_theta`` - if ``True``, runs :program:`theta` locally. Otherwise, :class:`Run` objects are returned which can be used e.g. to access the cfg file. Note that in case of
      ``run_theta=False``, the ``options`` parameter has no effect whatsoever.
    * ``lhclike_beta_signal`` - if not ``None``, it should be a floating point value for the LHC-like test statistic evalaution; it is the value of beta_signal tested; the
      beta_signal parameter is fixed to this value in the likelihood ratio calculation. If ``None``, the restrictions commonly used for a signal search are used:
      beta_signal=0 for the background-only  (null hypoythesis) case and a flat prior for beta_signal > 0.0 for the signal+background (alternative hypothesis) case.
    * ``pchi2`` - if ``True``, also calculates and returns the pseudo-chi2 value for the s+b fit for each toy
      
    The return value is a nested python dictionary. The first-level key is the signal process group id (see :ref:`what_is_signal`). The value depends on ``run_theta``:
    
    * in case ``run_theta`` is ``True``, the value is a dictionary with the keys "dnll" containing a list of delta-log-likelihood values. If
      ``pchi2`` is true, it also contains a key ``"pchi2"`` with the chi2 values.
    * if ``run_theta`` is ``False``, the value is an instance of the :class:`theta_interface.Run` class
    """
    if signal_process_groups is None: signal_process_groups = model.signal_process_groups
    if options is None: options = Options()
    result = {}
    for spid, signal_processes in signal_process_groups.iteritems():
        if lhclike_beta_signal is not None:
            pdnll = DeltaNllHypotest(model, signal_processes, nuisance_constraint, restrict_poi = 'beta_signal',
                           restrict_poi_value = lhclike_beta_signal, signal_prior_sb = 'flat', signal_prior_b = 'flat')
        else:
            pdnll = DeltaNllHypotest(model, signal_processes, nuisance_constraint, write_pchi2 = pchi2)
        r = Run(model, signal_processes, signal_prior = signal_prior, input = input, n = n,
             producers = [pdnll], nuisance_prior_toys = nuisance_prior_toys, seed = seed)
        if not run_theta:
            result[spid] = r
        else:
            r.run_theta(options)
            cols = ['dnll__nll_diff']
            if pchi2: cols += ['dnll__pchi2']
            data = r.get_products(cols)
            result[spid] = {'dnll': data['dnll__nll_diff']}
            if pchi2: result[spid]['pchi2'] = data['dnll__pchi2']
    return result


def pvalue_bkgtoys_runs(model, signal_process_groups = None, n_runs = 10, n = 10000, nuisance_constraint = None, nuisance_prior_toys = None, seed_min = 1, pchi2 = False):
    """
    Prepare Run instances for 'background-only' toys for p-value determination with ``pvalue``.
    
    For the parameters ``model``, ``signal_process_groups``, ``nuisance_constraint``, ``nuisance_prior_toys`` refer to :ref:`common_parameters`.
    
    Parameters:
    
    * ``n_runs`` is the number of ``Run`` instances to return
    * ``n`` is the number of toys per ``Run`` instance; so the total number of toys in the returned configuration will be ``n_runs * n``
    * ``seed_min`` is the minimum random seed to use. It is incremented by 1 for each returned ``Run``, so for ``seed_min = 1``, the random seeds used will be ``1..n_runs``.
    * ``pchi2`` - if ``True``, configure :program:`theta` to calculate the pseudo-chi2 value for the s+b fit for each toy
    
    Returns a dictionary with the signal process group id as key. The value is a list of :class:`theta_interface.Run` instances, which can be used
    for :ref:`distributed_running`, or to execute :program:`theta` locally with more control, e.g. running multiple :program:`theta` processes
    in parallel. See the documentation of :class:`theta_interface.Run` for more information.
    
    For distributed running, you have to:
    
    1. call ``pvalue_bkgtoys_runs``, get the config files, run :program:`theta` on all of them (maybe in a cluster) and copy the .cfg and the .db files to the "cache" directory (which is a subdirectory of the analysis workdir), see :ref:`distributed_running` for details.
    2. call :meth:`pvalue`, using the same ``bkgtoys_*`` parameters as in step 1., as only this ensures that the same .cfg files are created and the result from the cache created in step 1. will be used
   
.. note::

   If calling the method again with the same parameters and increased ``n_runs``, the created config files will be identical for the first previously created
   ones. This allows to increase the number of background-only toys without loosing the first ones. Note that this is *not* true for changing ``n``.
    
.. important::
    
    The random seed used for toys is always set explicitly to ``i_run + seed_min`` (with ``i_run = 0..n_run-1`` ).
    You have to be careful if calling this method more than once to use a different ``seed_min`` so that no overlapping seeds are used. In general,
    it is advisable to call this method only once per cited p-value, with ``n_runs`` set to a high enough value. You can look at the implementation of their
    ``discovery`` method for an example of how to do that correctly.
    """
    if signal_process_groups is None: signal_process_groups = model.signal_process_groups
    result = {}
    for i_run in range(n_runs):
        res = deltanll(model, 'toys:0.0', n, signal_process_groups = signal_process_groups, nuisance_constraint = nuisance_constraint,
            nuisance_prior_toys = nuisance_prior_toys, run_theta = False, seed = seed_min + i_run, pchi2 = pchi2)
        for spid in res:
            if not spid in result: result[spid] = []
            result[spid].append(res[spid])
    return result


def pvalue(model, input, n, signal_process_groups = None, nuisance_constraint = None, nuisance_prior_toys = None, options = None, bkgtoys_n_runs = 10, bkgtoys_n =  10000, bkgtoys_seed_min = 1):
    """
    Determine the p-value(s) for the dataset in 'input'
    
    For the parameters ``model``, ``input``, ``n``, ``signal_process_groups``, ``nuisance_constraint``, ``nuisance_prior_toys`` and ``options`` refer to :ref:`common_parameters`.
    Specifically, to get the "observed" p-value, use ``input="data"`` and ``n=1``.
    
    The remaining parameters ``bkgtoys_n_runs``, ``bkgtoys_n``, ``bkgtoys_seed_min`` are passed to :meth:`pvalue_bkgtoys_runs`, see the documentation
    there.
    
    Note that :program:`theta` will be executed locally. Use :meth:`pvalue_bkgtoys_runs` directly if you want to run :prorgam:`theta` distributed on a cluster.
    
    The return value is a dictionary where the key is the signal process group id. The value is a list of two-tuples ``(p, p_error)``.
    You can use See :func:`theta_auto.p_to_Z` to convert p-values to Z-values.
    
    For example, To get an ensemble of 1000 "expected" p-values for a signal strength ``beta_signal`` of 1.5 use::
    
      pvalue(model, input = 'toys:1.5', n = 1000)
      
    You can then also calculate the p-value for data, re-using the bvackground-only toys just generated with::
    
      pvalue(model, input = 'data', n = 1)
    """
    if options is None: options = Options()
    input_deltanlls = deltanll(model, input, n, signal_process_groups = signal_process_groups, nuisance_constraint = nuisance_constraint,
            nuisance_prior_toys = nuisance_prior_toys, options = options)
    bkg_runs = pvalue_bkgtoys_runs(model, signal_process_groups, n_runs = bkgtoys_n_runs, n = bkgtoys_n, nuisance_constraint = nuisance_constraint,
        nuisance_prior_toys = nuisance_prior_toys, seed_min = bkgtoys_seed_min)
    result = {}
    for spid in bkg_runs:
        result[spid] = []
        bkg_deltanlls = []
        for run in bkg_runs[spid]:
            run.run_theta(options)
            bkg_deltanlls += run.get_products(['dnll__nll_diff'])['dnll__nll_diff']
        bkg_deltanlls.sort()
        for dnll in input_deltanlls[spid]:
            # count how many background-only toys have a TS value >= dnll:
            n0 = len(bkg_deltanlls)
            n_above = n0 - bisect.bisect_left(bkg_deltanlls, dnll)
            result[spid].append(get_p(n_above, n0))
    return result



# debug_method is called using the test statistic values as argument in this order
#   debug_method('expected', ts_list)    -- after the toys for the "expected" significance using input_expected
#   debug_method('data', ts_single)   -- the ts value for data
#   debug_method('bkg', ts_list)    -- the ts value for background-only toys (using input="toys:0.0").
def discovery(model, spid = None, use_data = True, Z_error_max = 0.05, maxit = 100, n = 10000, input_expected = 'toys:1.0', n_expected = 1000,
   nuisance_constraint = None, nuisance_prior_toys_bkg = None, options = None, verbose = True, debug_method = None, chi2_trunc = 0.0):
    """
    Determine p-value / "N sigma" from tail distribution of background-only test statistic.

    The number of toys is be increased adaptively: at most ``maxit`` iterations are done, each with ``n`` backgrond-only toys.
    The procedure is stopped as soon as the (absolute) accuracy on all reported Z values is better than ``Z_error_max``.
    
    For ``nuisance_constraint`` and ``options``, refer to :ref:`common_parameters`. 
    
    Parameters:
    
    * ``spid`` - the signal process group id
    * ``use_data`` - if ``True``, also calculate the Z value for data
    * ``Z_error_max``, ``maxit`` define the stopping ctriteria, see above.
    * ``n`` number of background-only toys per iterations (there are ``maxit`` iterations maximum)
    * ``input_expected`` - a ``input``-like string which deinfes what is reported as "expected" Z-value
    * ``nuisance_prior_toys_bkg`` is like ``nuisance_prior_toys`` (see :ref:`common_parameters`), but only applied to the "background-only" toys.
    
    Returns a four-tuple (median expected significance, lower 1sigma expected, upper 1sigma expected, observed)
    each entry in the tuple is itself a two-tuple ``(Z, Z_error)`` where the Z_error is the uncertainty on ``Z`` from the limited number of background-only toys.
    
    In case ``use_data = False``, only the expected Z-values are computed and ``Z`` and ``Z_error`` for the observed Z-value in the return value are both set to ``None``.
    """
    if spid is None: spid = model.signal_process_groups.keys()[0]
    signal_process_groups = {spid : model.signal_process_groups[spid]}
    if options is None: options = Options()
    
    pchi2 = chi2_trunc > 0.0
    res = deltanll(model, signal_process_groups = signal_process_groups, nuisance_constraint = nuisance_constraint,
                   input = input_expected, n = n_expected, pchi2 = pchi2, options = options)[spid]
    ts_sorted = res['dnll']
    ts_sorted.sort()
    if chi2_trunc > 0.0:
        chi2s = res['pchi2']
        chi2s.sort()
        index = -int(len(chi2s) * chi2_trunc)
        if index == 0: chi2_threshols = float("inf")
        else: chi2_threshold = chi2s[index]
    if debug_method is not None: debug_method('expected', ts_sorted)
    expected = (ts_sorted[int(0.5 * len(ts_sorted))], ts_sorted[int(0.16 * len(ts_sorted))], ts_sorted[int(0.84 * len(ts_sorted))])
    del ts_sorted
    
    if use_data:
        res = deltanll(model, signal_process_groups = signal_process_groups, nuisance_constraint = nuisance_constraint, input = 'data', n = 1, pchi2 = pchi2, options = options)[spid]
        observed = res['dnll'][0]
        if chi2_trunc > 0.0 and res['pchi2'][0] > chi2_threshold:
            raise RuntimeError, "chi2 value on data too bad. No p-value on data available."
        if debug_method is not None: debug_method('data', observed)
        
    
    # (median [n, n0], -1sigma [n, n0], +1sigma [n, n0])
    expected_nn0 = ([0,0], [0,0], [0,0])
    # [n, n0] for observed p-value
    observed_nn0 = [0,0]
    observed_significance = None
    options.set('minimizer', 'strategy', 'fast')
    if verbose: print "making at most maxit=%d iterations of background-only toys, each with n=%d toys" % (maxit, n)
    for seed in range(1, maxit + 1):
        # only create only one run, so seed is never re-used.
        run = pvalue_bkgtoys_runs(model, signal_process_groups = signal_process_groups, n_runs = 1, n = n, nuisance_constraint = nuisance_constraint,
            nuisance_prior_toys = nuisance_prior_toys_bkg, seed_min = seed, pchi2 = pchi2)[spid][0]
        run.run_theta(options)
        if pchi2:
            res = run.get_products(['dnll__nll_diff', 'dnll__pchi2'])
            ts_bkgonly = [dnll for dnll, pchi2 in zip(res['dnll__nll_diff'], res['dnll__pchi2']) if pchi2 < chi2_threshold]
            if verbose: print "keeping %d / %d background-only toys according to chi2 threshold." % (len(ts_bkgonly), len(res['dnll__pchi2']))
        else:
            res = run.get_products(['dnll__nll_diff'])
            ts_bkgonly = res['dnll__nll_diff']
        if debug_method is not None: debug_method('bkg', ts_bkgonly)
        max_Z_error = 0.0
        expected_Z = [[0,0],[0,0],[0,0]]
        Z, Z_error = None, None
        for i in range(3):
            expected_nn0[i][1] += len(ts_bkgonly)
            expected_nn0[i][0] += count(lambda c: c >= expected[i], ts_bkgonly)
            expected_Z[i] = get_Z(*expected_nn0[i])
            max_Z_error = max(max_Z_error, expected_Z[i][1])
        exp_info = ''
        if expected_nn0[0][0] == 0:
            exp_info = ' (>~ %.3f)'  % get_Z(1, expected_nn0[0][1])[0]
        if use_data:
            observed_nn0[1] += len(ts_bkgonly)
            observed_nn0[0] += count(lambda c: c >= observed, ts_bkgonly)
            Z, Z_error = get_Z(*observed_nn0)
            obs_info = '' if observed_nn0[0]==0 else " (>~ %.3f)" % get_Z(1, observed_nn0[1])[0]
            max_Z_error = max(max_Z_error, Z_error)
            observed_significance = Z, Z_error
        if verbose:
            print "after %d iterations" % seed
            if use_data: print "    observed_significance = %.3f +- %.3f%s" % (Z, Z_error, obs_info)
            print "    expected significance (median, lower 1sigma, upper 1sigma): %.3f +-%.3f%s (%.3f--%.3f)" % (expected_Z[0][0], expected_Z[0][1], exp_info, expected_Z[1][0], expected_Z[2][0])
        if max_Z_error < Z_error_max:
            print "current max error on Z is %.3f, which is smaller than the provided threshold Z_error_max=%.3f; stopping iteration." % (max_Z_error, Z_error_max)
            break
    return tuple(expected_Z + [(Z, Z_error)])
