# -*- coding: utf-8 -*-

# This file contains statistical method connected to frequentist inference, such as p-values via toys, etc.

from theta_interface import *
import bisect


## Returns a new model, leaves the "model" parameter unchanged.
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
        if prior_nuisance['typ'] != 'gauss': raise RuntimeError, "only gaussian nuisance parameters are supported"
        rvobs = 'rvobs_%s' % p
        result.rvobs_distribution.set_distribution(rvobs, 'gauss', mean = p, width = prior_nuisance['width'], range = prior_nuisance['range'])
        result.distribution.set_distribution_parameters(p, width = inf)
        result.data_rvobsvalues[rvobs] = prior_nuisance['mean'] # usually 0.0
    return result


## Make toys to calculate the likelihood ratio test statistic
#
# This is a low-level intreface for caclulating p-values. For higher-level interfaces, see the pvalue_bkgtoys_runs, pvalue, and discovery methods below.
#
# If run_theta is True (default), returns a dictionary (spid) -> (list of negative log-likelihood ratio values).
# Otherwise, dictionary (spid) -> (Run instance) is returned which can be used to e.g. create all configfiles for parallel execution, etc.
#
# note that 'options' is ignored for run_theta = False
#
# seed is the random seed for toy data generation, in case of input = 'toys:...'. The default of None will use a procedure which is unlikely to produce collissions but
#  the result will not be exactly reproduced if eceuting again.
def deltanll(model, input, n, signal_process_groups = None, nuisance_constraint = None, nuisance_prior_toys = None, signal_prior = 'flat', options = None, run_theta = True, seed = None):
    if signal_process_groups is None: signal_process_groups = model.signal_process_groups
    if options is None: options = Options()
    result = {}
    for spid, signal_processes in signal_process_groups.iteritems():
        r = Run(model, signal_processes, signal_prior = signal_prior, input = input, n = n,
             producers = [DeltaNllHypotest(model, signal_processes, nuisance_constraint)], nuisance_prior_toys = nuisance_prior_toys, seed = seed)
        if not run_theta:
            result[spid] = r
        else:
            r.run_theta(options)
            data = r.get_products(['dnll__nll_diff'])
            result[spid] = data['dnll__nll_diff']
    return result


## \brief make 'background-only' toys for p-value determination
#
# Uses a likelihood ratio test statistic with beta_signal=0.0 as null hypothesis
#
# returns a dictionary (spid) -> (list of n_runs "Run" objects). Each Run object contains the config for the
# creation of n toys. The RUn objects can be used to:
# * get the config files, as preparation for running in parallel on a cluster
# * run locally, in sequence
# * run locally, in parallel
# * etc.
#
# Note that increasing n_runs later is easily possible and will not affect the config files / results already created. On the other hand, changing
# n is not advisable as it will lead to different configurations and nothing can be re-used.
#
# Important: the seed is always set explicitly to i_run + seed_min (with i_run = 0..n_run-1)
# You have to be careful if calling this method more than once to use a different seed_min so that no overlapping seeds are used. In general,
# it is advisable to call this method only once per cited p-value, with n_runs set high enough. See the "discovery" method for an example
# of how to do that correctly.
# that correctly
def pvalue_bkgtoys_runs(model, signal_process_groups = None, n_runs = 10, n = 10000, nuisance_constraint = None, nuisance_prior_toys = None, seed_min = 1):
    if signal_process_groups is None: signal_process_groups = model.signal_process_groups
    result = {}
    for i_run in range(n_runs):
        res = deltanll(model, 'toys:0.0', n, signal_process_groups = signal_process_groups, nuisance_constraint = nuisance_constraint,
            nuisance_prior_toys = nuisance_prior_toys, run_theta = False, seed = seed_min + i_run)
        for spid in res:
            if not spid in result: result[spid] = []
            result[spid].append(res[spid])
    return result


## \brief determine the p-value for the dataset in 'input'
#
# Set 'input' to 'data' (and n=1) to get the "observed" p-value.
# For the 'expected' p-values (tossing all nuisance parameters according to their prior) set 'input' to 'toys:1.0' and n to at least 1000 or so.
#
# It will call pvalue_bkgtoys_runs and run it locally, if necessary. Note that you can prevent running the background-only toys locally
# by calling pvalue_bkgtoys_runs with the exact same parameters, run the resulting theta configs whereever you like and copy
# the result back to the cache before calling this method.
#
# returns a dictionary (spid) -> (list p-values)
# each element in the list of p-values is a two-tuple (p0, p_error) with binomial error as calculated in utils.get_p
# Note that the list can have a length of less than n in case the minimization failed in some toys.
def pvalue(model, input, n, signal_process_groups = None, nuisance_constraint = None, nuisance_prior_toys = None, options = None, bkgtoys_n_runs = 10, bkgtoys_n =  10000, bkgtoys_seed_min = 1):
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


## \brief Determine p-value / "N sigma" from tail distribution of background-only test statistic
#
# The number of toys is be increased adaptively: at most maxit iterations are done, each with n backgrond-only toys.
# The procedure is stopped as soon as the accuracy on the Z value is better than Z_error_max.
#
# nuisance_prior_toys_bkg is the Distribution instance used as nuisance prior for the background-only toys.
# nuisance_contraint is the constraint terms applied for calculating the likelihood ratio test statistic.
#
# Returns a four-tuple (median expected significance, lower 1sigma expected, upper 1sigma expected, observed significance)
# each entry in the tuple is itself a two-tuple (Z, Z_error) where the Z_error is the uncertainty on Z from the limited number of background-only toys.
# If use_data was set to False, only the expected Z-values are computed and Z and Z_error are set to None in the return value.
#
# Options overwritten: data_source seed, minimizer strategy to 'fast' for toys
def discovery(model, spid = None, use_data = True, Z_error_max = 0.05, maxit = 100, n = 10000, input_expected = 'toys:1.0',
   nuisance_constraint = None, nuisance_prior_toys_bkg = None, options = None, verbose = True):
    if spid is None: spid = model.signal_process_groups.keys()[0]
    signal_process_groups = {spid : model.signal_process_groups[spid]}
    if options is None: options = Options()
    
    ts_sorted = deltanll(model, signal_process_groups = signal_process_groups, nuisance_constraint = nuisance_constraint, input = input_expected, n = n)[spid]
    ts_sorted.sort()
    expected = (ts_sorted[int(0.5 * len(ts_sorted))], ts_sorted[int(0.16 * len(ts_sorted))], ts_sorted[int(0.84 * len(ts_sorted))])
    del ts_sorted
    
    if use_data: observed = deltanll(model, signal_process_groups = signal_process_groups, nuisance_constraint = nuisance_constraint, input = 'data', n = 1, options = options)[spid][0]
    
    # (median [n, n0], -1sigma [n, n0], +1sigma [n, n0])
    expected_nn0 = ([0,0], [0,0], [0,0])
    # [n, n0] for observed p-value
    observed_nn0 = [0,0]
    observed_significance = None
    options.set('minimizer', 'strategy', 'fast')
    for seed in range(1, maxit + 1):
        run = pvalue_bkgtoys_runs(model, signal_process_groups = signal_process_groups, n_runs = 1, n = n, nuisance_constraint = nuisance_constraint,
            nuisance_prior_toys = nuisance_prior_toys_bkg, seed_min = seed)[spid][0]
        run.run_theta(options)
        ts_bkgonly = run.get_products(['dnll__nll_diff'])['dnll__nll_diff']
        max_Z_error = 0.0
        expected_Z = [[0,0],[0,0],[0,0]]
        Z, Z_error = None, None
        for i in range(3):
            expected_nn0[i][1] += len(ts_bkgonly)
            expected_nn0[i][0] += count(lambda c: c >= expected[i], ts_bkgonly)
            expected_Z[i] = get_Z(*expected_nn0[i])
            max_Z_error = max(max_Z_error, Z_error)
        if use_data:
            observed_nn0[1] += len(ts_bkgonly)
            observed_nn0[0] += count(lambda c: c >= observed, ts_bkgonly)
            Z, Z_error = get_Z(*observed_nn0)
            max_Z_error = max(max_Z_error, Z_error)
            observed_significance = Z, Z_error
        if verbose:
            print "after %d iterations" % seed
            if use_data: print "    observed_significance = %.3f +- %.3f" % (Z, Z_error)
            print "    expected significance (median, lower 1sigma, upper 1sigma): %.3f (%.3f--%.3f)" % (expected_Z[0][0], expected_Z[1][0], expected_Z[2][0])
        if max_Z_error < Z_error_max: break
    return tuple(expected_Z + [(Z, Z_error)])
