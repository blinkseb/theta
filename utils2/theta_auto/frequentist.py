# -*- coding: utf-8 -*-

# This file contains statistical method connected to frequentist inference, such as p-values via toys, etc.

from theta_interface import *


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


# Make toys to calculate the likelihood ratio test statistic
#
# returns a dictionary (spid) -> (list of negative log-likelihood ratio values)
def deltanll(model, input, n, signal_process_groups = None, nuisance_constraint = None, nuisance_prior_toys = None, signal_prior = 'flat', options = None):
    if signal_process_groups is None: signal_process_groups = model.signal_process_groups
    if options is None: options = Options()
    result = {}
    for spid, signal_processes in signal_process_groups.iteritems():
        r = Run(model, signal_processes, signal_prior = signal_prior, input = input, n = n,
             producers = [DeltaNllHypotest(model, signal_processes, nuisance_constraint)], nuisance_prior_toys = nuisance_prior_toys)
        r.run_theta(options)
        res = r.get_products()
        data = r.get_products(['dnll__nll_diff'])
        result[spid] = data['dnll__nll_diff']
    return result
        
## \brief Determine p-value / "N sigma" from tail distribution of background-only test statistic
#
# The number of toys will be increased adaptively: at most maxit iterations are done, each with n backgrond-only toys.
# The procedure is stopped as soon as the accuracy on the Z value is better than Z_error_max.
#
# If use_data is False, only the expected significance is calculated, not the observed one.
#
# nuisance_prior_toys_bkg is the Distribution instance used as nuisance prior for the background-only toys.
# nuisance_contraint is the constraint terms applied for calculating the likelihood ratio test statistic.
#
# Returns a four-tuple (median expected significance, lower 1sigma expected, upper 1sigma expected, observed significance)
# each entry in the tuple is itself a two-tuple (Z, Z_error) where the Z_error is the uncertainty on Z from the limited number of background-only toys.
# If use_data was set to False, Z and Z_error are None.
#
# Options overwritten: data_source seed, minimizer strategy to 'fast' for toys
def discovery(model, spid = None, use_data = True, Z_error_max = 0.05, maxit = 100, n = 10000, input_expected = 'toys:1.0', nuisance_constraint = None, nuisance_prior_toys_bkg = None, options = None):
    if spid is None: spid = model.signal_process_groups.keys()[0]
    signal_process_groups = {spid : model.signal_process_groups[spid]}
    if options is None: options = Options()
    
    ts_sorted = sorted(deltanll(model, signal_process_groups = signal_process_groups, nuisance_constraint = nuisance_constraint, input = input_expected, n = n)[spid])
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
        options.set('data_source', 'seed', "%d" % seed)
        ts_bkgonly = deltanll(model, signal_process_groups = signal_process_groups, input = 'toys:0.0', n = n, options = options, nuisance_constraint = nuisance_constraint,
            nuisance_prior_toys = nuisance_prior_toys_bkg)[spid]
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
            print "observed significance for process '%s' (with stat. error from limited toys): %f +- %f" % (spid, Z, Z_error)
        #print "max Z error: ", max_Z_error
        if max_Z_error < Z_error_max: break
    return tuple(expected_Z + [(Z, Z_error)])

