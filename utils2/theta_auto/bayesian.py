# -*- coding: utf-8 -*-
import config, os.path, math
from theta_interface import *
import Report
from utils import *


## \brief Compute bayesian posterior quantiles for a given ensemble
#
# Uses the theta plugin mcmc_quantiles.
#
# 
# Report: Writes a table containing the quantiles for each signal process to the report. For data, it will contain
# the best estimate and uncertainty of the determined limit (the uncertainty is only available for n > 1). For input!='data',
# the +-1sigma and +-2sigma (i.e., the centrgl 68% and central 84%) are given.
#
# Returns a dictionary (spid) --> (q) --> (list of results)
# where q is one element of quantiles (a float value). The list of results are the quantiles 
#
# q can also be the special string "accrate" to return the acceptance rate, if the parameter \c accrate was set to
# \c True
def bayesian_quantiles(model, input, n, quantiles = [0.95], signal_process_groups = None, nuisance_constraint = None, nuisance_prior_toys = None, signal_prior = 'flat',
   options = None, parameter = 'beta_signal', iterations = 10000, seed = 0):
    if signal_process_groups is None: signal_process_groups = model.signal_process_groups
    if options is None: options = Options()
    colnames = ['quant__quant%05d' % int(q*10000 + 0.5) for q in quantiles] + ['quant__accrate']
    result = {}
    for spid, signal_processes in signal_process_groups.iteritems():
        p = QuantilesProducer(model, signal_processes, nuisance_constraint, signal_prior, parameter = parameter, quantiles = quantiles, iterations = iterations, seed = seed)
        r = Run(model, signal_processes, signal_prior = signal_prior, input = input, n = n,
             producers = [p], nuisance_prior_toys = nuisance_prior_toys)
        r.run_theta(options)
        res = r.get_products(colnames)
        result[spid] = {}
        for i, q in enumerate(quantiles): result[spid][q] = res[colnames[i]]
        result[spid]['accrate'] = res['quant__accrate']
    return result
    
    
## \brief Compute bayesian posterior ratio
#
# Uses the theta plugin mcmc_posterior_ratio.
#
# beta_signal_values are the values for beta_signal to use in the nominator / denominator of the ratio.
#
# Returns a dictionary (spid) --> (list of results)
# The list of results are the ratios
def bayesian_nl_posterior_ratio(model, input, n, signal_prior_sb = 'fix:1.0', signal_prior_b = 'fix:0.0', signal_process_groups = None, nuisance_constraint = None,
   nuisance_prior_toys = None, options = None, iterations = 10000):
    if signal_process_groups is None: signal_process_groups = model.signal_process_groups
    if options is None: options = Options()
    result = {}
    for spid, signal_processes in signal_process_groups.iteritems():
        r = Run(model, signal_processes, signal_prior = 'flat', input = input, n = n,
             producers = [MCMCRatioProducer(model, signal_processes, nuisance_constraint, signal_prior_sb = signal_prior_sb, signal_prior_b = signal_prior_b, iterations = iterations)],
             nuisance_prior_toys = nuisance_prior_toys)
        r.run_theta(options)
        res = r.get_products(['mcmcratio__nl_posterior_sb', 'mcmcratio__nl_posterior_b'])
        result[spid] = map(lambda r: r[1] - r[0], zip(res['mcmcratio__nl_posterior_sb'], res['mcmcratio__nl_posterior_b']))
    return result

## \brief Calculate Bayesian limits.
#
# This is a high-level interface to calculate expected and observed limits and make a limit band plot:
# it will calculate expected and observed limits (using \ref bayesian_quantiles), make a "limit vs. mass" band plot
# (using \ref limit_band_plot) and write the result to the report (using \ref report_limit_band_plot).
#
# The 'what' parameter controls which limits are computed. Valid vaues are:
# * 'observed': compute observed limit on data only
# * 'expected': compute +-1sigma, +-2sigma limit bands for background only
# * 'all': both 'data' and 'expected'
#
# options are as for bayesian_quantiles, with the modification that
# * there should be n_toy, n_obs  instead of n;  'n' will be ignored
# * 'input' will be ignored
#
# returns a tuple of two plotutil.plotdata instances. The first contains expected limit (including the band) and the second the 'observed' limit.
# If 'what' is not 'all', one of the tuple entries is None.
def bayesian_limits(model, what = 'all', **options):
    if 'n' in options: del options['n']
    if 'input' in options: del options['input']
    n = options.get('n_toy', 1000)
    plot_expected, plot_observed = None, None
    if what in ('expected', 'all'):
        expected_limits = bayesian_quantiles(model, input = 'toys:0', n = n, **options)
        plot_expected = limit_band_plot(expected_limits, True)
    if what in ('observed', 'all'):
        assert model.has_data()
        n = options.get('n_data', 10)
        observed_limits = bayesian_quantiles(model, input = 'data', n = n, **options)
        plot_observed = limit_band_plot(observed_limits, False)
    # catch the case where the routines return None (e.g., if run_theta = False)
    report_limit_band_plot(plot_expected, plot_observed, 'Bayesian', 'bayesian')
    return (plot_expected, plot_observed)



## \brief Calculate the marginal posterior of the given parameters
#
# histogram_specs is a dictionary of (parameter name) -> tuple(int nbins, float xmin, float max) and determines for which
# parameters the posterior is computed on which range and binning. Note that the computational complexity is roughly
# proportional to nbins * mcmc_iterations. Therefore, use large values only if you have to / for the final result. It is suggested
# to use 30 bins as a start and mcmc_iterations = 10000.
#
# returns a dictionary (spid) --> (parameter name) -> (list of Histogram)
def bayesian_posteriors(model, input, n, histogram_specs, signal_process_groups = None, nuisance_constraint = None, nuisance_prior_toys = None, signal_prior = 'flat', options = None, smooth = True, iterations = 10000):
    if signal_process_groups is None: signal_process_groups = model.signal_process_groups
    if options is None: options = Options()
    result = {}
    parameters = sorted(histogram_specs.keys())
    colnames = ['post__posterior_%s' % p for p in parameters]
    for spid, signal_processes in signal_process_groups.iteritems():
        r = Run(model, signal_processes, signal_prior = signal_prior, input = input, n = n,
             producers = [PosteriorProducer(model, signal_processes, nuisance_constraint, signal_prior = signal_prior, histogram_specs = histogram_specs, smooth = smooth,
                         iterations = iterations)],
             nuisance_prior_toys = nuisance_prior_toys)
        r.run_theta(options)
        res = r.get_products(colnames)
        result[spid] = {}
        for i, p in enumerate(parameters): result[spid][p] = map(histogram_from_dbblob, res[colnames[i]])
    return result

