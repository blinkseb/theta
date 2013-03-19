# -*- coding: utf-8 -*-
import config, os.path, math
from theta_interface import *
import Report
from utils import *



def bayesian_quantiles(model, input, n, quantiles = [0.95], signal_process_groups = None, nuisance_constraint = None, nuisance_prior_toys = None, signal_prior = 'flat',
   options = None, parameter = 'beta_signal', iterations = 10000, seed = 0):
    """
    Compute bayesian posterior quantiles for a given ensemble
    
    For a documentation of the parameters, see :ref:`common_parameters`. The parameters specific to the method are:
    
    * ``quantiles``: a list a quantiles to compute.
    * ``iterations``: the number of iterations to use in the Markov-Chain Monte-Carlo algorithm
    * ``seed``: the random seed for the MCMC
    
    Returns a dictionary (spid) --> (q) --> (list of results)
    where ``q`` is one element of ``quantiles`` (i.e., a float value). The list of results are the quantiles.
    
    ``q`` can also be the special string "accrate" to return the acceptance rate of the Markov Chain. This is
    a useful diagnotical tool: acceptance rates below 10% or above 35% are usually suspicious.
    """
       
       
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
    
    

def bayesian_nl_posterior_ratio(model, input, n, signal_prior_sb = 'fix:1.0', signal_prior_b = 'fix:0.0', signal_process_groups = None, nuisance_constraint = None,
   nuisance_prior_toys = None, options = None, iterations = 10000):
    """
    Compute the negative logarithm of the Bayesian posterior ratio at different values of beta_signal.

    Returns a dictionary (spid) --> (list of results)
    The list of results are the ratios
    """
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


def bayesian_limits(model, what = 'all', **options):
    """
    Calculate Bayesian limits on the cross section beta_signal.
 
    This is a high-level interface to calculate expected and observed limits and make a limit band plot:
    it will calculate expected and observed limits (using :func:`bayesian_quantiles` ), make a "limit vs. mass" band plot
    and write the result to the global `report` object.

    The ``what`` parameter controls which limits are computed. Valid vaues are:
     * 'observed': compute observed limit on data only
     * 'expected': compute +-1sigma, +-2sigma limit bands for background only
     * 'all': both 'data' and 'expected'
 
    Further ``options`` are passed through to :func:`bayesian_quantiles` with the special case that
     * there should be "n_toy", "n_data"  instead of "n";  "n" will be ignored
     * "input" will be ignored
     
     
    The return value is a two-tuple of ``plotutil.plotdata`` instances. The first contains expected
    limit (including the 1sigma and 2sigma bands), the second the 'observed' limit.
    If ``what`` is different from 'all', the coorresponding tuple entry is `` None``.
     
    For example, to calculate expected limits using 2000 toys, and observed limits using 20 Markov Chains (instead of the default of 1000 / 10) use::
    
      expected, observed = bayesian_limits(model, 'all', n_toy = 2000, n_data = 20)
    """
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




def bayesian_posteriors(model, input, n, histogram_specs, signal_process_groups = None, nuisance_constraint = None, nuisance_prior_toys = None, signal_prior = 'flat', options = None, smooth = True, iterations = 10000):
    """
    Calculate the marginal posterior of the given parameters

    * ``histogram_specs`` is a dictionary of (parameter name) -> tuple(int nbins, float xmin, float max) and determines for which
      parameters the posterior is computed on which range and binning. A typical choice is ``histogram_specs = {'beta_signal': (100, 0.0, 3.0)}``. Note that the minimum and
      maximum range given should include the whole region where the posterior does not vanish.
    * ``smooth``: if ``True``, a whole histogram will be computed at each iterations, instead of just binning the chain elements. This results in a smooth
      posterior histogram, but has a runtime of nbins * iterations instead of just iterations.

    Returns a dictionary (spid) --> (parameter name) -> (list of Histograms)
    """
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



def bayesian_posterior_model_prediction(model, input, n, signal_process_groups = None, nuisance_constraint = None, nuisance_prior_toys = None, signal_prior = 'flat', options = None, iterations = 10000):
    """
    Get the mean and standard deviation of the predicted Poisson mean in each bin while running a Markov-Chain.
    This can be seen as the "posterior model prediction".
    
    While running the Markov Chain, the model prediction in each bin is evaluated. From the complete Markov-Chain, the mean and standard deviation
    in each bin is calculated and returned.
    
    See :ref:`common_parameters` for an explanation of the parameters of this method.
    
    The return value is a nested dictionary:
    
    * The first key is the signal process group id.
    * The second key is the observable name.
    
    The value is a list of ``n`` :class:`Histogram` instances with the values set to the mean and the uncertainties set to the standard deviation.
    
    .. note:: The posterior across different bins is correlated, but this correlation is not reported. There is currently no way to get these correlations. You should therefore
     not use the result for further statistical evaluation; it is only suitable for displaying the result.
     
    .. note:: The model prediction is evaluated at the current parameter values for each Markov Chain element, but without considering the minimization in the
      extra parameters in case of using the Barlow-Beeston light method for handling Monte-Carlo statistical uncertainties.
    """
    if signal_process_groups is None: signal_process_groups = model.signal_process_groups
    if options is None: options = Options()
    result = {}
    observables = model.get_observables()
    colnames = ['mp__%s_mean' % s for s in observables] + ['mp__%s_width' % s for s in observables]
    for spid, signal_processes in signal_process_groups.iteritems():
        r = Run(model, signal_processes, signal_prior = signal_prior, input = input, n = n,
             producers = [MCMCMeanPredictionProducer(model, signal_processes, nuisance_constraint, signal_prior = signal_prior, iterations = iterations)],
             nuisance_prior_toys = nuisance_prior_toys)
        r.run_theta(options)
        res = r.get_products(colnames)
        result[spid] = {}
        for o in observables:
            result[spid][o] = [histogram_from_dbblob(v,u) for v,u in zip(res['mp__%s_mean' % o], res['mp__%s_width' % o])]
    return result
