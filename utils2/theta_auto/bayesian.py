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
# where q is one element of quantiles. The list of results are the quantiles 
def bayesian_quantiles(model, input, n, quantiles = [0.95], signal_process_groups = None, nuisance_constraint = None, nuisance_prior_toys = None, signal_prior = 'flat', options = None, parameter = 'beta_signal'):
    if signal_process_groups is None: signal_process_groups = model.signal_process_groups
    if options is None: options = Options()
    colnames = ['quant__quant%05d' % int(q*10000 + 0.5) for q in quantiles]
    result = {}
    for spid, signal_processes in signal_process_groups.iteritems():
        r = Run(model, signal_processes, signal_prior = signal_prior, input = input, n = n,
             producers = [QuantilesProducer(model, signal_processes, nuisance_constraint, signal_prior, parameter = parameter, quantiles = quantiles)],
             nuisance_prior_toys = nuisance_prior_toys)
        r.run_theta(options)
        res = r.get_products(colnames)
        result[spid] = {}
        for i, q in enumerate(quantiles): result[spid][q] = res[colnames[i]]
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





"""

def posteriors(model, histogram_specs, input = 'data', n = 3, signal_prior = 'flat', nuisance_prior = '', signal_processes = None, mcmc_iterations = 10000, **options):
    if signal_processes is None: signal_processes = [[sp] for sp in model.signal_processes]
    signal_prior = signal_prior_dict(signal_prior)
    nuisance_prior = nuisance_prior_distribution(model, nuisance_prior)
    main = {'n-events': n, 'model': '@model', 'producers': ('@posteriors',), 'output_database': sqlite_database(), 'log-report': False}
    posteriors = {'type': 'mcmc_posterior_histo', 'name': 'posteriors', 'parameters': [],
       'override-parameter-distribution': product_distribution("@signal_prior", "@nuisance_prior"),
       'smooth': options.get('smooth', True), 'iterations': mcmc_iterations }
    for par in histogram_specs:
        posteriors['parameters'].append(par)
        nbins, xmin, xmax = histogram_specs[par]
        posteriors['histo_%s' % par] = {'range': [float(xmin), float(xmax)], 'nbins': nbins}
    toplevel_settings = {'signal_prior': signal_prior, 'posteriors': posteriors, 'main': main}
    options['load_root_plugins'] = False
    toplevel_settings.update(get_common_toplevel_settings(**options))
    cfg_names_to_run = []
    main['data_source'], toplevel_settings['model-distribution-signal'] = data_source_dict(model, input)
    cfg_names_to_run = []
    for sp in signal_processes:
        model_parameters = model.get_parameters(sp)
        toplevel_settings['nuisance_prior'] = nuisance_prior.get_cfg(model_parameters)
        name = write_cfg(model, sp, 'posteriors', input, additional_settings = toplevel_settings, **options)
        cfg_names_to_run.append(name)
    if 'run_theta' not in options or options['run_theta']:
        run_theta(cfg_names_to_run)
    else: return None
    
    cachedir = os.path.join(config.workdir, 'cache')
    plotsdir = os.path.join(config.workdir, 'plots')
    
    result = {}
    config.report.new_section('Posteriors %s' % input)
    result_table = Report.table()
    result_table.add_column('process', 'signal process')
    parameters = sorted([par for par in histogram_specs])
    for par in parameters:
        result_table.add_column('maximum posterior %s' % par)
    for name in cfg_names_to_run:
        method, sp, dummy = name.split('-',2)
        config.report.add_html('<h2>For signal "%s"</h2>' % sp)
        result_table.set_column('process', sp)
        sqlfile = os.path.join(cachedir, '%s.db' % name)
        cols = ['posteriors__posterior_%s' % par for par in parameters]
        data = sql(sqlfile, 'select %s from products' % ', '.join(cols))
        data = [map(plotdata_from_histoColumn, row) for row in data]
        result[sp] = {}
        i = 0
        for par in parameters:
            result[sp][par] = [row[i] for row in data]
            for pd in result[sp][par]: pd.as_function = True
            plot(result[sp][par], par, 'posterior density', os.path.join(plotsdir, '%s-%s.png' % (name, par)))
            config.report.add_html('<p>%s:<br/><img src="plots/%s-%s.png" /></p>' % (par, name, par))
            i += 1
            maxima = sorted(map(argmax, result[sp][par]))
            result_table.set_column('maximum posterior %s' % par, '%.3g' % maxima[int(0.5 * len(maxima))])
        result_table.add_row()
    config.report.add_p('input: %s, n: %d, signal prior: %s, nuisance prior: %s' % (input, n, str(signal_prior), str(nuisance_prior)))
    config.report.add_html(result_table.html())
    return result
"""
