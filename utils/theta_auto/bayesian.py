import config, os.path, math
#from Distribution import Distribution
from theta_interface import *
import Report

from utils import *


# cfg naming convention: <method>-<signal name>-<input>[-<input seed>]
# where <input> is either "data" for data, "toys:X" for toys with beta_signal = X,
# <input seed> is the random seed for the data_source in case of input "toys:X".
#
# signal_processes is a list of lists: the process groups to consider as (simultaneous!) signal. If none, each
# process in model.signal_processes is considered as (only) signal process ...
#
# options:
#   * 'toydata_seed' [default: 1] for another random seed used in toy data generation 
#   * 'run_theta' [default: True] whether to really run theta. Otherwise, only the cfg files are written. In case run_theta is False, None is returned
#   * 'mcmc_iterations' [default: 20000]
#   * 'mcmc_nproducers': number of independent mcmc producers to run per toy; used to derive uncertainties on the quantiles [default: 1, no uncertainty]
#   * 'mcmc_seed': the random sed to use for the mcmc_producer. The actual seed used for producer i is mcmc_seed + i
#      with i from 0 to mcmc_nproducers-1.
#
# returns:
#  a dictionary with signal process name as key. Value is an array with outcomes in the toy. For example
#  result['zp1000'] is the array of the quantiles for the signal of name 'zp1000'.
#  In case mcmc_nproducers is >= 2, there will also be a key (signal process)__uncertainty', e.g.
#  result['zp1000__uncertainty'] which contains the uncertainty for each toy
def bayesian_quantiles(model, input = 'toys:0', n = 1000, signal_prior = 'flat', nuisance_prior = '', quantile = 0.95, write_report = True, signal_processes = None, **options):
    if signal_processes is None: signal_processes = [[sp] for sp in model.signal_processes]
    mcmc_nproducers = 1
    if 'mcmc_nproducers' in options: mcmc_nproducers = options['mcmc_nproducers']
    # allow n==0 to simplify some scripts ...
    if n==0: return {}
    signal_prior = signal_prior_dict(signal_prior)
    nuisance_prior = nuisance_prior_distribution(model, nuisance_prior)
    main = {'n-events': n, 'model': '@model', 'producers': ['@bayes_interval%d' % d for d in range(mcmc_nproducers)], 'output_database': sqlite_database(), 'log-report': False}
    bayes_interval = {'type': 'mcmc_quantiles', 'name': 'bayes', 'parameter': 'beta_signal',
       'override-parameter-distribution': product_distribution("@signal_prior", "@nuisance_prior"), 'quantiles': [quantile], 'iterations': 20000 }
    if 'mcmc_iterations' in options: bayes_interval['iterations'] = options['mcmc_iterations']
    cfg_options = {'plugin_files': ('$THETA_DIR/lib/core-plugins.so',)}
    toplevel_settings = {'signal_prior': signal_prior, 'main': main, 'options': cfg_options}
    for i in range(mcmc_nproducers):
        if 'mcmc_seed' in options: bayes_interval['rnd_gen'] = {'seed': options['mcmc_seed'] + i }
        bayes_interval['name'] = 'bayes%d' % i
        toplevel_settings['bayes_interval%d' % i] = copy.deepcopy(bayes_interval)
    main['data_source'], toplevel_settings['model-distribution-signal'] = data_source_dict(model, input)
    cfg_names_to_run = []
    for sp in signal_processes:
        model_parameters = model.get_parameters(sp)
        toplevel_settings['nuisance_prior'] = nuisance_prior.get_cfg(model_parameters)
        name = write_cfg(model, sp, 'bayesian_quantiles', input, additional_settings = toplevel_settings)
        cfg_names_to_run.append(name)
    if 'run_theta' not in options or options['run_theta']:
        run_theta(cfg_names_to_run)
    else: return None
    cachedir = os.path.join(config.workdir, 'cache')
    plotsdir = os.path.join(config.workdir, 'plots')
    if not os.path.exists(plotsdir): os.mkdir(plotsdir)
    
    result = {}
    result_table = Report.table()
    result_table.add_column('process', 'signal process')
    if input=='data': header = '%f quantile' % quantile
    else: header = '%f quantile (median; central 68%%; central 95%%)' % quantile
    result_table.add_column('%f quantile' % quantile, header)
    for name in cfg_names_to_run:
        method, sp_id, dummy = name.split('-',2)
        result_table.set_column('process', sp_id)
        sqlfile = os.path.join(cachedir, '%s.db' % name)
        cols = ['bayes%d__quant%0.5d' % (d, quantile * 10000) for d in range(mcmc_nproducers)]
        data = sql(sqlfile, 'select %s from products' % ', '.join(cols))
        if mcmc_nproducers > 1:
            result[sp_id] = []
            result['%s__uncertainty' % sp_id] = []
            for row in data:
                mean, width = get_mean_width(row)
                result[sp_id].append(mean)
                result['%s__uncertainty' % sp_id].append(width / math.sqrt(mcmc_nproducers))
        else:
            result[sp_id] = [row[0] for row in data]
        outcomes = result[sp_id]
        n = len(outcomes)
        if n == 0:
           result_table.set_column('%f quantile' % quantile, 'N/A')
           continue
        if input == 'data':
            mean, width = get_mean_width(outcomes)
            result_table.set_column('%f quantile' % quantile, '%.5g +- %.3g' % (mean, width / math.sqrt(n)))
        else:
            sorted_res = sorted(outcomes)
            result_table.set_column('%f quantile' % quantile, '%.3g  (%.3g, %.3g; %.3g, %.3g)' % (sorted_res[n / 2], sorted_res[int(0.16 * n)], sorted_res[int(0.84 * n)],
                sorted_res[int(0.025 * n)], sorted_res[int(0.975 * n)]))
        result_table.add_row()
    if write_report:
        config.report.new_section("Bayesian Quantiles on ensemble '%s'" % input)
        config.report.add_p('input: %s, n: %d, signal prior: %s, nuisance prior: %s' % (input, n, str(signal_prior), str(nuisance_prior)))
        config.report.add_html(result_table.html())
    return result

# quantiles_toys is the result from bayesian_quantiles with input = 'toys:X'; this is used as 'expected' limit band
# quantiles_data is the result from bayesian_quantiles with input = 'data'; the median of this is used as solid line.
#   If quantiles_data is None, the quantiles_toys is used, i.e., the median expected
#   quantile is the quantile to use (usually 0.95)
#
# name is used to distinguish different calls to this function. This is used to make report section headings and plot filenames unique.
#
# options:
# - signalprocess_to_value: a dictionary mapping signal process name to values to be used as x axis for the band plot. As default, the first integer
#     in the signal process name is used for the x axis value.
#
# returns one plotutil.plotdata instance containing the 'observed' (or median expected) limit and 'expected' bands.
def limit_band_plot(quantiles_toys, quantiles_data = None, quantile = 0.95, name = '', **options):
    #expected results maps (process name) -> (median, band1, band2)
    # where band1 and band2 are tuples for the central 68 and 95%, resp.
    expected_results = {}
    observed_results = {}
    signal_processes = set()
    for sp in quantiles_toys:
        # ignore uncertainties or other special entries:
        if '__' in sp: continue
        signal_processes.add(sp)
        data = sorted(quantiles_toys[sp])
        n = len(data)
        if n!=0:
            expected_results[sp] = (data[n / 2], (data[int(0.16 * n)], data[int(0.84 * n)]), (data[int(0.025 * n)], data[int(0.975 * n)]))
        if quantiles_data is not None:
            data = sorted(quantiles_data[sp])
        if len(data) != 0:
            observed_results[sp] = data[len(data)/2]
    # map process names and x values:
    sp_to_x = {}
    x_to_sp = {}
    signal_processes = sorted(list(signal_processes))
    next_x = 0
    for sp in signal_processes:
        if 'signalprocess_to_value' in options and sp in options['signalprocess_to_value']: x = options['signalprocess_to_value'][sp]
        else: x = extract_number(sp)
        if x is None:
            print "WARNING: cannot find ordering value for signal process '%s', using %d" % (sp, next_x)
            x = next_x
            next_x += 1
        sp_to_x[sp] = x
        x_to_sp[x] = sp
    pd = plotdata()
    pd.color = '#000000'
    pd.as_function = True
    pd.x = sorted(list(x_to_sp.keys()))
    pd.y = []
    pd.bands = [([], [], '#00ff00'), ([], [], '#00aa00')]
    for x in pd.x:
        sp = x_to_sp[x]
        if sp in expected_results: median, band1, band2 = expected_results[sp]
        else: median, band1, band2 = 0, (0,0), (0,0)
        pd.y.append(observed_results[sp])
        pd.bands[1][0].append(band1[0])
        pd.bands[1][1].append(band1[1])
        pd.bands[0][0].append(band2[0])
        pd.bands[0][1].append(band2[1])
    plotsdir = os.path.join(config.workdir, 'plots')
    if len(pd.x) >=2:
        plot((pd,), 'signal process', '95% C.L. upper limit', os.path.join(plotsdir, 'limit_band_plot-%s.png' % name))
        config.report.new_section('Limit band Plot %s' % name)
        config.report.add_html('<p><img src="plots/limit_band_plot-%s.png" /></p>' % name)
    return pd

# options are as for bayesian_quantiles, with the modification that
# * additional option 'nodata': if exists and True, only calculate expected limits
# * there should be n_toy, n_obs  instead of n;  'n' will be ignored
# * 'input' will be ignored
def bayesian_limits(model, **options):
    if 'n' in options: del options['n']
    if 'input' in options: del options['input']
    n = 1000
    if 'n_toy' in options: n = options['n_toy']
    expected_limits = bayesian_quantiles(model, input = 'toys:0', n = n, **options)
    if model.has_data() and ('nodata' not in options or not options['nodata']):
        n = 10
        if 'n_data' in options: n = options['n_data']
        observed_limits = bayesian_quantiles(model, input = 'data', n = n, **options)
    else:
        observed_limits = expected_limits
    # e.g., if run_theta=False
    if expected_limits is None or observed_limits is None: return None
    return limit_band_plot(expected_limits, observed_limits)


#
# sp_group is a list of strings: the signal process names to scale with beta_signal
# options: same as bayesian_quantiles (in particular signal_processes!):
# * 'input' is ignored / replaced by settings according to beta_signal_range, beta_signal_n
def bayesian_limits_coveragetest(model, beta_signal_range = [0.0, 10.0], beta_signal_n = 10, **options): pass
    

# calculate the posteriors of all parameters
# model, input, n, signal_prior, nuisance_prior as as in bayesian_quantiles.
# histogram_specs
# * is a dictionary of (parameter name) -> (nbins, xmin, max)
#
# returns a dictionary (parameter name) -> (list of plotutil.plotdata)
def posteriors(model, histogram_specs, input = 'data', n = 3, signal_prior = 'flat', nuisance_prior = '', signal_processes = None, mcmc_iterations = 10000, **options):
    if signal_processes is None: signal_processes = [[sp] for sp in model.signal_processes]
    signal_prior = signal_prior_dict(signal_prior)
    nuisance_prior = nuisance_prior_distribution(model, nuisance_prior)
    main = {'n-events': n, 'model': '@model', 'producers': ('@posteriors',), 'output_database': sqlite_database(), 'log-report': False}
    posteriors = {'type': 'mcmc_posterior_histo', 'name': 'posteriors', 'parameters': [],
       'override-parameter-distribution': product_distribution("@signal_prior", "@nuisance_prior"),
       'smooth': True, 'iterations': mcmc_iterations }
    for par in histogram_specs:
        posteriors['parameters'].append(par)
        nbins, xmin, xmax = histogram_specs[par]
        posteriors['histo_%s' % par] = {'range': [float(xmin), float(xmax)], 'nbins': nbins}
    cfg_options = {'plugin_files': ('$THETA_DIR/lib/core-plugins.so',)}
    toplevel_settings = {'signal_prior': signal_prior, 'posteriors': posteriors, 'main': main, 'options': cfg_options}
    cfg_names_to_run = []
    main['data_source'], toplevel_settings['model-distribution-signal'] = data_source_dict(model, input)
    cfg_names_to_run = []
    for sp in signal_processes:
        model_parameters = model.get_parameters(sp)
        toplevel_settings['nuisance_prior'] = nuisance_prior.get_cfg(model_parameters)
        name = write_cfg(model, sp, 'posteriors', input, additional_settings = toplevel_settings)
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
        i = 0
        for par in parameters:
            result[par] = [row[i] for row in data]
            for pd in result[par]: pd.as_function = True
            plot(result[par], par, 'posterior density', os.path.join(plotsdir, '%s-%s.png' % (name, par)))
            config.report.add_html('<p>%s:<br/><img src="plots/%s-%s.png" /></p>' % (par, name, par))
            i += 1
            maxima = sorted(map(argmax, result[par]))
            result_table.set_column('maximum posterior %s' % par, '%.3g' % maxima[int(0.5 * len(maxima))])
        result_table.add_row()
    config.report.add_p('input: %s, n: %d, signal prior: %s, nuisance prior: %s' % (input, n, str(signal_prior), str(nuisance_prior)))
    config.report.add_html(result_table.html())
    return result

