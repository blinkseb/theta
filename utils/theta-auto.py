#!/usr/bin/env python

import os, os.path, datetime, re, sys, copy, traceback, shutil, hashlib

from theta_auto import *
inf = float("inf")

#TODO:
# * support for studies using +-1sigma toy / asimov data as input and re-run the method --> input='toys:scan-nuisance[-asimov]'?


## \page theta_auto_intro Introduction to theta-auto
# 
#  The theta-auto.py script provides a set of functions to automatically generate theta configuration files, run theta on them, and analyze the output.
#  This page assumes you have some basic working knowledge of theta, concerning the configuration file syntax, and running theta.
#
#  With the theta-auto scripts, some further restrictions apply regarding the statistical model, compared to the capabilities of using theta directly.
#  The prototype model is a model with multiple channels where in each channel, the predicted events yield is the sum of the background and signal processes,
#  where each process is affected by shape uncertainties
#  modeled with cubiclinear_histomorph and rate uncertainties modeled via multiplicative log-normal priors. While some simple extensions to this prototype model
#  are possible (such as using truncated Gaussian priors instead of log-normal ones for the uncertainties, or using a different kind of histogram morphing),
#  this naturally restricts the use case to limit setting, discovery and measurement for a cross-section type parameter. On the other hand, this restriction
#  allows the python scripts to make many assumptions.
#
#  The easiest way to get started is probably to specify the model as a datacard as documented in
#  https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideHiggsAnalysisCombinedLimit
#  If you have such a datacard called "datacard.txt", put it in the same directory as a script called "analysis.py" with the content
#  \code
#  model = higgs_datacard.build_model('datacard.txt')
#  result = discovery(model)
#  print result
#  report.write_html('/home/username/public_html/datacard/')
#  \endcode
#  Then, from this same directory, call the theta-auto.py script (without any arguments).
#
#  This example shows the basic work flow of an analysis:
#  <ol>
#   <li>Create the external file(s) required to specify the model. This can be done using the Higgs group datacard format; other formats are discussed later</li>
#   <li>Create the "analysis.py" script for this analysis. This usually consists of three simple steps: (i) read in the model, (ii) call a statistical
#      method such as limit setting, discovery, etc, (iii) write output.</li>
#   <li>Call theta-auto.py. This will create (in a subdirectory called after the analysis script, in this case "analysis") the theta configuration files,
#       call theta, and analyze the result.</li>
#  </ol>
#
# Each function returns the result as a return value. The type of returned object depends on the function, see the function documentation
# for details. Most functions also writes a summary of the result to a global object called "report" which you can instruct to write out
# its content as html in a specified directory; this should be done at the very end of the "analysis.py" script.
#
# Do not execute the "analysis.py" script directly. Instead pass the name of the script as first argument to "theta-auto.py". In this case, also the created
# subdirectory, which contains the automatically produced theta configuration files and other intermediate output, will change to the
# script name without trailing ".py".
#
# The created subdirectory contains intermediate results such as the generated theta configuration files and sqlite files with the result.
# You can (and from time to time should!) delete it, everything will be re-generated. However, it makes sense to keep the cached result because
# theta-auto does not run theta again if there is a result created by theta based on the same configuration file. While this does not
# help when making changes to the model, it can be very helpful if the model and the method are set up and you mainly work on
# the extraction of the result.
#
# The most important statistical methods are:
# <ol>
#  <li>\link theta_auto::ml_fit ml_fit \endlink: make a maximum likelihood fit for all model parameters</li>
#  <li>\link theta_auto::pl_intervals pl_intervals \endlink: calculate profile likelihood intervals for the signal cross section; this includes a maximum likelihood fit</li>
#  <li>\link theta_auto::discovery discovery\endlink: calculate the distribution of a test statistic (usually likelihood ratio) for the background-only hypothesis to get the p-value /
#      number of standard deviations "sigma".</li>
#  <li>\link theta_auto::bayesian_quantiles bayesian_quantiles\endlink: calculate marginal posterior quantiles for the signal
#        strength parameter = upper / lower Bayesian limits on the signal cross section</li>
#  <li>\link theta_auto::posteriors posteriors \endlink: calculate the marginal posteriors for the given parameters (which can be either nuisance parameters or the signal cross section)</li>
#  <li>\link theta_auto::cls_limits cls_limits \endlink: calculate CLs limits, using a certain test statistic (usually a variant of profile likelihood ratio)</li>
#  <li>\link theta_auto::ks_test ks_test \endlink: calculate the KS p-value by generating toys, running a maximum likelihood fit, and calculating the
#    KS test statistic. This is useful to get KS p-values
#    which both (i) are valid after fitting data, (ii) include the effect of systematic uncertainties (both shape and rate).</li>
# </ol>
#
# In addition, there are some more, higher-level methods which typically internally call one or more of the above. For example, the function bayesian_limits
# calculates the limits for different signal hypotheses (such as signals for different mass points) and constructs a limit band plot (= a plot
# with median, +-1sigma and +-2sigma expected and observed limit) as a function of the signal process.
#
#
# \section theta_auto_model The Model in theta-auto
#
# The Model class contains all relevant information of the model, including the observed data, and all the predictions including their dependence on the
# model parameters in the different channels, and which of the processes is to be considered as signal.
# The signal is always added to the prediction of the background yields with a multiplicative signal strength factor called "beta_signal" which corresponds
# to the cross section in units of the configured signal yield. All model parameters except "beta_signal" are called "nuisance parameters" throughout the documentation.
#
# Note that the python Model class in theta-auto is slightly different from the C++ theta::Model class you spefify in the %theta
# configuration files: theta::Model does not contain any information about data, nor about which process is considered as signal.
#
# The Model class also includes an instance of Distribution, model.distribution which is the prior distribution
# for all nuisance parameters. Toy data generation is always based on this prior. It is also used to add additional terms to the likelihood function, although
# this can be overridden by the parameter \c nuisance_constraint, which, if present, takes precedence over model.distribution; see \ref theta_auto_params below.
#
# An instance of the Model class is usually not constructed directly. Rather, use one of the functions which generate a model based on some external
# configuration like higgs_datacard.build_model or build_model_from_rootfile.
#
# It is possible to manipulate a Model instance with many functions. Important examples are to add a (log-normal) rate uncertainty for a certain process which
# can be done with
# \code
#  model.add_lognormal_uncertainty('ttbar_rate', math.log(1.12), 'ttbar', '*')
# \endcode
# which adds a 12% log-normal uncertainty controlled by the uncertainty 'ttbar_rate' on the process called 'ttbar', correlated across all channels ('*').
#
# Another important example is the combination of two Models which can be done via
# \code
# model1.combine(model2)
# \endcode
# which will add all channels of model2 to model1. In order for this to work correctly (=with correct correlations induced by shared uncertainties),
# model1 and model2 must have been built using the same convention for the names of nuisance parameters. Also, the signal processes
# should be the same for both models.
#
# For more information about Model manipulation, see the documentation of the Model class.
#
# \section theta_auto_signal Specifying what is 'signal'
#
# The Model and statistical methods in theta-auto are constructed such that it is easy to run them for situations in which there is more than one signal hypothesis, e.g.,
# a Zprime with different masses. In this case, the background model is always the same, and as signal model, one of the signal processes is chosen; the statistical
# procedure is applied for each of the signal processes in turn. For example, one might have the signal processes 'zp1000', 'zp2000', and 'zp3000' for Zprime models
# with different masses.
#
# In order for theta-auto to be able to correctly construct the background-only model, it must know what is signal. Use a call like
# \code
# model.set_signal_processes('zp*')
# \endcode
# to declare which processes exactly are the signal processes. Everything else is considered the background-only model.
#
# Sometimes, it makes sense to split the signal into several components which are conceptually the same process of interest
# (e.g., if the process was generated according to production process or decay or if the uncertainties for different signal components are different). For example,
# one might have split the production and processing of the Zprime signal sample according to production process into 'ggzp1000' and 'qqzp1000', etc., and in the example above
# with three masses one would have six signal processes in total.
# In such a case, one can specify more than one process name, all of which are scaled with the cross section scaling parameter "beta_signal". Such a group of processes to be
# used together at the same time is called "signal process group" and is specified via a list of strings. For different signal hypotheses, such as different masses,
# one has several such signal process groups. Therefore, a complete specification for which processes to run the statistical method on is a list of signal process groups.
# See also the documentation below of the \c signal_processes parameter. To run the statistical model on all three mass hypotheses from the above example, one would
# use
# \code
#  signal_processes = [['qqzp1000', 'ggzp1000'], ['qqzp2000', 'ggzp2000'], ['qqzp3000', 'ggzp3000']]
# \endcode
# For some purposes (such as a key in a dictionary return value or if used as part of a file name), identifying a signal process group via a list of
# strings is not possible. In this case, all process names of a signal process group are concatenated to for a "process group id". The process group
# ids in the above example are
# \code
#  'qqzp1000ggzp1000', 'qqzp2000ggzp2000', 'qqzp3000ggzp3000'
# \endcode
#
# \section theta_auto_params Common Parameters
#
# Almost all statistical methods allow to specify whether to run the method on the actually observed data or on toys. They also
# allow to specify alternative priors/likelihood terms (in theta, no distinction is made between these two). These
# parameters are:
# <ul>
#   <li>\c input: on which set of data to run the method. Can be either "data" or "toys:X" where X is the
#       signal scale factor; for example "toys:0" means no signal=background-only.</li>
#   <li>\c signal_prior: the signal prior to use for the likelihood function definition. Valid values are "fix:X", where X is the value to fix this parameter to, "flat" for
#     a flat prior on [0,inf] and "flat:[X,Y]" for a flat prior on the interval [X, Y]. For more complicated priors, you can directly specify the theta
#     configuration as dictionary.</li>
#   <li>\c nuisance_constraint: the signal prior to use for the likelihood function definition. Valid values are '' (=use model.distributions unmodified),
#     "shape:fix" (=fix shape-changing parameter, leave other parameter unchanged), "shape:free" (=use a flat prior for shape-changing parameters). Similarly, there are
#     "rate:fix", "rate:free", and combinations of "shape:..." and "rate:..." combinations, separated by semicolon, e.g., "shape:fixed;rate:free". For more complicated
#     situations, you can pass here an instance of the Distribution class which will be merged with the model.distribution (see Distribution.merge).</li>
#   <li>\c signal_processes: a list of signal process groups to considered together as signal. A signal process group is a list of strings of
#       the names of the signal processes. They all must have been declared as signal processes (via Model.set_signal_processes) in the model. The default (None)
#       will use all signals individually as a 'trivial' signal process group consisting of only this signal.</li>
# </ul>
#
# \section Some Internals
#
# This section describes some internals of theta-auto which are useful to know even if you do not want to modify theta-auto but only use it.
#
# The generated configuration file name follows the convention: <method>-<signal process group id>-<input>[-<id>]-<hash>.cfg
# where <method> is the python function name of the function generating the configuration file, <input>
# is the 'input' parameter as explained in the previous section, the optional <id> depends on the function, and <hash>
# is (part of a md5-)hash of the configuration file, in order to avoid clashes where otherwise the same values apply.
#
# Each configuration file creates one sqlite output file with the same name in the current directory; the suffix ".cfg"
# is replaced by ".db".


## \file theta-auto.py
#
# Main executable for the theta-auto scripts. See \ref theta_auto_intro for details.

## \namespace theta_auto
# \brief Python scripts of theta-auto
#


## \brief Perform a maximum likelihood fit
#
# Finds the parameter values for all model parameters at the maximum of the likelihood, using
# the theta plugin opf type 'mle'
#
# The parameters \c input, \c signal_prior, \c nuisance_constraint, and \c signal_processes are documented in more detail
# in \ref theta_auto_params.
#
# In the report, a table with the parameter values at the maximum is created, with one entry per signal process group.
#
# \param n The number of fits to run. Values &gtr; 1 only make sense for input="toys"
# \return A dictionary with the signal process group id as key (see \ref theta_auto_signal on the definition of signal process group id).
# The value is a dictionary with the parameter name as key and a list of length n with two-tuple of floats, (parameter value, uncertainty) as value.
# Note that in case the maximization does not converge, these toys will be omitted from the result and the lists of two-tuples have length
# smaller than \c n.
def ml_fit(model, input = 'data', signal_prior = 'flat', nuisance_constraint = 'shape:fix', signal_processes = None, n = 1, **options):
    if signal_processes is None: signal_processes = [[sp] for sp in model.signal_processes]
    beta_signal_value = 0.0
    nuisance_constraint = nuisance_prior_distribution(model, nuisance_constraint)
    signal_prior_spec = signal_prior
    signal_prior = signal_prior_dict(signal_prior)
    main = {'n-events': n, 'model': '@model', 'producers': ('@mle',), 'output_database': sqlite_database(), 'log-report': False}
    mle = {'type': 'mle', 'name': 'mle', 'parameters': None,
       'override-parameter-distribution': product_distribution("@signal_prior", "@nuisance_constraint"),
       'minimizer': minimizer()}
    cfg_options = {'plugin_files': ('$THETA_DIR/lib/core-plugins.so','$THETA_DIR/lib/root.so')}
    toplevel_settings = {'model-distribution-signal': delta_distribution(beta_signal = beta_signal_value), 'mle': mle, 'main': main, 'signal_prior':
        signal_prior_dict(signal_prior),  'options': cfg_options}
    main['data_source'], toplevel_settings['model-distribution-signal'] = data_source_dict(model, input)
    cfg_names_to_run = []
    for sp in signal_processes:
        model_parameters = sorted(list(model.get_parameters(sp)))
        mle['parameters'] = model_parameters
        if 'beta_signal' in model_parameters: mle['override-parameter-distribution'] = product_distribution("@signal_prior", "@nuisance_constraint")
        else: mle['override-parameter-distribution'] = "@nuisance_constraint"
        toplevel_settings['nuisance_constraint'] = nuisance_constraint.get_cfg(model_parameters)
        name = write_cfg(model, sp, 'ml_fit', input, additional_settings = toplevel_settings)
        cfg_names_to_run.append(name)
    if 'run_theta' not in options or options['run_theta']:
        run_theta(cfg_names_to_run)
    else: return None
    cachedir = os.path.join(config.workdir, 'cache')
    
    result = {}
    result_table = table()
    result_table.add_column('process', 'signal process')
    
    nuisance_parameters = sorted(list(model.get_parameters('')))
    for p in nuisance_parameters:
        suffix = ''
        if nuisance_constraint.get_distribution(p)['width'] == 0.0: suffix = ' (fixed)'
        result_table.add_column(p, '%s%s' % (p, suffix))
    suffix = ''
    if signal_prior_spec.startswith('fix:'): suffix = ' (fixed)'
    result_table.add_column('beta_signal', 'beta_signal%s' % suffix)
    for i in range(len(cfg_names_to_run)):
        sp = signal_processes[i]
        name = cfg_names_to_run[i]
        method, sp_id, dummy = name.split('-',2)
        result[sp_id] = {}
        result_table.set_column('process', sp_id)
        model_parameters = model.get_parameters(sp)
        sqlfile = os.path.join(cachedir, '%s.db' % name)
        cols = ['mle__%s, mle__%s_error' % (p, p) for p in model_parameters]
        data = sql(sqlfile, 'select %s from products' % ', '.join(cols))
        if len(data) == 0: raise RuntimeError, "no data in result file '%s'" % sqlfile
        i = 0
        for p in model_parameters:
            result[sp_id][p] = [(row[2*i], row[2*i+1]) for row in data]
            i += 1
            sorted_res = sorted([res[0] for res in result[sp_id][p]])
            n = len(sorted_res)
            if n >= 10:
                result_table.set_column(p, '%.3g (%.3g, %.3g)' % (sorted_res[int(0.5*n)], sorted_res[int(0.16*n)], sorted_res[int(0.84*n)]))
            else: result_table.set_column(p, '%.3g' % sorted_res[int(0.5*n)])
        for p in nuisance_parameters + ['beta_signal']:
            if p in model_parameters: continue
            result_table.set_column(p, 'n/a')
        result_table.add_row()
    config.report.new_section("Maximum Likelihood fit on ensemble '%s'" % input)
    config.report.add_p('The table entries give the median (and, if n>=10 the central 68%) of the parameter values at the found maximum of the likelihood function.')
    #config.report.add_p('input: %s, n: %d, signal prior: %s, nuisance prior: %s' % (input, n, str(signal_prior), str(nuisance_constraint)))
    config.report.add_html(result_table.html())
    return result
    
## \brief Perfoem a maximum likelihood fit and get the coefficient function values for all processes / channels
#
# 
#
# options: see ml_fit
#
# returns a dictionary
# (signal process) --> (observable name) --> (process name) --> (factor)
#
# Does not write anything to the report.
def ml_fit_coefficiencts(model, **options):
    result = {}
    res = ml_fit(model, **options)
    for sp in res:
        values = {}
        for param in res[sp]:
            values[param] = res[sp][param][0][0]
        result[sp] = {}
        for obs in model.observables:
            result[sp][obs] = {}
            for proc in model.get_processes(obs):
                # skip signal processes we are not interested in:
                if proc in model.signal_processes and proc not in sp: continue
                result[sp][obs][proc] = model.get_coeff(obs, proc).get_value(values)
    return result
            

## \brief Perform a KS-test on th whole model
#
# Perform a KS-test by (i) dicing toy data from the model (including nuisance parameters according to model.distribution),
# (ii) performing a maximum likelihood fit using nuisance_constraint and (iii) calculating the Kolmogorov-Smirnov test-statistic value comparing the
# prediction and data after the fit.
#
# This method always uses the background-only model. Therefore, it has no signal_priori or signal_processes parameters.
# If you want to include signal in the KS-test, you have to define it as background (via Model.set_signal_processes) before calling\
# this function.
#
# Returns the p-value for data. Does not write anything to the report.
def ks_test(model, n = 5000, nuisance_constraint = 'shape:fix', **options):
    nuisance_constraint = nuisance_prior_distribution(model, nuisance_constraint)
    main = {'n-events': n, 'model': '@model', 'producers': ('@mle',), 'output_database': sqlite_database(), 'log-report': False}
    mle = {'type': 'mle', 'name': 'mle', 'parameters': [], 'write_ks_ts': True,
       'override-parameter-distribution': product_distribution("@nuisance_constraint"),
       'minimizer': minimizer(need_error = False)}
    cfg_options = {'plugin_files': ('$THETA_DIR/lib/core-plugins.so','$THETA_DIR/lib/root.so')}
    toplevel_settings = {'mle': mle, 'main': main, 'options': cfg_options}
    cfg_names_to_run = []
    model_parameters = sorted(list(model.get_parameters('')))
    toplevel_settings['nuisance_constraint'] = nuisance_constraint.get_cfg(model_parameters)
    id = None
    if 'id' in options: id = options['id']
    # for toys:
    input = 'toys:0'
    main['data_source'], dummy = data_source_dict(model, input)
    name = write_cfg(model, '', 'ks_test', input, id=id, additional_settings = toplevel_settings)
    cfg_names_to_run.append(name)
    # for data:
    input = 'data'
    main['data_source'], dummy = data_source_dict(model, input)
    main['n-events'] = 1
    name = write_cfg(model, '', 'ks_test', input, id=id, additional_settings = toplevel_settings)
    cfg_names_to_run.append(name)
    if 'run_theta' not in options or options['run_theta']:
        run_theta(cfg_names_to_run)
    else: return None
    cachedir = os.path.join(config.workdir, 'cache')

    name_toys = cfg_names_to_run[0]
    sqlfile = os.path.join(cachedir, '%s.db' % name_toys)
    data = sql(sqlfile, 'select mle__ks_ts from products')
    if len(data) == 0: raise RuntimeError, "no data in result file '%s'" % sqlfile
    ks_values_bkg = [row[0] for row in data]
    
    name_data = cfg_names_to_run[1]
    sqlfile = os.path.join(cachedir, '%s.db' % name_toys)
    data = sql(sqlfile, 'select mle__ks_ts from products')
    if len(data) == 0: raise RuntimeError, "no data in result file '%s'" % sqlfile
    ks_value_data = data[0][0]
    p = len([x for x in ks_values_bkg if x >= ks_value_data]) * 1.0 / len(ks_values_bkg)
    return p
    

## \brief Perform a KS compatibility test for all individual channels of the given model
#
# Each channel is treated independently. For each channel in the model, a "restricted" model is built
# which contains each only that channel. This restricted model and all options are passed to to ks_test, see there
# for valid options.
#
# If fit_yield is True, one parameter is added to the model which is used to scale all processes
# at the same time (=the fraction between processes is not modified, at least not by this parameter).
# In this case, after the fit, the normalisation of simulation and data coincides by construction
# and the KS-test effectively compares the shapes only, not the rate.
# Note that even if fit_yield is False, there is still a maximum likelihood fit done which finds the
# parameter values at the maximum of the likelihood function, using nuisance_constraints (details see ks_test).
# Depending on the model, this can already mean that the yield is fitted.
#
# Returns a dictionary (channel name) --> (p-value).
#
# Does not write anything to the report.
def ks_test_individual_channels(model, fit_yield = False, **options):
    observables = model.observables.keys()
    result = {}
    for obs in observables:
        model_obs = copy.deepcopy(model)
        model_obs.restrict_to_observables([obs])
        if fit_yield:
            model_obs.distribution.set_distribution('scale_', 'gauss', 1.0, inf, (0, inf))
            procs = model_obs.get_processes(obs)
            for p in procs:
                model_obs.get_coeff(obs, p).add_factor('id', parameter = 'scale_')
        options['id'] = obs
        result[obs] = ks_test(model_obs, **options)
    return result
    

## \brief Calculate profile likelihood intervals
#
# Calculate profile likelihood intervals using the deltanll_intervals plugin in %theta.
#
# For the \c input, \c nuisance_constraint, \c signal_processes parameters, see \ref theta_auto_params.
#
# \param n The number of toys. A value n &gt; 1 only makes sense for input='toys'
# \param cls A list of confidence levels
# \param write_report If True, a summary will be written to the report. If None, the report will be written if and only if input is 'data'
#
# \return A dictionary with the signal process group id as key. The value is a dictionary with the confidence levels as keys; the values are lists of
#  two-tuples (lower interval border, upper interval border). In addition to the confidence levels, there is a special key 'mle' which contains a list
#  of beta_signal values at the maximum.
#
# For example, if the only signal process is called 'signal' and cl=[0.68, 0.84], the 1sigma intervals are
# \code
#  result['signal'][0.68][i] # i=0..n-1
# \endcode
# and the 2signma intervals are
# \code
#  result['signal'][0.84][i] # i=0..n-1
# \endcode
# and the beta_signal values at the maximum are 
# \code
#  result['signal']['mle'][i] # i=0..n-1
# \endcode
#
# In case the minimization fails, the lists have less than n entries.
#
# If write_report is true, it will only write out the result of the first fit. This usually makes only sense for data,
# not for toy MC.
def pl_intervals(model, input = 'toys:0', n = 100, nuisance_constraint = '', cls = [0.90], signal_processes = None, write_report = None, **options):
    if signal_processes is None: signal_processes = [[sp] for sp in model.signal_processes]
    if write_report is None: write_report = input == 'data'
    beta_signal_value = 0.0
    if input.startswith('toys:'): beta_signal_value = float(input[5:])
    elif input!='data': raise RuntimeError, "unexpected value for 'input': %s (expected 'toys:...' or 'data')" % input
    signal_prior_dict = {'type': 'flat_distribution', 'beta_signal': {'range': [1e-12, float("inf")], 'fix-sample-value': 1.0}}
    nuisance_constraint = nuisance_prior_distribution(model, nuisance_constraint)
    main = {'n-events': n, 'model': '@model', 'producers': ('@deltanll_interval',), 'output_database': sqlite_database(), 'log-report': False}
    deltanll_interval = {'type': 'deltanll_intervals', 'name': 'deltanll', 'parameter': 'beta_signal',
       'override-parameter-distribution': product_distribution("@signal_prior", "@nuisance_constraint"),
       'clevels': cls, 'minimizer': minimizer(need_error = False)}
    cfg_options = {'plugin_files': ('$THETA_DIR/lib/core-plugins.so','$THETA_DIR/lib/root.so')}
    toplevel_settings = {'signal_prior': signal_prior_dict, 
                'model-distribution-signal': delta_distribution(beta_signal = beta_signal_value), 'deltanll_interval': deltanll_interval, 'main': main,
                'options': cfg_options}
    cfg_names_to_run = []
    main['data_source'], toplevel_settings['model-distribution-signal'] = data_source_dict(model, input)
    for sp in signal_processes:
        model_parameters = model.get_parameters(sp)
        toplevel_settings['nuisance_constraint'] = nuisance_constraint.get_cfg(model_parameters)
        name = write_cfg(model, sp, 'deltanll_intervals', input, additional_settings = toplevel_settings)
        cfg_names_to_run.append(name)
    if 'run_theta' not in options or options['run_theta']:
        run_theta(cfg_names_to_run)
    else: return None
    
    result_table = Report.table()
    result_table.add_column('signal process')
    result_table.add_column('mle', 'mle')
    for cl in cls: result_table.add_column('cl%g' % cl, 'confidence level %g' % cl)
    cachedir = os.path.join(config.workdir, 'cache')
    col_suffixes= ['%05d' % (cl*10000) for cl in cls]
    result = {}
    for i in range(len(cfg_names_to_run)):
        name = cfg_names_to_run[i]
        method, sp_id, dummy = name.split('-',2)
        result_table.set_column('signal process', sp_id)
        result[sp_id] = {'mle': []}
        for cl in cls: result[sp_id][cl] = []
        model_parameters = model.get_parameters(sp)
        sqlfile = os.path.join(cachedir, '%s.db' % name)
        colnames = []
        for cs in col_suffixes:
            colnames.append('deltanll__lower%s' % cs)
            colnames.append('deltanll__upper%s' % cs)
        data = sql(sqlfile, 'select deltanll__maxl, %s from products' % (', '.join(colnames)))
        first_row = True
        for row in data:
            result[sp_id]['mle'].append(row[0])
            if first_row: result_table.set_column('mle', '%g' % row[0])
            for icol in range(len(cls)):
                interval = (row[2*icol+1], row[2*icol+2])
                result[sp_id][cls[icol]].append(interval)
                if first_row:
                    result_table.set_column('cl%g' % cls[icol], '(%g, %g)' % interval)
            first_row = False
        result_table.add_row()
    if write_report:
        config.report.new_section("deltanll intervals")
        config.report.add_html(result_table.html())
    return result

# runs deltanll_intervals and measures coverage as a function of the true beta_signal
# Returns: dictionary: (spid) --> (true beta signal) --> |  (cl) --> coverage
#                                                        |  'successrate' --> fit success rate
def pl_coveragetest(model, beta_signal_values = [0.2*i for i in range(10)], n = 1000, write_report = True, **deltanll_options):
    result = {}
    result_tables = {}
    for beta_signal in beta_signal_values:
        res = pl_intervals(model, input='toys:%g' % beta_signal, n = n, **deltanll_options)
        for spid in res:
            cls = [k for k in res[spid].keys() if type(k)==float]
            if spid not in result: result[spid] = {}
            if spid not in result_tables:
                 result_tables[spid] = Report.table()
                 result_tables[spid].add_column('beta_signal', 'true beta_signal')
                 for cl in cls: result_tables[spid].add_column('coverage %g' % cl, 'Coverage for cl=%g' % cl)
                 result_tables[spid].add_column('fit success fraction')
            result[spid][beta_signal] = {}
            result_tables[spid].set_column('beta_signal', '%g' % beta_signal)
            icl = 0
            for cl in cls:
                n_covered = len([1 for i in range(len(res[spid][cl])) if res[spid][cl][i][0] <= beta_signal and res[spid][cl][i][1] >= beta_signal])
                n_total = len(res[spid][cl])
                coverage = n_covered*1.0 / n_total
                result[spid][beta_signal][cl] = coverage
                successrate = (n_total*1.0 / n)
                result[spid][beta_signal]['successrate'] = successrate
                result_tables[spid].set_column('coverage %g' % cl, '%g' % coverage)
                result_tables[spid].set_column('fit success fraction', '%g' % successrate)
            result_tables[spid].add_row()
    if write_report:
        for spid in result_tables:
            config.report.new_section("deltanll interval coverage test for signal process '%s'" % spid)
            config.report.add_html(result_tables[spid].html())
    return result
    
def clean_workdir():
    shutil.rmtree(config.workdir, ignore_errors = True)
    setup_workdir()
    
def setup_workdir():
    if not os.path.exists(config.workdir): os.mkdir(config.workdir)
    if not os.path.exists(os.path.join(config.workdir, 'plots')): os.mkdir(os.path.join(config.workdir, 'plots'))
            
def main():
    scriptname = 'analysis.py'
    for arg in sys.argv[1:]:
        if '.py' in arg: scriptname = arg
    config.workdir = os.path.join(os.getcwd(), scriptname[:-3])
    config.workdir = os.path.realpath(config.workdir)    
    setup_workdir()
    config.theta_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), '..'))
    config.report = html_report(os.path.join(config.workdir, 'index.html'))
    variables = globals()
    variables['report'] = config.report
    utils.info("workdir is %s" % config.workdir)
    utils.info("cachedir is %s" % os.path.join(config.workdir, 'cache'))
    utils.info("executing script %s" % scriptname)
    try:
        execfile(scriptname, variables)
    except Exception as e:
        print "error while trying to execute analysis script %s:" % scriptname
        traceback.print_exc()
        sys.exit(1)
        
if __name__ == '__main__': main()



