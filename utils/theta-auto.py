#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, os.path, datetime, re, sys, copy, traceback, shutil, hashlib, tempfile

from theta_auto import *
inf = float("inf")

#TODO:
# * support for studies using +-1sigma toy / asimov data as input and re-run the method --> input='toys:scan-nuisance[-asimov]'?
# * support additional_nll_term everywhere:
#   - producers set additional-ll-term
#   - producers add additional terms in override-parameter-distribution
#  => make a base class Producer which contains takes care of these two!
# * switch Function to the one defined in FunctionBase


# base class for producers, providing consistent approach to override_parameter_distribution
# and additional_nll_term
class ProducerBase:
    def gef_cfg(self, parameters): pass
    
    # model_parameters is the list of parameters (not necesarily those of additional_nll_term)
    def get_cfg_base(self, model, signal_processes):
        result = {'name': self.name}
        if self.override_parameter_distribution is not None:
            parameters = set(model.get_parameters(signal_processes))
            if model.additional_nll_term is not None:
                parameters.update(model.additional_nll_term.get_parameters())
            result['override-parameter-distribution'] = self.override_parameter_distribution.get_cfg(parameters)
        if model.additional_nll_term is not None:
            result['additional-nll-term'] = model.additional_nll_term.get_cfg()
        return result
    
    def __init__(self, parameter_dist, name):
        self.override_parameter_distribution = parameter_dist
        self.name = name
    

class MleProducer(ProducerBase):
    def __init__(self, parameter_dist, name = 'mle'):
        ProducerBase.__init__(self, parameter_dist, name)
        
    # parameters_write is the list of parameters to write the mle for in the db. The default (None) means
    # to use all parameters.
    def get_cfg(self, model, signal_processes, parameters_write = None):
        model_parameters = model.get_parameters(signal_processes, True)
        if parameters_write is None: parameters_write = model_parameters
        result = {'type': 'mle', 'minimizer': minimizer(), 'parameters': list(parameters_write)}
        result.update(self.get_cfg_base(model, signal_processes))
        return result

class PliProducer(ProducerBase):
    def __init__(self, parameter_dist, cls, name = 'pli'):
        ProducerBase.__init__(self, parameter_dist, name)
        self.cls = cls[:]
        self.parameter = 'beta_signal'

    def get_cfg(self, model, signal_processes, parameters_write = None):
        result = {'type': 'deltanll_intervals', 'minimizer': minimizer(need_error = False), 'parameter': self.parameter, 'clevels': self.cls}
        result.update(self.get_cfg_base(model, signal_processes))
        return result


#the class for the configuration of the "main" path and common options.
class MainBase:
    def __init__(self, root_plugins = True):
        self.cfg_options = {'plugin_files': ['$THETA_DIR/lib/core-plugins.so']}
        if root_plugins: self.cfg_options['plugin_files'].append('$THETA_DIR/lib/root.so')
        
class Run(MainBase):
    def __init__(self, n_events, data_source_dict, model_dist, root_plugins = True):
        MainBase.__init__(self, root_plugins)
        self.n_events = n_events
        self.producers = []
        self.data_source_dict = data_source_dict
        self.model_dist = model_dist
        
    def add_producer(self, producer):
        self.producers.append(producer)
        
    def get_cfg(self, model, signal_processes):
        main = {'n-events': self.n_events, 'producers': [], 'log-report': False, 'output_database': sqlite_database(), 'model': "@model", 'data_source': self.data_source_dict}
        toplevel_settings = {'main': main, 'options': self.cfg_options, 'model': model.get_cfg(signal_processes)}
        model_parameters = model.get_parameters(signal_processes)
        toplevel_settings['model']['parameter-distribution'] = Distribution.merge(model.distribution, self.model_dist).get_cfg(model_parameters)
        for p in self.producers:
            toplevel_settings['p%s' % p.name] = p.get_cfg(model, signal_processes)
            main['producers'].append('@p%s' % p.name)
        return toplevel_settings


"""
# a statistical method 
class StatisticalMethod:
    # initialize the method to run the method on a dataset input specified by input_spec, using nuisance parameter prior
    # nuisance_prior and likelihood(/posterior) terms specified by likelihood_*
    def __init__(self, model, input_spec, n, model_nuisance_prior, likelihood_nuisance_prior='', likelihood_beta_signal_prior='flat'):
        #TODO: save / convert all!
"""


def model_signal_prior_dist(input_spec):
    if input_spec.startswith('toys:') or input_spec.startswith('toys-asimov:'):
        beta_signal_value = float(input_spec[input_spec.find(':') + 1:])
    else: beta_signal_value = 0.0
    result = Distribution()
    result.set_distribution('beta_signal', 'gauss', beta_signal_value, 0.0, [beta_signal_value, beta_signal_value])
    return result


def default_signal_processes(model, signal_processes):
    if signal_processes is not None: return signal_processes
    signal_processes = [[sp] for sp in model.signal_processes]
    if len(signal_processes)==0: signal_processes.append('')
    return signal_processes

def signal_prior_dist(spec):
    result = Distribution()
    if type(spec) == str:
        if spec.startswith('flat'):
            if spec.startswith('flat:'):
                res = re.match('flat:\[([^,]+),(.*)\]', spec)
                if res is None: raise RuntimeError, "signal_prior specification '%s' invalid (does not match range specification syntax)" % spec
                xmin, xmax = float(res.group(1)), float(res.group(2))
            else:
                if spec!='flat': raise RuntimeError, "signal_prior specification '%s' invalid" % spec
                xmin, xmax = 0.0, float("inf")
            value = 0.5 * (xmax - xmin)
            if value==float("inf"): value = 1.0
            result.set_distribution('beta_signal', 'gauss', value, float("inf"), [xmin, xmax])
            return result
        elif spec.startswith('fix:'):
            v = float(spec[4:])
            result.set_distribution('beta_signal', 'gauss', v, 0.0, [v,v])
            return result
        else: raise RuntimeError, "signal_prior specification '%s' unknown" % spec
    else:
        if not isinstance(spec,Distribution): raise RuntimeError, "signal_prior specification has to be a string ('flat' / 'fix:X') or a Distribution instance!"
        return spec


def write_cfg2(main, model, signal_processes, method, input, id = None, **options):
    all_parameters = model.get_parameters(signal_processes)
    if model.additional_nll_term is not None: all_parameters.update(model.additional_nll_term.get_parameters())
    all_parameters = sorted(list(all_parameters))
    theta_cfg = "parameters = " + settingvalue_to_cfg(all_parameters, 0, ['parameters']) + ";\n"
    obs = {}
    for o in model.observables:
        xmin, xmax, nbins = model.observables[o]
        obs[o] = {'range': [xmin, xmax], 'nbins': nbins}
    theta_cfg += "observables = " + settingvalue_to_cfg(obs, 0, ['observables']) + ";\n"
    cfg = main.get_cfg(model, signal_processes)
    cfg['main']['output_database']['filename'] = '@output_name';
    for s in cfg:
        theta_cfg += s + " = " + settingvalue_to_cfg(cfg[s]) + ";\n"
    m = hashlib.md5()
    m.update(theta_cfg)
    hash = m.hexdigest()[:10]
    if id is None:
        name = '%s-%s-%s-%s' % (method, ''.join(signal_processes), input, hash)
    else:
        name = '%s-%s-%s-%s-%s' % (method, ''.join(signal_processes), input, id, hash)
    f = open(os.path.join(global_config.workdir, name + '.cfg'), 'w')
    print >>f, theta_cfg
    print >>f, 'output_name = "%s.db";\n' % name
    f.close()
    return name

    

def ml_fit2(model, input = 'data', signal_prior = 'flat', nuisance_constraint = 'shape:fix', signal_processes = None, n = 1, **options):
    signal_processes = default_signal_processes(model, signal_processes)
    nuisance_constraint = nuisance_prior_distribution(model, nuisance_constraint)
    signal_prior_spec = signal_prior
    signal_prior = signal_prior_dist(signal_prior)
    model_signal_prior = model_signal_prior_dist(input)
    data_source_dict, model_dist_signal_dict = utils.data_source_dict(model, input)
    cfg_names_to_run = []
    for sp in signal_processes:
        main = Run(n, data_source_dict, model_signal_prior)
        mle = MleProducer(Distribution.merge(signal_prior, nuisance_constraint))
        main.add_producer(mle)
        name = write_cfg2(main, model, sp, 'ml_fit', input)
        cfg_names_to_run.append(name)
    run_theta_ = options.get('run_theta', True)
    if run_theta_: run_theta(cfg_names_to_run)
    else: return None

    cachedir = os.path.join(config.workdir, 'cache')    
    result = {}
    result_table = table()
    result_table.add_column('process', 'signal process')
    
    nuisance_parameters = sorted(list(model.get_parameters('', True)))
    for p in nuisance_parameters:
        suffix = ''
        if nuisance_constraint.get_distribution(p)['width'] == 0.0: suffix = ' (fixed)'
        result_table.add_column(p, '%s%s' % (p, suffix))
    suffix = ''
    if signal_prior_spec.startswith('fix:'): suffix = ' (fixed)'
    result_table.add_column('beta_signal', 'beta_signal%s' % suffix)
    for icfg in range(len(cfg_names_to_run)):
        sp = signal_processes[icfg]
        name = cfg_names_to_run[icfg]
        method, sp_id, dummy = name.split('-',2)
        result[sp_id] = {}
        result_table.set_column('process', sp_id)
        parameters = set(model.get_parameters(sp, True))
        sqlfile = os.path.join(cachedir, '%s.db' % name)
        cols = ['mle__%s, mle__%s_error' % (p, p) for p in parameters]
        data = sql(sqlfile, 'select %s from products' % ', '.join(cols))
        if len(data) == 0: raise RuntimeError, "no data in result file '%s'" % sqlfile
        i = 0
        for p in parameters:
            result[sp_id][p] = [(row[2*i], row[2*i+1]) for row in data]
            i += 1
            sorted_res = sorted([res[0] for res in result[sp_id][p]])
            n = len(sorted_res)
            if n >= 10:
                result_table.set_column(p, '%.3g (%.3g, %.3g)' % (sorted_res[int(0.5*n)], sorted_res[int(0.16*n)], sorted_res[int(0.84*n)]))
            else: result_table.set_column(p, '%.3g' % sorted_res[int(0.5*n)])
        for p in nuisance_parameters + ['beta_signal']:
            if p in parameters: continue
            result_table.set_column(p, 'n/a')
        result_table.add_row()
    config.report.new_section("Maximum Likelihood fit on ensemble '%s'" % input)
    config.report.add_p('The table entries give the median (and, if n>=10 the central 68%) of the parameter values at the found maximum of the likelihood function.')
    config.report.add_html(result_table.html())
    return result


       


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



## \brief Quantify the approximate impact of an unertainty on the result
#
# This will run
#    method(model, input = 'toys-asimov:1.0', n = 1, signal_processes = ..., nuisance_constraint = <see below>, **method_options)
# for each systematic uncertainty parameter twice, setting this parameter to +1sigma and -1sigma resp, 
# fixing all other parameters to their nominal value.
#
# It is possible to change this default of scanning -1sigma and +1sigma by setting sigma_factors.
#
# nuisance_constraint is set to a copy of model.distribution; if method_options contain
# 'nuisance_constraint', this will be considered as usual (i.e., you can set it to 'shape:fix', ...).
#
#
# returns a dictionary
# (spid) --> (parameter name) --> (sigma factor) --> (method result)
#
#
# Example: to see how the result changes if "switching on" the uncertainty in the toys and fixing the same parameter
# in the evaluation, you can use something along these line:
# for p in parameters:
#    fixp = Distribution()
#    fixp.set_distribution(p, 'gauss', 0.0, 0.0, [0., 0.])
#    individual_uncertainties(model, ..., nuisance_constraint = fixp, parameters = [p])
def individual_uncertainties(model, method, signal_processes = None, sigma_factors = [-1.0, 1.0], parameters = None,  **method_options):
    assert 'signal_processes' not in method_options
    if signal_processes is None: signal_processes = [[sp] for sp in model.signal_processes]
    dist_for_method = copy.deepcopy(model.distribution)
    if 'nuisance_constraint' in method_options:
        constr = nuisance_prior_distribution(model, method_options['nuisance_constraint'])
        dist_for_method = Distribution.merge(dist_for_method, constr)
        del method_options['nuisance_constraint']
    model_dist_orig = copy.deepcopy(model.distribution)
    result = {}
    if 'n' not in method_options: method_options['n'] = 1
    for sp in signal_processes:
        spid = ''.join(sp)
        result[spid] = {}
        if parameters is None: pars = model.get_parameters(sp)
        else: pars = parameters
        for par in pars:
            if par == 'beta_signal': continue
            result[spid][par] = {}
            # fix all parameters but par:
            model.distribution = copy.deepcopy(model_dist_orig)
            for par2 in model.get_parameters(sp):
                if par == par2 or par2 == 'beta_signal': continue
                model.distribution.set_distribution_parameters(par2, width=0.0)
            # find out mean and width of par:
            dist_pars = model.distribution.get_distribution(par)
            mean, width = dist_pars['mean'], dist_pars['width']
            for sf in sigma_factors:
                model.distribution.set_distribution_parameters(par, mean = mean + sf * width, width=0.0, range = (mean + sf * width, mean + sf * width))
                res = method(model, input = 'toys-asimov:1.0', nuisance_constraint = dist_for_method, signal_processes = [sp], **method_options)
                result[spid][par][sf] = res
    model.distribution = model_dist_orig
    return result


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
    
## \brief Perform a maximum likelihood fit and get the coefficient function values for all processes / channels
#
#
# options: see ml_fit
#
# returns a dictionary
# (signal process id) --> (observable name) --> (process name) --> (factor)
#
# Does not write anything to the report.
def ml_fit_coefficients(model, **options):
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
    sqlfile = os.path.join(cachedir, '%s.db' % name_data)
    data = sql(sqlfile, 'select mle__ks_ts from products')
    if len(data) == 0: raise RuntimeError, "no data in result file '%s'" % sqlfile
    ks_value_data = data[0][0]
    p = len([x for x in ks_values_bkg if x >= ks_value_data]) * 1.0 / len(ks_values_bkg)
    return p
    

# make a mle fit to data and compare the obtained chi2 value with the ones generated by toys according to input.
def chi2_test(model, signal_process_group, input = 'toys:1.0', n = 5000, signal_prior = 'flat', nuisance_constraint = '', **options):
    nuisance_constraint = nuisance_prior_distribution(model, nuisance_constraint)
    main = {'n-events': n, 'model': '@model', 'producers': ('@mle',), 'output_database': sqlite_database(), 'log-report': False}
    mle = {'type': 'mle', 'name': 'mle', 'parameters': [], 'write_pchi2': True,
       'override-parameter-distribution': product_distribution("@nuisance_constraint", "@signal_prior"),
       'minimizer': minimizer(need_error = False)}
    signal_prior = signal_prior_dict(signal_prior)
    cfg_options = {'plugin_files': ('$THETA_DIR/lib/core-plugins.so','$THETA_DIR/lib/root.so')}
    toplevel_settings = {'mle': mle, 'main': main, 'options': cfg_options, 'signal_prior': signal_prior}
    cfg_names_to_run = []
    model_parameters = sorted(list(model.get_parameters(signal_process_group)))
    toplevel_settings['nuisance_constraint'] = nuisance_constraint.get_cfg(model_parameters)
    # for toys:
    main['data_source'], toplevel_settings['model-distribution-signal'] = data_source_dict(model, input)
    name = write_cfg(model, signal_process_group, 'chi2_test', input, additional_settings = toplevel_settings)
    cfg_names_to_run.append(name)
    # for data:
    input = 'data'
    main['data_source'], toplevel_settings['model-distribution-signal'] = data_source_dict(model, input)
    main['n-events'] = 1
    name = write_cfg(model, signal_process_group, 'chi2_test', input, additional_settings = toplevel_settings)
    cfg_names_to_run.append(name)
    if 'run_theta' not in options or options['run_theta']:
        run_theta(cfg_names_to_run)
    else: return None
    cachedir = os.path.join(config.workdir, 'cache')

    name_toys = cfg_names_to_run[0]
    sqlfile = os.path.join(cachedir, '%s.db' % name_toys)
    data = sql(sqlfile, 'select mle__pchi2 from products')
    if len(data) == 0: raise RuntimeError, "no data in result file '%s'" % sqlfile
    chi2_values_bkg = [row[0] for row in data]
    
    name_data = cfg_names_to_run[1]
    sqlfile = os.path.join(cachedir, '%s.db' % name_data)
    data = sql(sqlfile, 'select mle__pchi2 from products')
    if len(data) == 0: raise RuntimeError, "no data in result file '%s'" % sqlfile
    chi2_value_data = data[0][0]
    print "chi2 in data: ", chi2_value_data
    p = len([x for x in chi2_values_bkg if x >= chi2_value_data]) * 1.0 / len(chi2_values_bkg)
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
def pl_intervals(model, input = 'toys:0', n = 100, signal_prior = 'flat', nuisance_constraint = '', cls = [0.90], signal_processes = None, write_report = None, **options):
    signal_processes = default_signal_processes(model, signal_processes)
    if write_report is None: write_report = input == 'data'
    nuisance_constraint = nuisance_prior_distribution(model, nuisance_constraint)
    signal_prior_spec = signal_prior
    signal_prior = signal_prior_dist(signal_prior)
    model_signal_prior = model_signal_prior_dist(input)
    data_source_dict, model_dist_signal_dict = utils.data_source_dict(model, input)
    cfg_names_to_run = []
    for sp in signal_processes:
        main = Run(n, data_source_dict, model_signal_prior)
        pl = PliProducer(Distribution.merge(signal_prior, nuisance_constraint), cls)
        main.add_producer(pl)
        name = write_cfg2(main, model, sp, 'pl_intervals', input)
        cfg_names_to_run.append(name)
    run_theta_ = options.get('run_theta', True)
    if run_theta_: run_theta(cfg_names_to_run)
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
            colnames.append('pli__lower%s' % cs)
            colnames.append('pli__upper%s' % cs)
        data = sql(sqlfile, 'select pli__maxl, %s from products' % (', '.join(colnames)))
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
    tmpdir = False
    for arg in sys.argv[1:]:
        if '.py' in arg: scriptname = arg
        if arg=='--tmpdir': tmpdir = True
    if tmpdir:
        config.workdir = tempfile.mkdtemp()
    else:
        config.workdir = os.path.join(os.getcwd(), scriptname[:-3])
        config.workdir = os.path.realpath(config.workdir)    
    setup_workdir()
    config.theta_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), '..'))
    config.report = html_report(os.path.join(config.workdir, 'index.html'))
    variables = globals()
    variables['report'] = config.report
    utils.info("executing script %s" % scriptname)
    try:
        execfile(scriptname, variables)
    except Exception as e:
        print "error while trying to execute analysis script %s:" % scriptname
        traceback.print_exc()
        sys.exit(1)
    utils.info("workdir is %s" % config.workdir)
        
if __name__ == '__main__': main()



