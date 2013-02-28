# -*- coding: utf-8 -*-

import hashlib, io, os, time, re, numpy, shutil, termios, threading
from Model import *
from utils import *
import config

# In general, each (plugin) class in theta (which in turn corresponds to a setting group with a "type" setting) corresponds to a class here.
#
# Conventions:
# * each class representing a concrete theta plugin should implement a constructor which takes all required arguments to build the class completely.
# * All classes are immutable in the sense that no change should be done after they are constructed. In case manipulation is needed,
#   use methods which return modified copies.
# * Note that any model instances, etc. passed to the modules here is assumed to be unmodified until theta is run.
#   This allows the modules to just save references instead of performing copies.
# * the theta configuration should be returned as dictionary by a method get_cfg(options) where options is an instance of Options class



# Conventions for theta config file generation:
# * there is only one model relevant to a config file, which is written to the path 'model', s.t. other modules
#   can refer to it via "@model"
# * the default nusiance parameter distribution (for toy generation and likelihood construction) is the one
#   from the model (model.distribution), but it can be overridden (even for a subset of parameters) for both toy generation and for
#   the likelihood construction in the producers.


from ConfigParser import SafeConfigParser


# Options to control certain details of the theta configuration generation which do not change the actual modules run etc. but generate
# somewhat different theta configurations. Examples are
# * number of threads to use
# * default minimizer options, such as whether or not to use mcmc_minimizer and how
# * whether or not to use llvm
# * some fine-tuning options for cls_limits, such as the number of toys for background-only and clb_cutoff, etc.
#
# The reason not to implement these as options for the methods are that the Options instance can be passed around more easily,
# even if new options are added, or old ones removed.
#
# Options affecting the actual statistical / physical _meaning_ of what should be done (and not
# just some details of how the theta config is generated)
# should NOT go here. Rather, they should be standard arguments for the respective methods.
#
class Options(SafeConfigParser):
    def __init__(self):
        SafeConfigParser.__init__(self)
        self.default_config = """
# global options concerning all / multiple modules or theta itself
[global]
debug = False
check_cache = True
        
# minimizer options
[minimizer]
always_mcmc = False
mcmc_iterations = 1000
bootstrap_mcmcpars = 0
# strategy can be 'fast' or 'robust'
strategy = fast
minuit_tolerance_factor = 1
       
[cls_limits]
clb_cutoff = 0.01
create_debuglog = True
        
[model]
use_llvm = False
        
[main]
n_threads = 1

[mcmc]
strategy = asimov_widths
stepsize_factor = None
"""
        self.readfp(io.BytesIO(self.default_config))
        
    def get_workdir(self):
        return os.path.realpath(config.workdir)

# each class representing a theta module should inherit from ModuleBase
class ModuleBase:
    def __init__(self):
        self.submodules = []
        self.required_plugins = frozenset(['core-plugins.so'])
    
    # derived classes should usually set self.required_plugins instead of re-implementing this function.
    # One exception is the case in which the plugins required depend on options.
    def get_required_plugins(self, options):
        plugins = set(self.required_plugins)
        for m in self.submodules:
            plugins.update(m.get_required_plugins(options))
        return plugins
        
    # declare module to be a submodule of the current module. This is used to track plugin
    # dependencies.
    def add_submodule(self, module):
        self.submodules.append(module)


class DataSource(ModuleBase):
    # input_string is either 'toys-asimov:XXX', 'toys:XXX', (where 'XXX' is a floating point value, the value of beta_signal), 'data',
    # or 'replay-toys:XXX' where XXX is the filename. It is also valid to directly specify the filename of the .db file.
    #
    # override_distribution is a Distribution which can override the default choice in the 'toys' case. If None, the default is to
    # * use model.distribution in case of input_string = 'toys:XXX'
    # * use the model.distrribution fixed to the most probable parameter values in case of input_string = 'toys-asimov:XXX'
    def __init__(self, input_string, model, signal_processes, override_distribution = None, name = 'source', seed = None):
        ModuleBase.__init__(self)
        flt_regex = r'([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)'
        if input_string.endswith('.db') and not input_string.startswith('replay-toys'):
            input_string = 'replay-toys:' + input_string
        self.name = name
        self.seed = seed
        self._id = input_string
        self.signal_processes = signal_processes
        self.model = model
        if input_string.startswith('toys'):
            self.mode = 'toys'
            self.beta_signal = None
            self.asimov = False
            match_toys = re.match('toys:%s' % flt_regex, input_string)
            if match_toys is not None:
                self.beta_signal = float(match_toys.group(1))
                if override_distribution is None:
                    self.distribution = model.distribution
                else:
                    self.distribution = Distribution.merge(model.distribution, override_distribution)
            else:
                match_toys_asimov = re.match('toys-asimov:%s' % flt_regex, input_string)
                if match_toys_asimov is None: raise RuntimeError, "invalid input specification '%s'" % input_string
                self.asimov = True
                self.beta_signal = float(match_toys_asimov.group(1))
                if override_distribution is None:
                    self.distribution = get_fixed_dist(model.distribution)
                else:
                    self.distribution = Distribution.merge(get_fixed_dist(model.distribution), override_distribution)
        elif input_string.startswith('replay-toys:'):
            self.mode = 'replay-toys'
            self.fname = input_string[len('replay-toys:'):]
            if not os.path.exists(self.fname): raise RuntimeError, "specified '%s' as data source, but file '%s' does not exist!" % (input_string, self.fname)
            self._id = 'replay-toys:' + self.fname[self.fname.rfind('-')+1:self.fname.rfind('.db')]
        else:
            if input_string != 'data': raise RuntimeError, "invalid input specification '%s'" % input_string
            self.mode = 'data'
            self.data_histos = {}
            for o in model.get_observables():
                self.data_histos[o] = model.get_data_histogram(o)
            self.data_rvobsvalues = model.get_data_rvobsvalues()
            
            
    def get_id(self): return self._id
    
    def get_cfg(self, options):
        result = {'name': self.name}
        if self.mode == 'toys':
            seed = self.seed
            if seed is None: seed = -1
            parameters = self.model.get_parameters(self.signal_processes)
            dist = {'type': 'product_distribution', 'distributions': [self.distribution.get_cfg(parameters)]}
            if 'beta_signal' in parameters: dist['distributions'].append({'type': 'delta_distribution', 'beta_signal' : self.beta_signal})
            result.update({'type': 'model_source', 'model': '@model', 'override-parameter-distribution': dist, 'rnd_gen': {'seed': seed}})
            if self.asimov:
                result['dice_poisson'] = False
                result['dice_template_uncertainties'] = False
                result['dice_rvobs'] = False
            return result
        elif self.mode == 'replay-toys':
            result.update({'type': 'replay_toys', 'input_database': {'type': 'sqlite_database_in', 'filename': self.fname}})
            return result
        else:
            result['type'] = 'histo_source'
            for o in self.data_histos: result[o] = self.data_histos[o].get_cfg(False)
            if len(self.data_rvobsvalues) > 0:
                result['rvobs-values'] = dict(self.data_rvobsvalues)
            return result


# a data source which adds many root_ntuple_sources. This is suitable for drawing toy datasets which are random subsets of
# events of a TTree which in turn is mainly useful for correlation studies ...
#
# Usage::
#
#  source = RootNtupleSource(model)
#  # Usually, the branch name is assumed to be the same as the observable name. If this is not the case, you have to explicitly state that:
#  source.define_observable_branchname(obsname = 'nnout', branchname = 'discriminator')
#  # add some files. Usually, there is one file per 'process'. You can specify optional total_nevents and weight_branchname to use
#  source.add_file('ttbar.root', mean_nevents = 107.2, relweight_branchname = 'weight2')
#
# The default for mean_nevents is the sum of 
#
class RootNtupleSource(ModuleBase):
    def __init__(self, model, name = 'source'):
        ModuleBase.__init__(self)
        self.required_plugins = frozenset(['root.so', 'core-plugins.so'])
        self.name = name
        # each entry in files is a dictionary with the keys 'filename', 'mean_nevents', etc.
        self.files = []
        # the 'observables' configuration for the root_ntuple_source:
        self.observables_cfg = {}
        for obs in model.get_observables():
            xmin, xmax, nbins = model.get_range_nbins(obs)
            self.observables_cfg[obs] = {'branchname': obs, 'nbins': nbins, 'range': (xmin, xmax)}
        
        
    def define_observable_branchname(self, obsname, branchname):
        if obsname not in self.observables_cfg: raise RuntimeError, "unknown observable '%s'" % obsname
        self.observables_cfg[obsname]['branchname'] = branchname
    
    # using None will use the defaults from the root_ntuple_source plugin.
    def add_file(self, filename, mean_nevents = None, relweight_branchname = None, treename = None, seed = None):
        d = {'filename': filename}
        if mean_nevents is not None: d['mean_nevents'] = float(mean_nevents)
        if relweight_branchname is not None: d['relweight_branchname'] = str(relweight_branchname)
        if treename is not None: d['treename'] = str(treename)
        if seed is not None: d['seed'] = int(seed)
        
        
    def get_cfg(self, options):
        result = {'name': self.name, 'type': 'add_sources'}
        result['sources'] = []
        for i, f in enumerate(self.files):
            s_cfg = {'name': '%s%d' % (self.name, i), 'type': 'root_ntuple_source', 'observables': self.observables_cfg}
            s_cfg['filename'] = f['filename']
            for key in ('treename', 'mean_nevents', 'relweight_branchname'):
                if key in f: s_cfg[key] = f[key]
            if 'seed' in f: s_cfg['rnd_gen'] = {'seed': f['seed']}
            result.sources.append(s_cfg)

# base class for producers, providing consistent approach to override_parameter_distribution and name
class ProducerBase(ModuleBase):
    def get_cfg_base(self, options):
        result = {'name': self.name}
        if self.override_distribution_cfg is not None:
            result['override-parameter-distribution'] = self.override_distribution_cfg
        return result
    
    
    # note that signal_prior is only respected if override_distribution is not None.
    # if model is None, override_distribution and signal_prior are ignored.
    def __init__(self, model, signal_processes, override_distribution, name, signal_prior):
        ModuleBase.__init__(self)
        assert type(name) == str
        self.override_distribution_cfg = None
        self.name = name
        if model is not None:
            signal_prior_cfg = _signal_prior_dict(signal_prior)
            parameters = set(model.get_parameters(signal_processes))
            if override_distribution is not None: dist = Distribution.merge(model.distribution, override_distribution)
            else: dist = model.distribution
            if 'beta_signal' in parameters:
                self.override_distribution_cfg = {'type': 'product_distribution', 'distributions': [dist.get_cfg(parameters), signal_prior_cfg]}
            else:
                self.override_distribution_cfg = dist.get_cfg(parameters)
    

class MleProducer(ProducerBase):
    # parameters_write is the list of parameters to write the mle for in the db. The default (None) means
    # to use all parameters.
    def __init__(self, model, signal_processes, override_distribution, signal_prior = 'flat', name = 'mle', need_error = True, with_covariance = False, parameters_write = None, ks = False, chi2 = False):
        ProducerBase.__init__(self, model, signal_processes, override_distribution, name, signal_prior)
        self.minimizer = Minimizer(need_error or with_covariance)
        self.add_submodule(self.minimizer)
        self.ks = ks
        self.chi2 = chi2
        self.with_covariance = with_covariance
        if parameters_write is None: self.parameters_write = sorted(list(model.get_parameters(signal_processes)))
        else: self.parameters_write = sorted(list(parameters_write))
        
    def get_cfg(self, options):
        result = {'type': 'mle', 'minimizer': self.minimizer.get_cfg(options), 'parameters': list(self.parameters_write)}
        if self.with_covariance: result['write_covariance'] = True
        if self.ks: result['write_ks_ts'] = True
        if self.chi2: result['write_pchi2'] = True
        result.update(self.get_cfg_base(options))
        return result
        
        
class QuantilesProducer(ProducerBase):
    def __init__(self, model, signal_processes, override_distribution, signal_prior = 'flat', name = 'quant', parameter = 'beta_signal', quantiles = [0.16, 0.5, 0.84], iterations = 10000, seed = 0):
        ProducerBase.__init__(self, model, signal_processes, override_distribution, name, signal_prior)
        self.parameter = parameter
        self.quantiles = quantiles
        self.iterations = iterations
        self.seed = seed
        
    def get_cfg(self, options):
        strategy = options.get('mcmc', 'strategy')
        stepsize_factor = options.get('mcmc', 'stepsize_factor')
        result = {'type': 'mcmc_quantiles', 'parameter': self.parameter, 'quantiles': self.quantiles}
        result.update(self.get_cfg_base(options))
        result['mcmc_strategy'] = {'type': strategy, 'name': self.name + "_mcs", 'iterations': self.iterations}
        if stepsize_factor != 'None':  result['mcmc_strategy']['factor'] = float(stepsize_factor)
        if self.seed != 0: result['mcmc_strategy']['rnd_gen'] = {'seed': self.seed}
        return result
        
        
        
class PosteriorProducer(ProducerBase):
    def __init__(self, model, signal_processes, override_distribution, histogram_specs, signal_prior = 'flat', name = 'post', iterations = 10000, smooth = True):
        ProducerBase.__init__(self, model, signal_processes, override_distribution, name, signal_prior)
        self.histogram_specs = histogram_specs
        self.iterations = iterations
        self.smooth = smooth
        
    def get_cfg(self, options):
        result = {'type': 'mcmc_posterior_histo', 'iterations': self.iterations, 'smooth': self.smooth, 'parameters': self.histogram_specs.keys()}
        for (p, (nbins, xmin, xmax)) in self.histogram_specs.iteritems():
            result['histo_%s' % p] = {'range': [xmin, xmax], 'nbins': nbins}
        result.update(self.get_cfg_base(options))
        return result


class PDWriter(ProducerBase):
    def __init__(self, name = 'pdw'):
        ProducerBase.__init__(self, model = None, signal_processes = None, override_distribution = None, name = name, signal_prior = None)
        
    def get_cfg(self, options):
        result = self.get_cfg_base(options)
        result['type'] = 'pseudodata_writer'
        result['write-data'] = True
        return result


class PliProducer(ProducerBase):
    def __init__(self, model, signal_processes, override_distribution = None, cls = [0.6827, 0.95], name = 'pli', parameter = 'beta_signal', signal_prior = None):
        ProducerBase.__init__(self, model, signal_processes, override_distribution, name, signal_prior)
        self.cls = [float(s) for s in cls]
        self.minimizer = Minimizer(need_error = False)
        self.add_submodule(self.minimizer)
        self.parameter = parameter

    def get_cfg(self, options):
        result = {'type': 'deltanll_intervals', 'minimizer': self.minimizer.get_cfg(options), 'parameter': self.parameter, 'clevels': self.cls}
        result.update(self.get_cfg_base(options))
        return result
        
        
class DeltaNllHypotest(ProducerBase):
    def __init__(self, model, signal_processes, override_distribution = None, name = 'dnll', restrict_poi = None, restrict_poi_value = None, signal_prior_sb = 'flat', signal_prior_b = 'fix:0.0', write_pchi2 = False):
        ProducerBase.__init__(self, model, signal_processes, override_distribution = None, name = name, signal_prior = None)
        self.minimizer = Minimizer(need_error = False)
        self.restrict_poi = restrict_poi
        self.write_pchi2 = write_pchi2
        self.restrict_poi_value = restrict_poi_value
        self.add_submodule(self.minimizer)
        if override_distribution is not None: dist = Distribution.merge(model.distribution, override_distribution)
        else: dist = model.distribution
        sb_parameters = set(model.get_parameters(signal_processes))
        b_parameters = set(model.get_parameters(''))
        means = dist.get_means()
        # the "signal parameters": those which the model depends only for the signal ...
        means_spar = {}
        for p in means:
            if p in sb_parameters and p not in b_parameters: means_spar[p] = means[p]
        dist_bkg = Distribution.merge(dist, get_fixed_dist_at_values(means_spar))
        self.sb_distribution_cfg = {'type': 'product_distribution', 'distributions': [dist.get_cfg(sb_parameters), _signal_prior_dict(signal_prior_sb)]}
        self.b_distribution_cfg = {'type': 'product_distribution', 'distributions': [dist_bkg.get_cfg(sb_parameters), _signal_prior_dict(signal_prior_b)]}
        
        
    def get_cfg(self, options):
        result = {'type': 'deltanll_hypotest', 'minimizer': self.minimizer.get_cfg(options),
           'background-only-distribution': self.b_distribution_cfg, 'signal-plus-background-distribution': self.sb_distribution_cfg}
        if self.write_pchi2: result['write_pchi2'] = True
        result.update(self.get_cfg_base(options))
        if 'override-parameter-distribution' in result: del result['override-parameter-distribution']
        if self.restrict_poi is not None: result['restrict_poi'] = self.restrict_poi
        if self.restrict_poi_value is not None: result['default_poi_value'] = self.restrict_poi_value
        return result
        
class MCMCRatioProducer(ProducerBase):
    def __init__(self, model, signal_processes, override_distribution = None, name = 'mcmcratio', signal_prior_sb = 'fix:1.0', signal_prior_b = 'fix:0.0', iterations = 10000):
        ProducerBase.__init__(self, model, signal_processes, override_distribution = None, name = name, signal_prior = None)
        if override_distribution is not None: dist = Distribution.merge(model.distribution, override_distribution)
        else: dist = model.distribution
        sb_parameters = set(model.get_parameters(signal_processes))
        b_parameters = set(model.get_parameters(''))
        means = dist.get_means()
        # the "signal parameters": those which the model depends only for the signal ...
        means_spar = {}
        for p in means:
            if p in sb_parameters and p not in b_parameters: means_spar[p] = means[p]
        dist_bkg = Distribution.merge(dist, get_fixed_dist_at_values(means_spar))
        self.sb_distribution_cfg = {'type': 'product_distribution', 'distributions': [dist.get_cfg(sb_parameters), _signal_prior_dict(signal_prior_sb)]}
        self.b_distribution_cfg = {'type': 'product_distribution', 'distributions': [dist_bkg.get_cfg(sb_parameters), _signal_prior_dict(signal_prior_b)]}
        self.iterations = iterations        
        
    def get_cfg(self, options):
        result = {'type': 'mcmc_posterior_ratio', 'background-only-distribution': self.b_distribution_cfg, 'signal-plus-background-distribution': self.sb_distribution_cfg, 'iterations': self.iterations}
        result.update(self.get_cfg_base(options))
        if 'override-parameter-distribution' in result: del result['override-parameter-distribution']
        return result
       
class NllScanProducer(ProducerBase):
    def __init__(self, model, signal_processes, override_distribution = None, signal_prior = 'flat', name = 'nllscan', parameter = 'beta_signal', range = [0.0, 3.0], npoints = 101, adaptive_startvalues = False):
        ProducerBase.__init__(self, model, signal_processes, override_distribution = override_distribution, name = name, signal_prior = signal_prior)
        self.parameter = parameter
        self.range = range
        self.npoints = npoints
        self.adaptive_startvalues = adaptive_startvalues
        self.minimizer = Minimizer(need_error = False)
        self.add_submodule(self.minimizer)

    def get_cfg(self, options):
        result = {'type': 'nll_scan', 'minimizer': self.minimizer.get_cfg(options), 'parameter': self.parameter,
               'parameter-values': {'start': self.range[0], 'stop': self.range[1], 'n-steps': self.npoints}, 'adaptive_startvalues': self.adaptive_startvalues}
        result.update(self.get_cfg_base(options))
        return result


class Minimizer(ModuleBase):
    def __init__(self, need_error):
        ModuleBase.__init__(self)
        self.need_error = need_error
        
    def get_required_plugins(self, options):
        plugins = set(['core-plugins.so'])
        strategy = options.get('minimizer', 'strategy')
        if strategy == 'simplex_vanilla':  plugins.add('simplex.so')
        elif strategy == 'lbfgs_vanilla': plugins.add('liblbfgs.so')
        else: plugins.add('root.so')
        return plugins
        
    
    def get_cfg(self, options):
        always_mcmc = options.getboolean('minimizer', 'always_mcmc')
        mcmc_iterations = options.getint('minimizer', 'mcmc_iterations')
        strategy = options.get('minimizer', 'strategy')
        tolerance_factor = options.getfloat('minimizer', 'minuit_tolerance_factor')
        assert strategy in ('fast', 'robust', 'minuit_vanilla', 'simplex_vanilla', 'lbfgs_vanilla')
        if strategy == 'minuit_vanilla':
            result = {'type': 'root_minuit'}
        elif strategy == 'simplex_vanilla':
            result = {'type': 'simplex_minimizer'}
            #result['max_iter'] = '100000'
            #result = {'type': "minimizer_chain", 'minimizers': [{'type': "mcmc_minimizer", 'iterations': 1000, 'bootstrap_mcmcpars': 2, 'name': "mcmc_min0"}],
            #    'last_minimizer': result}
            return result
        elif strategy == 'lbfgs_vanilla':
            result = {'type': 'lbfgs_minimizer'}
            return result
        else:
            minimizers = []
            #try, in this order: migrad, mcmc+migrad, simplex, mcmc+simplex, more mcmc+simplex
            if not always_mcmc: minimizers.append({'type': 'root_minuit'})
            minimizers.append({'type': 'mcmc_minimizer', 'name':'mcmc_min0', 'iterations': mcmc_iterations, 'after_minimizer': {'type': 'root_minuit'}})
            if not always_mcmc: minimizers.append({'type': 'root_minuit', 'method': 'simplex'})
            minimizers.append({'type': 'mcmc_minimizer', 'name':'mcmc_min1', 'iterations': mcmc_iterations, 'after_minimizer': {'type': 'root_minuit', 'method': 'simplex'}})
            if strategy == 'robust':
                minimizers.append({'type': 'mcmc_minimizer', 'name':'mcmc_min2', 'iterations': 50 * mcmc_iterations, 'after_minimizer': {'type': 'root_minuit', 'method': 'simplex'}})
                minimizers.append({'type': 'mcmc_minimizer', 'name':'mcmc_min3', 'iterations': 500 * mcmc_iterations, 'after_minimizer': {'type': 'root_minuit', 'method': 'simplex'}})
            bootstrap_mcmcpars = options.getint('minimizer', 'bootstrap_mcmcpars')
            if bootstrap_mcmcpars > 0:
                for m in minimizers:
                    if m['type'] == 'mcmc_minimizer': m['bootstrap_mcmcpars'] = bootstrap_mcmcpars
            result = {'type': 'minimizer_chain', 'minimizers': minimizers}
            if self.need_error: result['last_minimizer'] = {'type': 'root_minuit'}
            
        def apply_tolfactor(d):
            if type(d) in (list, tuple):
                for v in d: apply_tolfactor(v)
                return
            if type(d)==dict:
                if 'type' in d and d['type'] == 'root_minuit': d['tolerance_factor'] = tolerance_factor
                for key, value in d.iteritems():
                    apply_tolfactor(value)
        if tolerance_factor != 1: apply_tolfactor(result)
        return result


def _settingvalue_to_cfg(value, indent=0, current_path = [], value_is_toplevel_dict = False):
    if type(value) == numpy.float64: value = float(value)
    tval = type(value)
    if tval == array.array:
        return "(" + ",".join(map(lambda f: "%.5e" % f, value)) + ")"
    if tval == str: return '"%s"' % value
    if tval == bool: return 'true' if value else 'false'
    if tval == int: return '%d' % value
    if tval == float:
        if value == inf: return '"inf"'
        elif value == -inf: return '"-inf"'
        return '%.5e' % value
    if tval == list or tval == tuple:
        return "(" + ",".join([_settingvalue_to_cfg(value[i], indent + 4, current_path + [str(i)]) for i in range(len(value))]) + ')'
    if tval == dict:
        result = ''
        if not value_is_toplevel_dict: result += "{\n"
        new_indent = (indent + 4) if not value_is_toplevel_dict else 0
        # sort keys to make config files reproducible:
        for s in sorted(value.keys()):
            result += ' ' * new_indent + s + " = " + _settingvalue_to_cfg(value[s], new_indent, current_path + [s]) + ";\n"
        if not value_is_toplevel_dict: result += ' ' * indent + "}"
        return result
    raise RuntimeError, "Cannot convert type %s to theta cfg in path '%s'" % (type(value), '.'.join(current_path))




# return a config dictionary for the signal prior. s can be a string
# such as "flat", "fix:XXX", etc., or a dictionary in which case
# it is returned unmodified
#
# None is equivalent to 'flat'
def _signal_prior_dict(spec):
    if spec is None: spec = 'flat'
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
            signal_prior_dict = {'type': 'flat_distribution', 'beta_signal': {'range': [xmin, xmax], 'fix-sample-value': value}}
        elif spec.startswith('fix:'):
            v = float(spec[4:])
            signal_prior_dict = {'type': 'delta_distribution', 'beta_signal': v}
        else: raise RuntimeError, "signal_prior specification '%s' unknown" % spec
    else:
        if type(spec) != dict: raise RuntimeError, "signal_prior specification has to be a string ('flat' / 'fix:X') or a dictionary!"
        signal_prior_dict = spec
    return signal_prior_dict


class DbResult(object):
    def __init__(self, theta_db_fname = None):
        self.theta_db_fname = theta_db_fname

    # get the content of a table in the output database of theta.
    def get_results(self, table_name, columns, order_by = None):
        if self.theta_db_fname is None: return None
        if type(columns)!=str:
            columns = '"' + '", "'.join(columns) + '"'
        else: assert columns == '*'
        if order_by is not None:
            order_by = ' order by "%s"' % order_by
        else: order_by = ''
        data, columns = sql_singlefile(self.theta_db_fname, 'select %s from "%s"%s' % (columns, table_name, order_by), True)
        result = {}
        for i in range(len(columns)):
            result[columns[i][0]] = [row[i] for row in data]
        return result
        
    def get_products(self, columns = '*', fail_if_empty = True):
        """
        get the result (in the products table) of this invocation of theta as dictionary::
        
           column_name --> list of values
           
        Which columns are returned can be controlled by the `columns` argument which is either a list of column names
        or as special case the string '*'. The default is '*' which will return all columns.
    
        If `fail_if_empty` is `True` and no rows are found, some information from the log table is printed and a `RuntimeError` is raised.
        This can happen e.g. in case there was one attempt to fit on data which failed due to a minimizer which failed.
        """
        res = self.get_results('products', columns)
        if fail_if_empty and len(res[res.keys()[0]])==0:
            self.print_logtable_errors()
            raise RuntimeError, "No result is available, see errors above"
        return res
        
    def print_logtable_errors(self, limit = 10):
        if not self.theta_db_fname: return
        data = sql_singlefile(self.theta_db_fname, 'select "eventid", "message" from "log" where "severity"=0 limit %d' % int(limit))
        if len(data)==0: return
        print "There have been errors for some toys: eventid, message"
        for i in range(len(data)):
            print data[i][0], data[i][1]
        
    def get_products_nevents(self):
        data = sql_singlefile(self.theta_db_fname, 'select count(*) from products')
        return data[0][0]
    
    def get_db_fname(self):
        """
        Return the path to the db file created by theta. Is `None` if not available
        """
        return self.theta_db_fname


class MainBase(ModuleBase, DbResult):
    """
    Base class corresponding to a theta `Main` plugin. One instance of this class corresponds to
    one theta .cfg file and one execution of theta, and one db result file. This class contains management
    code for config file creation, theta execution, and (via `DbResult`) reading the resulting db file. It
    also takes care of the theta-auto cache management.
    
    Derived classes must implement a method `get_cfg(self, options)`.
    """
    
    def __init__(self, fname_fragment):
        """
        `name` is used in the created cnofiguration file name
        """
        ModuleBase.__init__(self)
        DbResult.__init__(self)
        self.result_available = False
        # full path + filename of theta configuration file, if created:
        self.theta_cfg_fname = None
        self.thread = None
        self.cfg_fname_fragment = fname_fragment
        
    
    # return a dictionary with settings for the theta configuration for the plugins, parameters, rvobservables, and observables.
    def common_main_cfg(self, model, signal_processes, options):
        plugins = self.get_required_plugins(options)
        result = {'options': {'plugin_files': ['$THETA_DIR/lib/'+s for s in sorted(list(plugins))]}}
        result['parameters'] = sorted(list(model.get_parameters(signal_processes)))
        rvobservables = sorted(list(model.rvobs_distribution.get_parameters()))
        if len(rvobservables) > 0:
            result['rvobservables'] = rvobservables
        result['observables'] = {}
        for obs in sorted(list(model.get_observables())):
            xmin, xmax, nbins = model.observables[obs]
            result['observables'][obs] = {'range': [xmin, xmax], 'nbins': nbins}
        return result
    
    
    def _check_cache(self, options):
        """
        Check in the cache directory whether result is already available there
        returns True if it is, false otherwise.
        Note that get_result_available will return `True` if the result is in the cache
        only after calling check_cache.
        """
        workdir = options.get_workdir()
        cache_dir = os.path.join(workdir, 'cache')
        cfgfile_full = self.get_configfile(options)
        cfgfile_cache = os.path.join(cache_dir, os.path.basename(cfgfile_full))
        dbfile_cache = cfgfile_cache.replace('.cfg', '.db')
        if os.path.exists(cfgfile_cache) and os.path.exists(dbfile_cache) and  open(cfgfile_cache, 'r').read() == open(cfgfile_full, 'r').read():
            self.theta_db_fname = dbfile_cache
            self.result_available = True
            return True
        return False
    
    
    def get_configfile(self, options):
        """
        Returns the path to the theta config file, creating it if necessary.
        """
        if self.theta_cfg_fname is not None: return self.theta_cfg_fname
        cfg_dict = self.get_cfg(options)
        cfg_string = _settingvalue_to_cfg(cfg_dict, value_is_toplevel_dict = True)            
        output_name =  self.cfg_fname_fragment + '-' + hashlib.md5(cfg_string).hexdigest()[:10]
        if '@output_name' in cfg_string: cfg_string += "output_name = \"%s.db\";\n" % output_name
        if '@debuglog_name' in cfg_string: cfg_string += 'debuglog_name = "%s-debuglog.txt";\n' % output_name
        workdir = options.get_workdir()
        assert os.path.exists(workdir)
        self.theta_cfg_fname = os.path.join(workdir, output_name + '.cfg')
        f = open(self.theta_cfg_fname, 'w')
        f.write(cfg_string)
        f.close()
        return self.theta_cfg_fname
    

    # can be run in a different thread
    def _exec(self, cmd):
        self.error = None
        ret = os.system(cmd)
        if ret != 0:
            self.error = "executing theta ('%s') failed with exit code %d" % (cmd, ret)
            if os.isatty(1):
                attr = termios.tcgetattr(1)
                attr[3] |= termios.ECHO
                termios.tcsetattr(1, termios.TCSANOW, attr)
                        
    # should always be run in the main thread, after _exec() returns
    def _cleanup_exec(self):
        cfgfile = self.theta_cfg_fname
        assert cfgfile is not None
        dbfile = os.path.basename(cfgfile).replace('.cfg', '.db')
        if self.error is None:
            cache_dir = os.path.join(os.path.dirname(cfgfile), 'cache')
            if not os.path.exists(cache_dir): os.mkdir(cache_dir)
            cfgfile_cache = os.path.join(cache_dir, os.path.basename(cfgfile))
            dbfile_cache = cfgfile_cache.replace('.cfg', '.db')
            shutil.move(dbfile, dbfile_cache)
            shutil.copy(cfgfile, cfgfile_cache)
            self.result_available = True
            self.theta_db_fname = dbfile_cache
        else:
            if os.path.exists(dbfile) and not self.debug: os.unlink(dbfile)
            error = self.error
            self.error = None
            raise RuntimeError, error
            
    def run_theta(self, options, in_background_thread = False):
        """
        Execute the :program:`theta` program on the current .cfg file, creating the .cfg file if necessary.
        
        Parameters:
        
        * `options` - an instance of :class:`Options` which is passed to the submodules to control some details of the configuration
        * `in_background_thread` - if `True`, the :program:`theta` this methods returns immediately. Otherwise, theta is executed in a background thread.
        
        Executing multiple instances of the theta program in parallel is usually done like this::
           
           mains = [] # the list of MainBase instances
           ...
           for m in mains: m.run_theta(options, in_background_thread = True)
           for m in mains: m.wait_for_result_available()
        
        The first for-loop spawns the threads which execute theta, so theta is executed in parallel. The second loop waits for all theta executions
        to terminate (either normally or abnormally).
        
        In case :program:`theta` exists with a failure, a RuntimeException will be thrown. This is done by
        this method in case `in_background_thread = False`. In case of `in_background_thread = True`, the exception
        will be raised in :meth:`Run.wait_for_result_availble`.
        """
        check_cache = options.getboolean('global', 'check_cache')
        if check_cache: self._check_cache(options)
        cfgfile = self.get_configfile(options)
        self.debug = options.getboolean('global', 'debug')
        if self.result_available: return
        workdir = options.get_workdir()
        cache_dir = os.path.join(workdir, 'cache')
        theta = os.path.realpath(os.path.join(config.theta_dir, 'bin', 'theta'))
        info("Running 'theta %s'" % cfgfile)
        cmd = theta + " " + cfgfile
        if self.debug: cmd += " --print-time"
        to_execute = lambda : self._exec(cmd)
        if in_background_thread:
            assert self.thread is None
            self.thread = threading.Thread(target = to_execute)
            self.thread.start()
        else:
            to_execute()
            self._cleanup_exec()

    def get_result_available(self):
        """
        return whether the result is available, either after check_cache, or after run_theta in a background thread.
        
        Will call :meth:`Run.wait_for_result_available` in case theta was executed in a in the background thread; this
        can lead ot an exception, see :meth:`Run.wait_for_result_available`.
        """
        if self.thread is not None:
            if not self.thread.is_alive():
                # do the cleanup
                self.wait_for_result_available()
        return self.result_available
    
    
    def wait_for_result_available(self):
        """
        This method is only useful after calling run_theta(in_background_thread = True)
        
        can throw an exception in case the thread exited with an exception
        """
        if self.result_available: return
        if self.thread is None: return
        self.thread.join()
        self.thread = None
        self._cleanup_exec()
    
class Run(MainBase):
    def __init__(self, model, signal_processes, signal_prior, input, n, producers, nuisance_prior_toys = None, seed = None):
        """
        signal_prior can be a string or a dictionary.
        `input` is a string (toys:XXX, toy-asimov:XXX, ...), see DataSource for documentation
        """
        fragment = ''.join([p.name for p in producers]) + '-' + input + ''.join(signal_processes)
        MainBase.__init__(self, fragment)
        self.model = model
        self.signal_processes = signal_processes
        self.n_events = n
        self.producers = producers
        self.nuisance_prior_toys = nuisance_prior_toys
        for p in producers: self.add_submodule(p)
        self.seed = seed
        if input is None: self.data_source = None
        else:
            if type(input) == str:
                self.data_source = DataSource(input, model, signal_processes, override_distribution = nuisance_prior_toys, seed = seed)
            else:
                assert type(input) in (RootNtupleSource, DataSource), "unexpected type for 'input' argument: %s" % type(input)
                self.data_source = input
        self.signal_prior_cfg = _signal_prior_dict(signal_prior)
        
    def get_cfg(self, options):
        result = self.common_main_cfg(self.model, self.signal_processes, options)
        main = {'n-events': self.n_events, 'producers': [], 'log-report': False, 'output_database': {'type': 'sqlite_database', 'filename': '@output_name'},
            'model': "@model", 'data_source': self.data_source.get_cfg(options)}
        n_threads = options.getint('main', 'n_threads')
        if n_threads > 1:
            main['type'] = 'run_mt'
            main['n_threads'] = n_threads
        for p in self.producers:
            result['p%s' % p.name] = p.get_cfg(options)
            main['producers'].append('@p%s' % p.name)
        result.update({'main': main, 'model': self.model.get_cfg(self.signal_processes, self.signal_prior_cfg, options)})
        use_llvm = options.getboolean('model', 'use_llvm')
        if use_llvm:
            print "using llvm. This is EXPERIMENTAL. Use at your own risk"
            result['model']['type'] = 'llvm_model'
            result['options']['plugin_files'].append('$THETA_DIR/lib/llvm-plugins.so')
        return result

class ClsMain(MainBase):
    def __init__(self, model, signal_processes, signal_prior, input, producers, nuisance_prior_toys = None, seed = None):
        """
        Note that `input` is used for the "observed" cls limit; only the first (pseudo-)dataset of this DataSource will be used.
        """
        MainBase.__init__(self, 'cls')
        assert len(producers)==1
        self.model = model
        self.signal_processes = signal_processes
        self.producers = producers
        self.nuisance_prior_toys = nuisance_prior_toys
        for p in producers: self.add_submodule(p)
        self.seed = seed
        if input is None: self.data_source = None
        else:
            if type(input) == str:
                self.data_source = DataSource(input, model, signal_processes, override_distribution = nuisance_prior_toys, seed = seed)
            else:
                assert type(input) in (RootNtupleSource, DataSource), "unexpected type for 'input' argument: %s" % type(input)
                self.data_source = input
        self.signal_prior_cfg = _signal_prior_dict(signal_prior)
        self.cls_options = {'ts_column': self.producers[0].name + '__nll_diff'}
        self.minimizer = Minimizer(need_error = True)
        self.add_submodule(self.minimizer)
        
    def set_cls_options(self, **args):
        allowed_keys = ['ts_column', 'expected_bands', 'clb_cutoff', 'tol_cls', 'truth_max', 'reltol_limit', 'frequentist_bootstrapping', 'write_debuglog', 'input_expected']
        for k in allowed_keys:
            if k in args:
                self.cls_options[k] = args[k]
                del args[k]
        if len(args) > 0: raise RuntimeError, "unrecognized cls options: %s" % str(args.keys())
    
    def get_cfg(self, options):
        assert 'ts_column' in self.cls_options
        result = self.common_main_cfg(self.model, self.signal_processes, options)
        seed = self.seed
        if seed is None: seed = -1
        main = {'type': 'cls_limits', 'producer': self.producers[0].get_cfg(options), 'output_database': {'type': 'sqlite_database', 'filename': '@output_name'},
            'truth_parameter': 'beta_signal', 'tol_cls': self.cls_options.get('tol_cls', 0.025),
            'clb_cutoff': self.cls_options.get('clb_cutoff', 0.02), 'model': '@model',
            'debuglog': '@debuglog_name', 'rnd_gen': {'seed': seed }, 'ts_column': self.cls_options['ts_column'],
            'minimizer': self.minimizer.get_cfg(options), 'expected_bands': self.cls_options.get('expected_bands', 2000)}
        if not self.cls_options.get('write_debuglog', True): del main['debuglog']
        if self.cls_options.get('frequentist_bootstrapping', False):  main['nuisancevalues-for-toys'] = 'datafit'
        if 'truth_max' in self.cls_options: main['truth_max'] = float(self.cls_options['truth_max'])
        if 'reltol_limit' in self.cls_options: main['reltol_limit'] = float(self.cls_options['reltol_limit'])
        if self.data_source is not None: main['data_source'] = self.data_source.get_cfg(options)
        if 'input_expected' in self.cls_options:
            print "using input = %s for expected limits" % self.cls_options['input_expected']
            ds = DataSource(self.cls_options['input_expected'], self.model, self.signal_processes,
                override_distribution = self.nuisance_prior_toys)
            main['data_source_expected'] = ds.get_cfg(options)
        result.update({'main': main, 'model': self.model.get_cfg(self.signal_processes, self.signal_prior_cfg, options)})
        return result
    
    
class AsymptoticClsMain(MainBase):
    def __init__(self, model, signal_processes, input, beta_signal_expected = 0.0):
        """
        `input` should be either 'data' or None.
        """
        MainBase.__init__(self, 'acls')
        self.minimizer = Minimizer(need_error = False)
        self.add_submodule(self.minimizer)
        self.model = model
        self.signal_processes = signal_processes
        self.beta_signal_expected = beta_signal_expected
        if input is None: self.input = None
        else:
            if input != 'data': raise RuntimeError, "input should be either 'data' or None, but is %s" % str(input)
            self.input = DataSource(input, self.model, self.signal_processes)
        self.signal_prior_cfg = _signal_prior_dict('flat')
        
    def get_cfg(self, options):
        result = self.common_main_cfg(self.model, self.signal_processes, options)
        main = {'type': 'asymptotic_cls', 'model': '@model', 'parameter': 'beta_signal', 'minimizer': self.minimizer.get_cfg(options),
            'output_database': {'type': 'sqlite_database', 'filename': '@output_name'}}
        if self.input is not None:
            main['data'] = self.input.get_cfg(options)
        if self.beta_signal_expected != 0.0:
            main['parameter_value_expected'] = self.beta_signal_expected
        result.update({'main': main, 'model': self.model.get_cfg(self.signal_processes, self.signal_prior_cfg, options)})
        return result


# create a Histogram instance from the blob data in a sqlite db
def histogram_from_dbblob(blob_data):
    a = array.array('d')
    a.fromstring(blob_data)
    return Histogram(xmin = a[0], xmax = a[1], values = a[3:-1])
    
def matrix_from_dbblob(blob_data):
    n = int(math.sqrt(len(blob_data) / 8 - 4) + 0.4)
    assert (n*n + 4) * 8 == len(blob_data)
    return numpy.ndarray(shape = (n, n), buffer = blob_data[3*8:-8])
