# -*- coding: utf-8 -*-

import hashlib, io, os, time, re, numpy, shutil
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
# workdir either relative to the current directory or an absolute path.
workdir = analysis
        
# minimizer options
[minimizer]
always_mcmc = True
mcmc_iterations = 1000
        
[data_source]
seed = 1
        
[cls_limits]
clb_cutoff = 0.01
create_debuglog = True
        
[model]
use_llvm = False
        
[main]
n_threads = 1
"""
        self.readfp(io.BytesIO(self.default_config))
        
    def get_workdir(self):
        workdir_rel = self.get('global', 'workdir')
        return os.path.realpath(workdir_rel)

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
        
    # decalre module to be a submodule of the current module. This is used to track plugin
    # dependencies.
    def add_submodule(self, module):
        self.submodules.append(module)


class DataSource(ModuleBase):
    # input_string is either 'toys-asimov:XXX', 'toys:XXX', (where 'XXX' is a floating point value, the value of beta_signal), 'data',
    # or 'replay-toys:XXX' where XXX is the filename.
    #
    # override_distribution is a Distribution which can override the default choice in the 'toys' case. If None, the default is to
    # * use model.distribution in case of input_string = 'toys:XXX'
    # * use the model.distrribution fixed to the most probable parameter values in case of input_string = 'toys-asimov:XXX'
    def __init__(self, input_string, model, signal_processes, override_distribution = None, name = 'source'):
        ModuleBase.__init__(self)
        flt_regex = r'([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)'
        self.name = name
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
            seed = options.getint('data_source', 'seed')
            parameters = self.model.get_parameters(self.signal_processes, True)
            dist = {'type': 'product_distribution', 'distributions': [self.distribution.get_cfg(parameters)]}
            if 'beta_signal' in parameters: dist['distributions'].append({'type': 'delta_distribution', 'beta_signal' : self.beta_signal})
            result.update({'type': 'model_source', 'model': '@model', 'override-parameter-distribution': dist, 'rnd_gen': {'seed': seed}})
            if self.asimov:
                result['dice_poisson'] = False
                result['dice_template_uncertainties'] = False
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


# base class for producers, providing consistent approach to override_parameter_distribution
# and additional_nll_term
class ProducerBase(ModuleBase):
    def get_cfg_base(self, options):
        result = {'name': self.name}
        if self.override_distribution_cfg is not None:
            result['override-parameter-distribution'] = self.override_distribution_cfg
        if self.additional_nll_cfg is not None:
            result['additional-nll-term'] = self.additional_nll_cfg
        return result
    
    
    # note that signal_prior is only respected is override_distribution
    # if model is None, override_distribution and signal_prior is ignored.
    def __init__(self, model, signal_processes, override_distribution, name, signal_prior):
        ModuleBase.__init__(self)
        assert type(name) == str
        self.override_distribution_cfg, self.additional_nll_cfg = None, None
        self.name = name
        if model is not None:
            signal_prior_cfg = _signal_prior_dict(signal_prior)
            parameters = set(model.get_parameters(signal_processes, True))
            if override_distribution is not None: dist = Distribution.merge(model.distribution, override_distribution)
            else: dist = model.distribution
            self.override_distribution_cfg = {'type': 'product_distribution', 'distributions': [dist.get_cfg(parameters), signal_prior_cfg]}
            if model.additional_nll_term is not None:
                self.additional_nll_cfg = model.additional_nll_term.get_cfg()
        
    

class MleProducer(ProducerBase):
    # parameters_write is the list of parameters to write the mle for in the db. The default (None) means
    # to use all parameters.
    def __init__(self, model, signal_processes, override_distribution, signal_prior = 'flat', name = 'mle', need_error = True, parameters_write = None):
        ProducerBase.__init__(self, model, signal_processes, override_distribution, name, signal_prior)
        self.minimizer = Minimizer(need_error)
        self.add_submodule(self.minimizer)
        if parameters_write is None: self.parameters_write = sorted(list(model.get_parameters(signal_processes, True)))
        else: self.parameters_write = sorted(list(parameters_write))
        
    def get_cfg(self, options):
        result = {'type': 'mle', 'minimizer': self.minimizer.get_cfg(options), 'parameters': list(self.parameters_write)}
        result.update(self.get_cfg_base(options))
        return result
        
        
class QuantilesProducer(ProducerBase):
    def __init__(self, model, signal_processes, override_distribution, signal_prior = 'flat', name = 'quant', parameter = 'beta_signal', quantiles = [0.16, 0.5, 0.84],
    iterations = 10000):
        ProducerBase.__init__(self, model, signal_processes, override_distribution, name, signal_prior)
        self.parameter = parameter
        self.quantiles = quantiles
        self.iterations = iterations
        
    def get_cfg(self, options):
        result = {'type': 'mcmc_quantiles', 'parameter': self.parameter, 'quantiles': self.quantiles, 'iterations': self.iterations}
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
    def __init__(self, model, signal_processes, override_distribution = None, name = 'dnll'):
        ProducerBase.__init__(self, model, signal_processes, override_distribution = None, name = name, signal_prior = None)
        self.minimizer = Minimizer(need_error = False)
        self.add_submodule(self.minimizer)
        if override_distribution is not None: dist = Distribution.merge(model.distribution, override_distribution)
        else: dist = model.distribution
        parameters = set(model.get_parameters(signal_processes, True))
        self.sb_distribution_cfg = {'type': 'product_distribution', 'distributions': [dist.get_cfg(parameters), _signal_prior_dict('flat')]}
        self.b_distribution_cfg = {'type': 'product_distribution', 'distributions': [dist.get_cfg(parameters), _signal_prior_dict('fix:0.0')]}
        
        
    def get_cfg(self, options):
        result = {'type': 'deltanll_hypotest', 'minimizer': self.minimizer.get_cfg(options), 'background-only-distribution': self.b_distribution_cfg, 'signal-plus-background-distribution': self.sb_distribution_cfg}
        result.update(self.get_cfg_base(options))
        return result
       
    
       
"""        
class NllScanProducer(ProducerBase):
    def __init__(self, parameter_dist, parameter = 'beta_signal', name = 'nll_scan', range = [0.0, 3.0], npoints = 101):
        ProducerBase.__init__(self, parameter_dist, name)
        self.parameter = parameter
        self.range = range
        self.npoints = npoints

    def get_cfg(self, model, signal_processes, parameters_write = None):
        result = {'type': 'nll_scan', 'minimizer': minimizer(need_error = False), 'parameter': self.parameter,
               'parameter-values': {'start': self.range[0], 'stop': self.range[1], 'n-steps': self.npoints}}
        result.update(self.get_cfg_base(model, signal_processes))
        return result
"""

class Minimizer(ModuleBase):
    def __init__(self, need_error):
        ModuleBase.__init__(self)
        self.need_error = need_error
        self.required_plugins = frozenset(['core-plugins.so', 'root.so'])
    
    def get_cfg(self, options):
        always_mcmc = options.getboolean('minimizer', 'always_mcmc')
        mcmc_iterations = options.getint('minimizer', 'mcmc_iterations')
        minimizers = []
        #try, in this order: migrad, mcmc+migrad, simplex, mcmc+simplex, more mcmc+simplex
        if not always_mcmc: minimizers.append({'type': 'root_minuit'})
        minimizers.append({'type': 'mcmc_minimizer', 'name':'mcmc_min0', 'iterations': mcmc_iterations, 'after_minimizer': {'type': 'root_minuit'}})
        if not always_mcmc: minimizers.append({'type': 'root_minuit', 'method': 'simplex'})
        minimizers.append({'type': 'mcmc_minimizer', 'name':'mcmc_min1', 'iterations': mcmc_iterations, 'after_minimizer': {'type': 'root_minuit', 'method': 'simplex'}})
        minimizers.append({'type': 'mcmc_minimizer', 'name':'mcmc_min1', 'iterations': 50 * mcmc_iterations, 'after_minimizer': {'type': 'root_minuit', 'method': 'simplex'}})
        minimizers.append({'type': 'mcmc_minimizer', 'name':'mcmc_min1', 'iterations': 500 * mcmc_iterations, 'after_minimizer': {'type': 'root_minuit', 'method': 'simplex'}})
        result = {'type': 'minimizer_chain', 'minimizers': minimizers}
        if self.need_error: result['last_minimizer'] = {'type': 'root_minuit'}
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
            signal_prior_dict = theta_interface.delta_distribution(beta_signal = v)
        else: raise RuntimeError, "signal_prior specification '%s' unknown" % spec
    else:
        if type(spec) != dict: raise RuntimeError, "signal_prior specification has to be a string ('flat' / 'fix:X') or a dictionary!"
        signal_prior_dict = spec
    return signal_prior_dict

        
# get_cfg returns a toplevel configuration
#
# corresponds to the "run" class in theta:
class Run:
    
    # signal_prior can be a string or a dictionary.
    # input is a string (toys:XXX, toy-asimov:XXX, ...), see DataSource for documentation
    def __init__(self, model, signal_processes, signal_prior, input, n, producers, nuisance_prior_toys = None):
        self.model = model
        self.signal_processes = signal_processes
        self.n_events = n
        self.producers = producers
        self.data_source = DataSource(input, model, signal_processes, override_distribution = nuisance_prior_toys)
        self.input_id = self.data_source.get_id()
        self.signal_prior_cfg = _signal_prior_dict(signal_prior)
        self.result_available = False
        # full path + filename of theta configuration file, if created:
        self.theta_cfg_fname = None
        self.theta_db_fname = None
        self.thread = None
    
    def _get_cfg(self, options):
        plugins = set()
        for p in self.producers:
            plugins.update(p.get_required_plugins(options))
        main = {'n-events': self.n_events, 'producers': [], 'log-report': False, 'output_database': {'type': 'sqlite_database', 'filename': '@output_name'},
                'model': "@model", 'data_source': self.data_source.get_cfg(options)}
        n_threads = options.getint('main', 'n_threads')
        if n_threads > 1:
            main['type'] = 'run_mt'
            main['n_threads'] = n_threads
        result = {'main': main, 'options': {'plugin_files': ['$THETA_DIR/lib/'+s for s in sorted(list(plugins))]},
                  'model': self.model.get_cfg2(self.signal_processes, self.signal_prior_cfg, options)}
        for p in self.producers:
            result['p%s' % p.name] = p.get_cfg(options)
            main['producers'].append('@p%s' % p.name)
        result['parameters'] = sorted(list(self.model.get_parameters(self.signal_processes, True)))
        rvobservables = sorted(list(self.model.rvobs_distribution.get_parameters()))
        if len(rvobservables) > 0:
            result['rvobservables'] = rvobservables
        result['observables'] = {}
        for obs in sorted(list(self.model.get_observables())):
            xmin, xmax, nbins = self.model.observables[obs]
            result['observables'][obs] = {'range': [xmin, xmax], 'nbins': nbins}
        return result
    
    # check in the cache directory whether result is already available there
    # returns True if it is, false otherwise.
    # Note that get_result_available will return True if the result is in the cache
    # only after calling check_cache
    def check_cache(self, options):
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
    
    # Returns the name of the config file created, including the full path.
    # creates the theta configuration file, if not done already.
    def get_configfile(self, options):
        if self.theta_cfg_fname is not None: return self.theta_cfg_fname
        cfg_dict = self._get_cfg(options)
        cfg_string = _settingvalue_to_cfg(cfg_dict, value_is_toplevel_dict = True)
        output_name = ''.join([p.name for p in self.producers]) + '-' + self.input_id + '-' + ''.join(self.signal_processes) + '-' + hashlib.md5(cfg_string).hexdigest()[:10]
        cfg_string += "output_name = \"%s.db\";\n" % output_name
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
            
    # runs theta
    def run_theta(self, options, in_background_thread = False):
        check_cache = options.getboolean('global', 'check_cache')
        if check_cache: self.check_cache(options)
        cfgfile = self.get_configfile(options)
        self.debug = options.getboolean('global', 'debug')
        if self.result_available: return
        workdir = options.get_workdir()
        cache_dir = os.path.join(workdir, 'cache')
        theta = os.path.realpath(os.path.join(config.theta_dir, 'bin', 'theta'))
        info("Running 'theta %s'" % cfgfile)
        to_execute = lambda : self._exec(theta + " " + cfgfile)
        if in_background_thread:
            assert self.thread is None
            self.thread = threading.Thread(target = to_execute)
            self.thread.start()
        else:
            to_execute()
            self._cleanup_exec()

    
    # check if a result is availbale, either after check_cache, or
    # after run_theta in a background thread.
    #
    # can throw an exception in case theta in the background thread exited with an exception
    def get_result_available(self):
        if self.thread is not None:
            if not self.thread.is_alive():
                # do the cleanup
                self.wait_for_result_available()
        return self.result_available
    
    # only useful after calling run_theta(in_background_thread = True)
    #
    # can throw an exception in case the thread exited with an exception
    def wait_for_result_available(self):
        if self.result_available: return
        if self.thread is None: return
        self.thread.join()
        self.thread = None
        self._cleanup_exec()
    
    # once the result is available, this is the path to the database filename
    def get_db_fname(self):
        return self.theta_db_fname
    
    # get the result (in the products table) of this invocation of theta as dictionary
    #  column_name --> list of values
    # which columns are returned can be controlled by the columns argument wich is either a list of column names
    # or as special case the string '*'. The defaul is '*' in which case all columns are returned.
    #
    #
    def get_products(self, columns = '*', fail_if_empty = True):
        if not self.result_available: return None
        assert self.theta_db_fname is not None
        if type(columns)!=str:
            columns = '"' + '", "'.join(columns) + '"'
        else: assert columns == '*'
        data, columns = sql_singlefile(self.theta_db_fname, 'select %s from products' % columns, True)
        if len(data)==0 and fail_empty:
            self.print_logtable_errors()
            raise RuntimeError, "No result is available, see errors above"
        result = {}
        for i in range(len(columns)):
            result[columns[i][0]] = [row[i] for row in data]
        return result
        
    def print_logtable_errors(self, limit = 10):
        if not self.result_available: return
        assert self.theta_db_fname is not None
        data = sql_singlefile(self.theta_db_fname, 'select "eventid", "message" from "log" where "severity"=0 limit %d' % int(limit))
        if len(data)==0: return
        print "There have been errors for some toys: eventid, message"
        for i in range(len(data)):
            print data[i][0], data[i][1]
        
    def get_products_nevents(self):
        data = sql_singlefile(self.theta_db_fname, 'select count(*) from products')
        return data[0][0]
        

# create a Histogram instance from the blob data in a sqlite db
def histogram_from_dbblob(blob_data):
    a = array.array('d')
    a.fromstring(blob_data)
    return Histogram(xmin = a[0], xmax = a[1], values = a[3:-1])

# the nuisance prior used for the 
"""
def toys(model, input, signal_prior, signal_processes, producers, n, nuisance_constraint = None, nuisance_prior_toys = None, options):
    data_source = DataSource(input, model, signal_processes, override_distribution = nuisance_prior_toys)
    r = Run(model, signal_processes, signal_prior, data_source, n, producers)
"""

def mle(model, input, signal_prior, signal_processes, n, nuisance_constraint = None, nuisance_prior_toys = None):
    r = Run(model, signal_processes, signal_prior, input, n, [MleProducer(model, signal_processes, nuisance_constraint)], nuisance_prior_toys)
    #if options is None: options = Options()
    return r
    #r.run_theta(options)
    #return r.get_products()




