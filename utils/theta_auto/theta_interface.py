import numpy, os.path, hashlib, shutil, copy

import config as global_config
import utils
from plotutil import *

def settingvalue_to_cfg(value, indent=0, errors = []):
    if type(value) == numpy.float64: value = float(value)
    if type(value) == str: return '"%s"' % value
    if type(value) == bool: return 'true' if value else 'false'
    if type(value) == int: return '%d' % value
    if type(value) == float:
        if value == float("inf"): return '"inf"'
        elif value == -float("inf"): return '"-inf"'
        return '%.5e' % value
    if type(value) == list or type(value) == tuple or type(value) == numpy.ndarray:
        return "(" + ",".join(map(lambda v: settingvalue_to_cfg(v, indent + 4), value)) + ')'
    if type(value) == dict:
        result = "{\n"
        # sort keys to make config files reproducible:
        for s in sorted(value.keys()):
            result += ' ' * (indent + 4) + s + " = " + settingvalue_to_cfg(value[s], indent + 4) + ";\n"
        return result + ' ' * indent + "}"
    errors.append("cannot convert type %s to cfg; writing '<type error>' to allow debugging" % type(value))
    return '<type error for type %s>' % type(value)

#histogram is a tuple (xmin, xmax, data). Returns a dictionary for cfg building
def get_histo_cfg(histogram):
    xmin, xmax, data = histogram
    return {'type': 'direct_data_histo', 'range': [xmin, xmax], 'nbins': len(data), 'data': data[:]}

# signal_processes is a list of process names (strings) to consider *simultaneously* as signal
# It can also be an empty list.
def write_cfg(model, signal_processes, method, input, id = None, additional_settings = {}):
    for sp in signal_processes:
        assert sp in model.signal_processes, "invalid signal process '%s'" % sp
    model_parameters = sorted(list(model.get_parameters(signal_processes)))
    additional_settings['main']['output_database']['filename'] = '@output_name'
    additional_settings = copy.deepcopy(additional_settings)
    config = ''
    config += "parameters = " + settingvalue_to_cfg(model_parameters) + ";\n"
    obs = {}
    for o in model.observables:
        xmin, xmax, nbins = model.observables[o]
        obs[o] = {'range': [xmin, xmax], 'nbins': nbins}
    config += "observables = " + settingvalue_to_cfg(obs) + ";\n"
    cfg_model = model.get_cfg(signal_processes)
    if 'beta_signal' in model_parameters:
        cfg_model['parameter-distribution'] = product_distribution("@model-distribution-signal", model.distribution.get_cfg(model_parameters))
    else:
        if 'model-distribution-signal' in additional_settings: del additional_settings['model-distribution-signal']
        cfg_model['parameter-distribution'] = model.distribution.get_cfg(model_parameters)
    config += "model = " + settingvalue_to_cfg(cfg_model) + ";\n"
    for s in additional_settings:
        config += s + " = " + settingvalue_to_cfg(additional_settings[s]) + ";\n"
    m = hashlib.md5()
    m.update(config)
    hash = m.hexdigest()[:10]
    if id is None:
        name = '%s-%s-%s-%s' % (method, ''.join(signal_processes), input, hash)
    else:
        name = '%s-%s-%s-%s-%s' % (method, ''.join(signal_processes), input, id, hash)
    f = open(os.path.join(global_config.workdir, name + '.cfg'), 'w')
    print >>f, config
    print >>f, 'output_name = "%s.db";\n' % name
    f.close()
    return name


def delta_distribution(**kwargs):
    kwargs['type'] = 'delta_distribution'
    return kwargs

def product_distribution(*args):
    result = {'type': 'product_distribution'}
    result['distributions'] = args[:]
    return result
    
def equidistant_deltas(parameter, range, n):
    return {'type': 'equidistant_deltas', 'parameter': parameter, 'range': range, 'n': n}
    
def sqlite_database(fname = ''):
    return {'type': 'sqlite_database', 'filename': fname}

# cfg_names is a list of filenames (without ".cfg") which are expected to be located in the working directory
# valid options:
#  * nodel: if True, do not delete the db file in case of failure
def run_theta(cfg_names, **options):
    cache_dir = os.path.join(global_config.workdir, 'cache')
    if not os.path.exists(cache_dir): os.mkdir(cache_dir)
    theta = os.path.realpath(os.path.join(global_config.theta_dir, 'bin', 'theta'))
    for name in cfg_names:
        cfgfile = name + '.cfg'
        dbfile = name + '.db'
        cfgfile_cache = os.path.join(cache_dir, cfgfile)
        dbfile_cache = os.path.join(cache_dir, dbfile)
        cfgfile_full = os.path.join(global_config.workdir, cfgfile)
        already_done = False
        if os.path.exists(cfgfile_cache) and os.path.exists(os.path.join(cache_dir, dbfile)):
            # compare the config files:
            already_done = open(cfgfile_cache, 'r').read() == open(cfgfile_full, 'r').read()
        if already_done:
            #info("not running \"theta %s\": found up-to-date output in cache" % cfgfile)
            continue
        utils.info("running 'theta %s'" % cfgfile)
        retval = os.system(theta + " " + cfgfile_full)
        if retval != 0:
            if os.path.exists(dbfile) and not options.get('nodel', False): os.unlink(dbfile)
            raise RuntimeError, "executing theta for cfg file '%s' failed with exit code %d" % (cfgfile, retval)
        # move to cache, also the config file ...
        shutil.move(dbfile, dbfile_cache)
        shutil.copy(cfgfile_full, cfgfile_cache)

