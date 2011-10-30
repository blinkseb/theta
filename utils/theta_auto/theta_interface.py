import numpy, os.path, hashlib, shutil, copy, threading

import config as global_config
import utils
from plotutil import *

def settingvalue_to_cfg(value, indent=0, current_path = []):
    if type(value) == numpy.float64: value = float(value)
    if type(value) == str: return '"%s"' % value
    if type(value) == bool: return 'true' if value else 'false'
    if type(value) == int: return '%d' % value
    if type(value) == float:
        if value == float("inf"): return '"inf"'
        elif value == -float("inf"): return '"-inf"'
        return '%.5e' % value
    if type(value) == list or type(value) == tuple or type(value) == numpy.ndarray:
        return "(" + ",".join([settingvalue_to_cfg(value[i], indent + 4, current_path + [str(i)]) for i in range(len(value))]) + ')'
    if type(value) == dict:
        result = "{\n"
        # sort keys to make config files reproducible:
        for s in sorted(value.keys()):
            result += ' ' * (indent + 4) + s + " = " + settingvalue_to_cfg(value[s], indent + 4, current_path + [s]) + ";\n"
        return result + ' ' * indent + "}"
    raise RuntimeError, "Cannot convert type %s to theta cfg in path '%s'" % (type(value), '.'.join(current_path))

#histogram is a tuple (xmin, xmax, data). Returns a dictionary for cfg building
def get_histo_cfg(histogram):
    xmin, xmax, data = histogram
    return {'type': 'direct_data_histo', 'range': [xmin, xmax], 'nbins': len(data), 'data': data[:]}

# signal_processes is a list of process names (strings) to consider *simultaneously* as signal
# It can also be an empty list.
#
# options are passed to model.get_cfg
def write_cfg(model, signal_processes, method, input, id = None, additional_settings = {}, **options):
    for sp in signal_processes:
        if sp == '__fakesp__': continue
        assert sp in model.signal_processes, "invalid signal process '%s'" % sp
    model_parameters = sorted(list(set(model.get_parameters(signal_processes))))
    all_parameters = set(model.get_parameters(signal_processes))
    if model.additional_nll_term is not None: all_parameters.update(model.additional_nll_term.get_parameters())
    all_parameters = sorted(list(all_parameters))
    additional_settings['main']['output_database']['filename'] = '@output_name'
    additional_settings = copy.deepcopy(additional_settings)
    config = ''
    config += "parameters = " + settingvalue_to_cfg(all_parameters, 0, ['parameters']) + ";\n"
    obs = {}
    for o in model.observables:
        xmin, xmax, nbins = model.observables[o]
        obs[o] = {'range': [xmin, xmax], 'nbins': nbins}
    config += "observables = " + settingvalue_to_cfg(obs, 0, ['observables']) + ";\n"
    cfg_model = model.get_cfg(signal_processes, **options)
    if 'beta_signal' in model_parameters:
        cfg_model['parameter-distribution'] = product_distribution("@model-distribution-signal", model.distribution.get_cfg(model_parameters))
    else:
        if 'model-distribution-signal' in additional_settings: del additional_settings['model-distribution-signal']
        cfg_model['parameter-distribution'] = model.distribution.get_cfg(model_parameters)
    config += "model = " + settingvalue_to_cfg(cfg_model, 0, ['model']) + ";\n"
    for s in additional_settings:
        config += s + " = " + settingvalue_to_cfg(additional_settings[s], 0, [s]) + ";\n"
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

def _run_theta_single(name, debug):
    cache_dir = os.path.join(global_config.workdir, 'cache')
    cfgfile = name + '.cfg'
    dbfile = name + '.db'
    cfgfile_cache = os.path.join(cache_dir, cfgfile)
    dbfile_cache = os.path.join(cache_dir, dbfile)
    cfgfile_full = os.path.join(global_config.workdir, cfgfile)
    already_done = False
    theta = os.path.realpath(os.path.join(global_config.theta_dir, 'bin', 'theta'))
    if os.path.exists(cfgfile_cache) and os.path.exists(os.path.join(cache_dir, dbfile)):
        # compare the config files:
        already_done = open(cfgfile_cache, 'r').read() == open(cfgfile_full, 'r').read()
    if already_done:
        utils.info("Skipping 'theta %s': found corresponding output file in cachedir" % cfgfile)
        return
    utils.info("Running 'theta %s'" % cfgfile)
    params = ""
    if debug: params += " --redirect-io=False"
    retval = os.system(theta + params + " " + cfgfile_full)
    if retval != 0:
        if os.path.exists(dbfile) and not debug: os.unlink(dbfile)
        raise RuntimeError, "executing theta for cfg file '%s' failed with exit code %d" % (cfgfile, retval)
    # move to cache, also the config file ...
    shutil.move(dbfile, dbfile_cache)
    shutil.copy(cfgfile_full, cfgfile_cache)

"""
def run_function_then_notify(f, condvar):
    ex, result = None, None
    try:
        result = f()
    except Exception, e:
        ex = e
    condvar.acquire()
    condvar.notify()
    condvar.release()
    if ex is not None: raise ex
    return result


# like map(func, l), but run in threads.
def parallel_map(func, l, n_threads):
    done = [False] * len(l)
    threads_current = []
    condvar = threading.Condition()
    while not all(done):
        while len(threads_current) < n_threads and not all(done):
            i=0
            while done[i]: i+=1
            done[i] = True
            t = threading.Thread(target = lambda: run_function_then_notify(lambda: func(l[i]), condvar))
            threads_current.append(t)
            t.start()
            print 'current threads: ', len(threads_current)
        condvar.acquire()
        condvar.wait(5.0)
        # remove non-alive threads:
        for i in range(len(threads_current)):
            if not threads_current[i].isAlive():
                threads_current[i].join()
                del threads_current[i]
"""

# cfg_names is a list of filenames (without ".cfg") which are expected to be located in the working directory
# valid options:
#  * debug: if True, do not delete the db file in case of failure and do not redirect theta stdio
#  * run_theta_parallel: number of concurrent theta processes to start (which run different config files). Default is 1
def run_theta(cfg_names, **options):
    cache_dir = os.path.join(global_config.workdir, 'cache')
    if not os.path.exists(cache_dir): os.mkdir(cache_dir)
    debug = options.get('debug', False)
    #run_theta_parallel = options.get('run_theta_parallel', 1)
    #parallel_map(lambda name: _run_theta_single(name, debug), cfg_names, run_theta_parallel)
    for name in cfg_names:
        _run_theta_single(name, debug)

