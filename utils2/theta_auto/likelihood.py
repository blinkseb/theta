# -*- coding: utf-8 -*-

from utils import *
from theta_interface import *
import itertools


# The result is a dictionary (spid) --> (parameter name) --> (list of results of length n)
# The "list of results" contains pairs of (value, error). In case with_error is False, error is always None.
# Note that the list can have less entries than n, if the minimization failed for some toys. If it fails for all toys,
# an RuntimeError is thrown.
def mle(model, input, n, with_error = True, signal_process_groups = None, nuisance_constraint = None, nuisance_prior_toys = None, signal_prior = 'flat', options = None):
    if signal_process_groups is None: signal_process_groups = model.signal_process_groups
    if options is None: options = Options()
    result = {}
    for spid, signal_processes in signal_process_groups.iteritems():
        r = Run(model, signal_processes, signal_prior = signal_prior, input = input, n = n,
             producers = [MleProducer(model, signal_processes, nuisance_constraint, signal_prior, need_error = with_error)], nuisance_prior_toys = nuisance_prior_toys)
        r.run_theta(options)
        res = r.get_products()
        parameters = [m.group(1) for m in map(lambda colname: re.match('mle__(.+)_error', colname), res.keys()) if m is not None]
        result[spid] = {}
        for p in parameters:
            if not with_error: result[spid][p] = list(itertools.izip_longest(res['mle__%s' % p], []))
            else: result[spid][p] = zip(res['mle__%s' % p], res['mle__%s_error' % p])
    return result

# The result is a dictionary (spid) -> (cl) -> (list of results)
# where cl is one floating point value passed in the cls argument. Each item in the list of results is a two-tuple (lower, upper)
# which define the both ends of this profile likelihood interval.
# The special cl=float(0.0) is always filled. This corresponds to the maximum lieklihood value of the parameter. In this case,
# each element in the list of results is one floating point value instead of a pair.
def pl_interval(model, input, n, cls = [cl_1sigma, cl_2sigma], signal_process_groups = None, nuisance_constraint = None, nuisance_prior_toys = None, signal_prior = 'flat', options = None, parameter = 'beta_signal'):
    if signal_process_groups is None: signal_process_groups = model.signal_process_groups
    if options is None: options = Options()
    result = {}
    colnames = ['pli__lower%05d' % int(cl*10000 + 0.5) for cl in cls] + ['pli__upper%05d' % int(cl*10000 + 0.5) for cl in cls] + ['pli__maxl']
    for spid, signal_processes in signal_process_groups.iteritems():
        r = Run(model, signal_processes, signal_prior = signal_prior, input = input, n = n,
             producers = [PliProducer(model, signal_processes, nuisance_constraint, cls = cls, parameter = parameter)],
             nuisance_prior_toys = nuisance_prior_toys)
        r.run_theta(options)
        res = r.get_products(colnames)
        result[spid] = {0.0: res['pli__maxl']}
        for i, cl in enumerate(cls):
            result[spid][cl] = zip(res[colnames[i]], res[colnames[len(cls) + i]])
    return result


# Get the approximate "Z value" (significance in sigma), based on Wilks' Theorem. Note that this only works if
# the conditions for Wilks' Theorem apply which is assume here for one degree of freedom.
#
# The result is a dictionary (spid) --> (list of z values)
#
# Note that Z-values can be <0 in principle, but this method yields always values >=0.0 by construction.
def zvalue_approx(model, input, n, signal_process_groups = None, nuisance_constraint = None, nuisance_prior_toys = None, signal_prior = 'flat', options = None):
    if signal_process_groups is None: signal_process_groups = model.signal_process_groups
    if options is None: options = Options()
    result = {}
    for spid, signal_processes in signal_process_groups.iteritems():
        r = Run(model, signal_processes = signal_processes, signal_prior = signal_prior, input = input, n=n,
             producers = [DeltaNllHypotest(model, signal_processes, nuisance_constraint)], nuisance_prior_toys = nuisance_prior_toys)
        r.run_theta(options)
        data = r.get_products(['dnll__nll_diff'])
        result[spid] = [math.sqrt(max(x, 0.0)) for x in data['dnll__nll_diff']]
    return result


