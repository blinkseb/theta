# -*- coding: utf-8 -*-
import scipy.special as sp
from theta_auto.test_model import *
from theta_auto.asymptotic_cls import *

N = 100000
xmin = 0
xmax = 5
nbins = 60
#betas = [0.01, 0.05, 0.5, 1.0]
betas = [0.05, 0.5]


spid = 's'

def colors():
    yield '#000000'
    yield '#888888'
    yield '#880000'
    yield '#008800'
    yield '#000088'
    yield '#888800'
    yield '#880088'
    yield '#008888'
    yield '#ff0000'
    yield '#00ff00'
    yield '#0000ff'
    yield '#ffff00'
    yield '#ff00ff'
    yield '#00ffff'

def count(pred, l):
    c = 0
    for el in l:
        if pred(el): c += 1
    return c

# Reference: http://arxiv.org/abs/1007.1727


# method 1: get sigma from the likelihood curvature of asimov data at that beta value:
def sigma1(model, beta_signal, spid):
    res = mle(model, input = 'toys-asimov:%g' % beta_signal, n=1, signal_prior = 'flat:[-inf, inf]')
    return res[spid]['beta_signal'][0][1]

# method 2 (usually preferred): get sigma from the value of the likelihood ratio in an asimov toy.
# This is called "sigma_A" in the reference, see eq. 30, 31 there:
def sigma2(model, beta_signal):
    res = deltanll(model, input = 'toys-asimov:%g' % beta_signal, n=1, signal_prior = 'flat:[-inf, inf]')
    spid = res.keys()[0]
    lr = res[spid][0]
    return math.sqrt(beta_signal**2 / (2*lr))

    
# returns a two-tuple (par_values, nll), where par_values are the paameter values, as dictionary, and nll is the negative log-likelihood.
def get_data_fitvalues(model, spid, signal_prior = 'flat'):
    signal_process_groups = {spid: model.signal_process_groups[spid]}
    pars = sorted(list(model.get_parameters(signal_process_groups[spid])))
    res = mle(model, 'data', 1, signal_process_groups = signal_process_groups, signal_prior = signal_prior, with_covariance = True)
    m = matrix_from_dbblob(res[spid]['__cov'][0])
    ibeta = pars.index('beta_signal')
    n = len(pars)
    nuisance_matrix = numpy.zeros((n-1, n-1))
    for i in range(n-1):
        io = 1 if i >= ibeta else 0
        for j in range(n-1):
            jo = 1 if j >= ibeta else 0
            nuisance_matrix[i][j] = m[i+io][j+jo]
    nuisance_vector = [res[spid][par][0][0] for par in pars if par!='beta_signal']
    result = {}
    for par in pars:
        result[par] = res[spid][par][0][0]
    nll = res[spid]['__nll'][0]
    return result, GaussDistribution(pars[0:ibeta] + pars[ibeta+1:], nuisance_vector, nuisance_matrix), nll
    
# return a tuple (nuisance parameter Distribution, beta_signal) for nuisance parameters from the fit to data
# with the given signal_prior value.
#def get_nuisance_dist_from_fit(model, spid, signal_prior):
#    vals, nll = get_data_fitvalues(model, spid, signal_prior)
#    beta_signal = vals['beta_signal']
#    return get_fixed_dist_at_values(vals), beta_signal
    
# Returns a new model, leaves the "model" parameter unchanged.
#
# makes a fit to data to get the parameter values to sample from. 
#
# * for each nuisance parameter:
#   - add a real-valued observable "rvobs_<parameter name>"
#   - add a Gaussian distribution for the rvobs with mean = nuisance parameter, parameter = rvobs_...
# * fix all nuisance parameters to the value fitted in data with a delta function.
# * add the default value to the real-valued observables in data (according to model.distribution)
#
# make all nuisance parameter priors flat, s.t. sampling nuisance parameters will be
#   around default_nuisance_values (use the ML estimate on data here ...). This ensures that
#   sampling always uses the values fitted in data, and all constraints come from the real-valued observables,
#   not from the prior.
"""
def frequentize_model(model, spid, signal_prior):
    signal_process_groups = {spid: model.signal_process_groups[spid]}
    parameters = model.distribution.get_parameters()
    par_values, dist, nll = get_data_fitvalues(model, spid, signal_prior)
    print "fitted nuisances for signal_prior = '%s': %s" % (signal_prior, str(par_values))
    print "  dist: " + str(dist)
    result = model.copy()
    for p in parameters:
        prior_nuisance = model.distribution.get_distribution(p)
        # have to use the conjugate distribution here. gauss is self-conjugate, so no problem here in most cases:
        assert prior_nuisance['typ'] == 'gauss'
        rvobs = 'rvobs_%s' % p
        dist.rename_parameter(p, rvobs)
        #result.rvobs_distribution.set_distribution(rvobs, 'gauss', mean = p, width = prior_nuisance['width'], range = prior_nuisance['range'])
        result.distribution.set_distribution_parameters(p, mean = par_values[p], width = inf)
        #result.distribution.set_distribution_parameters(p, mean = par_values[p], width = 0.0, range = (par_values[p], par_values[p]))
        result.data_rvobsvalues[rvobs] = prior_nuisance['mean'] # usually 0.0    
    result.rvobs_distribution = dist
    return result
"""

# construct original model
#model, mname = simple_counting(s = 100.0, n_obs = 1050.0, b = 1000.0, b_uncertainty = 100.0), '1cc'

def poisson_lr_p(pd_data, pd_model):
    assert len(pd_data)==len(pd_model)
    nll = 0.0
    ndof = 0
    for n,mu in zip(pd_data, pd_model):
        if n > 0.0: nll += n * math.log(n / mu)
        nll += mu - n
        ndof += 1
    return scipy.stats.chi2.sf(2 * nll, ndof-1)

def poisson(mus): return scipy.random.poisson(mus)


def test_zval_dist(model, n = 10000, xmin = 0.0, xmax = 4.0, nbins=50, plot_fname = None, signal_prior = 'flat:[0,inf]'):
    spid = model.signal_process_groups.keys()[0]
    tsvals = zvalue_approx(model, 'toys:0.0', 10000, signal_prior_sb = signal_prior)
    h = plotutil.plotdata()
    zs = tsvals[spid]['Z']
    #print "above 0: ", len([z for z in zs if z >= 0.0])
    #print "below 0: ", len([z for z in zs if z < 0.0])
    #sys.exit(0)
    h.histogram(tsvals[spid]['Z'], xmin, xmax, nbins, errors = True, include_uoflow = True)
    h_expected = plotutil.plotdata(color = '#ff0000')
    h_expected.x = h.x
    norm = len(tsvals[spid]['Z'])
    h_expected.y = map(lambda r: norm * (scipy.stats.norm.cdf(r[1]) - scipy.stats.norm.cdf(r[0])), zip([-inf] + h.x[1:], h.x[1:] + [inf]))
    """
    h_expected2 = plotutil.plotdata(color = '#0000ff')
    h_expected2.x = h.x
    s = sigma2(model, 1.0)
    print "sigma: %.3f" % s
    print "sigmas: ", [sigma2(model, b) for b in (0.0, 0.01, 0.1, 0.5, 1.0, 2.0)]
    dist = asymptotic_ts_dist(0.0, s)
    h_expected2.y = dist.histo_yvalues(0.0, h_expected2.x + [inf])
    """
    if plot_fname is not None:
        h.legend = 'Toys'
        pval = poisson_lr_p(h.y, h_expected.y)
        h_expected.legend = 'Gauss(0, 1): p=%.3g' % pval
        #h_expected2.legend = 'Asymptotic: p=%.3g' % poisson_lr_p(h.y, h_expected2.y)
        h.lw, h_expected.lw= 1, 1
        #plot([h_expected, h_expected2, h], 'Z', 'N', plot_fname, logy = True, ymin = 0.1)
        plot([h_expected, h], 'Z', 'N', plot_fname, logy = True, ymin = 0.1)
    return pval
    
"""
mus = [20., 10., 100.0]
ps = []
for i in range(1000):
    data = poisson(mus)
    p = poisson_lr_p(data, mus)
    ps.append(p)

pd_p = plotutil.plotdata()
pd_p.histogram(ps, 0.0, 1.0, 50)
plot([pd_p], 'p', 'N', 'ps.pdf', ymin = 0.0)

sys.exit(0)
""" 

sqrt2 = lambda x: math.sqrt(max([0.0, 2*x]))

# asymptotic (freq only) works more or less:
model1 = multichannel_counting(signals = [1000.0, 1000.0], n_obs = [11000.0, 11000.0], backgrounds = [10000.0, 10000.0], b_uncertainty1 = [0.1, 0.2], b_uncertainty2 = [0.3, 0.1])
model2 = multichannel_counting(signals = [1000.0, 1000.0], n_obs = [11000.0, 11000.0], backgrounds = [10000.0, 10000.0], b_uncertainty1 = [0.1, 0.2], b_uncertainty2 = [0.3, -0.25])
model3 = simple_counting_shape(100., 10100, 10000, 10300, 9900)


for model, name in zip((model1, model2, model3), ('model1', 'model2', 'model3')):
    #  1. model directly:
    pval = test_zval_dist(model, plot_fname = '%s.pdf' % name)
    print "%s: %.3g" % (name, pval)
    #assert pval < 1e-10


    # 2. frequentized version:
    model_freq = frequentize_model(model)
    pval = test_zval_dist(model_freq, plot_fname = '%s_freq.pdf' % name)
    print "%s (freq): %.3g" % (name, pval)
    #assert pval > 1e-3


sys.exit(0)


"""
sys.exit(0)


#model, mname = multichannel_counting(signals = [1000.0, 1000.0], n_obs = [11000.0, 11000.0], backgrounds = [10000.0, 10000.0], b_uncertainty1 = [0.05, 0.02], b_uncertainty2 = [0.02, 0.05]), '3'
#model, mname = multichannel_counting(signals = [1000.0, 1000.0], n_obs = [11000.0, 11000.0], backgrounds = [10000.0, 10000.0], b_uncertainty1 = [0.3, 0.3], b_uncertainty2 = [0.3, 0.3]), '30'


#model_freq = frequentize_model(model, 's', signal_prior = 'fix:1.0')

#dat = make_data(model_freq, input = 'toys-asimov:1.0', n=10)
#print dat
#sys.exit(0)

vals, dist, nll = get_data_fitvalues(model_freq, 's', 'flat')
print 's+b fit: ', vals, nll
vals, dist, nll = get_data_fitvalues(model_freq, 's', 'fix:0.0')
print 'b fit: ', vals, nll
sys.exit(0)
"""

# make fixed model around fitted nuisance parameter values:
    
#model_fixed = model.copy()
#dist, bs = get_nuisance_dist_from_fit(model, 's', signal_prior = 'fix:0.0')
#model_fixed.distribution = dist



pds = []
pds_cdf = []
binwidth = 1.0 * (xmax - xmin) / nbins

tsvals = {}
for beta in betas:
    model_freq = frequentize_model(model, 's', signal_prior = 'fix:%g' % beta)
    tsvals[('s', beta)] = deltanll(model_freq, input = 'toys:%g' % beta, n = N, lhclike_beta_signal = beta)['s']


# beta is used for both the "beta_signal_toy" value for toy data generation and the "beta_signal_ts" value
# which is used in the test statistic definition.
for beta in betas:
    # frequentize model around the found beta_signal value:
    model_freq = frequentize_model(model, 's', signal_prior = 'fix:%g' % beta)
    beta_signal_sigma = sigma1(model_freq, beta, 's')
    beta_signal_sigma2 = sigma2(model_freq, beta, 's')
    print "sigma_beta for beta=%g: %g (method2: %g)" % (beta, beta_signal_sigma, beta_signal_sigma2)
    pd = plotutil.plotdata()
    pd.histogram([ max(-x*2, 0.0) for x in tsvals[('s', beta)]], xmin, xmax, nbins)
    #pd.histogram([max(-x*2, 0.0) for x in tsvals[beta]['s']], xmin, xmax, nbins)
    for name, bs_sigma in zip(('method 1', 'method 2'), (beta_signal_sigma, beta_signal_sigma2)):
        print "starting ", name
        pd_as = plotutil.plotdata(as_function = True, legend = name)
        pd_as.x = numpy.linspace(0.0, xmax, nbins)
        dist = asymptotic_ts_dist(beta, bs_sigma)
        pd_as.y = numpy.ndarray(shape = len(pd_as.x))
        #pd_as.y = [dist.pdf(beta, x, binwidth, 20) * N for x in pd_as.x]   
        for i in range(len(pd_as.x)): pd_as.y[i] = dist.pdf(beta, pd_as.x[i], binwidth, 20) * N
        # shift the x values to coincide with the bin centers of the histogram:
        pd_as.x = [x + 0.5 * binwidth for x in pd_as.x]
        pds.append(pd_as)
    #ibin = int(-xmin / binwidth+1.5)
    #print beta, pd_as.y[0], sum(pd_as.y), sum(pd.y[:ibin]), sum(pd.y)
    pds.append(pd)
   
    """
    pd_cdf = plotutil.plotdata()
    pd_cdf.x = list(numpy.linspace(0.0, xmax, nbins))
    total = len(tsvals[('s', beta)])
    pd_cdf.y = [count(lambda x: -2*x > x0, tsvals[('s', beta)]) * 1.0 / total for x0 in pd_cdf.x]
    for name, bs_sigma in zip(('method 1', 'method 2'), (beta_signal_sigma, beta_signal_sigma2)):
        pd_as_cdf = plotutil.plotdata(as_function = True, legend = name)
        pd_as_cdf.x = pd_cdf.x
        dist = asymptotic_ts_dist(beta, bs_sigma)
        pd_as_cdf.y = [dist.cdf(beta, x) for x in pd_as_cdf.x]
        # also shift the x-values:
        pd_as_cdf.x = pd_as.x
        pds_cdf.append(pd_as_cdf)
    pds_cdf.append(pd_cdf)
    """

for pd, color in zip(pds, colors()): pd.color = color
#for pd, color in zip(pds_cdf, colors()): pd.color = color

plotutil.plot(pds, 'TS', 'N', 'ts-%s.pdf' % mname, logy = True, ymin = 1.0, xmax = xmax)
#plotutil.plot(pds_cdf, 'TS', 'N', 'ts_cdf.pdf', logy = True, ymin = 0.001, xmax = xmax)

print "exiting"

