# -*- coding: utf-8 -*-
import scipy.special as sp
import scipy.stats
import math
import subprocess

import os.path

import config as global_config

# lmbda is the non-centrality parameter defined as
# lmbda = sum_i=1^k  mu_i^2 / sigma_i^2
#
# x can be a list in which case the list is returned
def noncentral_chi2_pdf(x, k, lmbda):
    if type(x) == float: xlist = [x]
    else: xlist = x
    if lmbda==0:
        result = map(lambda x: scipy.stats.chi2.pdf(x, k), xlist)
    else:
        bessel = os.path.realpath(os.path.join(global_config.theta_dir, 'bin', 'bessel'))
        proc = subprocess.Popen([bessel], stdin=subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        head = '%g %d %d ' % (lmbda, k, len(x))
        inp = head + ' '.join([('%g' % x0) for x0 in xlist])
        out, err = proc.communicate(inp)
        assert err == ''
        result = map(lambda s: float(s), out.split())
        proc.wait()
    if type(x) == float: return result[0]
    return result

phi = lambda x: scipy.stats.norm.cdf(x)
phi_inverse = scipy.stats.norm.ppf

class half_chi2_dist: # delta(0) + 1/2 chi2
    def cdf(self, tsval):
        return 0.5 + 0.5 * scipy.stats.chi2.cdf(tsval, 1)
    
    def icdf(self, q):
        if q <= 0.5: return 0
        return scipy.stats.chi2.ppf(2 * q - 1, 1)


# represents an asymptotic/approximate test statistic distribution of the 'lhclike' test statistic.
# objects of this class are immutable; all parameters are fixed at the time of construction.
#
# IMPORTANT: tsval has to be -2*ln LR
#
# this implements formulae 63 to 66 of http://arxiv.org/abs/1007.1727 v2
class asymptotic_ts_dist:
    # beta_signal_toy is the value which one would use in toy data generation, if there were toys; it's called mu' in the Ref.
    def __init__(self, beta_signal_toy, sigma):
        self.truth = float(beta_signal_toy)
        self.sigma = float(sigma)

    # beta_signal here is the signal scale factor used in the TS definition (the "mu" in "t_mu" in the Ref.)
    def _pdf(self, beta_signal, tsval):
        if tsval<=0.0: return 0.0
        if tsval <= beta_signal**2 / self.sigma**2:
            return 0.5  / math.sqrt(2 * math.pi) / math.sqrt(tsval) * math.exp(-0.5 * (math.sqrt(tsval) - (beta_signal - self.truth) / self.sigma)**2)
        else:
            return 0.5  / math.sqrt(2 * math.pi) * self.sigma / beta_signal * \
               math.exp(-0.5 * (tsval - (beta_signal**2 - 2 * self.truth * beta_signal) / self.sigma**2)**2  /  (2 * beta_signal / self.sigma)**2)
        
    # returns a vector of y values appropriate for displaying a histogram with nbins bins. x_binborders should contains
    # all bin borders, i.e., it is of length nbins+1.
    # Note that the first bin should start exactly at 0.0.
    def histo_yvalues(self, beta_signal, x_binborders):
        assert x_binborders[0] == 0.0
        binranges = zip(x_binborders[:-1], x_binborders[1:])
        result = map(lambda r: self.cdf(beta_signal, r[1]) - self.cdf(beta_signal, r[0]), binranges)
        result[0] += phi((self.truth - beta_signal) / self.sigma)
        return result

    def cdf(self, beta_signal, tsval):
        if tsval <= beta_signal**2 / self.sigma**2:
            return phi(math.sqrt(tsval) - (beta_signal - self.truth) / self.sigma)
        else:
            return phi((tsval - (beta_signal**2 - 2*beta_signal * self.truth) / self.sigma**2) / (2 * beta_signal / self.sigma))

    def icdf(self, beta_signal, F):
        assert F > 0 and F < 1
        F_flip = phi(self.truth / self.sigma)
        if F < F_flip:
            return max([phi_inverse(F) + (beta_signal - self.truth) / self.sigma, 0.0])**2
        else: return 2 * phi_inverse(F) * beta_signal / self.sigma + (beta_signal**2 - 2 * beta_signal * self.truth) / self.sigma**2
        

