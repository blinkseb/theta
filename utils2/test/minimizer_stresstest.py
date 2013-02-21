# -*- coding: utf-8 -*-
from theta_auto import test_model
import time, pickle

def shifted_dist(dist, shift):
  dist = copy.deepcopy(dist)
  for p in dist.get_parameters():
    dist.set_distribution_parameters(p, mean = shift)
  return dist

  
#print "root version: ", ROOT.gROOT.GetVersionInt()

N = 2000

options = Options()
#options.set('minimizer', 'strategy', 'minuit_vanilla')
options.set('minimizer', 'strategy', 'lbfgs_vanilla')
#options.set('minimizer', 'minuit_tolerance_factor', '0.1')

print "using tolerance %s" % options.getfloat('minimizer', 'minuit_tolerance_factor')

data = {}
for shift in (-1.0,):#(0.0, -1.0):
    for n in range(3,4):#range(2, 4):
        model = test_model.bernstein_model(n, n_events_total = 10000)
        sp = 'proc%d' % n
        model.set_signal_processes(sp)
        #t0 = time.time()
        zvalues = zvalue_approx(model, 'toys:0.0', n = N, nuisance_constraint = shifted_dist(model.distribution, shift), options = options)[sp]['Z']
        n_1sigma = len([z for z in zvalues if z >= 1.0])
        n_2sigma = len([z for z in zvalues if z >= 2.0])
        p1 = get_p(n_1sigma, len(zvalues))
        p2 = get_p(n_2sigma, len(zvalues))
        p1_expected = (1 - cl_1sigma) / 2
        p2_expected = (1 - cl_2sigma) / 2
        print p1, p1_expected, (p1[0] - p1_expected) / p1[1]
        print p2, p2_expected, (p2[0] - p2_expected) / p2[1]
        print "failures: %.3f %%" % (100 - (len(zvalues) * 100. / N))
        
        #d = plotutil.plotdata()
        #d.histogram(zvalues, 0, 5, 50)
        #plot([d], '$Z$', '$N$', 'zvalues.pdf')
        #print res
        #sys.exit(0)
        #print "shift=%f, n=%d" % (shift, n)
        ##print "time: ", time.time() - t0
        #success_rate = len(res['']['proc0_unc']) * 1.0 / N
        #print "success rate: %f" % success_rate
  
# tolerance=1e-3 in 5.32 corresponds internally to tolerance=1e-6 in 5.28

