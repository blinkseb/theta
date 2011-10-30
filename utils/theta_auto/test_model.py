# this file defines some simple models useful for testing theta-auto

from Model import *

# returns a model with a counting experiment in one bin with the given backgruond with a log-normal uncertainty.
def simple_counting(s, n_obs, b=0.0, b_uncertainty=0.0):
    model = Model()
    model.set_data_histogram('obs', (0.0, 1.0, [float(n_obs)]))
    hf_s = HistogramFunction()
    hf_s.set_nominal_histo((0.0, 1.0, [float(s)]))
    model.set_histogram_function('obs', 's', hf_s)
    model.set_signal_processes('s')
    if b > 0:
        hf_b = HistogramFunction()
        hf_b.set_nominal_histo((0.0, 1.0, [float(b)]))
        model.set_histogram_function('obs', 'b', hf_b)
        if b_uncertainty > 0: model.add_lognormal_uncertainty('bunc', 1.0 * b_uncertainty / b, 'b')
    return model


