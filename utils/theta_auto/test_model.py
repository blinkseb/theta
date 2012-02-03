# -*- coding: utf-8 -*-
# this file defines some simple models useful for testing theta-auto

from Model import *


# returns a model with a counting experiment in one bin with the given backgruond with a log-normal uncertainty.
# b_uncertainty is the absolute uncertainty on b.
#
# If s2 is not None, will return a model with two signal processes, "s" and "s2"
def simple_counting(s, n_obs, b=0.0, b_uncertainty=0.0, s2 = None):
    model = Model()
    model.set_data_histogram('obs', (0.0, 1.0, [float(n_obs)]))
    hf_s = HistogramFunction()
    hf_s.set_nominal_histo((0.0, 1.0, [float(s)]))
    model.set_histogram_function('obs', 's', hf_s)
    if s2 is not None:
        hf_s2 = HistogramFunction()
        hf_s2.set_nominal_histo((0.0, 1.0, [float(s2)]))
        model.set_histogram_function('obs', 's2', hf_s2)
    model.set_signal_processes('s*')
    if b > 0:
        hf_b = HistogramFunction()
        hf_b.set_nominal_histo((0.0, 1.0, [float(b)]))
        model.set_histogram_function('obs', 'b', hf_b)
        if b_uncertainty > 0: model.add_lognormal_uncertainty('bunc', 1.0 * b_uncertainty / b, 'b')
    return model


# returns a model in one bin with the given background. b_plus and b_minus are the background yields at +1sigma and -1sigma
def simple_counting_shape(s, n_obs, b=0.0, b_plus=0.0, b_minus=0.0):
    model = Model()
    model.set_data_histogram('obs', (0.0, 1.0, [float(n_obs)]))
    hf_s = HistogramFunction()
    hf_s.set_nominal_histo((0.0, 1.0, [float(s)]))
    model.set_histogram_function('obs', 's', hf_s)
    model.set_signal_processes('s*')
    if b > 0:
        hf_b = HistogramFunction()
        hf_b.set_nominal_histo((0.0, 1.0, [float(b)]))
        plus_histo = (0.0, 1.0, [float(b_plus)])
        minus_histo = (0.0, 1.0, [float(b_minus)])
        hf_b.set_syst_histos('bunc', plus_histo, minus_histo)
        model.distribution.set_distribution('bunc', 'gauss', mean = 0.0, width = 1.0, range = [-float("inf"), float("inf")])
        model.set_histogram_function('obs', 'b', hf_b)
    return model


