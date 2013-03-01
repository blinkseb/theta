.. _examples:

********
Examples
********

This section contains many common examples. Most examples are shown using a simple counting experiment model constructed with the method
``test_model.simple_counting`` which has one signal process called ``s``, and an optional log-normal uncertainty on the background, modeled
via a nuisance parameter called ``"bunc"``. The model prediction for the Poisson mean is

.. math::

\mu = \beta_{\rm signal} \cdot s + e^{{\rm bunc}\cdot{\rm b_uncertainty}} \cdot b

where ``s``, `b`` and ``b_uncertainty`` are constants, passed as parameters to ``test_model.simple_counting``. The parameters of this model
are ``beta_signal`` and ``bunc``.


.. _examples_mle:

=======================
Maximum Likelihood Fits
=======================


.. _examples_mle_fit:

Perform a maximum likelihood fit
--------------------------------

See ``utils2/examples/mle/fit.py`` ::

 model = test_model.simple_counting(s = 10.0, n_obs = 12, b = 2.0)
 result = mle(model, input = 'data', n = 1)
 print result
 
This will print::

 {'s': {'__nll': [-17.818879797456002], 'beta_signal': [(0.9999999999999996, 0.34463876178582953)]}}
 
So the result is a python dictionary. As explained in :ref:`return_values`, the first-level key specifies
the signal group id, here "s". The dictionary `result['s']` has to keys: "__nll" contains the
negative log-likelihood values at the minimum. These values are not usually not useful on their own, but can be used to construct
likelihood ratios by making several maximum likelihood fit. The "beta_signal" contains a list of ``n = 1`` results, each
result is a two-tuple of the fitted parameter values and its uncertainty, so the result of the maximum likelihood fit
in this example can be summarized as

.. math::

 \hat{\beta}_{\rm signal} = 1.00 \pm 0.34


Making the fit more robust
--------------------------

In some cases, the fit might not converge. In this case, you can use the ``Options class`` (see :ref:`options_class`) for more control
about the minimizer settings. Also, in case you do not need the error, you can disable the error calculation.
So a robust maximum likelihood fit would be (see ``utils2/examples/mle/robust-fit.py``) ::

 model = test_model.bernstein_model(10) # a model of 10 bernstein polynomials, difficult to minimize
 mle(model, 'toys:0.0', 100) # will have many failures
 mle(model, 'toys:0.0', 100, with_error = False)
 options = Options()
 options.set('minimizer', 'strategy', 'robust')
 mle(model, 'toys:0.0', 100, with_error = False, options = options)

If you execute theta-auto on this python file, you'll note that the first theta execution will have many errors.
Specifying ``with_error = False`` for the ``mle`` method removes these errors already: this will not run MINUIT's migrad,
but you will not have an error estimate (if you need one, you can try using profile likelihood intervals via :meth:`theta_auto.pl_intervals` ).

The third call of ``mle`` in this example then uses different options for the minimizer.

.. warning:: If you find yourself in the situation where you the minimization fails, even the "robust" minimizer probably does not return
 the correct minimum. In this case, it might be a good idea to understand why the minimizer fails and
 to check that the result of the minimizer makes sense (e.g., by checking the pull distribution for toys).
  



