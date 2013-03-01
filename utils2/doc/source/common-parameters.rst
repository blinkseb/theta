.. _common_parameters:

*****************
Common Parameters
*****************

Many python methods in ``theta-auto`` have similar arguments:
 * ``model`` is  the statistical model to use; it should be a Model instance.
 * ``input`` is a string specifying on which (toy-)data a method should be run on. This is a string which is either
   - "data" to run on real data
   - "toys:X" where X is a floating point value. This uses toy data as input where beta_signal is set to X. Nuisance parameters are sampled according to their priors, using ``model.distribution``, possible overridden by ``nuisance_prior_toys`` (see below)
   - "toys-asimov:X" where X is a floating point value. This uses asimov toy data (i.e., toy data without Poisson smearing) where beta_signal is set to X. Nuisance parameters are fixed to the mean value from model.distribution (usually 0, but can be modified); this can be overridden by ``nuisance_prior_toys`` (see below)
   - "replay-toys:P" where P is the path to a .db output file which contains pseudo data from an earlier theta run. This is usually only needed for special cases in which e.g., the model for toy data and for the actual evaluation are very different. The matching between the channels from these both models is string-based; the binning must be consistent in range and number of bins.
 * ``n`` is the number of 'toys' to perform:  :program:`theta` does ``n`` passes of (i) getting data from the data source according to ``input`` and (ii) running the statistical
   method on this toy dataset. For ``input="data"` or ``input="toys-asimov:X"`` and a deterministic statistical method, you usually want to set ``n = 1``.
   For input="toys:X" or in case the statistical method is not deterministic (e.g., Markov-Chain Monte-Carlo), it makes sense to set it
   to a higher number (say ``n=100`` as a baseline) to study also the spread of the result.
 * ``signal_process_groups`` is an alternative for ``model.signal_process_groups``, i.e., a dictionary using signal process id as key and a list of
   process names (the signal process group) as value. The default ``None`` will use
   ``model.signal_process_groups``. See :ref:`what_is_signal` for details.
 * ``signal_prior`` is a string specifying the prior for beta_signal to use in the likelihood function. Valid values are
   - "fix:X", where X is the floating point value to fix beta_signal to
   - "flat" for a flat prior on the interval [0, infinity]
   - "flat:[X,Y]" where X, Y are floating point values (including "-inf", "inf") for a flat prior on the interval [X, Y]
 * ``nuisance_constraint`` the :class:Distribution instance to use in the likelihood function for the nuisance parameters. The default value ``None`` will use ``model.distribution``,
   ``nuisance_constraint`` does not have to define a prior for all model parameters:
   if ``nuisance_constraint`` does not define a prior for a certain parameter, the one from ``model.distribution`` is used. This makes it easier to "partially override" ``model.distribution``
   with a different constraint.
 * ``nuisance_prior_toys``: the prior distribution for nuisance parameters to use for toy data generation. This is only applicable in case of ``input="toys..."`` and will be ignored for
   other cases e.g., ``input="data"``. The default value ``None`` is to use ``model.distribution`` in case of ``input="toys:X"`` and a ``Distribution`` in which all parameters
   are fixed to the mean given in ``model.distribution`` in case of input="toys-asmiov:XXX". As for ``nuisance_prior``, you need only to include those parameters in
   ``nuisance_prior_toys`` for which you want to override this default.
 * ``options`` is an instance of the ``Options`` class. This parameter controls some technical details of the actual theta run, such as multithreading and minimizer settings; it should usually not affect the result. The default ``None`` will use a default-constructed ``Options`` object. This is documented in the next section.

.. _options_class:

=================
The Options class
=================

The options class inherits from python's `SafeConfigParser <http://docs.python.org/2/library/configparser.html>`_, which defines a simple ini-like
configuration file syntax with key/value settings.

The usual use within theta-auto is to default-construct an ``Options`` instance, which initializes a number of default values as discussed below.
You can then use the ``set`` method to set new parameters. So an example usage to run many maximum likelihood fits on pseudo data in multiple threads of theta would be::

 options = Options()
 options.set('main', 'n_threads', '20')
 mle(model, input = 'toys:1.0', n = 20000)
 
Note that the ``set`` method always requires a string, even if the setting is of a different type.

The default options are the result of parsing this ini-file configuration:

.. code-block:: none

 [global]
 debug = False
 check_cache = True
 
 [main]
 n_threads = 1
 
 [minimizer]
 strategy = fast
 mcmc_iterations = 1000
 always_mcmc = False
 bootstrap_mcmcpars = 0
 minuit_tolerance_factor = 1
       
 [cls_limits]
 write_debuglog = True
       
 [model]
 use_llvm = False
 
 [mcmc]
 strategy = asimov_widths
 stepsize_factor = None

All of these settings are discussed now, starting at tyhe top.

The option ``global debug`` controls how the :program:`theta` program is executed, adding some more verbosity and timing information. Future versions
might use this flag at other places in theta-auto as well.

The option ``global check_cache`` enables / disables cache checking in the method `theta_auto.Run.run_theta`: with check caching disabled, ``run_theta`` will
always execute :program:`theta` locally, even if a matching result .db file is found in the "analysis/cache" directory.

``main n_threads`` configures the number of threads used for theta. Note that this option has no effect for CLs limits.

The ``minimizer`` section controls aspects of the minimizer(s) used in theta. ``strategy`` selects a minimizer algorithm. The default value "fast"
should be a reasonable choice in most circumstances. In this setting, a chain of 4 minimizers is used: first root's minuit is tried. If it fails, a Markov Chain
is run. If this fails, root's simplex is run. If this fails, another Markov Chain is run. As a final step -- if uncertainties are required --
minuit is run at the very end to calculate the uncertainties. The default number of MCMC iterations for the ``fast`` strategy is controlled with the option ``mcmc_iterations`.

Other available values for ``minimizer strategy`` are "robust" and "minuit_vanilla". "robust" works the same as "fast", but adds to the chain
of minimizers two additional mcmc minimizers with a much larger number of iterations.
This can have substantial impact on the runtime and is therefore not the default. the strategy "minuit_vanilla" is present mainly for test purposes: it runs only minuit.

If ``minimizer always_mcmc`` is set to ``True``, only the Markov Chains are used in the minimizer chain, minuit is not run.

``minimizer minuit_tolerance_factor`` can be used to scale up the tolerance of minuit by th specified factor.

``cls_limits write_debuglog`` controls whether or not the text file with detailed debug information is created if running cls limits.

``model use_llvm`` controls whether llvm-compilation of the model is enabled. For a small set of models, this decreases execution time substantially (especially
for models with many one-bin channels). This option requires to compile with llvm. Note that llvm is somewhat experimental; you are advised to compare the results
without llvm.

The section ``mcmc`` controls Markov-Chain Monte-Carlo settings used by Bayesian methods (note that these settings do not affect the mcmc run in the minimizer):
``mcmc strategy`` controls how proposal steps in the Markov Chain are constructed. The default of "asimov_widths" uses a Gaussian proposal function with a diagonal
covariance matrix. This covariance matrix is obtained by inspecting the likelihood function for the Asimov dataset along the axes. This default should be a reasonable
choice for most cases. Other available strategies are the plugins of base type MCMCStrategy (see the `doxygen documentation <http://www-ekp.physik.uni-karlsruhe.de/~ott/theta/testing/html/classtheta_1_1MCMCStrategy.html>`_): "asimov_der_cov" will use the full covariance matrix (i.e., including off-diagonals) from
the asimov likelihood. This usually has even better convergence properties than the default "asimov_widths", but it takes longer to initialize, especially for models
with a large number of dimensions (say more than 100). "asimov_widths_1d" uses a Gaussian proposal function, but maklng a step only along one axes at a time.
This is a pretty robust algorithm, but it usually does not perform as well as the others. Finally, "asimov_mcmc_cov" calculates the covariance matrix
by iteratively running Markov Chains on the Asimov likelihood, calculating the covariance estimate from the preceding Markov Chain, until some convergence is
reached. This used to be the default in theta up to February 2013. It usually performs well unless the number of dimensions is large (larger than 100).

``mcmc stepsize_factor`` is a scale factor for the Gaussian proposal step size to be applied. Usually, there should be no need to change it, unless you have reason
to believe the default choice has poor performance which can be fixed by adjusting the Gaussian step size with a global factor.

