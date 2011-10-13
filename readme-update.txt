This file summarizes changes relevant to users, i.e., changes in config file
convention and output. It does not cover internal changes. New entries are added on
top.

within testing only:
--------------------
 * in r213:
   some options for 'cls_limits' changed because of architectural modifications and for more consistency:
     - instead of a list of producers in 'producers', specify a single producer in 'producer'
     - use 'data_source' instead of 'ts_values'; it is now no longer possible to specify a list of ts values directly
     - use 'expected_bands' instead of 'ts_values_background_bands'
     - the result table is has no 'index' and 'ts' column any more. Instead, there is the 'q' column
     See documentation of cls_limits for more details.

since r188:
-----------
 * The function plugin 'mult' has been removed; use 'multiply' instead which is much more general
 * normalize_to syntax specifying a list of doubles has been removed; use a single double instead

since r109:
-----------
* if using the 'range' setting in 'root_histogram', the range borders must coincide with bin borders,
  unless these specify underflow / overflow.
  [So far, any range setting was accepted and the bin range used was the one defined by the bins which contain the range border values,
  including these border bins. This is counter-intuitive as the range for the upper bin always extends beyond teh upper range value specified
  in the configuration file.]
* for flat_distribution, the 'fix-sample-value' setting must be supplied for parameters with an infinite range, so far thi setting was optional.
  This ensures there is a "default" for each parameter (these defaults have widespread use for bootstrapping the model, as starting
  point for minimizaiton, MCMC, etc.)
* the command line of theta changed: the name of the "main" setting cannot be specifiec any more, it is always "main".
  On the other hand, it is possible to specify more than one configuration file at once. In this case, they will be executed sequentially
  almost as if one would call theta for each of them.
  "almost" because there is one exception: the plugin files (configured via options.plugin_files) are not re-loaded between
  the runs for the configuration files. This will usually make no difference at all.
* the plugin 'model_source_norandom' has been dropped. Instead, use the 'model_source' plugin, set 'dice_poisson' to false
  and use delta distributions as 'override-parameter-distribution'.
* there used to be one global random number generator, configured via the (optional) setting 'main.seed'. Now,
  every plugin requiring random numbers has its own generator. With this approach, it is now
  possible to reproduce the exact same pseudo-data by using the same random seed in main.data_source. In earlier
  releases, a perfect pseudo-data reproduction would only work if using producers which do not consume random numbers.

  There are no required changes to the configuration files.
  It is recommended to drop 'main.seed' to avoid the warning about it being unused.
  To set the seed explicitely for plugins which use random numbers, you can add a setting 'rnd_gen' like this:
  rnd_gen = {
      seed = 123; // default of -1 means: use seed based on current (sub-second precision) time, hostname, and process id.
  };
  This affects the 'model_source', and the mcmc producers.
* The RndInfoTable schema has changed: it now contains seeds used for each module (instead of one per Run)
* The ProducerInfoTable has been dropped; it did not contain much useful information and was hardly used
* 'data-source' setting is not supported any more; use 'data_source'


