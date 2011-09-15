This file summarizes changes relevant to users, i.e., changes in config file
convention and output. It does not cover internal changes.


from June 2010 to trunk:
------------------------
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
      

from April 2010 to June 2010:
-----------------------------
* instead of the
  plugins = {
      files = ("some-so-file", "some-other-so-file");
  };
  setting group, you now have to use
  options = {
     plugin_files = ("some-so-file", "some-other-so-file");
  };
  The reason is that also other global options are supported now within this "options" setting group; although
  the only one so far is the n_threads parameter.
* you have to use a line
   name = "...";
  for each producer and for the data_source. This name will be used to construct a column name in
  the products table of the output.
* the run in the config file (main = {...}) now has "output_database" setting which specifies the
  database to which write the results. To emulate the behaviour, replace the
       result-file = "<filename>";
  setting by:
    output_database = {
        type = "sqlite_database";
        filename = "SigAndBkg.db";
    };
* instead of "data-source" in the run, one can now also use "data_source", with underscore, to avoid some problems
  with hyphens in column names.


from Feb 2010 to April 2010:
----------------------------
1. globally defined "parameters" setting is now only a list of names
2. what used to be "ranges" and "default" in the "parameters" setting is now at the model or producer level:
  a parameter-distribution setting in the model MUST be defined ) which contains
  a prior distribution for ALL parameters of the model (which can be flat).
  This also replaces the "constraints" setting in the model.
  (this is to minimize global information. It should make the process of pseudo data generation more transparent
  to the user, at least in more complex settings; see also 4. below.)
3. the component specification setting used to have a "coefficients" setting. For
  more consistency, this is replaced by a "coefficient-function"  setting which has
  to specify a function.
  Change "coefficients = (<list>);" to "coefficient-function = {type = "mult"; parameters = (<list>);};"
4. run now contains a data-source setting to specify where to take the data from.
  Use
  data-source = { type = "model_source"; model = "@..."; };
  for now where model is the model you want to draw the pseudo data from. The only model specified in the run
  is now simply called "model" and is now what was "model-producers" before.
  (the introduction of data-source can be used to write data-providing plugins, for example
  to run over real data or to run over special test pseudo data).
5. For deltanll_hypotest and mcmc_posterior_ratio, the signal-plus-background and background-only
  hypotheses are now specified as distributions which each have to specify ALL model parameters. Typically,
  you can use
  signal-plus-background-distribution = "@path-to-model.parameter-distribution";
  background-only-distribution = "@b-dist"; //some globally defined product_distribution
  

other changes:
* parameters used for pseudo data generation are now in the "events" table, not 
  in their own table any more
* for all producers, it is possible to specify additional likelihood functions to add
  (as additional prior, etc.)
* configuration links to nested setting paths are now possible ("@path1.path2.path3")

