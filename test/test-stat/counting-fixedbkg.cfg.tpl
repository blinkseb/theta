//note: this is not a valid configuration file, but a configuration file template.
//It is used by the test scripts which do some replacements, e.g., __THETA__ is replaced by a flating point value.

parameters = ("Theta", "mu");

observables = {
   o = {
      nbins = 1;
      range = [0.0, 1.0];
   };
};


flat-unit-histo = {
    type = "fixed_poly";
    observable = "o";
    normalize_to = 1.0;
    coefficients = [1.0];
};

//model definition:
counting = {
   o = {
      signal = {
          coefficient-function = { type = "mult"; parameters = ("Theta");};
          histogram = "@flat-unit-histo";
      };
      background = {
          coefficient-function = { type = "mult"; parameters = ("mu");};
          histogram = "@flat-unit-histo";
      };
   };
   
   parameter-distribution = {
      type = "flat_distribution";
      mu = {
         range = [__MU__, __MU__];
      };
      Theta = {
         //using 0.0 as lower range border would give a prediction of exactly 0 which would
         // have infinite likelihood in case of non-zero data. The delta_nll method
         // uses minimizers which cannot cope well with infinities. Therefore, this has currently
         // to be prevented by using a non-zero value here:
         range = (1E-12, "inf"); 
      	 fix-sample-value = __THETA__;
      };
   };
};

