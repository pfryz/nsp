# nsp
Software to accompany P. Fryzlewicz (2020) "Narrowest Significance Pursuit: inference for multiple change-points in linear models" and P. Fryzlewicz (2021) "Robust Narrowest Significance Pursuit: inference for multiple change-points in the median"

(A) "Narrowest Significance Pursuit: inference for multiple change-points in linear models"

**NOTE: the code for part (A) has now been turned into an R package called "nsp" and made available on CRAN [here](https://CRAN.R-project.org/package=nsp).**

**The entire content of this repository will be removed once (B) has also been incorporated into "nsp".**

**Please do not use the code described in part (A) here. Use the ["nsp" package](https://CRAN.R-project.org/package=nsp) instead.**

**(Legacy description starts here.)** The main files are
- NSP_for_Github_v*.R (the actual software), and
- NSP_simulations_and_data_examples_v*.R (as in the title).

To use NSP_for_Github_v*.R, do the following:

- Install the R package [lpSolve](https://CRAN.R-project.org/package=lpSolve).
- Save wiener_holder_norms.txt, tight_mres_norms_const.RData, tight_mres_norms_lin.RData to your R working directory.
- Source NSP_for_Github_v*.R into R.
- Read the descriptions within NSP_for_Github_v*.R.

Requires R package lpSolve. **(Legacy description ends here.)**

(B) "Robust Narrowest Significance Pursuit: inference for multiple change-points in the median"

The main file is
- RNSP_for_Github_v*.R (the actual software)

Note: some functions required in RNSP_for_Github_v*.R are defined in NSP_for_Github_v*.R.

Requires R package plyr.


This content will be turned into an R package and posted on CRAN at the earliest opportunity.

Questions/comments? - p.fryzlewicz@lse.ac.uk (Piotr Fryzlewicz)
