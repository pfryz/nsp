# nsp
Software to accompany P. Fryzlewicz (2020) "Narrowest Significance Pursuit: inference for multiple change-points in linear models".

The main functions are
- NSP_for_Github_faster.R (the actual software), and
- NSP_simulations_and_data_examples.R (as in the title).

To use NSP_for_Github_faster.R, do the following:

- Install the R package [lpSolve](https://CRAN.R-project.org/package=lpSolve).
- Save wiener_holder_norms.txt to your R working directory.
- Source both utils.R and NSP_for_Github_faster.R into R.
- Read the descriptions within NSP_for_Github_faster.R.

This content will be turned into an R package and posted on CRAN at the earliest opportunity.

NOTE: To reproduce the results from the paper exactly, use the older and slower code in NSP_for_Github_slow.R.

Questions/comments? - p.fryzlewicz@lse.ac.uk (Piotr Fryzlewicz)
