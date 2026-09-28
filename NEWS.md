# wdsmatch 0.3.0.9000

* Add `wm_match()` and `wm_fit()` for fixed-count matching on supplied
  Euclidean coordinates of any fixed dimension, with PATE and PATT targets,
  self-normalized donor fractions and polynomial or tensor-spline correction.
* Add stabilized residual imputation, common-score row/edge inference,
  single potential means and independently split PATE estimation. Each option
  records its own identification, centering and nuisance-rate requirements.
* Add fixed-contribution Gaussian multipliers and qualified corrections for
  parametric prediction, estimated weights and independent training. The
  Gaussian linear-model constructor permits actual fitted scalar coordinates
  under its explicit full-covariate and geometry assumptions.
* Add exact one-dimensional overlap constants and separate numerical
  boundary-rank and alpha-only integrations, retaining covariance and
  unresolved-component diagnostics.
* Include reproducible synthetic primary, first-stage and supplemental
  studies, explicit variance benchmarks, complete replication accounting,
  and numerical and simulation Monte Carlo error summaries.
* Preserve the legacy interfaces and their existing defaults. The new
  interfaces use exactly M distinct donors per query and fixed row-index
  tie ordering; their continuous-geometry theory does not cover unrestricted
  discrete ties or general dependent survey designs.

# wdsmatch 0.2.1

* Randomize boundary-distance donor selection independently for each recipient,
  replacing deterministic original-row priority. Matching remains with
  replacement and retains exactly M distinct donors per recipient.
* Add `tie.seed` and `tie.tolerance`, with a separate random stream that preserves
  the bootstrap stream and caller RNG state. Report tie and donor-reuse diagnostics.
* Preserve fixed original matches and weighted reuse during replication.
  No-tie estimates and replication draws are unchanged.
* Require R >= 3.6.0 for the specified reproducible rejection-sampling RNG.
* This is a candidate update; the published 0.2.0 artifact is retained.

# wdsmatch 0.2.0

This release corrects the numerical implementation. Point estimates, standard
errors, and confidence intervals can differ from 0.1.1. Analyses made with
0.1.1 should be rerun before numerical results are reused.

* Match on standardized propensity probabilities and arm-specific prognostic
  scores. Replace the old logit and joint-score polynomial with a complete
  quadratic basis in each arm's own two-dimensional double score.
* Correct PATE's replication residual contributions using the original weighted
  matching reuse coefficients. Use the corresponding fixed-reuse PATT formula.
* Apply survey-weighted propensity fitting for retrospective sampling and
  unweighted propensity fitting for prospective sampling, including replication
  refits. Use finite-checked, zero-initialized logistic fits.
* Use centered normal Wald intervals based on replication variance with divisor
  B, replacing the previous percentile intervals and B-1 variance convention.
* Respect disabled bias correction in both the point estimate and replication.
  This option retains matching discrepancies and does not provide the same
  asymptotic bias-correction guarantee.
* Hold supplied propensity/prognostic scores fixed in replication; corresponding
  standard errors are conditional on those supplied scores. PATT only requires
  the control-side prognostic model.
* Reject invalid inputs and unsuccessful numerical fits explicitly. Do not
  insert zero replicates or silently discard failed replicates.
* Retain the existing public arguments and result components, and add variance,
  confidence-level and procedure metadata. Print the requested confidence level
  rather than always labeling an interval as 95%.
* Add regression checks for the corrected calculations, improve the method and
  inference documentation, and add checks on Linux, Windows, and macOS.

# wdsmatch 0.1.1

* Add missing return-value documentation for the print and summary methods for
  CRAN resubmission. No numerical estimator changes from 0.1.0.

# wdsmatch 0.1.0

* Initial package version.
* `wdsmatchATE()`: Population average treatment effect estimation.
* `wdsmatchATT()`: Population average treatment effect on the treated.
* Supports retrospective and prospective sampling designs.
* Polynomial sieve bias correction with the original 10-term logit basis.
* Linearization-based multinomial bootstrap variance estimation.
* Bundled example dataset `survey_obs`.
