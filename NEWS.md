# wdsmatch 0.3.0

* Accept retained `wm_wdsm_fit` and `wm_model_fit` results in `wm_bootstrap`.
  Contribution replication binds declared complete public fitted inference.
  Default quadratic fixed-reuse replication refits the complete retained
  PS/PG/scaling/correction pipeline under each supplied count column, preserving
  the original point, graph, known weights, reuse and target denominator.
  Nondefault correction methods require an explicit compatible callback.
* Permit the protected `wm_model_inference` post-fit handoff for WDSM wrappers.
  Reject mismatched prepared source fits and noncanonical labeled count order;
  failed columns remain failures and are never deleted or redrawn.
* Add independent weighted estimator/reduction and full fitted-wrapper
  integration tests, an executable HTML tutorial, cross-platform package-check
  CI and package citation metadata. Existing point/analytic variance routines
  and frozen study results remain unchanged.
* Reorganize the landing README and study/application indexes. Document the
  principal producer's frozen source checkout and repair ECLS loading/count
  prerequisites. Historical source/check records retain their original scope.

# wdsmatch 0.3.0.9000

* Mark angular geometry components with zero empirical Monte Carlo variance
  as unresolved, including components with positive hits. Preserve their
  numerical estimates, MCSEs and covariance. Synchronize raw scalar inference
  references with Supplement S11.

* Add the common prepared-input, supplied-refit and WM application driver,
  deterministic selected39 ECLS decoding/construction and three case schemas.
  Complete saved-bundle reproduction retains all72 rows,12600 available draws
  and all nine unavailable ECLS PATE intervals; no participant/count/fitted
  objects are distributed. Add compact complete original-study reports with
  all608 cells, requested denominators, source reductions and adverse findings.

* Add the unchanged saved-graph application balance helper and explicit
  target-specific diagnostic documentation; no participant inputs or new
  package API are exported. Record the later focused source-tie validation
  separately from historical package checks and qualified formal results.
* Restrict the current Weighted Matching framework to supplied, known
  probability weights. Only scores, coordinate scales and outcome corrections
  are estimated. Earlier optional estimated-weight APIs and study artifacts
  are historical compatibility material, outside the current WM theory.

* Add `wm_scalar_replication()` for full-slope current-graph contribution
  variance and supplied-count/multiplier algebra under the qualified scalar
  model. It preserves raw points and full matching/nuisance covariance, uses
  total observed n for both targets, and generates no random draws or refits.
  Scalar fit objects now retain their original design for input binding.
  The row and explicit moment variances need not agree in finite samples.
* Add optional `wm_scalar_root_certificate()` for verified finite-data
  logistic-root containment, donor-set identity and point-arithmetic bounds.
  It leaves the fit unchanged and makes no asymptotic certification-success,
  variance-arithmetic or bootstrap claim. Exact rational enclosures use the
  optional `gmp` and `Rmpfr` dependencies.
* Add `wm_scalar_logistic_match()` for raw fitted-propensity matching and a
  separate explicit scalar variance under known smooth score-dependent weights
  and the documented bounded outcome/design model. Both arms' mean variation,
  finite-M drift and full covariance matrices are retained. Numerical-root
  certification and original-refit bootstrap validity are not supplied.
  Deterministic validation is available in `validation/check_scalar_logistic_match.R`.
* Add `wm_graph_transport()` and `wm_fitted_inference()` for the documented
  transport calculation and complete fitted variance in the qualified
  higher-dimensional branches. They preserve explicit application assumptions,
  quadrature-only error reporting, covariance blocks and unavailable intervals.
* Add the completed bounded iid WDSM calibration records and deterministic
  summary replay, retaining paired subset comparisons, Monte Carlo errors,
  small-sample departures and the historical outer-monitor receipt limitation.
  A separately labelled forward-replay configuration preserves frozen source
  pins after the fitted interface additions.
* Add `wm_wdsm_fit()` for the original fitted arm-specific double-score
  pipeline with arbitrary fixed M, known probability weights, and PATE/PATT.
  Point estimation is the default; explicit `inference = "full_x"` includes
  the complete fitted-pipeline correction under its documented conditions.
  PATT fits only the control prognostic model. This interface does not provide
  generic learned-score, estimated-weight or dependent-design inference.
* Add explicit `ps_weighting = "probability"` or `"unit"` to the fitted
  double-score interface, preserving the original retrospective/prospective
  PS-fitting conventions and their complete estimating equations. Paired
  original-source checks cover 96 point configurations across both designs;
  eight CorCor/M3 pairings also reproduce 64 same-count refit draws. This is
  implementation verification, not a simulation performance claim.
* Add separate fixed-map reciprocal variance and Gaussian replication APIs.
  These retain signed reciprocal covariance for distinct arm maps, preserve
  point estimates, and explicitly report unavailable or floored variances.
  Stabilized distinct-map PATE point calculations are available with
  `variance = FALSE`; the existing row/edge variance still requires a common map.
* Add `wm_bootstrap_refit()` for exact fixed-reuse multinomial replication
  with supplied counts and optional nuisance-prediction refitting, preserving
  the original WDSM normalization and variance convention. This procedure is
  separate from fixed-prediction Gaussian multipliers.
* Add independent original-survey WDSM and unit-weight matching reduction
  checks. Document that point-estimator reductions do not remove the
  joined-score centering or fitted-score inference requirements.
* Add `wm_match()` and `wm_fit()` for fixed-count matching on supplied
  Euclidean coordinates of any fixed dimension, with PATE and PATT targets,
  self-normalized donor fractions and polynomial or tensor-spline correction.
* Add stabilized residual imputation, common-score row/edge inference,
  single potential means and independently split PATE estimation. Each option
  records its own identification, centering and nuisance-rate requirements.
* Add fixed-contribution Gaussian multipliers and qualified corrections for
  known-weight parametric prediction and independent training. Historical
  estimated-weight helpers are retained separately for compatibility. The
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
