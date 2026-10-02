\name{wm_wdsm_fit}
\alias{wm_wdsm_fit}
\title{Fit Weighted Double Score Matching}
\usage{
wm_wdsm_fit(Y, Z, weights, ps_design, pg0_design,
  pg1_design = NULL, M = 3L, estimand = c("PATE", "PATT"),
  inference = c("none", "full_x"), controls = list(),
  conf.level = 0.95, ps_offset = NULL,
  pg0_offset = NULL, pg1_offset = NULL,
  ps_weighting = c("probability", "unit"),
  tie_rule = c("row_order", "source_random"),
  tie_seed = NULL, tie_tolerance = 64 * .Machine$double.eps)
}
\arguments{
  \item{Y}{Finite numeric outcome vector.}
  \item{Z}{Binary treatment vector, containing both zero and one.}
  \item{weights}{Finite, strictly positive known probability weights.
  Weights may depend on treatment and covariates, or on covariates alone.
  Their statistical meaning must identify the intended population estimand.
  Estimation uncertainty in the supplied weights is not included.}
  \item{ps_design}{Finite numeric matrix for the logistic propensity
  regression, with one row per observation. Include any intended intercept
  explicitly. Column names, if supplied, must be nonempty and unique.}
  \item{ps_weighting}{\code{"probability"} (default) fits the propensity
  regression with the supplied probability weights. \code{"unit"} fits it
  without those weights, reproducing the original WDSM prospective convention.
  This option affects propensity fitting only: matching donor fractions, target
  normalization and quadratic correction still use \code{weights}.
  It does not specify a sampling design or establish a model's correctness.}
  \item{pg0_design}{Finite numeric matrix for the control prognostic regression.
  Include any intended intercept explicitly.}
  \item{pg1_design}{Treated prognostic regression matrix, required for PATE.
  It may differ from \code{pg0_design}, including its number of columns.
  Leave \code{NULL} for PATT, which fits no treated prognostic model.}
  \item{M}{Positive integer number of distinct opposite-arm donors per query.
  Donors may be reused by different queries.}
  \item{tie_rule, tie_seed, tie_tolerance}{Matching policy forwarded to
  \code{\link{wm_match}}. The default retains exact distance/row-order ties.
  \code{"source_random"} reproduces the original WDSM squared-distance
  boundary policy, requires an explicit nonnegative integer seed, and permits
  only \code{inference = "none"}. Its default tolerance is 64 times machine
  epsilon; positive tolerance can change neighborhood geometry. The dedicated
  RNG stream is restored; external Box--Muller RNG is unsupported.}
  \item{estimand}{Population average effect (PATE), the default, or treated
  population average effect (PATT).}
  \item{inference}{\code{"none"}, the default, returns only point estimation.
  \code{"full_x"} computes the conditional full-X fitted-score variance
  described below. Selection of this method does not verify its assumptions.}
  \item{controls}{Named list overriding solver controls: \code{glm_maxit = 100L},
  \code{glm_epsilon = 1e-12}, \code{refine_maxit = 200L},
  \code{refine_epsilon = 1e-14}, and \code{rank_tolerance = 1e-10}.
  Only one refinement is permitted, after a root-accuracy failure. Refinement
  must use at least as many iterations and no looser tolerance. The normalized
  score and Newton acceptance thresholds below are not weakened by controls.}
  \item{conf.level}{Confidence level strictly between zero and one. The default
  is 0.95; no interval is computed when \code{inference = "none"}.}
  \item{ps_offset}{Optional finite numeric vector of known offsets on the
  logistic linear-predictor scale. \code{NULL} means zero.}
  \item{pg0_offset, pg1_offset}{Optional finite numeric vectors of known offsets
  in the respective linear prognostic regressions. \code{NULL} means zero.
  Omit \code{pg1_offset} for PATT. All offsets are treated as fixed, and their
  estimation uncertainty is not included.}
}
\description{
Fits propensity and arm-specific prognostic regressions, constructs separate
two-dimensional matching maps for the used arms, fits the weighted quadratic
bias correction, and computes the original self-normalized matching estimator.
}
\details{
The propensity regression uses the fitting weights selected by \code{ps_weighting}.
The same selection is used in its estimating equations, Jacobian, and influence
calculation. Prognostic
regressions use ordinary least squares within their respective treatment arms.
Each raw fitted score is centered and scaled using the pooled empirical mean
and variance (divisor \eqn{n}); PATE uses its own prognostic coordinate for each
arm. The correction is weighted least squares within each arm on an intercept,
the two standardized scores, both squared scores, and their interaction.
Regressions must be identifiable, pooled score variances positive, and fitted
propensities numerically nonsaturated.

Matching uses exact Euclidean distance on these two-dimensional coordinates,
with original row index breaking exact ties by default. The optional
\code{"source_random"} policy follows donor arm zero then one and ascending
recipient rows, uniformly selecting an oversized near-boundary group using
the supplied dedicated seed. Its squared-distance arithmetic, tolerance and
point-only inference scope are specified in \code{\link{wm_match}}.
For canonical applications supply the source's actual seed (for example,
\code{20260917L + M}); this wrapper does not choose it from M.
Each donor fraction is its weight
divided by the sum of weights of that query's selected donors. This preserves
the original weighted estimator and its arm-specific maps. The number of
regression columns can vary; the matching dimension remains two.

Both modes require a normalized propensity-score equation residual no larger
than \eqn{10^{-10}} and a maximum scaled Newton step no larger than
\eqn{\min(10^{-9},1/n)}. A finite, invertible stacked Jacobian is required.
The default initial and refinement controls are the controls used by the
reviewed bounded calibration implementation. Structural failure, rank failure,
or numerical saturation is not repaired by adding penalties, clipping scores,
changing models, or retrying with another dataset. Solver warnings and attempts
are retained; an exhausted nuisance fit raises a \code{wm_wdsm_fit_error}
condition containing its attempts.

\code{inference = "none"} skips the nuisance influence matrix and covariance
calculation. It retains the reviewed prediction and root-diagnostic path.
Standard errors, variances, and confidence limits are \code{NA}.
\code{"source_random"} cannot request full-X sampling inference and does not
inherit a discrete-tie or positive-tolerance sampling theorem. Its nested
fit supports empirical fixed-reuse replicate algebra with that explicit scope.

\code{inference = "full_x"} includes the joint influence of the chosen
propensity regression, arm prognostic regressions, pooled centering and scaling,
and quadratic correction. The full row/nuisance covariance is retained; the
reported estimator variance is the root-\eqn{n} variance divided by \eqn{n}.
This method uses zero graph transport under correct full-X population correction
means. Its sampling-law premises include the bounded branch, the
finite-moment d>=2 full-X alternative, and the separate complete Gaussian
route in \code{\link{wm_fitted_inference}}.
The finite-moment d>=2 alternative retains bounded smooth correct prediction functions, the stated
baseline/current score geometry, target identification and overlap, and a
complete regular joint root with empirical-L2-consistent influences. It requires
finite \eqn{2+\delta} moments of the conditional design mean of the complete
influence and uniformly bounded conditional \eqn{2+\delta} moments of the
full-X residual and centered full influence; an outcome moment alone is
insufficient. All PS, arm-PG, pooled centering/scaling and correction equations
are included. The full nuisance covariance may be singular, but the total
limiting variance must be positive. Actual-graph row laws and statistical-root
or numerical localization conditions remain application premises.
Correctness of a prognostic model label or numerical convergence alone does not
establish these conditions. They remain unverified model premises. The method
does not establish inference for arbitrary score dimensions, arbitrary
misspecified correction means, estimated weights or offsets, or dependent
cluster/stratum sampling.

The complete Gaussian route permits unbounded supplied weights and correct
full-X predictions under exact graph determination by the full design rows and
complete Gaussian noise projection, integrated uniform current-graph moment
laws, prediction derivative/Hessian envelopes and the complete regular influence.
The known-weight Gaussian selection application additionally retains its
double-score chart/minor and propensity-interaction restrictions. A misspecified
propensity minor requires separate verification. The numerical PS/PG roots and
the original graph must satisfy the statistical-root/graph-equivalence conditions;
these are not consequences of the wrapper's score or Newton diagnostics.
No bounded branch is extended merely by assuming Gaussian outcomes.

No bootstrap is run. The original count-refit procedure and the corrected
one-step multinomial procedure are distinct methods; their finite-sample
variances need not coincide. They are not newly connected to this wrapper.
The nested fixed-map \code{fit} object should not be passed to a fixed-map
bootstrap and interpreted as inference accounting for all estimated scores.
The existing \code{wdsmatchATE} and \code{wdsmatchATT} interfaces are unchanged.
}
\value{
A list of class \code{wm_wdsm_fit}, containing \code{estimate},
\code{raw_estimate}, \code{correction}, \code{n}, \code{M}, \code{estimand}, \code{ps_weighting},
\code{se}, \code{variance}, \code{root_n_variance}, and \code{conf.int}.
The \code{inference} component records the selected method, availability,
unavailability reason, scope, and variance components when requested.
\code{fit} retains the matching graph, donor loads and predictions;
\code{nuisance} retains fitted parameters, equations and diagnostics;
\code{assembly} retains the full inference assembly or is \code{NULL} in
point-only mode. \code{solver} records controls and every attempted root fit.
\code{assumptions_verified} is always \code{FALSE} and
\code{bootstrap_requested} is always \code{FALSE}.

A full-X interval is numerically available only if the assembled variance is
finite and strictly positive; availability is conditional on the documented
premises, not evidence that those premises hold in the data.
}
\seealso{\code{\link{wm_match}}, \code{\link{wm_fit}},
  \code{\link{wdsmatchATE}}, \code{\link{wdsmatchATT}}}
