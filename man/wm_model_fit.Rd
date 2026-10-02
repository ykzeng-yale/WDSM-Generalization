\name{wm_model_fit}
\alias{wm_model_fit}
\title{Weighted Matching with Finite Lists of Fitted Score Models}
\usage{
wm_model_fit(Y, Z, weights, ps_models = list(),
  pg0_models = list(), pg1_models = NULL, M = 3L,
  estimand = c("PATE", "PATT"), inference = c("none", "full_x"),
  controls = list(), conf.level = 0.95, max_stack_elements = 5e7,
  tie_rule = c("row_order", "source_random"), tie_seed = NULL,
  tie_tolerance = 64 * .Machine$double.eps)
}
\arguments{
  \item{Y}{Finite numeric outcome vector.}
  \item{Z}{Binary treatment vector containing both arms.}
  \item{weights}{Finite, strictly positive, supplied known probability weights.
  They may depend on treatment and covariates, or on covariates alone. Their
  interpretation must identify the intended population estimand. Estimation
  uncertainty in these weights is not included.}
  \item{ps_models}{A list of propensity model descriptors. Each descriptor is
  \code{list(design = D, offset = NULL, weighting = "probability")}, where
  \code{D} is a finite numeric matrix with one row per observation, including
  its intended intercept. The offset is a known finite vector on the logistic
  linear-predictor scale; \code{NULL} means zero. The optional fitting weight
  is \code{"probability"} or \code{"unit"}. This choice affects that PS fit
  and its estimating equation; target and donor weights remain \code{weights}.
  Each PS descriptor is fitted once and shared across used arms. An empty list
  is allowed when every used arm has a prognostic coordinate.}
  \item{pg0_models, pg1_models}{Lists of control and treated prognostic model
  descriptors, each of the form \code{list(design = D, offset = NULL)}.
  Regressions are ordinary least squares within the respective arm. Include
  intended intercepts in the designs. Known offsets are on the outcome scale.
  Each arm can have a different number of models and regression columns.
  PATT requires \code{pg1_models = NULL} and fits no treated prognostic model.
  An empty used-arm list is allowed if at least one PS coordinate is supplied.}
  \item{M}{Positive integer number of distinct opposite-arm donors per query.
  Donors may be reused across queries.}
  \item{estimand}{Population average effect, \code{"PATE"}, or treated
  population average effect, \code{"PATT"}.}
  \item{inference}{The default \code{"none"} returns point estimation only.
  \code{"full_x"} requests the conditional full-X fitted-score variance.
  Selecting it does not verify the required population assumptions.}
  \item{controls}{Supported root controls from \code{\link{wm_wdsm_fit}}:
  \code{glm_maxit}, \code{glm_epsilon}, \code{refine_maxit},
  \code{refine_epsilon} and \code{rank_tolerance}. An accuracy failure permits
  at most one refinement; structural failures are not repaired by changing models.}
  \item{conf.level}{Confidence level strictly between zero and one.}
  \item{max_stack_elements}{Positive finite cap on a conservative element
  estimate for the complete nuisance arrays. Exceeding it stops before fitting;
  no model or quadratic term is removed. It also applies to the legacy dispatch.
  Inputs, R copies, solver overhead and total resident memory are outside this cap.}
  \item{tie_rule, tie_seed, tie_tolerance}{Options passed to
  \code{\link{wm_match}}. The default breaks exact distance ties by row order.
  \code{"source_random"} requires an explicit seed and
  \code{inference = "none"}; it retains that function's point-only scope.}
}
\description{
Fits finite lists of logistic propensity and arm-specific linear prognostic
models, constructs arm-specific matching coordinates, fits the complete
weighted quadratic correction, and computes self-normalized weighted matching.
}
\details{
For J propensity models and K_z prognostic models in arm z, the matching
dimension is d_z = J + K_z. Every used arm needs at least one coordinate.
Model-list and design-column names, when supplied, must be nonempty and
unique. Unsupported descriptor fields fail. Models with equal fitted values
are not merged and dependent score coordinates are not deleted.

Propensity coordinates are fitted probabilities. Prognostic fits use ordinary
arm-specific least squares. All raw coordinates are centered and scaled using
pooled empirical means and variances with divisor n, without probability
weights. Each arm's correction uses its full intercept, linear and quadratic
basis, including every coordinate pair, with probability-weighted least
squares. Its basis size is (d_z+1)(d_z+2)/2. The full nuisance dimension also
includes all fitted regression coefficients, pooled centers and variances,
and correction coefficients. It is distinct from the matching dimensions.

Matching uses Euclidean distance on the standardized coordinates. The selected
donor fraction is its supplied weight divided by the sum of weights over that
query's M donors. The target denominator is the total supplied weight for PATE
and the treated supplied weight for PATT. This wrapper always fits a quadratic
correction, including in PS-only or PG-only configurations; it does not define
the raw propensity-score estimator by merely setting d_z to one.

When J=1 and every used K_z=1, this wrapper delegates to
\code{\link{wm_wdsm_fit}} with the same designs, offsets, PS fitting weight,
M, solver controls and ties. Underlying fit, nuisance and inference arithmetic
is unchanged. Call/class and model/allocation metadata differ. New allocation
admission applies before the dispatch.

The complete equation Jacobian includes both generated-basis and prediction
derivatives. Full-X inference retains every PS, PG, center, variance and
correction influence and all covariance blocks. Known weights have zero
parameter derivatives. Point-only mode omits influence and covariance arrays,
but retains the complete root accuracy checks. Normalized PS residuals must
not exceed 1e-10, and the maximum scaled Newton step must not exceed
min(1e-9,1/n). An exhausted root fit raises a \code{wm_wdsm_fit_error} with
its recorded attempts. Rank failures and saturated fits remain failures.

\code{"full_x"} sets graph transport to zero conditionally on correct
full-X population correction means. It uses the complete contribution and
nuisance covariance from \code{\link{wm_fitted_inference}}, with root-n
variance divided by total observed n for both estimands. Used dimensions must
either all be at least two or all equal one. Scalar maps retain their separate
fitted-graph conditions; mixed scalar/vector inference is unsupported. Point
estimation permits such mixed maps.

Independent observations, target identification, regular actual roots,
correction correctness and the applicable baseline/current geometry and
empirical row/slope conditions remain substantive premises. The bounded,
finite-moment and complete Gaussian alternatives retain their separate
requirements in \code{\link{wm_fitted_inference}}. Numerical convergence and
finite Gram rank do not verify those premises or full-dimensional score
support. A correct PS candidate alone does not establish full-X correction
correctness or automatic multiple-robust inference. The full conditional
variance calculation does not estimate a generic nonzero graph transport.

No bootstrap is run. A fixed-map bootstrap of the nested \code{fit} does not
account for this complete fitted pipeline. Known offsets and supplied weights
are held fixed; their estimation uncertainty and dependent cluster/stratum
inference are outside this wrapper's scope.
}
\value{
A list of class \code{wm_model_fit}, with point fields \code{estimate},
\code{raw_estimate}, \code{correction}, \code{n}, \code{M} and \code{estimand};
\code{se}, \code{variance}, \code{root_n_variance} and \code{conf.int}; and
nested \code{fit}, \code{nuisance}, \code{inference}, \code{assembly} and
\code{solver} records. \code{model_recipe} records model counts, matching
dimensions and shared PS use. \code{allocation} records the complete-stack
dimension and preflight element estimate. \code{assembly} is \code{NULL} in
point-only mode; interval fields are \code{NA}. Legacy specializations also
inherit \code{wm_wdsm_fit}. Numerical interval availability requires finite,
strictly positive variance. \code{assumptions_verified} and
\code{bootstrap_requested} are always \code{FALSE}.
}
\seealso{\code{\link{wm_wdsm_fit}}, \code{\link{wm_match}},
  \code{\link{wm_fit}}, \code{\link{wm_fitted_inference}}}
