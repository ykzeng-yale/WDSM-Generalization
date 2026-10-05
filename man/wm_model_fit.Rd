\name{wm_model_fit}
\alias{wm_model_fit}
\title{Weighted Matching with Finite Lists of Fitted Score Models}
\usage{
wm_model_fit(Y, Z, weights, ps_models = list(),
  pg0_models = list(), pg1_models = NULL, M = 3L,
  estimand = c("PATE", "PATT"), inference = c("none", "full_x"),
  controls = list(), conf.level = 0.95, max_stack_elements = 5e7,
  tie_rule = c("row_order", "source_random"), tie_seed = NULL,
  tie_tolerance = 64 * .Machine$double.eps, correction = NULL)
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
  Selecting it does not verify the required population assumptions.
  Both opt-in correction routes currently require \code{"none"}.
  They retain the complete score-root influence; the local-polynomial route
  additionally provides reference-derivative inputs.}
  \item{controls}{Supported root controls from \code{\link{wm_wdsm_fit}}:
  \code{glm_maxit}, \code{glm_epsilon}, \code{refine_maxit},
  \code{refine_epsilon} and \code{rank_tolerance}. An accuracy failure permits
  at most one refinement; structural failures are not repaired by changing models.}
  \item{conf.level}{Confidence level strictly between zero and one.}
  \item{max_stack_elements}{Positive finite cap on a conservative element
  estimate for the complete nuisance arrays. Exceeding it stops before fitting;
  no model or basis term is removed. It also applies to the legacy dispatch.
  For local polynomials this cap covers the complete score-only root, retained
  correction/derivative arrays and the largest serial local Gram workspace.
  For quadratic QR it covers the complete score-only root, retained correction
  bases and coordinate transforms, and serial regression workspaces.
  Inputs, R copies, solver overhead and total resident memory are outside this cap.}
  \item{correction}{NULL preserves the default complete quadratic correction
  and legacy dispatch. Alternatively, supply a list containing
  \code{method = "local_polynomial"}, \code{degree} and
  \code{bandwidth_exponent}, with optional \code{guard_exponent} (default 2).
  Degree and bandwidth exponent have length one, shared across used arms, or
  length two in control/treated order for PATE; PATT uses one. Degree is an
  integer at least one. The exponent must lie in the sufficient range below.
  The fixed positive guard exponent specifies \eqn{c_n=n^{-g}}. No arbitrary
  kernel, ridge, fold, adaptive bandwidth or basis-subset option is accepted.
  Alternatively, \code{list(method = "quadratic_qr")} uses invertible QR
  correction coordinates for the same complete quadratic WLS. It accepts no
  additional correction fields and requires \code{inference = "none"}.}
  \item{tie_rule, tie_seed, tie_tolerance}{Options passed to
  \code{\link{wm_match}}. The default breaks exact distance ties by row order.
  \code{"source_random"} requires an explicit seed and
  \code{inference = "none"}; it retains that function's point-only scope.}
}
\description{
Fits finite lists of logistic propensity and arm-specific linear prognostic
models, constructs arm-specific matching coordinates, fits the complete
weighted quadratic correction by default or an explicitly requested local
polynomial correction, and computes self-normalized weighted matching.
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
weights. With the default correction, each arm uses its full intercept, linear and quadratic
basis, including every coordinate pair, with probability-weighted least
squares. Its basis size is (d_z+1)(d_z+2)/2. The full nuisance dimension also
includes all fitted regression coefficients, pooled centers and variances,
and correction coefficients. It is distinct from the matching dimensions.

Matching uses Euclidean distance on the standardized coordinates. The selected
donor fraction is its supplied weight divided by the sum of weights over that
query's M donors. The target denominator is the total supplied weight for PATE
and the treated supplied weight for PATT. By default this wrapper fits a quadratic
correction, including in PS-only or PG-only configurations; it does not define
the raw propensity-score estimator by merely setting d_z to one.

With the default correction, when J=1 and every used K_z=1, this wrapper delegates to
\code{\link{wm_wdsm_fit}} with the same designs, offsets, PS fitting weight,
M, solver controls and ties. Underlying fit, nuisance and inference arithmetic
is unchanged. Call/class and model/allocation metadata differ. New allocation
admission applies before the dispatch.

For \code{quadratic_qr}, all original matching scores, distances, donor
fractions and target denominators are retained. QR transforms correction
coordinates only. An invertible linear coordinate change preserves the
complete degree-two polynomial space and hence its exact WLS predictions.
No correction column is removed. The complete score-only root excludes raw
quadratic coefficients; this matters when distinct fitted model outputs
coincide at their probability limit and those coefficients have no regular
full-rank population Jacobian.

The full score QR is rejected when its reciprocal condition is no greater
than \eqn{64\epsilon\max(n,d_z)}, with machine precision \eqn{\epsilon}.
The donor regression retains \code{controls$rank_tolerance}. Unresolved
rank, too few donor rows or nonfinite arithmetic stops without a ridge,
column deletion or bandwidth/model retry. Retained reconstruction,
orthogonality and normal-equation residual diagnostics are floating-point
checks, not a rigorous prediction-error bound or a verified root-n error
rate. Without explicit \code{quadratic_control},
\code{\link{wm_model_inference}} returns unavailable regular full-coefficient
inference for this option. Its separate \code{nested_contrast} preparation
retains the correction-projection contribution and complete joint covariance
under the scoped bounded outcome-correct nested model and numerical premises.
Pass that preparation's inference object to \code{\link{wm_bootstrap}} for its
untruncated Gaussian-mixture replication. The original point, supplied weights
and matching graph are retained; covariance-only local graphs vary with the
nuisance index. The route supports basic or no intervals, and any failed
requested slot makes all aggregate variance, standard-error and interval fields
unavailable. Numerical diagnostics do not verify the population or numerical
error-rate premises. This route does not substitute the regular
full-coefficient formula or establish original fixed-reuse bootstrap validity.

For the default quadratic route, the complete equation Jacobian includes both generated-basis and prediction
derivatives. Full-X inference retains every PS, PG, center, variance and
correction influence and all covariance blocks. Known weights have zero
parameter derivatives. Default-correction point-only mode omits influence and covariance arrays,
but retains the complete root accuracy checks. Normalized PS residuals must
not exceed 1e-10, and the maximum scaled Newton step must not exceed
min(1e-9,1/n). An exhausted root fit raises a \code{wm_wdsm_fit_error} with
its recorded attempts. Rank failures and saturated fits remain failures.

Every fitted propensity coordinate must lie strictly between
\code{.Machine$double.eps} and \code{1 - .Machine$double.eps}. This is a
fixed numerical acceptance rule, not statistical uniform overlap. It can
reject a finite regular root; one unbounded Gaussian propensity index is
enough for its pass probability to tend to zero along a consistent root
sequence. Exact-root statistical results therefore do not establish eventual
availability of this finite-data double-precision calculation. A numerical
realization of an unbounded result additionally needs a compatible precision,
threshold and root/graph-error regime. No probability clipping or silent
dataset/count-draw replacement repairs a failure. PS-free model lists do not
invoke this particular guard and retain their other numerical conditions.

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

With \code{correction = list(method = "local_polynomial", degree = q,
bandwidth_exponent = alpha, guard_exponent = g)}, the same score fits and
pooled standardization are used, but the quadratic correction is never fitted.
The complete finite score parameter contains every PS/PG coefficient and pooled
center/variance, with all Jacobian and influence blocks. Local regression
coefficients are auxiliary and do not enter this root. Its full influence is
retained even with \code{inference = "none"}; no false full-X declaration is
needed to obtain it.

For each used arm, the same full arm donor sample defines one local regression
function, evaluated at all original fitted scores. Let
\eqn{u_i=(\widehat S_{z,i}-s)/h}, let \eqn{r_q(u)} contain all monomials through
degree q with constant first, and use the fixed kernel
\deqn{K_d(u)=C_d(1-\|u\|^2)^3_+,\qquad
 C_d=\Gamma(d/2+4)/\{6\pi^{d/2}\}.}
The Gram and moment are
\deqn{G(s)=(nh^d)^{-1}\sum_{i:Z_i=z}W_i K_d(u_i)r_q(u_i)r_q(u_i)^{\mathsf T},
 \quad b(s)=(nh^d)^{-1}\sum_{i:Z_i=z}W_i K_d(u_i)r_q(u_i)Y_i.}
All weights remain in their original supplied units and n is the full observed
sample size. If \eqn{\lambda_{\min}(G)\ge c_n}, the correction is the intercept
of \eqn{G^{-1}b}; otherwise its declared pointwise completion is zero.
There is no ridge, term deletion or change of the requested bandwidth. The
absolute guard is affected by rescaling W or the kernel, even though a
successful unfloored intercept is invariant to common positive rescaling.
The kernel amplitude is fixed and disclosed, not estimated from the data.

The bandwidth is exactly \eqn{h=n^{-\alpha}}. For d<=2, require q>=1 and
\eqn{0<\alpha<1/4}. For d>2, require integer \eqn{q>(d^2-4)/4} and
\deqn{(d-2)/(2dq)<\alpha<2/\{d(d+2)\}.}
The guard uses a fixed g>0 and \eqn{c_n=n^{-g}}, default g=2. These are explicit
sufficient choices; no numerical success infers the required smoothness,
weighted conditional-mean identity, compact convex donor support, bounded
query density, uniform root-n score displacement or other population premises.
In particular, this constructive result requires bounded outcomes and known
positive weights bounded above and away from zero. It does not automatically
inherit the separate unbounded Gaussian/full-X inference alternatives.

The returned spatial gradient differentiates the fitted intercept itself:
\deqn{\partial_k\widehat m(s)=e_0^{\mathsf T}G(s)^{-1}
  \{\partial_k b(s)-\partial_kG(s)G(s)^{-1}b(s)\}.}
Both kernel and shifted monomial derivatives are included. Selecting a local
linear coefficient and dividing by h would be a different calculation.
The actual full-p score Jacobian B is contracted with this gradient to estimate
the reference derivative \eqn{B^{\mathsf T}\nabla m(S_z^0)}. This is not the
full derivative obtained by refitting the local regression after moving all
its donor coordinates. Under the stated local-polynomial and baseline matching
conditions, the constructive result supplies reference-derivative consistency
and correction point transfer. Using those inputs for inference additionally
requires the separate finite-edge, variance and fitted-law conditions.

Below the guard, the completed level is zero and the inference derivative is
NA, with retained failure/eigenvalue diagnostics and a summary warning. At
exact guard equality the level is computed but its derivative is withheld with
a boundary flag. Numerical eigenvalue/solve/overflow failures remain errors;
they do not trigger another degree, bandwidth or hidden regularization. Guard
success at finitely many evaluated rows does not verify the required uniform
shrinking-tube event. A nonfinite B/gradient contraction on an otherwise usable
row is separately marked by \code{reference_derivative_failure}; that row's
reference derivative is set to NA. Such arithmetic failure does not change the
LP level, point estimate or score-root influence. Required unavailable
derivatives suppress the optional common inference-input arguments while
preserving the completed point and the score root.

The local-polynomial option computes no variance or interval. When its
reference inputs are numerically available and the tie rule is supported,
\code{inference_inputs$arguments} can be combined with the caller's justified
\code{covariance_scope}, directional transport specifications and any required
critical controls in \code{\link{wm_fitted_inference}}. Actual matching-coordinate
arrays are separately in \code{inference_inputs$score_tangents}, in the complete
score-parameter order; their coordinate axes match the retained fit.
This does not establish generic multiple-robust inference, a general scalar
fitted-score theorem, bootstrap-refit validity or nominal coverage. Scalar
local-polynomial point correction and its transfer result alone do not extend
the common fitted-inference interface's scalar branches. The source_random
policy remains point-only and does not receive inference arguments.

No bootstrap is run. A fixed-map bootstrap of the nested \code{fit} does not
account for this complete fitted pipeline. Pass the whole retained wrapper to
\code{wm_bootstrap}: contribution replication uses its declared full-X or
explicitly prepared public inference, and \code{method="fixed_reuse"} with
counts and \code{refit=NULL} refits the default complete quadratic pipeline
with the same full-n counts, weights and controls while retaining the original
graph and reuse. Local-polynomial or quadratic-QR correction requires an
explicit compatible refit callback for that method. Preparation or successful
replication does not verify any sampling/refit-law premise. Known offsets and supplied weights
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

With local-polynomial correction, \code{nuisance} is the complete score-only
root, \code{correction_fit} retains controls, per-arm levels, gradients,
reference derivatives, eigenvalues, guard-failure/boundary and reference-
derivative-failure indicators and counts, and \code{inference_inputs} retains complete caller-use arguments and
actual score tangent arrays when numerically usable. Its numerical availability
is separate from verified assumptions or sampling inference. The top-level
\code{correction} still denotes the point estimate's correction amount.
The new correction path does not inherit \code{wm_wdsm_fit} or silently invoke
its quadratic correction. Default outputs and legacy dispatch are unchanged.

}
\seealso{\code{\link{wm_wdsm_fit}}, \code{\link{wm_match}},
  \code{\link{wm_fit}}, \code{\link{wm_bootstrap}}, \code{\link{wm_fitted_inference}}}
