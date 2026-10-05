\name{wm_fit}
\alias{wm_fit}
\title{Fit Nuisance Regressions for Matching on Fixed Known Maps}
\usage{
wm_fit(Y, Z, weights, scores0, scores1 = NULL, M = 3L,
  estimand = c("PATE", "PATT"),
  method = c("self_normalized", "stabilized"),
  regression = c("polynomial", "spline"), degree = 1L,
  spline_order = NULL, mesh_exponent = NULL,
  moment_order = NULL, support0 = NULL, support1 = NULL,
  rho_bounds = NULL, fold_id = NULL, strata = NULL,
  variance = TRUE, max_basis = 10000L,
  max_basis_elements = 1e7, qr_tol = 1e-10,
  tie_rule = c("row_order", "source_random"),
  tie_seed = NULL, tie_tolerance = 64 * .Machine$double.eps)
}
\arguments{
  \item{Y}{Finite numeric outcome vector.}
  \item{Z}{Binary treatment vector containing both zero and one.}
  \item{weights}{Known, finite, strictly positive individual weights.}
  \item{scores0}{Finite numeric matrix of fixed known control-direction
  matching coordinates, with one row per observation. A numeric vector is
  accepted for a scalar map.}
  \item{scores1}{Treated-direction coordinates for PATE. If \code{NULL},
  use \code{scores0}. Different dimensions are allowed for self-normalized
  PATE. Stabilized PATE with its row/edge variance requires identical coordinates;
  use \code{variance = FALSE} for a distinct-map point calculation and the
  separate \code{\link{wm_reciprocal_inference}} contract for inference.
  Omit this argument for PATT.}
  \item{M}{Positive integer number of distinct donors per query. Donors may
  be reused by different queries.}
  \item{estimand}{PATE or PATT. PATT fits only the control outcome regression
  and, for stabilized matching, the control conditional-weight mean.}
  \item{method}{Self-normalized donor fractions or stabilized residual
  imputation; see \code{\link{wm_match}}.}
  \item{tie_rule, tie_seed, tie_tolerance}{Matching policy forwarded to
  \code{\link{wm_match}}. The default preserves exact distance/row-order
  matching without RNG use. \code{"source_random"} requires an explicit
  nonnegative integer seed and \code{variance = FALSE}; it has no supported
  sampling-inference theorem. The tolerance is a nonnegative finite
  squared-distance boundary multiplier, defaulting to 64 times machine epsilon.
  The dedicated RNG stream is restored; external Box--Muller RNG is unsupported.}
  \item{regression}{Complete finite polynomial or tensor-product B-spline.}
  \item{degree}{Nonnegative total polynomial degree, including all
  interactions through that degree and an intercept. Supply one number for
  all used arms, or two numbers in control, treated order for PATE. Ignored
  for splines. The positive-degree polynomial constructor supports matching
  dimension at most 100; degree zero does not use that recursion limit.
  Basis-size and allocation guards still apply.}
  \item{spline_order}{Explicit integer B-spline order at least two.
  The piecewise polynomial degree is one less than this order. One number,
  or two in control, treated order for PATE. Required for splines.}
  \item{mesh_exponent}{Explicit positive mesh exponent \eqn{a}: a training
  arm of size \eqn{N} uses \eqn{\kappa=\lceil N^a\rceil} equal subintervals
  per coordinate. One number, or two in control, treated order for PATE.
  Required for splines and checked against the sufficient rate interval.}
  \item{moment_order}{Declared residual moment order \eqn{p>2}, or
  \code{Inf} to declare sufficiently many finite moments. One number, or
  two in control, treated order for PATE. Required for splines. This is a
  caller assumption, not an estimate from observed outcomes; it does not
  remove the bounded-outcome condition for common-score stabilized PATE.
  Stabilized PATT permits finite p greater than two under its separate
  moment, nuisance-rate and variance-transfer conditions.}
  \item{support0}{Known control-map support: a numeric \eqn{d_0} by two
  matrix of lower and upper endpoints. A two-element vector is accepted
  in dimension one. Required for splines, and every query and donor score
  must lie inside it. Empirical ranges are not substituted.}
  \item{support1}{Known treated-map support for spline PATE. Defaults to
  \code{support0} only when both maps have identical coordinate values.
  Required separately for distinct maps; omit for PATT.}
  \item{rho_bounds}{Required for stabilized matching: a positive finite
  lower/upper pair, or a matrix with one such row per used arm. The caller
  declares that these bounds strictly contain the true conditional means
  \eqn{E(W\mid S,Z=z)}. Fitted normalizers are clipped to this interval,
  with affected original row indices and counts recorded. Bounds and
  fitted normalizers are expressed in the original weight units.}
  \item{fold_id}{Optional fixed, data-independent fold labels with at
  least two nonempty folds. Required for stabilized splines. Each fit uses
  the other folds within its stratum, and matching uses only donors in the
  evaluation fold and stratum. Without folds, allowed fits use the same
  sample for fitting and matching.}
  \item{strata}{Optional nonmissing labels restricting both nuisance
  fitting and matching to the same stratum. These are matching strata,
  not a specification of a dependent survey sampling design. The spline
  wrapper uses the same declared arm support in every stratum; the
  theorem requires each stratum's donor law to have positive density
  throughout that arm rectangle. Stratum-specific supports or dimensions
  are not implemented here.}
  \item{variance}{Compute the contribution-based variance in
  \code{\link{wm_match}}. The mathematical assumptions remain caller
  obligations.}
  \item{max_basis}{Maximum number of basis columns in one arm fit.}
  \item{max_basis_elements}{Guard on requested design-matrix and basis-index
  allocation, in numeric elements. Increasing this limit does not ensure
  that an allocated model will be estimable.}
  \item{qr_tol}{Tolerance for the rank-revealing QR fit, strictly between
  zero and one. Singular or zero-column designs cause errors.}
}
\description{
Fits weighted least-squares outcome regressions and, for stabilized
matching, ordinary least-squares regressions of individual weight on the
same basis. Supplies their predictions to \code{\link{wm_match}}.
Score maps and individual weights are treated as fixed and known.
}
\details{
Polynomial fits use raw supplied coordinates and every monomial through
the declared total degree. Spline regression affinely maps the declared
support to the unit cube, uses simple equally spaced interior knots and
repeats each endpoint \code{spline_order} times. These fixed regression
coordinates do not change the matching distance. Tensor basis indices,
knots and support are retained in the fitted object. No penalty, ridge
repair, term deletion, empirical support estimation or automatic tuning
is performed.

For splines let \eqn{d} denote the dimension, \eqn{r} the order and
\eqn{t=\max(0,1/2-1/d)}. The mesh exponent must lie strictly between
\eqn{L=t/(r-1)} and
\eqn{U=\min\{(1-2t)/(d+2),(1-2/p)/d\}}.
For \code{moment_order = Inf}, the second upper bound is \eqn{1/d}.
Stabilized fits also require \eqn{a>1/(4r)} and \eqn{a<1/(2d)}.
These are sufficient smoothness/moment rate restrictions, not a
data-driven guarantee of correct models or adequate sample size.
The spline theorem additionally requires the declared known support,
positive bounded arm densities and weights, appropriate smoothness and
bounded conditional residual variance.

For self-normalized inference the strong graph-field residual-centering
assumption is required. Weighted score-mean centering alone is
insufficient for random multi-donor denominators. The finite polynomial
option requires a correct regular finite-dimensional outcome model.
Same-sample stabilized polynomial inference also requires a correct
regular conditional-weight mean model and bounds that strictly contain
that mean. Common-score stabilized PATE requires bounded outcomes; a
large declared moment order is not a substitute. Stabilized PATT permits
finite p greater than two when its residual and treated-outcome moments,
sampling-equivalence rates and uniform nuisance consistency hold. For its
concrete spline procedure, the declared smoothness, strict mesh window,
positive normalizer bounds and bounded conditional residual variance remain
required. These finite-moment PATT results also justify its fitted squared-row
variance and fixed-contribution Gaussian multipliers; they do not extend
the common-score PATE theorem to unbounded outcomes.

For stabilized splines, all nuisance fitting is outside the evaluation
fold, and all matching is inside it. Splitting only regression training
while retaining a full-sample matching graph is not this procedure.
The user must supply the data-independent folds. Original input row
positions are retained as training and evaluation identifiers.
The asymptotic construction uses a fixed finite number of folds and
strata with positive limiting cell proportions.

Insufficient evaluation donors, insufficient training rows, zero basis
columns and rank deficiency cause informative errors. This finite-sample
policy does not silently complete a failed fit with a zero function.
Estimated propensity/prognostic maps, estimated individual weights,
dependent survey designs and misspecified-model influence corrections
are not supplied by this wrapper.
}
\value{
An object inheriting from \code{wm_match}, with an additional
\code{wm_fit} class, \code{fit_call}, and \code{nuisance} list.
The latter contains \code{models}, one per fitted arm/stratum/evaluation
fold; all-row \code{predictions}; declared supports, moments and bounds;
allocation controls; total normalizer clipping count; and the
\code{inference_contract}. Each model records original
\code{training_rows} and \code{evaluation_rows}, the full basis
specification, fitted coefficients, rank, QR pivot and tolerance.
Normalizer models also retain clipping bounds, row indices, counts
and the pre-clipping prediction range.
}
\seealso{\code{\link{wm_match}}}
\examples{
x <- seq(0, 1, length.out = 16)
z <- rep(0:1, 8)
y <- 1 + 3 * x + 2 * z
fit <- wm_fit(y, z, weights = 1 + x, scores0 = x, M = 1)
fit$estimate
fit$nuisance$models[[1]]$training_rows
}
