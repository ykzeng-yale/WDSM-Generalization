\name{wm_match}
\alias{wm_match}
\title{Weighted Matching with Supplied Scores and Predictions}
\usage{
wm_match(Y, Z, weights, scores0, scores1 = NULL, M = 3L,
  estimand = c("PATE", "PATT"),
  method = c("self_normalized", "stabilized"),
  mean0 = NULL, mean1 = NULL, rho0 = NULL, rho1 = NULL,
  fold_id = NULL, strata = NULL, variance = TRUE,
  nuisance_influence = NULL, sensitivity = NULL,
  tie_rule = c("row_order", "source_random"),
  tie_seed = NULL, tie_tolerance = 64 * .Machine$double.eps)
}
\arguments{
  \item{Y}{Finite numeric outcome vector, with at least two rows.}
  \item{Z}{Binary treatment vector containing both 0 and 1.}
  \item{weights}{Strictly positive, finite individual weights, one per row.}
  \item{scores0}{Finite numeric matrix of matching coordinates for control donors,
    with one row per observation and at least one column.}
  \item{scores1}{Matching matrix for treated donors in PATE. Its dimension may
    differ for the self-normalized method. NULL reuses scores0. Omit for PATT.
    Stabilized PATE requires identical supplied matrices when requesting its
    row/edge variance; distinct maps are allowed with \code{variance = FALSE}
    and the separate \code{\link{wm_reciprocal_inference}} contract.}
  \item{M}{Positive integer number of distinct donors per query.}
  \item{estimand}{Population average treatment effect (PATE), or population
    average treatment effect on treated individuals (PATT). Interpretation of
    the weighted target requires separate identification assumptions.}
  \item{method}{Original within-set donor normalization, or the distinct
    stabilized residual-imputation rule.}
  \item{mean0, mean1}{Finite numeric prediction vectors of length n, evaluated
    on every row. PATE requires both or neither; PATT uses only mean0.
    No regression is fitted by this function. Omitting all relevant means
    requests a raw point estimate and requires variance = FALSE.}
  \item{rho0, rho1}{Strictly positive, finite length-n predictions of the
    within-arm conditional mean of the individual weight, in the original
    weight units. Required for every stabilized donor direction. Omit rho1
    for PATT and omit both for the self-normalized method.}
  \item{fold_id, strata}{Optional nonmissing row labels. Donors and queries
    must belong to the same fold and stratum. These are restrictions on the
    matching graph; strata does not request a survey-design variance.
    Each used cell must have at least M opposite-arm donors.}
  \item{variance}{Whether to compute contribution-based variance. TRUE requires
    all relevant supplied means and a valid nuisance/identification contract.}
  \item{nuisance_influence}{Optional finite numeric n-by-p matrix of joint
    nuisance influence values for self-normalized inference only. Supply
    together with sensitivity. The function does not establish its validity.}
  \item{sensitivity}{Finite numeric vector of length p giving the derivative
    of the already normalized effect estimate. Its product with each
    nuisance influence row is added to the normalized matching contribution,
    and the resulting rows are centered.}
  \item{tie_rule}{\code{"row_order"} (default) uses exact Euclidean distances
    and original row indices. \code{"source_random"} opts into the original
    WDSM randomized squared-distance boundary policy and requires
    \code{variance = FALSE}; sampling inference is unsupported.}
  \item{tie_seed}{Explicit nonnegative integer required for \code{"source_random"}.
    A dedicated Mersenne-Twister/Inversion/Rejection stream restores the
    caller's RNG kind and seed, including prior seed absence, on every exit.
    External Box--Muller normal RNG is rejected because its hidden cache
    cannot be restored. Ignored for \code{"row_order"}.}
  \item{tie_tolerance}{Nonnegative finite squared-distance boundary multiplier
    for \code{"source_random"}; default \eqn{64} times machine epsilon.
    The width is this value times \code{max(1, Mth squared distance)}.
    Positive tolerance can change neighborhood geometry. Ignored for
    \code{"row_order"}.}
}
\description{
An exact reference implementation of fixed-M weighted matching in any supplied
positive matching dimension. It returns the graph, incoming loads, point
estimate and the contribution arrays required for Gaussian multiplier
inference. It does not fit nuisance models or verify causal assumptions.
}
\details{
For iid observed rows with law \eqn{Q}, fixed individual weights
\eqn{W=w_Z(X)} define the target full-data law \eqn{dP/dQ=W/E_Q(W)}.
PATE is the mean causal effect under \eqn{P}; PATT is its mean among
treated individuals under \eqn{P}. Weights may depend on treatment and
covariates, with weights depending only on covariates as a special case.
Consistency, full-covariate ignorability (one-sided control-outcome
ignorability suffices for PATT), weighted outcome integrability, treatment
overlap and reduced-mean identification are substantive assumptions.
Recovering an external population with sampling or response weights requires
its recovery assumptions. Balancing weights define a target through the same
law; the label alone does not establish identification or residual centering.

The fixed-map inference theorem uses fixed M and fixed positive matching
dimensions, weights bounded above and away from zero, and positive arm
probabilities. In each direction the donor support is compact and convex
with nonempty interior; its density is bounded above and away from zero.
The query density is bounded and supported inside that donor support.
Required residual and weighted effect-heterogeneity moments have some order
\eqn{p>2}. The missing-arm means are Lipschitz, and corrected inference
requires a negligible root-n matching-bias remainder, consistent fitted
contribution rows and the stated nuisance rates. Positive limiting variance
is required for nondegenerate Wald inference. Finite input values and
successful numerical fitting do not verify these population conditions.
The row/edge variance for full-sample stabilized PATE additionally requires
bounded outcomes and a common score. Its point calculation permits distinct
maps with \code{variance = FALSE}. The separate
\code{\link{wm_reciprocal_inference}} API uses complete rows and reciprocal pairs
under its bounded fixed-map theorem, allowing distinct arm maps. Unknown
fitted maps or weights require a separately justified nuisance correction.

For this function's original PATE row variance, residuals must be centered conditional on treatment,
both matching maps and the individual weight. For original PATT, control
residuals around the control reduced mean must be centered conditional on
treatment, the control matching map and weight. The treated-side
decomposition also requires residuals around its treated conditional mean
given that matching map to be centered given weight and the same field.
Only the control mean is fitted in the PATT calculation; this does not
remove the treated-side theorem condition. A correct full-covariate mean
model is one sufficient route to the respective centering conditions.

The strata and evaluation-fold extensions use a fixed finite number of
cells with positive limiting cell and arm proportions and the preceding
conditions within each cell. Evaluation folds must be data-independent;
matching and any fold-specific inference use the same restricted graph.
Individual weights and matching strata do not specify a dependent
complex-survey variance estimator. Growing dimension, growing M and
unrestricted discrete ties are outside these inference claims.

Matching uses Euclidean distances on the supplied coordinates, with no automatic
standardization. Exact numerical distance ties are resolved by original row
index by default. Every query selects M distinct opposite-arm donors; a donor may be
reused by different queries. Calipers, variable match counts and matching
without replacement are not implemented. The continuous-score assumptions of
the statistical theory cannot be established by inspecting a finite matrix.
In particular, all-discrete coordinates with fixed-M row-priority matching
do not inherit the continuous-score variance theorem merely because ties
have a deterministic computational resolution.

The opt-in \code{"source_random"} policy computes squared coordinate
differences in supplied-column order, retaining the literal original arithmetic
in dimension two. Donor arm zero is processed before arm one, with ascending
original recipient rows and eligible donor rows; PATT uses only arm zero.
Fold/stratum restrictions preserve this global recipient order and the same
dedicated stream. Donors strictly below the Mth squared distance minus the
boundary width are retained. The remaining donors are uniformly selected
without replacement from that boundary, with a random draw only when the
boundary contains more donors than needed. Selected donors are ordered by
squared distance and original row. Public edges retain the usual query order
and report Euclidean distance (the square root of that squared distance).
Squared-distance overflow is rejected; supply a declared fixed rescaling.
Rescaling may change the positive-tolerance policy because of its absolute
\code{max(1, ...)} floor.

This optional policy supplies point calculations and, for self-normalized
matching, empirical fixed-reuse replication through \code{\link{wm_bootstrap_refit}}, with explicit unsupported
sampling-inference scope. It does not extend the continuous-score theorem
to discrete ties or positive boundary widths. Contribution, reciprocal and
fitted-score sampling inference retain their existing graph-policy guards.

For the self-normalized method, a selected donor receives fraction
\eqn{W_j/\sum_{k\in J(i)} W_k}. A missing outcome is its supplied query mean
plus the donor-fraction average of observed donor residuals. The incoming
load includes these fractions; it must not be divided by M again.

For stabilized matching, donor residuals are
\eqn{\eta_{z,j}=W_j(Y_j-m_{z,j})/\rho_{z,j}}, and a missing outcome is its query
mean plus the ordinary average of these M residuals. Its incoming load is
the sum of query weights divided by M. In this method the relevant outcome
regression is the weighted conditional mean, and rho is the conditional mean
of W. Omitting means sets them to zero only for the raw point calculation;
it does not remove the requirement to supply rho.

Original self-normalized inference requires strong residual centering in the
field determining the graphs and weights, together with the appropriate
prediction-rate and identification conditions. Weighted score centering alone
does not suffice. Stabilized PATT has a separate finite-moment result under
weighted score centering. The implemented stabilized PATE variance is the
reviewed common-score, bounded-outcome row-plus-unordered-edge formula. With
distinct supplied maps, request \code{variance = FALSE} and use the separate
\code{\link{wm_reciprocal_inference}} contract for reciprocal correction. A
diagonal-only squared-row variance is not substituted for it. Supplied
predictions and influences are reported as a user-supplied nuisance contract;
their numerical validity does not establish these statistical assumptions.

For flexible cross-fitted stabilized predictions, fit outside each evaluation
fold and supply fold_id so both queries and donors are inside that evaluation
fold. Splitting regression fits while keeping a full-sample graph is a
different procedure. Restricting by strata similarly requires the supplied
predictions and assumptions to apply within the relevant strata.

Weights are internally divided by their maximum, and supplied rho values are
divided by the same number. No weights are trimmed. Normalization that
underflows a positive input or produces a nonfinite rho value is rejected.
Common rescaling of weights, and of rho when used, preserves estimates and
normalized inference. Returned loads, denominator, rho predictions and
numerator identities are in these internal weight units.
}
\value{
An object of class \code{wm_match} containing:
\describe{
  \item{estimate, raw_estimate, correction}{The supplied-mean estimate, the
    corresponding zero-mean raw estimate, and their difference.}
  \item{root_n_variance, variance, se}{Estimated variance of the root-n error,
    sampling variance (root_n_variance/n), and standard error. All are NA
    when variance = FALSE.}
  \item{n, M, estimand, method}{Sample size and estimator specification.}
  \item{data}{Original-row outcome and treatment vectors Y and Z, retained for
    explicit finite-graph sensitivity calculations.}
  \item{numerator, denominator, gamma}{Direct weighted effect numerator,
    outer weight sum, and denominator/n, in internal weight units.}
  \item{weights, analysis_weights, weight_scale}{Original individual weights,
    weights divided by their maximum, and that maximum.}
  \item{imputed, raw_imputed}{n-by-2 potential-outcome arrays. Observed outcomes
    are retained. Unused treated potential outcomes for control rows in PATT
    are NA.}
  \item{predictions}{Supplied means and internally rescaled rho vectors.}
  \item{loads}{An n-by-2 incoming-load matrix and the conditional-mean-weight
    incoming loads used by the stabilized PATE row decomposition.}
  \item{contributions}{A list with normalized numeric vectors row, edge,
    actual and nuisance. actual contains the reindexed estimator rows;
    row and edge are the arrays used for inference. edge is empty except
    for stabilized PATE. Its row and edge arrays must not be centered
    separately. nuisance is the centered optional correction.}
  \item{numerator_identity}{Direct and independently reindexed numerators,
    their numerical difference, and centered actual and row-plus-edge sums.
    Equalities hold up to floating-point arithmetic. Optional nuisance
    corrections do not change the point-estimator identity.}
  \item{graph}{Directed edges with query, donor, arm, distance, share and
    outcome_share; per-query neighbors; unordered linked pairs with both
    direction indicators for stabilized PATE; cell/fold/stratum labels;
    supplied score matrices; and transformation and tie-rule records.
    The \code{"source_random"} graph additionally contains
    \code{tie_diagnostics} with the policy, seed, tolerance, stream order,
    fixed-original-graph flag and per-arm boundary/tie/reuse counts.
    share is the within-set weight fraction or 1/M. outcome_share multiplies
    the donor outcome directly and includes W/rho for stabilized matching.}
  \item{info}{Inference status, mathematical contract and weight-unit metadata.}
}
}
\examples{
Y <- c(1, 3, 8, 4, 7, 13)
Z <- c(0, 0, 0, 1, 1, 1)
W <- c(1, 2, 4, 3, 5, 7)
S <- matrix(c(0, 1, 2, 0.5, 1.5, 3), ncol = 1)
raw <- wm_match(Y, Z, W, S, M = 2, variance = FALSE)
raw$estimate
raw$graph$neighbors
}
