\name{wm_reciprocal_inference}
\alias{wm_reciprocal_inference}
\title{Fixed-Map Reciprocal Covariance Inference}
\usage{
wm_reciprocal_inference(object, conf.level = 0.95, variance_floor = NULL)
}
\arguments{
  \item{object}{A corrected, unadjusted \code{wm_match} or \code{wm_fit} result
    for PATE or PATT. A fit constructed with \code{variance = FALSE} is allowed.
    Known fixed matching maps and individual weights are substantive premises,
    not conclusions inferred from the object's class or stored metadata.}
  \item{conf.level}{Confidence level strictly between zero and one.}
  \item{variance_floor}{Optional strictly positive floor in the internally
    normalized \emph{numerator-variance} units, before division by
    \eqn{\gamma_n^2}. The default applies no floor.}
}
\description{
Computes a separate reciprocal-covariance variance estimate and Wald interval
using the complete reindexed contribution rows of a fixed matching graph.
It neither changes the fitted point estimate nor replaces its contribution
arrays. It does not fit models, rematch observations or run a bootstrap.
}
\details{
Write \eqn{w_i=W_i/\max_j W_j}, using the fit's stored analysis weights.
Let \eqn{\gamma_n=n^{-1}\sum_i w_i} for PATE and
\eqn{\gamma_n=n^{-1}\sum_i Z_iw_i} for PATT. Stored rho predictions already
use the same normalized weight units. If \eqn{a_i} is
\code{object$contributions$actual}, the diagonal numerator estimate is
\eqn{n^{-1}\sum_i(\gamma_n a_i)^2}. These are the actual complete rows;
the stabilized fit's separate row-plus-edge decomposition is not substituted.

For each PATE reciprocal pair, orient \eqn{i} as treated and \eqn{j} as
control. Both directed matching edges must be present. With original donor
normalization, its contribution to \eqn{\widehat R_n} is
\deqn{\frac{1}{n}\frac{w_i^2 w_j^2\epsilon_i\epsilon_j}
 {D_1(j)D_0(i)},}
where \eqn{\epsilon_i=Y_i-m_{Z_i}(S_{Z_i,i})} and \eqn{D_z(q)} is the
sum of normalized weights of query \eqn{q}'s selected arm-z donors.
For stabilized matching, put
\eqn{r_i=w_i\epsilon_i/\rho_{Z_i,i}}. The pair contribution is
\deqn{\frac{1}{nM^2}(w_i r_i)(w_j r_j).}
Every pair is counted once. For PATT there is only the control-donor direction,
the reciprocal table is empty and \eqn{\widehat R_n=0}.

The signed numerator estimate is
\deqn{\widehat v_n=n^{-1}\sum_i(\gamma_n a_i)^2-2\widehat R_n.}
The estimator variance is \eqn{\widehat v_n/(n\gamma_n^2)} when positive.
The signed unregularized values are always retained. Without an explicit
floor, nonpositive values produce \code{available = FALSE} and unavailable
standard errors and intervals. With a floor, the used estimate is
\eqn{\max(\widehat v_n,\text{variance_floor})}; its activation is recorded.
Regularization does not establish calibration. For an asymptotically inactive
regularization sequence under positive variance, specify a positive floor
tending to zero in these numerator units.

The theorem concerns iid observations, fixed positive dimensions and fixed M,
bounded outcomes and weights bounded above and away from zero. Each direction
has compact convex donor support with nonempty interior, density bounded above
and away from zero, and a bounded query density supported within it. The two
fixed maps may differ and need not have a joint density. The original
self-normalized PATE result requires
\eqn{E(\epsilon\mid Z=z,S_z,W)=0} in each arm. PATT needs this condition only
for control residuals about the control reduced mean, using the direct treated
outcome row. Stabilized matching instead requires own-score weighted centering
\eqn{E(W\epsilon\mid Z=z,S_z)=0} and the appropriate positive conditional
weight mean rho. Causal identification remains a separate assumption.

Variance plug-in consistency requires uniformly consistent evaluated means,
a consistent target center, and, when used, uniformly consistent rho bounded
away from zero. Feasible Wald inference additionally requires a negligible
root-n feasible-versus-oracle remainder and nondegenerate oracle numerator
variance. Uniform prediction consistency alone does not prove that remainder.
These premises give studentization; tracking the feasible estimator's actual
sampling variance requires an additional mean-square transfer or sufficient
uniform integrability. This function does not verify any population premise.
Estimated maps, estimated weights, dependent survey designs and first-stage
influence adjustments are outside this API's contract.

Nuisance-, training-, Gaussian-score-, weight- and prediction-adjusted objects
are rejected. Unadjusted \code{wm_fit} model metadata is allowed but does not
verify rates. Graphs with multiple fold or stratum cells are rejected.
Stored dimensions, graph indices, neighbor-edge agreement, selected-edge
distances, donor fractions, loads, point reconstruction and complete rows are
checked; numerical equalities use a componentwise absolute/relative tolerance
of \eqn{10^{-10}}. The existing graph is not reranked or rematched, so this
is not a certification that an arbitrary hand-built graph selects the nearest
donors. Discrete ties are not made eligible for the continuous-score theorem
by a deterministic computational tie rule.
}
\value{
A separate object of class \code{wm_reciprocal_inference}, never
\code{wm_match}, with:
\describe{
  \item{estimate, n, M, estimand, method}{The unchanged point estimate and
    matching specification.}
  \item{gamma, weight_scale}{Target denominator mean in normalized units and
    the original common weight scale.}
  \item{numerator_variance, used_numerator_variance}{Signed unregularized and
    optionally floored numerator variance estimates.}
  \item{reciprocal_numerator, diagonal_numerator_variance}{The reciprocal
    estimate \eqn{\widehat R_n} and complete-row diagonal estimate.}
  \item{unregularized_root_n_variance, unregularized_variance}{Signed values
    on the root-n and estimator scales, retained even when nonpositive.}
  \item{root_n_variance, variance, se, conf.int}{Used root-n variance,
    estimator variance, standard error and two-sided normal interval;
    unavailable values are \code{NA}.}
  \item{conf.level, available, variance_floor, floor_active}{Interval level,
    numerical availability and explicit regularization record. Availability
    does not assert that the statistical assumptions hold.}
  \item{reciprocal_pairs}{One row per reciprocal pair: treated/control
    indices, both directed edge indices, selected donor totals, own residuals,
    signed directed residual marks, their product and its contribution divided
    by n. It is empty for PATT.}
  \item{actual_contributions}{The original normalized complete rows.}
  \item{diagnostic, assumptions_verified, inference_contract, validation}{The
    positivity diagnostic and explicit limitations. \code{assumptions_verified}
    is always \code{FALSE}.}
}
}
\seealso{\code{\link{wm_match}}, \code{\link{wm_fit}}}
