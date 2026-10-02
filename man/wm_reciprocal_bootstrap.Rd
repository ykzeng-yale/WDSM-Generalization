\name{wm_reciprocal_bootstrap}
\alias{wm_reciprocal_bootstrap}
\title{Gaussian Replication with Reciprocal-Corrected Matching Variance}
\usage{
wm_reciprocal_bootstrap(object, B = 999L, seed = NULL,
  conf.level = 0.95, interval = c("normal", "basic", "none"))
}
\arguments{
\item{object}{An available \code{wm_reciprocal_inference} result.}
\item{B}{Positive integer number of draws.}
\item{seed}{Optional nonnegative integer seed, preserving the caller's RNG
state, including its prior absence. Explicit seeds with the Box--Muller
normal generator are rejected because its hidden cache cannot be restored.}
\item{conf.level}{Confidence level strictly between zero and one.}
\item{interval}{\code{"normal"} for a normal interval using the exact conditional
variance, \code{"basic"} for simulated quantiles, or \code{"none"}.}
}
\description{Generate centered Gaussian roots from the fixed-map variance
returned by \code{wm_reciprocal_inference()}. This is scalar Gaussian
replication, separate from row multipliers and original WDSM nuisance refits.}
\details{
The root is \eqn{T^*=\sqrt{\widehat V}G} with standard normal G and the corrected
root-n variance, including any explicitly requested floor. The estimator draw
is estimate plus \eqn{T^*/\sqrt n}. No observations are rematched or refitted.

Sampling validity requires the fixed-map sampling CLT, consistent reciprocal
variance, a positive variance lower bound and a negligible root-n feasible-point
remainder. The exact conditional Gaussian variance alone does not prove sampling
validity. Actual feasible sampling-variance agreement additionally requires
mean-square transfer or sufficient uniform integrability. Generic estimated PS/PG
maps, estimated weights and nuisance-refitting replication are not covered.
The normal interval does not depend on B. Floors are inherited and reported;
none is silently added. An unavailable upstream variance is rejected.
Asymptotic validity with a floor requires a deterministic positive sequence
tending to zero, or a separate proof that it is asymptotically inactive.
An arbitrary fixed positive floor is not automatically covered.
}
\value{A list with centered \code{root_n_draws}, shifted \code{draws},
\code{conditional_root_n_variance}, \code{conditional_variance},
\code{monte_carlo_root_n_variance} and \code{monte_carlo_variance} (divisor B-1,
NA for B=1), \code{conf.int}, interval and confidence level, estimate, n, B,
seed, inherited \code{variance_floor} and \code{floor_active}, method and
\code{inference_contract}.}
\seealso{\code{\link{wm_reciprocal_inference}}, \code{\link{wm_bootstrap}},
\code{\link{wm_bootstrap_refit}}}
