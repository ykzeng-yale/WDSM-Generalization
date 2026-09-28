\name{wm_bootstrap}
\alias{wm_bootstrap}
\title{Gaussian Contribution Bootstrap for Weighted Matching}
\usage{
wm_bootstrap(object, B = 999L, seed = NULL, conf.level = 0.95,
  interval = c("normal", "basic", "none"), chunk_size = 65536L)
}
\arguments{
\item{object}{A \code{wm_match} result with variance estimation enabled.}
\item{B}{Positive integer number of multiplier draws.}
\item{seed}{Optional nonnegative integer seed. An explicit seed restores the
caller's \code{.Random.seed}, including its prior absence. The Box--Muller
normal generator does not support this option because its hidden cached state
cannot be restored. Use \code{seed = NULL} with that generator or select another
normal generator before calling this function.}
\item{conf.level}{Confidence level strictly between zero and one.}
\item{interval}{\code{"normal"} uses the exact conditional Gaussian variance;
\code{"basic"} uses simulated root-n quantiles; \code{"none"} omits the interval.}
\item{chunk_size}{Positive integer maximum number of multiplier values
generated together. No contributions-by-B matrix is allocated.}
}
\description{
Applies independent Gaussian multipliers to normalized row contributions and,
where present, unordered-edge and independent-training contributions.
Graphs, predictions, weights,
and nuisance influence corrections remain fixed.
}
\details{
For normalized row contributions \eqn{r_i}, unordered-edge contributions
\eqn{e_j}, and optional independent-training contributions \eqn{t_a}, a root-n draw is
\deqn{T^* = (\sum_i r_i \xi_i + \sum_j e_j \zeta_j +
              \sum_a t_a \omega_a)/\sqrt{n},}
with mutually independent standard-normal multipliers. Its conditional variance
is \eqn{(\sum_i r_i^2 + \sum_j e_j^2 + \sum_a t_a^2)/n}.
The value \eqn{n} is always the evaluation sample size. Training coefficients
must already contain the required sample-size scaling: with \eqn{m} independent
training rows, the training adjustment supplies \eqn{n/m} times the centered
scalar training influence. A missing or \code{NULL}
\code{contributions$training} is an empty vector. Components are not centered separately.
The fitted contributions already contain any supplied nuisance influence
correction. The procedure inherits the fitted estimator's inference assumptions;
it does not refit nuisance models or rematch observations.

The normal interval uses the exact conditional variance, not Monte Carlo
quantiles, so increasing B does not change it. Its sampling validity still
depends on the fitted estimator's asymptotic inference contract. Monte Carlo
variance estimates use the usual sample variance with denominator B minus one.
They are unavailable when B equals one.

With \code{seed = NULL}, draws advance the caller's RNG stream. Within each
replicate, multipliers are drawn in row, edge, then training order. An absent
training component preserves the existing draw sequence. Changing \code{chunk_size}
preserves that order, apart from floating-point summation differences.
}
\value{
A list with \code{root_n_draws}, centered multiplier draws on the root-n scale;
\code{draws}, equal to the fitted estimate plus those draws divided by
\code{sqrt(n)}; exact \code{conditional_root_n_variance} and
\code{conditional_variance}; Monte Carlo sample variances
\code{monte_carlo_root_n_variance} and \code{monte_carlo_variance}; and
\code{conf.int}, \code{interval}, \code{conf.level}, \code{B}, \code{n},
\code{seed}, \code{estimate}, and \code{method}. The Monte Carlo variances are
\code{NA} for B equal to one. The interval is \code{NULL} for
\code{interval = "none"}.
}
\seealso{\code{\link{wm_match}}}
