\name{wm_logistic_weights}
\alias{wm_logistic_weights}
\title{Logistic Target Weights with Derivatives and Influence Estimates}
\usage{
wm_logistic_weights(Z, design, pi = 0.5, parameter_bound,
  maxit = 100L, tolerance = 1e-8, qr_tol = 1e-10)
}
\arguments{
  \item{Z}{Binary treatment vector containing both arms.}
  \item{design}{Finite numeric full-column-rank design matrix, in original
  row order. Include an intercept column explicitly if desired. The
  design map is treated as fixed and bounded in the population.}
  \item{pi}{Known desired treatment probability under the weighted target
  law, strictly between zero and one. It is not estimated from the data.}
  \item{parameter_bound}{Required positive finite scalar or vector of
  length equal to the number of design columns. Defines the declared
  compact box with endpoints plus/minus these bounds. The true parameter
  must be strictly inside this domain.}
  \item{maxit}{Maximum positive integer number of damped Newton iterations.}
  \item{tolerance}{Positive convergence tolerance below one. Both the
  empirical score and Newton update in linear-predictor units must be
  small.}
  \item{qr_tol}{Positive rank-checking tolerance below one.}
}
\description{
Fits a logistic treatment model on an explicitly bounded parameter
domain. Returns individual target weights, their analytic derivatives,
and same-sample logistic influence estimates for use with
\code{\link{wm_weight_adjust}}.
}
\details{
Writing \eqn{T_i} for a design row, the model is
\eqn{q_i=\operatorname{expit}(T_i^\top\alpha)}. The returned weights are
\eqn{W_i=\pi/q_i} for treated rows and
\eqn{W_i=(1-\pi)/(1-q_i)} for controls. Their derivative is
\eqn{\dot W_i=W_i(q_i-Z_i)T_i}. The empirical information is
\deqn{\widehat I=n^{-1}\sum_i q_i(1-q_i)T_iT_i^\top,}
and each influence estimate is
\eqn{\widehat U_i=\widehat I^{-1}T_i(Z_i-q_i)}.

The convex logistic likelihood is optimized using an analytic damped
Newton step, a stable log-likelihood calculation, and a line search
inside the declared box. A fit at its boundary does not receive the
regular interior influence formula. Singular information,
nonconvergence and numerical probability boundaries cause errors.
The Newton linear-predictor check guards against accepting a separated
fit merely because its score has become small. No probability clipping,
ridge, Firth adjustment, term deletion or automatic intercept is used.

The inference construction requires a correct logistic model, bounded
population design vector, positive definite population information and
a true parameter in the interior of the declared compact domain.
Finite observed design values do not establish a population bound.
The estimator is required to have a regular interior solution; the
routine deliberately reports other finite-sample outcomes as failures.

At the true parameter these weights preserve the observed covariate
distribution while setting the target treatment probability to
\code{pi}. Identification and matching residual-centering conditions
remain separate. If the fitted probabilities also change the matching
coordinates, this weight-only influence correction does not include
that coordinate-estimation effect.
}
\value{
A \code{wm_logistic_weights} list containing \code{weights},
\code{weight_derivative}, \code{nuisance_influence},
\code{parameters}, \code{propensity}, \code{empirical_information},
\code{score}, \code{iterations}, \code{converged}, the original
\code{design} and \code{Z}, the declared parameter bounds and numerical
controls, and an \code{assumptions} list.
Weights and derivatives use their original units.
}
\seealso{\code{\link{wm_weight_adjust}}}
\examples{
z <- c(rep(0, 14), rep(1, 6))
fit <- wm_logistic_weights(z, matrix(1, nrow = 20, ncol = 1),
                           pi = 0.5, parameter_bound = 5)
fit$parameters
colMeans(fit$nuisance_influence)
}
