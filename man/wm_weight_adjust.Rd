\name{wm_weight_adjust}
\alias{wm_weight_adjust}
\title{Account for Estimated Individual Weights on a Fixed Matching Graph}
\usage{
wm_weight_adjust(fit, weight_derivative, nuisance_influence)
}
\arguments{
  \item{fit}{A self-normalized \code{\link{wm_match}} or
  \code{\link{wm_fit}} object with variance estimation enabled, relevant
  fitted mean predictions, retained original data, and no existing
  nuisance correction.}
  \item{weight_derivative}{Finite numeric \eqn{n} by \eqn{p} matrix:
  derivatives of the original, unnormalized individual weights with
  respect to the fitted parameter. Rows must follow the original
  observations; columns follow the parameter order.}
  \item{nuisance_influence}{Finite numeric \eqn{n} by \eqn{p} matrix of
  influence estimates for the same parameter, normalized so that
  \deqn{\sqrt n(\widehat\alpha-\alpha_0)
    = n^{-1/2}\sum_i U_i + o_p(1).}
  If both matrices have column names, their names and order must agree.}
}
\description{
Differentiates the exact realized self-normalized matching estimator with
respect to individual weights, holding the matching graph and fitted
mean predictions fixed. Adds the resulting first-stage influence to the
normalized contribution rows and recomputes variance without changing
the point estimate.
}
\details{
For every selected donor, the derivative of its fraction includes the
derivatives of all weights in the matched set:
\deqn{\dot\lambda_{ij} =
  \frac{\dot W_j-\lambda_{ij}\sum_{k\in J(i)}\dot W_k}
       {\sum_{k\in J(i)}W_k}.}
The incoming derivative includes both recipient-weight and donor-fraction
derivatives. The random outer denominator is differentiated as well.
For PATT only control predictions are needed.

If \eqn{\widehat G} is the derivative of the normalized effect estimate
and \eqn{L_i} its baseline normalized contribution, the corrected row is
\eqn{L_i+\widehat G^\top\widehat U_i}. Rows are centered for the empirical
second moment. This retains the covariance between matching and weight
estimation; adding a separate nuisance variance alone omits that term.
The root-n variance is the mean squared corrected row; the variance of
the effect estimate divides it by \eqn{n}.

The function checks retained data, graph, means, point estimate and
contribution reconstruction before adjustment. It refuses an existing
nuisance correction to avoid double counting. Start from the unadjusted
fit to recompute a correction. A subsequent contribution bootstrap uses
the corrected rows already stored in the result.

Valid inference requires a proved influence expansion for the supplied
matrix, consistent influence estimates, bounded smooth positive weight
paths, the fixed-map matching assumptions, and an appropriate transfer
argument for fitted mean regressions. Supplying matrices of the right
shape does not prove those assumptions. Refitting the means during
differentiation is a different procedure. Estimated matching coordinates,
standardization, externally trained weights with different sample-size
scaling, and the stabilized donor rule require separate theory.
}
\value{
The original fitted object with its point estimate preserved and
updated \code{contributions$row}, \code{contributions$nuisance},
\code{root_n_variance}, \code{variance}, and \code{se}.
The \code{weight_adjustment} list records the sensitivity, an independent
direct-imputation derivative check, the supplied matrices, incoming-load
derivatives, their normalization scale, and the baseline, cross, and
nuisance terms in the variance decomposition.
The original \code{contributions$actual} still represents the centered
point-estimator numerator.
}
\seealso{\code{\link{wm_match}}, \code{\link{wm_fit}},
  \code{\link{wm_logistic_weights}}}
\examples{
x <- rep(seq(-1, 1, length.out = 6), 2)
z <- rep(0:1, each = 6)
weight_fit <- wm_logistic_weights(z, cbind(intercept = 1, x = x),
                                  parameter_bound = 5)
y <- 1 + x + 2 * z
fit <- wm_fit(y, z, weight_fit$weights, x, M = 1)
adjusted <- wm_weight_adjust(fit, weight_fit$weight_derivative,
                             weight_fit$nuisance_influence)
adjusted$estimate
}
