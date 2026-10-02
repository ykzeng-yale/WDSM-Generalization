\name{wm_scalar_logistic_match}
\alias{wm_scalar_logistic_match}
\title{Raw Scalar Matching After Weighted Logistic Propensity Fitting}
\usage{
wm_scalar_logistic_match(Y, Z, weights, design, M = 3L,
  estimand = c("PATE", "PATT"), weight_model = NULL,
  inference = TRUE, bandwidth = NULL, conf.level = 0.95,
  maxit = 100L, score_tolerance = 1e-7)
}
\arguments{
\item{Y}{Finite numeric observed outcomes.}
\item{Z}{Binary treatment indicators containing both arms.}
\item{weights}{Strictly positive supplied probability weights in their original units.}
\item{design}{Finite full-rank numeric matrix, with a leading column of ones
and fewer columns than observations. Every coefficient is estimated.}
\item{M}{Fixed positive integer donor count. Both arms must supply this many
donors for PATE; only the control arm must do so for PATT.}
\item{estimand}{Population average effect or population effect among the treated.}
\item{weight_model}{For inference, a list containing vectorized known functions
\code{w0}, \code{w1}, \code{dw0}, \code{dw1}, an ordered \code{domain} strictly
inside \eqn{(0,1)}, and an optional character \code{label}. Each function
returns one finite number per input probability; \code{w0} and \code{w1}
must be positive. The derivatives must be supplied, not inferred from the
observed weight vector. The weight functions and supplied weights use the
same units.}
\item{inference}{Whether to calculate the separate explicit scalar variance.
The raw matching fit is always kept separate from this calculation.}
\item{bandwidth}{Triangular-kernel bandwidth for variance-only own-arm mean
regression. \code{NULL} uses \eqn{n^{-1/5}}. A supplied value is recorded as
user supplied and does not verify an asymptotic bandwidth sequence.}
\item{conf.level}{Confidence level strictly between zero and one.}
\item{maxit, score_tolerance}{Controls for the existing ordinary weighted
logistic fitter. Its initialization uses zero coefficients and its convergence
and normalized score-residual diagnostics are retained.}
}
\description{
Fits all coefficients of a weighted logistic propensity model, calls
\code{\link{wm_match}} with the fitted probability as its scalar coordinate,
and optionally calculates an explicit model-conditional variance.
Matching uses replacement, the original normalized finite-donor weights,
and no outcome correction.}
\details{
The inference result is restricted to the bounded ordinary-logistic scalar
model. There must be a fixed nonsingular intercept-preserving representation
of the observed design as \eqn{(1,\operatorname{logit}(U),V^\top)^\top},
with bounded fixed-dimensional \eqn{V} having positive-definite covariance.
Conditional on \eqn{V}, \eqn{U} has a uniformly positive bounded density with
three uniformly bounded derivatives on a common compact interval inside
\eqn{(0,1)}. Treatment probability is \eqn{U}. Observed rows are iid under
the stated probability-selection law, with known fixed
\eqn{W=w_Z(U)}. The weight functions are positive and smooth, and
\eqn{w_z(u)-(z-u)w_z'(u)} is uniformly positive. This branch requires at
least three design columns; a lower-dimensional raw point can still be
computed without asserting its variance theorem.

Potential outcomes have smooth means \eqn{\mu_z(U)}, bounded smooth
interior-supported marks centered conditional on \eqn{U}, and bounded
centered errors independent of \eqn{(U,V,Z)}, with positive variances.
These conditions permit dependent covariates within the stated chart.
They are population assumptions, not facts certified by a successful fit
or finite sample checks. Unknown weight functions, arbitrary within-score
weight variation, estimated weights, dependent survey sampling and the
original refit bootstrap are outside this inference result.

Own-arm unweighted triangular-kernel regressions estimate the two means
only for the variance moments. They do not change the point, donor sets,
or raw fit's unavailable contribution variance. The calculation keeps
heterogeneity, matching--coefficient covariance, the finite-\eqn{M} weight
derivative terms, and the full nuisance sandwich matrices. Both PATE and
PATT divide their root-\eqn{n} variance by the full observed sample size.
Observed weights and all supplied weight functions are divided by the
same recorded constant for arithmetic; the final variance is invariant
to coherent common weight scaling.

For custom bandwidth sequences the sufficient conditions are
\eqn{h\to0}, \eqn{nh/\log n\to\infty}, and
\eqn{n^{-1/2}h^{-2}\to0}. A single supplied numeric bandwidth does not
establish these conditions. A mismatch between supplied \eqn{W_i} and
\eqn{w_{Z_i}(\widehat e_i)} is diagnostic only: estimated probabilities
need not reproduce the true weight arguments exactly.

The numerical root bridge is always reported as \code{"not_certified"}.
Ordinary solver convergence and a small floating-point score residual
do not prove equality to the exact matching graph or the theorem's
coefficient-accuracy sequence. A calculated interval is conditional on
the statistical and numerical-root requirements; no exact-root
certificate or scalar refit-bootstrap guarantee is supplied.

Malformed arguments, including an invalid weight-model record or bandwidth,
raise an error before fitting. After a point has been computed, runtime
inference failures retain that point and return an unavailable interval
with a reason. This includes callback evaluation, fitted-domain or branch-gap
failure, zero kernel denominators and singular nuisance matrices.
Nonpositive variances, floating-point underflow, nonfinite intervals and
rounded collapsed endpoints withhold the interval without flooring the
variance or changing the method. Available computed blocks and raw variance
calculations are retained. Warnings are recorded.
}
\value{
A list of class \code{wm_scalar_logistic_match}, with \code{status},
\code{estimate}, the unchanged raw \code{wm_match} object in \code{fit},
the propensity result and diagnostics, a separate \code{inference} record,
warnings, model label/domain, design-column names and bandwidth origin.
Successful results also retain the original numeric \code{fitting_design}
with an explicit \code{design_binding}, for the separate
\code{\link{wm_scalar_replication}} helper. This binding records input
identity; it does not verify population assumptions or the numerical root.
Propensity or matching failure returns a named failure record.

The inference record contains the root-\eqn{n} variance, variance divided
by \eqn{n}, standard error and interval when available; separate
\code{V0}, residual, heterogeneity, signed covariance and nuisance-quadratic
blocks; the full \code{b}, \code{C}, \code{J}, \code{B} and \code{Sigma};
variance-only predictions; and explicit scope/diagnostic fields.
The raw point fit is not relabelled as a contribution-based inference fit,
and it is not automatically accepted by a bootstrap function.
}
\seealso{\code{\link{wm_match}}, \code{\link{wm_fitted_inference}},
\code{\link{wm_scalar_replication}}}
