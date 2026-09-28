\name{wm_training_adjust}
\alias{wm_training_adjust}
\title{Uncertainty from Independently Trained Scores and Predictions}
\usage{
wm_training_adjust(fit, training_influence,
  mean_derivative0 = NULL, mean_derivative1 = NULL)
}
\arguments{
  \item{fit}{An unadjusted self-normalized \code{\link{wm_match}} object
  with variance estimation enabled, fitted mean predictions and retained
  original data. Its scores and prediction parameters must have been
  fitted using a sample independent of these evaluation observations.
  A \code{\link{wm_fit}} object is numerically accepted, but fitting
  regressions on the evaluation rows does not satisfy the independent
  training contract.}
  \item{training_influence}{Finite numeric \eqn{m} by \eqn{p} matrix,
  with \eqn{m\geq2}, of joint influence estimates from the independent
  training sample, normalized so that
  \deqn{\sqrt m(\widehat\theta-\theta_0)
    = m^{-1/2}\sum_{k=1}^m U_k^{\rm tr}+o_p(1).}
  Columns follow the same parameter order as the prediction derivatives.}
  \item{mean_derivative0, mean_derivative1}{Finite numeric \eqn{n} by
  \eqn{p} matrices of complete potential-outcome prediction derivatives
  on every evaluation row. Include the chain rule through any fitted
  scores or prediction basis. \code{NULL} means an identically zero
  prediction derivative. PATT uses only \code{mean_derivative0}.
  All named derivative and influence matrices must have the same
  nonempty, unique column names in the same order; unnamed columns are
  matched by position.}
}
\description{
Adds uncertainty from an independent training sample to the contribution
variance of self-normalized matching. Supports matching maps of any
fixed positive dimension under the independent-training theorem's
correct full-covariate prediction and uniform current-law assumptions.
Preserves the point estimate, evaluation contribution rows and graph.
}
\details{
The prediction sensitivity \eqn{\widehat b} uses the exact weighted
query--donor derivative imbalance described in
\code{\link{wm_prediction_adjust}}: \eqn{\widehat b_1-\widehat b_0}
for PATE and \eqn{-\widehat b_0} for PATT, with the actual effect
denominator. The scalar projection
\eqn{v_k=\widehat b^\top\widehat U_k^{\rm tr}} is centered over the
\eqn{m} training rows. The stored training multiplier coefficient is
\deqn{t_k=(n/m)(v_k-\overline v).}
The evaluation sample size remains \eqn{n}. The root-n variance becomes
\deqn{\widehat V_{\rm eval}+\frac{1}{n}\sum_{k=1}^m t_k^2
  =\widehat V_{\rm eval}
   +\frac{n}{m}\left\{\frac1m\sum_{k=1}^m(v_k-\overline v)^2\right\}.}
The influence covariance uses divisor \eqn{m}. Training and evaluation
influences are not stacked into one sample variance. There is no
same-sample cross term under the stipulated independence.

\code{\link{wm_bootstrap}} generates independent training multipliers
after the evaluation-row and edge multipliers. Its root-n draw includes
\eqn{n^{-1/2}\sum_k t_k\xi_k^{\rm tr}}, so its conditional variance
equals the stored variance formula. The point-estimator contribution
rows and \code{contributions$nuisance} are preserved.

Inference requires independent iid training and evaluation samples,
\eqn{n,m\to\infty} with \eqn{n/m} converging to a finite nonnegative
constant, a joint training influence expansion with mean-zero,
finite-variance influence, and empirically L2 consistent influence
estimates. All fitted components represented here must be trained
solely on the independent sample. Supplying a matrix with a different
number of rows does not establish independence. Ordinary cross-fitting
with overlapping training observations is a different procedure.

The outcome predictions must be correct for the full-covariate
conditional means at the true parameter. Their values and first two
parameter derivatives, and the residuals, are uniformly bounded
nearby. The target is identified under the weighted matching
assumptions. Weights are fixed known functions
\eqn{W=w_Z(X)}, bounded above and away from zero; \eqn{X} must include
any outcome-informative weight marks. Arm probabilities are positive,
dimensions and donor count are fixed, and distance ties use fixed row
order.

For every nearby deterministic parameter, each donor-score law must
have density bounded above and away from zero on a fixed compact convex
body, and its query density must be uniformly bounded and supported
there. Score maps have smooth bounded parameter derivatives.
For a convergence-determining family of bounded-Lipschitz tests of the
joint weight, prediction-derivative, mean and residual-variance marks,
the marked intensity densities must converge uniformly on compact
interior sets as the parameter approaches truth and be locally
equicontinuous. Bounded marks provide uniform integrability. A smooth
support-preserving conditional chart is one sufficient construction;
pointwise kernel continuity alone is insufficient. Fixed finite
restriction cells require these conditions within each cell and
positive limiting proportions.

These conditions are recorded in \code{training_adjustment$contract}.
The function does not verify population regularity or correct model
specification from finite arrays. Scalar maps are allowed by this
independent-training theorem; this is not a general same-sample fitted
scalar result. Estimated weights, stabilized matching, evaluation-sample
refitting, misspecified full-covariate means, arbitrary learners, and
overlapping cross-fit uncertainty are not covered.
Raw fits, point-only fits and existing nuisance adjustments are rejected.
Start from the unadjusted fit to replace a correction.
}
\value{
The original fitted object with \code{contributions$training} added
and \code{root_n_variance}, \code{variance}, and \code{se} updated.
The \code{training_adjustment} list records the sensitivity and
unsigned direction imbalances, prediction derivatives, training
influences and centered projections, evaluation and training sample
sizes and their ratio, dimensions, separate variance terms, and the
user-supplied unverified theoretical contract.
}
\seealso{\code{\link{wm_match}}, \code{\link{wm_prediction_adjust}},
  \code{\link{wm_bootstrap}}}
