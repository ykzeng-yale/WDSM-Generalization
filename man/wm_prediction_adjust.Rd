\name{wm_prediction_adjust}
\alias{wm_prediction_adjust}
\title{Qualified Parametric Prediction and Coordinate Uncertainty}
\usage{
wm_prediction_adjust(fit, nuisance_influence,
  mean_derivative0 = NULL, mean_derivative1 = NULL,
  changing_scores, weight_derivative = NULL)
}
\arguments{
  \item{fit}{An unadjusted self-normalized \code{\link{wm_match}} or
  \code{\link{wm_fit}} object with variance estimation enabled, fitted
  mean predictions and retained original data.}
  \item{nuisance_influence}{Finite numeric \eqn{n} by \eqn{p} matrix for
  one joint, same-sample parameter expansion
  \deqn{\sqrt n(\widehat\theta-\theta_0)
    = n^{-1/2}\sum_i U_i+o_p(1).}
  Rows follow the original observations. Include all fitted components
  in the same column order, retaining their covariance.}
  \item{mean_derivative0, mean_derivative1}{Finite numeric \eqn{n} by
  \eqn{p} matrices of the complete derivatives of the fitted potential
  outcome predictions with respect to that joint parameter, evaluated
  on every original row. Include the chain rule through any fitted
  coordinates or prediction basis. These are not just gradients with
  respect to score coordinates. \code{NULL} means an identically zero
  derivative for that prediction. PATT accepts only
  \code{mean_derivative0}.}
  \item{changing_scores}{Required logical vector indicating which
  individual matching maps depend on the fitted parameter. For PATE
  use length two, in potential-0 then potential-1 order; for PATT use
  length one for potential 0. Optional names must be \code{potential0}
  and then \code{potential1}. Fixed maps may have any positive dimension;
  changing maps must have dimension at least two.}
  \item{weight_derivative}{Optional finite numeric \eqn{n} by \eqn{p}
  derivative matrix retained for compatibility with the historical
  optional weight-adjustment route. Use \code{NULL} for the current
  supplied-weight WM framework. Nonzero derivatives invoke an excluded
  historical route and do not have the current paper's inference guarantee.
  The historical route rejects this argument if any changing map is two dimensional.
  All named derivative and influence matrices must have identical,
  nonempty, unique column names in the same order. Unnamed columns are
  matched by position.}
}
\description{
Computes prediction-imbalance sensitivity on the realized matching graph
and adds its joint first-stage influence to the normalized contribution
rows with supplied known weights. The point estimate and graph are preserved.
}
\details{
The current WM framework estimates scores, scales and predictions while
treating the supplied weights as known functions of treatment and covariates.
This helper retains a historical optional estimated-weight calculation for
compatibility; that route is outside the current paper and is not a supported
extension of its identification, variance or bootstrap results.

For direction \eqn{z}, let \eqn{D_z(i)} be the supplied prediction
derivative, \eqn{\lambda_{ij}} the self-normalized donor weight, and
\eqn{d_n} the estimator's outer weight sum. The unsigned imbalance is
\deqn{\widehat b_z =
  d_n^{-1}\sum_{i:Z_i=1-z}W_i
  \left\{D_z(i)-\sum_{j\in J_z(i)}\lambda_{ij}D_z(j)\right\}.}
For PATE the prediction sensitivity is \eqn{\widehat b_1-\widehat b_0}
and \eqn{d_n=\sum_i W_i}. For PATT it is \eqn{-\widehat b_0} and
\eqn{d_n=\sum_i Z_i W_i}. PATT needs only control predictions.

In the excluded historical route, if weight derivatives are supplied, the exact finite-graph derivative
\eqn{\widehat G} from \code{\link{wm_weight_adjust}} is computed with
means held fixed. It includes all donor-weight fraction derivatives and
the random outer denominator. The single combined sensitivity is
\eqn{\widehat g=\widehat G+\widehat b}. Corrected rows are centered
versions of \eqn{L_i+\widehat g^\top\widehat U_i}, where \eqn{L_i}
is the original normalized matching contribution. Their mean square
estimates the root-n variance. This includes matching--nuisance and
weight--prediction covariance without applying the influence twice.
A subsequent contribution bootstrap uses these corrected rows.

This is a conditional theoretical interface: finite matrices do not
verify an influence expansion, correct specification, or population
regularity. The returned \code{prediction_adjustment$contract} records
the requirements. In particular, inference uses iid observed rows,
fixed dimensions and donor count, positive bounded weights and arm
probabilities, compact convex baseline donor supports with densities
bounded above and away from zero, and bounded query densities supported
there. Weights obey \eqn{W=w_Z(X)}; if they involve additional
outcome-informative marks, those marks must be included in \eqn{X} for
the full-covariate conditional-mean requirement. Fixed finite strata or folds must satisfy the corresponding
conditions within cells and have positive limiting proportions.
Targets require the identification conditions of the matching theory.

The prediction family must be correctly specified for the
\emph{full-covariate} conditional outcome means at the true parameter.
Predictions, their first two parameter derivatives, and residuals are
uniformly bounded on a parameter neighborhood. The common fitted
parameter must have the joint expansion above, with a
\eqn{2+\delta} influence moment and empirically L2 consistent influence
estimates. A \code{wm_fit} object or derivative array alone does not
establish these requirements. For example, correctness of a regression
given a reduced score alone does not imply correctness given all
covariates.

Changing maps of dimension greater than two additionally require
uniformly bounded score derivatives and fixed-degree polynomial
squared-distance comparisons in the fixed-dimensional parameter,
including any fitted metric. Ties follow the fixed row-index rule.
Matching directly on logistic probabilities together with other
coordinates does not automatically satisfy the polynomial condition.
The excluded historical weight-adjustment route additionally requires smooth
weight paths with uniformly bounded first two parameter derivatives and
weights bounded and positive nearby.

Changing two-dimensional maps use the separate planar
\emph{distributional} transfer theorem with known weights: uniform
first-order score expansion, bounded derivatives, root-n local
displacement and local-index Lipschitz bounds, uniformly bounded current
donor densities at deterministic local parameters, and the source
centered-process modulus and marked local-law assumptions. Empirical
scales require the source bounded moment-rule and positive-variance
conditions. This branch does not assert pathwise equivalence of the
fitted-graph and baseline-graph residual sums.

Changing scalar maps, simultaneous planar score/weight estimation,
the stabilized donor rule, unrestricted misspecified outcome models,
arbitrary learners, and independent-training sample-size rescaling are
outside this helper's scope. Raw estimators, point-only fits and fits
with an existing nuisance correction are rejected. Start from the
unadjusted fit to replace a correction.
}
\value{
The original fitted object with updated \code{contributions$row},
\code{contributions$nuisance}, \code{root_n_variance}, \code{variance}
and \code{se}. The \code{prediction_adjustment} list records the combined,
prediction and weight sensitivities; unsigned direction imbalances;
supplied matrices; dimensions and theorem branches; the user-supplied,
unverified contract; and the baseline, cross and nuisance variance terms.
The original \code{contributions$actual} continues to represent the
centered point-estimator numerator.
}
\seealso{\code{\link{wm_match}}, \code{\link{wm_weight_adjust}},
  \code{\link{wm_bootstrap}}}
