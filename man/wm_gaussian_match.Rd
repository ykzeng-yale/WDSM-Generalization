\name{wm_gaussian_match}
\alias{wm_gaussian_match}
\title{Matching at an actual fitted Gaussian linear prediction parameter}
\usage{
wm_gaussian_match(Y, Z, weights, design0, design1 = design0,
  score_map, offset0 = NULL, offset1 = NULL, M = 3L,
  estimand = c("PATE", "PATT"), fold_id = NULL, strata = NULL,
  qr_tol = 1e-10)
}
\arguments{
\item{Y}{Finite numeric observed outcomes.}
\item{Z}{Binary treatment indicators containing both arms.}
\item{weights}{Fixed known positive individual weights \eqn{W=w_Z(X)} defining the target law.}
\item{design0, design1}{Numeric matrices with one row per observation and the same parameter columns, defining the two potential conditional means. Both are required for the joint outcome regression, including for PATT.}
\item{score_map}{Function of the named fitted parameter vector. Return a numeric matrix for a common matching map, or a list with exactly \code{scores0} and \code{scores1}; a NULL second map uses the first. The function's other inputs must be fixed design covariates, not outcomes.}
\item{offset0, offset1}{Known potential-mean offset vectors; NULL means zero.}
\item{M}{Fixed number of opposite-arm donors per query.}
\item{estimand}{PATE or PATT.}
\item{fold_id, strata}{Optional fixed matching restrictions, as in \code{wm_match}. These do not split the outcome regression into independent training samples.}
\item{qr_tol}{Positive rank tolerance below one; deficient designs are rejected without dropping columns.}
}
\description{
Fits a common finite-dimensional outcome parameter by ordinary least squares, evaluates the matching coordinates at the actual unrounded fit, and includes prediction-parameter uncertainty and its covariance with the matching contributions. This qualified Gaussian branch allows scalar matching coordinates.
}
\details{
The theorem requires correct full-covariate conditional means of the form offset plus design times parameter, conditionally independent Gaussian errors with variances bounded above and away from zero, bounded designs/offsets/true means, and a positive definite population Gram matrix. The fixed weights satisfy \eqn{W=w_Z(X)}: include every outcome-informative weight mark in \eqn{X}, and impose the full-covariate Gaussian model and design assumptions on that enlarged covariate vector. Conditional heteroskedasticity is permitted. Observation weights define the matching target and donor fractions; the outcome fit uses unweighted OLS.

Score maps must be uniformly Lipschitz in the fixed-dimensional parameter, with uniformly regular current donor and query densities on compact convex supports and convergence of current marked local-law expectations along deterministic parameter sequences. Score maps may depend on outcomes only through this fitted coefficient. These are substantive assumptions, not facts checked by supplying a function. Fixed known bounded positive weights, iid observations and causal identification remain required. No unrestricted fitted scalar-score result is asserted.

The nuisance influence is the inverse empirical OLS Gram matrix times observed design and fitted residual. Its matching sensitivity is the weighted query-minus-donor prediction-design imbalance. Centered contribution rows add their product, retaining the covariance term. \code{wm_bootstrap} then applies Gaussian multipliers to those corrected rows; it does not refit or rematch. A positive limiting variance is required for nondegenerate Wald inference.
}
\value{
A \code{wm_match} object with adjusted variance, standard error and contribution rows. \code{gaussian_adjustment} retains the fitted parameter, design, offsets, residuals, influence, sensitivity, covariance decomposition and explicit theorem contract. The supplied score function itself is not retained.
}
\seealso{\code{\link{wm_match}}, \code{\link{wm_prediction_adjust}}, \code{\link{wm_training_adjust}}, \code{\link{wm_bootstrap}}}
