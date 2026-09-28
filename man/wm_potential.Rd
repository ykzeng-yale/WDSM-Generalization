\name{wm_potential}
\alias{wm_potential}
\alias{wm_potential_fit}
\title{Stabilized Matching for One Potential-Outcome Mean}
\usage{
wm_potential(Y, Z, weights, scores, arm, M = 3L,
  mean, rho, fold_id = NULL, strata = NULL, variance = TRUE)

wm_potential_fit(Y, Z, weights, scores, arm, M = 3L,
  regression = c("polynomial", "spline"), degree = 1L,
  spline_order = NULL, mesh_exponent = NULL, moment_order = NULL,
  support = NULL, rho_bounds, fold_id = NULL, strata = NULL,
  variance = TRUE, max_basis = 10000L, max_basis_elements = 1e7,
  qr_tol = 1e-10)
}
\arguments{
\item{Y, Z, weights}{Finite outcome, binary treatment and positive individual-weight vectors.}
\item{scores}{Fixed supplied numeric matching matrix, with one row per observation. The fitting wrapper also accepts a vector for a scalar map.}
\item{arm}{Potential-outcome arm, zero or one.}
\item{M}{Fixed positive number of opposite-arm donors per query.}
\item{mean}{Supplied weighted mean predictions for the requested arm, at every row.}
\item{rho}{Supplied positive conditional-weight mean predictions for the requested arm, at every row.}
\item{fold_id, strata}{Matching restrictions. For the fitting wrapper, regression training excludes the evaluation fold and remains within its stratum. These are not survey-design variance specifications.}
\item{variance}{Whether to calculate contribution-based variance.}
\item{regression, degree, spline_order, mesh_exponent, moment_order}{Regression and explicit rate parameters, as in \code{\link{wm_fit}}.}
\item{support}{Known fixed support rectangle, as \code{support0} in \code{wm_fit}.}
\item{rho_bounds}{Declared positive lower and upper bounds strictly containing the true conditional-weight mean. Predictions are clipped to these bounds and clipping is recorded.}
\item{max_basis, max_basis_elements, qr_tol}{Basis allocation and full-rank least-squares controls, as in \code{wm_fit}.}
}
\description{Estimates a single weighted potential-outcome mean using observed outcomes in the requested arm and stabilized residual matching in the other arm. Each query uses exactly M distinct donors, with replacement across queries.}
\details{
The outer denominator is the sum of weights over all rows. The estimator imputes
\eqn{m(S_i)+M^{-1}\sum_{j\in J(i)} W_j\{Y_j-m(S_j)\}/\rho(S_j)}
for each query and retains each donor's observed outcome. Its normalized row contribution is
\eqn{[W_i\{m(S_i)-\widehat\theta\}+I(Z_i=a)\{\rho(S_i)+C_i\}\eta_i]/\overline W},
where \eqn{\eta_i=W_i\{Y_i-m(S_i)\}/\rho(S_i)} and
\eqn{C_i=M^{-1}\sum_{j:i\in J(j)}W_j}. There is no unordered-edge contribution.

The theorem requires iid rows, reduced-mean identification, weighted residual centering, positive bounded weights and overlap, fixed-dimensional continuous score geometry, finite moments of order p greater than two, and the relevant nuisance rates. The supplied predictions do not establish these assumptions. For spline fitting, declare the moment order and known support and satisfy the strict rate window; bounded conditional residual variance is additionally required for the regression rates. Stabilized splines require fixed data-independent folds with all matching inside the evaluation fold. Correct regular finite polynomial nuisance models have a separate same-sample justification. Estimated matching-map or estimated-weight uncertainty is not added automatically.

The internal fitted model labels use arm zero for the requested donor arm; the returned \code{nuisance$donor_arm} records the original treatment coding. The result's rho, loads and numerator use weights divided by the common \code{weight_scale}. \code{wm_bootstrap} uses the returned normalized row contributions directly.
}
\value{A \code{wm_potential} object inheriting from \code{wm_match}, containing the estimate, root-n and sampling variances, standard error, imputed outcomes, graph, incoming loads, contribution rows, and an exact numerator identity. The fitting wrapper also records regression models, training/evaluation rows and clipping.}
\seealso{\code{\link{wm_match}}, \code{\link{wm_fit}}, \code{\link{wm_split_pate}}, \code{\link{wm_bootstrap}}}
