#' Weighted Double Score Matching Estimator for Population Average Treatment Effect on the Treated
#'
#' Estimates the population average treatment effect on the treated (PATT)
#' using survey-weighted double score matching. Treated units are matched to
#' controls using the control-side double score consisting of propensity
#' probability and control prognostic score. Only the missing control outcome
#' requires imputation; treated outcomes remain observed. Aggregation uses
#' treated-side survey weights and Hajek normalization.
#'
#' @details
#' PATT uses only the control-side prognostic and bias-correction regressions;
#' a treated-outcome regression is not required. Treatment identification uses
#' one-sided unconfoundedness of the control potential outcome, together with
#' the required survey identification, positivity, and regularity assumptions.
#' The conditional residual-centering qualification in \code{\link{wdsmatchATE}}
#' also applies: a correct propensity score alone does not justify random
#' within-set donor normalization when weights vary within a score level.
#'
#' Matching, score fitting, arm-specific quadratic bias correction, strict
#' numerical checks, and the behavior of supplied scores follow
#' \code{\link{wdsmatchATE}}. The PATT replication expression retains the
#' original weighted control reuse coefficients and normalizes by the
#' replicate's treated-side survey-weight total. Variance uses divisor
#' \code{boots}; intervals are centered normal Wald intervals. All requested
#' replicates must succeed, otherwise an informative error is raised.
#'
#' Supplied scores remain fixed during replication, so their external
#' estimation uncertainty is excluded. With
#' \code{use.bias.correction = FALSE}, matching discrepancies remain and the
#' bias-corrected asymptotic justification does not automatically apply.
#' The weight-only interface and individual-unit multinomial replication do
#' not provide general variance estimation for arbitrary clustered or
#' stratified survey designs. See \code{\link{wdsmatchATE}} for full details.
#'
#' @inheritParams wdsmatchATE
#'
#' @return A list of class \code{wdsmatch} with the same components as
#'   \code{\link{wdsmatchATE}}, with \code{estimate} targeting PATT and
#'   \code{estimand} identifying PATT.
#'
#' @examples
#' data(survey_obs)
#' fit <- wdsmatchATT(
#'   Y = survey_obs$Y,
#'   X = survey_obs[, c("X1", "X2", "X3", "X4", "X5", "X6")],
#'   Z = survey_obs$Z,
#'   weights = survey_obs$survey_weight,
#'   M = 3,
#'   model.ps = Z ~ X1 + X2 + X3 + X4 + X5 + X6 + X1:X2,
#'   model.pg = Y ~ X1 + X2 + X3 + X4 + X5 + X6 + X1:X2,
#'   sampling = "retrospective",
#'   varest = FALSE
#' )
#' fit
#'
#' @export
wdsmatchATT <- function(Y, X, Z, weights, M = 5,
                        ps = NULL, pg = NULL,
                        model.ps = NULL, model.pg = NULL,
                        sampling = c("retrospective", "prospective"),
                        use.bias.correction = TRUE,
                        varest = TRUE, boots = 200, alpha = 0.05,
                        tie.seed = 20260917L,
                        tie.tolerance = 64 * .Machine$double.eps) {
  cl <- match.call()
  if (missing(weights)) stop("'weights' are required for WDSM.", call. = FALSE)
  wdsm_run(Y, X, Z, weights, M, ps, pg, model.ps, model.pg,
           match.arg(sampling), use.bias.correction, varest, boots, alpha,
           estimand = "PATT", call = cl, tie_seed = tie.seed, tie_tolerance = tie.tolerance)
}
