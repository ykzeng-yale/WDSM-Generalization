#' Weighted Double Score Matching Estimator for Population Average Treatment Effect
#'
#' Estimates the population average treatment effect (PATE) using weighted
#' double score matching with survey design weights. Matching uses propensity
#' probabilities and arm-specific prognostic scores. Missing potential outcomes
#' are imputed by survey-weighted averaging within match sets and aggregated
#' using Hajek normalization. Optional polynomial bias correction adjusts for
#' matching discrepancies. Inference uses fixed-reuse multinomial replication
#' with centered normal Wald intervals.
#'
#' @details
#' Under retrospective sampling, the propensity score is fitted with survey
#' weights; under prospective sampling it is fitted without survey weights.
#' Prognostic regressions are unweighted within each treatment arm in both
#' designs. These modeling choices require the corresponding identification
#' and sampling assumptions; weighting alone does not establish those assumptions.
#'
#' Matching uses unweighted pooled standardization of propensity probabilities
#' and the relevant prognostic score. The bias-correction regression uses survey
#' weights and a complete quadratic basis in that arm's own double score:
#' an intercept, both coordinates, their squares, and their interaction.
#' Distances are Euclidean and matching is with replacement. At the Mth-distance
#' boundary, the required donors are sampled uniformly without replacement,
#' separately for every recipient, while retaining all strictly closer donors.
#' A separate, reproducible random-number stream controls this selection and
#' leaves the caller's random-number state unchanged. Thus tie randomization
#' does not change the bootstrap count stream. Original matches remain fixed
#' in replication. Numerical boundary equivalence is controlled by
#' \code{tie.tolerance}; use zero for exact floating-point equality only.
#' This rule avoids reusing one fixed first-M subset throughout a tied score
#' cell; it does not establish inference guarantees for discrete scores.
#' Propensity fitting starts from zero
#' coefficients and checks convergence and finite, nonsaturated fitted probabilities.
#'
#' Multinomial replication draws counts for the original sample units and keeps
#' the original match sets and weighted reuse coefficients fixed. With bias
#' correction enabled, internally estimated nuisance scores are refitted using
#' the counts, with additional
#' survey weighting for retrospective propensity fitting and outcome bias
#' correction. Supplied scores remain fixed. Thus, when scores are supplied,
#' the resulting inference is conditional on them and excludes uncertainty
#' from their external estimation. The corresponding model formula is not used
#' to replace a supplied score during replication.
#'
#' With \code{B = boots} successful replicates, the variance estimate is
#' \code{mean((boot.estimates - mean(boot.estimates))^2)}, using divisor B.
#' The interval is \code{estimate + c(-1, 1) * qnorm(1 - alpha/2) * se}.
#' All requested replicates must succeed. An invalid fit, degenerate score,
#' unsuccessful replicate, or nonfinite calculation raises an error; failed
#' draws are not omitted or replaced by zeros. A successful numerical check
#' does not establish adequate statistical overlap.
#'
#' With \code{use.bias.correction = FALSE}, point estimation uses direct
#' matched imputation without polynomial outcome adjustment, and replication
#' uses the fixed-reuse expression without a bias-correction regression.
#' Unused nuisance scores are not refitted for these raw fixed-match replicates.
#' Matching discrepancies remain, so the asymptotic justification for the
#' bias-corrected estimator does not automatically apply to this option.
#'
#' Consistency claims for this self-normalized fixed-M estimator require
#' residual centering conditional on the individual weight and all matching-score
#' fields, in addition to causal identification, positivity, and the relevant
#' score-model and regularity assumptions. A correct propensity score alone does
#' not supply this centering when weights vary among subjects with the same
#' score. Weighted score centering alone is insufficient. These conditions do
#' not guarantee negligible finite-sample bias or nominal interval coverage in
#' every setting. This interface accepts
#' unit-level weights, not survey strata, clusters, or design-specific replicate
#' weights. The individual-unit multinomial procedure is not a general variance
#' estimator for arbitrary complex survey designs.
#'
#' @param Y Finite numeric vector of observed outcomes, with no missing values.
#' @param X Numeric matrix or data frame with finite numeric covariate columns,
#'   one row per unit, and unique names other than \code{Y} or \code{Z}.
#' @param Z Binary treatment indicator (1 = treated, 0 = control), with one
#'   entry per unit and both groups represented.
#' @param weights Finite, strictly positive numeric vector of survey design
#'   weights, with one entry per unit. Required.
#' @param M Positive integer number of nearest neighbors (default 5). Each
#'   required donor arm must contain at least \code{M} units.
#' @param ps Optional finite numeric vector of pre-estimated propensity
#'   probabilities strictly between zero and one. If \code{NULL} (default),
#'   estimate internally using \code{model.ps}. Supplied scores remain fixed
#'   during replication; see Details.
#' @param pg Optional matrix of pre-estimated prognostic scores, with columns
#'   \code{psi0} (control) and \code{psi1} (treated) in that order and one row
#'   per unit. PATT also accepts a one-column matrix containing only the
#'   control-side score; with two columns, only the first is used.
#'   If \code{NULL} (default),
#'   estimate the required arm-specific scores using \code{model.pg}. Supplied
#'   scores remain fixed during replication; see Details.
#' @param model.ps Formula for the propensity model, such as
#'   \code{Z ~ X1 + X2}. If \code{NULL}, use all columns of \code{X}.
#'   Used only when \code{ps} is \code{NULL}.
#' @param model.pg Formula for the prognostic regression, such as
#'   \code{Y ~ X1 + X2}. If \code{NULL}, use all columns of \code{X}.
#'   Used only when \code{pg} is \code{NULL}.
#' @param sampling Character: \code{"retrospective"} (default) for
#'   treatment-dependent sampling with survey-weighted propensity fitting, or
#'   \code{"prospective"} for treatment-independent sampling with unweighted
#'   propensity fitting. This choice also applies to replication refits.
#' @param use.bias.correction Logical: apply the arm-specific complete
#'   quadratic bias correction (default \code{TRUE}); see Details for the
#'   interpretation of \code{FALSE}.
#' @param varest Logical: compute replication variance and a centered normal
#'   Wald interval (default \code{TRUE}).
#' @param boots Integer number of multinomial replicates, at least 2 when
#'   variance is requested (default 200).
#' @param alpha Significance level strictly between zero and one (default
#'   0.05). The interval confidence level is \code{1 - alpha}.
#'
#' @param tie.seed Nonnegative integer seed for recipient-specific random
#'   tie selection (default 20260917), independent of the bootstrap seed.
#'   Use \code{set.seed()} separately to reproduce bootstrap draws.
#' @param tie.tolerance Nonnegative relative tolerance on squared distances.
#'   The default is \code{64 * .Machine$double.eps}, scaled by the larger
#'   of one and the Mth squared distance. Set zero to randomize exact ties only.
#'
#' @return A list of class \code{wdsmatch} with components:
#'   \item{estimate}{Point estimate of PATE.}
#'   \item{variance}{Replication variance with divisor B, or \code{NA} when
#'     \code{varest = FALSE}.}
#'   \item{se}{Replication standard error, or \code{NA} without inference.}
#'   \item{ci}{Centered normal Wald interval as \code{c(lower, upper)}, or
#'     two \code{NA} values without inference.}
#'   \item{boot.estimates}{All requested finite replication estimates, or
#'     \code{NULL} without inference.}
#'   \item{alpha}{Requested significance level.}
#'   \item{interval.type}{Interval-construction metadata.}
#'   \item{variance.divisor}{Replication-variance divisor metadata.}
#'   \item{n.boot}{Number of retained replicates.}
#'   \item{estimand}{Target estimand.}
#'   \item{sampling}{Sampling-design option.}
#'   \item{settings}{Estimation and replication settings.}
#'   \item{M}{Number of matches used.}
#'   \item{n}{Sample size.}
#'   \item{n.treated}{Number of treated units.}
#'   \item{n.control}{Number of control units.}
#'   \item{call}{The matched call.}
#'   \item{diagnostics}{Point-representation agreement and numerical
#'     propensity-fit diagnostics for the point estimate and replicates, plus
#'     tie counts, seed, tolerance, and donor-reuse diagnostics.}
#'
#' @examples
#' data(survey_obs)
#' fit <- wdsmatchATE(
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
wdsmatchATE <- function(Y, X, Z, weights, M = 5,
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
           estimand = "PATE", call = cl, tie_seed = tie.seed, tie_tolerance = tie.tolerance)
}
