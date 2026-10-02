# SPDX-License-Identifier: GPL-3.0-only
# Checked zero-start PS solver extracted from the supplied local revised
# WDSM case-study module. This is not an unmodified upstream-commit claim.
# WDSM/wdsmatch author credits: Yukang Zeng, Guangyu Tong, Jiaqi Tong,
# Haidong Lu, Bhramar Mukherjee, and Fan Li.
# Sourcing defines this function only; no shared tie module is sourced.

wdsm_case_fit_propensity <- function(formula, data, weights, maxit = 100L,
                              score_tolerance = 1e-7) {
  frame <- stats::model.frame(formula, data = data, na.action = stats::na.fail)
  response <- stats::model.response(frame)
  design <- stats::model.matrix(formula, data = frame)
  if (nrow(design) != nrow(data) || length(weights) != nrow(design) ||
      any(!is.finite(c(weights, design))) || any(weights < 0) ||
      !is.finite(sum(weights)) || sum(weights) <= 0 ||
      !all(response %in% c(0, 1)) ||
      length(unique(response[weights > 0])) != 2L)
    stop("Propensity fit requires finite inputs and both positive-weight treatment arms")
  fit <- do.call(stats::glm, list(formula = formula, data = data,
    family = stats::quasibinomial(link = "logit"), weights = weights,
    start = rep(0, ncol(design)), na.action = stats::na.fail,
    control = stats::glm.control(maxit = maxit, epsilon = 1e-10)))
  if (!isTRUE(fit$converged)) stop("Propensity-score fit did not converge from zero coefficients")
  if (fit$rank != ncol(design) || any(!is.finite(stats::coef(fit))))
    stop("Propensity-score fit has unidentified or nonfinite coefficients")
  probability <- as.numeric(stats::predict(fit, newdata = data, type = "response"))
  if (length(probability) != nrow(design) || any(!is.finite(probability)) ||
      any(probability <= .Machine$double.eps | probability >= 1 - .Machine$double.eps))
    stop("Propensity-score fit has nonfinite or numerically saturated probabilities")
  # Dimensionless numerical score residual, invariant to common weight scaling
  # and nonzero scaling of a design coordinate. This is an optimality check,
  # not a statistical assumption or a general separation diagnostic.
  scale <- pmax(colSums(weights * abs(design)), .Machine$double.eps * sum(weights))
  normalized_score <- drop(crossprod(design, weights * (response - probability))) / scale
  maximum_score <- max(abs(normalized_score))
  if (!is.finite(maximum_score) || maximum_score > score_tolerance)
    stop(sprintf("Propensity score residual %.3g exceeds numerical tolerance %.3g",
                 maximum_score, score_tolerance))
  attr(fit, "wdsm_ps_diagnostics") <- list(converged = TRUE,
    initialization = "zero_coefficients_input_weights", iterations = fit$iter,
    maxit = as.integer(maxit), normalized_score_max = maximum_score,
    normalized_score_tolerance = score_tolerance,
    probability_range = range(probability), coefficient_range = range(stats::coef(fit)),
    rank = fit$rank, positive_weight_count = sum(weights > 0), weight_range = range(weights))
  fit
}

