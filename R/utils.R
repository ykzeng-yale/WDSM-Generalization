# Internal numerical helpers shared by both estimands.
# Aligned with the independently verified September 2026 simulation formulas.

wdsm_fit_ps <- function(formula, data, weights, maxit = 100L,
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
  list(probability = probability, diagnostics = list(converged = TRUE,
    initialization = "zero_coefficients_raw_weights", iterations = fit$iter,
    maxit = as.integer(maxit), normalized_score_max = maximum_score,
    normalized_score_tolerance = score_tolerance,
    probability_range = range(probability), coefficient_range = range(stats::coef(fit)),
    rank = fit$rank, positive_weight_count = sum(weights > 0), weight_range = range(weights)))
}

wdsm_fit_regression <- function(formula, data, newdata, weights, label) {
  if (length(weights) != nrow(data) || any(!is.finite(weights)) ||
      any(weights < 0) || !is.finite(sum(weights)) || sum(weights) <= 0)
    stop(label, " requires finite nonnegative weights with positive total")
  fit <- do.call(stats::lm, list(formula = formula, data = data,
    weights = weights, na.action = stats::na.fail))
  if (fit$rank != length(stats::coef(fit)) || any(!is.finite(stats::coef(fit))))
    stop(label, " fit has unidentified or nonfinite coefficients")
  prediction <- as.numeric(stats::predict(fit, newdata = newdata))
  if (length(prediction) != nrow(newdata) || any(!is.finite(prediction)))
    stop(label, " produced nonfinite or incomplete predictions")
  prediction
}

wdsm_standardize <- function(scores, multiplicity) {
  scores <- as.matrix(scores)
  if (nrow(scores) != length(multiplicity) || any(!is.finite(c(scores, multiplicity))) ||
      any(multiplicity < 0) || sum(multiplicity) <= 0)
    stop("Invalid scores or multiplicities for standardization")
  center <- colSums(multiplicity * scores) / sum(multiplicity)
  centered <- sweep(scores, 2L, center, "-")
  scale <- sqrt(colSums(multiplicity * centered^2) / sum(multiplicity))
  if (any(!is.finite(scale) | scale <= 0)) stop("Degenerate matching score")
  sweep(centered, 2L, scale, "/")
}

wdsm_sieve <- function(score) {
  data.frame(ps = score[, 1L], psi = score[, 2L],
    ps_sq = score[, 1L]^2, ps_psi = score[, 1L] * score[, 2L], psi_sq = score[, 2L]^2)
}

wdsm_reuse <- function(A, weights, matches_0, matches_1 = NULL,
                             estimand = c("PATE", "PATT")) {
  estimand <- match.arg(estimand)
  n <- length(A)
  if (length(weights) != n || length(matches_0) != n ||
      !all(A %in% c(0, 1)) || any(!is.finite(weights) | weights <= 0) ||
      (estimand == "PATE" && length(matches_1) != n)) stop("Invalid original reuse inputs")
  reuse <- numeric(n)
  targets <- if (estimand == "PATE") seq_len(n) else which(A == 1)
  for (i in targets) {
    donors <- if (A[i] == 1) matches_0[[i]] else matches_1[[i]]
    if (!length(donors) || anyNA(donors) || any(donors != as.integer(donors)) ||
        any(donors < 1 | donors > n) || anyDuplicated(donors) || any(A[donors] == A[i]))
      stop("Missing or invalid opposite-arm donor set for unit ", i)
    donor_total <- sum(weights[donors])
    if (!is.finite(donor_total) || donor_total <= 0) stop("Invalid original donor weight total")
    reuse[donors] <- reuse[donors] + weights[i] * (weights[donors] / donor_total)
  }
  if (any(!is.finite(reuse))) stop("Nonfinite original weighted reuse")
  reuse
}

wdsm_replicate <- function(Y, A, weights, K, multiplicity, q0, q1 = NULL,
                                 estimand = c("PATE", "PATT")) {
  estimand <- match.arg(estimand)
  n <- length(Y)
  vectors <- list(A, weights, K, multiplicity, q0)
  if (estimand == "PATE") vectors <- c(vectors, list(q1))
  if (any(vapply(vectors, length, integer(1L)) != n) ||
      any(!is.finite(c(Y, unlist(vectors)))) || !all(A %in% c(0, 1)) ||
      any(weights <= 0) || any(K < 0) || any(multiplicity < 0))
    stop("Invalid linearized-replicate inputs")
  if (estimand == "PATE") {
    denominator <- sum(multiplicity * weights)
    q_observed <- ifelse(A == 1, q1, q0)
    numerator <- sum(multiplicity * weights * (q1 - q0) +
      multiplicity * (2*A - 1) * (weights + K) * (Y - q_observed))
  } else {
    treated <- which(A == 1)
    control <- which(A == 0)
    denominator <- sum(multiplicity[treated] * weights[treated])
    numerator <- sum(multiplicity[treated] * weights[treated] * (Y[treated] - q0[treated])) -
      sum(multiplicity[control] * K[control] * (Y[control] - q0[control]))
  }
  if (!is.finite(denominator) || denominator <= 0) stop("Invalid replicate normalization total")
  answer <- numerator / denominator
  if (!is.finite(answer)) stop("Nonfinite linearized replicate")
  answer
}
