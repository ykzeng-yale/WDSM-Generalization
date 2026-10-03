# Numerical representation of the same complete quadratic WLS correction.
# QR transforms correction coordinates only; matching uses the original scores.

.wm_quad_qr_allocation <- function(n, ps, pg, layout, max_stack_elements) {
  root <- .wm_model_score_allocation(n, ps, pg, layout, max_stack_elements)
  d <- as.double(layout$matching_dimensions)
  q <- (d + 1) * (d + 2) / 2
  if (any(!is.finite(q)) || any(q > .Machine$integer.max) ||
      any(n * q > .Machine$integer.max)) {
    stop("Complete quadratic QR allocation is not representable.", call. = FALSE)
  }
  retained <- n * sum(q + d + 5) + sum(q + d^2)
  workspace <- max(n * (8 * q + 6 * d + 12) + 8 * q^2 + 6 * d^2)
  planned <- root$planned_elements + retained + workspace
  if (!is.finite(planned) || planned > max_stack_elements) {
    stop("Complete score-root and quadratic QR allocation exceeds max_stack_elements; ",
         "no model or basis term is removed.", call. = FALSE)
  }
  list(parameter_dimension = root$parameter_dimension,
    raw_coordinates = root$raw_coordinates,
    basis_columns = stats::setNames(as.integer(q), names(layout$matching_dimensions)),
    score_root_elements = root$planned_elements,
    correction_retained_elements = retained,
    maximum_correction_workspace_elements = workspace,
    planned_elements = planned, max_stack_elements = max_stack_elements,
    contract = paste("Complete score-root and serial full quadratic QR workspaces.",
      "Element bound excludes inputs, R copies and solver overhead; not an RSS bound."))
}

.wm_quad_qr_arm <- function(Y, Z, W, scores, donor_arm, rank_tolerance) {
  n <- length(Y)
  d <- ncol(scores)
  q <- (as.double(d) + 1) * (as.double(d) + 2) / 2
  donor <- which(Z == donor_arm)
  if (n < d || length(donor) < q) {
    stop("Quadratic QR correction has too few donor rows for its complete basis.",
         call. = FALSE)
  }
  # LAPACK pivoting does not trim columns. Check the full triangular factor
  # explicitly because its reported rank is otherwise always the column count.
  decomposition <- qr(scores, LAPACK = TRUE)
  triangular <- qr.R(decomposition) / sqrt(n)
  reciprocal <- rcond(triangular)
  resolution <- 64 * .Machine$double.eps * max(n, d)
  if (!is.finite(reciprocal) || reciprocal <= resolution) {
    stop("Original score coordinates are singular or unresolved at machine precision; ",
         "quadratic_qr does not delete coordinates or add a ridge.", call. = FALSE)
  }
  transformed <- qr.Q(decomposition) * sqrt(n)
  colnames(transformed) <- paste0("correction_qr", seq_len(d))
  basis <- .wm_model_basis(transformed)
  coefficients <- .wm_ns_ols(basis[donor, , drop = FALSE], Y[donor],
    W[donor] / max(W[donor]), paste0("Quadratic QR correction arm ", donor_arm),
    rank_tolerance)
  mean <- as.vector(basis %*% coefficients)
  residual <- Y[donor] - mean[donor]
  w <- W[donor] / max(W[donor])
  moment <- as.vector(crossprod(basis[donor, , drop = FALSE], w * residual))
  scale <- max(1, max(abs(basis[donor, , drop = FALSE])) * sum(w * abs(residual)))
  reconstruction <- max(abs(scores[, decomposition$pivot, drop = FALSE] -
                              transformed %*% triangular)) /
    max(1, max(abs(scores)))
  orthogonality <- max(abs(crossprod(transformed) / n - diag(d)))
  if (any(!is.finite(c(mean, moment, scale, reconstruction, orthogonality)))) {
    stop("Quadratic QR prediction or diagnostic arithmetic is nonfinite.", call. = FALSE)
  }
  list(mean = mean, residual = Y - mean, donor_rows = donor,
    orthogonal_scores = transformed, basis = basis, coefficients = coefficients,
    score_pivot = decomposition$pivot, score_triangular = triangular,
    score_reciprocal_condition = reciprocal, score_resolution_guard = resolution,
    score_reconstruction_relative_error = reconstruction,
    score_orthogonality_error = orthogonality,
    normalized_correction_moment = max(abs(moment)) / scale,
    n = n, dimension = d, basis_columns = ncol(basis),
    ridge = 0, columns_removed = 0L, supplied_weight_units = "original raw W",
    assumptions_verified = FALSE,
    numerical_error_rate_verified = FALSE,
    contract = paste("Same complete probability-weighted quadratic WLS in invertible",
      "QR coordinates. Diagnostics are floating-point residuals, not a rigorous",
      "prediction-error bound or verification of a root-n numerical error rate."))
}

.wm_model_quadratic_qr_fit <- function(Y, Z, weights, ps, pg, layout,
    correction, M, estimand, inference, controls, conf.level, max_stack_elements,
    ties, call) {
  if (!is.list(correction) || is.data.frame(correction) ||
      !identical(names(correction), "method") ||
      !identical(correction$method, "quadratic_qr")) {
    stop("quadratic_qr accepts only correction=list(method='quadratic_qr'); ",
         "no ridge, rank trimming or automatic retry is supported.", call. = FALSE)
  }
  if (inference != "none") {
    stop("quadratic_qr currently requires inference='none'; its correction-projection ",
         "sampling law is not implemented by the regular full-coefficient inference.",
         call. = FALSE)
  }
  n <- length(Y)
  allocation <- .wm_quad_qr_allocation(n, ps, pg, layout, max_stack_elements)
  arguments <- list(Y = Y, Z = Z, weights = weights, ps_models = ps,
    pg0_models = pg[["0"]], pg1_models = if (estimand == "PATE") pg[["1"]] else NULL,
    estimand = estimand, max_stack_elements = max_stack_elements,
    score_only = TRUE, influence = TRUE)
  solved <- .wm_wdsm_fit_stack(arguments, .wm_model_stack, controls)
  stack <- solved$stack
  corrected <- list()
  for (arm in names(layout$matching_dimensions)) {
    corrected[[arm]] <- .wm_quad_qr_arm(stack$inputs$Y, stack$inputs$Z,
      stack$inputs$weights, stack[[paste0("scores", arm)]], as.integer(arm),
      controls$rank_tolerance)
  }
  fit <- wm_match(Y = stack$inputs$Y, Z = stack$inputs$Z,
    weights = stack$inputs$weights, scores0 = stack$scores0, scores1 = stack$scores1,
    M = M, estimand = estimand, method = "self_normalized",
    mean0 = corrected[["0"]]$mean,
    mean1 = if (estimand == "PATE") corrected[["1"]]$mean else NULL,
    variance = FALSE, tie_rule = ties$rule, tie_seed = ties$seed,
    tie_tolerance = ties$tolerance)
  reason <- paste("Quadratic QR point estimation is available; the correction",
    "projection law, its complete covariance and feasible replication are not",
    "implemented in the regular model-inference handoff.")
  inference_result <- list(method = "none", status = "point_estimate_only",
    available = FALSE, sampling_inference_available = FALSE,
    unavailable_reason = reason, root_n_variance = NA_real_,
    variance = NA_real_, se = NA_real_,
    conf.int = c(lower = NA_real_, upper = NA_real_), conf.level = conf.level,
    assumptions_verified = FALSE, original_refit_limit_agreement_declared = FALSE,
    contract = "Original quadratic-corrected point only; no sampling interval or bootstrap.")
  recipe <- list(J = length(ps), K = lengths(pg),
    matching_dimensions = layout$matching_dimensions, model_names = layout$model_names,
    shared_ps_fitted_once = TRUE, supplied_weights_known = TRUE,
    correction = "Complete weighted quadratic with invertible QR coordinates",
    legacy_double_score_reduction = FALSE, quadratic_correction_fitted = TRUE)
  structure(list(call = call, estimate = fit$estimate, raw_estimate = fit$raw_estimate,
    correction = fit$correction, n = fit$n, M = fit$M, estimand = estimand,
    ps_weighting = stack$ps_weighting, se = inference_result$se,
    variance = inference_result$variance, root_n_variance = inference_result$root_n_variance,
    conf.int = inference_result$conf.int, inference = inference_result,
    fit = fit, nuisance = stack, assembly = NULL, model_recipe = recipe,
    correction_fit = list(method = "quadratic_qr", arms = corrected,
      status = "computed", reference_derivatives_available = FALSE,
      numerical_error_rate_verified = FALSE, assumptions_verified = FALSE),
    inference_inputs = list(arguments = NULL, numerically_available = FALSE,
      score_root_available = TRUE, unavailable_reason = reason,
      assumptions_verified = FALSE),
    allocation = allocation,
    solver = list(controls = controls, attempts = solved$attempts,
      refined = solved$refined, normalized_ps_threshold = 1e-10,
      newton_threshold = min(1e-9, 1 / fit$n)),
    assumptions_verified = FALSE, bootstrap_requested = FALSE,
    scope = paste("Original known-W matching with the same complete quadratic WLS",
      "predictions in invertible correction coordinates. Original scores, metric",
      "and fractions are retained. Complete score-only IF excludes nonregular",
      "raw quadratic coefficients; it is not the full inference influence.",
      "Identification, population geometry, correction correctness, numerical",
      "error rates and the sampling/replication law remain separate premises.")),
    class = c("wm_model_fit", "list"))
}
