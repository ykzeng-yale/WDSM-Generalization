# Observable local-index inputs for the original quadratic correction law.
# This preparation/evaluation path computes no sampling interval or bootstrap.

.wm_qc_failure <- function(stage, reason, contrast = NULL) {
  list(status = paste0("quadratic_contrast_", stage, "_unavailable"),
    available = FALSE, numerically_available = FALSE,
    failure_stage = stage, unavailable_reason = reason, contrast = contrast,
    conditional_mean = NA_real_, conditional_variance = NA_real_,
    numerical_error_rate_verified = FALSE,
    sampling_inference_available = FALSE, assumptions_verified = FALSE)
}

# Full-column numerical QR. The resolution check is in floating-point units,
# not a statistical Gram threshold. A failure never removes a column.
.wm_qc_qr <- function(x) {
  if (!is.matrix(x) || !is.numeric(x) || ncol(x) < 1L ||
      nrow(x) < ncol(x) || any(!is.finite(x))) {
    stop("Complete QR input is nonfinite or has insufficient rows.", call. = FALSE)
  }
  factor <- qr(x, LAPACK = TRUE)
  R <- qr.R(factor)
  if (any(!is.finite(R))) stop("Complete QR factor is nonfinite.", call. = FALSE)
  reciprocal <- rcond(R)
  resolution <- 64 * .Machine$double.eps * max(dim(x))
  if (!is.finite(reciprocal) || reciprocal <= resolution) {
    stop("Complete QR rank is singular or unresolved at machine precision; no columns are removed.",
         call. = FALSE)
  }
  Q <- qr.Q(factor)
  reconstruction <- max(abs(x[, factor$pivot, drop = FALSE] - Q %*% R)) /
    max(1, max(abs(x)))
  orthogonality <- max(abs(crossprod(Q) - diag(ncol(x))))
  if (any(!is.finite(c(Q, reconstruction, orthogonality)))) {
    stop("Complete QR diagnostics are nonfinite.", call. = FALSE)
  }
  list(Q = Q, R = R, pivot = factor$pivot,
    reciprocal_condition = reciprocal, resolution = resolution,
    reconstruction_relative_error = reconstruction,
    orthogonality_error = orthogonality,
    numerical_error_rate_verified = FALSE)
}

# Exact finite logistic divided difference, evaluated without subtracting two
# probabilities or taking the h=0 directional limit for a nonzero contrast.
.wm_qc_logistic_contrast <- function(baseline, direction, step) {
  delta <- step * direction
  perturbed <- baseline + delta
  if (!is.finite(step) || step <= 0 ||
      any(!is.finite(c(baseline, direction, delta, perturbed))) ||
      any(delta == 0 & direction != 0)) {
    stop("Finite logistic path is nonfinite or its nonzero increment is unrepresentable.",
         call. = FALSE)
  }
  small <- stats::plogis(baseline)
  large <- stats::plogis(perturbed)
  divided <- numeric(length(delta))
  positive <- delta > 0
  negative <- delta < 0
  # For delta>0, expit(a+delta)-expit(a)
  # = expit(a+delta)*expit(-a)*[-expm1(-delta)].
  divided[positive] <- (direction[positive] *
    (-expm1(-delta[positive]) / delta[positive])) *
    large[positive] * stats::plogis(-baseline[positive])
  # For delta<0 use the corresponding sign-safe identity.
  divided[negative] <- (direction[negative] *
    (expm1(delta[negative]) / delta[negative])) *
    small[negative] * stats::plogis(-perturbed[negative])
  if (any(!is.finite(c(small, large, divided))) ||
      any(divided == 0 & direction != 0)) {
    stop("Finite logistic contrast is nonfinite or below representable precision.",
         call. = FALSE)
  }
  list(small_probability = small, large_probability = large,
    divided_probability = divided, logit_increment = delta)
}

.wm_qc_prepare <- function(object, control, maximum) {
  if (is.list(control) && "guard_exponent" %in% names(control)) {
    stop("guard_exponent is no longer supported; nested contrast evaluation uses no statistical Gram threshold.",
         call. = FALSE)
  }
  allowed <- c("method", "small_ps", "large_ps",
    "maximum_draw_entries", "maximum_distance_evaluations")
  required <- c("method", "small_ps", "large_ps")
  if (!is.list(control) || is.data.frame(control) || is.null(names(control)) ||
      anyNA(names(control)) || anyDuplicated(names(control)) ||
      !all(required %in% names(control)) || any(!names(control) %in% allowed) ||
      !identical(control$method, "nested_contrast")) {
    stop("quadratic_control requires method='nested_contrast', small_ps, large_ps and optional allocation caps only.",
         call. = FALSE)
  }
  replication <- list(maximum_draw_entries = 2000000,
    maximum_distance_evaluations = 100000000)
  for (name in names(replication)) {
    if (name %in% names(control)) replication[[name]] <-
      .wm_numeric_vector(control[[name]], 1L, name, positive = TRUE)
  }
  s <- object$nuisance; fit <- object$fit
  if (!inherits(s, ".wm_model_nuisance_stack") || !isTRUE(s$score_only) ||
      !identical(object$correction_fit$method, "quadratic_qr")) {
    stop("nested_contrast requires the retained quadratic_qr point and complete score-only root.",
         call. = FALSE)
  }
  n <- .wm_fit_integer(fit$n, "fit$n", 2L)
  parameters <- s$parameter_names
  p <- length(parameters)
  if (!is.character(parameters) || p < 1L || anyNA(parameters) ||
      any(!nzchar(parameters)) || anyDuplicated(parameters) ||
      !identical(names(s$parameter), parameters)) {
    stop("The complete retained parameter names/order are required.", call. = FALSE)
  }
  arms <- if (fit$estimand == "PATE") c("0", "1") else "0"
  ps <- s$inputs$ps_models; pg <- s$inputs$pg_models
  if (length(ps) != 2L || !identical(names(pg), arms) ||
      any(lengths(pg) != 1L) || any(s$layout$matching_dimensions != 3L)) {
    stop("nested_contrast currently requires exactly two retained PS models and one PG per used arm; no model is removed.",
         call. = FALSE)
  }
  select <- function(name, label) {
    if (!is.character(name) || length(name) != 1L || is.na(name) ||
        !name %in% names(ps)) stop(label, " must name a retained PS model.", call. = FALSE)
    match(name, names(ps))
  }
  small <- select(control$small_ps, "small_ps")
  large <- select(control$large_ps, "large_ps")
  if (small == large) stop("small_ps and large_ps must be different models.", call. = FALSE)
  if (any(vapply(ps, function(x) !identical(x$weighting, "probability"), logical(1L)))) {
    stop("Both nested PS models must use supplied probability weights.", call. = FALSE)
  }
  if (!identical(ps[[small]]$offset, ps[[large]]$offset)) {
    stop("The nested PS models require the same fixed known offset.", call. = FALSE)
  }
  D <- ps[[small]]$design; large_design <- ps[[large]]$design
  common <- match(colnames(D), colnames(large_design))
  if (anyNA(common) || !identical(unname(D), unname(large_design[, common, drop = FALSE]))) {
    stop("Small-model columns must occur identically in the declared large design.", call. = FALSE)
  }
  added <- setdiff(seq_len(ncol(large_design)), common)
  if (!length(added)) stop("The large PS requires a nonempty added coefficient block.", call. = FALSE)
  b <- length(added); q <- 10L
  # Full influence, projection/design workspaces, scores and fixed-M edge storage.
  planned <- n * (4 * p + 8 * b + 12 * q + 24 + 8 * length(arms) * fit$M) +
    4 * p^2 + 8 * q^2 + 8 * b^2 + 4 * ncol(D)^2
  if (any(!is.finite(c(planned, n * p, n * q, n * b))) ||
      any(c(n * p, n * q, n * b) > .Machine$integer.max) || planned > maximum) {
    return(.wm_qc_failure("allocation", paste("Complete contrast-input allocation exceeds",
      "max_influence_elements; no model, parameter or basis column is removed.")))
  }
  for (field in c("Y", "Z")) .wm_reciprocal_agree(s$inputs[[field]], fit$data[[field]], paste("Quadratic input", field))
  .wm_reciprocal_agree(s$inputs$weights, fit$weights, "Quadratic input supplied W")
  if (!identical(s$estimand, fit$estimand) ||
      any(.wm_numeric_vector(s$multiplicity, n, "multiplicity") != 1) ||
      any(.wm_numeric_vector(s$empirical_probability, n, "empirical probability") != 1 / n)) {
    stop("Original unit multiplicities and total-n empirical probabilities are required.", call. = FALSE)
  }
  phi <- .wm_fi_matrix(s$nuisance_influence, n, parameters, "Complete score influence")
  wd <- .wm_fi_matrix(s$weight_derivative, n, parameters, "Known-weight derivative")
  if (any(wd != 0)) stop("Only supplied known W with zero weight derivatives is supported.", call. = FALSE)
  if (!is.matrix(s$jacobian) || !identical(dim(s$jacobian), c(p, p)) ||
      !identical(rownames(s$jacobian), parameters) || !identical(colnames(s$jacobian), parameters)) {
    stop("The complete named retained Jacobian is required.", call. = FALSE)
  }
  psi <- .wm_fi_matrix(s$estimating_equations, n, parameters, "Complete retained equations")
  solved <- tryCatch(-t(solve(s$jacobian, t(psi))), error = identity)
  if (inherits(solved, "error")) return(.wm_qc_failure("influence", conditionMessage(solved)))
  if (any(!is.finite(solved))) return(.wm_qc_failure("influence", "Retained influence solve is nonfinite."))
  .wm_reciprocal_agree(phi, solved, "Retained equation/influence binding")
  phi <- sweep(phi, 2L, colMeans(phi), "-")
  # Scale before accumulation, as in the original retained score covariance.
  full_covariance <- crossprod(phi, phi / n)
  if (any(!is.finite(c(phi, full_covariance)))) {
    return(.wm_qc_failure("score_covariance", "Complete score covariance is nonfinite."))
  }
  Y <- s$inputs$Y; Z <- as.integer(s$inputs$Z); W <- s$inputs$weights
  w <- W / max(W)
  V <- large_design[, added, drop = FALSE]
  small_label <- s$layout$ps_names[small]
  large_label <- s$layout$ps_names[large]
  small_index <- s$blocks$ps[[small]]
  contrast_index <- s$blocks$ps[[large]][added]
  contrast_names <- parameters[contrast_index]
  H <- phi[, contrast_index, drop = FALSE]
  Omega <- crossprod(H, H / n)
  if (any(!is.finite(Omega)) || any(diag(Omega) <= 0)) return(.wm_qc_failure(
    "contrast_covariance", "Contrast covariance is nonfinite or its positive diagonal is unrepresentable."))
  contrast_qr <- tryCatch(.wm_qc_qr(H / sqrt(n)), error = identity)
  if (inherits(contrast_qr, "error")) return(.wm_qc_failure(
    "contrast_covariance", conditionMessage(contrast_qr)))
  alpha <- s$parameter[small_index]
  baseline <- as.vector(D %*% alpha) + ps[[small]]$offset
  probability <- stats::plogis(baseline)
  information <- tryCatch({
    multiplier <- sqrt(w) * sqrt(probability) * sqrt(stats::plogis(-baseline)) / sqrt(n)
    if (any(!is.finite(multiplier)) || any(multiplier <= 0)) {
      stop("Retained information row scales are nonfinite or unrepresentable.")
    }
    weighted_D <- D * multiplier
    weighted_V <- V * multiplier
    if (any(!is.finite(c(weighted_D, weighted_V)))) stop("Retained information inputs are nonfinite.")
    factor <- .wm_qc_qr(weighted_D)
    map <- matrix(0, ncol(D), ncol(V))
    map[factor$pivot, ] <- backsolve(factor$R, crossprod(factor$Q, weighted_V))
    if (any(!is.finite(map))) stop("Retained information QR solve is nonfinite.")
    list(map = map, factor = factor)
  }, error = identity)
  if (inherits(information, "error")) return(.wm_qc_failure("information", conditionMessage(information)))
  information_map <- information$map
  raw <- s$raw_scores
  mu <- derivatives <- list()
  # Reconstruct predictions at the retained coefficients to protect row bindings.
  expected_raw <- raw
  for (j in seq_along(ps)) expected_raw[, s$layout$ps_names[j]] <- stats::plogis(
    as.vector(ps[[j]]$design %*% s$parameter[s$blocks$ps[[j]]]) + ps[[j]]$offset)
  for (arm in arms) {
    spec <- pg[[arm]][[1L]]
    mu[[arm]] <- as.vector(spec$design %*% s$parameter[s$blocks$pg[[arm]][[1L]]]) + spec$offset
    expected_raw[, s$layout$pg_names[[arm]][1L]] <- mu[[arm]]
    derivatives[[arm]] <- spec$design
    .wm_reciprocal_agree(object$correction_fit$arms[[arm]]$mean,
      fit$predictions[[paste0("mean", arm)]], paste("Original correction prediction", arm))
  }
  .wm_reciprocal_agree(raw, expected_raw, "Retained primitive predictions")
  expected_scores <- sweep(sweep(raw, 2L, s$parameter[s$blocks$center], "-"),
                           2L, sqrt(s$parameter[s$blocks$variance]), "/")
  for (arm in arms) .wm_reciprocal_agree(
    unname(expected_scores[, s$layout$scores[[arm]], drop = FALSE]),
    fit$graph[[paste0("scores", arm)]], paste("Original score binding", arm))
  gamma <- mean(w * if (fit$estimand == "PATT") Z else 1)
  context <- list(n = n, p = p, arms = arms, parameters = parameters,
    contrast_names = contrast_names, phi = phi, H = H, Omega = Omega,
    contrast_qr = contrast_qr,
    Y = Y, Z = Z, W = W, w = w, gamma = gamma, M = fit$M,
    estimand = fit$estimand, estimate = fit$estimate,
    D = D, V = V, alpha = alpha, offset = ps[[small]]$offset,
    information_map = information_map, raw = raw, layout = s$layout,
    small_label = small_label, large_label = large_label,
    mu = mu, derivatives = derivatives, pg_blocks = s$blocks$pg,
    replication_control = replication)
  evaluator <- local({ retained <- context; function(contrast) .wm_qc_evaluate(retained, contrast) })
  list(status = "quadratic_contrast_inputs_ready", available = TRUE,
    numerically_available = TRUE, replication_available = TRUE, evaluate = evaluator,
    context_binding = context, original_fit_binding = fit,
    replication_control = replication,
    contrast_parameters = contrast_names, full_parameters = parameters,
    contrast_covariance = Omega, full_score_covariance = full_covariance,
    contrast_information_map = information_map, selected_models = control,
    dimension = b, matching_dimensions = s$layout$matching_dimensions,
    numerical_rank_contract = "Full-column QR, reciprocal condition above 64*machine epsilon*max(rows,columns); no statistical Gram threshold.",
    contrast_factorization = contrast_qr,
    information_factorization = information$factor,
    numerical_error_rate_verified = FALSE,
    n = n, M = fit$M, estimand = fit$estimand,
    allocation = list(planned_elements = planned, maximum = maximum,
      contract = "Retained input and serial evaluator element bound; solver/R overhead and copies are not an RSS bound."),
    models_refitted = FALSE, original_fit_preserved = TRUE,
    sampling_inference_available = FALSE, assumptions_verified = FALSE,
    independent_review_completed = FALSE,
    contract = paste("Observable candidate law inputs for exactly nested probability-weighted PS",
      "models and one arm PG. Full score root retained; all ten quadratic terms retained.",
      "Evaluator builds the literal untruncated local-index graphs for covariance only and never replaces the point.",
      "Complete QR diagnostics are not a certified numerical-error bound.",
      "Correct PS/PG, identification, bounded geometry, sampling/covariance limits and numerical",
      "prediction-error conditions remain required. No interval or bootstrap is computed."))
}

.wm_qc_evaluate <- function(context, contrast) {
  h <- .wm_numeric_vector(contrast, length(context$contrast_names), "contrast")
  if (!is.null(names(contrast)) && !identical(names(contrast), context$contrast_names)) {
    stop("Named contrast must retain the declared complete added-parameter order.", call. = FALSE)
  }
  magnitude <- max(abs(h))
  if (magnitude == 0) return(.wm_qc_failure("zero_contrast", "Exactly duplicated PS outputs do not identify the complete finite quadratic basis.", h))
  unit_norm <- sqrt(sum((h / magnitude)^2))
  omega <- (h / magnitude) / unit_norm
  step <- (magnitude / sqrt(context$n)) * unit_norm
  if (!is.finite(step) || step <= 0) return(.wm_qc_failure(
    "contrast_scale", "The finite contrast step is unrepresentable.", h))
  direction <- as.vector(context$V %*% omega) -
    as.vector(context$D %*% (context$information_map %*% omega))
  baseline_logit <- as.vector(context$D %*% context$alpha) + context$offset
  path <- tryCatch(.wm_qc_logistic_contrast(baseline_logit, direction, step), error = identity)
  if (inherits(path, "error")) return(.wm_qc_failure("logistic_path", conditionMessage(path), h))
  raw <- context$raw
  raw[, context$large_label] <- path$large_probability
  raw[, context$small_label] <- path$small_probability
  center <- colMeans(raw)
  variance <- colMeans(sweep(raw, 2L, center, "-")^2)
  if (any(!is.finite(variance)) || any(variance <= 0)) return(.wm_qc_failure("score_scale", "A perturbed pooled score variance is zero or nonfinite.", h))
  standardized <- sweep(sweep(raw, 2L, center, "-"), 2L, sqrt(variance), "/")
  scores <- lapply(context$arms, function(arm) standardized[, context$layout$scores[[arm]], drop = FALSE])
  names(scores) <- context$arms
  # Reuse the original matching implementation solely for its local-index graph.
  graph_fit <- tryCatch(wm_match(context$Y, context$Z, context$W, scores[["0"]],
    if (context$estimand == "PATE") scores[["1"]] else NULL,
    M = context$M, estimand = context$estimand, variance = FALSE), error = identity)
  if (inherits(graph_fit, "error")) return(.wm_qc_failure("graph", conditionMessage(graph_fit), h))
  edges <- graph_fit$graph$edges
  contrast_score <- path$divided_probability
  b <- stats::setNames(numeric(context$p), context$parameters)
  rows <- if (context$estimand == "PATE") {
    context$w * (context$mu[["1"]] - context$mu[["0"]] - context$estimate) / context$gamma
  } else context$Z * context$w * (context$Y - context$mu[["0"]] - context$estimate) / context$gamma
  projection <- list()
  for (arm in context$arms) {
    z <- as.integer(arm); donors <- which(context$Z == z)
    e <- edges[edges$arm == z, , drop = FALSE]
    # This finite raw-probability divided difference spans exactly the same
    # complete quadratic as the standardized three matching coordinates.
    pg_label <- context$layout$pg_names[[arm]][1L]
    basis <- .wm_model_basis(cbind(small_ps = standardized[, context$small_label],
      pg = standardized[, pg_label], contrast = contrast_score))
    # sign/gamma is already in ell. Total n, not arm count, scales the QR.
    linear <- function(f) {
      as.vector(crossprod(context$w[e$query] * e$share /
        (context$n * context$gamma),
        f[e$query, , drop = FALSE] - f[e$donor, , drop = FALSE])) * (2 * z - 1)
    }
    if (any(!is.finite(basis))) return(.wm_qc_failure("projection_basis", paste("Arm", arm, "complete basis is nonfinite."), h))
    ell <- linear(basis)
    solved <- tryCatch({
      multiplier <- sqrt(context$w[donors]) / sqrt(context$n)
      if (any(!is.finite(multiplier)) || any(multiplier <= 0)) stop("Donor row scales are unrepresentable.")
      factor <- .wm_qc_qr(basis[donors, , drop = FALSE] * multiplier)
      if (any(!is.finite(ell))) stop("Signed query/donor functional is nonfinite.")
      y <- forwardsolve(t(factor$R), ell[factor$pivot])
      coefficient <- numeric(ncol(basis))
      coefficient[factor$pivot] <- backsolve(factor$R, y)
      prediction <- as.vector(basis %*% coefficient)
      donor_projection <- as.vector(factor$Q %*% y)
      prediction[donors] <- donor_projection / multiplier
      derivative_projection <- as.vector(crossprod(context$derivatives[[arm]][donors, , drop = FALSE],
        multiplier * donor_projection))
      if (any(!is.finite(c(y, coefficient, prediction, derivative_projection)))) {
        stop("Complete QR projection arithmetic is nonfinite.")
      }
      list(factor = factor, y = y, coefficient = coefficient,
        prediction = prediction, derivative_projection = derivative_projection)
    }, error = identity)
    if (inherits(solved, "error")) return(.wm_qc_failure("projection_qr", paste("Arm", arm, conditionMessage(solved)), h))
    prediction <- solved$prediction
    block <- context$pg_blocks[[arm]][[1L]]
    b[block] <- b[block] + linear(context$derivatives[[arm]]) - solved$derivative_projection
    residual <- context$Y - context$mu[[arm]]
    if (any(!is.finite(c(ell, prediction, b, residual)))) {
      return(.wm_qc_failure("arithmetic", paste("Arm", arm, "projection arithmetic is nonfinite."), h))
    }
    incoming <- graph_fit$loads$incoming[, z + 1L]
    if (context$estimand == "PATE") rows[donors] <- rows[donors] +
      (2 * z - 1) * (context$w[donors] + incoming[donors]) * residual[donors] / context$gamma
    else rows[donors] <- rows[donors] - incoming[donors] * residual[donors] / context$gamma
    rows[donors] <- rows[donors] + context$w[donors] * residual[donors] * prediction[donors]
    projection[[arm]] <- list(basis = basis, weighted_qr = solved$factor,
      signed_functional = ell, orthogonal_functional = as.vector(solved$y),
      coefficient = solved$coefficient, prediction = prediction,
      derivative_projection = solved$derivative_projection, residual = residual,
      basis_contract = "Complete quadratic of standardized small PS, standardized arm PG and exact finite probability divided difference.")
  }
  centered_rows <- rows - mean(rows)
  complete_rows <- centered_rows + as.vector(context$phi %*% b)
  linear_root_variance <- sum((complete_rows / sqrt(context$n))^2)
  cross <- as.vector(crossprod(context$H, complete_rows / context$n))
  if (any(!is.finite(c(rows, complete_rows, cross, linear_root_variance)))) {
    return(.wm_qc_failure("arithmetic", "Complete row/covariance arithmetic is nonfinite.", h))
  }
  conditional <- tryCatch({
    factor <- context$contrast_qr
    scaled_rows <- complete_rows / sqrt(context$n)
    z <- as.vector(crossprod(factor$Q, scaled_rows))
    standardized_h <- forwardsolve(t(factor$R), h[factor$pivot])
    regression <- numeric(length(h))
    regression[factor$pivot] <- backsolve(factor$R, z)
    scaled_residual <- scaled_rows - as.vector(factor$Q %*% z)
    list(regression = regression, residual = scaled_residual * sqrt(context$n),
      mean = sum(z * standardized_h), variance = sum(scaled_residual^2))
  }, error = identity)
  if (inherits(conditional, "error")) return(.wm_qc_failure("conditional_qr", conditionMessage(conditional), h))
  regression <- conditional$regression
  conditional_residual <- conditional$residual
  conditional_mean <- conditional$mean
  conditional_variance <- conditional$variance
  if (any(!is.finite(c(regression, conditional_residual, conditional_mean, conditional_variance)))) {
    return(.wm_qc_failure("arithmetic", "Complete projection/covariance arithmetic is nonfinite.", h))
  }
  list(status = "quadratic_contrast_evaluated", available = TRUE,
    numerically_available = TRUE, contrast = stats::setNames(h, context$contrast_names),
    conditional_mean = conditional_mean, conditional_variance = conditional_variance,
    linear_root_variance = linear_root_variance, contrast_cross_covariance = cross,
    conditional_regression = regression, conditional_residual = conditional_residual,
    base_projection_rows = rows, complete_rows = complete_rows, sensitivity = b,
    full_parameter_names = context$parameters, full_influence = context$phi,
    graph = graph_fit$graph, incoming = graph_fit$loads$incoming,
    pooled_center = center, pooled_variance = variance, projections = projection,
    estimate = context$estimate, n = context$n, M = context$M, estimand = context$estimand,
    finite_probability_contrast = contrast_score,
    contrast_factorization = context$contrast_qr,
    numerical_error_rate_verified = FALSE,
    models_refitted = FALSE, original_point_preserved = TRUE,
    sampling_inference_available = FALSE, assumptions_verified = FALSE,
    contract = "Literal untruncated finite conditional-Gaussian inputs at this contrast; full QR diagnostics are not error certificates. No marginal variance, sampling interval, omitted direction, completed zero or deleted index.")
}
