# Constructive correction on the actual fitted scores. This is an opt-in
# correction for wm_model_fit, not another matching rule or automatic theorem.

.wm_lp_controls <- function(correction, dimensions, n, estimand) {
  required <- c("method", "degree", "bandwidth_exponent")
  allowed <- c(required, "guard_exponent")
  if (!is.list(correction) || is.data.frame(correction) ||
      is.null(names(correction)) || anyNA(names(correction)) ||
      any(!nzchar(names(correction))) || anyDuplicated(names(correction)) ||
      !all(required %in% names(correction)) || any(!names(correction) %in% allowed) ||
      !identical(correction$method, "local_polynomial")) {
    stop("correction must declare method='local_polynomial', degree and bandwidth_exponent, ",
         "with only optional guard_exponent; no ridge or arbitrary kernel is supported.", call. = FALSE)
  }
  exponent <- correction[["guard_exponent", exact = TRUE]]
  if (is.null(exponent)) exponent <- 2
  .wm_numeric_vector(exponent, 1L, "guard_exponent", positive = TRUE)
  guard <- n^(-exponent)
  if (!is.finite(guard) || guard <= 0) {
    stop("The declared Gram guard is not positive and representable.", call. = FALSE)
  }
  result <- vector("list", length(dimensions))
  names(result) <- names(dimensions)
  for (arm in names(dimensions)) {
    z <- as.integer(arm)
    d <- dimensions[[arm]]
    q <- .wm_fit_arm_argument(correction$degree, "correction degree", z, estimand)
    q <- .wm_fit_integer(q, "correction degree", 1L)
    alpha <- .wm_fit_arm_argument(correction$bandwidth_exponent,
                                  "correction bandwidth_exponent", z, estimand)
    .wm_numeric_vector(alpha, 1L, "correction bandwidth_exponent", positive = TRUE)
    lower <- if (d <= 2L) 0 else (d - 2) / (2 * d * q)
    upper <- if (d <= 2L) 1 / 4 else 2 / (d * (d + 2))
    if ((d > 2L && q <= (as.double(d)^2 - 4) / 4) ||
        alpha <= lower || alpha >= upper) {
      stop("Local-polynomial degree/bandwidth fail the stated sufficient correction-transfer window for arm ",
           arm, ". No degree or bandwidth is changed automatically.", call. = FALSE)
    }
    bandwidth <- n^(-alpha)
    kernel_constant <- exp(lgamma(d / 2 + 4) - lgamma(4) - (d / 2) * log(pi))
    if (any(!is.finite(c(bandwidth, kernel_constant))) ||
        bandwidth <= 0 || kernel_constant <= 0) {
      stop("Correction bandwidth or fixed kernel amplitude exceeds numerical range.", call. = FALSE)
    }
    result[[arm]] <- list(dimension = d, degree = q, bandwidth_exponent = alpha,
      bandwidth = bandwidth, bandwidth_interval = c(lower = lower, upper = upper),
      gram_guard = guard, guard_exponent = exponent, ridge = 0,
      kernel = "normalized_radial_triweight", kernel_constant = kernel_constant,
      numerical_rate_window_checked = TRUE, population_assumptions_verified = FALSE)
  }
  result
}

.wm_model_score_allocation <- function(n, ps, pg, layout, max_stack_elements) {
  .wm_numeric_vector(max_stack_elements, 1L, "max_stack_elements", positive = TRUE)
  if (n < 2L) stop("At least two rows are required.", call. = FALSE)
  primitive <- sum(vapply(ps, function(m) as.double(ncol(m$design)), numeric(1L))) +
    sum(vapply(pg, function(models) sum(vapply(models,
      function(m) as.double(ncol(m$design)), numeric(1L))), numeric(1L)))
  L <- length(layout$raw_names)
  dims <- as.double(layout$matching_dimensions)
  p <- primitive + 2 * L
  # Includes full equations/IF/zero weight derivative, all raw and standardized
  # tangents, generated coordinates and Jacobian/covariance workspaces.
  planned <- n * (3 * L * p + 4 * L + 4 * p + 2 * sum(dims)) + 4 * p^2
  largest <- c(p, n * p, p^2, n * L * p)
  if (any(!is.finite(c(largest, planned))) || any(largest > .Machine$integer.max)) {
    stop("Complete score-root allocation is not representable.", call. = FALSE)
  }
  if (planned > max_stack_elements) {
    stop("Complete score-root allocation exceeds max_stack_elements.", call. = FALSE)
  }
  list(parameter_dimension = as.integer(p), raw_coordinates = L,
    basis_columns = stats::setNames(integer(length(dims)), names(layout$matching_dimensions)),
    planned_elements = planned, max_stack_elements = max_stack_elements,
    contract = "Score-only root: no quadratic correction coefficients or basis. Conservative element guard, not an RSS bound.")
}

.wm_lp_allocation <- function(n, ps, pg, layout, controls, max_stack_elements) {
  root <- .wm_model_score_allocation(n, ps, pg, layout, max_stack_elements)
  p <- root$parameter_dimension
  dims <- as.double(layout$matching_dimensions)
  log_q <- vapply(controls, function(a)
    lchoose(as.double(a$dimension) + a$degree, a$degree), numeric(1L))
  if (any(!is.finite(log_q)) || any(log_q > log(.Machine$integer.max))) {
    stop("Requested complete local-polynomial basis is not representable.", call. = FALSE)
  }
  Q <- round(exp(log_q))
  # All arm results are retained, but evaluations and arms use one workspace
  # at a time. Conservative use of n also covers either donor-arm count.
  retained <- n * (length(dims) * (8 + 2 * p) + sum(dims) * (2 + p)) + sum(Q * dims)
  workspace <- max(n * (2 * dims + 5 * Q + 8) + 10 * Q^2 + Q * dims + 4 * Q)
  planned <- root$planned_elements + retained + workspace
  largest <- c(n * Q, Q^2, Q * dims, n * dims * p)
  if (any(!is.finite(c(planned, largest))) || any(largest > .Machine$integer.max)) {
    stop("Requested complete local-polynomial allocation is not representable.", call. = FALSE)
  }
  if (planned > max_stack_elements) {
    stop("Complete score-root and local-polynomial allocation exceeds max_stack_elements; ",
         "no requested basis term is removed.", call. = FALSE)
  }
  list(parameter_dimension = p, raw_coordinates = root$raw_coordinates,
    basis_columns = stats::setNames(as.integer(Q), names(controls)),
    score_root_elements = root$planned_elements, correction_retained_elements = retained,
    maximum_local_workspace_elements = workspace, planned_elements = planned,
    max_stack_elements = max_stack_elements,
    contract = "Full score-root/correction element guard with serial local workspaces; inputs, copies and solver/R overhead are not an RSS bound.")
}

.wm_lp_basis <- function(u, exponents) {
  r <- matrix(1, nrow(u), nrow(exponents))
  for (j in seq_len(nrow(exponents))) {
    for (k in seq_len(ncol(u))) {
      if (exponents[j, k] > 0L) r[, j] <- r[, j] * u[, k]^exponents[j, k]
    }
  }
  r
}

.wm_lp_basis_derivative <- function(u, exponents, coordinate, bandwidth) {
  dr <- matrix(0, nrow(u), nrow(exponents))
  for (j in which(exponents[, coordinate] > 0L)) {
    powers <- exponents[j, ]
    coefficient <- -powers[coordinate] / bandwidth
    powers[coordinate] <- powers[coordinate] - 1L
    value <- rep(coefficient, nrow(u))
    for (k in seq_len(ncol(u))) {
      if (powers[k] > 0L) value <- value * u[, k]^powers[k]
    }
    dr[, j] <- value
  }
  dr
}

.wm_lp_reference_derivative <- function(score_derivatives, gradient,
                                        gradient_available, parameter_names) {
  n <- nrow(gradient)
  d <- ncol(gradient)
  p <- length(parameter_names)
  B <- array(0, c(n, d, p), dimnames = list(NULL, NULL, parameter_names))
  D <- matrix(0, n, p, dimnames = list(NULL, parameter_names))
  for (k in seq_len(d)) {
    tangent <- score_derivatives[[k]]
    B[, k, ] <- tangent
    D <- D + tangent * gradient[, k]
  }
  # Finite B and spatial gradients do not guarantee a finite contraction.
  # Guard/boundary NAs are already disclosed; separately flag arithmetic
  # failure on otherwise usable spatial-gradient rows before handoff.
  finite <- rowSums(!is.finite(D)) == 0L
  failure <- gradient_available & !finite
  if (any(failure)) D[failure, ] <- NA_real_
  list(score_tangents = B, reference_mean_derivative = D,
    reference_derivative_failure = failure,
    reference_derivative_available = gradient_available & finite)
}

# Internal evaluation primitive; caller preflights the complete allocation.
# Every evaluation uses this same full arm donor sample, including own rows.
.wm_lp_arm <- function(Y, Z, W, scores, donor_arm, control, exponents,
                       evaluation_scores = scores) {
  n <- length(Y)
  d <- ncol(scores)
  Q <- nrow(exponents)
  donor <- which(Z == donor_arm)
  h <- control$bandwidth
  c_n <- control$gram_guard
  C_d <- control$kernel_constant
  normalizer <- n * h^d
  if (!length(donor) || !is.finite(normalizer) || normalizer <= 0) {
    stop("Local-polynomial donor sample or normalization is invalid.", call. = FALSE)
  }
  if (!is.matrix(evaluation_scores) || ncol(evaluation_scores) != d ||
      any(!is.finite(evaluation_scores))) {
    stop("Local-polynomial evaluation coordinates must be finite with the matching dimension.", call. = FALSE)
  }
  count <- nrow(evaluation_scores)
  value <- numeric(count)
  gradient <- matrix(NA_real_, count, d, dimnames = list(NULL, colnames(scores)))
  minimum <- rep(NA_real_, count)
  failed <- boundary <- logical(count)
  active_count <- integer(count)
  for (i in seq_len(count)) {
    u <- sweep(scores[donor, , drop = FALSE], 2L, evaluation_scores[i, ], "-") / h
    if (any(!is.finite(u))) stop("Local-polynomial coordinate arithmetic overflowed.", call. = FALSE)
    radius_squared <- rowSums(u^2)
    inside <- radius_squared < 1
    active_count[i] <- sum(inside)
    if (!any(inside)) {
      minimum[i] <- 0
      failed[i] <- TRUE
      next
    }
    u <- u[inside, , drop = FALSE]
    radial_gap <- 1 - radius_squared[inside]
    weight <- W[donor[inside]] / normalizer
    outcome <- Y[donor[inside]]
    kernel <- C_d * radial_gap^3
    kw <- weight * kernel
    r <- .wm_lp_basis(u, exponents)
    G <- crossprod(r, r * kw)
    b <- as.vector(crossprod(r, kw * outcome))
    if (any(!is.finite(c(G, b)))) {
      stop("Local-polynomial Gram or moment arithmetic overflowed.", call. = FALSE)
    }
    eigenvalues <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
    if (any(!is.finite(eigenvalues))) stop("Nonfinite local-polynomial Gram eigenvalues.", call. = FALSE)
    minimum[i] <- min(eigenvalues)
    failed[i] <- minimum[i] < c_n
    if (failed[i]) next
    # No ridge, rank trimming or bandwidth enlargement: the stated equations.
    R <- chol(G)
    solve_G <- function(rhs) backsolve(R, forwardsolve(t(R), rhs))
    beta <- as.vector(solve_G(b))
    if (any(!is.finite(beta))) stop("Nonfinite local-polynomial solve.", call. = FALSE)
    value[i] <- beta[1L]
    boundary[i] <- minimum[i] == c_n
    if (boundary[i]) next
    e0 <- numeric(Q)
    e0[1L] <- 1
    inverse_first_column <- as.vector(solve_G(e0))
    for (k in seq_len(d)) {
      dr <- .wm_lp_basis_derivative(u, exponents, k, h)
      dk <- 6 * C_d * u[, k] * radial_gap^2 / h
      dw <- weight * dk
      dG <- crossprod(r, r * dw) + crossprod(dr, r * kw) + crossprod(r, dr * kw)
      db <- as.vector(crossprod(r, dw * outcome) + crossprod(dr, kw * outcome))
      gradient[i, k] <- sum(inverse_first_column * (db - as.vector(dG %*% beta)))
    }
    if (any(!is.finite(c(value[i], gradient[i, ])))) {
      stop("Local-polynomial intercept or actual gradient arithmetic overflowed.", call. = FALSE)
    }
  }
  list(mean = value, gradient = gradient, minimum_gram_eigenvalue = minimum,
    guard_failure = failed, derivative_boundary = boundary,
    gradient_available = !failed & !boundary, local_donor_count = active_count,
    donor_rows = donor, evaluation_count = count, n = n, donor_arm = donor_arm,
    degree = control$degree, bandwidth = h, gram_guard = c_n,
    kernel = control$kernel, kernel_constant = C_d, ridge = 0,
    exponents = exponents, supplied_weight_units = "original raw W",
    normalization = "Original full n times h^d; no arm-size or weight normalization",
    derivative = "Actual spatial derivative of the fitted intercept, including kernel and basis derivatives",
    assumptions_verified = FALSE)
}

.wm_model_local_polynomial_fit <- function(Y, Z, weights, ps, pg, layout,
    correction, M, estimand, inference, controls, conf.level, max_stack_elements,
    ties, call) {
  if (inference != "none") {
    stop("local_polynomial currently returns point estimation and complete reference-inference inputs; ",
         "use inference='none' and supply the remaining premises to wm_fitted_inference explicitly.",
         call. = FALSE)
  }
  n <- length(Y)
  lp <- .wm_lp_controls(correction, layout$matching_dimensions, n, estimand)
  allocation <- .wm_lp_allocation(n, ps, pg, layout, lp, max_stack_elements)
  arguments <- list(Y = Y, Z = Z, weights = weights, ps_models = ps,
    pg0_models = pg[["0"]], pg1_models = if (estimand == "PATE") pg[["1"]] else NULL,
    estimand = estimand, max_stack_elements = max_stack_elements, score_only = TRUE,
    influence = TRUE)
  solved <- .wm_wdsm_fit_stack(arguments, .wm_model_stack, controls)
  stack <- solved$stack
  corrected <- reference <- tangents <- list()
  for (arm in names(lp)) {
    scores <- if (arm == "0") stack$scores0 else stack$scores1
    exponents <- .wm_fit_polynomial_exponents(ncol(scores), lp[[arm]]$degree,
                                              allocation$basis_columns[[arm]])
    corrected[[arm]] <- .wm_lp_arm(stack$inputs$Y, stack$inputs$Z,
      stack$inputs$weights, scores, as.integer(arm), lp[[arm]], exponents)
    derivative <- .wm_lp_reference_derivative(stack$score_derivatives[[arm]],
      corrected[[arm]]$gradient, corrected[[arm]]$gradient_available,
      stack$parameter_names)
    tangents[[arm]] <- derivative$score_tangents
    reference[[arm]] <- derivative$reference_mean_derivative
    corrected[[arm]]$reference_mean_derivative <- derivative$reference_mean_derivative
    corrected[[arm]]$reference_derivative_failure <- derivative$reference_derivative_failure
    corrected[[arm]]$reference_derivative_available <- derivative$reference_derivative_available
  }
  fit <- wm_match(Y = stack$inputs$Y, Z = stack$inputs$Z,
    weights = stack$inputs$weights, scores0 = stack$scores0, scores1 = stack$scores1,
    M = M, estimand = estimand, method = "self_normalized",
    mean0 = corrected[["0"]]$mean,
    mean1 = if (estimand == "PATE") corrected[["1"]]$mean else NULL,
    variance = FALSE, tie_rule = ties$rule, tie_seed = ties$seed,
    tie_tolerance = ties$tolerance)
  failures <- vapply(corrected, function(a) sum(a$guard_failure), integer(1L))
  boundaries <- vapply(corrected, function(a) sum(a$derivative_boundary), integer(1L))
  reference_failures <- vapply(corrected,
    function(a) sum(a$reference_derivative_failure), integer(1L))
  numerical <- sum(failures) == 0L && sum(boundaries) == 0L &&
    sum(reference_failures) == 0L
  if (!numerical) {
    warning("Local-polynomial correction used the declared zero completion at ", sum(failures),
      " evaluation(s), with ", sum(boundaries), " derivative boundary point(s) and ",
      sum(reference_failures), " nonfinite reference-derivative row(s); ",
      "reference-inference inputs are unavailable. Inspect correction_fit diagnostics.", call. = FALSE)
  }
  inference_ties <- identical(ties$rule, "row_order")
  input_arguments <- if (numerical && inference_ties) list(fit = fit,
    nuisance_influence = stack$nuisance_influence,
    mean_derivative0 = reference[["0"]],
    mean_derivative1 = if (estimand == "PATE") reference[["1"]] else NULL,
    weight_derivative = stack$weight_derivative) else NULL
  inference_result <- list(method = "none", status = "point_estimate_only",
    available = FALSE, sampling_inference_available = FALSE,
    unavailable_reason = "Sampling inference was not requested; complete population and covariance inputs remain required.",
    root_n_variance = NA_real_, variance = NA_real_, se = NA_real_,
    conf.int = c(lower = NA_real_, upper = NA_real_), conf.level = conf.level,
    assumptions_verified = FALSE, original_refit_limit_agreement_declared = FALSE,
    contract = "Point estimate only; no sampling variance, interval or bootstrap is asserted.")
  recipe <- list(J = length(ps), K = lengths(pg),
    matching_dimensions = layout$matching_dimensions, model_names = layout$model_names,
    shared_ps_fitted_once = TRUE, supplied_weights_known = TRUE,
    correction = "Constructive supplied-weight local polynomial",
    legacy_double_score_reduction = FALSE, quadratic_correction_fitted = FALSE)
  structure(list(call = call, estimate = fit$estimate, raw_estimate = fit$raw_estimate,
    correction = fit$correction, n = fit$n, M = fit$M, estimand = estimand,
    ps_weighting = stack$ps_weighting, se = inference_result$se,
    variance = inference_result$variance, root_n_variance = inference_result$root_n_variance,
    conf.int = inference_result$conf.int, inference = inference_result,
    fit = fit, nuisance = stack, assembly = NULL, model_recipe = recipe,
    correction_fit = list(method = "local_polynomial", controls = lp, arms = corrected,
      guard_failures = failures, derivative_boundaries = boundaries,
      reference_derivative_failures = reference_failures,
      status = if (numerical) "computed" else "completed_with_guard_or_derivative_failures",
      reference_derivatives_available = numerical, assumptions_verified = FALSE),
    inference_inputs = list(arguments = input_arguments, score_tangents = tangents,
      numerically_available = numerical, supported_tie_rule = inference_ties,
      assumptions_verified = FALSE,
      interpretation = paste("Complete score-root influence and plug-in reference derivative B'gradient.",
        "This is not the full derivative through refitting the local regression.",
        "Caller supplies covariance scope, transport and any critical controls under their own premises.",
        "The constructive correction requires bounded outcomes and bounded positive W;",
        "it does not inherit unbounded Gaussian/full-X alternatives.",
        "source_random remains point-only even when arrays are finite.")),
    allocation = allocation,
    solver = list(controls = controls, attempts = solved$attempts, refined = solved$refined,
      normalized_ps_threshold = 1e-10, newton_threshold = min(1e-9, 1 / fit$n)),
    assumptions_verified = FALSE, bootstrap_requested = FALSE,
    scope = paste("Original known-W self-normalized matching with one common armwise local-polynomial",
      "correction on actual fitted scores. Complete fixed-p score root excludes auxiliary local coefficients.",
      "The sufficient degree/bandwidth window is checked, but identification, weighted centering,",
      "bounded outcomes and bounded positive W, smoothness, score displacement, support",
      "and all sampling/covariance conditions",
      "remain unverified. No generic MR, scalar fitted-law, bootstrap or coverage claim.")),
    class = c("wm_model_fit", "list"))
}
