# Optional direct raw Gaussian average. No covariance support inverse is used.
# Complete centered covariance blocks follow the common fitted-row convention.
# Exact coefficient identities do not bound floating-point roundoff.

.wm_cp_direct_work <- function(M) {
  list(matrix_products_per_pair = 12 + 2 * as.double(M) * (M - 1),
    coefficient_terms_per_pair = (as.double(M) * (M + 1) / 2)^2 - M^2,
    coefficient_matrix_entries = 74 * as.double(M)^2)
}

.wm_cp_direct_psd <- function(x, label, tolerance, return_root = TRUE) {
  if (!is.matrix(x) || !is.numeric(x) || is.complex(x) ||
      nrow(x) != ncol(x) || !nrow(x) || any(!is.finite(x))) {
    stop(label, " must be a finite square numeric matrix.", call. = FALSE)
  }
  scale <- max(1, max(abs(x)))
  asymmetry <- max(abs(x - t(x)))
  if (asymmetry > tolerance * scale) {
    stop(label, " exceeds the declared symmetry precision guard.", call. = FALSE)
  }
  symmetric <- (x + t(x)) / 2
  eig <- eigen(symmetric, symmetric = TRUE, only.values = !return_root)
  if (min(eig$values) < -tolerance * scale) {
    stop(label, " exceeds the declared PSD precision guard.", call. = FALSE)
  }
  corrected <- pmax(eig$values, 0)
  root <- if (return_root) eig$vectors %*%
    diag(sqrt(corrected), nrow = length(corrected)) %*% t(eig$vectors) else NULL
  if (return_root && any(!is.finite(root))) stop(label, " square-root range failure.", call. = FALSE)
  list(root = root, raw_eigenvalues = eig$values,
    raw_asymmetry = asymmetry, negative_roundoff_correction = -pmin(eig$values, 0),
    positive_eigenvalues_discarded = 0L, support_threshold = NULL,
    roundoff_included = FALSE)
}

# Normalized Psi coefficients c00=1: the final prefactor contains det(Q)^(-1/2).
# Matrix-power coefficients retain every ordered word; no commutation shortcut.
.wm_cp_direct_coefficients <- function(P0, P1, M) {
  ell <- coefficients <- matrix(0, M, M)
  coefficients[1L, 1L] <- 1
  if (M == 1L) return(list(coefficients = coefficients, log_coefficients = ell))
  index <- function(a, b) a * M + b + 1L
  previous <- vector("list", M^2)
  previous[[1L]] <- diag(nrow(P0))
  for (degree in seq_len(2L * (M - 1L))) {
    current <- vector("list", M^2)
    for (a in seq.int(max(0L, degree - M + 1L), min(degree, M - 1L))) {
      b <- degree - a
      value <- matrix(0, nrow(P0), ncol(P0))
      if (a > 0L) value <- value + previous[[index(a - 1L, b)]] %*% P0
      if (b > 0L) value <- value + previous[[index(a, b - 1L)]] %*% P1
      current[[index(a, b)]] <- value
      ell[a + 1L, b + 1L] <- sum(diag(value)) / (2 * degree)
    }
    previous <- current
  }
  for (degree in seq_len(2L * (M - 1L))) {
    for (a in seq.int(max(0L, degree - M + 1L), min(degree, M - 1L))) {
      b <- degree - a; value <- 0
      for (i in 0:a) for (j in 0:b) {
        if (i + j) value <- value + (i + j) * ell[i + 1L, j + 1L] *
          coefficients[a - i + 1L, b - j + 1L]
      }
      coefficients[a + 1L, b + 1L] <- value / degree
    }
  }
  if (any(!is.finite(c(ell, coefficients)))) {
    stop("Direct coefficient recursion exceeded floating-point range.", call. = FALSE)
  }
  list(coefficients = coefficients, log_coefficients = ell)
}

.wm_cp_G_average <- function(E, density, Delta0, Delta1, Sigma, M, control) {
  M <- .wm_cp_integer(M, "M")
  if (!is.matrix(E) || !identical(dim(E), c(4L, 2L)) ||
      !is.numeric(E) || any(!is.finite(E)) ||
      !is.matrix(Delta0) || !is.matrix(Delta1) || nrow(Delta0) != 2L ||
      !identical(dim(Delta0), dim(Delta1)) || !ncol(Delta0) ||
      any(!is.finite(c(Delta0, Delta1))) ||
      !identical(dim(Sigma), c(ncol(Delta0), ncol(Delta0)))) {
    stop("Direct geometric dimensions or finite values are invalid.", call. = FALSE)
  }
  .wm_numeric_vector(density, 2L, "direct densities", positive = TRUE)
  work <- .wm_cp_direct_work(M)
  if (work$matrix_products_per_pair > control$maximum_analytic_matrix_products ||
      work$coefficient_terms_per_pair > control$maximum_analytic_coefficient_terms ||
      work$coefficient_matrix_entries > control$maximum_pair_entries) {
    stop("Direct finite-coefficient work budget exceeded.", call. = FALSE)
  }
  Delta <- rbind(Delta0, Delta1)
  Omega <- Delta %*% Sigma %*% t(Delta)
  shifted <- .wm_cp_direct_psd(Omega, "Four-shift Omega", control$analytic_precision_tolerance)
  J0 <- cbind(E[1:2, , drop = FALSE], shifted$root[1:2, , drop = FALSE])
  J1 <- cbind(E[3:4, , drop = FALSE], shifted$root[3:4, , drop = FALSE])
  Q0 <- pi * density[1L] * crossprod(J0)
  Q1 <- pi * density[2L] * crossprod(J1)
  Q <- Q0 + Q1 + diag(c(0, 0, rep(.5, 4L)))
  R <- chol(Q)
  reciprocal_condition <- rcond(Q)
  if (!is.finite(reciprocal_condition) ||
      reciprocal_condition < control$minimum_analytic_reciprocal_condition) {
    stop("Direct SPD solve exceeds the declared conditioning guard.", call. = FALSE)
  }
  inverse_R <- backsolve(R, diag(6L))
  P0 <- crossprod(inverse_R, Q0 %*% inverse_R)
  P1 <- crossprod(inverse_R, Q1 %*% inverse_R)
  value <- .wm_cp_direct_coefficients(P0, P1, M)
  minimum <- min(value$coefficients)
  if (minimum < -control$analytic_precision_tolerance *
      max(1, max(abs(value$coefficients)))) {
    stop("Direct normalized coefficients violate the precision/sign guard.", call. = FALSE)
  }
  coefficient_sum <- sum(value$coefficients)
  if (!is.finite(coefficient_sum) || coefficient_sum <= 0) {
    stop("Direct coefficient sum is unavailable.", call. = FALSE)
  }
  log_value <- log(pi / 4) - sum(log(diag(R))) + log(coefficient_sum)
  average <- exp(log_value)
  if (is.finite(log_value) && average == 0) {
    stop("Direct Gaussian average underflow; no certified absolute-error zero is substituted.",
         call. = FALSE)
  }
  bound <- M / (min(density) * control$projection_floor^2)
  if (!is.finite(average) || average < 0 ||
      average > bound * (1 + control$analytic_precision_tolerance)) {
    stop("Direct Gaussian average exceeded its numerical range/bound guard.", call. = FALSE)
  }
  list(value = average, underflow = is.finite(log_value) && average == 0,
    Omega = Omega, omega_square_root = shifted$root,
    omega_raw_eigenvalues = shifted$raw_eigenvalues,
    omega_negative_roundoff_correction = shifted$negative_roundoff_correction,
    omega_raw_asymmetry = shifted$raw_asymmetry,
    omega_cross_block_max = max(abs(Omega[1:2, 3:4, drop = FALSE])),
    reciprocal_condition = reciprocal_condition,
    normalized_coefficient_sum = coefficient_sum,
    normalized_coefficient_minimum = minimum,
    factor_columns = 4L, SPD_dimension = 6L,
    exact_arithmetic_geometric_error = 0,
    matrix_products_bound = work$matrix_products_per_pair,
    coefficient_terms = work$coefficient_terms_per_pair,
    roundoff_included = FALSE)
}

.wm_cp_direct_one <- function(inference, kernel, Sigma, C, label) {
  control <- kernel$control; n_pairs <- length(kernel$pairs)
  tolerance <- control$analytic_precision_tolerance
  covariance <- .wm_cp_direct_psd(Sigma, "Complete direct covariance", tolerance,
    return_root = FALSE)
  diagnostics <- data.frame(pair = seq_len(n_pairs),
    treated_row = vapply(kernel$pairs, function(x) x$i, 0L),
    control_row = vapply(kernel$pairs, function(x) x$j, 0L),
    status = rep("not_evaluated", n_pairs), average_G = rep(NA_real_, n_pairs),
    omega_minimum_eigenvalue = rep(NA_real_, n_pairs),
    omega_roundoff_correction = rep(NA_real_, n_pairs),
    omega_raw_asymmetry = rep(NA_real_, n_pairs),
    omega_cross_block_max = rep(NA_real_, n_pairs),
    reciprocal_condition = rep(NA_real_, n_pairs),
    coefficient_minimum = rep(NA_real_, n_pairs),
    G_underflow = rep(NA, n_pairs), error = rep("", n_pairs))
  average_R <- error <- 0
  for (block in .wm_cp_chunks(n_pairs, control$chunk_size)) for (k in block) {
    pair <- kernel$pairs[[k]]; first <- kernel$rows[[pair$first]]
    value <- tryCatch(.wm_cp_G_average(first$plane$E,
      c(first$marginal0$density, first$marginal1$density),
      pair$B0, pair$B1, Sigma, kernel$M, control), error = identity)
    if (inherits(value, "error")) {
      diagnostics$status[k] <- "numerical_evaluation_failed"
      diagnostics$error[k] <- conditionMessage(value)
      next
    }
    diagnostics$status[k] <- "complete"
    diagnostics$average_G[k] <- value$value
    diagnostics$omega_minimum_eigenvalue[k] <- min(value$omega_raw_eigenvalues)
    diagnostics$omega_roundoff_correction[k] <- max(value$omega_negative_roundoff_correction)
    diagnostics$omega_raw_asymmetry[k] <- value$omega_raw_asymmetry
    diagnostics$omega_cross_block_max[k] <- value$omega_cross_block_max
    diagnostics$reciprocal_condition[k] <- value$reciprocal_condition
    diagnostics$coefficient_minimum[k] <- value$normalized_coefficient_minimum
    diagnostics$G_underflow[k] <- value$underflow
    average_R <- average_R + pair$K * pair$coefficient * value$value
    error <- error + pair$K * pair$coefficient_error * value$value
  }
  complete <- all(diagnostics$status == "complete")
  average_R <- if (complete) average_R / kernel$normalization else NA_real_
  error <- if (complete) error / kernel$normalization else NA_real_
  gamma <- inference$analysis_target_weight_mean
  diagonal <- inference$diagonal_root_n_variance
  b <- inference$total_sensitivity
  cross <- 2 * sum(b * C); nuisance <- sum(b * as.vector(Sigma %*% b))
  subtraction <- 2 * average_R / gamma^2
  base <- diagonal - subtraction; raw <- base + cross + nuisance
  variance <- raw / inference$n; variance_error <- 2 * error / gamma^2
  finite <- complete && all(is.finite(c(average_R, error, cross, nuisance,
    subtraction, base, raw, variance, variance_error)))
  numerical <- finite && kernel$numerically_available && error <= control$kernel_error_tolerance
  list(status = if (numerical) "raw_direct_variance_computed" else "direct_numerically_unavailable",
    available = numerical, numerically_available = numerical,
    covariance_input = label,
    covariance_convention = "Existing common helper: P_n centered joint influence outer products and centered complete base-row cross products",
    Sigma = Sigma, C = C, total_sensitivity = b,
    diagonal_root_n_variance = diagonal, analysis_target_weight_mean = gamma,
    average_reciprocal = average_R, reciprocal_subtraction = subtraction,
    raw_base_root_n_variance = base, cross_term = cross, nuisance_variance = nuisance,
    raw_root_n_variance = raw, raw_variance = variance,
    negative_raw_root_n_variance = if (finite) raw < 0 else NA,
    positive_raw_root_n_variance = if (finite) raw > 0 else NA,
    reciprocal_quadrature_error_bound = error,
    root_n_variance_quadrature_error_bound = variance_error,
    variance_quadrature_error_bound = variance_error / inference$n,
    complete_raw_covariance_eigenvalues = covariance$raw_eigenvalues,
    complete_raw_covariance_asymmetry = covariance$raw_asymmetry,
    complete_raw_covariance_negative_roundoff = covariance$negative_roundoff_correction,
    complete_raw_covariance_roundoff_correction_applied = FALSE,
    pair_diagnostics = diagnostics, failed_pair_count = sum(diagnostics$status != "complete"),
    pair_count = n_pairs, factor_columns_per_pair = 4L, SPD_dimension = 6L,
    numerical_input_roundoff_corrections_disclosed = TRUE,
    roundoff_included_in_error_bound = FALSE,
    exact_operative_cone_root_n_variance = NA_real_,
    exact_operative_cone_variance = NA_real_,
    Gaussian_interval_authorized = FALSE, original_fit_changed = FALSE,
    population_assumptions_verified = FALSE)
}

.wm_cp_direct_results <- function(inference, kernel, control, blocks) {
  work <- .wm_cp_direct_work(kernel$M)
  products <- length(blocks) * as.double(length(kernel$pairs)) * work$matrix_products_per_pair
  terms <- length(blocks) * as.double(length(kernel$pairs)) * work$coefficient_terms_per_pair
  request <- list(mode = control$direct_variance,
    matrix_products_bound = products, coefficient_terms = terms,
    raw_only_authorizes_replication = FALSE,
    direct_quantity_is_exact_operative_cone_variance = FALSE,
    covariance_convention = "Common centered complete-row empirical covariance, not an uncentered finite-sample second moment")
  if (!is.finite(products) || !is.finite(terms) ||
      products > control$maximum_analytic_matrix_products ||
      terms > control$maximum_analytic_coefficient_terms ||
      work$coefficient_matrix_entries > control$maximum_pair_entries) {
    request$status <- "direct_work_budget_exceeded"
    request$available <- request$numerically_available <- FALSE
    request$unavailable_reason <- "Prospective all-pair/all-covariance direct coefficient budget exceeded; no pair or parameter is filtered."
    return(request)
  }
  request$status <- "direct_variance_results"
  for (name in names(blocks)) {
    value <- tryCatch(.wm_cp_direct_one(inference, kernel,
      blocks[[name]]$Sigma, blocks[[name]]$C, name), error = identity)
    request[[name]] <- if (inherits(value, "error")) {
      list(status = "direct_numerically_unavailable", available = FALSE,
        numerically_available = FALSE, unavailable_reason = conditionMessage(value),
        covariance_input = name, covariance_convention = request$covariance_convention,
        Sigma = blocks[[name]]$Sigma, C = blocks[[name]]$C,
        raw_root_n_variance = NA_real_, raw_variance = NA_real_,
        exact_operative_cone_root_n_variance = NA_real_,
        exact_operative_cone_variance = NA_real_, Gaussian_interval_authorized = FALSE)
    } else value
  }
  request$available <- request$numerically_available <- all(vapply(names(blocks),
    function(x) isTRUE(request[[x]]$numerically_available), FALSE))
  request
}

.wm_cp_complete_raw_direct <- function(inference, control) {
  kernel <- .wm_cp_kernel(inference$fit, inference$parameter_names, control)
  inference$diagonal_root_n_variance <- inference$V0
  inference$raw_C <- inference$C; inference$raw_Sigma <- inference$Sigma
  inference$raw_augmented_row_root_n_variance <- inference$root_n_variance
  inference$analysis_target_weight_mean <- inference$fit$gamma
  inference$reciprocal_units <- kernel$numerator_units
  inference$row_covariance_formula_error <- inference$covariance_formula_error
  inference$covariance_formula_error <- NA_real_
  inference$critical <- list(kernel = kernel)
  inference$critical_control <- control
  inference$direct_variance <- .wm_cp_direct_results(inference, kernel, control,
    list(raw_empirical = list(Sigma = inference$raw_Sigma, C = inference$raw_C)))
  inference$V0 <- inference$reciprocal_subtraction <- NA_real_
  inference$root_n_variance <- inference$variance <- inference$se <- NA_real_
  inference$conf.int <- stats::setNames(rep(NA_real_, 2L), c("lower", "upper"))
  inference$available <- inference$conditional_inference_available <- FALSE
  inference$numerically_available <- inference$direct_variance$numerically_available
  inference$status <- "raw_direct_variance_only"
  inference$unavailable_reason <- if (inference$numerically_available) {
    "Direct scalar variance computed; raw-only mode does not authorize conditional inference."
  } else "Direct scalar calculation is numerically unavailable; raw diagnostics retained."
  inference$interval <- "none: raw direct variance does not imply a Gaussian pivot"
  inference$variance_interpretation <- paste("Rate-free raw direct variance uses complete centered",
    "empirical blocks and the full Gaussian average. Its numerical availability and signed values",
    "are exposed in direct_variance, separately from unavailable conditional inference.",
    "No support threshold, pseudoinverse, bootstrap, Normal interval, operative cone expectation",
    "or actual sampling-variance UI conclusion is supplied.")
  inference
}
