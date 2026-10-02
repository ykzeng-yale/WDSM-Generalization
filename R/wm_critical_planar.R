# Whole-function critical-planar covariance estimation. Original fit is retained.
# Valid population guard/rate premises are caller declarations, never certified
# by these finite arrays. Analytic A bounds exclude floating-point roundoff.

.wm_cp_integer <- function(x, label, minimum = 1) {
  .wm_numeric_vector(x, 1L, label, positive = TRUE)
  if (x < minimum || x != floor(x) || x > .Machine$integer.max) {
    stop(label, " must be a representable integer >= ", minimum, ".", call. = FALSE)
  }
  as.integer(x)
}

.wm_cp_control <- function(control, fit, parameters) {
  control <- .wm_reciprocal_list(control, "critical_control")
  required <- c("tangents0", "tangents1", "auxiliary_rows", "bandwidth",
    "density_bounds", "marginal_mass_floor", "projection_floor",
    "tangent_cap", "residual_cap", "support")
  defaults <- list(quadrature_tolerance = .01 / fit$n,
    kernel_error_tolerance = 1 / sqrt(fit$n), max_nodes = 131073L,
    chunk_size = 2048L, maximum_distance_evaluations = 5000000L,
    maximum_pairs = 25000L, maximum_pair_entries = 2000000L,
    maximum_laplace_terms = 200000000L, maximum_hermite_order = 31L,
    maximum_geometric_evaluations = 20000000L,
    maximum_draw_entries = 2000000L)
  if (!all(required %in% names(control)) ||
      any(!names(control) %in% c(required, names(defaults)))) {
    stop("critical_control has missing or unsupported fields.", call. = FALSE)
  }
  for (name in names(defaults)) if (!name %in% names(control)) {
    control[[name]] <- defaults[[name]]
  }
  n <- fit$n; p <- length(parameters)
  for (name in c("tangents0", "tangents1")) {
    x <- control[[name]]
    if (!is.numeric(x) || is.complex(x) ||
        !identical(dim(x), c(as.integer(n), 2L, as.integer(p))) ||
        any(!is.finite(x)) || !identical(dimnames(x)[[3L]], parameters)) {
      stop(name, " must be finite n-by-2-by-p in the complete parameter order.",
           call. = FALSE)
    }
    .wm_fi_row_names(dimnames(x)[[1L]], n, name)
    if (!is.null(dimnames(x)[[2L]]) &&
        !identical(dimnames(x)[[2L]], c("1", "2"))) {
      stop(name, " coordinate labels must be NULL or positions '1','2'.", call. = FALSE)
    }
  }
  a <- control$auxiliary_rows
  if (!is.numeric(a) || is.complex(a) || !is.null(dim(a)) ||
      length(a) < 2L || any(!is.finite(a)) || any(a != floor(a)) ||
      any(a < 1 | a > n) || anyDuplicated(a) || is.unsorted(a, strictly = TRUE)) {
    stop("auxiliary_rows must be increasing distinct original row indices.", call. = FALSE)
  }
  control$auxiliary_rows <- as.integer(a)
  evaluation <- setdiff(seq_len(n), control$auxiliary_rows)
  if (length(evaluation) < 2L) stop("The evaluation fold requires >=2 rows.", call. = FALSE)
  for (name in c("bandwidth", "marginal_mass_floor", "projection_floor",
                 "tangent_cap", "residual_cap", "quadrature_tolerance",
                 "kernel_error_tolerance")) {
    .wm_numeric_vector(control[[name]], 1L, paste("critical_control", name), positive = TRUE)
  }
  if (control$projection_floor >= 1) {
    stop("projection_floor must be strictly below one.", call. = FALSE)
  }
  .wm_numeric_vector(control$density_bounds, 2L, "density_bounds", positive = TRUE)
  if (control$density_bounds[1L] >= control$density_bounds[2L]) {
    stop("density_bounds must be strictly increasing.", call. = FALSE)
  }
  if (!identical(control$support, "smooth_stack_threshold")) {
    stop("Declare support = 'smooth_stack_threshold' with its complete covariance-rate premises.",
         call. = FALSE)
  }
  for (name in c("max_nodes", "chunk_size", "maximum_distance_evaluations",
    "maximum_pairs", "maximum_pair_entries", "maximum_laplace_terms",
    "maximum_hermite_order", "maximum_geometric_evaluations",
    "maximum_draw_entries")) {
    control[[name]] <- .wm_cp_integer(control[[name]], name,
      if (name == "max_nodes") 3L else 1L)
  }
  if (2 * as.double(fit$M) - 1 > control$maximum_hermite_order) {
    stop("The unchanged M exceeds the declared Hermite-order budget.", call. = FALSE)
  }
  if ((2 * as.double(fit$M) - 1)^2 >
      min(control$maximum_pair_entries, control$maximum_geometric_evaluations)) {
    stop("The Hermite matrix/tensor exceeds prospective numerical budgets.", call. = FALSE)
  }
  if (as.double(p)^2 > control$maximum_pair_entries) {
    stop("The full nuisance covariance exceeds the declared matrix budget.", call. = FALSE)
  }
  control
}

.wm_cp_chunks <- function(n, size) {
  if (!n) return(list())
  lapply(seq.int(1L, n, by = size), function(first) {
    seq.int(first, min(n, first + as.double(size) - 1))
  })
}

.wm_cp_distance <- function(x, center) {
  sqrt(rowSums(sweep(x, 2L, center, "-")^2))
}

.wm_cp_logsum <- function(x) {
  maximum <- max(x)
  if (identical(maximum, -Inf)) return(-Inf)
  maximum + log(sum(exp(x - maximum)))
}

# Physicists' exp(-x^2) Hermite rule; weights normalized to sum one.
.wm_cp_hermite <- function(M) {
  q <- 2L * M - 1L
  J <- matrix(0, q, q)
  if (q > 1L) {
    k <- seq_len(q - 1L)
    J[cbind(k, k + 1L)] <- sqrt(k / 2)
    J[cbind(k + 1L, k)] <- sqrt(k / 2)
  }
  eig <- eigen(J, symmetric = TRUE)
  order <- order(eig$values)
  list(nodes = eig$values[order], weights = eig$vectors[1L, order]^2,
       order = q, exact_polynomial_degree = 2L * q - 1L,
       convention = "physicists exp(-x^2), normalized weights")
}

# The finite formula is exact in exact arithmetic; double roundoff is not bounded.
.wm_cp_G <- function(E, density, displacement0, displacement1, M, rule) {
  A0 <- E[1:2, , drop = FALSE]; A1 <- E[3:4, , drop = FALSE]
  lambda <- pi * density
  H <- lambda[1L] * crossprod(A0) + lambda[2L] * crossprod(A1)
  b <- lambda[1L] * as.vector(crossprod(A0, displacement0)) +
    lambda[2L] * as.vector(crossprod(A1, displacement1))
  R <- chol(H)
  mu <- -as.vector(backsolve(R, forwardsolve(t(R), b)))
  kappa <- lambda[1L] * sum((A0 %*% mu + displacement0)^2) +
    lambda[2L] * sum((A1 %*% mu + displacement1)^2)
  if (any(!is.finite(c(R, mu, kappa)))) {
    stop("Critical geometric solve exceeded floating-point range.", call. = FALSE)
  }
  logarithms <- numeric(rule$order^2)
  index <- 0L
  polynomial <- function(t) {
    if (M == 1L || t == 0) return(0)
    .wm_cp_logsum((0:(M - 1L)) * log(t) - lgamma(seq_len(M)))
  }
  for (a in seq_len(rule$order)) for (b in seq_len(rule$order)) {
    index <- index + 1L
    v <- mu + as.vector(backsolve(R, c(rule$nodes[a], rule$nodes[b])))
    t0 <- lambda[1L] * sum((A0 %*% v + displacement0)^2)
    t1 <- lambda[2L] * sum((A1 %*% v + displacement1)^2)
    if (any(!is.finite(c(t0, t1)))) stop("Critical polynomial range failure.", call. = FALSE)
    logarithms[index] <- log(rule$weights[a]) + log(rule$weights[b]) +
      polynomial(t0) + polynomial(t1)
  }
  log_value <- log(pi) - sum(log(diag(R))) - kappa + .wm_cp_logsum(logarithms)
  value <- exp(log_value)
  if (!is.finite(value)) stop("Critical geometric sum exceeded range.", call. = FALSE)
  list(value = value, underflow = is.finite(log_value) && value == 0,
       quadrature_error_exact_arithmetic = 0, roundoff_included = FALSE)
}

.wm_cp_marginal <- function(scores, center, donors, W, m, control) {
  h <- control$bandwidth; probability <- numeric(length(donors))
  for (block in .wm_cp_chunks(length(donors), control$chunk_size)) {
    probability[block] <- pmax(0, 1 -
      .wm_cp_distance(scores[donors[block], , drop = FALSE], center) / h) / (pi / 3)
  }
  mass <- sum(probability)
  raw_density <- mass / (m * h^2)
  rejected <- mass < control$marginal_mass_floor * m * h^2
  positive <- probability > 0
  if (rejected || !any(positive)) {
    law <- list(weights = W[1L], probability = 1)
  } else {
    law <- list(weights = W[donors[positive]], probability = probability[positive] / mass)
  }
  list(raw_density = raw_density,
    density = min(control$density_bounds[2L], max(control$density_bounds[1L], raw_density)),
    mass = mass, rejected = rejected,
    density_clipped = raw_density < control$density_bounds[1L] ||
      raw_density > control$density_bounds[2L], law = law)
}

.wm_cp_plane <- function(T, center, auxiliary, control) {
  l <- numeric(length(auxiliary)); first <- numeric(4); second <- matrix(0, 4, 4)
  for (block in .wm_cp_chunks(length(auxiliary), control$chunk_size)) {
    v <- sweep(T[auxiliary[block], , drop = FALSE], 2L, center, "-") / control$bandwidth
    weight <- pmax(0, 1 - sqrt(rowSums(v^2)))
    l[block] <- weight
    positive <- weight > 0
    if (any(positive)) {
      v <- v[positive, , drop = FALSE]; weight <- weight[positive]
      first <- first + colSums(v * weight)
      second <- second + crossprod(v, v * weight)
    }
  }
  mass <- sum(l)
  Q <- if (mass > 0) second / mass - tcrossprod(first / mass) else matrix(0, 4, 4)
  eig <- eigen((Q + t(Q)) / 2, symmetric = TRUE)
  E <- eig$vectors[, 1:2, drop = FALSE]
  projections <- c(min(svd(E[1:2, , drop = FALSE], nu = 0, nv = 0)$d),
    min(svd(E[3:4, , drop = FALSE], nu = 0, nv = 0)$d))
  retained <- mass > 0 && eig$values[2L] >= (3 / 20) / 2 &&
    all(projections >= control$projection_floor)
  list(E = E, covariance = Q, eigenvalues = eig$values, projections = projections,
    mass = mass, retained = retained)
}

.wm_cp_kernel <- function(fit, parameters, control) {
  n <- fit$n; p <- length(parameters); a <- control$auxiliary_rows
  evaluation <- setdiff(seq_len(n), a); m <- length(a); N <- length(evaluation)
  Z <- fit$data$Z; W <- fit$analysis_weights
  S0 <- fit$graph$scores0; S1 <- fit$graph$scores1; T <- cbind(S0, S1)
  raw_epsilon <- fit$data$Y - ifelse(Z == 1L, fit$predictions$mean1, fit$predictions$mean0)
  epsilon <- pmin(control$residual_cap, pmax(-control$residual_cap, raw_epsilon))
  raw_B0 <- control$tangents0; raw_B1 <- control$tangents1
  B0 <- raw_B0; B1 <- raw_B1
  B0[] <- pmin(control$tangent_cap, pmax(-control$tangent_cap, raw_B0))
  B1[] <- pmin(control$tangent_cap, pmax(-control$tangent_cap, raw_B1))
  treated <- evaluation[Z[evaluation] == 1L]; controls <- evaluation[Z[evaluation] == 0L]
  prospective_distances <- length(treated) * as.double(length(controls) + 3 * m)
  if (prospective_distances > control$maximum_distance_evaluations) {
    stop("Critical distance-work budget exceeded before construction.", call. = FALSE)
  }
  q <- 2L * fit$M - 1L
  rule <- .wm_cp_hermite(fit$M)
  rows <- vector("list", length(treated)); pairs <- list(); laplace <- 0
  pair_count <- 0L; raw_pair_count <- 0L
  A <- function(law, root_weight) {
    possible <- if (fit$M == 1L) control$max_nodes else
      min(control$max_nodes, floor((control$maximum_laplace_terms - laplace) /
        length(law$weights)))
    if (possible < 3 && fit$M != 1L) stop("Critical Laplace-work budget exhausted.", call. = FALSE)
    value <- .wm_gt_tagged_load(law$weights, law$probability,
      matrix(0, length(law$weights), 1L), root_weight, fit$M,
      control$quadrature_tolerance, max(3, possible))
    laplace <<- laplace + value$nodes * length(law$weights)
    list(value = value$value, error = value$error_bound[1L], nodes = value$nodes)
  }
  for (k in seq_along(treated)) {
    i <- treated[k]
    d0 <- .wm_cp_marginal(S0, S0[i, ], a[Z[a] == 0L], W, m, control)
    d1 <- .wm_cp_marginal(S1, S1[i, ], a[Z[a] == 1L], W, m, control)
    plane <- .wm_cp_plane(T, T[i, ], a, control)
    rows[[k]] <- list(row = i, marginal0 = d0, marginal1 = d1, plane = plane)
    A1 <- NULL
    for (block in .wm_cp_chunks(length(controls), control$chunk_size)) {
      K <- pmax(0, 1 - .wm_cp_distance(T[controls[block], , drop = FALSE], T[i, ]) /
        control$bandwidth)
      take <- which(K > 0)
      raw_pair_count <- raw_pair_count + length(take)
      if (!plane$retained || !length(take)) next
      if (pair_count + length(take) > control$maximum_pairs ||
          (pair_count + length(take)) * as.double(4 * p + 7) > control$maximum_pair_entries ||
          (pair_count + length(take)) * as.double(q)^2 > control$maximum_geometric_evaluations) {
        stop("Critical retained-pair/geometry budget exceeded; no pairs are trimmed.", call. = FALSE)
      }
      if (is.null(A1)) A1 <- A(d1$law, W[i])
      for (v in take) {
        j <- controls[block[v]]; A0 <- A(d0$law, W[j])
        factor <- W[i] * epsilon[i] * W[j] * epsilon[j]
        pair_count <- pair_count + 1L
        pairs[[pair_count]] <- list(i = i, j = j, first = k, K = K[v],
          coefficient = factor * A1$value * A0$value,
          coefficient_error = abs(factor) * (A1$error * abs(A0$value) +
            A0$error * abs(A1$value) + A1$error * A0$error),
          A1 = A1, A0 = A0,
          B0 = matrix(B0[j, , ], 2L, p) - matrix(B0[i, , ], 2L, p),
          B1 = matrix(B1[j, , ], 2L, p) - matrix(B1[i, , ], 2L, p))
      }
    }
  }
  normalization <- N * as.double(N - 1L) * control$bandwidth^2 * (pi / 3)
  G_bound <- fit$M / (control$density_bounds[1L] * control$projection_floor^2)
  global_error <- sum(vapply(pairs, function(x) x$K * x$coefficient_error, 0)) *
    G_bound / normalization
  global_envelope <- sum(vapply(pairs, function(x) x$K * abs(x$coefficient), 0)) *
    G_bound / normalization
  if (any(!is.finite(c(normalization, G_bound, global_error, global_envelope))) || normalization <= 0) {
    stop("Critical kernel bounds exceeded floating-point range.", call. = FALSE)
  }
  list(n = n, parameters = parameters, M = fit$M, control = control,
    scores0 = S0, scores1 = S1, Y = fit$data$Y, Z = Z, W = fit$weights,
    analysis_weights = W,
    numerator_units = "Original fit analysis weights; constant raw-W rescaling cancels in R/gamma^2",
    mean0 = fit$predictions$mean0, mean1 = fit$predictions$mean1,
    rows = rows, pairs = pairs, rule = rule, normalization = normalization,
    cache_binding = list(rows = rows, pairs = pairs, rule = rule, normalization = normalization),
    global_absolute_R_bound = global_envelope,
    global_quadrature_error_bound = global_error, roundoff_included = FALSE,
    numerically_available = global_error <= control$kernel_error_tolerance,
    diagnostics = list(auxiliary_rows = a, evaluation_rows = evaluation,
      auxiliary_size = m, evaluation_size = N, raw_local_pairs = raw_pair_count,
      retained_pairs = pair_count, distance_evaluations_bound = prospective_distances,
      laplace_terms = laplace, hermite_order = q, hermite_tensor_nodes = q^2,
      tangent0_cap_entries = sum(raw_B0 != B0), tangent1_cap_entries = sum(raw_B1 != B1),
      residual_cap_entries = sum(raw_epsilon != epsilon),
      rejected_planes = sum(!vapply(rows, function(x) x$plane$retained, FALSE)),
      marginal0_guard_count = sum(vapply(rows, function(x) x$marginal0$rejected, FALSE)),
      marginal1_guard_count = sum(vapply(rows, function(x) x$marginal1$rejected, FALSE)),
      density0_clip_count = sum(vapply(rows, function(x) x$marginal0$density_clipped, FALSE)),
      density1_clip_count = sum(vapply(rows, function(x) x$marginal1$density_clipped, FALSE))))
}

.wm_cp_evaluate <- function(kernel, u) {
  p <- length(kernel$parameters)
  .wm_fv_number(u, p, "critical index")
  value <- 0; error <- 0; underflow <- 0L
  for (block in .wm_cp_chunks(length(kernel$pairs), kernel$control$chunk_size)) {
    for (k in block) {
      pair <- kernel$pairs[[k]]; first <- kernel$rows[[pair$first]]
      G <- .wm_cp_G(first$plane$E,
        c(first$marginal0$density, first$marginal1$density),
        as.vector(pair$B0 %*% u), as.vector(pair$B1 %*% u), kernel$M, kernel$rule)
      value <- value + pair$K * pair$coefficient * G$value
      error <- error + pair$K * pair$coefficient_error * G$value
      underflow <- underflow + G$underflow
    }
  }
  value <- value / kernel$normalization; error <- error / kernel$normalization
  if (any(!is.finite(c(value, error)))) stop("Critical reciprocal sum range failure.", call. = FALSE)
  list(R = value, quadrature_error_bound = error,
    global_quadrature_error_bound = kernel$global_quadrature_error_bound,
    G_underflow_count = underflow, roundoff_included = FALSE)
}

.wm_cp_covariance <- function(inference) {
  raw <- inference$Sigma
  eig <- eigen((raw + t(raw)) / 2, symmetric = TRUE)
  threshold <- inference$n^(-1 / 3)
  keep <- eig$values > threshold
  E <- eig$vectors[, keep, drop = FALSE]
  values <- eig$values[keep]
  p <- nrow(raw)
  covariance <- inverse <- projector <- matrix(0, p, p)
  square_root <- matrix(0, p, p)
  if (length(values)) {
    projector <- tcrossprod(E)
    covariance <- tcrossprod(sweep(E, 2L, sqrt(values), "*"))
    inverse <- tcrossprod(sweep(E, 2L, 1 / sqrt(values), "*"))
    square_root <- E %*% diag(sqrt(values), nrow = length(values)) %*% t(E)
  }
  C <- as.vector(projector %*% inference$C)
  regression <- as.vector(inverse %*% C)
  q <- if (length(values)) sum(as.vector(crossprod(E, C))^2 / values) else 0
  if (any(!is.finite(c(covariance, inverse, square_root, C, regression, q))) || q < 0) {
    stop("Critical covariance-support arithmetic failed.", call. = FALSE)
  }
  names(C) <- names(regression) <- inference$parameter_names
  dimnames(covariance) <- dimnames(inverse) <- dimnames(square_root) <-
    dimnames(projector) <- dimnames(raw)
  list(raw_Sigma = raw, raw_C = inference$C, raw_eigenvalues = eig$values,
    threshold = threshold, selected = keep, rank = sum(keep), projector = projector,
    Sigma = covariance, Sigma_plus = inverse, square_root = square_root,
    C = C, regression = regression, q = q,
    support_mode = "smooth_stack_threshold", assumptions_verified = FALSE)
}

.wm_cp_schur <- function(object, u) {
  value <- .wm_cp_evaluate(object$critical$kernel, u)
  gamma <- object$analysis_target_weight_mean
  raw_V0 <- object$diagonal_root_n_variance - 2 * value$R / gamma^2
  raw_D <- raw_V0 - object$critical$support$q
  D <- max(raw_D, 0)
  error <- 2 * value$quadrature_error_bound / gamma^2
  if (any(!is.finite(c(raw_V0, raw_D, D, error)))) {
    stop("Critical Schur/cone arithmetic exceeded numerical range.", call. = FALSE)
  }
  list(R = value$R, raw_V0 = raw_V0, raw_D = raw_D,
    negative_raw_D = raw_D < 0, D = D, V0 = object$critical$support$q + D,
    cone_adjustment = D - raw_D,
    quadrature_error_bound = error,
    G_underflow_count = value$G_underflow_count, roundoff_included = FALSE)
}

.wm_cp_complete <- function(inference, control) {
  kernel <- .wm_cp_kernel(inference$fit, inference$parameter_names, control)
  support <- .wm_cp_covariance(inference)
  inference$diagonal_root_n_variance <- inference$V0
  inference$raw_C <- inference$C; inference$raw_Sigma <- inference$Sigma
  inference$raw_augmented_row_root_n_variance <- inference$root_n_variance
  inference$analysis_target_weight_mean <- inference$fit$gamma
  inference$reciprocal_units <- kernel$numerator_units
  inference$row_covariance_formula_error <- inference$covariance_formula_error
  inference$covariance_formula_error <- NA_real_
  inference$critical <- list(kernel = kernel, support = support)
  inference$critical_control <- control
  inference$critical$origin <- .wm_cp_schur(inference,
    stats::setNames(numeric(length(inference$parameter_names)), inference$parameter_names))
  inference$C <- support$C; inference$Sigma <- support$Sigma
  inference$cross_term <- 2 * sum(inference$total_sensitivity * support$C)
  inference$nuisance_variance <- sum(as.vector(crossprod(support$square_root,
    inference$total_sensitivity))^2)
  inference$V0 <- inference$critical$origin$V0
  inference$reciprocal_subtraction <- 2 * inference$critical$origin$R /
    inference$analysis_target_weight_mean^2
  inference$root_n_variance <- inference$variance <- inference$se <- NA_real_
  inference$conf.int <- stats::setNames(rep(NA_real_, 2L), c("lower", "upper"))
  inference$available <- inference$numerically_available <-
    inference$conditional_inference_available <- kernel$numerically_available
  inference$status <- if (kernel$numerically_available) "conditional_critical_planar" else
    "critical_numerical_error_unavailable"
  inference$unavailable_reason <- if (kernel$numerically_available) "" else
    "Propagated analytic A error exceeds kernel_error_tolerance."
  inference$interval <- "none: non-Gaussian law requires mixture quantiles"
  inference$variance_interpretation <- paste("Analytic conditional-law variance is an unevaluated",
    "Gaussian expectation; finite-B Monte Carlo summaries are supplied by wm_bootstrap.")
  inference
}

.wm_cp_check <- function(object) {
  if (!inherits(object, "wm_fitted_inference") ||
      !identical(object$status, "conditional_critical_planar") ||
      !identical(object$covariance_scope, "critical_planar") ||
      !isTRUE(object$available) || !isTRUE(object$numerically_available)) {
    stop("Require successful conditional critical_planar inference.", call. = FALSE)
  }
  state <- .wm_fv_state(object$fit)
  if (!state$pATE || !identical(state$dimensions, c(2L, 2L)) ||
      !identical(object$estimate, object$fit$estimate) || !identical(object$n, state$n)) {
    stop("Critical source point/graph dimensions changed.", call. = FALSE)
  }
  control <- .wm_cp_control(object$critical_control, object$fit, object$parameter_names)
  kernel <- object$critical$kernel
  if (!identical(control, kernel$control) || !identical(object$M, kernel$M) ||
      !identical(object$parameter_names, kernel$parameters)) {
    stop("Critical kernel/control binding changed.", call. = FALSE)
  }
  for (name in c("rows", "pairs", "rule", "normalization")) {
    if (!identical(kernel[[name]], kernel$cache_binding[[name]])) {
      stop("Critical constructed cache binding changed.", call. = FALSE)
    }
  }
  for (name in c("scores0", "scores1")) {
    if (!identical(kernel[[name]], object$fit$graph[[name]])) {
      stop("Critical matching-score binding changed.", call. = FALSE)
    }
  }
  for (name in c("Y", "Z")) if (!identical(kernel[[name]], object$fit$data[[name]])) {
    stop("Critical original-row data binding changed.", call. = FALSE)
  }
  if (!identical(kernel$W, object$fit$weights) ||
      !identical(kernel$analysis_weights, object$fit$analysis_weights) ||
      !identical(kernel$mean0, object$fit$predictions$mean0) ||
      !identical(kernel$mean1, object$fit$predictions$mean1)) {
    stop("Critical known-weight/prediction binding changed.", call. = FALSE)
  }
  expected <- .wm_wdsm_fitted_variance(object$fit, object$nuisance_influence,
    object$smooth_sensitivity, object$graph_sensitivity, "critical_planar", object$conf.level)
  .wm_reciprocal_agree(expected$base_rows, object$base_rows, "Critical complete base rows")
  .wm_reciprocal_agree(expected$Sigma, object$raw_Sigma, "Critical raw covariance")
  .wm_reciprocal_agree(expected$C, object$raw_C, "Critical raw cross block")
  .wm_reciprocal_agree(expected$V0, object$diagonal_root_n_variance, "Critical diagonal block")
  .wm_reciprocal_agree(expected$total_sensitivity, object$total_sensitivity, "Critical full slope")
  .wm_reciprocal_agree(expected$augmented_rows, object$augmented_rows, "Critical full augmented rows")
  .wm_reciprocal_agree(expected$nuisance_influence, object$nuisance_influence,
    "Critical complete centered influence")
  .wm_reciprocal_agree(object$analysis_target_weight_mean, object$fit$gamma,
    "Critical analysis target denominator")
  .wm_reciprocal_agree(object$raw_target_weight_mean, mean(object$fit$weights),
    "Critical raw target denominator")
  derivative <- object$smooth_derivative
  .wm_reciprocal_agree(derivative$sensitivity, object$smooth_sensitivity,
    "Critical smooth derivative binding")
  if (any(derivative$weight_sensitivity != 0) ||
      !identical(derivative$parameter_names, object$parameter_names) ||
      !identical(derivative$estimate, object$estimate)) {
    stop("Critical complete known-weight derivative binding changed.", call. = FALSE)
  }
  transports <- .wm_reciprocal_list(object$graph_transport, "critical graph transports")
  if (!identical(names(transports), c("potential0", "potential1"))) {
    stop("Critical direction bindings changed.", call. = FALSE)
  }
  j <- error <- matrix(0, 2L, length(object$parameter_names))
  for (k in 1:2) {
    value <- transports[[k]]
    if (!isTRUE(value$numerically_available) || !value$mode %in% c("zero", "estimate")) {
      stop("Critical graph transport is unavailable.", call. = FALSE)
    }
    j[k, ] <- .wm_fv_number(value$graph_drift, ncol(j), "critical graph drift")
    error[k, ] <- .wm_fv_number(value$quadrature_error_bound, ncol(j), "critical graph error")
    if (!identical(names(value$graph_drift), object$parameter_names) || any(error[k, ] < 0)) {
      stop("Critical transport parameter/error binding changed.", call. = FALSE)
    }
    if (value$mode == "zero") {
      if (!value$basis %in% c("fixed_map", "current_centering") || any(c(j[k, ], error[k, ]) != 0)) {
        stop("Critical zero transport declaration changed.", call. = FALSE)
      }
    } else if (!identical(value$n, object$n) || value$d != 2L || value$M != object$M ||
               value$donor_arm != k - 1L || !identical(value$parameter_names, object$parameter_names)) {
      stop("Critical estimated transport source binding changed.", call. = FALSE)
    }
  }
  .wm_reciprocal_agree(as.vector(crossprod(c(-1, 1), j)) / object$raw_target_weight_mean,
    object$graph_sensitivity, "Critical normalized full graph drift")
  .wm_reciprocal_agree(colSums(error) / object$raw_target_weight_mean,
    object$graph_sensitivity_quadrature_error_bound, "Critical normalized graph error")
  support <- .wm_cp_covariance(expected)
  for (name in c("Sigma", "Sigma_plus", "square_root", "C", "q", "projector", "regression")) {
    .wm_reciprocal_agree(support[[name]], object$critical$support[[name]],
      paste("Critical support", name))
  }
  invisible(object)
}

.wm_cp_bootstrap <- function(object, B, seed, conf.level, interval, chunk_size, counts) {
  if (!is.null(counts)) stop("critical_planar rejects counts; use full nuisance/Schur-mixture draws.", call. = FALSE)
  if (length(interval) > 1L && identical(interval, c("normal", "basic", "none"))) interval <- "basic"
  interval <- match.arg(interval, c("basic", "none"))
  B <- .wm_cp_integer(B, "B"); chunk_size <- .wm_cp_integer(chunk_size, "chunk_size")
  .wm_numeric_vector(conf.level, 1L, "conf.level")
  if (conf.level <= 0 || conf.level >= 1) stop("conf.level must lie in (0,1).", call. = FALSE)
  if (!is.null(seed)) {
    .wm_numeric_vector(seed, 1L, "seed")
    if (seed < 0 || seed > .Machine$integer.max || seed != floor(seed)) stop("Invalid seed.", call. = FALSE)
    if (identical(RNGkind()[2L], "Box-Muller")) stop("An explicit seed cannot preserve Box-Muller cache.", call. = FALSE)
  }
  .wm_cp_check(object)
  p <- length(object$parameter_names); c <- object$critical_control; kernel <- object$critical$kernel
  if (B * as.double(length(kernel$pairs)) * kernel$rule$order^2 > c$maximum_geometric_evaluations ||
      B * as.double(p + 12) > c$maximum_draw_entries) {
    stop("Critical all-B geometry/output budget exceeded before RNG; no draws are filtered.", call. = FALSE)
  }
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    previous <- if (had_seed) get(".Random.seed", .GlobalEnv) else NULL
    on.exit({
      if (had_seed) assign(".Random.seed", previous, .GlobalEnv)
      else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(as.integer(seed))
  }
  support <- object$critical$support
  beta <- object$total_sensitivity + support$regression
  beta_variance <- sum(as.vector(crossprod(support$square_root, beta))^2)
  roots <- rep(NA_real_, B); H <- matrix(NA_real_, p, B,
    dimnames = list(object$parameter_names, NULL))
  diagnostic <- data.frame(draw = seq_len(B), status = rep("not_evaluated", B),
    R = rep(NA_real_, B), raw_V0 = rep(NA_real_, B), raw_D = rep(NA_real_, B),
    D = rep(NA_real_, B), cone_adjustment = rep(NA_real_, B),
    negative_raw_D = rep(NA, B), quadrature_error_bound = rep(NA_real_, B),
    G_underflow_count = rep(NA_integer_, B), error = rep("", B))
  for (b in seq_len(B)) {
    multipliers <- numeric(p + 1L)
    for (block in .wm_cp_chunks(p + 1L, chunk_size)) multipliers[block] <- stats::rnorm(length(block))
    H[, b] <- as.vector(support$square_root %*% multipliers[seq_len(p)])
    value <- tryCatch(.wm_cp_schur(object, H[, b]), error = identity)
    if (inherits(value, "error")) {
      diagnostic$status[b] <- "numerical_evaluation_failed"
      diagnostic$error[b] <- conditionMessage(value)
      next
    }
    roots[b] <- sum(beta * H[, b]) + sqrt(value$D) * multipliers[p + 1L]
    for (name in c("R", "raw_V0", "raw_D", "D", "cone_adjustment", "negative_raw_D",
                   "quadrature_error_bound", "G_underflow_count")) diagnostic[[name]][b] <- value[[name]]
    diagnostic$status[b] <- if (is.finite(roots[b])) "complete" else "root_range_failed"
  }
  complete <- all(diagnostic$status == "complete")
  draws <- object$estimate + roots / sqrt(object$n)
  complete <- complete && all(is.finite(draws))
  mc <- if (complete && B > 1L) stats::var(roots) else NA_real_
  rb <- if (complete) mean(diagnostic$D) + beta_variance else NA_real_
  rb_se <- if (complete && B > 1L) stats::sd(diagnostic$D) / sqrt(B) else NA_real_
  if (complete && (any(!is.finite(c(beta_variance, rb))) ||
      (B > 1L && any(!is.finite(c(mc, rb_se)))))) {
    complete <- FALSE; mc <- rb <- rb_se <- NA_real_
  }
  alpha <- 1 - conf.level
  ci <- if (interval == "none") NULL else if (complete) object$estimate -
    as.numeric(stats::quantile(roots, c(1 - alpha / 2, alpha / 2), names = FALSE)) /
      sqrt(object$n) else rep(NA_real_, 2L)
  if (!is.null(ci)) names(ci) <- c("lower", "upper")
  list(root_n_draws = roots, draws = draws, nuisance_draws = H,
    draw_diagnostics = diagnostic, available = complete,
    conditional_root_n_variance = NA_real_, conditional_variance = NA_real_,
    conditional_variance_formula = "E*[Vhat0(H*)] + 2 b'C + b'Sigma b (unevaluated expectation)",
    monte_carlo_root_n_variance = mc, monte_carlo_variance = mc / object$n,
    rao_blackwell_root_n_variance_estimate = rb,
    rao_blackwell_root_n_variance_mcse = rb_se,
    rao_blackwell_variance_estimate = rb / object$n,
    rao_blackwell_variance_mcse = rb_se / object$n,
    conf.int = ci, interval = interval, conf.level = conf.level, B = B, n = object$n,
    seed = seed, estimate = object$estimate, estimand = object$estimand,
    covariance_scope = "critical_planar", source_inference = object,
    method = "Complete nuisance Gaussian and conditional Schur-mixture replication",
    draw_input = "generated_full_nuisance_schur_mixture",
    original_refit_bootstrap = FALSE, population_assumptions_verified = FALSE,
    variance_divisor = "B-1 for complete root-draw sample variance",
    replication_contract = paste("Budgeted regular rank-two baseline sheet sampling law, complete",
      "F1-F4, valid guard and vanishing numerical-error sequence, justified smooth-stack support",
      "rate, complete covariance/slope consistency and strict-CDF quantile premises remain.",
      "Raw signed kernel/Schur/cone diagnostics and every B slot are retained.",
      "Rao-Blackwell averages and root sample variances are Monte Carlo estimates, not exact",
      "finite-B conditional variance. No Normal calibration, rematching, refit or survivor filtering."))
}
