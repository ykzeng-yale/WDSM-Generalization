# Internal scalar_psm branch of wm_fitted_inference. No fit, RNG or rematching.
# Numerical integration bounds are real-arithmetic bounds; not roundoff certificates.

.wm_scalar_psm_controls <- function(control, n) {
  if (is.null(control)) control <- list()
  if (!is.list(control) || length(control)) {
    control <- .wm_reciprocal_list(control, "scalar_control")
  }
  defaults <- list(weight_lower_bound = NULL, bandwidth_constant = 1,
    density_constant = 0.01, outcome_clip_constant = 1,
    quadrature_error_constant = 0.01, max_nodes = 131072L,
    maximum_matrix_entries = 1000000L, maximum_laplace_terms = 2000000000,
    design = NULL, root_certificate = NULL)
  if (any(!names(control) %in% names(defaults))) {
    stop("Unsupported scalar_control field.", call. = FALSE)
  }
  for (name in names(control)) defaults[name] <- control[name]
  for (name in setdiff(names(defaults), c("design", "root_certificate"))) {
    value <- defaults[[name]]
    if (is.null(value) && name == "weight_lower_bound") next
    .wm_numeric_vector(value, 1L, paste("scalar_control", name), positive = TRUE)
    if (name %in% c("max_nodes", "maximum_matrix_entries", "maximum_laplace_terms") &&
        (value != floor(value) || value > .Machine$integer.max ||
         value < if (name == "max_nodes") 8 else 1)) {
      stop(name, " must be a positive representable integer (max_nodes >= 8).",
           call. = FALSE)
    }
  }
  defaults$bandwidth <- defaults$bandwidth_constant * n^(-1 / 12)
  defaults$density_threshold <- defaults$density_constant * n^(-1 / 256)
  defaults$outcome_clip <- defaults$outcome_clip_constant * n^(1 / 8)
  # This is a direct bound on b, after empirical amplification. Any prescribed
  # sequence tending to zero suffices; no density-threshold factor is needed.
  defaults$drift_error_target <- defaults$quadrature_error_constant / sqrt(n)
  if (any(!is.finite(unlist(defaults[c("bandwidth", "density_threshold",
      "outcome_clip", "drift_error_target")]))) ||
      any(unlist(defaults[c("bandwidth", "density_threshold",
        "outcome_clip", "drift_error_target")]) <= 0)) {
    stop("Scalar smoothing controls exceed numerical range.", call. = FALSE)
  }
  defaults
}

.wm_scalar_psm_bind <- function(object, control) {
  fit <- .wm_reciprocal_list(object$fit, "scalar fit")
  if (!inherits(fit, "wm_match") || !identical(fit$method, "self_normalized") ||
      !identical(fit$info$corrected, FALSE) ||
      !identical(fit$info$nuisance_correction, FALSE)) {
    stop("scalar_psm requires an unadjusted raw self_normalized point.", call. = FALSE)
  }
  n <- .wm_fit_integer(fit$n, "fit$n", 3L)
  M <- .wm_fit_integer(fit$M, "fit$M", 1L)
  if (!is.character(fit$estimand) || length(fit$estimand) != 1L ||
      is.na(fit$estimand) || !fit$estimand %in% c("PATE", "PATT")) {
    stop("Invalid scalar estimand.", call. = FALSE)
  }
  Y <- .wm_numeric_vector(fit$data$Y, n, "stored Y")
  Z <- .wm_numeric_vector(as.numeric(fit$data$Z), n, "stored Z")
  if (any(!Z %in% c(0, 1)) || length(unique(Z)) != 2L) {
    stop("Stored Z must contain both binary arms.", call. = FALSE)
  }
  pate <- identical(fit$estimand, "PATE")
  if (M > if (pate) min(table(Z)) else sum(Z == 0)) {
    stop("Stored M exceeds the required donor count.", call. = FALSE)
  }
  for (name in c("cell", "fold_id", "strata")) {
    value <- fit$graph[[name, exact = TRUE]]
    if (length(value) != n || anyNA(value) || length(unique(value)) != 1L) {
      stop("scalar_psm only accepts unrestricted full-sample graphs.", call. = FALSE)
    }
  }
  if (!identical(fit$graph$tie_rule, "Exact distance, then original row index")) {
    stop("Unsupported stored scalar tie rule.", call. = FALSE)
  }
  if (!is.list(fit$predictions)) stop("Missing raw prediction records.")
  mean0 <- .wm_numeric_vector(fit$predictions$mean0, n, "raw mean0")
  mean1 <- if (pate) .wm_numeric_vector(fit$predictions$mean1, n, "raw mean1") else NULL
  if (any(mean0 != 0) || (pate && any(mean1 != 0)) ||
      (!pate && !is.null(fit$predictions$mean1)) ||
      !is.null(fit$predictions$rho0) || !is.null(fit$predictions$rho1)) {
    stop("Raw scalar fit must not contain correction or rho predictions.", call. = FALSE)
  }
  W <- .wm_numeric_vector(fit$weights, n, "stored weights", TRUE)
  w <- .wm_numeric_vector(fit$analysis_weights, n, "analysis weights", TRUE)
  scale <- .wm_numeric_vector(fit$weight_scale, 1L, "weight scale", TRUE)
  if (!identical(as.numeric(scale), max(W)) ||
      !identical(as.numeric(w), as.numeric(W / scale))) {
    stop("Stored raw and analysis weights do not have the exact common scale.",
         call. = FALSE)
  }
  e <- .wm_numeric_vector(object$propensity$probability, n, "stored probability")
  if (any(e <= 0 | e >= 1)) stop("Stored probabilities must lie in (0,1).")
  if (!identical(object$estimate, fit$estimate)) stop("Stored raw point records differ.")
  incoming <- .wm_scalar_repl_graph(fit, e)
  # Validate the existing graph against the recorded distance/row-priority rule.
  # This check constructs no replacement graph and changes no donor identity.
  queries <- if (pate) seq_len(n) else which(Z == 1L)
  donor_cache <- list(which(Z == 0L), which(Z == 1L))
  for (i in queries) {
    donors <- donor_cache[[2L - Z[i]]]
    score <- if (Z[i] == 1L) fit$graph$scores0 else fit$graph$scores1
    distance <- .wm_distances(score, i, donors)
    selected <- match(fit$graph$neighbors[[i]], donors)
    d_selected <- distance[selected]
    j_selected <- donors[selected]
    if (!identical(order(d_selected, j_selected, method = "radix"), seq_len(M))) {
      stop("Stored donor ordering violates the distance/row-priority rule at row ", i,
           ".", call. = FALSE)
    }
    excluded <- rep(TRUE, length(donors))
    excluded[selected] <- FALSE
    if (any(distance[excluded] < d_selected[M] |
        (distance[excluded] == d_selected[M] & donors[excluded] < j_selected[M]))) {
      stop("Stored donors are not the recorded nearest donors at row ", i, ".",
           call. = FALSE)
    }
  }
  if (!identical(object$design_binding, "validated_input_matrix_v1") ||
      is.null(object$fitting_design)) {
    stop("Missing original fitting-design binding; a replacement design is insufficient.",
         call. = FALSE)
  }
  X <- object$fitting_design
  .wm_score_matrix(X, n, "stored fitting design")
  parameters <- object$design_columns
  if (!is.character(parameters) || anyNA(parameters) || any(!nzchar(parameters)) ||
      anyDuplicated(parameters) || !identical(colnames(X), parameters) ||
      length(parameters) != ncol(X) || ncol(X) < 2L || ncol(X) >= n ||
      any(X[, 1L] != 1) || qr(X)$rank != ncol(X)) {
    stop("Scalar PSM needs the bound full-rank design with p>=2 and a leading intercept.",
         call. = FALSE)
  }
  if (!is.null(control$design) && !identical(control$design, object$fitting_design)) {
    stop("Supplied design differs from the exact stored fitting design.", call. = FALSE)
  }
  # QR represents the stored logits; it does not fit or replace the propensity.
  logit <- stats::qlogis(e)
  coefficients <- qr.solve(X, logit)
  represented_logit <- drop(X %*% coefficients)
  if (any(!is.finite(represented_logit)) ||
      any(abs(represented_logit - logit) > 1e-10 * pmax(1, abs(logit)))) {
    stop("Stored probability logits are not in the bound logistic design space.")
  }
  diagnostics <- object$propensity$diagnostics
  tolerance <- .wm_numeric_vector(diagnostics$normalized_score_tolerance, 1L,
                                  "original normalized score tolerance", TRUE)
  score_scale <- pmax(colSums(w * abs(X)), .Machine$double.eps * sum(w))
  normalized_score <- drop(crossprod(X, w * (Z - e))) / score_scale
  if (!isTRUE(diagnostics$converged) ||
      !identical(diagnostics$initialization, "zero_coefficients_raw_weights") ||
      !is.finite(max(abs(normalized_score))) || max(abs(normalized_score)) > tolerance) {
    stop("Stored probabilities fail the original weighted-logistic numerical score check.")
  }
  if (!is.null(control$weight_lower_bound) &&
      control$weight_lower_bound > min(W)) {
    stop("Declared weight_lower_bound exceeds an observed raw weight.", call. = FALSE)
  }
  certificate <- control$root_certificate
  certificate_status <- list(numerical_root_status = "not_certified",
    donor_set_status = "not_certified", point_arithmetic_bound = NULL,
    bound_certificate = NULL, asymptotic_root_rate_verified = FALSE)
  if (!is.null(certificate)) {
    expected <- list(design = unname(X), Z = fit$data$Z, weights = W, probability = e,
      neighbors = fit$graph$neighbors,
      edges = fit$graph$edges[c("query", "donor", "arm")],
      Y = Y, estimate = fit$estimate, M = fit$M, estimand = fit$estimand)
    if (!inherits(certificate, "wm_scalar_root_certificate") ||
        !identical(certificate$inputs, expected)) {
      stop("Root certificate is not bound to these exact original inputs.", call. = FALSE)
    }
    certificate_status$numerical_root_status <- if (isTRUE(certificate$root_contained))
      "provided_finite_data_root_containment" else "not_certified"
    certificate_status$donor_set_status <- if (isTRUE(certificate$donor_sets_certified))
      "provided_finite_data_donor_certificate" else "not_certified"
    certificate_status$point_arithmetic_bound <- certificate$point_root_error_bound
    certificate_status$bound_certificate <- certificate
  }
  list(n = n, M = M, pate = pate, estimand = fit$estimand, Y = Y, Z = Z,
    W = W, w = w, scale = scale, e = e, X = X, parameters = parameters,
    incoming = incoming, gamma = mean(if (pate) w else Z * w),
    certificate = certificate_status,
    observed_logistic_check = list(logit_representation_max_error =
      max(abs(represented_logit - logit)), normalized_score_max = max(abs(normalized_score)),
      original_score_tolerance = tolerance, exact_root_verified = FALSE,
      stored_nearest_donors_validated = TRUE))
}

.wm_scalar_psm_kernel <- function(donor_score, evaluation_score, bandwidth) {
  u <- outer(donor_score, evaluation_score, function(x, s) (s - x) / bandwidth)
  active <- abs(u) < 1
  # Keep the matrix as pmax's first argument so row/column dimensions survive.
  v <- pmax(1 - u^2, 0)
  K <- (35 / 32) * v^3 / bandwidth
  derivative <- -(105 / 16) * u * v^2 / bandwidth^2
  K[!active] <- derivative[!active] <- 0
  if (any(!is.finite(c(K, derivative)))) stop("Kernel arithmetic exceeded numerical range.")
  list(K = K, derivative = derivative)
}

.wm_scalar_psm_smooth <- function(s, Y, Z, w, bandwidth, outcome_clip,
                                   maximum_entries, pate) {
  n <- length(s)
  a <- da <- weighted <- dweighted <- matrix(0, n, 2L,
    dimnames = list(NULL, c("arm0", "arm1")))
  mu <- dmu <- matrix(NA_real_, n, if (pate) 2L else 1L,
    dimnames = list(NULL, if (pate) c("mean0", "mean1") else "mean0"))
  Yclip <- pmax(-outcome_clip, pmin(outcome_clip, Y))
  for (z in 0:1) {
    donor <- which(Z == z)
    block <- max(1L, floor(maximum_entries / (2 * length(donor))))
    for (start in seq.int(1L, n, by = block)) {
      rows <- seq.int(start, min(n, start + block - 1L))
      kernel <- .wm_scalar_psm_kernel(s[donor], s[rows], bandwidth)
      den <- colSums(kernel$K)
      dp <- colSums(kernel$derivative)
      a[rows, z + 1L] <- den / n
      da[rows, z + 1L] <- dp / n
      weighted[rows, z + 1L] <- drop(crossprod(w[donor], kernel$K)) / n
      dweighted[rows, z + 1L] <- drop(crossprod(w[donor], kernel$derivative)) / n
      if (z == 0L || pate) {
        usable <- den > 0
        level <- drop(crossprod(Yclip[donor], kernel$K))
        prime <- drop(crossprod(Yclip[donor], kernel$derivative))
        mu[rows[usable], z + 1L] <- level[usable] / den[usable]
        dmu[rows[usable], z + 1L] <-
          (prime[usable] - mu[rows[usable], z + 1L] * dp[usable]) / den[usable]
      }
    }
  }
  if (any(!is.finite(c(a, da, weighted, dweighted)))) {
    stop("Scalar subdensity arithmetic exceeded numerical range.")
  }
  list(density = a, density_derivative = da, weighted_density = weighted,
    weighted_density_derivative = dweighted, mean = mu, mean_derivative = dmu,
    bandwidth = bandwidth, outcome_clip = outcome_clip,
    kernel = "triweight 35/32*(1-u^2)^3, abs(u)<1; compact C2")
}

.wm_scalar_psm_probabilities <- function(s, rows, donor, bandwidth) {
  k <- .wm_scalar_psm_kernel(s[donor], s[rows], bandwidth)
  den <- colSums(k$K)
  if (any(!is.finite(den) | den <= 0)) stop("An active donor kernel has zero mass.")
  p <- sweep(k$K, 2L, den, "/")
  dp <- sweep(k$derivative - sweep(p, 2L, colSums(k$derivative), "*"),
              2L, den, "/")
  if (any(!is.finite(c(p, dp)))) stop("Conditional mark probabilities are nonfinite.")
  list(probability = p, derivative = dp)
}

# Upper envelope calculations include a positive underflow floor, but ordinary
# binary64 roundoff is explicitly outside these real-arithmetic error bounds.
.wm_scalar_psm_positive_exp <- function(log_value) {
  if (anyNA(log_value) || any(log_value > log(.Machine$double.xmax))) {
    stop("Integration envelope exceeds numerical range.")
  }
  pmax(exp(log_value), .Machine$double.xmin)
}

.wm_scalar_psm_quad_bounds <- function(root_w, donor_min, donor_max, M,
                                       derivative_tv, cutoff, panels) {
  lower <- root_w + (M - 1) * donor_min
  upper <- root_w + (M - 1) * donor_max
  h <- cutoff / panels
  log_geom <- log(-expm1(-lower * cutoff)) - log(-expm1(-lower * h))
  # Order-eight Gauss-Legendre: (8!)^4 / {17 (16!)^3} times
  # panel_length^17 sup |f^(16)|. Mixture rates bound the derivative.
  log_constant <- 4 * lgamma(9) - log(17) - 3 * lgamma(17)
  error <- .wm_scalar_psm_positive_exp(log(root_w) + 16 * log(upper) +
                                        17 * log(h) + log_constant + log_geom)
  tail <- .wm_scalar_psm_positive_exp(log(root_w) - log(lower) - lower * cutoff)
  R <- (M - 1) * derivative_tv
  list(A = error, A_s = R * error, tail_A = tail, tail_A_s = R * tail,
       lambda_min = lower, lambda_max = upper)
}

.wm_scalar_psm_laplace_values <- function(donor_w, probability, derivative,
                                         root_w, M, cutoff, panels,
                                         maximum_entries) {
  q <- length(root_w)
  if (M == 1L) return(list(A = rep(1, q), A_s = numeric(q)))
  if (length(unique(donor_w)) == 1L) {
    return(list(A = root_w / (root_w + (M - 1) * donor_w[1L]),
                A_s = numeric(q)))
  }
  answer <- answer_s <- numeric(q)
  block <- max(1L, floor(maximum_entries / (length(donor_w) + 4 * q)))
  # These binary64 nodes/weights approximate the ideal Gauss rule. As with all
  # matrix arithmetic here, their roundoff is excluded from the analytic bound.
  positive_node <- c(0.18343464249564980494, 0.52553240991632898582,
                     0.79666647741362673959, 0.96028985649753623168)
  positive_weight <- c(0.36268378337836198297, 0.31370664587788728734,
                       0.22238103445337447054, 0.10122853629037625915)
  node <- c(-rev(positive_node), positive_node)
  weight <- c(rev(positive_weight), positive_weight)
  spacing <- cutoff / panels
  nodes <- 8 * panels
  for (start in seq.int(0, nodes - 1, by = block)) {
    index <- seq.int(start, min(nodes - 1, start + block - 1))
    within_panel <- index %% 8 + 1L
    t <- spacing * (index %/% 8 + (node[within_panel] + 1) / 2)
    exp_donor <- exp(-outer(donor_w, t))
    phi <- crossprod(probability, exp_donor)
    phi_s <- crossprod(derivative, exp_donor)
    exp_root <- exp(-outer(root_w, t))
    quadrature_weight <- spacing / 2 * weight[within_panel]
    integrand <- root_w * exp_root * phi^(M - 1)
    integrand_s <- (M - 1) * root_w * exp_root * phi^(M - 2) * phi_s
    answer <- answer + drop(integrand %*% quadrature_weight)
    answer_s <- answer_s + drop(integrand_s %*% quadrature_weight)
  }
  if (any(!is.finite(c(answer, answer_s)))) stop("Laplace arithmetic is nonfinite.")
  list(A = answer, A_s = answer_s)
}

.wm_scalar_psm_arm <- function(z, bound, smooth, active, control, cutoff,
                               remaining_terms, error_target) {
  n <- bound$n
  M <- bound$M
  X <- bound$X
  p <- ncol(X)
  donor <- which(bound$Z == z)
  rows <- which(active & bound$Z == z)
  query <- which(active & bound$Z != z)
  derivative <- bound$e * (1 - bound$e) * X
  m <- smooth$mean[, z + 1L]
  mp <- smooth$mean_derivative[, z + 1L]
  a <- smooth$density[, z + 1L]
  J <- Jp <- rep(NA_real_, n)
  J[active] <- smooth$weighted_density[active, 2L - z] / a[active]
  Jp[active] <- (smooth$weighted_density_derivative[active, 2L - z] -
    J[active] * smooth$density_derivative[active, z + 1L]) / a[active]
  if (any(!is.finite(c(m[active], mp[active], J[active], Jp[active])))) {
    stop("Active scalar conditional moments are nonfinite.")
  }
  residual <- bound$Y[rows] - m[rows]
  # These coefficients already include the final target denominator.
  coefficient_A <- M / (n * bound$gamma) * abs(derivative[rows, , drop = FALSE]) *
    abs(residual * Jp[rows] - J[rows] * mp[rows])
  coefficient_As <- M / (n * bound$gamma) * abs(derivative[rows, , drop = FALSE]) *
    abs(J[rows] * residual)
  if (any(!is.finite(c(coefficient_A, coefficient_As)))) {
    stop("Drift error amplification is nonfinite.")
  }
  zero <- stats::setNames(numeric(p), bound$parameters)
  A <- As <- rep(NA_real_, n)
  exact <- M == 1L || length(unique(bound$w[donor])) == 1L || !length(rows)
  method <- if (!length(rows)) "no_active_donor_rows" else if (M == 1L)
    "M1_exact" else if (exact) "constant_donor_weights_exact" else
    "shared_empirical_Laplace_composite_Gauss8"
  tv <- numeric(length(rows))
  panels <- 0L
  quadrature_bound <- tail_bound <- zero
  required_terms <- 0
  envelope <- NULL
  if (exact) {
    A[rows] <- if (M == 1L) 1 else
      bound$w[rows] / (bound$w[rows] + (M - 1) * bound$w[donor[1L]])
    As[rows] <- 0
  } else {
    block <- floor(control$maximum_matrix_entries / (2 * length(donor)))
    for (start in seq.int(1L, length(rows), by = block)) {
      index <- seq.int(start, min(length(rows), start + block - 1L))
      probabilities <- .wm_scalar_psm_probabilities(bound$e, rows[index], donor,
                                                   control$bandwidth)
      tv[index] <- colSums(abs(probabilities$derivative))
    }
    panels <- 1L
    tail_doublings <- 0L
    repeat {
      envelope <- .wm_scalar_psm_quad_bounds(bound$w[rows], min(bound$w[donor]),
        max(bound$w[donor]), M, tv, cutoff, panels)
      quadrature_bound <- colSums(coefficient_A * envelope$A +
                                   coefficient_As * envelope$A_s)
      tail_bound <- colSums(coefficient_A * envelope$tail_A +
                             coefficient_As * envelope$tail_A_s)
      if (any(!is.finite(c(quadrature_bound, tail_bound)))) {
        stop("Propagated integration bounds are nonfinite.")
      }
      if (max(tail_bound) > error_target / 2 && tail_doublings < 32L &&
          is.finite(2 * cutoff)) {
        cutoff <- 2 * cutoff
        tail_doublings <- tail_doublings + 1L
        next
      }
      if (max(quadrature_bound + tail_bound) <= error_target) break
      if (max(tail_bound) > error_target / 2 || panels >= floor(control$max_nodes / 8)) {
        return(list(available = FALSE, donor_arm = z, method = method,
          reason = if (max(tail_bound) > error_target / 2) "tail_error_budget" else
            "quadrature_node_budget", nodes = 8 * panels, cutoff = cutoff,
          quadrature_error_bound = quadrature_bound, tail_error_bound = tail_bound,
          drift_error_bound = quadrature_bound + tail_bound, error_target = error_target))
      }
      panels <- min(2 * panels, floor(control$max_nodes / 8))
    }
    # Two conditional Laplace matrix products, shared by A and partial_s A.
    required_terms <- 2 * length(donor) * length(rows) * (8 * panels)
    if (!is.finite(required_terms) || required_terms > remaining_terms) {
      return(list(available = FALSE, donor_arm = z, method = method,
        reason = "laplace_operation_budget", nodes = 8 * panels, cutoff = cutoff,
        required_laplace_terms = required_terms, remaining_laplace_terms = remaining_terms,
        quadrature_error_bound = quadrature_bound, tail_error_bound = tail_bound,
        drift_error_bound = quadrature_bound + tail_bound, error_target = error_target))
    }
    for (start in seq.int(1L, length(rows), by = block)) {
      index <- seq.int(start, min(length(rows), start + block - 1L))
      probabilities <- .wm_scalar_psm_probabilities(bound$e, rows[index], donor,
                                                   control$bandwidth)
      value <- .wm_scalar_psm_laplace_values(bound$w[donor],
        probabilities$probability, probabilities$derivative, bound$w[rows[index]],
        M, cutoff, panels, control$maximum_matrix_entries)
      A[rows[index]] <- value$A
      As[rows[index]] <- value$A_s
    }
  }
  query_term <- colSums(derivative[query, , drop = FALSE] *
                         (bound$w[query] * mp[query])) / n
  donor_term <- M * colSums(derivative[rows, , drop = FALSE] *
    (A[rows] * (residual * Jp[rows] - J[rows] * mp[rows]) +
       J[rows] * As[rows] * residual)) / n
  slope <- query_term + donor_term
  if (any(!is.finite(c(slope, A[rows], As[rows])))) {
    stop("Integrated scalar drift is nonfinite.")
  }
  list(available = TRUE, donor_arm = z, method = method, donor_rows = rows,
    query_rows = query, A = A, A_s = As, J = J, J_s = Jp,
    query_term = query_term, donor_term = donor_term, unnormalized_slope = slope,
    quadrature_error_bound = quadrature_bound, tail_error_bound = tail_bound,
    drift_error_bound = quadrature_bound + tail_bound, error_target = error_target,
    nodes = if (exact) 0L else 8 * panels, cutoff = if (exact) NULL else cutoff,
    conditional_derivative_total_variation = tv, per_row_error_envelope = envelope,
    laplace_terms = required_terms,
    partial_derivative_holds_root_weight_fixed = TRUE)
}

.wm_scalar_psm_inference <- function(object, scalar_control, conf.level) {
  object <- .wm_reciprocal_list(object, "scalar fit object")
  if (!inherits(object, "wm_scalar_logistic_match") ||
      !identical(object$status, "point_computed") ||
      !is.numeric(object$estimate) || length(object$estimate) != 1L ||
      !is.finite(object$estimate)) {
    stop("scalar_psm requires a completed wm_scalar_logistic_match object.", call. = FALSE)
  }
  n <- .wm_fit_integer(object$fit$n, "fit$n", 3L)
  control <- .wm_scalar_psm_controls(scalar_control, n)
  conf.level <- .wm_numeric_vector(conf.level, 1L, "conf.level")
  if (conf.level <= 0 || conf.level >= 1) stop("conf.level must lie in (0,1).")
  bound <- .wm_scalar_psm_bind(object, control)
  result <- list(estimate = object$estimate, fit = object$fit, source_object = object,
    n = n, M = bound$M, estimand = bound$estimand, covariance_scope = "scalar_psm",
    parameter_names = bound$parameters, dimensions = c(potential0 = 1L,
      if (bound$pate) c(potential1 = 1L)), controls = control,
    status = "scalar_nuisance_unavailable", available = FALSE,
    conditional_inference_available = FALSE, root_n_variance = NA_real_,
    variance = NA_real_, se = NA_real_, conf.int = c(lower = NA_real_, upper = NA_real_),
    conf.level = conf.level, assumptions_verified = FALSE, application_verified = FALSE,
    population_assumptions_verified = FALSE, numerical_certificate = bound$certificate,
    observed_logistic_check = bound$observed_logistic_check,
    contract = list(
      branch = "known smooth positive W(Z,X); raw fixed-M target-logistic PSM",
      population = paste("iid observed law; fixed full-rank bounded logistic design and M;",
        "correct target propensity, overlap, target identities and true own-(score,weight)",
        "centering (control only for PATT); compact zero-extended smooth conditional",
        "design densities, smooth means/known weights and the weighted-root derivative-span",
        "condition; uniform conditional outcome moment of order 2+eta; positive limit variance."),
      numerical_root = paste("The statistical theorem concerns the regular exact root.",
        "The stored numerical root's asymptotic bridge is not certified here;",
        "an optional bound input-matched certificate is only finite-data evidence."),
      weight_lower_bound = paste("Optional initial-cutoff bound, checked against observed raw W.",
        "Otherwise the empirical minimum is used. Neither verifies population boundedness."),
      variance_only = "Smoothing, clipping and trimming never change the raw point or donor graph.",
      integration = paste("Empirical Laplace interval plus tail bound propagated directly to b;",
        "excludes floating roundoff, smoothing bias, sampling and model error.",
        "A prescribed error constant/sqrt(n) tends to zero if the numerical budget passes;",
        "finite caps do not establish an asymptotic pass probability."),
      operations = paste("Scalar smoothing and conditional-mark construction use O(n^2) work;",
        "shared Laplace products use O(n^2 * nodes) work, prospectively capped.",
        "maximum_matrix_entries bounds kernel blocks, not total process memory."),
      replication = "Complete centered rows are returned; no refit bootstrap law is asserted.",
      estimated_weights_supported = FALSE, original_refit_limit_agreement_declared = FALSE),
    completed_components = list())
  class(result) <- "wm_fitted_inference"
  stage <- "scalar_smoothing"
  tryCatch({
    if (control$maximum_matrix_entries < 2 * max(table(bound$Z))) {
      stop("maximum_matrix_entries cannot hold one donor kernel/derivative column.")
    }
    normalized_lower <- if (is.null(control$weight_lower_bound)) min(bound$w) else
      control$weight_lower_bound / bound$scale
    cutoff <- 4 * log(n) / (bound$M * normalized_lower)
    if (!is.finite(cutoff) || !is.finite(normalized_lower) || normalized_lower <= 0) {
      stop("Normalized weight bound or Laplace cutoff is unrepresentable.")
    }
    smooth <- .wm_scalar_psm_smooth(bound$e, bound$Y, bound$Z, bound$w,
      control$bandwidth, control$outcome_clip, control$maximum_matrix_entries, bound$pate)
    active <- apply(smooth$density, 1L, min) >= 2 * control$density_threshold
    result$completed_components$smoothing <- smooth
    result$active <- active
    if (!any(active)) stop("No rows meet the two-arm density trim; no zero-drift fallback.")
    stage <- "finite_M_scalar_drift"
    arms <- if (bound$pate) 0:1 else 0L
    remaining_terms <- control$maximum_laplace_terms
    components <- list()
    for (z in arms) {
      component <- .wm_scalar_psm_arm(z, bound, smooth, active, control, cutoff,
        remaining_terms, control$drift_error_target / length(arms))
      components[[paste0("arm", z)]] <- component
      result$completed_components$drift <- components
      if (!isTRUE(component$available)) stop(component$reason)
      remaining_terms <- remaining_terms - component$laplace_terms
    }
    b <- if (bound$pate) (components$arm1$unnormalized_slope -
      components$arm0$unnormalized_slope) / bound$gamma else
      -components$arm0$unnormalized_slope / bound$gamma
    error_bound <- Reduce("+", lapply(components, `[[`, "drift_error_bound"))
    if (any(!is.finite(b)) || any(error_bound > control$drift_error_target)) {
      stop("Final drift or propagated error bound is invalid.")
    }
    result$total_sensitivity <- b
    result$drift_error_bound <- error_bound
    result$graph_sensitivity_quadrature_error_bound <-
      Reduce("+", lapply(components, `[[`, "quadrature_error_bound"))
    result$drift_tail_error_bound <- Reduce("+", lapply(components, `[[`, "tail_error_bound"))
    result$laplace_terms <- control$maximum_laplace_terms - remaining_terms
    stage <- "complete_rows_and_sandwich"
    means <- matrix(0, n, if (bound$pate) 2L else 1L,
      dimnames = list(NULL, if (bound$pate) c("mean0", "mean1") else "mean0"))
    means[active, ] <- pmax(-log(n), pmin(log(n), smooth$mean[active, , drop = FALSE]))
    if (any(!is.finite(means))) stop("Feasible row means are nonfinite.")
    w <- bound$w
    Z <- bound$Z
    Y <- bound$Y
    X <- bound$X
    e <- bound$e
    H <- crossprod(X, X * (w * e * (1 - e))) / n
    chol(H)
    psi <- X * (w * (Z - e))
    influence <- t(solve(H, t(psi)))
    influence <- sweep(influence, 2L, colMeans(influence), "-")
    if (bound$pate) {
      mu <- means[cbind(seq_len(n), Z + 1L)]
      own <- bound$incoming[cbind(seq_len(n), Z + 1L)]
      raw <- (w * (means[, 2L] - means[, 1L] - object$estimate) +
        (2 * Z - 1) * (w + own) * (Y - mu)) / bound$gamma
    } else {
      raw <- (Z * w * (Y - means[, 1L] - object$estimate) -
        (1 - Z) * bound$incoming[, 1L] * (Y - means[, 1L])) / bound$gamma
    }
    base <- raw - mean(raw)
    augmented <- base + drop(influence %*% b)
    augmented <- augmented - mean(augmented)
    V0 <- mean(base^2)
    C <- drop(crossprod(influence, base)) / n
    names(C) <- bound$parameters
    Sigma <- crossprod(influence) / n
    cross_term <- 2 * sum(b * C)
    nuisance_variance <- drop(crossprod(b, Sigma %*% b))
    V <- mean(augmented^2)
    if (any(!is.finite(c(H, influence, raw, augmented, V0, C, Sigma,
                         cross_term, nuisance_variance, V)))) {
      stop("Complete row or covariance arithmetic is nonfinite.")
    }
    .wm_scalar_repl_agree(V, V0 + cross_term + nuisance_variance,
                         "complete scalar covariance identity")
    result$feasible_means <- means
    result$base_rows <- base
    result$uncentered_base_rows <- raw
    result$augmented_rows <- augmented
    result$nuisance_influence <- influence
    result$H <- H
    result$V0 <- V0
    result$C <- C
    result$Sigma <- Sigma
    result$cross_term <- cross_term
    result$nuisance_variance <- nuisance_variance
    result$gamma <- bound$gamma
    interval <- .wm_scalar_interval(object$estimate, V, n, conf.level)
    for (name in setdiff(names(interval), "status")) result[[name]] <- interval[[name]]
    result$available <- identical(interval$status, "conditional_formula")
    result$conditional_inference_available <- result$available
    result$status <- if (result$available) "conditional_scalar_psm" else interval$status
    result
  }, error = function(error) {
    result$failure_stage <- stage
    result$reason <- conditionMessage(error)
    result
  })
}
