# Raw scalar matching with a separately qualified explicit variance.
# The moment formula concerns the reviewed bounded scalar logistic subclass.
# Numerical fitting checks below do not certify its exact-root bridge.

.wm_scalar_weight_model <- function(model) {
  required <- c("w0", "w1", "dw0", "dw1", "domain")
  if (!is.list(model) || !all(required %in% names(model)) ||
      anyDuplicated(names(model)) || any(!names(model) %in% c(required, "label"))) {
    stop("weight_model requires w0, w1, dw0, dw1 and domain, with optional label.",
         call. = FALSE)
  }
  for (name in required[1:4]) {
    if (!is.function(model[[name]])) stop(name, " must be a known function.", call. = FALSE)
  }
  domain <- model$domain
  if (!is.numeric(domain) || is.complex(domain) || !is.null(dim(domain)) ||
      length(domain) != 2L || anyNA(domain) || any(!is.finite(domain)) ||
      domain[1L] <= 0 || domain[2L] >= 1 || domain[1L] >= domain[2L]) {
    stop("weight_model$domain must be an ordered interval strictly inside (0, 1).",
         call. = FALSE)
  }
  if (!is.null(model$label) &&
      (!is.character(model$label) || length(model$label) != 1L ||
       is.na(model$label) || !nzchar(model$label))) {
    stop("weight_model$label must be one nonempty string.", call. = FALSE)
  }
  model
}

.wm_scalar_nw <- function(e, Y, Z, bandwidth, maximum_entries = 1000000L) {
  n <- length(e)
  out <- matrix(NA_real_, n, 2L, dimnames = list(NULL, c("mean0", "mean1")))
  denominators <- out
  for (z in 0:1) {
    donors <- which(Z == z)
    # Chunking bounds temporary storage; no n-by-n matrix is retained.
    block <- max(1L, floor(maximum_entries / length(donors)))
    for (start in seq.int(1L, n, by = block)) {
      rows <- seq.int(start, min(n, start + block - 1L))
      kernel <- 1 - abs(outer(e[donors], e[rows], "-")) / bandwidth
      kernel[kernel < 0] <- 0
      denominator <- colSums(kernel)
      if (any(!is.finite(denominator)) || any(denominator <= 0)) {
        stop("An own-arm triangular-kernel denominator is zero or nonfinite.",
             call. = FALSE)
      }
      prediction <- drop(crossprod(Y[donors], kernel)) / denominator
      if (any(!is.finite(prediction))) stop("Kernel predictions exceeded numerical range.",
                                           call. = FALSE)
      out[rows, z + 1L] <- prediction
      denominators[rows, z + 1L] <- denominator
    }
  }
  list(prediction = out, denominator_range = range(denominators),
       bandwidth = bandwidth, kernel = "triangular", weighted_regression = FALSE,
       point_estimator_uses_predictions = FALSE)
}

.wm_scalar_interval <- function(estimate, root_variance, n, conf.level) {
  raw_variance <- root_variance / n
  raw_se <- critical <- NA_real_
  raw_interval <- c(lower = NA_real_, upper = NA_real_)
  status <- "conditional_formula"
  if (root_variance <= 0) {
    status <- "nonpositive_root_variance"
  } else if (!is.finite(raw_variance) || raw_variance <= 0) {
    status <- "unrepresentable_sampling_variance"
  } else {
    raw_se <- sqrt(raw_variance)
    critical <- stats::qnorm((1 - conf.level) / 2, lower.tail = FALSE)
    if (!is.finite(raw_se) || raw_se <= 0 ||
        !is.finite(critical) || critical <= 0) {
      status <- "unrepresentable_interval_scale"
    } else {
      raw_interval <- estimate + c(lower = -1, upper = 1) * critical * raw_se
      if (any(!is.finite(raw_interval)) || raw_interval[1L] >= raw_interval[2L]) {
        status <- "unrepresentable_interval"
      }
    }
  }
  available <- identical(status, "conditional_formula")
  list(status = status, root_n_variance = root_variance,
       variance = if (available) raw_variance else NA_real_,
       se = if (available) raw_se else NA_real_,
       conf.int = if (available) raw_interval else c(lower = NA_real_, upper = NA_real_),
       conf.level = conf.level,
       computed_sampling_variance = raw_variance, computed_se = raw_se,
       computed_critical_value = critical, computed_conf.int = raw_interval)
}

.wm_scalar_moments <- function(fit, design, e, weight_model, bandwidth, conf.level,
                               bandwidth_rule = "user_supplied") {
  n <- fit$n
  Y <- fit$data$Y
  Z <- fit$data$Z
  if (ncol(design) < 3L) {
    stop("This scalar inference theorem requires at least three design columns, including the intercept.",
         call. = FALSE)
  }
  if (any(e < weight_model$domain[1L] | e > weight_model$domain[2L])) {
    stop("Fitted probabilities leave the declared weight-function domain.", call. = FALSE)
  }
  scale <- max(fit$weights)
  w <- fit$weights / scale
  sw <- sum(w)
  probability_weight <- w / sw
  known <- lapply(c("w0", "w1", "dw0", "dw1"), function(name) {
    value <- weight_model[[name]](e)
    .wm_numeric_vector(value, n, paste0("weight_model$", name),
                       positive = name %in% c("w0", "w1")) / scale
  })
  names(known) <- c("w0", "w1", "dw0", "dw1")
  if (any(!is.finite(unlist(known))) || any(known$w0 <= 0 | known$w1 <= 0)) {
    stop("Common normalization of the known weight functions is not representable.",
         call. = FALSE)
  }
  gap0 <- known$w0 + e * known$dw0
  gap1 <- known$w1 - (1 - e) * known$dw1
  if (any(gap0 <= 0 | gap1 <= 0)) {
    stop("The fitting-score branch gap is nonpositive at an evaluated probability.",
         call. = FALSE)
  }
  means <- .wm_scalar_nw(e, Y, Z, bandwidth)
  residual <- Y - means$prediction[cbind(seq_len(n), Z + 1L)]
  effect <- means$prediction[, 2L] - means$prediction[, 1L]
  f <- e * (1 - e)
  EP <- function(x) sum(probability_weight * x)
  EPvector <- function(x) drop(crossprod(design, probability_weight * x))
  moment <- function(z, F) {
    arm_probability <- if (z == 1L) e else 1 - e
    EP((Z == z) * F * residual^2 / arm_probability)
  }
  vector_moment <- function(z, F) {
    arm_probability <- if (z == 1L) e else 1 - e
    EPvector((Z == z) * F * residual / arm_probability)
  }
  J <- crossprod(design, design * (w * f)) / n
  B <- crossprod(design, design * (w * (Z - e))^2) / n
  if (any(!is.finite(c(J, B)))) stop("Nuisance matrices exceeded numerical range.",
                                     call. = FALSE)
  # Positive-definiteness is a numerical check, not a population certificate.
  chol(J)
  inverse <- solve(J)
  Sigma <- inverse %*% B %*% inverse
  a <- n / sw
  rho <- EP(Z)
  M <- fit$M
  w0 <- known$w0
  w1 <- known$w1
  F0 <- e + e^2 * (1 - e) * known$dw0 / w0 / M
  F1 <- (1 - e) - e * (1 - e)^2 * known$dw1 / w1 / M
  centered_effect <- effect - fit$estimate
  if (fit$estimand == "PATE") {
    H1 <- w1 / e + w1 * (1 - e)^2 / (2 * M * e) + w0 * (1 - e) / M
    H0 <- w0 / (1 - e) + w0 * e^2 / (2 * M * (1 - e)) + w1 * e / M
    residual_variance <- a * (moment(1L, H1) + moment(0L, H0))
    heterogeneity_variance <- a * EP((e * w1 + (1 - e) * w0) * centered_effect^2)
    b <- -vector_moment(0L, F0) - vector_moment(1L, F1)
    covariance_bracket <- vector_moment(1L, w1 * (1 - e)) +
      vector_moment(0L, w0 * e) + EPvector(f * (w1 - w0) * centered_effect)
  } else {
    H1 <- w1 * e
    H0 <- w1 * e / M + (1 + 1 / (2 * M)) * w0 * e^2 / (1 - e)
    residual_variance <- a / rho^2 * (moment(1L, H1) + moment(0L, H0))
    heterogeneity_variance <- a / rho^2 * EP(e * w1 * centered_effect^2)
    b <- -vector_moment(0L, F0) / rho
    covariance_bracket <- (vector_moment(1L, f * w1) +
      vector_moment(0L, e^2 * w0) + EPvector(f * w1 * centered_effect)) / rho
  }
  C <- drop(solve(J, covariance_bracket))
  V0 <- residual_variance + heterogeneity_variance
  cross_term <- 2 * sum(b * C)
  quadratic <- drop(crossprod(b, Sigma %*% b))
  root_variance <- V0 + cross_term + quadratic
  values <- c(V0, cross_term, quadratic, root_variance, b, C, Sigma, residual,
              effect, rho, a, F0, F1, H0, H1)
  if (any(!is.finite(values))) stop("The explicit variance calculation exceeded numerical range.",
                                    call. = FALSE)
  interval <- .wm_scalar_interval(fit$estimate, root_variance, n, conf.level)
  names(b) <- names(C) <- colnames(design)
  c(interval, list(
       blocks = list(V0 = V0, residual = residual_variance,
         heterogeneity = heterogeneity_variance, cross = cross_term,
         nuisance = quadratic, b = b, C = C, J = J, B = B, Sigma = Sigma),
       moments = list(a = a, rho = rho, weight_scale = scale,
         weight_units = "All weight functions and observed weights divided by weight_scale",
         normalizer = "Full observed n for both PATE and PATT"),
       predictions = c(means, list(bandwidth_rule = bandwidth_rule,
         user_bandwidth_rate_verified = FALSE,
         rate_requirements = c("h tends to zero", "n*h/log(n) tends to infinity",
                               "n^(-1/2)*h^(-2) tends to zero"))),
       diagnostics = list(branch_gap_range = range(gap0, gap1),
         branch_gap_is_uniformly_verified = FALSE,
         supplied_vs_fitted_weight_max_difference =
           max(abs(w - ifelse(Z == 1L, w1, w0))),
         weight_difference_is_not_a_rejection_rule = TRUE,
         numerical_J_rcond = rcond(J)),
       scope = "Known smooth score-dependent weights; reviewed bounded ordinary-logistic scalar model",
       model_contract = "bounded_scalar_logistic_known_score_weights",
       weight_model = list(label = weight_model$label, domain = weight_model$domain),
       population_assumptions_verified = FALSE,
       numerical_root_bridge = "not_certified",
       bootstrap_validity_claimed = FALSE))
}

# The raw point and explicit scalar variance stay separate.
wm_scalar_logistic_match <- function(Y, Z, weights, design, M = 3L,
                                    estimand = c("PATE", "PATT"),
                                    weight_model = NULL, inference = TRUE,
                                    bandwidth = NULL, conf.level = 0.95,
                                    maxit = 100L, score_tolerance = 1e-7) {
  estimand <- match.arg(estimand)
  if (!is.numeric(Y) || is.complex(Y) || !is.null(dim(Y)) || length(Y) < 3L) {
    stop("Y must be a numeric vector with at least three rows.", call. = FALSE)
  }
  n <- length(Y)
  Y <- .wm_numeric_vector(Y, n, "Y")
  weights <- .wm_numeric_vector(weights, n, "weights", positive = TRUE)
  if ((!is.numeric(Z) && !is.logical(Z)) || is.complex(Z) || !is.null(dim(Z)) ||
      length(Z) != n || anyNA(Z) || any(!Z %in% c(0, 1)) ||
      length(unique(Z)) != 2L) {
    stop("Z must be binary and contain both treatment arms.", call. = FALSE)
  }
  Z <- as.integer(Z)
  design_names <- colnames(design)
  design <- .wm_score_matrix(design, n, "design")
  if (ncol(design) >= n || any(design[, 1L] != 1) ||
      qr(design)$rank != ncol(design)) {
    stop("design must have a leading intercept and full rank, with fewer columns than rows.",
         call. = FALSE)
  }
  if (is.null(design_names)) design_names <- paste0("L", seq_len(ncol(design)))
  if (anyNA(design_names) || any(!nzchar(design_names)) || anyDuplicated(design_names)) {
    stop("Design column names must be nonempty and unique.", call. = FALSE)
  }
  colnames(design) <- design_names
  donor_limit <- if (estimand == "PATE") min(table(Z)) else sum(Z == 0L)
  if (!is.numeric(M) || is.complex(M) || length(M) != 1L || is.na(M) ||
      !is.finite(M) || M < 1 || M != floor(M) || M > donor_limit) {
    stop("M must be a positive integer no larger than each required donor arm.",
         call. = FALSE)
  }
  if (!is.logical(inference) || length(inference) != 1L || is.na(inference)) {
    stop("inference must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(conf.level) || is.complex(conf.level) || length(conf.level) != 1L ||
      !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("conf.level must lie strictly between zero and one.", call. = FALSE)
  }
  if (!is.numeric(maxit) || is.complex(maxit) || length(maxit) != 1L || !is.finite(maxit) ||
      maxit < 1 || maxit != floor(maxit) || maxit > .Machine$integer.max ||
      !is.numeric(score_tolerance) || is.complex(score_tolerance) || length(score_tolerance) != 1L ||
      !is.finite(score_tolerance) || score_tolerance <= 0) {
    stop("Supply finite positive solver controls and an integer maxit.", call. = FALSE)
  }
  bandwidth_rule <- if (is.null(bandwidth)) "n^(-1/5)" else "user_supplied"
  if (is.null(bandwidth)) bandwidth <- n^(-1 / 5)
  if (!is.numeric(bandwidth) || is.complex(bandwidth) || length(bandwidth) != 1L ||
      !is.finite(bandwidth) || bandwidth <= 0) {
    stop("bandwidth must be finite and positive.", call. = FALSE)
  }
  # Invalid arguments fail before fitting. Runtime inference failures below
  # preserve an already computed point; malformed arguments are not such failures.
  if (inference) weight_model <- .wm_scalar_weight_model(weight_model)
  model_record <- if (inference) list(label = weight_model$label,
                                      domain = weight_model$domain) else NULL
  warnings <- character()
  capture <- function(expr) withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  data <- data.frame(.response = Z)
  data$.design <- I(design)
  ps <- tryCatch(capture(wdsm_fit_ps(.response ~ 0 + .design, data, weights,
                                    maxit = as.integer(maxit),
                                    score_tolerance = score_tolerance)),
                 error = function(e) list(error = conditionMessage(e)))
  if (!is.null(ps$error)) {
    return(structure(list(status = "propensity_fit_failed", estimate = NA_real_,
      fit = NULL, propensity = ps, inference = NULL, warnings = warnings,
      design_columns = design_names, weight_model = model_record,
      bandwidth_rule = bandwidth_rule, numerical_root_bridge = "not_certified"),
      class = c("wm_scalar_logistic_match", "list")))
  }
  point <- tryCatch(capture(wm_match(Y, Z, weights,
      scores0 = matrix(ps$probability, ncol = 1L), M = as.integer(M),
      estimand = estimand, method = "self_normalized", variance = FALSE)),
      error = function(e) list(error = conditionMessage(e)))
  if (!is.null(point$error)) {
    return(structure(list(status = "matching_failed", estimate = NA_real_,
      fit = NULL, matching_error = point$error, propensity = ps, inference = NULL,
      warnings = warnings, design_columns = design_names,
      weight_model = model_record, bandwidth_rule = bandwidth_rule,
      numerical_root_bridge = "not_certified"),
      class = c("wm_scalar_logistic_match", "list")))
  }
  explicit <- if (inference) tryCatch(capture(.wm_scalar_moments(
      point, design, ps$probability, weight_model, bandwidth, conf.level, bandwidth_rule)),
      error = function(e) list(status = "unavailable", reason = conditionMessage(e),
        root_n_variance = NA_real_, variance = NA_real_, se = NA_real_,
        conf.int = c(lower = NA_real_, upper = NA_real_),
        population_assumptions_verified = FALSE,
        numerical_root_bridge = "not_certified", bootstrap_validity_claimed = FALSE)) else
    list(status = "not_requested")
  structure(list(status = "point_computed", estimate = point$estimate,
    fit = point, propensity = ps, inference = explicit, warnings = warnings,
    design_columns = design_names, fitting_design = design,
    design_binding = "validated_input_matrix_v1", weight_model = model_record,
    bandwidth_rule = bandwidth_rule,
    numerical_root_bridge = "not_certified",
    contract = "Raw fitted probability matching; separate model-conditional explicit variance; no exact-root or refit-bootstrap certificate"),
    class = c("wm_scalar_logistic_match", "list"))
}
