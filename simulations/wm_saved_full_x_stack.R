# Saved-array construction of a complete full-X prediction stack.
# Source this simulation helper; it is not an exported package API.
# No model fit, matching, bootstrap, data generation or RNG is performed.
#
# data is a list/data.frame with Y, Z and weights in original fit row order.
# fit is the original corrected self_normalized wm_match object. mean_models
# contains mean0 and, for PATE only, mean1; each has its recorded design and
# coefficient vector. family="PS" additionally requires ps$design,
# ps$probability and ps$weighting ("probability" or "unit"). This source
# adapter declares the original ZERO-offset logistic model; offsets are not
# introduced or inferred. Optional original
# ps$coefficients are checked; missing coefficients are never reconstructed.
# For family="full_X", raw_X has the original coordinates. The recorded mean
# designs must be the COMPLETE polynomial basis of its pooled standardized
# coordinates. polynomial_exponents supplies that basis's ordered multi-indices;
# NULL selects the source intercept/linear/complete-quadratic column order.
#
# The full PS stack contains every logistic coefficient, its pooled center and
# variance, and all raw-design mean coefficients. The full-X stack contains
# every coordinate center/variance and all RAW polynomial mean coefficients.
# For row equations Psi and J = d(P_n Psi)/d theta, influence rows are
# -Psi J^{-T}, evaluated without solving a new statistical estimating root.
# All blocks remain in the same-row covariance even when their mean derivative
# is zero. Supplied probability weights remain known; no weight parameter exists.
# raw_X denotes matched coordinates, not necessarily every observed design mark.
# W may depend on other observed marks; full-X residual centering concerns that
# full observed design. No W-measurable-given-raw_X or alpha-only variance is used.
#
# If S_j=(X_j-c_j)/sqrt(v_j), R(S)=R(X)T with
# T[b,a]=1(b<=a) prod_j choose(a_j,b_j)(-c_j)^(a_j-b_j)v_j^(-a_j/2).
# The complete basis makes T invertible. For any fixed multiplicities and arm
# weights, beta_RAW=T beta_STD is the same WLS prediction parameter because
# (T' G T)^-1 T' r = T^-1 G^-1 r. This is prediction-basis invariance only;
# matching geometry is NOT invariant to changing coordinate scales.
#
# Returned inputs feed existing wm_fitted_inference(); this helper does not call
# it. Full-X correctness, identification, iid sampling, regular root/geometry,
# actual-row/moment laws and numerical localization remain unverified premises.
# Original unbounded/dependent survey settings remain empirical approximations.

.wm_saved_stack_function <- function(wm, name) {
  if (!is.environment(wm) || !exists(name, envir = wm, inherits = FALSE) ||
      !is.function(get(name, envir = wm, inherits = FALSE))) {
    stop("The common package environment lacks ", name, ".", call. = FALSE)
  }
  get(name, envir = wm, inherits = FALSE)
}

.wm_saved_polynomial_basis <- function(x, exponents) {
  out <- matrix(1, nrow(x), nrow(exponents))
  for (a in seq_len(nrow(exponents))) for (j in seq_len(ncol(x))) {
    if (exponents[a, j] > 0) out[, a] <- out[, a] * x[, j]^exponents[a, j]
  }
  if (any(!is.finite(out))) stop("Polynomial basis exceeded numerical range.", call. = FALSE)
  colnames(out) <- rownames(exponents)
  out
}

.wm_saved_polynomial_transform <- function(exponents, center, variance) {
  r <- nrow(exponents); d <- ncol(exponents)
  out <- matrix(0, r, r, dimnames = list(rownames(exponents), rownames(exponents)))
  for (a in seq_len(r)) for (b in seq_len(r)) {
    if (all(exponents[b, ] <= exponents[a, ])) {
      value <- 1
      for (j in seq_len(d)) {
        value <- value * choose(exponents[a, j], exponents[b, j]) *
          (-center[j])^(exponents[a, j] - exponents[b, j]) *
          variance[j]^(-exponents[a, j] / 2)
      }
      out[b, a] <- value
    }
  }
  if (any(!is.finite(out))) stop("Polynomial coefficient transform overflowed.", call. = FALSE)
  out
}

wm_saved_full_x_stack <- function(data, fit, mean_models,
                                 family = c("PS", "full_X"), ps = NULL,
                                 raw_X = NULL, polynomial_exponents = NULL,
                                 wm = NULL) {
  family <- match.arg(family)
  if (!is.list(fit) || !inherits(fit, "wm_match")) {
    stop("Supply the recorded original wm_match fit.", call. = FALSE)
  }
  if (is.null(wm)) wm <- asNamespace("wdsmatch")
  vector_check <- .wm_saved_stack_function(wm, ".wm_ns_vector")
  design_check <- .wm_saved_stack_function(wm, ".wm_ns_design")
  agree <- .wm_saved_stack_function(wm, ".wm_reciprocal_agree")
  state <- .wm_saved_stack_function(wm, ".wm_fv_state")(fit)
  n <- state$n; pate <- state$pATE
  if (!is.list(data) || !all(c("Y", "Z", "weights") %in% names(data))) {
    stop("data must provide recorded Y, Z and weights in original row order.", call. = FALSE)
  }
  Y <- vector_check(data$Y, n, "data$Y")
  Z <- vector_check(data$Z, n, "data$Z")
  W <- vector_check(data$weights, n, "data$weights")
  if (!identical(Y, as.numeric(fit$data$Y)) ||
      !identical(Z, as.numeric(fit$data$Z)) ||
      !identical(W, as.numeric(fit$weights))) {
    stop("Recorded data do not identify the original matching fit.", call. = FALSE)
  }
  required_means <- if (pate) c("mean0", "mean1") else "mean0"
  if (!is.list(mean_models) || is.null(names(mean_models)) ||
      !setequal(names(mean_models), required_means) || anyDuplicated(names(mean_models))) {
    stop("Supply mean0 and, for PATE only, mean1 recorded models.", call. = FALSE)
  }
  if (family == "PS") {
    if (!is.null(raw_X) || !is.null(polynomial_exponents) || !is.list(ps) ||
        !all(c("design", "probability", "weighting") %in% names(ps)) ||
        any(!names(ps) %in% c("design", "probability", "weighting", "coefficients")) ||
        anyDuplicated(names(ps)) ||
        !is.character(ps$weighting) || length(ps$weighting) != 1L ||
        is.na(ps$weighting) || !ps$weighting %in% c("probability", "unit")) {
      stop("PS requires its saved design/probabilities and declared original fit weighting.", call. = FALSE)
    }
    Dps <- design_check(ps$design, n, "ps$design")
    probability <- vector_check(ps$probability, n, "ps$probability")
    if (any(probability <= 0 | probability >= 1)) {
      stop("Recorded PS probabilities must lie strictly inside (0,1).", call. = FALSE)
    }
    if (!is.null(ps$coefficients)) {
      alpha <- vector_check(ps$coefficients, ncol(Dps), "ps$coefficients")
      if (!is.null(names(ps$coefficients)) &&
          !identical(names(ps$coefficients), colnames(Dps))) {
        stop("PS coefficient names do not bind the original design order.", call. = FALSE)
      }
      agree(as.vector(stats::plogis(Dps %*% alpha)), probability,
            "Saved logistic coefficients and original probabilities")
    } else alpha <- rep(NA_real_, ncol(Dps))
    raw_coordinate <- matrix(probability, n, 1L, dimnames = list(NULL, "ps"))
    psw <- if (ps$weighting == "probability") W / max(W) else rep(1, n)
  } else {
    if (!is.null(ps)) stop("full_X has no logistic parameter block.", call. = FALSE)
    raw_coordinate <- design_check(raw_X, n, "raw_X")
    Dps <- NULL; alpha <- numeric(); psw <- NULL
  }
  coordinate_names <- colnames(raw_coordinate)
  d <- ncol(raw_coordinate)
  center <- colMeans(raw_coordinate)
  centered <- sweep(raw_coordinate, 2L, center, "-")
  variance <- colMeans(centered^2)
  if (any(!is.finite(variance)) || any(variance <= 0)) {
    stop("Recorded pooled coordinate variances are not positive.", call. = FALSE)
  }
  standardized <- sweep(centered, 2L, sqrt(variance), "/")
  if (ncol(fit$graph$scores0) != d ||
      (pate && ncol(fit$graph$scores1) != d)) {
    stop("Recorded family dimensions do not match the original graph.", call. = FALSE)
  }
  agree(as.numeric(standardized), as.numeric(fit$graph$scores0), "Original arm-0 standardized map")
  if (pate) agree(as.numeric(standardized), as.numeric(fit$graph$scores1),
                   "Original arm-1 standardized map")

  transform <- NULL; exponents <- NULL; raw_basis <- NULL; source_basis <- NULL
  if (family == "full_X") {
    if (is.null(polynomial_exponents)) {
      exponents <- rbind(rep(0, d), diag(d))
      basis_names <- c("intercept", coordinate_names)
      for (j in seq_len(d)) for (k in j:d) {
        a <- rep(0, d); a[j] <- a[j] + 1; a[k] <- a[k] + 1
        exponents <- rbind(exponents, a)
        basis_names <- c(basis_names, paste0(coordinate_names[j], ":", coordinate_names[k]))
      }
      dimnames(exponents) <- list(basis_names, coordinate_names)
    } else {
      exponents <- polynomial_exponents
      if (!is.matrix(exponents) || !is.numeric(exponents) || is.complex(exponents) ||
          ncol(exponents) != d || nrow(exponents) < 1L || any(!is.finite(exponents)) ||
          any(exponents < 0 | exponents != floor(exponents)) ||
          is.null(rownames(exponents)) || anyNA(rownames(exponents)) ||
          any(!nzchar(rownames(exponents))) || anyDuplicated(rownames(exponents)) ||
          !identical(colnames(exponents), coordinate_names)) {
        stop("Polynomial multi-indices must bind the recorded coordinate/basis order.", call. = FALSE)
      }
      degree <- max(rowSums(exponents))
      if (degree < 1 || nrow(unique(as.data.frame(exponents))) != nrow(exponents) ||
          !is.finite(choose(d + degree, degree)) ||
          nrow(exponents) != choose(d + degree, degree)) {
        stop("Use every distinct monomial through a fixed total polynomial degree.", call. = FALSE)
      }
    }
    raw_basis <- .wm_saved_polynomial_basis(raw_coordinate, exponents)
    source_basis <- .wm_saved_polynomial_basis(standardized, exponents)
    transform <- .wm_saved_polynomial_transform(exponents, center, variance)
    agree(as.numeric(raw_basis %*% transform), as.numeric(source_basis),
          "Complete polynomial affine identity")
  }

  mean_designs <- mean_parameters <- mean_predictions <- list()
  source_mean_coefficients <- list()
  for (label in required_means) {
    model <- mean_models[[label]]
    if (!is.list(model) || !all(c("design", "coefficients") %in% names(model))) {
      stop(label, " must contain its original design and coefficients.", call. = FALSE)
    }
    design <- design_check(model$design, n, paste0(label, "$design"))
    beta <- vector_check(model$coefficients, ncol(design), paste0(label, "$coefficients"))
    if (!is.null(names(model$coefficients)) &&
        !identical(names(model$coefficients), colnames(design))) {
      stop(label, " coefficient names do not bind its design order.", call. = FALSE)
    }
    source_mean_coefficients[[label]] <- beta
    recorded_mean <- fit$predictions[[label]]
    agree(as.vector(design %*% beta), as.numeric(recorded_mean),
          paste0(label, " original prediction binding"))
    if (family == "full_X") {
      if (!identical(colnames(design), colnames(source_basis))) {
        stop(label, " design columns do not bind the complete polynomial order.", call. = FALSE)
      }
      agree(as.numeric(design), as.numeric(source_basis), paste0(label, " standardized source basis"))
      beta <- as.vector(transform %*% beta)
      design <- raw_basis
      agree(as.vector(design %*% beta), as.numeric(recorded_mean),
            paste0(label, " raw polynomial prediction binding"))
    }
    mean_designs[[label]] <- design
    mean_parameters[[label]] <- beta
    mean_predictions[[label]] <- as.vector(design %*% beta)
  }

  parameter_names <- character(); blocks <- list()
  add_block <- function(prefix, labels) {
    index <- length(parameter_names) + seq_along(labels)
    parameter_names <<- c(parameter_names, paste0(prefix, ":", labels))
    index
  }
  if (family == "PS") blocks$ps <- add_block("ps", colnames(Dps))
  blocks$center <- add_block("center", coordinate_names)
  blocks$variance <- add_block("variance", coordinate_names)
  for (label in required_means) blocks[[label]] <- add_block(label, colnames(mean_designs[[label]]))
  p <- length(parameter_names)
  theta <- stats::setNames(numeric(p), parameter_names)
  if (family == "PS") theta[blocks$ps] <- alpha
  theta[blocks$center] <- center; theta[blocks$variance] <- variance
  equations <- matrix(0, n, p, dimnames = list(NULL, parameter_names))
  jacobian <- matrix(0, p, p, dimnames = list(parameter_names, parameter_names))
  score_derivative <- array(0, dim = c(n, d, p),
                            dimnames = list(NULL, coordinate_names, parameter_names))
  if (family == "PS") {
    derivative <- Dps * (probability * (1 - probability))
    equations[, blocks$ps] <- Dps * (psw * (Z - probability))
    jacobian[blocks$ps, blocks$ps] <- -crossprod(Dps, Dps * (psw * probability * (1 - probability))) / n
    jacobian[blocks$center, blocks$ps] <- colMeans(derivative)
    jacobian[blocks$variance, blocks$ps] <- 2 * colMeans(derivative * centered[, 1L])
    score_derivative[, 1L, blocks$ps] <- derivative / sqrt(variance[1L])
  }
  equations[, blocks$center] <- centered
  equations[, blocks$variance] <- sweep(centered^2, 2L, variance, "-")
  jacobian[blocks$center, blocks$center] <- -diag(d)
  jacobian[blocks$variance, blocks$variance] <- -diag(d)
  jacobian[blocks$variance, blocks$center] <- -2 * diag(colMeans(centered), nrow = d)
  for (j in seq_len(d)) {
    score_derivative[, j, blocks$center[j]] <- -1 / sqrt(variance[j])
    score_derivative[, j, blocks$variance[j]] <- -standardized[, j] / (2 * variance[j])
  }
  derivatives <- list(mean0 = NULL, mean1 = NULL)
  w <- W / max(W)
  for (label in required_means) {
    arm <- if (label == "mean0") 0L else 1L
    D <- mean_designs[[label]]; index <- blocks[[label]]
    theta[index] <- mean_parameters[[label]]
    iw <- w * as.numeric(Z == arm)
    equations[, index] <- D * (iw * (Y - mean_predictions[[label]]))
    jacobian[index, index] <- -crossprod(D, D * iw) / n
    derivatives[[label]] <- matrix(0, n, p, dimnames = list(NULL, parameter_names))
    derivatives[[label]][, index] <- D
  }
  if (any(!is.finite(c(equations, jacobian, score_derivative, unlist(derivatives))))) {
    stop("Complete stack evaluation exceeded numerical range.", call. = FALSE)
  }
  influence <- tryCatch(-t(solve(jacobian, t(equations))), error = function(e) {
    stop("Saved complete stack Jacobian is not invertible: ", conditionMessage(e), call. = FALSE)
  })
  colnames(influence) <- parameter_names
  if (any(!is.finite(influence))) stop("Complete influence calculation overflowed.", call. = FALSE)
  centered_influence <- sweep(influence, 2L, colMeans(influence), "-")
  Sigma <- crossprod(centered_influence) / n
  inference_arguments <- list(fit = fit, nuisance_influence = influence,
    mean_derivative0 = derivatives$mean0, mean_derivative1 = derivatives$mean1,
    covariance_scope = if (pate) "full_x" else "patt")
  if (!pate) inference_arguments$transport0 <- list(mode = "zero", basis = "current_centering")
  structure(list(family = family, n = n, estimand = fit$estimand, M = fit$M,
    parameter_names = parameter_names, parameter = theta,
    parameter_values_available = is.finite(theta), blocks = blocks,
    estimating_equations = equations, mean_equation = colMeans(equations),
    jacobian = jacobian, nuisance_influence = influence, nuisance_covariance = Sigma,
    mean0_derivative = derivatives$mean0, mean1_derivative = derivatives$mean1,
    score_derivative = score_derivative, coordinate_center = center,
    coordinate_variance = variance, polynomial_exponents = exponents,
    raw_polynomial_transform = transform, source_mean_coefficients = source_mean_coefficients,
    mean_models_raw = stats::setNames(lapply(required_means, function(label) list(
      design = mean_designs[[label]], coefficients = mean_parameters[[label]])), required_means),
    ps_weighting = if (family == "PS") ps$weighting else NULL,
    inference_arguments = inference_arguments,
    assumptions_verified = FALSE, application_verified = FALSE,
    numerical_root_verified = FALSE, original_refit_callback_verified = FALSE,
    contract = paste("Complete saved-stack arithmetic only; original point, graph, known weights",
      "and donor fractions are retained. No fit or graph is constructed.",
      "PS inputs declare the original zero-offset logistic design; finite probability",
      "bindings do not establish original design provenance or a statistical root.",
      "Full-X correctness, the original regular statistical root, numerical localization,",
      "iid/identification/geometry/moment and actual-row premises remain unverified.",
      "Original dependent/unbounded source designs are empirical approximations.",
      "No source full/partial callback, finite-B coverage, growing-B or variance-UI claim.")),
    class = c("wm_saved_full_x_stack", "list"))
}
