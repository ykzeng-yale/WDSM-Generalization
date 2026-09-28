# Estimated-weight uncertainty with a fixed realized matching graph.

.wm_weight_matrix <- function(x, n, name, p = NULL) {
  if (!is.matrix(x) || !is.numeric(x) || is.complex(x) ||
      nrow(x) != n || ncol(x) < 1L || any(!is.finite(x)) ||
      (!is.null(p) && ncol(x) != p)) {
    stop(name, " must be a finite numeric n-by-p matrix",
         if (!is.null(p)) paste0(" with p = ", p) else "",
         ".", call. = FALSE)
  }
  x
}

.wm_weight_agrees <- function(x, y, tolerance = 1e-8) {
  length(x) == length(y) && all(is.finite(x)) && all(is.finite(y)) &&
    all(abs(x - y) <= tolerance * (1 + pmax(abs(x), abs(y))))
}

#' Adjust fixed-graph matching inference for estimated individual weights
#' @export
wm_weight_adjust <- function(fit, weight_derivative, nuisance_influence) {
  if (!inherits(fit, "wm_match") || !identical(fit$method, "self_normalized")) {
    stop("wm_weight_adjust requires a self_normalized wm_match or wm_fit object.",
         call. = FALSE)
  }
  if (!isTRUE(fit$info$corrected) || length(fit$root_n_variance) != 1L ||
      !is.finite(fit$root_n_variance) || fit$root_n_variance < 0) {
    stop("Weight adjustment requires variance = TRUE and fitted mean predictions.",
         call. = FALSE)
  }
  if (isTRUE(fit$info$nuisance_correction) ||
      !is.null(fit$weight_adjustment) ||
      any(fit$contributions$nuisance != 0)) {
    stop("The fit already contains a nuisance correction; start from the ",
         "uncorrected fit to avoid double counting.", call. = FALSE)
  }
  n <- .wm_fit_integer(fit$n, "fit$n", 2L)
  if (is.null(fit$data) || is.null(fit$data$Y) || is.null(fit$data$Z)) {
    stop("The fit must retain clean original data in fit$data$Y and fit$data$Z; ",
         "refit with the current matching core.", call. = FALSE)
  }
  Y <- .wm_numeric_vector(fit$data$Y, n, "fit$data$Y")
  Z <- fit$data$Z
  if ((!is.numeric(Z) && !is.logical(Z)) || is.complex(Z) ||
      !is.null(dim(Z)) || length(Z) != n || anyNA(Z) ||
      any(!Z %in% c(0, 1)) || length(unique(Z)) != 2L) {
    stop("Retained fit$data$Z is not a valid binary treatment vector.",
         call. = FALSE)
  }
  Z <- as.integer(Z)
  W <- .wm_numeric_vector(fit$weights, n, "fit$weights", positive = TRUE)
  dW <- .wm_weight_matrix(weight_derivative, n, "weight_derivative")
  p <- ncol(dW)
  influence <- .wm_weight_matrix(nuisance_influence, n, "nuisance_influence", p)
  if (!is.null(colnames(dW)) && !is.null(colnames(influence)) &&
      !identical(colnames(dW), colnames(influence))) {
    stop("Derivative and influence column names must have the same parameter order.",
         call. = FALSE)
  }
  parameter_names <- if (!is.null(colnames(dW))) colnames(dW) else
    if (!is.null(colnames(influence))) colnames(influence) else
      paste0("parameter", seq_len(p))
  mu0 <- .wm_numeric_vector(fit$predictions$mean0, n, "fit mean0")
  pate <- identical(fit$estimand, "PATE")
  if (!pate && !identical(fit$estimand, "PATT")) {
    stop("Unsupported fit estimand.", call. = FALSE)
  }
  mu1 <- if (pate) .wm_numeric_vector(fit$predictions$mean1, n, "fit mean1") else
    numeric(n)
  M <- .wm_fit_integer(fit$M, "fit$M")
  edges <- fit$graph$edges
  queries <- if (pate) seq_len(n) else which(Z == 1L)
  if (!is.data.frame(edges) || !all(c("query", "donor") %in% names(edges)) ||
      nrow(edges) != length(queries) * M) {
    stop("The retained matching graph is incomplete or malformed.", call. = FALSE)
  }
  for (name in c("query", "donor")) {
    index <- edges[[name]]
    if (!is.numeric(index) || is.complex(index) || any(!is.finite(index)) ||
        any(index != floor(index)) || any(index < 1L | index > n)) {
      stop("The retained matching graph has invalid row indices.", call. = FALSE)
    }
  }
  if (any(!edges$query %in% queries) || any(Z[edges$query] == Z[edges$donor]) ||
      any(tabulate(edges$query, nbins = n)[queries] != M) ||
      anyDuplicated(edges[c("query", "donor")])) {
    stop("The retained graph must contain exactly M distinct opposite-arm ",
         "donors per query.", call. = FALSE)
  }
  scale <- max(W)
  w <- W / scale
  dw <- dW / scale
  if (any(w <= 0) || any(!is.finite(dw))) {
    stop("Weight or derivative normalization exceeds numerical range.",
         call. = FALSE)
  }
  # The scale is held constant while differentiating; its derivative would
  # cancel exactly from the normalized estimator.
  outer_weight <- if (pate) w else Z * w
  denominator <- sum(outer_weight)
  gamma <- denominator / n
  load <- numeric(n)
  dload <- matrix(0, nrow = n, ncol = p)
  contrast <- numeric(n)
  neighbor_derivative <- numeric(p)
  fraction_derivative_error <- 0
  for (i in queries) {
    donors <- edges$donor[edges$query == i]
    donor_sum <- sum(w[donors])
    lambda <- w[donors] / donor_sum
    dsum <- colSums(dw[donors, , drop = FALSE])
    dlambda <- (dw[donors, , drop = FALSE] - outer(lambda, dsum)) / donor_sum
    mu <- if (Z[i] == 1L) mu0 else mu1
    residual <- Y[donors] - mu[donors]
    missing <- mu[i] + sum(lambda * residual)
    sign <- 1L - 2L * Z[i]
    contrast[i] <- sign * (missing - Y[i])
    neighbor_derivative <- neighbor_derivative +
      w[i] * sign * as.vector(crossprod(residual, dlambda))
    load[donors] <- load[donors] + w[i] * lambda
    dload[donors, ] <- dload[donors, , drop = FALSE] +
      outer(lambda, dw[i, ]) + w[i] * dlambda
    fraction_derivative_error <- max(fraction_derivative_error,
                                      abs(colSums(dlambda)))
  }
  reconstructed <- sum(outer_weight * contrast) / denominator
  if (length(fit$estimate) != 1L ||
      !.wm_weight_agrees(reconstructed, fit$estimate)) {
    stop("Retained data, weights, graph or means are inconsistent with the ",
         "stored point estimate; refit before adjustment.", call. = FALSE)
  }
  tau <- fit$estimate
  if (pate) {
    own_mean <- ifelse(Z == 1L, mu1, mu0)
    epsilon <- Y - own_mean
    sign <- 2L * Z - 1L
    sensitivity <- colSums(dw * (mu1 - mu0 - tau) +
                              (dw + dload) * (sign * epsilon)) / denominator
    base <- (w * (mu1 - mu0 - tau) + sign * (w + load) * epsilon) / gamma
  } else {
    epsilon <- Y - mu0
    sensitivity <- colSums(dw * (Z * (epsilon - tau)) -
                              dload * ((1L - Z) * epsilon)) / denominator
    base <- (Z * w * (epsilon - tau) - (1L - Z) * load * epsilon) / gamma
  }
  direct_sensitivity <- (colSums(dw * (if (pate) 1 else Z) *
                                    (contrast - tau)) +
                           neighbor_derivative) / denominator
  if (!.wm_weight_agrees(sensitivity, direct_sensitivity) ||
      !.wm_weight_agrees(base, fit$contributions$row)) {
    stop("The retained contribution rows or finite-graph derivative identity ",
         "are inconsistent; refit before adjustment.", call. = FALSE)
  }
  base <- base - mean(base)
  component <- as.vector(influence %*% sensitivity)
  component <- component - mean(component)
  corrected <- base + component
  corrected <- corrected - mean(corrected)
  root_variance <- mean(corrected^2)
  if (any(!is.finite(c(sensitivity, component, corrected, root_variance)))) {
    stop("Weight-adjustment arithmetic exceeds numerical range.", call. = FALSE)
  }
  names(sensitivity) <- names(direct_sensitivity) <- parameter_names
  fit$weight_adjustment <- list(
    sensitivity = sensitivity, direct_sensitivity = direct_sensitivity,
    weight_derivative = dW, nuisance_influence = influence,
    incoming_derivative = dload, normalization_scale = scale,
    max_fraction_derivative_sum = fraction_derivative_error,
    baseline_root_n_variance = mean(base^2),
    cross_term = 2 * mean(base * component),
    nuisance_root_n_variance = mean(component^2),
    parameter_names = parameter_names,
    contract = paste(
      "Fixed known graph and means held fixed in differentiation.",
      "The caller must establish the supplied influence expansion, row alignment,",
      "smooth bounded weight path and regression-transfer/centering assumptions.",
      "This correction covers weight estimation only, not estimated coordinates."),
    incoming_derivative_units = "Derivatives use raw weights divided by normalization_scale")
  fit$contributions$row <- corrected
  fit$contributions$nuisance <- component
  fit$root_n_variance <- root_variance
  fit$variance <- root_variance / n
  fit$se <- sqrt(fit$variance)
  fit$info$nuisance_correction <- TRUE
  fit$info$inference_status <- "fixed_graph_weight_influence_contract"
  fit
}

.wm_logistic_state <- function(parameter, design, Z, qr_tol) {
  n <- nrow(design)
  eta <- as.vector(design %*% parameter)
  if (any(!is.finite(eta))) {
    stop("Logistic linear predictors exceed numerical range.", call. = FALSE)
  }
  probability <- stats::plogis(eta)
  if (any(probability <= 0 | probability >= 1)) {
    stop("Logistic probabilities reached a numerical boundary; separation or ",
         "inadequate positivity prevents regular weight inference.",
         call. = FALSE)
  }
  variance <- probability * (1 - probability)
  weighted_design <- design * sqrt(variance)
  if (qr(weighted_design, tol = qr_tol, LAPACK = FALSE)$rank < ncol(design)) {
    stop("Logistic information is rank deficient; separation or a singular ",
         "design prevents regular influence estimation.", call. = FALSE)
  }
  information <- crossprod(weighted_design) / n
  factor <- tryCatch(chol(information), error = function(e) NULL)
  if (is.null(factor) || any(!is.finite(factor))) {
    stop("Logistic information is not numerically positive definite.",
         call. = FALSE)
  }
  score <- as.vector(crossprod(design, Z - probability)) / n
  step <- as.vector(backsolve(factor, forwardsolve(t(factor), score)))
  list(probability = probability, eta = eta, information = information,
       information_factor = factor, score = score, step = step,
       objective = mean(pmax(eta, 0) - Z * eta + log1p(exp(-abs(eta)))))
}

#' Fit bounded-domain logistic target weights and their influence vectors
#' @export
wm_logistic_weights <- function(Z, design, pi = 0.5, parameter_bound,
                                maxit = 100L, tolerance = 1e-8,
                                qr_tol = 1e-10) {
  if (missing(parameter_bound)) {
    stop("Supply a declared positive parameter_bound for the compact logistic ",
         "parameter domain.", call. = FALSE)
  }
  if (!is.matrix(design) || !is.numeric(design) || is.complex(design) ||
      nrow(design) < 2L || ncol(design) < 1L || any(!is.finite(design))) {
    stop("design must be a finite numeric n-by-p matrix, including any chosen ",
         "intercept column.", call. = FALSE)
  }
  n <- nrow(design)
  p <- ncol(design)
  if ((!is.numeric(Z) && !is.logical(Z)) || is.complex(Z) ||
      !is.null(dim(Z)) || length(Z) != n || anyNA(Z) ||
      any(!Z %in% c(0, 1)) || length(unique(Z)) != 2L) {
    stop("Z must contain both treatment arms and one binary value per design row.",
         call. = FALSE)
  }
  Z <- as.integer(Z)
  if (!is.numeric(pi) || is.complex(pi) || length(pi) != 1L ||
      !is.finite(pi) || pi <= 0 || pi >= 1) {
    stop("pi must be a known number strictly between zero and one.", call. = FALSE)
  }
  if (!is.numeric(parameter_bound) || is.complex(parameter_bound) ||
      !(length(parameter_bound) %in% c(1L, p)) ||
      any(!is.finite(parameter_bound)) || any(parameter_bound <= 0)) {
    stop("parameter_bound must be positive and finite, with length one or p.",
         call. = FALSE)
  }
  bounds <- rep(parameter_bound, length.out = p)
  maxit <- .wm_fit_integer(maxit, "maxit")
  for (name in c("tolerance", "qr_tol")) {
    value <- get(name)
    if (!is.numeric(value) || is.complex(value) || length(value) != 1L ||
        !is.finite(value) || value <= 0 || value >= 1) {
      stop(name, " must be finite and strictly between zero and one.",
           call. = FALSE)
    }
  }
  if (p > n || qr(design, tol = qr_tol, LAPACK = FALSE)$rank < p) {
    stop("The logistic design must have full column rank; no columns are dropped.",
         call. = FALSE)
  }
  parameter_names <- colnames(design)
  if (is.null(parameter_names)) parameter_names <- paste0("parameter", seq_len(p))
  if (anyNA(parameter_names) || any(!nzchar(parameter_names)) ||
      anyDuplicated(parameter_names)) {
    stop("Logistic design column names must be nonempty and unique.",
         call. = FALSE)
  }
  parameter <- numeric(p)
  converged <- FALSE
  boundary_tolerance <- min(1e-7, sqrt(tolerance))
  state <- NULL
  for (iteration in seq_len(maxit)) {
    if (any((bounds - abs(parameter)) / bounds <= boundary_tolerance)) {
      stop("The logistic fit reached the declared parameter boundary; ",
           "noninterior or separated fits do not have this regular influence ",
           "estimate.", call. = FALSE)
    }
    state <- .wm_logistic_state(parameter, design, Z, qr_tol)
    predictor_step <- as.vector(design %*% state$step)
    if (any(!is.finite(c(state$score, state$step, predictor_step,
                         state$objective)))) {
      stop("Logistic fitting arithmetic exceeds numerical range.", call. = FALSE)
    }
    # A separated fit can have a tiny score and information while retaining
    # an order-one Newton predictor step. Both convergence checks are needed.
    if (max(abs(state$score)) <= tolerance &&
        max(abs(predictor_step)) <= tolerance * (1 + max(abs(state$eta)))) {
      converged <- TRUE
      break
    }
    positive <- state$step > 0
    negative <- state$step < 0
    box_step <- min(c(Inf, (bounds[positive] - parameter[positive]) /
                       state$step[positive],
                     (-bounds[negative] - parameter[negative]) /
                       state$step[negative]))
    fraction <- min(1, 0.995 * box_step)
    descent <- sum(state$score * state$step)
    if (!is.finite(fraction) || fraction <= 0 || !is.finite(descent) ||
        descent <= 0) {
      stop("Logistic Newton step failed; no regular interior solution was found.",
           call. = FALSE)
    }
    accepted <- FALSE
    for (line_search in seq_len(60L)) {
      candidate <- parameter + fraction * state$step
      eta <- as.vector(design %*% candidate)
      objective <- if (all(is.finite(eta)))
        mean(pmax(eta, 0) - Z * eta + log1p(exp(-abs(eta)))) else Inf
      if (is.finite(objective) && all(abs(candidate) < bounds) &&
          objective <= state$objective - 1e-4 * fraction * descent +
            10 * .Machine$double.eps * (1 + abs(state$objective))) {
        parameter <- candidate
        accepted <- TRUE
        break
      }
      fraction <- fraction / 2
    }
    if (!accepted) {
      stop("Logistic line search did not converge to a regular interior fit.",
           call. = FALSE)
    }
  }
  if (!converged) {
    stop("Logistic fit did not converge within maxit; separation, an inadequate ",
         "parameter domain, or numerical conditioning may be responsible.",
         call. = FALSE)
  }
  q <- state$probability
  weights <- ifelse(Z == 1L, pi / q, (1 - pi) / (1 - q))
  derivative <- design * (weights * (q - Z))
  score_rows <- design * (Z - q)
  influence <- t(backsolve(state$information_factor,
                           forwardsolve(t(state$information_factor),
                                        t(score_rows))))
  if (any(!is.finite(c(weights, derivative, influence)))) {
    stop("Logistic weights, derivatives or influence values exceed numerical range.",
         call. = FALSE)
  }
  names(parameter) <- names(state$score) <- parameter_names
  colnames(derivative) <- colnames(influence) <- parameter_names
  dimnames(state$information) <- list(parameter_names, parameter_names)
  structure(list(
    weights = weights, weight_derivative = derivative,
    nuisance_influence = influence, parameters = parameter,
    propensity = q, empirical_information = state$information,
    score = state$score, iterations = iteration, converged = TRUE,
    pi = pi, parameter_bound = stats::setNames(bounds, parameter_names),
    design = design, Z = Z, tolerance = tolerance, qr_tol = qr_tol,
    assumptions = list(
      model = "Correct Q(Z=1|X)=expit(T(X)'alpha), with the supplied fixed design map",
      design = "Population design vector is bounded; observed finite values alone do not verify this",
      parameter = "The true parameter is strictly inside the declared compact symmetric box",
      information = "Population information is positive definite",
      weights = "W1=pi/q and W0=(1-pi)/(1-q), with fixed known pi",
      influence = "Same-sample iid regular logistic MLE influence; no clipping, ridge or Firth correction",
      matching = "Identification, graph-field centering and fixed-map regression-transfer assumptions remain separate")),
    class = c("wm_logistic_weights", "list"))
}
