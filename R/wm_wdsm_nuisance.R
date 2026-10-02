# Source-derived package implementation; the original reviewed overlay is preserved privately.
# Extends the default PS equations with an explicit unit-fit-weight option.
# Audited source: wdsm_nuisance_stack_reference.R; SHA256 88ef7f11079605e50c395964347e2ad2d78a82c50ac5da5f012a0c7fdce4e65c.
# Standalone source-faithful WDSM nuisance stack. Not a package export.
# All inference columns use the single parameter order returned by this file.

.wm_ns_vector <- function(x, n, label, default = NULL) {
  if (is.null(x) && !is.null(default)) x <- rep(default, n)
  if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) ||
      length(x) != n || any(!is.finite(x))) {
    stop(label, " must be a finite numeric vector of length ", n, ".", call. = FALSE)
  }
  as.numeric(x)
}

.wm_ns_design <- function(x, n, label) {
  if (!is.matrix(x) || !is.numeric(x) || is.complex(x) ||
      nrow(x) != n || ncol(x) < 1L || any(!is.finite(x))) {
    stop(label, " must be a finite numeric matrix with n rows.", call. = FALSE)
  }
  nm <- colnames(x)
  if (is.null(nm)) nm <- paste0("x", seq_len(ncol(x)))
  if (anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm)) {
    stop(label, " column names must be nonempty and unique.", call. = FALSE)
  }
  storage.mode(x) <- "double"
  colnames(x) <- nm
  x
}

.wm_ns_basis <- function(s) {
  cbind(intercept = 1, ps = s[, 1L], pg = s[, 2L],
        ps_sq = s[, 1L]^2, ps_pg = s[, 1L] * s[, 2L],
        pg_sq = s[, 2L]^2)
}

.wm_ns_basis_derivatives <- function(s) {
  list(ps = cbind(0, 1, 0, 2 * s[, 1L], s[, 2L], 0),
       pg = cbind(0, 0, 1, 0, s[, 1L], 2 * s[, 2L]))
}

.wm_ns_ols <- function(x, y, w, label, tolerance) {
  if (sum(w > 0) < ncol(x)) stop(label, " has too few positive-weight rows.")
  fit <- stats::lm.wfit(x = x, y = y, w = w, tol = tolerance,
                       singular.ok = FALSE)
  if (fit$rank != ncol(x) || any(!is.finite(fit$coefficients))) {
    stop(label, " has unidentified or nonfinite coefficients.", call. = FALSE)
  }
  as.numeric(fit$coefficients)
}

# Evaluate the complete equation, Jacobian and prediction/score derivatives
# at an arbitrary admissible stack parameter, with data and counts fixed.
.wm_wdsm_nuisance_evaluate <- function(object,
                                               parameter = object$parameter) {
  if (!inherits(object, ".wm_wdsm_nuisance_stack") ||
      !is.list(object)) stop("Supply a nuisance-stack reference object.")
  d <- object$inputs
  n <- length(d$Y)
  p <- length(object$parameter_names)
  if (!is.null(names(parameter)) &&
      !identical(names(parameter), object$parameter_names)) {
    stop("Parameter names must use the stored order.", call. = FALSE)
  }
  theta <- .wm_ns_vector(parameter, p, "parameter")
  names(theta) <- object$parameter_names
  k <- object$blocks
  used <- object$used_arms
  pi <- object$empirical_probability
  w <- object$scaled_weights
  # Older saved stack objects used probability weights for PS fitting.
  psw <- if (is.null(object$scaled_ps_weights)) w else object$scaled_ps_weights
  J <- length(k$center)
  raw_names <- names(k$center)
  raw <- matrix(0, n, J, dimnames = list(NULL, raw_names))
  draw <- vector("list", J)
  names(draw) <- raw_names
  ps <- stats::plogis(as.vector(d$ps_design %*% theta[k$ps]) + d$ps_offset)
  if (any(ps <= 0 | ps >= 1)) stop("Propensity parameter gives saturated probabilities.")
  raw[, "ps"] <- ps
  draw$ps <- matrix(0, n, p, dimnames = list(NULL, object$parameter_names))
  draw$ps[, k$ps] <- d$ps_design * (ps * (1 - ps))
  for (z in used) {
    label <- paste0("pg", z)
    raw[, label] <- as.vector(d[[paste0(label, "_design")]] %*%
                               theta[k[[label]]]) + d[[paste0(label, "_offset")]]
    draw[[label]] <- matrix(0, n, p, dimnames = list(NULL, object$parameter_names))
    draw[[label]][, k[[label]]] <- d[[paste0(label, "_design")]]
  }
  center <- theta[k$center]
  variance <- theta[k$variance]
  if (any(variance <= 0)) stop("All pooled variance parameters must be positive.")
  centered <- sweep(raw, 2L, center, "-")
  standardized <- sweep(centered, 2L, sqrt(variance), "/")
  ds <- vector("list", J)
  names(ds) <- raw_names
  for (j in seq_len(J)) {
    ds[[j]] <- draw[[j]] / sqrt(variance[j])
    ds[[j]][, k$center[j]] <- ds[[j]][, k$center[j]] - 1 / sqrt(variance[j])
    ds[[j]][, k$variance[j]] <- ds[[j]][, k$variance[j]] -
      standardized[, j] / (2 * variance[j])
  }

  equation <- matrix(0, n, p, dimnames = list(NULL, object$parameter_names))
  jacobian <- matrix(0, p, p,
                     dimnames = list(object$parameter_names, object$parameter_names))
  equation[, k$ps] <- d$ps_design * (psw * (d$Z - ps))
  jacobian[k$ps, k$ps] <- -crossprod(d$ps_design,
                                    d$ps_design * (pi * psw * ps * (1 - ps)))
  for (z in used) {
    label <- paste0("pg", z)
    x <- d[[paste0(label, "_design")]]
    iz <- as.numeric(d$Z == z)
    equation[, k[[label]]] <- x * (iz * (d$Y - raw[, label]))
    jacobian[k[[label]], k[[label]]] <- -crossprod(x, x * (pi * iz))
  }
  equation[, k$center] <- centered
  equation[, k$variance] <- sweep(centered^2, 2L, variance, "-")
  for (j in seq_len(J)) {
    jacobian[k$center[j], ] <- colSums(draw[[j]] * pi)
    jacobian[k$center[j], k$center[j]] <- -1
    jacobian[k$variance[j], ] <- colSums(draw[[j]] * (2 * pi * centered[, j]))
    jacobian[k$variance[j], k$center[j]] <- -2 * sum(pi * centered[, j])
    jacobian[k$variance[j], k$variance[j]] <- -1
  }

  scores <- means <- mean_derivative <- score_derivative <- basis <- list()
  for (z in used) {
    arm <- as.character(z)
    label <- paste0("pg", z)
    bc <- k[[paste0("bc", z)]]
    s <- standardized[, c("ps", label), drop = FALSE]
    colnames(s) <- c("ps", "pg")
    b <- .wm_ns_basis(s)
    db <- .wm_ns_basis_derivatives(s)
    beta <- theta[bc]
    q <- as.vector(b %*% beta)
    dq <- ds$ps * as.vector(db$ps %*% beta) +
      ds[[label]] * as.vector(db$pg %*% beta)
    dq[, bc] <- dq[, bc, drop = FALSE] + b
    iw <- as.numeric(d$Z == z) * w
    residual <- d$Y - q
    equation[, bc] <- b * (iw * residual)
    # Differentiate both the generated basis and the generated residual.
    jacobian[bc, ] <-
      crossprod(db$ps * (pi * iw * residual), ds$ps) +
      crossprod(db$pg * (pi * iw * residual), ds[[label]]) -
      crossprod(b * (pi * iw), dq)
    scores[[arm]] <- s
    means[[arm]] <- q
    mean_derivative[[arm]] <- dq
    score_derivative[[arm]] <- list(ps = ds$ps, pg = ds[[label]])
    basis[[arm]] <- b
  }
  if (any(!is.finite(c(raw, standardized, equation, jacobian,
                       unlist(mean_derivative), unlist(score_derivative))))) {
    stop("Nuisance evaluation exceeded numerical range.", call. = FALSE)
  }
  list(parameter = theta, raw_scores = raw, center = center, variance = variance,
       scores0 = scores[["0"]], scores1 = scores[["1"]],
       mean0 = means[["0"]], mean1 = means[["1"]],
       mean0_derivative = mean_derivative[["0"]],
       mean1_derivative = mean_derivative[["1"]],
       score_derivatives = score_derivative, basis = basis,
       estimating_equations = equation,
       mean_equation = as.vector(crossprod(pi, equation)), jacobian = jacobian)
}

# Input design matrices must include their intended intercept columns.
# Multiplicities may be nonnegative real numbers for deterministic reweighting;
# they are not asserted to have a multinomial distribution.
.wm_wdsm_nuisance_stack <- function(Y, Z, weights, ps_design,
                                             pg0_design, pg1_design = NULL,
                                             estimand = c("PATE", "PATT"),
                                             ps_offset = NULL,
                                             pg0_offset = NULL, pg1_offset = NULL,
                                             multiplicity = NULL,
                                             glm_maxit = 100L,
                                             glm_epsilon = 1e-12,
                                             rank_tolerance = 1e-10,
                                             ps_weighting = c("probability", "unit")) {
  estimand <- match.arg(estimand)
  ps_weighting <- match.arg(ps_weighting)
  n <- length(Y)
  if (n < 2L) stop("At least two rows are required.")
  Y <- .wm_ns_vector(Y, n, "Y")
  if (is.logical(Z) && is.null(dim(Z))) Z <- as.numeric(Z)
  Z <- .wm_ns_vector(Z, n, "Z")
  weights <- .wm_ns_vector(weights, n, "weights")
  if (!all(Z %in% c(0, 1)) || any(weights <= 0)) {
    stop("Z must be binary and weights strictly positive.")
  }
  m <- .wm_ns_vector(multiplicity, n, "multiplicity", 1)
  if (any(m < 0) || !is.finite(sum(m)) || sum(m) <= 0 ||
      any(vapply(0:1, function(z) sum(m[Z == z]) <= 0, logical(1L)))) {
    stop("Multiplicities need a finite positive total and positive mass in both arms.")
  }
  if (!is.numeric(glm_maxit) || is.complex(glm_maxit) || length(glm_maxit) != 1L ||
      !is.finite(glm_maxit) || glm_maxit < 1 || glm_maxit != floor(glm_maxit) ||
      glm_maxit > .Machine$integer.max ||
      !is.numeric(glm_epsilon) || is.complex(glm_epsilon) || length(glm_epsilon) != 1L ||
      !is.finite(glm_epsilon) || glm_epsilon <= 0 || glm_epsilon >= 1 ||
      !is.numeric(rank_tolerance) || is.complex(rank_tolerance) || length(rank_tolerance) != 1L ||
      !is.finite(rank_tolerance) || rank_tolerance <= 0 || rank_tolerance >= 1) {
    stop("Invalid solver controls.")
  }
  used <- if (estimand == "PATE") 0:1 else 0L
  if (estimand == "PATT" && (!is.null(pg1_design) || !is.null(pg1_offset))) {
    stop("PATT does not fit an unused treated prognostic model.")
  }
  d <- list(Y = Y, Z = Z, weights = weights,
            ps_design = .wm_ns_design(ps_design, n, "ps_design"),
            pg0_design = .wm_ns_design(pg0_design, n, "pg0_design"),
            ps_offset = .wm_ns_vector(ps_offset, n, "ps_offset", 0),
            pg0_offset = .wm_ns_vector(pg0_offset, n, "pg0_offset", 0))
  if (estimand == "PATE") {
    d$pg1_design <- .wm_ns_design(pg1_design, n, "pg1_design")
    d$pg1_offset <- .wm_ns_vector(pg1_offset, n, "pg1_offset", 0)
  }
  prob <- m / sum(m)
  count_weight <- n * prob
  weight_scale <- max(weights)
  w <- weights / weight_scale
  if (any(!is.finite(w)) || any(w <= 0)) {
    stop("Common weight scaling underflowed; supply representable weight ratios.")
  }
  # PS fitting weights are separate from the effect and correction weights.
  psw <- if (ps_weighting == "probability") w else rep(1, n)
  parameter_names <- character()
  k <- list()
  add_block <- function(prefix, labels) {
    index <- length(parameter_names) + seq_along(labels)
    parameter_names <<- c(parameter_names, paste0(prefix, ":", labels))
    index
  }
  k$ps <- add_block("ps", colnames(d$ps_design))
  for (z in used) {
    label <- paste0("pg", z)
    k[[label]] <- add_block(label, colnames(d[[paste0(label, "_design")]]))
  }
  raw_names <- c("ps", paste0("pg", used))
  k$center <- stats::setNames(add_block("center", raw_names), raw_names)
  k$variance <- stats::setNames(add_block("variance", raw_names), raw_names)
  for (z in used) {
    k[[paste0("bc", z)]] <- add_block(paste0("bc", z),
      c("intercept", "ps", "pg", "ps_sq", "ps_pg", "pg_sq"))
  }
  theta <- stats::setNames(numeric(length(parameter_names)), parameter_names)
  psfit <- stats::glm.fit(x = d$ps_design, y = Z,
    weights = count_weight * psw, offset = d$ps_offset,
    start = rep(0, ncol(d$ps_design)),
    family = stats::quasibinomial(link = "logit"),
    control = stats::glm.control(epsilon = glm_epsilon, maxit = as.integer(glm_maxit)))
  if (!isTRUE(psfit$converged) || psfit$rank != ncol(d$ps_design) ||
      any(!is.finite(psfit$coefficients))) {
    stop("Weighted propensity fit failed to select a finite full-rank root.")
  }
  theta[k$ps] <- psfit$coefficients
  raw <- matrix(0, n, length(raw_names), dimnames = list(NULL, raw_names))
  raw[, "ps"] <- stats::plogis(as.vector(d$ps_design %*% theta[k$ps]) + d$ps_offset)
  if (any(raw[, "ps"] <= .Machine$double.eps |
          raw[, "ps"] >= 1 - .Machine$double.eps)) {
    stop("Propensity fit has numerically saturated probabilities.")
  }
  for (z in used) {
    label <- paste0("pg", z)
    x <- d[[paste0(label, "_design")]]
    offset <- d[[paste0(label, "_offset")]]
    theta[k[[label]]] <- .wm_ns_ols(x, Y - offset,
      count_weight * as.numeric(Z == z), label, rank_tolerance)
    raw[, label] <- as.vector(x %*% theta[k[[label]]]) + offset
  }
  theta[k$center] <- colSums(raw * prob)
  centered <- sweep(raw, 2L, theta[k$center], "-")
  theta[k$variance] <- colSums(centered^2 * prob)
  if (any(!is.finite(theta[k$variance])) || any(theta[k$variance] <= 0)) {
    stop("A pooled matching-score variance is zero or nonfinite.")
  }
  standardized <- sweep(centered, 2L, sqrt(theta[k$variance]), "/")
  for (z in used) {
    b <- .wm_ns_basis(standardized[, c("ps", paste0("pg", z)), drop = FALSE])
    theta[k[[paste0("bc", z)]]] <- .wm_ns_ols(b, Y,
      count_weight * w * as.numeric(Z == z), paste0("bc", z), rank_tolerance)
  }
  object <- structure(list(parameter = theta, parameter_names = parameter_names,
    blocks = k, used_arms = used, estimand = estimand, inputs = d,
    multiplicity = m, empirical_probability = prob,
    weight_scale = weight_scale, scaled_weights = w,
    ps_weighting = ps_weighting, scaled_ps_weights = psw),
    class = c(".wm_wdsm_nuisance_stack", "list"))
  evaluated <- .wm_wdsm_nuisance_evaluate(object)
  A <- evaluated$jacobian
  if (qr(A, tol = rank_tolerance)$rank != ncol(A)) {
    stop("The full stacked estimating-equation Jacobian is numerically singular.")
  }
  influence <- -t(solve(A, t(evaluated$estimating_equations)))
  colnames(influence) <- parameter_names
  influence_mean <- as.vector(crossprod(prob, influence))
  influence_centered <- sweep(influence, 2L, influence_mean, "-")
  Sigma <- crossprod(influence_centered, influence_centered * prob)
  dimnames(Sigma) <- list(parameter_names, parameter_names)
  if (any(!is.finite(c(influence, Sigma)))) stop("Stacked sandwich calculation overflowed.")
  score_scale <- colSums(abs(d$ps_design) * (prob * psw))
  ps_equation <- evaluated$mean_equation[k$ps]
  normalized_ps_residual <- max(abs(ps_equation) /
    pmax(score_scale, .Machine$double.eps * sum(prob * psw)))
  if (!is.finite(normalized_ps_residual) || normalized_ps_residual > 1e-7) {
    stop("Weighted propensity root fails the normalized score check.")
  }
  for (name in names(evaluated)) object[[name]] <- evaluated[[name]]
  object$n <- n
  object$nuisance_influence <- influence
  object$nuisance_covariance <- Sigma
  object$weight_derivative <- matrix(0, n, length(theta),
                                     dimnames = list(NULL, parameter_names))
  object$diagnostics <- list(glm_iterations = psfit$iter,
    normalized_ps_score = normalized_ps_residual,
    max_mean_equation = max(abs(evaluated$mean_equation)),
    max_influence_mean = max(abs(influence_mean)),
    jacobian_reciprocal_condition = rcond(A),
    glm_epsilon = glm_epsilon, rank_tolerance = rank_tolerance)
  object$assumptions_verified <- FALSE
  object$contract <- paste("Known probability weights. Logistic PS with declared fitting weights, unweighted",
    "arm OLS, pooled count-weighted centers/variances, weighted quadratic BC.",
    "Influence is -J^{-1}psi with every generated-score cross block retained.",
    "For nonunit multiplicities this is the influence of the weighted empirical",
    "measure, not an iid-row sampling variance assertion. Model centering, geometry,",
    "identification, regular roots and graph transport require separate verification.")
  object
}
