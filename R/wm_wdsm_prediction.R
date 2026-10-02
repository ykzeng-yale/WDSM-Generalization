# Source-derived package implementation; the original reviewed overlay is preserved privately.
# Extends the default PS equations with an explicit unit-fit-weight option.
# Audited source: wdsm_calibration_prediction_stack.R; SHA256 f547ce5a1fc207f5bc861894082927a50ecd377b4e2612bc50a493d369d1ae74.
# Count-refit constructor copied from the pinned reference constructor:
# 88ef7f11079605e50c395964347e2ad2d78a82c50ac5da5f012a0c7fdce4e65c.
# Omits the unused nuisance IF/covariance and its diagnostics.
# Uses the package nuisance helpers/evaluator, including declared PS fit weights.
# Newton/root checks are performed by .wm_wdsm_fit_stack before predictions are used.

.wm_wdsm_prediction_stack <- function(Y, Z, weights, ps_design,
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
  score_scale <- colSums(abs(d$ps_design) * (prob * psw))
  ps_equation <- evaluated$mean_equation[k$ps]
  normalized_ps_residual <- max(abs(ps_equation) /
    pmax(score_scale, .Machine$double.eps * sum(prob * psw)))
  if (!is.finite(normalized_ps_residual) || normalized_ps_residual > 1e-7) {
    stop("Weighted propensity root fails the normalized score check.")
  }
  for (name in names(evaluated)) object[[name]] <- evaluated[[name]]
  object$n <- n
  object$weight_derivative <- matrix(0, n, length(theta),
                                     dimnames = list(NULL, parameter_names))
  object$diagnostics <- list(glm_iterations = psfit$iter,
    normalized_ps_score = normalized_ps_residual,
    max_mean_equation = max(abs(evaluated$mean_equation)),
    jacobian_reciprocal_condition = rcond(A),
    glm_epsilon = glm_epsilon, rank_tolerance = rank_tolerance)
  object$assumptions_verified <- FALSE
  object$contract <- paste("Count-refit predictions and complete estimating equations/Jacobian only.",
    "No nuisance influence matrix or sampling sandwich is calculated.",
    "Shared helpers are from pinned wdsm_nuisance_stack_reference.R.")
  object
}
