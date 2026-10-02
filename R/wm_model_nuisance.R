# Finite lists of logistic and arm-specific OLS models.
# Existing .wm_ns_* helpers and primitive solver controls are reused unchanged.
# Matching geometry, correction correctness and fitted-law premises are not
# inferred from successful roots or a model label.

.wm_model_specs <- function(models, n, kind, label) {
  if (is.null(models)) models <- list()
  if (!is.list(models) || is.data.frame(models)) {
    stop(label, " must be a list of model descriptors.", call. = FALSE)
  }
  labels <- names(models)
  if (length(models) && !is.null(labels) &&
      (anyNA(labels) || any(!nzchar(labels)) || anyDuplicated(labels))) {
    stop(label, " model names must be nonempty and unique.", call. = FALSE)
  }
  if (is.null(labels)) labels <- if (length(models))
    paste0("model", seq_along(models)) else character()
  allowed <- if (kind == "ps") c("design", "offset", "weighting") else
    c("design", "offset")
  result <- vector("list", length(models))
  names(result) <- labels
  for (j in seq_along(models)) {
    spec <- models[[j]]
    if (!is.list(spec) || is.data.frame(spec) || is.null(names(spec)) ||
        anyNA(names(spec)) || any(!nzchar(names(spec))) ||
        anyDuplicated(names(spec)) || any(!names(spec) %in% allowed) ||
        !"design" %in% names(spec)) {
      stop(label, "[[", j, "]] requires a design and only supported fields.",
           call. = FALSE)
    }
    design <- .wm_ns_design(spec$design, n, paste0(label, "[[", j, "]] design"))
    offset <- .wm_ns_vector(spec$offset, n,
                            paste0(label, "[[", j, "]] offset"), 0)
    entry <- list(design = design, offset = offset)
    if (kind == "ps") {
      weighting <- if (is.null(spec$weighting)) "probability" else spec$weighting
      if (!is.character(weighting) || length(weighting) != 1L ||
          is.na(weighting) || !weighting %in% c("probability", "unit")) {
        stop(label, "[[", j, "]] weighting must be 'probability' or 'unit'.",
             call. = FALSE)
      }
      entry$weighting <- weighting
    }
    result[[j]] <- entry
  }
  result
}

# No nominally identical models are merged. One PS descriptor is genuinely
# shared: it is fitted once and selected in every required donor-arm map.
.wm_model_layout <- function(ps, pg, used) {
  ps_names <- if (length(ps)) paste0("ps", seq_along(ps)) else character()
  pg_names <- lapply(used, function(z) {
    models <- pg[[as.character(z)]]
    if (length(models)) paste0("pg", z, "_", seq_along(models)) else character()
  })
  names(pg_names) <- as.character(used)
  raw_names <- c(ps_names, unlist(pg_names, use.names = FALSE))
  scores <- lapply(used, function(z) c(ps_names, pg_names[[as.character(z)]]))
  names(scores) <- as.character(used)
  if (any(lengths(scores) < 1L)) {
    stop("Every required arm needs at least one PS or arm-PG coordinate.",
         call. = FALSE)
  }
  list(ps_names = ps_names, pg_names = pg_names, raw_names = raw_names,
       scores = scores, matching_dimensions = lengths(scores),
       model_names = list(ps = names(ps), pg = lapply(pg, names)))
}

# Conservative element budget for the full retained stack, not a basis trim.
# Existing input matrices and R/solver overhead still require external resources.
.wm_model_allocation <- function(n, ps, pg, layout, max_stack_elements) {
  if (!is.numeric(max_stack_elements) || is.complex(max_stack_elements) ||
      !is.null(dim(max_stack_elements)) || length(max_stack_elements) != 1L ||
      !is.finite(max_stack_elements) || max_stack_elements < 1) {
    stop("max_stack_elements must be a positive finite allocation limit.", call. = FALSE)
  }
  if (n < 2L) stop("At least two rows are required.", call. = FALSE)
  dims <- as.double(layout$matching_dimensions)
  basis <- (dims + 1) * (dims + 2) / 2
  primitive <- sum(vapply(ps, function(x) as.double(ncol(x$design)), numeric(1L))) +
    sum(vapply(pg, function(models)
      sum(vapply(models, function(x) as.double(ncol(x$design)), numeric(1L))), numeric(1L)))
  L <- length(layout$raw_names)
  p <- primitive + 2 * L + sum(basis)
  # Raw/score derivatives, equations, prediction derivatives, generated bases
  # and their derivatives, and complete Jacobian/sandwich matrices are counted.
  # Include the retained zero weight derivative and complete influence rows,
  # including the single-coordinate case where raw derivative slack is smaller.
  planned <- n * (2 * L * p + 2 * L + 3 * p + length(dims) * p +
                   sum(basis * (dims + 2))) + 3 * p^2
  largest <- c(p, basis, n * p, p^2, n * basis)
  if (any(!is.finite(c(largest, planned))) || any(largest > .Machine$integer.max)) {
    stop("The complete nuisance/basis allocation is not representable.", call. = FALSE)
  }
  if (planned > max_stack_elements) {
    stop("Complete nuisance allocation exceeds max_stack_elements; increase the explicit limit only with adequate resources.",
         call. = FALSE)
  }
  list(parameter_dimension = as.integer(p), raw_coordinates = L,
    basis_columns = stats::setNames(as.integer(basis), names(layout$matching_dimensions)),
    planned_elements = planned, max_stack_elements = max_stack_elements,
    contract = "Conservative retained-array element guard; input matrices, copies and R/solver overhead are not an RSS bound.")
}

# Complete quadratic, lexicographic pairs a <= b, with no reduced subset.
.wm_model_pairs <- function(d) {
  q <- (as.double(d) + 1) * (as.double(d) + 2) / 2
  if (d < 1L || !is.finite(q) || q > .Machine$integer.max) {
    stop("The complete quadratic basis is not representable.", call. = FALSE)
  }
  a <- rep(seq_len(d), times = d:1L)
  b <- unlist(lapply(seq_len(d), function(j) seq.int(j, d)), use.names = FALSE)
  cbind(a = a, b = b)
}

.wm_model_basis <- function(s) {
  d <- ncol(s)
  pairs <- .wm_model_pairs(d)
  nm <- colnames(s)
  if (is.null(nm)) nm <- paste0("score", seq_len(d))
  b <- cbind(intercept = 1, s,
             s[, pairs[, 1L], drop = FALSE] * s[, pairs[, 2L], drop = FALSE])
  colnames(b) <- c("intercept", nm,
    paste0("quadratic:", nm[pairs[, 1L]], ":", nm[pairs[, 2L]]))
  b
}

.wm_model_basis_derivatives <- function(s) {
  n <- nrow(s)
  d <- ncol(s)
  pairs <- .wm_model_pairs(d)
  q <- 1L + d + nrow(pairs)
  output <- vector("list", d)
  names(output) <- colnames(s)
  for (h in seq_len(d)) {
    db <- matrix(0, n, q)
    db[, 1L + h] <- 1
    for (j in seq_len(nrow(pairs))) {
      a <- pairs[j, 1L]
      b <- pairs[j, 2L]
      if (a == h) db[, 1L + d + j] <- db[, 1L + d + j] + s[, b]
      if (b == h) db[, 1L + d + j] <- db[, 1L + d + j] + s[, a]
    }
    output[[h]] <- db
  }
  output
}

# At arbitrary admissible parameters, retain every generated-coordinate block.
.wm_model_nuisance_evaluate <- function(object, parameter = object$parameter) {
  if (inherits(object, ".wm_wdsm_nuisance_stack")) {
    return(.wm_wdsm_nuisance_evaluate(object, parameter))
  }
  if (!inherits(object, ".wm_model_nuisance_stack") || !is.list(object)) {
    stop("Supply a finite candidate-model nuisance stack.", call. = FALSE)
  }
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
  layout <- object$layout
  pi <- object$empirical_probability
  w <- object$scaled_weights
  raw <- matrix(0, n, length(layout$raw_names),
                dimnames = list(NULL, layout$raw_names))
  draw <- ds <- vector("list", ncol(raw))
  names(draw) <- names(ds) <- layout$raw_names
  equation <- matrix(0, n, p, dimnames = list(NULL, object$parameter_names))
  jacobian <- matrix(0, p, p,
    dimnames = list(object$parameter_names, object$parameter_names))
  for (j in seq_along(d$ps_models)) {
    label <- layout$ps_names[j]
    spec <- d$ps_models[[j]]
    index <- k$ps[[j]]
    psw <- object$scaled_ps_weights[[j]]
    value <- stats::plogis(as.vector(spec$design %*% theta[index]) + spec$offset)
    if (any(value <= 0 | value >= 1)) {
      stop("Propensity parameter gives saturated probabilities.")
    }
    raw[, label] <- value
    draw[[label]] <- matrix(0, n, p, dimnames = list(NULL, object$parameter_names))
    draw[[label]][, index] <- spec$design * (value * (1 - value))
    equation[, index] <- spec$design * (psw * (d$Z - value))
    jacobian[index, index] <- -crossprod(spec$design,
      spec$design * (pi * psw * value * (1 - value)))
  }
  for (z in object$used_arms) {
    arm <- as.character(z)
    for (j in seq_along(d$pg_models[[arm]])) {
      label <- layout$pg_names[[arm]][j]
      spec <- d$pg_models[[arm]][[j]]
      index <- k$pg[[arm]][[j]]
      value <- as.vector(spec$design %*% theta[index]) + spec$offset
      raw[, label] <- value
      draw[[label]] <- matrix(0, n, p,
                            dimnames = list(NULL, object$parameter_names))
      draw[[label]][, index] <- spec$design
      iz <- as.numeric(d$Z == z)
      equation[, index] <- spec$design * (iz * (d$Y - value))
      jacobian[index, index] <- -crossprod(spec$design, spec$design * (pi * iz))
    }
  }
  center <- theta[k$center]
  variance <- theta[k$variance]
  if (any(variance <= 0)) stop("All pooled variance parameters must be positive.")
  centered <- sweep(raw, 2L, center, "-")
  standardized <- sweep(centered, 2L, sqrt(variance), "/")
  equation[, k$center] <- centered
  equation[, k$variance] <- sweep(centered^2, 2L, variance, "-")
  for (j in seq_len(ncol(raw))) {
    label <- layout$raw_names[j]
    ds[[label]] <- draw[[label]] / sqrt(variance[j])
    ds[[label]][, k$center[j]] <- ds[[label]][, k$center[j]] - 1 / sqrt(variance[j])
    ds[[label]][, k$variance[j]] <- ds[[label]][, k$variance[j]] -
      standardized[, j] / (2 * variance[j])
    jacobian[k$center[j], ] <- colSums(draw[[label]] * pi)
    jacobian[k$center[j], k$center[j]] <- -1
    jacobian[k$variance[j], ] <- colSums(draw[[label]] * (2 * pi * centered[, j]))
    jacobian[k$variance[j], k$center[j]] <- -2 * sum(pi * centered[, j])
    jacobian[k$variance[j], k$variance[j]] <- -1
  }
  scores <- means <- mean_derivative <- score_derivatives <- basis <- list()
  for (z in object$used_arms) {
    arm <- as.character(z)
    index <- k$correction[[arm]]
    labels <- layout$scores[[arm]]
    s <- standardized[, labels, drop = FALSE]
    b <- .wm_model_basis(s)
    db <- .wm_model_basis_derivatives(s)
    beta <- theta[index]
    value <- as.vector(b %*% beta)
    dq <- matrix(0, n, p, dimnames = list(NULL, object$parameter_names))
    iw <- as.numeric(d$Z == z) * w
    residual <- d$Y - value
    for (h in seq_along(labels)) {
      tangent <- ds[[labels[h]]]
      dq <- dq + tangent * as.vector(db[[h]] %*% beta)
      jacobian[index, ] <- jacobian[index, , drop = FALSE] +
        crossprod(db[[h]] * (pi * iw * residual), tangent)
    }
    dq[, index] <- dq[, index, drop = FALSE] + b
    equation[, index] <- b * (iw * residual)
    jacobian[index, ] <- jacobian[index, , drop = FALSE] -
      crossprod(b * (pi * iw), dq)
    scores[[arm]] <- s
    means[[arm]] <- value
    mean_derivative[[arm]] <- dq
    score_derivatives[[arm]] <- ds[labels]
    basis[[arm]] <- b
  }
  if (any(!is.finite(c(raw, standardized, equation, jacobian,
                        unlist(mean_derivative), unlist(score_derivatives))))) {
    stop("Nuisance evaluation exceeded numerical range.", call. = FALSE)
  }
  list(parameter = theta, raw_scores = raw, center = center, variance = variance,
    scores0 = scores[["0"]], scores1 = scores[["1"]],
    mean0 = means[["0"]], mean1 = means[["1"]],
    mean0_derivative = mean_derivative[["0"]],
    mean1_derivative = mean_derivative[["1"]],
    score_derivatives = score_derivatives, basis = basis,
    estimating_equations = equation,
    mean_equation = as.vector(crossprod(pi, equation)), jacobian = jacobian)
}

.wm_model_stack <- function(Y, Z, weights, ps_models = list(),
                             pg0_models = list(), pg1_models = NULL,
                             estimand = c("PATE", "PATT"), multiplicity = NULL,
                             glm_maxit = 100L, glm_epsilon = 1e-12,
                             rank_tolerance = 1e-10, max_stack_elements = 5e7,
                             influence = TRUE) {
  estimand <- match.arg(estimand)
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
  solver <- .wm_wdsm_controls(list(glm_maxit = glm_maxit,
    glm_epsilon = glm_epsilon, refine_maxit = glm_maxit,
    refine_epsilon = glm_epsilon, rank_tolerance = rank_tolerance))
  if (!is.logical(influence) || length(influence) != 1L || is.na(influence)) {
    stop("influence must be TRUE or FALSE.")
  }
  used <- if (estimand == "PATE") 0:1 else 0L
  if (estimand == "PATT" && !is.null(pg1_models)) {
    stop("PATT does not fit an unused treated prognostic model.")
  }
  ps <- .wm_model_specs(ps_models, n, "ps", "ps_models")
  pg <- list("0" = .wm_model_specs(pg0_models, n, "pg", "pg0_models"))
  if (estimand == "PATE") {
    pg[["1"]] <- .wm_model_specs(pg1_models, n, "pg", "pg1_models")
  }
  layout <- .wm_model_layout(ps, pg, used)
  allocation <- .wm_model_allocation(n, ps, pg, layout, max_stack_elements)
  # Exact existing double-score recipe uses existing constructor/evaluator.
  legacy <- length(ps) == 1L && all(lengths(pg) == 1L)
  if (legacy) {
    fitter <- if (influence) .wm_wdsm_nuisance_stack else .wm_wdsm_prediction_stack
    return(fitter(Y, Z, weights, ps[[1L]]$design, pg[["0"]][[1L]]$design,
      pg1_design = if (estimand == "PATE") pg[["1"]][[1L]]$design else NULL,
      estimand = estimand, ps_offset = ps[[1L]]$offset,
      pg0_offset = pg[["0"]][[1L]]$offset,
      pg1_offset = if (estimand == "PATE") pg[["1"]][[1L]]$offset else NULL,
      multiplicity = m, glm_maxit = solver$glm_maxit,
      glm_epsilon = solver$glm_epsilon, rank_tolerance = solver$rank_tolerance,
      ps_weighting = ps[[1L]]$weighting))
  }
  d <- list(Y = Y, Z = Z, weights = weights, ps_models = ps, pg_models = pg)
  prob <- m / sum(m)
  count_weight <- n * prob
  weight_scale <- max(weights)
  w <- weights / weight_scale
  if (any(!is.finite(w)) || any(w <= 0)) {
    stop("Common weight scaling underflowed; supply representable weight ratios.")
  }
  psw <- lapply(ps, function(spec)
    if (spec$weighting == "probability") w else rep(1, n))
  parameter_names <- character()
  add_block <- function(prefix, labels) {
    if (length(parameter_names) + length(labels) > .Machine$integer.max) {
      stop("The complete nuisance dimension is not representable.")
    }
    index <- length(parameter_names) + seq_along(labels)
    parameter_names <<- c(parameter_names, paste0(prefix, ":", labels))
    index
  }
  k <- list(ps = list(), pg = list(), correction = list())
  for (j in seq_along(ps)) {
    k$ps[[j]] <- add_block(layout$ps_names[j], colnames(ps[[j]]$design))
  }
  for (z in used) {
    arm <- as.character(z)
    k$pg[[arm]] <- list()
    for (j in seq_along(pg[[arm]])) {
      k$pg[[arm]][[j]] <- add_block(layout$pg_names[[arm]][j],
                                   colnames(pg[[arm]][[j]]$design))
    }
  }
  k$center <- stats::setNames(add_block("center", layout$raw_names), layout$raw_names)
  k$variance <- stats::setNames(add_block("variance", layout$raw_names), layout$raw_names)
  for (z in used) {
    arm <- as.character(z)
    dummy <- matrix(0, 1L, length(layout$scores[[arm]]),
                    dimnames = list(NULL, layout$scores[[arm]]))
    k$correction[[arm]] <- add_block(paste0("bc", z), colnames(.wm_model_basis(dummy)))
  }
  if (anyDuplicated(parameter_names)) stop("Complete parameter labels must be unique.")
  theta <- stats::setNames(numeric(length(parameter_names)), parameter_names)
  raw <- matrix(0, n, length(layout$raw_names),
                dimnames = list(NULL, layout$raw_names))
  iterations <- stats::setNames(integer(length(ps)), names(ps))
  for (j in seq_along(ps)) {
    spec <- ps[[j]]
    fitted <- stats::glm.fit(x = spec$design, y = Z,
      weights = count_weight * psw[[j]], offset = spec$offset,
      start = rep(0, ncol(spec$design)), family = stats::quasibinomial(link = "logit"),
      control = stats::glm.control(epsilon = solver$glm_epsilon,
                                    maxit = as.integer(solver$glm_maxit)))
    if (!isTRUE(fitted$converged) || fitted$rank != ncol(spec$design) ||
        any(!is.finite(fitted$coefficients))) {
      stop("Weighted propensity fit failed to select a finite full-rank root.")
    }
    theta[k$ps[[j]]] <- fitted$coefficients
    raw[, layout$ps_names[j]] <- stats::plogis(
      as.vector(spec$design %*% fitted$coefficients) + spec$offset)
    iterations[j] <- fitted$iter
  }
  if (length(ps) && any(raw[, layout$ps_names, drop = FALSE] <= .Machine$double.eps |
      raw[, layout$ps_names, drop = FALSE] >= 1 - .Machine$double.eps)) {
    stop("Propensity fit has numerically saturated probabilities.")
  }
  for (z in used) {
    arm <- as.character(z)
    for (j in seq_along(pg[[arm]])) {
      spec <- pg[[arm]][[j]]
      index <- k$pg[[arm]][[j]]
      theta[index] <- .wm_ns_ols(spec$design, Y - spec$offset,
        count_weight * as.numeric(Z == z), layout$pg_names[[arm]][j], solver$rank_tolerance)
      raw[, layout$pg_names[[arm]][j]] <- as.vector(spec$design %*% theta[index]) + spec$offset
    }
  }
  theta[k$center] <- colSums(raw * prob)
  centered <- sweep(raw, 2L, theta[k$center], "-")
  theta[k$variance] <- colSums(centered^2 * prob)
  if (any(!is.finite(theta[k$variance])) || any(theta[k$variance] <= 0)) {
    stop("A pooled matching-score variance is zero or nonfinite.")
  }
  standardized <- sweep(centered, 2L, sqrt(theta[k$variance]), "/")
  for (z in used) {
    arm <- as.character(z)
    b <- .wm_model_basis(standardized[, layout$scores[[arm]], drop = FALSE])
    theta[k$correction[[arm]]] <- .wm_ns_ols(b, Y,
      count_weight * w * as.numeric(Z == z), paste0("bc", z), solver$rank_tolerance)
  }
  object <- structure(list(parameter = theta, parameter_names = parameter_names,
    blocks = k, used_arms = used, estimand = estimand, inputs = d, layout = layout,
    multiplicity = m, empirical_probability = prob, weight_scale = weight_scale,
    scaled_weights = w, ps_weighting = vapply(ps, function(x) x$weighting, character(1L)),
    scaled_ps_weights = psw), class = c(".wm_model_nuisance_stack", "list"))
  evaluated <- .wm_model_nuisance_evaluate(object)
  A <- evaluated$jacobian
  if (qr(A, tol = solver$rank_tolerance)$rank != ncol(A)) {
    stop("The full stacked estimating-equation Jacobian is numerically singular.")
  }
  normalized <- stats::setNames(numeric(length(ps)), names(ps))
  for (j in seq_along(ps)) {
    scale <- colSums(abs(ps[[j]]$design) * (prob * psw[[j]]))
    normalized[j] <- max(abs(evaluated$mean_equation[k$ps[[j]]]) /
      pmax(scale, .Machine$double.eps * sum(prob * psw[[j]])))
  }
  normalized_max <- if (length(normalized)) max(normalized) else 0
  if (!is.finite(normalized_max) || normalized_max > 1e-7) {
    stop("Weighted propensity root fails the normalized score check.")
  }
  for (name in names(evaluated)) object[[name]] <- evaluated[[name]]
  object$n <- n
  object$allocation <- allocation
  object$weight_derivative <- matrix(0, n, length(theta),
                                    dimnames = list(NULL, parameter_names))
  if (influence) {
    complete <- -t(solve(A, t(evaluated$estimating_equations)))
    colnames(complete) <- parameter_names
    influence_mean <- as.vector(crossprod(prob, complete))
    centered_influence <- sweep(complete, 2L, influence_mean, "-")
    Sigma <- crossprod(centered_influence, centered_influence * prob)
    dimnames(Sigma) <- list(parameter_names, parameter_names)
    if (any(!is.finite(c(complete, Sigma)))) stop("Stacked sandwich calculation overflowed.")
    object$nuisance_influence <- complete
    object$nuisance_covariance <- Sigma
  }
  object$diagnostics <- list(glm_iterations = iterations,
    normalized_ps_score = normalized_max, normalized_ps_scores = normalized,
    max_mean_equation = max(abs(evaluated$mean_equation)),
    jacobian_reciprocal_condition = rcond(A), glm_epsilon = solver$glm_epsilon,
    rank_tolerance = solver$rank_tolerance)
  if (influence) object$diagnostics$max_influence_mean <- max(abs(influence_mean))
  object$assumptions_verified <- FALSE
  object$contract <- paste("Known supplied weights; each declared logistic PS fitted once",
    "and shared across arms, unweighted arm OLS, pooled empirical centers/variances,",
    "weighted complete quadratic correction. All coefficients and generated-coordinate",
    "cross blocks retained.", if (influence) "Empirical influence is -A^{-1}psi; centered empirical covariance is not a verification of population assumptions." else
      "Prediction-only construction omits nuisance influence and covariance.", "For nonunit",
    "multiplicities this is empirical-measure arithmetic, not iid sampling inference.")
  object
}

.wm_model_nuisance_stack <- function(...) .wm_model_stack(..., influence = TRUE)
.wm_model_prediction_stack <- function(...) .wm_model_stack(..., influence = FALSE)
