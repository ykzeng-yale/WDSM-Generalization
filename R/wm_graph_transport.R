# Private additive package candidate; not installed or included in formal sources.
# Exact numerical reference body below, with two identifier-only substitutions.
# Reference: simulations/wm_graph_transport_reference.R
# SHA256: ff62fd5a5c5868acc68f75821e0aa95e7aeec60aa822bf50dc891204199f2af3

# Standalone reference for the reviewed incoming-load graph transport.
# Requires explicit raw chart tangents and a justified interior cutoff.
# This does not verify the statistical assumptions or fit a nuisance model.

.wm_gt_finite <- function(x, label) {
  if (!is.numeric(x) || is.complex(x) || any(!is.finite(x))) {
    stop(label, " must be finite and numeric.", call. = FALSE)
  }
  invisible(x)
}

# Conditional empirical tagged-weight expectation and its spatial derivative.
# Simpson error bounds below exclude floating-point roundoff.
.wm_gt_tagged_load <- function(weight, probability, probability_gradient,
                                     root_weight, M, tolerance = 1e-7,
                                     max_nodes = 131073L) {
  .wm_gt_finite(weight, "weight")
  .wm_gt_finite(probability, "probability")
  .wm_gt_finite(probability_gradient, "probability_gradient")
  .wm_gt_finite(root_weight, "root_weight")
  .wm_gt_finite(M, "M")
  if (!is.null(dim(weight)) || !is.null(dim(probability)) ||
      !length(weight) || length(weight) != length(probability) ||
      any(weight <= 0) || any(probability < 0) ||
      abs(sum(probability) - 1) > 1e-10 ||
      !is.matrix(probability_gradient) ||
      nrow(probability_gradient) != length(weight) ||
      ncol(probability_gradient) < 1L || length(root_weight) != 1L ||
      root_weight <= 0 || length(M) != 1L || M < 1 || M != floor(M)) {
    stop("Invalid positive weight law, spatial derivative, root weight, or M.")
  }
  gradient_size <- colSums(abs(probability_gradient))
  if (any(abs(colSums(probability_gradient)) > 1e-10 * pmax(1, gradient_size))) {
    stop("Each probability-gradient column must sum to zero.")
  }
  if (any(probability == 0 & rowSums(abs(probability_gradient)) != 0)) {
    stop("A zero-probability interior atom must have zero spatial derivative.")
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance <= 0 ||
      !is.numeric(max_nodes) || length(max_nodes) != 1L ||
      !is.finite(max_nodes) || max_nodes < 3 || max_nodes != floor(max_nodes)) {
    stop("tolerance and max_nodes must be positive finite controls.")
  }
  d <- ncol(probability_gradient)
  if (M == 1) {
    return(list(value = 1, gradient = numeric(d), error_bound = rep(0, d + 1L),
                nodes = 0L, method = "exact M=1", roundoff_included = FALSE))
  }
  positive <- probability > 0
  weight <- weight[positive]
  probability <- probability[positive]
  probability_gradient <- probability_gradient[positive, , drop = FALSE]
  # Project only the already accepted roundoff-sized simplex discrepancies.
  mass <- sum(probability)
  probability <- probability / mass
  probability_gradient <- (probability_gradient -
    outer(probability, colSums(probability_gradient))) / mass
  gradient_size <- colSums(abs(probability_gradient))
  weight_scale <- max(root_weight, weight)
  root_weight <- root_weight / weight_scale
  weight <- weight / weight_scale
  if (root_weight <= 0 || any(weight <= 0)) stop("Weight dynamic range underflow.")
  if (all(weight == weight[1L])) {
    # The conditional law is a point mass; its derivative integrates to zero.
    exact <- root_weight / (root_weight + (M - 1) * weight[1L])
    if (!is.finite((M - 1) * weight[1L])) stop("M exceeds numerical limits.")
    return(list(value = exact,
                gradient = numeric(d), error_bound = rep(0, d + 1L),
                nodes = 0L, method = "exact point mass", roundoff_included = FALSE))
  }
  m <- M - 1
  a <- root_weight + m * min(weight)
  largest <- root_weight + m * max(weight)
  coefficient <- (root_weight / a) * c(1, m * gradient_size)
  derivative_bound <- coefficient * (largest / a)^4
  if (any(!is.finite(c(a, largest, coefficient, derivative_bound)))) {
    stop("Weight or derivative range exceeds numerical limits.")
  }
  upper <- max(1, log(2) + log(max(coefficient)) - log(tolerance))
  panels <- 1
  repeat {
    step <- upper / (2 * panels)
    # Sum of the fourth-derivative bounds over Simpson panels.
    simpson_error <- derivative_bound * step^5 /
      (90 * (-expm1(-2 * step))) * (-expm1(-upper))
    tail_error <- coefficient * exp(-upper)
    error_bound <- simpson_error + tail_error
    if (max(error_bound) <= tolerance) break
    panels <- 2 * panels
    if (2 * panels + 1 > max_nodes) {
      stop("Quadrature node budget cannot meet the analytic error bound.")
    }
  }
  count <- 2 * panels + 1
  integral <- numeric(d + 1L)
  # Bound temporary matrices independently of the number of quadrature nodes.
  chunk <- max(1L, min(count, floor(2e6 / length(weight))))
  for (start in seq.int(1, count, by = chunk)) {
    index <- seq.int(start, min(count, start + chunk - 1))
    u <- (index - 1) * step
    laplace_matrix <- exp(-outer((weight - min(weight)) / a, u))
    local_laplace <- as.vector(crossprod(probability, laplace_matrix))
    local_gradient <- crossprod(probability_gradient, laplace_matrix)
    multiplier <- (root_weight / a) * exp(-u)
    integrand <- rbind(multiplier * local_laplace^m,
      sweep(m * local_gradient, 2L,
            multiplier * local_laplace^(m - 1), `*`))
    simpson_weight <- ifelse(index == 1 | index == count, 1,
      ifelse((index - 1) %% 2 == 1, 4, 2))
    integral <- integral + as.vector(integrand %*% simpson_weight) * step / 3
  }
  if (any(!is.finite(integral))) stop("Nonfinite Laplace quadrature result.")
  list(value = integral[1L], gradient = integral[-1L], error_bound = error_bound,
       nodes = count, upper = upper, step = step,
       method = "scaled Laplace integral, composite Simpson with analytic bound",
       roundoff_included = FALSE)
}

.wm_gt_kernel <- function(coordinates, evaluation, bandwidth) {
  u <- sweep(coordinates, 2L, evaluation, function(x, y) (y - x) / bandwidth)
  inside <- abs(u) < 1
  k <- (35 / 32) * pmax(0, 1 - u^2)^3 * inside
  dk <- -(210 / 32) * u * pmax(0, 1 - u^2)^2 * inside / bandwidth
  n <- nrow(u)
  d <- ncol(u)
  value <- apply(k, 1L, prod)
  gradient <- matrix(0, n, d)
  for (j in seq_len(d)) {
    other <- if (d == 1L) rep(1, n) else apply(k[, -j, drop = FALSE], 1L, prod)
    gradient[, j] <- dk[, j] * other
  }
  list(value = value, gradient = gradient)
}

# Returns unnormalized j_z in ORIGINAL raw-weight units.
# tangents[i,,k] is the effective raw-chart tangent v_k(X_i), not merely
# the derivative of a standardized or transformed matching coordinate.
.wm_gt_compute <- function(raw_scores, Z, W, residual,
                                         tangents, cutoff, donor_arm, M,
                                         bandwidth = NULL, density_floor = NULL,
                                         quadrature_tolerance = 1e-7,
                                         max_nodes = 131073L) {
  for (name in c("raw_scores", "Z", "W", "residual", "tangents", "cutoff")) {
    .wm_gt_finite(get(name), name)
  }
  if (!is.matrix(raw_scores) || nrow(raw_scores) < 2L || ncol(raw_scores) < 1L) {
    stop("raw_scores must be an n-by-d matrix.")
  }
  n <- nrow(raw_scores)
  d <- ncol(raw_scores)
  if (length(Z) != n || !is.null(dim(Z)) || any(!Z %in% 0:1) ||
      length(W) != n || !is.null(dim(W)) || any(W <= 0) ||
      length(residual) != n || !is.null(dim(residual)) ||
      length(cutoff) != n || !is.null(dim(cutoff)) || any(cutoff < 0 | cutoff > 1) ||
      length(dim(tangents)) != 3L ||
      !identical(dim(tangents)[1:2], c(n, d)) || dim(tangents)[3] < 1L ||
      !is.numeric(donor_arm) || length(donor_arm) != 1L ||
      !is.finite(donor_arm) || !donor_arm %in% 0:1 ||
      !is.numeric(M) || length(M) != 1L || !is.finite(M) || M < 1 || M != floor(M)) {
    stop("Invalid design, weight, residual, tangent, cutoff, arm or M input.")
  }
  donor <- which(Z == donor_arm)
  query <- which(Z != donor_arm)
  if (length(donor) < M || !length(query)) stop("Insufficient arm sizes.")
  if (is.null(bandwidth)) bandwidth <- n^(-1 / (2 * (d + 4)))
  if (is.null(density_floor)) density_floor <- 0.01 * n^(-1/4)
  for (name in c("bandwidth", "density_floor", "quadrature_tolerance")) {
    control <- get(name)
    if (!is.numeric(control) || length(control) != 1L || !is.finite(control) || control <= 0) {
      stop(name, " must be a positive finite scalar.")
    }
  }
  p <- dim(tangents)[3]
  active <- donor[cutoff[donor] > 0]
  loads <- numeric(length(active))
  gradients <- errors <- matrix(0, length(active), d)
  density <- numeric(length(active))
  nodes <- numeric(length(active))
  j <- j_error <- numeric(p)
  for (position in seq_along(active)) {
    i <- active[position]
    kernel <- .wm_gt_kernel(raw_scores, raw_scores[i, ], bandwidth)
    kd <- kernel$value[donor]
    gd <- kernel$gradient[donor, , drop = FALSE]
    total <- sum(kd)
    density[position] <- total / (n * bandwidth^d)
    if (!is.finite(density[position]) || density[position] <= density_floor) {
      stop("Donor density is below the declared floor at row ", i, ".")
    }
    probability <- kd / total
    probability_gradient <- (gd - outer(probability, colSums(gd))) / total
    qmass <- sum(W[query] * kernel$value[query])
    ratio <- qmass / total
    ratio_gradient <- (colSums(kernel$gradient[query, , drop = FALSE] * W[query]) -
      ratio * colSums(gd)) / total
    tagged <- .wm_gt_tagged_load(W[donor], probability, probability_gradient,
      W[i], M, tolerance = quadrature_tolerance, max_nodes = max_nodes)
    loads[position] <- M * ratio * tagged$value
    gradients[position, ] <- M * (ratio_gradient * tagged$value + ratio * tagged$gradient)
    errors[position, ] <- M * (abs(ratio_gradient) * tagged$error_bound[1L] +
      abs(ratio) * tagged$error_bound[-1L])
    nodes[position] <- tagged$nodes
    tangent <- matrix(tangents[i, , ], d, p)
    j <- j + residual[i] * cutoff[i] * as.vector(crossprod(tangent, gradients[position, ])) / n
    j_error <- j_error + abs(residual[i] * cutoff[i]) *
      as.vector(crossprod(abs(tangent), errors[position, ])) / n
  }
  if (any(!is.finite(c(loads, gradients, j, j_error)))) stop("Nonfinite graph transport result.")
  names(j) <- names(j_error) <- dimnames(tangents)[[3L]]
  list(graph_drift = j, quadrature_error_bound = j_error,
       evaluated_donor_rows = active, incoming_load = loads,
       incoming_load_gradient = gradients, gradient_error_bound = errors,
       donor_subdensity = density, quadrature_nodes = nodes,
       n = n, d = d, M = M, donor_arm = donor_arm, bandwidth = bandwidth,
       density_floor = density_floor, weight_units = "original raw W",
       assumptions_verified = FALSE, roundoff_included = FALSE,
       contract = paste("Unnormalized graph drift j_z from the specified raw chart.",
         "Requires own-(score,weight) centering, interior signed support and marked-density smoothness.",
         "Cutoff must equal one on that signed support; it is not arbitrary trimming.",
         "Sampling error, bias and floating-point roundoff are outside the quadrature bound.",
         "Use signed arm contrast divided by the target-weight mean; add smooth sensitivity separately."))
}

#' Estimate the incoming graph-transport component in a supplied raw chart
#' @export
wm_graph_transport <- function(raw_scores, Z, W, residual, tangents, cutoff,
                               donor_arm, M, bandwidth = NULL,
                               density_floor = NULL,
                               quadrature_tolerance = 1e-7,
                               max_nodes = 131073L) {
  # Parameter names identify one common nuisance parameter order, not score axes.
  # The numerical core below separately validates all dimensions and values.
  parameter_names <- if (length(dim(tangents)) == 3L) dimnames(tangents)[[3L]] else NULL
  if (!is.null(parameter_names) &&
      (anyNA(parameter_names) || any(!nzchar(parameter_names)) ||
       anyDuplicated(parameter_names))) {
    stop("Tangent parameter names must be nonempty and unique.", call. = FALSE)
  }
  result <- .wm_gt_compute(raw_scores, Z, W, residual, tangents, cutoff,
    donor_arm, M, bandwidth, density_floor, quadrature_tolerance, max_nodes)
  result$parameter_names <- parameter_names
  result$parameter_order <- if (is.null(parameter_names)) "positional" else "named"
  result$status <- if (length(result$evaluated_donor_rows))
    "component_computed" else "no_active_donors"
  result$numerically_available <- TRUE
  result$sampling_inference_available <- FALSE
  result$bound_scope <- paste("Deterministic quadrature error for the supplied",
    "finite empirical kernel law in exact arithmetic; excludes sampling error,",
    "kernel bias, nuisance error and floating-point roundoff.")
  result$normalization <- paste("Unnormalized j_z in original raw-W units.",
    "PATE graph sensitivity is (j_1-j_0)/mean(W);",
    "PATT is -j_0/mean(Z*W). Add the complete fixed-edge smooth ratio",
    "derivative exactly once. No extra M, n or target normalization is included.")
  result$scope <- paste("This component calculation allows any fixed positive",
    "matching dimension and finite M. Applying same-sample fitted-map inference",
    "requires a separately justified point theorem and covariance branch;",
    "the general fitted-map companion theorem has used dimensions at least two.",
    "Scalar fitted inference requires its separate model-specific theorem.",
    "No variance, interval, model fit or bootstrap is produced here.")
  class(result) <- c("wm_graph_transport", "list")
  result
}
