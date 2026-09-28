# Universal Euclidean overlap constants. All general-dimensional results below
# are numerical estimates of the exact finite angular integral.

.wm_log_ball_union <- function(radius1, radius2, separation, d) {
  n <- length(radius1)
  if (length(radius2) != n || length(separation) != n ||
      any(!is.finite(c(radius1, radius2, separation))) ||
      any(c(radius1, radius2, separation) < 0)) {
    stop("Ball radii and separations must be finite nonnegative vectors of equal length.",
         call. = FALSE)
  }
  log_unit_ball <- d * log(pi) / 2 - lgamma(d / 2 + 1)
  answer <- rep(-Inf, n)
  positive <- pmax(radius1, radius2) > 0
  if (!any(positive)) return(answer)
  scale <- pmax(radius1[positive], radius2[positive])
  r <- radius1[positive] / scale
  s <- radius2[positive] / scale
  delta <- separation[positive] / scale
  r_power <- r^d
  s_power <- s^d
  union_scaled <- pmax(r_power, s_power)
  disjoint <- delta >= r + s
  union_scaled[disjoint] <- r_power[disjoint] + s_power[disjoint]
  partial <- !disjoint & delta > abs(r - s)
  if (any(partial)) {
    rr <- r[partial]
    ss <- s[partial]
    dd <- delta[partial]
    # Signed distances to the common chord/hyperplane, in units of scale.
    x <- (dd + (rr - ss) * (rr + ss) / dd) / 2
    y <- dd - x
    cap_fraction <- function(signed_ratio) {
      signed_ratio <- pmax(-1, pmin(1, signed_ratio))
      small_cap <- stats::pbeta(signed_ratio^2, 1 / 2, (d + 1) / 2,
                               lower.tail = FALSE) / 2
      ifelse(signed_ratio >= 0, small_cap, 1 - small_cap)
    }
    intersection <- r_power[partial] * cap_fraction(x / rr) +
      s_power[partial] * cap_fraction(y / ss)
    # Correct only round-off at the geometric bounds, never a statistical
    # estimate or its Jensen bound.
    tolerance <- 100 * .Machine$double.eps
    upper <- pmin(r_power[partial], s_power[partial])
    if (any(intersection < -tolerance | intersection > upper + tolerance)) {
      stop("Ball intersection arithmetic violated its geometric bounds.", call. = FALSE)
    }
    intersection <- pmax(0, pmin(upper, intersection))
    union_scaled[partial] <- r_power[partial] + s_power[partial] - intersection
  }
  answer[positive] <- log_unit_ball + d * log(scale) + log(union_scaled)
  answer
}

.wm_ball_union_volume <- function(radius1, radius2, separation, d) {
  exp(.wm_log_ball_union(radius1, radius2, separation, d))
}

.wm_angular_component <- function(normal_draws, M, d, overlap) {
  private <- M - 1L - overlap
  point_count <- 2 * M - overlap
  sphere_dimension <- d * point_count
  x <- normal_draws[, seq_len(sphere_dimension), drop = FALSE]
  norm <- sqrt(rowSums(x * x))
  if (any(!is.finite(norm)) || any(norm <= 0)) {
    stop("Uniform-sphere generation encountered an invalid Gaussian norm.", call. = FALSE)
  }
  x <- x / norm
  point <- function(index) x[, (index - 1L) * d + seq_len(d), drop = FALSE]
  norm_at <- function(v) sqrt(rowSums(v * v))
  y <- point(1L)
  t <- point(2L)
  radius_y <- norm_at(y)
  radius_t <- norm_at(t)
  if (overlap > 0L) {
    for (index in seq_len(overlap) + 2L) {
      donor <- point(index)
      radius_y <- pmax(radius_y, norm_at(donor - y))
      radius_t <- pmax(radius_t, norm_at(donor - t))
    }
  }
  if (private > 0L) {
    for (index in seq_len(private)) {
      radius_y <- pmax(radius_y, norm_at(point(2L + overlap + index) - y))
      radius_t <- pmax(radius_t, norm_at(point(2L + overlap + private + index) - t))
    }
  }
  include <- rep(TRUE, nrow(x))
  if (private > 0L) {
    for (index in seq_len(private)) {
      include <- include &
        norm_at(point(2L + overlap + index) - t) > radius_t &
        norm_at(point(2L + overlap + private + index) - y) > radius_y
    }
  }
  values <- numeric(nrow(x))
  log_area <- log(2) + sphere_dimension * log(pi) / 2 - lgamma(sphere_dimension / 2)
  log_prefactor <- lgamma(point_count) + log_area - log(d) -
    lgamma(overlap + 1) - 2 * lgamma(private + 1)
  if (any(include)) {
    log_union <- .wm_log_ball_union(radius_y[include], radius_t[include],
                                    norm_at(y - t)[include], d)
    log_value <- log_prefactor - point_count * log_union
    if (any(!is.finite(log_value)) || any(log_value > log(.Machine$double.xmax))) {
      stop("Angular integrand exceeds numerical range; use a specialized integration method.",
           call. = FALSE)
    }
    values[include] <- exp(log_value)
    if (any(values[include] == 0)) {
      stop("An included angular integrand underflowed; use a specialized integration method.",
           call. = FALSE)
    }
  }
  values
}

.wm_geometry_mc <- function(M, d, draws, chunk_size) {
  count <- 0L
  average <- numeric(M)
  cross_deviation <- matrix(0, M, M)
  alpha_average <- alpha_deviation <- 0
  nonzero <- numeric(M)
  max_value <- numeric(M)
  max_dimension <- 2 * M * d
  # This cap keeps accumulated squares and cross-products representable.
  value_limit <- sqrt(.Machine$double.xmax) / (4 * M * sqrt(draws))
  while (count < draws) {
    batch <- min(chunk_size, draws - count)
    normals <- matrix(stats::rnorm(batch * max_dimension), nrow = batch, byrow = TRUE)
    values <- matrix(0, batch, M)
    for (overlap in 0:(M - 1L)) {
      values[, overlap + 1L] <- .wm_angular_component(normals, M, d, overlap)
    }
    if (any(values > value_limit)) {
      stop("Angular sample moments exceed numerical range; use specialized integration.",
           call. = FALSE)
    }
    nonzero <- nonzero + colSums(values > 0)
    max_value <- pmax(max_value, apply(values, 2L, max))
    batch_average <- colMeans(values)
    centered <- sweep(values, 2L, batch_average, "-")
    change <- batch_average - average
    updated <- count + batch
    cross_deviation <- cross_deviation + crossprod(centered) +
      tcrossprod(change) * (count * as.double(batch) / updated)
    average <- average + change * (batch / updated)
    alpha_values <- rowSums(values)
    alpha_batch_average <- mean(alpha_values)
    alpha_change <- alpha_batch_average - alpha_average
    alpha_deviation <- alpha_deviation + sum((alpha_values - alpha_batch_average)^2) +
      alpha_change^2 * (count * as.double(batch) / updated)
    alpha_average <- alpha_average + alpha_change * (batch / updated)
    count <- updated
  }
  covariance <- cross_deviation / (draws - 1) / draws
  list(beta = average, beta_mcse = sqrt(pmax(0, diag(covariance))),
       covariance = covariance, alpha_mcse = sqrt(max(0, alpha_deviation / (draws - 1) / draws)),
       nonzero = nonzero, max_value = max_value)
}

#' Universal Euclidean matching overlap constants
#'
#' @export
wm_geometry <- function(M, d, draws = 4096L, seed = NULL, chunk_size = 256L) {
  positive_integer <- function(x) {
    is.numeric(x) && !is.complex(x) && length(x) == 1L && is.finite(x) &&
      x >= 1 && x <= .Machine$integer.max && x == floor(x)
  }
  for (name in c("M", "d", "draws", "chunk_size")) {
    if (!positive_integer(get(name))) {
      stop(name, " must be one positive integer.", call. = FALSE)
    }
  }
  if (!is.null(seed) && (!is.numeric(seed) || is.complex(seed) || length(seed) != 1L ||
                        !is.finite(seed) || seed < 0 || seed > .Machine$integer.max || seed != floor(seed))) {
    stop("seed must be NULL or one nonnegative integer.", call. = FALSE)
  }
  M <- as.integer(M)
  d <- as.integer(d)
  draws <- as.integer(draws)
  requested_chunk <- as.integer(chunk_size)
  if (M > 1000000L) {
    stop("M exceeds the one-million-component output allocation cap.", call. = FALSE)
  }
  overlap <- seq_len(M) - 1L
  names_overlap <- as.character(overlap)
  definition <- "beta_l = Gamma(2M-l) |S^{d(2M-l)-1}| E[I(Theta) L(Theta)^(-(2M-l))] / (d l! ((M-1-l)!)^2); alpha = sum_l beta_l"
  if (d == 1L) {
    beta <- 2 * (overlap + 1)
    beta[M] <- 3 * M / 2
    names(beta) <- names_overlap
    return(structure(list(M = M, d = d, overlap = overlap, beta = beta,
      alpha = M^2 + M / 2, beta_mcse = stats::setNames(numeric(M), names_overlap),
      alpha_mcse = 0, covariance = NULL, method = "exact_1d_closed_form",
      precision_status = "exact_formula_evaluated_in_floating_point",
      component_precision = rep("exact_formula", M), draws = 0L,
      requested_draws = draws, seed = seed, chunk_size = 0L,
      requested_chunk_size = requested_chunk, definition = definition,
      diagnostics = list(jensen_lower_bound = M^2, below_jensen = FALSE,
                         rng_used = FALSE)), class = c("wm_geometry", "list")))
  }
  if (draws < 2L) stop("General-dimensional Monte Carlo requires draws >= 2.", call. = FALSE)
  max_dimension <- 2 * as.double(M) * d
  total_coordinates <- as.double(draws) * d * (3 * as.double(M)^2 + M) / 2
  if (M > 128L || max_dimension > 1000000 || total_coordinates > 50000000) {
    stop("Requested integration exceeds the bounded work cap (M <= 128, at most 50 million component coordinates); reduce M, d or draws.",
         call. = FALSE)
  }
  chunk_size <- as.integer(min(requested_chunk, draws, floor(1000000 / max_dimension)))
  if (chunk_size < 1L) stop("Sphere dimension exceeds the allocation cap.", call. = FALSE)
  if (!is.null(seed)) {
    if (identical(RNGkind()[2L], "Box-Muller")) {
      stop("An explicit seed cannot preserve the Box-Muller cache; use seed = NULL or another normal RNG.",
           call. = FALSE)
    }
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    previous_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
    on.exit({
      if (had_seed) assign(".Random.seed", previous_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(as.integer(seed))
  }
  result <- .wm_geometry_mc(M, d, draws, chunk_size)
  names(result$beta) <- names(result$beta_mcse) <- names_overlap
  dimnames(result$covariance) <- list(names_overlap, names_overlap)
  alpha <- sum(result$beta)
  zero_hits <- result$nonzero == 0
  component_precision <- ifelse(zero_hits, "unresolved_zero_hits", "monte_carlo_estimate")
  point_count <- 2 * M - overlap
  sphere_dimension <- d * point_count
  normalization <- data.frame(overlap = overlap, point_count = point_count,
    sphere_dimension = sphere_dimension,
    log_sphere_area = log(2) + sphere_dimension * log(pi) / 2 - lgamma(sphere_dimension / 2),
    log_prefactor = lgamma(point_count) - log(d) - lgamma(overlap + 1) -
      2 * lgamma(M - overlap) + log(2) + sphere_dimension * log(pi) / 2 - lgamma(sphere_dimension / 2))
  structure(list(M = M, d = d, overlap = overlap, beta = result$beta,
    alpha = alpha, beta_mcse = result$beta_mcse, alpha_mcse = result$alpha_mcse,
    covariance = result$covariance, method = "uniform_sphere_monte_carlo",
    precision_status = if (any(zero_hits)) "unresolved_components" else "monte_carlo_estimate",
    component_precision = component_precision, draws = draws, requested_draws = draws,
    seed = seed, chunk_size = chunk_size, requested_chunk_size = requested_chunk,
    definition = definition, normalization = normalization,
    diagnostics = list(nonzero_draws = stats::setNames(result$nonzero, names_overlap),
      maximum_draw_value = stats::setNames(result$max_value, names_overlap),
      jensen_lower_bound = M^2, below_jensen = alpha < M^2,
      rng_used = TRUE, coupled_components = TRUE,
      gaussian_allocation_cap = 1000000, coordinate_work_cap = 50000000,
      note = "MCSE is empirical, not a rigorous error bound. Zero-hit component MCSEs are uninformative. Estimates are not projected onto the Jensen bound.")),
    class = c("wm_geometry", "list"))
}
