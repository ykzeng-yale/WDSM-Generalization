# The Poisson intersection count below is not the complete-set overlap beta.
.wm_alpha_integrand_bound <- function(M, d) {
  exp((d + 1) * log(2) + 2 * lgamma(d / 2 + 1) - log(d) - lgamma(d)) *
    M * (2 * M - 1)
}

.wm_alpha_integrand <- function(radius1, radius2, separation, M, d) {
  log_union <- .wm_log_ball_union(radius1, radius2, separation, d)
  log_ball <- d * log(pi) / 2 - lgamma(d / 2 + 1)
  x_fraction <- exp(log_ball + d * log(radius1) - log_union)
  y_fraction <- exp(log_ball + d * log(radius2) - log_union)
  tolerance <- 1000 * .Machine$double.eps
  intersection <- x_fraction + y_fraction - 1
  if (any(!is.finite(c(x_fraction, y_fraction, intersection))) ||
      any(x_fraction > 1 + tolerance | y_fraction > 1 + tolerance |
          intersection < -tolerance))
    stop("Reduced-alpha volume fractions violated their geometric bounds.", call. = FALSE)
  # Correct only floating-point geometry error at the simplex boundary.
  x_fraction <- pmin(1, x_fraction)
  y_fraction <- pmin(1, y_fraction)
  p0 <- pmax(0, x_fraction + y_fraction - 1)
  p1 <- 1 - y_fraction
  p2 <- 1 - x_fraction
  total <- p0 + p1 + p2
  p0 <- p0 / total
  p1 <- p1 / total
  p2 <- p2 / total
  powers <- function(p) {
    result <- matrix(1, length(p), M)
    if (M > 1L) for (j in 2:M) result[, j] <- result[, j - 1L] * p
    result
  }
  shared <- powers(p0)
  first <- powers(p1)
  second <- powers(p2)
  polynomial <- numeric(length(radius1))
  for (ell in 0:(M - 1L)) {
    for (a in 0:(M - 1L - ell)) for (b in 0:(M - 1L - ell)) {
      N <- ell + a + b
      coefficient <- exp(lgamma(N + 2) - lgamma(ell + 1) -
                           lgamma(a + 1) - lgamma(b + 1))
      polynomial <- polynomial +
        coefficient * shared[, ell + 1L] * first[, a + 1L] * second[, b + 1L]
    }
  }
  log_prefactor <- log(2) + d * log(pi) - lgamma(d) - log(d) - 2 * log_union
  result <- exp(log_prefactor) * polynomial
  bound <- .wm_alpha_integrand_bound(M, d)
  if (any(!is.finite(result)) || any(result <= 0) ||
      any(result > bound * (1 + 1e-10)))
    stop("Reduced-alpha integrand exceeded its finite numerical bounds.", call. = FALSE)
  result
}

.wm_alpha_mc <- function(M, d, draws, chunk_size) {
  count <- 0L
  average <- deviation <- 0
  minimum <- Inf
  maximum <- 0
  while (count < draws) {
    batch <- min(chunk_size, draws - count)
    normals <- matrix(stats::rnorm(batch * 2 * d), nrow = batch, byrow = TRUE)
    norm <- sqrt(rowSums(normals^2))
    if (any(!is.finite(norm)) || any(norm <= 0))
      stop("Reduced-alpha sphere generation encountered an invalid norm.", call. = FALSE)
    normalized <- normals / norm
    x <- normalized[, seq_len(d), drop = FALSE]
    y <- normalized[, d + seq_len(d), drop = FALSE]
    value <- .wm_alpha_integrand(sqrt(rowSums(x^2)), sqrt(rowSums(y^2)),
                                 sqrt(rowSums((x - y)^2)), M, d)
    batch_mean <- mean(value)
    change <- batch_mean - average
    updated <- count + batch
    deviation <- deviation + sum((value - batch_mean)^2) +
      change^2 * (count * as.double(batch) / updated)
    average <- average + change * (batch / updated)
    minimum <- min(minimum, value)
    maximum <- max(maximum, value)
    count <- updated
  }
  list(alpha = average, alpha_mcse = sqrt(max(0, deviation / (draws - 1) / draws)),
       minimum = minimum, maximum = maximum)
}

#' Reduced Poisson-count integration of the universal alpha constant
#' @export
wm_geometry_alpha <- function(M, d, draws = 4096L, seed = NULL, chunk_size = 256L) {
  positive_integer <- function(x) is.numeric(x) && !is.complex(x) &&
    length(x) == 1L && is.finite(x) && x >= 1 && x == floor(x) &&
    x <= .Machine$integer.max
  for (name in c("M", "d", "draws", "chunk_size")) {
    if (!positive_integer(get(name)))
      stop(name, " must be one positive integer.", call. = FALSE)
  }
  if (!is.null(seed) && (!is.numeric(seed) || is.complex(seed) ||
      length(seed) != 1L || !is.finite(seed) || seed < 0 ||
      seed > .Machine$integer.max || seed != floor(seed)))
    stop("seed must be NULL or one nonnegative integer.", call. = FALSE)
  M <- as.integer(M)
  d <- as.integer(d)
  draws <- as.integer(draws)
  requested_chunk <- as.integer(chunk_size)
  definition <- paste("alpha = E[|S^(2d-1)|/(d U^2) * sum_{l=0}^{M-1}",
    "sum_{a,b=0}^{M-1-l} (l+a+b+1)! p0^l p1^a p2^b/(l!a!b!)];",
    "p0=intersection/union, p1=first-exclusive/union, p2=second-exclusive/union")
  if (d == 1L) {
    return(structure(list(M = M, d = d, alpha = as.double(M)^2 + M / 2,
      alpha_mcse = 0, method = "exact_1d_closed_form",
      precision_status = "exact_formula_evaluated_in_floating_point", draws = 0L,
      requested_draws = draws, seed = seed, chunk_size = 0L,
      requested_chunk_size = requested_chunk, definition = definition,
      scope = "alpha only; no overlap beta components",
      diagnostics = list(jensen_lower_bound = as.double(M)^2, below_jensen = FALSE,
                         rng_used = FALSE)), class = c("wm_geometry_alpha", "list")))
  }
  terms <- as.double(M) * (as.double(M) + 1) * (2 * as.double(M) + 1) / 6
  if (draws < 2L) stop("Reduced-alpha Monte Carlo requires draws >= 2.", call. = FALSE)
  if (M > 128L || d > 128L || draws * terms > 50000000 ||
      as.double(draws) * 2 * d > 50000000)
    stop("Reduced-alpha integration exceeds its bounded work cap (M,d <= 128; ",
         "at most 50 million count terms and Gaussian coordinates).", call. = FALSE)
  chunk_size <- as.integer(min(requested_chunk, draws,
                               floor(1000000 / max(2 * d, 3 * M))))
  if (!is.null(seed)) {
    if (identical(RNGkind()[2L], "Box-Muller"))
      stop("An explicit seed cannot preserve the Box-Muller cache; use seed = NULL or another normal RNG.",
           call. = FALSE)
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    previous_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
    on.exit({
      if (had_seed) assign(".Random.seed", previous_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    })
    set.seed(as.integer(seed))
  }
  result <- .wm_alpha_mc(M, d, draws, chunk_size)
  bound <- .wm_alpha_integrand_bound(M, d)
  structure(list(M = M, d = d, alpha = result$alpha, alpha_mcse = result$alpha_mcse,
    method = "reduced_poisson_count_sphere_monte_carlo",
    precision_status = if (result$alpha_mcse > 0) "monte_carlo_estimate" else
      "unresolved_zero_empirical_variance",
    draws = draws, requested_draws = draws, seed = seed, chunk_size = chunk_size,
    requested_chunk_size = requested_chunk, definition = definition,
    scope = "alpha only; no overlap beta components",
    diagnostics = list(integrand_upper_bound = bound,
      true_mcse_upper_bound = bound / (2 * sqrt(draws)),
      hoeffding_95_half_width = bound * sqrt(log(40) / (2 * draws)),
      minimum_draw_value = result$minimum, maximum_draw_value = result$maximum,
      jensen_lower_bound = as.double(M)^2, below_jensen = result$alpha < as.double(M)^2,
      sphere_dimension = 2 * d, count_terms = terms, rng_used = TRUE,
      gaussian_coordinate_work = as.double(draws) * 2 * d,
      count_term_work = draws * terms,
      note = paste("The direct alpha integrand is positive and bounded; no beta components are returned.",
        "The MCSE bound and Hoeffding width follow from the mathematical integrand bound",
        "and do not account for floating-point error. The estimate is not projected onto Jensen's bound."))),
    class = c("wm_geometry_alpha", "list"))
}
