# Source first_stage.R for configuration/integer helpers. No RNG is used here.

.wm_fs_geometry <- function(M, d, geometry) {
  if (d == 1L) {
    exact <- M^2 + M / 2
    if (!is.null(geometry) &&
        (!isTRUE(all.equal(geometry$M, M)) || !isTRUE(all.equal(geometry$d, d)) ||
         length(geometry$alpha) != 1L || !is.finite(geometry$alpha) ||
         abs(geometry$alpha - exact) > 1e-12 * exact))
      stop("supplied geometry disagrees with the exact one-dimensional constant")
    return(list(alpha = exact, alpha_mcse = 0, reported_mcse = 0,
                status = "exact_1d", draws = 0L, seed = NA_integer_,
                below_jensen = FALSE))
  }
  if (is.null(geometry)) return(list(alpha = NA_real_, alpha_mcse = NA_real_,
    reported_mcse = NA_real_, status = "geometry_not_provided",
    draws = NA_integer_, seed = NA_integer_, below_jensen = NA))
  alpha_only <- inherits(geometry, "wm_geometry_alpha")
  if ((!inherits(geometry, "wm_geometry") && !alpha_only) ||
      !isTRUE(all.equal(geometry$M, M)) || !isTRUE(all.equal(geometry$d, d)) ||
      length(geometry$alpha) != 1L || !is.finite(geometry$alpha) || geometry$alpha <= 0 ||
      length(geometry$alpha_mcse) != 1L || !is.finite(geometry$alpha_mcse) ||
      geometry$alpha_mcse < 0)
    stop("geometry must be a matching finite wm_geometry or wm_geometry_alpha result")
  if (alpha_only) {
    bound <- geometry$diagnostics$integrand_upper_bound
    maximum <- geometry$diagnostics$maximum_draw_value
    if (length(bound) != 1L || !is.finite(bound) || bound <= 0 ||
        length(maximum) != 1L || !is.finite(maximum) || maximum <= 0 ||
        maximum > bound * (1 + 1e-10) || geometry$alpha > maximum * (1 + 1e-10))
      stop("alpha-only geometry has inconsistent bounded-integrand diagnostics")
  }
  zero_hits <- any(geometry$component_precision == "unresolved_zero_hits") ||
    (!is.null(geometry$diagnostics$nonzero_draws) &&
     any(geometry$diagnostics$nonzero_draws == 0))
  meaningful <- !zero_hits && geometry$alpha_mcse > 0 &&
    identical(geometry$precision_status, "monte_carlo_estimate")
  list(alpha = geometry$alpha,
       alpha_mcse = if (meaningful) geometry$alpha_mcse else NA_real_,
       reported_mcse = geometry$alpha_mcse,
       status = if (meaningful) "monte_carlo_estimate" else "unresolved_geometry_precision",
       draws = geometry$draws, seed = if (is.null(geometry$seed)) NA_integer_ else geometry$seed,
       below_jensen = geometry$alpha < M^2,
       method = geometry$method,
       integrand_upper_bound = if (alpha_only) geometry$diagnostics$integrand_upper_bound else NA_real_,
       true_mcse_upper_bound = if (alpha_only) geometry$diagnostics$true_mcse_upper_bound else NA_real_,
       maximum_draw_value = if (alpha_only) geometry$diagnostics$maximum_draw_value else NA_real_)
}

wm_first_stage_weight_integrals <- function() {
  q <- function(t) stats::plogis(3 * t)
  J <- function(t) 1 / q(t) + 1 / (1 - q(t))
  functions <- list(
    I00 = function(t) q(t) * (1 - q(t)),
    I01 = function(t) t * q(t) * (1 - q(t)),
    I11 = function(t) t^2 * q(t) * (1 - q(t)),
    gA0 = function(t) (q(t) - 0.5) * t,
    gA1 = function(t) (q(t) - 0.5) * t^2,
    gT0 = function(t) (q(t) - 1) * t,
    gT1 = function(t) (q(t) - 1) * t^2,
    J = J, t2J = function(t) t^2 * J(t),
    inverse_q = function(t) 1 / q(t),
    inverse_control = function(t) 1 / (1 - q(t)),
    t2_inverse_q = function(t) t^2 / q(t))
  results <- lapply(functions, function(fun) stats::integrate(fun, -0.5, 0.5,
    rel.tol = 1e-11, abs.tol = 1e-12, subdivisions = 200L, stop.on.error = TRUE))
  value <- vapply(results, function(x) x$value, numeric(1L))
  error <- vapply(results, function(x) x$abs.error, numeric(1L))
  I <- matrix(c(value[["I00"]], value[["I01"]], value[["I01"]], value[["I11"]]), 2L)
  dimnames(I) <- list(c("intercept", "t"), c("intercept", "t"))
  factor <- chol(I)
  gA <- c(intercept = value[["gA0"]], t = value[["gA1"]])
  gT <- c(intercept = value[["gT0"]], t = value[["gT1"]])
  inverse_times <- function(g) as.vector(backsolve(factor, forwardsolve(t(factor), g)))
  reduction <- c(PATE = sum(gA * inverse_times(gA)), PATT = sum(gT * inverse_times(gT)))
  if (any(!is.finite(c(value, error, reduction))) || any(reduction < 0))
    stop("invalid integrated first-stage weight benchmark")
  list(values = value, absolute_errors = error, information = I,
       sensitivity = list(PATE = gA, PATT = gT),
       covariance = list(PATE = -inverse_times(gA), PATT = -inverse_times(gT)),
       variance_reduction = reduction,
       quadrature = list(relative_tolerance = 1e-11, absolute_tolerance = 1e-12,
                         subdivisions = vapply(results, function(x) x$subdivisions, integer(1L))),
       scope = "Deterministic quadrature error estimates, separate from geometry Monte Carlo uncertainty")
}

wm_first_stage_benchmark <- function(branch, d, M, n = 200L, m = n,
                                      geometry = NULL) {
  branch <- match.arg(branch, c("same_sample", "independent_training", "estimated_weights", "gaussian_same_sample"))
  d <- .wm_fs_integer(d, "d")
  M <- .wm_fs_integer(M, "M")
  n <- .wm_fs_integer(n, "n", 2L)
  m <- .wm_fs_integer(m, "m", 2L)
  if (branch != "independent_training" && m != n) stop("this branch requires m = n")
  if (branch == "same_sample" && d == 1L) stop("same-sample changing scalar scores are unsupported")
  g <- .wm_fs_geometry(M, d, geometry)
  sigma2 <- 1 / 16
  r <- 0.6
  integral <- NULL
  if (branch == "estimated_weights") {
    integral <- wm_first_stage_weight_integrals()
    v <- integral$values
    constant <- c(PATE = 0.25 * v[["t2J"]] +
                    sigma2 * (0.75 + 0.25 / M) * v[["J"]],
                  PATT = v[["t2_inverse_q"]] +
                    sigma2 * (1 + 1 / M) * v[["inverse_q"]])
    coefficient <- c(PATE = sigma2 * 0.25 / M^2 * v[["J"]],
                     PATT = sigma2 / M^2 * v[["inverse_control"]])
    change <- -integral$variance_reduction
    target <- 1.5
    constants <- list(sigma2 = sigma2, target_pi = 0.5, true_parameter = c(0, 3))
  } else {
    constant <- c(PATE = if (branch == "gaussian_same_sample") sigma2 * (26/9 + 10/(9*M)) else sigma2 * (3 + 1 / M),
                  PATT = 2 * sigma2 * (1 + 1 / M))
    coefficient <- c(PATE = if (branch == "gaussian_same_sample") sigma2 * 10/(9*M^2) else sigma2 / M^2,
                     PATT = 2 * sigma2 / M^2)
    magnitude <- (10 / 3) * r^2 * sigma2
    change <- rep(if (branch %in% c("same_sample", "gaussian_same_sample")) -magnitude else (n / m) * magnitude, 2L)
    names(change) <- c("PATE", "PATT")
    target <- 1
    constants <- list(r = r, sigma2 = sigma2, true_parameter = 0,
      population_D2 = 1 / 30, prediction_sensitivity = -r / 3,
      parameter_variance = 30 * sigma2, same_sample_covariance = 10 * r * sigma2,
      independent_training_multiplier = n / m)
  }
  rows <- list()
  for (estimand in c("PATE", "PATT"))
    for (comparison in c("known_first_stage", "fitted_naive", "fitted_adjusted")) {
      sampling_change <- if (comparison == "known_first_stage") 0 else change[[estimand]]
      reported_change <- if (comparison == "fitted_adjusted") change[[estimand]] else 0
      a <- coefficient[[estimand]]
      sampling_constant <- constant[[estimand]] + sampling_change
      reported_constant <- constant[[estimand]] + reported_change
      rows[[length(rows) + 1L]] <- data.frame(branch = branch, d = d, M = M, n = n, m = m,
        estimand = estimand, comparison = comparison, target = target,
        probability_limit = target, sampling_variance_constant = sampling_constant,
        reported_variance_constant = reported_constant, alpha_coefficient = a,
        first_stage_sampling_change = sampling_change,
        benchmark_root_n_variance = sampling_constant + a * g$alpha,
        reported_root_n_variance_limit = reported_constant + a * g$alpha,
        benchmark_geometry_mcse = a * g$alpha_mcse,
        geometry_reported_mcse = a * g$reported_mcse,
        alpha = g$alpha, alpha_mcse = g$alpha_mcse, geometry_status = g$status,
        geometry_draws = g$draws, geometry_seed = g$seed, geometry_below_jensen = g$below_jensen,
        geometry_method = if (is.null(g$method)) g$status else g$method,
        geometry_integrand_upper_bound = if (is.null(g$integrand_upper_bound)) NA_real_ else g$integrand_upper_bound,
        geometry_true_mcse_upper_bound = if (is.null(g$true_mcse_upper_bound)) NA_real_ else g$true_mcse_upper_bound,
        geometry_maximum_draw_value = if (is.null(g$maximum_draw_value)) NA_real_ else g$maximum_draw_value,
        inference_status = if (comparison == "fitted_naive")
          "first_stage_omission_diagnostic" else "within_declared_parametric_scope",
        stringsAsFactors = FALSE)
    }
  table <- do.call(rbind, rows)
  if (any(table$benchmark_root_n_variance < 0, na.rm = TRUE) ||
      any(table$reported_root_n_variance_limit < 0, na.rm = TRUE))
    stop("supplied geometry produced a negative variance benchmark")
  list(table = table, constants = constants, weight_integrals = integral,
       geometry = g, interpretation = paste(
         "Sampling variance and reported variance limits differ for fitted_naive.",
         "Propagated geometry uncertainty is independent of quadrature error;",
         "unresolved Monte Carlo precision is retained, not replaced by zero."))
}
