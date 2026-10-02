# Deterministic formula and guard checks. This is not coverage simulation.
source("simulations/fitted_wm_components.R")
checks <- 0L
near <- function(x, y, tol = 1e-10) {
  checks <<- checks + 1L
  if (length(x) != length(y) || any(!is.finite(c(x, y))) ||
      max(abs(x - y)) > tol * max(1, abs(x), abs(y))) {
    stop("Numeric check failed: ", checks)
  }
}
reject <- function(expr) {
  checks <<- checks + 1L
  if (!inherits(try(force(expr), silent = TRUE), "try-error")) {
    stop("Missing guard: ", checks)
  }
}
# Independent query-level estimator, keeping the supplied donor lists fixed.
frozen_point <- function(fit, W, mean0, mean1) {
  Z <- fit$data$Z
  Y <- fit$data$Y
  queries <- if (fit$estimand == "PATE") seq_along(Y) else which(Z == 1L)
  pieces <- vapply(queries, function(i) {
    donors <- fit$graph$edges$donor[fit$graph$edges$query == i]
    mu <- if (Z[i] == 1L) mean0 else mean1
    missing <- mu[i] + sum(W[donors] * (Y[donors] - mu[donors])) /
      sum(W[donors])
    W[i] * (2 * Z[i] - 1) * (Y[i] - missing)
  }, numeric(1))
  sum(pieces) / sum(W[queries])
}
n <- 18L
i <- seq_len(n)
Z <- as.integer(i %% 2)
Y <- sin(i / 2) + i / 7 + Z * (1 + cos(i / 5))
W <- .7 + (i %% 7) / 6
S0 <- cbind(i / n, sin(i))
S1 <- cbind(cos(i / 3), i^2 / n^2, sin(i / 5))
mu0 <- i / 9 + cos(i / 4)
mu1 <- 1 + i / 8 - sin(i / 3)
dw <- cbind(location = W * sin(i / 4), scale = W * cos(i / 3))
d0 <- cbind(location = cos(i / 5), scale = i / 13)
d1 <- cbind(location = i / 11, scale = sin(i / 6))
ell <- cbind(location = sin(i / 3), scale = cos(i / 4))
configs <- 0L
max_gradient_error <- 0
for (estimand in c("PATE", "PATT")) {
  for (M in c(1L, 3L, 5L)) {
    make_fit <- function(weights = W, variance = TRUE) {
      wdsmatch::wm_match(Y = Y, Z = Z, weights = weights,
        scores0 = S0, scores1 = if (estimand == "PATE") S1 else NULL,
        estimand = estimand, M = M, mean0 = mu0,
        mean1 = if (estimand == "PATE") mu1 else NULL,
        variance = variance)
    }
    fit <- make_fit()
    g <- wm_graph_gradient_reference(fit, dw, d0,
      if (estimand == "PATE") d1 else NULL)
    step <- 1e-5
    numeric_gradient <- vapply(seq_len(ncol(dw)), function(k) {
      (frozen_point(fit, W + step * dw[, k], mu0 + step * d0[, k],
                    mu1 + step * d1[, k]) -
       frozen_point(fit, W - step * dw[, k], mu0 - step * d0[, k],
                    mu1 - step * d1[, k])) / (2 * step)
    }, numeric(1))
    max_gradient_error <- max(max_gradient_error,
                              abs(g$sensitivity - numeric_gradient))
    near(g$sensitivity, numeric_gradient, 2e-8)
    near(g$sensitivity, g$weight_sensitivity + g$prediction_sensitivity)
    near(g$max_fraction_derivative_sum, 0)
    zero_if <- matrix(0, n, 2L, dimnames = list(NULL, colnames(dw)))
    old_weight <- wdsmatch::wm_weight_adjust(fit,
      weight_derivative = dw, nuisance_influence = zero_if)
    near(g$weight_sensitivity, old_weight$weight_adjustment$sensitivity)
    near(wm_graph_gradient_reference(fit, matrix(W, n, 1L))$sensitivity, 0)
    near(wm_graph_gradient_reference(fit,
      mean0_derivative = matrix(1, n, 1L))$sensitivity, 0)
    for (factor in c(1e-100, 1e100)) {
      scaled <- make_fit(W * factor, variance = FALSE)
      gs <- wm_graph_gradient_reference(scaled, dw * factor, d0,
        if (estimand == "PATE") d1 else NULL)
      near(gs$sensitivity, g$sensitivity)
    }
    graph <- c(location = .3, scale = -.2)
    scope <- if (estimand == "PATE") "distinct_rarity" else "patt"
    inference <- wm_fitted_variance_reference(fit, ell,
      g$sensitivity, graph, scope)
    base <- .wm_fv_state(fit)$rows
    base <- base - mean(base)
    centered_if <- sweep(ell, 2L, colMeans(ell), "-")
    b <- g$sensitivity + graph
    expected <- mean((base + as.vector(centered_if %*% b))^2)
    near(inference$root_n_variance, expected)
    near(inference$root_n_variance,
      inference$V0 + 2 * sum(b * inference$C) +
        as.numeric(crossprod(b, inference$Sigma %*% b)))
    near(inference$variance, expected / n)
    near(inference$estimate, fit$estimate)
    shifted <- wm_fitted_variance_reference(fit,
      sweep(ell, 2L, c(10, -3), "+"), g$sensitivity, graph, scope)
    near(shifted$root_n_variance, inference$root_n_variance)
    configs <- configs + 1L
  }
}

# Enumerate all 4^4 ordered draws. Their equal probabilities induce exactly
# Multinomial(4, 1/4) probabilities on count vectors, without random sampling.
small <- wdsmatch::wm_match(
  Y = c(0, 4, 2, 7), Z = c(0, 1, 0, 1), weights = c(1.1, .9, 1.4, 1.2),
  scores0 = cbind(c(0, 1, 3, 4), c(0, .2, -.4, .9)), M = 1L,
  mean0 = c(0, .3, .4, .5), mean1 = c(1, 1.2, 1.3, 1.5),
  variance = FALSE)
small_if <- cbind(alpha = c(-1, 2, 0, 3), beta = c(2, 1, -1, 0))
small_smooth <- c(alpha = .2, beta = -.1)
small_graph <- c(alpha = .4, beta = .3)
small_inf <- wm_fitted_variance_reference(small, small_if,
  small_smooth, small_graph, "common_field")
high_confidence <- wm_fitted_variance_reference(small, small_if,
  small_smooth, small_graph, "common_field",
  conf.level = 1 - .Machine$double.eps / 2)
near(as.numeric(all(is.finite(high_confidence$conf.int))), 1)
ordered <- as.matrix(expand.grid(rep(list(seq_len(4L)), 4L)))
counts <- apply(ordered, 1L, tabulate, nbins = 4L)
boot <- wm_fitted_count_reference(small_inf, counts)
near(mean(boot$root_n_draws), 0)
near(mean(boot$root_n_draws^2), small_inf$root_n_variance)
near(boot$conditional_variance, small_inf$root_n_variance / 4)
near(boot$monte_carlo_root_n_variance, small_inf$root_n_variance)
near(boot$draws, small_inf$estimate +
  as.vector(crossprod(small_inf$augmented_rows, counts - 1)) / 4)
empty_arm <- colSums(counts[c(1L, 3L), , drop = FALSE]) == 0 |
  colSums(counts[c(2L, 4L), , drop = FALSE]) == 0
near(sum(empty_arm), 32)
near(as.numeric(all(is.finite(boot$draws[empty_arm]))), 1)
single <- wm_fitted_count_reference(small_inf, counts[, 1L, drop = FALSE])
near(single$conditional_root_n_variance, small_inf$root_n_variance)
near(single$monte_carlo_root_n_variance, 0)

reject(wm_graph_gradient_reference(small))
reject(wm_graph_gradient_reference(small, matrix(0, 3L, 2L)))
duplicated <- matrix(0, 4L, 2L, dimnames = list(NULL, c("x", "x")))
reject(wm_graph_gradient_reference(small, duplicated))
reject(wm_fitted_variance_reference(small, small_if,
  small_smooth, small_graph, "patt"))
reject(wm_fitted_variance_reference(small, small_if,
  small_smooth, small_graph, "full_x"))
reject(wm_fitted_variance_reference(small, small_if,
  small_smooth, small_graph, "unspecified"))
reject(wm_fitted_variance_reference(small, small_if,
  rev(small_smooth), small_graph, "common_field"))
scalar <- wdsmatch::wm_match(Y = c(0, 4, 2, 7), Z = c(0, 1, 0, 1),
  weights = rep(1, 4), scores0 = matrix(1:4, 4L, 1L), M = 1L,
  mean0 = rep(0, 4), mean1 = rep(1, 4), variance = FALSE)
reject(wm_fitted_variance_reference(scalar, small_if,
  small_smooth, small_graph, "common_field"))
stale <- small
stale$estimate <- stale$estimate + .5
reject(wm_graph_gradient_reference(stale, matrix(1, 4, 1)))
bad_counts <- counts
bad_counts[1:2, 1] <- c(.5, 3.5)
reject(wm_fitted_count_reference(small_inf, bad_counts))
bad_counts <- counts
bad_counts[1, 1] <- bad_counts[1, 1] + 1
reject(wm_fitted_count_reference(small_inf, bad_counts))
tampered <- small_inf
tampered$root_n_variance <- tampered$root_n_variance + 1
reject(wm_fitted_count_reference(tampered, counts))
bad_if <- small_if
bad_if[1, 1] <- Inf
reject(wm_fitted_variance_reference(small, bad_if,
  small_smooth, small_graph, "common_field"))

result <- list(passed = TRUE, deterministic_checks = checks,
  gradient_configurations = configs, max_gradient_error = max_gradient_error,
  exact_ordered_multinomial_assignments = nrow(ordered),
  simulation_replications = 0L,
  wdsmatch_version = as.character(utils::packageVersion("wdsmatch")))
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  dir.create(dirname(args[1L]), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(result, args[1L], auto_unbox = TRUE, pretty = TRUE)
}
print(result)
