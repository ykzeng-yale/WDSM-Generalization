gaussian_fixture <- function() {
  n <- 18L
  s <- seq(.025, .975, length.out = n)
  mark <- rep(c(-1, 1, 1), length.out = n)
  z <- rep(c(0L, 1L, 1L, 0L, 1L, 0L), 3L)
  d <- mark * s * (1 - s)
  e <- rep(c(-.04, .03, .01), 6L)
  e <- e - d * sum(d * e) / sum(d^2)
  list(Y = s + z + .21 * d + e, Z = z,
    W = 1 + .5 * (mark == 1) + .2 * z, S = s,
    D = matrix(d, ncol = 1L, dimnames = list(NULL, "theta")))
}

test_that("the actual unrounded OLS parameter determines scalar maps and predictions", {
  a <- gaussian_fixture()
  map <- function(theta) matrix(a$S + theta[[1L]] * a$D[, 1L], ncol = 1L)
  f <- wm_gaussian_match(a$Y, a$Z, a$W, a$D, score_map = map,
    offset0 = a$S, offset1 = a$S + 1, M = 2L)
  expect_equal(unname(f$gaussian_adjustment$parameters), .21, tolerance = 1e-12)
  expect_equal(f$graph$scores0[, 1L], a$S + .21 * a$D[, 1L], tolerance = 1e-12)
  expect_equal(f$predictions$mean0, a$S + .21 * a$D[, 1L], tolerance = 1e-12)
  expect_equal(f$gaussian_adjustment$estimating_equation, 0, tolerance = 1e-12)
  expected_residual <- a$Y - a$S - a$Z - .21 * a$D[, 1L]
  expect_equal(f$gaussian_adjustment$nuisance_influence[, 1L],
    a$D[, 1L] * expected_residual / mean(a$D^2), tolerance = 1e-12)
  baseline <- wm_match(a$Y, a$Z, a$W, matrix(a$S, ncol = 1L), M = 2L,
    mean0 = a$S, mean1 = a$S + 1)
  expect_false(identical(f$graph$edges[c("query", "donor")],
                         baseline$graph$edges[c("query", "donor")]))
  expect_identical(f$info$inference_status, "conditional_gaussian_linear_score_contract")
})

test_that("current-graph finite differences recover Gaussian prediction sensitivity", {
  a <- gaussian_fixture()
  for (estimand in c("PATE", "PATT")) {
    f <- wm_gaussian_match(a$Y, a$Z, a$W, a$D,
      score_map = function(theta) matrix(a$S + theta[[1L]] * a$D[, 1L], ncol = 1L),
      offset0 = a$S, offset1 = a$S + 1, M = 2L, estimand = estimand)
    point_on_fixed_graph <- function(theta) {
      mu0 <- a$S + theta * a$D[, 1L]
      mu1 <- mu0 + 1
      query <- if (estimand == "PATE") seq_along(a$Y) else which(a$Z == 1L)
      total <- 0
      for (i in query) {
        donors <- f$graph$edges$donor[f$graph$edges$query == i]
        mu <- if (a$Z[i] == 1L) mu0 else mu1
        missing <- mu[i] + sum(a$W[donors] * (a$Y[donors] - mu[donors])) / sum(a$W[donors])
        total <- total + a$W[i] * (2 * a$Z[i] - 1) * (a$Y[i] - missing)
      }
      total / sum(a$W[query])
    }
    theta <- f$gaussian_adjustment$parameters[[1L]]
    expect_equal(f$estimate, point_on_fixed_graph(theta), tolerance = 1e-12)
    numerical <- (point_on_fixed_graph(theta + 1e-5) - point_on_fixed_graph(theta - 1e-5)) / 2e-5
    expect_equal(unname(f$gaussian_adjustment$sensitivity), numerical, tolerance = 1e-8)
    g <- f$gaussian_adjustment
    expect_gt(abs(g$cross_term), 1e-9)
    expect_equal(f$root_n_variance,
      g$baseline_root_n_variance + g$cross_term + g$nuisance_root_n_variance,
      tolerance = 1e-12)
    expect_equal(sum(f$contributions$row), 0, tolerance = 1e-12)
    boot <- wm_bootstrap(f, B = 7L, seed = 718L)
    expect_equal(boot$conditional_root_n_variance, f$root_n_variance, tolerance = 1e-12)
  }
})

test_that("Gaussian OLS is independent of a common scaling of matching weights", {
  a <- gaussian_fixture()
  arguments <- list(Y = a$Y, Z = a$Z, weights = a$W, design0 = a$D,
    score_map = function(theta) matrix(a$S + theta[[1L]] * a$D[, 1L], ncol = 1L),
    offset0 = a$S, offset1 = a$S + 1, M = 2L)
  f <- do.call(wm_gaussian_match, arguments)
  arguments$weights <- 1000 * a$W
  scaled <- do.call(wm_gaussian_match, arguments)
  expect_equal(scaled$gaussian_adjustment$parameters, f$gaussian_adjustment$parameters)
  expect_equal(scaled$estimate, f$estimate, tolerance = 1e-12)
  expect_equal(scaled$root_n_variance, f$root_n_variance, tolerance = 1e-12)
})

test_that("different potential designs fit one joint outcome parameter for PATT", {
  a <- gaussian_fixture()
  D0 <- cbind(intercept = 1, effect = 0)
  D0 <- D0[rep(1L, length(a$Y)), , drop = FALSE]
  D1 <- D0
  D1[, "effect"] <- 1
  f <- wm_gaussian_match(a$Y, a$Z, a$W, D0, D1,
    score_map = function(theta) matrix(a$S, ncol = 1L), estimand = "PATT", M = 2L)
  expect_equal(unname(f$gaussian_adjustment$parameters),
    c(mean(a$Y[a$Z == 0L]), mean(a$Y[a$Z == 1L]) - mean(a$Y[a$Z == 0L])), tolerance = 1e-12)
  expect_equal(f$gaussian_adjustment$sensitivity, c(intercept = 0, effect = 0), tolerance = 1e-12)
})

test_that("Gaussian branch rejects unidentified or malformed inputs explicitly", {
  a <- gaussian_fixture()
  map <- function(theta) matrix(a$S, ncol = 1L)
  expect_error(wm_gaussian_match(a$Y, a$Z, a$W, cbind(first = a$D[, 1], second = a$D[, 1]), score_map = map), "rank deficient")
  expect_error(wm_gaussian_match(a$Y, a$Z, a$W, cbind(a$D, a$D), score_map = map), "names.*unique")
  expect_error(wm_gaussian_match(a$Y, a$Z, a$W, a$D, score_map = function(theta) a$S), "return")
  expect_error(wm_gaussian_match(a$Y, a$Z, a$W, a$D, score_map = function(theta) list(scores0 = matrix(a$S), extra = 1)), "return")
  expect_error(wm_gaussian_match(a$Y, a$Z, a$W, a$D, score_map = map, qr_tol = 0), "qr_tol")
  expect_error(wm_gaussian_match(a$Y, a$Z, a$W, a$D, score_map = map, offset1 = 1), "offset1")
  expect_error(wm_gaussian_match(a$Y, rep(1L, length(a$Y)), a$W, a$D, score_map = map), "both arms")
})
