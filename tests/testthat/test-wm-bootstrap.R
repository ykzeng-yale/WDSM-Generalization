wm_bootstrap_fixture <- function(row = c(1, 2, -1, 3),
                                 edge = c(-2, 0.5), estimate = 1.25,
                                 training = NULL) {
  n <- length(row)
  v <- (sum(row^2) + sum(edge^2) + sum(training^2)) / n
  structure(list(n = n, estimate = estimate, root_n_variance = v,
    variance = v / n, se = sqrt(v / n),
    contributions = list(row = row, edge = edge, actual = rep(1e3, n),
                         training = training)),
    class = c("wm_match", "list"))
}

test_that("Gaussian draws use the uncentered row and independent edge components", {
  fit <- wm_bootstrap_fixture()
  before <- fit
  set.seed(513)
  expected <- replicate(7L, {
    sum(fit$contributions$row * rnorm(fit$n)) +
      sum(fit$contributions$edge * rnorm(length(fit$contributions$edge)))
  }) / sqrt(fit$n)
  after_expected <- .Random.seed
  set.seed(513)
  out <- wm_bootstrap(fit, B = 7, interval = "none", chunk_size = 2)
  expect_equal(out$root_n_draws, expected, tolerance = 1e-14)
  expect_identical(.Random.seed, after_expected)
  expect_equal(out$draws, fit$estimate + expected / sqrt(fit$n))
  expect_equal(out$conditional_root_n_variance, 19.25 / 4)
  expect_equal(out$conditional_variance, 19.25 / 16)
  expect_equal(out$monte_carlo_root_n_variance, var(expected))
  expect_equal(out$monte_carlo_variance, var(expected) / fit$n)
  expect_null(out$conf.int)
  expect_identical(fit, before)
})

test_that("independent training multipliers follow rows and edges with evaluation scaling", {
  scalar_influence <- c(-2, 0.5, 1.5)
  n <- 4L
  m <- length(scalar_influence)
  training <- (n / m) * scalar_influence
  fit <- wm_bootstrap_fixture(training = training)
  before <- fit
  set.seed(815)
  expected <- replicate(9L, {
    sum(fit$contributions$row * rnorm(n)) +
      sum(fit$contributions$edge * rnorm(length(fit$contributions$edge))) +
      sum(training * rnorm(m))
  }) / sqrt(n)
  after_expected <- .Random.seed
  set.seed(815)
  out <- wm_bootstrap(fit, B = 9, interval = "none", chunk_size = 2)
  expect_equal(out$root_n_draws, expected, tolerance = 1e-14)
  expect_identical(.Random.seed, after_expected)
  expect_equal(out$draws, fit$estimate + expected / sqrt(n))
  base_variance <- (sum(fit$contributions$row^2) +
                      sum(fit$contributions$edge^2)) / n
  expected_variance <- base_variance + (n / m) * mean(scalar_influence^2)
  expect_equal(out$conditional_root_n_variance, expected_variance)
  expect_equal(out$conditional_variance, expected_variance / n)
  expect_identical(out$n, fit$n)
  expect_identical(fit, before)
  small <- wm_bootstrap(fit, B = 9, seed = 815, chunk_size = 1)
  large <- wm_bootstrap(fit, B = 9, seed = 815, chunk_size = 1000)
  expect_equal(small$root_n_draws, large$root_n_draws, tolerance = 1e-14)
})

test_that("absent, NULL and empty training vectors preserve the original stream", {
  absent <- wm_bootstrap_fixture()
  absent$contributions$training <- NULL
  explicit_null <- absent
  explicit_null$contributions["training"] <- list(NULL)
  empty <- absent
  empty$contributions$training <- numeric()
  expected <- wm_bootstrap(absent, B = 5, seed = 203)
  expect_identical(wm_bootstrap(explicit_null, B = 5, seed = 203), expected)
  expect_identical(wm_bootstrap(empty, B = 5, seed = 203), expected)
})

test_that("training contributions are validated against the total fitted variance", {
  fit <- wm_bootstrap_fixture()
  set.seed(93)
  before <- .Random.seed
  for (bad in list(Inf, NA_real_, TRUE, "1", 1 + 1i, matrix(1, 2, 1), list(1))) {
    invalid <- fit
    invalid$contributions$training <- bad
    expect_error(wm_bootstrap(invalid), "training contributions must be")
  }
  omitted_variance <- fit
  omitted_variance$contributions$training <- c(-2, 2)
  expect_error(wm_bootstrap(omitted_variance), "does not agree")
  expect_identical(.Random.seed, before)
})

test_that("normal intervals use exact variance and scales follow the root-n law", {
  fit <- wm_bootstrap_fixture()
  first <- wm_bootstrap(fit, B = 1, seed = 51, conf.level = 0.9)
  many <- wm_bootstrap(fit, B = 19, seed = 52, conf.level = 0.9)
  expected <- fit$estimate + c(-1, 1) * qnorm(0.95) * fit$se
  expect_equal(unname(first$conf.int), expected)
  expect_identical(first$conf.int, many$conf.int)
  expect_true(is.na(first$monte_carlo_variance))
  expect_true(is.na(first$monte_carlo_root_n_variance))

  enlarged <- wm_bootstrap_fixture(fit$contributions$row * 3,
                                  fit$contributions$edge * 3, fit$estimate * 3)
  same_seed <- wm_bootstrap(fit, B = 17, seed = 20)
  transformed <- wm_bootstrap(enlarged, B = 17, seed = 20)
  expect_equal(transformed$root_n_draws, 3 * same_seed$root_n_draws)
  expect_equal(transformed$draws, 3 * same_seed$draws)
  expect_equal(transformed$conditional_variance, 9 * same_seed$conditional_variance)
})

test_that("chunking preserves multiplier order and row-only and zero laws work", {
  fit <- wm_bootstrap_fixture(row = seq(-2, 3, length.out = 37), edge = numeric())
  small <- wm_bootstrap(fit, B = 11, seed = 123, chunk_size = 1)
  large <- wm_bootstrap(fit, B = 11, seed = 123, chunk_size = 1000)
  expect_equal(small$root_n_draws, large$root_n_draws, tolerance = 1e-14)
  expect_equal(small$conditional_root_n_variance, mean(fit$contributions$row^2))
  zero <- wm_bootstrap(wm_bootstrap_fixture(rep(0, 4), numeric(), 2),
                       B = 3, seed = 2)
  expect_identical(zero$root_n_draws, numeric(3))
  expect_equal(zero$draws, rep(2, 3))
  expect_equal(unname(zero$conf.int), c(2, 2))
  expect_equal(zero$monte_carlo_variance, 0)
})

test_that("explicit seeds restore existing and absent caller RNG states", {
  prior_kind <- RNGkind()
  had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  prior_seed <- if (had_seed) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    do.call(RNGkind, as.list(prior_kind))
    if (had_seed) assign(".Random.seed", prior_seed, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG", normal.kind = "Inversion")
  set.seed(182)
  before <- .Random.seed
  kind <- RNGkind()
  out <- wm_bootstrap(wm_bootstrap_fixture(), B = 8, seed = 75)
  expect_identical(.Random.seed, before)
  expect_identical(RNGkind(), kind)
  expect_identical(wm_bootstrap(wm_bootstrap_fixture(), B = 8, seed = 75), out)
  rm(".Random.seed", envir = .GlobalEnv)
  invisible(wm_bootstrap(wm_bootstrap_fixture(), B = 2, seed = 5))
  expect_false(exists(".Random.seed", .GlobalEnv, inherits = FALSE))
  expect_identical(RNGkind(), kind)
})

test_that("Box-Muller explicit seeds fail before disturbing its hidden cache", {
  kind <- RNGkind()
  had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  previous <- if (had_seed) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    do.call(RNGkind, as.list(kind))
    if (had_seed) assign(".Random.seed", previous, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", normal.kind = "Box-Muller")
  set.seed(67)
  expected <- rnorm(2)
  set.seed(67)
  expect_equal(rnorm(1), expected[1])
  before <- .Random.seed
  expect_error(wm_bootstrap(wm_bootstrap_fixture(), seed = 7), "Box-Muller")
  expect_identical(.Random.seed, before)
  expect_equal(rnorm(1), expected[2])
  expect_length(wm_bootstrap(wm_bootstrap_fixture(), B = 2)$root_n_draws, 2)
})

test_that("basic intervals use root-n quantiles with the correct sign", {
  fit <- wm_bootstrap_fixture()
  out <- wm_bootstrap(fit, B = 21, seed = 700, interval = "basic", conf.level = 0.8)
  expected <- fit$estimate - quantile(out$root_n_draws, c(0.9, 0.1)) / sqrt(fit$n)
  expect_equal(unname(out$conf.int), unname(expected))
})

test_that("invalid requests and point-only or corrupted fits fail before drawing", {
  fit <- wm_bootstrap_fixture()
  set.seed(41)
  before <- .Random.seed
  for (bad in list(0, -1, NA, Inf, 1.2, c(1, 2), TRUE, "4", 1 + 1i)) {
    expect_error(wm_bootstrap(fit, B = bad), "positive integer")
    expect_error(wm_bootstrap(fit, chunk_size = bad), "positive integer")
  }
  for (bad in list(-1, NA, Inf, 0.5, c(1, 2), TRUE, "4", 1 + 1i))
    expect_error(wm_bootstrap(fit, seed = bad), "nonnegative integer")
  for (bad in list(0, 1, NA, Inf, c(0.9, 0.95)))
    expect_error(wm_bootstrap(fit, conf.level = bad), "strictly between")
  expect_error(wm_bootstrap(unclass(fit)), "wm_match result")
  point_only <- fit
  point_only$root_n_variance <- NA_real_
  expect_error(wm_bootstrap(point_only), "enable variance estimation")
  corrupt <- fit
  corrupt$contributions$edge[1] <- Inf
  expect_error(wm_bootstrap(corrupt), "finite row and edge")
  corrupt <- fit
  corrupt$root_n_variance <- fit$root_n_variance * 2
  expect_error(wm_bootstrap(corrupt), "does not agree")
  expect_identical(.Random.seed, before)
})
