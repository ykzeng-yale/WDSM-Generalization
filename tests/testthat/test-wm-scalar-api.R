# Reuse the original24-row synthetic software fixture; no Monte Carlo study.
wm_scalar_api_fixture <- function() {
  readRDS(testthat::test_path("fixtures", "scalar_replication_24.rds"))
}

test_that("raw scalar logistic options preserve the public matching point", {
  a <- wm_scalar_api_fixture()
  expect_true(a$provenance$synthetic)
  for (target in c("PATE", "PATT")) {
    original <- a$fits[[target]]
    d <- original$fit$data
    fit <- wm_scalar_logistic_match(d$Y, d$Z, original$fit$weights, a$design,
      M = 3, estimand = target, inference = FALSE)
    direct <- wm_match(d$Y, d$Z, original$fit$weights,
      matrix(fit$propensity$probability), M = 3, estimand = target,
      variance = FALSE)
    expect_identical(fit$status, "point_computed")
    expect_identical(fit$fit, direct)
    expect_identical(fit$inference$status, "not_requested")
    expect_identical(fit$numerical_root_bridge, "not_certified")
    expect_identical(fit$fitting_design, a$design)
    expect_error(wm_scalar_logistic_match(d$Y, d$Z, original$fit$weights,
      a$design), "weight_model requires")
    expect_error(wm_scalar_logistic_match(d$Y, d$Z, original$fit$weights,
      a$design, M = 0, inference = FALSE), "M must be")
    expect_error(wm_scalar_logistic_match(d$Y, d$Z, original$fit$weights,
      a$design[, -1, drop = FALSE], inference = FALSE), "leading intercept")
    expect_error(wm_scalar_logistic_match(d$Y, d$Z, original$fit$weights,
      cbind(a$design, a$design[, 2]), inference = FALSE), "full rank")
  }
})

test_that("scalar replication binds rows and distinguishes analytic from draw variance", {
  a <- wm_scalar_api_fixture()
  for (target in c("PATE", "PATT")) {
    fit <- a$fits[[target]]; before <- fit; n <- fit$fit$n
    counts <- matrix(1L, n, 1L)
    r <- wm_scalar_replication(fit, counts = counts, conf.level = .9)
    expect_identical(fit, before)
    expect_identical(r$status, "row_variance_computed")
    expect_equal(r$variance, r$root_n_variance/n)
    expect_equal(r$root_n_draws, 0)
    expect_equal(r$draws, fit$estimate)
    expect_true(is.na(r$empirical_draw_variance_B))
    expect_false(r$comparison$finite_sample_equality_expected)
    expect_error(wm_scalar_replication(fit, counts = counts,
      multipliers = counts), "not both")
    bad <- counts; bad[1:2, 1] <- c(.5, 1.5)
    expect_error(wm_scalar_replication(fit, counts = bad), "nonnegative integers")
    bad <- counts; bad[1, 1] <- -1
    expect_error(wm_scalar_replication(fit, counts = bad), "nonnegative integers")
    bad <- counts; bad[1, 1] <- 0
    expect_error(wm_scalar_replication(fit, counts = bad), "column sums n")
    expect_error(wm_scalar_replication(fit, design = a$design[n:1, ]),
      "differs from the stored fitting-design")
    bad <- fit; bad$fit$loads$incoming[1, 1] <- bad$fit$loads$incoming[1, 1]+1
    expect_error(wm_scalar_replication(bad), "disagrees with the bound inputs")
  }
})

test_that("optional scalar root checks retain scope and explicit failure states", {
  skip_if_not_installed("gmp")
  skip_if_not_installed("Rmpfr")
  a <- wm_scalar_api_fixture(); fit <- a$fits$PATE; before <- fit
  expect_error(wm_scalar_root_certificate(fit, a$design, precision = 63),
    "precision must be")
  expect_error(wm_scalar_root_certificate(fit, a$design, radius = 0),
    "radius must be")
  r <- wm_scalar_root_certificate(fit, a$design)
  expect_identical(fit, before)
  expect_identical(r$status, "root_and_donor_sets_certified")
  expect_false(r$population_assumptions_verified)
  expect_false(r$asymptotic_success_rate_claimed)
  expect_false(r$variance_arithmetic_certified)
  expect_false(r$bootstrap_validity_claimed)
  far <- wm_scalar_root_certificate(fit, a$design, center = c(10, 0, 0),
    radius = 1e-8)
  expect_identical(far$status, "not_certified")
})
