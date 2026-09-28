test_that("one-dimensional overlap constants use the exact all-M formula", {
  set.seed(519)
  before <- .Random.seed
  for (M in c(1L, 2L, 3L, 9L, 10000L)) {
    out <- wm_geometry(M, 1L, draws = 2, seed = 57)
    expected <- c(if (M > 1) 2 * seq_len(M - 1L) else numeric(), 1.5 * M)
    expect_equal(unname(out$beta), expected)
    expect_equal(out$alpha, M^2 + M / 2)
    expect_equal(sum(out$beta), out$alpha)
    expect_equal(out$beta_mcse, setNames(numeric(M), as.character(0:(M - 1L))))
    expect_identical(out$alpha_mcse, 0)
    expect_identical(out$draws, 0L)
    expect_null(out$covariance)
    expect_identical(out$method, "exact_1d_closed_form")
  }
  expect_identical(.Random.seed, before)
})

test_that("ball union volumes agree with independent interval, disk and sphere formulas", {
  union <- wdsmatch:::.wm_ball_union_volume
  expect_equal(union(c(2, 2, 2, 0), c(1, 1, 1, 0), c(0, 2, 4, 0), 1),
               c(4, 5, 6, 0))
  expect_equal(union(1, 1, 0, 2), pi)
  expect_equal(union(1, 1, 2, 2), 2 * pi)
  expect_equal(union(1, 1, 1, 2), 4 * pi / 3 + sqrt(3) / 2)
  expect_equal(union(1, 1, 1, 3), 9 * pi / 4)
  # This unequal-radius lens includes a cap larger than a half-ball.
  r <- 2; s <- 1; distance <- 1.25
  intersection <- r^2 * acos((distance^2 + r^2 - s^2) / (2 * distance * r)) +
    s^2 * acos((distance^2 + s^2 - r^2) / (2 * distance * s)) -
    sqrt((-distance + r + s) * (distance + r - s) *
           (distance - r + s) * (distance + r + s)) / 2
  expected <- pi * (r^2 + s^2) - intersection
  expect_equal(union(r, s, distance, 2), expected)
  expect_equal(union(s, r, distance, 2), expected)
  expect_equal(wdsmatch:::.wm_log_ball_union(r * 1e100, s * 1e100, distance * 1e100, 3),
               wdsmatch:::.wm_log_ball_union(r, s, distance, 3) + 3 * log(1e100))
  expect_error(union(-1, 1, 1, 2), "nonnegative")
})

test_that("angular factors and private-donor exclusions have exact fixtures", {
  component <- wdsmatch:::.wm_angular_component
  normals1 <- rbind(c(1, 0), c(1, 1), c(1, -1))
  expect_equal(component(normals1, 1L, 1L, 0L), c(pi / 2, pi, pi / 4))
  normals2 <- rbind(c(1, 0, 0, 0), c(1, 0, 1, 0), c(1, 0, -1, 0))
  expect_equal(component(normals2, 1L, 2L, 0L), c(1, 4, 1))
  # For M=2 and no shared donors, [-1,1,-2,2]/sqrt(10) has union
  # length 4/sqrt(10), sphere area 2*pi^2 and Gamma(4)=6.
  private <- rbind(c(-1, 1, -2, 2), c(-1, 1, 0, 2))
  expect_equal(component(private, 2L, 1L, 0L), c(75 * pi^2 / 16, 0))
})

test_that("coupled MC moments include within-draw component covariance", {
  M <- 2L; d <- 2L; draws <- 73L
  set.seed(404)
  normals <- matrix(rnorm(draws * 2 * M * d), nrow = draws, byrow = TRUE)
  values <- cbind(wdsmatch:::.wm_angular_component(normals, M, d, 0L),
                  wdsmatch:::.wm_angular_component(normals, M, d, 1L))
  covariance <- stats::cov(values) / draws
  out <- wm_geometry(M, d, draws = draws, seed = 404, chunk_size = 17L)
  expect_equal(unname(out$beta), colMeans(values), tolerance = 1e-12)
  expect_equal(unname(out$covariance), covariance, tolerance = 1e-12)
  expect_equal(unname(out$beta_mcse), sqrt(diag(covariance)), tolerance = 1e-12)
  expect_equal(out$alpha_mcse, sqrt(var(rowSums(values)) / draws), tolerance = 1e-12)
  expect_equal(out$alpha_mcse^2, sum(out$covariance), tolerance = 1e-12)
  expect_equal(out$alpha, sum(out$beta))
  expect_gt(abs(sum(covariance) - sum(diag(covariance))), 1e-10)
  expect_identical(out$method, "uniform_sphere_monte_carlo")
  expect_false(grepl("exact", out$precision_status))
  expect_identical(out$diagnostics$below_jensen, out$alpha < M^2)
  expect_equal(out$normalization$sphere_dimension, c(8, 6))
  expected_area <- 2 * pi^(out$normalization$sphere_dimension / 2) /
    gamma(out$normalization$sphere_dimension / 2)
  expect_equal(exp(out$normalization$log_sphere_area), expected_area)
})

test_that("chunking preserves coupled draws and explicit seeds restore RNG state", {
  set.seed(771)
  before <- .Random.seed
  small <- wm_geometry(2, 2, draws = 53, seed = 36, chunk_size = 1)
  expect_identical(.Random.seed, before)
  large <- wm_geometry(2, 2, draws = 53, seed = 36, chunk_size = 1000)
  expect_identical(.Random.seed, before)
  expect_equal(small$beta, large$beta, tolerance = 1e-12)
  expect_equal(small$covariance, large$covariance, tolerance = 1e-12)
  expect_equal(small$alpha_mcse, large$alpha_mcse, tolerance = 1e-12)
  expect_equal(large$chunk_size, 53)
  no_seed <- wm_geometry(1, 2, draws = 5, seed = NULL)
  expect_false(identical(.Random.seed, before))
  expect_true(no_seed$diagnostics$rng_used)
})

test_that("invalid and oversized integration requests fail before drawing", {
  set.seed(81)
  before <- .Random.seed
  for (bad in list(0, -1, NA, Inf, 1.5, TRUE, "2", c(1, 2), 1 + 1i)) {
    expect_error(wm_geometry(bad, 2), "positive integer")
    expect_error(wm_geometry(1, bad), "positive integer")
    expect_error(wm_geometry(1, 2, draws = bad), "positive integer")
    expect_error(wm_geometry(1, 2, chunk_size = bad), "positive integer")
  }
  expect_error(wm_geometry(2, 2, draws = 1), "draws >= 2")
  expect_error(wm_geometry(129, 2, draws = 2), "work cap")
  expect_error(wm_geometry(10, 10, draws = 1000000), "work cap")
  expect_error(wm_geometry(1000001, 1), "allocation cap")
  expect_error(wm_geometry(1, 2, seed = -1), "nonnegative integer")
  expect_identical(.Random.seed, before)
})
