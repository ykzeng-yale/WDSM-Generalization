test_that("boundary volumes and surface additions have independent geometric fixtures", {
  parts <- wdsmatch:::.wm_overlap_ball_parts
  one <- parts(2, 1, 2, 1)
  expect_equal(exp(unlist(one)), c(I = 1, A = 3, B = 1, U = 5,
                                   s = 1, t = 1, p = 1, q = 1))
  disk <- parts(1, 1, 1, 2)
  intersection <- 2 * pi / 3 - sqrt(3) / 2
  expect_equal(exp(unlist(disk)), c(I = intersection, A = pi - intersection,
    B = pi - intersection, U = 2 * pi - intersection,
    s = 2 * pi / 3, t = 2 * pi / 3, p = 4 * pi / 3, q = 4 * pi / 3))
  sphere <- parts(1, 1, 1, 3)
  expect_equal(exp(sphere$I), 5 * pi / 12)
  expect_equal(exp(sphere$s), pi)
  expect_equal(exp(sphere$p), 3 * pi)
  expect_equal(exp(parts(1, 2, 0, 3)$s), 4 * pi)
  expect_equal(exp(parts(2, 1, 0, 3)$s), 0)
  a <- unlist(parts(2, 1, 1.25, 3))
  b <- unlist(parts(2e100, 1e100, 1.25e100, 3))
  expect_equal(b - a, c(rep(3 * log(1e100), 4), rep(2 * log(1e100), 4)),
               ignore_attr = TRUE)
  expect_error(parts(-1, 1, 1, 2), "nonnegative")
})

test_that("positive boundary polynomials retain each rank case and factorial choice", {
  polynomial <- wdsmatch:::.wm_overlap_log_polynomial
  p <- as.list(log(c(I = 2, A = 3, B = 5, U = 10, s = 7, t = 11, p = 13, q = 17)))
  # Independently expanded: F=30; P_R=105+130; distinct=357+715+442;
  # common boundary leaves one private point on each side, giving 15.
  expect_equal(vapply(1:4, function(f) exp(polynomial(p, 3, 1, f)), numeric(1)),
               c(30, 235, 1514, 15))
  expect_equal(vapply(1:4, function(f) exp(polynomial(p, 3, 2, f)), numeric(1)),
               c(4, 28, 154, 4))
  expect_equal(vapply(1:4, function(f) exp(polynomial(p, 2, 0, f)), numeric(1)),
               c(15, 65, 221, 0))
  zero_parts <- lapply(p, function(x) -Inf)
  zero_parts$U <- 0
  expect_equal(exp(polynomial(zero_parts, 1, 0, 1)), 1)
  expect_equal(exp(polynomial(zero_parts, 1, 0, 2)), 0)
  expect_false(anyNA(wdsmatch:::.wm_overlap_log_add(c(-Inf, log(2)), c(-Inf, log(3)))))
})

test_that("M1 boundary integrands reduce to the full-position integral", {
  for (d in 1:3) {
    set.seed(317 + d)
    normal <- matrix(rnorm(47 * 2 * d), 47, 2 * d)
    reduced <- wdsmatch:::.wm_overlap_angular_terms(normal, 1L, d)
    full <- wdsmatch:::.wm_angular_component(normal, 1L, d, 0L)
    expect_equal(reduced[, 1], full, tolerance = 1e-12)
    expect_equal(reduced[, 2:4, drop = FALSE], matrix(0, 47, 3))
  }
  active <- wdsmatch:::.wm_overlap_active(3L, 1L)
  expect_false(active[2L, "root_root"])
  expect_false(active[1L, "shared_nonroot"])
})

test_that("all boundary and overlap covariance is retained with common Gaussian draws", {
  M <- 2L; d <- 2L; draws <- 79L
  set.seed(819)
  normal <- matrix(rnorm(draws * 6L), draws, byrow = TRUE)
  terms <- wdsmatch:::.wm_overlap_angular_terms(normal, M, d)
  beta <- cbind(rowSums(terms[, c(1, 3, 5, 7)]),
                 rowSums(terms[, c(2, 4, 6, 8)]))
  result <- wm_geometry_overlap(M, d, draws = draws, seed = 819, chunk_size = 13)
  expect_equal(unname(result$term_estimates), matrix(colMeans(terms), M, 4), tolerance = 1e-12)
  expect_equal(unname(result$term_covariance), stats::cov(terms) / draws, tolerance = 1e-12)
  expect_equal(unname(result$beta), colMeans(beta), tolerance = 1e-12)
  expect_equal(unname(result$covariance), stats::cov(beta) / draws, tolerance = 1e-12)
  expect_equal(unname(result$beta_mcse), sqrt(diag(stats::cov(beta)) / draws), tolerance = 1e-12)
  expect_equal(result$alpha_mcse, sqrt(stats::var(rowSums(terms)) / draws), tolerance = 1e-12)
  expect_equal(result$alpha_mcse^2, sum(result$covariance), tolerance = 1e-12)
  expect_gt(abs(sum(result$term_covariance) - sum(diag(result$term_covariance))), 1e-10)
  expect_equal(result$diagnostics$term_nonzero_draws, matrix(colSums(terms > 0), M, 4),
               ignore_attr = TRUE)
  expect_identical(result$method, "boundary_rank_monte_carlo")
  expect_equal(result$normalization$sphere_dimension, c(4, 5, 6, 6))
  expect_true(result$diagnostics$coupled_boundary_terms)
  expect_identical(result$diagnostics$below_jensen, result$alpha < M^2)
})

test_that("exact public reduction, internal 1d integration and chunking are explicit", {
  set.seed(720)
  previous <- .Random.seed
  exact <- wm_geometry_overlap(3, 1, seed = 12)
  expect_equal(unname(exact$beta), c(2, 4, 4.5))
  expect_identical(exact$draws, 0L)
  expect_identical(.Random.seed, previous)
  small <- wm_geometry_overlap(3, 2, draws = 59, seed = 13, chunk_size = 1)
  large <- wm_geometry_overlap(3, 2, draws = 59, seed = 13, chunk_size = 1000)
  expect_identical(.Random.seed, previous)
  expect_equal(small$beta, large$beta, tolerance = 1e-12)
  expect_equal(small$term_covariance, large$term_covariance, tolerance = 1e-12)
  expect_equal(large$chunk_size, 59)
  numerical <- wdsmatch:::.wm_geometry_overlap_dispatch(2, 1, 7, 19, 3, exact_1d = FALSE)
  expect_identical(numerical$method, "boundary_rank_monte_carlo")
  expect_identical(numerical$draws, 7L)
  expect_identical(.Random.seed, previous)
  wm_geometry_overlap(1, 2, draws = 3)
  expect_false(identical(.Random.seed, previous))
})

test_that("invalid and excessive requests fail before consuming random numbers", {
  set.seed(921)
  previous <- .Random.seed
  for (bad in list(0, NA, Inf, 1.5, "2", TRUE, 1 + 1i, c(1, 2))) {
    expect_error(wm_geometry_overlap(bad, 2), "positive integer")
    expect_error(wm_geometry_overlap(1, bad), "positive integer")
    expect_error(wm_geometry_overlap(1, 2, draws = bad), "positive integer")
    expect_error(wm_geometry_overlap(1, 2, chunk_size = bad), "positive integer")
  }
  expect_error(wm_geometry_overlap(1, 2, draws = 1), "draws >= 2")
  expect_error(wm_geometry_overlap(65, 2, draws = 2), "work cap")
  expect_error(wm_geometry_overlap(1, 129, draws = 2), "work cap")
  expect_error(wm_geometry_overlap(3, 5, draws = 1000000), "work cap")
  expect_error(wm_geometry_overlap(1, 2, seed = -1), "nonnegative integer")
  expect_identical(.Random.seed, previous)
})
