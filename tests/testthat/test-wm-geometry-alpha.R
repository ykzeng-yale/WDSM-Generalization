test_that("alpha-only exact geometry has the correct one-dimensional identity", {
  set.seed(31)
  before <- .Random.seed
  for (M in c(1L, 2L, 3L, 5L, 200L)) {
    result <- wm_geometry_alpha(M, 1, draws = 1)
    expect_s3_class(result, "wm_geometry_alpha")
    expect_equal(result$alpha, M^2 + M / 2)
    expect_equal(result$alpha_mcse, 0)
    expect_identical(result$draws, 0L)
    expect_false(result$diagnostics$rng_used)
    expect_null(result$beta)
    expect_false(inherits(result, "wm_geometry"))
  }
  expect_identical(.Random.seed, before)
})

test_that("the reduced alpha integral reproduces the original M1 draws", {
  for (d in c(2L, 3L, 5L)) {
    original <- wm_geometry(1, d, draws = 128L, seed = 492L, chunk_size = 32L)
    reduced <- wm_geometry_alpha(1, d, draws = 128L, seed = 492L, chunk_size = 32L)
    expect_equal(reduced$alpha, original$alpha, tolerance = 1e-12)
    expect_equal(reduced$alpha_mcse, original$alpha_mcse, tolerance = 1e-12)
  }
})

test_that("direct alpha draws preserve RNG and batch invariance and finite bounds", {
  set.seed(953L)
  before <- .Random.seed
  one <- wm_geometry_alpha(3, 5, draws = 257L, seed = 412L, chunk_size = 1L)
  expect_identical(.Random.seed, before)
  many <- wm_geometry_alpha(3, 5, draws = 257L, seed = 412L, chunk_size = 64L)
  expect_identical(.Random.seed, before)
  expect_equal(one$alpha, many$alpha, tolerance = 1e-12)
  expect_equal(one$alpha_mcse, many$alpha_mcse, tolerance = 1e-12)
  expect_gt(one$diagnostics$minimum_draw_value, 0)
  expect_lte(one$diagnostics$maximum_draw_value, one$diagnostics$integrand_upper_bound)
  expect_equal(one$diagnostics$count_terms, 14)
  expect_equal(one$diagnostics$sphere_dimension, 10)
  expect_equal(one$diagnostics$true_mcse_upper_bound,
               one$diagnostics$integrand_upper_bound / (2 * sqrt(257)))
  expect_equal(one$diagnostics$hoeffding_95_half_width,
               one$diagnostics$integrand_upper_bound * sqrt(log(40) / (2 * 257)))
  expect_true("wm_geometry_alpha" %in% getNamespaceExports("WeightedMatching"))
})

test_that("alpha geometry refuses invalid or unbounded requests", {
  expect_error(wm_geometry_alpha(0, 2), "positive integer")
  expect_error(wm_geometry_alpha(3, 2.5), "positive integer")
  expect_error(wm_geometry_alpha(3, 2, draws = 1), "draws >= 2")
  expect_error(wm_geometry_alpha(3, 2, seed = -1), "nonnegative")
  expect_error(wm_geometry_alpha(129, 2), "work cap")
  expect_error(wm_geometry_alpha(3, 129), "work cap")
  expect_error(wm_geometry_alpha(128, 2, draws = 100L), "work cap")
  expect_error(wm_geometry_alpha(3, 2, chunk_size = 0), "positive integer")
})
