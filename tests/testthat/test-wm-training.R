wm_training_fixture <- function(estimand = "PATE") {
  x <- c(0, 2, 0.2, 1.8)
  wm_match(c(2, 5, 1, 2), c(1, 1, 0, 0), c(1, 3, 2, 4),
           matrix(x, ncol = 1), M = 1L, estimand = estimand,
           mean0 = 0.1 + 0.2 * x,
           mean1 = if (estimand == "PATE") -0.3 + 0.5 * x else NULL)
}

test_that("independent training uses the evaluation-to-training ratio exactly", {
  D0 <- matrix(c(2, 5, 7, 11), ncol = 1)
  D1 <- matrix(c(-1, 3, 4, 8), ncol = 1)
  U <- matrix(c(4, -2, 1), ncol = 1)
  for (estimand in c("PATE", "PATT")) {
    fit <- wm_training_fixture(estimand)
    adjusted <- wm_training_adjust(fit, U, D0,
                                    if (estimand == "PATE") D1 else NULL)
    b <- if (estimand == "PATE") 53 / 10 else 23 / 4
    expect_equal(unname(adjusted$training_adjustment$sensitivity), b,
                 tolerance = 1e-12)
    expect_equal(adjusted$contributions$training, c(4 * b, -4 * b, 0),
                 tolerance = 1e-12)
    expect_equal(adjusted$training_adjustment$training_root_n_variance, 8 * b^2,
                 tolerance = 1e-12)
    expect_equal(adjusted$root_n_variance, fit$root_n_variance + 8 * b^2,
                 tolerance = 1e-12)
    expect_identical(adjusted$contributions$row, fit$contributions$row)
    expect_identical(adjusted$contributions$nuisance, fit$contributions$nuisance)
    expect_identical(adjusted$contributions$actual, fit$contributions$actual)
    expect_identical(adjusted$estimate, fit$estimate)
    expect_identical(adjusted$graph, fit$graph)
    expect_identical(adjusted$n, fit$n)
    expect_equal(adjusted$training_adjustment$evaluation_n, 4)
    expect_equal(adjusted$training_adjustment$training_n, 3)
    expect_equal(adjusted$training_adjustment$evaluation_training_ratio, 4 / 3)
    expect_equal(adjusted$variance, adjusted$root_n_variance / 4)
    expect_null(adjusted$weight_adjustment)
    expect_null(adjusted$prediction_adjustment)
    boot <- wm_bootstrap(adjusted, B = 2L, seed = 17L)
    expect_equal(boot$conditional_root_n_variance, adjusted$root_n_variance)
    expect_equal(boot$conditional_variance, adjusted$variance)
  }
})

test_that("training centering and sample sizes do not mix with evaluation rows", {
  fit <- wm_training_fixture()
  D0 <- cbind(a = c(2, 5, 7, 11), b = c(1, 3, -2, 4))
  D1 <- cbind(a = c(-1, 3, 4, 8), b = c(2, -3, 1, 0))
  U <- cbind(a = c(4, -2, 1), b = c(-1, 3, 2))
  adjusted <- wm_training_adjust(fit, U, D0, D1)
  b <- adjusted$training_adjustment$sensitivity
  projection <- as.vector(U %*% b)
  centered <- projection - mean(projection)
  expected <- (4 / 3) * mean(centered^2)
  expect_equal(adjusted$training_adjustment$training_root_n_variance, expected)
  expect_equal(adjusted$contributions$training, (4 / 3) * centered)
  expect_equal(mean(adjusted$contributions$training), 0, tolerance = 1e-12)
  expect_equal(adjusted$root_n_variance,
               mean(fit$contributions$row^2) + expected)
  shifted <- wm_training_adjust(fit, sweep(U, 2, c(100, -50), "+"), D0, D1)
  expect_equal(shifted$contributions$training, adjusted$contributions$training,
               tolerance = 1e-10)
  doubled <- wm_training_adjust(fit, rbind(U, U), D0, D1)
  expect_equal(doubled$training_adjustment$training_n, 6)
  expect_equal(doubled$training_adjustment$evaluation_training_ratio, 2 / 3)
  expect_equal(doubled$training_adjustment$training_root_n_variance,
               expected / 2, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(adjusted$root_n_variance,
                                stats::var(c(fit$contributions$row, projection)))))
})

test_that("independently trained scalar maps and zero derivatives are supported", {
  fit <- wm_training_fixture()
  adjusted <- wm_training_adjust(fit, matrix(c(4, -2, 1), ncol = 1))
  expect_equal(unname(adjusted$training_adjustment$dimensions), c(1L, 1L))
  expect_equal(adjusted$contributions$training, numeric(3))
  expect_equal(adjusted$root_n_variance, fit$root_n_variance, tolerance = 1e-12)
  expect_match(adjusted$training_adjustment$contract$independence,
                "independent samples", fixed = TRUE)
  expect_match(adjusted$training_adjustment$contract$weights,
                "W=w_Z(X)", fixed = TRUE)
  expect_identical(adjusted$training_adjustment$contract$status,
                   "user_supplied_independent_training_assumptions_not_verified_by_arrays")
})

test_that("training adjustment rejects existing corrections and unsupported fits", {
  fit <- wm_training_fixture()
  U <- matrix(c(4, -2, 1), ncol = 1)
  adjusted <- wm_training_adjust(fit, U)
  expect_error(wm_training_adjust(adjusted, U), "double counting")
  prediction <- wm_prediction_adjust(fit, matrix(0, 4, 1),
                                      changing_scores = c(FALSE, FALSE))
  expect_error(wm_training_adjust(prediction, U), "double counting")
  weights <- wm_weight_adjust(fit, matrix(0, 4, 1), matrix(0, 4, 1))
  expect_error(wm_training_adjust(weights, U), "double counting")
  bad <- fit
  bad$method <- "stabilized"
  expect_error(wm_training_adjust(bad, U), "self_normalized")
  bad <- fit
  bad$root_n_variance <- NA_real_
  expect_error(wm_training_adjust(bad, U), "variance = TRUE")
  bad <- fit
  bad$info$corrected <- FALSE
  expect_error(wm_training_adjust(bad, U), "fitted mean predictions")
  bad <- fit
  bad$data <- NULL
  expect_error(wm_training_adjust(bad, U), "retain clean original data")
  bad <- fit
  bad$data$Y[1] <- bad$data$Y[1] + 10
  expect_error(wm_training_adjust(bad, U), "inconsistent")
  expect_error(wm_training_adjust(wm_training_fixture("PATT"), U,
                                   mean_derivative1 = matrix(0, 4, 1)),
               "only mean_derivative0")
})

test_that("training derivative columns and unequal sample row counts are validated", {
  fit <- wm_training_fixture()
  U <- cbind(a = c(4, -2, 1), b = c(-1, 3, 2))
  D <- cbind(a = c(2, 5, 7, 11), b = c(1, 3, -2, 4))
  expect_error(wm_training_adjust(fit, matrix(1, 1, 2)), "m >= 2")
  expect_error(wm_training_adjust(fit, as.vector(U)), "m-by-p")
  expect_error(wm_training_adjust(fit, U, mean_derivative0 = U), "n-by-p")
  expect_error(wm_training_adjust(fit, U, mean_derivative0 = D[, 1, drop = FALSE]),
               "p = 2")
  expect_error(wm_training_adjust(fit, U, mean_derivative0 = D[, 2:1]),
               "parameter order")
  invalid <- U
  invalid[1, 1] <- NA_real_
  expect_error(wm_training_adjust(fit, invalid), "finite")
  colnames(invalid) <- c("a", "a")
  invalid[1, 1] <- 0
  expect_error(wm_training_adjust(fit, invalid), "unique")
  valid <- wm_training_adjust(fit, U, D, -D)
  expect_identical(valid$training_adjustment$parameter_names, colnames(U))
  expect_true("wm_training_adjust" %in% getNamespaceExports("wdsmatch"))
})
