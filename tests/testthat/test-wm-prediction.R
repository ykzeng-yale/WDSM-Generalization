wm_prediction_fixture <- function(estimand = "PATE", dimension = 1L,
                                   M = 1L) {
  x <- c(0, 2, 0.2, 1.8)
  scores <- vapply(seq_len(dimension), function(k) x^k, numeric(4))
  wm_match(c(2, 5, 1, 2), c(1, 1, 0, 0), c(1, 3, 2, 4), scores,
           M = M, estimand = estimand, mean0 = 0.1 + 0.2 * x,
           mean1 = if (estimand == "PATE") -0.3 + 0.5 * x else NULL)
}

wm_prediction_direct <- function(fit, mean0, mean1, weights = fit$weights) {
  y <- fit$data$Y
  z <- fit$data$Z
  query <- if (fit$estimand == "PATE") seq_along(y) else which(z == 1)
  effect <- numeric(length(query))
  for (k in seq_along(query)) {
    i <- query[k]
    donors <- fit$graph$edges$donor[fit$graph$edges$query == i]
    mu <- if (z[i] == 1) mean0 else mean1
    missing <- mu[i] + sum(weights[donors] * (y[donors] - mu[donors])) /
      sum(weights[donors])
    effect[k] <- if (z[i] == 1) y[i] - missing else missing - y[i]
  }
  sum(weights[query] * effect) / sum(weights[query])
}

test_that("prediction imbalance has the PATE and PATT signs and denominators", {
  D0 <- matrix(c(2, 5, 7, 11), ncol = 1)
  D1 <- matrix(c(-1, 3, 4, 8), ncol = 1)
  U <- matrix(c(-1, 2, 0, -1), ncol = 1)
  pate <- wm_prediction_fixture()
  adjusted <- wm_prediction_adjust(pate, U, D0, D1,
                                    changing_scores = c(FALSE, FALSE))
  expect_equal(unname(adjusted$prediction_adjustment$potential0_imbalance),
               -23 / 10, tolerance = 1e-12)
  expect_equal(unname(adjusted$prediction_adjustment$potential1_imbalance),
               30 / 10, tolerance = 1e-12)
  expect_equal(unname(adjusted$prediction_adjustment$sensitivity),
               53 / 10, tolerance = 1e-12)
  patt <- wm_prediction_fixture("PATT")
  adjusted_att <- wm_prediction_adjust(patt, U, D0, changing_scores = FALSE)
  expect_equal(unname(adjusted_att$prediction_adjustment$sensitivity),
               23 / 4, tolerance = 1e-12)
  expect_null(adjusted_att$prediction_adjustment$potential1_imbalance)
  expect_identical(adjusted$estimate, pate$estimate)
  expect_identical(adjusted$graph, pate$graph)
  expect_identical(adjusted$contributions$actual, pate$contributions$actual)
  boot <- wm_bootstrap(adjusted, B = 2L, seed = 91L)
  expect_equal(boot$conditional_root_n_variance, adjusted$root_n_variance)
})

test_that("joint prediction and weight sensitivity matches frozen-graph differences", {
  for (estimand in c("PATE", "PATT")) {
    fit <- wm_prediction_fixture(estimand, dimension = 3L, M = 2L)
    D0 <- cbind(a = c(2, 5, 7, 11), b = c(1, 3, -2, 4))
    D1 <- if (estimand == "PATE")
      cbind(a = c(-1, 3, 4, 8), b = c(2, -3, 1, 0)) else NULL
    dW <- cbind(a = c(0.2, -0.1, 0.3, 0.1),
                 b = c(-0.3, 0.1, 0.2, -0.1))
    U <- cbind(a = fit$contributions$row, b = c(1, -2, 0, 1))
    changing <- if (estimand == "PATE") c(TRUE, TRUE) else TRUE
    adjusted <- wm_prediction_adjust(fit, U, D0, D1, changing, dW)
    numeric_derivative <- vapply(seq_len(2), function(k) {
      h <- 1e-5
      m1plus <- if (is.null(D1)) NULL else fit$predictions$mean1 + h * D1[, k]
      m1minus <- if (is.null(D1)) NULL else fit$predictions$mean1 - h * D1[, k]
      (wm_prediction_direct(fit, fit$predictions$mean0 + h * D0[, k],
                            m1plus, fit$weights + h * dW[, k]) -
         wm_prediction_direct(fit, fit$predictions$mean0 - h * D0[, k],
                              m1minus, fit$weights - h * dW[, k])) / (2 * h)
    }, numeric(1))
    expect_equal(unname(adjusted$prediction_adjustment$sensitivity),
                 numeric_derivative, tolerance = 1e-8)
    metadata <- adjusted$prediction_adjustment
    expect_equal(metadata$sensitivity,
                 metadata$weight_sensitivity + metadata$prediction_sensitivity,
                 tolerance = 1e-12)
    component <- as.vector(U %*% metadata$sensitivity)
    component <- component - mean(component)
    base <- fit$contributions$row - mean(fit$contributions$row)
    expect_equal(adjusted$contributions$nuisance, component, tolerance = 1e-12)
    expect_equal(adjusted$root_n_variance, mean((base + component)^2),
                 tolerance = 1e-12)
    expect_equal(adjusted$root_n_variance,
                 metadata$baseline_root_n_variance + metadata$cross_term +
                   metadata$nuisance_root_n_variance, tolerance = 1e-12)
    expect_gt(abs(metadata$cross_term), 1e-6)
    prediction_component <- as.vector(U %*% metadata$prediction_sensitivity)
    weight_component <- as.vector(U %*% metadata$weight_sensitivity)
    expect_gt(abs(stats::cov(prediction_component, weight_component)), 1e-8)
    expect_equal(adjusted$variance, adjusted$root_n_variance / fit$n)
    expect_null(adjusted$weight_adjustment)
  }
})

test_that("zero prediction derivatives and fixed scalar maps preserve inference", {
  fit <- wm_prediction_fixture()
  U <- matrix(c(-1, 2, 0, -1), ncol = 1)
  adjusted <- wm_prediction_adjust(fit, U, changing_scores = c(FALSE, FALSE))
  expect_equal(unname(adjusted$prediction_adjustment$sensitivity), 0)
  expect_equal(adjusted$root_n_variance, fit$root_n_variance, tolerance = 1e-12)
  expect_equal(adjusted$contributions$nuisance, rep(0, 4))
  expect_identical(unname(adjusted$prediction_adjustment$branch),
                   c("fixed_map", "fixed_map"))
  planar <- wm_prediction_adjust(wm_prediction_fixture(dimension = 2L), U,
                                  changing_scores = c(TRUE, TRUE))
  expect_identical(unname(planar$prediction_adjustment$branch),
                   rep("planar_distributional_transfer", 2))
  expect_match(planar$prediction_adjustment$contract$planar,
                "not fitted/baseline graph pathwise equivalence", fixed = TRUE)
  expect_identical(planar$prediction_adjustment$contract$status,
                   "user_supplied_assumptions_not_verified_by_arrays")
})

test_that("prediction adjustment enforces dimension and inference restrictions", {
  fit <- wm_prediction_fixture()
  U <- matrix(0, 4, 1)
  expect_error(wm_prediction_adjust(fit, U), "Supply changing_scores")
  expect_error(wm_prediction_adjust(fit, U, changing_scores = TRUE), "length 2")
  expect_error(wm_prediction_adjust(fit, U, changing_scores = c(0, 0)), "logical")
  expect_error(wm_prediction_adjust(fit, U, changing_scores = c(FALSE, NA)), "logical")
  expect_error(wm_prediction_adjust(fit, U, changing_scores = c(TRUE, FALSE)),
               "changing scalar")
  expect_error(wm_prediction_adjust(fit, U,
                                     changing_scores = c(potential1 = FALSE,
                                                         potential0 = FALSE)),
               "Named changing_scores")
  planar <- wm_prediction_fixture(dimension = 2L)
  expect_error(wm_prediction_adjust(planar, U, changing_scores = c(TRUE, FALSE),
                                     weight_derivative = U), "fixed known weights")
  fixed <- wm_prediction_adjust(planar, U, changing_scores = c(FALSE, FALSE),
                                 weight_derivative = U)
  expect_equal(unname(fixed$prediction_adjustment$sensitivity), 0)
  mixed <- wm_match(fit$data$Y, fit$data$Z, fit$weights,
                    fit$graph$scores0, scores1 = cbind(x = c(0, 2, 0.2, 1.8),
                                                       y = c(0, 4, 0.04, 3.24),
                                                       t = c(0, 8, 0.008, 5.832)),
                    M = 1L, mean0 = fit$predictions$mean0,
                    mean1 = fit$predictions$mean1)
  mixed_adjusted <- wm_prediction_adjust(mixed, U,
                                         changing_scores = c(FALSE, TRUE),
                                         weight_derivative = U)
  expect_identical(unname(mixed_adjusted$prediction_adjustment$branch),
                   c("fixed_map", "higher_dimensional_graph_transfer"))
  att <- wm_prediction_fixture("PATT")
  expect_error(wm_prediction_adjust(att, U, mean_derivative1 = U,
                                     changing_scores = FALSE), "only mean_derivative0")
  bad <- fit
  bad$method <- "stabilized"
  expect_error(wm_prediction_adjust(bad, U, changing_scores = c(FALSE, FALSE)),
               "self_normalized")
  bad <- fit
  bad$root_n_variance <- NA_real_
  expect_error(wm_prediction_adjust(bad, U, changing_scores = c(FALSE, FALSE)),
               "variance = TRUE")
  bad <- fit
  bad$info$corrected <- FALSE
  expect_error(wm_prediction_adjust(bad, U, changing_scores = c(FALSE, FALSE)),
               "fitted mean predictions")
  adjusted <- wm_prediction_adjust(fit, U, changing_scores = c(FALSE, FALSE))
  expect_error(wm_prediction_adjust(adjusted, U, changing_scores = c(FALSE, FALSE)),
               "double counting")
  weight_adjusted <- wm_weight_adjust(fit, U, U)
  expect_error(wm_prediction_adjust(weight_adjusted, U,
                                     changing_scores = c(FALSE, FALSE)),
               "double counting")
})

test_that("joint derivative matrices are validated and namespace export is present", {
  fit <- wm_prediction_fixture()
  U <- cbind(alpha = c(1, 0, -1, 0), beta = c(0, 1, 0, -1))
  expect_error(wm_prediction_adjust(fit, U, mean_derivative0 = U[, 2:1],
                                     changing_scores = c(FALSE, FALSE)),
               "parameter order")
  expect_error(wm_prediction_adjust(fit, U, mean_derivative0 = matrix(0, 4, 1),
                                     changing_scores = c(FALSE, FALSE)), "p = 2")
  expect_error(wm_prediction_adjust(fit, U[-1, , drop = FALSE],
                                     changing_scores = c(FALSE, FALSE)), "n-by-p")
  repeated <- U
  colnames(repeated) <- c("a", "a")
  expect_error(wm_prediction_adjust(fit, repeated,
                                     changing_scores = c(FALSE, FALSE)), "unique")
  invalid <- U
  invalid[1, 1] <- Inf
  expect_error(wm_prediction_adjust(fit, invalid,
                                     changing_scores = c(FALSE, FALSE)), "finite")
  valid <- wm_prediction_adjust(fit, U, U, -U,
                                 changing_scores = c(FALSE, FALSE))
  expect_identical(valid$prediction_adjustment$parameter_names, colnames(U))
  expect_true("wm_prediction_adjust" %in% getNamespaceExports("WeightedMatching"))
})
