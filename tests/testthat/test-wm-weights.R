wm_weights_fixture <- function(estimand = "PATE", M = 2L, weights = NULL) {
  x <- c(0, 0.25, 0.6, 1, 0.1, 0.3, 0.7, 0.9)
  z <- c(rep(1, 4), rep(0, 4))
  y <- c(-1, -1, 1, 1, 0.2, -0.3, 0.4, 0.1)
  if (is.null(weights)) weights <- c(1, 2, 1.5, 3, 2.5, 1, 4, 2)
  m0 <- 0.1 + 0.2 * x
  m1 <- -0.2 + 0.3 * x
  wm_match(y, z, weights, matrix(x, ncol = 1), M = M, estimand = estimand,
           mean0 = m0, mean1 = if (estimand == "PATE") m1 else NULL)
}

wm_weights_direct <- function(fit, weights) {
  y <- fit$data$Y
  z <- fit$data$Z
  effect <- numeric(length(y))
  query <- if (fit$estimand == "PATE") seq_along(y) else which(z == 1)
  for (i in query) {
    donors <- fit$graph$edges$donor[fit$graph$edges$query == i]
    mean <- if (z[i] == 1) fit$predictions$mean0 else fit$predictions$mean1
    missing <- mean[i] +
      sum(weights[donors] * (y[donors] - mean[donors])) / sum(weights[donors])
    effect[i] <- if (z[i] == 1) y[i] - missing else missing - y[i]
  }
  sum(weights[query] * effect[query]) / sum(weights[query])
}

test_that("finite-graph weight derivatives agree with direct finite differences", {
  for (estimand in c("PATE", "PATT")) {
    for (m in c(1L, 3L)) {
      fit <- wm_weights_fixture(estimand, m)
      derivative <- cbind(a = seq(-0.4, 0.5, length.out = 8),
                           b = c(0.2, -0.1, 0.3, -0.2, 0.1, 0.4, -0.3, 0.2))
      influence <- cbind(a = seq(-1, 1, length.out = 8),
                          b = sin(seq_len(8)))
      adjusted <- wm_weight_adjust(fit, derivative, influence)
      finite_difference <- vapply(seq_len(2), function(k) {
        h <- 1e-5
        (wm_weights_direct(fit, fit$weights + h * derivative[, k]) -
           wm_weights_direct(fit, fit$weights - h * derivative[, k])) / (2 * h)
      }, numeric(1))
      expect_equal(unname(adjusted$weight_adjustment$sensitivity),
                   finite_difference, tolerance = 1e-8)
      expect_equal(adjusted$weight_adjustment$sensitivity,
                   adjusted$weight_adjustment$direct_sensitivity,
                   tolerance = 1e-12)
      expect_identical(adjusted$estimate, fit$estimate)
      expect_identical(adjusted$graph, fit$graph)
      expect_lt(adjusted$weight_adjustment$max_fraction_derivative_sum, 1e-12)
    }
  }
})

test_that("eight-row mark averages reproduce exact derivative benchmarks", {
  x <- c(0, 0.25, 0.6, 1, 0.1, 0.3, 0.7, 0.9)
  z <- c(rep(1, 4), rep(0, 4))
  marks <- as.matrix(expand.grid(rep(list(c(-1, 1)), 4)))
  for (m in 1:4) {
    average <- c(PATE = 0, PATT = 0)
    for (estimand in names(average)) {
      for (row in seq_len(nrow(marks))) {
        y <- c(marks[row, ], rep(0, 4))
        fit <- wm_match(y, z, rep(1, 8), matrix(x, ncol = 1),
                         M = m, estimand = estimand, mean0 = rep(0, 8),
                         mean1 = if (estimand == "PATE") rep(0, 8) else NULL)
        adjusted <- wm_weight_adjust(fit, matrix(y, ncol = 1),
                                      matrix(0, nrow = 8, ncol = 1))
        average[estimand] <- average[estimand] +
          unname(adjusted$weight_adjustment$sensitivity) / nrow(marks)
      }
    }
    expect_equal(unname(average["PATE"]), 1 - 1 / (2 * m) - 1 / 8,
                 tolerance = 1e-12)
    expect_equal(unname(average["PATT"]), 3 / 4, tolerance = 1e-12)
  }
})

test_that("common scaling has zero sensitivity and adjusted rows retain covariance", {
  for (estimand in c("PATE", "PATT")) {
    fit <- wm_weights_fixture(estimand)
    scale <- wm_weight_adjust(fit, matrix(fit$weights, ncol = 1),
                               matrix(seq_len(8), ncol = 1))
    expect_equal(unname(scale$weight_adjustment$sensitivity), 0,
                 tolerance = 1e-12)
    expect_equal(scale$root_n_variance, fit$root_n_variance, tolerance = 1e-12)
    derivative <- matrix(fit$data$Y, ncol = 1)
    adjusted <- wm_weight_adjust(fit, derivative,
                                  matrix(fit$contributions$row, ncol = 1))
    decomposition <- adjusted$weight_adjustment
    expect_gt(abs(decomposition$cross_term), 1e-6)
    expect_equal(adjusted$root_n_variance,
                 decomposition$baseline_root_n_variance +
                   decomposition$cross_term +
                   decomposition$nuisance_root_n_variance, tolerance = 1e-12)
    expect_equal(mean(adjusted$contributions$row), 0, tolerance = 1e-12)
    expect_equal(adjusted$variance, adjusted$root_n_variance / 8)
    scaled_fit <- wm_weights_fixture(estimand, weights = 9 * fit$weights)
    scaled_adjusted <- wm_weight_adjust(scaled_fit, 9 * derivative,
                                         matrix(fit$contributions$row, ncol = 1))
    expect_equal(scaled_adjusted$weight_adjustment$sensitivity,
                 adjusted$weight_adjustment$sensitivity, tolerance = 1e-12)
    expect_equal(scaled_adjusted$root_n_variance, adjusted$root_n_variance,
                 tolerance = 1e-12)
  }
})

test_that("weight adjustment rejects incomplete data and double correction", {
  fit <- wm_weights_fixture()
  derivative <- matrix(seq_len(8), ncol = 1)
  influence <- matrix(rep(c(-1, 1), 4), ncol = 1)
  bad <- fit
  bad$data <- NULL
  expect_error(wm_weight_adjust(bad, derivative, influence), "retain clean original data")
  bad <- fit
  bad$data$Y[1] <- bad$data$Y[1] + 1
  expect_error(wm_weight_adjust(bad, derivative, influence), "inconsistent")
  bad <- fit
  bad$root_n_variance <- NA_real_
  expect_error(wm_weight_adjust(bad, derivative, influence), "variance = TRUE")
  bad <- fit
  bad$method <- "stabilized"
  expect_error(wm_weight_adjust(bad, derivative, influence), "self_normalized")
  expect_error(wm_weight_adjust(fit, derivative, matrix(1, 8, 2)), "p = 1")
  expect_error(wm_weight_adjust(fit, derivative[-1, , drop = FALSE], influence),
               "n-by-p")
  adjusted <- wm_weight_adjust(fit, derivative, influence)
  expect_error(wm_weight_adjust(adjusted, derivative, influence), "double counting")
  named_derivative <- derivative
  named_influence <- influence
  colnames(named_derivative) <- "a"
  colnames(named_influence) <- "b"
  expect_error(wm_weight_adjust(fit, named_derivative, named_influence),
               "parameter order")
})

test_that("intercept logistic weights have analytic coefficients and influence", {
  z <- c(rep(0, 14), rep(1, 6))
  fit <- wm_logistic_weights(z, matrix(1, 20, 1), pi = 0.4,
                              parameter_bound = 5)
  q <- mean(z)
  expect_s3_class(fit, "wm_logistic_weights")
  expect_equal(unname(fit$parameters), stats::qlogis(q), tolerance = 1e-8)
  expect_equal(fit$propensity, rep(q, 20), tolerance = 1e-9)
  expected_weights <- ifelse(z == 1, 0.4 / q, 0.6 / (1 - q))
  expect_equal(fit$weights, expected_weights, tolerance = 1e-8)
  expect_equal(unname(fit$weight_derivative[, 1]),
               expected_weights * (q - z), tolerance = 1e-8)
  expect_equal(unname(fit$nuisance_influence[, 1]),
               (z - q) / (q * (1 - q)), tolerance = 1e-8)
  expect_lt(max(abs(fit$score)), 1e-8)
})

test_that("multivariate logistic fitting solves the score and derivative equations", {
  x <- rep(seq(-1, 1, length.out = 5), each = 5)
  successes <- rep(c(1, 2, 2, 3, 4), each = 5)
  z <- as.integer(rep(1:5, 5) <= successes)
  design <- cbind(intercept = 1, x = x)
  fit <- wm_logistic_weights(z, design, pi = 0.5, parameter_bound = c(5, 5))
  reference <- stats::glm.fit(design, z, family = stats::binomial(),
                              control = stats::glm.control(epsilon = 1e-12,
                                                           maxit = 100))
  expect_equal(unname(fit$parameters), unname(reference$coefficients),
               tolerance = 1e-7)
  expect_lt(max(abs(crossprod(design, z - fit$propensity) / length(z))), 1e-8)
  expect_equal(fit$nuisance_influence %*% fit$empirical_information,
               design * (z - fit$propensity), tolerance = 1e-10)
  weights_at <- function(parameter) {
    q <- stats::plogis(as.vector(design %*% parameter))
    ifelse(z == 1, 0.5 / q, 0.5 / (1 - q))
  }
  for (k in 1:2) {
    direction <- numeric(2)
    direction[k] <- 1e-6
    numeric_derivative <- (weights_at(fit$parameters + direction) -
                             weights_at(fit$parameters - direction)) / 2e-6
    expect_equal(unname(fit$weight_derivative[, k]), numeric_derivative,
                 tolerance = 1e-7)
  }
})

test_that("logistic failures do not silently change the model", {
  z <- rep(0:1, 5)
  design <- cbind(intercept = 1, x = seq(-1, 1, length.out = 10))
  expect_error(wm_logistic_weights(z, design), "parameter_bound")
  expect_error(wm_logistic_weights(z, cbind(design, copy = design[, 2]),
                                    parameter_bound = 5), "full column rank")
  expect_error(wm_logistic_weights(z, design, pi = 1, parameter_bound = 5),
               "strictly between")
  expect_error(wm_logistic_weights(z, design, parameter_bound = 0), "positive")
  separated <- cbind(intercept = 1, x = c(-2, -1, -0.5, 0.5, 1, 2))
  expect_error(wm_logistic_weights(c(0, 0, 0, 1, 1, 1), separated,
                                    parameter_bound = 8),
               "boundary|separat|converge")
  expect_error(wm_logistic_weights(c(rep(0, 14), rep(1, 6)),
                                    matrix(1, 20, 1), parameter_bound = 0.1),
               "boundary")
  expect_error(wm_logistic_weights(z, design, parameter_bound = 5, maxit = 1),
               "did not converge")
})
