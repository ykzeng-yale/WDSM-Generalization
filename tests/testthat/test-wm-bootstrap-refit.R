test_that("multinomial refits obey the exact centered ratio and prediction-imbalance identity", {
  set.seed(28119)
  n <- 48L
  Z <- rep(0:1, each = n / 2)
  X <- matrix(runif(3 * n), n, 3)
  W <- 0.5 + X[, 3] + Z / 2
  q0 <- X[, 1]^2 + X[, 2]
  q1 <- 1 + X[, 1] - X[, 3]
  Y <- ifelse(Z == 1, q1, q0) + rnorm(n, sd = 0.2)
  counts <- rmultinom(12L, n, rep(1 / n, n))
  for (estimand in c("PATE", "PATT")) {
    pate <- estimand == "PATE"
    fit <- wm_match(Y, Z, W, X[, 1:2],
      scores1 = if (pate) X[, c(1, 3)] else NULL, M = 3,
      mean0 = q0, mean1 = if (pate) q1 else NULL, estimand = estimand)
    refit <- function(m) {
      list(mean0 = q0 + X[, 3] * sum((m - 1) * X[, 1]) / n,
           mean1 = q1 + X[, 2] * sum((m - 1) * X[, 3]) / n)
    }
    expected <- vapply(seq_len(ncol(counts)), function(b) {
      m <- counts[, b]
      prediction <- refit(m)
      change <- list(prediction$mean0 - q0, prediction$mean1 - q1)
      remainder <- 0
      for (arm in if (pate) 0:1 else 0L) {
        queries <- which(Z != arm)
        donors <- which(Z == arm)
        # Prediction imbalance, evaluated separately on queries and donors.
        imbalance <- sum(m[queries] * fit$analysis_weights[queries] * change[[arm + 1L]][queries]) -
          sum(m[donors] * fit$loads$incoming[donors, arm + 1L] * change[[arm + 1L]][donors])
        remainder <- remainder + (2 * arm - 1) * imbalance
      }
      denominator <- sum(m * fit$analysis_weights * if (pate) 1 else Z)
      fit$estimate + (sum((m - 1) * fit$contributions$actual * fit$gamma) + remainder) / denominator
    }, numeric(1))
    before <- .Random.seed
    out <- wm_bootstrap_refit(fit, counts, refit, conf.level = 0.9)
    expect_identical(.Random.seed, before)
    expect_equal(out$draws, expected, tolerance = 1e-12)
    expect_equal(out$variance, mean((expected - mean(expected))^2))
    expect_equal(unname(out$conf.int), fit$estimate + c(-1, 1) * qnorm(0.95) * sqrt(out$variance))
    expect_identical(out$variance_divisor, "B")
    expect_equal(out$root_n_draws, sqrt(n) * (expected - fit$estimate))
    frozen <- wm_bootstrap_refit(fit, counts)
    expected_frozen <- vapply(seq_len(ncol(counts)), function(b) {
      m <- counts[, b]
      fit$estimate + sum((m - 1) * fit$contributions$actual * fit$gamma) /
        sum(m * fit$analysis_weights * if (pate) 1 else Z)
    }, numeric(1))
    expect_equal(frozen$draws, expected_frozen, tolerance = 1e-12)
    expect_gt(max(abs(out$draws - frozen$draws)), 1e-7)
  }
})

test_that("raw replication and constant-effect replication retain the original point", {
  n <- 12L
  Z <- rep(0:1, each = 6)
  S <- matrix(seq_len(n) / n, n, 1)
  Y <- S[, 1] + Z
  counts <- matrix(1, n, 2)
  for (estimand in c("PATE", "PATT")) {
    fit <- wm_match(Y, Z, 1 + S[, 1], S, estimand = estimand, variance = FALSE)
    out <- wm_bootstrap_refit(fit, counts)
    expect_equal(out$draws, rep(fit$estimate, 2), tolerance = 1e-12)
    expect_equal(out$variance, 0)
    expect_error(wm_bootstrap_refit(fit, counts, function(m) NULL), "Raw matching")
  }
  fit <- wm_match(Y, Z, 1 + S[, 1], S, mean0 = S[, 1], mean1 = S[, 1] + 1)
  counts[, 2] <- rep(c(0, 2), 6)
  expect_equal(wm_bootstrap_refit(fit, counts)$draws, c(1, 1), tolerance = 1e-12)
})

test_that("failed refits and malformed replication requests are not silently retained", {
  n <- 12L
  Z <- rep(0:1, each = 6)
  S <- matrix(seq_len(n) / n, n, 1)
  fit <- wm_match(S[, 1] + Z, Z, rep(1, n), S, mean0 = S[, 1], mean1 = S[, 1] + 1)
  counts <- matrix(1, n, 3)
  expect_error(wm_bootstrap_refit(fit, counts[, 1, drop = FALSE]), "B >= 2")
  invalid <- counts
  invalid[, 2] <- rep(c(0, 2), each = 6)
  expect_error(wm_bootstrap_refit(fit, invalid), "Both arms")
  invalid <- counts
  invalid[1, 1] <- NA_real_
  expect_error(wm_bootstrap_refit(fit, invalid), "count matrix")
  invalid <- fit
  invalid$method <- "stabilized"
  expect_error(wm_bootstrap_refit(invalid, counts), "self_normalized")
  calls <- 0L
  failed <- function(m) {
    calls <<- calls + 1L
    if (calls == 2L) stop("deliberate fit failure")
    list(mean0 = S[, 1], mean1 = S[, 1] + 1)
  }
  expect_error(wm_bootstrap_refit(fit, counts, failed), "replicate 2/3 failed: deliberate fit failure")
  expect_identical(calls, 2L)
  expect_error(wm_bootstrap_refit(fit, counts, function(m) list(mean0 = NA_real_)), "mean0")
  invalid <- fit
  invalid$loads$incoming <- invalid$loads$incoming + 1
  invalid$data$Y <- invalid$data$Y + seq_len(n)
  expect_error(wm_bootstrap_refit(invalid, counts), "do not reconstruct")
})
