test_that("wdsmatchATT returns correct structure", {
  set.seed(123)
  n <- 200
  X <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  Z <- rbinom(n, 1, plogis(0.5 * X$X1))
  Y <- 1 + X$X1 + Z * 0.5 + rnorm(n)
  w <- runif(n, 0.5, 2)

  fit <- wdsmatchATT(Y, X, Z, weights = w, M = 3, varest = FALSE)

  expect_s3_class(fit, "wdsmatch")
  expect_true(is.numeric(fit$estimate))
  expect_equal(fit$M, 3)
  expect_equal(fit$n, n)
})

test_that("wdsmatchATT returns finite estimate for simple DGP", {
  set.seed(42)
  n <- 500
  X <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  Z <- rbinom(n, 1, 0.5)
  Y <- X$X1 + Z * 1.0 + rnorm(n)
  w <- rep(1, n)

  fit <- wdsmatchATT(Y, X, Z, weights = w, M = 5, varest = FALSE)
  expect_true(is.finite(fit$estimate))
  expect_true(abs(fit$estimate) < 5)
})

test_that("wdsmatchATT with bootstrap produces valid CI", {
  set.seed(99)
  n <- 200
  X <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  Z <- rbinom(n, 1, 0.5)
  Y <- X$X1 + Z * 0.8 + rnorm(n)
  w <- runif(n, 1, 3)

  fit <- wdsmatchATT(Y, X, Z, weights = w, M = 3,
                     varest = TRUE, boots = 50, alpha = 0.05)

  expect_true(!is.na(fit$se))
  expect_true(fit$se > 0)
  expect_true(fit$ci[1] < fit$ci[2])
})

test_that("wdsmatchATT respects sampling argument", {
  set.seed(7)
  n <- 300
  X <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  Z <- rbinom(n, 1, 0.4)
  Y <- X$X1 + Z + rnorm(n)
  w <- runif(n, 0.5, 3)

  fit_retro <- wdsmatchATT(Y, X, Z, w, M = 3, sampling = "retrospective", varest = FALSE)
  fit_prosp <- wdsmatchATT(Y, X, Z, w, M = 3, sampling = "prospective", varest = FALSE)

  expect_true(is.numeric(fit_retro$estimate))
  expect_true(is.numeric(fit_prosp$estimate))
})

test_that("wdsmatchATT works with M = 1, 3, 5", {
  set.seed(55)
  n <- 300
  X <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  Z <- rbinom(n, 1, 0.5)
  Y <- X$X1 + Z + rnorm(n)
  w <- rep(1, n)

  fit1 <- wdsmatchATT(Y, X, Z, w, M = 1, varest = FALSE)
  fit5 <- wdsmatchATT(Y, X, Z, w, M = 5, varest = FALSE)

  expect_true(all(sapply(list(fit1, fit5), function(f) is.finite(f$estimate))))
})

test_that("wdsmatchATT validates missing weights", {
  expect_error(wdsmatchATT(1:5, data.frame(X1 = 1:5), c(0,1,0,1,0)),
               "weights")
})

test_that("wdsmatchATT is reproducible with same seed", {
  run_once <- function() {
    set.seed(88)
    n <- 200; X <- data.frame(X1 = rnorm(n))
    Z <- rbinom(n, 1, 0.5); Y <- X$X1 + Z + rnorm(n)
    wdsmatchATT(Y, X, Z, rep(1, n), M = 3, varest = FALSE)$estimate
  }
  expect_equal(run_once(), run_once())
})

test_that("wdsmatchATT handles unequal treatment groups", {
  set.seed(66)
  n <- 300
  X <- data.frame(X1 = rnorm(n))
  Z <- rbinom(n, 1, 0.2)
  Y <- X$X1 + Z * 0.5 + rnorm(n)
  w <- runif(n, 1, 5)

  fit <- wdsmatchATT(Y, X, Z, w, M = 3, varest = FALSE)
  expect_true(is.finite(fit$estimate))
})

test_that("wdsmatchATT without bias correction differs from with", {
  set.seed(77)
  n <- 300
  X <- data.frame(X1 = rnorm(n))
  Z <- rbinom(n, 1, 0.5)
  Y <- X$X1 + Z + rnorm(n)
  w <- rep(1, n)

  fit_bc <- wdsmatchATT(Y, X, Z, w, M = 3, use.bias.correction = TRUE, varest = FALSE)
  fit_nobc <- wdsmatchATT(Y, X, Z, w, M = 3, use.bias.correction = FALSE, varest = FALSE)
  expect_true(fit_bc$estimate != fit_nobc$estimate)
})
