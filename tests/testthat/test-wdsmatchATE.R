test_that("wdsmatchATE returns correct structure", {
  set.seed(123)
  n <- 200
  X <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  Z <- rbinom(n, 1, plogis(0.5 * X$X1))
  Y <- 1 + X$X1 + Z * 0.5 + rnorm(n)
  w <- runif(n, 0.5, 2)

  fit <- wdsmatchATE(Y, X, Z, weights = w, M = 3, varest = FALSE)

  expect_s3_class(fit, "wdsmatch")
  expect_true(is.numeric(fit$estimate))
  expect_equal(fit$M, 3)
  expect_equal(fit$n, n)
  expect_equal(fit$n.treated + fit$n.control, n)
  expect_true(!is.null(fit$call))
})

test_that("wdsmatchATE returns finite estimate for simple DGP", {
  set.seed(42)
  n <- 500
  X <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  Z <- rbinom(n, 1, 0.5)
  Y <- X$X1 + Z * 1.0 + rnorm(n)
  w <- rep(1, n)

  fit <- wdsmatchATE(Y, X, Z, weights = w, M = 5, varest = FALSE)
  expect_true(is.finite(fit$estimate))
  expect_true(abs(fit$estimate) < 5)
})

test_that("wdsmatchATE with bootstrap produces valid CI", {
  set.seed(99)
  n <- 200
  X <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  Z <- rbinom(n, 1, 0.5)
  Y <- X$X1 + Z * 0.8 + rnorm(n)
  w <- runif(n, 1, 3)

  fit <- wdsmatchATE(Y, X, Z, weights = w, M = 3,
                     varest = TRUE, boots = 50, alpha = 0.05)

  expect_true(!is.na(fit$se))
  expect_true(fit$se > 0)
  expect_true(fit$ci[1] < fit$ci[2])
  expect_true(fit$ci[1] < fit$estimate)
  expect_true(fit$ci[2] > fit$estimate)
  expect_length(fit$boot.estimates, 50L)
})

test_that("wdsmatchATE validates inputs correctly", {
  expect_error(wdsmatchATE(1:5, data.frame(X1 = 1:5), c(0,1,0,1,0)),
               "weights")
  expect_error(wdsmatchATE("a", data.frame(X1 = 1), 0, weights = 1),
               "numeric")
  expect_error(wdsmatchATE(1:3, data.frame(X1 = 1:3), c(0, 2, 1), weights = rep(1, 3)),
               "binary")
  expect_error(wdsmatchATE(1:3, data.frame(X1 = 1:4), c(0,1,0), weights = rep(1, 3)))
})

test_that("wdsmatchATE respects sampling argument", {
  set.seed(7)
  n <- 300
  X <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  Z <- rbinom(n, 1, 0.4)
  Y <- X$X1 + Z + rnorm(n)
  w <- runif(n, 0.5, 3)

  fit_retro <- wdsmatchATE(Y, X, Z, w, M = 3, sampling = "retrospective", varest = FALSE)
  fit_prosp <- wdsmatchATE(Y, X, Z, w, M = 3, sampling = "prospective", varest = FALSE)

  expect_true(is.numeric(fit_retro$estimate))
  expect_true(is.numeric(fit_prosp$estimate))
})

test_that("wdsmatchATE accepts pre-computed scores", {
  set.seed(11)
  n <- 200
  X <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  Z <- rbinom(n, 1, 0.5)
  Y <- X$X1 + Z * 0.5 + rnorm(n)
  w <- rep(1, n)
  ps_pre <- plogis(0.3 * X$X1)
  pg_pre <- cbind(psi0 = X$X1, psi1 = X$X1 + 0.5)

  fit <- wdsmatchATE(Y, X, Z, w, M = 3, ps = ps_pre, pg = pg_pre, varest = FALSE)
  expect_true(is.numeric(fit$estimate))
})

test_that("wdsmatchATE with custom model formulas works", {
  set.seed(22)
  n <- 200
  X <- data.frame(X1 = rnorm(n), X2 = rnorm(n), X3 = rnorm(n))
  Z <- rbinom(n, 1, plogis(X$X1 + X$X1 * X$X2))
  Y <- X$X1 + X$X2 + Z + rnorm(n)
  w <- rep(1, n)

  fit <- wdsmatchATE(Y, X, Z, w, M = 3,
                     model.ps = Z ~ X1 + X2 + X1:X2,
                     model.pg = Y ~ X1 + X2 + X3,
                     varest = FALSE)
  expect_true(is.numeric(fit$estimate))
})

test_that("wdsmatchATE without bias correction differs from with", {
  set.seed(33)
  n <- 300
  X <- data.frame(X1 = rnorm(n))
  Z <- rbinom(n, 1, 0.5)
  Y <- X$X1 + Z + rnorm(n)
  w <- rep(1, n)

  fit_bc <- wdsmatchATE(Y, X, Z, w, M = 3, use.bias.correction = TRUE, varest = FALSE)
  fit_nobc <- wdsmatchATE(Y, X, Z, w, M = 3, use.bias.correction = FALSE, varest = FALSE)
  expect_true(fit_bc$estimate != fit_nobc$estimate)
})

test_that("wdsmatchATE works with M = 1, 3, 5", {
  set.seed(44)
  n <- 300
  X <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  Z <- rbinom(n, 1, 0.5)
  Y <- X$X1 + Z + rnorm(n)
  w <- rep(1, n)

  fit1 <- wdsmatchATE(Y, X, Z, w, M = 1, varest = FALSE)
  fit3 <- wdsmatchATE(Y, X, Z, w, M = 3, varest = FALSE)
  fit5 <- wdsmatchATE(Y, X, Z, w, M = 5, varest = FALSE)

  expect_true(all(sapply(list(fit1, fit3, fit5), function(f) is.finite(f$estimate))))
})

test_that("wdsmatchATE is reproducible with same seed", {
  run_once <- function() {
    set.seed(77)
    n <- 200; X <- data.frame(X1 = rnorm(n))
    Z <- rbinom(n, 1, 0.5); Y <- X$X1 + Z + rnorm(n)
    wdsmatchATE(Y, X, Z, rep(1, n), M = 3, varest = FALSE)$estimate
  }
  expect_equal(run_once(), run_once())
})

test_that("wdsmatchATE handles unequal treatment groups", {
  set.seed(55)
  n <- 300
  X <- data.frame(X1 = rnorm(n))
  Z <- rbinom(n, 1, 0.2)
  Y <- X$X1 + Z * 0.5 + rnorm(n)
  w <- runif(n, 1, 5)

  fit <- wdsmatchATE(Y, X, Z, w, M = 3, varest = FALSE)
  expect_true(is.finite(fit$estimate))
  expect_true(fit$n.treated < fit$n.control)
})

test_that("wdsmatchATE with large survey weights works", {
  set.seed(66)
  n <- 200
  X <- data.frame(X1 = rnorm(n))
  Z <- rbinom(n, 1, 0.5)
  Y <- X$X1 + Z + rnorm(n)
  w <- rexp(n, rate = 0.01)

  fit <- wdsmatchATE(Y, X, Z, w, M = 3, varest = FALSE)
  expect_true(is.finite(fit$estimate))
})
