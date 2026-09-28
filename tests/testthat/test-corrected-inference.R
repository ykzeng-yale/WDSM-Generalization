test_that("fixed original reuse equals an independent six-unit imputation oracle", {
  A <- c(0, 0, 0, 1, 1, 1)
  Y <- c(7, -1, 6, 5, 10, 9)
  w <- c(1, 2, 4, 3, 5, 6)
  q0 <- 1:6
  q1 <- 5:10
  matches0 <- list(NULL, NULL, NULL, c(1L, 2L), c(2L, 3L), c(1L, 3L))
  matches1 <- list(c(4L, 5L), c(5L, 6L), c(4L, 6L), NULL, NULL, NULL)
  # These fractions are hand-calculated from the donor weights, independently
  # of the reuse implementation. K has survey-weight units for both estimands.
  raw_control <- c(11/5, 11/3, 122/15)
  raw_treated <- c(41/24, 135/88, 124/33)
  for (estimand in c("PATE", "PATT")) {
    expected_K <- c(raw_control, if (estimand == "PATE") raw_treated else rep(0, 3))
    K <- wdsmatch:::wdsm_reuse(A, w, matches0, matches1, estimand)
    expect_equal(K, expected_K, tolerance = 1e-14)
    expect_equal(sum(K[A == 0]), sum(w[A == 1]), tolerance = 1e-14)
    for (bias_correction in c(FALSE, TRUE)) {
      mu0 <- if (bias_correction) q0 else numeric(6)
      mu1 <- if (bias_correction) q1 else numeric(6)
      difference <- numeric(6)
      for (i in which(A == 1)) {
        j <- matches0[[i]]
        difference[i] <- Y[i] - (mu0[i] + sum(w[j] * (Y[j] - mu0[j])) / sum(w[j]))
      }
      if (estimand == "PATE") for (i in which(A == 0)) {
        j <- matches1[[i]]
        difference[i] <- mu1[i] + sum(w[j] * (Y[j] - mu1[j])) / sum(w[j]) - Y[i]
      }
      target <- if (estimand == "PATE") 1:6 else which(A == 1)
      expected_point <- sum(w[target] * difference[target]) / sum(w[target])
      actual <- wdsmatch:::wdsm_replicate(Y, A, w, K, rep(1, 6), mu0, mu1, estimand)
      expect_equal(actual, expected_point, tolerance = 1e-14)
      # A second, nonuniform count vector catches donor-denominator
      # renormalization and omitted residual contributions.
      m <- c(2, 0, 1, 1, 0, 2)
      if (estimand == "PATE") {
        oracle <- (sum(m*w*(mu1-mu0)) +
          sum(m[A == 1]*(w[A == 1]+expected_K[A == 1])*(Y[A == 1]-mu1[A == 1])) -
          sum(m[A == 0]*(w[A == 0]+expected_K[A == 0])*(Y[A == 0]-mu0[A == 0]))) / sum(m*w)
      } else {
        oracle <- (sum(m[A == 1]*w[A == 1]*(Y[A == 1]-mu0[A == 1])) -
          sum(m[A == 0]*expected_K[A == 0]*(Y[A == 0]-mu0[A == 0]))) / sum(m[A == 1]*w[A == 1])
      }
      expect_equal(wdsmatch:::wdsm_replicate(Y, A, w, K, m, mu0, mu1, estimand), oracle,
        tolerance = 1e-14)
      expect_equal(wdsmatch:::wdsm_replicate(Y, A, 100*w, 100*K, m, mu0, mu1, estimand), oracle,
        tolerance = 1e-13)
    }
  }
})

test_that("uncorrected public estimates and every replicate match an independent oracle", {
  set.seed(98111)
  n <- 120L
  X <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  A <- rep(0:1, each = n/2)
  Y <- 0.4*X$X1^2 + A*(1+X$X2) + rnorm(n)
  w <- exp(0.3*X$X1 - 0.2*X$X2)
  ps <- plogis(0.7*X$X1 - 0.3*X$X2)
  pg <- cbind(X$X1 + 0.2*X$X2, X$X1 - 0.4*X$X2)
  B <- 6L
  alpha <- 0.10
  seed <- 14731
  # Independent matching uses the supplied probabilities without a logit,
  # pooled unweighted component moments, and stable original-row tie breaking.
  standardize <- function(x) {
    centered <- sweep(x, 2L, colMeans(x), "-")
    sweep(centered, 2L, sqrt(colMeans(centered^2)), "/")
  }
  D0 <- standardize(cbind(ps, pg[, 1]))
  D1 <- standardize(cbind(ps, pg[, 2]))
  for (estimand in c("PATE", "PATT")) for (M in c(1L, 3L, 5L)) {
    K <- numeric(n)
    direct <- numeric(n)
    target <- if (estimand == "PATE") seq_len(n) else which(A == 1)
    for (i in target) {
      donors <- which(A != A[i])
      D <- if (A[i] == 1) D0 else D1
      distance <- rowSums(sweep(D[donors, , drop = FALSE], 2L, D[i, ], "-")^2)
      j <- donors[order(distance, method = "radix")[seq_len(M)]]
      K[j] <- K[j] + w[i]*w[j]/sum(w[j])
      direct[i] <- (2*A[i]-1) * (Y[i] - sum(w[j]*Y[j])/sum(w[j]))
    }
    expected_point <- sum(w[target]*direct[target])/sum(w[target])
    set.seed(seed)
    counts <- rmultinom(B, n, rep(1/n, n))
    expected_replicates <- apply(counts, 2L, function(m) {
      if (estimand == "PATE") sum(m*(2*A-1)*(w+K)*Y)/sum(m*w) else
        (sum(m[A == 1]*w[A == 1]*Y[A == 1]) - sum(m[A == 0]*K[A == 0]*Y[A == 0])) /
          sum(m[A == 1]*w[A == 1])
    })
    for (sampling in c("retrospective", "prospective")) {
      fun <- if (estimand == "PATE") wdsmatchATE else wdsmatchATT
      set.seed(seed)
      fit <- fun(Y, X, A, w, M = M, ps = ps, pg = pg, sampling = sampling,
        use.bias.correction = FALSE, boots = B, alpha = alpha)
      expect_equal(fit$estimate, expected_point, tolerance = 1e-12)
      expect_equal(fit$boot.estimates, expected_replicates, tolerance = 1e-12)
      expected_se <- sqrt(mean((expected_replicates - mean(expected_replicates))^2))
      expect_equal(fit$se, expected_se, tolerance = 1e-12)
      expect_equal(fit$ci, expected_point + c(-1, 1)*qnorm(1-alpha/2)*expected_se, tolerance = 1e-12)
      expect_length(fit$boot.estimates, B)
      expect_s3_class(fit, "wdsmatch")
    }
  }
})

test_that("supplied scores stay fixed in conditional inference", {
  set.seed(67715)
  n <- 160
  X <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  A <- rep(0:1, each = n/2)
  Y <- rnorm(n) + A
  w <- exp(0.2 * X$X1)
  ps <- plogis(X$X1)
  pg <- cbind(X$X2, X$X1 + X$X2)
  # X affects no fit when all scores are supplied and bias correction is off.
  # The old implementation silently fitted unrelated default formulas here.
  for (fun in list(wdsmatchATE, wdsmatchATT)) {
    set.seed(829)
    fit1 <- fun(Y, X, A, w, ps = ps, pg = pg, use.bias.correction = FALSE, boots = 5)
    set.seed(829)
    fit2 <- fun(Y, -X, A, w, ps = ps, pg = pg, use.bias.correction = FALSE, boots = 5)
    expect_identical(fit1$estimate, fit2$estimate)
    expect_identical(fit1$boot.estimates, fit2$boot.estimates)
    expect_identical(fit1$ci, fit2$ci)
  }
})

test_that("PATT uses only the control prognostic and control correction model", {
  set.seed(51713)
  n <- 160
  X <- data.frame(X1 = rnorm(n), X2 = rnorm(n), X3 = rnorm(n))
  A <- rep(0:1, each = n/2)
  Y <- X$X1 - 0.4*X$X2 + rnorm(n) + A
  w <- exp(0.2*X$X3)
  ps <- plogis(0.2*X$X1 + X$X3)
  pg0 <- X$X1 - 0.4*X$X2
  pg1 <- cbind(pg0, X$X1 + X$X2)
  pg2 <- cbind(pg0, X$X1^2 + X$X3)
  set.seed(829)
  fit1 <- wdsmatchATT(Y, X, A, w, ps = ps, pg = pg1, boots = 5)
  set.seed(829)
  fit2 <- wdsmatchATT(Y, X, A, w, ps = ps, pg = pg2, boots = 5)
  expect_identical(fit1$estimate, fit2$estimate)
  expect_identical(fit1$boot.estimates, fit2$boot.estimates)

  # One treated observation is sufficient for the point estimator when the
  # propensity is supplied and the control models are identified. Fitting an
  # unnecessary treated prognostic model would make this example fail.
  single <- rep(0L, n)
  single[n] <- 1L
  fit <- wdsmatchATT(Y, X, single, w, ps = ps, M = 3, varest = FALSE)
  changed <- Y
  changed[n] <- changed[n] + 3
  shifted <- wdsmatchATT(changed, X, single, w, ps = ps, M = 3, varest = FALSE)
  expect_true(is.finite(fit$estimate))
  expect_equal(shifted$estimate - fit$estimate, 3, tolerance = 1e-12)
})

test_that("an empty-arm bootstrap draw fails with its index instead of becoming zero", {
  Y <- c(7, -1, 6, 5, 10, 9)
  Z <- c(0, 0, 0, 1, 1, 1)
  w <- c(1, 2, 4, 3, 5, 6)
  counts <- cbind(rep(1, 6), c(2, 2, 2, 0, 0, 0))
  for (estimand in c("PATE", "PATT")) {
    K <- c(11/5, 11/3, 122/15, if (estimand == "PATE") c(41/24, 135/88, 124/33) else rep(0, 3))
    pt <- list(Y = Y, Z = Z, weights = w, K = K, estimate = 1, estimand = estimand)
    expect_error(wdsmatch:::wdsm_bootstrap(pt, boots = 2L,
      use.bias.correction = FALSE, bootstrap_counts = counts),
      "bootstrap replicate 2/2 failed: both treatment arms")
  }
})

test_that("dot formulas use exactly X and reproduce explicit formulas in every design", {
  set.seed(511729)
  n <- 180L
  X <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  Z <- rbinom(n, 1, plogis(0.35*X$X1 - 0.45*X$X2))
  Y <- 3*Z + X$X1 + 0.4*X$X2^2 + rnorm(n)
  w <- exp(0.3*X$X1 + 0.25*Z)
  expanded_ps <- wdsmatch:::wdsm_validate_formula(Z ~ ., "Z", X)
  expect_setequal(all.vars(expanded_ps), c("Z", "X1", "X2"))
  expect_false("Y" %in% all.vars(expanded_ps))
  for (fun in list(wdsmatchATE, wdsmatchATT)) for (sampling in c("retrospective", "prospective")) {
    set.seed(14377)
    dots <- fun(Y, X, Z, w, M = 3, model.ps = Z ~ ., model.pg = Y ~ .,
      sampling = sampling, boots = 4)
    set.seed(14377)
    explicit <- fun(Y, X, Z, w, M = 3, model.ps = Z ~ X1 + X2, model.pg = Y ~ X1 + X2,
      sampling = sampling, boots = 4)
    expect_identical(dots$estimate, explicit$estimate)
    expect_identical(dots$boot.estimates, explicit$boot.estimates)
    expect_identical(dots$se, explicit$se)
    expect_identical(dots$ci, explicit$ci)
  }
})

test_that("mixed supplied scores ignore only the corresponding unused formula", {
  set.seed(661729)
  n <- 180L
  X <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  Z <- rbinom(n, 1, plogis(0.35*X$X1 - 0.45*X$X2))
  Y <- 2*Z + X$X1 + 0.6*X$X2 + rnorm(n)
  w <- exp(0.3*X$X1 + 0.25*Z)
  ps <- plogis(0.4*X$X1 - 0.2*X$X2)
  pg <- cbind(X$X1 + 0.7*X$X2, X$X1 - 0.3*X$X2)
  for (fun in list(wdsmatchATE, wdsmatchATT)) for (sampling in c("retrospective", "prospective")) {
    set.seed(5537)
    ps_fit <- fun(Y, X, Z, w, M = 3, ps = ps, model.pg = Y ~ X1 + X2,
      sampling = sampling, boots = 4)
    set.seed(5537)
    ps_unused <- fun(Y, X, Z, w, M = 3, ps = ps,
      model.ps = Z ~ nonexistent_predictor, model.pg = Y ~ X1 + X2,
      sampling = sampling, boots = 4)
    expect_identical(ps_fit$estimate, ps_unused$estimate)
    expect_identical(ps_fit$boot.estimates, ps_unused$boot.estimates)
    expect_identical(ps_fit$ci, ps_unused$ci)
    expect_true(ps_fit$settings$supplied.ps.fixed)
    expect_false(ps_fit$settings$supplied.pg.fixed)
    expect_true(all(vapply(ps_fit$diagnostics$bootstrap.ps,
      function(x) identical(x$source, "supplied_fixed"), logical(1))))

    set.seed(5537)
    pg_fit <- fun(Y, X, Z, w, M = 3, pg = pg, model.ps = Z ~ X1 + X2,
      sampling = sampling, boots = 4)
    set.seed(5537)
    pg_unused <- fun(Y, X, Z, w, M = 3, pg = pg,
      model.pg = Y ~ nonexistent_predictor, model.ps = Z ~ X1 + X2,
      sampling = sampling, boots = 4)
    expect_identical(pg_fit$estimate, pg_unused$estimate)
    expect_identical(pg_fit$boot.estimates, pg_unused$boot.estimates)
    expect_identical(pg_fit$ci, pg_unused$ci)
    expect_false(pg_fit$settings$supplied.ps.fixed)
    expect_true(pg_fit$settings$supplied.pg.fixed)
    expect_true(all(vapply(pg_fit$diagnostics$bootstrap.ps,
      function(x) isTRUE(x$converged), logical(1))))
  }
  for (sampling in c("retrospective", "prospective")) {
    set.seed(77111)
    single <- wdsmatchATT(Y, X, Z, w, pg = pg[, 1, drop = FALSE],
      model.ps = Z ~ X1 + X2, sampling = sampling, boots = 4)
    set.seed(77111)
    double <- wdsmatchATT(Y, X, Z, w, pg = pg,
      model.ps = Z ~ X1 + X2, sampling = sampling, boots = 4)
    expect_identical(single$estimate, double$estimate)
    expect_identical(single$boot.estimates, double$boot.estimates)
    expect_identical(single$ci, double$ci)
  }
})
