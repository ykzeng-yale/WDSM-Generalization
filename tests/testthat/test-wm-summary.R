test_that("supplied-map summaries retain statistics and graph counts", {
  Y <- c(1, 3, 8, 4, 7, 13); Z <- rep(0:1, each = 3)
  S0 <- matrix(c(0, 1, 2, .5, 1.5, 3))
  fit <- wm_match(Y, Z, 1:6, S0, scores1 = cbind(S0, 0:5), M = 2,
    mean0 = c(0, 1, 3, .4, 2, 5), mean1 = c(2, 4, 6, 3, 5, 9))
  before <- fit
  set.seed(82); rng <- .Random.seed
  s <- summary(fit)
  printed <- capture.output(returned <- print(fit))
  expect_identical(returned, before)
  expect_identical(fit, before)
  expect_identical(.Random.seed, rng)
  expect_s3_class(s, "summary_wm")
  expect_identical(s$statistics, c(estimate = fit$estimate,
    raw_estimate = fit$raw_estimate, correction = fit$correction, se = fit$se,
    variance = fit$variance, root_n_variance = fit$root_n_variance))
  expect_identical(s$dimensions, c(arm0 = 1L, arm1 = 2L))
  expect_identical(s$donor_usage$donor_arm, 0:1)
  expect_equal(s$donor_usage$queries, c(3, 3))
  expect_equal(s$donor_usage$selected_pairs, c(6, 6))
  # Treated-direction rows (S0, row-1) select rows4/5 for all three queries.
  expect_equal(s$donor_usage$donors_used, c(3, 2))
  expect_equal(s$donor_usage$maximum_reuse, c(3, 3))
  expect_null(s$conf.int)
  expect_identical(s$contract, fit$info$contract)
  expect_lt(length(printed), 30L)
  expect_true(any(grepl("arm0=1, arm1=2", printed, fixed = TRUE)))
})

test_that("point-only summaries and inherited fitted-map summaries are explicit", {
  x <- seq(0, 1, length.out = 16); z <- rep(0:1, 8)
  fit <- wm_fit(1 + 3*x + 2*z, z, 1+x, x, M = 1, estimand = "PATT")
  s <- summary(fit)
  expect_identical(s$call, fit$fit_call)
  expect_identical(s$dimensions, c(arm0 = 1L))
  expect_identical(s$estimand, "PATT")
  expect_identical(s$donor_usage$selected_pairs, 8L)
  point <- wm_match(1 + 3*x + 2*z, z, 1+x, matrix(x), M = 1,
    variance = FALSE)
  out <- capture.output(print(point))
  expect_identical(summary(point)$inference_status, "point_estimate_only")
  expect_true(is.na(summary(point)$statistics[["se"]]))
  expect_null(summary(point)$conf.int)
  expect_true(any(grepl("Sampling SE/interval unavailable", out, fixed = TRUE)))
  for (bad in list(0, 23, NA_real_, Inf, 1.5, 1+1i, "4", matrix(4))) {
    expect_error(print(point, digits = bad), "digits must be")
  }
})

test_that("fitted wrapper summaries preserve the selected inference and interval", {
  set.seed(7132)
  n <- 320L; x <- matrix(runif(2*n, -1, 1), n, 2)
  D <- cbind(intercept = 1, x1 = x[, 1], x2 = x[, 2])
  z <- rbinom(n, 1, plogis(-.2+.4*x[, 1]-.3*x[, 2]))
  y <- 1+x[, 1]+.4*x[, 2]+z*(.7+.3*x[, 2])+rnorm(n, sd = .4)
  w <- exp(.2*z+.3*x[, 1])
  for (target in c("PATE", "PATT")) for (inference in c("none", "full_x")) {
    pate <- target == "PATE"
    fit <- wm_wdsm_fit(y, z, w, D, D, if (pate) D else NULL,
      estimand = target, inference = inference, conf.level = .9)
    s <- summary(fit); before <- fit
    expect_identical(s$statistics[["estimate"]], fit$estimate)
    expect_identical(s$statistics[["se"]], fit$se)
    expect_identical(s$conf.int, fit$conf.int)
    expect_identical(s$conf.level, .9)
    expect_identical(s$inference_status, fit$inference$status)
    expect_identical(s$contract, fit$inference$contract)
    text <- capture.output(returned <- print(fit))
    expect_identical(returned, before)
    expect_identical(fit, before)
    expect_identical(any(grepl("90% interval:", text, fixed = TRUE)),
      inference == "full_x")
    # The single-model public route delegates exactly to the WDSM pipeline.
    models <- list(list(design = D))
    model <- wm_model_fit(y, z, w, ps_models = models, pg0_models = models,
      pg1_models = if (pate) models else NULL, estimand = target,
      inference = inference, conf.level = .9)
    expect_identical(summary(model)$statistics, s$statistics)
    expect_identical(summary(model)$conf.int, s$conf.int)
    expect_lt(length(capture.output(print(model))), 30L)
  }
})

test_that("potential and independent-block summaries do not mislabel their graphs", {
  y <- c(1, 5, 7, 11); z <- c(0, 0, 1, 1); S <- matrix(c(0, 3, 1, 2))
  f <- wm_potential(y, z, 1:4, S, arm = 1, M = 1,
    mean = rep(0, 4), rho = rep(3.5, 4))
  s <- summary(f)
  expect_identical(s$dimensions, c(arm1 = 1L))
  expect_identical(s$donor_usage$donor_arm, 1L)
  expect_true(is.na(s$statistics[["raw_estimate"]]))
  expect_lt(length(capture.output(print(f))), 30L)
  split <- wm_split_pate(rep(y, 2), rep(z, 2), rep(1:4, 2),
    rbind(S, S), rbind(S, S), rep(0:1, each = 4),
    mean0 = rep(0, 8), mean1 = rep(0, 8), rho0 = rep(1.5, 8),
    rho1 = rep(3.5, 8), M = 1)
  expect_null(summary(split)$dimensions)
  expect_null(summary(split)$donor_usage)
  expect_identical(summary(split)$statistics[["estimate"]], split$estimate)
  expect_lt(length(capture.output(print(split))), 30L)
})
