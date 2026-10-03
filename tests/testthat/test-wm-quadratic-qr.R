# Deterministic regression algebra and one existing finite-model fixture.
# These are software checks, not a repeated-sampling or coverage experiment.
wm_qr_reference_basis <- function(s) {
  b <- cbind(intercept = 1, s)
  for (a in seq_len(ncol(s))) for (bb in a:ncol(s)) {
    b <- cbind(b, s[, a] * s[, bb])
  }
  b
}

test_that("QR quadratic predictions retain the complete near-coincident polynomial space", {
  x <- as.matrix(expand.grid(t = seq(-1, 1, length.out = 5L),
                            r = seq(-.9, .9, length.out = 5L),
                            v = seq(-.8, .8, length.out = 5L)))
  n <- nrow(x)
  z <- rep(0:1, length.out = n)
  w <- exp(.3 * z + .2 * x[, 1] - .1 * x[, 2])
  y <- sin(1.3 * x[, 1]) + x[, 2] * x[, 3] + .6 * z + .02 * sin(seq_len(n))
  s <- cbind(ps1 = x[, 1], ps2 = x[, 1] + 1e-6 * x[, 3], pg = x[, 2])
  # A separate explicit contrast gives a well-conditioned basis of exactly the
  # same polynomial space. It neither calls QR-coordinate nor package-basis code.
  reference <- wm_qr_reference_basis(cbind(s[, 1], s[, 3],
                                           (s[, 2] - s[, 1]) / 1e-6))
  raw <- wm_qr_reference_basis(s)
  result <- list()
  for (arm in 0:1) {
    donors <- which(z == arm)
    beta <- stats::lm.wfit(reference[donors, ], y[donors], w[donors],
                          singular.ok = FALSE)$coefficients
    expected <- as.vector(reference %*% beta)
    out <- wdsmatch:::.wm_quad_qr_arm(y, z, w, s, arm, 1e-10)
    result[[as.character(arm)]] <- out
    expect_equal(out$mean, expected, tolerance = 2e-8)
    expect_identical(out$basis_columns, 10L)
    expect_identical(out$columns_removed, 0L)
    expect_equal(out$ridge, 0)
    expect_identical(out$donor_rows, donors)
    expect_lt(out$score_reconstruction_relative_error, 1e-12)
    expect_lt(out$normalized_correction_moment, 1e-12)
    expect_false(out$numerical_error_rate_verified)
    # Raw coefficient rank rejection is the actual source-level problem.
    expect_error(wdsmatch:::.wm_ns_ols(raw[donors, ], y[donors], w[donors],
                                     "raw", 1e-10), "singular|unidentified")
    rescaled <- wdsmatch:::.wm_quad_qr_arm(y, z, 11 * w, s, arm, 1e-10)
    expect_equal(rescaled$mean, out$mean, tolerance = 1e-12)
  }
  # Check the original donor fractions and target denominators independently.
  for (target in c("PATE", "PATT")) for (M in c(1L, 3L)) {
    fit <- wm_match(y, z, w, s, if (target == "PATE") s else NULL,
      M = M, estimand = target, mean0 = result[["0"]]$mean,
      mean1 = if (target == "PATE") result[["1"]]$mean else NULL, variance = FALSE)
    y0 <- y1 <- y
    for (i in which(z == 1L | target == "PATE")) {
      donor <- fit$graph$neighbors[[i]]
      mu <- result[[as.character(1L - z[i])]]$mean
      value <- mu[i] + sum(w[donor] * (y[donor] - mu[donor])) / sum(w[donor])
      if (z[i] == 1L) y0[i] <- value else y1[i] <- value
    }
    target_weight <- w * if (target == "PATT") z else 1
    expect_equal(fit$estimate, sum(target_weight * (y1 - y0)) / sum(target_weight),
                 tolerance = 1e-12)
  }
  duplicate <- cbind(x[, 1], x[, 1], x[, 2])
  expect_error(wdsmatch:::.wm_quad_qr_arm(y, z, w, duplicate, 0L, 1e-10),
               "singular or unresolved")
})

test_that("common API preserves the original quadratic point and graph", {
  # Reuse the pre-existing model-fit fixture and its actual model descriptors.
  set.seed(29317)
  n <- 320L
  x <- cbind(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
  D <- cbind(intercept = 1, x)
  G <- cbind(D, product = x[, 1] * x[, 2])
  Z <- rbinom(n, 1, plogis(-.2 + .5 * x[, 1] - .3 * x[, 2]))
  Y <- 1 + x[, 1] + .4 * x[, 2] + Z * (.8 + x[, 2]) + rnorm(n, sd = .4)
  W <- exp(.5 * Z + .4 * x[, 1] - .2 * x[, 2])
  arguments <- list(Y = Y, Z = Z, weights = W,
    ps_models = list(joint = list(design = D, weighting = "probability"),
      second = list(design = D[, c("intercept", "x2"), drop = FALSE], weighting = "unit")),
    pg0_models = list(product = list(design = G)),
    pg1_models = list(linear = list(design = D)),
    M = 3L, estimand = "PATE", inference = "none")
  original <- do.call(wm_model_fit, arguments)
  out <- do.call(wm_model_fit, c(arguments, list(correction = list(method = "quadratic_qr"))))
  expect_identical(out$fit$graph$neighbors, original$fit$graph$neighbors)
  expect_equal(out$fit$graph$scores0, original$fit$graph$scores0, tolerance = 1e-13)
  expect_equal(out$fit$graph$scores1, original$fit$graph$scores1, tolerance = 1e-13)
  expect_equal(out$fit$loads, original$fit$loads, tolerance = 1e-13)
  expect_equal(out$estimate, original$estimate, tolerance = 1e-10)
  expect_equal(out$fit$predictions, original$fit$predictions, tolerance = 1e-9)
  expect_identical(out$fit$weights, W)
  expect_true(out$nuisance$score_only)
  expect_length(out$nuisance$blocks$correction, 0L)
  expect_identical(unname(out$allocation$basis_columns), c(10L, 10L))
  expect_false(out$inference$available)
  expect_true(all(is.na(out$conf.int)))
  expect_false(out$bootstrap_requested)
  after <- wm_model_inference(out, covariance_scope = "full_x")
  expect_identical(after$fit, out$fit)
  expect_identical(after$nuisance, out$nuisance)
  expect_identical(after$inference$status, "model_correction_law_unavailable")
  expect_false(after$inference_handoff$dispatched)
  expect_match(after$inference$unavailable_reason, "projection contribution")
  expect_match(after$inference_handoff$correction_contract, "quadratic")

  # PATT uses only its control correction and original treated denominator.
  arguments$estimand <- "PATT"
  arguments$pg1_models <- NULL
  patt <- do.call(wm_model_fit, c(arguments, list(correction = list(method = "quadratic_qr"))))
  expect_identical(names(patt$correction_fit$arms), "0")
  expect_null(patt$nuisance$scores1)
  expect_identical(patt$fit$weights, W)
  expect_identical(wm_model_inference(patt, "patt")$inference$status,
                   "model_correction_law_unavailable")

  # Rejected options fail before any model fit or matching operation.
  local_mocked_bindings(.wm_wdsm_fit_stack = function(...) stop("UNEXPECTED_FIT"),
                        .package = "wdsmatch")
  arguments$inference <- "full_x"
  expect_error(do.call(wm_model_fit, c(arguments,
    list(correction = list(method = "quadratic_qr")))), "requires inference='none'")
  arguments$inference <- "none"
  expect_error(do.call(wm_model_fit, c(arguments,
    list(correction = list(method = "quadratic_qr", ridge = 1)))), "accepts only")
  arguments$max_stack_elements <- 1
  expect_error(do.call(wm_model_fit, c(arguments,
    list(correction = list(method = "quadratic_qr")))), "max_stack_elements")
})
