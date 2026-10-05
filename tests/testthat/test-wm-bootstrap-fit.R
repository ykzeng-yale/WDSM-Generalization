# Bounded implementation tests. These do not establish population premises,
# asymptotic graph transfer or repeated-sampling coverage.
wm_bootstrap_fit_fixture <- function() {
  set.seed(29317)
  n <- 320L
  x <- cbind(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
  D <- cbind(intercept = 1, x)
  Z <- rbinom(n, 1, plogis(-.2 + .5*x[, 1] - .3*x[, 2]))
  Y <- 1 + x[, 1] + .4*x[, 2] + Z*(.8 + x[, 2]) + rnorm(n, sd = .4)
  W <- exp(.5*Z + .4*x[, 1] - .2*x[, 2])
  list(n = n, D = D, Z = Z, Y = Y, W = W,
    counts = cbind(rep(1, n), rep(c(0, 2), n/2), rep(c(2, 0), n/2)))
}

# Literal full same-count prediction recipe, without candidate stack/basis
# helpers: weighted PS, multiplicity-only PG, count pooling, weighted quadratic.
wm_bootstrap_literal_prediction <- function(a, m, target, ps_weighting = "probability") {
  n <- a$n; D <- a$D; Z <- a$Z; W <- a$W/max(a$W)
  beta <- stats::glm.fit(D, Z,
    weights = m * if (ps_weighting == "probability") W else 1,
    start = rep(0, ncol(D)), family = stats::quasibinomial(),
    control = stats::glm.control(epsilon = 1e-12, maxit = 100L))$coefficients
  raw <- cbind(ps = plogis(drop(D %*% beta)))
  arms <- if (target == "PATE") 0:1 else 0L
  for (arm in arms) {
    pg <- stats::lm.wfit(D, a$Y, w = m*(Z == arm), tol = 1e-10,
      singular.ok = FALSE)$coefficients
    raw <- cbind(raw, drop(D %*% pg))
  }
  prob <- m/sum(m)
  center <- sweep(raw, 2L, colSums(raw*prob), "-")
  standardized <- sweep(center, 2L, sqrt(colSums(center^2*prob)), "/")
  means <- list()
  for (k in seq_along(arms)) {
    score <- standardized[, c(1L, k + 1L)]
    B <- cbind(1, score[, 1], score[, 2], score[, 1]^2,
      score[, 1]*score[, 2], score[, 2]^2)
    bc <- stats::lm.wfit(B, a$Y, w = m*W*(Z == arms[k]),
      tol = 1e-10, singular.ok = FALSE)$coefficients
    means[[paste0("mean", arms[k])]] <- drop(B %*% bc)
  }
  means
}

# Literal original fixed-reuse ratio, independently reconstructing the incoming
# loads from each stored donor set rather than using candidate load helpers.
wm_bootstrap_literal_ratio <- function(fit, m, predictions) {
  n <- fit$n; Z <- fit$data$Z; Y <- fit$data$Y
  W <- fit$weights/max(fit$weights); K <- numeric(n)
  pate <- fit$estimand == "PATE"
  for (i in if (pate) seq_len(n) else which(Z == 1L)) {
    j <- fit$graph$neighbors[[i]]
    K[j] <- K[j] + W[i]*W[j]/sum(W[j])
  }
  q0 <- predictions$mean0
  if (pate) {
    q1 <- predictions$mean1
    (sum(m*W*(q1-q0)) + sum(m*(2*Z-1)*(W+K)*
      (Y-ifelse(Z == 1, q1, q0))))/sum(m*W)
  } else {
    (sum(m*Z*W*(Y-q0)) - sum(m*(1-Z)*K*(Y-q0)))/sum(m*Z*W)
  }
}

test_that("fitted wrappers use bound complete contribution inference", {
  a <- wm_bootstrap_fit_fixture()
  for (target in c("PATE", "PATT")) {
    pate <- target == "PATE"
    f <- wm_wdsm_fit(a$Y, a$Z, a$W, a$D, a$D,
      if (pate) a$D else NULL, estimand = target, inference = "full_x")
    before <- f
    b <- wm_bootstrap(f, B = 3L, seed = 91L)
    h <- wm_model_inference(f, if (pate) "full_x" else "patt",
      transport0 = if (pate) NULL else list(mode = "zero", basis = "current_centering"))
    direct <- wm_bootstrap(h$inference, B = 3L, seed = 91L)
    expect_identical(f, before)
    expect_identical(b$draws, direct$draws)
    expect_equal(b$conditional_variance, f$variance, tolerance = 1e-12)
    expect_equal(b$conditional_root_n_variance, f$root_n_variance, tolerance = 1e-12)
    expect_identical(b$n, a$n)
    expect_false(b$models_refitted)
    expect_false(b$original_refit_bootstrap)
    expect_false(b$population_assumptions_verified)
    expect_true(b$original_fit_preserved)
    counted <- wm_bootstrap(f, counts = a$counts)
    expect_identical(counted$B, ncol(a$counts))
    expect_equal(counted$root_n_draws, drop(crossprod(a$counts-1,
      h$inference$augmented_rows))/sqrt(a$n), tolerance = 1e-12)
  }
})

test_that("default fitted-wrapper refits reproduce the full same-count recipe", {
  a <- wm_bootstrap_fit_fixture()
  for (target in c("PATE", "PATT")) for (weighting in c("probability", "unit")) {
    pate <- target == "PATE"
    f <- wm_wdsm_fit(a$Y, a$Z, a$W, a$D, a$D,
      if (pate) a$D else NULL, estimand = target, ps_weighting = weighting)
    before <- f
    b <- wm_bootstrap(f, counts = a$counts, method = "fixed_reuse")
    expected <- vapply(seq_len(ncol(a$counts)), function(j) {
      m <- a$counts[, j]
      wm_bootstrap_literal_ratio(f$fit, m,
        wm_bootstrap_literal_prediction(a, m, target, weighting))
    }, numeric(1))
    expect_identical(f, before)
    expect_equal(b$draws, expected, tolerance = 1e-11)
    expect_equal(b$draws[1L], f$estimate, tolerance = 1e-11)
    expect_equal(b$variance, mean((expected - mean(expected))^2), tolerance = 1e-11)
    expect_identical(b$variance_divisor, "B")
    expect_identical(b$refit_mode, "complete_retained_pipeline")
    expect_identical(b$predictions, "refitted")
    expect_false(b$population_assumptions_verified)
    # Explicit callback and adapter are identical statistics; callback receives
    # all-n multiplicities, including rows with zero count.
    callback <- function(m) wm_bootstrap_literal_prediction(a, m, target, weighting)
    supplied <- wm_bootstrap(f, counts = a$counts, method = "fixed_reuse", refit = callback)
    expect_equal(supplied$draws, expected, tolerance = 1e-11)
    expect_identical(supplied$refit_mode, "supplied_callback")
  }
})

test_that("model-list wrappers retain generic and legacy count-refit routes", {
  a <- wm_bootstrap_fit_fixture()
  for (target in c("PATE", "PATT")) for (generic in c(FALSE, TRUE)) {
    pate <- target == "PATE"
    models <- list(joint = list(design = a$D))
    if (generic) models$second <- list(design = a$D[, c("intercept", "x2"), drop = FALSE],
      weighting = "unit")
    f <- wm_model_fit(a$Y, a$Z, a$W, ps_models = models,
      pg0_models = list(list(design = a$D)),
      pg1_models = if (pate) list(list(design = a$D)) else NULL,
      estimand = target, inference = "full_x")
    before <- f
    counts <- a$counts
    if (generic) {
      counts <- matrix(1, a$n, 3L)
      counts[1:2, 2L] <- c(0, 2)
      counts[1:2, 3L] <- c(2, 0)
    }
    out <- wm_bootstrap(f, counts = counts, method = "fixed_reuse")
    # Existing guarded constructor is a separate lower-level reference for the
    # generic p stack; independent explicit recipe is tested for legacy above.
    expected <- vapply(seq_len(ncol(counts)), function(j) {
      current <- wdsmatch:::.wm_model_prediction_stack(a$Y, a$Z, a$W,
        ps_models = models, pg0_models = list(list(design = a$D)),
        pg1_models = if (pate) list(list(design = a$D)) else NULL,
        estimand = target, multiplicity = counts[, j])
      wm_bootstrap_literal_ratio(f$fit, counts[, j], current)
    }, numeric(1))
    expect_identical(f, before)
    expect_equal(out$draws, expected, tolerance = 1e-11)
    expect_equal(wm_bootstrap(f, B = 2L, seed = 31)$conditional_variance,
      f$variance, tolerance = 1e-11)
  }
})

test_that("wrapper options reject unsupported scopes and count failures", {
  a <- wm_bootstrap_fit_fixture()
  f <- wm_wdsm_fit(a$Y, a$Z, a$W, a$D, a$D, a$D)
  expect_error(wm_bootstrap(f, B = 2L), "point-only", fixed = TRUE)
  expect_error(wm_bootstrap(f, method = "fixed_reuse"), "supplied n-by-B", fixed = TRUE)
  expect_error(wm_bootstrap(f, counts = a$counts, method = "fixed_reuse", interval = "basic"), "basic is unsupported", fixed = TRUE)
  expect_error(wm_bootstrap(f, counts = a$counts, method = "fixed_reuse", seed = 1), "frozen counts", fixed = TRUE)
  expect_error(wm_bootstrap(f, counts = a$counts, method = "fixed_reuse", B = 4), "count columns", fixed = TRUE)
  wrong <- a$counts; rownames(wrong) <- as.character(rev(seq_len(a$n)))
  expect_error(wm_bootstrap_refit(f$fit, wrong), "original observation indices", fixed = TRUE)
  expect_error(wm_bootstrap(f, counts = wrong, method = "fixed_reuse"), "original observation indices", fixed = TRUE)
  invalid <- a$counts
  invalid[, 2L] <- 0; invalid[which(a$Z == 1)[1L], 2L] <- a$n
  expect_error(wm_bootstrap(f, counts = invalid, method = "fixed_reuse"), "Both arms", fixed = TRUE)
  # Both arms are present but the requested root has only two positive rows.
  sparse <- a$counts; sparse[, 2L] <- 0
  sparse[c(which(a$Z == 0)[1L], which(a$Z == 1)[1L]), 2L] <- a$n/2
  expect_error(wm_bootstrap(f, counts = sparse, method = "fixed_reuse"), "replicate 2/3 failed", fixed = TRUE)
  calls <- 0L
  fail <- function(m) {
    calls <<- calls + 1L
    if (calls == 2L) stop("deliberate complete-root failure")
    list(mean0 = f$fit$predictions$mean0, mean1 = f$fit$predictions$mean1)
  }
  expect_error(wm_bootstrap(f, counts = a$counts, method = "fixed_reuse", refit = fail), "replicate 2/3 failed", fixed = TRUE)
  expect_identical(calls, 2L)
  special <- f; special$correction_fit <- list(method = "quadratic_qr")
  expect_error(wm_bootstrap(special, counts = a$counts, method = "fixed_reuse"), "explicit compatible refit callback", fixed = TRUE)
  bad <- f; bad$nuisance$inputs$weights[1L] <- 2*bad$nuisance$inputs$weights[1L]
  expect_error(wm_bootstrap(bad, counts = a$counts, method = "fixed_reuse"), "Nuisance supplied W", fixed = TRUE)
})

test_that("cached inference from another valid source cannot replace the wrapper", {
  a <- wm_bootstrap_fit_fixture()
  for (target in c("PATE", "PATT")) {
    pate <- target == "PATE"
    f <- wm_wdsm_fit(a$Y, a$Z, a$W, a$D, a$D,
      if (pate) a$D else NULL, estimand = target, inference = "full_x")
    other <- wm_wdsm_fit(a$Y + .2*a$Z, a$Z, a$W, a$D, a$D,
      if (pate) a$D else NULL, estimand = target, inference = "full_x")
    prepare <- function(x) wm_model_inference(x, if (pate) "full_x" else "patt",
      transport0 = if (pate) NULL else list(mode = "zero", basis = "current_centering"))
    current <- prepare(f); foreign <- prepare(other)
    # Both public preparations independently pass their own source validation.
    expect_equal(wm_bootstrap(current, B = 2L, seed = 31L)$estimate, f$estimate)
    expect_equal(wm_bootstrap(foreign, B = 2L, seed = 31L)$estimate, other$estimate)
    expect_gt(abs(foreign$estimate-current$estimate), .1)
    current$inference <- foreign$inference
    expect_error(wm_bootstrap(current, B = 2L, seed = 31L),
      "exact original fitted-wrapper source fit", fixed = TRUE)
    # The same exact-fit binding precedes all public advanced dispatches.
    for (scope in c("critical_planar", "model_quadratic_contrast_inputs_ready")) {
      replaced <- current
      if (scope == "critical_planar") replaced$inference$covariance_scope <- scope else
        replaced$inference$status <- scope
      expect_error(wm_bootstrap(replaced, B = 2L, seed = 31L),
        "exact original fitted-wrapper source fit", fixed = TRUE)
    }
  }
})
