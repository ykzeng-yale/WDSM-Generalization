# Finite-data reference checks only. These literal fixtures do not certify
# continuous population geometry, own-field centering, fitted roots or coverage.
wm_joint_fixture <- function(M = 1L, d = 3L, unit = FALSE) {
  x <- c(.02, .13, .29, .48, .69, .91, .07, .21, .36, .56, .77, .96)
  u <- c(.2, .8, .4, .1, .9, .5, .6, .3, .7, .4, .2, .8)
  v <- c(.7, .1, .9, .3, .6, .2, .5, .8, .1, .9, .4, .6)
  S0 <- cbind(x, u, v)
  S1 <- cbind(x^2, v, u)
  if (d == 6L) {
    S0 <- cbind(S0, x * u, u * v, x + v)
    S1 <- cbind(S1, x * v, x * u, u + v)
  }
  n <- length(x)
  list(data = list(Y = c(1, -2, 4, 0, 3, 2, 4, 1, -1, 5, 2, 6),
    Z = rep(0:1, each = 6L), weights = if (unit) rep(1, n) else 1 + seq_len(n) / 10,
    scores0 = S0, scores1 = S1, M = M, estimand = "PATE",
    mean0 = x, mean1 = 1 + x^2, variance = FALSE),
    influence = cbind(p = x - .4, q = u - .5),
    D0 = cbind(p = x, q = v), D1 = cbind(p = u, q = x^2))
}

wm_joint_inference <- function(a, scope = "regular_joint_d_gt2", fit = NULL) {
  if (is.null(fit)) fit <- do.call(wm_match, a$data)
  transport <- list(mode = "zero", basis = "fixed_map")
  arguments <- list(fit = fit, nuisance_influence = a$influence,
    mean_derivative0 = a$D0, mean_derivative1 = a$D1, covariance_scope = scope)
  if (scope != "full_x") arguments[c("transport0", "transport1")] <-
    list(transport, transport)
  do.call(wm_fitted_inference, arguments)
}

# Independent query imputation/reindexing reference: raw weights and complete
# donor totals; no stored load, fraction, contribution or variance arrays used.
wm_joint_direct <- function(a, neighbors) {
  dat <- a$data; W <- dat$weights; Z <- dat$Z; n <- length(Z)
  eps <- dat$Y - ifelse(Z == 1L, dat$mean1, dat$mean0)
  load <- total <- contrast <- numeric(n)
  derivative <- numeric(ncol(a$influence))
  for (i in seq_len(n)) {
    donor <- neighbors[[i]]
    total[i] <- sum(W[donor])
    fraction <- W[donor] / total[i]
    load[donor] <- load[donor] + W[i] * fraction
    mu <- if (Z[i] == 1L) dat$mean0 else dat$mean1
    D <- if (Z[i] == 1L) a$D0 else a$D1
    contrast[i] <- (2 * Z[i] - 1) *
      (dat$Y[i] - mu[i] - sum(fraction * eps[donor]))
    derivative <- derivative + W[i] * (1 - 2 * Z[i]) *
      (D[i, ] - colSums(D[donor, , drop = FALSE] * fraction))
  }
  estimate <- sum(W * contrast) / sum(W)
  gamma <- mean(W)
  base <- (W * (dat$mean1 - dat$mean0 - estimate) +
    (2 * Z - 1) * (W + load) * eps) / gamma
  base <- base - mean(base)
  influence <- sweep(a$influence, 2L, colMeans(a$influence), "-")
  b <- derivative / sum(W)
  C <- as.vector(crossprod(influence, base)) / n
  Sigma <- crossprod(influence) / n
  R <- 0; pairs <- 0L
  for (i in which(Z == 1L)) for (j in which(Z == 0L)) {
    if (j %in% neighbors[[i]] && i %in% neighbors[[j]]) {
      pairs <- pairs + 1L
      R <- R + W[i]^2 * W[j]^2 * eps[i] * eps[j] / (total[i] * total[j])
    }
  }
  R <- R / n
  subtraction <- 2 * R / gamma^2
  diagonal <- mean(base^2)
  cross <- 2 * sum(b * C)
  nuisance <- as.numeric(crossprod(b, Sigma %*% b))
  list(estimate = estimate, base = base, b = b, C = C, Sigma = Sigma,
    diagonal = diagonal, R = R, subtraction = subtraction, pairs = pairs,
    V0 = diagonal - subtraction, cross = cross, nuisance = nuisance,
    V = diagonal - subtraction + cross + nuisance)
}

test_that("qualified fitted reciprocal covariance agrees with complete direct sums", {
  for (M in c(1L, 3L)) for (d in c(3L, 6L)) for (unit in c(FALSE, TRUE)) {
    a <- wm_joint_fixture(M, d, unit)
    fit <- do.call(wm_match, a$data); before <- fit
    out <- wm_joint_inference(a, fit = fit)
    ref <- wm_joint_direct(a, fit$graph$neighbors)
    expect_identical(fit, before)
    expect_identical(out$fit, fit)
    expect_equal(out$estimate, ref$estimate, tolerance = 1e-12)
    expect_equal(out$base_rows, ref$base, tolerance = 1e-12)
    expect_equal(unname(out$total_sensitivity), unname(ref$b), tolerance = 1e-12)
    expect_equal(unname(out$C), ref$C, tolerance = 1e-12)
    expect_equal(unname(out$Sigma), unname(ref$Sigma), tolerance = 1e-12)
    expect_equal(out$reciprocal_numerator * max(a$data$weights)^2, ref$R,
                 tolerance = 1e-12)
    expect_equal(nrow(out$reciprocal_pairs), ref$pairs)
    expect_equal(out$diagonal_root_n_variance, ref$diagonal, tolerance = 1e-12)
    expect_equal(out$reciprocal_subtraction, ref$subtraction, tolerance = 1e-12)
    expect_equal(out$V0, ref$V0, tolerance = 1e-12)
    expect_equal(out$cross_term, ref$cross, tolerance = 1e-12)
    expect_equal(out$nuisance_variance, ref$nuisance, tolerance = 1e-12)
    expect_equal(out$unregularized_root_n_variance, ref$V, tolerance = 1e-12)
    expect_identical(out$available, ref$V > 0)
    if (ref$V > 0) expect_equal(out$variance, ref$V / length(a$data$Z))
    else expect_true(is.na(out$se))
    scaled <- a; scaled$data$weights <- 29 * scaled$data$weights
    scale_out <- wm_joint_inference(scaled)
    expect_equal(scale_out$estimate, out$estimate, tolerance = 1e-12)
    expect_equal(scale_out$V0, out$V0, tolerance = 1e-12)
    expect_equal(scale_out$unregularized_root_n_variance,
                 out$unregularized_root_n_variance, tolerance = 1e-12)
    expect_false(out$assumptions_verified)
    expect_false(out$application_verified)
    expect_false(out$floor_active)
  }
})

# Exact two-row identities deliberately exercise finite signed covariance.
# Declaring three retained coordinates does not claim regular 3d population law.
wm_joint_signed <- function(negative_R = FALSE, influence_scale = 4) {
  fit <- wm_match(c(1, 1), c(0, 1), c(1, 1),
    scores0 = cbind(c(0, 1), c(0, 2), c(0, 3)), M = 1,
    mean0 = c(0, 2), mean1 = if (negative_R) c(2, 2) else c(2, 0),
    variance = FALSE)
  transport <- list(mode = "zero", basis = "fixed_map")
  wm_fitted_inference(fit, cbind(p = c(-influence_scale, influence_scale)),
    mean_derivative0 = cbind(p = c(1, 0)),
    transport0 = transport, transport1 = transport,
    covariance_scope = "regular_joint_d_gt2")
}

test_that("signed base is retained and final contrast governs availability", {
  out <- wm_joint_signed()
  expect_equal(out$reciprocal_numerator, .5)
  expect_equal(out$diagonal_root_n_variance, 0)
  expect_equal(out$V0, -1)
  expect_equal(out$total_sensitivity, c(p = .5))
  expect_equal(out$nuisance_variance, 4)
  expect_equal(out$root_n_variance, 3)
  expect_true(out$available)
  expect_equal(out$variance, 1.5)
  negative <- wm_joint_signed(TRUE)
  expect_equal(negative$reciprocal_numerator, -.5)
  expect_equal(negative$reciprocal_subtraction, -1)
  expect_equal(negative$V0, 2)
  expect_equal(negative$root_n_variance, 2)
  failed <- wm_joint_signed(influence_scale = 0)
  expect_equal(failed$V0, -1)
  expect_equal(failed$unregularized_root_n_variance, -1)
  expect_false(failed$available)
  expect_true(is.na(failed$root_n_variance))
  expect_true(is.na(failed$variance))
  expect_true(is.na(failed$se))
  expect_true(all(is.na(failed$conf.int)))
  expect_error(wm_bootstrap(failed, seed = 5), "successful")
})

test_that("corrected fitted Gaussian dispatch uses full V and preserves RNG", {
  out <- wm_joint_signed()
  set.seed(119); previous <- .Random.seed
  set.seed(45); normal <- rnorm(11); .Random.seed <- previous
  boot <- wm_bootstrap(out, B = 11, seed = 45, chunk_size = 3)
  expect_identical(.Random.seed, previous)
  expect_equal(boot$root_n_draws, sqrt(3) * normal)
  expect_equal(boot$draws, out$estimate + boot$root_n_draws / sqrt(2))
  expect_equal(boot$conditional_root_n_variance, 3)
  expect_equal(boot$conditional_variance, 1.5)
  expect_equal(boot$monte_carlo_root_n_variance, var(boot$root_n_draws))
  expect_equal(boot$conf.int, out$conf.int)
  expect_identical(boot$source_inference, out)
  expect_identical(boot$draw_input, "generated_gaussian_corrected_variance")
  expect_false(boot$original_refit_bootstrap)
  single <- wm_bootstrap(out, B = 1, seed = 45)
  expect_true(is.na(single$monte_carlo_variance))
  expect_equal(single$conf.int, boot$conf.int)
  basic <- wm_bootstrap(out, B = 11, seed = 45, interval = "basic")
  expect_equal(unname(basic$conf.int), out$estimate -
    unname(quantile(boot$root_n_draws, c(.975, .025))) / sqrt(out$n))
  expect_null(wm_bootstrap(out, B = 1, seed = 45, interval = "none")$conf.int)
  expect_error(wm_bootstrap(out, counts = matrix(1, 2, 2)), "rejects supplied counts")
  expect_error(wm_bootstrap(out, refit = function(...) stop("never")), "fixed_reuse")
  expect_identical(.Random.seed, previous)
})

test_that("new scope rejects incompatible dimensions and supplied-weight uncertainty", {
  for (d in c(1L, 2L)) {
    a <- wm_joint_fixture(); a$data$scores0 <- a$data$scores0[, seq_len(d), drop = FALSE]
    a$data$scores1 <- a$data$scores1[, seq_len(d), drop = FALSE]
    expect_error(wm_joint_inference(a), "equal used dimensions d > 2")
  }
  a <- wm_joint_fixture(); a$data$scores1 <- cbind(a$data$scores1, extra = seq_len(12))
  expect_error(wm_joint_inference(a), "equal used dimensions d > 2")
  a <- wm_joint_fixture(); fit <- do.call(wm_match, a$data)
  expect_error(wm_fitted_inference(fit, a$influence,
    weight_derivative = matrix(1e-12, 12, 2, dimnames = list(NULL, c("p", "q"))),
    covariance_scope = "regular_joint_d_gt2"), "known")
  expect_error(wm_fitted_inference(fit, a$influence,
    covariance_scope = "regular_joint_d_gt2"), "transport0")
  a$data$estimand <- "PATT"; a$data[c("scores1", "mean1")] <- NULL
  fit <- do.call(wm_match, a$data)
  expect_error(wm_fitted_inference(fit, a$influence,
    covariance_scope = "regular_joint_d_gt2"), "PATT")
})

test_that("corrected Gaussian seed restores absence and rejects Box-Muller cache changes", {
  out <- wm_joint_signed()
  kind <- RNGkind()
  had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  previous <- if (had_seed) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    do.call(RNGkind, as.list(kind))
    if (had_seed) assign(".Random.seed", previous, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", normal.kind = "Inversion")
  rm(".Random.seed", envir = .GlobalEnv)
  invisible(wm_bootstrap(out, B = 2, seed = 5))
  expect_false(exists(".Random.seed", .GlobalEnv, inherits = FALSE))
  RNGkind("Mersenne-Twister", normal.kind = "Box-Muller")
  set.seed(67); expected <- rnorm(2)
  set.seed(67); expect_equal(rnorm(1), expected[1])
  before <- .Random.seed
  expect_error(wm_bootstrap(out, B = 2, seed = 5), "Box-Muller")
  expect_identical(.Random.seed, before)
  expect_equal(rnorm(1), expected[2])
})

test_that("corrected bootstrap rejects stale finite provenance before drawing", {
  out <- wm_joint_signed(); set.seed(41); previous <- .Random.seed
  changes <- list(
    function(x) {x$V0 <- x$V0 + 1; x},
    function(x) {x$reciprocal_subtraction <- x$reciprocal_subtraction + 1; x},
    function(x) {x$reciprocal_numerator <- x$reciprocal_numerator + 1; x},
    function(x) {x$reciprocal_pairs$product[1] <- 99; x},
    function(x) {x$reciprocal_pairs$control_donor_total[1] <- 99; x},
    function(x) {x$augmented_rows[1] <- 99; x},
    function(x) {x$total_sensitivity[1] <- 99; x},
    function(x) {x$Sigma[1, 1] <- 99; x},
    function(x) {x$fit$data$Y[1] <- 99; x},
    function(x) {x$fit$graph$edges$share[1] <- .5; x},
    function(x) {x$parameter_names <- "other"; x},
    function(x) {x$floor_active <- TRUE; x},
    function(x) {x$analysis_target_weight_mean <- 99; x},
    function(x) {x$unregularized_root_n_variance <- 99; x},
    function(x) {x$smooth_derivative$weight_sensitivity[1] <- 1;
      x$smooth_derivative$prediction_sensitivity[1] <-
        x$smooth_derivative$prediction_sensitivity[1] - 1; x})
  for (change in changes) expect_error(wm_bootstrap(change(out), B = 2, seed = 5))
  for (B in list(0, -1, NA, Inf, 1.2)) expect_error(wm_bootstrap(out, B = B, seed = 5))
  expect_identical(.Random.seed, previous)
})


test_that("private full-slope counts reject the new scope even when observed R is zero", {
  fit <- wm_match(c(0, 1), c(0, 1), c(1, 1),
    scores0 = cbind(c(0, 1), c(0, 2), c(0, 3)), M = 1,
    mean0 = c(0, 0), mean1 = c(0, 0), variance = FALSE)
  reference <- WeightedMatching:::.wm_wdsm_fitted_variance(fit,
    cbind(p = c(-4, 4)), c(p = 0), c(p = 0), "regular_joint_d_gt2")
  expect_equal(reference$reciprocal_numerator, 0)
  expect_true(reference$available)
  expect_error(WeightedMatching:::.wm_wdsm_fitted_count(reference, matrix(1, 2, 2)),
               "cannot use ordinary full-slope count")
})


# Currency changes multiply Y, correction means and their derivatives, while
# the same literal graph, weights and dimensionless nuisance influence remain.
# This is a finite-arithmetic check, not a fitted-model or coverage experiment.
test_that("complete fitted bootstrap is equivariant to outcome currency units", {
  for (estimand in c("PATE", "PATT")) {
    a <- wm_joint_fixture(M = 3L)
    a$data$estimand <- estimand
    a$data$scores0 <- a$data$scores0[, 1L, drop = FALSE]
    a$data$scores1 <- a$data$scores1[, 1L, drop = FALSE]
    if (estimand == "PATT") a$data[c("scores1", "mean1")] <- NULL
    infer <- function(x) {
      fit <- do.call(wm_match, x$data)
      wm_fitted_inference(fit, x$influence, mean_derivative0 = x$D0,
        mean_derivative1 = if (estimand == "PATE") x$D1 else NULL,
        covariance_scope = if (estimand == "PATE") "full_x" else "patt",
        transport0 = if (estimand == "PATT")
          list(mode = "zero", basis = "current_centering") else NULL)
    }
    reference <- infer(a)
    counts <- cbind(rep(1, reference$n), rep(c(2, 0), reference$n / 2))
    reference_boot <- wm_bootstrap(reference, counts = counts)
    for (currency in c(.01, 100, 1e4)) {
      changed <- a
      changed$data$Y <- currency * a$data$Y
      changed$data$mean0 <- currency * a$data$mean0
      if (estimand == "PATE") changed$data$mean1 <- currency * a$data$mean1
      changed$D0 <- currency * a$D0
      changed$D1 <- currency * a$D1
      out <- infer(changed)
      boot <- wm_bootstrap(out, counts = counts)
      expect_identical(out$fit$graph, reference$fit$graph)
      expect_identical(out$fit$weights, reference$fit$weights)
      expect_equal(out$estimate / currency, reference$estimate, tolerance = 1e-12)
      expect_equal(out$base_rows / currency, reference$base_rows, tolerance = 1e-12)
      expect_equal(out$augmented_rows / currency, reference$augmented_rows,
                   tolerance = 1e-12)
      expect_equal(out$total_sensitivity / currency, reference$total_sensitivity,
                   tolerance = 1e-12)
      expect_equal(out$C / currency, reference$C, tolerance = 1e-12)
      expect_equal(out$Sigma, reference$Sigma, tolerance = 1e-12)
      for (field in c("V0", "cross_term", "nuisance_variance", "root_n_variance",
                      "variance")) {
        expect_equal(out[[field]] / currency^2, reference[[field]], tolerance = 1e-12)
      }
      expect_equal(boot$root_n_draws / currency, reference_boot$root_n_draws,
                   tolerance = 1e-12)
      expect_equal(boot$draws / currency, reference_boot$draws, tolerance = 1e-12)
      expect_equal(boot$conditional_variance / currency^2,
                   reference_boot$conditional_variance, tolerance = 1e-12)
      expect_equal(boot$conf.int / currency, reference_boot$conf.int, tolerance = 1e-12)
      expect_identical(boot$source_inference, out)

      # Model a few last-bit differences between equivalent covariance
      # evaluations. No point, row, slope, covariance or variance is changed.
      block_scale <- max(abs(c(out$V0, out$cross_term, out$nuisance_variance)))
      roundoff <- out
      roundoff$covariance_formula_error <- out$covariance_formula_error +
        8 * .Machine$double.eps * block_scale
      if (currency == 1e4) expect_gt(roundoff$covariance_formula_error -
        out$covariance_formula_error, 1e-10)
      retained <- roundoff$covariance_formula_error
      roundoff_boot <- wm_bootstrap(roundoff, counts = counts)
      expect_identical(roundoff_boot$root_n_draws, boot$root_n_draws)
      expect_identical(roundoff_boot$conditional_variance, boot$conditional_variance)
      expect_identical(roundoff_boot$source_inference$covariance_formula_error, retained)

      corrupted <- out
      corrupted$covariance_formula_error <- out$covariance_formula_error + 1e-6 * block_scale
      expect_error(wm_bootstrap(corrupted, counts = counts), "covariance_formula_error")
      corrupted <- out
      corrupted$V0 <- out$V0 + 1e-3 * block_scale
      expect_error(wm_bootstrap(corrupted, counts = counts), "Complete fitted V0")
    }
    for (invalid in c(-1, Inf, NA_real_)) {
      corrupted <- reference; corrupted$covariance_formula_error <- invalid
      expect_error(wm_bootstrap(corrupted, counts = counts), "covariance_formula_error")
    }
  }
})


# These existing literal matching arrays exercise units and centering only;
# no nuisance model is fitted and no population calibration is asserted.
test_that("complete fitted centering validation respects nuisance parameter units", {
  for (estimand in c("PATE", "PATT")) {
    a <- wm_joint_fixture(M = 3L)
    a$data$estimand <- estimand
    a$data$scores0 <- a$data$scores0[, 1:2, drop = FALSE]
    a$data$scores1 <- a$data$scores1[, 1:2, drop = FALSE]
    # Antisymmetric retained marks have zero sum and two zero entries; the
    # near-zero positions expose an absolute-unit validation error directly.
    v <- a$influence[1:5, "p"]
    a$influence[, "p"] <- c(v, 0, -v, 0)
    if (estimand == "PATT") {
      a$data$scores1 <- NULL; a$data$mean1 <- NULL
    }
    infer <- function(x) {
      wm_fitted_inference(do.call(wm_match, x$data), x$influence,
        mean_derivative0 = x$D0,
        mean_derivative1 = if (estimand == "PATE") x$D1 else NULL,
        covariance_scope = if (estimand == "PATE") "full_x" else "patt",
        transport0 = if (estimand == "PATT")
          list(mode = "zero", basis = "current_centering") else NULL)
    }
    reference <- infer(a)
    counts <- cbind(rep(1, reference$n), rep(c(2, 0), reference$n / 2))
    reference_boot <- wm_bootstrap(reference, counts = counts)
    units <- c(p = -1e10, q = 1e-8)
    for (currency in c(.01, 1e4)) {
      b <- a
      b$data$Y <- currency * a$data$Y
      b$data$mean0 <- currency * a$data$mean0
      if (estimand == "PATE") b$data$mean1 <- currency * a$data$mean1
      b$influence <- sweep(a$influence, 2L, units, "*")
      b$D0 <- currency * sweep(a$D0, 2L, units, "/")
      b$D1 <- currency * sweep(a$D1, 2L, units, "/")
      out <- infer(b)
      boot <- wm_bootstrap(out, counts = counts)
      expect_identical(out$fit$graph, reference$fit$graph)
      expect_identical(out$fit$weights, reference$fit$weights)
      expect_equal(out$estimate / currency, reference$estimate, tolerance = 1e-12)
      expect_equal(sweep(out$nuisance_influence, 2L, units, "/"),
                   reference$nuisance_influence, tolerance = 1e-12)
      expect_equal(out$C / (units * currency), reference$C, tolerance = 1e-12)
      expect_equal(out$Sigma / outer(units, units), reference$Sigma, tolerance = 1e-12)
      expect_equal(out$total_sensitivity * units / currency,
                   reference$total_sensitivity, tolerance = 1e-12)
      expect_equal(out$augmented_rows / currency, reference$augmented_rows,
                   tolerance = 1e-12)
      expect_equal(out$variance / currency^2, reference$variance, tolerance = 1e-12)
      expect_equal(boot$root_n_draws / currency, reference_boot$root_n_draws,
                   tolerance = 1e-12)
      expect_equal(boot$conf.int / currency, reference_boot$conf.int, tolerance = 1e-12)

      scale <- max(abs(out$nuisance_influence[, "p"]))
      roundoff <- out
      roundoff$nuisance_influence[, "p"] <- out$nuisance_influence[, "p"] +
        8 * .Machine$double.eps * scale
      recentered <- sweep(roundoff$nuisance_influence, 2L,
                          colMeans(roundoff$nuisance_influence), "-")
      # The exact old elementwise validator rejects the same few-ULP column
      # shift, despite unchanged scientific rows and complete covariance.
      expect_error(.wm_reciprocal_agree(roundoff$nuisance_influence, recentered,
        "Centered joint influence"), "Centered joint influence")
      accepted <- wm_bootstrap(roundoff, counts = counts)
      expect_identical(accepted$root_n_draws, boot$root_n_draws)
      expect_identical(accepted$conditional_variance, boot$conditional_variance)
      expect_identical(accepted$source_inference, roundoff)
      expect_identical(out$nuisance_influence, boot$source_inference$nuisance_influence)

      corrupted <- out
      corrupted$nuisance_influence[, "p"] <- out$nuisance_influence[, "p"] + 1e-6 * scale
      expect_error(wm_bootstrap(corrupted, counts = counts), "Centered joint influence")
      corrupted <- out; corrupted$Sigma["p", "p"] <- 1.001 * out$Sigma["p", "p"]
      expect_error(wm_bootstrap(corrupted, counts = counts), "Complete joint covariance")
    }
    # Exact zero influence coordinates remain admissible with their complete
    # zero covariance row/column, rather than a unit floor or parameter drop.
    zero <- a
    zero$influence <- cbind(a$influence, zero = 0)
    zero$D0 <- cbind(a$D0, zero = 0); zero$D1 <- cbind(a$D1, zero = 0)
    zero_out <- infer(zero)
    zero_boot <- wm_bootstrap(zero_out, counts = counts)
    expect_identical(zero_out$nuisance_influence[, "zero"], rep(0, reference$n))
    expect_equal(zero_boot$root_n_draws, reference_boot$root_n_draws, tolerance = 1e-12)
    expect_equal(zero_out$Sigma["zero", ], c(p = 0, q = 0, zero = 0))
  }
})
