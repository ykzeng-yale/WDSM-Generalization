# One finite fixture fit checks the new stack's arithmetic. It is not a
# repeated-sampling study or a verification of population model/geometry laws.
test_that("finite model lists retain complete shared and arm-specific derivatives", {
  # Reuse the existing test-wm-wdsm-ps-weighting.R fixture exactly.
  set.seed(29317)
  n <- 320L
  x <- cbind(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
  D <- cbind(intercept = 1, x)
  G <- cbind(D, product = x[, 1] * x[, 2])
  Z <- rbinom(n, 1, plogis(-0.2 + 0.5 * x[, 1] - 0.3 * x[, 2]))
  Y <- 1 + x[, 1] + 0.4 * x[, 2] + Z * (0.8 + x[, 2]) + rnorm(n, sd = 0.4)
  W <- exp(0.5 * Z + 0.4 * x[, 1] - 0.2 * x[, 2])
  f <- wm_model_fit(Y, Z, W,
    ps_models = list(joint = list(design = D, weighting = "probability"),
                     second = list(design = D[, c("intercept", "x2"), drop = FALSE],
                                   weighting = "unit")),
    pg0_models = list(product = list(design = G)),
    pg1_models = list(linear = list(design = D)),
    M = 3L, estimand = "PATE", inference = "full_x")
  s <- f$nuisance
  expect_s3_class(f, "wm_model_fit")
  expect_length(s$parameter, 40L)
  expect_identical(unname(s$layout$matching_dimensions), c(3L, 3L))
  expect_identical(unname(f$allocation$basis_columns), c(10L, 10L))
  expect_length(s$blocks$ps, 2L)
  expect_identical(sum(lengths(s$blocks$ps)), 5L)
  expect_identical(s$scores0[, 1:2], s$scores1[, 1:2])
  expect_identical(s$inputs$weights, W)
  expect_false(f$assumptions_verified)
  expect_false(f$bootstrap_requested)
  expect_false(f$inference$original_refit_limit_agreement_declared)

  # Primitive rows and polynomial values are reconstructed without the
  # candidate's basis/derivative helpers, at the one retained fitted parameter.
  th <- s$parameter
  k <- s$blocks
  psi <- matrix(0, n, length(th))
  raw <- NULL
  for (j in seq_along(s$inputs$ps_models)) {
    m <- s$inputs$ps_models[[j]]
    r <- plogis(drop(m$design %*% th[k$ps[[j]]]) + m$offset)
    w <- if (m$weighting == "probability") W / max(W) else rep(1, n)
    psi[, k$ps[[j]]] <- m$design * (w * (Z - r))
    raw <- cbind(raw, r)
  }
  for (z in s$used_arms) {
    arm <- as.character(z)
    for (j in seq_along(s$inputs$pg_models[[arm]])) {
      m <- s$inputs$pg_models[[arm]][[j]]
      ix <- k$pg[[arm]][[j]]
      r <- drop(m$design %*% th[ix]) + m$offset
      psi[, ix] <- m$design * ((Z == z) * (Y - r))
      raw <- cbind(raw, r)
    }
  }
  colnames(raw) <- s$layout$raw_names
  center <- sweep(raw, 2L, th[k$center], "-")
  standardized <- sweep(center, 2L, sqrt(th[k$variance]), "/")
  psi[, k$center] <- center
  psi[, k$variance] <- sweep(center^2, 2L, th[k$variance], "-")
  for (z in s$used_arms) {
    arm <- as.character(z)
    score <- standardized[, s$layout$scores[[arm]], drop = FALSE]
    b <- cbind(1, score)
    for (a in seq_len(ncol(score))) for (bb in a:ncol(score)) {
      b <- cbind(b, score[, a] * score[, bb])
    }
    expect_equal(unname(b), unname(s$basis[[arm]]), tolerance = 1e-12)
    r <- drop(b %*% th[k$correction[[arm]]])
    psi[, k$correction[[arm]]] <- b * ((Z == z) * W / max(W) * (Y - r))
  }
  expect_equal(unname(s$estimating_equations), psi, tolerance = 1e-11)
  expect_lt(max(abs((-s$nuisance_influence %*% t(s$jacobian) - psi) /
                     (1 + abs(psi)))), 1e-9)
  ell <- sweep(s$nuisance_influence, 2L, colMeans(s$nuisance_influence), "-")
  expect_equal(unname(s$nuisance_covariance), unname(crossprod(ell) / n),
               tolerance = 1e-9)

  # Parameter perturbations evaluate the same object; they do not refit models.
  evaluator <- WeightedMatching:::.wm_model_nuisance_evaluate
  numeric_A <- matrix(0, length(th), length(th))
  derivative_error <- 0
  for (j in seq_along(th)) {
    h <- 1e-6 * (1 + abs(th[j]))
    if (j %in% k$variance) h <- min(h, th[j] * 1e-5)
    plus <- minus <- th
    plus[j] <- plus[j] + h
    minus[j] <- minus[j] - h
    ep <- evaluator(s, plus)
    em <- evaluator(s, minus)
    numeric_A[, j] <- (ep$mean_equation - em$mean_equation) / (2 * h)
    for (z in s$used_arms) {
      arm <- as.character(z)
      sf <- paste0("scores", z)
      mf <- paste0("mean", z)
      numerical <- (ep[[sf]] - em[[sf]]) / (2 * h)
      analytic <- do.call(cbind, lapply(s$score_derivatives[[arm]], function(v) v[, j]))
      derivative_error <- max(derivative_error,
        abs(numerical - analytic) / (1 + abs(analytic)))
      numerical <- (ep[[mf]] - em[[mf]]) / (2 * h)
      analytic <- s[[paste0(mf, "_derivative")]][, j]
      derivative_error <- max(derivative_error,
        abs(numerical - analytic) / (1 + abs(analytic)))
    }
  }
  expect_lt(max(abs(numeric_A - s$jacobian) / (1 + abs(s$jacobian))), 5e-5)
  expect_lt(derivative_error, 5e-5)

  expect_true(all(lengths(f$fit$graph$neighbors) == 3L))
  y0 <- y1 <- Y
  for (i in seq_len(n)) {
    donors <- f$fit$graph$neighbors[[i]]
    mean_z <- if (Z[i] == 1L) s$mean0 else s$mean1
    imputed <- mean_z[i] + sum(W[donors] * (Y[donors] - mean_z[donors])) / sum(W[donors])
    if (Z[i] == 1L) y0[i] <- imputed else y1[i] <- imputed
  }
  expect_equal(f$estimate, sum(W * (y1 - y0)) / sum(W), tolerance = 1e-11)
  a <- f$inference
  b <- a$total_sensitivity
  expect_length(b, length(th))
  expect_true(is.matrix(a$Sigma))
  expect_equal(a$root_n_variance,
    a$V0 + 2 * sum(b * a$C) + drop(t(b) %*% a$Sigma %*% b), tolerance = 1e-9)
  expect_equal(a$variance, a$root_n_variance / n, tolerance = 1e-12)
})

test_that("model layouts and resource guards retain the complete requested stack", {
  specs <- WeightedMatching:::.wm_model_specs
  layout <- WeightedMatching:::.wm_model_layout
  budget <- WeightedMatching:::.wm_model_allocation
  D <- cbind(intercept = 1, x = seq(-1, 1, length.out = 8L))
  ps <- specs(list(one = list(design = D)), 8L, "ps", "ps_models")
  pg <- specs(list(one = list(design = D)), 8L, "pg", "pg0_models")
  expect_identical(unname(layout(ps, list("0" = list()), 0L)$matching_dimensions), 1L)
  expect_identical(unname(layout(list(), list("0" = pg), 0L)$matching_dimensions), 1L)
  expect_identical(unname(layout(ps, list("0" = pg, "1" = c(pg, pg)), 0:1)$matching_dimensions),
                   c(2L, 3L))
  expect_error(layout(list(), list("0" = list()), 0L), "at least one")
  expect_error(specs(list(list(design = D, unsupported = TRUE)), 8L, "ps", "ps_models"),
               "supported fields")
  expect_error(specs(list(list(design = D, weighting = "estimated")), 8L, "ps", "ps_models"),
               "probability.*unit")
  duplicated <- WeightedMatching:::.wm_model_basis(cbind(D[, 2L], D[, 2L]))
  expect_identical(ncol(duplicated), 6L)
  expect_lt(qr(duplicated)$rank, ncol(duplicated))

  # Allocation-only regression check: no large matrix or nuisance fit is made.
  one_ps <- list(list(design = matrix(0, 1L, 32L)))
  empty_pg <- list("0" = list(), "1" = list())
  l <- layout(one_ps, empty_pg, 0:1)
  admitted <- budget(10000L, one_ps, empty_pg, l, 3100000)
  expect_identical(admitted$parameter_dimension, 40L)
  expect_equal(admitted$planned_elements, 3004800)
  expect_gte(admitted$planned_elements, 6 * 10000 * 40 + 2 * 40^2)
  expect_error(budget(10000L, one_ps, empty_pg, l, 2500000), "max_stack_elements")

  # These public checks fail before any regression or graph construction.
  y <- seq_len(8L)
  z <- rep(0:1, each = 4L)
  expect_error(wm_model_fit(y, z, rep(1, 8L), ps_models = ps,
    pg1_models = list(), estimand = "PATT"), "unused treated")
  expect_error(wm_model_fit(y, z, rep(1, 8L), ps_models = ps,
    inference = "full_x", tie_rule = "source_random", tie_seed = 7L), "point-only")
  expect_error(wm_model_fit(y, z, rep(1, 8L), pg0_models = pg,
    pg1_models = list(one = list(design = D), two = list(design = D)),
    inference = "full_x"), "mixed scalar/vector")
})
