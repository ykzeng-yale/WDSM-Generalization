# Independent finite-sample estimator checks. These do not establish population
# identification, residual centering or confidence-interval coverage.
wm_direct_imputation_reference <- function(Y, Z, W, S0, S1, M,
                                            estimand, mean0, mean1) {
  n <- length(Y)
  potential <- matrix(NA_real_, n, 2L)
  potential[cbind(seq_len(n), Z + 1L)] <- Y
  neighbors <- vector("list", n)
  queries <- if (estimand == "PATE") seq_len(n) else which(Z == 1L)
  for (i in queries) {
    donor_arm <- 1L - Z[i]
    pool <- which(Z == donor_arm)
    coordinates <- if (donor_arm == 0L) S0 else S1
    # Exhaustively sort squared distances; no package graph helpers are used.
    squared <- vapply(pool, function(j)
      sum((coordinates[i, ] - coordinates[j, ])^2), numeric(1L))
    donors <- pool[order(squared, pool)[seq_len(M)]]
    neighbors[[i]] <- donors
    mu <- if (donor_arm == 0L) mean0 else mean1
    potential[i, donor_arm + 1L] <- mu[i] +
      sum(W[donors] * (Y[donors] - mu[donors])) / sum(W[donors])
  }
  contrasts <- if (estimand == "PATE") potential[, 2L] - potential[, 1L] else
    Y - potential[, 1L]
  list(estimate = sum(W[queries] * contrasts[queries]) / sum(W[queries]),
       denominator = sum(W[queries]), potential = potential,
       neighbors = neighbors)
}

test_that("weighted estimates use direct query denominators and each donor map", {
  set.seed(100503L)
  n <- 64L
  Z <- rep(0:1, n / 2L)
  X <- matrix(runif(n * 6L, -1, 1), n, 6L)
  # Known treatment-dependent weights test both outer and donor normalization.
  W <- exp(0.4 * Z + 0.35 * X[, 1L] - 0.2 * X[, 2L])
  Y <- 1 + X[, 2L]^2 + Z * (1.2 + 0.2 * X[, 1L]) +
    0.2 * sin(seq_len(n))
  for (d in c(1L, 2L, 6L)) for (M in c(1L, 3L)) {
    S0 <- X[, seq_len(d), drop = FALSE]
    # Distinct arm maps can differ in coordinates and dimension.
    S1 <- X[, seq_len(if (d == 6L) 2L else d + 1L), drop = FALSE]
    mean0 <- 0.3 + S0[, 1L]^2
    mean1 <- 1.2 - 0.2 * S1[, 1L]
    for (estimand in c("PATE", "PATT")) {
      args <- list(Y = Y, Z = Z, weights = W, scores0 = S0, M = M,
                   estimand = estimand, mean0 = mean0, variance = FALSE)
      if (estimand == "PATE") {
        args$scores1 <- S1
        args$mean1 <- mean1
      }
      fit <- do.call(wm_match, args)
      expected <- wm_direct_imputation_reference(Y, Z, W, S0, S1, M,
        estimand, mean0, mean1)
      expect_identical(fit$graph$neighbors, expected$neighbors)
      expect_equal(fit$estimate, expected$estimate, tolerance = 1e-12)
      expect_equal(fit$denominator * fit$weight_scale, expected$denominator)
      expect_equal(unname(fit$imputed), expected$potential, tolerance = 1e-12)
      expect_true(is.na(fit$variance))

      # A common weight-unit change preserves all donor fractions and effects.
      args$weights <- 17 * W
      rescaled <- do.call(wm_match, args)
      expect_identical(rescaled$graph$neighbors, expected$neighbors)
      expect_equal(rescaled$estimate, expected$estimate, tolerance = 1e-12)

      # Unit weights reduce to the ordinary arithmetic donor average and n/n1
      # outer denominators; this reference uses the estimator definition.
      args$weights <- rep(1, n)
      unweighted <- do.call(wm_match, args)
      ordinary <- wm_direct_imputation_reference(Y, Z, rep(1, n), S0, S1,
        M, estimand, mean0, mean1)
      expect_equal(unweighted$estimate, ordinary$estimate, tolerance = 1e-12)
      expect_equal(unweighted$denominator,
                   if (estimand == "PATE") n else sum(Z))
      expect_equal(unweighted$graph$edges$share,
                   rep(1 / M, nrow(unweighted$graph$edges)))
    }
  }
})
