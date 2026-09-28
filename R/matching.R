wdsm_point <- function(Y, Z, weights, scores, M, estimand,
    tie_seed = 20260917L, tie_tolerance = 64 * .Machine$double.eps) {
  n <- length(Y)
  treated <- which(Z == 1)
  control <- which(Z == 0)
  allocation <- wdsm_make_matches(Z, scores$D0, scores$D1, M, estimand,
    tie_seed = tie_seed, tie_tolerance = tie_tolerance)
  matches_0 <- allocation$matches_0
  matches_1 <- allocation$matches_1
  K <- wdsm_reuse(Z, weights, matches_0, matches_1, estimand)
  contrast <- rep(NA_real_, n)
  for (i in treated) {
    donors <- matches_0[[i]]
    imputed <- scores$q0[i] +
      sum(weights[donors] * (Y[donors] - scores$q0[donors])) / sum(weights[donors])
    contrast[i] <- Y[i] - imputed
  }
  if (estimand == "PATE") for (i in control) {
    donors <- matches_1[[i]]
    imputed <- scores$q1[i] +
      sum(weights[donors] * (Y[donors] - scores$q1[donors])) / sum(weights[donors])
    contrast[i] <- imputed - Y[i]
  }
  targets <- if (estimand == "PATE") seq_len(n) else treated
  estimate <- sum(weights[targets] * contrast[targets]) / sum(weights[targets])
  check <- wdsm_replicate(Y, Z, weights, K, rep(1, n), scores$q0, scores$q1, estimand)
  identity_error <- abs(estimate - check)
  if (!is.finite(identity_error) || identity_error > 1e-9 * (1 + abs(estimate)))
    stop("Direct imputation and original-reuse representations disagree", call. = FALSE)
  list(estimate = estimate, Y = Y, Z = Z, weights = weights,
       K = K, matches_0 = matches_0, matches_1 = matches_1,
       scores = scores, estimand = estimand, identity_error = identity_error,
       tie_diagnostics = allocation$tie_diagnostics)
}
