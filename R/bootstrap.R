wdsm_bootstrap <- function(point_est, boots = 200, alpha = 0.05,
                            model.ps = NULL, model.pg = NULL,
                            sampling = "retrospective", ps = NULL, pg = NULL,
                            use.bias.correction = TRUE, bootstrap_counts = NULL) {
  Y <- point_est$Y
  Z <- point_est$Z
  w <- point_est$weights
  n <- length(Y)
  if (!is.null(bootstrap_counts)) {
    if (!is.matrix(bootstrap_counts) || !is.numeric(bootstrap_counts) ||
        !identical(dim(bootstrap_counts), c(as.integer(n), as.integer(boots))) ||
        any(!is.finite(bootstrap_counts)) || any(bootstrap_counts < 0) ||
        any(bootstrap_counts != floor(bootstrap_counts)) ||
        any(colSums(bootstrap_counts) != n))
      stop("bootstrap_counts must contain n-by-boots multinomial counts", call. = FALSE)
  }
  estimates <- rep(NA_real_, boots)
  diagnostics <- vector("list", boots)
  for (b in seq_len(boots)) {
    estimates[b] <- tryCatch({
      multiplicity <- if (is.null(bootstrap_counts)) {
        as.numeric(stats::rmultinom(1L, n, rep(1 / n, n)))
      } else bootstrap_counts[, b]
      if (sum(multiplicity[Z == 0]) == 0 || sum(multiplicity[Z == 1]) == 0)
        stop("both treatment arms must have positive bootstrap multiplicity")
      if (use.bias.correction) {
        nuisance <- estimate_scores(Y, point_est$X_internal, Z, w, ps, pg,
                                    model.ps, model.pg, sampling, point_est$estimand,
                                    multiplicity, use.bias.correction = TRUE)
        diagnostics[[b]] <- nuisance$ps_diagnostics
        q0 <- nuisance$q0
        q1 <- nuisance$q1
      } else {
        # Raw matching has no augmentation; refitting unused scores cannot
        # change a fixed-match raw replicate.
        q0 <- q1 <- rep(0, n)
        diagnostics[[b]] <- list(source = "not_refitted_raw_fixed_match")
      }
      wdsm_replicate(Y, Z, w, point_est$K, multiplicity, q0, q1, point_est$estimand)
    }, error = function(e) stop(sprintf("%s %s bootstrap replicate %d/%d failed: %s",
              sampling, point_est$estimand, b, boots, conditionMessage(e)), call. = FALSE))
  }
  if (any(!is.finite(estimates))) stop("All requested bootstrap draws must be finite")
  variance <- mean((estimates - mean(estimates))^2)
  se <- sqrt(variance)
  ci <- point_est$estimate + c(-1, 1) * stats::qnorm(1 - alpha / 2) * se
  if (any(!is.finite(c(variance, se, ci)))) stop("Nonfinite replication uncertainty")
  list(variance = variance, se = se, ci = ci, boot_estimates = estimates,
       ps_diagnostics = diagnostics)
}
