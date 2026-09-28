# Fit only nuisance components required by the estimand. Supplied scores stay fixed.
estimate_scores <- function(Y, X, Z, sw = NULL, ps = NULL, pg = NULL,
                            model.ps = NULL, model.pg = NULL,
                            sampling = "retrospective", estimand = "PATE",
                            multiplicity = rep(1, length(Y)),
                            use.bias.correction = TRUE) {
  n <- length(Y)
  df <- data.frame(Z = Z, Y = Y, X, check.names = FALSE)
  if (is.null(sw)) sw <- rep(1, n)
  if (is.null(ps)) {
    ps_weights <- multiplicity * if (sampling == "retrospective") sw else 1
    fit <- wdsm_fit_ps(model.ps, df, ps_weights)
    ps <- fit$probability
    diagnostics <- fit$diagnostics
  } else {
    diagnostics <- list(source = "supplied_fixed", probability_range = range(ps))
  }
  result <- list(ps = ps, e_score = ps, ps_diagnostics = diagnostics,
                 psi0 = NULL, psi1 = NULL, q0 = rep(0, n), q1 = rep(0, n),
                 D0 = NULL, D1 = NULL)
  for (arm in if (estimand == "PATE") 0:1 else 0L) {
    idx <- which(Z == arm)
    prognostic <- if (is.null(pg)) {
      wdsm_fit_regression(model.pg, df[idx, , drop = FALSE], df,
                          multiplicity[idx], paste("Arm", arm, "prognostic"))
    } else pg[, arm + 1L]
    score <- wdsm_standardize(cbind(ps, prognostic), multiplicity)
    if (use.bias.correction) {
      sieve <- wdsm_sieve(score)
      regression_data <- data.frame(Y = Y[idx], sieve[idx, , drop = FALSE])
      result[[paste0("q", arm)]] <- wdsm_fit_regression(
        Y ~ ., regression_data, sieve, multiplicity[idx] * sw[idx],
        paste("Arm", arm, "bias correction"))
    }
    result[[paste0("psi", arm)]] <- as.numeric(prognostic)
    result[[paste0("D", arm)]] <- score
  }
  result
}
