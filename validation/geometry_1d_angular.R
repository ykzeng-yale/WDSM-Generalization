# Independent reduction benchmark: the angular integrator uses d-ball unions
# and private-donor exclusions, while the reference uses exact 1d query cells.
# Run explicitly, e.g. Rscript validation/geometry_1d_angular.R 20000 output.csv.
# This script is not invoked by package loading or routine unit tests.
args <- commandArgs(trailingOnly = TRUE)
draws <- if (length(args)) as.numeric(args[[1L]]) else 20000
if (length(draws) != 1L || !is.finite(draws) || draws != floor(draws) ||
    draws < 2 || draws > 100000) stop("draws must be an integer from 2 to 100000")
output <- if (length(args) >= 2L) args[[2L]] else "geometry_1d_angular.csv"
records <- list()
for (M in 1:3) {
  seed <- 20260928L + M
  set.seed(seed)
  numerical <- WeightedMatching:::.wm_geometry_mc(M, 1L, as.integer(draws), 256L)
  exact <- WeightedMatching::wm_geometry(M, 1L)
  precision <- ifelse(numerical$nonzero == 0, "unresolved_zero_hits",
                      ifelse(numerical$beta_mcse == 0,
                             "unresolved_zero_empirical_variance", "monte_carlo_estimate"))
  standardized_error <- rep(NA_real_, M)
  positive_mcse <- numerical$beta_mcse > 0
  standardized_error[positive_mcse] <-
    (numerical$beta[positive_mcse] - exact$beta[positive_mcse]) /
    numerical$beta_mcse[positive_mcse]
  records[[M]] <- data.frame(M = M, d = 1L, overlap = 0:(M - 1L),
    seed = seed, draws = draws, exact = unname(exact$beta),
    estimate = numerical$beta, mcse = numerical$beta_mcse,
    nonzero_draws = numerical$nonzero, precision_status = precision,
    standardized_error = standardized_error,
    status = "Numerical reduction check; not a proof or estimator simulation")
}
result <- do.call(rbind, records)
utils::write.csv(result, output, row.names = FALSE)
print(result)
