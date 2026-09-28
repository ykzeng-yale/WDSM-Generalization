# Base-R deterministic benchmark-join checks; no integration or simulations.
# Rscript code-release/validation/check_join_benchmarks.R
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this validation with Rscript")
script <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
simulation_dir <- file.path(dirname(script), "..", "simulations")
source(file.path(simulation_dir, "benchmarks.R"))
source(file.path(simulation_dir, "join_benchmarks.R"))
equal <- function(actual, expected, tolerance = 1e-11) {
  result <- all.equal(unname(actual), unname(expected), tolerance = tolerance,
                     check.attributes = FALSE)
  if (!isTRUE(result)) stop(paste(result, collapse = "; "), call. = FALSE)
}
fails <- function(expression, pattern) {
  message <- tryCatch({force(expression); NULL}, error = conditionMessage)
  if (is.null(message) || !grepl(pattern, message)) stop("Expected error: ", pattern)
}
fixture <- function(design = "strong", method = "self_normalized", estimand = "PATE",
                    correction = "oracle", d = 1, M = 1, n = 200,
                    scenario = paste(design, d, M, n, sep = "_")) {
  unsupported <- design == "weak" && method == "self_normalized"
  data.frame(scenario = scenario, n = n, d = d, M = M, design = design,
    method = method, estimand = estimand, correction = correction,
    target = if (design == "strong") 1.5 else 0,
    requested = 40, completed = 38, failed = 2,
    inference_status = if (correction == "raw") "point_estimate_only" else
      if (unsupported) "assumption_failure_diagnostic" else "within_declared_theorem_scope",
    empirical_variance = 0.03, empirical_variance_mcse_asymptotic = 0.003,
    mean_root_n_variance = if (correction == "raw") NA_real_ else 6,
    mean_root_n_variance_mcse = if (correction == "raw") NA_real_ else 0.2,
    failure_rate = 0.05, coverage = if (correction == "raw") NA_real_ else 0.92,
    stringsAsFactors = FALSE)
}
geometry <- function(M = 1L, d = 2L, beta = 1.6, covariance = matrix(1e-4, 1, 1)) {
  structure(list(M = M, d = d, beta = beta, covariance = covariance,
    beta_mcse = if (is.null(covariance)) rep(NA_real_, M) else sqrt(diag(covariance)),
    method = "uniform_sphere_monte_carlo", precision_status = "monte_carlo_estimate",
    component_precision = rep("monte_carlo_estimate", M), draws = 10000,
    seed = 928102, diagnostics = list(nonzero_draws = rep(8000, M))),
    class = c("wm_geometry", "list"))
}

# Hand-derived rational M=1 limits independently check row selection and scaling.
known <- c(29695 / 4704, 1565 / 192, 42265 / 4704, 8555 / 768)
rows <- do.call(rbind, lapply(seq_len(4), function(k) {
  method <- if (k <= 2) "self_normalized" else "stabilized"
  estimand <- if (k %% 2 == 1) "PATE" else "PATT"
  x <- fixture(method = method, estimand = estimand)
  x$empirical_variance <- (known[k] + 0.4) / x$n
  x$empirical_variance_mcse_asymptotic <- 0.1 / x$n
  x$mean_root_n_variance <- known[k] + 0.2
  x
}))
joined <- wm_join_benchmarks(rows)
stopifnot(identical(joined[names(rows)], rows))
equal(joined$benchmark_root_n_variance, known)
equal(joined$benchmark_geometry_mcse, rep(0, 4))
equal(joined$empirical_root_n_variance, known + 0.4)
equal(joined$empirical_minus_benchmark_root_n, rep(0.4, 4))
equal(joined$empirical_minus_benchmark_mcse_independent, rep(0.1, 4))
equal(joined$empirical_benchmark_standardized_discrepancy, rep(4, 4))
equal(joined$fitted_benchmark_standardized_discrepancy, rep(1, 4))
stopifnot(all(joined$geometry_precision_goal_met), all(joined$failed == 2),
          all(joined$geometry_draws == 0))

weak <- rbind(fixture("weak", "self_normalized", "PATE"),
              fixture("weak", "self_normalized", "PATT"),
              fixture("weak", "stabilized", "PATE"),
              fixture("weak", "stabilized", "PATT"))
w <- wm_join_benchmarks(weak)
equal(w$benchmark_probability_limit, c(-1/24, 1/6, 0, 0))
equal(w$benchmark_root_n_variance, c(NA_real_, NA_real_, 167/96, 221/96))
stopifnot(all(is.na(w$empirical_benchmark_standardized_discrepancy[1:2])),
          all(w$benchmark_scope[1:2] == "original_weak_variance_not_established"))
raw <- wm_join_benchmarks(fixture("weak", "self_normalized", "PATE", "raw"))
equal(raw$benchmark_probability_limit, -1/24)
stopifnot(is.na(raw$benchmark_root_n_variance), is.na(raw$benchmark_geometry_mcse),
          raw$benchmark_scope == "raw_point_only", raw$failed == 2)
polynomial <- wm_join_benchmarks(fixture(correction = "polynomial"))
equal(polynomial$benchmark_root_n_variance, known[1])

# In M=1 strong original PATE, the beta coefficient is exactly 235/294.
g <- geometry()
input <- fixture(d = 2)
j <- wm_join_benchmarks(input, g)
c1 <- 235/294
expected_value <- known[1] + c1 * (1.6 - 1.5)
expected_geometry_se <- c1 * 0.01
equal(j$benchmark_root_n_variance, expected_value)
equal(j$benchmark_geometry_mcse, expected_geometry_se)
equal(j$empirical_minus_benchmark_mcse_independent,
      sqrt((200 * 0.003)^2 + expected_geometry_se^2))
equal(j$empirical_benchmark_standardized_discrepancy,
      (6 - expected_value) / sqrt(0.6^2 + expected_geometry_se^2))
stopifnot(j$geometry_precision_goal_met)

# Constant coefficients for weak stabilized PATE make the full covariance
# contraction c^2*sum(Cov), including all off-diagonal terms, independently visible.
covariance <- matrix(c(.0009, .0001, .0002, .0001, .0016, .0003,
                      .0002, .0003, .0025), 3, 3)
g3 <- geometry(3L, 2L, c(3, 4, 4), covariance)
j3 <- wm_join_benchmarks(fixture("weak", "stabilized", d = 2, M = 3), g3)
equal(j3$benchmark_geometry_mcse, (119 / 3456) * sqrt(sum(covariance)))

missing <- wm_join_benchmarks(input)
stopifnot(is.na(missing$benchmark_root_n_variance),
          missing$comparison_status == "geometry_not_provided", missing$failed == 2)
failed <- wm_join_benchmarks(input, list(M1_d2 = NULL))
stopifnot(failed$geometry_status == "geometry_run_failed", is.na(failed$benchmark_root_n_variance))
unresolved <- g
unresolved$beta <- 0
unresolved$covariance[,] <- 0
unresolved$precision_status <- "unresolved_components"
unresolved$component_precision <- "unresolved_zero_hits"
unresolved$diagnostics$nonzero_draws <- 0
u <- wm_join_benchmarks(input, unresolved)
stopifnot(is.finite(u$benchmark_root_n_variance), u$geometry_reported_mcse == 0,
          is.na(u$benchmark_geometry_mcse), is.na(u$empirical_minus_benchmark_mcse_independent),
          is.na(u$empirical_benchmark_standardized_discrepancy), u$geometry_below_jensen,
          u$comparison_status == "geometry_unresolved")
zero_variance <- g
zero_variance$covariance[,] <- 0
z <- wm_join_benchmarks(input, zero_variance)
stopifnot(z$comparison_status == "geometry_unresolved", is.na(z$benchmark_geometry_mcse))
unknown_cov <- g
unknown_cov$covariance <- NULL
stopifnot(wm_join_benchmarks(input, unknown_cov)$comparison_status == "geometry_unresolved")
coarse <- g
coarse$covariance[,] <- 4
cj <- wm_join_benchmarks(input, coarse)
stopifnot(cj$comparison_status == "geometry_precision_goal_not_met",
          is.finite(cj$benchmark_geometry_mcse),
          is.finite(cj$empirical_minus_benchmark_mcse_independent),
          is.na(cj$empirical_benchmark_standardized_discrepancy))

all_failed <- fixture()
all_failed$completed <- 0; all_failed$failed <- all_failed$requested
all_failed$empirical_variance <- all_failed$empirical_variance_mcse_asymptotic <- NA_real_
all_failed$mean_root_n_variance <- all_failed$mean_root_n_variance_mcse <- NA_real_
af <- wm_join_benchmarks(all_failed)
stopifnot(af$failed == 40, af$completed == 0,
          is.finite(af$benchmark_root_n_variance), is.na(af$empirical_root_n_variance),
          is.na(af$fitted_benchmark_standardized_discrepancy))

# Corrupt input is rejected rather than silently dropped, pooled or relabeled.
fails(wm_join_benchmarks(rbind(input, input)), "duplicate summary")
bad <- input; bad$target <- 0
fails(wm_join_benchmarks(bad), "target disagrees")
bad <- input; bad$inference_status <- "point_estimate_only"
fails(wm_join_benchmarks(bad), "inference label")
bad <- input; bad$failed <- 0
fails(wm_join_benchmarks(bad), "counts disagree")
bad <- rbind(input, transform(input, method = "stabilized", n = 201))
fails(wm_join_benchmarks(bad), "inconsistent scenario")
fails(wm_join_benchmarks(input, list(M3_d2 = g)), "key disagrees")
fails(wm_join_benchmarks(input, list(g, g)), "duplicate geometry")
bad_geometry <- g; bad_geometry$covariance[,] <- -1
fails(wm_join_benchmarks(input, bad_geometry), "positive semidefinite")
fails(wm_join_benchmarks(joined), "already contains")
cat("Benchmark join deterministic checks passed. No simulation or geometry integration was run.\n")
