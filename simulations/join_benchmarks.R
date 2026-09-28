#!/usr/bin/env Rscript
# Rscript simulations/join_benchmarks.R summary.csv output.csv [geometry.rds ...]
# Source benchmarks.R first when using wm_join_benchmarks() as a function.
# Geometry and estimator simulations must use independent random streams.

.wm_jb_geometry_list <- function(geometry) {
  if (is.null(geometry)) return(list())
  if (inherits(geometry, "wm_geometry") ||
      (is.list(geometry) && all(c("M", "d", "beta") %in% names(geometry)))) {
    geometry <- list(geometry)
  }
  if (!is.list(geometry)) stop("geometry must be a geometry result or list of results")
  result <- list()
  for (i in seq_along(geometry)) {
    g <- geometry[[i]]
    label <- if (is.null(names(geometry))) "" else names(geometry)[i]
    if (is.na(label)) stop("geometry keys cannot be missing")
    if (is.null(g)) {
      parsed <- regmatches(label, regexec("^M([1-9][0-9]*)_d([1-9][0-9]*)$", label))[[1L]]
      if (length(parsed) != 3L) stop("failed geometry entries require an M{M}_d{d} key")
      g <- list(M = as.numeric(parsed[2L]), d = as.numeric(parsed[3L]), failed = TRUE)
    }
    if (!is.list(g)) stop("invalid geometry result")
    M <- .wm_bm_integer(g$M)
    d <- .wm_bm_integer(g$d)
    key <- sprintf("M%d_d%d", M, d)
    if (nzchar(label) && !identical(label, key)) stop("geometry key disagrees with M and d")
    if (!is.null(result[[key]])) stop("duplicate geometry key: ", key)
    if (isTRUE(g$failed)) {
      result[[key]] <- g
      next
    }
    checked <- .wm_bm_geometry(M, g$beta, g$covariance)
    if (d == 1L) {
      if (!isTRUE(all.equal(unname(checked$beta), unname(wm_beta_1d(M)),
                           tolerance = 1e-12))) stop("incorrect one-dimensional geometry")
    }
    result[[key]] <- g
  }
  result
}

.wm_jb_geometry <- function(M, d, geometry) {
  key <- sprintf("M%d_d%d", M, d)
  if (d == 1L) return(list(beta = wm_beta_1d(M), covariance = matrix(0, M, M),
    available = TRUE, resolved = TRUE, status = "exact_1d_closed_form",
    precision_status = "exact_formula_evaluated_in_floating_point",
    draws = 0, seed = NA_real_, key = key))
  g <- geometry[[key]]
  if (is.null(g) || isTRUE(g$failed)) {
    status <- if (is.null(g)) "not_provided" else "geometry_run_failed"
    return(list(beta = wm_beta_1d(M), covariance = NULL,
      available = FALSE, resolved = FALSE, status = status,
      precision_status = status, draws = NA_real_, seed = NA_real_, key = key))
  }
  status <- "unresolved_geometry_precision"
  resolved <- FALSE
  # A zero-hit or numerically zero-variance component cannot establish precision,
  # even when the API returned a numerical zero empirical standard error.
  if (!is.null(g$covariance) &&
      identical(g$precision_status, "monte_carlo_estimate") &&
      is.character(g$component_precision) && length(g$component_precision) == M &&
      !anyNA(g$component_precision) && all(g$component_precision == "monte_carlo_estimate") &&
      all(diag(g$covariance) > 0)) {
    hits <- g$diagnostics$nonzero_draws
    if (is.numeric(hits) && length(hits) == M && all(is.finite(hits)) && all(hits > 0)) {
      resolved <- TRUE
      status <- "monte_carlo_precision_estimated"
    }
  }
  scalar_or_na <- function(x) if (is.numeric(x) && !is.complex(x) &&
    length(x) == 1L && is.finite(x)) as.numeric(x) else NA_real_
  list(beta = g$beta, covariance = g$covariance, available = TRUE,
    resolved = resolved, status = status,
    precision_status = if (is.character(g$precision_status) &&
      length(g$precision_status) == 1L && !is.na(g$precision_status))
      g$precision_status else "not_declared",
    draws = scalar_or_na(g$draws), seed = scalar_or_na(g$seed), key = key)
}

wm_join_benchmarks <- function(summary, geometry = NULL, geometry_relative_goal = 0.01) {
  if (!exists("wm_benchmark_strong", mode = "function") ||
      !exists("wm_benchmark_weak", mode = "function")) stop("Source benchmarks.R first")
  geometry_relative_goal <- .wm_bm_scalar(geometry_relative_goal,
    "geometry_relative_goal", 0, 1, TRUE, TRUE)
  required <- c("scenario", "n", "d", "M", "design", "method", "estimand", "correction",
    "target", "requested", "completed", "failed", "inference_status",
    "empirical_variance", "empirical_variance_mcse_asymptotic",
    "mean_root_n_variance", "mean_root_n_variance_mcse")
  if (!is.data.frame(summary) || nrow(summary) == 0L ||
      !all(required %in% names(summary)) || anyDuplicated(names(summary)))
    stop("summary must contain the required unique columns and at least one row")
  key_columns <- c("scenario", "method", "estimand", "correction")
  if (anyDuplicated(summary[key_columns])) stop("duplicate summary key")
  for (field in c("scenario", "design", "method", "estimand", "correction", "inference_status")) {
    if (!is.character(summary[[field]]) || anyNA(summary[[field]]) ||
        any(!nzchar(summary[[field]]))) stop("invalid summary field: ", field)
  }
  for (field in c("n", "d", "M", "requested", "completed", "failed")) {
    x <- summary[[field]]
    if (!is.numeric(x) || is.complex(x) || any(!is.finite(x)) ||
        any(x != floor(x)) || any(x > .Machine$integer.max) ||
        any(x < if (field %in% c("completed", "failed")) 0 else 1))
      stop("invalid summary integer: ", field)
  }
  if (any(summary$completed + summary$failed != summary$requested))
    stop("requested/completed/failed counts disagree")
  for (scenario in unique(summary$scenario)) {
    block <- summary[summary$scenario == scenario, , drop = FALSE]
    for (field in c("n", "d", "M", "design")) {
      if (length(unique(block[[field]])) != 1L) stop("inconsistent scenario configuration")
    }
  }
  if (any(!summary$design %in% c("strong", "weak")) ||
      any(!summary$method %in% c("self_normalized", "stabilized")) ||
      any(!summary$estimand %in% c("PATE", "PATT")) ||
      any(!summary$correction %in% c("oracle", "polynomial", "raw")))
    stop("unsupported design, method, estimand or correction")
  for (field in c("target", "empirical_variance", "empirical_variance_mcse_asymptotic",
                  "mean_root_n_variance", "mean_root_n_variance_mcse")) {
    x <- summary[[field]]
    if (!is.numeric(x) || is.complex(x) || any(!is.na(x) & !is.finite(x)) ||
        (field != "target" && any(x < 0, na.rm = TRUE))) stop("invalid numeric summary field: ", field)
  }
  if (anyNA(summary$target)) stop("target must be finite")
  geometry <- .wm_jb_geometry_list(geometry)
  cache <- list()
  additions <- vector("list", nrow(summary))
  for (i in seq_len(nrow(summary))) {
    row <- summary[i, , drop = FALSE]
    M <- .wm_bm_integer(row$M); d <- .wm_bm_integer(row$d)
    g <- .wm_jb_geometry(M, d, geometry)
    cache_key <- paste(row$design, g$key, sep = ":")
    if (is.null(cache[[cache_key]])) {
      fun <- if (row$design == "strong") wm_benchmark_strong else wm_benchmark_weak
      cache[[cache_key]] <- fun(M, g$beta, g$covariance)
    }
    benchmark <- cache[[cache_key]]
    b <- benchmark$table[benchmark$table$method == row$method &
      benchmark$table$estimand == row$estimand, , drop = FALSE]
    if (nrow(b) != 1L || abs(row$target - b$target) > 1e-12)
      stop("summary target disagrees with the fixed design benchmark")
    if ("probability_limit" %in% names(row) &&
        (!is.finite(row$probability_limit) || abs(row$probability_limit - b$probability_limit) > 1e-12))
      stop("summary probability limit disagrees with the benchmark")
    raw <- row$correction == "raw"
    supported <- !raw && b$variance_scope == "bias_corrected"
    expected_inference <- if (raw) "point_estimate_only" else
      if (!supported) "assumption_failure_diagnostic" else "within_declared_theorem_scope"
    if (row$inference_status != expected_inference) stop("summary inference label disagrees with benchmark scope")
    scope <- if (raw) "raw_point_only" else if (supported) "bias_corrected" else
      "original_weak_variance_not_established"
    value <- if (supported && g$available) b$root_n_variance else NA_real_
    reported_se <- if (supported && g$available) b$geometry_mcse else NA_real_
    meaningful <- supported && g$available && g$resolved && is.finite(reported_se) &&
      (d == 1L || reported_se > 0)
    geometry_se <- if (meaningful) reported_se else NA_real_
    relative_se <- if (meaningful && value > 0) geometry_se / value else NA_real_
    precision_met <- meaningful && is.finite(relative_se) && relative_se <= geometry_relative_goal
    comparison_status <- if (!supported) "not_applicable" else
      if (!g$available) "geometry_not_provided" else if (!meaningful) "geometry_unresolved" else
        if (!precision_met) "geometry_precision_goal_not_met" else "monte_carlo_diagnostic_available"
    empirical <- row$n * row$empirical_variance
    empirical_se <- row$n * row$empirical_variance_mcse_asymptotic
    fitted <- row$mean_root_n_variance
    fitted_se <- row$mean_root_n_variance_mcse
    combined <- function(se) if (is.finite(se) && is.finite(geometry_se))
      sqrt(se^2 + geometry_se^2) else NA_real_
    empirical_combined <- combined(empirical_se)
    fitted_combined <- combined(fitted_se)
    empirical_difference <- empirical - value
    fitted_difference <- fitted - value
    standardize <- function(difference, se) if (precision_met && is.finite(difference) &&
      is.finite(se) && se > 0) difference / se else NA_real_
    additions[[i]] <- data.frame(
      benchmark_scope = scope, benchmark_probability_limit = b$probability_limit,
      benchmark_target_bias = b$bias, benchmark_root_n_variance = value,
      benchmark_geometry_mcse = geometry_se, geometry_reported_mcse = reported_se,
      geometry_key = g$key, geometry_status = g$status,
      geometry_precision_status = g$precision_status, geometry_draws = g$draws,
      geometry_seed = g$seed, geometry_relative_mcse = relative_se,
      geometry_relative_goal = geometry_relative_goal,
      geometry_precision_goal_met = if (supported) precision_met else NA,
      geometry_below_jensen = if (g$available) !benchmark$alpha_jensen_bound_satisfied else NA,
      comparison_status = comparison_status,
      empirical_root_n_variance = empirical,
      empirical_root_n_variance_mcse_asymptotic = empirical_se,
      empirical_minus_benchmark_root_n = empirical_difference,
      empirical_minus_benchmark_mcse_independent = empirical_combined,
      empirical_benchmark_standardized_discrepancy = standardize(empirical_difference, empirical_combined),
      fitted_minus_benchmark_root_n = fitted_difference,
      fitted_minus_benchmark_mcse_independent = fitted_combined,
      fitted_benchmark_standardized_discrepancy = standardize(fitted_difference, fitted_combined),
      empirical_to_benchmark_ratio = if (is.finite(value) && value > 0) empirical / value else NA_real_,
      fitted_to_benchmark_ratio = if (is.finite(value) && value > 0) fitted / value else NA_real_,
      comparison_interpretation = "Independent-stream Monte Carlo diagnostics; finite-n departures from an asymptotic limit are possible; no hypothesis-test or coverage certification",
      stringsAsFactors = FALSE)
  }
  added <- do.call(rbind, additions)
  if (any(names(added) %in% names(summary))) stop("summary already contains benchmark output columns")
  cbind(summary, added)
}

.wm_join_benchmarks_main <- function(args) {
  if (length(args) < 2L) stop("usage: summary.csv output.csv [geometry.rds ...]")
  output <- args[2L]
  sidecar <- paste0(output, ".metadata.rds")
  if (file.exists(output) || file.exists(sidecar)) stop("output or metadata already exists")
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(script_arg) != 1L) stop("Run this interface with Rscript")
  script <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
  benchmark_file <- file.path(dirname(script), "benchmarks.R")
  source(benchmark_file)
  geometry_files <- if (length(args) > 2L) args[-c(1L, 2L)] else character()
  if (anyDuplicated(normalizePath(geometry_files))) stop("duplicate geometry input file")
  geometry <- list()
  for (file in geometry_files) {
    block <- .wm_jb_geometry_list(readRDS(file))
    if (any(names(block) %in% names(geometry))) stop("duplicate geometry key across files")
    geometry <- c(geometry, block)
  }
  summary <- utils::read.csv(args[1L], stringsAsFactors = FALSE)
  joined <- wm_join_benchmarks(summary, geometry)
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  metadata <- list(created = Sys.time(), session = utils::sessionInfo(),
    input_md5 = tools::md5sum(c(args[1L], geometry_files)),
    code_md5 = tools::md5sum(c(script, benchmark_file)),
    geometry_relative_goal = 0.01,
    independence_contract = "Geometry integration streams are independent of estimator-simulation streams",
    interpretation = "Descriptive Monte Carlo discrepancies from asymptotic benchmarks, not hypothesis tests")
  saveRDS(metadata, sidecar)
  utils::write.csv(joined, output, row.names = FALSE, na = "NA")
  cat("Joined", nrow(joined), "summary rows; requested/completed/failed counts are unchanged.\n")
}

if (sys.nframe() == 0L) .wm_join_benchmarks_main(commandArgs(trailingOnly = TRUE))
