#!/usr/bin/env Rscript
# Fixed validation grid; no adaptive tuning, projection, or estimator simulation.
# Usage: Rscript validation/check_overlap_boundary.R draws new_output_directory
# Suggested bounded pilot: 20000 draws; maximum accepted here: 250000.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: draws new_output_directory")
draws <- suppressWarnings(as.numeric(args[[1L]]))
if (length(draws) != 1L || !is.finite(draws) || draws < 2 ||
    draws > 250000 || draws != floor(draws)) stop("draws must be an integer from 2 to 250000")
draws <- as.integer(draws)
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = TRUE)
if (length(script) != 1L) stop("Run with Rscript")
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
output <- args[[2L]]
if (file.exists(output) || dir.exists(output)) stop("Use a new output directory; existing results are never overwritten")
if (!dir.create(output, recursive = TRUE, showWarnings = FALSE)) stop("Unable to create output directory")
output <- normalizePath(output, mustWork = TRUE)
source_paths <- file.path(root, c("R/wm_geometry.R", "R/wm_geometry_alpha.R",
  "R/wm_geometry_overlap.R", "simulations/benchmarks.R"))
angle_path <- file.path(root, "validation", "geometry_angle_reference.R")
env <- new.env(parent = globalenv())
for (path in source_paths) sys.source(path, envir = env)
angle <- new.env(parent = globalenv())
angle$env <- env
sys.source(angle_path, envir = angle)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
config <- rbind(data.frame(M = 1:3, d = 1L),
  expand.grid(M = c(1L, 3L), d = c(2L, 3L, 5L)), data.frame(M = 2L, d = 2L))
config$draws <- draws
config$seed <- 928800L + seq_len(nrow(config))
config$case_id <- sprintf("M%d_d%d", config$M, config$d)
utils::write.csv(config, file.path(output, "config.csv"), row.names = FALSE)
code_paths <- c(script, source_paths, angle_path)
metadata <- list(status = "started", started = Sys.time(), config = config,
  code_md5 = tools::md5sum(code_paths),
  config_md5 = tools::md5sum(file.path(output, "config.csv")),
  immutable_sha256_manifest = Sys.getenv("WDSM_RUN_MANIFEST", unset = NA_character_),
  session = utils::sessionInfo(), host = Sys.info(),
  resources = Sys.getenv(c("SLURM_JOB_ID", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE",
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"), unset = NA_character_),
  precision_goal = "Propagated geometry MCSE <= 1% of positive benchmark variance; diagnostic goal, never a value-selection rule.",
  comparison_policy = paste("Exact 1d component checks; independent two-angle alpha quadrature;",
    "separately seeded direct-alpha MC; original full-position MC only at predefined low-dimensional cases.",
    "Standardized discrepancies are numerical diagnostics, not rigorous error guarantees.",
    "No draws or methods are changed in response to results."))
saveRDS(metadata, file.path(output, "metadata.rds"))
records <- comparisons <- benchmark_rows <- list()
objects <- stats::setNames(vector("list", nrow(config)), config$case_id)
started <- proc.time()[["elapsed"]]
append_table <- function(x, name) {
  path <- file.path(output, name)
  present <- file.exists(path)
  utils::write.table(x, path, sep = ",", row.names = FALSE, append = present,
    col.names = !present, na = "NA", qmethod = "double")
}
standardized <- function(difference, error) {
  result <- rep(NA_real_, length(difference))
  take <- is.finite(error) & error > 0
  result[take] <- difference[take] / error[take]
  result
}
for (index in seq_len(nrow(config))) {
  cfg <- config[index, ]
  case_started <- proc.time()[["elapsed"]]
  result <- tryCatch(env$.wm_geometry_overlap_dispatch(cfg$M, cfg$d, cfg$draws,
    cfg$seed, 1024L, exact_1d = FALSE), error = function(e) e)
  ok <- inherits(result, "wm_geometry")
  status <- data.frame(case_id = cfg$case_id, M = cfg$M, d = cfg$d, draws = draws,
    seed = cfg$seed, status = if (ok) "completed" else "error",
    elapsed_seconds = proc.time()[["elapsed"]] - case_started,
    error = if (ok) "" else conditionMessage(result))
  records[[index]] <- status
  append_table(status, "runs.csv")
  if (!ok) next
  objects[[cfg$case_id]] <- result
  saveRDS(result, file.path(output, paste0(cfg$case_id, ".rds")))
  append_table(data.frame(case_id = cfg$case_id, M = cfg$M, d = cfg$d,
    overlap = result$overlap, beta = unname(result$beta),
    mcse = unname(result$beta_mcse), component_precision = result$component_precision,
    nonzero_draws = unname(result$diagnostics$nonzero_draws),
    maximum_draw_value = unname(result$diagnostics$maximum_draw_value)), "beta.csv")
  term <- expand.grid(overlap = result$overlap, family = colnames(result$term_estimates),
                      stringsAsFactors = FALSE)
  append_table(data.frame(case_id = cfg$case_id, M = cfg$M, d = cfg$d, term,
    estimate = as.vector(result$term_estimates), mcse = as.vector(result$term_mcse),
    precision = as.vector(result$term_precision),
    nonzero_draws = as.vector(result$diagnostics$term_nonzero_draws),
    maximum_draw_value = as.vector(result$diagnostics$term_maximum_draw_value)), "boundary_terms.csv")
  if (cfg$d == 1L) {
    exact <- env$wm_geometry(cfg$M, 1L)
    delta <- unname(result$beta - exact$beta)
    z <- standardized(delta, unname(result$beta_mcse))
    append_table(data.frame(case_id = cfg$case_id, comparison = "exact_1d_beta",
      component = as.character(result$overlap), estimate = unname(result$beta),
      reference = unname(exact$beta), combined_mcse = unname(result$beta_mcse),
      standardized_difference = z, reference_error_scope = "exact_closed_form",
      diagnostic = ifelse(is.na(z), "unresolved_zero_empirical_mcse",
        ifelse(abs(z) > 6, "review_discrepancy", "within_six_empirical_mcse"))), "comparisons.csv")
    reference <- list(alpha = exact$alpha, outer_absolute_error = 0,
                       maximum_observed_inner_error = 0, evaluations = 0L,
                       error_scope = "Exact one-dimensional formula")
  } else {
    reference <- angle$alpha_angles(cfg$M, cfg$d)
  }
  saveRDS(reference, file.path(output, paste0(cfg$case_id, "_alpha_reference.rds")))
  z <- standardized(result$alpha - reference$alpha, result$alpha_mcse)
  append_table(data.frame(case_id = cfg$case_id, M = cfg$M, d = cfg$d,
    alpha = result$alpha, alpha_mcse = result$alpha_mcse,
    alpha_reference = reference$alpha, standardized_difference = z,
    outer_reference_error = reference$outer_absolute_error,
    largest_observed_inner_error = reference$maximum_observed_inner_error,
    reference_error_scope = reference$error_scope,
    precision_status = result$precision_status,
    diagnostic = if (is.na(z)) "unresolved_zero_empirical_mcse" else
      if (abs(z) > 6) "review_discrepancy" else "within_six_empirical_mcse"), "alpha.csv")
  if (cfg$d > 1L) {
    independent <- env$wm_geometry_alpha(cfg$M, cfg$d, draws = draws,
                                        seed = cfg$seed + 10000L, chunk_size = 1024L)
    saveRDS(independent, file.path(output, paste0(cfg$case_id, "_alpha_mc.rds")))
    error <- sqrt(result$alpha_mcse^2 + independent$alpha_mcse^2)
    z <- standardized(result$alpha - independent$alpha, error)
    append_table(data.frame(case_id = cfg$case_id, comparison = "independently_seeded_alpha_mc",
      component = "sum_beta", estimate = result$alpha, reference = independent$alpha,
      combined_mcse = error, standardized_difference = z,
      reference_error_scope = "Independent empirical Monte Carlo errors added in variance",
      diagnostic = if (is.na(z)) "unresolved_zero_empirical_mcse" else
        if (abs(z) > 6) "review_discrepancy" else "within_six_empirical_mcse"), "comparisons.csv")
  }
  if ((cfg$M == 1L && cfg$d %in% c(2L, 3L)) || (cfg$M == 2L && cfg$d == 2L)) {
    full <- env$wm_geometry(cfg$M, cfg$d, draws = min(draws, 20000L),
                            seed = cfg$seed + 20000L, chunk_size = 1024L)
    saveRDS(full, file.path(output, paste0(cfg$case_id, "_full_position.rds")))
    error <- sqrt(result$beta_mcse^2 + full$beta_mcse^2)
    z <- standardized(result$beta - full$beta, error)
    resolved <- full$beta > 0 & full$beta_mcse / full$beta <= 0.05 &
      !grepl("^unresolved", full$component_precision)
    append_table(data.frame(case_id = cfg$case_id, comparison = "independently_seeded_full_position",
      component = as.character(result$overlap), estimate = unname(result$beta),
      reference = unname(full$beta), combined_mcse = unname(error),
      standardized_difference = unname(z),
      reference_error_scope = "Empirical MCSE; full-position reference requires relative MCSE <= 5%",
      diagnostic = ifelse(!resolved, "full_position_reference_unresolved",
        ifelse(is.na(z), "unresolved_zero_empirical_mcse",
          ifelse(abs(z) > 6, "review_discrepancy", "within_six_empirical_mcse")))), "comparisons.csv")
  }
  for (design in c("strong", "weak")) {
    benchmark <- if (design == "strong") env$wm_benchmark_strong else env$wm_benchmark_weak
    table <- benchmark(cfg$M, result$beta, result$covariance)$table
    positive <- is.finite(table$root_n_variance) & table$root_n_variance > 0
    relative <- rep(NA_real_, nrow(table))
    relative[positive] <- table$geometry_mcse[positive] / table$root_n_variance[positive]
    unresolved <- grepl("^unresolved", result$precision_status)
    table$geometry_relative_mcse <- relative
    table$geometry_precision_goal <- ifelse(!positive, "variance_not_established",
      ifelse(unresolved, "unresolved_geometry",
        ifelse(relative <= 0.01, "within_1_percent_goal", "additional_precision_needed")))
    append_table(data.frame(case_id = cfg$case_id, M = cfg$M, d = cfg$d,
      design = design, table), "benchmark_precision.csv")
  }
  cat("Completed", cfg$case_id, "; alpha =", format(result$alpha, digits = 10),
      "+/-", format(result$alpha_mcse, digits = 5), ";", result$precision_status, "\n")
}
saveRDS(objects, file.path(output, "geometry.rds"))
metadata$completed <- Sys.time()
metadata$elapsed_seconds <- proc.time()[["elapsed"]] - started
metadata$runs <- do.call(rbind, records)
metadata$status <- if (any(metadata$runs$status != "completed")) "completed_with_errors" else "completed"
metadata$code_md5_at_completion <- tools::md5sum(code_paths)
metadata$code_unchanged <- identical(metadata$code_md5, metadata$code_md5_at_completion)
saveRDS(metadata, file.path(output, "metadata.rds"))
if (!metadata$code_unchanged) stop("Validation code changed during the run; inspect provenance")
if (metadata$status != "completed") stop("Some cases failed; preserve partial results and inspect runs.csv")
cat("Boundary-rank validation records saved. Precision goals were reported, not enforced by masking or projection.\n")
