#!/usr/bin/env Rscript
# Bounded universal-constant integration, not an estimator simulation.
# Usage: Rscript simulations/geometry_constants.R config.csv new_output_directory [source_package]
# Default grid: simulations/geometry_constants.csv. The third argument defaults
# to the package source enclosing this script; no installed package is required.
# No projection onto alpha >= M^2 and no precision-based deletion are applied.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) {
  stop("usage: config.csv new_output_directory [source_package]", call. = FALSE)
}
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this driver with Rscript", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source_root <- if (length(args) == 3L) args[[3L]] else
  file.path(dirname(script_path), "..")
source_root <- normalizePath(source_root, mustWork = TRUE)
geometry_source <- normalizePath(file.path(source_root, "R", "wm_geometry.R"), mustWork = TRUE)
config_path <- normalizePath(args[[1L]], mustWork = TRUE)
config <- utils::read.csv(config_path, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("M", "d", "draws", "seed")
if (!identical(names(config), required) || nrow(config) < 1L || nrow(config) > 32L) {
  stop("Configuration must have 1 to 32 rows and exactly M,d,draws,seed columns.", call. = FALSE)
}
integer_column <- function(x, minimum, maximum = .Machine$integer.max) {
  is.numeric(x) && !is.complex(x) && all(is.finite(x)) &&
    all(x >= minimum & x <= maximum & x == floor(x))
}
if (!integer_column(config$M, 1, 128) || !integer_column(config$d, 1) ||
    !integer_column(config$draws, 1, 1000000) || !integer_column(config$seed, 0) ||
    any(config$d > 1 & config$draws < 2)) {
  stop("Invalid integer configuration; M <= 128, draws <= 1000000, and numerical draws >= 2 are required.",
       call. = FALSE)
}
case_id <- sprintf("M%d_d%d", config$M, config$d)
if (anyDuplicated(case_id) || anyDuplicated(config$seed)) {
  stop("Each (M,d) and each seed must be unique in a run.", call. = FALSE)
}
coordinate_work <- ifelse(config$d == 1, 0,
  config$draws * as.double(config$d) * (3 * as.double(config$M)^2 + config$M) / 2)
max_dimension <- 2 * as.double(config$M) * config$d
if (any(!is.finite(coordinate_work)) || any(coordinate_work > 50000000) ||
    any(max_dimension[config$d > 1] > 1000000) || sum(coordinate_work) > 100000000) {
  stop("Configuration exceeds 50 million component coordinates per call or 100 million per run.",
       call. = FALSE)
}
output <- args[[2L]]
if (file.exists(output) || dir.exists(output)) {
  stop("Output path exists; use a new directory and preserve completed or partial results.", call. = FALSE)
}
if (!dir.create(output, recursive = TRUE, showWarnings = FALSE)) {
  stop("Unable to create output directory.", call. = FALSE)
}
output <- normalizePath(output, mustWork = TRUE)
if (!file.copy(config_path, file.path(output, "config.csv"), overwrite = FALSE)) {
  stop("Unable to retain the exact input configuration.", call. = FALSE)
}
env <- new.env(parent = globalenv())
sys.source(geometry_source, envir = env)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
chunk_size <- 2048L
manifest_link <- Sys.getenv("WDSM_RUN_MANIFEST", unset = NA_character_)
manifest_exists <- !is.na(manifest_link) && nzchar(manifest_link) &&
  file.exists(manifest_link) && !dir.exists(manifest_link)
code_paths <- c(driver = script_path, geometry = geometry_source)
started <- Sys.time()
elapsed_start <- proc.time()[["elapsed"]]
metadata <- list(
  schema_version = 1L, task = "universal_geometry_constants",
  status = "started", started = started, config = config, case_id = case_id,
  code_md5 = tools::md5sum(code_paths), config_file_md5 = tools::md5sum(config_path),
  retained_config_md5 = tools::md5sum(file.path(output, "config.csv")),
  immutable_sha256_manifest = manifest_link,
  linked_manifest_md5 = if (manifest_exists) tools::md5sum(manifest_link) else NULL,
  manifest_link_status = if (manifest_exists) "local_file_linked_and_hashed" else
    if (is.na(manifest_link) || !nzchar(manifest_link)) "not_supplied" else "external_link_retained",
  source_root = source_root, source_mode = "explicit_source_file",
  session = utils::sessionInfo(), rng_kind = RNGkind(), host = Sys.info(),
  resource_environment = Sys.getenv(c("WDSM_RUN_MANIFEST", "SLURM_JOB_ID", "SLURM_ARRAY_JOB_ID",
    "SLURM_ARRAY_TASK_ID", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE", "SLURM_MEM_PER_CPU",
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS"),
    unset = NA_character_),
  bounds = list(max_config_rows = 32L, max_M = 128L, max_requested_draws = 1000000L,
    max_component_coordinates_per_call = 50000000, max_component_coordinates_per_run = 100000000,
    requested_chunk_size = chunk_size, expected_component_coordinates = coordinate_work,
    total_expected_component_coordinates = sum(coordinate_work)),
  precision_goal = list(relative_to_positive_limit_variance = 0.01,
    expression = "sqrt(t(a) %*% Cov(beta_hat) %*% a) <= 0.01 * V",
    meaning = "Goal for propagated geometry MCSE in a specified positive asymptotic variance benchmark V with coefficient vector a.",
    status = "not_evaluated_without_benchmark_coefficients",
    policy = "A precision goal, not a guarantee, rejection rule, mask, or clipping/projection instruction."),
  output_schema = list(geometry_rds = "Named list keyed M{M}_d{d} of unmodified successful wm_geometry objects; failed cases remain NULL.",
    individual_rds = "One unmodified wm_geometry object per successful case, named M{M}_d{d}.rds.",
    covariance = "Covariance of beta estimates, including cross-component covariance; exact 1d CSV entries are zero, while original RDS covariance is NULL.",
    statuses = "API precision status is retained separately from additional zero-empirical-variance diagnostics."))
saveRDS(metadata, file.path(output, "metadata.rds"))

objects <- stats::setNames(vector("list", nrow(config)), case_id)
run_records <- list()
append_csv <- function(record, name) {
  path <- file.path(output, name)
  present <- file.exists(path)
  utils::write.table(record, path, sep = ",", row.names = FALSE, col.names = !present,
                     append = present, na = "NA", qmethod = "double")
}
for (index in seq_len(nrow(config))) {
  cfg <- config[index, , drop = FALSE]
  id <- case_id[[index]]
  warnings <- character()
  case_started <- Sys.time()
  clock_started <- proc.time()[["elapsed"]]
  result <- tryCatch(withCallingHandlers(
    env$wm_geometry(cfg$M, cfg$d, draws = cfg$draws, seed = cfg$seed, chunk_size = chunk_size),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  elapsed <- proc.time()[["elapsed"]] - clock_started
  success <- inherits(result, "wm_geometry")
  record <- data.frame(case_id = id, M = cfg$M, d = cfg$d, requested_draws = cfg$draws,
    seed = cfg$seed, status = if (success) "completed" else "error",
    method = if (success) result$method else NA_character_,
    precision_status = if (success) result$precision_status else "unavailable_error",
    started = format(case_started, tz = "UTC", usetz = TRUE), elapsed_seconds = elapsed,
    component_coordinates = coordinate_work[index], warning = paste(warnings, collapse = " | "),
    error = if (success) "" else conditionMessage(result), stringsAsFactors = FALSE)
  run_records[[index]] <- record
  append_csv(record, "runs.csv")
  if (!success) {
    saveRDS(list(case_id = id, config = cfg, message = conditionMessage(result),
                 warnings = warnings, elapsed_seconds = elapsed), file.path(output, paste0(id, ".error.rds")))
    next
  }
  objects[[id]] <- result
  saveRDS(result, file.path(output, paste0(id, ".rds")))
  exact <- identical(result$method, "exact_1d_closed_form")
  nonzero <- if (exact) rep(NA_real_, cfg$M) else unname(result$diagnostics$nonzero_draws)
  maxima <- if (exact) rep(NA_real_, cfg$M) else unname(result$diagnostics$maximum_draw_value)
  diagnostic <- if (exact) rep("exact_formula", cfg$M) else
    ifelse(nonzero == 0, "unresolved_zero_hits",
           ifelse(result$beta_mcse == 0, "unresolved_zero_empirical_variance", "monte_carlo_estimate"))
  beta <- data.frame(case_id = id, M = cfg$M, d = cfg$d, overlap = result$overlap,
    beta = unname(result$beta), beta_mcse = unname(result$beta_mcse),
    method = result$method, precision_status = result$precision_status,
    component_precision = result$component_precision, diagnostic_status = unname(diagnostic),
    requested_draws = result$requested_draws, draws = result$draws, seed = result$seed,
    chunk_size = result$chunk_size, nonzero_draws = nonzero, maximum_draw_value = maxima,
    rng_used = result$diagnostics$rng_used, stringsAsFactors = FALSE)
  alpha_diagnostic <- if (exact) "exact_formula" else
    if (any(grepl("^unresolved", diagnostic))) "unresolved_components" else
    if (result$alpha_mcse == 0) "unresolved_zero_empirical_variance" else "monte_carlo_estimate"
  covariance <- if (exact) matrix(0, cfg$M, cfg$M) else result$covariance
  # Record both equivalent MCSE calculations, without changing either result.
  covariance_alpha_mcse <- sqrt(max(0, sum(covariance)))
  alpha <- data.frame(case_id = id, M = cfg$M, d = cfg$d, alpha = result$alpha,
    alpha_mcse = result$alpha_mcse, method = result$method,
    alpha_mcse_from_covariance = covariance_alpha_mcse,
    precision_status = result$precision_status, diagnostic_status = alpha_diagnostic,
    requested_draws = result$requested_draws, draws = result$draws, seed = result$seed,
    chunk_size = result$chunk_size, jensen_lower_bound = result$diagnostics$jensen_lower_bound,
    below_jensen = result$diagnostics$below_jensen, rng_used = result$diagnostics$rng_used,
    stringsAsFactors = FALSE)
  pairs <- expand.grid(overlap_row = result$overlap, overlap_column = result$overlap)
  cov_record <- data.frame(case_id = id, M = cfg$M, d = cfg$d, pairs,
    covariance = as.vector(covariance), covariance_kind = "covariance_of_beta_estimates",
    method = result$method, precision_status = result$precision_status,
    row_diagnostic = diagnostic[pairs$overlap_row + 1L],
    column_diagnostic = diagnostic[pairs$overlap_column + 1L],
    draws = result$draws, seed = result$seed, stringsAsFactors = FALSE)
  append_csv(beta, "beta.csv")
  append_csv(alpha, "alpha.csv")
  append_csv(cov_record, "covariance.csv")
  cat("Completed", id, "using", result$method, "in", format(elapsed, digits = 5), "seconds;",
      result$precision_status, "\n")
}
saveRDS(objects, file.path(output, "geometry.rds"))
metadata$completed <- Sys.time()
metadata$elapsed_seconds <- proc.time()[["elapsed"]] - elapsed_start
metadata$case_results <- do.call(rbind, run_records)
metadata$status <- if (any(metadata$case_results$status != "completed")) "completed_with_errors" else "completed"
metadata$code_md5_at_completion <- tools::md5sum(code_paths)
metadata$config_md5_at_completion <- tools::md5sum(config_path)
metadata$input_files_unchanged <- identical(metadata$code_md5, metadata$code_md5_at_completion) &&
  identical(metadata$config_file_md5, metadata$config_md5_at_completion)
metadata$session_at_completion <- utils::sessionInfo()
saveRDS(metadata, file.path(output, "metadata.rds"))
if (!metadata$input_files_unchanged) stop("Source or configuration changed during integration; inspect retained provenance.", call. = FALSE)
if (metadata$status != "completed") stop("One or more integrations failed; retain all partial results and inspect runs.csv.", call. = FALSE)
cat("Universal-constant integration complete. No estimator simulation or precision-goal acceptance was performed.\n")
