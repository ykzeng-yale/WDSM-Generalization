#!/usr/bin/env Rscript
# Install with the complete audited aggregate exports in results/.
# Replays benchmark arithmetic from audited aggregates; never generates subjects,
# fits a model, integrates geometry, or draws random numbers.

.wm_fsr_close <- function(x, y, tolerance = 2e-10) {
  isTRUE(all.equal(x, y, tolerance = tolerance, check.attributes = FALSE))
}
.wm_fsr_read <- function(directory, name, fields = NULL) {
  x <- utils::read.csv(file.path(directory, name), stringsAsFactors = FALSE)
  if (!is.null(fields) && !identical(names(x), fields)) stop("Unexpected schema: ", name)
  x
}
.wm_fsr_alpha_fields <- c("key", "M", "d", "source_class", "method", "precision_status",
  "draws", "requested_draws", "seed", "chunk_size", "alpha", "alpha_mcse",
  "integrand_upper_bound", "true_mcse_upper_bound", "hoeffding_95_half_width",
  "minimum_draw_value", "maximum_draw_value", "jensen_lower_bound", "below_jensen",
  "sphere_dimension", "count_terms", "gaussian_coordinate_work", "count_term_work",
  "rng_used", "process_elapsed_seconds", "source_rds_md5", "run_script_md5",
  "geometry_api_md5", "geometry_common_md5", "manifest_md5")
.wm_fsr_selection_fields <- c("key", "M", "d", "used_by_first_stage", "source_kind",
  "source_csv", "source_class", "method", "precision_status", "draws", "seed",
  "alpha", "alpha_mcse", "source_rds_md5", "assembled_rds_md5")

wm_results_first_stage_geometry <- function(results_directory, boundary, used_keys) {
  directory <- file.path(results_directory, "first_stage")
  selection <- .wm_fsr_read(directory, "geometry_selection.csv", .wm_fsr_selection_fields)
  alpha <- .wm_fsr_read(directory, "alpha_only.csv", .wm_fsr_alpha_fields)
  expected <- as.vector(outer(c(1L, 3L), c(1L, 2L, 3L, 5L),
    function(M, d) sprintf("M%d_d%d", M, d)))
  if (nrow(selection) != 8L || anyNA(selection) || anyDuplicated(selection$key) ||
      !setequal(selection$key, expected) || !setequal(names(boundary), expected) ||
      !is.logical(selection$used_by_first_stage) ||
      !setequal(selection$key[selection$used_by_first_stage], used_keys) ||
      any(selection$key != sprintf("M%d_d%d", selection$M, selection$d)) ||
      any(!grepl("^[0-9a-f]{32}$", selection$source_rds_md5)) ||
      any(!grepl("^[0-9a-f]{32}$", selection$assembled_rds_md5)) ||
      length(unique(selection$assembled_rds_md5)) != 1L)
    stop("Incomplete or inconsistent actual first-stage geometry selection")
  if (nrow(alpha) != 1L || anyNA(alpha) || alpha$key != "M3_d5" ||
      alpha$M != 3L || alpha$d != 5L || alpha$source_class != "wm_geometry_alpha" ||
      alpha$method != "reduced_poisson_count_sphere_monte_carlo" ||
      alpha$precision_status != "monte_carlo_estimate" || alpha$draws != 200000L ||
      alpha$requested_draws != alpha$draws || alpha$seed != 935205L ||
      alpha$chunk_size != 4096L || !is.logical(alpha$below_jensen) ||
      !is.logical(alpha$rng_used) || !isTRUE(alpha$rng_used))
    stop("Alpha-only export differs from the predeclared first-stage run")
  numeric_fields <- names(alpha)[vapply(alpha, is.numeric, logical(1L))]
  if (any(!is.finite(unlist(alpha[numeric_fields], use.names = FALSE))) ||
      any(!vapply(alpha[c("source_rds_md5", "run_script_md5", "geometry_api_md5",
        "geometry_common_md5", "manifest_md5")], function(x)
          length(x) == 1L && grepl("^[0-9a-f]{32}$", x), logical(1L))))
    stop("Nonfinite alpha diagnostics or invalid source fingerprints")
  bound <- 2^(alpha$d + 1) * gamma(alpha$d / 2 + 1)^2 /
    (alpha$d * gamma(alpha$d)) * alpha$M * (2 * alpha$M - 1)
  terms <- sum(vapply(0:(alpha$M - 1L), function(l) (alpha$M - l)^2, numeric(1L)))
  if (!.wm_fsr_close(alpha$integrand_upper_bound, bound) ||
      !.wm_fsr_close(alpha$true_mcse_upper_bound, bound / (2 * sqrt(alpha$draws))) ||
      !.wm_fsr_close(alpha$hoeffding_95_half_width, bound * sqrt(log(40) / (2 * alpha$draws))) ||
      alpha$minimum_draw_value <= 0 || alpha$minimum_draw_value > alpha$alpha ||
      alpha$alpha > alpha$maximum_draw_value || alpha$maximum_draw_value > bound * (1 + 1e-10) ||
      alpha$alpha_mcse <= 0 ||
      alpha$alpha_mcse > bound / (2 * sqrt(alpha$draws - 1)) * (1 + 1e-10) ||
      alpha$jensen_lower_bound != alpha$M^2 ||
      alpha$below_jensen != (alpha$alpha < alpha$M^2) ||
      alpha$sphere_dimension != 2 * alpha$d || alpha$count_terms != terms ||
      alpha$gaussian_coordinate_work != 2 * alpha$d * alpha$draws ||
      alpha$count_term_work != terms * alpha$draws || alpha$process_elapsed_seconds < 0)
    stop("Inconsistent bounded-integrand alpha diagnostics")
  diagnostics <- as.list(alpha[c("integrand_upper_bound", "true_mcse_upper_bound",
    "hoeffding_95_half_width", "minimum_draw_value", "maximum_draw_value",
    "jensen_lower_bound", "below_jensen", "sphere_dimension", "count_terms",
    "gaussian_coordinate_work", "count_term_work", "rng_used")])
  direct <- structure(c(as.list(alpha[c("M", "d", "alpha", "alpha_mcse", "method",
    "precision_status", "draws", "requested_draws", "seed", "chunk_size")]),
    list(diagnostics = diagnostics)), class = c("wm_geometry_alpha", "list"))
  result <- boundary
  for (i in seq_len(nrow(selection))) {
    r <- selection[i, , drop = FALSE]
    is_direct <- r$key == "M3_d5"
    if (r$source_kind != (if (is_direct) "alpha_only" else "boundary") ||
        r$source_csv != (if (is_direct) "first_stage/alpha_only.csv" else "geometry/geometry.csv") ||
        r$source_class != (if (is_direct) "wm_geometry_alpha" else "wm_geometry"))
      stop("Geometry selection does not preserve the declared separate sources")
    g <- if (is_direct) direct else boundary[[r$key]]
    for (f in c("M", "d", "method", "precision_status", "draws", "seed", "alpha", "alpha_mcse"))
      if (!.wm_fsr_close(r[[f]], g[[f]])) stop("Geometry index differs: ", r$key, "/", f)
    if (is_direct && r$source_rds_md5 != alpha$source_rds_md5)
      stop("Alpha source fingerprint differs from the selection manifest")
    result[[r$key]] <- g
  }
  if (length(unique(selection$source_rds_md5[selection$source_kind == "boundary"])) != 1L ||
      !is.null(result$M3_d5$beta) || !is.null(result$M3_d5$covariance))
    stop("Boundary sources differ or alpha-only object contains invented beta moments")
  result
}

wm_results_reproduce_first_stage <- function(package_directory, output_directory = NULL) {
  root <- normalizePath(package_directory, mustWork = TRUE)
  saved <- file.path(root, "results")
  h <- new.env(parent = globalenv())
  sys.source(file.path(saved, "reproduce_benchmarks.R"), envir = h)
  for (file in c("first_stage.R", "first_stage_benchmarks.R", "first_stage_join.R",
                 "extension_production.R"))
    sys.source(file.path(root, "simulations", file), envir = h)
  hashes <- .wm_fsr_read(saved, "artifact_hashes.csv", c("path", "bytes", "sha256", "md5"))
  needed <- c("first_stage/summary.csv", "first_stage/summary_joined.csv",
    "first_stage/config.csv", "first_stage/task_plan.csv", "first_stage/alpha_only.csv",
    "first_stage/geometry_selection.csv", "first_stage/recovery_summary.csv",
    "first_stage/provenance.json", "geometry/geometry.csv", "geometry/beta.csv",
    "geometry/beta_covariance.csv", "reproduce_benchmarks.R", "reproduce_first_stage.R")
  if (anyNA(hashes) || anyDuplicated(hashes$path) || !all(needed %in% hashes$path) ||
      any(grepl("(^/|(^|/)\\.\\.(/|$))", hashes$path)) ||
      any(!grepl("^[0-9a-f]{32}$", hashes$md5)) ||
      any(!grepl("^[0-9a-f]{64}$", hashes$sha256))) stop("Invalid public artifact manifest")
  paths <- file.path(saved, hashes$path)
  if (any(!file.exists(paths)) || !identical(unname(tools::md5sum(paths)), hashes$md5) ||
      any(file.info(paths)$size != hashes$bytes)) stop("Public artifact hash/size mismatch")
  directory <- file.path(saved, "first_stage")
  config <- h$wm_first_stage_config(.wm_fsr_read(directory, "config.csv"))
  expected_plan <- h$wm_extension_plan("first_stage", config)
  plan <- .wm_fsr_read(directory, "task_plan.csv")
  if (!identical(names(plan), names(expected_plan)) || !.wm_fsr_close(plan, expected_plan, 0))
    stop("Saved task plan differs from the complete prescribed 446-task grid")
  summary <- .wm_fsr_read(directory, "summary.csv")
  expected_grid <- expand.grid(scenario = config$scenario,
    estimand = c("PATE", "PATT"),
    comparison = c("known_first_stage", "fitted_naive", "fitted_adjusted"),
    stringsAsFactors = FALSE)
  group_key <- function(x) paste(x$scenario, x$estimand, x$comparison, sep = "/")
  if (nrow(summary) != 234L || anyDuplicated(group_key(summary)) ||
      !setequal(group_key(summary), group_key(expected_grid)) ||
      !all(c("requested", "completed", "failed", "independent_replications",
             "records_per_replication") %in% names(summary)))
    stop("Saved first-stage summary is not the complete 39-cell grid")
  if (anyNA(summary[c("requested", "completed", "failed")]) ||
      any(summary$completed + summary$failed != summary$requested) ||
      any(summary$independent_replications != 1000L) ||
      any(summary$records_per_replication != 6L))
    stop("Saved first-stage summary lost planned replication/failure accounting")
  recovery <- .wm_fsr_read(directory, "recovery_summary.csv",
    c("scenario", "estimand", "comparison", "action", "selected_status", "count"))
  if (!nrow(recovery) || nrow(recovery) > 234L * 6L || anyNA(recovery) ||
      anyDuplicated(recovery[setdiff(names(recovery), "count")]) ||
      any(!group_key(recovery) %in% group_key(expected_grid)) ||
      any(!recovery$action %in% c("retained_original", "guard_replay", "missing_dataset")) ||
      any(!recovery$selected_status %in% c("ok", "error")) ||
      !is.numeric(recovery$count) || any(!is.finite(recovery$count)) ||
      any(recovery$count <= 0 | recovery$count != floor(recovery$count)))
    stop("Invalid aggregate-only recovery accounting")
  for (i in seq_len(nrow(summary))) {
    r <- recovery[group_key(recovery) == group_key(summary[i, , drop = FALSE]), , drop = FALSE]
    if (sum(r$count[r$selected_status == "ok"]) != summary$completed[i] ||
        sum(r$count[r$selected_status == "error"]) != summary$failed[i])
      stop("Recovery counts do not reconcile with retained aggregate successes/failures")
  }
  index <- match(summary$scenario, config$scenario)
  for (f in c("branch", "n", "m", "d", "M"))
    if (!.wm_fsr_close(summary[[f]], config[[f]][index], 0)) stop("Summary/config mismatch: ", f)
  if (!.wm_fsr_close(summary$requested, config$replications[index], 0))
    stop("Summary does not retain every planned replication")
  used <- unique(sprintf("M%d_d%d", config$M, config$d))
  geometry <- wm_results_first_stage_geometry(saved,
    h$wm_results_geometry(file.path(saved, "geometry")), used)
  joined <- h$wm_first_stage_join(summary, geometry)
  reference <- .wm_fsr_read(directory, "summary_joined.csv")
  if (nrow(reference) != 234L || !identical(names(joined), names(reference)))
    stop("Saved first-stage joined schema differs")
  for (f in names(reference)) if (!.wm_fsr_close(joined[[f]], reference[[f]]))
    stop("First-stage benchmark reproduction differs: ", f)
  # All rows using one key share the same alpha estimate, including across n.
  keys <- sprintf("M%d_d%d", joined$M, joined$d)
  covariance <- outer(joined$alpha_coefficient * joined$alpha_mcse,
                      joined$alpha_coefficient * joined$alpha_mcse)
  covariance[outer(keys, keys, "!=")] <- 0
  dimnames(covariance) <- rep(list(paste(joined$scenario, joined$estimand,
    joined$comparison, sep = "/")), 2L)
  if (!.wm_fsr_close(sqrt(diag(covariance)), joined$benchmark_geometry_mcse))
    stop("Shared geometry covariance disagrees with row MCSEs")
  if (!is.null(output_directory)) {
    if (file.exists(output_directory)) stop("Use a new output directory")
    dir.create(output_directory, recursive = TRUE)
    utils::write.csv(joined, file.path(output_directory, "first_stage_joined.csv"), row.names = FALSE)
    saveRDS(list(geometry = geometry, benchmark_geometry_covariance = covariance),
      file.path(output_directory, "first_stage_geometry.rds"))
  }
  cat("PASS: artifact checks and all 234 saved first-stage benchmark joins agree.\n")
  cat("This replays benchmark arithmetic; it does not independently reconstruct sampling summaries.\n")
  invisible(list(joined = joined, geometry = geometry,
    benchmark_geometry_covariance = covariance))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 2L) stop("usage: Rscript reproduce_first_stage.R [package_directory [new_output_directory]]")
  script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(script) != 1L) stop("Run this interface with Rscript")
  root <- if (length(args)) args[1L] else file.path(dirname(normalizePath(script)), "..")
  wm_results_reproduce_first_stage(root, if (length(args) == 2L) args[2L] else NULL)
}
