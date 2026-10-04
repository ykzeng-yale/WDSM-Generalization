# Archived saved-array calculation source; execution status and scope are documented in results/README.md.
# Definitions: code-release/simulations/paired_comparison_metrics.R and
# simulations/original_wdsm/original_survey_reporting.R. No fit or RNG.
# CLI: companion_records.csv NEW_output_directory [original_report_metrics.csv]
ows_read_saved_summary_csv <- function(path, type = c("companion", "original_metrics")) {
  type <- match.arg(type)
  character_fields <- if (type == "companion") c("original_run_id", "overlap", "sampling_design",
    "estimand", "scenario", "method", "point_source", "point_status", "interval_status", "status", "error",
    "saved_refit_point_status", "saved_refit_interval_status", "saved_refit_error", "saved_refit_variance_divisor",
    "pairing_id", "input_data_sha256", "population_sha256", "input_artifact_sha256", "input_receipt_sha256",
    "sample_completion_sha256", "case_completion_sha256", "case_artifact_sha256", "case_receipt_sha256",
    "source_plan_sha256", "config_sha256", "variance_method", "inference_scope") else
    c("sampling_design", "overlap", "estimand", "scenario", "method", "coverage_rule", "point_metric_scope",
      "relative_bias_status", "reference_scenario", "reference_method", "relative_efficiency_status")
  # Explicit text classes keep quoted empty error strings as "", including
  # an all-success file. Numeric/logical columns retain standard CSV parsing.
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
    colClasses = stats::setNames(rep("character", length(character_fields)), character_fields), na.strings = "NA")
}

ows_summarize_saved_augmentation <- function(records, original_metrics = NULL) {
  required <- c("original_run_id", "overlap", "sampling_design", "estimand", "scenario", "M", "method",
    "point_source", "replicate", "target", "estimate", "variance", "lower", "upper", "status", "point_status",
    "interval_status", "error", "covered_requested", "aligned_with_recorded_point", "saved_refit_estimate",
    "saved_refit_variance", "saved_refit_lower", "saved_refit_upper", "saved_refit_point_status",
    "saved_refit_interval_status", "saved_refit_error", "saved_refit_covered_requested", "saved_refit_B",
    "saved_refit_variance_divisor", "pairing_id", "input_data_sha256", "population_sha256", "config_sha256")
  if (!is.data.frame(records) || !all(required %in% names(records)) || anyDuplicated(names(records)))
    stop("Invalid companion record schema.")
  keys <- c("sampling_design", "overlap", "estimand", "method")
  key <- function(x, fields) do.call(paste, c(unname(x[fields]), sep = "|"))
  plan <- expand.grid(sampling_design = c("retrospective", "prospective"),
    overlap = c("GoodOverlap", "PoorOverlap"), estimand = c("PATE", "PATT"),
    method = c("WM_PS_analytic_M3", "WM_X6_analytic_M3"), stringsAsFactors = FALSE)
  stopifnot(nrow(records) == 16000L, !anyNA(records[c(keys, "scenario", "M", "replicate")]),
    all(records$original_run_id == "original-survey-formal-20260930-v2"),
    all(records$scenario == "CorCor"), all(records$M == 3L),
    all(is.finite(records$replicate)), all(records$replicate %in% 1:1000),
    !anyDuplicated(key(records, c(keys, "replicate"))),
    all(key(records, keys) %in% key(plan, keys)),
    all(records$point_source == sub("_analytic", "", records$method, fixed = TRUE)))
  for (field in c("pairing_id", "input_data_sha256", "population_sha256", "config_sha256")) {
    if (!is.character(records[[field]]) || anyNA(records[[field]]) ||
        any(!grepl("^[0-9a-f]{64}$", records[[field]])))
      stop("Complete authenticated input IDs are required: ", field)
  }
  for (field in c("covered_requested", "saved_refit_covered_requested", "aligned_with_recorded_point"))
    stopifnot(is.logical(records[[field]]), !anyNA(records[[field]]))
  for (field in c("status", "point_status", "interval_status", "saved_refit_point_status", "saved_refit_interval_status"))
    stopifnot(is.character(records[[field]]), !anyNA(records[[field]]))
  stopifnot(all(records$status %in% c("ok", "failed")),
    all(records$point_status %in% c("ok", "failed", "unavailable")),
    all(records$saved_refit_point_status %in% c("ok", "failed", "unavailable")),
    all(records$interval_status %in% c("ok", "failed", "unavailable")),
    all(records$saved_refit_interval_status %in% c("ok", "failed", "unavailable")),
    all(is.finite(records$target)), !anyNA(records$error), !any(records$error == "Not evaluated"))
  # The same original sample must retain its ID across both estimands/families.
  sample_groups <- split(seq_len(nrow(records)), key(records, c("sampling_design", "overlap", "replicate")))
  stopifnot(length(sample_groups) == 4000L)
  for (ix in sample_groups) {
    stopifnot(length(ix) == 4L, length(unique(records$pairing_id[ix])) == 1L,
      length(unique(records$input_data_sha256[ix])) == 1L, length(unique(records$population_sha256[ix])) == 1L)
  }
  one_per_sample <- vapply(sample_groups, `[`, integer(1), 1L)
  stopifnot(!anyDuplicated(records$pairing_id[one_per_sample]))
  cell_groups <- split(seq_len(nrow(records)), factor(key(records, keys), levels = key(plan, keys)))
  for (ix in cell_groups) stopifnot(length(ix) == 1000L, identical(sort(as.integer(records$replicate[ix])), 1:1000),
    length(unique(records$target[ix])) == 1L)
  group_targets <- split(records$target, key(records, c("sampling_design", "overlap", "estimand")))
  stopifnot(all(vapply(group_targets, function(x) length(unique(x)) == 1L, logical(1))))
  status_point <- records$point_status == "ok"
  status_refit_point <- records$saved_refit_point_status == "ok"
  stopifnot(all(is.finite(records$estimate[status_point])), all(is.finite(records$saved_refit_estimate[status_refit_point])),
    identical(status_point, status_refit_point),
    identical(records$estimate[status_point], records$saved_refit_estimate[status_refit_point]))
  interval_ok <- function(status, point, variance, lower, upper) {
    value <- status == "ok"
    stopifnot(all(point[value]), all(is.finite(variance[value])), all(variance[value] >= 0),
      all(is.finite(lower[value])), all(is.finite(upper[value])), all(lower[value] <= upper[value]))
    value
  }
  a <- interval_ok(records$interval_status, status_point, records$variance, records$lower, records$upper)
  b <- interval_ok(records$saved_refit_interval_status, status_refit_point, records$saved_refit_variance,
    records$saved_refit_lower, records$saved_refit_upper)
  stopifnot(all(records$variance[a] > 0), all(records$aligned_with_recorded_point[a]),
    identical(records$status == "ok", status_point & a), all(records$saved_refit_B[b] == 200L),
    all(records$saved_refit_variance_divisor[b] == "B"))
  hit_a <- hit_b <- rep(FALSE, nrow(records))
  hit_a[a] <- records$lower[a] <= records$target[a] & records$target[a] <= records$upper[a]
  hit_b[b] <- records$saved_refit_lower[b] <= records$target[b] & records$target[b] <= records$saved_refit_upper[b]
  stopifnot(identical(hit_a, records$covered_requested), identical(hit_b, records$saved_refit_covered_requested))
  pair_variance <- a & b & records$saved_refit_variance > 0
  pair_variance[is.na(pair_variance)] <- FALSE
  # A row's finite number never overrides its failed source status.
  ratio <- rep(NA_real_, nrow(records))
  ratio[pair_variance] <- records$variance[pair_variance] / records$saved_refit_variance[pair_variance]
  stopifnot(all(is.finite(ratio[pair_variance])))
  paired_records <- records[c(keys, "replicate", "point_source", "pairing_id", "input_data_sha256",
    "target", "estimate", "variance", "saved_refit_variance", "point_status", "interval_status",
    "saved_refit_point_status", "saved_refit_interval_status", "error", "saved_refit_error")]
  paired_records$analytic_available <- a; paired_records$refit_available <- b
  paired_records$analytic_hit <- hit_a; paired_records$refit_hit <- hit_b
  paired_records$coverage_difference <- as.integer(hit_a) - as.integer(hit_b)
  paired_records$common_positive_variance_pair <- pair_variance
  paired_records$analytic_refit_variance_ratio <- ratio
  finite <- function(x) if (length(x) == 1L && is.finite(x)) unname(x) else NA_real_
  mean_or_na <- function(x) if (length(x)) finite(mean(x)) else NA_real_
  mcse <- function(x) if (length(x) > 1L) finite(stats::sd(x) / sqrt(length(x))) else NA_real_
  rows <- comparisons <- diagnostics <- vector("list", nrow(plan))
  for (j in seq_len(nrow(plan))) {
    ix <- cell_groups[[j]]; ix <- ix[order(records$replicate[ix])]
    x <- records[ix, , drop = FALSE]; R <- 1000L; truth <- x$target[1L]
    point <- status_point[ix]; complete_points <- all(point)
    estimate <- x$estimate
    ev <- if (complete_points) finite(stats::var(estimate)) else NA_real_
    bias <- if (complete_points) finite(mean(estimate - truth)) else NA_real_
    q <- if (complete_points) R / (R - 1) * (estimate - mean(estimate))^2 else rep(NA_real_, R)
    common <- pair_variance[ix]; all_common <- all(a[ix] & b[ix])
    modes <- lapply(c("analytic", "original_refit"), function(mode) {
      analytic <- mode == "analytic"
      avail <- if (analytic) a[ix] else b[ix]; hit <- if (analytic) hit_a[ix] else hit_b[ix]
      variance <- if (analytic) x$variance else x$saved_refit_variance
      lower <- if (analytic) x$lower else x$saved_refit_lower
      upper <- if (analytic) x$upper else x$saved_refit_upper
      p <- mean(hit); av <- mean(avail)
      calibration <- if (complete_points && all(avail) && is.finite(ev) && ev > 0) finite(mean(variance) / ev) else NA_real_
      data.frame(plan[j, ], inference = mode, requested = R, received = R, target = truth,
        point_available = sum(point), point_failed = sum(!point), interval_available = sum(avail), interval_failed = sum(!avail),
        point_metric_scope = if (complete_points) "all_requested" else "withheld_incomplete_points",
        relative_bias_status = if (!complete_points) "withheld_incomplete_points" else if (truth == 0) "undefined_zero_target" else "defined",
        mean_estimate = if (complete_points) finite(mean(estimate)) else NA_real_, signed_bias = bias,
        signed_relative_bias_percent = if (truth != 0) finite(100 * bias / abs(truth)) else NA_real_,
        signed_relative_bias_mcse = if (truth != 0) finite(100 * sqrt(ev / R) / abs(truth)) else NA_real_,
        empirical_variance = ev, empirical_variance_divisor = "R_minus_1",
        empirical_variance_mcse = if (complete_points) mcse(q) else NA_real_,
        operational_coverage = p, coverage_denominator = R, operational_coverage_mcse = sqrt(p * (1 - p) / R),
        interval_availability = av, availability_mcse = sqrt(av * (1 - av) / R),
        mean_reported_variance_all_requested = if (all(avail)) finite(mean(variance)) else NA_real_,
        variance_calibration_ratio = calibration,
        variance_calibration_mcse_delta = if (is.finite(calibration)) mcse((variance - calibration * q) / ev) else NA_real_,
        calibration_status = if (!complete_points) "withheld_incomplete_points" else if (!all(avail))
          "withheld_incomplete_intervals" else if (!is.finite(calibration)) "undefined_zero_or_nonfinite_arithmetic" else "all_requested",
        stringsAsFactors = FALSE)
    })
    rows[[j]] <- do.call(rbind, modes)
    delta <- as.integer(hit_a[ix]) - as.integer(hit_b[ix])
    variance_ratio <- if (all_common && mean(x$saved_refit_variance) > 0)
      finite(mean(x$variance) / mean(x$saved_refit_variance)) else NA_real_
    if (all_common && mean(x$saved_refit_variance) > 0 &&
        (!is.finite(variance_ratio) || variance_ratio <= 0)) stop("Nonrepresentable complete paired variance ratio.")
    comparisons[[j]] <- data.frame(plan[j, ], requested_pairs = R, common_original_ids = R,
      common_point_pairs = sum(point), common_available_intervals = sum(a[ix] & b[ix]),
      analytic_only_available = sum(a[ix] & !b[ix]), refit_only_available = sum(!a[ix] & b[ix]),
      neither_interval_available = sum(!a[ix] & !b[ix]), common_positive_variance_pairs = sum(common),
      unavailable_variance_pairs = sum(!(a[ix] & b[ix])), nonpositive_or_unavailable_individual_ratio_pairs = sum(!common),
      both_cover = sum(hit_a[ix] & hit_b[ix]),
      analytic_only_covers = sum(hit_a[ix] & !hit_b[ix]), refit_only_covers = sum(!hit_a[ix] & hit_b[ix]),
      neither_covers = sum(!hit_a[ix] & !hit_b[ix]),
      paired_operational_coverage_difference = mean(delta), paired_operational_coverage_difference_mcse = mcse(delta),
      analytic_refit_ratio_of_mean_variances = variance_ratio,
      ratio_of_mean_variances_mcse_delta = if (is.finite(variance_ratio))
        mcse((x$variance - variance_ratio * x$saved_refit_variance) / mean(x$saved_refit_variance)) else NA_real_,
      mean_paired_variance_difference = if (all_common) finite(mean(x$variance - x$saved_refit_variance)) else NA_real_,
      paired_variance_difference_mcse = if (all_common) mcse(x$variance - x$saved_refit_variance) else NA_real_,
      variance_comparison_scope = if (!all_common) "withheld_incomplete_variance_pairs" else
        if (!is.finite(variance_ratio)) "undefined_zero_or_nonfinite_mean_ratio" else "all_requested",
      stringsAsFactors = FALSE)
    # Optional diagnostic values cannot replace the all-requested primary rows.
    diagnostics[[j]] <- data.frame(plan[j, ], scope = "conditional_on_common_available_positive_variance_pairs",
      requested = R, denominator = sum(common), unavailable_pairs = sum(!common),
      analytic_mean_variance = mean_or_na(x$variance[common]), refit_mean_variance = mean_or_na(x$saved_refit_variance[common]),
      ratio_of_means = if (any(common)) finite(mean(x$variance[common]) / mean(x$saved_refit_variance[common])) else NA_real_,
      mean_individual_variance_ratio = mean_or_na(ratio[ix][common]),
      mean_individual_ratio_mcse = mcse(ratio[ix][common]), stringsAsFactors = FALSE)
  }
  metrics <- do.call(rbind, rows); paired <- do.call(rbind, comparisons); diagnostic <- do.call(rbind, diagnostics)
  complete <- metrics$point_metric_scope == "all_requested"
  point_fields <- c("mean_estimate", "signed_bias", "empirical_variance", "empirical_variance_mcse")
  stopifnot(all(is.finite(as.matrix(metrics[complete, point_fields, drop = FALSE]))),
    all(is.finite(metrics$signed_relative_bias_percent[complete & metrics$target != 0])),
    all(is.finite(metrics$signed_relative_bias_mcse[complete & metrics$target != 0])))
  complete_calibration <- complete & metrics$interval_available == 1000L & metrics$empirical_variance > 0
  complete_calibration[is.na(complete_calibration)] <- FALSE
  stopifnot(all(is.finite(metrics$variance_calibration_ratio[complete_calibration])),
    all(is.finite(metrics$variance_calibration_mcse_delta[complete_calibration])))
  complete_pairs <- paired$common_available_intervals == 1000L & paired$analytic_refit_ratio_of_mean_variances > 0
  complete_pairs[is.na(complete_pairs)] <- FALSE
  pair_fields <- c("analytic_refit_ratio_of_mean_variances", "ratio_of_mean_variances_mcse_delta",
    "mean_paired_variance_difference", "paired_variance_difference_mcse")
  stopifnot(all(is.finite(as.matrix(paired[complete_pairs, pair_fields, drop = FALSE]))))
  stopifnot(all(is.finite(paired$mean_paired_variance_difference[paired$common_available_intervals == 1000L])),
    all(is.finite(paired$paired_variance_difference_mcse[paired$common_available_intervals == 1000L])))
  metrics$reference_method <- "PSM_M1"
  metrics$reference_scenario <- "CorCor"
  metrics$reference_empirical_variance <- metrics$relative_efficiency <- metrics$relative_efficiency_mcse <- NA_real_
  metrics$relative_efficiency_status <- "not_computed_reference_absent_from_companion_records"
  if (!is.null(original_metrics)) {
    required_reference <- c(keys, "scenario", "requested", "truth", "point_metric_scope", "mean_estimate",
      "empirical_variance", "reference_method", "reference_scenario", "reference_empirical_variance",
      "relative_efficiency", "relative_efficiency_mcse", "relative_efficiency_status")
    stopifnot(is.data.frame(original_metrics), all(required_reference %in% names(original_metrics)))
    close <- function(a, b) (is.na(a) && is.na(b)) || (is.finite(a) && is.finite(b) && abs(a - b) <= 1e-10 * max(1, abs(a), abs(b)))
    for (j in seq_len(nrow(metrics))) {
      m <- metrics[j, ]
      match_group <- original_metrics$sampling_design == m$sampling_design & original_metrics$overlap == m$overlap &
        original_metrics$estimand == m$estimand & original_metrics$scenario == "CorCor"
      source <- original_metrics[match_group & original_metrics$method == sub("_analytic", "", m$method, fixed = TRUE), ]
      baseline <- original_metrics[match_group & original_metrics$method == "PSM_M1", ]
      stopifnot(nrow(source) == 1L, nrow(baseline) == 1L, source$requested == 1000L, baseline$requested == 1000L,
        identical(source$reference_method, "PSM_M1"), identical(source$reference_scenario, "CorCor"),
        identical(source$point_metric_scope, m$point_metric_scope), close(source$truth, m$target),
        close(source$mean_estimate, m$mean_estimate), close(source$empirical_variance, m$empirical_variance))
      if (is.finite(source$reference_empirical_variance)) stopifnot(close(source$reference_empirical_variance, baseline$empirical_variance))
      if (identical(source$relative_efficiency_status, "defined")) stopifnot(
        close(source$relative_efficiency, baseline$empirical_variance / source$empirical_variance))
      metrics$reference_empirical_variance[j] <- source$reference_empirical_variance
      metrics$relative_efficiency[j] <- source$relative_efficiency
      metrics$relative_efficiency_mcse[j] <- source$relative_efficiency_mcse
      metrics$relative_efficiency_status[j] <- paste0("original_point_report:", source$relative_efficiency_status)
    }
  }
  list(metrics = metrics, paired_comparisons = paired, paired_records = paired_records,
    conditional_variance_diagnostics = diagnostic,
    scope = paste("Original finite-population/dependent-design empirical comparison; all 1000 requested samples per cell retained.",
      "Statistical failures contribute unavailable intervals, not deleted rows. Primary point/calibration metrics are withheld when incomplete.",
      "No cluster/stratum inference theorem, nominal-coverage guarantee or analytic/refit equality is asserted.",
      "Relative efficiency, when supplied, is copied from the original PSM_M1-baseline point report after binding checks."))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  stopifnot(length(args) %in% c(2L, 3L), requireNamespace("digest", quietly = TRUE),
    requireNamespace("jsonlite", quietly = TRUE), !exists(".Random.seed", .GlobalEnv, inherits = FALSE))
  input_path <- normalizePath(args[1L], mustWork = TRUE)
  out <- args[2L]
  stopifnot(!file.exists(out), dir.exists(dirname(out)))
  records <- ows_read_saved_summary_csv(input_path, "companion")
  reference_path <- if (length(args) == 3L) normalizePath(args[3L], mustWork = TRUE) else NULL
  original_metrics <- if (!is.null(reference_path)) ows_read_saved_summary_csv(reference_path, "original_metrics") else NULL
  result <- ows_summarize_saved_augmentation(records, original_metrics)
  stopifnot(!exists(".Random.seed", .GlobalEnv, inherits = FALSE), dir.create(out, recursive = FALSE, showWarnings = FALSE))
  for (name in c("metrics", "paired_comparisons", "paired_records", "conditional_variance_diagnostics"))
    utils::write.csv(result[[name]], file.path(out, paste0(name, ".csv")), row.names = FALSE, na = "NA")
  script <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])
  jsonlite::write_json(list(status = "COMPLETE", source = input_path,
    input_sha256 = digest::digest(file = input_path, algo = "sha256"),
    summary_script_sha256 = digest::digest(file = script, algo = "sha256"),
    original_point_report = reference_path,
    original_point_report_sha256 = if (!is.null(reference_path)) digest::digest(file = reference_path, algo = "sha256") else NULL,
    original_configuration_hashes = sort(unique(records$config_sha256)),
    requested_rows = 16000L, cells = 16L, requested_per_cell = 1000L,
    source_point_definitions = "paired_comparison_metrics.R:118-165; original_survey_reporting.R:146-147",
    scope = result$scope), file.path(out, "completion.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
  cat("Summarized 16000 saved companion rows; requested denominator1000 in all16 cells.\n")
}
