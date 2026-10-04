# Archived summary source. Requires caller-supplied saved inputs and contract.
# 54000-row comparison. No estimator, fitted influence, graph or bootstrap call.
args <- commandArgs(trailingOnly = TRUE); stopifnot(length(args) == 3L)
root <- normalizePath(args[1L], mustWork = TRUE)
out <- normalizePath(args[2L], mustWork = TRUE)
cfg <- jsonlite::read_json(normalizePath(args[3L], mustWork = TRUE), simplifyVector = FALSE)
sha <- function(p) digest::digest(file = p, algo = "sha256", serialize = FALSE)
rng <- function() list(kind = RNGkind(), present = exists(".Random.seed", .GlobalEnv, inherits = FALSE),
  seed = get0(".Random.seed", .GlobalEnv, inherits = FALSE))
initial_rng <- rng()
binding <- jsonlite::read_json(file.path(out, "binding.json"), simplifyVector = FALSE)
stopifnot(identical(binding$status, "COMPLETE_SOURCE_POINTS_BOUND"), binding$rows == 54000L,
  binding$cases == 3000L, binding$source_metric_rows == 72L, binding$source_audit_rows == 72L)
for (f in names(binding$output_pins)) stopifnot(identical(sha(file.path(out, f)), binding$output_pins[[f]]))
metric_source <- file.path(root, cfg$paired_metrics_source)
stopifnot(identical(sha(metric_source), cfg$source_pins[[cfg$paired_metrics_source]]))
metric_env <- new.env(parent = globalenv()); sys.source(metric_source, metric_env)
read_csv <- function(p, text) {
  header <- names(read.csv(p, nrows = 0L, check.names = FALSE, colClasses = "character"))
  classes <- setNames(rep(NA_character_, length(header)), header)
  classes[intersect(text, header)] <- "character"
  read.csv(p, check.names = FALSE, stringsAsFactors = FALSE, colClasses = classes, na.strings = "NA")
}
text <- c("response", "estimand", "family", "source_method", "point_status", "analytic_status", "error",
  "saved_refit_status", "saved_refit_error", "saved_refit_divisor", "scope", "pairing_id")
x <- read_csv(file.path(out, "bound_records.csv"), text)
# The pilot predates compact warning receipts. Read its already accepted result
# once for metadata only; no nuisance/inference function is called again.
pilot_path <- file.path(root, cfg$pilot_result)
stopifnot(identical(sha(pilot_path), cfg$source_pins[[cfg$pilot_result]]))
pilot <- readRDS(pilot_path)
pilot_ix <- which(x$population_scenario == 1L & x$replicate == 1L)
stopifnot(length(pilot_ix) == 18L, all(is.na(x$analytic_warning_count[pilot_ix])))
for (i in pilot_ix) {
  family <- x$family[i]; q <- x$multiplier[i]
  assembly <- pilot$stacks[[family]]
  inference <- pilot$assemblies[[paste0("WM_", family, "_analytic_M3_q", q)]]
  stopifnot(is.list(assembly), is.list(inference), isTRUE(assembly$ok), isTRUE(inference$ok))
  captured <- "warnings" %in% names(assembly) && "warnings" %in% names(inference) &&
    is.character(assembly$warnings) && is.character(inference$warnings)
  x$analytic_warning_count[i] <- if (captured) length(assembly$warnings) + length(inference$warnings) else NA_real_
}
rm(pilot); gc(verbose = FALSE)
original <- read_csv(file.path(out, "original_metrics72.csv"), c("response", "method", "estimand", "scenario",
  "coverage_rule", "point_metric_scope", "relative_bias_status", "reference_scenario", "reference_method", "relative_efficiency_status"))
original_audit <- read_csv(file.path(out, "original_audit72.csv"), c("response", "method", "estimand", "scenario", "coverage_rule"))
groups <- c("population_scenario", "multiplier", "response", "estimand", "family")
key <- function(d, cols) do.call(paste, c(unname(d[cols]), list(sep = "|")))
stopifnot(nrow(x) == 54000L, all(x$point_status == "ok"), all(is.finite(x$estimate)),
  all(is.finite(x$target)), all(x$saved_refit_B == 200L), all(x$saved_refit_divisor == "B"),
  !anyNA(x$error), !anyNA(x$saved_refit_error), !anyDuplicated(key(x, c(groups, "replicate"))))
plan <- expand.grid(population_scenario = 1:3, multiplier = 1:6, family = c("PS", "DSM", "X6"),
  stringsAsFactors = FALSE)
plan$response <- "No"; plan$estimand <- "PATT"
stopifnot(setequal(key(x, groups), key(plan, groups)))
make_inference <- function(mode) {
  analytic <- mode == "analytic"
  y <- x[c(groups, "replicate", "target", "estimate", "point_status", "pairing_id")]
  y$scenario <- "same_saved_point"; y$method <- mode
  y$variance <- if (analytic) x$analytic_variance else x$saved_refit_variance
  y$lower <- if (analytic) x$analytic_lower else x$saved_refit_lower
  y$upper <- if (analytic) x$analytic_upper else x$saved_refit_upper
  component <- if (analytic) x$analytic_status else x$saved_refit_status
  stopifnot(all(component %in% c("ok", "failed", "failed_interval", "unavailable")))
  y$interval_status <- ifelse(component == "ok", "ok", "failed_interval")
  y$status <- ifelse(component == "ok", "ok", "failed")
  y$error <- if (analytic) x$error else x$saved_refit_error
  warning_count <- if (analytic) x$analytic_warning_count else x$saved_refit_warning_count
  y$warning_count_known <- is.finite(warning_count)
  # Generic helper requires a nonnegative number. This internal placeholder
  # affects only its warning diagnostic, which is withheld below if unknown;
  # it never changes points, intervals, coverage or variance calculations.
  y$warning_count <- ifelse(y$warning_count_known, warning_count, 0L)
  y
}
records <- rbind(make_inference("analytic"), make_inference("original_refit"))
planned <- rbind(transform(plan, method = "analytic"), transform(plan, method = "original_refit"))
planned$scenario <- "same_saved_point"; planned$requested <- 1000L
planned$interval_required <- TRUE; planned$coverage_rule <- "inclusive"
base <- metric_env$wm_paired_comparison_metrics(records, planned, groups,
  reference_scenario = "same_saved_point", reference_method = "original_refit")
stopifnot(base$records_complete, base$point_metrics_complete, base$reporting_inputs_valid,
  nrow(base$metrics) == 108L, nrow(base$audit) == 108L)

# Copy canonical original point summaries after complete source correspondence.
# No baseline-point input or new RE/paired-RE-MCSE calculation is requested.
point_fields <- c("truth", "point_metric_scope", "relative_bias_status", "mean_estimate", "signed_bias",
  "absolute_bias", "signed_relative_bias_percent", "empirical_variance", "empirical_sd", "bias_mcse",
  "signed_relative_bias_mcse", "mse", "rmse", "rmse_mcse", "variance_mcse", "reference_scenario",
  "reference_method", "reference_empirical_variance", "relative_efficiency", "relative_efficiency_mcse",
  "relative_efficiency_status")
close <- function(a, b) length(a) == length(b) && all((is.na(a) & is.na(b)) |
  (is.finite(a) & is.finite(b) & abs(a - b) <= 1e-10 * pmax(1, abs(a), abs(b))))
mcse <- function(z) stats::sd(z) / sqrt(length(z))
metrics <- paired <- vector("list", nrow(plan))
x_group <- key(x, groups); base_group <- key(base$metrics, groups); plan_group <- key(plan, groups)
records_group <- key(records, groups)
base$audit$warning_known_replicates <- base$audit$warning_unknown_replicates <- integer(nrow(base$audit))
for (j in seq_len(nrow(plan))) {
  ix <- which(x_group == plan_group[j])
  d <- x[ix, , drop = FALSE]; d <- d[order(d$replicate), , drop = FALSE]
  stopifnot(nrow(d) == 1000L, identical(as.integer(d$replicate), 1:1000))
  method <- paste0("WM_", plan$family[j], "_M3")
  common <- original$population_scenario == plan$population_scenario[j] &
    original$multiplier == plan$multiplier[j] & original$response == "No" & original$estimand == "PATT" &
    original$scenario == "source_specification"
  source <- original[common & original$method == method, , drop = FALSE]
  reference <- original[common & original$method == "Lenis_U_PS_U_OM_OW_U", , drop = FALSE]
  stopifnot(nrow(source) == 1L, nrow(reference) == 1L, source$requested == 1000L,
    reference$requested == 1000L, source$point_metric_scope == "all_requested",
    reference$point_metric_scope == "all_requested", source$reference_scenario == "source_specification",
    source$reference_method == "Lenis_U_PS_U_OM_OW_U",
    close(source$truth, d$target[1L]), close(source$mean_estimate, mean(d$estimate)),
    close(source$empirical_variance, stats::var(d$estimate)),
    close(source$reference_empirical_variance, reference$empirical_variance),
    source$relative_efficiency_status == "defined",
    close(source$relative_efficiency, reference$empirical_variance / source$empirical_variance),
    is.finite(source$relative_efficiency_mcse), source$relative_efficiency_mcse >= 0)
  ev <- source$empirical_variance; q <- 1000 / 999 * (d$estimate - mean(d$estimate))^2
  stopifnot(close(mean(q), ev), is.finite(ev), ev > 0)
  avail <- hit <- list(); rows <- list()
  for (mode in c("analytic", "original_refit")) {
    i <- which(base_group == plan_group[j] & base$metrics$method == mode)
    stopifnot(length(i) == 1L)
    m <- base$metrics[i, , drop = FALSE]; a <- base$audit[i, , drop = FALSE]
    if (mode == "original_refit") {
      # Same authenticated saved procedure: these summaries must reproduce its
      # canonical report. The fixed1e-10 scaled tolerance in close() covers only
      # CSV serialization/arithmetic; it is not chosen from observed calibration.
      for (field in c("operational_coverage", "interval_availability", "mean_reported_variance_valid_intervals"))
        if (!close(m[[field]], source[[field]]))
          stop("Original refit summary differs from authenticated canonical report: ", field, call. = FALSE)
    }
    # Self-reference RE from the generic summary is never presented as Lenis RE.
    keep <- c(groups, "requested", "coverage_rule", "interval_availability", "availability_mcse",
      "operational_coverage", "operational_coverage_mcse", "coverage_exact95_lower", "coverage_exact95_upper",
      "operational_denominator", "conditional_coverage", "conditional_coverage_mcse", "conditional_coverage_denominator",
      "mean_valid_ci_length", "mean_reported_variance_valid_intervals", "warning_rate")
    row <- m[keep]; row$inference <- mode; row$source_method <- method
    warning_known <- records$warning_count_known[records_group == plan_group[j] & records$method == mode]
    row$warning_known_replicates <- sum(warning_known); row$warning_unknown_replicates <- sum(!warning_known)
    base$audit$warning_known_replicates[i] <- sum(warning_known)
    base$audit$warning_unknown_replicates[i] <- sum(!warning_known)
    if (!all(warning_known)) {
      row$warning_rate <- NA_real_
      base$audit$warning_replicates[i] <- NA_integer_
    }
    row[point_fields] <- source[point_fields]
    for (k in c("received", "missing", "evaluated", "not_run", "point_available", "point_unavailable",
      "interval_available", "interval_unavailable", "invalid_points", "invalid_intervals", "invalid_status")) row[[k]] <- a[[k]]
    analytic <- mode == "analytic"
    v <- if (analytic) d$analytic_variance else d$saved_refit_variance
    lo <- if (analytic) d$analytic_lower else d$saved_refit_lower
    hi <- if (analytic) d$analytic_upper else d$saved_refit_upper
    status <- if (analytic) d$analytic_status else d$saved_refit_status
    valid <- status == "ok" & is.finite(v) & v >= 0 & is.finite(lo) & is.finite(hi) & lo <= hi
    stopifnot(sum(valid) == a$interval_available)
    covered <- valid & lo <= d$target & d$target <= hi; covered[!valid] <- FALSE
    avail[[mode]] <- valid; hit[[mode]] <- covered
    complete <- all(valid)
    row$mean_reported_variance_all_requested <- if (complete) mean(v) else NA_real_
    row$mean_reported_variance_mcse <- if (complete) mcse(v) else NA_real_
    row$variance_calibration_ratio <- if (complete) mean(v) / ev else NA_real_
    row$variance_calibration_mcse <- if (complete) mcse((v - row$variance_calibration_ratio * q) / ev) else NA_real_
    row$primary_variance_status <- if (complete) "all_requested" else "withheld_incomplete_intervals"
    row$point_summary_provenance <- "Copied authenticated original point report after54000row binding"
    row$empirical_variance_divisor <- "R_minus_1"; row$saved_refit_B <- 200L; row$saved_refit_divisor <- "B"
    row$scope <- cfg$statistical_scope
    if (complete) stopifnot(all(is.finite(unlist(row[c("mean_reported_variance_all_requested",
      "mean_reported_variance_mcse", "variance_calibration_ratio", "variance_calibration_mcse")]))))
    rows[[mode]] <- row
  }
  metrics[[j]] <- do.call(rbind, rows)
  a <- avail$analytic; b <- avail$original_refit; ha <- hit$analytic; hb <- hit$original_refit
  delta <- as.integer(ha) - as.integer(hb); complete <- all(a & b)
  va <- d$analytic_variance; vb <- d$saved_refit_variance
  ratio <- if (complete && is.finite(mean(vb)) && mean(vb) > 0) mean(va) / mean(vb) else NA_real_
  pr <- data.frame(plan[j, ], source_method = method, requested_pairs = 1000L, common_original_ids = 1000L,
    common_point_pairs = 1000L, common_available_intervals = sum(a & b), analytic_only_available = sum(a & !b),
    refit_only_available = sum(!a & b), neither_interval_available = sum(!a & !b),
    both_cover = sum(ha & hb), analytic_only_covers = sum(ha & !hb), refit_only_covers = sum(!ha & hb),
    neither_covers = sum(!ha & !hb), paired_coverage_difference = mean(delta), paired_coverage_difference_mcse = mcse(delta),
    analytic_refit_ratio_of_mean_variances = ratio,
    ratio_of_mean_variances_mcse = if (is.finite(ratio)) mcse((va - ratio * vb) / mean(vb)) else NA_real_,
    mean_paired_variance_difference = if (complete) mean(va - vb) else NA_real_,
    paired_variance_difference_mcse = if (complete) mcse(va - vb) else NA_real_,
    variance_comparison_scope = if (!complete) "withheld_incomplete_variance_pairs" else
      if (!is.finite(ratio)) "undefined_zero_or_nonfinite_mean_ratio" else "all_requested",
    saved_refit_B = 200L, saved_refit_divisor = "B", stringsAsFactors = FALSE)
  if (complete) stopifnot(all(is.finite(unlist(pr[c("mean_paired_variance_difference", "paired_variance_difference_mcse")]))))
  if (complete && mean(vb) > 0) stopifnot(is.finite(ratio), is.finite(pr$ratio_of_mean_variances_mcse))
  paired[[j]] <- pr
}
metrics <- do.call(rbind, metrics); paired <- do.call(rbind, paired)
stopifnot(nrow(metrics) == 108L, nrow(paired) == 54L, identical(initial_rng, rng()))
records$warning_count[!records$warning_count_known] <- NA_real_
failure <- records$status != "ok" | !records$warning_count_known |
  (records$warning_count_known & records$warning_count > 0)
outputs <- list(metrics = metrics, paired_comparisons = paired, cell_audit = base$audit,
  conditional_point_diagnostics = base$point_success_diagnostics,
  failure_audit = records[failure, , drop = FALSE])
pins <- list()
for (name in names(outputs)) {
  path <- file.path(out, paste0(name, ".csv")); stopifnot(!file.exists(path))
  write.csv(outputs[[name]], path, row.names = FALSE, na = "NA"); pins[[basename(path)]] <- sha(path)
}
for (f in names(binding$output_pins)) stopifnot(identical(sha(file.path(out, f)), binding$output_pins[[f]]))
stopifnot(identical(sha(metric_source), cfg$source_pins[[cfg$paired_metrics_source]]), identical(initial_rng, rng()))
jsonlite::write_json(list(status = "COMPLETE_SAVED_SUMMARY", rows = 54000L, cases = 3000L,
  inference_summary_rows = 108L, paired_summary_rows = 54L, original_metric_rows_copied = 72L,
  records_complete = base$records_complete, reporting_inputs_valid = base$reporting_inputs_valid,
  RNG_unchanged = TRUE, no_estimator_inference_fitting_matching_bootstrap = TRUE,
  runtime = list(R = R.version.string, R_home = R.home(),
    digest = as.character(utils::packageVersion("digest")), jsonlite = as.character(utils::packageVersion("jsonlite")),
    threads = Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS"))),
  output_pins = pins, source_point_copy_fields = point_fields,
  variance_comparison = "Ratio of all-requested mean variances, not mean of individual ratios",
  scope = cfg$statistical_scope), file.path(out, "summary_completion.json"),
  auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = NA)
