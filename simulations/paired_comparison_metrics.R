# Paired Monte Carlo reporting with explicit failure and missing-run accounting.
# This adds no estimator, variance theorem, or data-generating process.

wm_paired_comparison_metrics <- function(records, plan, group_cols,
    reference_scenario, reference_method) {
  for (x in list(records, plan))
    if (!is.data.frame(x) || is.null(names(x)) || anyNA(names(x)) ||
        any(!nzchar(names(x))) || anyDuplicated(names(x))) stop("Invalid or duplicate column names")
  keys <- c(group_cols, "scenario", "method")
  fields <- c("replicate", "target", "estimate", "variance", "lower", "upper",
              "status", "point_status", "interval_status", "warning_count", "pairing_id")
  if (!is.character(group_cols) || !length(group_cols) || anyNA(keys) ||
      anyDuplicated(keys) || any(!nzchar(keys)) ||
      !is.data.frame(plan) || !nrow(plan) || !is.data.frame(records) ||
      !all(c(keys, "requested", "interval_required", "coverage_rule") %in% names(plan)) ||
      !all(c(keys, fields) %in% names(records))) stop("Invalid comparison schema")
  if (!is.character(reference_scenario) || length(reference_scenario) != 1L ||
      !is.character(reference_method) || length(reference_method) != 1L ||
      anyNA(c(reference_scenario, reference_method)) ||
      any(!nzchar(trimws(c(reference_scenario, reference_method))))) stop("Invalid reference")
  for (k in keys) if (anyNA(plan[[k]]) || anyNA(records[[k]])) stop("Missing cell key")
  numeric_fields <- c("replicate", "target", "estimate", "variance", "lower", "upper", "warning_count")
  if (!all(vapply(records[numeric_fields], function(x) is.numeric(x) && !is.complex(x), logical(1))))
    stop("Performance fields must be real numeric columns")
  if (!is.numeric(plan$requested) || is.complex(plan$requested) ||
      any(!is.finite(plan$requested) | plan$requested < 2 | plan$requested != floor(plan$requested)))
    stop("requested must be an integer >=2")
  if (!is.logical(plan$interval_required) || anyNA(plan$interval_required) ||
      anyNA(plan$coverage_rule) || any(!plan$coverage_rule %in% c("strict", "inclusive")))
    stop("Interval requirement and coverage endpoint rule must be explicit")
  for (k in c("status", "point_status", "interval_status"))
    if (!is.character(records[[k]]) || anyNA(records[[k]]) || any(!nzchar(records[[k]])))
      stop("Missing or noncharacter outcome status")
  if (any(!records$status %in% c("ok", "failed", "not_run"))) stop("Unknown execution status")
  if (any(!records$point_status %in% c("ok", "unavailable", "failed", "failed_point")) ||
      any(!records$interval_status %in% c("ok", "unavailable", "failed", "failed_interval", "not_requested")))
    stop("Unknown or nonterminal component status")
  if (!is.character(records$pairing_id) || anyNA(records$pairing_id) || any(!nzchar(records$pairing_id)))
    stop("Actual dataset pairing identifiers are required")
  if ("error" %in% names(records) && any(records$status != "not_run" &
      !is.na(records$error) & records$error == "Not evaluated")) stop("Unevaluated placeholder is not a method failure")
  if (any(!is.finite(records$warning_count) | records$warning_count < 0 |
          records$warning_count != floor(records$warning_count))) stop("Invalid warning count")
  if (any(!is.finite(records$target))) stop("Every received record must retain a finite target")
  key <- function(x, cols) {
    if (!nrow(x)) return(character())
    do.call(paste, c(unname(lapply(x[cols], function(z) {
      z <- as.character(z); paste0(nchar(z, type = "bytes"), ":", z)
    })), list(sep = "|")))
  }
  pk <- key(plan, keys); rk <- key(records, keys)
  if (anyDuplicated(pk)) stop("Duplicate planned cell")
  if (any(!rk %in% pk)) stop("Unplanned cell")
  cell <- match(rk, pk)
  if (any(!is.finite(records$replicate) | records$replicate < 1 |
          records$replicate != floor(records$replicate)) ||
      any(records$replicate > plan$requested[cell])) stop("Unplanned replicate")
  if (anyDuplicated(key(records, c(keys, "replicate")))) stop("Duplicate replicate")
  groups <- key(plan, group_cols)
  for (g in unique(groups)) {
    ix <- which(groups == g)
    targets <- unique(records$target[cell %in% ix])
    if (length(targets) > 1L) stop("Inconsistent target within comparison group")
  }
  ref_plan <- plan; ref_plan$scenario <- reference_scenario; ref_plan$method <- reference_method
  ref_index <- match(key(ref_plan, keys), pk)
  has_ref <- !is.na(ref_index)
  if (any(plan$requested[has_ref] != plan$requested[ref_index[has_ref]]))
    stop("Paired methods must have the same prescribed replicate IDs")
  # Split once: formal Lenis has1.8million rows; repeated full-table filtering
  # for each method cell would make aggregation needlessly quadratic.
  row_index <- split(seq_len(nrow(records)), factor(cell, levels = seq_len(nrow(plan))))
  audit <- plan
  for (nm in c("received", "missing", "evaluated", "not_run", "point_available", "point_unavailable",
               "invalid_points", "interval_available", "interval_unavailable",
               "invalid_intervals", "invalid_status", "warning_replicates")) audit[[nm]] <- 0L
  audit$reference_present <- has_ref
  audit$calibration_expected <- audit$paired_point_expected <- audit$paired_coverage_expected <- FALSE
  metric_rows <- diagnostic_rows <- vector("list", nrow(plan))
  paired <- vector("list", nrow(plan))
  mean_or_na <- function(x) if (length(x)) mean(x) else NA_real_
  coverage_summary <- function(hit) {
    if (!length(hit)) return(c(rate = NA_real_, mcse = NA_real_, low = NA_real_, high = NA_real_))
    p <- mean(hit); ci <- stats::binom.test(sum(hit), length(hit))$conf.int
    c(rate = p, mcse = sqrt(p * (1 - p) / length(hit)), low = ci[1L], high = ci[2L])
  }
  for (i in seq_len(nrow(plan))) {
    x <- records[row_index[[i]], , drop = FALSE]
    x <- x[order(x$replicate), , drop = FALSE]
    R <- plan$requested[i]; n <- nrow(x)
    evaluated <- x$status != "not_run"
    point <- evaluated & x$point_status == "ok" & is.finite(x$estimate)
    invalid_point <- x$point_status == "ok" & !is.finite(x$estimate)
    interval_numeric <- is.finite(x$lower) & is.finite(x$upper) & x$lower <= x$upper &
      is.finite(x$variance) & x$variance >= 0 & point
    interval <- x$interval_status == "ok" & interval_numeric
    invalid_interval <- x$interval_status == "ok" & !interval_numeric
    hit <- rep(FALSE, n)
    if (plan$coverage_rule[i] == "strict")
      hit[interval] <- x$lower[interval] < x$target[interval] & x$target[interval] < x$upper[interval]
    else hit[interval] <- x$lower[interval] <= x$target[interval] & x$target[interval] <= x$upper[interval]
    if ("source_coverage" %in% names(x)) {
      declared <- interval & !is.na(x$source_coverage)
      if (any(x$source_coverage[declared] != as.numeric(hit[declared])))
        stop("Source coverage disagrees with the declared endpoint rule")
    }
    audit$received[i] <- n; audit$missing[i] <- R - n
    audit$evaluated[i] <- sum(evaluated); audit$not_run[i] <- sum(!evaluated)
    audit$point_available[i] <- sum(point); audit$point_unavailable[i] <- sum(evaluated & !point)
    audit$invalid_points[i] <- sum(invalid_point)
    audit$interval_available[i] <- sum(interval); audit$interval_unavailable[i] <- sum(evaluated & !interval)
    audit$invalid_intervals[i] <- sum(invalid_interval)
    audit$invalid_status[i] <- sum((x$status == "ok" & (!point | (plan$interval_required[i] & !interval))) |
      (x$status == "failed" & point & (!plan$interval_required[i] | interval)) |
      (!evaluated & (x$point_status != "unavailable" |
        !x$interval_status %in% c("unavailable", "not_requested"))) |
      (plan$interval_required[i] & x$interval_status == "not_requested"))
    audit$warning_replicates[i] <- sum(evaluated & x$warning_count > 0)
    complete <- n == R && all(evaluated)
    complete_points <- complete && all(point)
    complete_intervals <- complete && all(interval)
    target <- if (n) x$target[1L] else NA_real_
    estimates <- x$estimate[point]
    error <- estimates - target
    ev <- if (length(estimates) >= 2L) stats::var(estimates) else NA_real_
    audit$calibration_expected[i] <- complete_points && complete_intervals && is.finite(ev) && ev > 0
    bias <- mean_or_na(error); mse <- mean_or_na(error^2); rmse <- sqrt(mse)
    values <- c(mean_estimate = mean_or_na(estimates), signed_bias = bias, absolute_bias = abs(bias),
      signed_relative_bias_percent = if (is.finite(target) && target != 0) 100 * bias / abs(target) else NA_real_,
      empirical_variance = ev, empirical_sd = sqrt(ev), bias_mcse = sqrt(ev / length(estimates)),
      signed_relative_bias_mcse = if (is.finite(target) && target != 0) 100 * sqrt(ev / length(estimates)) / abs(target) else NA_real_,
      mse = mse, rmse = rmse,
      rmse_mcse = if (length(estimates) >= 2L && is.finite(rmse) && rmse > 0)
        stats::sd(error^2) / (2 * sqrt(length(estimates)) * rmse) else NA_real_)
    if (!length(estimates)) values[] <- NA_real_
    # The only successful-point-only summaries are in this separate diagnostic table.
    diagnostic_rows[[i]] <- data.frame(plan[i, , drop = FALSE], truth = target,
      denominator = length(estimates), scope = "conditional_on_available_point",
      as.list(values), stringsAsFactors = FALSE)
    if (!complete_points) values[] <- NA_real_
    q <- if (complete_points) R / (R - 1) * (x$estimate - mean(x$estimate))^2 else rep(NA_real_, n)
    cv <- if (complete && plan$interval_required[i]) coverage_summary(hit) else coverage_summary(logical())
    av <- if (complete && plan$interval_required[i]) coverage_summary(interval) else coverage_summary(logical())
    conditional <- if (plan$interval_required[i]) coverage_summary(hit[interval]) else coverage_summary(logical())
    metric_rows[[i]] <- data.frame(plan[i, , drop = FALSE], truth = target,
      point_metric_scope = if (complete_points) "all_requested" else "withheld_incomplete_points",
      relative_bias_status = if (!complete_points) "withheld_incomplete_points" else if (target == 0)
        "undefined_zero_target" else "defined",
      as.list(values), variance_mcse = if (complete_points) stats::sd(q) / sqrt(R) else NA_real_,
      interval_availability = unname(av["rate"]), availability_mcse = unname(av["mcse"]),
      operational_coverage = unname(cv["rate"]), operational_coverage_mcse = unname(cv["mcse"]),
      coverage_exact95_lower = unname(cv["low"]), coverage_exact95_upper = unname(cv["high"]),
      operational_denominator = if (complete && plan$interval_required[i]) R else NA_integer_,
      conditional_coverage = unname(conditional["rate"]), conditional_coverage_mcse = unname(conditional["mcse"]),
      conditional_coverage_denominator = if (plan$interval_required[i]) sum(interval) else 0L,
      mean_valid_ci_length = mean_or_na(x$upper[interval] - x$lower[interval]),
      mean_reported_variance_valid_intervals = mean_or_na(x$variance[interval]),
      variance_calibration_ratio = if (complete_points && complete_intervals && ev > 0)
        mean(x$variance) / ev else NA_real_,
      warning_rate = if (complete) mean(x$warning_count > 0) else NA_real_,
      stringsAsFactors = FALSE)
    paired[[i]] <- list(complete = complete, complete_points = complete_points,
      ids = x$replicate, pairing_id = x$pairing_id,
      estimate = x$estimate, squared_error = (x$estimate - target)^2,
      q = q, hit = hit, variance = if (complete_points) ev else NA_real_)
  }
  metrics <- do.call(rbind, metric_rows); diagnostics <- do.call(rbind, diagnostic_rows)
  metrics$reference_scenario <- reference_scenario; metrics$reference_method <- reference_method
  for (nm in c("reference_empirical_variance", "relative_efficiency", "relative_efficiency_mcse",
               "paired_bias_difference", "paired_bias_difference_mcse", "paired_mse_difference",
               "paired_mse_difference_mcse", "paired_operational_coverage_difference",
               "paired_operational_coverage_difference_mcse")) metrics[[nm]] <- NA_real_
  metrics$relative_efficiency_status <- "missing_reference"
  for (i in seq_len(nrow(plan))) {
    if (!has_ref[i]) next
    j <- ref_index[i]; a <- paired[[i]]; b <- paired[[j]]; R <- plan$requested[i]
    common <- intersect(a$ids, b$ids)
    if (!identical(a$pairing_id[match(common, a$ids)], b$pairing_id[match(common, b$ids)]))
      stop("Dataset pairing identifiers disagree")
    metrics$relative_efficiency_status[i] <- "withheld_incomplete_paired_points"
    if (a$complete && b$complete && plan$interval_required[i] && plan$interval_required[j]) {
      audit$paired_coverage_expected[i] <- TRUE
      stopifnot(identical(a$ids, b$ids))
      delta <- as.numeric(a$hit) - as.numeric(b$hit)
      metrics$paired_operational_coverage_difference[i] <- mean(delta)
      metrics$paired_operational_coverage_difference_mcse[i] <- stats::sd(delta) / sqrt(R)
    }
    if (!a$complete_points || !b$complete_points) next
    audit$paired_point_expected[i] <- TRUE
    stopifnot(identical(a$ids, b$ids))
    delta <- a$estimate - b$estimate
    squared_delta <- a$squared_error - b$squared_error
    metrics$paired_bias_difference[i] <- mean(delta)
    metrics$paired_bias_difference_mcse[i] <- stats::sd(delta) / sqrt(R)
    metrics$paired_mse_difference[i] <- mean(squared_delta)
    metrics$paired_mse_difference_mcse[i] <- stats::sd(squared_delta) / sqrt(R)
    metrics$reference_empirical_variance[i] <- b$variance
    metrics$relative_efficiency_status[i] <- "undefined_zero_or_nonfinite_variance"
    if (is.finite(a$variance) && a$variance > 0 && is.finite(b$variance) && b$variance > 0) {
      ratio <- b$variance / a$variance
      influence <- (b$q - ratio * a$q) / a$variance
      mcse <- stats::sd(influence) / sqrt(R)
      if (is.finite(ratio) && ratio > 0 && is.finite(mcse) && mcse >= 0) {
        metrics$relative_efficiency[i] <- ratio
        metrics$relative_efficiency_mcse[i] <- mcse
        metrics$relative_efficiency_status[i] <- "defined"
      } else metrics$relative_efficiency_status[i] <- "undefined_nonrepresentable_ratio_or_mcse"
    }
  }
  # Arithmetic overflow must remain visible, not turn into a plausible table.
  performance <- c("mean_estimate", "signed_bias", "absolute_bias", "empirical_variance",
    "empirical_sd", "bias_mcse", "mse", "rmse", "variance_mcse")
  audit$derived_point_metrics_finite <- metrics$point_metric_scope != "all_requested" |
    apply(metrics[performance], 1L, function(x) all(is.finite(x)))
  audit$derived_point_metrics_finite <- audit$derived_point_metrics_finite &
    (metrics$point_metric_scope != "all_requested" | metrics$truth == 0 |
      (is.finite(metrics$signed_relative_bias_percent) & is.finite(metrics$signed_relative_bias_mcse))) &
    (metrics$point_metric_scope != "all_requested" | metrics$rmse == 0 | is.finite(metrics$rmse_mcse))
  audit$derived_interval_metrics_finite <- audit$interval_available == 0L |
    (is.finite(metrics$mean_valid_ci_length) & is.finite(metrics$mean_reported_variance_valid_intervals))
  audit$derived_interval_metrics_finite <- audit$derived_interval_metrics_finite &
    (!audit$calibration_expected | is.finite(metrics$variance_calibration_ratio))
  paired_fields <- c("paired_bias_difference", "paired_bias_difference_mcse",
    "paired_mse_difference", "paired_mse_difference_mcse", "reference_empirical_variance")
  audit$derived_paired_metrics_finite <- (!audit$paired_point_expected |
    apply(metrics[paired_fields], 1L, function(x) all(is.finite(x)))) &
    (!audit$paired_coverage_expected | (is.finite(metrics$paired_operational_coverage_difference) &
      is.finite(metrics$paired_operational_coverage_difference_mcse))) &
    metrics$relative_efficiency_status != "undefined_nonrepresentable_ratio_or_mcse"
  diag_basic <- c("mean_estimate", "signed_bias", "absolute_bias", "mse", "rmse")
  diag_se <- c("empirical_variance", "empirical_sd", "bias_mcse")
  audit$derived_diagnostic_metrics_finite <- (diagnostics$denominator == 0L |
    apply(diagnostics[diag_basic], 1L, function(x) all(is.finite(x)))) &
    (diagnostics$denominator < 2L | apply(diagnostics[diag_se], 1L, function(x) all(is.finite(x)))) &
    (diagnostics$denominator == 0L | diagnostics$truth == 0 |
      is.finite(diagnostics$signed_relative_bias_percent)) &
    (diagnostics$denominator < 2L | diagnostics$truth == 0 |
      is.finite(diagnostics$signed_relative_bias_mcse)) &
    (diagnostics$denominator < 2L | diagnostics$rmse == 0 | is.finite(diagnostics$rmse_mcse))
  rownames(audit) <- rownames(metrics) <- rownames(diagnostics) <- NULL
  list(records_complete = all(audit$missing == 0L & audit$not_run == 0L),
    point_metrics_complete = all(metrics$point_metric_scope == "all_requested"),
    reporting_inputs_valid = all(audit$invalid_points == 0L & audit$invalid_intervals == 0L &
      audit$invalid_status == 0L & audit$derived_point_metrics_finite & audit$derived_interval_metrics_finite &
      audit$derived_paired_metrics_finite & audit$derived_diagnostic_metrics_finite),
    audit = audit, metrics = metrics,
    point_success_diagnostics = diagnostics,
    scope = "Monte Carlo reporting; all method failures retained; not a theory or publication approval")
}
