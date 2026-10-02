# Publication metrics for a planned, paired-method simulation study.
# Matches SW_DSM's signed-RB and fixed CorCor / PSM_M1 reference conventions.
# This is reporting code, not a causal estimator or an inference procedure.

wm_comparison_metrics <- function(records, plan,
    group_cols = c("estimand", "design", "overlap", "n"),
    reference_scenario = "CorCor", reference_method = "PSM_M1") {
  # plan: one row per group/scenario/method, requested (replicates 1:requested),
  # and interval_required. records: those keys, replicate, target, estimate,
  # lower, upper, variance (estimator variance, not root-n variance), status.
  keys <- c(group_cols, "scenario", "method")
  if (!is.character(group_cols) || !length(group_cols) || anyNA(group_cols) ||
      anyDuplicated(keys) || any(!nzchar(keys))) stop("Invalid grouping columns")
  if (!is.data.frame(records) || !is.data.frame(plan) || !nrow(plan) ||
      !all(c(keys, "requested", "interval_required") %in% names(plan)) ||
      !all(c(keys, "replicate", "target", "estimate", "lower", "upper",
             "variance", "status") %in% names(records))) stop("Incomplete comparison schema")
  if (length(reference_scenario) != 1L || is.na(reference_scenario) ||
      length(reference_method) != 1L || is.na(reference_method)) stop("Invalid reference")
  for (nm in keys) {
    if (anyNA(plan[[nm]]) || anyNA(records[[nm]])) stop("Missing grouping value")
  }
  # Length prefixes prevent distinct character keys from colliding.
  key <- function(x, cols) {
    if (!nrow(x)) return(character())
    do.call(paste, c(lapply(x[cols], function(v) {
      v <- as.character(v); paste0(nchar(v, type = "bytes"), ":", v)
    }), sep = "|"))
  }
  pk <- key(plan, keys); rk <- key(records, keys)
  if (anyDuplicated(pk)) stop("Duplicate planned cell")
  if (any(!rk %in% pk)) stop("Unplanned records")
  numeric_cols <- c("replicate", "target", "estimate", "lower", "upper", "variance")
  if (!all(vapply(records[numeric_cols], function(x) is.numeric(x) && !is.complex(x), logical(1))))
    stop("Performance fields must be real numeric columns")
  if (!is.numeric(plan$requested) || is.complex(plan$requested) ||
      any(!is.finite(plan$requested) | plan$requested < 2 |
            plan$requested != floor(plan$requested))) stop("requested must be an integer >= 2")
  if (!is.logical(plan$interval_required) || anyNA(plan$interval_required))
    stop("interval_required must be explicit TRUE/FALSE")
  if (any(!is.finite(records$replicate) | records$replicate < 1 |
            records$replicate != floor(records$replicate)) ||
      any(records$replicate > plan$requested[match(rk, pk)])) stop("Unplanned replicate")
  if (anyDuplicated(key(records, c(keys, "replicate")))) stop("Duplicate replicate")
  if (anyNA(records$status) || any(!records$status %in% c("ok", "failed")))
    stop("Unknown status")
  if (any(!is.finite(records$target))) stop("Every record must retain its finite target")
  # A scenario changes working models, not the target within a design group.
  for (g in unique(key(records, group_cols))) {
    if (length(unique(records$target[key(records, group_cols) == g])) != 1L)
      stop("Inconsistent target within comparison group")
  }
  audit <- plan
  audit$received <- audit$completed <- audit$failed <- audit$invalid <- audit$intervals <- 0L
  audit$missing <- plan$requested
  audit$reference_present <- FALSE
  pg <- key(plan, group_cols)
  for (i in seq_len(nrow(plan))) {
    x <- records[rk == pk[i], , drop = FALSE]
    ok <- x$status == "ok"
    ci <- is.finite(x$lower) & is.finite(x$upper) & x$lower <= x$upper
    # NA variance is allowed (e.g. a point-only comparator); Inf/negative is not.
    bad_var <- !is.na(x$variance) & (!is.finite(x$variance) | x$variance < 0)
    bad_ci <- (!is.na(x$lower) | !is.na(x$upper)) & !ci
    bad <- ok & (!is.finite(x$estimate) | bad_var | bad_ci |
                   (plan$interval_required[i] & !ci))
    audit$received[i] <- nrow(x)
    audit$completed[i] <- sum(ok & !bad)
    audit$failed[i] <- sum(!ok)
    audit$invalid[i] <- sum(bad)
    audit$intervals[i] <- sum(ok & !bad & ci)
    audit$missing[i] <- plan$requested[i] - nrow(x)
    audit$reference_present[i] <- sum(pg == pg[i] & plan$scenario == reference_scenario &
                                      plan$method == reference_method) == 1L
  }
  complete <- all(audit$completed == audit$requested & audit$reference_present)
  if (!complete) return(list(records_complete = all(audit$completed == audit$requested),
    publication_ready = FALSE, audit = audit, metrics = NULL))
  metrics <- do.call(rbind, lapply(seq_len(nrow(plan)), function(i) {
    x <- records[rk == pk[i], , drop = FALSE]
    target <- x$target[1L]; bias <- mean(x$estimate) - target
    ev <- stats::var(x$estimate)
    ci <- is.finite(x$lower) & is.finite(x$upper)
    hits <- x$lower[ci] <= target & target <= x$upper[ci]
    coverage <- if (length(hits)) mean(hits) else NA_real_
    cp <- if (length(hits)) stats::binom.test(sum(hits), length(hits))$conf.int else c(NA_real_, NA_real_)
    rv <- is.finite(x$variance)
    mean_var <- if (any(rv)) mean(x$variance[rv]) else NA_real_
    data.frame(plan[i, , drop = FALSE], truth = target,
      mean_estimate = mean(x$estimate), signed_bias = bias, absolute_bias = abs(bias),
      signed_relative_bias_percent = if (target != 0) 100 * bias / abs(target) else NA_real_,
      absolute_relative_bias_percent = if (target != 0) 100 * abs(bias) / abs(target) else NA_real_,
      empirical_variance = ev, mc_sd = sqrt(ev), bias_mcse = sqrt(ev / nrow(x)),
      rmse = sqrt(mean((x$estimate - target)^2)),
      intervals = sum(ci), coverage = coverage,
      coverage_mcse = if (length(hits)) sqrt(coverage * (1 - coverage) / length(hits)) else NA_real_,
      coverage_exact95_lower = cp[1L], coverage_exact95_upper = cp[2L],
      mean_ci_length = if (any(ci)) mean(x$upper[ci] - x$lower[ci]) else NA_real_,
      coverage_denominator = if (sum(ci) == nrow(x)) "all_requested" else "available_optional_intervals",
      variance_reports = sum(rv), mean_reported_variance = mean_var,
      mean_reported_se = if (any(rv)) mean(sqrt(x$variance[rv])) else NA_real_,
      variance_calibration_ratio = if (is.finite(mean_var) && ev > 0) mean_var / ev else NA_real_,
      relative_bias_status = if (target == 0) "undefined_zero_target" else "defined",
      stringsAsFactors = FALSE)
  }))
  metrics$reference_scenario <- reference_scenario
  metrics$reference_method <- reference_method
  metrics$reference_mc_variance <- metrics$relative_efficiency <- NA_real_
  metrics$relative_efficiency_status <- "undefined_zero_or_nonfinite_variance"
  for (i in seq_len(nrow(metrics))) {
    j <- which(pg == pg[i] & plan$scenario == reference_scenario & plan$method == reference_method)
    vref <- metrics$empirical_variance[j]
    metrics$reference_mc_variance[i] <- vref
    if (is.finite(vref) && vref > 0 && is.finite(metrics$empirical_variance[i]) &&
        metrics$empirical_variance[i] > 0) {
      metrics$relative_efficiency[i] <- vref / metrics$empirical_variance[i]
      metrics$relative_efficiency_status[i] <- if (is.finite(metrics$relative_efficiency[i]))
        "defined" else "undefined_nonfinite_ratio"
    }
  }
  # Match the original complete-publication gate: undefined RE blocks a table.
  # A zero target only makes relative bias undefined; it does not block other metrics.
  finite_fields <- c("mean_estimate", "signed_bias", "absolute_bias", "empirical_variance",
    "mc_sd", "bias_mcse", "rmse")
  valid <- apply(metrics[finite_fields], 1L, function(x) all(is.finite(x)))
  valid <- valid & (metrics$truth == 0 |
    (is.finite(metrics$signed_relative_bias_percent) & is.finite(metrics$absolute_relative_bias_percent)))
  valid <- valid & (metrics$intervals == 0 | is.finite(metrics$mean_ci_length))
  valid <- valid & (metrics$variance_reports == 0 |
    (is.finite(metrics$mean_reported_variance) & is.finite(metrics$mean_reported_se) &
       is.finite(metrics$variance_calibration_ratio)))
  audit$relative_efficiency_status <- metrics$relative_efficiency_status
  audit$derived_metrics_finite <- valid
  if (any(!valid | metrics$relative_efficiency_status != "defined"))
    return(list(records_complete = TRUE, publication_ready = FALSE, audit = audit, metrics = NULL))
  rownames(metrics) <- NULL
  list(records_complete = TRUE, publication_ready = TRUE, audit = audit, metrics = metrics)
}
