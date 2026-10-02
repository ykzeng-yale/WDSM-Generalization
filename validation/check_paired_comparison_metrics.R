# Deterministic reporting fixtures; no scientific samples or model fits.
source("simulations/paired_comparison_metrics.R")
checks <- 0L
check <- function(ok) {
  checks <<- checks + 1L
  if (!isTRUE(ok)) stop("Fixture check ", checks, " failed")
}
near <- function(x, y) {
  result <- all.equal(x, y, tolerance = 1e-12, check.attributes = FALSE)
  if (!isTRUE(result)) stop("Fixture check ", checks + 1L, ": ", paste(result, collapse = "; "))
  check(TRUE)
}
reject <- function(expr, text) {
  msg <- tryCatch({ force(expr); "" }, error = conditionMessage)
  check(grepl(text, msg, fixed = TRUE))
}
plan <- data.frame(population = 1L, multiplier = 1L, response = "No", scenario = "source",
  method = c("reference", "candidate", "scaled"), requested = 4L,
  interval_required = TRUE, coverage_rule = c("strict", "inclusive", "inclusive"))
estimates <- list(c(1, 3, 5, 7), c(2, 2, 6, 8), c(4, 8, 12, 16))
records <- do.call(rbind, lapply(seq_len(nrow(plan)), function(i)
  data.frame(plan[rep(i, 4L), c("population", "multiplier", "response", "scenario", "method")],
    replicate = 1:4, pairing_id = paste0("fixed-dataset-", 1:4), target = 4,
    estimate = estimates[[i]], lower = 3, upper = 5, variance = 1,
    status = "ok", point_status = "ok", interval_status = "ok", warning_count = 0L)))
run <- function(x = records, p = plan) wm_paired_comparison_metrics(x, p,
  c("population", "multiplier", "response"), "source", "reference")
a <- run(); m <- a$metrics
check(a$records_complete && a$point_metrics_complete && a$reporting_inputs_valid)
near(m$empirical_variance, c(20/3, 9, 80/3))
near(m$signed_bias, c(0, .5, 6))
near(m$absolute_bias, c(0, .5, 6))
near(m$signed_relative_bias_percent, c(0, 12.5, 150))
near(m$relative_efficiency, c(1, 20/27, 1/4))
near(m$relative_efficiency_mcse[c(1, 3)], c(0, 0))
# Explicit paired ratio derivative on the four hand-chosen values.
q_ref <- c(12, 4/3, 4/3, 12)
q_candidate <- c(25/3, 25/3, 3, 49/3)
phi <- (q_ref - (20/27) * q_candidate) / 9
near(m$relative_efficiency_mcse[2], sqrt(sum((phi - mean(phi))^2) / 12))
near(m$paired_bias_difference[2], .5)
# Paired estimate differences (1,-1,1,1) have sample variance1.
near(m$paired_bias_difference_mcse[2], .5)
near(m$paired_mse_difference[2], 2)
# Squared-error differences (-5,3,3,7) have centered sum of squares76.
near(m$paired_mse_difference_mcse[2], sqrt(19/3))
near(m$rmse, c(sqrt(5), sqrt(7), sqrt(56)))
near(m$bias_mcse, sqrt(c(20/3, 9, 80/3) / 4))
near(m$operational_coverage, c(1, 1, 1))
check(all(m$coverage_exact95_lower < 1))
near(m$variance_calibration_ratio, c(3/20, 1/9, 3/80))
near(run(records[rev(seq_len(nrow(records))), ])$metrics, m)

# Pair covariance matters even when all marginal summaries are unchanged.
permuted <- records; permuted$estimate[5:8] <- c(2, 6, 8, 2)
mp <- run(permuted)$metrics
near(mp$relative_efficiency, m$relative_efficiency)
check(abs(mp$relative_efficiency_mcse[2] - m$relative_efficiency_mcse[2]) > 1e-3)

# Interval failure keeps a valid point, counts operational noncoverage, and
# exposes the conditional coverage denominator without selecting point results.
x <- records; x$status[6] <- "failed"; x$interval_status[6] <- "failed_interval"
x$variance[6] <- x$lower[6] <- x$upper[6] <- NA_real_
f <- run(x)
check(f$records_complete && f$point_metrics_complete && f$reporting_inputs_valid)
near(f$metrics$relative_efficiency, m$relative_efficiency)
near(f$metrics$operational_coverage[2], .75)
near(f$metrics$conditional_coverage[2], 1)
near(f$metrics$conditional_coverage_denominator[2], 3)
near(f$metrics$interval_availability[2], .75)
near(f$metrics$paired_operational_coverage_difference[2], -.25)
near(f$metrics$paired_operational_coverage_difference_mcse[2], .25)
check(is.na(f$metrics$variance_calibration_ratio[2]))
near(f$metrics$operational_coverage_mcse[2], sqrt(3)/8)

# Missing and explicitly unrun rows are not failed replications.
missing <- run(records[-6, ])
check(!missing$records_complete && missing$audit$missing[2] == 1L)
check(is.na(missing$metrics$operational_coverage[2]))
check(is.na(missing$metrics$signed_bias[2]))
check(is.na(missing$metrics$relative_efficiency[2]))
near(missing$point_success_diagnostics$denominator[2], 3)
unrun <- x; unrun$status[6] <- "not_run"; unrun$point_status[6] <- "unavailable"
unrun$estimate[6] <- NA_real_; unrun$interval_status[6] <- "unavailable"
u <- run(unrun)
check(!u$records_complete && u$audit$not_run[2] == 1L && u$audit$missing[2] == 0L)
near(u$audit$point_unavailable[2], 0)
check(is.na(u$metrics$operational_coverage[2]))

# An observed point failure withholds unconditional point summaries, while
# complete-record operational coverage still includes that failed replication.
x$point_status[6] <- "unavailable"; x$estimate[6] <- NA_real_
pf <- run(x)
check(pf$records_complete && !pf$point_metrics_complete && pf$reporting_inputs_valid)
check(is.na(pf$metrics$rmse[2]) && is.na(pf$metrics$paired_bias_difference[2]))
near(pf$metrics$operational_coverage[2], .75)
near(pf$point_success_diagnostics$denominator[2], 3)
all_failed <- records
all_failed$status <- "failed"; all_failed$point_status <- all_failed$interval_status <- "unavailable"
all_failed$estimate <- all_failed$variance <- all_failed$lower <- all_failed$upper <- NA_real_
af <- run(all_failed)
check(af$records_complete && !af$point_metrics_complete && af$reporting_inputs_valid)
near(af$metrics$operational_coverage, rep(0, 3))
check(all(af$metrics$coverage_exact95_upper > 0))
check(all(is.na(af$metrics$conditional_coverage)))
empty <- run(records[FALSE, ])
check(!empty$records_complete && empty$reporting_inputs_valid)
check(all(empty$audit$missing == 4) && all(is.na(empty$metrics$operational_coverage)))

# Endpoint convention is explicit and any native coverage field is checked.
endpoint <- records; endpoint$lower <- 4; endpoint$upper <- 5
ep <- run(endpoint)
near(ep$metrics$operational_coverage, c(0, 1, 1))
endpoint$source_coverage <- c(rep(0, 4), rep(1, 8))
near(run(endpoint)$metrics$operational_coverage, c(0, 1, 1))
endpoint$source_coverage[1] <- 1
reject(run(endpoint), "Source coverage disagrees")
negative <- records; negative$target <- -4
near(run(negative)$metrics$signed_relative_bias_percent[1], 200)
zero <- records; zero$target <- 0
z <- run(zero)
check(all(is.na(z$metrics$signed_relative_bias_percent)))
check(all(z$metrics$relative_bias_status == "undefined_zero_target"))
flat <- records; flat$estimate[5:8] <- 4
flat_result <- run(flat)
check(is.na(flat_result$metrics$relative_efficiency[2]))
check(is.na(flat_result$metrics$rmse_mcse[2]) && flat_result$metrics$rmse[2] == 0)
check(flat_result$reporting_inputs_valid)
warnings <- records; warnings$warning_count[6] <- 2L
near(run(warnings)$metrics$signed_bias, m$signed_bias)
near(run(warnings)$metrics$warning_rate[2], .25)

# Claimed successes with invalid arithmetic remain visible and invalidate
# reporting; finite points are not thrown away because their intervals failed.
bad <- records; bad$upper[6] <- NA_real_
b <- run(bad)
check(!b$reporting_inputs_valid && b$audit$invalid_intervals[2] == 1L)
near(b$metrics$signed_bias[2], .5)
bad <- records; bad$estimate[6] <- Inf
check(!run(bad)$reporting_inputs_valid)
bad <- records; bad$variance[6] <- -1
check(!run(bad)$reporting_inputs_valid)
bad <- records; bad$estimate[6] <- 1e300
check(!run(bad)$reporting_inputs_valid)
bad <- records; bad$lower[6] <- -1e308; bad$upper[6] <- 1e308
check(!run(bad)$reporting_inputs_valid)

# Pairing means actual dataset identity, not a shared row or seed label.
wrong <- records; wrong$pairing_id[6] <- "different-realized-data"
reject(run(wrong), "Dataset pairing identifiers disagree")
wrong <- records; wrong$status[6] <- "pending"
reject(run(wrong), "Unknown execution status")
wrong <- records; wrong$error <- ""; wrong$error[6] <- "Not evaluated"
reject(run(wrong), "Unevaluated placeholder")
reject(run(rbind(records, records[1, ])), "Duplicate replicate")
reject(run(records, rbind(plan, plan[1, ])), "Duplicate planned cell")
wrong <- records; wrong$target[1] <- 5
reject(run(wrong), "Inconsistent target")
wrong <- records; wrong$replicate[1] <- 5
reject(run(wrong), "Unplanned replicate")
wrong <- records; wrong$method[1] <- "unplanned"
reject(run(wrong), "Unplanned cell")
wrong <- plan; wrong$requested[2] <- 5
reject(run(records, wrong), "same prescribed replicate IDs")

# Independent-review counterexamples: derived overflow, representational
# underflow, terminal statuses, and duplicate schema names.
tiny <- records; tiny$target <- 0; tiny$estimate <- rep(c(-2,-1,1,2) * 1e-50, 3)
tiny$variance <- 1e300
overflow <- run(tiny)
check(!overflow$reporting_inputs_valid && is.infinite(overflow$metrics$variance_calibration_ratio[2]))
large <- records; large$target <- 0
large$estimate <- sqrt(1.2e154) * c(0,0,1,1, 1,1,0,0, 0,0,1,1)
paired_overflow <- run(large)
check(!paired_overflow$reporting_inputs_valid)
check(!paired_overflow$audit$derived_paired_metrics_finite[2])
underflow <- records; underflow$target <- 0
underflow$estimate <- c(c(-2,-1,1,2)*1e-150, c(-2,-1,1,2)*1e50, c(-2,-1,1,2))
under <- run(underflow)
check(is.na(under$metrics$relative_efficiency[2]))
check(under$metrics$relative_efficiency_status[2] == "undefined_nonrepresentable_ratio_or_mcse")
check(!under$reporting_inputs_valid)
pending <- records; pending$status[6] <- "failed"
pending$point_status[6] <- pending$interval_status[6] <- "pending"
pending$estimate[6] <- pending$lower[6] <- pending$upper[6] <- pending$variance[6] <- NA_real_
reject(run(pending), "Unknown or nonterminal component status")
contradiction <- records; contradiction$status[6] <- "failed"
check(!run(contradiction)$reporting_inputs_valid)
point_only <- plan; point_only$interval_required[2] <- FALSE
point_records <- records; point_records$interval_status[5:8] <- "not_requested"
point_records$lower[5:8] <- point_records$upper[5:8] <- point_records$variance[5:8] <- NA_real_
po <- run(point_records, point_only)
check(po$reporting_inputs_valid && po$point_metrics_complete)
check(is.na(po$metrics$operational_coverage[2]))
point_records$status[6] <- "failed"
check(!run(point_records, point_only)$reporting_inputs_valid)
duplicate <- cbind(records, second = Inf); names(duplicate)[ncol(duplicate)] <- "estimate"
reject(run(duplicate), "Invalid or duplicate column names")
duplicate <- cbind(plan, second = 4L); names(duplicate)[ncol(duplicate)] <- "requested"
reject(run(records, duplicate), "Invalid or duplicate column names")
wrong <- records; names(wrong)[1] <- ""
reject(run(wrong), "Invalid or duplicate column names")
wrong <- records; names(wrong)[1] <- NA_character_
reject(run(wrong), "Invalid or duplicate column names")
wrong <- records; wrong$point_status <- 1
reject(run(wrong), "Missing or noncharacter outcome status")
reject(wm_paired_comparison_metrics(records, plan, c("population", "multiplier", "response"),
  "source", " "), "Invalid reference")
cat("PASS", checks, "deterministic paired-metric, failure-accounting and malformed-input checks\n")
