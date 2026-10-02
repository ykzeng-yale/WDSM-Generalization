# Deterministic reporting checks, not Monte Carlo sampling experiments.
source("simulations/comparison_metrics.R")
equal <- function(a, b) stopifnot(isTRUE(all.equal(a, b, tolerance = 1e-12, check.attributes = FALSE)))
reject <- function(expr, pattern) {
  msg <- tryCatch({force(expr); NA_character_}, error = conditionMessage)
  stopifnot(!is.na(msg), grepl(pattern, msg, fixed = TRUE))
}
plan <- data.frame(estimand = "PATE", design = "iid_response", overlap = "good", n = 200L,
  scenario = c("CorCor", "CorCor", "MisMis", "MisMis"),
  method = c("PSM_M1", "WM_M3", "PSM_M1", "WM_M3"), requested = 4L,
  interval_required = TRUE)
est <- list(c(-5, -3, -1, 1), c(-3, -2, -1, 0), c(-8, -4, 0, 4), c(-2, -1, 0, 1))
records <- do.call(rbind, lapply(1:4, function(i) data.frame(plan[rep(i, 4), 1:6],
  replicate = 1:4, target = -2, estimate = est[[i]], lower = est[[i]] - 1,
  upper = est[[i]] + 1, variance = 1, status = "ok")))
ans <- wm_comparison_metrics(records, plan); m <- ans$metrics
stopifnot(ans$publication_ready, all(ans$audit$completed == 4))
# Hand moments: baseline centered deviations (-3,-1,1,3), sum squares 20.
equal(m$empirical_variance, c(20/3, 5/3, 80/3, 5/3))
equal(m$signed_relative_bias_percent, c(0, 25, 0, 75))
equal(m$absolute_relative_bias_percent, c(0, 25, 0, 75))
# Mean absolute errors would be positive for row 1; aggregate absolute bias is 0.
equal(m$absolute_bias[1], 0)
equal(m$relative_efficiency, c(1, 4, 1/4, 4))
equal(m$reference_mc_variance, rep(20/3, 4))
equal(m$variance_calibration_ratio, c(3/20, 3/5, 3/80, 3/5))
equal(m$coverage, c(1/2, 3/4, 0, 1/2))
equal(m$coverage_mcse[2], sqrt(3)/8)
equal(m$bias_mcse[1], sqrt(5/3))
equal(m$rmse, c(sqrt(5), sqrt(1.5), sqrt(20), sqrt(3.5)))
equal(m$mean_reported_se, rep(1, 4))
stopifnot(all(m$coverage_denominator == "all_requested"))
stopifnot(m$coverage_exact95_upper[3] > 0) # zero observed coverage != no MC uncertainty
# Target sign and target zero: absolute denominator, no invented percentage.
positive <- records; positive$target <- 2
mp <- wm_comparison_metrics(positive, plan)$metrics
equal(mp$signed_relative_bias_percent[1], -200)
zero <- records; zero$target <- 0
mz <- wm_comparison_metrics(zero, plan)$metrics
stopifnot(all(is.na(mz$signed_relative_bias_percent)), all(mz$relative_bias_status == "undefined_zero_target"))
# Do not silently publish successful-only metrics or count a missing CI as noncoverage.
failed <- records; failed$status[2] <- "failed"; failed$estimate[2] <- NA_real_
f <- wm_comparison_metrics(failed, plan)
stopifnot(!f$publication_ready, is.null(f$metrics), f$audit$failed[1] == 1, f$audit$completed[1] == 3)
missing <- wm_comparison_metrics(records[-2, ], plan)
stopifnot(!missing$publication_ready, missing$audit$missing[1] == 1)
no_ci <- records; no_ci$upper[2] <- NA_real_
bad <- wm_comparison_metrics(no_ci, plan)
stopifnot(!bad$publication_ready, bad$audit$invalid[1] == 1)
badvar <- records; badvar$variance[2] <- -1
stopifnot(!wm_comparison_metrics(badvar, plan)$publication_ready)
no_reference <- plan; no_reference$method[1] <- "OTHER_M1"
nr <- records; nr$method[nr$scenario == "CorCor" & nr$method == "PSM_M1"] <- "OTHER_M1"
stopifnot(!wm_comparison_metrics(nr, no_reference)$publication_ready)
# Entirely point-only comparator: interval absence is declared in advance.
pp <- plan; pp$interval_required[2] <- FALSE
pr <- records; pr$lower[5:8] <- pr$upper[5:8] <- pr$variance[5:8] <- NA_real_
pm <- wm_comparison_metrics(pr, pp)$metrics
stopifnot(pm$intervals[2] == 0, is.na(pm$coverage[2]), is.na(pm$variance_calibration_ratio[2]))
flat <- records; flat$estimate[5:8] <- 1
fm <- wm_comparison_metrics(flat, plan)
stopifnot(fm$records_complete, !fm$publication_ready, is.null(fm$metrics),
          fm$audit$relative_efficiency_status[2] == "undefined_zero_or_nonfinite_variance")
huge <- records; huge$estimate[1] <- 1e300
hm <- wm_comparison_metrics(huge, plan)
stopifnot(hm$records_complete, !hm$publication_ready, !hm$audit$derived_metrics_finite[1])
reject(wm_comparison_metrics(rbind(records, records[1, ]), plan), "Duplicate replicate")
reject(wm_comparison_metrics(records, rbind(plan, plan[1, ])), "Duplicate planned cell")
wrong <- records; wrong$target[1] <- -3
reject(wm_comparison_metrics(wrong, plan), "Inconsistent target")
wrong <- records; wrong$replicate[1] <- 5
reject(wm_comparison_metrics(wrong, plan), "Unplanned replicate")
wrong <- records; wrong$method[1] <- "not_planned"
reject(wm_comparison_metrics(wrong, plan), "Unplanned records")
equal(wm_comparison_metrics(records[rev(seq_len(nrow(records))), ], plan)$metrics, m)
empty <- wm_comparison_metrics(records[FALSE, ], plan)
stopifnot(!empty$publication_ready, all(empty$audit$missing == 4))
cat("PASS: hand-derived signed/absolute RB, fixed CorCor PSM RE, variance calibration, coverage denominators, strict publication gate, and malformed-record rejection.\n")
