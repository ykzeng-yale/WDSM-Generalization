#!/usr/bin/env Rscript
# Source first_stage.R and first_stage_benchmarks.R before calling this helper.
# CLI: summary.csv output.csv [geometry.rds ...]
wm_first_stage_join <- function(summary, geometry = list(), geometry_relative_goal = .01) {
  required <- c("scenario", "branch", "n", "m", "d", "M", "estimand", "comparison", "target",
    "inference_status", "requested", "completed", "failed", "empirical_root_n_variance",
    "empirical_root_n_variance_mcse_asymptotic", "mean_root_n_variance", "mean_root_n_variance_mcse")
  if (!is.data.frame(summary) || !nrow(summary) || !all(required %in% names(summary)) ||
      anyNA(summary[c("scenario", "branch", "estimand", "comparison")]) ||
      anyDuplicated(summary[c("scenario", "estimand", "comparison")])) stop("invalid or duplicate first-stage summary")
  if (!is.numeric(geometry_relative_goal) || length(geometry_relative_goal) != 1 ||
      !is.finite(geometry_relative_goal) || geometry_relative_goal <= 0) stop("invalid geometry precision goal")
  if (inherits(geometry, "wm_geometry") || inherits(geometry, "wm_geometry_alpha"))
    geometry <- setNames(list(geometry), sprintf("M%d_d%d", geometry$M, geometry$d))
  if (!is.list(geometry) || (length(geometry) && (is.null(names(geometry)) || anyDuplicated(names(geometry)))))
    stop("geometry must be a named list or individual geometry object")
  for (s in unique(summary$scenario)) {
    r <- summary[summary$scenario == s, , drop = FALSE]
    for (field in c("branch", "n", "m", "d", "M", "target", "requested"))
      if (anyNA(r[[field]]) || length(unique(r[[field]])) != 1L) stop("inconsistent scenario summary")
    if (nrow(r) != 6L || !setequal(paste(r$estimand, r$comparison), as.vector(outer(c("PATE", "PATT"),
      c("known_first_stage", "fitted_naive", "fitted_adjusted"), paste)))) stop("incomplete summary comparison grid")
  }
  for (f in c("requested", "completed", "failed"))
    if (!is.numeric(summary[[f]]) || any(!is.finite(summary[[f]])) || any(summary[[f]] < 0) ||
        any(summary[[f]] != floor(summary[[f]]))) stop("invalid summary counts")
  if (any(summary$requested != summary$completed + summary$failed) || any(summary$requested < 1)) stop("inconsistent failure counts")
  for (f in c("empirical_root_n_variance", "empirical_root_n_variance_mcse_asymptotic", "mean_root_n_variance", "mean_root_n_variance_mcse"))
    if (!is.numeric(summary[[f]]) || any(!is.na(summary[[f]]) & (!is.finite(summary[[f]]) | summary[[f]] < 0)))
      stop("invalid variance summary or MCSE")
  added <- vector("list", nrow(summary))
  for (i in seq_len(nrow(summary))) {
    r <- summary[i, ]; key <- sprintf("M%d_d%d", r$M, r$d)
    geo <- geometry[[key]]
    result <- wm_first_stage_benchmark(r$branch, r$d, r$M, r$n, r$m, geo)
    b <- result$table[result$table$estimand == r$estimand & result$table$comparison == r$comparison, , drop = FALSE]
    if (nrow(b) != 1L || !isTRUE(all.equal(r$target, b$target)) || r$inference_status != b$inference_status)
      stop("summary target, comparison or inference scope disagrees with benchmark")
    if (key %in% names(geometry) && is.null(geo) && r$d != 1L) b$geometry_status <- "geometry_run_failed"
    retain <- setdiff(names(b), c("branch", "d", "M", "n", "m", "estimand", "comparison", "target", "inference_status"))
    b <- b[retain]
    resolved <- is.finite(b$benchmark_geometry_mcse) && is.finite(b$benchmark_root_n_variance) &&
      is.finite(b$reported_root_n_variance_limit) &&
      b$benchmark_geometry_mcse <= geometry_relative_goal * min(b$benchmark_root_n_variance, b$reported_root_n_variance_limit) &&
      !isTRUE(b$geometry_below_jensen)
    b$geometry_precision_goal <- geometry_relative_goal
    b$benchmark_comparison_status <- if (resolved) "precision_goal_met" else "geometry_unresolved_or_goal_not_met"
    b$empirical_minus_sampling_benchmark <- r$empirical_root_n_variance - b$benchmark_root_n_variance
    b$mean_reported_minus_reported_limit <- r$mean_root_n_variance - b$reported_root_n_variance_limit
    b$mean_reported_minus_sampling_benchmark <- r$mean_root_n_variance - b$benchmark_root_n_variance
    emp_se <- sqrt(r$empirical_root_n_variance_mcse_asymptotic^2 + b$benchmark_geometry_mcse^2)
    report_se <- sqrt(r$mean_root_n_variance_mcse^2 + b$benchmark_geometry_mcse^2)
    b$empirical_discrepancy_combined_mcse <- emp_se
    b$reported_discrepancy_combined_mcse <- report_se
    b$empirical_discrepancy_z <- if (resolved && is.finite(emp_se) && emp_se > 0)
      b$empirical_minus_sampling_benchmark / emp_se else NA_real_
    b$reported_limit_discrepancy_z <- if (resolved && is.finite(report_se) && report_se > 0)
      b$mean_reported_minus_reported_limit / report_se else NA_real_
    t <- result$table
    # Alpha cancels exactly even when the numeric geometric constant is missing.
    constants <- t[t$estimand == r$estimand, ]
    delta <- constants$reported_variance_constant[constants$comparison == "fitted_adjusted"] -
      constants$reported_variance_constant[constants$comparison == "fitted_naive"]
    b$adjusted_minus_naive_reported_limit <- delta
    b$quadrature_maximum_reported_absolute_error <- if (is.null(result$weight_integrals)) 0 else
      max(result$weight_integrals$absolute_errors)
    added[[i]] <- b
  }
  extra <- do.call(rbind, added); rownames(extra) <- NULL
  if (any(names(extra) %in% names(summary))) stop("summary already contains benchmark fields")
  out <- cbind(summary, extra)
  attr(out, "interpretation") <- paste("Actual sampling variance and naive reported variance limits are distinct.",
    "Geometry MCSE is propagated separately; combined discrepancy MCSE assumes an independently simulated geometric constant.",
    "Quadrature error is a separate deterministic integration diagnostic, not a propagated probabilistic bound.",
    "Z scores are descriptive finite-simulation discrepancies, not tests of the asymptotic theorem.",
    "Shared data, paired comparisons and shared geometry induce dependence across rows.")
  out
}

.wm_first_stage_join_main <- function(args) {
  if (length(args) < 2L) stop("usage: summary.csv output.csv [geometry.rds ...]")
  output <- args[2]; metadata <- paste0(output, ".metadata.rds")
  if (any(file.exists(c(output, metadata)))) stop("joined output or metadata already exists")
  script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = TRUE)
  dir <- dirname(script)
  source(file.path(dir, "first_stage.R")); source(file.path(dir, "first_stage_benchmarks.R"))
  inputs <- if (length(args) > 2L) args[-c(1, 2)] else character()
  geometry <- list()
  for (file in inputs) {
    g <- readRDS(file)
    if (inherits(g, "wm_geometry") || inherits(g, "wm_geometry_alpha")) g <- setNames(list(g), sprintf("M%d_d%d", g$M, g$d))
    if (!is.list(g) || is.null(names(g)) || anyDuplicated(names(g)) || any(names(g) %in% names(geometry)))
      stop("invalid or duplicate geometry input keys")
    geometry <- c(geometry, g)
  }
  joined <- wm_first_stage_join(utils::read.csv(args[1], stringsAsFactors = FALSE), geometry)
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(joined, output, row.names = FALSE, na = "NA")
  saveRDS(list(created = Sys.time(), session = utils::sessionInfo(),
    input_md5 = tools::md5sum(c(args[1], inputs)), output_md5 = tools::md5sum(output),
    source_md5 = tools::md5sum(c(script, file.path(dir, c("first_stage.R", "first_stage_benchmarks.R")))),
    interpretation = attr(joined, "interpretation")), metadata)
  cat("Joined", nrow(joined), "first-stage summary rows; geometric uncertainty is separate.\n")
}
if (sys.nframe() == 0L) .wm_first_stage_join_main(commandArgs(trailingOnly = TRUE))
