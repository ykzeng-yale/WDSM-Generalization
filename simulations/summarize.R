#!/usr/bin/env Rscript
# Usage: Rscript simulations/summarize.R summary.csv result1.csv [result2.csv ...]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("provide output summary and one or more result CSVs")
if (file.exists(args[1L])) stop("summary output already exists")
inputs <- args[-1L]
if (anyDuplicated(normalizePath(inputs))) stop("duplicate input file")
records <- do.call(rbind, lapply(inputs, utils::read.csv, stringsAsFactors = FALSE))
keys <- c("scenario", "replication", "method", "estimand", "correction")
if (anyDuplicated(records[keys])) stop("duplicate scenario/replication/method record")
if (any(!records$status %in% c("ok", "error"))) stop("unknown record status")
config_fields <- c("n", "d", "M", "design", "seed")
if (any(!is.finite(records$target))) stop("nonfinite target")
for (scenario in unique(records$scenario)) {
  block <- records[records$scenario == scenario, , drop = FALSE]
  for (field in config_fields) {
    if (anyNA(block[[field]]) || length(unique(block[[field]])) != 1L)
      stop("inconsistent scenario configuration: ", scenario, " / ", field)
  }
}
# Completeness is checked against each batch's manifest, including failed records.
for (file in inputs) {
  m <- readRDS(paste0(file, ".metadata.rds"))
  r <- utils::read.csv(file, stringsAsFactors = FALSE)
  expected <- expand.grid(replication = seq.int(m$first, m$last),
                          method = c("self_normalized", "stabilized"),
                          estimand = c("PATE", "PATT"),
                          correction = c("oracle", "polynomial", "raw"),
                          stringsAsFactors = FALSE)
  present <- r[c("replication", "method", "estimand", "correction")]
  for (field in config_fields) {
    if (is.null(m$config[[field]]) || length(m$config[[field]]) != 1L ||
        anyNA(r[[field]]) || any(r[[field]] != m$config[[field]]))
      stop("record/manifest configuration mismatch: ", file, " / ", field)
  }
  if (is.null(m$completed) || nrow(r) != m$requested_records ||
      nrow(merge(expected, present)) != nrow(expected) ||
      any(r$scenario != m$config$scenario)) stop("incomplete batch: ", file)
}
groups <- unique(records[c("scenario", "method", "estimand", "correction")])
mcse_mean <- function(x) if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_
rows <- vector("list", nrow(groups))
for (i in seq_len(nrow(groups))) {
  take <- rep(TRUE, nrow(records))
  for (name in names(groups)) take <- take & records[[name]] == groups[[name]][i]
  r <- records[take, , drop = FALSE]
  ok <- r$status == "ok"
  x <- r$estimate[ok]
  if (any(!is.finite(x)) || length(unique(r$target)) != 1L) stop("invalid successful result or target")
  error <- x - r$target[1L]
  if (length(unique(r$inference_status)) != 1L) stop("inconsistent inference label")
  probability_limit <- r$target[1L]
  if (r$design[1L] == "weak" && r$method[1L] == "self_normalized") {
    k <- 0:r$M[1L]
    v0 <- sum(stats::dbinom(k, r$M[1L], 1/3) * 2*k/(r$M[1L]+k))
    v1 <- sum(stats::dbinom(k, r$M[1L], 1/4) * 3*k/(r$M[1L]+2*k))
    probability_limit <- if (r$estimand[1L] == "PATE") (v1-v0)/2 else .5-v0
  }
  R <- length(x)
  variance <- if (R > 1L) stats::var(x) else NA_real_
  sd_empirical <- sqrt(variance)
  centered_square <- if (R) (x - mean(x))^2 else numeric()
  variance_mcse <- if (R > 1L) mcse_mean(centered_square) * R / (R - 1) else NA_real_
  mse <- if (R) mean(error^2) else NA_real_
  rmse <- sqrt(mse)
  v <- r$sampling_variance[ok]
  has_variance <- R > 0L && all(is.finite(v))
  mean_v <- if (has_variance) mean(v) else NA_real_
  ratio <- if (has_variance && is.finite(variance) && variance > 0) mean_v / variance else NA_real_
  ratio_if <- if (is.finite(ratio)) (v - mean_v) / variance -
    mean_v / variance^2 * (centered_square * R / (R - 1) - variance) else numeric()
  ci <- ok & is.finite(r$lower) & is.finite(r$upper)
  coverage <- if (any(ci)) mean(r$lower[ci] <= r$target[ci] & r$upper[ci] >= r$target[ci]) else NA_real_
  width <- r$upper[ci] - r$lower[ci]
  rows[[i]] <- cbind(groups[i, , drop = FALSE], data.frame(
    n = r$n[1L], d = r$d[1L], M = r$M[1L], design = r$design[1L],
    requested = nrow(r), completed = sum(ok), failed = sum(!ok),
    inference_status = r$inference_status[1L], probability_limit = probability_limit,
    bias_from_limit = if (R) mean(x-probability_limit) else NA_real_,
    bias_from_limit_mcse = mcse_mean(x-probability_limit),
    failure_rate = mean(!ok), target = r$target[1L],
    bias = if (R) mean(error) else NA_real_, bias_mcse = mcse_mean(error),
    empirical_sd = sd_empirical,
    empirical_sd_mcse_delta = if (is.finite(sd_empirical) && sd_empirical > 0) variance_mcse / (2 * sd_empirical) else NA_real_,
    empirical_variance = variance, empirical_variance_mcse_asymptotic = variance_mcse,
    rmse = rmse, rmse_mcse_delta = if (is.finite(rmse) && rmse > 0) mcse_mean(error^2) / (2 * rmse) else NA_real_,
    mean_root_n_variance = if (has_variance) mean(r$root_n_variance[ok]) else NA_real_,
    mean_root_n_variance_mcse = if (has_variance) mcse_mean(r$root_n_variance[ok]) else NA_real_,
    mean_diagonal_root_n_variance = if (has_variance) mean(r$diagonal_root_n_variance[ok]) else NA_real_,
    diagonal_minus_reported_root_n = if (has_variance)
      mean(r$diagonal_root_n_variance[ok]-r$root_n_variance[ok]) else NA_real_,
    diagonal_minus_reported_root_n_mcse = if (has_variance)
      mcse_mean(r$diagonal_root_n_variance[ok]-r$root_n_variance[ok]) else NA_real_,
    mean_estimated_variance = mean_v,
    mean_estimated_variance_mcse = if (has_variance) mcse_mean(v) else NA_real_,
    variance_ratio = ratio, variance_ratio_mcse_delta = mcse_mean(ratio_if),
    intervals = sum(ci), coverage = coverage,
    coverage_mcse = if (any(ci)) sqrt(coverage * (1 - coverage) / sum(ci)) else NA_real_,
    mean_interval_width = if (length(width)) mean(width) else NA_real_,
    mean_interval_width_mcse = mcse_mean(width),
    total_elapsed_seconds = sum(r$elapsed_seconds),
    warning_records = sum(!is.na(r$warning) & nzchar(r$warning)),
    clipping_records = sum(r$clipped > 0), stringsAsFactors = FALSE))
}
dir.create(dirname(args[1L]), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(do.call(rbind, rows), args[1L], row.names = FALSE, na = "NA")
cat("Summarized", nrow(records), "records. Failure counts are retained; numerical summaries condition on successful records.\n")
