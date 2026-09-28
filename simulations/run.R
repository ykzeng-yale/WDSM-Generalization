#!/usr/bin/env Rscript
# Usage: Rscript simulations/run.R config.csv scenario first last output.csv [source_package]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5L || length(args) > 6L)
  stop("usage: config.csv scenario first last output.csv [source_package]")
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "dgp.R"))
if (length(args) == 6L) {
  source_files <- sort(list.files(file.path(args[6L], "R"), "\\.R$", full.names = TRUE))
  for (file in source_files) source(file)
} else library(wdsmatch)
config <- utils::read.csv(args[1L], stringsAsFactors = FALSE)
required <- c("scenario", "n", "d", "M", "design", "seed")
if (!all(required %in% names(config)) || anyDuplicated(config$scenario)) stop("invalid configuration")
cfg <- config[config$scenario == args[2L], , drop = FALSE]
if (nrow(cfg) != 1L) stop("select exactly one known scenario")
range_values <- suppressWarnings(as.numeric(args[3:4]))
if (any(!is.finite(range_values)) || any(range_values != floor(range_values)) ||
    range_values[1L] < 1L || range_values[2L] < range_values[1L] || range_values[2L] > 100000L)
  stop("invalid replication range")
first <- as.integer(range_values[1L]); last <- as.integer(range_values[2L])
output <- args[5L]
if (file.exists(output) || file.exists(paste0(output, ".metadata.rds")))
  stop("output or metadata already exists; do not overwrite completed or partial results")
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
stream <- wm_sim_stream(cfg$seed, first)
metadata <- list(config = cfg, first = first, last = last,
                 session = utils::sessionInfo(), started = Sys.time(),
                 seed_method = "L'Ecuyer-CMRG; replication r is nextRNGStream^(r-1)(base stream)",
                 requested_records = (last - first + 1L) * 12L,
                 config_file_md5 = tools::md5sum(args[1L]),
                 code_md5 = tools::md5sum(c(list.files(script_dir, "\\.R$", full.names = TRUE),
                   if (length(args) == 6L) source_files else character())),
                 package_description = if (length(args) == 5L) utils::packageDescription("wdsmatch") else NULL,
                 immutable_sha256_manifest = Sys.getenv("WDSM_RUN_MANIFEST", unset = NA_character_),
                 host = unname(Sys.info()["nodename"]),
                 allocation = Sys.getenv(c("SLURM_JOB_ID", "SLURM_ARRAY_JOB_ID", "SLURM_ARRAY_TASK_ID",
                   "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                   "VECLIB_MAXIMUM_THREADS"), unset = NA_character_))
saveRDS(metadata, paste0(output, ".metadata.rds"))
append_records <- function(records) {
  present <- file.exists(output)
  utils::write.table(do.call(rbind, records), output, sep = ",", row.names = FALSE,
                     col.names = !present, append = present, na = "NA", qmethod = "double")
}
for (replication in seq.int(first, last)) {
  assign(".Random.seed", stream, .GlobalEnv)
  dat <- wm_simulate(cfg$n, cfg$d, cfg$design)
  stream <- parallel::nextRNGStream(stream)
  records <- list()
  for (estimand in c("PATE", "PATT")) for (method in c("self_normalized", "stabilized")) {
    for (correction in c("oracle", "polynomial")) {
      captured_warnings <- character()
      started <- proc.time()[["elapsed"]]
      fit <- tryCatch(withCallingHandlers({
        if (correction == "oracle") {
          a <- list(Y = dat$Y, Z = dat$Z, weights = dat$weights, scores0 = dat$scores,
                    M = cfg$M, estimand = estimand, method = method, mean0 = dat$mean0)
          if (estimand == "PATE") a$mean1 <- dat$mean1
          if (method == "stabilized") {
            a$rho0 <- dat$rho0
            if (estimand == "PATE") a$rho1 <- dat$rho1
          }
          do.call(wm_match, a)
        } else {
          a <- list(Y = dat$Y, Z = dat$Z, weights = dat$weights, scores0 = dat$scores,
                    M = cfg$M, estimand = estimand, method = method, degree = dat$degree)
          if (method == "stabilized") a$rho_bounds <- c(.25, 4)
          do.call(wm_fit, a)
        }
      }, warning = function(w) {captured_warnings <<- c(captured_warnings, conditionMessage(w)); invokeRestart("muffleWarning")}),
      error = function(e) e)
      elapsed <- proc.time()[["elapsed"]] - started
      ok <- inherits(fit, "wm_match")
      estimate <- if (ok) fit$estimate else NA_real_
      se <- if (ok) fit$se else NA_real_
      record <- data.frame(scenario = cfg$scenario, replication = replication,
        n = cfg$n, d = cfg$d, M = cfg$M, design = cfg$design, seed = cfg$seed,
        method = method, estimand = estimand, correction = correction,
        target = unname(dat$target[estimand]), status = if (ok) "ok" else "error",
        inference_status = if (cfg$design == "weak" && method == "self_normalized")
          "assumption_failure_diagnostic" else "within_declared_theorem_scope",
        estimate = estimate, sampling_variance = if (ok) fit$variance else NA_real_,
        root_n_variance = if (ok) fit$root_n_variance else NA_real_,
        diagonal_root_n_variance = if (ok) mean(fit$contributions$actual^2) else NA_real_,
        lower = estimate - qnorm(.975) * se, upper = estimate + qnorm(.975) * se,
        elapsed_seconds = elapsed,
        clipped = if (ok && !is.null(fit$nuisance)) fit$nuisance$rho_clipped_count else 0L,
        warning = paste(captured_warnings, collapse = " | "),
        error = if (ok) "" else conditionMessage(fit), stringsAsFactors = FALSE)
      records[[length(records) + 1L]] <- record
      if (correction == "oracle") {
        raw <- record; raw$correction <- "raw"; raw$inference_status <- "point_estimate_only"
        raw$estimate <- if (ok) fit$raw_estimate else NA_real_
        raw$sampling_variance <- raw$root_n_variance <- raw$diagonal_root_n_variance <- NA_real_
        raw$lower <- raw$upper <- NA_real_
        raw$elapsed_seconds <- 0
        records[[length(records) + 1L]] <- raw
      }
    }
  }
  append_records(records)
}
metadata$completed <- Sys.time()
metadata$records_written <- nrow(utils::read.csv(output))
stopifnot(metadata$records_written == metadata$requested_records)
saveRDS(metadata, paste0(output, ".metadata.rds"))
cat("Completed", metadata$records_written, "records, including every estimator failure.\n")
