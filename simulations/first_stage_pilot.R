#!/usr/bin/env Rscript
# Bounded, pre-specified first-stage pilot: new_output_directory [wall_seconds]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 2L) stop("usage: new_output_directory [wall_seconds]")
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
source(file.path(root, "simulations", "first_stage.R"))
budget <- if (length(args) == 2L) suppressWarnings(as.numeric(args[2L])) else 600
budget <- .wm_fs_integer(budget, "wall_seconds", 1L, 3600L)
output <- args[1L]
if (file.exists(output)) stop("Use a new output directory")
if (!dir.create(output, recursive = TRUE)) stop("Cannot create output directory")
output <- normalizePath(output, mustWork = TRUE)
config <- wm_first_stage_pilot()
config_path <- file.path(output, "config.csv")
utils::write.csv(config, config_path, row.names = FALSE)
started <- proc.time()[["elapsed"]]
rscript <- file.path(R.home("bin"), "Rscript")
run <- function(script, args, log) {
  remaining <- floor(budget - (proc.time()[["elapsed"]] - started))
  if (remaining < 1L) stop("Pilot wall-time budget exhausted; retain partial records")
  status <- system2(rscript, c("--vanilla", shQuote(script), shQuote(args)),
                   stdout = log, stderr = log, timeout = remaining)
  if (status != 0L) stop("Pilot step failed; inspect ", log, " and retain partial records")
}
run(file.path(root, "validation", "check_first_stage.R"), character(),
    file.path(output, "validation.log"))
for (i in seq_len(nrow(config))) {
  record <- file.path(output, paste0(config$scenario[i], ".csv"))
  run(file.path(root, "simulations", "first_stage_run.R"),
      c(config_path, config$scenario[i], 1L, config$replications[i], record, root),
      file.path(output, paste0(config$scenario[i], ".log")))
  dat <- utils::read.csv(record, stringsAsFactors = FALSE)
  meta <- readRDS(paste0(record, ".metadata.rds"))
  stopifnot(nrow(dat) == config$replications[i] * 6L,
            meta$records_written == meta$requested_records,
            !is.null(meta$completed))
  cat(config$scenario[i], ":", nrow(dat), "records;", sum(dat$status != "ok"), "failures\n")
}
cat("Pilot completed in", proc.time()[["elapsed"]] - started,
    "seconds. Forty replications do not certify interval coverage.\n")
