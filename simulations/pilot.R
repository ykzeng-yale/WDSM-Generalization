#!/usr/bin/env Rscript
# Sequential bounded pilot. Production arrays use run.R directly.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("usage: config.csv output_directory replications wall_seconds")
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
root <- normalizePath(file.path(dirname(script), ".."))
config_path <- normalizePath(args[1L])
config <- utils::read.csv(config_path, stringsAsFactors = FALSE)
replications <- suppressWarnings(as.numeric(args[3L])); wall <- suppressWarnings(as.numeric(args[4L]))
if (!is.finite(replications) || replications < 1 || replications != floor(replications) ||
    !is.finite(wall) || wall < 1) stop("invalid replication or wall-time limit")
out <- args[2L]
if (dir.exists(out) || file.exists(out)) stop("use a new output directory")
dir.create(out, recursive = TRUE)
out <- normalizePath(out)
rscript <- file.path(R.home("bin"), "Rscript")
started <- proc.time()[["elapsed"]]
files <- character(nrow(config))
for (i in seq_len(nrow(config))) {
  remaining <- floor(wall - (proc.time()[["elapsed"]] - started))
  if (remaining < 1) stop("pilot wall-time budget exhausted; retain partial results")
  files[i] <- file.path(out, paste0(config$scenario[i], ".csv"))
  log <- file.path(out, paste0(config$scenario[i], ".log"))
  status <- system2(rscript, c("--vanilla", shQuote(file.path(root,"simulations","run.R")),
    shQuote(config_path), shQuote(config$scenario[i]), "1", as.character(replications),
    shQuote(files[i]), shQuote(root)), stdout = log, stderr = log, timeout = remaining)
  if (status != 0L) stop("pilot scenario failed: ", config$scenario[i], "; inspect ", log)
  cat("Completed", config$scenario[i], "in", round(proc.time()[["elapsed"]]-started,2), "elapsed seconds total\n")
}
remaining <- floor(wall - (proc.time()[["elapsed"]] - started))
if (remaining < 1) stop("budget exhausted before summarization; retain results")
status <- system2(rscript, c("--vanilla", shQuote(file.path(root,"simulations","summarize.R")),
  shQuote(file.path(out,"summary.csv")), shQuote(files)), stdout = file.path(out,"summary.log"),
  stderr = file.path(out,"summary.log"), timeout = remaining)
if (status != 0L) stop("pilot summarization failed")
cat("Pilot complete. Total elapsed seconds:", proc.time()[["elapsed"]]-started, "\n")
