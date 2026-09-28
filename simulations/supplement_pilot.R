#!/usr/bin/env Rscript
# new_output_directory [wall_seconds]; pilot uses production IDs 1 through 5.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 2L) stop("usage: new_output_directory [wall_seconds]")
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
source(file.path(root, "simulations", "first_stage.R"))
source(file.path(root, "simulations", "supplement.R"))
budget <- .wm_fs_integer(if (length(args) == 2L) as.numeric(args[2L]) else 600, "wall_seconds", 1L, 3600L)
output <- args[1L]
if (file.exists(output) || !dir.create(output, recursive = TRUE)) stop("Use a new output directory")
output <- normalizePath(output)
config <- wm_supplement_config(); config_path <- file.path(output, "config.csv")
utils::write.csv(config, config_path, row.names = FALSE)
started <- proc.time()[["elapsed"]]
for (i in seq_len(nrow(config))) {
  record <- file.path(output, paste0(config$scenario[i], ".csv"))
  log <- file.path(output, paste0(config$scenario[i], ".log"))
  remaining <- floor(budget - (proc.time()[["elapsed"]] - started))
  if (remaining < 1L) stop("Pilot budget exhausted; retain partial records")
  status <- system2(file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(file.path(root, "simulations", "supplement_run.R")),
      shQuote(c(config_path, config$scenario[i], 1L, 5L, record, root))),
    stdout = log, stderr = log, timeout = remaining)
  if (status != 0L) stop("Pilot child failed; retain partial records and inspect ", log)
  dat <- utils::read.csv(record, stringsAsFactors = FALSE)
  cat(config$scenario[i], nrow(dat), "records", sum(dat$status != "ok"), "failures\n")
}
cat("Pilot elapsed", proc.time()[["elapsed"]] - started, "seconds. No coverage claim from five replications.\n")
