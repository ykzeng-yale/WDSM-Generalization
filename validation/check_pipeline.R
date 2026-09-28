# Deterministic pipeline checks; this is not a performance/coverage experiment.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
root <- normalizePath(file.path(dirname(script), ".."))
work <- tempfile("wdsm-pipeline-")
dir.create(work)
# This script retains its small diagnostic directory if a check fails.
run <- file.path(root, "simulations", "run.R")
summarize <- file.path(root, "simulations", "summarize.R")
rscript <- file.path(R.home("bin"), "Rscript")
invoke <- function(script, args, expected = 0L, pattern = NULL) {
  log <- tempfile("command-", tmpdir = work)
  status <- system2(rscript, c("--vanilla", shQuote(script), shQuote(args)),
                     stdout = log, stderr = log, timeout = 120)
  output <- paste(readLines(log, warn = FALSE), collapse = "\n")
  if ((expected == 0L && status != 0L) || (expected != 0L && status == 0L) ||
      (!is.null(pattern) && !grepl(pattern, output))) stop(output)
  invisible(output)
}
config <- file.path(work, "config.csv")
write.csv(data.frame(scenario = "check", n = 40L, d = 2L, M = 1L,
                     design = "strong", seed = 928009L), config, row.names = FALSE)
full <- file.path(work, "full.csv")
one <- file.path(work, "one.csv")
two <- file.path(work, "two.csv")
invoke(run, c(config, "check", 1, 2, full, root))
invoke(run, c(config, "check", 1, 1, one, root))
invoke(run, c(config, "check", 2, 2, two, root))
a <- read.csv(full); b <- rbind(read.csv(one), read.csv(two))
stopifnot(nrow(a) == 24L, all(a$status == "ok"))
fields <- setdiff(names(a), "elapsed_seconds")
stopifnot(isTRUE(all.equal(a[fields], b[fields], check.attributes = FALSE)))
m <- readRDS(paste0(full, ".metadata.rds"))
stopifnot(length(m$code_md5) > 3L, length(m$config_file_md5) == 1L,
          length(m$host) == 1L, !is.null(m$completed))
summary <- file.path(work, "summary.csv")
invoke(summarize, c(summary, one, two))
s <- read.csv(summary)
stopifnot(nrow(s) == 12L, all(s$requested == 2L), all(s$failed == 0L),
          all(is.na(s$coverage[s$correction == "raw"])))
invoke(summarize, c(file.path(work, "duplicate.csv"), full, one), expected = 1,
       pattern = "duplicate scenario")
invoke(run, c(config, "check", 1.9, 2, file.path(work, "fractional.csv"), root),
       expected = 1, pattern = "invalid replication")
orphan <- file.path(work, "orphan.csv")
saveRDS(list(partial = TRUE), paste0(orphan, ".metadata.rds"))
invoke(run, c(config, "check", 1, 1, orphan, root), expected = 1,
       pattern = "metadata already exists")
# Conflicting configurations must not pool under one scenario name.
altered <- read.csv(two); altered$n <- altered$n + 1L
write.csv(altered, two, row.names = FALSE)
meta <- readRDS(paste0(two, ".metadata.rds")); meta$config$n <- meta$config$n + 1L
saveRDS(meta, paste0(two, ".metadata.rds"))
invoke(summarize, c(file.path(work, "mixed.csv"), one, two), expected = 1,
       pattern = "inconsistent scenario configuration")
# An impossible donor count produces every requested error row and no CI.
failed_config <- file.path(work, "failed-config.csv")
write.csv(data.frame(scenario = "failure", n = 20L, d = 2L, M = 50L,
                     design = "weak", seed = 928010L), failed_config, row.names = FALSE)
failed <- file.path(work, "failed.csv")
invoke(run, c(failed_config, "failure", 1, 1, failed, root))
f <- read.csv(failed)
stopifnot(nrow(f) == 12L, all(f$status == "error"), all(is.na(f$estimate)))
failed_summary <- file.path(work, "failed-summary.csv")
invoke(summarize, c(failed_summary, failed))
fs <- read.csv(failed_summary)
stopifnot(all(fs$failed == 1L), all(fs$failure_rate == 1), all(is.na(fs$bias)))
cat("PASS: indexed-stream batching, exact record counts, provenance, failure retention, duplicate/configuration guards, strict endpoints and orphan-manifest preservation.\n")
unlink(work, recursive = TRUE)
