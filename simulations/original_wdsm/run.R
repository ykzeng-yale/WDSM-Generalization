#!/usr/bin/env Rscript
# A direct sequential entry point; preflight is read-only and consumes no RNG.
args <- commandArgs(TRUE)
if (length(args) != 7L || !args[1L] %in% c("preflight", "run"))
  stop("Usage: Rscript run.R preflight|run UPSTREAM POPULATION_RDS NEW_OUTPUT_DIR GoodOverlap|PoorOverlap retrospective|prospective REPLICATES_COMMA_SEPARATED")
upstream <- normalizePath(args[2L], mustWork = TRUE)
population_path <- normalizePath(args[3L], mustWork = TRUE)
output <- args[4L]; overlap <- args[5L]; design <- args[6L]
ids <- strsplit(args[7L], ",", fixed = TRUE)[[1L]]
if (!length(ids) || any(!grepl("^[0-9]+$", ids))) stop("Invalid replicate IDs.")
replicates <- as.integer(ids)
if (anyNA(replicates) || any(!replicates %in% 1:1000) || anyDuplicated(replicates))
  stop("Replicate IDs must be distinct integers from 1 through 1000.")
if (!overlap %in% c("GoodOverlap", "PoorOverlap") ||
    !design %in% c("retrospective", "prospective")) stop("Unsupported original setting.")
if (isTRUE(file.info(population_path)$isdir)) stop("Population input must be an RDS file.")
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("Call as an Rscript file.")
directory <- dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE))
had_rng <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
if (had_rng) rng_before <- get(".Random.seed", .GlobalEnv)
source(file.path(directory, "load.R"))
output <- original_wdsm_output(output, c(file.path(directory, "..", ".."),
  upstream, dirname(population_path)))
loaded <- original_wdsm_load(directory, upstream)
if (args[1L] == "preflight") {
  if (!identical(had_rng, exists(".Random.seed", .GlobalEnv, inherits = FALSE)) ||
      (had_rng && !identical(rng_before, get(".Random.seed", .GlobalEnv))))
    stop("Preflight changed RNG state.")
  cat("Preflight passed: statistical sources, upstream, package dependencies, selected paths and IDs checked.\n")
  cat("Population contents and the run branch were not executed; no data, counts, fits, files or RNG draws created.\n")
} else {
  population <- readRDS(population_path)
  required <- c(paste0("X", 1:6), "A", "Y_0", "Y_1", "Y_obs", "Cluster", "Strata")
  if (!is.data.frame(population) || nrow(population) != 1000000L ||
      !all(required %in% names(population))) stop("Invalid original population schema.")
  truth <- c(PATE = mean(population$Y_1) - mean(population$Y_0),
    PATT = mean(population$Y_1[population$A == 1]) - mean(population$Y_0[population$A == 1]))
  population_sha <- digest::digest(file = population_path, algo = "sha256", serialize = FALSE)
  if (!dir.create(output, showWarnings = FALSE)) stop("Cannot reserve new output directory.")
  for (replicate in replicates) {
    input <- loaded$engine$ows_generate_sample(population, overlap, design, replicate,
      file.path(directory, "sampler.R"), population_sha, truth, 200L)
    loaded$engine$ows_run_sample(input, file.path(output, sprintf("replicate-%04d", replicate)),
      loaded$wm, loaded$original, loaded$upstream)
  }
}
