#!/usr/bin/env Rscript
# Explicit generation action only; preflight parses/hashes without evaluation.
args <- commandArgs(TRUE)
if (length(args) != 2L || !args[1L] %in% c("preflight", "generate"))
  stop("Usage: Rscript generate_populations.R preflight|generate NEW_POPULATION_DIR")
output <- args[2L]
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("Call as an Rscript file.")
directory <- dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE))
source(file.path(directory, "load.R"))
output <- original_wdsm_output(output, file.path(directory, "..", ".."))
for (package in c("digest", "jsonlite"))
  if (!requireNamespace(package, quietly = TRUE)) stop("Missing package: ", package)
provenance <- jsonlite::read_json(file.path(directory, "provenance.json"))
for (entry in provenance$statistical_modules) {
  if (!identical(digest::digest(file = file.path(directory, entry$file),
                 algo = "sha256", serialize = FALSE), entry$sha256))
    stop("Transferred statistical source differs: ", entry$file)
}
for (overlap in c("GoodOverlap", "PoorOverlap"))
  invisible(parse(file.path(directory, "generators", paste0(overlap, ".R"))))
if (args[1L] == "preflight") {
  cat("Generator sources parsed and matched provenance; generation/output creation was not executed.\n")
} else {
  if (!dir.create(output, showWarnings = FALSE)) stop("Cannot reserve output directory.")
  for (overlap in c("GoodOverlap", "PoorOverlap")) {
    env <- new.env(parent = globalenv())
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    sys.source(file.path(directory, "generators", paste0(overlap, ".R")), env)
    raw <- env$data_sets[[1L]]
    # Literal column reconstruction used in the original upstream driver.
    population <- data.frame(X1 = raw$x1, X2 = raw$x2, X3 = raw$x3, X4 = raw$x4,
      X5 = raw$x5, X6 = raw$x6, A = raw$z, Y_0 = raw$y0, Y_1 = raw$y1, Y_obs = raw$y,
      Cluster = raw$Cluster, Strata = raw$Strata, S = rep(1L, nrow(raw)))
    stopifnot(nrow(population) == 1000000L, length(unique(population$Strata)) == 10L,
      all(table(population$Strata) == 100000L))
    saveRDS(population, file.path(output, paste0(overlap, "_population.rds")),
      version = 2L, compress = FALSE)
  }
}
