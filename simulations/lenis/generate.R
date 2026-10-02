# Explicit generation of all 18 literal source populations; never runs on source().
lenis_generate_populations <- function(source_dir, output, export_dir) {
  setup <- lenis_export_load(export_dir, "generate")
  sources <- lenis_export_sources(source_dir, setup$export_dir)
  out <- lenis_export_output(output, c(source_dir, setup$code_root))
  source_dir <- dirname(sources[[1L]])
  sha <- lenis_export_sha
  manifest <- jsonlite::read_json(file.path(setup$export_dir, "external_sources.json"))$files
  started <- proc.time()[[3L]]
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rounding"))
  stopifnot(dir.create(out, recursive = FALSE, showWarnings = FALSE))
  started <- proc.time()[[3L]]; env <- new.env(parent = globalenv()); evaluated <- list()
  for (e in as.list(parse(file.path(source_dir, "biosts-17039-File015.R")))) {
    if (is.call(e) && identical(e[[1L]], as.name("library"))) next
    eval(e, env); evaluated[[length(evaluated) + 1L]] <- e
    if (is.call(e) && identical(e[[1L]], as.name("for")) &&
        identical(e[[2L]], as.name("scenario"))) break
  }
  generation_seconds <- proc.time()[[3L]] - started
  stopifnot(length(env$data_sets) == 18L,
    all(vapply(env$data_sets, nrow, integer(1)) == 1000000L))
  writeLines(unlist(lapply(evaluated, deparse)), file.path(out, "evaluated_population_source.R"))
  saveRDS(.Random.seed, file.path(out, "population_final_rng.rds"), version = 2)
  truth <- do.call(rbind, lapply(seq_along(env$data_sets), function(i) {
    d <- env$data_sets[[i]]
    data.frame(index = i, scenario = unique(d$Scenario), multiplier = unique(d$Multiplier),
      N = nrow(d), PATT = mean(d$y1[d$z == 1]) - mean(d$y0[d$z == 1]),
      PATE = mean(d$y1) - mean(d$y0))
  }))
  write.csv(truth, file.path(out, "population_truth.csv"), row.names = FALSE)
  checks <- list(); cache_paths <- character()
  for (scenario in 1:3) {
    indices <- scenario + 3L * (0:5)
    first <- env$data_sets[[indices[1L]]]
    common_names <- setdiff(names(first), c("y", "y1", "Multiplier"))
    Y1 <- matrix(NA_real_, nrow(first), 6L, dimnames = list(NULL, as.character(1:6)))
    for (j in 1:6) {
      d <- env$data_sets[[indices[j]]]
      common <- all(vapply(common_names, function(k) identical(first[[k]], d[[k]]), logical(1)))
      observed <- identical(d$y, d$z * d$y1 + (1 - d$z) * d$y0)
      controls <- identical(d$y[d$z == 0], first$y[first$z == 0])
      stopifnot(common, observed, controls)
      checks[[length(checks) + 1L]] <- data.frame(scenario, multiplier = j,
        common_columns_identical = common, observed_outcome_exact = observed,
        control_outcomes_identical = controls)
      Y1[, j] <- d$y1
    }
    cache <- list(common = first[common_names], y1 = Y1, original_columns = names(first),
      scenario = scenario, truth = truth[indices, ],
      representation = "Reconstruct y = z*y1 + (1-z)*y0; exact source arithmetic checked")
    path <- file.path(out, paste0("scenario_", scenario, ".rds"))
    saveRDS(cache, path, version = 2, compress = "gzip"); cache_paths <- c(cache_paths, path)
    cat("CACHED source scenario", scenario, "six exact multiplier populations\n"); flush.console()
  }
  write.csv(do.call(rbind, checks), file.path(out, "multiplier_identity_checks.csv"), row.names = FALSE)
  receipt <- list(status = "PASS", source_specification_not_historical_numeric_publication = TRUE,
    source_population_seed = 1357L, population_scenarios = 3L, multipliers = 6L,
    N_per_population = 1000000L, generation_seconds = generation_seconds,
    RNGkind = RNGkind(), source_hashes = manifest,
    cache_hashes = as.list(setNames(vapply(cache_paths, sha, character(1)), basename(cache_paths))),
    scope = "Literal source population generation and multiplier identity checks only",
    limitations = c("Original printed Table 1 population identity remains unresolved",
      "Historical replicate list and package versions are unavailable"),
    session = capture.output(sessionInfo()))
  jsonlite::write_json(receipt, file.path(out, "receipt.json"), auto_unbox = TRUE, pretty = TRUE, digits = 16)
  invisible(receipt)
}

if (sys.nframe() == 0L) {
  file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  stopifnot(length(file) == 1L)
  export_dir <- dirname(normalizePath(sub("^--file=", "", file), mustWork = TRUE))
  source(file.path(export_dir, "entrypoint.R"))
  args <- commandArgs(TRUE)
  stopifnot(length(args) == 3L, args[1L] %in% c("preflight", "generate"))
  if (args[1L] == "preflight") {
    setup <- lenis_export_load(export_dir, "generate")
    lenis_export_sources(args[2L], export_dir)
    lenis_export_output(args[3L], c(args[2L], setup$code_root))
    cat("PASS preflight only: no output creation, RNG, populations, samples, counts or fits\n")
  } else lenis_generate_populations(args[2L], args[3L], export_dir)
}
