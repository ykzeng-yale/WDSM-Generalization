# Path/dependency setup only. Sourcing this file creates no datasets or counts.
lenis_export_sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)

lenis_export_directory <- function(path) {
  stopifnot(is.character(path), length(path) == 1L, !is.na(path), nzchar(path), dir.exists(path))
  normalizePath(path, mustWork = TRUE)
}

lenis_export_output <- function(path, inputs) {
  stopifnot(is.character(path), length(path) == 1L, !is.na(path), nzchar(path), !file.exists(path))
  link <- Sys.readlink(path)
  stopifnot(is.na(link) || !nzchar(link))
  parent <- normalizePath(dirname(path), mustWork = TRUE)
  stopifnot(dir.exists(parent), file.access(parent, 2L) == 0L,
    !basename(path) %in% c(".", "..", ""))
  output <- file.path(parent, basename(path))
  for (input in inputs) {
    input <- normalizePath(input, mustWork = TRUE)
    stopifnot(output != input, !startsWith(output, paste0(input, "/")))
    if (dir.exists(input)) stopifnot(!startsWith(input, paste0(output, "/")))
  }
  output
}

lenis_export_load <- function(export_dir, mode = c("run", "generate", "report")) {
  mode <- match.arg(mode); export_dir <- lenis_export_directory(export_dir)
  packages <- switch(mode, run = c("digest", "jsonlite", "sampling", "survey", "MatchIt",
    "data.table", "sandwich", "quickmatch", "wdsmatch"),
    generate = c("digest", "jsonlite", "sampling", "data.table"), report = c("digest", "jsonlite"))
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) stop("Missing dependency: ", package, call. = FALSE)
  }
  functions <- new.env(parent = globalenv())
  for (name in c("lenis_comparison_engine.R", "lenis_comparison_batch.R",
                 "lenis_reporting.R", "lenis_canonical_counts.R"))
    sys.source(file.path(export_dir, name), functions)
  # Literal source sampling expressions resolve these names from globalenv().
  if (mode != "report") suppressPackageStartupMessages({
    library(sampling); library(data.table)
  })
  wm <- if (mode == "run") asNamespace("wdsmatch") else NULL
  if (mode == "run") stopifnot(all(vapply(c("wm_match", "wm_wdsm_fit",
    ".wm_wdsm_fit_stack", ".wm_wdsm_prediction_stack", ".wm_wdsm_controls"),
    exists, logical(1), envir = wm, inherits = FALSE)))
  plan <- jsonlite::read_json(file.path(export_dir, "design.json"), simplifyVector = TRUE)
  stopifnot(identical(as.integer(plan$population_scenarios), 1:3),
    identical(as.integer(plan$effect_multipliers), 1:6),
    identical(plan$response_models, c("No", "MAR", "MARX", "MART")),
    identical(as.integer(plan$M), c(1L, 3L, 5L)), plan$B == 200L,
    plan$population_N == 1000000L, plan$selected_n == 5000L, plan$formal_requested_R == 1000L)
  list(export_dir = export_dir, code_root = normalizePath(file.path(export_dir, "..", ".."), mustWork = TRUE),
    functions = functions, wm = wm, plan = plan,
    packages = setNames(lapply(packages, function(package) list(
      version = as.character(utils::packageVersion(package)),
      path = getNamespaceInfo(asNamespace(package), "path"))), packages))
}

lenis_export_sources <- function(source_dir, export_dir) {
  source_dir <- lenis_export_directory(source_dir)
  pins <- jsonlite::read_json(file.path(export_dir, "external_sources.json"), simplifyVector = TRUE)
  paths <- setNames(file.path(source_dir, names(pins$files)), names(pins$files))
  for (name in names(paths)) {
    stopifnot(file.exists(paths[[name]]), lenis_export_sha(paths[[name]]) == pins$files[[name]]$sha256)
  }
  # Parse without evaluation: in particular, never source File013's setwd/load/loop.
  src <- parse(paths[["biosts-17039-File013.R"]])
  helper <- new.env(parent = baseenv())
  sys.source(file.path(export_dir, "lenis_comparison_batch.R"), helper)
  for (name in c("size1", "size2", "size3", "s", "sample"))
    stopifnot(!is.null(helper$lenis_find_assignment(src, name)))
  population <- parse(paths[["biosts-17039-File015.R"]])
  stopifnot(any(vapply(as.list(population), function(e) is.call(e) &&
    identical(e[[1L]], as.name("for")) && identical(e[[2L]], as.name("scenario")), logical(1))))
  paths
}

lenis_export_cache <- function(population_dir, scenarios, population_source) {
  population_dir <- lenis_export_directory(population_dir)
  receipt <- jsonlite::read_json(file.path(population_dir, "receipt.json"))
  paths <- file.path(population_dir, paste0("scenario_", scenarios, ".rds"))
  stopifnot(receipt$status == "PASS", receipt$source_population_seed == 1357L,
    receipt$population_scenarios == 3L, receipt$multipliers == 6L, receipt$N_per_population == 1000000L,
    receipt$source_hashes[["biosts-17039-File015.R"]]$sha256 == lenis_export_sha(population_source))
  for (path in paths) stopifnot(file.exists(path),
    lenis_export_sha(path) == receipt$cache_hashes[[basename(path)]])
  list(paths = paths, receipt = receipt)
}

lenis_export_job <- function(path, plan) {
  path <- normalizePath(path, mustWork = TRUE)
  job <- jsonlite::read_json(path, simplifyVector = TRUE)
  stopifnot(job$kind == "formal", length(job$replicates) > 0L, !anyDuplicated(job$replicates),
    all(is.finite(job$replicates)), all(job$replicates == floor(job$replicates)),
    all(job$replicates %in% seq_len(plan$formal_requested_R)),
    length(job$population_scenarios) > 0L, !anyDuplicated(job$population_scenarios),
    all(job$population_scenarios %in% plan$population_scenarios),
    length(job$max_seconds) == 1L, is.finite(job$max_seconds), job$max_seconds > 0, job$max_seconds <= 7200)
  job
}

lenis_export_cli_directory <- function() {
  file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  stopifnot(length(file) == 1L)
  dirname(normalizePath(sub("^--file=", "", file), mustWork = TRUE))
}
