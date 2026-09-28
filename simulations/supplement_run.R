#!/usr/bin/env Rscript
# config.csv scenario first last output.csv [source_package]
.wm_supplement_main <- function(args) {
  if (length(args) < 5L || length(args) > 6L)
    stop("usage: config.csv scenario first last output.csv [source_package]")
  script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = TRUE)
  directory <- dirname(script)
  helpers <- file.path(directory, c("first_stage.R", "first_stage_benchmarks.R", "benchmarks.R", "supplement.R"))
  for (file in helpers) source(file)
  config <- utils::read.csv(args[1L], stringsAsFactors = FALSE)
  if (!is.data.frame(config) || !nrow(config) || anyDuplicated(config$scenario) ||
      !identical(names(config), names(wm_supplement_config()))) stop("Invalid supplemental configuration")
  cfg <- config[config$scenario == args[2L], , drop = FALSE]
  if (nrow(cfg) != 1L) stop("Select exactly one supplemental scenario")
  declared <- wm_supplement_config(cfg$replications)
  if (!isTRUE(all.equal(cfg, declared[declared$scenario == cfg$scenario, ], check.attributes = FALSE)))
    stop("Use the declared supplemental configuration")
  bounds <- suppressWarnings(as.numeric(args[3:4]))
  first <- .wm_fs_integer(bounds[1L], "first", 1L, cfg$replications)
  last <- .wm_fs_integer(bounds[2L], "last", first, cfg$replications)
  output <- args[5L]; meta_path <- paste0(output, ".metadata.rds")
  diag_path <- paste0(output, ".diagnostics.rds")
  if (any(file.exists(c(output, meta_path, diag_path)))) stop("Existing output is never overwritten")
  source_files <- character()
  if (length(args) == 6L) {
    root <- normalizePath(args[6L], mustWork = TRUE)
    source_files <- sort(list.files(file.path(root, "R"), "\\.R$", full.names = TRUE))
    if (!length(source_files)) stop("No package sources")
    for (file in source_files) source(file)
  } else library(wdsmatch)
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  records_per_rep <- if (cfg$branch == "split") 3L else 8L
  metadata <- list(schema_version = "supplement_v1", config = cfg, first = first, last = last,
    records_per_replication = records_per_rep, requested_records = (last-first+1L)*records_per_rep,
    records_written = 0L, session = utils::sessionInfo(), started = Sys.time(),
    code_md5 = tools::md5sum(c(script, helpers, source_files)),
    design_md5 = tools::md5sum(file.path(directory, "supplement_design.md")),
    config_file_md5 = tools::md5sum(args[1L]),
    immutable_sha256_manifest = Sys.getenv("WDSM_RUN_MANIFEST", unset = NA_character_),
    host = unname(Sys.info()["nodename"]),
    scope = "Fixed pre-specified configurations; all failures retained; component rows are paired, not independent datasets")
  saveRDS(metadata, meta_path); diagnostics <- list()
  for (replication in seq.int(first, last)) {
    result <- wm_supplement_run_rep(cfg, replication)
    stopifnot(nrow(result$records) == records_per_rep)
    present <- file.exists(output)
    utils::write.table(result$records, output, sep = ",", row.names = FALSE,
      col.names = !present, append = present, na = "NA", qmethod = "double")
    diagnostics[[as.character(replication)]] <- result$diagnostics
    saveRDS(diagnostics, diag_path)
    metadata$records_written <- metadata$records_written + records_per_rep
    metadata$columns <- names(result$records); metadata$last_completed_replication <- replication
    saveRDS(metadata, meta_path)
  }
  metadata$completed <- Sys.time()
  metadata$code_unchanged <- identical(metadata$code_md5, tools::md5sum(c(script, helpers, source_files)))
  saveRDS(metadata, meta_path)
  stopifnot(metadata$records_written == metadata$requested_records, metadata$code_unchanged)
  cat("Completed", metadata$records_written, "supplemental records, retaining every failure.\n")
}
if (sys.nframe() == 0L) .wm_supplement_main(commandArgs(trailingOnly = TRUE))
