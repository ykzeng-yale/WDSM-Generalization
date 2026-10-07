#!/usr/bin/env Rscript
# Rscript first_stage_run.R config.csv scenario first last output.csv [source_package]
# Rscript first_stage_run.R --write-pilot new_config.csv
.wm_first_stage_main <- function(args) {
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(script_arg) != 1L) stop("Run first_stage_run.R with Rscript")
  script <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
  directory <- dirname(script)
  helper <- file.path(directory, "first_stage.R")
  source(helper)
  if (length(args) == 2L && identical(args[1L], "--write-pilot")) {
    if (file.exists(args[2L])) stop("pilot configuration already exists")
    dir.create(dirname(args[2L]), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(wm_first_stage_pilot(), args[2L], row.names = FALSE)
    cat("Wrote the pre-specified three-scenario, 40-replication pilot configuration.\n")
    return(invisible(NULL))
  }
  if (length(args) < 5L || length(args) > 6L)
    stop("usage: config.csv scenario first last output.csv [source_package], or --write-pilot new_config.csv")
  config <- wm_first_stage_config(utils::read.csv(args[1L], stringsAsFactors = FALSE))
  cfg <- config[config$scenario == args[2L], , drop = FALSE]
  if (nrow(cfg) != 1L) stop("select exactly one known first-stage scenario")
  bounds <- suppressWarnings(as.numeric(args[3:4]))
  first <- .wm_fs_integer(bounds[1L], "first replication", 1L, cfg$replications)
  last <- .wm_fs_integer(bounds[2L], "last replication", first, cfg$replications)
  output <- args[5L]
  metadata_path <- paste0(output, ".metadata.rds")
  diagnostics_path <- paste0(output, ".diagnostics.rds")
  if (any(file.exists(c(output, metadata_path, diagnostics_path))))
    stop("output or companion metadata/diagnostics already exists; partial runs are not overwritten")
  source_files <- character()
  if (length(args) == 6L) {
    source_root <- normalizePath(args[6L], mustWork = TRUE)
    source_files <- sort(list.files(file.path(source_root, "R"), "\\.R$", full.names = TRUE))
    if (!length(source_files)) stop("source_package contains no R source files")
    for (file in source_files) source(file)
  } else {
    library(WeightedMatching)
    source_root <- NULL
  }
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  stream <- wm_first_stage_stream(cfg$seed, first)
  metadata <- list(schema_version = "first_stage_v1", config = cfg,
    first = first, last = last, requested_records = (last - first + 1L) * 6L,
    records_written = 0L, session = utils::sessionInfo(), started = Sys.time(),
    seed_method = paste("L'Ecuyer-CMRG/Inversion/Rejection; replication r is",
      "nextRNGStream^(r-1)(base stream). Independent training uses nextRNGSubStream",
      "of its replication stream; complete states are recorded in each CSV row."),
    pairing = "PATE/PATT x known_first_stage/fitted_naive/fitted_adjusted; fitted variants share point and graph",
    output_schema = paste("One row per requested estimand/comparison, including failures.",
      "n is evaluation size; m is independent training size only in that branch.",
      "parameter1 is theta or logistic intercept; parameter2 is logistic slope.",
      "first_stage_seconds and generation_seconds repeat across paired rows;",
      "elapsed_seconds is comparison-specific and must not be confused with total replication time."),
    truth = list(score = list(theta = 0, tau = 1, r = 0.6, sigma2 = 1 / 16,
      population_D2 = 1 / 30, parameter_bound = 0.9, accepted_bound = 0.8),
      weights = list(parameter = c(intercept = 0, t = 3), pi = 0.5,
        tau = 1.5, sigma2 = 1 / 16, parameter_bound = 5)),
    config_file_md5 = tools::md5sum(args[1L]),
    code_md5 = tools::md5sum(c(script, helper,
      file.path(directory, "first_stage_benchmarks.R"), source_files)),
    design_md5 = tools::md5sum(file.path(directory, "first_stage_design.md")),
    source_root = source_root,
    package_description = if (length(args) == 5L) utils::packageDescription("WeightedMatching") else NULL,
    immutable_sha256_manifest = Sys.getenv("WDSM_RUN_MANIFEST", unset = NA_character_),
    host = unname(Sys.info()["nodename"]),
    allocation = Sys.getenv(c("SLURM_JOB_ID", "SLURM_ARRAY_JOB_ID", "SLURM_ARRAY_TASK_ID",
      "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE", "OMP_NUM_THREADS",
      "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"), unset = NA_character_))
  saveRDS(metadata, metadata_path)
  diagnostics <- list()
  for (replication in seq.int(first, last)) {
    result <- wm_first_stage_run_rep(cfg, replication, stream = stream)
    stopifnot(nrow(result$records) == 6L)
    present <- file.exists(output)
    utils::write.table(result$records, output, sep = ",", row.names = FALSE,
      col.names = !present, append = present, na = "NA", qmethod = "double")
    diagnostics[[as.character(replication)]] <- result$diagnostics
    saveRDS(diagnostics, diagnostics_path)
    metadata$records_written <- metadata$records_written + nrow(result$records)
    metadata$columns <- names(result$records)
    metadata$last_completed_replication <- replication
    saveRDS(metadata, metadata_path)
    stream <- parallel::nextRNGStream(stream)
  }
  stopifnot(metadata$records_written == metadata$requested_records)
  metadata$completed <- Sys.time()
  saveRDS(metadata, metadata_path)
  cat("Completed", metadata$records_written, "first-stage records, retaining every fit failure.\n")
  invisible(metadata)
}

if (sys.nframe() == 0L) .wm_first_stage_main(commandArgs(trailingOnly = TRUE))
