# Read saved outputs and apply the existing paired reporting functions.
lenis_report_preflight <- function(input_path, output, export_dir) {
  setup <- lenis_export_load(export_dir, "report")
  input_path <- normalizePath(input_path, mustWork = TRUE)
  input <- jsonlite::read_json(input_path, simplifyVector = TRUE)
  stopifnot(is.character(input$runs), length(input$runs) > 0L, !anyNA(input$runs),
    all(nzchar(input$runs)), !anyDuplicated(input$runs),
    length(input$requested) == 1L, is.finite(input$requested),
    input$requested >= 2L, input$requested <= setup$plan$formal_requested_R,
    input$requested == floor(input$requested))
  runs <- vapply(input$runs, lenis_export_directory, character(1))
  stopifnot(!anyDuplicated(runs))
  for (run in runs) stopifnot(all(file.exists(file.path(run,
    c("manifest.json", "receipt.json", "records.csv")))))
  out <- lenis_export_output(output, c(input_path, runs, setup$code_root))
  list(setup = setup, runs = runs, requested = input$requested, output = out)
}

lenis_report_saved <- function(input_path, output, export_dir) {
  checked <- lenis_report_preflight(input_path, output, export_dir)
  setup <- checked$setup; out <- checked$output
  reporting <- setup$functions
  sys.source(file.path(setup$export_dir, "..", "paired_comparison_metrics.R"), reporting)
  source_names <- c("entrypoint.R", "generate.R", "run.R", "report.R",
    "lenis_comparison_engine.R", "lenis_comparison_batch.R", "lenis_reporting.R",
    "lenis_canonical_counts.R", "design.json", "external_sources.json")
  files <- setNames(file.path(setup$export_dir, source_names), source_names)
  files <- c(files, paired_comparison_metrics.R = file.path(setup$export_dir, "..", "paired_comparison_metrics.R"))
  pins <- as.list(setNames(vapply(files, lenis_export_sha, character(1)), names(files)))
  saved <- lapply(checked$runs, function(run) {
    manifest <- jsonlite::read_json(file.path(run, "manifest.json"), simplifyVector = TRUE)
    stopifnot(identical(manifest$job, jsonlite::read_json(file.path(run, "frozen", "job.json"),
      simplifyVector = TRUE)))
    reporting$lenis_saved_reporting_records(run, setup$plan, pins)
  })
  records <- do.call(rbind, lapply(saved, `[[`, "records"))
  metrics <- reporting$wm_paired_comparison_metrics(records,
    reporting$lenis_reporting_plan(setup$plan, requested = checked$requested),
    group_cols = c("population_scenario", "multiplier", "response", "estimand"),
    reference_scenario = "source_specification", reference_method = "Lenis_U_PS_U_OM_OW_U")
  stopifnot(dir.create(out, recursive = FALSE, showWarnings = FALSE))
  saveRDS(metrics, file.path(out, "reporting_result.rds"), version = 2)
  for (name in names(metrics)) if (is.data.frame(metrics[[name]]))
    write.csv(metrics[[name]], file.path(out, paste0(name, ".csv")), row.names = FALSE)
  write.csv(records, file.path(out, "records_with_pairing.csv"), row.names = FALSE)
  saveRDS(lapply(saved, function(x) x[c("input_pins", "receipt", "source")]),
    file.path(out, "saved_input_audits.rds"), version = 2)
  valid <- isTRUE(metrics$records_complete) && isTRUE(metrics$reporting_inputs_valid)
  receipt <- list(status = if (valid) "SAVED_RECORDS_AGGREGATED_REVIEW_REQUIRED" else
    "REPORTING_REJECTED_RESULTS_RETAINED", records = nrow(records), requested_R = checked$requested,
    records_complete = metrics$records_complete, reporting_inputs_valid = metrics$reporting_inputs_valid,
    point_metrics_complete = metrics$point_metrics_complete, reporting_scope = metrics$scope,
    no_RNG_or_fits = TRUE,
    inference_scope = "External comparison with supplied fitted response weights frozen in WM intervals")
  jsonlite::write_json(receipt, file.path(out, "receipt.json"), auto_unbox = TRUE, pretty = TRUE)
  if (!valid) stop("Reporting rejected; exported diagnostics and failure records retained", call. = FALSE)
  invisible(receipt)
}

if (sys.nframe() == 0L) {
  file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  stopifnot(length(file) == 1L)
  export_dir <- dirname(normalizePath(sub("^--file=", "", file), mustWork = TRUE))
  source(file.path(export_dir, "entrypoint.R"))
  args <- commandArgs(TRUE)
  stopifnot(length(args) == 3L, args[1L] %in% c("preflight", "report"))
  if (args[1L] == "preflight") {
    lenis_report_preflight(args[2L], args[3L], export_dir)
    cat("PASS preflight only: no output creation, RNG, counts, fits or metrics\n")
  } else lenis_report_saved(args[2L], args[3L], export_dir)
}
