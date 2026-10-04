#!/usr/bin/env Rscript
# Reconstruct summary metrics from accepted synthetic estimator records only.
# No wdsmatch installation, model fitting, matching, counts, or RNG is used.
wm_compact_schema <- function(study) {
  if (identical(study, "original_wdsm")) {
    fields <- strsplit(paste("sampling_design overlap estimand scenario replicate method n",
      "target estimate raw_estimate variance lower upper status point_status interval_status",
      "warning_count error pairing_id bootstrap_requested bootstrap_draws_available",
      "variance_method inference_scope"), " ")[[1L]]
    integers <- c("replicate", "n", "warning_count", "bootstrap_requested", "bootstrap_draws_available")
    numbers <- c("target", "estimate", "raw_estimate", "variance", "lower", "upper")
  } else if (identical(study, "lenis")) {
    fields <- strsplit(paste("population_scenario multiplier response replicate scenario estimand",
      "selected_n respondent_n method target estimate raw_estimate variance lower upper",
      "source_coverage status point_status interval_status warning_count error inference_scope pairing_id"), " ")[[1L]]
    integers <- c("population_scenario", "multiplier", "replicate", "selected_n", "respondent_n", "warning_count")
    numbers <- c("target", "estimate", "raw_estimate", "variance", "lower", "upper", "source_coverage")
  } else stop("study must be original_wdsm or lenis")
  setNames(ifelse(fields %in% integers, "integer", ifelse(fields %in% numbers, "numeric", "character")), fields)
}

wm_compact_reaggregate <- function(study, output, code_root, archive = NULL) {
  stopifnot(requireNamespace("jsonlite", quietly = TRUE), requireNamespace("digest", quietly = TRUE))
  schema <- wm_compact_schema(study)
  code_root <- normalizePath(code_root, mustWork = TRUE)
  if (is.null(archive)) archive <- file.path(code_root, "results", study, "records")
  archive <- normalizePath(archive, mustWork = TRUE)
  if (length(output) != 1L || !is.character(output) || is.na(output) || !nzchar(output) ||
      file.exists(output) || dir.exists(output)) stop("Supply a new output directory")
  parent <- normalizePath(dirname(output), mustWork = TRUE)
  output <- file.path(parent, basename(output))
  inside <- function(path, root) identical(path, root) || startsWith(path, paste0(root, .Platform$file.sep))
  if (inside(output, code_root) || inside(output, archive))
    stop("Place generated summaries outside the source checkout and record archive")
  sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
  hash_raw <- function(bytes) digest::digest(bytes, algo = "sha256", serialize = FALSE)
  manifest_path <- file.path(archive, "manifest.json")
  manifest <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)
  stopifnot(identical(manifest$format, "wm_compact_csv_archive_v1"), identical(manifest$study, study),
    identical(vapply(manifest$schema, `[[`, character(1), "name"), names(schema)),
    identical(vapply(manifest$schema, `[[`, character(1), "type"), unname(schema)),
    identical(manifest$missing_value_token, "NA"), manifest$requested_replicates_per_cell == 1000,
    manifest$max_part_decoded_bytes == 32 * 1024^2, length(manifest$parts) > 0L)
  expected_sources <- c("simulations/paired_comparison_metrics.R",
    if (study == "original_wdsm") "simulations/original_wdsm/original_survey_reporting.R" else
      c("simulations/lenis/lenis_reporting.R", "simulations/lenis/design.json"))
  stopifnot(identical(vapply(manifest$reporting_sources, `[[`, character(1), "file"), expected_sources))
  for (pin in manifest$reporting_sources)
    if (!identical(sha(file.path(code_root, pin$file)), pin$sha256)) stop("Reporting source hash mismatch: ", pin$file)
  expected_rows <- if (study == "original_wdsm") 608000L else 1800000L
  expected_cells <- if (study == "original_wdsm") 608L else 1800L
  groups <- if (study == "original_wdsm") c("sampling_design", "overlap", "estimand") else
    c("population_scenario", "multiplier", "response", "estimand")
  stopifnot(manifest$rows == expected_rows, manifest$planned_cells == expected_cells,
    identical(unlist(manifest$unique_record_key, use.names = FALSE), c(groups, "scenario", "method", "replicate")),
    identical(unlist(manifest$pairing_group, use.names = FALSE), c(groups, "replicate")))
  # Authenticate both gzip files and their decoded bytes. Reconstruct the exact
  # original CSV stream in a temporary file to verify its whole-source digest.
  temporary <- tempfile("wm-compact-source-", fileext = ".csv")
  on.exit(unlink(temporary), add = TRUE)
  reconstructed <- file(temporary, "wb")
  on.exit(try(close(reconstructed), silent = TRUE), add = TRUE)
  frames <- vector("list", length(manifest$parts))
  next_row <- 1L
  header <- NULL
  for (i in seq_along(manifest$parts)) {
    part <- manifest$parts[[i]]
    if (!identical(part$file, sprintf("records-%03d.csv.gz", i))) stop("Invalid or unordered part name")
    path <- file.path(archive, part$file)
    stopifnot(part$rows > 0, part$first_source_row == next_row,
      part$last_source_row == next_row + part$rows - 1L,
      part$decoded_bytes <= manifest$max_part_decoded_bytes,
      file.info(path)$size == part$compressed_bytes,
      identical(sha(path), part$compressed_sha256))
    connection <- gzfile(path, "rb")
    bytes <- tryCatch(readBin(connection, "raw", n = part$decoded_bytes + 1L),
      finally = close(connection))
    stopifnot(length(bytes) == part$decoded_bytes, identical(hash_raw(bytes), part$decoded_sha256),
      manifest$header_bytes >= 1L, manifest$header_bytes < length(bytes))
    current_header <- bytes[seq_len(manifest$header_bytes)]
    stopifnot(identical(hash_raw(current_header), manifest$header_sha256),
      tail(current_header, 1L) == as.raw(10L))
    if (is.null(header)) header <- current_header else stopifnot(identical(header, current_header))
    writeBin(if (i == 1L) bytes else bytes[-seq_len(manifest$header_bytes)], reconstructed)
    # Explicit classes preserve pairing hashes and distinguish empty character
    # fields from the literal missing-value token NA. Parse whole CSV records,
    # including quoted multiline fields; never split records on newline alone.
    frame <- read.csv(text = rawToChar(bytes), check.names = FALSE, colClasses = unname(schema),
      stringsAsFactors = FALSE, na.strings = "NA", row.names = NULL, fill = FALSE,
      comment.char = "", strip.white = FALSE)
    stopifnot(identical(names(frame), names(schema)), nrow(frame) == part$rows)
    frames[[i]] <- frame
    next_row <- part$last_source_row + 1L
  }
  close(reconstructed)
  stopifnot(next_row - 1L == expected_rows, file.info(temporary)$size == manifest$source_csv_bytes,
    identical(sha(temporary), manifest$source_csv_sha256))
  records <- do.call(rbind, frames)
  rm(frames, bytes, frame)
  stopifnot(nrow(records) == expected_rows, !anyNA(records$pairing_id),
    all(grepl("^[0-9a-f]{64}$", records$pairing_id)))
  for (field in c("status", "point_status", "interval_status")) {
    counts <- table(records[[field]], useNA = "ifany")
    wanted <- unlist(manifest$status_counts[[field]])
    stopifnot(setequal(names(counts), names(wanted)),
      identical(as.numeric(counts[names(wanted)]), as.numeric(wanted)))
  }
  nonempty_error <- !is.na(records$error) & nzchar(records$error)
  expected_error_counts <- manifest$error_counts
  if (length(expected_error_counts)) {
    wanted <- setNames(vapply(expected_error_counts, `[[`, numeric(1), "rows"),
      vapply(expected_error_counts, `[[`, character(1), "error"))
    counts <- table(records$error[nonempty_error])
    stopifnot(setequal(names(counts), names(wanted)),
      identical(as.numeric(counts[names(wanted)]), as.numeric(wanted)))
  } else stopifnot(!any(nonempty_error))
  environment <- new.env(parent = globalenv())
  sys.source(file.path(code_root, expected_sources[1L]), environment)
  sys.source(file.path(code_root, expected_sources[2L]), environment)
  if (study == "original_wdsm") {
    report <- environment$ows_report(records, 1000L)
  } else {
    config <- jsonlite::read_json(file.path(code_root, "simulations/lenis/design.json"), simplifyVector = TRUE)
    report <- environment$wm_paired_comparison_metrics(records,
      environment$lenis_reporting_plan(config, requested = 1000L), group_cols = groups,
      reference_scenario = "source_specification", reference_method = "Lenis_U_PS_U_OM_OW_U")
  }
  stopifnot(nrow(report$metrics) == expected_cells, nrow(report$audit) == expected_cells)
  if (!dir.create(output, recursive = FALSE, showWarnings = FALSE)) stop("Cannot reserve output directory")
  for (name in names(report)) if (is.data.frame(report[[name]]))
    write.csv(report[[name]], file.path(output, paste0(name, ".csv")), row.names = FALSE)
  failures <- records[records$status != "ok" | records$point_status != "ok" |
    records$interval_status != "ok" | nonempty_error, , drop = FALSE]
  write.csv(failures, file.path(output, "failure_records.csv"), row.names = FALSE)
  saveRDS(report, file.path(output, "report.rds"), version = 2L)
  receipt <- list(study = study, rows = nrow(records), planned_cells = expected_cells,
    requested_replicates_per_cell = 1000L, manifest_sha256 = sha(manifest_path),
    source_csv_sha256 = manifest$source_csv_sha256, source_stream_sha256_verified = TRUE,
    compressed_and_decoded_part_hashes_verified = TRUE,
    reporting_inputs_valid = report$reporting_inputs_valid, records_complete = report$records_complete,
    point_metrics_complete = report$point_metrics_complete, failed_records_retained = nrow(failures),
    no_RNG_fits_matching_or_counts = TRUE,
    scope = "Saved synthetic record aggregation; not a new simulation or model validation")
  jsonlite::write_json(receipt, file.path(output, "aggregation_receipt.json"), auto_unbox = TRUE, pretty = TRUE)
  if (!isTRUE(report$reporting_inputs_valid) || !isTRUE(report$records_complete))
    stop("Reporting inputs rejected; complete diagnostics were retained in the output directory")
  invisible(receipt)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  if (!length(args) || identical(args, "--help")) {
    cat("Usage: Rscript simulations/reaggregate_compact_archive.R original_wdsm|lenis NEW_OUTPUT_DIR [ARCHIVE_DIR]\n",
      "Requires jsonlite and digest. Reaggregates saved records only; no fitting or simulation.\n",
      "Output parent must exist and output must be new and outside this checkout/archive.\n", sep = "")
  } else {
    if (!length(args) %in% 2:3) stop("Expected study, new output directory, and optional archive directory")
    script <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    stopifnot(length(script) == 1L)
    code_root <- dirname(dirname(normalizePath(sub("^--file=", "", script), mustWork = TRUE)))
    result <- wm_compact_reaggregate(args[1L], args[2L], code_root,
      if (length(args) == 3L) args[3L] else NULL)
    cat("Records:", result$rows, "Cells:", result$planned_cells,
      "Failed records retained:", result$failed_records_retained, "\n")
  }
}
