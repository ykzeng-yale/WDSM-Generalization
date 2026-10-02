#!/usr/bin/env Rscript
# Report supplied estimator-record CSV files without simulation or fitting.
args <- commandArgs(TRUE)
if (length(args) < 2L) stop("Usage: Rscript report.R NEW_OUTPUT_DIR RECORDS_CSV [RECORDS_CSV ...]")
output <- args[1L]
inputs <- normalizePath(args[-1L], mustWork = TRUE)
if (anyDuplicated(inputs)) stop("Duplicate record inputs.")
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("Call as an Rscript file.")
directory <- dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE))
source(file.path(directory, "load.R"))
output <- original_wdsm_output(output, c(file.path(directory, "..", ".."),
  unique(dirname(inputs))))
source(file.path(directory, "..", "paired_comparison_metrics.R"))
source(file.path(directory, "original_survey_reporting.R"))
records <- do.call(rbind, lapply(inputs, read.csv, stringsAsFactors = FALSE,
                               check.names = FALSE))
report <- ows_report(records, 1000L)
if (!dir.create(output, showWarnings = FALSE)) stop("Cannot reserve output directory.")
for (name in names(report)) if (is.data.frame(report[[name]]))
  write.csv(report[[name]], file.path(output, paste0(name, ".csv")), row.names = FALSE)
saveRDS(report, file.path(output, "report.rds"), version = 2L)
cat("Reporting inputs valid:", report$reporting_inputs_valid,
    "Complete planned records:", report$records_complete, "\n")
cat("Reporting/accounting is not statistical acceptance; partial/missing/failed rows remain visible.\n")
