#!/usr/bin/env Rscript
# SPDX-License-Identifier: GPL-3.0-only
# Sourcing defines the CLI only. Explicit Rscript invocation runs one phase.
meps_cli <- function(args, script_path) {
  help <- paste(
    "MEPS 2009 portable analysis: run exactly one phase per invocation.",
    "Required: --prepared-dir PATH --output-dir PATH --pair white_asian|white_hispanic",
    "          --phase preflight|points|balance|counts|refits|original-summary|contribution|summary",
    "Optional: --release-dir PATH (default: repository containing this script)",
    "          --columns 1:200 or 1,2,3 (refits only; completed columns are reused)",
    "          --methods PS_PATE,DSM_PATT,... (contribution only)",
    "Recipes are fixed at M=3, B=200 and the declared cohort seeds.",
    "Keep every output/checkpoint private. See README.md before launching fitting phases.",
    sep = "\n")
  if (!length(args) || identical(args, "--help")) {
    cat(help, "\n")
    return(invisible(NULL))
  }
  if (length(args) %% 2L) stop("Arguments must be --name value pairs; use --help")
  keys <- args[seq.int(1L, length(args), by = 2L)]
  values <- args[seq.int(2L, length(args), by = 2L)]
  allowed <- c("--prepared-dir", "--output-dir", "--pair", "--phase",
               "--release-dir", "--columns", "--methods")
  if (any(!keys %in% allowed) || anyDuplicated(keys) || any(!nzchar(values)))
    stop("Unknown, repeated or empty CLI argument")
  options <- stats::setNames(as.list(values), keys)
  required <- c("--prepared-dir", "--output-dir", "--pair", "--phase")
  if (!all(required %in% names(options))) stop("Missing required arguments; use --help")
  if (!is.null(options[["--columns"]]) && options[["--phase"]] != "refits")
    stop("--columns is only available for the refits phase")
  if (!is.null(options[["--methods"]]) && options[["--phase"]] != "contribution")
    stop("--methods is only available for the contribution phase")
  script_path <- normalizePath(script_path, mustWork = TRUE)
  release <- options[["--release-dir"]]
  if (is.null(release)) release <- dirname(dirname(dirname(script_path)))
  functions <- new.env(parent = globalenv())
  sys.source(file.path(release, "applications", "meps2009", "meps_analysis.R"), envir = functions)
  columns <- seq_len(200L)
  if (!is.null(options[["--columns"]])) {
    value <- options[["--columns"]]
    if (grepl("^[0-9]+:[0-9]+$", value)) {
      ends <- as.integer(strsplit(value, ":", fixed = TRUE)[[1L]])
      if (anyNA(ends) || ends[1L] > ends[2L] || any(!ends %in% seq_len(200L)))
        stop("Invalid --columns range")
      columns <- seq.int(ends[1L], ends[2L])
    } else if (grepl("^[0-9]+(,[0-9]+)*$", value)) {
      columns <- as.integer(strsplit(value, ",", fixed = TRUE)[[1L]])
    } else stop("Use integer comma-list or ascending range for --columns")
  }
  methods <- functions$meps_method_names()
  if (!is.null(options[["--methods"]]))
    methods <- strsplit(options[["--methods"]], ",", fixed = TRUE)[[1L]]
  result <- functions$meps_run(prepared_dir = options[["--prepared-dir"]],
    output_dir = options[["--output-dir"]], pair = options[["--pair"]],
    release_dir = release, phase = options[["--phase"]], columns = columns, methods = methods)
  cat("Completed phase:", options[["--phase"]], "for", options[["--pair"]], "\n")
  invisible(result)
}

if (sys.nframe() == 0L) {
  script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(script) != 1L) stop("Cannot identify the Rscript entrypoint")
  meps_cli(commandArgs(TRUE), script)
}
