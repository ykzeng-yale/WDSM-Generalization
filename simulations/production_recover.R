#!/usr/bin/env Rscript
# Preserve a timed-out production attempt; run only complete-prefix missing IDs.
# Usage: plan|run|collect frozen_package production_directory recovery_directory [task_id]
# Stage this script OUTSIDE the frozen package so its source hash is unchanged.
.wm_recovery_grid <- function(r, task, allow_prefix = FALSE) {
  keys <- c("replication", "method", "estimand", "correction")
  if (!nrow(r) || !all(c(keys, "scenario", "status", "n", "d", "M", "design", "seed") %in% names(r)) ||
      anyNA(r[keys]) || anyDuplicated(r[keys])) stop("Invalid or duplicated attempted record grid")
  for (f in c("scenario", "n", "d", "M", "design", "seed"))
    if (anyNA(r[[f]]) || any(r[[f]] != task[[f]])) stop("Attempt/configuration mismatch: ", f)
  if (any(!r$status %in% c("ok", "error"))) stop("Unknown estimator status")
  ids <- sort(unique(r$replication))
  last <- if (allow_prefix) max(ids) else task$last
  if (!is.numeric(ids) || any(!is.finite(ids)) || any(ids != floor(ids)) ||
      last > task$last || !identical(as.integer(ids), seq.int(as.integer(task$first), as.integer(last))))
    stop("Attempt is not an exact contiguous prefix of requested replication IDs")
  expected <- expand.grid(replication = seq.int(task$first, last),
    method = c("self_normalized", "stabilized"), estimand = c("PATE", "PATT"),
    correction = c("oracle", "polynomial", "raw"), stringsAsFactors = FALSE)
  key <- function(x) do.call(paste, c(x[keys], sep = "\r"))
  if (nrow(r) != nrow(expected) || !setequal(key(r), key(expected)))
    stop("An incomplete replication exists; automatic prefix recovery refuses it")
  last
}

.wm_recovery_main <- function(args) {
  if (length(args) < 4L || !args[1] %in% c("plan", "run", "collect") ||
      (args[1] == "run" && length(args) != 5L) ||
      (args[1] != "run" && length(args) != 4L))
    stop("usage: plan|run|collect frozen_package production_directory recovery_directory [task_id]")
  mode <- args[1]; root <- normalizePath(args[2], mustWork = TRUE)
  production <- normalizePath(args[3], mustWork = TRUE); recovery <- args[4]
  source(file.path(root, "simulations", "production.R"))
  frozen <- readRDS(file.path(production, "production_plan.metadata.rds"))
  if (!identical(.wm_production_hashes(root), frozen$code_md5)) stop("Frozen estimator/runner source changed")
  plan <- frozen$plan
  if (!identical(unname(tools::md5sum(file.path(root, "simulations", "production_config.csv"))), frozen$config_md5))
    stop("Frozen configuration changed")
  original <- file.path(production, "results", plan$result_file)
  original_meta <- paste0(original, ".metadata.rds")
  task_meta <- paste0(original, ".task.rds")
  originals <- c(original, original_meta, task_meta)
  if (any(!file.exists(originals))) stop("Missing original production artifacts")
  audit_file <- file.path(recovery, "recovery_plan.rds")
  if (mode == "plan") {
    if (file.exists(recovery)) stop("Use a new recovery directory")
    # Audit all attempts before creating derivative files or executing work.
    audit <- vector("list", nrow(plan))
    for (i in seq_len(nrow(plan))) {
      task <- plan[i, , drop = FALSE]; tm <- readRDS(task_meta[i]); m <- readRDS(original_meta[i])
      if (!identical(tm$task, task) || !identical(tm$code_md5, frozen$code_md5) ||
          !isTRUE(tm$code_unchanged)) stop("Attempt has mismatched or changed source: ", i)
      if (!identical(m$config$scenario, task$scenario) || m$first != task$first || m$last != task$last)
        stop("Attempt metadata range mismatch: ", i)
      failed <- !identical(tm$status, "completed")
      if (failed && (!identical(tm$status, "failed") || tm$exit_status != 124L))
        stop("Recovery only supports explicitly timed-out child processes: ", i)
      r <- utils::read.csv(original[i], stringsAsFactors = FALSE)
      retained_last <- .wm_recovery_grid(r, task, allow_prefix = failed)
      if (!failed && (is.null(m$completed) || m$records_written != task$requested_records))
        stop("Completed task lacks completed record metadata: ", i)
      audit[[i]] <- data.frame(task_id = task$task_id, scenario = task$scenario,
        original_status = tm$status, original_exit_status = tm$exit_status,
        first = task$first, retained_last = retained_last, last = task$last,
        retained_records = nrow(r), missing_replications = task$last - retained_last,
        original_file = original[i], stringsAsFactors = FALSE)
    }
    table <- do.call(rbind, audit)
    dir.create(recovery, recursive = TRUE)
    dir.create(file.path(recovery, "prefix")); dir.create(file.path(recovery, "suffix"))
    recovery <- normalizePath(recovery, mustWork = TRUE)
    table$prefix_file <- ""; table$suffix_file <- ""
    for (i in which(table$original_status == "failed")) {
      prefix <- file.path(recovery, "prefix", plan$result_file[i])
      stopifnot(file.copy(original[i], prefix, overwrite = FALSE))
      m <- readRDS(original_meta[i])
      m$last <- table$retained_last[i]
      m$requested_records <- m$records_written <- table$retained_records[i]
      m$completed <- Sys.time()
      m$recovery <- list(kind = "verified_complete_prefix_of_timed_out_attempt",
        original_file = original[i], original_md5 = tools::md5sum(original[i]),
        original_metadata_md5 = tools::md5sum(original_meta[i]),
        original_task_metadata_md5 = tools::md5sum(task_meta[i]),
        original_requested_last = plan$last[i], recovered = Sys.time(),
        statement = "Every saved replication has the entire 12-record grid, including all estimator failures. No saved record was discarded or recomputed.")
      saveRDS(m, paste0(prefix, ".metadata.rds")); table$prefix_file[i] <- prefix
      if (table$missing_replications[i] > 0)
        table$suffix_file[i] <- file.path(recovery, "suffix", sprintf("task%04d_r%04d-%04d.csv",
          i, table$retained_last[i] + 1L, plan$last[i]))
    }
    script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    audit <- list(table = table, created = Sys.time(), source_md5 = tools::md5sum(script),
      original_md5 = tools::md5sum(originals), code_md5 = frozen$code_md5,
      prefix_md5 = tools::md5sum(c(table$prefix_file[nzchar(table$prefix_file)],
        paste0(table$prefix_file[nzchar(table$prefix_file)], ".metadata.rds"))),
      plan = plan, missing_total = sum(table$missing_replications),
      requested_records = frozen$requested_records, session = utils::sessionInfo())
    saveRDS(audit, file.path(recovery, "recovery_plan.rds"))
    utils::write.csv(table, file.path(recovery, "recovery_plan.csv"), row.names = FALSE)
    cat("Verified", sum(table$retained_records), "saved records; missing", audit$missing_total,
      "replications across", sum(table$missing_replications > 0), "suffix tasks.\n")
    return(invisible(audit))
  }
  recovery <- normalizePath(recovery, mustWork = TRUE); audit <- readRDS(audit_file)
  if (!identical(audit$plan, plan) || !identical(audit$code_md5, frozen$code_md5) ||
      !identical(tools::md5sum(names(audit$original_md5)), audit$original_md5) ||
      !identical(tools::md5sum(names(audit$prefix_md5)), audit$prefix_md5))
    stop("Original attempts or derived prefixes changed after recovery audit")
  script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (!identical(unname(tools::md5sum(script)), unname(audit$source_md5))) stop("Recovery source changed after audit")
  table <- audit$table; rscript <- file.path(R.home("bin"), "Rscript")
  if (mode == "run") {
    id <- suppressWarnings(as.numeric(args[5]))
    if (!is.finite(id) || id != floor(id) || !id %in% table$task_id || table$missing_replications[id] < 1)
      stop("Choose one declared missing suffix task ID")
    output <- table$suffix_file[id]; log <- paste0(output, ".log"); report <- paste0(output, ".recovery.rds")
    if (any(file.exists(c(output, paste0(output, ".metadata.rds"), log, report)))) stop("Recovery attempt already exists")
    started <- Sys.time()
    status <- system2(rscript, c("--vanilla", shQuote(file.path(root, "simulations", "run.R")),
      shQuote(file.path(root, "simulations", "production_config.csv")), shQuote(table$scenario[id]),
      table$retained_last[id] + 1L, table$last[id], shQuote(output), shQuote(root)),
      stdout = log, stderr = log, timeout = 1700)
    final <- list(task_id = id, started = started, finished = Sys.time(), exit_status = status,
      original_prefix_last = table$retained_last[id], first = table$retained_last[id] + 1L,
      last = table$last[id], code_unchanged = identical(.wm_production_hashes(root), frozen$code_md5),
      recovery_plan_md5 = tools::md5sum(audit_file), job = Sys.getenv("SLURM_JOB_ID"))
    saveRDS(final, report)
    if (status != 0L || !final$code_unchanged) stop("Recovery suffix failed; preserve it for separate audit")
    task <- plan[id, , drop = FALSE]; task$first <- final$first
    r <- utils::read.csv(output, stringsAsFactors = FALSE)
    .wm_recovery_grid(r, task)
    cat("Completed missing suffix for original task", id, ":", nrow(r), "records.\n")
    return(invisible(final))
  }
  inputs <- original[table$original_status == "completed"]
  inputs <- c(inputs, table$prefix_file[nzchar(table$prefix_file)], table$suffix_file[nzchar(table$suffix_file)])
  expected_records <- 0
  for (i in seq_len(nrow(plan))) {
    paths <- if (table$original_status[i] == "completed") original[i] else
      c(table$prefix_file[i], table$suffix_file[i][nzchar(table$suffix_file[i])])
    if (any(!file.exists(paths))) stop("Missing recovery suffix for task ", i)
    if (nzchar(table$suffix_file[i])) {
      rr <- readRDS(paste0(table$suffix_file[i], ".recovery.rds"))
      if (rr$exit_status != 0L || !rr$code_unchanged) stop("Incomplete suffix for task ", i)
    }
    r <- do.call(rbind, lapply(paths, utils::read.csv, stringsAsFactors = FALSE))
    .wm_recovery_grid(r, plan[i, , drop = FALSE]); expected_records <- expected_records + nrow(r)
  }
  stopifnot(expected_records == frozen$requested_records)
  output <- file.path(recovery, "summary.csv"); log <- file.path(recovery, "summary.log")
  if (any(file.exists(c(output, log)))) stop("Recovered summary attempt already exists")
  status <- system2(rscript, c("--vanilla", shQuote(file.path(root, "simulations", "summarize.R")),
    shQuote(output), shQuote(inputs)), stdout = log, stderr = log, timeout = 600)
  if (status != 0L) stop("Recovered summary failed; inspect log")
  saveRDS(list(completed = Sys.time(), input_md5 = tools::md5sum(inputs),
    metadata_md5 = tools::md5sum(paste0(inputs, ".metadata.rds")),
    recovery_plan_md5 = tools::md5sum(audit_file), record_count = expected_records,
    source_code_md5 = frozen$code_md5, summary_md5 = tools::md5sum(output)),
    file.path(recovery, "completed.rds"))
  cat("Reconciled", expected_records, "unique requested records; preserved all original attempts.\n")
}
if (sys.nframe() == 0L) .wm_recovery_main(commandArgs(trailingOnly = TRUE))
