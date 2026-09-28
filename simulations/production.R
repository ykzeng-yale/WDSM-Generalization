#!/usr/bin/env Rscript
# Sourceable production task plan; CLI modes prepare, run, collect.
# No Slurm submission is performed by this script.

wm_production_plan <- function(config) {
  required <- c("scenario", "n", "d", "M", "design", "seed")
  if (!is.data.frame(config) || !identical(names(config), required) || nrow(config) != 48L ||
      anyNA(config) || anyDuplicated(config$scenario) || anyDuplicated(config$seed))
    stop("Production requires the frozen 48-scenario configuration with unique names/seeds")
  for (field in c("n", "d", "M", "seed")) {
    x <- config[[field]]
    if (!is.numeric(x) || is.complex(x) || any(!is.finite(x)) ||
        any(x != floor(x)) || any(x < 1 | x > .Machine$integer.max)) stop("Invalid config integer")
  }
  if (any(!config$n %in% c(200, 800, 2000)) || any(!config$d %in% c(1, 2, 3, 5)) ||
      any(!config$M %in% c(1, 3)) || any(!config$design %in% c("strong", "weak")) ||
      anyDuplicated(config[c("n", "d", "M", "design")])) stop("Configuration is not the full primary grid")
  names_expected <- with(config, sprintf("%s_d%d_m%d_n%d", design, d, M, n))
  if (!identical(config$scenario, names_expected)) stop("Scenario names must describe their exact configuration")
  timing <- which(config$n == 2000 & config$M == 3 & config$d %in% c(1, 5))
  blocks <- list()
  add <- function(i, first, last, stage) {
    blocks[[length(blocks) + 1L]] <<- data.frame(config[i, , drop = FALSE],
      first = first, last = last, replications = last - first + 1L,
      requested_records = 12L * (last - first + 1L), stage = stage,
      stringsAsFactors = FALSE, row.names = NULL)
  }
  # These observations are part of the 1,000; the timing stage never reruns them.
  for (i in timing) add(i, 1L, 10L, "timing")
  for (i in seq_len(nrow(config))) {
    start <- if (i %in% timing) 11L else 1L
    size <- if (config$n[i] == 200L) 1000L else if (config$n[i] == 800L) 250L else 100L
    for (first in seq.int(start, 1000L, by = size)) add(i, first, min(1000L, first + size - 1L), "main")
  }
  plan <- do.call(rbind, blocks)
  plan$task_id <- seq_len(nrow(plan))
  plan <- plan[c("task_id", setdiff(names(plan), "task_id"))]
  plan$result_file <- sprintf("task%04d_%s_r%04d-%04d.csv", plan$task_id,
                              plan$scenario, plan$first, plan$last)
  plan
}

.wm_production_hashes <- function(root) {
  relative <- sort(c(file.path("R", list.files(file.path(root, "R"), "\\.R$")),
    file.path("simulations", list.files(file.path(root, "simulations"), "\\.R$"))))
  stats::setNames(unname(tools::md5sum(file.path(root, relative))), relative)
}

.wm_production_main <- function(args) {
  if (length(args) < 3L || !args[1L] %in% c("prepare", "run", "collect"))
    stop("usage: prepare config.csv run_directory | run config.csv run_directory task_id [wall_seconds] | collect config.csv run_directory")
  mode <- args[1L]
  if ((mode %in% c("prepare", "collect") && length(args) != 3L) ||
      (mode == "run" && !length(args) %in% c(4L, 5L))) stop("Invalid mode arguments")
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(script_arg) != 1L) stop("Run CLI with Rscript")
  script <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
  root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
  config_path <- normalizePath(args[2L], mustWork = TRUE)
  config <- utils::read.csv(config_path, stringsAsFactors = FALSE)
  plan <- wm_production_plan(config)
  out <- args[3L]
  plan_file <- file.path(out, "production_plan.csv")
  plan_meta <- file.path(out, "production_plan.metadata.rds")
  if (mode == "prepare") {
    if (file.exists(plan_file) || file.exists(plan_meta)) stop("Production plan already exists")
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(out)) stop("Cannot create run directory")
    utils::write.csv(plan, plan_file, row.names = FALSE)
    saveRDS(list(config = config, plan = plan, created = Sys.time(),
      code_md5 = .wm_production_hashes(root), config_md5 = unname(tools::md5sum(config_path)),
      plan_md5 = unname(tools::md5sum(plan_file)), requested_replications = 48000L,
      requested_records = 576000L, tasks = nrow(plan), timing_tasks = 1:4,
      main_tasks = 5:244, session = utils::sessionInfo()), plan_meta)
    cat("Prepared 244 tasks: timing IDs1-4, then main IDs5-244; 48,000 datasets, 576,000 records.\n")
    return(invisible(plan))
  }
  out <- normalizePath(out, mustWork = TRUE)
  frozen <- readRDS(plan_meta)
  if (!identical(config, frozen$config) || !identical(plan, frozen$plan) ||
      !identical(unname(tools::md5sum(config_path)), frozen$config_md5) ||
      !identical(unname(tools::md5sum(plan_file)), frozen$plan_md5) ||
      !identical(.wm_production_hashes(root), frozen$code_md5))
    stop("Code, configuration or production plan changed since preparation")
  rscript <- file.path(R.home("bin"), "Rscript")
  results <- file.path(out, "results", plan$result_file)
  task_meta <- paste0(results, ".task.rds")
  if (mode == "collect") {
    for (i in seq_len(nrow(plan))) {
      if (!file.exists(results[i]) || !file.exists(task_meta[i])) stop("Missing production task: ", i)
      m <- readRDS(task_meta[i])
      if (!identical(m$status, "completed") || !isTRUE(m$code_unchanged) ||
          !identical(m$task, plan[i, , drop = FALSE]) ||
          !identical(m$code_md5, frozen$code_md5)) stop("Uncompleted or mismatched task: ", i)
    }
    summary <- file.path(out, "summary.csv")
    if (file.exists(summary)) stop("Summary already exists")
    log <- file.path(out, "summary.log")
    if (file.exists(log)) stop("Summary log already exists; preserve prior attempt")
    status <- system2(rscript, c("--vanilla", shQuote(file.path(root, "simulations", "summarize.R")),
      shQuote(summary), shQuote(results)), stdout = log, stderr = log, timeout = 600)
    if (status != 0L) stop("Production summarization failed; inspect retained log")
    cat("Collected all 244 completed tasks. Benchmark joining remains separate.\n")
    return(invisible(summary))
  }
  task_id <- suppressWarnings(as.numeric(args[4L]))
  wall <- if (length(args) == 5L) suppressWarnings(as.numeric(args[5L])) else 1700
  if (!is.finite(task_id) || task_id != floor(task_id) || task_id < 1 || task_id > nrow(plan) ||
      !is.finite(wall) || wall < 1 || wall > 1700) stop("Invalid task ID or bounded wall time")
  manifest <- Sys.getenv("WDSM_RUN_MANIFEST", unset = "")
  if (!nzchar(manifest) || !file.exists(manifest) || dir.exists(manifest))
    stop("Set WDSM_RUN_MANIFEST to the existing immutable snapshot manifest")
  task_id <- as.integer(task_id)
  task <- plan[task_id, , drop = FALSE]
  log <- file.path(out, "logs", sprintf("task%04d.log", task_id))
  guards <- c(results[task_id], paste0(results[task_id], ".metadata.rds"), task_meta[task_id], log)
  if (any(file.exists(guards))) stop("Task output exists; preserve attempts and do not rerun completed replications")
  for (dir in c(dirname(results[task_id]), dirname(log))) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(dir)) stop("Cannot create task output directory")
  }
  m <- list(task = task, status = "started", started = Sys.time(), wall_seconds = wall,
    code_md5 = frozen$code_md5, config_md5 = frozen$config_md5, manifest = normalizePath(manifest),
    manifest_md5 = unname(tools::md5sum(manifest)), session = utils::sessionInfo(),
    slurm = Sys.getenv(c("SLURM_JOB_ID", "SLURM_ARRAY_JOB_ID", "SLURM_ARRAY_TASK_ID", "SLURM_JOB_ACCOUNT")))
  saveRDS(m, task_meta[task_id])
  started <- proc.time()[["elapsed"]]
  status <- system2(rscript, c("--vanilla", shQuote(file.path(root, "simulations", "run.R")),
    shQuote(config_path), shQuote(task$scenario), as.character(task$first), as.character(task$last),
    shQuote(results[task_id]), shQuote(root)), stdout = log, stderr = log, timeout = wall)
  m$elapsed_seconds <- proc.time()[["elapsed"]] - started
  m$exit_status <- status
  m$code_unchanged <- identical(.wm_production_hashes(root), frozen$code_md5) &&
    identical(unname(tools::md5sum(config_path)), frozen$config_md5)
  m$status <- if (status == 0L && m$code_unchanged) "completed" else "failed"
  m$finished <- Sys.time()
  saveRDS(m, task_meta[task_id])
  if (m$status != "completed") stop("Production task failed; retain outputs and inspect ", log)
  cat("Completed task", task_id, "with", task$replications, "replications in", m$elapsed_seconds, "seconds.\n")
}

if (sys.nframe() == 0L) .wm_production_main(commandArgs(trailingOnly = TRUE))
