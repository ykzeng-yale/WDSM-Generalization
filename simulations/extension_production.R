#!/usr/bin/env Rscript
# CLI: prepare|run family config.csv run_directory [pilot_directory|task_id]
# All scheduler submission remains explicit and external to this script.
wm_extension_plan <- function(family, config) {
  family <- match.arg(family, c("first_stage", "supplement"))
  if (family == "first_stage") {
    config <- wm_first_stage_config(config)
    if (nrow(config) != 39L || any(config$replications != 1000L) ||
        any(!config$n %in% c(200,800,2000)) || any(config$n != config$m) ||
        anyDuplicated(config[c("branch","n","d","M")])) stop("Expected complete 39-cell first-stage grid")
    expected <- data.frame()
    for (n in c(200L,800L,2000L)) for (branch in c("same_sample","independent_training","estimated_weights","gaussian_same_sample")) {
      dimensions <- if (branch == "same_sample") c(2L,3L,5L) else if (branch == "gaussian_same_sample") 1L else c(1L,2L,3L,5L)
      for (d in dimensions) for (M in if (branch == "gaussian_same_sample") c(1L,3L) else 3L)
        expected <- rbind(expected, data.frame(scenario = sprintf("%s_d%d_m%d_n%d",branch,d,M,n),
          branch=branch,n=n,m=n,d=d,M=M,seed=936001L+nrow(expected),replications=1000L))
    }
    if (!isTRUE(all.equal(config, expected, check.attributes=FALSE))) stop("First-stage grid or seeds changed")
    timing <- which(config$n == 2000L & config$M == 3L &
      (config$d == 5L | config$branch == "gaussian_same_sample"))
  } else {
    if (!isTRUE(all.equal(config, wm_supplement_config(), check.attributes=FALSE)))
      stop("Supplemental configuration differs from pre-specified design")
    timing <- seq_len(nrow(config))
  }
  blocks <- list()
  add <- function(i, first, last, stage) {
    block <- config[i,,drop=FALSE]
    block$first <- first; block$last <- last; block$stage <- stage
    block$requested_records <- (last-first+1L) * if (family == "first_stage") 6L else if (block$branch == "split") 3L else 8L
    blocks[[length(blocks)+1L]] <<- block
  }
  for (i in timing) add(i,1L,5L,if (family == "first_stage") "timing" else "imported_pilot")
  for (i in seq_len(nrow(config))) {
    size <- if (config$n[i] == 200L) 250L else if (config$n[i] == 800L) 100L else 50L
    for (first in seq.int(if (i %in% timing) 6L else 1L,config$replications[i],by=size))
      add(i,first,min(config$replications[i],first+size-1L),"main")
  }
  plan <- do.call(rbind,blocks); rownames(plan) <- NULL
  plan$task_id <- seq_len(nrow(plan))
  plan$result_file <- sprintf("task%04d_%s_r%04d-%04d.csv",plan$task_id,plan$scenario,plan$first,plan$last)
  plan
}

.wm_extension_batch <- function(path, task, cfg, family, code_hashes) {
  m <- readRDS(paste0(path,".metadata.rds"))
  if (is.null(m$completed) || !isTRUE(all.equal(m$config,cfg,check.attributes=FALSE)) ||
      !isTRUE(m$first == task$first) || !isTRUE(m$last == task$last) ||
      !isTRUE(m$last_completed_replication == task$last) ||
      !isTRUE(m$records_written == task$requested_records) ||
      !isTRUE(m$requested_records == task$requested_records) ||
      !identical(.wm_fss_hash(m$code_md5,"simulation"), code_hashes)) stop("Incomplete or mismatched batch metadata")
  r <- utils::read.csv(path,stringsAsFactors=FALSE,
    colClasses=if (family == "first_stage") c(training_rng_state="character",warning="character",error="character") else c(warning="character",error="character"))
  if (nrow(r) != task$requested_records || !setequal(unique(r$replication),seq.int(task$first,task$last)) ||
      any(r$scenario != cfg$scenario) || anyDuplicated(r[if (family == "first_stage")
        c("replication","estimand","comparison") else c("replication","estimand","method","correction")]))
    stop("Incomplete or duplicated batch grid")
  if (family == "first_stage") wm_first_stage_validate_records(r) else {
    per <- if (cfg$branch == "split") 3L else 8L
    if (any(table(r$replication) != per)) stop("Incomplete supplemental replication")
    for (field in c("branch","n","d","M")) if (any(r[[field]] != cfg[[field]])) stop("Supplemental configuration mismatch")
    expected <- if (cfg$branch == "split") data.frame(estimand=c("PATE","potential0","potential1"),method="stabilized",correction="oracle") else
      expand.grid(estimand=c("PATE","PATT"),method=c("self_normalized","stabilized"),correction=c("oracle","feasible"))
    key <- function(x) paste(x$estimand,x$method,x$correction)
    for (i in unique(r$replication)) if (!setequal(key(r[r$replication == i,]),key(expected))) stop("Invalid supplemental comparisons")
    if (any(!r$status %in% c("ok","error"))) stop("Invalid supplemental status")
  }
  invisible(m)
}

.wm_extension_main <- function(args) {
  if (length(args) < 4L || length(args) > 5L || !args[1] %in% c("prepare","run"))
    stop("usage: prepare|run first_stage|supplement config.csv run_directory [pilot_directory|task_id]")
  mode <- args[1]; family <- match.arg(args[2],c("first_stage","supplement"))
  script <- normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)),mustWork=TRUE)
  root <- normalizePath(file.path(dirname(script),".."))
  for (name in c("first_stage.R","first_stage_benchmarks.R","first_stage_summarize.R","supplement.R","production.R"))
    source(file.path(root,"simulations",name))
  config_path <- normalizePath(args[3],mustWork=TRUE)
  config <- utils::read.csv(config_path,stringsAsFactors=FALSE);plan <- wm_extension_plan(family,config)
  output <- args[4];plan_path <- file.path(output,"plan.csv");meta_path <- file.path(output,"plan.rds")
  runner <- file.path(root,"simulations",paste0(family,"_run.R"))
  helpers <- if (family == "first_stage") c("first_stage.R","first_stage_benchmarks.R") else
    c("first_stage.R","first_stage_benchmarks.R","benchmarks.R","supplement.R")
  # Order/path differences are immaterial for imported pilot hashes; contents are not.
  batch_hashes <- .wm_fss_hash(tools::md5sum(c(runner,file.path(root,"simulations",helpers),
    sort(list.files(file.path(root,"R"),"\\.R$",full.names=TRUE)))),"simulation")
  full_hashes <- .wm_production_hashes(root)
  if (mode == "prepare") {
    if (dir.exists(output) || file.exists(output)) stop("Use a new production output directory")
    if (family == "supplement" && length(args) != 5L) stop("Supplement requires retained pilot directory")
    dir.create(output,recursive=TRUE);dir.create(file.path(output,"results"));dir.create(file.path(output,"logs"))
    utils::write.csv(plan,plan_path,row.names=FALSE)
    frozen <- list(family=family,config=config,plan=plan,code_md5=full_hashes,batch_hashes=batch_hashes,
      config_md5=unname(tools::md5sum(config_path)),plan_md5=unname(tools::md5sum(plan_path)),
      created=Sys.time(),session=utils::sessionInfo())
    saveRDS(frozen,meta_path)
    if (family == "supplement") for (i in which(plan$stage == "imported_pilot")) {
      task <- plan[i,,drop=FALSE];cfg <- config[config$scenario == task$scenario,,drop=FALSE]
      original <- file.path(args[5],paste0(task$scenario,".csv"))
      .wm_extension_batch(original,task,cfg,family,batch_hashes)
      target <- file.path(output,"results",task$result_file)
      source_paths <- paste0(original,c("",".metadata.rds",".diagnostics.rds"))
      target_paths <- paste0(target,c("",".metadata.rds",".diagnostics.rds"))
      if (!all(file.copy(source_paths,target_paths,overwrite=FALSE)) ||
          !identical(unname(tools::md5sum(source_paths)),unname(tools::md5sum(target_paths)))) stop("Pilot copy verification failed")
      saveRDS(list(task=task,status="completed",imported=TRUE,code_unchanged=TRUE,code_md5=full_hashes,
        source_files=normalizePath(source_paths),source_md5=tools::md5sum(source_paths)),paste0(target,".task.rds"))
    }
    cat("Prepared",nrow(plan),"tasks and",sum(plan$requested_records),"records. Stages:\n");print(table(plan$stage))
    return(invisible(plan))
  }
  if (length(args) != 5L) stop("run requires task_id")
  frozen <- readRDS(meta_path)
  if (!identical(frozen$family,family) || !identical(frozen$config,config) || !identical(frozen$plan,plan) ||
      !identical(frozen$code_md5,full_hashes) || !identical(frozen$batch_hashes,batch_hashes) ||
      !identical(frozen$config_md5,unname(tools::md5sum(config_path))) ||
      !identical(frozen$plan_md5,unname(tools::md5sum(plan_path)))) stop("Frozen code, configuration or plan changed")
  id <- .wm_fs_integer(as.numeric(args[5]),"task_id",1L,nrow(plan));task <- plan[id,,drop=FALSE]
  if (task$stage == "imported_pilot") stop("Retained pilot observations must not be rerun")
  target <- file.path(output,"results",task$result_file);log <- file.path(output,"logs",sprintf("task%04d.log",id))
  if (any(file.exists(c(target,paste0(target,c(".metadata.rds",".task.rds",".diagnostics.rds")),log))))
    stop("Output exists; preserve previous attempts")
  manifest <- Sys.getenv("WDSM_RUN_MANIFEST",unset="")
  if (!nzchar(manifest) || !file.exists(manifest)) stop("Set immutable snapshot manifest")
  m <- list(task=task,status="started",code_md5=full_hashes,manifest_md5=tools::md5sum(manifest),
    started=Sys.time(),session=utils::sessionInfo(),slurm=Sys.getenv(c("SLURM_JOB_ID","SLURM_ARRAY_JOB_ID","SLURM_ARRAY_TASK_ID","SLURM_JOB_ACCOUNT")))
  saveRDS(m,paste0(target,".task.rds"));started <- proc.time()[["elapsed"]]
  status <- system2(file.path(R.home("bin"),"Rscript"),c("--vanilla",shQuote(runner),
    shQuote(c(config_path,task$scenario,task$first,task$last,target,root))),stdout=log,stderr=log,timeout=1700)
  m$elapsed_seconds <- proc.time()[["elapsed"]]-started;m$exit_status <- status
  m$code_unchanged <- identical(full_hashes,.wm_production_hashes(root)) && identical(frozen$config_md5,unname(tools::md5sum(config_path)))
  cfg <- config[config$scenario == task$scenario,,drop=FALSE]
  m$validation_error <- if (status == 0L && m$code_unchanged) tryCatch({
    .wm_extension_batch(target,task,cfg,family,batch_hashes);NULL
  },error=conditionMessage) else "Child execution or code freeze failed"
  m$status <- if (is.null(m$validation_error)) "completed" else "failed";m$finished <- Sys.time()
  saveRDS(m,paste0(target,".task.rds"))
  if (m$status != "completed") stop("Child or batch validation failed; retain partial records: ",m$validation_error)
  cat("Completed",family,"task",id,"in",m$elapsed_seconds,"seconds.\n")
}
if (sys.nframe() == 0L) .wm_extension_main(commandArgs(trailingOnly=TRUE))
