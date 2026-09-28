#!/usr/bin/env Rscript
# family frozen_package production_directory new_summary.csv
# Collector integrity repair: software_sources, 2026-09-28. Postprocessing only.
.wm_extension_collect_hashes <- function(x) {
  if (!is.character(x) || !length(x) || is.null(names(x)) ||
      anyNA(x) || anyNA(names(x)) || any(!nzchar(names(x))) ||
      anyDuplicated(names(x)) || any(!grepl("^[[:xdigit:]]{32}$", x)))
    stop("Invalid frozen source hash mapping")
  # Production files can be retrieved on a host with a different collation.
  # Compare path-to-hash mappings, retaining every relative path and hash.
  x[order(names(x), method = "radix")]
}

.wm_extension_collect_plan <- function(family, frozen, plan_path) {
  if (!is.list(frozen) || !identical(frozen$family, family)) stop("Frozen production family differs")
  # A CSV hash alone does not bind the plan embedded in its RDS sidecar. Rebuild
  # the complete declared grid before accepting either representation.
  expected <- wm_extension_plan(family, frozen$config)
  if (!identical(frozen$plan, expected)) stop("Frozen RDS plan differs from reconstructed complete plan")
  if (!identical(frozen$plan_md5, unname(tools::md5sum(plan_path)))) stop("Frozen plan CSV hash differs")
  csv <- utils::read.csv(plan_path, stringsAsFactors = FALSE)
  if (!identical(names(csv), names(expected)) ||
      !isTRUE(all.equal(csv, expected, check.attributes = FALSE)))
    stop("Plan CSV differs from reconstructed complete plan")
  expected
}

.wm_extension_collect_import <- function(path, task, metadata) {
  if (task$stage != "imported_pilot") {
    if (isTRUE(metadata$imported)) stop("Unexpected imported-pilot task flag")
    return(invisible(NULL))
  }
  copied <- paste0(path, c("", ".metadata.rds", ".diagnostics.rds"))
  saved <- metadata$source_md5
  if (!isTRUE(metadata$imported) || !is.character(saved) || length(saved) != 3L ||
      anyNA(saved) || any(!grepl("^[[:xdigit:]]{32}$", saved)) || any(!file.exists(copied)))
    stop("Retained pilot files or original copy hashes are missing")
  current <- tools::md5sum(copied)
  if (!identical(unname(current), unname(saved))) stop("Retained pilot copy differs from its original hashes")
  invisible(current)
}

.wm_extension_collect_main <- function(args) {
  if (length(args)!=4L) stop("usage: first_stage|supplement frozen_package production_directory new_summary.csv")
  family<-match.arg(args[1],c("first_stage","supplement"))
  root<-normalizePath(args[2],mustWork=TRUE);run<-normalizePath(args[3],mustWork=TRUE);output<-args[4]
  script<-normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)))
  # Postprocessing may be updated separately; estimator source stays frozen.
  for (name in c("first_stage.R","first_stage_summarize.R","supplement.R","supplement_summarize.R","production.R","extension_production.R"))
    source(file.path(dirname(script),name))
  frozen<-readRDS(file.path(run,"plan.rds"))
  plan<-.wm_extension_collect_plan(family,frozen,file.path(run,"plan.csv"))
  if (!identical(.wm_extension_collect_hashes(frozen$code_md5),
      .wm_extension_collect_hashes(.wm_production_hashes(root)))) stop("Frozen production source differs")
  batch_hashes <- .wm_fss_hash(frozen$batch_hashes,"frozen simulation")
  if (any(file.exists(c(output,paste0(output,".metadata.rds"))))) stop("Output exists")
  inputs<-file.path(run,"results",plan$result_file)
  tasks<-paste0(inputs,".task.rds")
  for (i in seq_len(nrow(plan))) {
    task<-plan[i,,drop=FALSE];m<-readRDS(tasks[i])
    if (!identical(m$status,"completed") || !isTRUE(m$code_unchanged) || !identical(m$task,task) ||
        !identical(m$code_md5,frozen$code_md5)) stop("Task is incomplete or mismatched: ",i)
    .wm_extension_collect_import(inputs[i],task,m)
    cfg<-frozen$config[frozen$config$scenario==task$scenario,,drop=FALSE]
    .wm_extension_batch(inputs[i],task,cfg,family,batch_hashes)
  }
  records<-if (family=="first_stage") wm_first_stage_read_batches(inputs) else wm_supplement_read_batches(inputs)
  if (nrow(records)!=sum(plan$requested_records)) stop("Production record total differs")
  summary<-if (family=="first_stage") wm_first_stage_summarize(records) else wm_supplement_summarize(records)
  dir.create(dirname(output),recursive=TRUE,showWarnings=FALSE)
  utils::write.csv(summary,output,row.names=FALSE)
  saveRDS(list(completed=Sys.time(),family=family,records=nrow(records),groups=nrow(summary),
    frozen_plan_md5=tools::md5sum(file.path(run,c("plan.csv","plan.rds"))),
    task_metadata_md5=tools::md5sum(tasks),output_md5=tools::md5sum(output),
    plan_reconstructed=TRUE,retained_pilot_tasks_verified=sum(plan$stage=="imported_pilot"),
    retained_pilot_hashes_verified=TRUE,
    provenance=attr(records,"batch_provenance"),postprocessing_md5=tools::md5sum(unique(c(script,
      list.files(dirname(script),"\\.R$",full.names=TRUE)))),session=utils::sessionInfo()),paste0(output,".metadata.rds"))
  cat("Verified all",nrow(plan),"tasks and",nrow(records),"records; produced",nrow(summary),"summary groups.\n")
}
if (sys.nframe()==0L) .wm_extension_collect_main(commandArgs(trailingOnly=TRUE))
