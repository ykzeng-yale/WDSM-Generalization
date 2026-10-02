#!/usr/bin/env Rscript
# Exclusive, resumeless execution: streams/data/counts are persisted before fits.
# Run only through an approved resource-owning supervisor. This script never
# launches workers, retries a dataset/count, or interprets calibration results.

wm_cal_runner <- function(args = commandArgs(trailingOnly = TRUE)) {
  usage <- paste("Usage: run_bounded_wdsm_calibration.R --config FILE --output-root EXISTING_DIRECTORY",
    "[--shard-index I --shard-count S] [--preflight-only]")
  options <- list(config=NULL,output_root=NULL,shard_index=1L,shard_count=1L,preflight_only=FALSE)
  seen <- character(); position <- 1L
  while (position <= length(args)) {
    flag <- args[position]
    if (!flag %in% c("--config","--output-root","--shard-index","--shard-count","--preflight-only") ||
        flag %in% seen) stop(usage,"; unknown or duplicate option: ",flag)
    seen <- c(seen,flag)
    key <- gsub("-","_",substring(flag,3L),fixed=TRUE)
    if (flag == "--preflight-only") {
      options[[key]] <- TRUE; position <- position+1L; next
    }
    if (position == length(args) || startsWith(args[position+1L],"--")) stop(usage,"; missing option value")
    value <- args[position+1L]
    if (flag %in% c("--shard-index","--shard-count")) {
      if (!grepl("^[1-9][0-9]*$",value) || !is.finite(as.numeric(value)) ||
          as.numeric(value) > .Machine$integer.max) stop("Shard index/count must be positive decimal integers.")
      value <- as.integer(value)
    }
    options[[key]] <- value; position <- position+2L
  }
  if (is.null(options$config) || is.null(options$output_root)) stop(usage)
  if (options$shard_index > options$shard_count) stop("Shard index exceeds shard count.")
  preflight_only <- options$preflight_only
  shard_index <- options$shard_index; shard_count <- options$shard_count
  config_path <- normalizePath(options$config, mustWork = TRUE)
  output_root <- normalizePath(options$output_root, mustWork = TRUE)
  if (!dir.exists(output_root)) stop("Output root must be an existing supervisor-owned directory.")
  if (!requireNamespace("jsonlite", quietly = TRUE) ||
      !requireNamespace("digest", quietly = TRUE) ||
      !requireNamespace("wdsmatch", quietly = TRUE)) stop("Missing pinned runtime dependencies.")
  sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
  cfg <- jsonlite::read_json(config_path, simplifyVector = TRUE)
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) != 1L) stop("Runner must be called as an Rscript file.")
  script <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
  code_root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
  required <- c("bounded_wdsm_design_reference.R", "wdsm_nuisance_stack_reference.R",
    "fitted_wm_components.R", "wdsm_fitted_pipeline_reference.R",
    "wdsm_calibration_prediction_stack.R", "bounded_wdsm_calibration.R",
    "run_bounded_wdsm_calibration.R")
  expected_names <- paste0("simulations/", required)
  if (!identical(names(cfg$source_sha256), expected_names)) stop("Source pin grid is incomplete or reordered.")
  source_paths <- file.path(code_root, expected_names)
  actual_pins <- vapply(source_paths, sha, character(1L), USE.NAMES = FALSE)
  if (!identical(unname(unlist(cfg$source_sha256)), actual_pins)) stop("Source checksum mismatch.")
  package_root <- find.package("wdsmatch")
  pkg_paths <- c("DESCRIPTION", "NAMESPACE", "R/wdsmatch", "R/wdsmatch.rdb", "R/wdsmatch.rdx")
  if (!identical(names(cfg$package_sha256), pkg_paths) ||
      !identical(unname(unlist(cfg$package_sha256)),
        vapply(file.path(package_root,pkg_paths),sha,character(1L),USE.NAMES=FALSE)) ||
      !identical(cfg$package_version,as.character(utils::packageVersion("wdsmatch"))))
    stop("Installed wdsmatch source/version mismatch.")
  package_sources <- sort(list.files(file.path(code_root,"R"),pattern="\\.R$",full.names=TRUE))
  package_relative <- paste0("R/",basename(package_sources))
  if (!identical(names(cfg$package_source_sha256),package_relative) ||
      !identical(unname(unlist(cfg$package_source_sha256)),
        vapply(package_sources,sha,character(1L),USE.NAMES=FALSE))) stop("Package checkout source mismatch.")
  ns <- asNamespace("wdsmatch"); source_env <- new.env(parent=ns)
  for (path in package_sources) sys.source(path,source_env)
  function_names <- Filter(function(name)is.function(source_env[[name]]),ls(source_env,all.names=TRUE))
  for (name in function_names) {
    if (!exists(name,ns,inherits=FALSE) ||
        !identical(deparse(formals(source_env[[name]])),deparse(formals(get(name,ns)))) ||
        !identical(deparse(body(source_env[[name]])),deparse(body(get(name,ns)))))
      stop("Installed/source function mismatch: ",name)
  }
  if (!identical(cfg$schema_version,1L) || !cfg$mode %in% c("runtime_pilot","calibration") ||
      !identical(cfg$index,"source_good") || cfg$radius != .25 ||
      !identical(as.integer(cfg$sample_sizes),c(500L,2000L))) stop("Unreviewed design configuration.")
  if (!identical(cfg$models,c("CorCor","MisCor")) ||
      !identical(cfg$estimands,c("PATE","PATT")) ||
      !identical(as.integer(cfg$M),c(1L,3L))) stop("Unreviewed analysis grid.")
  pilot <- cfg$mode == "runtime_pilot"
  R <- if (pilot) 4L else 1000L
  B <- if (pilot) 39L else 399L
  seed <- if (pilot) 29092026L else 29092027L
  if (cfg$datasets_per_n != R || cfg$count_draws != B || cfg$master_seed != seed ||
      cfg$refit_modulus != (if (pilot) 1L else 4L) || cfg$refit_remainder != 1L ||
      cfg$max_candidate_multiplier != 16L) stop("Unreviewed allocation, seed or candidate budget.")
  if (!is.character(cfg$reviewed_plan_sha256) || length(cfg$reviewed_plan_sha256) != 1L ||
      !grepl("^[a-f0-9]{64}$",cfg$reviewed_plan_sha256)) stop("Missing reviewed plan identity.")
  records_root <- file.path(output_root,"records")
  if (file.exists(records_root) || !dir.create(records_root)) stop("Records root must be new.")
  exclusive <- function(path, writer) {
    if (file.exists(path)) stop("Refusing to overwrite ",path)
    temporary <- paste0(path,".partial")
    if (file.exists(temporary)) stop("Existing partial write: ",temporary)
    writer(temporary)
    if (file.exists(path) || !file.rename(temporary,path)) stop("Exclusive finalization failed: ",path)
    invisible(path)
  }
  save_rds <- function(path,x) exclusive(path,function(p) saveRDS(x,p,compress=FALSE,version=3))
  save_csv <- function(path,x) exclusive(path,function(p) utils::write.csv(x,p,row.names=FALSE,na="NA"))
  save_json <- function(path,x) exclusive(path,function(p)
    jsonlite::write_json(x,p,auto_unbox=TRUE,pretty=TRUE,digits=17,na="null"))
  config_copy <- file.path(records_root,"config.json")
  exclusive(config_copy,function(p) {
    if (!file.copy(config_path,p,overwrite=FALSE)) stop("Config copy failed.")
  })
  provenance <- list(config_sha256=sha(config_path), source_sha256=cfg$source_sha256,
    package_sha256=cfg$package_sha256,package_version=cfg$package_version,
    package_source_sha256=cfg$package_source_sha256,installed_functions_verified=function_names,
    preflight_only=preflight_only,shard_index=shard_index,shard_count=shard_count,
    parent_config_sha256=cfg$parent_config_sha256,
    R_version=R.version.string,platform=R.version$platform,
    dependency_versions=list(jsonlite=as.character(utils::packageVersion("jsonlite")),
      digest=as.character(utils::packageVersion("digest"))),
    started_utc=format(Sys.time(),tz="UTC",usetz=TRUE),pid=Sys.getpid(),
    mode=cfg$mode,calibration_interpretation_permitted=!pilot,
    runtime_note="Supervisor owns wall, aggregate memory, and output limits; no internal workers.")
  save_json(file.path(records_root,"provenance.json"),provenance)
  for (path in source_paths[-length(source_paths)]) source(path,local=environment())
  # nextRNGStream allocates one independent stream for data and one for counts.
  # Future data and count allocations are fixed before observing any fit.
  RNGkind("L'Ecuyer-CMRG","Inversion","Rejection"); set.seed(seed)
  stream <- .Random.seed; allocations <- list(); j <- 0L
  for (n in cfg$sample_sizes) for (id in seq_len(R)) {
    j <- j+1L; data_stream <- stream; stream <- parallel::nextRNGStream(stream)
    count_stream <- stream; stream <- parallel::nextRNGStream(stream)
    allocations[[j]] <- list(dataset_id=sprintf("%s-n%d-id%04d",cfg$mode,n,id),
      n=as.integer(n),replication_id=id,refit=pilot || (id-1L) %% 4L == 0L,
      data_seed=data_stream,count_seed=count_stream)
  }
  save_rds(file.path(records_root,"streams.rds"),allocations)
  grid <- do.call(rbind,lapply(allocations,function(a) {
    g <- rbind(expand.grid(branch=c("actual_fitted","oracle_means_fitted_map"),
      model=cfg$models,estimand=cfg$estimands,M=cfg$M,stringsAsFactors=FALSE),
      expand.grid(branch="known_map_oracle",model="shared",estimand=cfg$estimands,
        M=cfg$M,stringsAsFactors=FALSE))
    cbind(dataset_id=a$dataset_id,n=a$n,replication_id=a$replication_id,
      requested_count_draws=if (a$refit) B else 0L,g)
  }))
  # Assign after the entire scientific grid and both streams have been fixed.
  # Per-n contiguous blocks balance costly n and original-refit allocations.
  assignment <- data.frame(global_index=seq_along(allocations),
    dataset_id=vapply(allocations,function(a)a$dataset_id,character(1L)),
    n=vapply(allocations,function(a)a$n,integer(1L)),
    replication_id=vapply(allocations,function(a)a$replication_id,integer(1L)),
    refit=vapply(allocations,function(a)a$refit,logical(1L)),stringsAsFactors=FALSE)
  assignment$owner_shard <- as.integer(1+floor((as.double(assignment$replication_id)-1)*shard_count/R))
  selected_indices <- which(assignment$owner_shard == shard_index)
  selected <- assignment[selected_indices,,drop=FALSE]
  global_count <- length(allocations)
  grid$global_index <- match(grid$dataset_id,assignment$dataset_id)
  grid$owner_shard <- assignment$owner_shard[grid$global_index]
  save_csv(file.path(records_root,"global_assignment.csv"),assignment)
  save_csv(file.path(records_root,"selected_ids.csv"),selected)
  save_csv(file.path(records_root,"global_requested_grid.csv"),grid)
  grid <- grid[grid$owner_shard == shard_index,,drop=FALSE]
  save_csv(file.path(records_root,"requested_grid.csv"),grid)
  save_json(file.path(records_root,"shard_assignment.json"),list(
    shard_index=shard_index,shard_count=shard_count,
    rule="owner_shard = 1 + floor((replication_id - 1) * shard_count / datasets_per_n)",
    global_dataset_count=global_count,selected_dataset_count=nrow(selected),
    selected_global_indices=selected$global_index,selected_dataset_ids=selected$dataset_id,
    global_assignment_sha256=sha(file.path(records_root,"global_assignment.csv")),
    global_grid_sha256=sha(file.path(records_root,"global_requested_grid.csv")),
    full_streams_sha256=sha(file.path(records_root,"streams.rds")),
    note="Each global dataset has exactly one owner for this fixed shard count; an empty shard is an explicit no-op."))
  allocations <- allocations[selected_indices]
  truth <- wm_bounded_wdsm_truth_reference("source_good",radius=.25)
  save_rds(file.path(records_root,"target.rds"),truth)
  if (preflight_only) {
    save_json(file.path(records_root,"preflight_completion.json"),
      list(status="preflight_only",datasets_generated=0L,counts_generated=0L,
        fits_executed=0L,requested_dataset_count=length(allocations),
        global_dataset_count=global_count,shard_index=shard_index,shard_count=shard_count,
        empty_shard=length(allocations)==0L,
        package_functions_verified=length(function_names),scientific_validation=FALSE))
    return(invisible(NULL))
  }
  statuses <- list()
  for (j in seq_along(allocations)) {
    a <- allocations[[j]]; directory <- file.path(records_root,a$dataset_id)
    if (!dir.create(directory)) stop("Dataset directory already exists.")
    began <- proc.time()[[3L]]; event_id <- 0L
    event <- function(stage,tag="",draw=NA_integer_) {
      event_id <<- event_id+1L
      save_json(file.path(directory,sprintf("event-%05d.json",event_id)),
        list(stage=stage,tag=tag,draw=draw,elapsed_seconds=proc.time()[[3L]]-began,
          utc=format(Sys.time(),tz="UTC",usetz=TRUE)))
    }
    save_rds(file.path(directory,"allocation.rds"),a)
    # Neither generator failure nor fitting can alter the count allocation.
    assign(".Random.seed",a$data_seed,envir=.GlobalEnv)
    data <- tryCatch(wm_bounded_wdsm_design_reference(n=a$n,index="source_good",radius=.25,
      max_candidates=16L*a$n),error=identity)
    data_end <- .Random.seed
    save_rds(file.path(directory,"dataset.rds"),list(data=data,seed_before=a$data_seed,seed_after=data_end))
    assign(".Random.seed",a$count_seed,envir=.GlobalEnv)
    counts <- if (a$refit) stats::rmultinom(B,a$n,rep(1/a$n,a$n)) else matrix(integer(),a$n,0L)
    save_rds(file.path(directory,"counts.rds"),list(counts=counts,seed_before=a$count_seed,seed_after=.Random.seed))
    event("data_and_counts_saved")
    result <- tryCatch({
      if (inherits(data,"error")) stop(data)
      wm_cal_analyze_dataset(data,counts,
        save_object=function(name,x) save_rds(file.path(directory,name),x),
        save_table=function(name,x) save_csv(file.path(directory,name),x),event=event)
    },error=identity)
    okay <- !inherits(result,"error")
    status <- list(dataset_id=a$dataset_id,status=if (okay) "completed" else "failed",
      error=if (okay) "" else conditionMessage(result),elapsed_seconds=proc.time()[[3L]]-began,
      data_sha256=sha(file.path(directory,"dataset.rds")),counts_sha256=sha(file.path(directory,"counts.rds")),
      arm_counts=if (!inherits(data,"error")) as.list(table(data$Z)) else NULL)
    if (!okay) {
      # A terminal error preserves all existing partial outputs and explicitly
      # marks every unproduced requested row/draw as unavailable. No replacement.
      wanted <- grid[grid$dataset_id == a$dataset_id,,drop=FALSE]
      suffix <- ifelse(wanted$branch == "known_map_oracle","-point.csv",
        ifelse(wanted$branch == "actual_fitted","-actual-point.csv","-oracle-fitted-point.csv"))
      prefix <- ifelse(wanted$model == "shared","known",wanted$model)
      wanted$record_file <- paste0(prefix,"-",wanted$estimand,"-",wanted$M,suffix)
      wanted$produced <- file.exists(file.path(directory,wanted$record_file))
      wanted$terminal_error <- status$error
      save_csv(file.path(directory,"terminal_record_accounting.csv"),wanted)
      if (a$refit) {
        missing <- expand.grid(model=cfg$models,estimand=cfg$estimands,M=cfg$M,draw=seq_len(B),stringsAsFactors=FALSE)
        missing$record_file <- sprintf("%s-%s-%d-count%03d.csv",missing$model,missing$estimand,missing$M,missing$draw)
        missing$produced <- file.exists(file.path(directory,missing$record_file))
        missing$terminal_error <- status$error
        save_csv(file.path(directory,"terminal_draw_accounting.csv"),missing)
      }
    }
    save_json(file.path(directory,"completion.json"),status);statuses[[j]] <- status
    cat(a$dataset_id,status$status,sprintf("%.2f seconds",status$elapsed_seconds),"\n")
    rm(data,counts,result);invisible(gc())
  }
  if (!identical(actual_pins,vapply(source_paths,sha,character(1L),USE.NAMES=FALSE)) ||
      !identical(sha(config_path),provenance$config_sha256) ||
      !identical(unname(unlist(cfg$package_source_sha256)),vapply(package_sources,sha,character(1L),USE.NAMES=FALSE)) ||
      !identical(unname(unlist(cfg$package_sha256)),vapply(file.path(package_root,pkg_paths),sha,character(1L),USE.NAMES=FALSE))) stop("Inputs changed during execution.")
  save_json(file.path(records_root,"completion.json"),list(status=if (length(allocations)) "requested_grid_processed" else "empty_shard_processed",
    shard_index=shard_index,shard_count=shard_count,global_datasets_requested=global_count,
    datasets_requested=length(allocations),datasets_terminal_errors=sum(vapply(statuses,function(x)x$status!="completed",logical(1L))),
    completed_utc=format(Sys.time(),tz="UTC",usetz=TRUE),dataset_status=statuses,
    scientific_validation=FALSE,mode=cfg$mode))
  invisible(statuses)
}

if (sys.nframe() == 0L) wm_cal_runner()
