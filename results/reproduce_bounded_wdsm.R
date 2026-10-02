#!/usr/bin/env Rscript
# Deterministic reconstruction from synthetic saved estimator/count-set records.
# No raw-data generation, model fitting, count generation or raw archive access.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)) stop("Run without arguments from the repository root.")
root <- normalizePath(".", mustWork=TRUE)
bundle <- file.path(root,"results/bounded_wdsm")
manifest <- utils::read.csv(file.path(bundle,"artifact_hashes.csv"),stringsAsFactors=FALSE)
paths <- file.path(root,manifest$file)
stopifnot(!anyDuplicated(manifest$file),all(file.exists(paths)),
          identical(unname(tools::md5sum(paths)),manifest$md5),
          all(file.info(paths)$size==manifest$bytes))
read <- function(name) utils::read.csv(file.path(bundle,name),stringsAsFactors=FALSE,
  na.strings="NA",check.names=FALSE)
records <- read("record_rows.csv.gz")
counts <- read("count_summaries.csv")
execution <- read("dataset_execution.csv")
grid <- read("requested_grid.csv.gz")
helper <- file.path(root,"simulations/summarize_bounded_wdsm_calibration.R")
had_rng <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
if(had_rng) rng <- get(".Random.seed",envir=.GlobalEnv)
source(helper,local=TRUE)
stopifnot(nrow(records)==60000L,nrow(counts)==4000L,nrow(execution)==2000L,
          nrow(grid)==40000L,!anyDuplicated(.wm_ca_key(records)),
          !anyDuplicated(.wm_ca_key(grid)),!anyDuplicated(execution$dataset_id),
          all(execution$status=="completed"),
          all(records$point_status=="completed"),all(records$interval_status=="completed"),
          all(counts$B==399),all(counts$successful_draws==399),
          all(counts$missing_draws==0),all(counts$terminal_unproduced_draws==0),
          all(counts$all_refits_succeeded),all(counts$all_one_step_available),
          sum(counts$B)==1596000)
original <- records[records$branch %in% c("actual_fitted","oracle_means_fitted_map","known_map_oracle"),]
stopifnot(setequal(.wm_ca_key(original),.wm_ca_key(grid)))
refit <- records[records$branch=="original_refit",]
stopifnot(all((refit$replication_id-1)%%4==0),
          length(unique(refit$dataset_id))==500L)
collected <- list(mode="calibration",records=records,count_summaries=counts,
  dataset_execution=execution,execution_complete=TRUE,assessment_ready=TRUE,
  draw_status=NULL)
replayed <- wm_cal_summarize(collected)
comparisons <- list(summary=c("branch","model","estimand","M","n"),
  paired_comparisons=c("comparison","model","estimand","M","n"),
  count_pair_summary=c("model","estimand","M","n"),
  membership=c("dataset_id","branch","model","estimand","M","n"),
  paired_membership=c("dataset_id","comparison","model","estimand","M","n"))
for (name in names(comparisons)) {
  ext <- if (name %in% c("membership","paired_membership")) ".csv.gz" else ".csv"
  saved <- read(paste0(name,ext)); fresh <- replayed[[name]]; fields <- comparisons[[name]]
  stopifnot(setequal(names(saved),names(fresh)))
  skey <- .wm_ca_key(saved,fields); fkey <- .wm_ca_key(fresh,fields)
  stopifnot(!anyDuplicated(skey),!anyDuplicated(fkey),setequal(skey,fkey))
  fresh <- fresh[match(skey,fkey),names(saved),drop=FALSE]
  for(field in names(saved)) {
    a <- fresh[[field]]; b <- saved[[field]]
    if(is.numeric(a)&&is.numeric(b)) {
      stopifnot(identical(is.na(a),is.na(b)),identical(is.finite(a),is.finite(b)))
      ok <- is.finite(a)&is.finite(b)
      if(any(abs(a[ok]-b[ok])>1e-9*(1+abs(b[ok])))) stop(name," numeric mismatch: ",field)
    } else if(!identical(as.character(a),as.character(b))) stop(name," mismatch: ",field)
  }
  cat(name,": ",nrow(saved)," rows, all columns agree\n",sep="")
}
stopifnot(nrow(replayed$failures)==0L,nrow(read("failures.csv"))==0L,
  identical(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE),had_rng))
if(had_rng) stopifnot(identical(get(".Random.seed",envir=.GlobalEnv),rng))
cat("PASS: outer summaries and paired memberships reconstructed; RNG unchanged.\n",
    "Raw refits/count draws were not rerun or independently reconstructed.\n",sep="")
