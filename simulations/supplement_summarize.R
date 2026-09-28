#!/usr/bin/env Rscript
# Source first_stage.R, first_stage_summarize.R and supplement.R first.
# CLI: output.csv input1.csv [input2.csv ...]; --pilot as final arg permits subsets.
wm_supplement_validate <- function(r) {
  fields <- c("scenario","branch","replication","n","analysis_n","d","M","fraction1","base_seed",
    "method","estimand","correction","target","status","estimate","root_n_variance",
    "sampling_variance","lower","upper","clipped","rng_state","warning","error")
  if (!is.data.frame(r) || !nrow(r) || !all(fields %in% names(r))) stop("Missing supplemental columns")
  key <- c("scenario","replication","method","estimand","correction")
  if (anyNA(r[key]) || anyDuplicated(r[key])) stop("Duplicated or missing supplemental keys")
  config <- wm_supplement_config()
  for (scenario in unique(r$scenario)) {
    cfg <- config[config$scenario == scenario,,drop=FALSE];x <- r[r$scenario == scenario,,drop=FALSE]
    if (nrow(cfg) != 1L) stop("Unknown supplemental scenario")
    for (field in c("branch","n","d","M"))
      if (anyNA(x[[field]]) || any(x[[field]] != cfg[[field]])) stop("Supplemental configuration differs")
    if (anyNA(x$base_seed) || any(x$base_seed != cfg$seed) ||
        any(!is.finite(x$replication)) || any(x$replication != floor(x$replication)) ||
        any(x$replication < 1 | x$replication > cfg$replications)) stop("Invalid supplemental seed or replication")
    split <- cfg$branch == "split"
    if (split && (anyNA(x$fraction1) || any(abs(x$fraction1-cfg$fraction1)>1e-14))) stop("Block fraction differs")
    if (!split && any(!is.na(x$fraction1))) stop("Unexpected block fraction")
    expected <- if (split) data.frame(estimand=c("PATE","potential0","potential1"),method="stabilized",correction="oracle") else
      expand.grid(estimand=c("PATE","PATT"),method=c("self_normalized","stabilized"),correction=c("oracle","feasible"))
    pairkey <- function(y) paste(y$estimand,y$method,y$correction)
    for (i in unique(x$replication)) {
      z <- x[x$replication == i,,drop=FALSE]
      if (nrow(z) != nrow(expected) || !setequal(pairkey(z),pairkey(expected))) stop("Incomplete supplemental replication")
      if (anyNA(z$rng_state) || length(unique(z$rng_state)) != 1L || any(!nzchar(z$rng_state))) stop("Paired RNG states differ")
    }
    size <- ifelse(x$estimand == "potential1",floor(cfg$n*cfg$fraction1),
      ifelse(x$estimand == "potential0",cfg$n-floor(cfg$n*cfg$fraction1),cfg$n))
    if (!split) size <- rep(cfg$n,nrow(x))
    if (anyNA(x$analysis_n) || any(x$analysis_n != size)) stop("Analysis sample size differs")
    target <- if (cfg$branch == "strata") ifelse(x$estimand == "PATE",393/281,41/27) else
      ifelse(x$estimand == "potential1",3,1.5)
    if (anyNA(x$target) || any(abs(x$target-target)>1e-12)) stop("Supplemental target differs")
  }
  if (anyNA(r$status) || any(!r$status %in% c("ok","error"))) stop("Invalid status")
  ok <- r$status == "ok";numeric <- c("estimate","root_n_variance","sampling_variance","lower","upper")
  if (any(!is.finite(as.matrix(r[ok,numeric,drop=FALSE]))) || any(r$root_n_variance[ok]<0) ||
      any(r$sampling_variance[ok]<0) || any(!is.na(as.matrix(r[!ok,numeric,drop=FALSE])))) stop("Invalid successful/failed numeric result")
  if (any(!ok & (is.na(r$error) | !nzchar(r$error)))) stop("Failed record lacks reason")
  if (any(abs(r$root_n_variance[ok]-r$analysis_n[ok]*r$sampling_variance[ok])>1e-10*pmax(1,r$root_n_variance[ok])))
    stop("Component variance scaling differs")
  half <- stats::qnorm(.975)*sqrt(r$sampling_variance[ok])
  if (any(abs(r$lower[ok]-(r$estimate[ok]-half))>1e-10*pmax(1,abs(r$estimate[ok]))) ||
      any(abs(r$upper[ok]-(r$estimate[ok]+half))>1e-10*pmax(1,abs(r$estimate[ok])))) stop("Interval differs")
  if (any(!is.finite(r$clipped)) || any(r$clipped<0) || any(r$clipped!=floor(r$clipped))) stop("Invalid clipping count")
  for (id in unique(paste(r$scenario[r$branch == "split"],r$replication[r$branch == "split"]))) {
    x <- r[paste(r$scenario,r$replication) == id,,drop=FALSE]
    if (all(x$status == "ok")) {
      a <- x[x$estimand == "PATE",];b <- x[x$estimand == "potential0",];c <- x[x$estimand == "potential1",]
      if (abs(a$estimate-(c$estimate-b$estimate))>1e-10 ||
          abs(a$sampling_variance-(c$sampling_variance+b$sampling_variance))>1e-10) stop("Split component identity differs")
    }
  }
  invisible(r)
}

wm_supplement_read_batches <- function(inputs, require_complete=TRUE) {
  inputs <- normalizePath(inputs,mustWork=TRUE)
  if (!length(inputs) || anyDuplicated(inputs)) stop("Missing or duplicate input paths")
  rows <- list();hash <- design <- NULL
  for (path in inputs) {
    r <- utils::read.csv(path,stringsAsFactors=FALSE,colClasses=c(warning="character",error="character",fraction1="numeric"))
    wm_supplement_validate(r);m <- readRDS(paste0(path,".metadata.rds"))
    if (length(unique(r$scenario)) != 1L) stop("Batch contains multiple scenarios")
    if (!identical(m$schema_version,"supplement_v1") || is.null(m$completed) || !isTRUE(m$code_unchanged)) stop("Incomplete batch")
    cfg <- wm_supplement_config();cfg <- cfg[cfg$scenario == unique(r$scenario),,drop=FALSE]
    if (nrow(cfg) != 1L || !isTRUE(all.equal(m$config,cfg,check.attributes=FALSE)) ||
        !isTRUE(m$requested_records==nrow(r)) || !isTRUE(m$records_written==nrow(r)) ||
        !identical(m$columns,names(r)) || !isTRUE(m$last_completed_replication==m$last) ||
        !setequal(unique(r$replication),seq.int(m$first,m$last))) stop("Batch metadata differs")
    current <- .wm_fss_hash(m$code_md5,"code");current_design <- .wm_fss_hash(m$design_md5,"design")
    if (is.null(hash)) {hash<-current;design<-current_design}
    if (!identical(hash,current) || !identical(design,current_design)) stop("Mixed simulation code or design hashes")
    rows[[length(rows)+1L]] <- r
  }
  r <- do.call(rbind,rows);rownames(r)<-NULL;wm_supplement_validate(r)
  if (require_complete) {
    config <- wm_supplement_config()
    if (!setequal(unique(r$scenario),config$scenario)) stop("Missing supplemental cells")
    for (s in config$scenario) if (!setequal(unique(r$replication[r$scenario == s]),1:200)) stop("Missing planned replication")
  }
  attr(r,"batch_provenance") <- list(inputs=inputs,input_md5=tools::md5sum(inputs),
    metadata_md5=tools::md5sum(paste0(inputs,".metadata.rds")),code_md5=hash,design_md5=design,require_complete=require_complete)
  r
}

wm_supplement_summarize <- function(records) {
  wm_supplement_validate(records)
  avg <- function(x) if (length(x)) mean(x) else NA_real_
  mcse <- function(x) if (length(x)>1L) stats::sd(x)/sqrt(length(x)) else NA_real_
  key <- c("scenario","method","estimand","correction");groups<-unique(records[key]);rows<-list()
  for (i in seq_len(nrow(groups))) {
    use <- .wm_fss_key(records,key) == .wm_fss_key(groups[i,,drop=FALSE],key)
    r <- records[use,,drop=FALSE];x <- r[r$status == "ok",,drop=FALSE];R<-nrow(x);n<-r$analysis_n[1]
    v <- if (R>1L) stats::var(x$estimate) else NA_real_
    sq <- if (R>1L) (x$estimate-mean(x$estimate))^2*R/(R-1) else numeric()
    meanv <- avg(x$sampling_variance);ratio<-if (is.finite(v) && v>0) meanv/v else NA_real_
    ratio_if <- if (is.finite(ratio)) (x$sampling_variance-meanv)/v-meanv/v^2*(sq-v) else numeric()
    hit <- x$lower<=x$target & x$upper>=x$target;coverage<-avg(hit)
    cp <- if (R) stats::binom.test(sum(hit),R)$conf.int else c(NA_real_,NA_real_)
    pair <- records[records$scenario==r$scenario[1] & records$method==r$method[1] & records$estimand==r$estimand[1],]
    a<-pair[pair$correction=="feasible",];b<-pair[pair$correction=="oracle",]
    joint<-logical();delta<-numeric()
    if (nrow(a)) {b<-b[match(a$replication,b$replication),];joint<-a$status=="ok" & b$status=="ok";delta<-a$estimate[joint]-b$estimate[joint]}
    rows[[i]] <- cbind(groups[i,],data.frame(branch=r$branch[1],n=r$n[1],analysis_n=n,d=r$d[1],M=r$M[1],
      fraction1=r$fraction1[1],target=r$target[1],requested=nrow(r),completed=R,failed=nrow(r)-R,
      failure_rate=mean(r$status!="ok"),mean_estimate=avg(x$estimate),bias=avg(x$estimate-x$target),bias_mcse=mcse(x$estimate),
      rmse=if (R) sqrt(mean((x$estimate-x$target)^2)) else NA_real_,
      empirical_root_n_variance=n*v,empirical_root_n_variance_mcse=n*mcse(sq),
      mean_root_n_variance=avg(x$root_n_variance),mean_root_n_variance_mcse=mcse(x$root_n_variance),
      variance_ratio=ratio,variance_ratio_mcse=mcse(ratio_if),coverage=coverage,
      coverage_mcse=if (R) sqrt(coverage*(1-coverage)/R) else NA_real_,coverage_exact95_lower=cp[1],coverage_exact95_upper=cp[2],
      mean_interval_width=avg(x$upper-x$lower),mean_interval_width_mcse=mcse(x$upper-x$lower),
      clipped_records=sum(r$clipped>0),warning_records=sum(!is.na(r$warning)&nzchar(r$warning)),
      paired_requested=if (nrow(a)) nrow(a) else NA_integer_,paired_completed=if (nrow(a)) sum(joint) else NA_integer_,
      feasible_minus_oracle_point=avg(delta),feasible_minus_oracle_point_mcse=mcse(delta)))
  }
  out<-do.call(rbind,rows);rownames(out)<-NULL;out
}

if (sys.nframe() == 0L) {
  args<-commandArgs(trailingOnly=TRUE);pilot<-length(args)>0L && tail(args,1L)=="--pilot"
  if (pilot) args<-head(args,-1L)
  if (length(args)<2L) stop("usage: output.csv input1.csv [...] [--pilot]")
  script<-normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)))
  helpers<-file.path(dirname(script),c("first_stage.R","first_stage_summarize.R","supplement.R"))
  for (file in helpers) source(file)
  if (file.exists(args[1]) || file.exists(paste0(args[1],".metadata.rds"))) stop("Output exists")
  records<-wm_supplement_read_batches(args[-1],require_complete=!pilot)
  result<-wm_supplement_summarize(records)
  dir.create(dirname(args[1]),recursive=TRUE,showWarnings=FALSE)
  utils::write.csv(result,args[1],row.names=FALSE)
  saveRDS(list(provenance=attr(records,"batch_provenance"),code_md5=tools::md5sum(c(script,helpers)),
    output_md5=tools::md5sum(args[1]),session=utils::sessionInfo(),completed=Sys.time()),paste0(args[1],".metadata.rds"))
  cat("Summarized",nrow(records),"supplemental records in",nrow(result),"groups; all failures retained.\n")
}
