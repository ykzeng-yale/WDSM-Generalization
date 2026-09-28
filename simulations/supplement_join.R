#!/usr/bin/env Rscript
# Source benchmarks.R, first_stage.R, first_stage_benchmarks.R and supplement.R.
wm_supplement_join <- function(summary, geometry=list(), relative_goal=.01) {
  required<-c("scenario","method","estimand","correction","n","analysis_n","d","M","target",
    "requested","completed","failed","empirical_root_n_variance","empirical_root_n_variance_mcse",
    "mean_root_n_variance","mean_root_n_variance_mcse")
  if (!is.data.frame(summary) || !nrow(summary) || !all(required %in% names(summary)) ||
      anyDuplicated(summary[c("scenario","method","estimand","correction")])) stop("Invalid supplemental summary")
  if (!is.numeric(relative_goal) || length(relative_goal)!=1L || !is.finite(relative_goal) || relative_goal<=0) stop("Invalid geometry goal")
  if (inherits(geometry,"wm_geometry") || inherits(geometry,"wm_geometry_alpha"))
    geometry<-setNames(list(geometry),sprintf("M%d_d%d",geometry$M,geometry$d))
  if (!is.list(geometry) || (length(geometry) && (is.null(names(geometry)) || anyDuplicated(names(geometry))))) stop("Named geometry list required")
  config<-wm_supplement_config();addition<-list()
  for (s in unique(summary$scenario)) {
    cfg<-config[config$scenario==s,,drop=FALSE];block<-summary[summary$scenario==s,,drop=FALSE]
    if (nrow(cfg)!=1L) stop("Unknown supplemental cell")
    expected<-if (cfg$branch=="split") data.frame(estimand=c("PATE","potential0","potential1"),method="stabilized",correction="oracle") else
      expand.grid(estimand=c("PATE","PATT"),method=c("self_normalized","stabilized"),correction=c("oracle","feasible"))
    key<-function(x) paste(x$estimand,x$method,x$correction)
    if (nrow(block)!=nrow(expected) || !setequal(key(block),key(expected))) stop("Incomplete supplemental comparison grid")
  }
  for (i in seq_len(nrow(summary))) {
    row<-summary[i,];cfg<-config[config$scenario==row$scenario,,drop=FALSE]
    if (nrow(cfg)!=1L || anyNA(row[c("requested","completed","failed","target")]) ||
        row$requested!=row$completed+row$failed || row$requested<1 || min(row$completed,row$failed)<0) stop("Invalid cell or counts")
    for (field in c("n","d","M")) if (!isTRUE(row[[field]]==cfg[[field]])) stop("Cell configuration differs")
    size<-if (row$estimand=="potential0") cfg$n-floor(cfg$n*cfg$fraction1) else
      if (row$estimand=="potential1") floor(cfg$n*cfg$fraction1) else cfg$n
    if (!isTRUE(row$analysis_n==size)) stop("Benchmark analysis sample size differs")
    d<-if (cfg$branch=="split") 2L else cfg$d;key<-sprintf("M%d_d%d",cfg$M,d);geo<-geometry[[key]]
    table<-wm_supplement_benchmarks(cfg,geo)
    b<-table[table$method==row$method & table$estimand==row$estimand,,drop=FALSE]
    if (nrow(b)!=1L || abs(b$target-row$target)>1e-12) stop("Benchmark scope or target differs")
    exact<-cfg$branch=="strata" || (cfg$branch=="spline" && cfg$d==1L) || row$estimand=="potential0"
    g<-if (exact) list(status="exact_1d",below_jensen=FALSE) else .wm_fs_geometry(cfg$M,d,geo)
    relative<-b$geometry_mcse/b$root_n_variance
    available<-is.finite(b$root_n_variance) && is.finite(b$geometry_mcse)
    resolved<-available && g$status %in% c("exact_1d","monte_carlo_estimate") &&
      !isTRUE(g$below_jensen) && is.finite(relative) && relative<=relative_goal
    empirical_difference<-row$empirical_root_n_variance-b$root_n_variance
    empirical_mcse<-sqrt(row$empirical_root_n_variance_mcse^2+b$geometry_mcse^2)
    fitted_difference<-row$mean_root_n_variance-b$root_n_variance
    fitted_mcse<-sqrt(row$mean_root_n_variance_mcse^2+b$geometry_mcse^2)
    addition[[i]]<-data.frame(benchmark_root_n_variance=b$root_n_variance,benchmark_geometry_mcse=b$geometry_mcse,
      geometry_key=if (exact) sprintf("M%d_d1",cfg$M) else key,geometry_status=if (available) g$status else "benchmark_unavailable",
      relative_geometry_mcse=relative,geometry_goal_met=resolved,geometry_below_jensen=g$below_jensen,
      empirical_to_benchmark=row$empirical_root_n_variance/b$root_n_variance,
      fitted_to_benchmark=row$mean_root_n_variance/b$root_n_variance,
      empirical_minus_benchmark=empirical_difference,empirical_combined_mcse=empirical_mcse,
      empirical_discrepancy_z=if (resolved && is.finite(empirical_mcse) && empirical_mcse>0) empirical_difference/empirical_mcse else NA_real_,
      fitted_minus_benchmark=fitted_difference,fitted_combined_mcse=fitted_mcse,
      fitted_discrepancy_z=if (resolved && is.finite(fitted_mcse) && fitted_mcse>0) fitted_difference/fitted_mcse else NA_real_)
  }
  cbind(summary,do.call(rbind,addition))
}

if (sys.nframe()==0L) {
  args<-commandArgs(trailingOnly=TRUE)
  if (length(args)<2L) stop("usage: summary.csv output.csv [geometry.rds ...]")
  script<-normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)))
  helpers<-file.path(dirname(script),c("first_stage.R","first_stage_benchmarks.R","benchmarks.R","supplement.R"))
  for (file in helpers) source(file)
  sidecar<-paste0(args[2],".metadata.rds")
  if (any(file.exists(c(args[2],sidecar)))) stop("Output or companion metadata exists")
  geometry<-list()
  if (length(args)>2L) for (path in args[-c(1,2)]) {
    value<-readRDS(path)
    if (inherits(value,"wm_geometry") || inherits(value,"wm_geometry_alpha")) value<-setNames(list(value),sprintf("M%d_d%d",value$M,value$d))
    if (length(intersect(names(value),names(geometry)))) stop("Duplicate geometry key")
    geometry<-c(geometry,value)
  }
  result<-wm_supplement_join(utils::read.csv(args[1],stringsAsFactors=FALSE),geometry)
  dir.create(dirname(args[2]),recursive=TRUE,showWarnings=FALSE)
  utils::write.csv(result,args[2],row.names=FALSE)
  saveRDS(list(input_md5=tools::md5sum(args[-2]),code_md5=tools::md5sum(c(script,helpers)),
    output_md5=tools::md5sum(args[2]),relative_geometry_goal=.01,
    interpretation="Geometry and estimator Monte Carlo errors are separate; discrepancies from asymptotic limits are descriptive, not finite-n theorem tests.",
    session=utils::sessionInfo(),completed=Sys.time()),sidecar)
  cat("Joined",nrow(result),"supplemental summary groups without changing counts.\n")
}
