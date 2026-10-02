# Source-study record assembly. No fitting, sampling, RNG, or file mutation.
ows_method_ids <- function() c(paste0("WDSM_M",c(1L,3L,5L)),
  as.vector(outer(c("PS","DSM","X6"),c(1L,3L,5L),function(f,m)paste0("WM_",f,"_M",m))),
  paste0("WM_DSM_analytic_M",c(1L,3L,5L)),"PSM_M1","PGM_M1","DSM_M1","SWPSM_full")

ows_reporting_plan <- function(requested=1000L) {
  stopifnot(length(requested)==1L,is.finite(requested),requested>=2L,requested==floor(requested))
  plan <- expand.grid(sampling_design=c("retrospective","prospective"),
    overlap=c("GoodOverlap","PoorOverlap"),estimand=c("PATE","PATT"),
    scenario=c("CorCor","CorMis","MisCor","MisMis"),method=ows_method_ids(),stringsAsFactors=FALSE)
  plan$requested <- as.integer(requested); plan$interval_required <- TRUE
  plan$coverage_rule <- "inclusive"
  plan
}

ows_warning_text <- function(x) {
  if (!is.list(x)) return(character())
  own <- if (is.character(x$warnings)) x$warnings else character()
  unique(c(own,unlist(lapply(x,ows_warning_text),use.names=FALSE)))
}

ows_reduction_table <- function(weighted,tolerance=1e-10) {
  stopifnot(length(tolerance)==1L,is.finite(tolerance),tolerance>0)
  canonical <- function(x,Z,M,pATE) {
    stopifnot(is.list(x),length(x)==length(Z))
    lapply(seq_along(Z),function(i) {
      v <- x[[i]]
      if(!pATE&&Z[i]==0L) {stopifnot(length(v)==0L); return(integer())}
      stopifnot(is.numeric(v),length(v)==M,all(is.finite(v)),all(v==floor(v)),
        all(v>=1L&v<=length(Z)),!anyDuplicated(v),all(Z[v]!=Z[i]))
      sort(as.integer(v))
    })
  }
  compare <- function(a,b) {
    stopifnot(length(a)==length(b))
    common <- is.finite(a)&is.finite(b)
    list(paired=sum(common),different_availability=sum(is.finite(a)!=is.finite(b)),
      max_abs=if(any(common))max(abs(a[common]-b[common])) else NA_real_,
      exceeds=which(common & abs(a-b)>tolerance*(1+pmax(abs(a),abs(b)))))
  }
  output <- lapply(c(1L,3L,5L),function(M) {
    oldkey <- paste0("WDSM_M",M); newkey <- paste0("WM_DSM_M",M)
    old <- weighted$points[[oldkey]]; new <- weighted$points[[newkey]]
    a <- weighted$summaries[[oldkey]]; b <- weighted$summaries[[newkey]]
    graph <- NA; raw <- list(paired=0L,max_abs=NA_real_,exceeds=integer())
    if (isTRUE(old$ok)&&isTRUE(new$ok)) {
      Z <- new$value$data$Z
      oldgraph <- lapply(seq_len(weighted$n),function(i)
        if(Z[i]==1L)old$value$graph$matches_0[[i]] else
          if(weighted$estimand=="PATE")old$value$graph$matches_1[[i]] else integer())
      graph <- identical(canonical(oldgraph,Z,M,weighted$estimand=="PATE"),
        canonical(new$value$graph$neighbors,Z,M,weighted$estimand=="PATE"))
      raw <- compare(old$value$raw_estimate,new$value$raw_estimate)
    }
    point <- compare(a$estimate,b$estimate); draws <- compare(a$draws,b$draws)
    variance <- compare(a$variance,b$variance); ci <- compare(a$conf.int,b$conf.int)
    full <- isTRUE(graph)&&point$paired==1L&&raw$paired==1L&&draws$paired==weighted$B&&
      variance$paired==1L&&ci$paired==2L&&
      !length(c(raw$exceeds,point$exceeds,draws$exceeds,variance$exceeds,ci$exceeds))
    data.frame(M=M,original_point_status=a$point_status,modern_point_status=b$point_status,
      original_interval_status=a$interval_status,modern_interval_status=b$interval_status,
      donor_sets_identical=graph,raw_point_max_abs=raw$max_abs,point_max_abs=point$max_abs,
      paired_draws=draws$paired,draw_availability_mismatches=draws$different_availability,
      draw_max_abs=draws$max_abs,draws_exceeding_tolerance=length(draws$exceeds),
      variance_max_abs=variance$max_abs,ci_max_abs=ci$max_abs,tolerance=tolerance,
      result=if(full)"complete_numeric_agreement" else
        if(point$paired==0L)"unavailable" else "requires_investigation",stringsAsFactors=FALSE)
  })
  do.call(rbind,output)
}

ows_case_records <- function(data,weighted,analytic,native,overlap,replicate,target,population_sha256) {
  stopifnot(overlap %in% c("GoodOverlap","PoorOverlap"),
    length(replicate)==1L,is.finite(replicate),replicate>=1L,replicate==floor(replicate),
    length(target)==1L,is.finite(target),length(population_sha256)==1L,
    grepl("^[0-9a-f]{64}$",population_sha256),identical(nrow(data),weighted$n),
    is.numeric(weighted$B),length(weighted$B)==1L,is.finite(weighted$B),
    weighted$B>=2L,weighted$B==floor(weighted$B),
    all(c("PSM","PGM","DSM","SWPSM")%in%names(native)))
  # Bind pairing to actual observed rows, not only a nominal replication label.
  fields <- c("id","Y","A",paste0("X",1:6),"survey_weight")
  stopifnot(all(fields%in%names(data)),!anyDuplicated(data$id))
  actual_hash <- digest::digest(data,algo="sha256")
  stopifnot(identical(weighted$input_data_sha256,actual_hash),
    all(vapply(native,function(x)identical(x$input_data_sha256,actual_hash),logical(1))))
  pairing <- digest::digest(list(population_sha256=population_sha256,
    sampling_design=weighted$sampling_design,data=data[fields]),algo="sha256")
  rows <- lapply(ows_method_ids(),function(method) {
    raw <- NA_real_; draws <- 0L; requested <- 0L
    if (method%in%names(weighted$summaries)) {
      s <- weighted$summaries[[method]]; p <- weighted$points[[method]]
      family <- if(startsWith(method,"WDSM"))"original" else sub("_M[0-9]+$","",sub("^WM_","",method))
      bases <- if(family=="PS")weighted$bases[c("PS","PS_means")] else
        if(family=="X6")weighted$bases["X6_means"] else weighted$bases[family]
      warnings <- ows_warning_text(list(p,bases,lapply(weighted$refit_diagnostics,function(x)x$fits[[family]])))
      if(p$ok)raw <- p$value$raw_estimate
      variance_method <- "fixed_reuse_prediction_refit_B_divisor"
      error <- s$interval_error; scope <- s$scope
      requested <- s$B; draws <- sum(is.finite(s$draws)&s$errors=="")
    } else if (method%in%names(analytic$rows)) {
      s <- analytic$rows[[method]]
      p <- weighted$points[[s$point_source]]
      if(p$ok)raw <- p$value$raw_estimate
      warnings <- ows_warning_text(list(p,weighted$bases$DSM,
        analytic$nuisance,analytic$assemblies[[method]]))
      variance_method <- s$variance_method; error <- s$interval_error; scope <- s$scope
    } else {
      family <- sub("_(M1|full)$","",method); x <- native[[family]]
      interval <- x$interval_status=="ok"
      s <- list(estimate=x$estimate,point_status=x$point_status,interval_status=x$interval_status,
        variance=if(interval)x$result$value$variance else NA_real_,
        conf.int=if(interval)c(x$result$value$ci_lower,x$result$value$ci_upper) else c(NA_real_,NA_real_))
      warnings <- x$warning_inventory; error <- x$result$error
      variance_method <- if(family=="SWPSM")"full_matching_subclass_sandwich" else "vendored_dsmatch_bootstrap_B_minus_1"
      scope <- "Unchanged original comparator procedure; algorithm-specific geometry and inference"
      if(family!="SWPSM") {
        requested <- weighted$B
        if(!is.null(x$bootstrap))stopifnot(is.numeric(x$bootstrap$requested),
          length(x$bootstrap$requested)==1L,is.finite(x$bootstrap$requested),
          x$bootstrap$requested==weighted$B)
        if(interval)stopifnot(!is.null(x$bootstrap),
          length(x$bootstrap$draws)==weighted$B,all(is.finite(x$bootstrap$draws)),
          length(x$bootstrap$calls)==weighted$B+1L,
          all(vapply(x$bootstrap$calls,function(z)isTRUE(z$ok),logical(1))))
        draws <- sum(is.finite(x$bootstrap$draws))
      }
    }
    if(s$point_status=="ok")stopifnot(length(s$estimate)==1L,is.finite(s$estimate))
    if(s$interval_status=="ok")stopifnot(s$point_status=="ok",is.finite(s$variance),s$variance>=0,
      length(s$conf.int)==2L,all(is.finite(s$conf.int)),s$conf.int[1L]<=s$conf.int[2L])
    data.frame(sampling_design=weighted$sampling_design,overlap=overlap,
      estimand=weighted$estimand,scenario=weighted$specification,replicate=as.integer(replicate),
      method=method,n=nrow(data),target=target,estimate=s$estimate,raw_estimate=raw,
      variance=s$variance,lower=s$conf.int[1L],upper=s$conf.int[2L],
      status=if(s$point_status=="ok"&&s$interval_status=="ok")"ok" else "failed",
      point_status=s$point_status,interval_status=s$interval_status,
      warning_count=length(unique(warnings)),error=error,pairing_id=pairing,
      bootstrap_requested=requested,bootstrap_draws_available=draws,
      variance_method=variance_method,inference_scope=scope,stringsAsFactors=FALSE,row.names=NULL)
  })
  records <- do.call(rbind,rows)
  stopifnot(nrow(records)==19L,!anyDuplicated(records$method))
  records
}

ows_report <- function(records,requested=1000L) wm_paired_comparison_metrics(records,
  ows_reporting_plan(requested),c("sampling_design","overlap","estimand"),"CorCor","PSM_M1")
