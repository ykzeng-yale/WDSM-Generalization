# Private original-study orchestration. Statistical procedures are the separately
# reviewed weighted/native/analytic modules; this file does not redefine them.

.ows_sha_file <- function(path) digest::digest(file=path,algo="sha256")
.ows_json <- function(object,path) {
  stopifnot(!file.exists(path))
  pending <- paste0(path,".pending")
  stopifnot(!file.exists(pending))
  jsonlite::write_json(object,pending,pretty=TRUE,auto_unbox=TRUE,digits=NA,null="null")
  stopifnot(file.rename(pending,path))
  invisible(.ows_sha_file(path))
}
.ows_ref <- function(receipt,root) {
  prefix <- paste0(normalizePath(root,mustWork=TRUE),"/")
  stopifnot(startsWith(receipt$directory,prefix))
  list(directory=substring(receipt$directory,nchar(prefix)+1L),
    file_sha256=receipt$file_sha256,receipt_sha256=receipt$receipt_sha256)
}

ows_generate_sample <- function(population,overlap,sampling_design,replicate,
                                sampler_path,population_sha256,truth,B=200L) {
  stopifnot(is.data.frame(population),nrow(population)==1000000L,
    overlap %in% c("GoodOverlap","PoorOverlap"),
    sampling_design %in% c("retrospective","prospective"),
    length(replicate)==1L,is.finite(replicate),replicate>=1L,replicate<=1000L,
    replicate==floor(replicate),identical(B,200L),
    identical(names(truth),c("PATE","PATT")),length(truth)==2L,all(is.finite(truth)),
    length(population_sha256)==1L,grepl("^[0-9a-f]{64}$",population_sha256))
  env <- new.env(parent=globalenv())
  env$config <- list(sample_seed_base=2025L,n_sample=5000L)
  env$population <- population; env$design <- sampling_design
  sys.source(sampler_path,env)
  stopifnot(is.function(env$draw_sample))
  RNGkind("Mersenne-Twister","Inversion","Rejection")
  sample_seed <- as.integer(2025L+replicate); set.seed(sample_seed)
  before <- .Random.seed
  data <- env$draw_sample(as.integer(replicate))
  after <- .Random.seed
  stopifnot(is.data.frame(data),all(c("id","Y_obs","A","survey_weight",paste0("X",1:6))%in%names(data)),
    !anyDuplicated(data$id),all(is.finite(data$survey_weight)),all(data$survey_weight>0))
  data$Y <- data$Y_obs
  count_seed <- as.integer(2000000L+100L*replicate+13L); set.seed(count_seed)
  count_before <- .Random.seed; n <- nrow(data)
  counts <- vapply(seq_len(B),function(b)as.vector(stats::rmultinom(1L,n,rep(1/n,n))),integer(n))
  list(data=data,counts=list(counts=counts,seed=count_seed,rng_before=count_before,rng_after=.Random.seed,
    column_sha256=vapply(seq_len(B),function(b)digest::digest(counts[,b],algo="sha256"),character(1))),
    context=list(overlap=overlap,sampling_design=sampling_design,replicate=as.integer(replicate),
      population_sha256=population_sha256,truth=truth,sample_seed=sample_seed,
      sample_rng_before=before,sample_rng_after=after,sampler_sha256=.ows_sha_file(sampler_path),
      sample_n=n,treated=sum(data$A==1L),input_data_sha256=digest::digest(data,algo="sha256"),
      rng_kind=RNGkind(),B=B,alpha=.05,
      scope="Source-design empirical comparison; unchanged main WM assumptions; no design-variance or hidden-confounding extension"))
}

ows_run_case <- function(input,specification,estimand,root,input_reference,native_cache,
                         native_receipts,wm,original,
                         weighted_function=ows_weighted_batch,analytic_function=ows_analytic_batch) {
  stopifnot(specification %in% c("CorCor","CorMis","MisCor","MisMis"),
    estimand %in% c("PATE","PATT"),is.environment(native_receipts))
  context <- input$context; data <- input$data; counts <- input$counts$counts
  key <- paste(estimand,specification,sep="_")
  dest <- file.path(root,"cases",key)
  if(!dir.create(dest,showWarnings=FALSE))stop("Case directory already exists or parent is unavailable")
  stage <- "weighted_batch"; checkpoints <- list(); refs <- dependencies <- list()
  committed_artifact <- NULL; retired <- list()
  tryCatch({
    weighted <- weighted_function(data,specification,estimand,context$sampling_design,
      context$replicate,counts,wm,original)
    # These task-created checkpoints survive any later error. Only after exact
    # whole-bundle disk restoration may their duplicate files be retired.
    checkpoint <- function(x,name) {
      path <- file.path(dest,paste0(name,".checkpoint.rds"))
      stopifnot(!file.exists(path)); saveRDS(x,path,version=2L,compress=FALSE)
      list(path=path,sha256=.ows_sha_file(path))
    }
    checkpoints$weighted <- checkpoint(weighted,"weighted")
    stopifnot(identical(weighted$input_data_sha256,context$input_data_sha256),
      identical(weighted$counts_sha256,digest::digest(counts,algo="sha256")),
      identical(weighted$specification,specification),identical(weighted$estimand,estimand),
      identical(weighted$sampling_design,context$sampling_design),
      identical(weighted$n,nrow(data)),identical(weighted$B,ncol(counts)))
    stage <- "analytic_batch"
    analytic <- analytic_function(data,weighted,wm)
    checkpoints$analytic <- checkpoint(analytic,"analytic")
    stage <- "native_comparators"
    native <- list()
    for(method in c("PSM","PGM","DSM","SWPSM")) {
      returned <- native_cache$get(method,specification)
      native[[method]] <- returned$value
      cache_key <- paste(estimand,returned$dependency$key,sep="_")
      if(!returned$dependency$cache_hit) {
        stopifnot(!exists(cache_key,envir=native_receipts,inherits=FALSE))
        receipt <- ows_write_artifact(returned$value,file.path(root,"native",cache_key),returned$dependency)
        assign(cache_key,.ows_ref(receipt,root),envir=native_receipts)
      }
      stopifnot(exists(cache_key,envir=native_receipts,inherits=FALSE))
      refs[[method]] <- get(cache_key,envir=native_receipts,inherits=FALSE)
      dependencies[[method]] <- returned$dependency
    }
    stage <- "records_and_reduction"
    records <- ows_case_records(data,weighted,analytic,native,context$overlap,
      context$replicate,context$truth[[estimand]],context$population_sha256)
    reductions <- ows_reduction_table(weighted)
    stage <- "case_artifact"
    bundle <- list(weighted=weighted,analytic=analytic,records=records,reductions=reductions,
      input_reference=input_reference,native_references=refs,native_dependencies=dependencies)
    receipt <- ows_write_artifact(bundle,file.path(dest,"result"),
      list(case=key,input_data_sha256=context$input_data_sha256,
        counts_sha256=weighted$counts_sha256,scope=context$scope))
    committed_artifact <- .ows_ref(receipt,root)
    # The reviewed writer has checked strict identity of this entire bundle,
    # including weighted and analytic objects and all attributes, on disk.
    stage <- "record_export"
    write.csv(records,file.path(dest,"records.csv"),row.names=FALSE)
    write.csv(reductions,file.path(dest,"reductions.csv"),row.names=FALSE)
    stage <- "retire_duplicate_checkpoints"
    for(name in names(checkpoints)) {
      x <- checkpoints[[name]]
      stopifnot(identical(.ows_sha_file(x$path),x$sha256))
      stopifnot(unlink(x$path)==0L,!file.exists(x$path))
      retired[[name]] <- list(file=basename(x$path),sha256=x$sha256,
        retained_in="result/artifact.rds",exact_bundle_roundtrip=TRUE)
    }
    result <- list(status="CASE_COMPLETE_WITH_FAILURES_RETAINED",case=key,records=nrow(records),
      point_failures=sum(records$point_status!="ok"),interval_failures=sum(records$interval_status!="ok"),
      reduction_status=as.list(setNames(reductions$result,paste0("M",reductions$M))),
      artifact=.ows_ref(receipt,root),input_reference=input_reference,native_references=refs,
      records_sha256=.ows_sha_file(file.path(dest,"records.csv")),
      reductions_sha256=.ows_sha_file(file.path(dest,"reductions.csv")),
      retired_duplicate_checkpoints=retired)
    .ows_json(result,file.path(dest,"completion.json"))
    list(records=records,reductions=reductions,receipt=result,
      completion_sha256=.ows_sha_file(file.path(dest,"completion.json")))
  },error=function(e) {
    report <- tryCatch({.ows_json(list(status="CASE_INTEGRATION_FAILED",case=key,stage=stage,
      error=conditionMessage(e),completed_checkpoints=checkpoints,native_references=refs,
      committed_artifact=committed_artifact,retired_duplicate_checkpoints=retired),
      file.path(dest,"failure.json")); "saved"},error=function(err)conditionMessage(err))
    stop(paste0("Case ",key," failed at ",stage,": ",conditionMessage(e),
      "; failure receipt: ",report),call.=FALSE)
  })
}

ows_run_sample <- function(input,directory,wm,original,repo) {
  context <- input$context; data <- input$data; counts <- input$counts$counts
  stopifnot(context$overlap%in%c("GoodOverlap","PoorOverlap"),
    context$sampling_design%in%c("retrospective","prospective"),
    identical(context$B,200L),identical(context$alpha,.05),
    identical(dim(counts),c(nrow(data),200L)),all(is.finite(counts)),
    all(counts>=0&counts==floor(counts)),all(colSums(counts)==nrow(data)),
    identical(context$input_data_sha256,digest::digest(data,algo="sha256")),
    identical(input$counts$column_sha256,
      vapply(seq_len(ncol(counts)),function(b)digest::digest(counts[,b],algo="sha256"),character(1))),
    !file.exists(directory),dir.create(directory,showWarnings=FALSE))
  root <- normalizePath(directory,mustWork=TRUE)
  native_receipts <- new.env(hash=TRUE,parent=emptyenv())
  results <- list(); ledgers <- list(); stage <- "sample_initialization"; input_reference <- NULL
  tryCatch({
    stopifnot(dir.create(file.path(root,"cases")),dir.create(file.path(root,"native")))
    stage <- "input_artifact"
    input_receipt <- ows_write_artifact(input,file.path(root,"input"),context[c("overlap","sampling_design","replicate","population_sha256")])
    input_reference <- .ows_ref(input_receipt,root)
    stage <- "started_receipt"
    .ows_json(list(status="SAMPLE_STARTED",input=input_reference,context=context[c("overlap","sampling_design","replicate","B")]),
      file.path(root,"started.json"))
    for(estimand in c("PATE","PATT")) {
      native_cache <- ows_native_cache(data,estimand,context$sampling_design,context$replicate,
        context$B,original,repo,context$alpha)
      for(spec in c("CorCor","CorMis","MisCor","MisMis")) {
        key <- paste(estimand,spec,sep="_")
        stage <- key
        results[[key]] <- ows_run_case(input,spec,estimand,root,input_reference,native_cache,
          native_receipts,wm,original)
        cat("Completed",context$overlap,context$sampling_design,"rep",context$replicate,key,"\n")
        invisible(gc(FALSE))
      }
      stopifnot(native_cache$size()==10L,length(native_cache$ledger())==16L)
      ledgers[[estimand]] <- native_cache$ledger()
    }
    stage <- "sample_summary"
    records <- do.call(rbind,lapply(results,`[[`,"records")); rownames(records) <- NULL
    stopifnot(nrow(records)==152L,length(unique(records$pairing_id))==1L,
      !anyDuplicated(records[c("estimand","scenario","method")]))
    write.csv(records,file.path(root,"records.csv"),row.names=FALSE)
    ledger_receipt <- ows_write_artifact(ledgers,file.path(root,"native_ledger"))
    completion <- list(status="SAMPLE_COMPLETE_WITH_FAILURES_RETAINED",case_count=length(results),records=nrow(records),
      point_failures=sum(records$point_status!="ok"),interval_failures=sum(records$interval_status!="ok"),
      input=input_reference,ledger=.ows_ref(ledger_receipt,root),
      records_sha256=.ows_sha_file(file.path(root,"records.csv")),
      cases=lapply(results,function(x)list(completion_sha256=x$completion_sha256,receipt=x$receipt)),
      native_unique_executions=length(ls(native_receipts,all.names=TRUE)),
      native_requested_procedures=32L,scientific_performance_reviewed=FALSE)
    stopifnot(completion$native_unique_executions==20L)
    .ows_json(completion,file.path(root,"completion.json"))
    completion
  },error=function(e) {
    report <- tryCatch({.ows_json(list(status="SAMPLE_INTEGRATION_FAILED",stage=stage,
      error=conditionMessage(e),completed_cases=names(results),input=input_reference),
      file.path(root,"failure.json")); "saved"},error=function(err)conditionMessage(err))
    stop(paste0("Sample integration failed: ",conditionMessage(e),"; failure receipt: ",report),call.=FALSE)
  })
}
