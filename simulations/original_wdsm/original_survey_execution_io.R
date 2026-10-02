# Private execution helpers for the original-survey study.
# Source the reviewed native adapter and lossless-storage helper first.
# Caching never crosses data/estimand/design/replication/runtime instances.

ows_native_cache <- function(data, estimand, sampling_design, replicate, B,
                             original, repo, alpha=.05,
                             adapter=ows_native_comparator) {
  stopifnot(is.data.frame(data), estimand %in% c("PATE","PATT"),
    sampling_design %in% c("retrospective","prospective"),
    length(replicate)==1L,is.finite(replicate),replicate>=1L,replicate<=1000L,
    replicate==floor(replicate),length(B)==1L,is.finite(B),B>=2L,B==floor(B),
    length(alpha)==1L,is.finite(alpha),alpha>0,alpha<1,is.function(adapter))
  # This cache is private to one immutable call context. Source/runtime pins
  # are an additional responsibility of the enclosing formal-run manifest.
  data_sha <- digest::digest(data,algo="sha256")
  force(original); force(repo); force(adapter)
  values <- new.env(hash=TRUE,parent=emptyenv())
  events <- list()
  get <- function(method,specification) {
    stopifnot(length(method)==1L,method %in% c("PSM","PGM","DSM","SWPSM"),
      length(specification)==1L,specification %in% c("CorCor","CorMis","MisCor","MisMis"))
    ps <- specification %in% c("CorCor","CorMis")
    pg <- specification %in% c("CorCor","MisCor")
    # The pinned source consumes only PS for PSM/SWPSM, only PG for PGM,
    # and both for DSM. SWPSM's outcome model always includes the interaction.
    canonical <- switch(method,
      PSM=if(ps)"CorCor" else "MisCor",
      SWPSM=if(ps)"CorCor" else "MisCor",
      PGM=if(pg)"CorCor" else "CorMis", DSM=specification)
    key <- paste(method,canonical,sep="_")
    seed <- as.integer(2000000L+100L*replicate+
      10L*c(DSM=2L,PGM=3L,PSM=4L,SWPSM=5L)[[method]]+
      if(method=="SWPSM")0L else 1L)
    hit <- exists(key,envir=values,inherits=FALSE)
    if (!hit) {
      RNGkind("Mersenne-Twister","Inversion","Rejection")
      set.seed(seed)
      before <- .Random.seed
      value <- adapter(data,method,estimand,sampling_design,canonical,B,original,repo,alpha)
      stopifnot(identical(value$input_data_sha256,data_sha),
        identical(value$rng_before,before),identical(value$rng_after,.Random.seed))
      # Failed prescribed fits are cacheable too: no retry/reseed to hide them.
      assign(key,list(value=value,executed_specification=canonical,seed=seed),envir=values)
    }
    saved <- base::get(key,envir=values,inherits=FALSE)
    dependency <- list(key=key,method=method,requested_specification=specification,
      executed_specification=saved$executed_specification,cache_hit=hit,
      input_data_sha256=data_sha,estimand=estimand,sampling_design=sampling_design,
      replicate=as.integer(replicate),B=as.integer(B),alpha=alpha,seed=saved$seed,
      effective_inputs=switch(method,PGM="PG specification",DSM="PS and PG specifications",
        "PS specification; source outcome rule unchanged"))
    events[[length(events)+1L]] <<- dependency
    list(value=saved$value,dependency=dependency)
  }
  list(get=get,ledger=function()events,
    size=function()length(ls(envir=values,all.names=TRUE)))
}

ows_write_artifact <- function(object, directory, metadata=list()) {
  stopifnot(is.character(directory),length(directory)==1L,!is.na(directory),
    nzchar(directory),is.list(metadata))
  # Exclusive directory reservation prevents overwriting an earlier success or
  # a partial failure. Failed directories are evidence, not automatic retries.
  if (file.exists(directory)||dir.exists(directory)||
      !dir.create(directory,recursive=FALSE,showWarnings=FALSE))
    stop("Artifact directory must be a new child of an existing directory")
  directory <- normalizePath(directory,mustWork=TRUE)
  pending <- file.path(directory,"artifact.pending.rds")
  path <- file.path(directory,"artifact.rds")
  sha <- function(p)digest::digest(file=p,algo="sha256")
  tryCatch({
    packed <- ows_pack_arrays(object)
    saveRDS(packed,pending,version=2L,compress="gzip")
    external_sha <- sha(pending)
    restored <- ows_unpack_arrays(readRDS(pending))
    if (!identical(object,restored,attrib.as.set=FALSE))
      stop("Lossless artifact failed strict disk roundtrip")
    stopifnot(identical(sha(pending),external_sha),file.rename(pending,path),
      identical(sha(path),external_sha))
    receipt <- list(status="ARTIFACT_COMMITTED",format=packed$format,
      file="artifact.rds",file_sha256=external_sha,file_bytes=unname(file.info(path)$size),
      exact_roundtrip=TRUE,attribute_order_preserved=TRUE,
      references=packed$references,unique_blocks=packed$unique_blocks,metadata=metadata)
    receipt_tmp <- file.path(directory,"receipt.pending.json")
    jsonlite::write_json(receipt,receipt_tmp,auto_unbox=TRUE,pretty=TRUE,digits=NA,null="null")
    stopifnot(file.rename(receipt_tmp,file.path(directory,"receipt.json")))
    # Bind this SHA to the enclosing run manifest/inventory. The internal block
    # checks alone are not authentication of a scientific output artifact.
    list(directory=directory,file_sha256=external_sha,
      receipt_sha256=sha(file.path(directory,"receipt.json")),receipt=receipt)
  },error=function(e) {
    # Preserve the unmodified object as well as any partial packed bytes.
    # A failed fallback is itself exposed in the failure receipt.
    retained <- tryCatch({saveRDS(object,file.path(directory,"original-on-failure.rds"),
      version=2L,compress="gzip"); TRUE},error=function(err)conditionMessage(err))
    failure <- list(status="ARTIFACT_WRITE_FAILED",error=conditionMessage(e),
      original_retained=retained,metadata=metadata)
    failure_receipt <- tryCatch({
      jsonlite::write_json(failure,file.path(directory,"failure.json"),
        auto_unbox=TRUE,pretty=TRUE,digits=NA,null="null"); TRUE
    },error=function(err)conditionMessage(err))
    stop(paste0("Artifact write failed: ",conditionMessage(e),
      "; original retention: ",if(isTRUE(retained))"saved" else retained,
      "; failure receipt: ",if(isTRUE(failure_receipt))"saved" else failure_receipt),call.=FALSE)
  })
}

ows_read_artifact <- function(directory, file_sha256, receipt_sha256) {
  # The expected hashes come from the enclosing run's inventory, not merely
  # from the artifact's self-reported receipt.
  stopifnot(is.character(file_sha256),length(file_sha256)==1L,
    !is.na(file_sha256),grepl("^[0-9a-f]{64}$",file_sha256),
    is.character(receipt_sha256),length(receipt_sha256)==1L,
    !is.na(receipt_sha256),grepl("^[0-9a-f]{64}$",receipt_sha256),
    !file.exists(file.path(directory,"failure.json")))
  receipt_path <- file.path(directory,"receipt.json")
  stopifnot(identical(digest::digest(file=receipt_path,algo="sha256"),receipt_sha256))
  receipt <- jsonlite::read_json(receipt_path,simplifyVector=TRUE)
  stopifnot(identical(receipt$status,"ARTIFACT_COMMITTED"),
    identical(receipt$file,"artifact.rds"),identical(receipt$file_sha256,file_sha256),
    isTRUE(receipt$exact_roundtrip),isTRUE(receipt$attribute_order_preserved))
  path <- file.path(directory,"artifact.rds")
  stopifnot(identical(digest::digest(file=path,algo="sha256"),file_sha256),
    file.info(path)$size==receipt$file_bytes)
  ows_unpack_arrays(readRDS(path))
}
