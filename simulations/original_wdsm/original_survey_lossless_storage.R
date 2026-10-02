# Private lossless storage for repeated numeric arrays in saved study artifacts.
# This changes serialization only: never rounding, fitting, RNG, or field removal.
.ows_storage_plain <- function(x) {
  if (!typeof(x)%in%c("NULL","logical","integer","double","complex","character","raw","list"))
    stop("Storage supports only plain list/atomic artifacts, including their attributes")
  if (typeof(x)=="list")for(i in seq_along(x)).ows_storage_plain(x[[i]])
  for(attribute in attributes(x)).ows_storage_plain(attribute)
  invisible(NULL)
}

ows_pack_arrays <- function(object, minimum_length=128L) {
  stopifnot(length(minimum_length)==1L,is.finite(minimum_length),minimum_length>=1L,
    minimum_length==floor(minimum_length))
  .ows_storage_plain(object)
  blocks <- new.env(hash=TRUE,parent=emptyenv()); references <- 0L
  encode <- function(x) {
    if (inherits(x,".ows_atomic_reference")) stop("Reserved storage reference class in input")
    if (typeof(x)%in%c("environment","closure","externalptr","weakref"))
      stop("Storage requires data artifacts without functions or external references")
    if (is.atomic(x)&&length(x)>=minimum_length) {
      key <- digest::digest(x,algo="sha256",serializeVersion=2L)
      if (exists(key,envir=blocks,inherits=FALSE)) {
        if(!identical(get(key,envir=blocks,inherits=FALSE),x))stop("Storage hash collision")
      } else assign(key,x,envir=blocks)
      references <<- references+1L
      return(structure(list(key=key),class=".ows_atomic_reference"))
    }
    if (is.list(x)) {
      y <- lapply(x,encode); attributes(y) <- attributes(x); return(y)
    }
    x
  }
  payload <- encode(object)
  # R objects can be identical despite different compact row-name or ALTREP
  # representations. Authenticate the stored payload bytes, not a presumed
  # canonical serialization of the reconstructed object.
  payload_blob <- serialize(payload,NULL,version=2L)
  keys <- sort(ls(blocks,all.names=TRUE))
  structure(list(format="ows_lossless_atomic_arrays_v2",payload_blob=payload_blob,
    blocks=mget(keys,envir=blocks,inherits=FALSE),references=references,
    unique_blocks=length(keys),minimum_length=minimum_length,
    payload_sha256=digest::digest(payload_blob,algo="sha256",serialize=FALSE)),class="ows_packed_arrays")
}

ows_unpack_arrays <- function(object) {
  .ows_storage_plain(object)
  stopifnot(inherits(object,"ows_packed_arrays"),
    identical(object$format,"ows_lossless_atomic_arrays_v2"),is.list(object$blocks),
    length(object$blocks)==object$unique_blocks,!anyDuplicated(names(object$blocks)),
    is.raw(object$payload_blob),
    identical(digest::digest(object$payload_blob,algo="sha256",serialize=FALSE),object$payload_sha256))
  for(key in names(object$blocks))stopifnot(identical(digest::digest(object$blocks[[key]],algo="sha256",serializeVersion=2L),key))
  decode <- function(x) {
    if (inherits(x,".ows_atomic_reference")) {
      stopifnot(identical(names(x),"key"),is.character(x$key),length(x$key)==1L,
        !is.na(x$key),x$key%in%names(object$blocks))
      return(object$blocks[[x$key]])
    }
    if(is.list(x)) {
      y <- lapply(x,decode); attributes(y) <- attributes(x); return(y)
    }
    x
  }
  payload <- unserialize(object$payload_blob)
  .ows_storage_plain(payload)
  decode(payload)
}
