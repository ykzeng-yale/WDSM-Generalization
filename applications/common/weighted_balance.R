# Pure saved-graph diagnostic functions. Sourcing executes no statistic or RNG.
wb_moments <- function(X,w) {
  stopifnot(is.matrix(X),is.numeric(X),length(w)==nrow(X),all(is.finite(X)),
    all(is.finite(w)&w>=0),sum(w)>0)
  # Normalized probabilities give scale-invariant means/covariance.
  p <- (w/max(w)); if(any(w>0&p==0))stop("Weight normalization underflow")
  p <- p/sum(p); mu <- colSums(X*p)
  positive <- which(w>0)
  constant <- vapply(seq_len(ncol(X)),function(j)all(X[positive,j]==X[positive[1L],j]),logical(1))
  mu[constant] <- X[positive[1L],constant]
  C <- sweep(X,2L,mu,"-")
  population_covariance <- crossprod(C*sqrt(p))
  df <- 1-sum(p^2)
  covariance <- if(df>0)population_covariance/df else matrix(NA_real_,ncol(X),ncol(X),dimnames=dimnames(population_covariance))
  list(mean=mu,covariance=covariance,population_covariance=population_covariance,
    weight_total=sum(w),positive_units=sum(w>0),ESS=1/sum(p^2),finite_weight_df=df)
}
wb_reference <- function(X,Z,W,binary,estimand) {
  stopifnot(estimand %in% c("PATE","PATT"),length(Z)==nrow(X),length(W)==nrow(X),
    all(Z %in% 0:1),all(W>0&is.finite(W)),length(binary)==ncol(X),is.logical(binary))
  if(any(binary))stopifnot(all(X[,binary,drop=FALSE] %in% 0:1))
  before <- lapply(0:1,function(z)wb_moments(X[Z==z,,drop=FALSE],W[Z==z]))
  names(before) <- c("control","treated")
  variances <- lapply(before,function(v) {
    s <- diag(v$covariance);names(s)<-colnames(X)
    # Binary indicator standardization follows documented Bernoulli variance.
    s[binary]<-v$mean[binary]*(1-v$mean[binary]);s
  })
  scale <- sqrt(if(estimand=="PATE")(variances[[1L]]+variances[[2L]])/2 else variances[[2L]])
  target <- wb_moments(X,if(estimand=="PATE")W else W*Z)$mean
  list(estimand=estimand,before=before,fixed_SD=scale,binary=binary,target_mean=target,
    mean_difference_before=before$treated$mean-before$control$mean,
    denominator=if(estimand=="PATE")"pre_match_W_weighted_pooled_arm_SD"else"pre_match_W_weighted_treated_SD")
}
wb_diagnostics <- function(data,schema) {
  # The manifest defines diagnostic columns only; it never changes metric X.
  X <- matrix(numeric(),nrow(data),0L);map <- list();binary <- logical()
  for(i in seq_along(schema)) {
    item <- schema[[i]];v <- data[[item$source_variable]]
    if(is.null(v)||!is.numeric(v)||any(!is.finite(v)))stop("Invalid diagnostic source: ",item$source_variable)
    if(item$kind=="binary_indicator") {
      allowed <- as.numeric(unlist(item$allowed_codes,use.names=FALSE))
      if(any(!v %in% allowed))stop("Unexpected category retained as failure: ",item$source_variable)
      value <- as.numeric(v==as.numeric(item$category));is_binary<-TRUE
    } else if(item$kind=="binary_numeric") {
      stopifnot(all(v %in% 0:1));value<-v;is_binary<-TRUE
    } else {
      stopifnot(item$kind=="continuous");value<-v;is_binary<-FALSE
    }
    X<-cbind(X,value);colnames(X)[ncol(X)]<-item$diagnostic
    map[[i]]<-data.frame(diagnostic=item$diagnostic,source_variable=item$source_variable,
      kind=item$kind,category=if(is.null(item$category))NA_character_ else as.character(item$category),label=item$label)
    binary<-c(binary,is_binary)
  }
  stopifnot(!anyDuplicated(colnames(X)))
  list(X=X,binary=binary,mapping=do.call(rbind,map),scope="Manifest diagnostics only; original matching metric unchanged")
}
wb_saved_graph <- function(point,X,Z,W,reference) {
  stopifnot(inherits(point,"wm_match"),identical(point$method,"self_normalized"),
    identical(point$estimand,reference$estimand),point$n==nrow(X),
    identical(as.integer(point$data$Z),as.integer(Z)),identical(as.numeric(point$weights),as.numeric(W)),
    is.finite(point$weight_scale),point$weight_scale>0,
    identical(as.numeric(point$analysis_weights),as.numeric(W/point$weight_scale)),
    is.data.frame(point$graph$edges),nrow(point$graph$edges)>0,
    is.list(point$graph$neighbors),length(point$graph$neighbors)==nrow(X))
  K <- point$loads$incoming
  stopifnot(is.matrix(K),identical(dim(K),c(nrow(X),2L)),
    identical(colnames(K),c("arm0","arm1")),all(is.finite(K)&K>=0))
  # Incoming stored in W/max(W) units; return to the supplied known-W units.
  K <- K*point$weight_scale
  stopifnot(all(K[Z==1L,1L]==0),all(K[Z==0L,2L]==0))
  if(reference$estimand=="PATT")stopifnot(all(K[,2L]==0))
  own <- K[cbind(seq_len(nrow(X)),as.integer(Z)+1L)]
  after_w <- if(reference$estimand=="PATE")W+own else ifelse(Z==1L,W,K[,1L])
  post <- lapply(0:1,function(z)wb_moments(X[Z==z,,drop=FALSE],after_w[Z==z]))
  names(post)<-c("control","treated")
  diff <- post$treated$mean-post$control$mean
  sd <- reference$fixed_SD
  valid <- is.finite(sd)&sd>0
  status <- ifelse(!is.finite(sd),"undefined_pre_match_SD",
    ifelse(sd==0,"zero_pre_match_SD_SMD_undefined","finite_fixed_scale"))
  standardize <- function(delta) {out<-rep(NA_real_,length(delta));out[valid]<-delta[valid]/sd[valid];out}
  before_smd<-standardize(reference$mean_difference_before);after_smd<-standardize(diff)
  expected_mass <- sum(if(reference$estimand=="PATE")W else W*Z)
  observed_mass <- vapply(post,`[[`,numeric(1),"weight_total")
  table <- data.frame(variable=colnames(X),kind=ifelse(reference$binary,"binary_indicator","continuous"),
    fixed_pre_match_SD=sd,SD_status=status,
    mean_control_before=reference$before$control$mean,mean_treated_before=reference$before$treated$mean,
    difference_before=reference$mean_difference_before,SMD_before=before_smd,abs_SMD_before=abs(before_smd),
    mean_control_after=post$control$mean,mean_treated_after=post$treated$mean,
    difference_after=diff,SMD_after=after_smd,abs_SMD_after=abs(after_smd),target_mean=reference$target_mean,
    deviation_control_after_from_target=post$control$mean-reference$target_mean,
    deviation_treated_after_from_target=post$treated$mean-reference$target_mean,
    estimand=reference$estimand,denominator=reference$denominator,
    before_label="same_original_sample_W_weighted_arm_means",
    after_label=if(reference$estimand=="PATE")"saved_graph_W_plus_K_each_arm"else"saved_graph_treated_W_control_K0",
    stringsAsFactors=FALSE)
  list(table=table,before=reference$before,after=post,
    mass=list(expected_each_arm=expected_mass,actual_each_arm=observed_mass,
      residual_each_arm=observed_mass-expected_mass,
      relative_residual_each_arm=(observed_mass-expected_mass)/expected_mass,
      mathematical_identity="PATE eacharm=sumW;PATT eacharm=treatedW",
      scope="Retained arithmetic integrity evidence; no tolerance-based scientific acceptance"),
    method=point$method,M=point$M,estimand=point$estimand,
    scope="Descriptive saved raw matching distribution; outcome-model correction may affect effects, not these graph weights",
    scientific_acceptance=FALSE)
}
