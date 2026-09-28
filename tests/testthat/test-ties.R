test_that("boundary allocation retains closer donors and spreads exact ties", {
  A <- c(rep(0, 20), rep(1, 60))
  D <- cbind(c(0, rep(1,19),rep(0,60)), 0)
  set.seed(96); before <- .Random.seed; kind <- RNGkind()
  fit <- wdsmatch:::wdsm_make_matches(A,D,D,M=3,estimand="PATT",tie_seed=2024)
  expect_identical(.Random.seed,before)
  expect_identical(RNGkind(),kind)
  choices <- fit$matches_0[A==1]
  expect_true(all(vapply(choices, function(j) length(j)==3 && !anyDuplicated(j) &&
    1L %in% j && all(A[j]==0), logical(1))))
  expect_gt(length(unique(vapply(choices,paste,collapse=",",character(1)))),1L)
  expect_gt(length(unique(unlist(choices))),3L)
  expect_identical(fit$matches_0,wdsmatch:::wdsm_make_matches(A,D,D,M=3,
    estimand="PATT",tie_seed=2024)$matches_0)
  expect_false(identical(fit$matches_0,wdsmatch:::wdsm_make_matches(A,D,D,M=3,
    estimand="PATT",tie_seed=2025)$matches_0))
  expect_identical(fit$tie_diagnostics$arms[[1]]$randomized_recipients,60L)
  expect_error(wdsmatch:::wdsm_make_matches(A,D,D,M=30,estimand="PATT"))
  expect_identical(.Random.seed,before)
})

test_that("tie helper restores a missing seed and nondefault RNG kind", {
  previous_kind <- RNGkind(); had <- exists(".Random.seed",.GlobalEnv)
  previous <- if(had) get(".Random.seed",.GlobalEnv) else NULL
  on.exit({do.call(RNGkind,as.list(previous_kind));if(had)assign(".Random.seed",previous,.GlobalEnv)
    else if(exists(".Random.seed",.GlobalEnv))rm(".Random.seed",envir=.GlobalEnv)},add=TRUE)
  RNGkind("L'Ecuyer-CMRG");set.seed(19);before <- .Random.seed;kind <- RNGkind()
  A<-rep(0:1,each=10);D<-matrix(0,20,2)
  invisible(wdsmatch:::wdsm_make_matches(A,D,D,M=5))
  expect_identical(RNGkind(),kind);expect_identical(.Random.seed,before)
  rm(".Random.seed",envir=.GlobalEnv)
  invisible(wdsmatch:::wdsm_make_matches(A,D,D,M=5))
  expect_false(exists(".Random.seed",.GlobalEnv,inherits=FALSE))
  expect_identical(RNGkind(),kind)
})

test_that("public PATE and PATT fixed-match replicates agree with a discrete oracle", {
  set.seed(3301)
  grid<-expand.grid(x1=c(-1,0,1),x2=c(-1,0,1))
  X<-grid[rep(seq_len(9),each=20),]
  A<-rep(rep(0:1,each=10),9);n<-length(A)
  Y<-X$x1+X$x2^2+A+rnorm(n);w<-exp(runif(n,-.5,.5))
  ps<-plogis(.3*X$x1+.2*X$x2)
  pg<-cbind(X$x1-.4*X$x2,X$x1+.4*X$x2)
  for(estimand in c("PATE","PATT")) for(M in c(1L,3L,5L)) {
    scores<-wdsmatch:::estimate_scores(Y,X,A,w,ps,pg,NULL,NULL,
      "prospective",estimand,use.bias.correction=FALSE)
    point<-wdsmatch:::wdsm_point(Y,A,w,scores,M,estimand,tie_seed=2024)
    target<-if(estimand=="PATE")seq_len(n) else which(A==1)
    K<-numeric(n);contrasts<-numeric(n)
    for(i in target) {
      j<-if(A[i]==1)point$matches_0[[i]] else point$matches_1[[i]]
      expect_length(j,M)
      K[j]<-K[j]+w[i]*w[j]/sum(w[j])
      contrasts[i]<-(2*A[i]-1)*(Y[i]-sum(w[j]*Y[j])/sum(w[j]))
    }
    oracle<-sum(w[target]*contrasts[target])/sum(w[target])
    set.seed(662);counts<-rmultinom(8,n,rep(1/n,n));after<-.Random.seed
    draws<-apply(counts,2,function(m)if(estimand=="PATE")
      sum(m*(2*A-1)*(w+K)*Y)/sum(m*w) else
      (sum(m[A==1]*w[A==1]*Y[A==1])-sum(m[A==0]*K[A==0]*Y[A==0]))/sum(m[A==1]*w[A==1]))
    for(design in c("prospective","retrospective")) {
      fun<-if(estimand=="PATE")wdsmatchATE else wdsmatchATT
      set.seed(662)
      result<-fun(Y,X,A,w,M=M,ps=ps,pg=pg,sampling=design,
        use.bias.correction=FALSE,boots=8,tie.seed=2024)
      expect_equal(result$estimate,oracle,tolerance=1e-12)
      expect_equal(result$boot.estimates,draws,tolerance=1e-12)
      expect_identical(.Random.seed,after)
      expect_true(result$diagnostics$ties$independent_recipient_selection)
      expect_equal(result$se,sqrt(mean((draws-mean(draws))^2)),tolerance=1e-12)
    }
  }
})

test_that("invalid tie options fail through the public interface", {
  data(survey_obs)
  args<-list(Y=survey_obs$Y,X=survey_obs[,paste0("X",1:6)],
    Z=survey_obs$Z,weights=survey_obs$survey_weight,varest=FALSE)
  for(seed in list(-1,NA,Inf,1.5,c(1,2)))
    expect_error(do.call(wdsmatchATE,c(args,list(tie.seed=seed))),"Invalid matching")
  for(tol in list(-1,NA,Inf,c(0,1)))
    expect_error(do.call(wdsmatchATT,c(args,list(tie.tolerance=tol))),"Invalid matching")
})
