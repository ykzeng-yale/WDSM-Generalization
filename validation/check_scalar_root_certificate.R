# Deterministic certificate checks on one small fit and exact constructed data.
source('R/utils.R')
source('R/wm_core.R')
source('R/wm_scalar_logistic_match.R')
source('R/wm_scalar_root_certificate.R')
had_rng <- exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE)
if (had_rng) rng <- get('.Random.seed',envir=.GlobalEnv)
start <- proc.time()[['elapsed']]
# One deterministic 24-row fit; no random generation or Monte Carlo.
T <- seq(-1.1,1.1,length.out=24); V <- rep(c(-.8,.2,.6,-.3),6)
L <- cbind(intercept=1,index=T,auxiliary=V); Z <- rep(c(0,1,1,0,1,0),4)
Y <- .5+.8*Z+T^2-.3*V+rep(c(-.2,.1,.3),8)
W <- ifelse(Z==1,1.5+.1*plogis(T),2+.2*plogis(T))
saved <- list(fit=wm_scalar_logistic_match(Y,Z,W,L,M=3,inference=FALSE),design=L)
checks <- character()
check <- function(value, name) {
  if (!isTRUE(value)) stop(name)
  checks <<- c(checks, name)
}
fails <- function(expr, name) check(inherits(try(expr,silent=TRUE),'try-error'), name)
A <- .wm_cert_arithmetic(64L); Q <- A$Q
x <- Q(c(-17,-8,-1,0,1,8,17),c(3,3,3,1,3,3,3))
for (bits in c(2L,64L,128L)) {
  a <- .wm_cert_arithmetic(bits)
  check(all(a$roundq(x,FALSE)<=x & a$roundq(x,TRUE)>=x),paste('rational outward',bits))
  check(all(a$roundq(-x,FALSE)==-a$roundq(x,TRUE)),paste('negative duality',bits))
}
binary <- Q(c(0,.Machine$double.xmin,2^-1074,.Machine$double.xmax,-.Machine$double.xmax))
check(all(A$roundq(binary,FALSE)==binary & A$roundq(binary,TRUE)==binary),'exact binary endpoints')
iv <- A$pack(x); printed <- A$record(iv)
check(all(Q(printed$lower)<=iv$lo & Q(printed$upper)>=iv$hi),'directed double enclosure')
check(all(A$pack(x)$lo <= x),'lower enclosures')
check(all(A$mul(A$exact(x),A$exact(rev(x)))$lo<=x*rev(x)),'signed interval multiplication lower')
check(all(A$mul(A$exact(x),A$exact(rev(x)))$hi>=x*rev(x)),'signed interval multiplication upper')
check(all(A$total(A$pack(x))$lo<=sum(x) & A$total(A$pack(x))$hi>=sum(x)),'signed interval sum')
fails(A$div(A$exact(1),A$pack(Q(-1),Q(1))),'zero denominator rejected')
fails(A$exp(A$exact(1001)),'exponential argument guard')
# Independent exact Taylor-series reference: exp(t), 0<=t<=1.
exp_ref <- function(t) {
  term <- total <- Q(1)
  for (k in seq_len(100L)) { term <- term*t/k; total <- total+term }
  next_term <- term*t/101
  list(lo=total,hi=total+next_term/(1-t/102))
}
B <- .wm_cert_arithmetic(128L)
for (t in list(Q(0),Q(1,2),Q(1))) {
  ref <- exp_ref(t)
  e <- B$exp(B$exact(t))
  check(e$lo<=ref$lo && e$hi>=ref$hi,paste('Taylor positive exp enclosure',as.character(t)))
  en <- B$exp(B$exact(-t))
  check(en$lo<=1/ref$hi && en$hi>=1/ref$lo,paste('Taylor negative exp enclosure',as.character(t)))
}
f <- saved$fit; L <- saved$design
c0 <- wm_scalar_root_certificate(f,L)
check(identical(c0$status,'root_and_donor_sets_certified'),'saved fit root and graph')
check(c0$point_root_error_bound$upper < 1e-14,'point arithmetic bound')
check(c0$stored_probability_error_upper$upper < 1e-10,'stored score error bound')
check(!c0$population_assumptions_verified && !c0$asymptotic_success_rate_claimed &&
      !c0$variance_arithmetic_certified && !c0$bootstrap_validity_claimed,'scope nonclaims')
check(identical(c0$inputs$weights,f$fit$weights) && identical(c0$inputs$probability,f$propensity$probability),
      'actual input binding')
# Independent exact raw-point evaluation; all fractions retained, no interval helpers.
exact_point <- function(fit) {
  w <- Q(fit$weights); y <- Q(fit$data$Y); z <- fit$data$Z
  query <- if (fit$estimand=='PATE') seq_along(z) else which(z==1)
  num <- Q(0)
  for (i in query) {
    j <- fit$graph$neighbors[[i]]
    num <- num+w[i]*(2*z[i]-1)*(y[i]-sum(w[j]*y[j])/sum(w[j]))
  }
  num/sum(w[query])
}
truth <- exact_point(f$fit)
check(Q(c0$exact_current_graph_point$lower_exact)<=truth &&
      Q(c0$exact_current_graph_point$upper_exact)>=truth,'independent exact rational PATE')
check(Q(c0$point_root_error_bound$upper_exact)>=abs(Q(f$estimate)-truth),'independent arithmetic discrepancy')
# Reuse stored probabilities: these calls construct graphs, never refit models.
make_object <- function(Y,Z,W,L,e,M=3L,estimand='PATE') {
  fit <- wm_match(Y,Z,W,matrix(e,ncol=1L),M=M,estimand=estimand,variance=FALSE)
  structure(list(status='point_computed',fit=fit,estimate=fit$estimate,
    propensity=list(probability=e)),class=c('wm_scalar_logistic_match','list'))
}
att <- make_object(f$fit$data$Y,f$fit$data$Z,f$fit$weights,L,f$propensity$probability,estimand='PATT')
ca <- wm_scalar_root_certificate(att,L)
check(ca$donor_sets_certified,'PATT graph certificate')
ta <- exact_point(att$fit)
check(Q(ca$exact_current_graph_point$lower_exact)<=ta && Q(ca$exact_current_graph_point$upper_exact)>=ta,
      'independent exact rational PATT')
check(identical(ca$inputs$estimand,'PATT'),'PATT label')
bad_center <- wm_scalar_root_certificate(f,L,center=c(10,0,0),radius=1e-8)
check(!bad_center$root_contained && identical(bad_center$status,'not_certified'),'far center inconclusive')
small_ball <- wm_scalar_root_certificate(f,L,radius=1e-20)
check(!small_ball$root_contained,'insufficient root ball inconclusive')
extreme <- wm_scalar_root_certificate(f,L,center=c(10000,0,0))
check(!extreme$root_contained && !is.null(extreme$reason),'saturated proposal retained as failure')
bad <- f; bad$fit$graph$neighbors[[1]][1] <- 1L
fails(wm_scalar_root_certificate(bad,L),'same-arm donor rejected')
bad <- f; bad$fit$graph$edges$donor[1] <- 1L
fails(wm_scalar_root_certificate(bad,L),'edge-neighbor discrepancy rejected')
for (nm in c('query','donor','arm')) {
  bad <- f; bad$fit$graph$edges[[nm]][1] <- bad$fit$graph$edges[[nm]][1] + .25
  fails(wm_scalar_root_certificate(bad,L),paste('fractional edge rejected',nm))
}
bad <- f; bad$propensity$probability[1] <- bad$propensity$probability[1]+.001
fails(wm_scalar_root_certificate(bad,L),'different stored score rejected')
bad <- f; bad$fit$info$corrected <- TRUE
fails(wm_scalar_root_certificate(bad,L),'corrected fit rejected')
fails(wm_scalar_root_certificate(f,L,precision=63),'invalid precision rejected')
fails(wm_scalar_root_certificate(f,L,radius=0),'invalid radius rejected')
# Exact balanced data have beta*=0. Boundary ties remain inconclusive.
v <- rep(-3:3,2); lt <- cbind(1,v,v^2); zt <- rep(0:1,each=7)
tie <- make_object(v+zt,zt,rep(1,14),lt,rep(.5,14),M=1L)
ct <- wm_scalar_root_certificate(tie,lt,center=c(0,0,0))
check(ct$root_contained && !ct$donor_sets_certified,'boundary ties explicitly unresolved')
all <- make_object(v+zt,zt,rep(1,14),lt,rep(.5,14),M=7L)
call <- wm_scalar_root_certificate(all,lt,center=c(0,0,0))
check(call$root_contained && call$donor_sets_certified && is.null(call$minimum_stored_gap),
      'all donors certify without strict gap')
check(identical(f,saved$fit),'input object unchanged')
check(identical(had_rng,exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE)) &&
      (!had_rng || identical(rng,get('.Random.seed',envir=.GlobalEnv))),'RNG unchanged')
result <- list(status='PASS',checks=checks,count=length(checks),elapsed=proc.time()[['elapsed']]-start,
  new_fits=1L,random_draws=0L,certificates=list(PATE=c0,PATT=ca,ties=ct,all_donors=call))
cat(length(checks),'checks PASS; elapsed',result$elapsed,'seconds;1 small fit;0 random draws\n')
