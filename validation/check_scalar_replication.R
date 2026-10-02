# One saved 24-row fixture, no fits, kernel regression, RNG or new graph.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1L] else ".", mustWork = TRUE)
source(file.path(root, "R", "utils.R"))
source(file.path(root, "R", "wm_core.R"))
source(file.path(root, "R", "wm_scalar_logistic_match.R"))
source(file.path(root, "R", "wm_scalar_replication.R"))
checks <- character()
check <- function(ok, label) {
  if (!isTRUE(ok)) stop(label)
  checks <<- c(checks, label)
}
near <- function(x, y, label, tol = 2e-10) {
  check(identical(dim(x), dim(y)) && length(x) == length(y) &&
    all(is.finite(c(x, y))) && all(abs(x - y) <= tol * pmax(1, abs(x), abs(y))), label)
}
fails <- function(expr, label) check(inherits(try(expr, silent = TRUE), "try-error"), label)
had_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
if (had_rng) rng_before <- get(".Random.seed", envir = .GlobalEnv)
forbidden <- function(...) stop("Forbidden statistical fitting or RNG call.")
wdsm_fit_ps <- wm_match <- .wm_scalar_nw <- forbidden
rnorm <- runif <- sample <- rmultinom <- forbidden
start <- proc.time()[["elapsed"]]
fixture <- readRDS(file.path(root, "validation", "fixtures", "scalar_replication_24.rds"))
saved <- fixture$fits
design <- fixture$design
check(isTRUE(fixture$provenance$synthetic) &&
  identical(dim(design), c(24L, 3L)) &&
  identical(fixture$provenance$original_fits_sha256,
    "58485ec76a4cc7c5427a5f49df29552d22476c3b9cd5cbd21de8c297ccc54336") &&
  identical(design, saved$PATE$fitting_design) &&
  identical(design, saved$PATT$fitting_design), "portable synthetic fixture provenance")
for (target in c("PATE", "PATT")) {
  x <- saved[[target]]
  old <- x
  old$fitting_design <- old$design_binding <- NULL
  old$inference$model_contract <- NULL
  unbound <- wm_scalar_replication(old, design)
  check(unbound$status == "unavailable" && unbound$design_binding == "unbound" &&
    identical(unbound$fit, old$fit), paste(target, "unbound legacy preserves raw fit"))
  check(identical(old$fit$data, saved$PATE$fit$data) &&
    identical(old$fit$weights, saved$PATE$fit$weights) &&
    identical(old$propensity$probability, saved$PATE$propensity$probability) &&
    identical(old$design_columns, colnames(design)), paste(target, "original design provenance"))
  # The fixture retains its original numerical values. Removing only the
  # documented metadata above tests the unavailable legacy-object branch.
  original <- serialize(x, NULL)
  a <- wm_scalar_replication(x)
  check(a$status == "row_variance_computed", paste(target, "computed row variance"))
  check(identical(a$fit, x$fit) && identical(a$source_inference, x$inference) &&
    identical(serialize(x, NULL), original), paste(target, "original object fully preserved"))
  check(identical(a$numerical_root_bridge, x$numerical_root_bridge) &&
    !a$population_assumptions_verified && !a$supplied_draw_law_verified &&
    !a$original_refit_bootstrap, paste(target, "qualified scope retained"))
  n <- x$fit$n; Y <- x$fit$data$Y; Z <- x$fit$data$Z
  w <- x$fit$analysis_weights; e <- x$propensity$probability
  mu <- x$inference$predictions$prediction; p <- ncol(design)
  gamma <- if (target == "PATE") sum(w)/n else sum(Z*w)/n
  K <- matrix(0, n, 2L)
  query <- if (target == "PATE") seq_len(n) else which(Z == 1)
  corrected_numerator <- 0
  for (i in query) {
    donors <- x$fit$graph$neighbors[[i]]; z <- 1-Z[i]
    share <- w[donors]/sum(w[donors])
    K[donors,z+1] <- K[donors,z+1] + w[i]*share
    impute <- mu[i,z+1] + sum(share*(Y[donors]-mu[donors,z+1]))
    corrected_numerator <- corrected_numerator + w[i]*(2*Z[i]-1)*(Y[i]-impute)
  }
  near(unname(K), unname(x$fit$loads$incoming), paste(target, "direct weighted incoming reconstruction"))
  row <- vapply(seq_len(n), function(i) {
    if (target == "PATE") {
      (w[i]*(mu[i,2]-mu[i,1]-x$estimate) +
         (2*Z[i]-1)*(w[i]+K[i,Z[i]+1])*(Y[i]-mu[i,Z[i]+1]))/gamma
    } else {
      (Z[i]*w[i]*(Y[i]-mu[i,1]-x$estimate) -
         (1-Z[i])*K[i,1]*(Y[i]-mu[i,1]))/gamma
    }
  }, numeric(1))
  near(row, a$rows$contribution, paste(target, "independent complete row"))
  near(mean(row), corrected_numerator/(n*gamma)-x$estimate,
       paste(target, "row mean preserves raw point rather than re-centering point"))
  J <- Reduce("+", lapply(seq_len(n), function(i)
    w[i]*e[i]*(1-e[i])*tcrossprod(design[i,])))/n
  phi <- t(vapply(seq_len(n), function(i)
    solve(J, w[i]*design[i,]*(Z[i]-e[i])), numeric(p)))
  near(unname(phi), unname(a$rows$influence), paste(target, "independent full influence"))
  near(unname(crossprod(phi)/n), unname(x$inference$blocks$Sigma),
       paste(target, "uncentered influence sandwich"))
  phic <- sweep(phi,2,colMeans(phi),"-"); Lc <- row-mean(row)
  near(unname(crossprod(phic)/n),
       unname(x$inference$blocks$Sigma-tcrossprod(colMeans(phi))),
       paste(target, "centered sandwich subtraction"))
  b <- unname(x$inference$blocks$b)
  aug <- vapply(seq_len(n), function(i) Lc[i]+sum(b*phic[i,]), numeric(1))
  aug <- aug-mean(aug)
  near(aug,a$rows$augmented,paste(target,"plus full signed b"))
  near(a$root_n_variance, mean(Lc^2)+2*sum(b*colMeans(phic*Lc))+
    drop(crossprod(b,crossprod(phic)%*%b))/n,paste(target,"complete centered variance identity"))
  near(a$conditional_point_variance,a$root_n_variance/n,paste(target,"full-n point variance"))
  check(!a$comparison$finite_sample_equality_expected &&
    identical(a$comparison$explicit_root_n_variance,x$inference$root_n_variance),
    paste(target,"separate explicit variance kept"))
  i <- which(Z==0)[1]; j <- which(Z==1)[1]
  count <- matrix(1,n,4); count[i,2] <- 2; count[j,2] <- 0
  count[i,3] <- 0; count[j,3] <- 2; count[,4] <- 0; count[i,4] <- n
  cfit <- wm_scalar_replication(x, counts=count)
  check(cfit$status=="row_variance_computed",paste(target,"collapsed single-arm count admitted"))
  expected <- c(0,(aug[i]-aug[j])/sqrt(n),(aug[j]-aug[i])/sqrt(n),sqrt(n)*aug[i])
  near(cfit$root_n_draws,expected,paste(target,"direct finite-count roots"))
  near(cfit$draws,x$estimate+expected/sqrt(n),paste(target,"direct full-n point draws"))
  near(drop(crossprod(aug,(diag(n)-matrix(1/n,n,n))%*%aug))/n,
    a$root_n_variance,paste(target,"exact multinomial covariance algebra"))
  near(cfit$empirical_draw_variance_B,mean((cfit$draws-mean(cfit$draws))^2),
    paste(target,"finite-draw divisor B only diagnostic"))
  one <- wm_scalar_replication(x,counts=count[,1,drop=FALSE])
  check(one$B==1L && is.na(one$empirical_draw_variance_B) &&
    one$root_n_variance==a$root_n_variance,paste(target,"one column retains analytic variance"))
  mult <- cbind(rep(0,n),rep(1,n),diag(n)[,i],seq_len(n)/n)
  gfit <- wm_scalar_replication(x,multipliers=mult)
  near(gfit$root_n_draws,drop(crossprod(aug,mult))/sqrt(n),
    paste(target,"supplied multiplier algebra without randomness"))
  negative <- x; negative$inference$root_n_variance <- -1
  negative$inference$status <- "nonpositive_root_variance"
  nf <- wm_scalar_replication(negative)
  near(nf$root_n_variance,a$root_n_variance,paste(target,"nonpositive explicit variance does not select row algorithm"))
  scaled <- x; scaled$fit$weights <- 7*scaled$fit$weights
  scaled$fit$weight_scale <- 7*scaled$fit$weight_scale
  scaled$inference$moments$weight_scale <- 7*scaled$inference$moments$weight_scale
  sf <- wm_scalar_replication(scaled)
  near(sf$rows$augmented,a$rows$augmented,paste(target,"coherent common weight-scale invariance"))
  noinf <- x; noinf$inference <- list(status="not_requested")
  check(wm_scalar_replication(noinf)$status=="unavailable",paste(target,"no implicit missing-mean fit"))
  wrong <- design; wrong[1,2] <- wrong[1,2]+1
  fails(wm_scalar_replication(x,wrong),paste(target,"changed design rejected"))
  fails(wm_scalar_replication(x,design[n:1,,drop=FALSE]),paste(target,"row order rejected"))
  wrong <- design; colnames(wrong)[2] <- "changed"
  fails(wm_scalar_replication(x,wrong),paste(target,"column name rejected"))
  wrongfit <- x; wrongfit$fit$loads$incoming[1,1] <- wrongfit$fit$loads$incoming[1,1]+1
  fails(wm_scalar_replication(wrongfit),paste(target,"changed load rejected"))
  wrongfit <- x; wrongfit$inference$blocks$J[1,1] <- 2*wrongfit$inference$blocks$J[1,1]
  fails(wm_scalar_replication(wrongfit),paste(target,"changed J rejected"))
  fails(wm_scalar_replication(x,counts=count,multipliers=mult),paste(target,"ambiguous draw inputs rejected"))
  bad <- count; bad[1,1] <- -1
  fails(wm_scalar_replication(x,counts=bad),paste(target,"negative counts rejected"))
  bad <- count; bad[1,1] <- 0.5
  fails(wm_scalar_replication(x,counts=bad),paste(target,"fractional counts rejected"))
  bad <- count; bad[1,1] <- 2
  fails(wm_scalar_replication(x,counts=bad),paste(target,"wrong count sum rejected"))
  fails(wm_scalar_replication(x,counts=count[-1,,drop=FALSE]),paste(target,"wrong n rejected"))
  fails(wm_scalar_replication(x,multipliers=mult+1i),paste(target,"complex draws rejected"))
  bad <- mult; bad[1,1] <- Inf
  fails(wm_scalar_replication(x,multipliers=bad),paste(target,"nonfinite draws rejected"))
  if (target == "PATT") {
    changed <- x; changed$inference$predictions$prediction[,2] <- 1000
    near(wm_scalar_replication(changed)$rows$contribution,a$rows$contribution,"PATT direct row does not use treated prediction")
  }
  check(identical(serialize(x,NULL),original),paste(target,"all calls preserved original"))
}
check(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)==had_rng &&
  (!had_rng || identical(rng_before,get(".Random.seed",envir=.GlobalEnv))),"RNG unchanged")
parsed <- tools::parse_Rd(file.path(root,"man","wm_scalar_replication.Rd"))
check(length(parsed)>0,"Rd parsed")
receipt <- list(status="PASS",checks=length(checks),labels=checks,
  elapsed_seconds=proc.time()[["elapsed"]]-start,statistical_fits=0,
  kernel_regressions=0,new_datasets=0,RNG_used=FALSE,
  saved_fixture_rows=24,saved_estimands=c("PATE","PATT"),
  source_data="Portable saved synthetic fixture; original numerical values with documented input-binding metadata",
  scope="Finite-data algebra and guards only; no asymptotic validation")
if (length(args) >= 2L) saveRDS(list(receipt = receipt, last_result = a), args[2L])
cat(length(checks),"checks PASS; zero fits, kernel regressions, RNG or datasets.\n")
