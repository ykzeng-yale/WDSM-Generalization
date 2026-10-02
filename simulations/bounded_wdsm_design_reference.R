# Source-form bounded iid selection design. Not the historical clustered DGP.
# No simulation is run when this file is sourced.

.wm_bd_constants <- function(index) {
  index <- match.arg(index, c("source_good", "source_poor"))
  list(a = log(c(1.10,1.25,1.50,1.75,2,2.50)),
       b = c(2.5,-2,1.75,-1.25,1.5,1.1), a7 = log(1.10),
       a0 = if (index == "source_good") log(35/80) else log(20/80),
       delta = if (index == "source_good") .6 else 2,
       c0 = log(.005/.995),
       selection = log(c(1.05,1.10,1.15,1.10,1.05,1.10)),
       treatment_selection = log(.9), index = index)
}

# Deterministic Simpson truth calculation with an analytic truncation-free
# integration error bound. The reported bound excludes floating-point roundoff.
wm_bounded_wdsm_truth_reference <- function(index = "source_good", radius = .25,
                                           tolerance = 1e-11,
                                           max_intervals = 65536L) {
  k <- .wm_bd_constants(index)
  if (!is.numeric(radius) || length(radius) != 1L || !is.finite(radius) ||
      radius < 0 || radius >= .5 || !is.numeric(tolerance) ||
      length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0 ||
      !is.numeric(max_intervals) || length(max_intervals) != 1L ||
      !is.finite(max_intervals) || max_intervals < 2 || max_intervals != floor(max_intervals)) {
    stop("Invalid radius or integration controls.")
  }
  kap <- sqrt(1 + k$a7^2)
  center <- k$a0 + k$delta*k$a7
  width <- k$delta*radius*kap
  B <- abs(width)
  # Bounds for theta derivatives of expit(center+width*cos(theta)).
  # Polynomial absolute-coefficient sums bound expit's first 4 derivatives.
  derivative <- c(1, 2*B, 6*B^2+2*B,
    26*B^3+18*B^2+2*B, 150*B^4+156*B^3+42*B^2+2*B)
  rho_factor <- c(1,1,2,4,8)
  first_factor <- (1+3^(0:4))/4
  fourth <- (2/pi)*c(
    sum(choose(4,0:4)*rev(rho_factor)*derivative),
    sum(choose(4,0:4)*rev(first_factor)*derivative))
  intervals <- 2L
  repeat {
    h <- pi/intervals
    error <- pi*h^4*fourth/180
    if (max(error) <= tolerance) break
    intervals <- 2*intervals
    if (intervals > max_intervals) stop("Truth integral exceeds interval budget.")
  }
  theta <- (0:intervals)*h
  simpson <- rep(c(2,4), length.out=intervals+1)
  simpson[c(1,intervals+1)] <- 1
  common <- (2/pi)*sin(theta)^2*stats::plogis(center+width*cos(theta))
  rho <- sum(simpson*common)*h/3
  first <- sum(simpson*common*cos(theta))*h/3
  lower <- stats::plogis(center-B)
  multiplier <- .3*radius*k$a7/kap
  patt <- 1.3+multiplier*first/rho
  patt_error <- abs(multiplier)*(error[2]/lower + abs(first)*error[1]/(lower*rho))
  list(PATE=1.3, PATT=patt, treatment_probability=rho,
       projected_first_moment=first, integration_error_bound=error,
       PATT_integration_error_bound=patt_error, intervals=intervals,
       radius=radius,index=k$index,roundoff_included=FALSE,
       target_law="iid population P before selection; PATT conditions on population Z=1")
}

# Exact covariate reconstruction, also useful for deterministic law checks.
.wm_bd_covariates <- function(indices, sign, radial, extras, constants) {
  n <- nrow(indices)
  X <- matrix(0,n,6)
  X[,1] <- sign*radial
  X[,2] <- sign*indices[,3]/radial
  X[,5:6] <- extras
  a <- constants$a; b <- constants$b
  rhs <- cbind(indices[,1]-as.vector(X[,c(1,2,5,6),drop=FALSE] %*% a[c(1,2,5,6)]),
               indices[,2]-as.vector(X[,c(1,2,5,6),drop=FALSE] %*% b[c(1,2,5,6)]))
  X[,3:4] <- t(solve(rbind(a[3:4],b[3:4]), t(rhs)))
  colnames(X) <- paste0("X",1:6)
  X
}

wm_bounded_wdsm_design_reference <- function(n, index = "source_good", radius = .25,
                                            max_candidates = 10000000L) {
  k <- .wm_bd_constants(index)
  if (!is.numeric(n) || length(n)!=1L || !is.finite(n) || n<2 || n!=floor(n) ||
      !is.numeric(radius) || length(radius)!=1L || !is.finite(radius) || radius<=0 || radius>=.5 ||
      !is.numeric(max_candidates) || length(max_candidates)!=1L ||
      !is.finite(max_candidates) || max_candidates<n || max_candidates!=floor(max_candidates)) {
    stop("Invalid sample size, fixed radius, or candidate budget.")
  }
  # Interval bounds on the exact reconstruction give a deterministic rejection
  # envelope. Dividing pi by this constant changes efficiency, not Q or W.
  bound <- c(2, 1+radius, 0, 0, .5, .5)
  rhs_bound <- c(radius+sum(abs(k$a[c(1,2,5,6)])*bound[c(1,2,5,6)]),
                 radius+sum(abs(k$b[c(1,2,5,6)])*bound[c(1,2,5,6)]))
  bound[3:4] <- as.vector(abs(solve(rbind(k$a[3:4],k$b[3:4]))) %*% rhs_bound)
  projection <- solve(t(rbind(k$a[3:4],k$b[3:4])),k$selection[3:4])
  remainder <- k$selection-projection[1]*k$a-projection[2]*k$b
  selection_index_bound <- min(sum(abs(k$selection)*bound),
    radius*sum(abs(projection))+sum(abs(remainder)*bound))
  pi_upper <- stats::plogis(k$c0+selection_index_bound+max(0,k$treatment_selection))
  kept_X <- matrix(numeric(),0,6)
  kept_Z <- kept_e <- kept_pi <- numeric()
  generated <- 0
  while (nrow(kept_X)<n) {
    size <- min(max(256L, 2*(n-nrow(kept_X))),max_candidates-generated)
    if (size<1) stop("Candidate budget exhausted before completing iid selection.")
    normal <- matrix(stats::rnorm(size*4),size,4)
    norms <- sqrt(rowSums(normal^2))
    if (any(!is.finite(norms)) || any(norms==0)) stop("Invalid spherical proposal.")
    U <- radius*normal[,1:3,drop=FALSE]/norms
    U[,3] <- U[,3]+1
    sign <- ifelse(stats::runif(size)<.5,-1,1)
    radial <- stats::runif(size,1,2)
    extras <- matrix(stats::runif(size*2,-.5,.5),size,2)
    X <- .wm_bd_covariates(U,sign,radial,extras,k)
    e <- stats::plogis(k$a0+k$delta*(U[,1]+k$a7*U[,3]))
    Z <- stats::rbinom(size,1,e)
    pi_select <- stats::plogis(k$c0+k$treatment_selection*Z+as.vector(X %*% k$selection))
    if (any(pi_select>pi_upper*(1+1e-12))) stop("Selection envelope violated.")
    accepted <- which(stats::runif(size)<pi_select/pi_upper)
    keep <- head(accepted,n-nrow(kept_X))
    kept_X <- rbind(kept_X,X[keep,,drop=FALSE])
    kept_Z <- c(kept_Z,Z[keep]); kept_e <- c(kept_e,e[keep])
    kept_pi <- c(kept_pi,pi_select[keep]); generated <- generated+size
  }
  truncation <- 2
  lower <- stats::pnorm(-truncation); upper <- stats::pnorm(truncation)
  scale <- sqrt(1-2*truncation*stats::dnorm(truncation)/(upper-lower))
  noise <- matrix(stats::qnorm(stats::runif(2*n,lower,upper))/scale,n,2)
  L <- as.vector(kept_X %*% k$b)
  T <- kept_X[,1]*kept_X[,2]
  mu0 <- .3*L+.75*T; mu1 <- 1+.5*L+1.05*T
  Y0 <- mu0+noise[,1]; Y1 <- mu1+noise[,1]+noise[,2]
  Xmain <- cbind(intercept=1,kept_X)
  Xcorrect <- cbind(Xmain,interaction=T)
  list(Y=ifelse(kept_Z==1,Y1,Y0),Z=kept_Z,weights=1/kept_pi,
       X=kept_X,ps_main=Xmain,ps_correct=Xcorrect,
       pg_main=Xmain,pg_correct=Xcorrect,mu0=mu0,mu1=mu1,
       population_propensity=kept_e,inclusion_probability=kept_pi,
       truth=wm_bounded_wdsm_truth_reference(k$index,radius),
       radius=radius,index=k$index,candidates_generated=generated,
       selection_envelope=pi_upper,covariate_abs_bound=bound,
       theory_assumptions_automatically_verified=FALSE,
       design="Bounded iid selection tilt of original structural equations; no clusters or fixed finite population")
}
