args <- commandArgs(trailingOnly = TRUE)
base <- if (length(args)) args[1L] else "."
source(file.path(base, "simulations", "wm_graph_transport_reference.R"))
checks <- 0L
assert <- function(value) {
  if (!isTRUE(value)) stop("Check ", checks + 1L, " failed.")
  checks <<- checks + 1L
}
exact_tagged <- function(w, p, g, root, M) {
  d <- ncol(g)
  if (M == 1) return(c(1, rep(0, d)))
  tuple <- as.matrix(expand.grid(rep(list(seq_along(w)), M - 1L)))
  out <- numeric(d + 1L)
  for (row in seq_len(nrow(tuple))) {
    ids <- tuple[row, ]
    fraction <- root / (root + sum(w[ids]))
    out[1L] <- out[1L] + fraction * prod(p[ids])
    for (r in seq_along(ids)) {
      out[-1L] <- out[-1L] + fraction * g[ids[r], ] * prod(p[ids[-r]])
    }
  }
  out
}
w <- c(0.8, 1.7, 2.4)
p <- c(.2, .35, .45)
g <- cbind(c(.4, -.1, -.3), c(-.2, .7, -.5))
max_error <- 0
for (M in c(1, 2, 3, 5)) for (root in c(.6, 1.2, 2.4)) {
  fit <- wm_tagged_load_reference(w, p, g, root, M, tolerance = 1e-8)
  oracle <- exact_tagged(w, p, g, root, M)
  error <- abs(c(fit$value, fit$gradient) - oracle)
  max_error <- max(max_error, error)
  assert(all(error <= fit$error_bound + 2e-14))
  assert(all(fit$error_bound <= 1e-8))
  scaled <- wm_tagged_load_reference(12*w, p, g, 12*root, M, tolerance = 1e-8)
  assert(max(abs(c(fit$value,fit$gradient) - c(scaled$value,scaled$gradient))) < 2e-14)
}
mass <- wm_tagged_load_reference(rep(2,3), p, g, 3, 5)
assert(abs(mass$value - 3/11) < 1e-15 && all(mass$gradient == 0))
assert(inherits(try(wm_tagged_load_reference(w,p,g,1,3,tolerance=1e-12,max_nodes=3),
                    silent=TRUE), "try-error"))

# Independently differentiate exact local-law enumeration in the raw coordinate.
grid <- expand.grid(x = seq(-.4,.4,length.out=4), y=seq(-.3,.3,length.out=4))
R <- as.matrix(grid)
n <- nrow(R)
Z <- rep(c(0,1), n/2)
W <- c(.8,1.3,1.7)[1L + (seq_len(n) %% 3L)]
residual <- sin(seq_len(n))
tangent <- array(cos(seq_len(n*2*2)), c(n,2,2),
                 dimnames=list(NULL,NULL,c("alpha","beta")))
cutoff <- rep(.75,n)
bw <- .95
kernel_value <- function(point) {
  u <- sweep(R,2L,point,`-`) / bw
  apply((35/32)*pmax(0,1-u^2)^3*(abs(u)<1),1L,prod)
}
load_exact <- function(point, donor_arm, root, M) {
  K <- kernel_value(point)
  donor <- which(Z == donor_arm)
  query <- which(Z != donor_arm)
  # M=3 allows direct finite empirical enumeration, independent of integration.
  A <- exact_tagged(W[donor], K[donor]/sum(K[donor]),
                   matrix(0,length(donor),1),root,M)[1L]
  M * sum(W[query]*K[query]) / sum(K[donor]) * A
}
max_gradient_error <- 0
for (arm in 0:1) {
  fit <- wm_graph_transport_reference(R,Z,W,residual,tangent,cutoff,arm,3,
    bandwidth=bw,density_floor=1e-6,quadrature_tolerance=1e-9)
  grad <- matrix(0,length(fit$evaluated_donor_rows),2)
  for (k in seq_along(fit$evaluated_donor_rows)) {
    i <- fit$evaluated_donor_rows[k]
    assert(abs(fit$incoming_load[k] - load_exact(R[i,],arm,W[i],3)) < 1e-8)
    for (d in 1:2) {
      shift <- c(0,0); shift[d] <- 1e-5
      grad[k,d] <- (load_exact(R[i,]+shift,arm,W[i],3) -
                     load_exact(R[i,]-shift,arm,W[i],3)) / 2e-5
    }
  }
  max_gradient_error <- max(max_gradient_error, abs(grad-fit$incoming_load_gradient))
  assert(max(abs(grad-fit$incoming_load_gradient)) < 1e-7)
  expected <- numeric(2)
  for (k in seq_along(fit$evaluated_donor_rows)) {
    i <- fit$evaluated_donor_rows[k]
    expected <- expected + residual[i]*cutoff[i]*
      as.vector(crossprod(matrix(tangent[i,,],2,2),grad[k,])) / n
  }
  assert(max(abs(expected-fit$graph_drift)) < 1e-7)
  scaled <- wm_graph_transport_reference(R,Z,10*W,residual,tangent,cutoff,arm,3,
    bandwidth=bw,density_floor=1e-6,quadrature_tolerance=1e-9)
  assert(max(abs(scaled$graph_drift-10*fit$graph_drift)) < 1e-10)
}
assert(inherits(try(wm_graph_transport_reference(R,Z,W,residual,tangent,cutoff,0,3,
   bandwidth=bw,density_floor=100),silent=TRUE),"try-error"))
zero <- wm_graph_transport_reference(R,Z,W,residual,tangent,rep(0,n),0,3)
assert(all(zero$graph_drift == 0) && length(zero$evaluated_donor_rows) == 0)
# Unit weights remove M from the first incoming-load moment in every dimension.
unit_reference <- NULL
for (M in c(1,3,5)) {
  unit <- wm_graph_transport_reference(R,Z,rep(1,n),residual,tangent,cutoff,0,M,
                                      bandwidth=bw,density_floor=1e-6)
  if (is.null(unit_reference)) unit_reference <- unit
  assert(max(abs(unit$incoming_load-unit_reference$incoming_load)) < 1e-14)
  assert(max(abs(unit$incoming_load_gradient-unit_reference$incoming_load_gradient)) < 1e-14)
  assert(all(unit$quadrature_nodes == 0) && all(unit$quadrature_error_bound == 0))
}
cat(sprintf("PASS: %d deterministic checks; tagged integral error %.3g; raw-gradient error %.3g. No Monte Carlo.\n",
            checks,max_error,max_gradient_error))
