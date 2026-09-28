# Bounded base-R checks before the complete package test suite.
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1L]] else "."
for (f in sort(list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE))) source(f)
assert_close <- function(x, y, tolerance = 1e-10) {
  stopifnot(isTRUE(all.equal(x, y, tolerance = tolerance, check.attributes = FALSE)))
}
x <- matrix(c(0, 3, 1, 2), ncol = 1)
y <- c(1, 5, 7, 11); z <- c(0, 0, 1, 1); w <- 1:4
f <- wm_match(y, z, w, x, M = 1, mean0 = rep(0, 4), mean1 = rep(0, 4))
stopifnot(identical(f$graph$neighbors, list(3L, 4L, 1L, 2L)))
assert_close(f$estimate, 6)
assert_close(f$contributions$row, c(-4, -16.8, 4, 16.8))
assert_close(f$root_n_variance, 149.12)
g <- wm_match(y, z, w, x, M = 1, mean0 = rep(0, 4), estimand = "PATT")
assert_close(g$estimate, 6)
assert_close(g$contributions$row, c(-12, -80, 12, 80) / 7)
assert_close(g$root_n_variance, 3272 / 49)
set.seed(928)
rng_before <- .Random.seed
b <- wm_bootstrap(f, B = 41, seed = 17, chunk_size = 3)
stopifnot(identical(rng_before, .Random.seed))
assert_close(b$conditional_root_n_variance, 149.12)
assert_close(b$conditional_variance, 37.28)
assert_close(b$conf.int, 6 + c(-1, 1) * qnorm(.975) * sqrt(37.28))
checked <- 0L
for (d in c(1L, 2L, 3L, 6L)) for (M in c(1L, 3L)) {
  n <- 80L
  S <- matrix(runif(n * d), n, d)
  Z <- rep(0:1, n / 2)
  W <- 1 + runif(n)
  mu0 <- rowSums(S) / d
  mu1 <- mu0 + 1 + S[, 1]
  Y <- ifelse(Z == 1, mu1, mu0) + runif(n, -.5, .5)
  for (estimand in c("PATE", "PATT")) for (method in c("self_normalized", "stabilized")) {
    a <- list(Y = Y, Z = Z, weights = W, scores0 = S, M = M,
              mean0 = mu0, estimand = estimand, method = method)
    if (estimand == "PATE") a$mean1 <- mu1
    if (method == "stabilized") {
      a$rho0 <- rep(1.5, n)
      if (estimand == "PATE") a$rho1 <- rep(1.5, n)
    }
    fit <- do.call(wm_match, a)
    assert_close(fit$numerator_identity$difference, 0)
    assert_close(fit$numerator_identity$row_edge_sum, 0)
    assert_close(fit$root_n_variance,
                 (sum(fit$contributions$row^2) + sum(fit$contributions$edge^2)) / n)
    boot <- wm_bootstrap(fit, B = 3, seed = 92)
    assert_close(boot$conditional_variance, fit$variance)
    a$weights <- a$weights * 1e100
    if (method == "stabilized") {
      a$rho0 <- a$rho0 * 1e100
      if (estimand == "PATE") a$rho1 <- a$rho1 * 1e100
    }
    scaled <- do.call(wm_match, a)
    assert_close(fit$estimate, scaled$estimate)
    assert_close(fit$root_n_variance, scaled$root_n_variance)
    checked <- checked + 1L
  }
}
cat("PASS: hand-calculated PATE/PATT fixtures, exact Gaussian law/RNG restoration,", checked,
    "dimensional/method/estimand cases with identities and scale invariance.\n")
print(sessionInfo())
