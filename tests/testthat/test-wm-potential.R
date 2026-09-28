test_that("single potential means agree with hand-directed imputation", {
  Y <- c(1, 5, 7, 11); Z <- c(0, 0, 1, 1); W <- 1:4
  S <- matrix(c(0, 3, 1, 2))
  f0 <- wm_potential(Y, Z, W, S, arm = 0, M = 1, mean = rep(0, 4), rho = rep(1.5, 4))
  expect_equal(f0$imputed, c(1, 5, 2/3, 20/3))
  expect_equal(f0$estimate, 119/30)
  expect_equal(f0$contributions$row, c(-29, 862, -357, -476)/75)
  expect_equal(f0$root_n_variance, sum(c(-29, 862, -357, -476)^2)/(4*75^2))
  expect_equal(f0$numerator_identity$difference, 0, tolerance = 1e-12)
  expect_equal(sum(f0$contributions$row), 0, tolerance = 1e-12)
  expect_length(f0$contributions$edge, 0)
  expect_true(all(f0$graph$edges$arm == 0))
  f1 <- wm_potential(Y, Z, W, S, arm = 1, M = 1, mean = rep(0, 4), rho = rep(3.5, 4))
  expect_equal(f1$imputed, c(6, 88/7, 7, 11))
  expect_equal(f1$estimate, 673/70)
  expect_true(all(f1$graph$edges$arm == 1))
  expect_equal(f1$root_n_variance, mean(f1$contributions$row^2))
  expect_equal(wm_bootstrap(f1, B = 2, seed = 91)$conditional_variance, f1$variance)
})

test_that("query outcomes are irrelevant to a one-direction potential mean", {
  Y <- c(1, 2, 1e200, -1e200); Z <- c(0, 0, 1, 1)
  S <- matrix(c(0, 1, .1, .9))
  fit <- wm_potential(Y, Z, rep(1, 4), S, 0, M = 1,
                     mean = rep(0, 4), rho = rep(1, 4))
  expect_equal(fit$estimate, 1.5)
  expect_equal(fit$imputed, c(1, 2, 1, 2))
  expect_true(is.finite(fit$root_n_variance))
  fitted <- wm_potential_fit(Y, Z, rep(1, 4), S, 0, M = 1, degree = 0,
                            rho_bounds = c(.5, 2))
  expect_equal(fitted$estimate, 1.5)
  expect_true(is.finite(fitted$root_n_variance))
  expect_equal(fit$data$Y, Y)
  expect_error(wm_potential(c(1, 2, Inf, 3), Z, rep(1, 4), S, 0,
                           mean = rep(0, 4), rho = rep(1, 4)), "finite numeric")
})

test_that("potential fitting uses only the requested donor arm", {
  x <- rep(seq(0, 1, length.out = 6), 2)
  Z <- rep(0:1, each = 6)
  Y <- ifelse(Z == 0, 1 + 2*x, 8 - x)
  fit <- wm_potential_fit(Y, Z, 1 + x, matrix(x), arm = 1, M = 2,
                         rho_bounds = c(.5, 3))
  expect_equal(fit$predictions$mean, 8 - x, tolerance = 1e-12)
  expect_equal(fit$estimate, sum((1+x)*(8-x))/sum(1+x), tolerance = 1e-12)
  expect_equal(fit$nuisance$donor_arm, 1L)
  expect_equal(fit$nuisance$models[[1]]$training_rows, which(Z == 1))
  expect_null(fit$nuisance$inference_contract$bounded_outcomes)
})

test_that("independent blocks rescale by their own sample sizes", {
  Y <- c(1,5,7,11, 1,5,7,11,1,7)
  Z <- c(0,0,1,1, 0,0,1,1,0,1)
  W <- c(1:4,1:4,2,3)
  S <- matrix(c(0,3,1,2, 0,3,1,2,4,4.1))
  blocks <- c(rep(0,4),rep(1,6))
  fit <- wm_split_pate(Y,Z,W,S,cbind(S,2*S),blocks,
                       mean0=rep(0,10),mean1=rep(0,10),
                       rho0=rep(1.5,10),rho1=rep(3.5,10),M=1)
  f0 <- wm_potential(Y[1:4],Z[1:4],W[1:4],S[1:4,,drop=FALSE],0,
                     mean=rep(0,4),rho=rep(1.5,4),M=1)
  f1 <- wm_potential(Y[5:10],Z[5:10],W[5:10],cbind(S,2*S)[5:10,,drop=FALSE],1,
                     mean=rep(0,6),rho=rep(3.5,6),M=1)
  expect_equal(fit$estimate,f1$estimate-f0$estimate)
  expect_equal(fit$variance,f1$variance+f0$variance)
  expect_equal(fit$root_n_variance,10*(f1$root_n_variance/6+f0$root_n_variance/4))
  expect_equal(fit$contributions$row,c(-10/4*f0$contributions$row,10/6*f1$contributions$row))
  expect_equal(fit$components$potential0$global_rows,1:4)
  expect_equal(fit$components$potential1$global_rows,5:10)
  expect_equal(wm_bootstrap(fit,B=3,seed=91)$conditional_variance,fit$variance)
})

test_that("potential means preserve weight scale and reject invalid coding", {
  S <- matrix(c(0,3,1,2)); Y <- c(1,5,7,11); Z <- c(0,0,1,1)
  a <- list(Y=Y,Z=Z,weights=1:4,scores=S,arm=0,M=1,mean=rep(0,4),rho=rep(1.5,4))
  f <- do.call(wm_potential,a)
  a$weights <- a$weights*1e100; a$rho <- a$rho*1e100
  fs <- do.call(wm_potential,a)
  expect_equal(fs$estimate,f$estimate)
  expect_equal(fs$root_n_variance,f$root_n_variance)
  a$arm <- 2
  expect_error(do.call(wm_potential,a),"arm must")
  a$arm <- 0; a$Z[1] <- 2
  expect_error(do.call(wm_potential,a),"Z must")
})
