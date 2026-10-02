# These fixtures check graph algebra and API contracts, not repeated-sampling
# calibration. The direct reference does not use stored outcome shares/rows.
reciprocal_fixture <- function(M = 2L, method = "self_normalized",
                               estimand = "PATE", unit = FALSE) {
  n <- 12L
  x <- c(0.02, 0.13, 0.29, 0.48, 0.69, 0.91,
         0.07, 0.21, 0.36, 0.56, 0.77, 0.96)
  a <- list(Y = c(1, -2, 4, 0, 3, 2, 4, 1, -1, 5, 2, 6),
    Z = rep(0:1, each = 6L), weights = if (unit) rep(1, n) else 1 + (1:n) / 10,
    scores0 = cbind(x, sin(x)), scores1 = cbind(x^2, cos(x), x),
    M = M, estimand = estimand, method = method,
    mean0 = x, mean1 = 1 + x^2, variance = FALSE)
  if (method == "stabilized") {
    a$rho0 <- if (unit) rep(1, n) else 1.5 + x / 3
    a$rho1 <- if (unit) rep(1, n) else 1.7 + x / 5
  }
  if (estimand == "PATT") a[c("scores1", "mean1", "rho1")] <- NULL
  a
}

reciprocal_direct_reference <- function(a, neighbors) {
  n <- length(a$Y); pate <- a$estimand == "PATE"
  w <- a$weights
  mu1 <- if (pate) a$mean1 else rep(0, n)
  own_mu <- ifelse(a$Z == 1, mu1, a$mean0)
  eps <- a$Y - own_mu
  m <- cbind(a$mean0, mu1)
  loads <- numeric(n); query <- if (pate) seq_len(n) else which(a$Z == 1)
  contrasts <- numeric(n)
  for (i in query) {
    j <- neighbors[[i]]; arm <- 1 - a$Z[i]
    if (a$method == "self_normalized") {
      avg <- sum(w[j] * eps[j]) / sum(w[j])
      loads[j] <- loads[j] + w[i] * w[j] / sum(w[j])
    } else {
      rho <- if (arm == 0) a$rho0 else a$rho1
      avg <- sum(w[j] * eps[j] / rho[j]) / a$M
      loads[j] <- loads[j] + w[i] / a$M
    }
    contrasts[i] <- (2 * a$Z[i] - 1) * (a$Y[i] - m[i, arm + 1] - avg)
  }
  outer <- if (pate) w else w * a$Z
  estimate <- sum(outer * contrasts) / sum(outer)
  gamma <- mean(outer)
  if (a$method == "self_normalized") {
    eta <- if (pate) w * (mu1 - a$mean0 - estimate) +
      (2 * a$Z - 1) * (w + loads) * eps else
      a$Z * w * (a$Y - a$mean0 - estimate) - (1 - a$Z) * loads * eps
  } else {
    r <- numeric(n)
    r[a$Z == 0] <- w[a$Z == 0] * eps[a$Z == 0] / a$rho0[a$Z == 0]
    if (pate) {
      r[a$Z == 1] <- w[a$Z == 1] * eps[a$Z == 1] / a$rho1[a$Z == 1]
      rho <- ifelse(a$Z == 1, a$rho1, a$rho0)
      eta <- w * (mu1 - a$mean0 - estimate) +
        (2 * a$Z - 1) * (rho + loads) * r
    } else eta <- a$Z * w * (a$Y - a$mean0 - estimate) -
      (1 - a$Z) * loads * r
  }
  total <- 0; pairs <- 0L
  if (pate) for (i in which(a$Z == 1)) for (j in which(a$Z == 0)) {
    if (j %in% neighbors[[i]] && i %in% neighbors[[j]]) {
      pairs <- pairs + 1L
      total <- total + if (a$method == "self_normalized") {
        w[i]^2 * w[j]^2 * eps[i] * eps[j] /
          (sum(w[neighbors[[j]]]) * sum(w[neighbors[[i]]]))
      } else w[i] * r[i] * w[j] * r[j] / a$M^2
    }
  }
  list(estimate = estimate, root_variance = (sum(eta^2) - 2 * total) / n / gamma^2,
       diagonal = mean(eta^2) / gamma^2, reciprocal = total / n / gamma^2,
       pairs = pairs)
}

test_that("distinct-map reciprocal variance agrees with independent direct sums", {
  for (M in c(1L, 3L)) for (method in c("self_normalized", "stabilized"))
    for (estimand in c("PATE", "PATT")) for (unit in c(FALSE, TRUE)) {
      a <- reciprocal_fixture(M, method, estimand, unit)
      fit <- do.call(wm_match, a)
      before <- fit
      ans <- wm_reciprocal_inference(fit)
      ref <- reciprocal_direct_reference(a, fit$graph$neighbors)
      expect_identical(fit, before)
      expect_s3_class(ans, "wm_reciprocal_inference")
      expect_false(inherits(ans, "wm_match"))
      expect_equal(ans$estimate, ref$estimate, tolerance = 1e-12)
      expect_equal(ans$numerator_variance / ans$gamma^2, ref$root_variance, tolerance = 1e-12)
      expect_equal(ans$diagonal_numerator_variance / ans$gamma^2, ref$diagonal, tolerance = 1e-12)
      expect_equal(ans$reciprocal_numerator / ans$gamma^2, ref$reciprocal, tolerance = 1e-12)
      expect_equal(nrow(ans$reciprocal_pairs), ref$pairs)
      expect_equal(ans$variance, ans$root_n_variance / length(a$Y))
      expect_equal(ans$se^2, ans$variance)
      if (estimand == "PATT") expect_identical(ans$reciprocal_numerator, 0)
    }
})

test_that("unit weights, weight scale, and treatment relabeling have their exact reductions", {
  a <- reciprocal_fixture(M = 3L, unit = TRUE)
  original <- wm_reciprocal_inference(do.call(wm_match, a))
  b <- a; b$method <- "stabilized"; b$rho0 <- b$rho1 <- rep(1, length(a$Y))
  stabilized <- wm_reciprocal_inference(do.call(wm_match, b))
  expect_equal(original$estimate, stabilized$estimate)
  expect_equal(original$root_n_variance, stabilized$root_n_variance)
  for (method in c("self_normalized", "stabilized")) {
    a <- reciprocal_fixture(M = 3L, method = method)
    fit <- wm_reciprocal_inference(do.call(wm_match, a))
    b <- a; b$weights <- 29 * b$weights
    if (method == "stabilized") { b$rho0 <- 29 * b$rho0; b$rho1 <- 29 * b$rho1 }
    scaled <- wm_reciprocal_inference(do.call(wm_match, b))
    expect_equal(scaled$estimate, fit$estimate)
    expect_equal(scaled$root_n_variance, fit$root_n_variance)
    b <- a; b$Z <- 1 - a$Z; b$scores0 <- a$scores1; b$scores1 <- a$scores0
    b$mean0 <- a$mean1; b$mean1 <- a$mean0
    if (method == "stabilized") { b$rho0 <- a$rho1; b$rho1 <- a$rho0 }
    flipped <- wm_reciprocal_inference(do.call(wm_match, b))
    expect_equal(flipped$estimate, -fit$estimate)
    expect_equal(flipped$root_n_variance, fit$root_n_variance)
  }
})

test_that("negative reciprocal variances are retained and explicit floors are recorded", {
  fit <- wm_match(c(1, 1), c(0, 1), c(1, 1), matrix(c(0, 1)), M = 1,
    mean0 = c(0, 2), mean1 = c(2, 0), variance = FALSE)
  ans <- wm_reciprocal_inference(fit)
  expect_equal(ans$diagonal_numerator_variance, 0)
  expect_equal(ans$reciprocal_numerator, 0.5)
  expect_equal(ans$numerator_variance, -1)
  expect_false(ans$available)
  expect_true(is.na(ans$se))
  expect_error(wm_reciprocal_bootstrap(ans), "unavailable")
  floored <- wm_reciprocal_inference(fit, variance_floor = 0.1)
  expect_equal(floored$numerator_variance, -1)
  expect_equal(floored$used_numerator_variance, 0.1)
  expect_equal(floored$variance, 0.05)
  expect_true(floored$floor_active)
  expect_true(floored$available)
})

test_that("stale graphs, adjusted fits, and restricted graphs cannot enter the new contract", {
  fit <- do.call(wm_match, reciprocal_fixture())
  bad <- fit; bad$graph$edges$share[1] <- 0
  expect_error(wm_reciprocal_inference(bad))
  bad <- fit; bad$contributions$actual[1] <- bad$contributions$actual[1] + 1
  expect_error(wm_reciprocal_inference(bad))
  bad <- fit; bad$estimate <- bad$estimate + 1
  expect_error(wm_reciprocal_inference(bad))
  bad <- fit; bad$graph$edges$donor[1] <- bad$graph$edges$query[1]
  expect_error(wm_reciprocal_inference(bad))
  bad <- fit; bad$graph$neighbors[[1]][1] <- 1L
  expect_error(wm_reciprocal_inference(bad))
  bad <- fit; bad$info$nuisance_correction <- TRUE
  expect_error(wm_reciprocal_inference(bad))
  bad <- fit; bad$contributions$training <- 0
  expect_error(wm_reciprocal_inference(bad))
  bad <- fit; bad$graph$cell[1] <- 2L
  expect_error(wm_reciprocal_inference(bad))
  a <- reciprocal_fixture(); a$mean0 <- a$mean1 <- NULL
  expect_error(wm_reciprocal_inference(do.call(wm_match, a)))
})

test_that("Gaussian reciprocal replication has the specified draws and preserves RNG", {
  fit <- do.call(wm_match, reciprocal_fixture())
  inference <- wm_reciprocal_inference(fit)
  set.seed(119); previous <- .Random.seed
  set.seed(45); normal <- rnorm(11); .Random.seed <- previous
  boot <- wm_reciprocal_bootstrap(inference, B = 11, seed = 45)
  expect_identical(.Random.seed, previous)
  expect_equal(boot$root_n_draws, sqrt(inference$root_n_variance) * normal)
  expect_equal(boot$draws, inference$estimate + boot$root_n_draws / sqrt(fit$n))
  expect_equal(boot$conditional_variance, inference$variance)
  expect_equal(boot$monte_carlo_variance, var(boot$draws))
  expect_equal(boot$conf.int, inference$conf.int)
  basic <- wm_reciprocal_bootstrap(inference, B = 11, seed = 45, interval = "basic")
  expect_equal(unname(basic$conf.int), inference$estimate -
    unname(quantile(boot$root_n_draws, c(0.975, 0.025))) / sqrt(fit$n))
  single <- wm_reciprocal_bootstrap(inference, B = 1, seed = 45)
  expect_true(is.na(single$monte_carlo_variance))
  expect_equal(single$conf.int, boot$conf.int)
  expect_null(wm_reciprocal_bootstrap(inference, B = 1, seed = 45, interval = "none")$conf.int)
  expect_error(wm_bootstrap(inference), "wm_match")
  bad <- inference; bad$variance <- bad$variance * 2
  expect_error(wm_reciprocal_bootstrap(bad, seed = 45), "disagree")
  expect_identical(.Random.seed, previous)
  for (B in list(0, -1, 1.1, NA_real_, Inf, 1i))
    expect_error(wm_reciprocal_bootstrap(inference, B = B))
})

test_that("stabilized PATT ignores treated rho values and extreme accepted levels stay finite", {
  a <- reciprocal_fixture(method = "stabilized", estimand = "PATT")
  a$rho0[a$Z == 1] <- 1e-310
  fit <- do.call(wm_match, a)
  ans <- wm_reciprocal_inference(fit, conf.level = 1 - .Machine$double.eps / 2)
  expect_true(all(is.finite(ans$conf.int)))
  expect_true(all(is.finite(wm_reciprocal_bootstrap(ans, B = 1, seed = 1,
    conf.level = 1 - .Machine$double.eps / 2)$conf.int)))
})

test_that("distinct-map stabilized mean fits permit point calculation without old inference", {
  a <- reciprocal_fixture(method = "stabilized")
  a[c("mean0", "mean1", "rho0", "rho1")] <- NULL
  a$rho_bounds <- c(0.5, 3)
  fit <- do.call(wm_fit, a)
  expect_s3_class(fit, "wm_fit")
  expect_true(is.finite(wm_reciprocal_inference(fit)$numerator_variance))
  a$variance <- TRUE
  expect_error(do.call(wm_fit, a), "same supplied score")
})

test_that("reciprocal replication restores absent RNG state and rejects hidden normal caches", {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  saved_seed <- if (had_seed) .Random.seed else NULL
  saved_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(saved_kind))
    if (had_seed) assign(".Random.seed", saved_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  inference <- wm_reciprocal_inference(do.call(wm_match, reciprocal_fixture()))
  RNGkind(normal.kind = "Inversion")
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    rm(".Random.seed", envir = .GlobalEnv)
  wm_reciprocal_bootstrap(inference, B = 1, seed = 19)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  RNGkind(normal.kind = "Box-Muller")
  before <- .Random.seed
  expect_error(wm_reciprocal_bootstrap(inference, B = 1, seed = 19), "Box-Muller")
  expect_identical(.Random.seed, before)
})
