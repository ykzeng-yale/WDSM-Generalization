# Fixed compatibility fixtures; no Monte Carlo study or nuisance calibration.
wm_tie_fixture <- function() {
  list(Y = c(1, 3, 8, 4, 7, 13, 2, 9, 5, 12, 6, 10),
       Z = rep(0:1, 6), weights = c(1, 2, 4, 3, 5, 7, 2, 1, 3, 4, 8, 6),
       scores0 = cbind(rep(c(0, 0, 1, 1, 2, 2), 2), rep(c(0, 1), 6)),
       scores1 = cbind(rep(c(2, 1, 0), 4), rep(c(1, 1, 0, 0), 3)),
       mean0 = seq_len(12) / 4, mean1 = 1 + seq_len(12) / 3,
       M = 3L, variance = FALSE)
}

wm_tie_restore_rng <- function() {
  kind <- RNGkind()
  had <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  seed <- if (had) get(".Random.seed", .GlobalEnv, inherits = FALSE) else NULL
  function() {
    do.call(RNGkind, as.list(kind))
    if (had) assign(".Random.seed", seed, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }
}

test_that("the default exact-distance policy remains deterministic without RNG use", {
  restore <- wm_tie_restore_rng(); on.exit(restore(), add = TRUE)
  args <- wm_tie_fixture()
  set.seed(938L); before <- .Random.seed; kind <- RNGkind()
  default <- do.call(wm_match, args)
  explicit <- do.call(wm_match, c(args, list(tie_rule = "row_order")))
  expect_identical(default, explicit)
  expect_identical(.Random.seed, before)
  expect_identical(RNGkind(), kind)
  expect_identical(default$graph$tie_rule, "Exact distance, then original row index")
  # Distances beyond squared-distance range still use the frozen robust norm.
  args$scores0 <- args$scores0 * 2^660
  args$scores1 <- args$scores1 * 2^660
  expect_identical(do.call(wm_match, args)$graph$neighbors, default$graph$neighbors)
  rm(".Random.seed", envir = .GlobalEnv)
  invisible(do.call(wm_match, args))
  expect_false(exists(".Random.seed", .GlobalEnv, inherits = FALSE))
})

test_that("source randomized donors, fractions and reindexed points agree for both targets", {
  args <- wm_tie_fixture()
  source <- wdsmatch:::wdsm_make_matches
  for (estimand in c("PATE", "PATT")) for (M in c(1L, 3L, 5L)) {
    args$estimand <- estimand; args$M <- M
    args$mean1 <- if (estimand == "PATE") 1 + seq_len(12) / 3 else NULL
    args$scores1 <- if (estimand == "PATE") wm_tie_fixture()$scores1 else NULL
    seed <- 20260917L + M
    oracle <- source(args$Z, args$scores0, args$scores1, M = M,
                     estimand = estimand, tie_seed = seed)
    fit <- do.call(wm_match, c(args, list(tie_rule = "source_random", tie_seed = seed)))
    queries <- if (estimand == "PATE") seq_along(args$Y) else which(args$Z == 1L)
    expected <- vector("list", length(args$Y))
    incoming <- matrix(0, length(args$Y), 2)
    contrast <- raw_contrast <- numeric(length(args$Y))
    for (i in queries) {
      arm <- 1L - args$Z[i]
      j <- oracle[[paste0("matches_", arm)]][[i]]
      expected[[i]] <- j
      fractions <- args$weights[j] / sum(args$weights[j])
      edge <- fit$graph$edges[fit$graph$edges$query == i, ]
      expect_equal(edge$share, fractions, tolerance = 1e-14)
      mu <- if (arm == 0L) args$mean0 else args$mean1
      missing <- mu[i] + sum(fractions * (args$Y[j] - mu[j]))
      contrast[i] <- (2 * args$Z[i] - 1) * (args$Y[i] - missing)
      raw_contrast[i] <- (2 * args$Z[i] - 1) *
        (args$Y[i] - sum(fractions * args$Y[j]))
      incoming[j, arm + 1L] <- incoming[j, arm + 1L] +
        args$weights[i] / max(args$weights) * fractions
    }
    expect_identical(fit$graph$neighbors, expected)
    expect_identical(fit$graph$edges$query, rep(queries, each = M))
    expect_identical(fit$graph$tie_diagnostics$arms, oracle$tie_diagnostics$arms)
    expect_equal(unname(fit$loads$incoming), incoming, tolerance = 1e-14)
    expect_equal(fit$estimate, sum(args$weights[queries] * contrast[queries]) /
      sum(args$weights[queries]), tolerance = 1e-13)
    expect_equal(fit$raw_estimate, sum(args$weights[queries] * raw_contrast[queries]) /
      sum(args$weights[queries]), tolerance = 1e-13)
  }
})

test_that("arbitrary map dimensions and restricted cells retain source stream ordering", {
  args <- wm_tie_fixture()
  # Complete ties have known geometry in every dimension, including unequal maps.
  args$scores0 <- matrix(0, 12, 1)
  args$scores1 <- matrix(0, 12, 4)
  common <- wdsmatch:::wdsm_make_matches(args$Z, matrix(0, 12, 2), matrix(0, 12, 2),
    M = 3L, tie_seed = 20260920L)
  fit <- do.call(wm_match, c(args, list(tie_rule = "source_random", tie_seed = 20260920L)))
  expected <- lapply(seq_len(12), function(i)
    common[[paste0("matches_", 1L - args$Z[i])]][[i]])
  expect_identical(fit$graph$neighbors, expected)
  expect_identical(dim(fit$graph$scores0), c(12L, 1L))
  expect_identical(dim(fit$graph$scores1), c(12L, 4L))
  # An informative third coordinate must participate in Euclidean matching.
  coordinates <- cbind(c(0, 1, 2, 0.25, 1.4, 3), 0, c(0, 0, 0, 10, 0, 0))
  higher <- wm_match(c(1, 3, 8, 4, 7, 13), c(0, 0, 0, 1, 1, 1), rep(1, 6),
    coordinates, M = 1L, variance = FALSE, tie_rule = "source_random", tie_seed = 7L)
  expect_identical(higher$graph$neighbors, list(5L, 5L, 5L, 1L, 2L, 3L))
  args$M <- 1L; args$strata <- rep(c("a", "b"), each = 6)
  args$fold_id <- rep("evaluation", 12)
  restricted <- do.call(wm_match, c(args, list(tie_rule = "source_random", tie_seed = 11L)))
  restore <- wm_tie_restore_rng(); on.exit(restore(), add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(11L)
  expected <- vector("list", 12)
  for (arm in 0:1) for (i in which(args$Z != arm)) {
    eligible <- which(args$Z == arm & args$strata == args$strata[i])
    expected[[i]] <- eligible[sample.int(length(eligible), 1L)]
  }
  expect_identical(restricted$graph$neighbors, expected)
  expect_true(all(args$strata[restricted$graph$edges$query] ==
                    args$strata[restricted$graph$edges$donor]))
})

test_that("near boundaries use squared tolerance and draw only on oversized boundaries", {
  Z <- c(0, 1, 0, 1, 0, 1)
  # Two control donors differ in squared distance by less than the source width.
  S <- cbind(c(1, 0, 1 + 2 * .Machine$double.eps, 0, 4, 0), 0)
  for (tol in c(0, 64 * .Machine$double.eps)) {
    source <- wdsmatch:::wdsm_make_matches(Z, S, M = 1, estimand = "PATT",
      tie_seed = 7L, tie_tolerance = tol)
    fit <- wm_match(seq_len(6), Z, rep(1, 6), S, M = 1,
      estimand = "PATT", variance = FALSE, tie_rule = "source_random",
      tie_seed = 7L, tie_tolerance = tol)
    expect_identical(fit$graph$neighbors, source$matches_0)
    expect_identical(fit$graph$tie_diagnostics$arms, source$tie_diagnostics$arms)
    expect_identical(fit$graph$tie_diagnostics$arms[[1]]$randomized_recipients,
                     if (tol == 0) 0L else 3L)
  }
  # Selecting all eligible donors consumes no sample draw despite exact ties.
  fit <- wm_match(seq_len(6), Z, rep(1, 6), matrix(0, 6, 3), M = 3,
    variance = FALSE, tie_rule = "source_random", tie_seed = 7L)
  expect_identical(vapply(fit$graph$tie_diagnostics$arms, `[[`, integer(1),
                          "randomized_recipients"), c(`0` = 0L, `1` = 0L))
})

test_that("source seeds and caller RNG are preserved on success and failures", {
  restore <- wm_tie_restore_rng(); on.exit(restore(), add = TRUE)
  args <- wm_tie_fixture()
  args$scores0 <- matrix(0, 12, 2); args$scores1 <- matrix(0, 12, 2)
  opts <- list(tie_rule = "source_random", tie_seed = 19L)
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); set.seed(84L)
  before <- .Random.seed; kind <- RNGkind()
  first <- do.call(wm_match, c(args, opts))
  expect_identical(.Random.seed, before); expect_identical(RNGkind(), kind)
  expect_identical(first, do.call(wm_match, c(args, opts)))
  opts$tie_seed <- 20L
  expect_false(identical(first$graph$neighbors, do.call(wm_match, c(args, opts))$graph$neighbors))
  args$M <- 7L
  expect_error(do.call(wm_match, c(args, opts)), "Insufficient donors")
  expect_identical(.Random.seed, before); expect_identical(RNGkind(), kind)
  args$M <- 3L; args$scores0[1, 1] <- 2^660
  expect_error(do.call(wm_match, c(args, opts)), "Squared distances exceed")
  expect_identical(.Random.seed, before); expect_identical(RNGkind(), kind)
  args$scores0[1, 1] <- 0
  rm(".Random.seed", envir = .GlobalEnv)
  invisible(do.call(wm_match, c(args, opts)))
  expect_false(exists(".Random.seed", .GlobalEnv, inherits = FALSE))
  expect_identical(RNGkind(), kind)
  RNGkind("Mersenne-Twister", "Box-Muller", "Rejection"); set.seed(44L)
  invisible(rnorm(1)); before <- .Random.seed
  expect_error(do.call(wm_match, c(args, opts)), "Box-Muller cache")
  expect_identical(.Random.seed, before)
  cached <- rnorm(1)
  set.seed(44L); expect_identical(cached, rnorm(2)[2])
})

test_that("source policy has an explicit unsupported sampling-inference contract", {
  args <- wm_tie_fixture()
  opts <- list(tie_rule = "source_random", tie_seed = 20260920L)
  fit <- do.call(wm_match, c(args, opts))
  args$variance <- TRUE
  expect_error(do.call(wm_match, c(args, opts)), "point-only")
  expect_error(do.call(wm_fit, c(args[names(args) != "mean0" & names(args) != "mean1"], opts)), "point-only")
  expect_error(wm_reciprocal_inference(fit), "Unsupported graph metric or tie convention")
  expect_error(wm_bootstrap(fit, B = 2L), "enable variance estimation")
  counts <- matrix(1, 12, 2)
  replicated <- wm_bootstrap_refit(fit, counts)
  expect_identical(wm_bootstrap(fit, counts = counts, method = "fixed_reuse"), replicated)
  expect_equal(replicated$draws, rep(fit$estimate, 2))
  expect_false(replicated$sampling_inference_available)
  expect_match(replicated$inference_contract, "no supported sampling-inference theorem")
  expect_false(fit$graph$tie_diagnostics$sampling_inference_supported)
  for (seed in list(NULL, -1, NA_real_, Inf, 1.5, c(1, 2))) {
    invalid <- opts; invalid$tie_seed <- seed
    expect_error(do.call(wm_match, c(wm_tie_fixture(), invalid)), "tie_seed")
  }
  for (tol in list(-1, NA_real_, Inf, c(0, 1))) {
    invalid <- opts; invalid$tie_tolerance <- tol
    expect_error(do.call(wm_match, c(wm_tie_fixture(), invalid)), "tie_tolerance")
  }
})

test_that("general nuisance and WDSM wrappers forward source policy without changing fits", {
  data(survey_obs)
  D <- cbind(intercept = 1, as.matrix(survey_obs[, c("X1", "X2", "X3")]))
  opts <- list(tie_rule = "source_random", tie_seed = 20260920L)
  for (estimand in c("PATE", "PATT")) {
    args <- list(Y = survey_obs$Y, Z = survey_obs$Z, weights = survey_obs$survey_weight,
      ps_design = D, pg0_design = D, pg1_design = if (estimand == "PATE") D else NULL,
      estimand = estimand, M = 3L, inference = "none")
    a <- do.call(wm_wdsm_fit, args)
    b <- do.call(wm_wdsm_fit, c(args, opts))
    expect_identical(a$nuisance$parameter, b$nuisance$parameter)
    expect_identical(a$nuisance$mean0, b$nuisance$mean0)
    expect_identical(a$nuisance$mean1, b$nuisance$mean1)
    source <- wdsmatch:::wdsm_make_matches(args$Z, b$nuisance$scores0,
      b$nuisance$scores1, M = 3L, estimand = estimand, tie_seed = opts$tie_seed)
    rows <- if (estimand == "PATE") seq_len(b$n) else which(args$Z == 1L)
    for (i in rows) {
      expect_identical(b$fit$graph$neighbors[[i]],
        source[[paste0("matches_", 1L - args$Z[i])]][[i]])
    }
    expect_false(b$inference$sampling_inference_available)
    expect_match(b$inference$unavailable_reason, "unsupported for source_random")
    args$inference <- "full_x"
    expect_error(do.call(wm_wdsm_fit, c(args, opts)), "point-only")
    general <- wm_fit(args$Y, args$Z, args$weights, b$nuisance$scores0,
      scores1 = if (estimand == "PATE") b$nuisance$scores1 else NULL,
      M = 3L, estimand = estimand, degree = 1, variance = FALSE,
      tie_rule = "source_random", tie_seed = opts$tie_seed)
    expect_identical(general$graph$neighbors, b$fit$graph$neighbors)
    expect_match(general$nuisance$inference_contract$tie_policy, "no sampling inference theorem")
  }
})
