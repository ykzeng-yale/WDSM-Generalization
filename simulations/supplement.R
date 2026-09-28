# Sourceable supplemental fixtures; require package functions, benchmarks.R and
# first_stage.R and first_stage_benchmarks.R for stream/capture/geometry helpers.
wm_supplement_config <- function(replications = 200L) {
  data.frame(scenario = c("spline_d1_n800", "spline_d1_n2000", "spline_d2_n800", "spline_d2_n2000",
    "strata_n800", "strata_n2000", "split_half_n2000", "split_third_n2000"),
    branch = c(rep("spline", 4), rep("strata", 2), rep("split", 2)),
    n = c(800L, 2000L, 800L, 2000L, 800L, 2000L, 2000L, 2000L),
    d = c(1L, 1L, 2L, 2L, 1L, 1L, 2L, 2L), M = 3L,
    fraction1 = c(rep(NA_real_, 6), .5, 1/3), seed = 933101:933108,
    replications = .wm_fs_integer(replications, "replications", 1L, 1000L), stringsAsFactors = FALSE)
}

wm_supplement_generate <- function(cfg) {
  n <- cfg$n; d <- cfg$d
  S <- matrix(stats::runif(n * d), n, d)
  U <- stats::rbinom(n, 1, .5)
  L <- if (cfg$branch == "strata") 1L + stats::rbinom(n, 1, .35) else rep(1L, n)
  q <- if (cfg$branch == "strata") c(.3, .6)[L] else rep(.4, n)
  Z <- stats::rbinom(n, 1, q)
  W <- ifelse(Z == 1L, 1 + 2 * U, .5 + U)
  if (cfg$branch == "spline") {
    mean0 <- 1 + sin(2 * pi * S[, 1]) + (.3 / d) * rowSums(cos(2 * pi * S))
    mean1 <- mean0 + 1 + S[, 1]
  } else if (cfg$branch == "strata") {
    mean0 <- 1 + S[, 1] + S[, 1]^2 / 2
    mean1 <- mean0 + c(.5, 1.5)[L] + S[, 1]
  } else {
    mean0 <- 1 + S[, 1]
    mean1 <- 2 + S[, 1] + S[, 2]
  }
  variance <- ifelse(Z == 1L, ifelse(U == 1L, 1.5, .5), ifelse(U == 1L, 1, .25))
  Y <- ifelse(Z == 1L, mean1, mean0) + sqrt(variance) * (2 * stats::rbinom(n, 1, .5) - 1)
  list(Y = Y, Z = Z, weights = W, scores = S, mean0 = mean0, mean1 = mean1,
    rho0 = rep(1, n), rho1 = rep(2, n), strata = if (cfg$branch == "strata") L else NULL,
    folds = if (cfg$branch == "spline") rep(1:2, length.out = n) else NULL,
    block = if (cfg$branch == "split") as.integer(seq_len(n) <= floor(n * cfg$fraction1)) else NULL)
}

wm_supplement_benchmarks <- function(cfg, geometry = NULL) {
  if (cfg$branch == "strata") {
    beta <- wm_beta_1d(cfg$M)
    prob <- c(.65, .35); q <- c(.3, .6); local_tau <- c(1, 2)
    local <- lapply(q, function(value) wm_benchmark_strong(cfg$M, beta, q = value))
    out <- local[[1]]$table
    for (i in seq_len(nrow(out))) {
      pate <- out$estimand[i] == "PATE"
      gamma <- if (pate) 1 + q else 2 * q
      nu <- if (pate) 1.25 + 3.75 * q else 5 * q
      target <- sum(prob * gamma * local_tau) / sum(prob * gamma)
      out$target[i] <- out$probability_limit[i] <- target
      out$root_n_variance[i] <- sum(prob * (gamma^2 * vapply(local, function(x) x$table$root_n_variance[i], numeric(1)) +
        nu * (local_tau - target)^2)) / sum(prob * gamma)^2
      out$geometry_mcse[i] <- 0
    }
    return(out[c("method", "estimand", "target", "root_n_variance", "geometry_mcse")])
  }
  if (cfg$branch == "spline") {
    alpha <- .wm_fs_geometry(cfg$M, cfg$d, geometry)
    beta <- if (cfg$d == 1L) wm_beta_1d(cfg$M) else if (inherits(geometry, "wm_geometry")) geometry$beta else NULL
    covariance <- if (cfg$d == 1L) matrix(0, cfg$M, cfg$M) else if (inherits(geometry, "wm_geometry")) geometry$covariance else NULL
    if (is.null(beta)) {
      symbolic <- wm_benchmark_strong(cfg$M, wm_beta_1d(cfg$M))
      out <- symbolic$table
      out$root_n_variance <- out$geometry_mcse <- NA_real_
      for (i in which(out$method == "stabilized")) {
        coefficient <- symbolic$beta_coefficients[i, ]
        stopifnot(max(coefficient) == min(coefficient))
        out$root_n_variance[i] <- out$variance_intercept[i] + coefficient[1] * alpha$alpha
        out$geometry_mcse[i] <- coefficient[1] * alpha$alpha_mcse
      }
    } else out <- wm_benchmark_strong(cfg$M, beta, covariance)$table
    if (!alpha$status %in% c("exact_1d", "monte_carlo_estimate")) out$geometry_mcse <- NA_real_
    return(out[c("method", "estimand", "target", "root_n_variance", "geometry_mcse")])
  }
  M <- cfg$M; n1 <- floor(cfg$n * cfg$fraction1); n0 <- cfg$n - n1
  alpha1 <- M^2 + M / 2
  g <- .wm_fs_geometry(M, 2L, geometry)
  alpha2 <- g$alpha; error2 <- g$alpha_mcse
  V0 <- (25/49) * (11/48 + 407/160 + 37/(16*M) + 37*alpha1/(30*M^2))
  V1 <- (25/49) * (11/24 + 7 + 21/(16*M) + 63*alpha2/(40*M^2))
  error1 <- (25/49) * 63/(40*M^2) * error2
  data.frame(method = "stabilized", estimand = c("PATE", "potential0", "potential1"),
    target = c(1.5, 1.5, 3), root_n_variance = c(cfg$n/n1*V1 + cfg$n/n0*V0, V0, V1),
    geometry_mcse = c(cfg$n/n1*error1, 0, error1))
}

wm_supplement_run_rep <- function(cfg, replication) {
  declared <- wm_supplement_config(cfg$replications)
  selected <- declared[declared$scenario == cfg$scenario, , drop = FALSE]
  if (nrow(cfg) != 1L || nrow(selected) != 1L ||
      !isTRUE(all.equal(cfg, selected, check.attributes = FALSE))) stop("Use a declared supplemental configuration")
  replication <- .wm_fs_integer(replication, "replication", 1L, cfg$replications)
  stream <- wm_first_stage_stream(cfg$seed, replication)
  old_kind <- RNGkind(); had <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  old_seed <- if (had) get(".Random.seed", .GlobalEnv) else NULL
  if (old_kind[2] == "Box-Muller") stop("Use a noncached normal RNG")
  on.exit({do.call(RNGkind, as.list(old_kind)); if (had) assign(".Random.seed", old_seed, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)})
  assign(".Random.seed", stream, .GlobalEnv)
  dat <- wm_supplement_generate(cfg)
  records <- list(); details <- list()
  add <- function(captured, estimand, method, correction, analysis_n = cfg$n) {
    fit <- captured$value; ok <- captured$ok
    target <- if (estimand == "potential1") 3 else if (cfg$branch == "strata")
      if (estimand == "PATE") 393/281 else 41/27 else 1.5
    records[[length(records)+1L]] <<- data.frame(scenario = cfg$scenario, branch = cfg$branch,
      replication = replication, n = cfg$n, analysis_n = analysis_n, d = cfg$d, M = cfg$M,
      fraction1 = cfg$fraction1, base_seed = cfg$seed, rng_state = paste(stream, collapse = ":"),
      method = method, estimand = estimand, correction = correction, target = target,
      status = if (ok) "ok" else "error", estimate = if (ok) fit$estimate else NA_real_,
      root_n_variance = if (ok) fit$root_n_variance else NA_real_,
      sampling_variance = if (ok) fit$variance else NA_real_,
      lower = if (ok) fit$estimate - stats::qnorm(.975)*fit$se else NA_real_,
      upper = if (ok) fit$estimate + stats::qnorm(.975)*fit$se else NA_real_,
      clipped = if (ok && !is.null(fit$nuisance)) fit$nuisance$rho_clipped_count else 0L,
      elapsed_seconds = captured$elapsed, warning = captured$warning, error = captured$error,
      stringsAsFactors = FALSE)
    details[[paste(estimand, method, correction, sep = "/")]] <<- if (ok) list(
      nuisance = fit$nuisance, numerator_identity = fit$numerator_identity,
      analysis_n = analysis_n, graph_restrictions = list(folds = dat$folds, strata = dat$strata),
      block_sizes = if (cfg$branch == "split") table(dat$block) else NULL) else list(error = captured$error)
  }
  if (cfg$branch == "split") {
    cap <- .wm_fs_capture(function() wm_split_pate(dat$Y, dat$Z, dat$weights,
      dat$scores[, 1, drop = FALSE], dat$scores, block_id = dat$block,
      mean0 = dat$mean0, mean1 = dat$mean1, rho0 = dat$rho0, rho1 = dat$rho1, M = cfg$M))
    add(cap, "PATE", "stabilized", "oracle")
    for (z in 0:1) {
      component <- cap
      if (cap$ok) component$value <- cap$value$components[[z+1L]]
      component$elapsed <- 0
      add(component, paste0("potential", z), "stabilized", "oracle", sum(dat$block == z))
    }
  } else for (estimand in c("PATE", "PATT")) for (method in c("self_normalized", "stabilized"))
    for (correction in c("oracle", "feasible")) {
      cap <- .wm_fs_capture(function() {
        arguments <- list(Y = dat$Y, Z = dat$Z, weights = dat$weights, scores0 = dat$scores,
          M = cfg$M, estimand = estimand, method = method, strata = dat$strata, fold_id = dat$folds)
        if (correction == "oracle") {
          arguments$mean0 <- dat$mean0
          if (estimand == "PATE") arguments$mean1 <- dat$mean1
          if (method == "stabilized") {
            arguments$rho0 <- dat$rho0
            if (estimand == "PATE") arguments$rho1 <- dat$rho1
          }
          do.call(wm_match, arguments)
        } else {
          if (cfg$branch == "spline") {
            arguments$regression <- "spline"; arguments$spline_order <- 3L
            arguments$mesh_exponent <- 1/6; arguments$moment_order <- Inf
            arguments$support0 <- cbind(rep(0, cfg$d), rep(1, cfg$d))
            if (estimand == "PATE") arguments$support1 <- arguments$support0
          } else {arguments$regression <- "polynomial"; arguments$degree <- 2L}
          if (method == "stabilized") arguments$rho_bounds <- c(.25, 4)
          do.call(wm_fit, arguments)
        }
      })
      add(cap, estimand, method, correction)
    }
  list(records = do.call(rbind, records), diagnostics = details)
}
