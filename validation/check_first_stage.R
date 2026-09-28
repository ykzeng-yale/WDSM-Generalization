# Deterministic first-stage algebra/pipeline checks, not coverage experiments.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this script with Rscript")
script <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
env <- new.env(parent = globalenv())
for (file in sort(list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)))
  sys.source(file, envir = env)
sys.source(file.path(root, "simulations", "first_stage.R"), envir = env)
sys.source(file.path(root, "simulations", "first_stage_benchmarks.R"), envir = env)
equal <- function(x, y, tolerance = 1e-10) {
  comparison <- all.equal(x, y, tolerance = tolerance, check.attributes = FALSE)
  if (!isTRUE(comparison)) stop(paste(comparison, collapse = "; "))
}
expect_error <- function(fun, pattern) {
  message <- tryCatch({ fun(); NULL }, error = conditionMessage)
  stopifnot(is.character(message), grepl(pattern, message))
}

# Exact score-model integrals and independent/same-sample corrections.
integral <- function(fun) stats::integrate(fun, 0, 1, rel.tol = 1e-12)$value
equal(integral(function(s) s^2 * (1 - s)^2), 1 / 30)
equal(integral(function(s) s * (1 - s)), 1 / 6)
same <- env$wm_first_stage_benchmark("same_sample", 2L, 3L)
independent <- env$wm_first_stage_benchmark("independent_training", 1L, 3L, n = 80L, m = 40L)
equal(same$constants$prediction_sensitivity, -0.2)
equal(same$constants$parameter_variance, 30 / 16)
equal(same$constants$same_sample_covariance, 6 / 16)
equal(unique(same$table$first_stage_sampling_change[same$table$comparison != "known_first_stage"]),
      -0.075)
equal(unique(independent$table$first_stage_sampling_change[
  independent$table$comparison != "known_first_stage"]), 0.15)
stopifnot(all(is.na(same$table$benchmark_root_n_variance)),
          all(independent$table$alpha == 10.5),
          all(independent$table$benchmark_geometry_mcse == 0))
naive <- independent$table[independent$table$comparison == "fitted_naive", ]
equal(naive$benchmark_root_n_variance - naive$reported_root_n_variance_limit, rep(0.15, 2))
expect_error(function() env$wm_first_stage_benchmark("same_sample", 1L, 3L),
             "scalar")

# Synthetic geometry object tests propagation/status only, not an alpha truth.
geometry <- structure(list(M = 3L, d = 2L, alpha = 12, alpha_mcse = 0.4,
  precision_status = "monte_carlo_estimate",
  component_precision = rep("monte_carlo_estimate", 3L), draws = 100L, seed = 31L,
  diagnostics = list(nonzero_draws = c(20L, 70L, 100L))), class = "wm_geometry")
reported <- env$wm_first_stage_benchmark("same_sample", 2L, 3L, geometry = geometry)$table
equal(reported$benchmark_geometry_mcse, 0.4 * reported$alpha_coefficient)
geometry$component_precision[1L] <- "unresolved_zero_hits"
unresolved <- env$wm_first_stage_benchmark("same_sample", 2L, 3L, geometry = geometry)$table
stopifnot(all(is.na(unresolved$benchmark_geometry_mcse)),
          all(unresolved$geometry_status == "unresolved_geometry_precision"),
          all(is.finite(unresolved$benchmark_root_n_variance)))
geometry$component_precision[] <- "monte_carlo_estimate"
geometry$alpha <- 8
below <- env$wm_first_stage_benchmark("same_sample", 2L, 3L, geometry = geometry)$table
stopifnot(all(below$geometry_below_jensen), all(below$alpha == 8))
alpha_geometry <- env$wm_geometry_alpha(3L, 2L, draws = 64L, seed = 519L)
alpha_benchmark <- env$wm_first_stage_benchmark("same_sample", 2L, 3L,
                                               geometry = alpha_geometry)$table
equal(alpha_benchmark$alpha, rep(alpha_geometry$alpha, 6L))
equal(alpha_benchmark$benchmark_geometry_mcse,
      alpha_geometry$alpha_mcse * alpha_benchmark$alpha_coefficient)
equal(alpha_benchmark$geometry_integrand_upper_bound,
      rep(alpha_geometry$diagnostics$integrand_upper_bound, 6L))
bad_alpha <- alpha_geometry
bad_alpha$diagnostics$maximum_draw_value <- 2 * bad_alpha$diagnostics$integrand_upper_bound
expect_error(function() env$wm_first_stage_benchmark("same_sample", 2L, 3L,
                                                      geometry = bad_alpha),
             "bounded-integrand")
bad_alpha <- alpha_geometry
bad_alpha$d <- 3L
expect_error(function() env$wm_first_stage_benchmark("same_sample", 2L, 3L,
                                                      geometry = bad_alpha), "matching finite")
expect_error(function() env$wm_first_stage_benchmark("same_sample", 2L, 3L,
                                                      geometry = unclass(alpha_geometry)),
             "wm_geometry_alpha")

# Logistic sensitivities agree with finite differences of the normalized target.
weight <- env$wm_first_stage_weight_integrals()
equal(weight$information[1L, 2L], 0, 1e-11)
equal(weight$sensitivity$PATE[2L], 0, 1e-11)
equal(weight$sensitivity$PATT[2L], -1 / 24, 1e-11)
equal(weight$sensitivity$PATE[1L], weight$sensitivity$PATT[1L], 1e-11)
equal(weight$values[["inverse_q"]], weight$values[["inverse_control"]], 1e-11)
stopifnot(all(weight$absolute_errors < 1e-9), all(eigen(weight$information)$values > 0))
target <- function(parameter, estimand) {
  integrand <- function(t, outcome) {
    q0 <- stats::plogis(3 * t)
    qfit <- stats::plogis(parameter[1L] + parameter[2L] * t)
    factor <- if (estimand == "PATE")
      0.5 * q0 / qfit + 0.5 * (1 - q0) / (1 - qfit) else 0.5 * q0 / qfit
    factor * if (outcome) 1.5 + t else 1
  }
  numerator <- stats::integrate(function(t) integrand(t, TRUE), -0.5, 0.5,
                                 rel.tol = 1e-12)$value
  denominator <- stats::integrate(function(t) integrand(t, FALSE), -0.5, 0.5,
                                   rel.tol = 1e-12)$value
  numerator / denominator
}
for (estimand in c("PATE", "PATT")) {
  derivative <- vapply(1:2, function(k) {
    direction <- numeric(2)
    direction[k] <- 1e-5
    (target(c(0, 3) + direction, estimand) -
       target(c(0, 3) - direction, estimand)) / 2e-5
  }, numeric(1))
  equal(derivative, unname(weight$sensitivity[[estimand]]), 1e-8)
  equal(as.vector(weight$information %*% weight$covariance[[estimand]]),
        -unname(weight$sensitivity[[estimand]]))
}
wb <- env$wm_first_stage_benchmark("estimated_weights", 1L, 3L)$table
for (estimand in c("PATE", "PATT")) {
  rows <- wb[wb$estimand == estimand, ]
  known <- rows$benchmark_root_n_variance[rows$comparison == "known_first_stage"]
  fitted <- rows$benchmark_root_n_variance[rows$comparison == "fitted_adjusted"]
  equal(known - fitted, unname(weight$variance_reduction[estimand]))
}

# Shared-data pairing and total-n versus training-m uncertainty.
config <- env$wm_first_stage_pilot()
stopifnot(nrow(config) == 3L, all(config$replications == 40L))
config$n <- config$m <- 80L
config$M <- 1L
config$replications <- 3L
config$m[2L] <- 40L
timing_fields <- c("elapsed_seconds", "first_stage_seconds", "generation_seconds")
set.seed(9021)
saved <- .Random.seed
for (case in 1:3) {
  a <- env$wm_first_stage_run_rep(config[case, ], 1L, retain_fits = TRUE)
  b <- env$wm_first_stage_run_rep(config[case, ], 1L)
  stopifnot(identical(.Random.seed, saved), nrow(a$records) == 6L)
  fields <- setdiff(names(a$records), timing_fields)
  equal(a$records[fields], b$records[fields], 0)
  for (estimand in c("PATE", "PATT")) {
    naive_fit <- a$fits[[paste(estimand, "fitted_naive", sep = "/")]]
    adjusted_fit <- a$fits[[paste(estimand, "fitted_adjusted", sep = "/")]]
    if (!is.null(naive_fit)) {
      stopifnot(!is.null(adjusted_fit), identical(naive_fit$estimate, adjusted_fit$estimate),
                identical(naive_fit$graph, adjusted_fit$graph))
      if (case == 2L) {
        equal(adjusted_fit$training_adjustment$evaluation_training_ratio, 2)
        stopifnot(length(adjusted_fit$contributions$training) == 40L,
                  identical(adjusted_fit$contributions$row, naive_fit$contributions$row))
      }
    }
  }
  if (case == 2L)
    stopifnot(a$records$evaluation_rng_state[1L] != a$records$training_rng_state[1L])
}
# Direct fit guard and deliberately forced first-stage failure preservation.
dat <- env$wm_first_stage_generate(40L, 2L, "same_sample")
dat$Y <- dat$scores[, 1L] + dat$Z + 0.85 * dat$D
expect_error(function() env$wm_first_stage_score_fit(dat), "interior rule")
forced <- (function() {
  old <- env$wm_first_stage_generate
  on.exit(assign("wm_first_stage_generate", old, envir = env))
  assign("wm_first_stage_generate", function(...) {
    dat <- old(...)
    dat$Y <- dat$scores[, 1L] + dat$Z + 0.85 * dat$D
    dat
  }, envir = env)
  env$wm_first_stage_run_rep(config[1L, ], 1L)$records
})()
stopifnot(sum(forced$status == "ok") == 2L,
          all(forced$status[forced$comparison != "known_first_stage"] == "error"),
          all(forced$first_stage_status == "error"))

# CLI batching, failure capture, immutable outputs and companion provenance.
work <- tempfile("wdsm-first-stage-")
dir.create(work)
driver <- file.path(root, "simulations", "first_stage_run.R")
rscript <- file.path(R.home("bin"), "Rscript")
invoke <- function(args, expect_success = TRUE, pattern = NULL) {
  log <- tempfile("command-", tmpdir = work)
  status <- system2(rscript, c("--vanilla", shQuote(driver), shQuote(args)),
                     stdout = log, stderr = log, timeout = 120)
  output <- paste(readLines(log, warn = FALSE), collapse = "\n")
  if ((expect_success && status != 0L) || (!expect_success && status == 0L) ||
      (!is.null(pattern) && !grepl(pattern, output))) stop(output)
}
config_path <- file.path(work, "config.csv")
utils::write.csv(config, config_path, row.names = FALSE)
full <- file.path(work, "full.csv")
one <- file.path(work, "one.csv")
two <- file.path(work, "two.csv")
invoke(c(config_path, config$scenario[1L], 1, 2, full, root))
invoke(c(config_path, config$scenario[1L], 1, 1, one, root))
invoke(c(config_path, config$scenario[1L], 2, 2, two, root))
a <- utils::read.csv(full)
b <- rbind(utils::read.csv(one), utils::read.csv(two))
equal(a[setdiff(names(a), timing_fields)], b[setdiff(names(b), timing_fields)], 0)
stopifnot(nrow(a) == 12L)
manifest <- readRDS(paste0(full, ".metadata.rds"))
diagnostics <- readRDS(paste0(full, ".diagnostics.rds"))
stopifnot(identical(manifest$schema_version, "first_stage_v1"),
          manifest$records_written == 12L, length(diagnostics) == 2L,
          length(manifest$code_md5) > 3L, !is.null(manifest$completed))
invoke(c(config_path, config$scenario[1L], 1, 1, full, root), FALSE, "already exists")
orphan <- file.path(work, "orphan.csv")
saveRDS(list(partial = TRUE), paste0(orphan, ".diagnostics.rds"))
invoke(c(config_path, config$scenario[1L], 1, 1, orphan, root), FALSE, "already exists")
invoke(c(config_path, config$scenario[1L], 1.5, 2, file.path(work, "bad.csv"), root),
        FALSE, "integer")
bad <- config[1L, ]
bad$M <- bad$n + 1L
bad_path <- file.path(work, "failed-config.csv")
utils::write.csv(bad, bad_path, row.names = FALSE)
failed <- file.path(work, "failed.csv")
invoke(c(bad_path, bad$scenario, 1, 1, failed, root))
failure <- utils::read.csv(failed)
stopifnot(nrow(failure) == 6L, all(failure$status == "error"),
          all(is.na(failure$estimate)), all(is.na(failure$lower)))
pilot_path <- file.path(work, "pilot.csv")
invoke(c("--write-pilot", pilot_path))
equal(utils::read.csv(pilot_path), env$wm_first_stage_pilot())
invoke(c("--write-pilot", pilot_path), FALSE, "already exists")
cat("PASS: first-stage integrals, variance corrections, geometry precision, pairing, independent sample scaling, RNG restoration, CLI batching, provenance and failure retention.\n")
unlink(work, recursive = TRUE)
