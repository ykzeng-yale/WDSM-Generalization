# Deterministic base-R fixtures only: no estimator, simulation or geometry run.
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = TRUE)
dir <- file.path(dirname(script), "..", "simulations")
for (file in c("first_stage.R", "first_stage_benchmarks.R", "first_stage_summarize.R", "first_stage_join.R")) source(file.path(dir, file))
equal <- function(x, y, tolerance = 1e-10) {
  check <- all.equal(x, y, tolerance = tolerance, check.attributes = FALSE)
  if (!isTRUE(check)) stop(paste(check, collapse = "; "))
}
fails <- function(expr, pattern) {
  error <- tryCatch({force(expr); NULL}, error = conditionMessage)
  if (is.null(error) || !grepl(pattern, error)) stop("expected error containing: ", pattern, "; got: ", error)
}
fixture <- function() {
  x <- expand.grid(replication = 1:4, estimand = c("PATE", "PATT"),
    comparison = c("known_first_stage", "fitted_naive", "fitted_adjusted"), stringsAsFactors = FALSE)
  x$schema_version <- 1L; x$scenario <- "manual_independent"; x$branch <- "independent_training"
  x$n <- 10L; x$m <- 20L; x$d <- 1L; x$M <- 1L; x$base_seed <- 12L
  x$pairing_id <- paste(x$scenario, x$replication, x$estimand, sep = "/")
  x$target <- 1; x$status <- "ok"; x$first_stage_status <- "ok"
  x$inference_status <- ifelse(x$comparison == "fitted_naive", "first_stage_omission_diagnostic", "within_declared_parametric_scope")
  x$estimate <- 1 + c(-.3, -.1, .1, .3)[x$replication]
  x$root_n_variance <- 1 + x$replication / 10 + .2 * (x$comparison == "fitted_adjusted")
  x$sampling_variance <- x$root_n_variance / x$n
  x$lower <- x$estimate - qnorm(.975) * sqrt(x$sampling_variance)
  x$upper <- x$estimate + qnorm(.975) * sqrt(x$sampling_variance)
  x$evaluation_rng_state <- paste0("evaluation/", x$replication)
  x$training_rng_state <- paste0("training/", x$replication)
  x$parameter1 <- .01; x$parameter2 <- NA_real_; x$sigma2 <- 1/16; x$r <- .6
  x$evaluation_training_ratio <- .5; x$warning <- ""; x$error <- ""
  x
}
x <- fixture(); s <- wm_first_stage_summarize(x)
stopifnot(nrow(s) == 6L, all(s$requested == 4), all(s$independent_replications == 4), all(s$records_per_replication == 6))
equal(s$bias, rep(0, 6)); equal(s$bias_mcse, rep(sqrt(1/60), 6))
equal(s$empirical_root_n_variance, rep(2/3, 6))
equal(s$mean_root_n_variance[s$comparison == "fitted_naive"], rep(1.25, 2))
equal(s$adjusted_minus_naive_root_n_variance, rep(.2, 6))
equal(s$adjusted_minus_naive_root_n_variance_mcse, rep(0, 6))
stopifnot(all(s$coverage == 1), all(s$coverage_mcse == 0))
equal(s$coverage_exact95_lower, rep(.025^(1/4), 6)); equal(s$coverage_exact95_upper, rep(1, 6))
# Ratio IF uses the paired per-replication fit variance and squared centered point.
b <- s[s$estimand == "PATE" & s$comparison == "fitted_naive", ]
a <- x[x$estimand == "PATE" & x$comparison == "fitted_naive", ]
v <- 1/15; vbar <- .125; square <- c(.09, .01, .01, .09) * 4/3
expected_if <- (a$sampling_variance-vbar)/v-vbar/v^2*(square-v)
equal(b$variance_ratio_mcse_delta, sd(expected_if)/2)
# One adjusted-only failure changes the joint-success count, not naive results.
y <- x; bad <- y$replication == 1 & y$estimand == "PATE" & y$comparison == "fitted_adjusted"
y$status[bad] <- "error"; y$error[bad] <- "intentional adjustment failure"
for (name in c("estimate", "root_n_variance", "sampling_variance", "lower", "upper")) y[[name]][bad] <- NA_real_
f <- wm_first_stage_summarize(y)
stopifnot(f$completed[f$estimand == "PATE" & f$comparison == "fitted_adjusted"] == 3,
  all(f$paired_variance_completed[f$estimand == "PATE"] == 3), all(f$paired_variance_failed[f$estimand == "PATE"] == 1))
fails(wm_first_stage_summarize(x[-1, ]), "incomplete")
fails(wm_first_stage_summarize(rbind(x, x[1, ])), "duplicate")
y <- x
for (field in c("estimate", "lower", "upper")) y[[field]][y$comparison == "fitted_adjusted"] <- y[[field]][y$comparison == "fitted_adjusted"] + .001
fails(wm_first_stage_summarize(y), "point identity")
y <- x; y$training_rng_state[1] <- "corrupted"; fails(wm_first_stage_summarize(y), "RNG states")
y <- x; y$n[1] <- 11; fails(wm_first_stage_summarize(y), "configuration")
y <- x; y$root_n_variance[1] <- Inf; fails(wm_first_stage_summarize(y), "numeric")
y <- x; y$status[1] <- "error"; fails(wm_first_stage_summarize(y), "failed record contains")
# Same records split into complete declared batches reproduce the full summary.
tmp <- tempfile("first-stage-summary-"); dir.create(tmp)
cfg <- data.frame(scenario = "manual_independent", branch = "independent_training", n = 10L, m = 20L,
  d = 1L, M = 1L, seed = 12L, replications = 4L)
hash <- function(name, letter = "a") setNames(paste(rep(letter, 32), collapse = ""), name)
write_batch <- function(path, rows, first, last) {
  utils::write.csv(rows, path, row.names = FALSE)
  metadata <- list(schema_version = "first_stage_v1", config = cfg, first = first, last = last,
    requested_records = 6L*(last-first+1L), records_written = nrow(rows), columns = names(rows),
    last_completed_replication = last, completed = as.POSIXct("2026-09-28", tz = "UTC"),
    code_md5 = hash("/remote/a.R"), design_md5 = hash("/remote/design.md", "b"),
    config_file_md5 = hash("/remote/config.csv", "c"))
  saveRDS(metadata, paste0(path, ".metadata.rds")); path
}
p1 <- write_batch(file.path(tmp, "first.csv"), x[x$replication <= 2, ], 1L, 2L)
p2 <- write_batch(file.path(tmp, "second.csv"), x[x$replication >= 3, ], 3L, 4L)
r <- wm_first_stage_read_batches(c(p1,p2)); split_summary <- wm_first_stage_summarize(r)
order_rows <- function(z) z[order(z$scenario, z$estimand, z$comparison), ]
equal(order_rows(split_summary), order_rows(s))
fails(wm_first_stage_read_batches(c(p1,p1)), "duplicate")
fails(wm_first_stage_read_batches(p1), "planned replications")
stopifnot(nrow(wm_first_stage_read_batches(p1, require_complete = FALSE)) == 12L)
m <- readRDS(paste0(p2, ".metadata.rds")); old <- m
m$completed <- NULL; saveRDS(m, paste0(p2, ".metadata.rds"))
fails(wm_first_stage_read_batches(c(p1,p2)), "incomplete")
m <- old; m$code_md5 <- hash("/remote/a.R", "d"); saveRDS(m, paste0(p2, ".metadata.rds"))
fails(wm_first_stage_read_batches(c(p1,p2)), "hashes differ")
m <- old; m$config$n <- 11L; saveRDS(m, paste0(p2, ".metadata.rds"))
fails(wm_first_stage_read_batches(c(p1,p2)), "configuration mismatch")
m <- old; m$records_written <- 11L; saveRDS(m, paste0(p2, ".metadata.rds"))
fails(wm_first_stage_read_batches(c(p1,p2)), "record grid")
saveRDS(old, paste0(p2, ".metadata.rds"))
# Exact d=1 benchmark: V0,A=11/32; V0,T=7/16, training addition=(n/m)*.075=.0375.
j <- wm_first_stage_join(s)
equal(j$benchmark_root_n_variance[j$estimand == "PATE" & j$comparison == "known_first_stage"], 11/32)
equal(j$benchmark_root_n_variance[j$estimand == "PATT" & j$comparison == "known_first_stage"], 7/16)
equal(j$benchmark_root_n_variance[j$estimand == "PATE" & j$comparison == "fitted_naive"], 11/32+.0375)
equal(j$reported_root_n_variance_limit[j$estimand == "PATE" & j$comparison == "fitted_naive"], 11/32)
equal(j$adjusted_minus_naive_reported_limit, rep(.0375, 6))
stopifnot(all(j$geometry_status == "exact_1d"), all(j$benchmark_geometry_mcse == 0))
equal(j[names(s)], s)
# Missing and unresolved high-dimensional geometry never fabricates precision.
s2 <- s; s2$d <- 2L; missing <- wm_first_stage_join(s2)
stopifnot(all(is.na(missing$benchmark_root_n_variance)), all(is.na(missing$empirical_discrepancy_z)))
equal(missing$adjusted_minus_naive_reported_limit, rep(.0375, 6))
g <- structure(list(M = 1L, d = 2L, alpha = 1.6, alpha_mcse = .02,
  precision_status = "monte_carlo_estimate", component_precision = "monte_carlo_estimate",
  draws = 1000L, seed = 1L, diagnostics = list(nonzero_draws = 800L)), class = c("wm_geometry", "list"))
j2 <- wm_first_stage_join(s2, g)
equal(j2$benchmark_geometry_mcse[j2$estimand == "PATE"], rep(.02/16, 3))
equal(j2$empirical_discrepancy_combined_mcse[j2$estimand == "PATE"],
  sqrt(s2$empirical_root_n_variance_mcse_asymptotic[s2$estimand == "PATE"]^2 + (.02/16)^2))
g$alpha_mcse <- 0; unresolved <- wm_first_stage_join(s2, g)
stopifnot(all(is.na(unresolved$benchmark_geometry_mcse)), all(is.na(unresolved$empirical_discrepancy_z)))
failed <- wm_first_stage_join(s2, setNames(list(NULL), "M1_d2"))
stopifnot(all(failed$geometry_status == "geometry_run_failed"))
y <- s; y$target <- 2; fails(wm_first_stage_join(y), "target")
y <- s; y$failed[1] <- 1; fails(wm_first_stage_join(y), "failure counts")
fails(wm_first_stage_join(s[-1, ]), "comparison grid")
# CSV import must preserve all-empty same-sample training/error columns.
y <- x; y$scenario <- "manual_same"; y$branch <- "same_sample"
y$m <- y$n; y$d <- 2L; y$training_rng_state <- ""; y$evaluation_training_ratio <- NA_real_
y$pairing_id <- paste(y$scenario, y$replication, y$estimand, sep = "/")
cfg$scenario <- "manual_same"; cfg$branch <- "same_sample"; cfg$m <- cfg$n; cfg$d <- 2L
p3 <- write_batch(file.path(tmp, "same.csv"), y, 1L, 4L)
same <- wm_first_stage_read_batches(p3)
stopifnot(all(same$training_rng_state == ""), all(same$error == ""), all(same$warning == ""))
stopifnot(nrow(wm_first_stage_summarize(same)) == 6L)
unlink(tmp, recursive = TRUE)
cat("First-stage summary and benchmark-join deterministic fixtures passed.\n")
