script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
for (file in sort(list.files(file.path(root, "R"), "\\.R$", full.names = TRUE))) source(file)
for (file in c("first_stage.R", "first_stage_benchmarks.R", "benchmarks.R", "first_stage_summarize.R"))
  source(file.path(root, "simulations", file))
cfg <- data.frame(scenario = "gaussian_check", branch = "gaussian_same_sample", n = 200L,
  m = 200L, d = 1L, M = 3L, seed = 653L, replications = 2L)
set.seed(654); old <- .Random.seed
x <- do.call(rbind, lapply(1:2, function(i) wm_first_stage_run_rep(cfg, i)$records))
stopifnot(identical(old, .Random.seed), nrow(x) == 12L, all(x$status == "ok"))
wm_first_stage_validate_records(x)
benchmark <- wm_first_stage_benchmark("gaussian_same_sample", 1L, 3L)$table
independent <- wm_benchmark_strong(3L, wm_beta_1d(3L), q = .5,
  weight0 = c(1,1), weight1 = c(2,2), variance0 = c(1,1)/16,
  variance1 = c(1,1)/16, effect_variance = 0, effect_mean = 1)$table
for (estimand in c("PATE", "PATT")) {
  expected <- independent$root_n_variance[independent$method == "self_normalized" & independent$estimand == estimand]
  b <- benchmark[benchmark$estimand == estimand, ]
  stopifnot(abs(b$benchmark_root_n_variance[b$comparison == "known_first_stage"]-expected) < 1e-12,
    abs(b$benchmark_root_n_variance[b$comparison == "fitted_adjusted"]-(expected-.075)) < 1e-12,
    abs(b$reported_root_n_variance_limit[b$comparison == "fitted_naive"]-expected) < 1e-12)
}
summary <- wm_first_stage_summarize(x)
stopifnot(nrow(summary) == 6L, all(summary$completed == 2L))
cat("PASS: Gaussian actual fitted-score constructor, paired graph/point identity, RNG restoration, independent mark-variance benchmark and first-stage covariance reduction.\n")
