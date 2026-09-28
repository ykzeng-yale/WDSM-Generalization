script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
for (file in sort(list.files(file.path(root, "R"), "\\.R$", full.names = TRUE))) source(file)
for (file in c("first_stage.R", "first_stage_benchmarks.R", "benchmarks.R", "supplement.R"))
  source(file.path(root, "simulations", file))
equal <- function(x, y) stopifnot(isTRUE(all.equal(x, y, tolerance = 1e-10, check.attributes = FALSE)))
config <- wm_supplement_config()
stopifnot(nrow(config) == 8L, !anyDuplicated(config$seed),
  sum(ifelse(config$branch == "split", 3, 8)*5) == 270L)
set.seed(651); old <- .Random.seed
for (i in c(1L, 3L, 5L, 7L, 8L)) {
  result <- wm_supplement_run_rep(config[i, ], 1L)
  stopifnot(identical(old, .Random.seed))
  x <- result$records
  stopifnot(all(x$status == "ok"), !anyDuplicated(x[c("method", "estimand", "correction")]))
  equal(x$root_n_variance, x$analysis_n*x$sampling_variance)
  if (config$branch[i] == "split") {
    a <- x[x$estimand == "PATE", ]; b <- x[x$estimand == "potential0", ]; c <- x[x$estimand == "potential1", ]
    equal(a$estimate, c$estimate-b$estimate)
    equal(a$sampling_variance, b$sampling_variance+c$sampling_variance)
    stopifnot(b$analysis_n+c$analysis_n == a$analysis_n)
  } else for (entry in result$diagnostics) if (!is.null(entry$nuisance)) {
    for (model in entry$nuisance$models) {
      stopifnot(model$mean$rank == length(model$mean$coefficients))
      if (config$branch[i] == "spline") {
        stopifnot(!length(intersect(model$training_rows, model$evaluation_rows)))
        stopifnot(all(entry$graph_restrictions$folds[model$training_rows] != model$fold))
      }
      if (config$branch[i] == "strata")
        stopifnot(all(entry$graph_restrictions$strata[model$training_rows] == model$stratum))
    }
  }
}
stratified <- wm_supplement_benchmarks(config[5, ])
equal(stratified$target[stratified$estimand == "PATE"], rep(393/281, 2))
equal(stratified$target[stratified$estimand == "PATT"], rep(41/27, 2))
stopifnot(all(stratified$geometry_mcse == 0))
missing <- wm_supplement_benchmarks(config[3, ])
stopifnot(all(is.na(missing$root_n_variance)))
# One-direction coefficient follows from raw mark enumeration, not a full-PATE reuse.
g <- wm_geometry_alpha(M = 3, d = 2, draws = 1000, seed = 652)
split <- wm_supplement_benchmarks(config[7, ], g)
rho <- c(1, 2); second <- c(1.25, 5); t <- c(37/32, 7/4); p <- c(.6, .4)
V <- vapply(1:2, function(z) {
  j <- 3-z; ratio <- p[j]/p[z]; alpha <- if (z == 1) 10.5 else g$alpha
  (sum(p*second)*c(1/12, 1/6)[z]+p[z]*t[z]*(rho[z]^2+
    2*ratio*rho[j]*rho[z]+ratio*second[j]/3+alpha*ratio^2*rho[j]^2/9))/sum(p*rho)^2
}, numeric(1))
equal(split$root_n_variance, c(2*sum(V), V))
alpha_only <- wm_supplement_benchmarks(config[3, ], g)
stopifnot(all(is.na(alpha_only$root_n_variance[alpha_only$method == "self_normalized"])),
  all(is.finite(alpha_only$root_n_variance[alpha_only$method == "stabilized"])))
wrong <- g; wrong$d <- 3L
stopifnot(inherits(try(wm_supplement_benchmarks(config[3, ], wrong), silent = TRUE), "try-error"))
cat("PASS: supplemental RNG restoration, complete paired records, support/rank/fold/stratum guards, unequal block scaling and analytic mark benchmarks.\n")
