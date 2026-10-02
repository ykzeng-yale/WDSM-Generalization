# Run from code-release. Deterministic formulas and aligned external calls;
# this is not a Monte Carlo coverage or bootstrap-calibration experiment.
source("simulations/ai06_variance_reference.R")
checks <- 0L
same <- function(actual, expected, label, tolerance = 1e-10) {
  if (!isTRUE(all.equal(actual, expected, tolerance = tolerance,
                       check.attributes = FALSE))) stop(label)
  checks <<- checks + 1L
}
reject <- function(expr, label) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  if (is.null(error)) stop(label)
  checks <<- checks + 1L
}
Y <- c(0, 2, 1, 7, 4, 8)
Z <- rep(0:1, 3)
S <- matrix(c(0, 1, 3, 5, 6, 9), ncol = 1)
a <- ai06_variance_reference(Y, Z, S, estimand = "PATE")
same(a$donors[, 1], c(2, 1, 2, 5, 4, 5), "ATE hand donor list")
same(a$reuse_count, c(1, 2, 0, 1, 2, 0), "ATE integer reuse")
same(a$auxiliary_donors[, 1], c(3, 4, 1, 2, 3, 4), "Other same-arm donors")
same(a$conditional_outcome_variance, c(.5, 12.5, .5, 12.5, 4.5, .5), "AI06 J=1 factors")
same(a$raw_contrast, c(2, 2, 1, 3, 3, 4), "ATE hand contrasts")
same(a$estimate, 5 / 2, "ATE estimate")
same(a$conditional_root_n_variance, 103 / 3, "ATE conditional root-n variance")
same(a$marginal_root_n_variance, 89 / 4, "ATE marginal root-n variance")
same(a$marginal_variance, 89 / 24, "ATE estimator scale")
t <- ai06_variance_reference(Y, Z, S, estimand = "PATT")
same(t$donors[, 1], c(1, 5, 5), "ATT hand donor list")
same(t$estimate, 3, "ATT estimate")
same(t$conditional_root_target_count_variance, 44 / 3, "ATT source scale conditional")
same(t$marginal_root_target_count_variance, 11 / 3, "ATT source scale marginal")
same(t$conditional_root_n_variance, 88 / 3, "ATT converted root-n conditional")
same(t$marginal_variance, 11 / 9, "ATT estimator scale")

# J changes the auxiliary variance graph, not the M-donor point estimator.
j2 <- ai06_variance_reference(Y, Z, S, J = 2)
same(j2$estimate, a$estimate, "J must not change the point estimate")
same(j2$donors, a$donors, "J must not change the imputation graph")
same(j2$conditional_outcome_variance,
     (2 / 3) * c(6.25, 30.25, 1, 4, 12.25, 12.25), "AI06 J=2 squared-average formula")
shift <- ai06_variance_reference(Y + 100, Z, S)
same(shift$marginal_variance, a$marginal_variance, "Outcome-location invariance")
scale <- ai06_variance_reference(-3 * Y, Z, S)
same(scale$estimate, -3 * a$estimate, "Outcome-scale point equivariance")
same(scale$marginal_variance, 9 * a$marginal_variance, "Outcome-scale variance equivariance")
swap <- ai06_variance_reference(Y, 1 - Z, S)
same(swap$estimate, -a$estimate, "Treatment relabeling ATE point")
same(swap$marginal_variance, a$marginal_variance, "Treatment relabeling ATE variance")
reject(ai06_variance_reference(Y, Z, S, J = 3), "Insufficient auxiliary arm accepted")
reject(ai06_variance_reference(Y, Z, S, M = 4), "Insufficient opposite arm accepted")
reject(ai06_variance_reference(Y, Z, S, J = 1.5), "Noninteger J accepted")
reject(ai06_variance_reference(Y, Z, S, weights = rep(1, 6)), "Undeclared weights accepted")

# Review regression: squared raw distances silently collapsed tiny coordinates.
tiny_Y <- c(1, 4, -2, 3, 8, 2, 7, 1, 9)
tiny_Z <- c(rep(0, 5), rep(1, 4))
tiny_S <- matrix(c(0, 3, 7, 12, 20, 1, 5, 10, 18), ncol = 1)
ordinary <- ai06_variance_reference(tiny_Y, tiny_Z, tiny_S, M = 2, J = 2)
same(ordinary$estimate, 37 / 18, "Independent nine-row reference point")
for (factor in c(1e-200, 1e200)) {
  rescaled <- ai06_variance_reference(tiny_Y, tiny_Z, factor * tiny_S, M = 2, J = 2)
  # Floating-point scaling may order an internal tie differently; both members
  # remain selected. Compare donor sets and the estimator, not that internal order.
  same(t(apply(rescaled$donors, 1, sort)), t(apply(ordinary$donors, 1, sort)),
       "Extreme common score scaling must preserve donor sets")
  same(rescaled$estimate, ordinary$estimate, "Extreme scaling point estimate")
  same(rescaled$auxiliary_donors, ordinary$auxiliary_donors, "Extreme scaling auxiliary graph")
  same(rescaled$marginal_variance, ordinary$marginal_variance, "Extreme scaling variance")
}

stopifnot(requireNamespace("Matching", quietly = TRUE),
          requireNamespace("wdsmatch", quietly = TRUE))
external_cases <- 0L
max_point_error <- max_variance_error <- 0
index <- seq_len(48)
all_scores <- outer(index, sqrt(c(2, 3, 5, 7, 11)), function(i, x) sin(i * x))
Z <- rep(0:1, 24)
Y <- .2 * all_scores[, 1] + Z * (1 + .4 * all_scores[, 1]) + cos(index * .71)
for (dimension in c(1L, 2L, 5L)) for (M in c(1L, 3L)) {
  S <- all_scores[, seq_len(dimension), drop = FALSE]
  for (estimand in c("PATE", "PATT")) {
    ref <- ai06_variance_reference(Y, Z, S, M = M, J = 1L, estimand = estimand)
    wm <- wdsmatch::wm_match(Y, Z, rep(1, length(Y)), S,
      M = M, estimand = estimand, variance = FALSE)
    same(wm$estimate, ref$estimate, "WM raw point reduction")
    wm_donors <- matrix(wm$graph$edges$donor, ncol = M, byrow = TRUE)
    same(wm_donors, ref$donors, "WM exact donor ordering")
    for (conditional in c(FALSE, TRUE)) {
      native <- Matching::Match(Y = Y, Tr = Z, X = S,
        estimand = if (estimand == "PATE") "ATE" else "ATT", M = M,
        replace = TRUE, ties = TRUE, distance.tolerance = 0,
        Weight = 3, Weight.matrix = diag(apply(S, 2, var), dimension),
        BiasAdjust = FALSE, Var.calc = 1L, sample = conditional)
      expected <- if (conditional) ref$conditional_variance else ref$marginal_variance
      same(native$version, "standard", "Matching must retain variance calculation")
      same(as.numeric(native$weights), rep(1 / M, length(ref$query) * M),
           "Native graph must retain exactly M equal-weight donors per query")
      native_pairs <- sort(paste(native$MatchLoopC[, 1], native$MatchLoopC[, 2], sep = ":"))
      reference_pairs <- sort(paste(rep(ref$query, each = M), as.vector(t(ref$donors)), sep = ":"))
      same(native_pairs, reference_pairs, "Native exact directed donor sets")
      same(as.numeric(native$est), ref$estimate, "Matching raw point")
      same(as.numeric(native$se)^2, expected, "Matching aligned AI variance", 1e-8)
      max_point_error <- max(max_point_error, abs(as.numeric(native$est) - ref$estimate))
      max_variance_error <- max(max_variance_error, abs(as.numeric(native$se)^2 - expected))
      external_cases <- external_cases + 1L
    }
  }
}
result <- list(passed = TRUE, deterministic_checks = checks,
  external_variance_configurations = external_cases,
  matching_version = as.character(utils::packageVersion("Matching")),
  wdsmatch_version = as.character(utils::packageVersion("wdsmatch")),
  native_variance_settings = "standard version, ties=TRUE, no boundary ties, Var.calc=1",
  max_absolute_point_error = max_point_error,
  max_absolute_variance_error = max_variance_error,
  simulation_replications = 0L)
print(result)
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  if (length(args) != 1L) stop("At most one result JSON path is accepted.")
  jsonlite::write_json(result, args[[1]], auto_unbox = TRUE, pretty = TRUE, digits = 17)
}
