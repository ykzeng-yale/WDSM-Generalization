# Unit-weight reduction, with expected values built without any package helper.
# Run: Rscript code-release/validation/check_unweighted_reductions.R OUTPUT_DIR [PREDECESSOR_DIR]
# OUTPUT_DIR must be new. PREDECESSOR_DIR optionally supplies original wdsmatch R/.
# These are deterministic algebra/implementation tests, not coverage experiments.
# Fitted-map cases check the realized graph and estimator; they do not certify
# inference that omits score-estimation uncertainty. No bootstrap SE equality
# with predecessor implementations is asserted.
a <- commandArgs(TRUE)
if (!length(a) || length(a) > 2L) stop("Supply a new output directory and optional predecessor source directory")
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
package <- normalizePath(file.path(dirname(script), ".."))
if (file.exists(a[1L])) stop("Output directory already exists")
dir.create(a[1L], recursive = TRUE)
out <- normalizePath(a[1L])

# Independent full-distance reference: deliberately uses exhaustive sorting,
# direct imputation and integer reuse counts, rather than package matching code.
reference <- function(Y, Z, S0, S1, M, estimand, mu0 = NULL, mu1 = NULL) {
  n <- length(Y)
  if (is.null(mu0)) mu0 <- rep(0, n)
  if (is.null(mu1)) mu1 <- rep(0, n)
  maps <- list(as.matrix(S0), as.matrix(S1))
  means <- list(mu0, mu1)
  queries <- if (estimand == "PATE") seq_len(n) else which(Z == 1L)
  neighbors <- vector("list", n)
  counts <- integer(n)
  imputed <- raw <- matrix(NA_real_, n, 2L)
  imputed[cbind(seq_len(n), Z + 1L)] <- Y
  raw[cbind(seq_len(n), Z + 1L)] <- Y
  edge_rows <- vector("list", length(queries))
  gaps <- numeric(length(queries))
  correction_sum <- 0
  for (k in seq_along(queries)) {
    i <- queries[k]
    arm <- 1L - Z[i]
    pool <- which(Z == arm)
    S <- maps[[arm + 1L]]
    squared <- vapply(pool, function(j) sum((S[i, ] - S[j, ])^2), numeric(1))
    ordering <- order(squared, pool)
    donors <- pool[ordering[seq_len(M)]]
    neighbors[[i]] <- donors
    for (j in donors) counts[j] <- counts[j] + 1L
    mu <- means[[arm + 1L]]
    raw[i, arm + 1L] <- sum(Y[donors]) / M
    imputed[i, arm + 1L] <- mu[i] + sum(Y[donors] - mu[donors]) / M
    correction_sum <- correction_sum + (1 - 2 * Z[i]) * (mu[i] - mean(mu[donors]))
    edge_rows[[k]] <- data.frame(query = i, donor = donors, arm = arm,
                                distance = sqrt(squared[ordering[seq_len(M)]]))
    gaps[k] <- if (length(pool) > M) squared[ordering[M + 1L]] - squared[ordering[M]] else Inf
  }
  denom <- if (estimand == "PATE") n else sum(Z)
  tau <- if (estimand == "PATE") mean(imputed[, 2L] - imputed[, 1L]) else
    mean(Y[Z == 1L] - imputed[Z == 1L, 1L])
  tau_raw <- if (estimand == "PATE") mean(raw[, 2L] - raw[, 1L]) else
    mean(Y[Z == 1L] - raw[Z == 1L, 1L])
  reuse <- counts / M
  coefficients <- if (estimand == "PATE") (2 * Z - 1) * (1 + reuse) / n else
    (Z - (1 - Z) * reuse) / sum(Z)
  row <- if (estimand == "PATE") {
    mu1 - mu0 - tau + (2 * Z - 1) * (1 + reuse) * (Y - ifelse(Z == 1L, mu1, mu0))
  } else {
    (Z * (Y - mu0 - tau) - (1 - Z) * reuse * (Y - mu0)) / mean(Z)
  }
  loads <- matrix(0, n, 2L)
  loads[cbind(seq_len(n), Z + 1L)] <- reuse
  list(estimate = tau, raw_estimate = tau_raw, correction = correction_sum / denom,
       neighbors = neighbors, edges = do.call(rbind, edge_rows), loads = loads,
       imputed = imputed, raw_imputed = raw, coefficients = coefficients,
       row = row, root_n_variance = mean(row^2), variance = sum(row^2) / n^2,
       minimum_boundary_squared_gap = min(gaps), denominator = denom)
}
standardize <- function(S) {
  S <- as.matrix(S)
  centered <- sweep(S, 2L, colMeans(S), "-")
  sweep(centered, 2L, sqrt(colMeans(centered^2)), "/")
}
quadratic_predictions <- function(S, Y, Z, arm) {
  S <- as.matrix(S)
  # Explicit intercept, linear coordinates, then each unordered quadratic term.
  design <- cbind(1, S)
  for (j in seq_len(ncol(S))) for (k in j:ncol(S)) design <- cbind(design, S[, j] * S[, k])
  fit <- stats::lm.fit(design[Z == arm, , drop = FALSE], Y[Z == arm])
  if (fit$rank != ncol(design) || anyNA(fit$coefficients)) stop("Reference quadratic is rank deficient")
  drop(design %*% fit$coefficients)
}

sut <- new.env(parent = globalenv())
source_files <- sort(list.files(file.path(package, "R"), "\\.R$", full.names = TRUE))
for (f in source_files) sys.source(f, sut)
predecessor <- NULL
if (length(a) == 2L) {
  predecessor <- new.env(parent = globalenv())
  old_files <- sort(list.files(file.path(normalizePath(a[2L]), "R"), "\\.R$", full.names = TRUE))
  if (!length(old_files)) stop("No predecessor R sources")
  for (f in old_files) sys.source(f, predecessor)
  source_files <- c(source_files, old_files)
}
write.csv(data.frame(path = c(script, source_files), md5 = unname(tools::md5sum(c(script, source_files)))),
          file.path(out, "source_md5.csv"), row.names = FALSE)

set.seed(930128)
n <- 180L
X <- as.data.frame(matrix(runif(4L * n, -1, 1), n, 4L))
names(X) <- paste0("x", seq_len(4L))
ps <- plogis(0.2 + 0.55 * X$x1 - 0.4 * X$x2)
Z <- as.integer(runif(n) < ps)
mu0 <- 0.3 + X$x2 + 0.7 * X$x3 - 0.3 * X$x4
mu1 <- 1.1 + 0.4 * X$x1 - 0.5 * X$x2 + 0.6 * X$x4
Y <- ifelse(Z == 1L, mu1, mu0) + 0.4 * (2 * rbinom(n, 1L, 0.5) - 1)
D <- data.frame(Y = Y, Z = Z, X)
ps_fit <- stats::glm(Z ~ x1 + x2 + x3 + x4, data = D,
  family = stats::quasibinomial(), weights = rep(1, n), start = rep(0, 5L),
  control = stats::glm.control(maxit = 100L, epsilon = 1e-10))
fitted_ps <- as.numeric(stats::predict(ps_fit, newdata = D, type = "response"))
fitted_pg <- sapply(0:1, function(z) {
  f <- stats::lm(Y ~ x1 + x2 + x3 + x4, data = D[Z == z, ], weights = rep(1, sum(Z == z)))
  as.numeric(stats::predict(f, newdata = D))
})
maps <- list(
  known_PSM = list(S0 = matrix(ps), S1 = matrix(ps), state = "known", family = "PSM"),
  known_DSM = list(S0 = standardize(cbind(ps, mu0)), S1 = standardize(cbind(ps, mu1)), state = "known", family = "DSM"),
  known_full = list(S0 = as.matrix(X), S1 = as.matrix(X), state = "known", family = "full"),
  fitted_PSM = list(S0 = matrix(fitted_ps), S1 = matrix(fitted_ps), state = "fitted", family = "PSM"),
  fitted_DSM = list(S0 = standardize(cbind(fitted_ps, fitted_pg[, 1L])), S1 = standardize(cbind(fitted_ps, fitted_pg[, 2L])), state = "fitted", family = "DSM"),
  fitted_scaling_full = list(S0 = standardize(X), S1 = standardize(X), state = "empirical_scaling", family = "full")
)
saveRDS(list(seed = 930128L, X = X, Y = Y, Z = Z, ps = ps, pg = cbind(mu0, mu1),
             fitted_ps = fitted_ps, fitted_pg = fitted_pg, maps = maps,
             rng_after_fixture = .Random.seed), file.path(out, "fixture.rds"))

checks <- list()
append_check <- function(case, check, pass, error = NA_real_, detail = "") {
  checks[[length(checks) + 1L]] <<- data.frame(case = case, check = check, pass = pass,
                                            maximum_absolute_error = error, detail = detail)
}
near <- function(case, name, observed, expected, tolerance = 2e-11) {
  same_missing <- identical(is.na(observed), is.na(expected))
  ok_shape <- identical(dim(observed), dim(expected)) && length(observed) == length(expected)
  diff <- if (same_missing && ok_shape) max(c(0, abs(observed - expected)), na.rm = TRUE) else Inf
  bound <- tolerance * max(c(1, abs(expected)), na.rm = TRUE)
  append_check(case, name, same_missing && ok_shape && is.finite(diff) && diff <= bound, diff)
}
exact <- function(case, name, observed, expected) append_check(case, name, identical(observed, expected))
run_case <- function(id, expr) tryCatch(force(expr), error = function(e)
  append_check(id, "execution", FALSE, detail = conditionMessage(e)))
compare_fit <- function(id, got, expected, has_variance) {
  for (field in c("estimate", "raw_estimate", "correction", "denominator")) near(id, field, got[[field]], expected[[field]])
  exact(id, "donor_neighbors", got$graph$neighbors, expected$neighbors)
  for (field in c("query", "donor", "arm")) exact(id, paste0("edge_", field), got$graph$edges[[field]], expected$edges[[field]])
  near(id, "edge_distance", got$graph$edges$distance, expected$edges$distance)
  near(id, "edge_share", got$graph$edges$share, rep(1 / got$M, nrow(expected$edges)))
  near(id, "outcome_share", got$graph$edges$outcome_share, rep(1 / got$M, nrow(expected$edges)))
  near(id, "incoming_K_over_M", unname(got$loads$incoming), expected$loads)
  near(id, "direct_imputation", unname(got$imputed), expected$imputed)
  near(id, "raw_imputation", unname(got$raw_imputed), expected$raw_imputed)
  near(id, "raw_outcome_coefficients", sum(expected$coefficients * got$data$Y), got$raw_estimate)
  near(id, "normalized_row_contributions", got$contributions$row, expected$row)
  near(id, "actual_contributions", got$contributions$actual, expected$row)
  near(id, "zero_edge_contribution", got$contributions$edge, rep(0, length(got$contributions$edge)))
  if (has_variance) {
    near(id, "root_n_variance_arithmetic", got$root_n_variance, expected$root_n_variance)
    near(id, "sampling_variance_arithmetic", got$variance, expected$variance)
    near(id, "standard_error_arithmetic", got$se, sqrt(expected$variance))
  } else exact(id, "raw_no_inference", is.na(got$variance), TRUE)
}
references <- list()
for (map_name in names(maps)) for (estimand in c("PATE", "PATT")) for (M in c(1L, 3L)) {
  mapping <- maps[[map_name]]
  for (correction in c("raw", "oracle_BC", "quadratic_BC")) {
    id <- paste(map_name, estimand, M, correction, sep = "/")
    run_case(id, {
      m0 <- if (correction == "raw") NULL else if (correction == "oracle_BC") mu0 else quadratic_predictions(mapping$S0, Y, Z, 0L)
      m1 <- if (correction == "raw") NULL else if (correction == "oracle_BC") mu1 else quadratic_predictions(mapping$S1, Y, Z, 1L)
      ref <- reference(Y, Z, mapping$S0, mapping$S1, M, estimand, m0, m1)
      references[[id]] <- ref
      args <- list(Y = Y, Z = Z, weights = rep(1, n), scores0 = mapping$S0,
                   M = M, estimand = estimand, mean0 = m0, variance = correction != "raw")
      if (estimand == "PATE") args <- c(args, list(scores1 = mapping$S1, mean1 = m1))
      got <- do.call(sut$wm_match, args)
      compare_fit(id, got, ref, correction != "raw")
      if (estimand == "PATT" || mapping$family != "DSM") {
        args$method <- "stabilized"
        args$rho0 <- rep(1, n)
        if (estimand == "PATE") args$rho1 <- rep(1, n)
        stabilized <- do.call(sut$wm_match, args)
        compare_fit(paste0(id, "/stabilized"), stabilized, ref, correction != "raw")
      }
      if (correction == "quadratic_BC") {
        fit_args <- list(Y = Y, Z = Z, weights = rep(1, n), scores0 = mapping$S0,
                         M = M, estimand = estimand, degree = 2L, variance = TRUE)
        if (estimand == "PATE") fit_args$scores1 <- mapping$S1
        fitted <- do.call(sut$wm_fit, fit_args)
        compare_fit(paste0(id, "/wm_fit"), fitted, ref, TRUE)
        near(paste0(id, "/wm_fit"), "independent_mean0", fitted$predictions$mean0, m0)
        if (estimand == "PATE") near(paste0(id, "/wm_fit"), "independent_mean1", fitted$predictions$mean1, m1)
      }
    })
  }
}
# Exact ties: the reference explicitly encodes the core's original-row priority.
for (estimand in c("PATE", "PATT")) for (M in c(1L, 3L)) {
  id <- paste("exact_tie", estimand, M, sep = "/")
  run_case(id, {
    ty <- seq_len(12L) / 3
    tz <- rep(0:1, each = 6L)
    ts <- matrix(0, 12L, 2L)
    got <- sut$wm_match(ty, tz, rep(1, 12L), ts, M = M, estimand = estimand, variance = FALSE)
    compare_fit(id, got, reference(ty, tz, ts, ts, M, estimand), FALSE)
  })
}
run_case("distinct_map_stabilized_PATE_boundary", {
  message <- tryCatch({sut$wm_match(Y, Z, rep(1, n), maps$known_DSM$S0,
    maps$known_DSM$S1, mean0 = mu0, mean1 = mu1, rho0 = rep(1, n),
    rho1 = rep(1, n), method = "stabilized"); "NO ERROR"}, error = conditionMessage)
  append_check("distinct_map_stabilized_PATE_boundary", "explicit_scope_rejection",
               grepl("requires the same supplied score matrix", message, fixed = TRUE), detail = message)
})

# Only the independent expectations above are reused for predecessor public calls.
# Its multinomial bootstrap is intentionally disabled: SE conventions differ.
if (!is.null(predecessor)) for (state in c("known", "fitted")) for (estimand in c("PATE", "PATT")) for (M in c(1L, 3L)) for (bc in c(FALSE, TRUE)) {
  id <- paste("predecessor", state, estimand, M, if (bc) "quadratic_BC" else "raw", sep = "/")
  run_case(id, {
    args <- list(Y = Y, X = X, Z = Z, weights = rep(1, n), M = M,
                 use.bias.correction = bc, varest = FALSE)
    if (state == "known") args <- c(args, list(ps = ps, pg = cbind(mu0, mu1))) else
      args <- c(args, list(model.ps = Z ~ x1 + x2 + x3 + x4, model.pg = Y ~ x1 + x2 + x3 + x4))
    public <- if (estimand == "PATE") predecessor$wdsmatchATE else predecessor$wdsmatchATT
    old <- do.call(public, args)
    ties_inactive <- all(vapply(old$diagnostics$ties$arms, function(x)
      x$boundary_tie_recipients == 0L && x$randomized_recipients == 0L, logical(1)))
    exact(id, "legacy_boundary_ties_inactive", ties_inactive, TRUE)
    ref_id <- paste(paste0(state, "_DSM"), estimand, M, if (bc) "quadratic_BC" else "raw", sep = "/")
    near(id, "public_point_estimate", old$estimate, references[[ref_id]]$estimate)
  })
}
saveRDS(references, file.path(out, "independent_references.rds"))
results <- do.call(rbind, checks)
write.csv(results, file.path(out, "checks.csv"), row.names = FALSE)
write.csv(data.frame(case = names(references), minimum_boundary_squared_gap = vapply(references,
  function(x) x$minimum_boundary_squared_gap, numeric(1))), file.path(out, "boundary_gaps.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "session.txt"))
summary <- list(passed = all(results$pass), checks = nrow(results), cases = length(unique(results$case)),
                failed_checks = sum(!results$pass), seed = 930128L, n = n, n0 = sum(Z == 0L), n1 = sum(Z == 1L),
                predecessor_public_executed = !is.null(predecessor),
                statistical_scope = "Finite-sample graph/point/contribution algebra; no coverage or general fitted-map variance claim")
saveRDS(summary, file.path(out, "summary.rds"))
cat(if (summary$passed) "PASS" else "FAIL", ":", summary$checks, "checks in", summary$cases, "cases;",
    summary$failed_checks, "failed checks.\n")
if (!summary$passed) print(results[!results$pass, ], row.names = FALSE)
if (!summary$passed) quit(status = 1L)
