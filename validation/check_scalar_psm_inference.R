# Deterministic saved-array validation; no fitting, DGP, RNG or new simulation.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1L] else ".", mustWork = TRUE)
for (f in list.files(file.path(root, "R"), full.names = TRUE)) source(f)
labels <- character()
check <- function(ok, label) {
  if (!isTRUE(ok)) stop(label, call. = FALSE)
  labels <<- c(labels, label)
}
near <- function(x, y, label, tolerance = 2e-10) {
  check(identical(dim(x), dim(y)) && length(x) == length(y) &&
    all(is.finite(c(x, y))) &&
    all(abs(x - y) <= tolerance * pmax(1, abs(x), abs(y))), label)
}
fails <- function(expr, label) check(inherits(try(expr, silent = TRUE), "try-error"), label)
had_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
if (had_rng) rng_before <- get(".Random.seed", envir = .GlobalEnv)
forbidden <- function(...) stop("Forbidden fitting, rematching or RNG call.")
wdsm_fit_ps <- wm_match <- .wm_scalar_nw <- forbidden
rnorm <- runif <- sample <- rmultinom <- forbidden
start <- proc.time()[["elapsed"]]
fixture <- readRDS(file.path(root, "validation", "fixtures", "scalar_replication_24.rds"))
check(isTRUE(fixture$provenance$synthetic) && fixture$provenance$M == 3L &&
  identical(fixture$provenance$original_fits_sha256,
    "58485ec76a4cc7c5427a5f49df29552d22476c3b9cd5cbd21de8c297ccc54336"),
  "existing portable synthetic fixture provenance")

# Literal finite-sum oracle for two independent donor marks (M=3), and its
# derivative with the root argument fixed. It does not call the Laplace helper.
mark_oracle <- function(w, probability, derivative, root_weight, M) {
  if (M == 1L) return(c(A = 1, A_s = 0))
  if (M == 2L) return(c(A = sum(root_weight * probability / (root_weight + w)),
    A_s = sum(root_weight * derivative / (root_weight + w))))
  stopifnot(M == 3L)
  A <- A_s <- 0
  for (j in seq_along(w)) for (k in seq_along(w)) {
    fraction <- root_weight / (root_weight + w[j] + w[k])
    A <- A + fraction * probability[j] * probability[k]
    A_s <- A_s + fraction * (derivative[j] * probability[k] +
      probability[j] * derivative[k])
  }
  c(A = A, A_s = A_s)
}
root_weights <- c(0.7, 1.1)
mark_weights <- c(0.4, 1.3, 2.2)
probability <- c(0.2, 0.3, 0.5)
derivative <- c(0.11, -0.07, -0.04)
P <- matrix(rep(probability, 2), 3L)
DP <- matrix(rep(derivative, 2), 3L)
for (M in 1:3) {
  exact <- vapply(root_weights, function(w)
    mark_oracle(mark_weights, probability, derivative, w, M), numeric(2))
  value <- .wm_scalar_psm_laplace_values(mark_weights, P, DP, root_weights,
    M, cutoff = 40, panels = 64L, maximum_entries = 1000L)
  near(value$A, exact[1L, ], paste("M", M, "independent finite mark A oracle"))
  near(value$A_s, exact[2L, ], paste("M", M, "independent fixed-root partial derivative oracle"))
  if (M == 1L) check(identical(value$A, c(1, 1)) &&
    identical(value$A_s, c(0, 0)), "M1 exact constants")
}
constant <- .wm_scalar_psm_laplace_values(rep(1.25, 3L), P, DP,
  root_weights, 3L, cutoff = 1, panels = 1L, maximum_entries = 100L)
near(constant$A, root_weights / (root_weights + 2.5), "constant donor weight exact A")
check(identical(constant$A_s, c(0, 0)), "constant donor weight exact zero derivative")
check(all(abs(vapply(root_weights, function(w)
  mark_oracle(mark_weights, probability, derivative, w, 3L)[2L], numeric(1))) > 1e-4),
  "nonconstant fixed-root derivative check is informative")

results <- list()
oracle_results <- list()
for (target in c("PATE", "PATT")) {
  x <- fixture$fits[[target]]
  unchanged <- serialize(x, NULL)
  result <- wm_fitted_inference(x, covariance_scope = "scalar_psm")
  check(result$available && result$status == "conditional_scalar_psm",
    paste(target, "M3 nonconstant weights available under default controls"))
  check(identical(result$fit, x$fit) && identical(result$source_object, x) &&
    identical(serialize(x, NULL), unchanged), paste(target, "original point and object preserved"))
  check(!result$population_assumptions_verified && !result$assumptions_verified &&
    !result$numerical_certificate$asymptotic_root_rate_verified &&
    result$observed_logistic_check$stored_nearest_donors_validated,
    paste(target, "population and numerical-root status remain distinct"))
  n <- x$fit$n
  Y <- x$fit$data$Y
  Z <- x$fit$data$Z
  w <- x$fit$analysis_weights
  e <- x$propensity$probability
  X <- x$fitting_design
  p <- ncol(X)
  gamma <- mean(if (target == "PATE") w else Z * w)
  h <- result$controls$bandwidth
  active <- result$active
  d <- e * (1 - e) * X
  b_arm <- list()
  oracle_arm <- list()
  for (z in if (target == "PATE") 0:1 else 0L) {
    donor <- which(Z == z)
    component <- result$completed_components$drift[[paste0("arm", z)]]
    A_exact <- As_exact <- rep(NA_real_, n)
    for (i in component$donor_rows) {
      u <- (e[i] - e[donor]) / h
      kernel <- ifelse(abs(u) < 1, 35 / 32 * (1 - u^2)^3 / h, 0)
      prime <- ifelse(abs(u) < 1, -105 / 16 * u * (1 - u^2)^2 / h^2, 0)
      prob <- kernel / sum(kernel)
      dp <- (prime - prob * sum(prime)) / sum(kernel)
      literal <- mark_oracle(w[donor], prob, dp, w[i], x$fit$M)
      A_exact[i] <- literal[1L]
      As_exact[i] <- literal[2L]
      at <- match(i, component$donor_rows)
      envelope <- component$per_row_error_envelope
      check(abs(component$A[i] - A_exact[i]) <= envelope$A[at] +
        envelope$tail_A[at] + 1e-13, paste(target, z, i, "A real-arithmetic envelope"))
      check(abs(component$A_s[i] - As_exact[i]) <= envelope$A_s[at] +
        envelope$tail_A_s[at] + 1e-13, paste(target, z, i, "A_s signed-mixture envelope"))
    }
    sm <- result$completed_components$smoothing
    m <- sm$mean[, z + 1L]
    mp <- sm$mean_derivative[, z + 1L]
    donor_sum <- query_sum <- numeric(p)
    for (i in seq_len(n)) if (active[i]) {
      if (Z[i] != z) query_sum <- query_sum + w[i] * mp[i] * d[i, ]
      else donor_sum <- donor_sum + x$fit$M * d[i, ] *
        (A_exact[i] * (Y[i] - m[i]) * component$J_s[i] +
          component$J[i] * As_exact[i] * (Y[i] - m[i]) -
          component$J[i] * A_exact[i] * mp[i])
    }
    b_arm[[paste0("arm", z)]] <- (query_sum + donor_sum) / n
    oracle_arm[[paste0("arm", z)]] <- list(A = A_exact, A_s = As_exact)
  }
  b_oracle <- if (target == "PATE") (b_arm$arm1 - b_arm$arm0) / gamma else
    -b_arm$arm0 / gamma
  check(all(abs(result$total_sensitivity - b_oracle) <= result$drift_error_bound + 1e-12),
    paste(target, "complete signed drift within propagated oracle bound"))
  check(all(result$drift_error_bound <= result$controls$drift_error_target) &&
    result$laplace_terms > 0 && result$laplace_terms <= result$controls$maximum_laplace_terms,
    paste(target, "available finite-M computation meets declared error and work caps"))
  K <- matrix(0, n, 2L)
  query <- if (target == "PATE") seq_len(n) else which(Z == 1L)
  for (i in query) {
    donors <- x$fit$graph$neighbors[[i]]
    shares <- w[donors] / sum(w[donors])
    K[donors, 2L - Z[i]] <- K[donors, 2L - Z[i]] + w[i] * shares
  }
  mu <- result$feasible_means
  row <- vapply(seq_len(n), function(i) {
    if (target == "PATE") (w[i] * (mu[i, 2L] - mu[i, 1L] - x$estimate) +
      (2 * Z[i] - 1) * (w[i] + K[i, Z[i] + 1L]) * (Y[i] - mu[i, Z[i] + 1L])) / gamma
    else (Z[i] * w[i] * (Y[i] - mu[i, 1L] - x$estimate) -
      (1 - Z[i]) * K[i, 1L] * (Y[i] - mu[i, 1L])) / gamma
  }, numeric(1))
  near(row, result$uncentered_base_rows, paste(target, "independent complete raw-point row"))
  H <- Reduce("+", lapply(seq_len(n), function(i)
    w[i] * e[i] * (1 - e[i]) * tcrossprod(X[i, ]))) / n
  ell <- t(vapply(seq_len(n), function(i) solve(H, w[i] * (Z[i] - e[i]) * X[i, ]), numeric(p)))
  ell <- sweep(ell, 2L, colMeans(ell), "-")
  near(unname(ell), unname(result$nuisance_influence), paste(target, "independent sandwich influence"))
  base <- row - mean(row)
  augmented <- base + drop(ell %*% result$total_sensitivity)
  augmented <- augmented - mean(augmented)
  near(augmented, result$augmented_rows, paste(target, "complete plus-signed row augmentation"))
  near(result$root_n_variance, mean(base^2) +
    2 * sum(result$total_sensitivity * colMeans(ell * base)) +
    mean(drop(ell %*% result$total_sensitivity)^2), paste(target, "full covariance identity"))
  near(result$variance, result$root_n_variance / n, paste(target, "total-n sampling variance"))
  near(mean(result$conf.int), x$estimate, paste(target, "interval centered at original raw point"))
  if (target == "PATT") check(ncol(result$feasible_means) == 1L &&
    length(result$completed_components$drift) == 1L &&
    all(is.finite(row[Z == 1L])) && any(abs(row[Z == 1L]) > 1e-8),
    "PATT control mean only and complete treated-outcome rows")
  scaled <- x
  scaled$fit$weights <- 7 * scaled$fit$weights
  scaled$fit$weight_scale <- max(scaled$fit$weights)
  scaled$fit$analysis_weights <- scaled$fit$weights / scaled$fit$weight_scale
  sf <- wm_fitted_inference(scaled, covariance_scope = "scalar_psm")
  check(sf$available, paste(target, "scaled inference available"))
  near(sf$total_sensitivity, result$total_sensitivity, paste(target, "common weight scale drift invariance"))
  near(sf$augmented_rows, result$augmented_rows, paste(target, "common weight scale row invariance"))
  near(sf$root_n_variance, result$root_n_variance, paste(target, "common weight scale variance invariance"))
  chunk <- wm_fitted_inference(x, covariance_scope = "scalar_psm",
    scalar_control = list(maximum_matrix_entries = 100L))
  check(chunk$available, paste(target, "small bounded matrix chunks available"))
  near(chunk$augmented_rows, result$augmented_rows, paste(target, "chunk size preserves formula"))
  for (control in list(list(density_constant = 100), list(maximum_laplace_terms = 1),
      list(max_nodes = 8L, quadrature_error_constant = 1e-30))) {
    failed <- wm_fitted_inference(x, covariance_scope = "scalar_psm", scalar_control = control)
    check(!failed$available && identical(failed$estimate, x$estimate) &&
      identical(failed$fit, x$fit) && is.na(failed$variance) && !is.null(failed$reason),
      paste(target, paste(names(control), collapse = "/"), "runtime failure preserves raw point"))
  }
  mismatch <- X
  mismatch[1L, 2L] <- mismatch[1L, 2L] + 0.1
  fails(wm_fitted_inference(x, covariance_scope = "scalar_psm",
    scalar_control = list(design = mismatch)), paste(target, "changed design rejected"))
  changed <- x
  changed$fit$weights[1L] <- 2 * changed$fit$weights[1L]
  fails(wm_fitted_inference(changed, covariance_scope = "scalar_psm"),
    paste(target, "incoherent changed weights rejected"))
  changed <- x
  changed$fit$loads$incoming[1L, 1L] <- changed$fit$loads$incoming[1L, 1L] + 0.1
  fails(wm_fitted_inference(changed, covariance_scope = "scalar_psm"),
    paste(target, "changed graph loads rejected"))
  fails(wm_fitted_inference(x, covariance_scope = "scalar_psm", nuisance_influence = X),
    paste(target, "separately supplied influence rejected"))
  fails(wm_fitted_inference(x$fit, nuisance_influence = X,
    covariance_scope = if (target == "PATE") "full_x" else "patt"),
    paste(target, "old d>=2 scalar guard preserved"))
  if (length(args) >= 3L) {
    certificates <- readRDS(args[3L])$certificates
    certificate <- certificates[[target]]
    certified <- wm_fitted_inference(x, covariance_scope = "scalar_psm",
      scalar_control = list(root_certificate = certificate))
    check(certified$available &&
      certified$numerical_certificate$numerical_root_status ==
        "provided_finite_data_root_containment" &&
      certified$numerical_certificate$donor_set_status ==
        "provided_finite_data_donor_certificate" &&
      !certified$population_assumptions_verified &&
      !certified$numerical_certificate$asymptotic_root_rate_verified,
      paste(target, "existing exact finite-data certificate remains separate"))
    wrong <- certificate
    wrong$inputs$design[1L, 2L] <- wrong$inputs$design[1L, 2L] + 0.1
    fails(wm_fitted_inference(x, covariance_scope = "scalar_psm",
      scalar_control = list(root_certificate = wrong)),
      paste(target, "certificate for different bound design rejected"))
    inconclusive <- certificate
    inconclusive$root_contained <- inconclusive$donor_sets_certified <- FALSE
    partial <- wm_fitted_inference(x, covariance_scope = "scalar_psm",
      scalar_control = list(root_certificate = inconclusive))
    check(partial$available &&
      partial$numerical_certificate$numerical_root_status == "not_certified" &&
      partial$numerical_certificate$donor_set_status == "not_certified",
      paste(target, "inconclusive finite certificate does not change point or conditional formula"))
  }
  results[[target]] <- result
  oracle_results[[target]] <- list(b = b_oracle, arms = oracle_arm)
}

# A finite engineering reparameterization of stored logits tests p=2 algebra.
# It is not presented as a newly fitted or historically bound statistical model.
x2 <- fixture$fits$PATE
x2$fitting_design <- cbind(intercept = 1, saved_logit = stats::qlogis(x2$propensity$probability))
x2$design_columns <- colnames(x2$fitting_design)
r2 <- wm_fitted_inference(x2, covariance_scope = "scalar_psm")
check(r2$available && ncol(r2$nuisance_influence) == 2L &&
  length(r2$total_sensitivity) == 2L, "p2 stored-logit finite-array inference algebra")

# Coordinated wrong donor graph can satisfy old edge/load/point reconstruction;
# the new nearest-donor validation must reject it. Preserve every original array.
bad <- fixture$fits$PATE
n <- bad$fit$n
Z <- bad$fit$data$Z
w <- bad$fit$analysis_weights
Y <- bad$fit$data$Y
e <- bad$propensity$probability
i <- 1L
donors <- which(Z != Z[i])
excluded <- setdiff(donors, bad$fit$graph$neighbors[[i]])
j <- excluded[which.max(abs(e[excluded] - e[i]))]
replacement <- c(bad$fit$graph$neighbors[[i]][-bad$fit$M], j)
replacement <- replacement[order(abs(e[replacement] - e[i]), replacement)]
bad$fit$graph$neighbors[[i]] <- replacement
edges <- bad$fit$graph$edges
at <- which(edges$query == i)
edges$donor[at] <- replacement
edges$share[at] <- edges$outcome_share[at] <- w[replacement] / sum(w[replacement])
bad$fit$graph$edges <- edges
incoming <- matrix(0, n, 2L)
missing <- numeric(n)
for (r in seq_len(n)) {
  js <- bad$fit$graph$neighbors[[r]]
  share <- w[js] / sum(w[js])
  incoming[js, 2L - Z[r]] <- incoming[js, 2L - Z[r]] + w[r] * share
  missing[r] <- sum(share * Y[js])
}
bad$fit$loads$incoming <- incoming
bad$estimate <- bad$fit$estimate <- bad$fit$raw_estimate <-
  sum(w * (2 * Z - 1) * (Y - missing)) / sum(w)
check(!inherits(try(.wm_scalar_repl_graph(bad$fit, e), silent = TRUE), "try-error"),
  "coordinated wrong graph passes old finite reconstruction")
fails(wm_fitted_inference(bad, covariance_scope = "scalar_psm"),
  "coordinated wrong graph rejected by nearest donor binding")

check(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE) == had_rng &&
  (!had_rng || identical(rng_before, get(".Random.seed", envir = .GlobalEnv))), "RNG unchanged")
check(length(tools::parse_Rd(file.path(root, "man", "wm_fitted_inference.Rd"))) > 0,
  "current fitted inference documentation parses")
receipt <- list(status = "PASS", checks = length(labels), labels = labels,
  elapsed_seconds = proc.time()[["elapsed"]] - start,
  propensity_refits = 0, outcome_model_refits = 0, new_datasets = 0, RNG_used = FALSE,
  variance_nuisance_smoothers_computed = TRUE,
  saved_fixture_rows = 24, saved_M = 3,
  scope = "Independent finite-array formulas, numerical bounds, availability and guards; not sampling calibration",
  max_drift_error_bound = vapply(results, function(x) max(x$drift_error_bound), numeric(1)),
  default_laplace_terms = vapply(results, `[[`, numeric(1), "laplace_terms"))
if (length(args) >= 2L) saveRDS(list(receipt = receipt, results = results,
  oracle_results = oracle_results), args[2L])
cat(length(labels), "checks PASS; zero propensity/outcome model refits, RNG or new datasets.\n")
