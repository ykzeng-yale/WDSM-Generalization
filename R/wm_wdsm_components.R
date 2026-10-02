# Private package-overlay copy. See SOURCE_MAPPING.md.
# Audited source: fitted_wm_components.R; SHA256 84f106202539be24d9211f296454ea099b230f503bb1f056887a2b6b2f7201bc.
# Reference implementation of the reviewed fitted-map variance components.
# Not an exported package API or an automatic verification of a fitted model.

.wm_fv_number <- function(x, expected_length, label) {
  if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) ||
      length(x) != expected_length || any(!is.finite(x))) {
    stop(label, " must be a finite numeric vector of length ", expected_length, ".",
         call. = FALSE)
  }
  x
}

.wm_fv_state <- function(fit) {
  if (!is.list(fit) || !identical(fit$method, "self_normalized")) {
    stop("Only original self_normalized matching is supported.", call. = FALSE)
  }
  # Reuse strict finite-graph reconstruction without computing an unused
  # fixed-map interval or reciprocal variance.
  checked <- .wm_reciprocal_validated_state(fit)
  list(n = checked$n, pATE = identical(checked$estimand, "PATE"),
       estimate = checked$estimate, gamma = checked$gamma,
       scale = checked$weight_scale, rows = checked$actual,
       Y = fit$data$Y, Z = fit$data$Z, w = fit$analysis_weights,
       mu0 = fit$predictions$mean0, mu1 = fit$predictions$mean1,
       edges = fit$graph$edges, estimand = checked$estimand,
       dimensions = c(ncol(fit$graph$scores0),
         if (identical(checked$estimand, "PATE")) ncol(fit$graph$scores1)),
       reciprocal_state = checked)
}

.wm_fv_matrices <- function(values, n) {
  supplied <- values[!vapply(values, is.null, logical(1L))]
  if (!length(supplied)) stop("Supply at least one derivative matrix.", call. = FALSE)
  p <- ncol(supplied[[1L]])
  if (is.null(p) || p < 1L) stop("Derivatives must be n-by-p matrices.", call. = FALSE)
  parameter_names <- NULL
  for (name in names(supplied)) {
    x <- supplied[[name]]
    if (!is.matrix(x) || !is.numeric(x) || is.complex(x) ||
        !identical(dim(x), c(as.integer(n), as.integer(p))) ||
        any(!is.finite(x))) {
      stop(name, " must be a finite numeric n-by-p matrix.", call. = FALSE)
    }
    nm <- colnames(x)
    if (!is.null(nm)) {
      if (anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm) ||
          (!is.null(parameter_names) && !identical(nm, parameter_names))) {
        stop("Derivative column names must be unique and use one parameter order.",
             call. = FALSE)
      }
      parameter_names <- nm
    }
  }
  if (is.null(parameter_names)) parameter_names <- paste0("parameter", seq_len(p))
  for (name in names(values)) {
    if (is.null(values[[name]])) values[[name]] <- matrix(0, n, p)
    colnames(values[[name]]) <- parameter_names
  }
  list(values = values, p = p, names = parameter_names)
}

# Exact derivative of the corrected estimator with its graph held fixed.
# Weight derivatives are of the original raw weights, not normalized weights.
.wm_wdsm_graph_gradient <- function(fit, weight_derivative = NULL,
                                        mean0_derivative = NULL,
                                        mean1_derivative = NULL) {
  s <- .wm_fv_state(fit)
  if (!s$pATE && !is.null(mean1_derivative)) {
    stop("PATT has no arm-1 prediction derivative.", call. = FALSE)
  }
  d <- .wm_fv_matrices(list(weight = weight_derivative,
                            mean0 = mean0_derivative,
                            mean1 = mean1_derivative), s$n)
  dw <- d$values$weight / s$scale
  if (any(!is.finite(dw))) stop("Weight derivative normalization overflow.")
  queries <- if (s$pATE) seq_len(s$n) else which(s$Z == 1L)
  denominator <- sum(s$w[queries])
  weight_part <- prediction_part <- numeric(d$p)
  reconstructed_numerator <- 0
  fraction_error <- 0
  for (i in queries) {
    take <- s$edges$query == i
    donors <- s$edges$donor[take]
    donor_sum <- sum(s$w[donors])
    fraction <- s$w[donors] / donor_sum
    dfraction <- (dw[donors, , drop = FALSE] -
      outer(fraction, colSums(dw[donors, , drop = FALSE]))) / donor_sum
    if (s$Z[i] == 1L) {
      mu <- s$mu0
      dm <- d$values$mean0
    } else {
      mu <- s$mu1
      dm <- d$values$mean1
    }
    residual <- s$Y[donors] - mu[donors]
    missing <- mu[i] + sum(fraction * residual)
    sign <- 1 - 2 * s$Z[i]
    contrast <- sign * (missing - s$Y[i])
    weight_part <- weight_part + dw[i, ] * (contrast - s$estimate) +
      s$w[i] * sign * as.vector(crossprod(residual, dfraction))
    prediction_part <- prediction_part + s$w[i] * sign *
      (dm[i, ] - colSums(dm[donors, , drop = FALSE] * fraction))
    reconstructed_numerator <- reconstructed_numerator + s$w[i] * contrast
    fraction_error <- max(fraction_error, abs(colSums(dfraction)))
  }
  weight_part <- weight_part / denominator
  prediction_part <- prediction_part / denominator
  total <- weight_part + prediction_part
  if (any(!is.finite(c(total, weight_part, prediction_part,
                       reconstructed_numerator, fraction_error)))) {
    stop("Gradient arithmetic exceeded numerical range.", call. = FALSE)
  }
  if (abs(reconstructed_numerator / denominator - s$estimate) >
      1e-10 * max(1, abs(s$estimate))) {
    stop("The independent query reconstruction disagrees with the fit.")
  }
  names(total) <- names(weight_part) <- names(prediction_part) <- d$names
  list(sensitivity = total, weight_sensitivity = weight_part,
       prediction_sensitivity = prediction_part,
       max_fraction_derivative_sum = fraction_error,
       parameter_names = d$names, estimate = s$estimate,
       weight_derivative_units = "Original raw weights; one constant normalization cancels",
       contract = paste("Exact fixed-edge smooth-mark derivative only.",
         "Prediction, donor-weight and target-denominator derivatives are included.",
         "Graph-transport sensitivity is not estimated by holding the graph fixed."))
}

# Assemble the new theorem's diagonal Gaussian block with a separately
# supplied full graph-transport sensitivity. Both slopes use estimator units.
.wm_wdsm_fitted_variance <- function(fit, nuisance_influence,
                                         smooth_sensitivity,
                                         graph_sensitivity,
                                         covariance_scope,
                                         conf.level = 0.95) {
  s <- .wm_fv_state(fit)
  if (!is.matrix(nuisance_influence) || !is.numeric(nuisance_influence) ||
      is.complex(nuisance_influence) || nrow(nuisance_influence) != s$n ||
      ncol(nuisance_influence) < 1L || any(!is.finite(nuisance_influence))) {
    stop("nuisance_influence must be a finite numeric n-by-p matrix.")
  }
  p <- ncol(nuisance_influence)
  smooth <- .wm_fv_number(smooth_sensitivity, p, "smooth_sensitivity")
  graph <- .wm_fv_number(graph_sensitivity, p, "graph_sensitivity")
  nm <- colnames(nuisance_influence)
  if (!is.null(nm) && (anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm))) {
    stop("Influence column names must be unique and nonempty.")
  }
  for (x in list(smooth, graph)) {
    if (!is.null(names(x)) &&
        (anyNA(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x)) ||
         (!is.null(nm) && !identical(names(x), nm)))) {
      stop("Sensitivity names must match influence parameter order.")
    }
    if (is.null(nm) && !is.null(names(x))) nm <- names(x)
  }
  if (is.null(nm)) nm <- paste0("parameter", seq_len(p))
  if (!is.character(covariance_scope) || length(covariance_scope) != 1L ||
      is.na(covariance_scope) ||
      !covariance_scope %in% c("distinct_rarity", "full_x", "common_field", "patt",
                              "regular_joint_d_gt2")) {
    stop("Declare one supported covariance_scope explicitly.")
  }
  reciprocal_scope <- identical(covariance_scope, "regular_joint_d_gt2")
  if (reciprocal_scope && (!s$pATE || length(s$dimensions) != 2L ||
      any(s$dimensions <= 2L) || s$dimensions[1L] != s$dimensions[2L])) {
    stop("regular_joint_d_gt2 requires PATE with equal used dimensions d > 2.",
         call. = FALSE)
  }
  scalar <- all(s$dimensions == 1L)
  if (!all(s$dimensions >= 2L) &&
      !(scalar && covariance_scope %in% c("full_x", "patt"))) {
    stop("Use all dimensions >= 2, or all scalar dimensions with full_x/PATT patt; mixed and other scalar scopes are unsupported.")
  }
  if ((s$pATE && covariance_scope == "patt") ||
      (!s$pATE && covariance_scope != "patt")) {
    stop("PATT uses covariance_scope = 'patt'; PATE needs a PATE branch.")
  }
  if ((covariance_scope == "full_x" || scalar) && any(graph != 0)) {
    stop("The full-X centering branch has zero graph transport.")
  }
  conf.level <- .wm_fv_number(conf.level, 1L, "conf.level")
  if (conf.level <= 0 || conf.level >= 1) stop("conf.level must be in (0, 1).")
  base <- s$rows - mean(s$rows)
  influence <- sweep(nuisance_influence, 2L, colMeans(nuisance_influence), "-")
  total_sensitivity <- smooth + graph
  names(total_sensitivity) <- nm
  nuisance_rows <- as.vector(influence %*% total_sensitivity)
  augmented <- base + nuisance_rows
  augmented <- augmented - mean(augmented)
  V0 <- mean(base^2)
  C <- as.vector(crossprod(influence, base)) / s$n
  Sigma <- crossprod(influence) / s$n
  cross <- 2 * sum(total_sensitivity * C)
  nuisance_variance <- mean(nuisance_rows^2)
  root_variance <- mean(augmented^2)
  if (any(!is.finite(c(base, influence, total_sensitivity, V0, C, Sigma,
                       cross, nuisance_variance, root_variance)))) {
    stop("Fitted-variance arithmetic exceeded numerical range.")
  }
  if (reciprocal_scope) {
    return(.wm_regular_joint_fitted_variance(s, base, influence, smooth, graph,
      total_sensitivity, nm, augmented, V0, C, Sigma, cross,
      nuisance_variance, root_variance, conf.level))
  }
  discrepancy <- abs(root_variance - (V0 + cross + nuisance_variance))
  if (discrepancy > 1e-10 * max(1, V0, abs(cross), nuisance_variance)) {
    stop("Augmented rows disagree with the covariance-block formula.")
  }
  variance <- root_variance / s$n
  se <- sqrt(variance)
  available <- root_variance > 0 && variance > 0
  ci <- if (available) s$estimate + c(-1, 1) *
    stats::qnorm((1 - conf.level) / 2, lower.tail = FALSE) * se else rep(NA_real_, 2L)
  if (available && any(!is.finite(ci))) stop("Confidence interval overflow.")
  names(ci) <- c("lower", "upper")
  names(C) <- nm
  dimnames(Sigma) <- list(nm, nm)
  structure(list(estimate = s$estimate, n = s$n, estimand = s$estimand,
    smooth_sensitivity = smooth, graph_sensitivity = graph,
    total_sensitivity = total_sensitivity, parameter_names = nm,
    V0 = V0, C = C, Sigma = Sigma, cross_term = cross,
    nuisance_variance = nuisance_variance, root_n_variance = root_variance,
    variance = variance, se = if (available) se else NA_real_,
    conf.int = ci, conf.level = conf.level, available = available,
    base_rows = base, augmented_rows = augmented,
    nuisance_influence = influence, covariance_scope = covariance_scope,
    covariance_formula_error = discrepancy,
    reciprocal_subtraction = 0, assumptions_verified = FALSE,
    contract = if (scalar) paste("Original corrected all-scalar estimator with known positive",
      "W=w_Z(X) and correct full-X predictions. The design-measurable",
      "monotone-chart bounded-mark law, or the exact outcome-fitted graph-root law with",
      "continuous block-score submersion, quantitative finite-piecewise marked charts",
      "and finite/conditional 2+delta moments, retains its bounded weight/prediction premises.",
      "Alternatively the applicable noncompact design-fitted or complete Gaussian law",
      "must supply integrated current-graph square/cross/maximum and prediction-gradient",
      "laws, prediction envelopes and complete regular influence. Gaussian graph determination",
      "is exact given all design rows and the complete noise projection. Complete joint influence,",
      "finite-M prediction slope, empirical L2 consistency and zero graph transport",
      "retain every covariance block, including singular full prediction covariance.",
      "Population, exact/numerical graph-root and callback conditions remain unverified.") else
      paste("Original corrected estimator; qualified iid theorem with fixed",
      "used dimensions >= 2, verified identification/centering/current geometry,",
      "joint regular influence and consistent full graph plus smooth sensitivity.",
      "Bounded alternatives retain their premises. Correct-full-X finite-moment/current-graph",
      "or complete Gaussian alternatives require their integrated moments, prediction",
      "envelopes, actual row/slope laws and exact Gaussian graph determination where applicable.",
      "The declared covariance branch is a caller premise, not inferred from data.",
      "Uses the proved diagonal block; does not justify unrestricted reciprocal",
      "geometry, general scalar maps, or the original nuisance-refit bootstrap.")),
    class = c(".wm_wdsm_fitted_variance", "list"))
}

# Complete contrast variance for the qualified equal-d>2 baseline-sheet branch.
# A nonpositive empirical oracle block never gates a positive full contrast.
.wm_regular_joint_fitted_variance <- function(s, base, influence, smooth, graph,
    total_sensitivity, nm, augmented, diagonal, C, Sigma, cross,
    nuisance_variance, augmented_variance, conf.level) {
  pairs <- .wm_reciprocal_pair_state(s$reciprocal_state)
  subtraction <- 2 * pairs$reciprocal / s$gamma^2
  V0 <- diagonal - subtraction
  raw_root <- augmented_variance - subtraction
  raw_variance <- raw_root / s$n
  if (any(!is.finite(c(subtraction, V0, raw_root, raw_variance)))) {
    stop("Reciprocal fitted-variance arithmetic exceeded numerical range.",
         call. = FALSE)
  }
  discrepancy <- abs(raw_root - (V0 + cross + nuisance_variance))
  if (discrepancy > 1e-10 * max(1, diagonal, abs(V0), abs(subtraction),
                                abs(cross), nuisance_variance)) {
    stop("Corrected variance disagrees with its complete covariance blocks.",
         call. = FALSE)
  }
  available <- raw_root > 0 && raw_variance > 0
  root_variance <- if (available) raw_root else NA_real_
  variance <- if (available) raw_variance else NA_real_
  se <- if (available) sqrt(variance) else NA_real_
  ci <- if (available) s$estimate + c(-1, 1) *
    stats::qnorm((1 - conf.level) / 2, lower.tail = FALSE) * se else rep(NA_real_, 2L)
  if (available && any(!is.finite(c(root_variance, variance, se, ci)))) {
    stop("Reciprocal fitted confidence interval exceeded numerical range.", call. = FALSE)
  }
  names(ci) <- c("lower", "upper")
  names(C) <- nm
  dimnames(Sigma) <- list(nm, nm)
  structure(list(estimate = s$estimate, n = s$n, estimand = s$estimand,
    smooth_sensitivity = smooth, graph_sensitivity = graph,
    total_sensitivity = total_sensitivity, parameter_names = nm,
    V0 = V0, C = C, Sigma = Sigma, cross_term = cross,
    nuisance_variance = nuisance_variance, root_n_variance = root_variance,
    variance = variance, se = se, conf.int = ci, conf.level = conf.level,
    available = available, base_rows = base, augmented_rows = augmented,
    nuisance_influence = influence, covariance_scope = "regular_joint_d_gt2",
    covariance_formula_error = discrepancy, reciprocal_subtraction = subtraction,
    reciprocal_numerator = pairs$reciprocal, reciprocal_pairs = pairs$pairs,
    diagonal_root_n_variance = diagonal,
    augmented_row_root_n_variance = augmented_variance,
    unregularized_root_n_variance = raw_root, unregularized_variance = raw_variance,
    analysis_target_weight_mean = s$gamma, weight_scale = s$scale,
    reciprocal_units = "Analysis-weight numerator; subtraction is root-n variance",
    variance_floor = NULL, floor_active = FALSE,
    diagnostic = if (available) "Positive complete reciprocal-corrected contrast variance; no floor."
      else "Nonpositive complete reciprocal-corrected contrast variance; inference unavailable; no floor.",
    assumptions_verified = FALSE,
    contract = paste("Original corrected known-weight PATE, equal fixed d>2, bounded iid",
      "own-field framework and regular baseline joint sheets. Full nuisance slope and",
      "signed reciprocal covariance retained. Population premises remain unverified;",
      "no rematching/refit bootstrap or actual sampling-variance convergence claim.")),
    class = c(".wm_wdsm_fitted_variance", "list"))
}

# Corrected one-step multinomial replication, not the original refit algorithm.
# Counts are supplied to make exact paired/algebra comparisons reproducible.
.wm_wdsm_fitted_count <- function(inference, counts) {
  if (!inherits(inference, ".wm_wdsm_fitted_variance") ||
      !is.list(inference)) stop("Supply a fitted-variance reference result.")
  if (identical(inference$covariance_scope, "regular_joint_d_gt2")) {
    stop("regular_joint_d_gt2 cannot use ordinary full-slope count replication.", call. = FALSE)
  }
  n <- inference$n
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 2 || n != floor(n)) stop("Invalid stored sample size.")
  rows <- .wm_fv_number(inference$augmented_rows, n, "augmented_rows")
  estimate <- .wm_fv_number(inference$estimate, 1L, "estimate")
  stored <- .wm_fv_number(inference$root_n_variance, 1L, "root_n_variance")
  rows <- rows - mean(rows)
  exact_root <- mean(rows^2)
  if (!is.finite(exact_root) || abs(stored - exact_root) >
      1e-10 * max(1, exact_root)) {
    stop("Stored variance disagrees with the augmented rows.")
  }
  if (!is.matrix(counts) || !is.numeric(counts) || is.complex(counts) ||
      nrow(counts) != n || ncol(counts) < 1L || any(!is.finite(counts)) ||
      any(counts < 0) || any(counts != floor(counts)) ||
      any(colSums(counts) != n)) {
    stop("counts must be n-by-B nonnegative integers with column sums n.")
  }
  roots <- as.vector(crossprod(rows, counts - 1)) / sqrt(n)
  draws <- estimate + roots / sqrt(n)
  mc_root <- mean((roots - mean(roots))^2)
  if (any(!is.finite(c(roots, draws, mc_root)))) stop("Replication overflow.")
  list(estimate = estimate, n = n, B = ncol(counts), draws = draws,
       root_n_draws = roots, conditional_root_n_variance = exact_root,
       conditional_variance = exact_root / n,
       monte_carlo_root_n_variance = mc_root,
       monte_carlo_variance = mc_root / n, variance_divisor = "B",
       method = "Full-slope one-step multinomial replication",
       contract = paste("Conditional variance assumes Multinomial(n, 1/n) counts.",
         "Supplied counts are checked algebraically; their random generation is not verified.",
         "All full-slope nuisance corrections are in the augmented rows.",
         "No graph or nuisance refit is performed; empty-arm count draws remain defined.",
         "Sampling validity requires the fitted theorem and empirical row Lindeberg conditions."))
}
