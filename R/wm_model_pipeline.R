# Finite-model assembly. Exact arithmetic body of the existing WDSM
# pipeline is retained; only stack class and scope descriptions are generalized.
# Uses unchanged wm_match, graph-gradient and fitted-variance helpers.
# This conditional full-X branch never certifies model correctness/geometry.
.wm_model_fitted_pipeline <- function(stack, M = 3L,
                                    conf.level = 0.95,
                                    tie_rule = c("row_order", "source_random"),
                                    tie_seed = NULL, tie_tolerance = 64 * .Machine$double.eps) {
  ties <- .wm_tie_options(tie_rule, tie_seed, tie_tolerance)
  if (ties$rule == "source_random") {
    stop("source_random is point-only; full-X sampling assembly is unsupported.", call. = FALSE)
  }
  if (!inherits(stack, ".wm_model_nuisance_stack") ||
      !is.list(stack)) stop("Supply a finite candidate-model nuisance stack.")
  n <- stack$n
  if (length(stack$multiplicity) != n || anyNA(stack$multiplicity) ||
      any(stack$multiplicity != 1)) {
    stop("Sampling assembly requires original unit multiplicities.")
  }
  if (!is.matrix(stack$weight_derivative) ||
      any(!is.finite(stack$weight_derivative)) || any(stack$weight_derivative != 0)) {
    stop("This branch requires known weights with zero parameter derivative.")
  }
  d <- stack$inputs
  fit <- wm_match(Y = d$Y, Z = d$Z, weights = d$weights,
    scores0 = stack$scores0, scores1 = stack$scores1,
    M = M, estimand = stack$estimand, method = "self_normalized",
    mean0 = stack$mean0, mean1 = stack$mean1, variance = FALSE,
    tie_rule = ties$rule, tie_seed = ties$seed, tie_tolerance = ties$tolerance)
  gradient <- .wm_wdsm_graph_gradient(fit,
    weight_derivative = stack$weight_derivative,
    mean0_derivative = stack$mean0_derivative,
    mean1_derivative = stack$mean1_derivative)
  if (any(gradient$weight_sensitivity != 0)) {
    stop("Known-weight assembly produced a nonzero weight sensitivity.")
  }
  zero_graph <- stats::setNames(numeric(length(gradient$sensitivity)),
                          gradient$parameter_names)
  inference <- .wm_wdsm_fitted_variance(fit,
    nuisance_influence = stack$nuisance_influence,
    smooth_sensitivity = gradient$sensitivity,
    graph_sensitivity = zero_graph,
    covariance_scope = if (stack$estimand == "PATE") "full_x" else "patt",
    conf.level = conf.level)

  # Independently collect coefficients in the original fixed-K replicate
  # numerator. K and w here have the SAME fixed common normalization.
  Z <- d$Z
  w <- fit$analysis_weights
  K <- fit$loads$incoming[cbind(seq_len(n), Z + 1L)]
  q0 <- stack$mean0_derivative
  a <- colSums(q0 * ((1 - Z) * K - Z * w))
  if (stack$estimand == "PATE") {
    a <- a + colSums(stack$mean1_derivative * ((1 - Z) * w - Z * K))
  }
  a <- a / fit$denominator
  names(a) <- gradient$parameter_names
  slope_error <- max(abs(a - inference$total_sensitivity))
  if (!is.finite(slope_error) ||
      slope_error > 1e-10 * max(1, abs(a), abs(inference$total_sensitivity))) {
    stop("Original fixed-reuse slope disagrees with the query prediction gradient.")
  }
  inference$method <- "full_x"
  inference$sampling_inference_available <- inference$available
  inference$status <- if (inference$available) "conditional_full_x" else
    "nonpositive_variance"
  inference$unavailable_reason <- if (inference$available) "" else
    "The assembled variance is not strictly positive."
  inference$assumptions_verified <- FALSE
  inference$original_refit_limit_agreement_declared <- FALSE
  inference$contract <- paste("Known supplied weights; finite declared PS/arm-PG model lists and pooled",
    "scaling/quadratic-correction pipeline. This conditional full-X method sets",
    "graph transport to zero under full-X-correct correction means and the",
    "stated iid, identification, branch-specific support/overlap, geometry and joint regular-root",
    "conditions. Those premises are not verified by this calculation.",
    "Bounded alternatives retain their premises. The complete Gaussian alternative",
    "permits unbounded supplied weights and correct full-X predictions under exact",
    "design/noise-projection graph determination, integrated current-graph moment laws,",
    "prediction derivative/Hessian envelopes, actual slope and complete influence conditions.",
    "Source geometry and numerical statistical-root/graph equivalence remain application premises.",
    "All row/nuisance covariance is retained. No bootstrap claim is made.")
  donors0 <- lapply(seq_len(n), function(i)
    if (Z[i] == 1L) fit$graph$neighbors[[i]] else integer())
  donors1 <- if (stack$estimand == "PATE") lapply(seq_len(n), function(i)
    if (Z[i] == 0L) fit$graph$neighbors[[i]] else integer()) else NULL
  structure(list(estimate = fit$estimate, estimand = stack$estimand, n = n,
    M = fit$M, stack = stack, fit = fit, gradient = gradient,
    inference = inference, original_refit_slope = a,
    original_reuse_normalized = K, original_weight_scale = fit$weight_scale,
    original_donors = list(matches_0 = donors0, matches_1 = donors1),
    slope_identity_error = slope_error,
    original_refit_limit_agreement_declared = FALSE,
    assumptions_verified = FALSE, contract = inference$contract,
    matching_rule = "Exact Euclidean distance on supplied pooled-standardized scores; row-index ties",
    scope = paste("Original self-normalized PATE/PATT with the supplied finite arm maps,",
      "known weights and full-X-correct population correction means.",
      "Finite slope identity does not itself prove sampling or bootstrap validity.")),
    class = c(".wm_model_fitted_pipeline", "list"))
}
