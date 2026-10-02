# Qualified source-form WDSM integration. Source the nuisance-stack and
# fitted-component reference files first. This is not a package export.

.wm_wp_contract <- function(contract) {
  if (!is.list(contract) || !is.character(contract$mode) ||
      length(contract$mode) != 1L || is.na(contract$mode) ||
      !contract$mode %in% c("algebra_only", "declared_full_x") ||
      !is.character(contract$justification) ||
      length(contract$justification) != 1L || is.na(contract$justification) ||
      !nzchar(trimws(contract$justification))) {
    stop("Supply a contract with mode and a nonempty justification.")
  }
  required <- c("iid_bounded_sampling", "known_probability_weights",
    "target_identification_overlap", "correct_full_x_means",
    "current_score_geometry", "regular_joint_nuisance_influence",
    "positive_limit_variance")
  if (contract$mode == "declared_full_x") {
    a <- contract$assumptions
    if (!is.logical(a) || !is.null(dim(a)) || anyNA(a) ||
        is.null(names(a)) || anyDuplicated(names(a)) ||
        !setequal(names(a), required) || !all(a)) {
      stop("The declared full-X contract requires all named theorem premises.")
    }
    if (!is.logical(contract$conditional_refit_expansion) ||
        length(contract$conditional_refit_expansion) != 1L ||
        is.na(contract$conditional_refit_expansion)) {
      stop("Declare conditional_refit_expansion as TRUE or FALSE separately.")
    }
  } else if (!is.null(contract$assumptions) ||
             !is.null(contract$conditional_refit_expansion)) {
    stop("The algebra-only contract does not declare sampling/refit premises.")
  }
  contract
}

# stack is the unit-multiplicity output of wm_wdsm_nuisance_stack_reference.
# No model label or goodness-of-fit diagnostic certifies full-X correctness.
wm_wdsm_fitted_pipeline_reference <- function(stack, M = 3L, contract,
                                               conf.level = 0.95) {
  contract <- .wm_wp_contract(contract)
  if (!inherits(stack, "wm_wdsm_nuisance_stack_reference") ||
      !is.list(stack)) stop("Supply the reviewed nuisance-stack reference output.")
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
  fit <- wdsmatch::wm_match(Y = d$Y, Z = d$Z, weights = d$weights,
    scores0 = stack$scores0, scores1 = stack$scores1,
    M = M, estimand = stack$estimand, method = "self_normalized",
    mean0 = stack$mean0, mean1 = stack$mean1, variance = FALSE)
  gradient <- wm_graph_gradient_reference(fit,
    weight_derivative = stack$weight_derivative,
    mean0_derivative = stack$mean0_derivative,
    mean1_derivative = stack$mean1_derivative)
  if (any(gradient$weight_sensitivity != 0)) {
    stop("Known-weight assembly produced a nonzero weight sensitivity.")
  }
  zero_graph <- setNames(numeric(length(gradient$sensitivity)),
                          gradient$parameter_names)
  inference <- wm_fitted_variance_reference(fit,
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
  inference$reference_contract <- contract
  inference$sampling_contract_declared <- contract$mode == "declared_full_x"
  inference$sampling_inference_available <-
    inference$sampling_contract_declared && inference$available
  if (!inference$sampling_contract_declared) {
    inference$available <- FALSE
    inference$se <- NA_real_
    inference$conf.int[] <- NA_real_
  }
  inference$contract <- paste("Known-weight source-form PS/arm-PG/pooled-scaling/",
    "quadratic-BC pipeline. Graph drift is set to zero only under the supplied",
    "full-X branch; its truth is not tested. Complete row/nuisance covariance is",
    "retained. Algebra-only mode withholds sampling intervals and standard errors.")
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
    original_refit_limit_agreement_declared =
      contract$mode == "declared_full_x" && isTRUE(contract$conditional_refit_expansion),
    assumptions_verified = FALSE, contract = contract,
    matching_rule = "Exact Euclidean distance on supplied pooled-standardized scores; row-index ties",
    scope = paste("Original self-normalized PATE/PATT with two-dimensional arm maps,",
      "known weights and full-X-correct population correction means.",
      "Finite slope identity does not itself prove sampling or bootstrap validity.")),
    class = c("wm_wdsm_fitted_pipeline_reference", "list"))
}

# Supplied count columns are used without drawing random numbers or refitting.
wm_wdsm_fitted_count_reference <- function(object, counts) {
  if (!inherits(object, "wm_wdsm_fitted_pipeline_reference")) {
    stop("Supply a fitted WDSM pipeline reference object.")
  }
  result <- wm_fitted_count_reference(object$inference, counts)
  result$sampling_contract_declared <- object$inference$sampling_contract_declared
  result$original_refit_limit_agreement_declared <-
    object$original_refit_limit_agreement_declared
  result$assumptions_verified <- FALSE
  result$reference_contract <- object$contract
  result
}
