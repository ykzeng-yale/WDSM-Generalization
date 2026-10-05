# Public fitted-wrapper handoff. Statistical routines and source fits remain
# unchanged: contribution preparation binds the complete retained system, and
# original replication refits predictions with the same full-n count vector.

.wm_wrapper_prediction_refit <- function(object) {
  if (!is.null(object$correction_fit)) {
    stop("Automatic fixed_reuse refitting supports the default complete quadratic correction only; supply an explicit compatible refit callback for this correction.",
         call. = FALSE)
  }
  bound <- .wm_mi_bind(object)
  if (!is.null(bound$unavailable)) {
    stop("The retained fitted pipeline cannot be bound for automatic refitting: ",
         bound$unavailable, call. = FALSE)
  }
  stack <- object$nuisance
  d <- stack$inputs
  controls <- .wm_wdsm_controls(object$solver$controls)
  pate <- identical(object$estimand, "PATE")
  if (inherits(stack, ".wm_wdsm_nuisance_stack")) {
    arguments <- list(Y = d$Y, Z = d$Z, weights = d$weights,
      ps_design = d$ps_design, pg0_design = d$pg0_design,
      pg1_design = if (pate) d$pg1_design else NULL,
      estimand = object$estimand, ps_offset = d$ps_offset,
      pg0_offset = d$pg0_offset, pg1_offset = if (pate) d$pg1_offset else NULL,
      ps_weighting = stack$ps_weighting)
    fitter <- .wm_wdsm_prediction_stack
  } else if (inherits(stack, ".wm_model_nuisance_stack")) {
    maximum <- object$allocation$max_stack_elements
    if (is.null(maximum)) maximum <- 5e7
    arguments <- list(Y = d$Y, Z = d$Z, weights = d$weights,
      ps_models = d$ps_models, pg0_models = d$pg_models[["0"]],
      pg1_models = if (pate) d$pg_models[["1"]] else NULL,
      estimand = object$estimand, max_stack_elements = maximum)
    fitter <- .wm_model_prediction_stack
  } else {
    stop("The retained prediction pipeline is unsupported for automatic refitting.",
         call. = FALSE)
  }
  function(multiplicity) {
    current <- arguments
    current$multiplicity <- multiplicity
    solved <- .wm_wdsm_fit_stack(current, fitter, controls)$stack
    list(mean0 = solved$mean0, mean1 = if (pate) solved$mean1 else NULL)
  }
}

.wm_fitted_wrapper_bootstrap <- function(object, B, B_missing, seed, conf.level,
                                         interval, chunk_size, counts, method, refit) {
  .wm_mi_structure(object)
  # Keep B truly omitted when the caller omitted it: supplied count columns
  # determine their number through the unchanged downstream validation.
  arguments <- list(seed = seed, conf.level = conf.level, interval = interval,
    chunk_size = chunk_size, counts = counts, method = method)
  if (!B_missing) arguments$B <- B
  if (identical(method, "fixed_reuse")) {
    automatic <- is.null(refit)
    if (automatic) refit <- .wm_wrapper_prediction_refit(object)
    if (!is.function(refit)) {
      stop("refit must be NULL or a function of one multiplicity vector.", call. = FALSE)
    }
    arguments$object <- object$fit
    arguments$refit <- refit
    result <- do.call(wm_bootstrap, arguments)
    result$refit_mode <- if (automatic) "complete_retained_pipeline" else "supplied_callback"
    result$inference_contract <- paste(result$inference_contract,
      if (automatic) paste("The fitted-wrapper adapter refits each declared PS model,",
        "arm PG model, pooled center/variance and complete quadratic prediction",
        "using the same full-n counts and retained solver controls. The original",
        "supplied weights, matching graph and incoming loads are unchanged.",
        "This execution does not verify the sampling/refit-law or all-draw numerical premises.") else
        "The supplied callback's complete-system or compatible prediction-only premises remain caller requirements.")
  } else {
    if (!is.null(refit)) {
      stop("refit is supported only with method = 'fixed_reuse'.", call. = FALSE)
    }
    if (inherits(object$inference, "wm_fitted_inference")) {
      inference <- object$inference
      if (!identical(inference$fit, object$fit)) {
        stop("Prepared inference must retain the exact original fitted-wrapper source fit; a replaced or stale source is unsupported.",
             call. = FALSE)
      }
    } else {
      if (!identical(object$inference$method, "full_x")) {
        stop("Contribution replication of a fitted wrapper requires inference = 'full_x' or an explicit wm_model_inference handoff; point-only fits do not supply a sampling law.",
             call. = FALSE)
      }
      pate <- identical(object$estimand, "PATE")
      prepared <- wm_model_inference(object,
        covariance_scope = if (pate) "full_x" else "patt",
        transport0 = if (pate) NULL else list(mode = "zero", basis = "current_centering"),
        conf.level = conf.level)
      inference <- prepared$inference
    }
    arguments$object <- inference
    result <- do.call(wm_bootstrap, arguments)
    result$models_refitted <- FALSE
    result$original_refit_bootstrap <- FALSE
  }
  result$fitted_wrapper_class <- class(object)
  result$original_fit_preserved <- TRUE
  result$supplied_weights_known <- TRUE
  result$population_assumptions_verified <- FALSE
  result
}
