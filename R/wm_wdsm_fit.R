# Public double-score wrapper. Does not change legacy interfaces.

.wm_wdsm_controls <- function(controls) {
  defaults <- list(glm_maxit = 100L, glm_epsilon = 1e-12,
    refine_maxit = 200L, refine_epsilon = 1e-14, rank_tolerance = 1e-10)
  if (!is.list(controls) || (length(controls) &&
      (is.null(names(controls)) || anyNA(names(controls)) ||
       anyDuplicated(names(controls)) ||
       any(!names(controls) %in% names(defaults))))) {
    stop("controls must be a named list of supported solver controls.", call. = FALSE)
  }
  for (name in names(controls)) defaults[[name]] <- controls[[name]]
  number <- function(x) is.numeric(x) && !is.complex(x) &&
    is.null(dim(x)) && length(x) == 1L && is.finite(x)
  for (name in c("glm_maxit", "refine_maxit")) {
    x <- defaults[[name]]
    if (!number(x) || x < 1 || x != floor(x) || x > .Machine$integer.max) {
      stop(name, " must be a positive representable integer.", call. = FALSE)
    }
  }
  for (name in c("glm_epsilon", "refine_epsilon", "rank_tolerance")) {
    x <- defaults[[name]]
    if (!number(x) || x <= 0 || x >= 1) {
      stop(name, " must be strictly between zero and one.", call. = FALSE)
    }
  }
  if (defaults$refine_maxit < defaults$glm_maxit ||
      defaults$refine_epsilon > defaults$glm_epsilon) {
    stop("Refinement must have at least as many iterations and no looser tolerance.",
         call. = FALSE)
  }
  defaults
}

# Gate arithmetic is adapted from .wm_cal_attempt without statistical changes.
# Arguments replace simulation model labels; warning/error diagnostics survive.
.wm_wdsm_attempt <- function(arguments, fitter, solver) {
  began <- proc.time()[[3L]]
  warnings <- character()
  value <- tryCatch(withCallingHandlers(do.call(fitter, c(arguments, solver)),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = identity)
  elapsed <- proc.time()[[3L]] - began
  if (inherits(value, "error")) {
    message <- conditionMessage(value)
    accuracy <- identical(message, "Weighted propensity root fails the normalized score check.")
    return(list(stack = NULL, record = list(status = "failed", error = message,
      accuracy_only = accuracy, controls = solver, elapsed_seconds = elapsed,
      warnings = warnings)))
  }
  step <- tryCatch(as.vector(solve(value$jacobian, value$mean_equation)), error = identity)
  inverse_ok <- !inherits(step, "error") && all(is.finite(step)) &&
    is.finite(value$diagnostics$jacobian_reciprocal_condition) &&
    value$diagnostics$jacobian_reciprocal_condition > 0
  residual <- if (inverse_ok) max(abs(step)/(1 + abs(value$parameter))) else Inf
  threshold <- min(1e-9, 1/length(arguments$Y))
  score <- value$diagnostics$normalized_ps_score
  okay <- inverse_ok && is.finite(score) && score <= 1e-10 && residual <= threshold
  record <- list(status = if (okay) "completed" else "failed",
    error = if (okay) "" else if (!inverse_ok) "Nonfinite inverse-Jacobian diagnostic" else
      "Root accuracy criterion failed", accuracy_only = !okay && inverse_ok && is.finite(score),
    controls = solver, elapsed_seconds = elapsed, warnings = warnings,
    parameter = value$parameter, mean_equation = value$mean_equation,
    newton_step = if (inverse_ok) step else NULL,
    newton_residual = residual, newton_threshold = threshold,
    diagnostics = value$diagnostics)
  list(stack = if (okay) value else NULL, record = record)
}

.wm_wdsm_fit_stack <- function(arguments, fitter, controls) {
  solvers <- list(
    list(glm_maxit = controls$glm_maxit, glm_epsilon = controls$glm_epsilon,
         rank_tolerance = controls$rank_tolerance),
    list(glm_maxit = controls$refine_maxit, glm_epsilon = controls$refine_epsilon,
         rank_tolerance = controls$rank_tolerance))
  attempts <- list()
  for (i in seq_along(solvers)) {
    result <- .wm_wdsm_attempt(arguments, fitter, solvers[[i]])
    attempts[[i]] <- result$record
    if (!is.null(result$stack)) {
      return(list(stack = result$stack, attempts = attempts, refined = i == 2L))
    }
    if (!isTRUE(result$record$accuracy_only)) break
  }
  condition <- structure(list(
    message = paste0("WDSM nuisance fitting failed: ", attempts[[length(attempts)]]$error),
    call = NULL, attempts = attempts, refined = length(attempts) == 2L),
    class = c("wm_wdsm_fit_error", "error", "condition"))
  stop(condition)
}

#' Fit weighted double score matching with optional full-X inference
#'
#' Matching tie options are forwarded to \code{wm_match()}.
#' \code{tie_rule = "source_random"} requires an explicit \code{tie_seed}
#' and \code{inference = "none"}; sampling inference for this policy is unsupported.
#'
#' @export
wm_wdsm_fit <- function(Y, Z, weights, ps_design, pg0_design,
                        pg1_design = NULL, M = 3L,
                        estimand = c("PATE", "PATT"),
                        inference = c("none", "full_x"), controls = list(),
                        conf.level = 0.95, ps_offset = NULL,
                        pg0_offset = NULL, pg1_offset = NULL,
                        ps_weighting = c("probability", "unit"),
                        tie_rule = c("row_order", "source_random"),
                        tie_seed = NULL, tie_tolerance = 64 * .Machine$double.eps) {
  call <- match.call()
  estimand <- match.arg(estimand)
  inference <- match.arg(inference)
  ps_weighting <- match.arg(ps_weighting)
  ties <- .wm_tie_options(tie_rule, tie_seed, tie_tolerance)
  if (ties$rule == "source_random" && inference != "none") {
    stop("source_random is point-only; use inference = 'none'. Sampling inference is unsupported.",
         call. = FALSE)
  }
  controls <- .wm_wdsm_controls(controls)
  if (!is.numeric(M) || is.complex(M) || !is.null(dim(M)) ||
      length(M) != 1L || !is.finite(M) || M < 1 || M != floor(M) ||
      M > .Machine$integer.max) {
    stop("M must be a positive representable integer.", call. = FALSE)
  }
  conf.level <- .wm_fv_number(conf.level, 1L, "conf.level")
  if (conf.level <= 0 || conf.level >= 1) {
    stop("conf.level must be strictly between zero and one.", call. = FALSE)
  }
  if (estimand == "PATT" && (!is.null(pg1_design) || !is.null(pg1_offset))) {
    stop("PATT does not fit an unused treated prognostic model.", call. = FALSE)
  }
  arguments <- list(Y = Y, Z = Z, weights = weights,
    ps_design = ps_design, pg0_design = pg0_design, pg1_design = pg1_design,
    estimand = estimand, ps_offset = ps_offset,
    pg0_offset = pg0_offset, pg1_offset = pg1_offset,
    ps_weighting = ps_weighting)
  fitter <- if (inference == "none") .wm_wdsm_prediction_stack else
    .wm_wdsm_nuisance_stack
  solved <- .wm_wdsm_fit_stack(arguments, fitter, controls)
  stack <- solved$stack
  if (inference == "full_x") {
    assembly <- .wm_wdsm_fitted_pipeline(stack, M = M, conf.level = conf.level,
      tie_rule = ties$rule, tie_seed = ties$seed, tie_tolerance = ties$tolerance)
    fit <- assembly$fit
    inference_result <- assembly$inference
  } else {
    fit <- wm_match(Y = stack$inputs$Y, Z = stack$inputs$Z,
      weights = stack$inputs$weights, scores0 = stack$scores0,
      scores1 = stack$scores1, M = M, estimand = estimand,
      method = "self_normalized", mean0 = stack$mean0,
      mean1 = stack$mean1, variance = FALSE,
      tie_rule = ties$rule, tie_seed = ties$seed, tie_tolerance = ties$tolerance)
    assembly <- NULL
    inference_result <- list(method = "none", status = "point_estimate_only",
      available = FALSE, sampling_inference_available = FALSE,
      unavailable_reason = "Sampling inference was not requested.",
      root_n_variance = NA_real_, variance = NA_real_, se = NA_real_,
      conf.int = c(lower = NA_real_, upper = NA_real_), conf.level = conf.level,
      assumptions_verified = FALSE, original_refit_limit_agreement_declared = FALSE,
      contract = "Point estimate only; no sampling variance or interval is asserted.")
    if (ties$rule == "source_random") {
      inference_result$unavailable_reason <- "Sampling inference is unsupported for source_random."
      inference_result$contract <- fit$info$contract
    }
  }
  scope <- paste("Two-dimensional arm-specific fitted PS/PG matching with known",
    "strictly positive probability weights. Logistic PS with declared fitting weights, unweighted arm",
    "OLS prognostic models, pooled sample centering/scaling, weighted quadratic",
    "correction, and original self-normalized donor fractions.",
    "Model identification, full-X correction correctness and theorem premises are",
    "unverified. Individual weights alone do not specify cluster or stratum dependence.")
  structure(list(call = call, estimate = fit$estimate,
    raw_estimate = fit$raw_estimate, correction = fit$correction,
    n = fit$n, M = fit$M, estimand = estimand, ps_weighting = ps_weighting,
    se = inference_result$se, variance = inference_result$variance,
    root_n_variance = inference_result$root_n_variance,
    conf.int = inference_result$conf.int, inference = inference_result,
    fit = fit, nuisance = stack, assembly = assembly,
    solver = list(controls = controls, attempts = solved$attempts,
                  refined = solved$refined,
                  normalized_ps_threshold = 1e-10,
                  newton_threshold = min(1e-9, 1/fit$n)),
    assumptions_verified = FALSE, bootstrap_requested = FALSE,
    scope = scope), class = c("wm_wdsm_fit", "list"))
}
