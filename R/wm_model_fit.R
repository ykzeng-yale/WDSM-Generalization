# Common finite-model wrapper; legacy public functions remain available.
# Model descriptors contain explicit design matrices (including their intended
# intercepts), offsets, and for PS only a declared probability/unit fit weight.

#' Weighted Matching with Finite Lists of Fitted Score Models
#'
#' Fits finite lists of logistic propensity and arm-specific linear prognostic
#' models and pooled coordinate standardization, with complete weighted quadratic
#' corrections by default or an explicit local-polynomial option, then applies the original self-normalized matching estimator.
#'
#' @param Y Finite numeric outcome vector.
#' @param Z Binary treatment vector containing both arms.
#' @param weights Finite, strictly positive supplied known probability weights,
#'   which may depend on treatment and covariates. Their interpretation must
#'   identify the intended estimand. Weight-estimation uncertainty is excluded.
#' @param ps_models List of descriptors with a finite numeric \code{design}
#'   matrix, optional fixed logistic-scale \code{offset}, and optional
#'   \code{weighting = "probability"} or \code{"unit"}. Include intended
#'   intercept columns. Each PS model is fitted once and shared across used arms.
#' @param pg0_models,pg1_models Lists of descriptors with \code{design} and optional
#'   fixed outcome-scale \code{offset}. Fits are ordinary arm-specific OLS.
#'   PATT requires \code{pg1_models = NULL}. Empty lists are allowed if every
#'   used arm retains at least one PS or PG coordinate.
#' @param M Positive integer number of distinct opposite-arm donors per query;
#'   donors may be reused across queries.
#' @param estimand Population average effect, PATE, or treated-population effect,
#'   PATT.
#' @param inference \code{"none"} returns point estimation only;
#'   \code{"full_x"} requests conditional full-X fitted-score inference.
#' @param controls Root controls supported by \code{wm_wdsm_fit()}.
#' @param conf.level Confidence level strictly between zero and one.
#' @param max_stack_elements Conservative element cap for complete nuisance
#'   arrays, including the legacy dispatch. Inputs, copies, solver overhead and
#'   total resident memory are outside this cap. No model/basis term is trimmed.
#' @param correction NULL retains the exact default quadratic/legacy path.
#'   An opt-in list with method='local_polynomial', explicit degree and
#'   bandwidth_exponent (one value or control/treated values), and optional
#'   positive guard_exponent (default 2) uses the constructive correction.
#'   This option requires inference='none', but retains the complete score-root
#'   IF and caller-completed reference-inference inputs without a quadratic fit.
#'   Alternatively, list(method='quadratic_qr') solves the same complete
#'   quadratic WLS in invertible QR correction coordinates. It requires
#'   inference='none', retains the score-only root and all correction columns,
#'   and leaves the matching scores unchanged.
#' @param tie_rule,tie_seed,tie_tolerance Matching options from \code{wm_match()}.
#'   \code{"source_random"} requires an explicit seed and point-only inference.
#'
#' @details Matching dimensions are \eqn{d_z=J+K_z}, distinct from the full
#'   nuisance dimension. Model/column names must be nonempty and unique when
#'   supplied. Models or dependent output coordinates are not silently merged
#'   or deleted. PS coordinates are probabilities. Pooled centers and variances
#'   use unweighted empirical moments with divisor n. Corrections include the
#'   complete intercept/linear/quadratic basis, fitted with probability weights.
#'
#'   The donor fraction remains supplied weight divided by that query's donor
#'   weight total. PATE uses the total weight denominator and PATT the treated
#'   total. PS-only configurations remain quadratically corrected and are not
#'   automatically the raw propensity-score estimator.
#'
#'   For one PS and one PG in every used arm, the wrapper delegates to
#'   \code{wm_wdsm_fit()} with identical statistical arguments. Underlying
#'   arithmetic is preserved; call/class and recipe/allocation metadata differ.
#'
#'   The full Jacobian differentiates both generated bases and predictions.
#'   Full-X inference retains all nuisance/covariance blocks and sets graph
#'   transport to zero conditionally on correct full-X correction means. Used
#'   dimensions must all be at least two or all equal one, with the separate
#'   scalar conditions retained. Mixed scalar/vector inference is unsupported.
#'   Point estimation permits mixed maps. Independent observations, identification,
#'   regular roots, correction correctness, support/geometry and row/slope laws
#'   remain unverified premises. See \code{wm_fitted_inference()} for the
#'   separate bounded, finite-moment and complete Gaussian alternatives.
#'   A correct PS candidate does not automatically establish these conditions.
#'
#'   Default-correction point-only mode skips IF/covariance arrays but retains root checks:
#'   normalized PS residual at most 1e-10 and maximum scaled Newton step at
#'   most min(1e-9,1/n). Accuracy-only refinement and structured
#'   \code{wm_wdsm_fit_error} conditions follow \code{wm_wdsm_fit()}.
#'   No bootstrap is run; fixed-map replication of the nested fit does not
#'   account for the complete fitted pipeline. Weights and offsets are known.
#'   Their estimation uncertainty and cluster/stratum inference are excluded.
#'
#'   The local-polynomial option uses a fixed normalized radial C2 triweight
#'   kernel on the unit ball, all monomials through the requested degree,
#'   bandwidth n^(-alpha), raw-W total-n Gram normalization and guard n^(-g).
#'   Its derivative differentiates the fitted intercept, including the kernel
#'   and shifted basis. No ridge, basis deletion or bandwidth retry is used.
#'   Guard failures retain the stated zero-level completion and NA derivatives
#'   with explicit diagnostics. Successful evaluated rows do not verify the
#'   required uniform neighborhood event. The bounded-Y/bounded-W correction
#'   theorem does not inherit unbounded Gaussian full-X alternatives.
#'   Local coefficients are auxiliary and do not enter the finite score root;
#'   B'gradient estimates the reference derivative, not an LP-refit derivative.
#'
#'   Quadratic QR retains a complete degree-two polynomial space under an
#'   invertible coordinate transform; it does not fit a different regression
#'   or use QR distances for matching. Unresolved score rank or unidentified
#'   donor regression stops without a ridge, deletion or retry. Numerical
#'   diagnostics do not prove a root-n prediction-error rate. Without explicit
#'   quadratic_control, wm_model_inference returns unavailable regular
#'   full-coefficient inference. Its separate nested_contrast preparation can
#'   be passed to wm_bootstrap under the scoped bounded outcome-correct nested
#'   model and numerical premises. The original point and graph are retained;
#'   covariance-only local graphs vary with the nuisance index. That route
#'   supports basic or no intervals, and any failed requested slot makes all
#'   aggregate variance, standard-error and interval fields unavailable.
#'
#' @return A \code{wm_model_fit} list with point/interval fields and nested
#'   \code{fit}, \code{nuisance}, \code{inference}, \code{assembly},
#'   \code{solver}, \code{model_recipe} and \code{allocation} records.
#'   Point-only interval fields are NA and assembly is NULL. Legacy
#'   specializations also inherit \code{wm_wdsm_fit}. Interval availability
#'   requires finite positive variance. Assumptions are never automatically
#'   verified, and no bootstrap is requested by this wrapper. The local-polynomial
#'   option adds correction_fit diagnostics and inference_inputs containing
#'   complete caller-use arguments and separately bound actual score tangents.
#' @seealso \code{wm_wdsm_fit()}, \code{wm_match()}, \code{wm_fit()},
#'   \code{wm_fitted_inference()}
#' @export

wm_model_fit <- function(Y, Z, weights, ps_models = list(),
                          pg0_models = list(), pg1_models = NULL, M = 3L,
                          estimand = c("PATE", "PATT"),
                          inference = c("none", "full_x"), controls = list(),
                          conf.level = 0.95, max_stack_elements = 5e7,
                          tie_rule = c("row_order", "source_random"),
                          tie_seed = NULL,
                          tie_tolerance = 64 * .Machine$double.eps,
                          correction = NULL) {
  call <- match.call()
  estimand <- match.arg(estimand)
  inference <- match.arg(inference)
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
  if (estimand == "PATT" && !is.null(pg1_models)) {
    stop("PATT does not fit an unused treated prognostic model.", call. = FALSE)
  }
  n <- length(Y)
  ps <- .wm_model_specs(ps_models, n, "ps", "ps_models")
  pg <- list("0" = .wm_model_specs(pg0_models, n, "pg", "pg0_models"))
  if (estimand == "PATE") {
    pg[["1"]] <- .wm_model_specs(pg1_models, n, "pg", "pg1_models")
  }
  used <- if (estimand == "PATE") 0:1 else 0L
  layout <- .wm_model_layout(ps, pg, used)
  if (!is.null(correction)) {
    if (is.list(correction) &&
        identical(correction[["method", exact = TRUE]], "quadratic_qr")) {
      return(.wm_model_quadratic_qr_fit(Y, Z, weights, ps, pg, layout,
        correction, M, estimand, inference, controls, conf.level, max_stack_elements,
        ties, call))
    }
    return(.wm_model_local_polynomial_fit(Y, Z, weights, ps, pg, layout,
      correction, M, estimand, inference, controls, conf.level, max_stack_elements,
      ties, call))
  }
  allocation <- .wm_model_allocation(n, ps, pg, layout, max_stack_elements)
  recipe <- list(J = length(ps), K = lengths(pg),
    matching_dimensions = layout$matching_dimensions,
    model_names = layout$model_names, shared_ps_fitted_once = TRUE,
    supplied_weights_known = TRUE, correction = "Complete weighted quadratic",
    legacy_double_score_reduction = length(ps) == 1L && all(lengths(pg) == 1L))
  # Exact source specialization delegates to the existing reviewed wrapper.
  # Core fit/nuisance/inference/solver payloads retain that path's arithmetic.
  if (recipe$legacy_double_score_reduction) {
    answer <- wm_wdsm_fit(Y, Z, weights,
      ps_design = ps[[1L]]$design, pg0_design = pg[["0"]][[1L]]$design,
      pg1_design = if (estimand == "PATE") pg[["1"]][[1L]]$design else NULL,
      M = M, estimand = estimand, inference = inference, controls = controls,
      conf.level = conf.level, ps_offset = ps[[1L]]$offset,
      pg0_offset = pg[["0"]][[1L]]$offset,
      pg1_offset = if (estimand == "PATE") pg[["1"]][[1L]]$offset else NULL,
      ps_weighting = ps[[1L]]$weighting,
      tie_rule = ties$rule, tie_seed = ties$seed, tie_tolerance = ties$tolerance)
    answer$call <- call
    answer$model_recipe <- recipe
    answer$allocation <- allocation
    class(answer) <- c("wm_model_fit", class(answer))
    return(answer)
  }
  dimensions <- layout$matching_dimensions
  if (inference == "full_x" &&
      !(all(dimensions >= 2L) || all(dimensions == 1L))) {
    stop("full_x uses the existing all-vector or all-scalar inference branches; mixed scalar/vector maps are unsupported.",
         call. = FALSE)
  }
  arguments <- list(Y = Y, Z = Z, weights = weights, ps_models = ps,
    pg0_models = pg[["0"]],
    pg1_models = if (estimand == "PATE") pg[["1"]] else NULL,
    estimand = estimand, max_stack_elements = max_stack_elements)
  fitter <- if (inference == "none") .wm_model_prediction_stack else
    .wm_model_nuisance_stack
  solved <- .wm_wdsm_fit_stack(arguments, fitter, controls)
  stack <- solved$stack
  if (inference == "full_x") {
    assembly <- .wm_model_fitted_pipeline(stack, M = M, conf.level = conf.level,
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
  scope <- paste("Fixed finite candidate logistic PS and arm OLS lists with known",
    "supplied positive weights. Shared PS coefficients fitted once; no duplicate",
    "or dependent outputs deleted. Pooled ordinary empirical centers/variances,",
    "complete weighted quadratic correction and original fixed-M self-normalized",
    "PATE/PATT fractions. All generated-coordinate Jacobian/IF covariance blocks",
    "retained. Fitted model identification, root regularity, prediction correctness,",
    "weighted model-union compatibility, support and transport premises remain",
    "unverified. full_x is conditional on full-X-correct population predictions,",
    "not automatic validity under a correct PS candidate. No original refit",
    "bootstrap or dependent cluster/stratum inference is asserted.")
  structure(list(call = call, estimate = fit$estimate,
    raw_estimate = fit$raw_estimate, correction = fit$correction,
    n = fit$n, M = fit$M, estimand = estimand, ps_weighting = stack$ps_weighting,
    se = inference_result$se, variance = inference_result$variance,
    root_n_variance = inference_result$root_n_variance,
    conf.int = inference_result$conf.int, inference = inference_result,
    fit = fit, nuisance = stack, assembly = assembly, model_recipe = recipe,
    allocation = allocation,
    solver = list(controls = controls, attempts = solved$attempts,
      refined = solved$refined, normalized_ps_threshold = 1e-10,
      newton_threshold = min(1e-9, 1/fit$n)),
    assumptions_verified = FALSE, bootstrap_requested = FALSE,
    scope = scope), class = c("wm_model_fit", "list"))
}
