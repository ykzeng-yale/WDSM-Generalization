# Post-fit binding only. This file never fits a model, rebuilds a matching graph,
# or changes a supplied individual weight.

.wm_mi_unavailable <- function(object, scope, level, stage, reason) {
  structure(list(estimate = object$fit$estimate, n = object$fit$n,
    M = object$fit$M, estimand = object$fit$estimand, fit = object$fit,
    covariance_scope = scope, status = paste0("model_", stage, "_unavailable"),
    available = FALSE, numerically_available = FALSE,
    conditional_inference_available = FALSE, sampling_inference_available = FALSE,
    failure_stage = stage, unavailable_reason = reason,
    root_n_variance = NA_real_, variance = NA_real_, se = NA_real_,
    conf.int = c(lower = NA_real_, upper = NA_real_), conf.level = level,
    assumptions_verified = FALSE, application_verified = FALSE,
    original_refit_limit_agreement_declared = FALSE,
    contract = paste("The original fitted point is retained. No fit, graph, weight,",
      "correction or derivative is repaired or substituted; no sampling interval",
      "or original-refit bootstrap claim is made.")),
    class = c("wm_fitted_inference", "list"))
}

.wm_mi_finish <- function(object, result, request, call, dispatched, parameters = NULL) {
  previous <- object$inference_handoff
  original_inference <- if (is.null(previous)) object$inference else previous$original_inference
  original_assembly <- if (is.null(previous)) object$assembly else previous$original_assembly
  object$inference <- result
  for (field in c("root_n_variance", "variance", "se", "conf.int")) {
    object[[field]] <- result[[field]]
  }
  # A historical full-X assembly is not an active assembly for another scope.
  # Preserve it once as provenance, without nesting repeated handoffs.
  object$assembly <- NULL
  object$inference_handoff <- list(call = call, request = request,
    dispatched = dispatched, original_inference = original_inference,
    original_assembly = original_assembly, parameter_names = parameters,
    fit_preserved = TRUE, nuisance_parameters_preserved = TRUE, models_refitted = FALSE,
    influence_completion = object$nuisance$influence_completion,
    graph_rebuilt = FALSE, individual_weights_estimated = FALSE,
    assumptions_verified = FALSE,
    correction_contract = if (!is.null(object$correction_fit)) paste(
      "Local-polynomial score-only IF and B'gradient reference derivative.",
      "This is not differentiation through refitting local coefficients.",
      "The adopted correction transfer and derivative consistency use bounded",
      "outcomes and positive bounded known W plus their smoothness/support/rate",
      "premises; this handoff does not inherit unbounded Gaussian/full-X alternatives.",
      "Variance-row, slope, transport, covariance-kernel and joint sampling laws",
      "retain their separate premises and are not verified by finite arrays.") else paste(
      "Complete retained quadratic/prediction root and its actual derivatives.",
      "The selected low-level sampling and covariance premises remain unverified."),
    contract = paste("Post-fit dispatch to wm_fitted_inference with protected original",
      "fit, complete IF, reference derivatives and exact-zero weight derivative.",
      "Covariance/transport choices are caller declarations. Numerical availability",
      "does not verify population assumptions; no bootstrap is run."))
  object$assumptions_verified <- FALSE
  object
}

.wm_mi_structure <- function(object) {
  object <- .wm_reciprocal_list(object, "object")
  if (!inherits(object, "wm_model_fit")) {
    stop("object must be a retained wm_model_fit result.", call. = FALSE)
  }
  fit <- .wm_reciprocal_list(object$fit, "object$fit")
  if (!inherits(fit, "wm_match") || !identical(fit$method, "self_normalized")) {
    stop("The original self_normalized matching fit is required.", call. = FALSE)
  }
  if (!identical(object$estimand, fit$estimand) ||
      !fit$estimand %in% c("PATE", "PATT")) {
    stop("Model and original fit estimands disagree.", call. = FALSE)
  }
  for (name in c("n", "M", "estimate", "raw_estimate", "correction")) {
    .wm_reciprocal_agree(object[[name]], fit[[name]], paste("Model", name))
  }
  invisible(TRUE)
}

.wm_mi_bind <- function(object) {
  fit <- object$fit
  s <- .wm_reciprocal_list(object$nuisance, "object$nuisance")
  if (!inherits(s, ".wm_model_nuisance_stack") &&
      !inherits(s, ".wm_wdsm_nuisance_stack")) {
    stop("A retained complete model nuisance stack is required.", call. = FALSE)
  }
  n <- .wm_fit_integer(fit$n, "fit$n", 2L)
  pate <- identical(fit$estimand, "PATE")
  if (!identical(s$estimand, fit$estimand)) stop("Nuisance estimand disagrees with fit.", call. = FALSE)
  m <- .wm_numeric_vector(s$multiplicity, n, "nuisance multiplicity")
  if (any(m != 1)) stop("Post-fit sampling handoff requires original unit multiplicities.", call. = FALSE)
  for (name in c("Y", "Z")) {
    .wm_reciprocal_agree(s$inputs[[name]], fit$data[[name]], paste("Nuisance", name))
  }
  .wm_reciprocal_agree(s$inputs$weights, fit$weights, "Nuisance supplied W")
  parameters <- s$parameter_names
  if (!is.character(parameters) || !length(parameters) || anyNA(parameters) ||
      any(!nzchar(parameters)) || anyDuplicated(parameters) ||
      !identical(names(s$parameter), parameters)) {
    stop("Nuisance parameters must retain their complete names/order.", call. = FALSE)
  }
  arms <- if (pate) c("0", "1") else "0"
  lp <- !is.null(object$correction_fit)
  if (lp && (!identical(object$correction_fit$method, "local_polynomial") ||
             !isTRUE(s$score_only))) {
    stop("Local correction and score-only nuisance layout disagree.", call. = FALSE)
  }
  if (!lp && isTRUE(s$score_only)) stop("Score-only roots require their retained LP correction.", call. = FALSE)
  D <- list()
  for (arm in arms) {
    score <- .wm_score_matrix(s[[paste0("scores", arm)]], n, paste("Nuisance scores", arm))
    observed <- fit$graph[[paste0("scores", arm)]]
    if (!identical(dim(score), dim(observed))) stop("Nuisance and graph dimensions disagree.", call. = FALSE)
    .wm_reciprocal_agree(score, observed, paste("Nuisance score binding", arm))
    mu <- if (lp) object$correction_fit$arms[[arm]]$mean else s[[paste0("mean", arm)]]
    .wm_reciprocal_agree(mu, fit$predictions[[paste0("mean", arm)]], paste("Prediction binding", arm))
    D[arm] <- list(if (lp) object$correction_fit$arms[[arm]]$reference_mean_derivative else
      s[[paste0("mean", arm, "_derivative")]])
  }
  influence <- s$nuisance_influence
  values <- c(if (!is.null(influence)) list(nuisance_influence = influence),
              list(weight_derivative = s$weight_derivative),
              stats::setNames(D, paste0("mean_derivative", arms)))
  for (name in names(values)) {
    x <- values[[name]]
    if (is.null(x)) return(list(unavailable = paste("Missing retained", name, "is not replaced by zero or refitted.")))
    if (!is.matrix(x) || !is.numeric(x) || is.complex(x) ||
        !identical(dim(x), c(as.integer(n), as.integer(length(parameters)))) ||
        !identical(colnames(x), parameters)) {
      stop(name, " must retain the complete n-by-p parameter layout.", call. = FALSE)
    }
    .wm_fi_row_names(rownames(x), n, name)
    if (any(!is.finite(x))) return(list(unavailable = paste("Nonfinite retained", name, "prevents inference; the original point is retained.")))
  }
  if (any(s$weight_derivative != 0)) stop("Only known W with an exact-zero weight derivative is supported.", call. = FALSE)
  list(arguments = list(fit = fit, nuisance_influence = influence,
    mean_derivative0 = D[["0"]], mean_derivative1 = if (pate) D[["1"]] else NULL,
    weight_derivative = s$weight_derivative), parameters = parameters, arms = arms)
}

.wm_mi_complete_influence <- function(stack, n, parameters, maximum) {
  p <- length(parameters)
  planned <- 8 * as.double(n) * p + 8 * as.double(p)^2 + 4 * p
  if (!is.finite(planned) || as.double(n) * p > .Machine$integer.max ||
      as.double(p)^2 > .Machine$integer.max || planned > maximum) {
    return(list(unavailable = "Complete retained-equation IF/covariance exceeds max_influence_elements; no columns are removed and no fit is performed."))
  }
  psi <- stack$estimating_equations; J <- stack$jacobian
  if (is.null(psi) || is.null(J)) {
    return(list(unavailable = "The retained fit has neither a complete IF nor both complete equations/Jacobian; no refit is performed."))
  }
  if (!is.matrix(psi) || !is.numeric(psi) || is.complex(psi) ||
      !identical(dim(psi), c(as.integer(n), as.integer(p))) ||
      !identical(colnames(psi), parameters)) {
    stop("Retained estimating equations must have the complete original n-by-p layout.", call. = FALSE)
  }
  .wm_fi_row_names(rownames(psi), n, "Retained estimating equations")
  if (!is.matrix(J) || !is.numeric(J) || is.complex(J) ||
      !identical(dim(J), c(as.integer(p), as.integer(p))) ||
      !identical(rownames(J), parameters) || !identical(colnames(J), parameters)) {
    stop("Retained Jacobian must have the complete named p-by-p layout.", call. = FALSE)
  }
  probability <- .wm_numeric_vector(stack$empirical_probability, n, "Retained empirical probability")
  if (any(probability != 1 / n)) stop("IF completion requires original uniform empirical probabilities.", call. = FALSE)
  if (any(!is.finite(psi)) || any(!is.finite(J))) {
    return(list(unavailable = "Retained equations/Jacobian are nonfinite; no repair or refit is performed."))
  }
  # Exactly the existing complete-stack arithmetic at the stored parameter.
  # This solves for influence rows, never for new model coefficients.
  completed <- tryCatch({
    influence <- -t(solve(J, t(psi)))
    colnames(influence) <- parameters
    influence_mean <- as.vector(crossprod(probability, influence))
    centered <- sweep(influence, 2L, influence_mean, "-")
    Sigma <- crossprod(centered, centered * probability)
    dimnames(Sigma) <- list(parameters, parameters)
    if (any(!is.finite(influence)) || any(!is.finite(Sigma))) {
      stop("Retained-equation influence arithmetic exceeded numerical range.")
    }
    list(influence = influence, covariance = Sigma)
  }, error = identity)
  if (inherits(completed, "error")) {
    return(list(unavailable = paste("Retained-equation influence solve failed:",
      conditionMessage(completed), "No retry, parameter change or refit is performed.")))
  }
  completed$record <- list(method = "-t(solve(retained_J,t(retained_psi)))",
    original_unit_multiplicities = TRUE, original_uniform_probability = TRUE,
    parameter_names = parameters, planned_elements = planned,
    max_influence_elements = maximum, models_refitted = FALSE,
    contract = "Complete IF and centered Sigma at the retained root, with no coefficient update. Conservative element guard; not an RSS bound or a root/population certificate.")
  completed
}

.wm_mi_tangents <- function(object, arm, parameters) {
  fit_scores <- object$fit$graph[[paste0("scores", arm)]]
  n <- object$fit$n; d <- ncol(fit_scores); p <- length(parameters)
  derivatives <- object$nuisance$score_derivatives[[arm]]
  if (!is.list(derivatives) || length(derivatives) != d) {
    stop("Complete actual score derivatives are missing for the fitted sentinel.", call. = FALSE)
  }
  B <- array(0, c(as.integer(n), as.integer(d), as.integer(p)),
    dimnames = list(NULL, colnames(fit_scores), parameters))
  for (k in seq_len(d)) {
    B[, k, ] <- .wm_fi_matrix(derivatives[[k]], n, parameters,
                              paste("Actual score derivative", arm, k))
  }
  B
}

.wm_mi_transport <- function(spec, object, arm, parameters) {
  if (is.null(spec)) return(NULL)
  spec <- .wm_reciprocal_list(spec, paste0("transport", arm))
  if (identical(spec[["tangents", exact = TRUE]], "fitted")) {
    if (!identical(spec[["mode", exact = TRUE]], "estimate") ||
        !identical(spec[["representation", exact = TRUE]], "score_measurable_weight")) {
      stop("The fitted tangent sentinel requires explicit score_measurable_weight estimate mode; it is not a raw-chart derivative.", call. = FALSE)
    }
    if (!is.null(spec[["raw_to_matching", exact = TRUE]])) {
      stop("The fitted tangent sentinel uses actual matching coordinates and no raw_to_matching map.", call. = FALSE)
    }
    scores <- object$fit$graph[[paste0("scores", arm)]]
    if (!is.null(spec[["raw_scores", exact = TRUE]])) {
      stop("Omit raw_scores with the fitted tangent sentinel; actual graph coordinates are bound internally.", call. = FALSE)
    }
    spec$raw_scores <- scores
    spec$tangents <- .wm_mi_tangents(object, arm, parameters)
  }
  spec
}

.wm_mi_dispatch <- function(arguments) {
  # Existing low-level transport/critical numerical failures already return
  # unavailable objects. Only these explicit remaining arithmetic failures
  # become an unavailable model handoff; API/binding/invariant errors propagate.
  numerical <- c("Weight derivative normalization overflow.",
    "Gradient arithmetic exceeded numerical range.",
    "Fitted-variance arithmetic exceeded numerical range.",
    "Confidence interval overflow.",
    "Reciprocal fitted-variance arithmetic exceeded numerical range.",
    "Reciprocal fitted confidence interval exceeded numerical range.",
    "Normalized transport arithmetic exceeded numerical range.")
  tryCatch(do.call(wm_fitted_inference, arguments), error = function(e) {
    if (!conditionMessage(e) %in% numerical) stop(e)
    e
  })
}

#' Add explicitly qualified inference to a retained finite-model matching fit
#'
#' No model, local correction or graph is refitted. Scientific scope and transport
#' premises are supplied explicitly and are not verified from the retained arrays.
#' @export
wm_model_inference <- function(object, covariance_scope,
                               transport0 = NULL, transport1 = NULL,
                               critical_control = NULL, conf.level = 0.95,
                               max_influence_elements = 5e7) {
  call <- match.call()
  .wm_mi_structure(object)
  allowed <- c("distinct_rarity", "full_x", "common_field", "patt",
               "regular_joint_d_gt2", "critical_planar")
  if (missing(covariance_scope) || !is.character(covariance_scope) ||
      length(covariance_scope) != 1L || is.na(covariance_scope) ||
      !covariance_scope %in% allowed) {
    stop("Declare an existing fitted covariance scope explicitly; scalar_psm is not a generic model handoff.", call. = FALSE)
  }
  pate <- identical(object$fit$estimand, "PATE")
  if ((pate && covariance_scope == "patt") || (!pate && covariance_scope != "patt")) {
    stop("PATE and PATT require their respective covariance scopes.", call. = FALSE)
  }
  if (!pate && !is.null(transport1)) stop("PATT uses only transport0.", call. = FALSE)
  if (!is.null(critical_control) && covariance_scope != "critical_planar") {
    stop("critical_control is only supported with critical_planar.", call. = FALSE)
  }
  if (covariance_scope == "full_x" && (!is.null(transport0) || !is.null(transport1))) {
    stop("full_x supplies its declared zero transport; omit transport specifications.", call. = FALSE)
  }
  conf.level <- .wm_numeric_vector(conf.level, 1L, "conf.level")
  if (conf.level <= 0 || conf.level >= 1) stop("conf.level must lie in (0,1).", call. = FALSE)
  max_influence_elements <- .wm_numeric_vector(max_influence_elements, 1L,
    "max_influence_elements", positive = TRUE)
  request <- list(covariance_scope = covariance_scope, transport0 = transport0,
    transport1 = transport1, critical_control = critical_control, conf.level = conf.level,
    max_influence_elements = max_influence_elements)
  unavailable <- function(stage, reason) .wm_mi_finish(object,
    .wm_mi_unavailable(object, covariance_scope, conf.level, stage, reason),
    request, call, FALSE)
  if (!identical(object$fit$graph$tie_rule, "Exact distance, then original row index")) {
    return(unavailable("tie_rule", "Sampling handoff is unavailable for source_random or an unsupported tie policy; the point is retained."))
  }
  if (!is.null(object$correction_fit) &&
      !isTRUE(object$correction_fit$reference_derivatives_available)) {
    return(unavailable("correction", "Retained LP guard, boundary or reference-derivative failure prevents inference; no tuning, row deletion or refit is performed."))
  }
  if (!is.null(object$correction_fit)) {
    arms <- if (pate) c("0", "1") else "0"
    for (arm in arms) {
      a <- object$correction_fit$arms[[arm]]
      for (field in c("guard_failure", "derivative_boundary", "reference_derivative_failure",
                      "reference_derivative_available")) {
        flag <- a[[field]]
        if (!is.logical(flag) || length(flag) != object$fit$n || anyNA(flag)) {
          stop("LP diagnostics must retain all original-row flags.", call. = FALSE)
        }
      }
      if (any(a$guard_failure | a$derivative_boundary | a$reference_derivative_failure) ||
          !all(a$reference_derivative_available)) {
        return(unavailable("correction", "At least one original LP evaluation lacks a usable reference derivative; no row is removed or refitted."))
      }
    }
  }
  bound <- .wm_mi_bind(object)
  if (!is.null(bound$unavailable)) return(unavailable("inputs", bound$unavailable))
  parameters <- bound$parameters
  if (is.null(bound$arguments$nuisance_influence)) {
    completed <- .wm_mi_complete_influence(object$nuisance, object$fit$n,
      parameters, max_influence_elements)
    if (!is.null(completed$unavailable)) return(unavailable("influence", completed$unavailable))
    object$nuisance$nuisance_influence <- completed$influence
    object$nuisance$nuisance_covariance <- completed$covariance
    object$nuisance$influence_completion <- completed$record
    bound$arguments$nuisance_influence <- completed$influence
  }
  transport0 <- .wm_mi_transport(transport0, object, "0", parameters)
  if (pate) transport1 <- .wm_mi_transport(transport1, object, "1", parameters)
  if (!is.null(critical_control)) {
    critical_control <- .wm_reciprocal_list(critical_control, "critical_control")
    for (arm in c("0", "1")) {
      key <- paste0("tangents", arm)
      if (identical(critical_control[[key, exact = TRUE]], "fitted")) {
        critical_control[[key]] <- .wm_mi_tangents(object, arm, parameters)
      }
    }
  }
  # The common function validates transport/critical controls exactly once.
  # Do not prevalidate raw_to_matching here: caller callbacks can have state.
  # The dispatch handler propagates API/binding errors unchanged.
  arguments <- c(bound$arguments, list(transport0 = transport0, transport1 = transport1,
    covariance_scope = covariance_scope, critical_control = critical_control,
    conf.level = conf.level))
  result <- .wm_mi_dispatch(arguments)
  if (inherits(result, "error")) {
    result <- .wm_mi_unavailable(object, covariance_scope, conf.level,
      "dispatch", conditionMessage(result))
  }
  .wm_mi_finish(object, result, request, call, TRUE, parameters)
}
