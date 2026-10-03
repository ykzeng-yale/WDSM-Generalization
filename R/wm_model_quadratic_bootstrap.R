# Common-API marginal replication for the scoped nested quadratic law.
# All original point/root objects are retained; only local-index graphs vary.

.wm_qc_bootstrap_check <- function(object, conf.level) {
  if (!inherits(object, "wm_fitted_inference") ||
      !identical(object$status, "model_quadratic_contrast_inputs_ready") ||
      !isTRUE(object$conditional_law_inputs_available) ||
      !isTRUE(object$replication_available) ||
      !isTRUE(object$numerically_available) || isTRUE(object$available) ||
      isTRUE(object$sampling_inference_available)) {
    stop("Require prepared nested quadratic conditional-law inputs.", call. = FALSE)
  }
  inputs <- .wm_reciprocal_list(object$contrast_inputs, "quadratic contrast inputs")
  if (!isTRUE(inputs$available) || !isTRUE(inputs$numerically_available) ||
      !isTRUE(inputs$replication_available) ||
      !is.function(inputs$evaluate) ||
      !identical(formals(inputs$evaluate), formals(function(contrast) NULL)) ||
      !identical(body(inputs$evaluate), quote(.wm_qc_evaluate(retained, contrast)))) {
    stop("The retained quadratic evaluator binding is invalid.", call. = FALSE)
  }
  env <- environment(inputs$evaluate)
  if (is.null(env) || !exists("retained", envir = env, inherits = FALSE)) {
    stop("The quadratic evaluator lacks its retained context.", call. = FALSE)
  }
  context <- get("retained", envir = env, inherits = FALSE)
  if (!identical(context, inputs$context_binding)) {
    stop("The quadratic evaluator/context binding changed.", call. = FALSE)
  }
  if (!identical(object$fit, inputs$original_fit_binding)) {
    stop("The original quadratic fit binding changed.", call. = FALSE)
  }
  state <- .wm_reciprocal_validated_state(object$fit, conf.level)
  n <- as.integer(state$n); parameters <- inputs$full_parameters
  if (!is.character(parameters) || !length(parameters) || anyNA(parameters) ||
      any(!nzchar(parameters)) || anyDuplicated(parameters) ||
      !identical(parameters, context$parameters)) {
    stop("The complete quadratic parameter order changed.", call. = FALSE)
  }
  expected_scope <- if (object$fit$estimand == "PATE") "full_x" else "patt"
  if (!identical(object$covariance_scope, expected_scope) ||
      !identical(object$fit$method, "self_normalized")) {
    stop("Nested quadratic replication requires its original PATE/full_x or PATT/patt point.", call. = FALSE)
  }
  for (name in c("n", "M", "estimand", "estimate")) {
    if (!identical(context[[name]], object[[name]]) ||
        !identical(context[[name]], object$fit[[name]]) ||
        (name != "estimate" && !identical(inputs[[name]], context[[name]]))) {
      stop("The quadratic original point/count/estimand binding changed.", call. = FALSE)
    }
  }
  for (name in c("Y", "Z")) .wm_reciprocal_agree(context[[name]], object$fit$data[[name]],
    paste("Quadratic source", name))
  .wm_reciprocal_agree(context$W, object$fit$weights, "Quadratic raw W")
  .wm_reciprocal_agree(context$w, object$fit$analysis_weights, "Quadratic analysis W")
  .wm_reciprocal_agree(context$gamma, object$fit$gamma, "Quadratic target denominator")
  phi <- .wm_fi_matrix(context$phi, n, parameters, "Complete quadratic influence")
  .wm_reciprocal_agree(phi, sweep(phi, 2L, colMeans(phi), "-"), "Centered quadratic influence")
  indices <- match(inputs$contrast_parameters, parameters)
  if (anyNA(indices) || !length(indices) || anyDuplicated(indices) ||
      !identical(inputs$contrast_parameters, context$contrast_names) ||
      inputs$dimension != length(indices)) stop("Quadratic contrast order changed.", call. = FALSE)
  for (name in c("contrast_covariance", "full_score_covariance")) {
    axes <- if (name == "contrast_covariance") inputs$contrast_parameters else parameters
    value <- inputs[[name]]
    if (!is.matrix(value) || !identical(dim(value), c(length(axes), length(axes))) ||
        !identical(rownames(value), axes) || !identical(colnames(value), axes)) {
      stop("Quadratic covariance parameter axes changed.", call. = FALSE)
    }
  }
  .wm_reciprocal_agree(context$H, phi[, indices, drop = FALSE], "Quadratic full-root contrast")
  .wm_reciprocal_agree(context$Omega, crossprod(context$H, context$H / n), "Quadratic contrast covariance")
  .wm_reciprocal_agree(inputs$contrast_covariance, context$Omega, "Quadratic covariance binding")
  .wm_reciprocal_agree(inputs$full_score_covariance, crossprod(phi, phi / n), "Quadratic full covariance")
  .wm_reciprocal_agree(inputs$contrast_information_map, context$information_map, "Quadratic information map")
  if (!identical(inputs$replication_control, context$replication_control) ||
      !identical(names(context$replication_control),
        c("maximum_draw_entries", "maximum_distance_evaluations"))) {
    stop("The declared quadratic allocation controls changed.", call. = FALSE)
  }
  factor <- context$contrast_qr
  d <- length(indices)
  if (!identical(factor, inputs$contrast_factorization) ||
      !is.matrix(factor$Q) || !identical(dim(factor$Q), c(n, d)) ||
      !is.matrix(factor$R) || !identical(dim(factor$R), c(d, d)) ||
      !identical(sort(factor$pivot), seq_len(d)) ||
      any(!is.finite(c(factor$Q, factor$R))) || any(diag(factor$R) == 0)) {
    stop("The retained full contrast QR binding changed.", call. = FALSE)
  }
  .wm_reciprocal_agree(factor$Q %*% factor$R,
    context$H[, factor$pivot, drop = FALSE] / sqrt(n), "Quadratic contrast QR reconstruction")
  .wm_reciprocal_agree(crossprod(factor$Q), diag(d), "Quadratic contrast QR orthogonality")
  for (name in c("maximum_draw_entries", "maximum_distance_evaluations")) {
    .wm_numeric_vector(context$replication_control[[name]], 1L, name, positive = TRUE)
  }
  arms <- if (object$fit$estimand == "PATE") c("0", "1") else "0"
  if (!identical(context$arms, arms) ||
      !identical(inputs$matching_dimensions, context$layout$matching_dimensions) ||
      any(inputs$matching_dimensions != 3L)) stop("Quadratic arm maps changed.", call. = FALSE)
  centered <- sweep(context$raw, 2L, colMeans(context$raw), "-")
  standardized <- sweep(centered, 2L, sqrt(colMeans(centered^2)), "/")
  for (arm in arms) {
    .wm_reciprocal_agree(unname(standardized[, context$layout$scores[[arm]], drop = FALSE]),
      object$fit$graph[[paste0("scores", arm)]], paste("Quadratic original score", arm))
    .wm_reciprocal_agree(context$mu[[arm]], context$raw[, context$layout$pg_names[[arm]][1L]],
      paste("Quadratic original PG", arm))
  }
  context
}

.wm_qc_bootstrap <- function(object, B, seed, conf.level, interval, chunk_size, counts) {
  if (!is.null(counts)) stop("Nested quadratic replication rejects counts; use its Gaussian contrast law.", call. = FALSE)
  if (length(interval) > 1L && identical(interval, c("normal", "basic", "none"))) interval <- "basic"
  interval <- match.arg(interval, c("basic", "none"))
  B <- .wm_cp_integer(B, "B"); chunk_size <- .wm_cp_integer(chunk_size, "chunk_size")
  .wm_numeric_vector(conf.level, 1L, "conf.level")
  if (conf.level <= 0 || conf.level >= 1) stop("conf.level must lie in (0,1).", call. = FALSE)
  if (!is.null(seed)) {
    .wm_numeric_vector(seed, 1L, "seed")
    if (seed < 0 || seed > .Machine$integer.max || seed != floor(seed)) stop("Invalid seed.", call. = FALSE)
    if (identical(RNGkind()[2L], "Box-Muller")) stop("An explicit seed cannot preserve Box-Muller cache.", call. = FALSE)
  }
  # Original state, full covariance, cached full QR and allocations bind before RNG.
  context <- .wm_qc_bootstrap_check(object, conf.level)
  d <- length(context$contrast_names); n <- context$n
  distance_count <- sum(vapply(context$arms, function(arm)
    sum(context$Z == as.integer(arm)) * as.double(sum(context$Z != as.integer(arm))), numeric(1L)))
  planned <- B * as.double(2 * d + 18L) + 4 * as.double(d)^2
  if (!is.finite(planned) || B * as.double(d) > .Machine$integer.max ||
      planned > context$replication_control$maximum_draw_entries ||
      B * distance_count > context$replication_control$maximum_distance_evaluations) {
    stop("Quadratic all-B output/distance budget exceeded before RNG; no draws are filtered.", call. = FALSE)
  }
  factor <- context$contrast_qr
  roots <- eta <- rep(NA_real_, B)
  H <- zeta <- matrix(NA_real_, d, B)
  rownames(H) <- context$contrast_names
  rownames(zeta) <- paste0("qr_normal", seq_len(d))
  diagnostic <- data.frame(draw = seq_len(B), status = rep("not_evaluated", B),
    conditional_mean = rep(NA_real_, B), conditional_variance = rep(NA_real_, B),
    failure_stage = rep("", B), error = rep("", B), stringsAsFactors = FALSE)
  rng_kind <- RNGkind()
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    previous <- if (had_seed) get(".Random.seed", .GlobalEnv) else NULL
    on.exit({
      if (had_seed) assign(".Random.seed", previous, .GlobalEnv)
      else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(as.integer(seed))
  }
  for (draw in seq_len(B)) {
    normals <- numeric(d + 1L)
    for (block in .wm_cp_chunks(d + 1L, chunk_size)) normals[block] <- stats::rnorm(length(block))
    zeta[, draw] <- normals[seq_len(d)]; eta[draw] <- normals[d + 1L]
    # H/sqrt(n) [,pivot] = Q R, hence H* [pivot] = R' zeta.
    H[factor$pivot, draw] <- as.vector(crossprod(factor$R, zeta[, draw]))
    if (any(!is.finite(c(H[, draw], eta[draw])))) {
      diagnostic$status[draw] <- "contrast_range_failed"
      diagnostic$failure_stage[draw] <- "arithmetic"
      diagnostic$error[draw] <- "Generated contrast or innovation is nonfinite."
      next
    }
    value <- tryCatch(object$contrast_inputs$evaluate(H[, draw]), error = identity)
    if (inherits(value, "error")) {
      diagnostic$status[draw] <- "numerical_evaluation_failed"
      diagnostic$failure_stage[draw] <- "evaluation"
      diagnostic$error[draw] <- conditionMessage(value)
      next
    }
    if (!isTRUE(value$available)) {
      diagnostic$status[draw] <- "numerical_evaluation_failed"
      diagnostic$failure_stage[draw] <- if (is.null(value$failure_stage)) "evaluation" else value$failure_stage
      diagnostic$error[draw] <- if (is.null(value$unavailable_reason)) "Unavailable evaluated inputs." else value$unavailable_reason
      next
    }
    m <- value$conditional_mean; v <- value$conditional_variance
    if (!is.numeric(m) || is.complex(m) || length(m) != 1L || !is.finite(m) ||
        !is.numeric(v) || is.complex(v) || length(v) != 1L || !is.finite(v) || v < 0) {
      diagnostic$status[draw] <- "conditional_input_failed"
      diagnostic$failure_stage[draw] <- "arithmetic"
      diagnostic$error[draw] <- "Conditional mean or variance is invalid."
      next
    }
    diagnostic$conditional_mean[draw] <- m; diagnostic$conditional_variance[draw] <- v
    roots[draw] <- m + sqrt(v) * eta[draw]
    diagnostic$status[draw] <- if (is.finite(roots[draw])) "complete" else "root_range_failed"
    if (!is.finite(roots[draw])) {
      diagnostic$failure_stage[draw] <- "arithmetic"
      diagnostic$error[draw] <- "Root draw is nonfinite."
    }
  }
  draws <- object$estimate + roots / sqrt(n)
  shifted_failure <- diagnostic$status == "complete" & !is.finite(draws)
  diagnostic$status[shifted_failure] <- "shifted_draw_range_failed"
  diagnostic$failure_stage[shifted_failure] <- "arithmetic"
  diagnostic$error[shifted_failure] <- "Shifted estimate draw is nonfinite."
  complete <- all(diagnostic$status == "complete") && all(is.finite(c(roots, draws)))
  mc <- if (complete && B > 1L) stats::var(roots) else NA_real_
  arithmetic_ok <- complete && (B == 1L || (is.finite(mc) && mc >= 0))
  available <- arithmetic_ok && B > 1L
  alpha <- 1 - conf.level
  ci <- if (interval == "none") NULL else if (available) object$estimate -
    as.numeric(stats::quantile(roots, c(1 - alpha / 2, alpha / 2), names = FALSE)) / sqrt(n) else
    rep(NA_real_, 2L)
  if (available && any(!is.finite(ci))) {
    available <- FALSE; arithmetic_ok <- FALSE
  }
  if (!available) {
    mc <- NA_real_
    if (!is.null(ci)) ci[] <- NA_real_
  }
  if (!is.null(ci)) names(ci) <- c("lower", "upper")
  status <- if (!arithmetic_ok) "numerical_inference_unavailable" else if (B == 1L)
    "insufficient_replicates" else "complete"
  list(root_n_draws = roots, draws = draws, contrast_draws = H,
    standard_contrast_normals = zeta, scalar_innovations = eta,
    draw_diagnostics = diagnostic, available = available, status = status,
    sampling_inference_available = available, draws_available = complete,
    unavailable_reason = if (!complete) "At least one requested slot failed; see draw_diagnostics." else
      if (!arithmetic_ok) "All-slot variance or interval arithmetic is nonfinite." else
      if (B == 1L) "B=1 supplies draws only." else NULL,
    conditional_root_n_variance = NA_real_, conditional_variance = NA_real_,
    conditional_variance_formula = "Unevaluated full finite-data mixture variance; its existence/convergence is not asserted.",
    monte_carlo_root_n_variance = mc, monte_carlo_variance = mc / n,
    root_n_variance = mc, variance = mc / n,
    se = if (available) sqrt(mc / n) else NA_real_,
    conf.int = ci, interval = interval, conf.level = conf.level, B = B, n = n,
    seed = seed, rng_kind = rng_kind, estimate = object$estimate, estimand = object$estimand,
    covariance_scope = object$covariance_scope, source_inference = object,
    method = "Nested quadratic contrast Gaussian-mixture replication",
    draw_input = "generated_untruncated_nested_contrast_mixture", original_refit_bootstrap = FALSE,
    original_fit_preserved = TRUE, models_refitted = FALSE,
    numerical_error_rate_verified = FALSE, assumptions_verified = FALSE,
    population_assumptions_verified = FALSE, variance_divisor = "B-1",
    allocation = list(planned_draw_entries = planned, all_B_distance_evaluations = B * distance_count,
      controls = context$replication_control,
      contract = "Computational allocation/work caps before RNG; not statistical selection, an RSS bound or an asymptotic success guarantee."),
    replication_contract = paste("Scoped bounded outcome-correct nested weighted logistic PS",
      "and arm-specific linear PG law; full score/correction covariance, known W and original point.",
      "Within each slot all contrast normals precede the independent scalar innovation.",
      "Every requested slot is retained and attempted; no radial/Gram cutoff, zero completion,",
      "column deletion, rejection, redraw or survivor summaries. Numerical failures invalidate",
      "all aggregate inference. B=1 supplies draws only. All-B B-1 sample variance includes",
      "conditional-mean variation; basic intervals use uncentered root quantiles.",
      "Limit-law quantiles require B_n increasing and continuous unique limiting quantiles;",
      "empirical variance requires predeclared deterministic B_n increasing at most polynomially.",
      "Same-draw numerical RMS and root-n original-point numerical error conditions, population",
      "and sampling/input premises remain unverified. No full finite-data Var* convergence,",
      "normal interval, finite-sample coverage or original fixed-reuse validity claim."))
}
