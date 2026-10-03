# Private qualified package candidate. Uses reviewed package helpers unchanged.

.wm_fi_row_names <- function(x, n, label) {
  if (!is.null(x) && !identical(x, as.character(seq_len(n)))) {
    stop(label, " row labels must be original row indices 1,...,n in order.",
         call. = FALSE)
  }
}

.wm_fi_matrix <- function(x, n, parameter_names, label, zero = FALSE) {
  p <- length(parameter_names)
  if (is.null(x) && zero) {
    return(matrix(0, n, p, dimnames = list(NULL, parameter_names)))
  }
  x <- .wm_weight_matrix(x, n, label, p)
  if (!identical(colnames(x), parameter_names)) {
    stop(label, " columns must use the exact influence parameter names/order.",
         call. = FALSE)
  }
  .wm_fi_row_names(rownames(x), n, label)
  x
}

.wm_fi_transport_spec <- function(spec, scores, n, parameter_names, label) {
  spec <- .wm_reciprocal_list(spec, label)
  mode <- spec[["mode", exact = TRUE]]
  if (!is.character(mode) || length(mode) != 1L || is.na(mode) ||
      !mode %in% c("zero", "estimate")) {
    stop(label, " must declare mode = 'zero' or 'estimate'.", call. = FALSE)
  }
  if (identical(mode, "zero")) {
    if (!setequal(names(spec), c("mode", "basis")) ||
        !is.character(spec$basis) || length(spec$basis) != 1L ||
        is.na(spec$basis) ||
        !spec$basis %in% c("fixed_map", "current_centering")) {
      stop(label, " zero mode requires only mode and basis = 'fixed_map' or ",
           "'current_centering'.", call. = FALSE)
    }
    return(spec)
  }
  required <- c("mode", "raw_scores", "tangents", "cutoff")
  optional <- c("raw_to_matching", "bandwidth", "density_floor",
                "quadrature_tolerance", "max_nodes", "representation")
  if (!all(required %in% names(spec)) ||
      any(!names(spec) %in% c(required, optional))) {
    stop(label, " estimate mode has missing or unsupported fields; data, weights, ",
         "residuals, donor arm and M must come from fit.", call. = FALSE)
  }
  representation <- spec[["representation", exact = TRUE]]
  if (is.null(representation)) representation <- "conditional_weight_chart"
  if (!is.character(representation) || length(representation) != 1L ||
      is.na(representation) ||
      !representation %in% c("conditional_weight_chart", "score_measurable_weight")) {
    stop(label, " representation must be conditional_weight_chart or score_measurable_weight.",
         call. = FALSE)
  }
  score_weight <- identical(representation, "score_measurable_weight")
  raw <- spec$raw_scores
  .wm_score_matrix(raw, n, paste(label, "raw_scores"))
  d <- ncol(scores)
  if (ncol(raw) != d) stop(label, " raw and matching dimensions differ.", call. = FALSE)
  .wm_fi_row_names(rownames(raw), n, paste(label, "raw_scores"))
  tangent <- spec$tangents
  if (!is.numeric(tangent) || is.complex(tangent) ||
      !identical(dim(tangent), c(as.integer(n), as.integer(d),
                                 as.integer(length(parameter_names)))) ||
      any(!is.finite(tangent))) {
    stop(label, " tangents must be a finite n-by-d-by-p array.", call. = FALSE)
  }
  if (!identical(dimnames(tangent)[[3L]], parameter_names)) {
    stop(label, " tangent parameter names/order must match influence columns.",
         call. = FALSE)
  }
  .wm_fi_row_names(dimnames(tangent)[[1L]], n, paste(label, "tangents"))
  raw_names <- colnames(raw)
  if (!is.null(raw_names) &&
      (anyNA(raw_names) || any(!nzchar(raw_names)) || anyDuplicated(raw_names))) {
    stop(label, " raw coordinate names must be unique and nonempty.", call. = FALSE)
  }
  if (score_weight && !identical(raw_names, colnames(scores))) {
    stop(label, " score_measurable_weight coordinate names/order must match fit scores.",
         call. = FALSE)
  }
  if (!identical(dimnames(tangent)[[2L]], raw_names)) {
    stop(label, " tangent coordinate names/order must match raw_scores (or both be unnamed).",
         call. = FALSE)
  }
  .wm_numeric_vector(spec$cutoff, n, paste(label, "cutoff"))
  if (any(spec$cutoff < 0 | spec$cutoff > 1)) {
    stop(label, " cutoff values must lie in [0,1].", call. = FALSE)
  }
  .wm_fi_row_names(names(spec$cutoff), n, paste(label, "cutoff"))
  map <- spec[["raw_to_matching", exact = TRUE]]
  if (score_weight && !is.null(map)) {
    stop(label, " score_measurable_weight requires identity matching coordinates; ",
         "omit raw_to_matching or use NULL.", call. = FALSE)
  }
  if (!is.null(map) && !is.function(map)) {
    stop(label, " raw_to_matching must be NULL or a function.", call. = FALSE)
  }
  mapped <- if (is.null(map)) raw else tryCatch(map(raw), error = function(e) {
    stop(label, " raw_to_matching failed: ", conditionMessage(e), call. = FALSE)
  })
  mapped <- .wm_score_matrix(mapped, n, paste(label, "mapped raw_scores"))
  if (!identical(dim(mapped), dim(scores))) {
    stop(label, " mapped raw_scores have the wrong matching dimension.", call. = FALSE)
  }
  .wm_reciprocal_agree(mapped, scores, paste(label, "raw-to-matching data binding"))
  defaults <- list(bandwidth = NULL, density_floor = NULL,
                   quadrature_tolerance = 1e-7, max_nodes = 131073L)
  controls <- lapply(names(defaults), function(name) {
    if (name %in% names(spec)) spec[[name, exact = TRUE]] else defaults[[name]]
  })
  names(controls) <- names(defaults)
  for (name in names(controls)) {
    value <- controls[[name]]
    if (is.null(value) && name %in% c("bandwidth", "density_floor")) next
    .wm_numeric_vector(value, 1L, paste(label, name), positive = TRUE)
    if (name == "max_nodes" && (value < 3 || value != floor(value))) {
      stop(label, " max_nodes must be an integer at least three.", call. = FALSE)
    }
  }
  result <- list(mode = mode, raw_scores = raw, tangents = tangent,
       cutoff = spec$cutoff, controls = controls,
       map_binding = list(method = if (is.null(map)) "identity" else "supplied_function",
         max_absolute_difference = max(abs(mapped - scores)),
         componentwise_tolerance = 1e-10,
         interpretation = paste("Observed coordinates bound to this fit in original row order;",
           "population invertibility, chart validity and tangent interpretation remain unverified.")))
  if (score_weight) {
    result$representation <- representation
    result$map_binding$interpretation <- paste("Actual matching coordinates bound to this fit",
      "in original row and coordinate order. Tangents must be the complete actual score",
      "derivative; population score measurability and derivative validity remain unverified.")
  }
  result
}

.wm_fi_zero_transport <- function(parameter_names, basis) {
  zero <- stats::setNames(numeric(length(parameter_names)), parameter_names)
  list(mode = "zero", basis = basis, graph_drift = zero,
       quadrature_error_bound = zero, status = "zero_by_declared_premise",
       numerically_available = TRUE, assumptions_verified = FALSE,
       contract = "No kernel estimated; zero transport follows only under the declared premise.")
}

#' Qualified known-weight fitted-map inference with complete nuisance slopes
#'
#' The PATE regular_joint_d_gt2 scope additionally retains signed reciprocal
#' covariance under the bounded equal-d>2 regular baseline joint-sheet theorem.
#' Its common bootstrap dispatch uses corrected-variance Gaussian calibration;
#' ordinary augmented-row/count multipliers are unsupported for this scope.
#' @export
wm_fitted_inference <- function(fit, nuisance_influence = NULL,
                                mean_derivative0 = NULL,
                                mean_derivative1 = NULL,
                                weight_derivative = NULL,
                                transport0 = NULL, transport1 = NULL,
                                covariance_scope, conf.level = 0.95,
                                scalar_control = NULL, critical_control = NULL,
                                scalar_model = NULL) {
  if (!is.null(scalar_model) &&
      (missing(covariance_scope) || !identical(covariance_scope, "scalar_psm"))) {
    stop("scalar_model is only supported with covariance_scope = 'scalar_psm'.",
         call. = FALSE)
  }
  if (!is.null(critical_control) &&
      (missing(covariance_scope) || !identical(covariance_scope, "critical_planar"))) {
    stop("critical_control is only supported with critical_planar.", call. = FALSE)
  }
  if (!missing(covariance_scope) && identical(covariance_scope, "scalar_psm")) {
    if (!is.null(nuisance_influence) || !is.null(mean_derivative0) ||
        !is.null(mean_derivative1) || !is.null(weight_derivative) ||
        !is.null(transport0) || !is.null(transport1)) {
      stop("scalar_psm derives its own influence and drift; omit supplied derivatives and transport.",
           call. = FALSE)
    }
    return(.wm_scalar_psm_inference(fit, scalar_control, conf.level, scalar_model))
  }
  if (!is.null(scalar_control)) {
    stop("scalar_control is only supported with covariance_scope = 'scalar_psm'.",
         call. = FALSE)
  }
  fit <- .wm_reciprocal_list(fit, "fit")
  if (!inherits(fit, "wm_match") || !identical(fit$method, "self_normalized")) {
    stop("Supply an unadjusted self_normalized wm_match or wm_fit object.", call. = FALSE)
  }
  n <- .wm_fit_integer(fit$n, "fit$n", 2L)
  if (!is.character(fit$estimand) || length(fit$estimand) != 1L ||
      is.na(fit$estimand) || !fit$estimand %in% c("PATE", "PATT")) {
    stop("Unsupported estimand.", call. = FALSE)
  }
  pate <- identical(fit$estimand, "PATE")
  if (missing(covariance_scope) || !is.character(covariance_scope) ||
      length(covariance_scope) != 1L || is.na(covariance_scope) ||
      !covariance_scope %in% c("distinct_rarity", "full_x", "common_field", "patt",
                               "regular_joint_d_gt2", "critical_planar") ||
      (pate && covariance_scope == "patt") || (!pate && covariance_scope != "patt")) {
    stop("Declare PATE distinct_rarity/full_x/common_field/regular_joint_d_gt2/critical_planar, or PATT patt explicitly.",
         call. = FALSE)
  }
  if (!pate && (!is.null(mean_derivative1) || !is.null(transport1))) {
    stop("PATT uses only control prediction derivatives and transport0.", call. = FALSE)
  }
  conf.level <- .wm_numeric_vector(conf.level, 1L, "conf.level")
  if (conf.level <= 0 || conf.level >= 1) stop("conf.level must lie in (0,1).", call. = FALSE)
  influence <- .wm_weight_matrix(nuisance_influence, n, "nuisance_influence")
  parameters <- colnames(influence)
  if (is.null(parameters) || anyNA(parameters) || any(!nzchar(parameters)) ||
      anyDuplicated(parameters)) {
    stop("nuisance_influence must have unique, nonempty parameter names.", call. = FALSE)
  }
  .wm_fi_row_names(rownames(influence), n, "nuisance_influence")
  D0 <- .wm_fi_matrix(mean_derivative0, n, parameters, "mean_derivative0", TRUE)
  D1 <- if (pate) .wm_fi_matrix(mean_derivative1, n, parameters, "mean_derivative1", TRUE) else NULL
  DW <- .wm_fi_matrix(weight_derivative, n, parameters, "weight_derivative", TRUE)
  if (any(DW != 0)) {
    stop("wm_fitted_inference treats supplied weights as known; estimated-weight ",
         "inference is outside this WM scope. Use NULL or exact zero weight_derivative.",
         call. = FALSE)
  }
  score0 <- .wm_score_matrix(fit$graph$scores0, n, "fit scores0")
  score1 <- if (pate) .wm_score_matrix(fit$graph$scores1, n, "fit scores1") else NULL
  dimensions <- c(potential0 = ncol(score0), if (pate) c(potential1 = ncol(score1)))
  if (covariance_scope == "regular_joint_d_gt2" &&
      (!pate || length(dimensions) != 2L || any(dimensions <= 2L) ||
       dimensions[1L] != dimensions[2L])) {
    stop("regular_joint_d_gt2 requires PATE with equal used dimensions d > 2.",
         call. = FALSE)
  }
  if (covariance_scope == "critical_planar" &&
      (!pate || !identical(unname(dimensions), c(2L, 2L)))) {
    stop("critical_planar requires corrected PATE with both used dimensions exactly two.", call. = FALSE)
  }
  if (covariance_scope == "critical_planar") {
    critical_control <- .wm_cp_control(critical_control, fit, parameters)
  }
  scalar <- all(dimensions == 1L)
  if (!all(dimensions >= 2L) &&
      !(scalar && covariance_scope %in% c("full_x", "patt"))) {
    stop("Use all dimensions >= 2, or all scalar dimensions with full_x/PATT patt; mixed and other scalar scopes are unsupported.",
         call. = FALSE)
  }
  if (scalar && !pate) {
    declaration <- .wm_reciprocal_list(transport0, "transport0")
    if (!setequal(names(declaration), c("mode", "basis")) ||
        !identical(declaration$mode, "zero") ||
        !identical(declaration$basis, "current_centering")) {
      stop("Scalar PATT requires transport0 = list(mode = 'zero', basis = 'current_centering') and full-X predictions.",
           call. = FALSE)
    }
  }
  if (covariance_scope == "common_field" && dimensions[1L] != dimensions[2L]) {
    stop("The stated common_field branch requires equal used dimensions.", call. = FALSE)
  }
  if (covariance_scope == "full_x") {
    if (!is.null(transport0) || !is.null(transport1)) {
      stop("full_x supplies zero graph transport; omit transport specifications.", call. = FALSE)
    }
    specs <- list(potential0 = list(mode = "zero", basis = "full_x"),
                  potential1 = list(mode = "zero", basis = "full_x"))
  } else {
    specs <- list(potential0 = .wm_fi_transport_spec(transport0, fit$graph$scores0, n,
                                                    parameters, "transport0"))
    if (pate) specs$potential1 <- .wm_fi_transport_spec(transport1, fit$graph$scores1, n,
                                                      parameters, "transport1")
  }
  # This reviewed helper reconstructs the complete finite fit and rejects prior
  # adjustments. It never substitutes its validator's fixed-map interval.
  smooth <- .wm_wdsm_graph_gradient(fit, DW, D0, D1)
  W <- fit$weights
  Z <- fit$data$Z
  gamma_raw <- mean(if (pate) W else Z * W)
  if (!is.finite(gamma_raw) || gamma_raw <= 0) {
    stop("The raw target-weight mean is not positive and finite.", call. = FALSE)
  }
  .wm_reciprocal_agree(gamma_raw / fit$weight_scale, fit$gamma,
                       "Raw target-weight normalization")
  contract <- list(
    status = "user_supplied_population_assumptions_not_verified",
    target = paste("Iid observed rows; identified baseline PATE/PATT; original corrected",
      "self-normalized full-sample rule with fixed M, dimensions and row-priority ties."),
    fitted_law = paste("All used d_z>=2. The bounded own-field alternative requires bounded outcomes, residuals, predictions and",
      "joint influence; own-(baseline score,baseline weight) centering; compact regular",
      "baseline geometry, current upper donor-density control and the proved compact-index",
      "modulus/current-row laws. Baseline predictions equal the own-score population",
      "functions in the point theorem and satisfy its target identities. Predictions",
      "have bounded first two parameter derivatives.",
      "One joint regular root-n expansion includes scores, means and fitted scales;",
      "influence estimates are empirically L2 consistent. This is distributional transfer,",
      "not same-realization equality of fitted and baseline graphs.",
      "Alternatively, PATE full_x with correct full-X predictions, or PATT patt with",
      "the correct full-X control prediction and zero current-centering transport,",
      "permits unbounded outcomes and complete",
      "influences. Predictions and their first two parameter derivatives remain bounded.",
      "Each used direction retains compact convex full-dimensional donor support,",
      "positive bounded baseline donor density, bounded supported query density,",
      "bounded C2 parameter derivatives and uniform current donor-density upper control",
      "on every root-n parameter ball. Known positive bounded W=w_Z(X) is retained.",
      "For B=(X,Z), m(B)=E(ell|B) and some 0<delta<=2, E||m(B)||^(2+delta) is finite;",
      "the conditional (2+delta) moments of full-X residuals and ell-m(B) are uniformly bounded.",
      "The regular complete score/prediction/centering/scale stack, empirical L2 influence",
      "consistency and actual square/cross/maximum row laws retain all covariance blocks;",
      "the full nuisance covariance may be singular. These premises are not verified by",
      "finite arrays, a fitted-model label or numerical convergence. Other own-field and",
      "estimated-transport routes retain their bounded assumptions.",
      "A separate correct-full-X finite-moment current-graph route permits unbounded",
      "supplied weights and predictions under the applicable sampling theorem's",
      "joint oracle/nuisance law, integrated square/cross/maximum row laws, prediction derivative/Hessian envelopes,",
      "complete regular influence and empirical L2 consistency. The complete Gaussian",
      "route additionally requires independent conditional Gaussian innovations with",
      "positive bounded variances, an influence m(B)+A(B)epsilon with global second",
      "moments, and exact graph determination by all design rows and the complete",
      "Gaussian noise projection. Uniform current-graph moments and the actual",
      "prediction-slope limit are required; an influence approximation in rankings",
      "or finite solver convergence does not establish these conditions."),
    covariance = switch(covariance_scope,
      distinct_rarity = paste("Baseline maps (A,B0),(A,B1), after allowed fixed regular",
        "coordinate changes, have bounded full-dimensional joint density with",
        "dim(A)<max(d0,d1); the actual local-index reciprocal-rarity theorem applies."),
      full_x = paste("Baseline residuals are centered given full covariates and treatment;",
        "the fixed baseline weight is a covariate/treatment function. This implies exact",
        "current centering and zero graph transport, without a joint-map density condition."),
      common_field = paste("Both current maps generate the same score field for every",
        "deterministic parameter in a neighborhood, e.g. invertible transforms of one",
        "common map. The stated branch uses equal d>=2. Sample agreement is not proof."),
      patt = "Only control-direction conditions are needed; no PATE reciprocal branch is imposed."),
    transport = paste("Every zero mode is an explicit fixed-map/current-centering declaration.",
      "Estimated directions additionally require the effective raw signed chart, compact",
      "interior signed support and cutoff, integrably C3 ordinary marked densities,",
      "positive interior donor density, uniform generated-input control and a vanishing",
      "propagated quadrature error. Coordinate binding verifies only observed alignment."),
    derivatives = paste("NULL mean derivative arguments explicitly declare identically",
      "zero derivatives in the joint parameter; they are never inferred from the fitted",
      "prediction values or a successful model fit. All supplied prediction derivatives",
      "are complete and in the common parameter order."),
    smooth_weights = paste("Positive supplied weights define the target and are treated",
      "as known throughout score, prediction and scale estimation. NULL or exact zero",
      "weight_derivative is required. Each sampling branch retains its stated weight",
      "bounds or integrated moments; estimated-weight uncertainty is outside this scope."),
    inference = paste("Positive total limiting variance is required. All same-row C and",
      "Sigma cross-covariances are retained; no fixed-map reciprocal subtraction is used.",
      "No generic scalar maps, nonzero-reciprocal fitted branch, dependent-design inference,",
      "arbitrary learner, original-refit bootstrap or actual sampling-variance convergence claim."))
  if (scalar) {
    contract$fitted_law <- paste("Design-fitted alternative: all used matching dimensions equal one; iid rows with",
      "bounded Y and known positive bounded W=w_Z(X); correct full-X predictions",
      "with bounded first two parameter derivatives; one regular root-n joint",
      "score/prediction/scale influence, bounded and empirically L2 consistent.",
      "The actual graph parameter and its influence component are measurable in",
      "the full observed design B=(X,Z); prediction parameters may use outcomes.",
      "Each direction has a compact monotone scalar chart reconstructing X,",
      "a shared arm reference measure with joint arm densities bounded above/below",
      "and bounded first-coordinate density derivatives, and uniformly bounded C2",
      "maps with positive bounded coordinate derivative. Same direction maps apply",
      "to queries and donors. Root marks may be bounded Borel; charts may differ",
      "between directions. These premises are not verified from finite arrays.",
      "Alternatively, the actual scalar graph may use outcomes under the exact-root",
      "full-X sampling theorem. Its closed finite-dimensional graph estimating equation",
      "has a nonsingular population Jacobian, bounded parameter Jacobian/Hessian,",
      "finite 2+delta score moments and the",
      "accessible smooth continuous-coordinate block-score domination/submersion",
      "conditions. The selected root is exact, or graph-equivalent; ordinary numerical",
      "root error alone does not certify scalar graph transfer. Each direction has the",
      "finite-piecewise Lipschitz-moving density/mark chart, comparable query/donor",
      "densities and within-piece Wasserstein-1 Lipschitz bounded joint marks.",
      "This alternative retains known positive bounded W=w_Z(X), correct bounded",
      "full-X predictions and bounded C2 prediction derivatives, the complete regular",
      "joint influence and empirical L2 consistency. Residuals and full influences",
      "have finite 2+delta moments; centered innovations have uniformly bounded",
      "conditional 2+delta moments. The full prediction covariance may be singular.",
      "A separate noncompact design-fitted correct-full-X route permits unbounded",
      "supplied weights and predictions when its scalar current-graph square/cross/maximum",
      "and prediction-gradient laws, moment/cap-removal conditions and complete regular",
      "influence establish the applicable joint sampling law. The complete Gaussian route requires independent",
      "conditional Gaussian innovations with positive bounded variances, global second",
      "moments for the complete m(B)+A(B)epsilon influence, and exact graph determination",
      "by all design rows and the complete Gaussian noise projection. Its uniform",
      "integrated current-graph moments, square-integrable prediction derivative/Hessian",
      "envelopes, actual slope limit and empirical L2 influence consistency remain.",
      "These routes do not imply arbitrary scalar-map or misspecified-score geometry.",
      "At least one complete alternative must hold; none is verified from arrays.")
    contract$covariance <- paste("Under the applicable stated scalar sampling-law alternative,",
      "correct full-X means and known W=w_Z(X) give zero graph drift and zero",
      "reciprocal subtraction on the actual graph through the corresponding proved law.",
      "PATT requires only the control full-X mean. No joint density or independence",
      "of the two directional coordinates is imposed; all full nuisance cross blocks remain.")
    contract$transport <- paste("Scalar PATE full_x supplies zero full-X graph transport;",
      "scalar PATT requires explicit zero/current_centering with correct full-X",
      "prediction. No estimated scalar transport or mixed dimensions are admitted.")
    contract$inference <- paste("Positive total limiting variance is required. Under the stated scalar premises,",
      "empirical augmented-row variance consistency and row Lindeberg are derived.",
      "The complete finite-M prediction slope and all same-row C/Sigma blocks remain.",
      "These contracts do not certify population assumptions, an original numerical",
      "scalar graph, an actual full/partial refit callback, growing-B moments or raw",
      "sampling-variance convergence. Outcome-fitted maps outside the declared exact-root",
      "class, dependent designs and unrestricted learners remain unsupported.")
  }
  if (covariance_scope == "regular_joint_d_gt2") {
    contract$fitted_law <- paste("Equal fixed used d>2, original corrected PATE, bounded iid rows,",
      "known positive bounded W=w_Z(X), baseline target identities and own-score/weight",
      "residual centering. The complete F1-F4 regular root-n stack includes scores,",
      "predictions and scales, bounded joint influence, empirical L2 influence consistency,",
      "bounded smooth prediction derivatives, baseline/current marginal geometry,",
      "current residual projection/transport, compact-index moment/modulus and actual-fit",
      "laws. Full sensitivity estimation retains its separate additional conditions.",
      "This branch does not inherit the unbounded full-X or Gaussian alternatives.")
    contract$covariance <- paste("Finite compact embedded C2 baseline joint-score sheets with",
      "piecewise smooth Hausdorff-null boundaries, dimension r>=d and full-rank score",
      "projections. Retain regular marginal donor/query geometry and full donor",
      "subdensities h_{z|z} continuous off Lebesgue-null sets; bounded arm sheet",
      "densities, finite regular atlases with uniformly bounded coordinate Jacobians",
      "and local inverses, and uniform component small-ball bounds are required.",
      "Distinct d-sheets have Hausdorff-d-null intersections. Coincident",
      "positive-measure sheets are merged with their arm measures and root-mark",
      "mixtures only when a finite smooth decomposition is admissible; other",
      "positive-measure coincidences are outside this branch. The full marginal",
      "donor-weight kernel, including every sheet and preimage, retains its stated",
      "a.e. Wasserstein-1 continuity. Root-n movement is smaller than the neighbor",
      "scale; current joined sheets need not stay coincident. Finite arrays and",
      "scope labels do not verify these population conditions.")
    contract$inference <- paste("Positive complete limiting contrast variance is required.",
      "Signed actual-graph reciprocal covariance is subtracted from the diagonal base",
      "block, and full matching/nuisance C and Sigma terms remain. A nonpositive finite",
      "base block does not gate a positive full contrast. No hidden floor is applied.",
      "Common bootstrap dispatch uses explicit corrected-variance Gaussian marginal",
      "replication; ordinary row/count multipliers and original-refit/rematching validity",
      "are not established. Actual sampling-variance convergence needs separate UI/L2.")
  }
  if (covariance_scope == "critical_planar") {
    contract$fitted_law <- paste("Original corrected PATE, equal fixed d=2, iid bounded known-W own-field",
      "framework, complete F1-F4 and baseline S7 finite regular joint sheets.",
      "The complete regular root-n stack and same-sample compact-index laws retain all nuisance",
      "scores, prediction and scaling parameters, covariance and full transport/smooth slope.",
      "This branch does not inherit the unbounded full-X or Gaussian alternatives.")
    contract$covariance <- paste("Whole critical reciprocal function estimated by macro-scale",
      "radial pairing, guarded local rank-two PCA and full marginal donor-weight laws.",
      "Both endpoint tangents are actual matching-coordinate derivatives.",
      "Finite valid density/mass/projection/mark guard budgets, a deterministic positive-fraction",
      "split and h=n^-a, 0<a<1/2, are caller premises. Empty/rejected auxiliary summaries do",
      "not repair a failed original fit. Numerical error must vanish along the study sequence.")
    contract$inference <- paste("Generally non-Gaussian complete nuisance/Schur-mixture law.",
      "The twice-smooth finite-p stack, its derivative envelopes and bounded/fourth-moment",
      "influence rate justify the disclosed n^-1/3 covariance-support threshold; generic",
      "empirical-L2 consistency alone does not. The fixed nuisance block/range and raw signed",
      "kernel/Schur/cone diagnostics are retained. No Normal interval or analytic finite-B",
      "variance is asserted. Complete-input consistency and strict-CDF quantile premises",
      "remain; actual sampling-variance convergence requires separate UI.")
    if (identical(critical_control$direct_variance, "raw_only")) {
      contract$inference <- paste("Direct signed scalar variance only, using the complete",
        "centered empirical influence covariance and cross block already used by the common",
        "fitted-row helper. This differs finitely from an uncentered second moment.",
        "Generic empirical-L2 covariance consistency plus separate complete diagonal,",
        "normalizer, C, full-slope, compact-uniform whole-kernel and global-envelope",
        "consistency suffice for the raw variance limit. No covariance support rate, rank",
        "recovery, pseudoinverse, mixture replication, Normal pivot or original refit",
        "bootstrap is authorized. Actual sampling-variance convergence needs separate UI.")
    } else if (identical(critical_control$direct_variance, "with_mixture")) {
      contract$direct_variance <- paste("Optional raw-empirical and support-aligned signed",
        "direct Gaussian averages retain all donor/root marks and the complete nuisance",
        "blocks. These raw variances need not equal the finite operative cone-positive-part",
        "expectation, which remains unevaluated/NA; finite-B Monte Carlo summaries are",
        "separate. The existing mixture support and quantile premises remain unchanged.")
    }
  }
  if (any(vapply(specs, function(spec) {
    identical(spec[["representation", exact = TRUE]], "score_measurable_weight")
  }, logical(1)))) {
    contract$transport <- paste("Every zero mode is an explicit fixed-map/current-centering declaration.",
      "Conditional-weight-chart estimates retain the effective raw signed chart, ordinary",
      "marked-density and vanishing quadrature-error premises. Score-measurable-weight",
      "estimates instead require donor W=omega_z(S_z^0) with a positive bounded smooth",
      "omega_z, actual complete-p score tangents, three smooth score subdensities,",
      "uniform generated-score/tangent control and residual-level consistency.",
      "Each estimated direction retains own-score/weight centering and compact interior",
      "signed support with a justified cutoff. Score-weight estimates use denominator floors",
      "with zero active-floor derivatives and preserve the exact constant-donor-weight",
      "log-gradient reduction. Their zero quadrature bound excludes statistical estimation",
      "error and does not certify the sampling, covariance or bootstrap premises.",
      "Coordinate binding verifies only observed alignment; supplied individual W is unchanged.")
  }
  transports <- list()
  for (direction in names(specs)) {
    spec <- specs[[direction]]
    arm <- if (direction == "potential0") 0L else 1L
    if (spec$mode == "zero") {
      transports[[direction]] <- .wm_fi_zero_transport(parameters, spec$basis)
      next
    }
    mu <- if (arm == 0L) fit$predictions$mean0 else fit$predictions$mean1
    arguments <- c(list(raw_scores = spec$raw_scores, Z = Z, W = W,
      residual = fit$data$Y - mu, tangents = spec$tangents,
      cutoff = spec$cutoff, donor_arm = arm, M = fit$M), spec$controls)
    if (!is.null(spec[["representation", exact = TRUE]])) {
      arguments$representation <- spec$representation
    }
    value <- tryCatch(do.call(wm_graph_transport, arguments), error = identity)
    if (inherits(value, "error")) {
      return(structure(list(estimate = fit$estimate, n = n, M = fit$M,
        estimand = fit$estimand, fit = fit, status = "transport_unavailable",
        available = FALSE, numerically_available = FALSE,
        conditional_inference_available = FALSE,
        failure_stage = direction, unavailable_reason = conditionMessage(value),
        root_n_variance = NA_real_, variance = NA_real_, se = NA_real_,
        conf.int = stats::setNames(rep(NA_real_, 2L), c("lower", "upper")),
        conf.level = conf.level, parameter_names = parameters,
        dimensions = dimensions, covariance_scope = covariance_scope,
        smooth_derivative = smooth, graph_transport = transports,
        assumptions_verified = FALSE, application_verified = FALSE,
        contract = contract), class = c("wm_fitted_inference", "list")))
    }
    value$map_binding <- spec$map_binding
    value$mode <- "estimate"
    transports[[direction]] <- value
  }
  graph <- -transports$potential0$graph_drift / gamma_raw
  graph_error <- transports$potential0$quadrature_error_bound / gamma_raw
  if (pate) {
    graph <- graph + transports$potential1$graph_drift / gamma_raw
    graph_error <- graph_error + transports$potential1$quadrature_error_bound / gamma_raw
  }
  if (any(!is.finite(c(graph, graph_error)))) {
    stop("Normalized transport arithmetic exceeded numerical range.", call. = FALSE)
  }
  inference <- .wm_wdsm_fitted_variance(fit, influence, smooth$sensitivity,
    graph, covariance_scope, conf.level)
  inference$fit <- fit
  inference$M <- fit$M
  inference$dimensions <- dimensions
  inference$smooth_derivative <- smooth
  inference$graph_transport <- transports
  inference$graph_sensitivity_quadrature_error_bound <- graph_error
  inference$raw_target_weight_mean <- gamma_raw
  inference$numerically_available <- inference$available
  inference$conditional_inference_available <- inference$available
  inference$status <- if (inference$available) paste0("conditional_", covariance_scope) else
    "nonpositive_variance"
  inference$unavailable_reason <- if (inference$available) "" else
    "The assembled root-n or estimator-scale variance is not strictly positive."
  inference$application_verified <- FALSE
  inference$contract <- contract
  inference$original_refit_limit_agreement_declared <- FALSE
  if (covariance_scope == "critical_planar") {
    value <- tryCatch(.wm_cp_complete(inference, critical_control), error = identity)
    if (inherits(value, "error")) {
      inference$available <- inference$numerically_available <-
        inference$conditional_inference_available <- FALSE
      inference$status <- "critical_kernel_unavailable"
      inference$failure_stage <- "critical_planar_kernel"
      inference$unavailable_reason <- conditionMessage(value)
      inference$root_n_variance <- inference$variance <- inference$se <- NA_real_
      inference$conf.int <- stats::setNames(rep(NA_real_, 2L), c("lower", "upper"))
      inference$critical_control <- critical_control
      if (!is.null(critical_control$direct_variance)) {
        inference$direct_variance <- list(mode = critical_control$direct_variance,
          status = "critical_kernel_unavailable", available = FALSE,
          numerically_available = FALSE, unavailable_reason = conditionMessage(value),
          raw_empirical = list(raw_root_n_variance = NA_real_, raw_variance = NA_real_,
            available = FALSE, numerically_available = FALSE))
      }
    } else inference <- value
  }
  class(inference) <- c("wm_fitted_inference", "list")
  inference
}
