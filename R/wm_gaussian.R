# Actual fitted-coordinate inference in the qualified Gaussian linear branch.
# The score function is evaluated at the unrounded same-sample OLS coefficient.

#' Fit weighted matching with a Gaussian linear prediction/score model
#' @export
wm_gaussian_match <- function(Y, Z, weights, design0, design1 = design0,
                              score_map, offset0 = NULL, offset1 = NULL,
                              M = 3L, estimand = c("PATE", "PATT"),
                              fold_id = NULL, strata = NULL, qr_tol = 1e-10) {
  estimand <- match.arg(estimand)
  if (!is.numeric(Y) || is.complex(Y) || !is.null(dim(Y)) || length(Y) < 2L)
    stop("Y must be a finite numeric vector with at least two rows.", call. = FALSE)
  n <- length(Y)
  Y <- .wm_numeric_vector(Y, n, "Y")
  .wm_potential_arm(Z, 0L, n)
  Z <- as.integer(Z)
  weights <- .wm_numeric_vector(weights, n, "weights", positive = TRUE)
  D0 <- .wm_weight_matrix(design0, n, "design0")
  D1 <- .wm_weight_matrix(design1, n, "design1", ncol(D0))
  p <- ncol(D0)
  if (p >= n) stop("The joint outcome design must have fewer columns than observations.", call. = FALSE)
  if (!is.null(colnames(D0)) && !is.null(colnames(D1)) &&
      !identical(colnames(D0), colnames(D1)))
    stop("Potential-outcome designs must use the same parameter-column order.", call. = FALSE)
  labels <- if (!is.null(colnames(D0))) colnames(D0) else colnames(D1)
  if (is.null(labels)) labels <- paste0("parameter", seq_len(p))
  if (anyNA(labels) || any(!nzchar(labels)) || anyDuplicated(labels))
    stop("Parameter names must be nonempty and unique.", call. = FALSE)
  colnames(D0) <- colnames(D1) <- labels
  offset0 <- if (is.null(offset0)) numeric(n) else .wm_numeric_vector(offset0, n, "offset0")
  offset1 <- if (is.null(offset1)) numeric(n) else .wm_numeric_vector(offset1, n, "offset1")
  if (!is.function(score_map)) stop("score_map must be a function of the fitted parameter vector.", call. = FALSE)
  if (!is.numeric(qr_tol) || is.complex(qr_tol) || length(qr_tol) != 1L ||
      !is.finite(qr_tol) || qr_tol <= 0 || qr_tol >= 1)
    stop("qr_tol must lie strictly between zero and one.", call. = FALSE)
  design <- D0
  design[Z == 1L, ] <- D1[Z == 1L, , drop = FALSE]
  observed_offset <- ifelse(Z == 1L, offset1, offset0)
  model <- qr(design, tol = qr_tol, LAPACK = FALSE)
  if (model$rank != p)
    stop("The joint outcome design is rank deficient; no columns are dropped.", call. = FALSE)
  parameter <- as.vector(qr.coef(model, Y - observed_offset))
  names(parameter) <- labels
  mean0 <- as.vector(offset0 + D0 %*% parameter)
  mean1 <- as.vector(offset1 + D1 %*% parameter)
  residual <- Y - ifelse(Z == 1L, mean1, mean0)
  gram <- crossprod(design) / n
  factor <- tryCatch(chol(gram), error = function(e) NULL)
  if (is.null(factor))
    stop("The outcome Gram matrix cannot be factored; no regularization is used.", call. = FALSE)
  influence <- t(backsolve(factor, forwardsolve(t(factor), t(design)))) * residual
  colnames(influence) <- labels
  if (any(!is.finite(c(parameter, mean0, mean1, residual, influence))))
    stop("Gaussian-model arithmetic exceeds numerical range.", call. = FALSE)
  maps <- score_map(parameter)
  if (is.matrix(maps)) maps <- list(scores0 = maps, scores1 = NULL)
  if (!is.list(maps) || !identical(sort(names(maps)), c("scores0", "scores1")))
    stop("score_map must return a numeric matrix or a list with scores0 and scores1.", call. = FALSE)
  S0 <- .wm_score_matrix(maps$scores0, n, "score_map scores0")
  S1 <- if (is.null(maps$scores1)) S0 else .wm_score_matrix(maps$scores1, n, "score_map scores1")
  fit <- wm_match(Y, Z, weights, S0,
    scores1 = if (estimand == "PATE") S1 else NULL, M = M,
    estimand = estimand, method = "self_normalized", mean0 = mean0,
    mean1 = if (estimand == "PATE") mean1 else NULL,
    fold_id = fold_id, strata = strata, variance = TRUE)
  imbalance <- .wm_prediction_imbalance(fit, D0, D1)
  b <- imbalance$sensitivity
  names(b) <- labels
  base <- fit$contributions$row
  component <- as.vector(influence %*% b)
  component <- component - mean(component)
  row <- base + component
  row <- row - mean(row)
  root_variance <- mean(row^2)
  if (any(!is.finite(c(b, component, row, root_variance))))
    stop("Gaussian prediction adjustment exceeds numerical range.", call. = FALSE)
  fit$contributions$row <- row
  fit$contributions$nuisance <- component
  fit$root_n_variance <- root_variance
  fit$variance <- root_variance / n
  fit$se <- sqrt(fit$variance)
  fit$gaussian_adjustment <- list(
    parameters = parameter, design0 = D0, design1 = D1,
    offset0 = offset0, offset1 = offset1, gram = gram,
    residual = residual, nuisance_influence = influence,
    sensitivity = b, potential0_imbalance = imbalance$potential0,
    potential1_imbalance = if (estimand == "PATE") imbalance$potential1 else NULL,
    baseline_root_n_variance = mean(base^2), cross_term = 2 * mean(base * component),
    nuisance_root_n_variance = mean(component^2),
    estimating_equation = as.vector(crossprod(design, residual) / n),
    dimensions = if (estimand == "PATE") c(potential0 = ncol(S0), potential1 = ncol(S1)) else
      c(potential0 = ncol(S0)),
    contract = list(
      status = "model_and_geometry_assumptions_not_verified_by_data",
      outcome_model = paste("Correct full-X finite linear conditional-mean family",
        "offset_z(X) + design_z(X)'theta; conditionally independent Gaussian errors",
        "given all (X,Z), with conditional variances bounded above and away from zero.",
        "Known offsets, design vectors and true means are bounded; population Gram is positive definite."),
      fitting = paste("The same-sample unweighted OLS coefficient is computed exactly",
        "up to floating-point arithmetic; score_map is evaluated at this unrounded coefficient.",
        "Conditional heteroskedasticity is allowed; observation weights enter the matching target, not OLS."),
      score_map = paste("A fixed-dimensional function of design covariates and theta only;",
        "no dependence on outcomes except through the supplied theta. Uniform parameter Lipschitz",
        "continuity and regular current donor/query densities on compact convex supports are required.",
        "The current marked local-law expectations must converge along every deterministic theta_n to theta0."),
      sampling = paste("Iid observed rows, fixed M and parameter/score dimensions, fixed known",
        "weights W=w_Z(X) bounded above and away from zero, identified causal target.",
        "X must include every weight mark that can inform outcomes; the full-X Gaussian model",
        "and design assumptions apply to this enlarged X. Any fixed finite",
        "matching strata/folds require positive cell/arm proportions and the geometry in each cell."),
      exclusions = paste("No generic Gaussianity test or causal-identification check is performed.",
        "Arbitrary non-Gaussian outcome fitting, estimated individual weights, misspecified full-X",
        "predictions, data-adaptive score families and dependent survey designs are outside this contract.")))
  fit$info$nuisance_correction <- TRUE
  fit$info$inference_status <- "conditional_gaussian_linear_score_contract"
  fit
}
