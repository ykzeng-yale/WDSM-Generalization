# Independently trained parametric scores and predictions in any dimension.

#' Account for independent training of matching scores and predictions
#' @export
wm_training_adjust <- function(fit, training_influence,
                               mean_derivative0 = NULL,
                               mean_derivative1 = NULL) {
  if (!inherits(fit, "wm_match") || !identical(fit$method, "self_normalized")) {
    stop("wm_training_adjust requires a self_normalized wm_match or wm_fit object.",
         call. = FALSE)
  }
  if (isTRUE(fit$info$nuisance_correction) || !is.null(fit$weight_adjustment) ||
      !is.null(fit$prediction_adjustment) || !is.null(fit$training_adjustment) ||
      length(fit$contributions$training) > 0L) {
    stop("The fit already contains a nuisance correction; start from the ",
         "uncorrected fit to avoid double counting.", call. = FALSE)
  }
  n <- .wm_fit_integer(fit$n, "fit$n", 2L)
  pate <- identical(fit$estimand, "PATE")
  if (!pate && !identical(fit$estimand, "PATT")) {
    stop("Unsupported fit estimand.", call. = FALSE)
  }
  if (!pate && !is.null(mean_derivative1)) {
    stop("PATT uses only mean_derivative0; mean_derivative1 must be NULL.",
         call. = FALSE)
  }
  if (!is.matrix(training_influence) || nrow(training_influence) < 2L) {
    stop("training_influence must be a finite numeric m-by-p matrix with m >= 2.",
         call. = FALSE)
  }
  m <- nrow(training_influence)
  influence <- .wm_weight_matrix(training_influence, m, "training_influence")
  p <- ncol(influence)
  supplied <- list(training_influence = influence,
                   mean_derivative0 = mean_derivative0,
                   mean_derivative1 = mean_derivative1)
  named_order <- NULL
  for (name in names(supplied)) {
    if (is.null(supplied[[name]])) next
    rows <- if (name == "training_influence") m else n
    supplied[[name]] <- .wm_weight_matrix(supplied[[name]], rows, name, p)
    order <- colnames(supplied[[name]])
    if (!is.null(order)) {
      if (anyNA(order) || any(!nzchar(order)) || anyDuplicated(order)) {
        stop("Parameter column names must be nonempty and unique.", call. = FALSE)
      }
      if (!is.null(named_order) && !identical(order, named_order)) {
        stop("All named derivative and influence matrices must have the same ",
             "parameter order.", call. = FALSE)
      }
      named_order <- order
    }
  }
  parameter_names <- if (is.null(named_order)) paste0("parameter", seq_len(p)) else
    named_order
  zero <- matrix(0, nrow = n, ncol = p,
                 dimnames = list(NULL, parameter_names))
  D0 <- if (is.null(mean_derivative0)) zero else supplied$mean_derivative0
  D1 <- if (is.null(mean_derivative1)) zero else supplied$mean_derivative1
  # A zero weight derivative invokes the clean-data, graph, mean and row
  # reconstruction checks without treating training rows as evaluation rows.
  wm_weight_adjust(fit, zero, zero)
  if (length(fit$contributions$edge) != 0L) {
    stop("A self-normalized fit must not contain edge multiplier contributions.",
         call. = FALSE)
  }
  score0 <- .wm_score_matrix(fit$graph$scores0, n, "retained scores0")
  dimensions <- c(potential0 = ncol(score0))
  if (pate) {
    score1 <- .wm_score_matrix(fit$graph$scores1, n, "retained scores1")
    dimensions <- c(dimensions, potential1 = ncol(score1))
  }
  imbalance <- .wm_prediction_imbalance(fit, D0, D1)
  sensitivity <- imbalance$sensitivity
  training_projection <- as.vector(influence %*% sensitivity)
  training_projection <- training_projection - mean(training_projection)
  ratio <- n / m
  training_contribution <- ratio * training_projection
  baseline_variance <- mean(fit$contributions$row^2)
  training_variance <- sum(training_contribution^2) / n
  root_variance <- baseline_variance + training_variance
  if (any(!is.finite(c(sensitivity, training_projection, training_contribution,
                       baseline_variance, training_variance, root_variance)))) {
    stop("Training-adjustment arithmetic exceeds numerical range.", call. = FALSE)
  }
  names(sensitivity) <- names(imbalance$potential0) <-
    names(imbalance$potential1) <- parameter_names
  fit$training_adjustment <- list(
    sensitivity = sensitivity, potential0_imbalance = imbalance$potential0,
    potential1_imbalance = if (pate) imbalance$potential1 else NULL,
    mean_derivative0 = D0, mean_derivative1 = if (pate) D1 else NULL,
    training_influence = influence, centered_training_projection = training_projection,
    parameter_names = parameter_names, evaluation_n = n, training_n = m,
    evaluation_training_ratio = ratio, dimensions = dimensions,
    baseline_root_n_variance = baseline_variance,
    training_root_n_variance = training_variance,
    contract = list(
      status = "user_supplied_independent_training_assumptions_not_verified_by_arrays",
      independence = paste("The m training observations and n evaluation observations",
                           "are independent samples. All fitted score/prediction components",
                           "represented here are trained solely on the independent sample.",
                           "Ordinary cross-fitting or supplying different matrix rows",
                           "does not establish this independence."),
      sampling = paste("Iid evaluation rows; fixed M, score and parameter dimensions;",
                       "identified target; fixed row-index ties; m and n tend to infinity",
                       "with n/m converging to a finite nonnegative constant.",
                       "Any fixed finite restriction cells satisfy the same conditions",
                       "and have positive limiting proportions."),
      weights = paste("Fixed known weights W=w_Z(X), bounded above and away from zero;",
                      "X includes all outcome-informative weight marks; positive arm probabilities."),
      predictions = paste("Correct full-X conditional means at the true parameter;",
                          "predictions, residuals and first two parameter derivatives",
                          "uniformly bounded nearby; complete prediction chain-rule derivatives."),
      geometry = paste("For every nearby deterministic parameter, donor scores have",
                       "densities bounded above and away from zero on a fixed compact convex",
                       "body, and query densities are uniformly bounded and supported there.",
                       "Score maps have smooth bounded parameter derivatives."),
      marked_law = paste("For a convergence-determining family of bounded-Lipschitz tests",
                         "of the joint weight, prediction-derivative, mean and residual-variance",
                         "marks, the marked intensity densities converge uniformly on compact",
                         "interior sets as the parameter approaches truth and are locally",
                         "equicontinuous. Mark bounds give uniform integrability.",
                         "A smooth support-preserving conditional chart is one sufficient",
                         "construction; mere pointwise kernel continuity is insufficient."),
      influence = paste("Iid training rows with mean-zero finite-variance influence;",
                        "a joint training expansion sqrt(m)(theta_hat-theta0)=",
                        "m^(-1/2) sum(U_train)+o_p(1), consistent centered training",
                        "influence covariance and empirical L2-consistent influence estimates."),
      exclusions = paste("No estimated-weight correction, stabilized rule, same-sample",
                         "or overlapping cross-fit influence adjustment, misspecified full-X",
                         "means, or unrestricted fitted-scalar theorem is asserted.")))
  fit$contributions$training <- training_contribution
  fit$root_n_variance <- root_variance
  fit$variance <- root_variance / n
  fit$se <- sqrt(fit$variance)
  fit$info$nuisance_correction <- TRUE
  fit$info$inference_status <- "independent_training_prediction_influence_contract"
  fit
}
