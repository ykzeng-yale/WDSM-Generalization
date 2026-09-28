# Qualified same-sample prediction and coordinate uncertainty.

.wm_prediction_imbalance <- function(fit, D0, D1) {
  n <- fit$n
  p <- ncol(D0)
  pate <- identical(fit$estimand, "PATE")
  Z <- as.integer(fit$data$Z)
  w <- fit$weights / max(fit$weights)
  denominator <- sum(if (pate) w else Z * w)
  edges <- fit$graph$edges
  b0 <- b1 <- numeric(p)
  queries <- if (pate) seq_len(n) else which(Z == 1L)
  for (i in queries) {
    donors <- edges$donor[edges$query == i]
    lambda <- w[donors] / sum(w[donors])
    D <- if (Z[i] == 1L) D0 else D1
    imbalance <- D[i, ] - colSums(D[donors, , drop = FALSE] * lambda)
    if (Z[i] == 1L) b0 <- b0 + w[i] * imbalance else
      b1 <- b1 + w[i] * imbalance
  }
  b0 <- b0 / denominator
  b1 <- b1 / denominator
  list(potential0 = b0, potential1 = b1, sensitivity = b1 - b0)
}

#' Account for a qualified parametric prediction and coordinate fit
#' @export
wm_prediction_adjust <- function(fit, nuisance_influence,
                                 mean_derivative0 = NULL,
                                 mean_derivative1 = NULL,
                                 changing_scores, weight_derivative = NULL) {
  if (!inherits(fit, "wm_match") || !identical(fit$method, "self_normalized")) {
    stop("wm_prediction_adjust requires a self_normalized wm_match or wm_fit object.",
         call. = FALSE)
  }
  if (!is.null(fit$prediction_adjustment) ||
      isTRUE(fit$info$nuisance_correction) || !is.null(fit$weight_adjustment)) {
    stop("The fit already contains a nuisance correction; start from the ",
         "uncorrected fit to avoid double counting.", call. = FALSE)
  }
  n <- .wm_fit_integer(fit$n, "fit$n", 2L)
  pate <- identical(fit$estimand, "PATE")
  if (!pate && !identical(fit$estimand, "PATT")) {
    stop("Unsupported fit estimand.", call. = FALSE)
  }
  directions <- if (pate) c("potential0", "potential1") else "potential0"
  if (missing(changing_scores) || !is.logical(changing_scores) ||
      !is.null(dim(changing_scores)) || length(changing_scores) != length(directions) ||
      anyNA(changing_scores)) {
    stop("Supply changing_scores as a logical vector of length ",
         length(directions), " in potential-0 then potential-1 order.",
         call. = FALSE)
  }
  if (!is.null(names(changing_scores)) &&
      !identical(names(changing_scores), directions)) {
    stop("Named changing_scores must use potential0, then potential1 for PATE, ",
         "or potential0 for PATT.", call. = FALSE)
  }
  names(changing_scores) <- directions
  if (!pate && !is.null(mean_derivative1)) {
    stop("PATT uses only mean_derivative0; mean_derivative1 must be NULL.",
         call. = FALSE)
  }
  score0 <- .wm_score_matrix(fit$graph$scores0, n, "retained scores0")
  dimensions <- c(potential0 = ncol(score0))
  if (pate) {
    score1 <- .wm_score_matrix(fit$graph$scores1, n, "retained scores1")
    dimensions <- c(dimensions, potential1 = ncol(score1))
  }
  if (any(changing_scores & dimensions == 1L)) {
    stop("Same-sample changing scalar coordinates are outside this helper's ",
         "proved scope; unchanged scalar maps are allowed.", call. = FALSE)
  }
  planar <- any(changing_scores & dimensions == 2L)
  higher <- any(changing_scores & dimensions > 2L)
  if (planar && !is.null(weight_derivative)) {
    stop("Changing planar coordinates require fixed known weights; omit ",
         "weight_derivative. Simultaneous planar weight estimation is unsupported.",
         call. = FALSE)
  }
  influence <- .wm_weight_matrix(nuisance_influence, n, "nuisance_influence")
  p <- ncol(influence)
  supplied <- list(nuisance_influence = influence,
                   mean_derivative0 = mean_derivative0,
                   mean_derivative1 = mean_derivative1,
                   weight_derivative = weight_derivative)
  named_order <- NULL
  for (name in names(supplied)) {
    if (is.null(supplied[[name]])) next
    supplied[[name]] <- .wm_weight_matrix(supplied[[name]], n, name, p)
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
  dW <- if (is.null(weight_derivative)) zero else supplied$weight_derivative
  # Reuse the independent finite-graph reconstruction and weight derivative.
  # Zero influence prevents a correction before the combined sensitivity exists.
  checked <- wm_weight_adjust(fit, dW, zero)
  G <- checked$weight_adjustment$sensitivity
  base <- checked$contributions$row
  imbalance <- .wm_prediction_imbalance(fit, D0, D1)
  b0 <- imbalance$potential0
  b1 <- imbalance$potential1
  prediction_sensitivity <- imbalance$sensitivity
  sensitivity <- as.vector(G) + prediction_sensitivity
  component <- as.vector(influence %*% sensitivity)
  component <- component - mean(component)
  corrected <- base + component
  corrected <- corrected - mean(corrected)
  root_variance <- mean(corrected^2)
  if (any(!is.finite(c(b0, b1, sensitivity, component, corrected, root_variance)))) {
    stop("Prediction-adjustment arithmetic exceeds numerical range.", call. = FALSE)
  }
  names(b0) <- names(b1) <- names(prediction_sensitivity) <-
    names(sensitivity) <- parameter_names
  branch <- rep("fixed_map", length(directions))
  branch[changing_scores & dimensions == 2L] <- "planar_distributional_transfer"
  branch[changing_scores & dimensions > 2L] <- "higher_dimensional_graph_transfer"
  names(branch) <- directions
  fit$prediction_adjustment <- list(
    sensitivity = sensitivity, prediction_sensitivity = prediction_sensitivity,
    weight_sensitivity = G, potential0_imbalance = b0,
    potential1_imbalance = if (pate) b1 else NULL,
    mean_derivative0 = D0, mean_derivative1 = if (pate) D1 else NULL,
    weight_derivative = weight_derivative, nuisance_influence = influence,
    parameter_names = parameter_names, changing_scores = changing_scores,
    dimensions = dimensions, branch = branch,
    baseline_root_n_variance = mean(base^2),
    cross_term = 2 * mean(base * component),
    nuisance_root_n_variance = mean(component^2),
    weight_derivative_check = if (is.null(weight_derivative)) NULL else
      checked$weight_adjustment[c("direct_sensitivity", "max_fraction_derivative_sum")],
    contract = list(
      status = "user_supplied_assumptions_not_verified_by_arrays",
      sampling = paste("Iid observed rows; fixed M, score and parameter dimensions;",
                       "identified target; fixed row-index tie rule; any fixed finite",
                       "strata/folds satisfy positive proportions and the conditions within cells."),
      predictions = paste("Correct full-X conditional means at the true parameter;",
                          "predictions, first two parameter derivatives and residuals",
                          "uniformly bounded on a parameter neighborhood.",
                          "Supplied derivative matrices include the complete chain rule."),
      geometry = paste("Baseline score supports compact and convex; donor densities",
                       "bounded above and away from zero; query densities bounded",
                       "and supported there; positive arm probabilities; weights",
                       "bounded above and away from zero. W=w_Z(X), with X enlarged",
                       "to contain any outcome-informative weight marks."),
      influence = paste("One same-sample joint root-n asymptotic linear expansion",
                        "for all fitted components; influence has a 2+delta moment",
                        "and the supplied estimates are empirically L2 consistent."),
      planar = if (planar) paste(
        "Known weights; uniform first-order score expansion with bounded derivatives,",
        "root-n local displacement and local-index Lipschitz bounds; uniform upper",
        "bound on each deterministic current donor density and the source planar",
        "centered-process modulus and marked-law assumptions. Any empirical scales",
        "obey the source bounded moment-rule and positive-variance conditions.",
        "This is distributional transfer, not fitted/baseline graph pathwise equivalence.") else NULL,
      higher_dimensional = if (higher) paste(
        "Score maps have uniformly bounded parameter derivatives; squared-distance",
        "comparisons are fixed-degree polynomials in the fixed-dimensional parameter,",
        "including any fitted metric, with ties resolved by fixed row order.",
        "Probability-scale logistic coordinates do not automatically meet this condition.") else NULL,
      weights = if (!is.null(weight_derivative)) paste(
        "Weights are smooth functions of the same joint parameter with uniformly",
        "bounded first two derivatives and remain positive and bounded nearby.",
        "Every changing matching direction has dimension greater than two.") else
          "Weights are fixed and known.",
      exclusions = paste("No changing scalar coordinates, stabilized rule, misspecified",
                         "full-X prediction, arbitrary learner, independent-training",
                         "rescaling, or unrestricted generated-score guarantee.")))
  fit$contributions$row <- corrected
  fit$contributions$nuisance <- component
  fit$root_n_variance <- root_variance
  fit$variance <- root_variance / n
  fit$se <- sqrt(fit$variance)
  fit$info$nuisance_correction <- TRUE
  fit$info$inference_status <- "bounded_full_x_prediction_influence_contract"
  fit
}
