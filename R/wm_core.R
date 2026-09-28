# Exact fixed-M reference implementation. No nuisance models are fitted here.

.wm_numeric_vector <- function(x, n, name, positive = FALSE) {
  if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) ||
      length(x) != n || anyNA(x) || any(!is.finite(x))) {
    stop(name, " must be a finite numeric vector of length ", n, ".", call. = FALSE)
  }
  if (positive && any(x <= 0)) {
    stop(name, " must be strictly positive.", call. = FALSE)
  }
  as.numeric(x)
}

.wm_score_matrix <- function(x, n, name) {
  if (!is.matrix(x) || !is.numeric(x) || is.complex(x) ||
      nrow(x) != n || ncol(x) < 1L || anyNA(x) || any(!is.finite(x))) {
    stop(name, " must be a finite numeric matrix with ", n,
         " rows and at least one column.", call. = FALSE)
  }
  storage.mode(x) <- "double"
  unname(x)
}

.wm_cell_labels <- function(x, n, name) {
  if (is.null(x)) return(rep.int("all", n))
  if ((!is.atomic(x) && !is.factor(x)) || !is.null(dim(x)) ||
      length(x) != n || anyNA(x) ||
      (is.numeric(x) && any(!is.finite(x)))) {
    stop(name, " must contain one nonmissing cell label per row.", call. = FALSE)
  }
  x
}

.wm_sum_at <- function(index, value, n) {
  ans <- numeric(n)
  if (length(index)) {
    sums <- rowsum(value, index, reorder = FALSE)
    ans[as.integer(rownames(sums))] <- sums[, 1L]
  }
  ans
}

.wm_distances <- function(scores, query, donors) {
  delta <- sweep(scores[donors, , drop = FALSE], 2L, scores[query, ], "-")
  if (any(!is.finite(delta))) {
    stop("Score differences exceed numerical range; supply a fixed rescaling.",
         call. = FALSE)
  }
  # Columnwise pmax avoids one R function call per candidate donor while
  # preserving the scaled Euclidean norm used by the frozen reference.
  largest <- abs(delta[, 1L])
  if (ncol(delta) > 1L) {
    for (column in 2:ncol(delta)) largest <- pmax(largest, abs(delta[, column]))
  }
  distance <- largest * sqrt(rowSums((delta / ifelse(largest > 0, largest, 1))^2))
  if (any(!is.finite(distance))) {
    stop("Euclidean distances exceed numerical range; supply a fixed rescaling.",
         call. = FALSE)
  }
  distance
}

.wm_nearest_positions <- function(distance, donors, M) {
  # Select the exact Mth value, then resolve its boundary by original row.
  # Candidate donors are cached in increasing original-row order.
  if (M == length(distance)) {
    selected <- seq_along(distance)
  } else {
    cutoff <- sort.int(distance, partial = M)[M]
    closer <- which(distance < cutoff)
    boundary <- which(distance == cutoff)
    selected <- c(closer, boundary[seq_len(M - length(closer))])
  }
  selected[order(distance[selected], donors[selected], method = "radix")]
}

#' Weighted matching with supplied scores and nuisance predictions
#'
#' @details
#' For iid observed rows with law \eqn{Q}, fixed individual weights
#' \eqn{W=w_Z(X)} define the target full-data law \eqn{dP/dQ=W/E_Q(W)}.
#' PATE is the mean causal effect under \eqn{P}; PATT is its mean among
#' treated individuals under \eqn{P}. Weights may depend on treatment and
#' covariates, with weights depending only on covariates as a special case.
#' Consistency, full-covariate ignorability (one-sided control-outcome
#' ignorability suffices for PATT), weighted outcome integrability, treatment
#' overlap and reduced-mean identification are substantive assumptions.
#' Recovering an external population with sampling or response weights requires
#' its recovery assumptions. Balancing weights define a target through the same
#' law; the label alone does not establish identification or residual centering.
#'
#' The fixed-map inference theorem uses fixed M and fixed positive matching
#' dimensions, weights bounded above and away from zero, and positive arm
#' probabilities. In each direction the donor support is compact and convex
#' with nonempty interior; its density is bounded above and away from zero.
#' The query density is bounded and supported inside that donor support.
#' Required residual and weighted effect-heterogeneity moments have some order
#' \eqn{p>2}. The missing-arm means are Lipschitz, and corrected inference
#' requires a negligible root-n matching-bias remainder, consistent fitted
#' contribution rows and the stated nuisance rates. Positive limiting variance
#' is required for nondegenerate Wald inference. Finite input values and
#' successful numerical fitting do not verify these population conditions.
#' Full-sample stabilized PATE additionally requires the bounded-outcome and
#' common-score assumptions stated below; finite moments alone do not replace
#' them. Unknown fitted maps or weights require a separately justified
#' nuisance correction.
#'
#' For original PATE, residuals must be centered conditional on treatment,
#' both matching maps and the individual weight. For original PATT, control
#' residuals around the control reduced mean must be centered conditional on
#' treatment, the control matching map and weight. The treated-side
#' decomposition also requires residuals around its treated conditional mean
#' given that matching map to be centered given weight and the same field.
#' Only the control mean is fitted in the PATT calculation; this does not
#' remove the treated-side theorem condition. A correct full-covariate mean
#' model is one sufficient route to the respective centering conditions.
#'
#' The strata and evaluation-fold extensions use a fixed finite number of
#' cells with positive limiting cell and arm proportions and the preceding
#' conditions within each cell. Evaluation folds must be data-independent;
#' matching and any fold-specific inference use the same restricted graph.
#' Individual weights and matching strata do not specify a dependent
#' complex-survey variance estimator. Growing dimension, growing M and
#' unrestricted discrete ties are outside these inference claims.
#'
#' @export
wm_match <- function(Y, Z, weights, scores0, scores1 = NULL, M = 3L,
                     estimand = c("PATE", "PATT"),
                     method = c("self_normalized", "stabilized"),
                     mean0 = NULL, mean1 = NULL, rho0 = NULL, rho1 = NULL,
                     fold_id = NULL, strata = NULL, variance = TRUE,
                     nuisance_influence = NULL, sensitivity = NULL) {
  estimand <- match.arg(estimand)
  method <- match.arg(method)
  if (!is.numeric(Y) || is.complex(Y) || !is.null(dim(Y)) || length(Y) < 2L) {
    stop("Y must be a finite numeric vector with at least two rows.", call. = FALSE)
  }
  n <- length(Y)
  Y <- .wm_numeric_vector(Y, n, "Y")
  if ((!is.numeric(Z) && !is.logical(Z)) || is.complex(Z) || !is.null(dim(Z)) ||
      length(Z) != n || anyNA(Z) || any(!Z %in% c(0, 1)) ||
      length(unique(Z)) != 2L) {
    stop("Z must have length n, contain only 0 and 1, and include both arms.",
         call. = FALSE)
  }
  Z <- as.integer(Z)
  weights <- .wm_numeric_vector(weights, n, "weights", positive = TRUE)
  if (!is.numeric(M) || is.complex(M) || length(M) != 1L || is.na(M) || !is.finite(M) ||
      M < 1 || M != floor(M) || M > .Machine$integer.max) {
    stop("M must be one positive integer.", call. = FALSE)
  }
  M <- as.integer(M)
  if (!is.logical(variance) || length(variance) != 1L || is.na(variance)) {
    stop("variance must be TRUE or FALSE.", call. = FALSE)
  }
  scores0 <- .wm_score_matrix(scores0, n, "scores0")
  if (estimand == "PATT" &&
      (!is.null(scores1) || !is.null(mean1) || !is.null(rho1))) {
    stop("PATT uses only scores0, mean0 and rho0; omit arm-1 inputs.", call. = FALSE)
  }
  scores1 <- if (is.null(scores1)) scores0 else .wm_score_matrix(scores1, n, "scores1")
  if (method == "stabilized" && estimand == "PATE" &&
      !identical(scores0, scores1)) {
    stop("Stabilized PATE requires the same supplied score matrix in both arms.",
         call. = FALSE)
  }
  corrected <- !is.null(mean0)
  if (estimand == "PATE" && (is.null(mean0) != is.null(mean1))) {
    stop("Supply both mean0 and mean1 for PATE, or omit both for a raw point estimate.",
         call. = FALSE)
  }
  if (variance && !corrected) {
    stop("variance = TRUE requires the relevant supplied mean predictions.",
         call. = FALSE)
  }
  mean0 <- if (is.null(mean0)) numeric(n) else .wm_numeric_vector(mean0, n, "mean0")
  mean1 <- if (is.null(mean1)) numeric(n) else .wm_numeric_vector(mean1, n, "mean1")
  if (method == "self_normalized" && (!is.null(rho0) || !is.null(rho1))) {
    stop("rho predictions are used only by method = 'stabilized'.", call. = FALSE)
  }
  if (method == "stabilized") {
    if (is.null(rho0) || (estimand == "PATE" && is.null(rho1))) {
      stop("Stabilized matching requires rho predictions for every donor direction.",
           call. = FALSE)
    }
    rho0 <- .wm_numeric_vector(rho0, n, "rho0", positive = TRUE)
    if (estimand == "PATE") rho1 <- .wm_numeric_vector(rho1, n, "rho1", positive = TRUE)
  }
  has_influence <- !is.null(nuisance_influence) || !is.null(sensitivity)
  if (has_influence) {
    if (method != "self_normalized" || !variance) {
      stop("Nuisance influence correction requires self_normalized inference.", call. = FALSE)
    }
    if (is.null(nuisance_influence) || is.null(sensitivity)) {
      stop("Supply nuisance_influence and sensitivity together.", call. = FALSE)
    }
    nuisance_influence <- .wm_score_matrix(nuisance_influence, n, "nuisance_influence")
    sensitivity <- .wm_numeric_vector(sensitivity, ncol(nuisance_influence), "sensitivity")
  }

  # Common normalization changes neither the target nor any normalized contribution.
  weight_scale <- max(weights)
  w <- weights / weight_scale
  if (any(!is.finite(w)) || any(w <= 0)) {
    stop("Common weight normalization underflowed; the weight range is too wide.", call. = FALSE)
  }
  if (method == "stabilized") {
    rho0 <- rho0 / weight_scale
    if (estimand == "PATE") rho1 <- rho1 / weight_scale
    if (any(!is.finite(rho0)) || any(rho0 <= 0) ||
        (estimand == "PATE" && (any(!is.finite(rho1)) || any(rho1 <= 0)))) {
      stop("Normalized rho predictions must remain finite and strictly positive.", call. = FALSE)
    }
  }
  folds <- .wm_cell_labels(fold_id, n, "fold_id")
  restrictions <- .wm_cell_labels(strata, n, "strata")
  # Encode equality before formatting: raw labels can contain separators, and
  # formatting distinct floating-point labels can collapse their identities.
  fold_code <- match(folds, unique(folds))
  stratum_code <- match(restrictions, unique(restrictions))
  cell_key <- paste(fold_code, stratum_code, sep = ":")
  cell <- match(cell_key, unique(cell_key))
  query_rows <- if (estimand == "PATE") seq_len(n) else which(Z == 1L)
  neighbors <- vector("list", n)
  # Integer cell IDs preserve label equality; split retains original row order.
  donor_key <- 2 * (cell - 1) + Z + 1
  donor_cache <- split(seq_len(n), donor_key)
  edge_count <- as.double(length(query_rows)) * M
  for (k in seq_along(query_rows)) {
    i <- query_rows[k]
    arm <- 1L - Z[i]
    donors <- donor_cache[[as.character(2 * (cell[i] - 1) + arm + 1)]]
    if (length(donors) < M) {
      stop("Insufficient donors for query row ", i, ": arm ", arm,
           " in its fold/stratum cell has ", length(donors),
           " donors but M = ", M, ".", call. = FALSE)
    }
    score <- if (arm == 0L) scores0 else scores1
    distance <- .wm_distances(score, i, donors)
    selected <- .wm_nearest_positions(distance, donors, M)
    j <- donors[selected]
    neighbors[[i]] <- j
    share <- if (method == "self_normalized") w[j] / sum(w[j]) else rep(1 / M, M)
    rho <- if (arm == 0L) rho0 else rho1
    outcome_share <- if (method == "self_normalized") share else share * w[j] / rho[j]
    if (k == 1L) {
      # Delay allocation until the first query validates M and its distances.
      edge_query <- rep(query_rows, each = M)
      edge_donor <- integer(edge_count)
      edge_arm <- 1L - Z[edge_query]
      edge_distance <- edge_share <- edge_outcome_share <- numeric(edge_count)
      block_offset <- seq_len(M)
    }
    block <- (as.double(k) - 1) * M + block_offset
    edge_donor[block] <- j
    edge_distance[block] <- distance[selected]
    edge_share[block] <- share
    edge_outcome_share[block] <- outcome_share
  }
  edges <- data.frame(query = edge_query, donor = edge_donor, arm = edge_arm,
                      distance = edge_distance, share = edge_share,
                      outcome_share = edge_outcome_share)
  rownames(edges) <- NULL
  if (any(!is.finite(edges$outcome_share))) {
    stop("Donor outcome coefficients exceed numerical range.", call. = FALSE)
  }
  incoming <- matrix(0, n, 2L, dimnames = list(NULL, c("arm0", "arm1")))
  imputed <- matrix(NA_real_, n, 2L, dimnames = list(NULL, c("Y0", "Y1")))
  raw_imputed <- imputed
  imputed[cbind(seq_len(n), Z + 1L)] <- Y
  raw_imputed[cbind(seq_len(n), Z + 1L)] <- Y
  for (arm in if (estimand == "PATE") 0:1 else 0L) {
    take <- edges$arm == arm
    e <- edges[take, , drop = FALSE]
    incoming[, arm + 1L] <- .wm_sum_at(e$donor, w[e$query] * e$share, n)
    mu <- if (arm == 0L) mean0 else mean1
    fitted_missing <- .wm_sum_at(e$query, e$outcome_share * (Y[e$donor] - mu[e$donor]), n)
    raw_missing <- .wm_sum_at(e$query, e$outcome_share * Y[e$donor], n)
    q <- unique(e$query)
    imputed[q, arm + 1L] <- mu[q] + fitted_missing[q]
    raw_imputed[q, arm + 1L] <- raw_missing[q]
  }
  outer <- if (estimand == "PATE") w else w * Z
  denominator <- sum(outer)
  gamma <- denominator / n
  if (!is.finite(gamma) || gamma <= 0) {
    stop("The normalized outer-weight mean is outside numerical range.", call. = FALSE)
  }
  contrast <- if (estimand == "PATE") imputed[, 2L] - imputed[, 1L] else Y - imputed[, 1L]
  raw_contrast <- if (estimand == "PATE") raw_imputed[, 2L] - raw_imputed[, 1L] else Y - raw_imputed[, 1L]
  numerator <- sum(outer * contrast)
  estimate <- numerator / denominator
  raw_estimate <- sum(outer * raw_contrast) / denominator
  sign <- 2L * Z - 1L
  own_mean <- ifelse(Z == 1L, mean1, mean0)
  residual <- Y - own_mean
  own_load <- incoming[cbind(seq_len(n), Z + 1L)]
  unordered <- data.frame(i = integer(), j = integer(), i_to_j = logical(), j_to_i = logical())
  edge_numerator <- numeric()
  mean_incoming <- numeric(n)
  if (method == "self_normalized") {
    uncentered <- if (estimand == "PATE") {
      w * (mean1 - mean0) + sign * (w + own_load) * residual
    } else {
      Z * w * (Y - mean0) - (1L - Z) * own_load * (Y - mean0)
    }
    actual_numerator <- uncentered - outer * estimate
    row_numerator <- actual_numerator
  } else {
    if (estimand == "PATE") {
      own_rho <- ifelse(Z == 1L, rho1, rho0)
      eta <- w * residual / own_rho
      u <- w - own_rho
      uncentered <- w * (mean1 - mean0) + sign * (own_rho + own_load) * eta
      actual_numerator <- uncentered - w * estimate
      mean_incoming <- .wm_sum_at(edges$donor, own_rho[edges$query] / M, n)
      row_numerator <- w * (mean1 - mean0 - estimate) + sign * (own_rho + mean_incoming) * eta
      lo <- pmin(edges$query, edges$donor)
      hi <- pmax(edges$query, edges$donor)
      pair_key <- paste(lo, hi, sep = ":")
      unique_key <- unique(pair_key)
      first <- match(unique_key, pair_key)
      pair_index <- match(pair_key, unique_key)
      edge_numerator <- .wm_sum_at(pair_index, sign[edges$donor] * eta[edges$donor] * u[edges$query] / M,
                                   length(unique_key))
      unordered <- data.frame(i = lo[first], j = hi[first],
                              i_to_j = .wm_sum_at(pair_index, as.numeric(edges$query == lo), length(unique_key)) > 0,
                              j_to_i = .wm_sum_at(pair_index, as.numeric(edges$query == hi), length(unique_key)) > 0)
    } else {
      eta <- numeric(n)
      controls <- Z == 0L
      eta[controls] <- w[controls] * (Y[controls] - mean0[controls]) / rho0[controls]
      uncentered <- Z * w * (Y - mean0) - (1L - Z) * own_load * eta
      actual_numerator <- uncentered - outer * estimate
      row_numerator <- actual_numerator
    }
  }
  row <- row_numerator / gamma
  edge <- edge_numerator / gamma
  actual <- actual_numerator / gamma
  nuisance_component <- numeric(n)
  if (has_influence) {
    nuisance_component <- as.vector(nuisance_influence %*% sensitivity)
    nuisance_component <- nuisance_component - mean(nuisance_component)
    row <- row + nuisance_component
    row <- row - mean(row)
  }
  if (any(!is.finite(c(estimate, raw_estimate, row, edge, actual)))) {
    stop("Estimator or contribution arithmetic exceeded numerical range.", call. = FALSE)
  }
  root_n_variance <- if (variance) (sum(row^2) + sum(edge^2)) / n else NA_real_
  if (variance && !is.finite(root_n_variance)) {
    stop("Contribution variance exceeded numerical range.", call. = FALSE)
  }
  identity <- list(direct = numerator, reconstructed = sum(uncentered),
                   difference = numerator - sum(uncentered),
                   centered_sum = sum(actual_numerator),
                   row_edge_sum = sum(row_numerator) + sum(edge_numerator))
  contract <- if (method == "self_normalized") {
    "Original donor normalization: identification and strong graph-field residual centering; supplied means and any nuisance correction need their separate rate/influence contract."
  } else if (estimand == "PATE") {
    "Common-score stabilized PATE: weighted centering, bounded-outcome theorem and the supplied mean/rho nuisance-rate contract; row-plus-edge variance."
  } else {
    "Stabilized PATT: weighted centering, finite-moment theorem and the supplied mean/rho nuisance-rate contract."
  }
  structure(list(estimate = estimate, raw_estimate = raw_estimate,
                 correction = estimate - raw_estimate, n = n, M = M,
                 estimand = estimand, method = method,
                 root_n_variance = root_n_variance, variance = root_n_variance / n,
                 se = sqrt(root_n_variance / n), numerator = numerator,
                 denominator = denominator, gamma = gamma,
                 weights = weights, analysis_weights = w, weight_scale = weight_scale,
                 data = list(Y = Y, Z = Z),
                 imputed = imputed, raw_imputed = raw_imputed,
                 predictions = list(mean0 = mean0, mean1 = if (estimand == "PATE") mean1 else NULL,
                                    rho0 = rho0, rho1 = rho1),
                 loads = list(incoming = incoming, mean_incoming = mean_incoming),
                 contributions = list(row = row, edge = edge, actual = actual,
                                      nuisance = nuisance_component),
                 numerator_identity = identity,
                 graph = list(edges = edges, neighbors = neighbors, unordered = unordered,
                              cell = cell, fold_id = folds, strata = restrictions,
                              scores0 = scores0, scores1 = if (estimand == "PATE") scores1 else NULL,
                              transform = "Supplied coordinates; no automatic scaling",
                              tie_rule = "Exact distance, then original row index"),
                 info = list(inference_status = if (variance) "user_supplied_nuisance_contract" else "point_estimate_only",
                             contract = contract, corrected = corrected,
                             nuisance_correction = has_influence,
                             weight_units = "Denominator, incoming loads, rho predictions and numerator identities use weights / weight_scale.",
                             restriction = "fold_id and strata restrict donor graphs; strata is not a survey-design variance specification.")),
            class = c("wm_match", "list"))
}
