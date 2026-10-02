# Fixed-map reciprocal covariance inference. This does not alter a matching fit.

.wm_reciprocal_list <- function(x, label) {
  if (!is.list(x) || is.null(names(x)) || anyNA(names(x)) ||
      any(!nzchar(names(x))) || anyDuplicated(names(x))) {
    stop(label, " must be a named list with unique fields.", call. = FALSE)
  }
  x
}

.wm_reciprocal_agree <- function(x, expected, label) {
  if (!is.numeric(x) || is.complex(x) || length(x) != length(expected) ||
      anyNA(x) || any(!is.finite(x)) || any(!is.finite(expected)) ||
      any(abs(x - expected) > 1e-10 * pmax(1, abs(x), abs(expected)))) {
    stop(label, " is inconsistent with the stored graph, weights or predictions.",
         call. = FALSE)
  }
  invisible(TRUE)
}

# Validate the finite graph and point without computing an inference interval.
.wm_reciprocal_validated_state <- function(object, conf.level = 0.95,
                                             variance_floor = NULL) {
  object <- .wm_reciprocal_list(object, "object")
  if (!inherits(object, "wm_match") ||
      any(!class(object) %in% c("wm_match", "wm_fit", "list"))) {
    stop("object must be an unadjusted wm_match or wm_fit result.", call. = FALSE)
  }
  get <- function(x, name) x[[name, exact = TRUE]]
  n <- .wm_numeric_vector(get(object, "n"), 1L, "object n")
  M <- .wm_numeric_vector(get(object, "M"), 1L, "object M")
  if (n < 2 || n != floor(n) || n > .Machine$integer.max ||
      M < 1 || M != floor(M) || M > n) {
    stop("object n and M must be valid integer counts.", call. = FALSE)
  }
  n <- as.integer(n)
  M <- as.integer(M)
  estimand <- get(object, "estimand")
  method <- get(object, "method")
  if (!identical(estimand, "PATE") && !identical(estimand, "PATT")) {
    stop("Only PATE and PATT matching fits are supported.", call. = FALSE)
  }
  if (!identical(method, "self_normalized") && !identical(method, "stabilized")) {
    stop("Unsupported matching method.", call. = FALSE)
  }
  pate <- identical(estimand, "PATE")
  stabilized <- identical(method, "stabilized")
  conf.level <- .wm_numeric_vector(conf.level, 1L, "conf.level")
  if (conf.level <= 0 || conf.level >= 1) {
    stop("conf.level must lie strictly between zero and one.", call. = FALSE)
  }
  if (!is.null(variance_floor)) {
    variance_floor <- .wm_numeric_vector(variance_floor, 1L, "variance_floor", TRUE)
  }
  info <- .wm_reciprocal_list(get(object, "info"), "object info")
  status <- get(info, "inference_status")
  if (!identical(get(info, "corrected"), TRUE) ||
      !identical(get(info, "nuisance_correction"), FALSE) ||
      !is.character(status) || length(status) != 1L || is.na(status) ||
      !status %in% c("point_estimate_only", "user_supplied_nuisance_contract")) {
    stop("A corrected fit without nuisance or first-stage adjustments is required.",
         call. = FALSE)
  }
  adjustment_names <- grep("_adjustment$", names(object), value = TRUE)
  if (any(vapply(adjustment_names, function(key) !is.null(get(object, key)), logical(1)))) {
    stop("Adjusted fits cannot use fixed-map reciprocal inference.", call. = FALSE)
  }
  contributions <- .wm_reciprocal_list(get(object, "contributions"), "object contributions")
  if (!is.null(get(contributions, "training"))) {
    stop("Training contributions are outside this inference contract.", call. = FALSE)
  }
  nuisance <- .wm_numeric_vector(get(contributions, "nuisance"), n, "nuisance contributions")
  if (any(nuisance != 0)) stop("Nuisance-adjusted contributions are unsupported.", call. = FALSE)
  actual <- .wm_numeric_vector(get(contributions, "actual"), n, "actual contributions")
  dat <- .wm_reciprocal_list(get(object, "data"), "object data")
  Y <- .wm_numeric_vector(get(dat, "Y"), n, "object outcomes")
  Z <- .wm_numeric_vector(get(dat, "Z"), n, "object treatment")
  if (any(!Z %in% c(0, 1)) || length(unique(Z)) != 2L) {
    stop("The stored data must include both binary treatment arms.", call. = FALSE)
  }
  Z <- as.integer(Z)
  w <- .wm_numeric_vector(get(object, "analysis_weights"), n, "analysis weights", TRUE)
  weights <- .wm_numeric_vector(get(object, "weights"), n, "original weights", TRUE)
  weight_scale <- .wm_numeric_vector(get(object, "weight_scale"), 1L, "weight scale", TRUE)
  .wm_reciprocal_agree(weight_scale, max(weights), "Weight scale")
  .wm_reciprocal_agree(w, weights / weight_scale, "Analysis weights")
  prediction <- .wm_reciprocal_list(get(object, "predictions"), "object predictions")
  mu0 <- .wm_numeric_vector(get(prediction, "mean0"), n, "mean0")
  mu1 <- if (pate) .wm_numeric_vector(get(prediction, "mean1"), n, "mean1") else numeric(n)
  if (!pate && (!is.null(get(prediction, "mean1")) || !is.null(get(prediction, "rho1")))) {
    stop("PATT must not contain arm-1 predictions.", call. = FALSE)
  }
  rho0 <- rho1 <- NULL
  if (stabilized) {
    rho0 <- .wm_numeric_vector(get(prediction, "rho0"), n, "normalized rho0", TRUE)
    if (pate) rho1 <- .wm_numeric_vector(get(prediction, "rho1"), n, "normalized rho1", TRUE)
  } else if (!is.null(get(prediction, "rho0")) || !is.null(get(prediction, "rho1"))) {
    stop("Self-normalized fits must not contain rho predictions.", call. = FALSE)
  }
  graph <- .wm_reciprocal_list(get(object, "graph"), "object graph")
  for (label in c("cell", "fold_id", "strata")) {
    value <- get(graph, label)
    if (is.null(value) || length(value) != n || !is.null(dim(value)) ||
        !is.atomic(value) || anyNA(value) || length(unique(value)) != 1L ||
        (is.numeric(value) && any(!is.finite(value)))) {
      stop("Only unrestricted full-sample graphs are supported; invalid ", label, ".",
           call. = FALSE)
    }
  }
  S0 <- .wm_score_matrix(get(graph, "scores0"), n, "graph scores0")
  S1 <- if (pate) .wm_score_matrix(get(graph, "scores1"), n, "graph scores1") else NULL
  if (!pate && !is.null(get(graph, "scores1"))) stop("PATT has no arm-1 graph.", call. = FALSE)
  if (!identical(get(graph, "tie_rule"), "Exact distance, then original row index") ||
      !identical(get(graph, "transform"), "Supplied coordinates; no automatic scaling")) {
    stop("Unsupported graph metric or tie convention.", call. = FALSE)
  }
  edges <- get(graph, "edges")
  columns <- c("query", "donor", "arm", "distance", "share", "outcome_share")
  queries <- if (pate) seq_len(n) else which(Z == 1L)
  edge_count <- as.double(length(queries)) * M
  if (!is.data.frame(edges) || !identical(names(edges), columns) || nrow(edges) != edge_count) {
    stop("The graph must contain exactly M ordered edges per eligible query.", call. = FALSE)
  }
  for (column in columns) {
    .wm_numeric_vector(get(edges, column), edge_count, paste("edge", column))
  }
  q <- edges[["query"]]
  j <- edges[["donor"]]
  if (any(q != rep(queries, each = M)) || any(j != floor(j)) ||
      any(j < 1 | j > n) || any(edges[["arm"]] != 1L - Z[q]) ||
      any(Z[j] == Z[q]) || any(edges[["distance"]] < 0) ||
      any(edges[["share"]] <= 0) || any(edges[["outcome_share"]] <= 0)) {
    stop("Invalid graph indices, treatment arms, distances or donor fractions.", call. = FALSE)
  }
  neighbors <- get(graph, "neighbors")
  if (!is.list(neighbors) || length(neighbors) != n) {
    stop("The graph must retain one neighbor-list entry per observation.", call. = FALSE)
  }
  donor_total <- numeric(n)
  query_position <- match(seq_len(n), queries)
  expected_share <- expected_outcome_share <- numeric(edge_count)
  for (i in seq_len(n)) {
    block <- if (pate) (as.double(i) - 1) * M + seq_len(M) else {
      position <- query_position[i]
      if (is.na(position)) integer() else (as.double(position) - 1) * M + seq_len(M)
    }
    if (!length(block)) {
      if (!is.null(neighbors[[i]])) stop("A nonquery row has stored neighbors.", call. = FALSE)
      next
    }
    donors <- j[block]
    stored_neighbors <- neighbors[[i]]
    if (!is.numeric(stored_neighbors) || is.complex(stored_neighbors) ||
        !is.null(dim(stored_neighbors)) || length(stored_neighbors) != M ||
        anyNA(stored_neighbors) || any(stored_neighbors != donors) || anyDuplicated(donors)) {
      stop("Neighbor lists disagree with the directed edges.", call. = FALSE)
    }
    score <- if (Z[i] == 1L) S0 else S1
    distance <- .wm_distances(score, i, donors)
    .wm_reciprocal_agree(edges[["distance"]][block], distance, "Selected-edge distances")
    if (!identical(order(distance, donors, method = "radix"), seq_len(M))) {
      stop("Selected donors are not in stored distance/row order.", call. = FALSE)
    }
    donor_total[i] <- sum(w[donors])
    expected_share[block] <- if (stabilized) rep(1 / M, M) else w[donors] / donor_total[i]
    expected_outcome_share[block] <- if (stabilized) {
      rho <- if (Z[i] == 1L) rho0 else rho1
      expected_share[block] * w[donors] / rho[donors]
    } else expected_share[block]
  }
  .wm_reciprocal_agree(edges[["share"]], expected_share, "Donor fractions")
  .wm_reciprocal_agree(edges[["outcome_share"]], expected_outcome_share, "Donor outcome coefficients")
  incoming <- matrix(0, n, 2L)
  for (arm in if (pate) 0:1 else 0L) {
    take <- edges[["arm"]] == arm
    incoming[, arm + 1L] <- .wm_sum_at(j[take], w[q[take]] * expected_share[take], n)
  }
  loads <- .wm_reciprocal_list(get(object, "loads"), "object loads")
  stored_incoming <- get(loads, "incoming")
  if (!is.matrix(stored_incoming) || !identical(dim(stored_incoming), c(n, 2L))) {
    stop("Incoming loads must be an n-by-2 matrix.", call. = FALSE)
  }
  .wm_reciprocal_agree(stored_incoming, incoming, "Incoming loads")
  K <- incoming[cbind(seq_len(n), Z + 1L)]
  own_mean <- ifelse(Z == 1L, mu1, mu0)
  residual <- Y - own_mean
  own_rho <- if (stabilized) {
    if (pate) ifelse(Z == 1L, rho1, rho0) else rho0
  } else rep(1, n)
  r <- residual
  if (stabilized) {
    r <- numeric(n)
    used_donors <- if (pate) seq_len(n) else which(Z == 0L)
    r[used_donors] <- w[used_donors] * residual[used_donors] / own_rho[used_donors]
  }
  outer <- if (pate) w else Z * w
  gamma <- mean(outer)
  estimate <- .wm_numeric_vector(get(object, "estimate"), 1L, "point estimate")
  uncentered <- if (pate) {
    w * (mu1 - mu0) + (2L * Z - 1L) *
      (if (stabilized) own_rho + K else w + K) * r
  } else Z * w * (Y - mu0) - (1L - Z) * K * r
  expected_actual <- (uncentered - outer * estimate) / gamma
  .wm_reciprocal_agree(get(object, "gamma"), gamma, "Target denominator mean")
  .wm_reciprocal_agree(get(object, "denominator"), sum(outer), "Target denominator")
  .wm_reciprocal_agree(actual, expected_actual, "Actual contribution rows")
  # Reject stale decomposition arrays as well, even though inference uses actual.
  mean_incoming <- if (stabilized && pate) {
    .wm_sum_at(j, own_rho[q] / M, n)
  } else numeric(n)
  .wm_reciprocal_agree(get(loads, "mean_incoming"), mean_incoming,
                       "Mean-weight incoming loads")
  stored_row <- .wm_numeric_vector(get(contributions, "row"), n, "row contributions")
  expected_row <- if (stabilized && pate) {
    (w * (mu1 - mu0 - estimate) +
       (2L * Z - 1L) * (own_rho + mean_incoming) * r) / gamma
  } else expected_actual
  .wm_reciprocal_agree(stored_row, expected_row, "Decomposition rows")
  expected_edge <- numeric()
  unordered <- get(graph, "unordered")
  expected_unordered <- data.frame(i = integer(), j = integer(),
                                   i_to_j = logical(), j_to_i = logical())
  if (stabilized && pate) {
    lo <- pmin(q, j)
    hi <- pmax(q, j)
    pair_key <- paste(lo, hi, sep = ":")
    unique_key <- unique(pair_key)
    first <- match(unique_key, pair_key)
    pair_index <- match(pair_key, unique_key)
    expected_edge <- .wm_sum_at(pair_index,
      (2L * Z[j] - 1L) * r[j] * (w[q] - own_rho[q]) / M,
      length(unique_key)) / gamma
    expected_unordered <- data.frame(i = lo[first], j = hi[first],
      i_to_j = .wm_sum_at(pair_index, as.numeric(q == lo), length(unique_key)) > 0,
      j_to_i = .wm_sum_at(pair_index, as.numeric(q == hi), length(unique_key)) > 0)
  }
  stored_edge <- .wm_numeric_vector(get(contributions, "edge"),
                                    length(expected_edge), "edge contributions")
  .wm_reciprocal_agree(stored_edge, expected_edge, "Decomposition edges")
  if (!is.data.frame(unordered) || !identical(names(unordered), names(expected_unordered)) ||
      nrow(unordered) != nrow(expected_unordered) ||
      !all(vapply(names(expected_unordered), function(column) {
        identical(unname(unordered[[column]]), unname(expected_unordered[[column]]))
      }, logical(1)))) {
    stop("Stored unordered edge structure is inconsistent with directed edges.", call. = FALSE)
  }
  # A second reconstruction uses the query form rather than incoming loads.
  missing_residual <- .wm_sum_at(q, expected_outcome_share * residual[j], n)
  query_contrast <- if (pate) (2L * Z - 1L) *
    (Y - ifelse(Z == 1L, mu0, mu1) - missing_residual) else Y - mu0 - missing_residual
  direct_numerator <- sum(outer * query_contrast)
  .wm_reciprocal_agree(get(object, "numerator"), direct_numerator, "Point numerator")
  .wm_reciprocal_agree(estimate, direct_numerator / sum(outer), "Point estimate")
  list(n = n,
       M = M,
       estimand = estimand,
       method = method,
       gamma = gamma,
       weight_scale = weight_scale,
       estimate = estimate,
       actual = actual,
       q = q,
       j = j,
       Z = Z,
       residual = residual,
       w = w,
       expected_outcome_share = expected_outcome_share,
       donor_total = donor_total,
       conf.level = conf.level,
       variance_floor = variance_floor)
}

# Pure actual-graph pair arithmetic; no interval, floor or availability gate.
# The same once-oriented products serve fixed-map and qualified fitted inference.
.wm_reciprocal_pair_state <- function(state) {
  n <- state$n
  q <- state$q
  j <- state$j
  Z <- state$Z
  w <- state$w
  residual <- state$residual
  expected_outcome_share <- state$expected_outcome_share
  donor_total <- state$donor_total
  # Each reciprocal pair appears once, oriented treated -> control.
  forward <- which(Z[q] == 1L)
  key <- paste(q, j, sep = ":")
  reverse <- match(paste(j[forward], q[forward], sep = ":"), key)
  keep <- !is.na(reverse)
  forward <- forward[keep]
  reverse <- reverse[keep]
  edge_mark <- w[q] * expected_outcome_share * residual[j]
  pair_value <- edge_mark[forward] * edge_mark[reverse]
  pairs <- data.frame(treated = q[forward], control = j[forward],
    treated_to_control_edge = forward, control_to_treated_edge = reverse,
    control_donor_total = donor_total[q[forward]],
    treated_donor_total = donor_total[j[forward]],
    treated_residual = residual[q[forward]], control_residual = residual[j[forward]],
    treated_to_control_mark = edge_mark[forward],
    control_to_treated_mark = edge_mark[reverse],
    product = pair_value, numerator_contribution = pair_value / n)
  reciprocal <- sum(pair_value) / n
  if (any(!is.finite(c(edge_mark, pair_value, reciprocal)))) {
    stop("Reciprocal variance arithmetic exceeded numerical range.", call. = FALSE)
  }
  list(edge_mark = edge_mark, pair_value = pair_value,
       pairs = pairs, reciprocal = reciprocal)
}

#' Reciprocal covariance inference for a corrected fixed-map matching fit
#' @export
wm_reciprocal_inference <- function(object, conf.level = 0.95,
                                    variance_floor = NULL) {
  state <- .wm_reciprocal_validated_state(object, conf.level, variance_floor)
  n <- state[["n"]]
  M <- state[["M"]]
  estimand <- state[["estimand"]]
  method <- state[["method"]]
  gamma <- state[["gamma"]]
  weight_scale <- state[["weight_scale"]]
  estimate <- state[["estimate"]]
  actual <- state[["actual"]]
  q <- state[["q"]]
  j <- state[["j"]]
  Z <- state[["Z"]]
  residual <- state[["residual"]]
  w <- state[["w"]]
  expected_outcome_share <- state[["expected_outcome_share"]]
  donor_total <- state[["donor_total"]]
  conf.level <- state[["conf.level"]]
  variance_floor <- state[["variance_floor"]]
  pair_state <- .wm_reciprocal_pair_state(state)
  edge_mark <- pair_state$edge_mark
  pair_value <- pair_state$pair_value
  pairs <- pair_state$pairs
  reciprocal <- pair_state$reciprocal
  diagonal <- mean((actual * gamma)^2)
  signed <- diagonal - 2 * reciprocal
  raw_root <- signed / gamma^2
  raw_variance <- raw_root / n
  if (any(!is.finite(c(edge_mark, pair_value, reciprocal, diagonal, signed,
                       raw_root, raw_variance)))) {
    stop("Reciprocal variance arithmetic exceeded numerical range.", call. = FALSE)
  }
  floor_active <- !is.null(variance_floor) && signed < variance_floor
  used <- if (floor_active) variance_floor else signed
  available <- used > 0
  root_variance <- if (available) used / gamma^2 else NA_real_
  variance <- if (available) root_variance / n else NA_real_
  se <- if (available) sqrt(variance) else NA_real_
  ci <- if (available) estimate + c(-1, 1) *
    stats::qnorm((1 - conf.level) / 2, lower.tail = FALSE) * se else rep(NA_real_, 2L)
  names(ci) <- c("lower", "upper")
  if (available && (any(!is.finite(c(root_variance, variance, se, ci))) ||
                    root_variance <= 0 || variance <= 0 || se <= 0)) {
    stop("The requested inference scale or confidence interval exceeds numerical range.", call. = FALSE)
  }
  structure(list(estimate = estimate, n = n, M = M, estimand = estimand,
    method = method, gamma = gamma, weight_scale = weight_scale,
    numerator_variance = signed, used_numerator_variance = used,
    reciprocal_numerator = reciprocal, diagonal_numerator_variance = diagonal,
    unregularized_root_n_variance = raw_root, unregularized_variance = raw_variance,
    root_n_variance = root_variance, variance = variance, se = se, conf.int = ci,
    conf.level = conf.level, available = available, variance_floor = variance_floor,
    floor_active = floor_active, reciprocal_pairs = pairs,
    actual_contributions = actual,
    diagnostic = if (signed <= 0) {
      if (floor_active) "Nonpositive unregularized variance; explicit floor applied."
      else "Nonpositive unregularized variance; standard error and interval unavailable."
    } else if (floor_active) "Positive unregularized variance below the explicit floor."
    else "Positive unregularized variance; no floor applied.",
    assumptions_verified = FALSE,
    inference_contract = paste(
      "Bounded iid rows, known fixed maps and positive bounded fixed weights,",
      "own-field residual centering, continuous-score geometry, nondegenerate",
      "oracle variance and appropriate nuisance rates are unverified premises.",
      "Feasible Wald validity requires root-n oracle equivalence; tracking the",
      "actual feasible sampling variance additionally requires L2 transfer or",
      "another sufficient uniform-integrability argument."),
    validation = paste(
      "Stored graph structure, selected-edge distances, donor fractions, loads,",
      "point reconstruction and actual rows checked at relative/absolute 1e-10",
      "arithmetic tolerance. No rematching or nearest-neighbor reranking performed.")),
    class = c("wm_reciprocal_inference", "list"))
}
