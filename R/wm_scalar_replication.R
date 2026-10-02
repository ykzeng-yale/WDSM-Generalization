# Separate current-graph contribution replication for the qualified scalar model.
# No fitting, donor search, random-number generation or modification of the point.

.wm_scalar_repl_agree <- function(value, expected, label) {
  if (!is.numeric(value) || is.complex(value) ||
      !identical(dim(value), dim(expected)) || length(value) != length(expected) ||
      anyNA(value) || any(!is.finite(value)) || any(!is.finite(expected)) ||
      any(abs(value - expected) > 1e-10 * pmax(1, abs(value), abs(expected)))) {
    stop(label, " disagrees with the bound inputs.", call. = FALSE)
  }
  invisible(TRUE)
}

.wm_scalar_repl_graph <- function(fit, probability) {
  n <- fit$n; M <- fit$M; Z <- fit$data$Z; Y <- fit$data$Y
  w <- fit$analysis_weights
  pate <- identical(fit$estimand, "PATE")
  graph <- fit$graph
  if (!is.list(graph) || length(graph$cell) != n || anyNA(graph$cell) ||
      length(unique(graph$cell)) != 1L ||
      !is.matrix(graph$scores0) || !identical(dim(graph$scores0), c(n, 1L)) ||
      !identical(as.numeric(graph$scores0), probability) ||
      (pate && (!is.matrix(graph$scores1) ||
        !identical(dim(graph$scores1), c(n, 1L)) ||
        !identical(as.numeric(graph$scores1), probability))) ||
      (!pate && !is.null(graph$scores1))) {
    stop("Stored graph must use the same unrestricted scalar probability.", call. = FALSE)
  }
  query <- if (pate) seq_len(n) else which(Z == 1L)
  neighbors <- graph$neighbors
  if (!is.list(neighbors) || length(neighbors) != n) stop("Invalid stored neighbor list.")
  for (i in query) {
    j <- neighbors[[i]]
    if (!is.numeric(j) || is.complex(j) || length(j) != M || anyNA(j) ||
        any(!is.finite(j)) || any(j != floor(j) | j < 1 | j > n) ||
        anyDuplicated(j) || any(Z[j] == Z[i])) stop("Invalid donor list at row ", i, ".")
  }
  if (!pate && any(lengths(neighbors[Z == 0L]) != 0L)) stop("PATT cannot contain control queries.")
  edges <- graph$edges
  keys <- c("query", "donor", "arm", "share", "outcome_share")
  if (!is.data.frame(edges) || !all(keys %in% names(edges)) ||
      anyDuplicated(names(edges)) ||
      !all(vapply(edges[keys], function(x) is.numeric(x) && !is.complex(x) &&
        !anyNA(x) && all(is.finite(x)), logical(1)))) stop("Invalid stored edges.")
  q <- rep(query, each = M)
  j <- unlist(neighbors[query], use.names = FALSE)
  if (!identical(as.numeric(edges$query), as.numeric(q)) ||
      !identical(as.numeric(edges$donor), as.numeric(j)) ||
      !identical(as.numeric(edges$arm), as.numeric(1L - Z[q]))) {
    stop("Stored edges and donor lists disagree.", call. = FALSE)
  }
  share <- w[j] / stats::ave(w[j], q, FUN = sum)
  .wm_scalar_repl_agree(edges$share, share, "Donor shares")
  .wm_scalar_repl_agree(edges$outcome_share, share, "Raw donor coefficients")
  incoming <- matrix(0, n, 2L)
  for (z in if (pate) 0:1 else 0L) {
    take <- edges$arm == z
    incoming[, z + 1L] <- .wm_sum_at(j[take], w[q[take]] * share[take], n)
  }
  .wm_scalar_repl_agree(fit$loads$incoming, incoming, "Incoming loads")
  outer <- if (pate) w else Z * w
  .wm_scalar_repl_agree(fit$gamma, mean(outer), "Target normalizer")
  .wm_scalar_repl_agree(fit$denominator, sum(outer), "Target denominator")
  missing <- .wm_sum_at(q, share * Y[j], n)
  contrast <- (2 * Z[query] - 1) * (Y[query] - missing[query])
  point <- sum(w[query] * contrast) / sum(w[query])
  .wm_scalar_repl_agree(fit$estimate, point, "Raw point")
  .wm_scalar_repl_agree(fit$raw_estimate, point, "Raw point record")
  incoming
}

wm_scalar_replication <- function(object, design = NULL, counts = NULL,
                                  multipliers = NULL, conf.level = 0.95) {
  if (!inherits(object, "wm_scalar_logistic_match") || !is.list(object) ||
      !identical(object$status, "point_computed") || !inherits(object$fit, "wm_match")) {
    stop("Supply a computed wm_scalar_logistic_match object.", call. = FALSE)
  }
  fit <- object$fit
  n <- .wm_numeric_vector(fit$n, 1L, "n")
  M <- .wm_numeric_vector(fit$M, 1L, "M")
  if (n < 3 || n != floor(n) || n > .Machine$integer.max ||
      M < 1 || M != floor(M) || M > n ||
      !fit$estimand %in% c("PATE", "PATT") ||
      !identical(fit$method, "self_normalized") ||
      !identical(fit$info$corrected, FALSE) ||
      !identical(fit$info$nuisance_correction, FALSE)) {
    stop("Require an unadjusted raw scalar self-normalized matching point.", call. = FALSE)
  }
  n <- as.integer(n)
  conf.level <- .wm_numeric_vector(conf.level, 1L, "conf.level")
  if (conf.level <= 0 || conf.level >= 1) stop("conf.level must lie between zero and one.")
  if (!is.null(counts) && !is.null(multipliers)) stop("Supply counts or multipliers, not both.")
  draw_input <- if (!is.null(counts)) "multinomial_counts" else
    if (!is.null(multipliers)) "gaussian_multipliers" else "none"
  draws <- if (!is.null(counts)) counts else multipliers
  if (!is.null(draws)) {
    draws <- .wm_score_matrix(draws, n, "replication input")
    if (!is.null(counts) && (any(draws < 0 | draws != floor(draws)) ||
        any(colSums(draws) != n))) stop("Counts must be nonnegative integers with column sums n.")
  }
  estimate <- .wm_numeric_vector(object$estimate, 1L, "estimate")
  .wm_scalar_repl_agree(fit$estimate, estimate, "Point records")
  base <- list(status = "unavailable", reason = NULL, estimate = estimate,
    fit = fit, source_inference = object$inference, n = n, M = as.integer(M),
    estimand = fit$estimand, root_n_variance = NA_real_, variance = NA_real_,
    se = NA_real_, conf.int = c(lower = NA_real_, upper = NA_real_),
    conf.level = conf.level, interval_available = FALSE, row_variance_finite = FALSE,
    draw_input = draw_input, B = if (is.null(draws)) 0L else ncol(draws),
    supplied_draws = draws, design_binding = "unbound",
    numerical_root_bridge = object$numerical_root_bridge,
    population_assumptions_verified = FALSE, supplied_draw_law_verified = FALSE,
    original_refit_bootstrap = FALSE,
    scope = paste("Full signed-drift contribution replication under the bounded ordinary-logistic",
      "scalar known-weight model and its current-graph moment theorem.",
      "Population and numerical-root premises remain required; no fitting or rematching."))
  finish <- function(x) structure(x, class = c("wm_scalar_replication", "list"))
  if (is.null(object$fitting_design) ||
      !identical(object$design_binding, "validated_input_matrix_v1")) {
    base$reason <- "No stored fitting-design binding; supplied design or J/B agreement cannot establish historical row identity."
    return(finish(base))
  }
  columns <- object$design_columns
  bound <- object$fitting_design
  if (!is.character(columns) || anyNA(columns) || any(!nzchar(columns)) ||
      anyDuplicated(columns) || !identical(colnames(bound), columns)) stop("Invalid stored design columns.")
  bound <- .wm_score_matrix(bound, n, "stored fitting design")
  p <- ncol(bound)
  if (p < 3L || p >= n || length(columns) != p || any(bound[, 1L] != 1) ||
      qr(bound)$rank != p) stop("The bound scalar inference design needs at least three full-rank columns and a leading intercept.")
  if (!is.null(design)) {
    provided_names <- colnames(design)
    checked <- .wm_score_matrix(design, n, "design")
    if ((!is.null(provided_names) && !identical(provided_names, columns)) ||
        !identical(checked, bound)) stop("Supplied design differs from the stored fitting-design binding.")
  }
  design <- bound
  colnames(design) <- columns
  base$design <- design
  base$design_binding <- object$design_binding
  Y <- .wm_numeric_vector(fit$data$Y, n, "stored Y")
  Z <- .wm_numeric_vector(as.numeric(fit$data$Z), n, "stored Z")
  if (any(!Z %in% c(0, 1)) || length(unique(Z)) != 2L) stop("Stored treatment must contain both binary arms.")
  W <- .wm_numeric_vector(fit$weights, n, "stored weights", TRUE)
  w <- .wm_numeric_vector(fit$analysis_weights, n, "analysis weights", TRUE)
  scale <- .wm_numeric_vector(fit$weight_scale, 1L, "weight scale", TRUE)
  .wm_scalar_repl_agree(scale, max(W), "Weight scale")
  .wm_scalar_repl_agree(w, W / scale, "Analysis weights")
  e <- .wm_numeric_vector(object$propensity$probability, n, "stored probability")
  if (any(e <= 0 | e >= 1)) stop("Stored probabilities must lie strictly between zero and one.")
  incoming <- .wm_scalar_repl_graph(fit, e)
  inf <- object$inference
  blocks <- inf$blocks
  if (!is.list(inf) || !is.list(blocks) ||
      !all(c("V0", "b", "C", "J", "B", "Sigma") %in% names(blocks)) ||
      is.null(inf$predictions$prediction)) {
    base$reason <- "The original scalar inference record lacks predictions or full nuisance blocks; no model will be refitted."
    return(finish(base))
  }
  if (!identical(inf$model_contract, "bounded_scalar_logistic_known_score_weights")) stop("Unsupported scalar inference contract.")
  .wm_scalar_repl_agree(inf$moments$weight_scale, scale, "Inference weight scale")
  mu <- .wm_score_matrix(inf$predictions$prediction, n, "variance predictions")
  if (ncol(mu) != 2L) stop("Variance predictions must have two columns, mean0 and mean1.")
  if (!identical(colnames(inf$predictions$prediction), c("mean0", "mean1"))) stop("Variance prediction column order is invalid.")
  if (!identical(names(blocks$b), columns) || !identical(names(blocks$C), columns)) stop("Nuisance vector coordinates do not match the bound design.")
  b <- .wm_numeric_vector(blocks$b, p, "full graph drift b")
  .wm_numeric_vector(blocks$C, p, "explicit covariance C")
  .wm_numeric_vector(blocks$V0, 1L, "explicit V0")
  J <- crossprod(design, design * (w * e * (1 - e))) / n
  psi <- design * (w * (Z - e))
  B <- crossprod(psi) / n
  .wm_scalar_repl_agree(blocks$J, J, "Logistic information J")
  .wm_scalar_repl_agree(blocks$B, B, "Logistic score covariance B")
  calculate <- function() {
    chol(J)
    phi <- t(solve(J, t(psi)))
    .wm_scalar_repl_agree(blocks$Sigma, crossprod(phi) / n, "Explicit nuisance sandwich")
    K <- incoming[cbind(seq_len(n), Z + 1L)]
    gamma <- if (fit$estimand == "PATE") mean(w) else mean(Z * w)
    L <- if (fit$estimand == "PATE") {
      (w * (mu[, 2L] - mu[, 1L] - estimate) +
         (2 * Z - 1) * (w + K) * (Y - mu[cbind(seq_len(n), Z + 1L)])) / gamma
    } else {
      (Z * w * (Y - mu[, 1L] - estimate) -
         (1 - Z) * incoming[, 1L] * (Y - mu[, 1L])) / gamma
    }
    Lc <- L - mean(L)
    phic <- sweep(phi, 2L, colMeans(phi), "-")
    augmented <- Lc + drop(phic %*% b)
    augmented_center <- mean(augmented)
    augmented <- augmented - augmented_center
    V0_row <- mean(Lc^2)
    C_row <- drop(crossprod(phic, Lc)) / n
    Sigma_row <- crossprod(phic) / n
    cross <- 2 * sum(b * C_row)
    quadratic <- drop(crossprod(b, Sigma_row %*% b))
    V <- mean(augmented^2)
    if (any(!is.finite(c(phi, L, augmented, V0_row, C_row, Sigma_row,
                        cross, quadratic, V)))) stop("Row arithmetic exceeded numerical range.")
    .wm_scalar_repl_agree(V, V0_row + cross + quadratic, "Augmented-row variance identity")
    explicit <- .wm_numeric_vector(inf$root_n_variance, 1L, "explicit root variance")
    interval <- .wm_scalar_interval(estimate, V, n, conf.level)
    root <- point <- NULL
    empirical <- NA_real_
    if (!is.null(draws)) {
      root <- drop(crossprod(augmented, if (draw_input == "multinomial_counts") draws - 1 else draws)) / sqrt(n)
      point <- estimate + root / sqrt(n)
      if (any(!is.finite(c(root, point)))) stop("Replication arithmetic exceeded numerical range.")
      if (length(point) >= 2L) empirical <- mean((point - mean(point))^2)
      if (length(point) >= 2L && !is.finite(empirical)) stop("Empirical draw variance exceeded numerical range.")
    }
    out <- base
    out$status <- "row_variance_computed"
    out$reason <- NULL
    for (name in setdiff(names(interval), "status")) out[[name]] <- interval[[name]]
    out$interval_status <- interval$status
    out$interval_available <- identical(interval$status, "conditional_formula")
    out$row_variance_finite <- is.finite(V) && V >= 0
    out$rows <- list(contribution = L, contribution_centered = Lc, influence = phi,
      influence_centered = phic, augmented = augmented, b = stats::setNames(b, columns),
      original_row_mean = mean(L), original_influence_mean = colMeans(phi),
      final_roundoff_center = augmented_center)
    out$blocks <- list(V0_row = V0_row, C_row = stats::setNames(C_row, columns),
      Sigma_row = Sigma_row, cross = cross, nuisance = quadratic,
      expansion_root_n_variance = V0_row + cross + quadratic)
    out$comparison <- list(explicit_root_n_variance = explicit,
      row_root_n_variance = V, row_minus_explicit = V - explicit,
      finite_sample_equality_expected = FALSE)
    out$normalization <- list(gamma = gamma, weight_scale = scale,
      sample_size = n, convention = "Full observed n for PATE and PATT; all loads and weights share one scale")
    out$root_n_draws <- root
    out$draws <- point
    out$empirical_draw_variance_B <- empirical
    out$empirical_variance_divisor <- if (length(point) >= 2L) "B" else "not_calculated"
    out$conditional_root_variance <- V
    out$conditional_point_variance <- V / n
    out$conditional_law_requirement <- switch(draw_input,
      gaussian_multipliers = "Independent standard Gaussian entries; supplied values alone do not verify this law",
      multinomial_counts = "Each count column is Multinomial(n, rep(1/n,n)); supplied values alone do not verify this law",
      none = "Independent standard Gaussian multipliers or full-n uniform multinomial counts have this conditional variance; no draws were supplied")
    finish(out)
  }
  tryCatch(calculate(), error = function(error) {
    base$reason <- conditionMessage(error)
    finish(base)
  })
}
