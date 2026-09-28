# One-direction stabilized potential means use row contributions only.
.wm_potential_arm <- function(Z, arm, n) {
  if (!is.numeric(arm) || is.complex(arm) || length(arm) != 1L ||
      is.na(arm) || !arm %in% c(0, 1)) stop("arm must be 0 or 1", call. = FALSE)
  if ((!is.numeric(Z) && !is.logical(Z)) || is.complex(Z) ||
      !is.null(dim(Z)) || length(Z) != n || anyNA(Z) ||
      any(!Z %in% c(0, 1)) || length(unique(Z)) != 2L)
    stop("Z must contain both arms and only 0 and 1", call. = FALSE)
  as.integer(Z == arm)
}

.wm_potential_convert <- function(fit, Y, Z, arm, variance) {
  n <- fit$n
  w <- fit$analysis_weights
  donor <- as.integer(Z == arm)
  mu <- fit$predictions$mean0
  rho <- fit$predictions$rho0
  C <- fit$loads$incoming[, 1L]
  gamma <- mean(w)
  imputed <- fit$imputed[, 1L]
  numerator <- sum(w * imputed)
  estimate <- numerator / sum(w)
  eta <- numeric(n)
  own <- donor == 1L
  eta[own] <- w[own] * (Y[own] - mu[own]) / rho[own]
  uncentered <- w * mu + donor * (rho + C) * eta
  row_numerator <- uncentered - w * estimate
  row <- row_numerator / gamma
  if (any(!is.finite(c(estimate, row)))) stop("potential-mean arithmetic exceeded numerical range")
  root_n_variance <- if (variance) mean(row^2) else NA_real_
  if (variance && !is.finite(root_n_variance)) stop("potential-mean variance exceeded numerical range")
  graph <- fit$graph
  graph$edges$arm <- rep.int(as.integer(arm), nrow(graph$edges))
  graph$scores <- graph$scores0
  graph$scores0 <- graph$scores1 <- NULL
  graph$donor_arm <- as.integer(arm)
  ans <- list(estimate = estimate, n = n, M = fit$M,
    estimand = paste0("potential_mean_", arm), arm = as.integer(arm), method = "stabilized",
    root_n_variance = root_n_variance, variance = root_n_variance / n,
    se = sqrt(root_n_variance / n), numerator = numerator, denominator = sum(w), gamma = gamma,
    weights = fit$weights, analysis_weights = w, weight_scale = fit$weight_scale,
    data = list(Y = Y, Z = as.integer(Z)), imputed = imputed,
    predictions = list(mean = mu, rho = rho), loads = list(incoming = C),
    graph = graph, contributions = list(row = row, edge = numeric(), actual = row, nuisance = numeric(n)),
    numerator_identity = list(direct = numerator, reconstructed = sum(uncentered),
      difference = numerator - sum(uncentered), centered_sum = sum(row_numerator)),
    info = list(inference_status = if (variance) "user_supplied_nuisance_contract" else "point_estimate_only",
      contract = "One-direction stabilized potential mean: reduced-mean identification, weighted centering, finite p>2 moments, fixed-map geometry and the mean/rho nuisance-rate contract.",
      weight_units = "rho, incoming load, numerator and denominator use weights / weight_scale.",
      nuisance_correction = FALSE))
  if (!is.null(fit$nuisance)) {
    ans$nuisance <- fit$nuisance
    ans$nuisance$donor_arm <- as.integer(arm)
    ans$nuisance$internal_arm_coding <- "Nuisance model arm 0 denotes the requested potential-outcome arm."
    ans$nuisance$inference_contract$bounded_outcomes <- NULL
    ans$nuisance$inference_contract$moment_scope <- "The potential-mean feasible branch allows the declared finite p>2 moments; spline rates additionally require bounded conditional residual variance."
  }
  structure(ans, class = c("wm_potential", "wm_match", "list"))
}

#' Stabilized one-direction potential mean
#' @export
wm_potential <- function(Y, Z, weights, scores, arm, M = 3L,
                         mean, rho, fold_id = NULL, strata = NULL,
                         variance = TRUE) {
  donor <- .wm_potential_arm(Z, arm, length(Y))
  Y <- .wm_numeric_vector(Y, length(Y), "Y")
  graph_Y <- Y
  graph_Y[donor == 0L] <- 0
  fit <- wm_match(graph_Y, 1L - donor, weights, scores, M = M,
                  estimand = "PATT", method = "stabilized", mean0 = mean, rho0 = rho,
                  fold_id = fold_id, strata = strata, variance = variance)
  .wm_potential_convert(fit, Y, Z, arm, variance)
}

#' Fit a stabilized one-direction potential mean
#' @export
wm_potential_fit <- function(Y, Z, weights, scores, arm, M = 3L,
                             regression = c("polynomial", "spline"), degree = 1L,
                             spline_order = NULL, mesh_exponent = NULL,
                             moment_order = NULL, support = NULL, rho_bounds,
                             fold_id = NULL, strata = NULL, variance = TRUE,
                             max_basis = 10000L, max_basis_elements = 1e7,
                             qr_tol = 1e-10) {
  donor <- .wm_potential_arm(Z, arm, length(Y))
  Y <- .wm_numeric_vector(Y, length(Y), "Y")
  graph_Y <- Y
  graph_Y[donor == 0L] <- 0
  fit <- wm_fit(graph_Y, 1L - donor, weights, scores, M = M,
                estimand = "PATT", method = "stabilized", regression = regression,
                degree = degree, spline_order = spline_order, mesh_exponent = mesh_exponent,
                moment_order = moment_order, support0 = support, rho_bounds = rho_bounds,
                fold_id = fold_id, strata = strata, variance = variance,
                max_basis = max_basis, max_basis_elements = max_basis_elements, qr_tol = qr_tol)
  .wm_potential_convert(fit, Y, Z, arm, variance)
}

#' Independent-block stabilized PATE with potentially different matching maps
#' @export
wm_split_pate <- function(Y, Z, weights, scores0, scores1, block_id,
                          mean0, mean1, rho0, rho1, M = 3L,
                          strata = NULL, variance = TRUE) {
  n <- length(Y)
  .wm_potential_arm(Z, 0, n)
  Y <- .wm_numeric_vector(Y, n, "Y")
  weights <- .wm_numeric_vector(weights, n, "weights", positive = TRUE)
  if ((!is.numeric(block_id) && !is.logical(block_id)) || is.complex(block_id) ||
      !is.null(dim(block_id)) || length(block_id) != n || anyNA(block_id) ||
      any(!block_id %in% c(0, 1)) || length(unique(block_id)) != 2L)
    stop("block_id must specify both independent evaluation blocks, coded 0 and 1", call. = FALSE)
  scores0 <- .wm_score_matrix(scores0, n, "scores0")
  scores1 <- .wm_score_matrix(scores1, n, "scores1")
  mean0 <- .wm_numeric_vector(mean0, n, "mean0")
  mean1 <- .wm_numeric_vector(mean1, n, "mean1")
  rho0 <- .wm_numeric_vector(rho0, n, "rho0", positive = TRUE)
  rho1 <- .wm_numeric_vector(rho1, n, "rho1", positive = TRUE)
  labels <- .wm_cell_labels(strata, n, "strata")
  components <- vector("list", 2L)
  rows <- numeric(n)
  for (arm in 0:1) {
    index <- which(block_id == arm)
    S <- if (arm == 0L) scores0 else scores1
    mu <- if (arm == 0L) mean0 else mean1
    rho <- if (arm == 0L) rho0 else rho1
    fit <- wm_potential(Y[index], Z[index], weights[index], S[index, , drop = FALSE],
                         arm = arm, M = M, mean = mu[index], rho = rho[index],
                         strata = labels[index], variance = variance)
    fit$global_rows <- index
    components[[arm + 1L]] <- fit
    rows[index] <- (2L * arm - 1L) * n / length(index) * fit$contributions$row
  }
  estimate <- components[[2L]]$estimate - components[[1L]]$estimate
  root_n_variance <- if (variance) mean(rows^2) else NA_real_
  if (!is.finite(estimate) || (variance && !is.finite(root_n_variance)))
    stop("split-estimator arithmetic exceeded numerical range")
  names(components) <- c("potential0", "potential1")
  structure(list(estimate = estimate, n = n, M = M,
    estimand = "PATE_independent_blocks", method = "stabilized",
    root_n_variance = root_n_variance, variance = root_n_variance / n,
    se = sqrt(root_n_variance / n), block_id = as.integer(block_id), components = components,
    contributions = list(row = rows, edge = numeric(), actual = rows, nuisance = numeric(n)),
    info = list(inference_status = if (variance) "user_supplied_independent_block_contract" else "point_estimate_only",
      contract = "Fixed data-independent evaluation blocks, both arms in each block, independent block observations and externally supplied or separately trained nuisance predictions satisfying each potential-mean contract. Full-data nuisance training requires a separate argument.",
      precision = "Each arm-specific potential mean uses its own block and denominator; this is a different estimator from full-sample PATE.")),
    class = c("wm_split_pate", "wm_match", "list"))
}
