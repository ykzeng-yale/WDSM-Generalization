# Fixed-array component/preflight checks only: no nuisance fit, simulated data,
# bootstrap or population-theorem/coverage claim.
wm_score_weight_fixture <- function() {
  x <- c(.02, .13, .29, .48, .69, .91, .07, .21, .36, .56, .77, .96)
  u <- c(.2, .8, .4, .1, .9, .5, .6, .3, .7, .4, .2, .8)
  S <- cbind(x = x, u = u)
  B <- array(0, c(12L, 2L, 4L),
    dimnames = list(NULL, colnames(S), c("ps", "pg", "center", "scale")))
  B[, , 1L] <- cbind(1 + x, u)
  B[, , 2L] <- cbind(x * u, 1 - x)
  B[, , 3L] <- cbind(-1, .3 + u)
  B[, , 4L] <- cbind(x^2, -u)
  list(raw_scores = S, Z = rep(0:1, each = 6L),
    W = c(.7, 1.1, 1.4, .9, 1.8, 1.2, 1.3, .8, 1.9, 1.5, 1, 1.7),
    residual = c(1, -2, .4, -.7, 2, .5, -.3, 1.5, -1, .8, 1.2, -.9),
    tangents = B, cutoff = c(1, .7, 0, 1, .5, 1, .6, 1, 1, 0, .8, 1),
    donor_arm = 0L, M = 3L, bandwidth = 1.3, density_floor = 1e-8)
}

# Independent scalar KDE sums, differentiated by central finite differences.
# No candidate kernel, analytic-gradient or contraction helper is used here.
wm_score_weight_numeric_reference <- function(a) {
  n <- nrow(a$raw_scores)
  d <- ncol(a$raw_scores)
  donor <- which(a$Z == a$donor_arm)
  query <- which(a$Z != a$donor_arm)
  active <- donor[a$cutoff[donor] > 0]
  density <- function(s) {
    K <- numeric(n)
    for (i in seq_len(n)) {
      u <- (s - a$raw_scores[i, ]) / a$bandwidth
      K[i] <- if (any(abs(u) >= 1)) 0 else prod((35 / 32) * (1 - u^2)^3)
    }
    c(h = sum(K[donor]), q = sum(K[query] * a$W[query]),
      w = sum(K[donor] * a$W[donor])) / (n * a$bandwidth^d)
  }
  scalar <- function(s) {
    f <- density(s)
    hp <- max(f["h"], a$density_floor)
    wp <- max(f["w"], a$density_floor)
    c(kappa = unname(f["q"] / hp), log_weight = unname(log(wp) - log(hp)))
  }
  grad <- field <- matrix(0, length(active), d)
  j <- numeric(dim(a$tangents)[3L])
  constant <- all(a$W[donor] == a$W[donor[1L]])
  for (r in seq_along(active)) {
    i <- active[r]
    s <- a$raw_scores[i, ]
    kappa <- scalar(s)["kappa"]
    log_grad <- numeric(d)
    for (v in seq_len(d)) {
      plus <- minus <- s
      plus[v] <- plus[v] + 1e-5
      minus[v] <- minus[v] - 1e-5
      delta <- (scalar(plus) - scalar(minus)) / 2e-5
      grad[r, v] <- delta["kappa"]
      if (a$M > 1L && !constant) log_grad[v] <- delta["log_weight"]
    }
    field[r, ] <- grad[r, ] - (1 - 1 / a$M) * kappa * log_grad
    for (p in seq_along(j)) {
      j[p] <- j[p] + a$residual[i] * a$cutoff[i] *
        sum(a$tangents[i, , p] * field[r, ]) / n
    }
  }
  list(gradient = grad, field = field, drift = j)
}

test_that("score-weight drift uses complete p tangents and total-n density gradients", {
  for (arm in 0:1) for (M in c(1L, 3L)) {
    a <- wm_score_weight_fixture()
    a$donor_arm <- arm
    a$M <- M
    before <- a
    out <- do.call(wm_graph_transport, c(a, list(representation = "score_measurable_weight")))
    ref <- wm_score_weight_numeric_reference(a)
    expect_identical(a, before)
    expect_length(out$graph_drift, 4L)
    expect_identical(names(out$graph_drift), dimnames(a$tangents)[[3L]])
    expect_equal(out$incoming_load_gradient, ref$gradient, tolerance = 2e-8)
    expect_equal(out$transport_field, ref$field, tolerance = 2e-8)
    expect_equal(unname(out$graph_drift), ref$drift, tolerance = 2e-8)
    expect_identical(out$quadrature_error_bound,
      stats::setNames(numeric(4L), dimnames(a$tangents)[[3L]]))
    expect_true(all(out$quadrature_nodes == 0))
    expect_false(any(out$donor_density_floor_active))
    expect_false(any(out$donor_weighted_density_floor_active))
    expect_match(out$bound_scope, "does not bound statistical estimation error", fixed = TRUE)
    expect_identical(out$representation, "score_measurable_weight")
    expect_false(out$sampling_inference_available)
    if (M == 1L) {
      expect_identical(out$transport_field, out$incoming_load_gradient)
      expect_identical(out$log_weight_gradient_method, "not_used_M1")
    }
  }
})

test_that("new floors differentiate the maximum and disclose exact constant-weight completion", {
  a <- wm_score_weight_fixture()
  a$density_floor <- 100
  out <- do.call(wm_graph_transport, c(a, list(representation = "score_measurable_weight")))
  ref <- wm_score_weight_numeric_reference(a)
  expect_true(all(out$donor_density_floor_active))
  expect_true(all(out$donor_weighted_density_floor_active))
  expect_equal(out$transport_field, ref$field, tolerance = 2e-8)
  expect_equal(unname(out$graph_drift), ref$drift, tolerance = 2e-8)
  expect_true(all(out$log_donor_weight_gradient == 0))

  # Isolate a one-active-floor event for constant donor weights c != 1.
  a$W[] <- 100
  a$density_floor <- 1e-8
  baseline <- do.call(wm_graph_transport,
    c(a, list(representation = "score_measurable_weight")))
  a$density_floor <- 2 * max(baseline$donor_subdensity)
  out <- do.call(wm_graph_transport, c(a, list(representation = "score_measurable_weight")))
  expect_true(all(out$donor_density_floor_active))
  expect_false(any(out$donor_weighted_density_floor_active))
  expect_true(out$constant_donor_weights)
  expect_identical(out$log_weight_gradient_method, "exact_constant_donor_weights")
  expect_identical(out$transport_field, out$incoming_load_gradient)
  expect_true(all(out$log_donor_weight_gradient == 0))

  # W=1 reduction is exact, with and without active floors, for multiple d.
  for (d in c(1L, 2L, 3L)) for (floor in c(1e-8, 100)) {
    a <- wm_score_weight_fixture()
    a$raw_scores <- cbind(a$raw_scores, v = a$raw_scores[, 1]^2)[, seq_len(d), drop = FALSE]
    a$tangents <- array(rep(seq_len(12L) / 12, d * 4L), c(12L, d, 4L),
      dimnames = list(NULL, colnames(a$raw_scores), c("ps", "pg", "center", "scale")))
    a$W[] <- 1
    a$density_floor <- floor
    out <- do.call(wm_graph_transport, c(a, list(representation = "score_measurable_weight")))
    expect_identical(out$transport_field, out$incoming_load_gradient)
    expect_identical(out$donor_subdensity, out$donor_weighted_subdensity)
  }
})

test_that("old conditional-weight-chart dispatch and default objects are unchanged", {
  a <- wm_score_weight_fixture()
  a$W[] <- 1
  before <- do.call(WeightedMatching:::.wm_gt_compute, a)
  implicit <- do.call(wm_graph_transport, a)
  explicit <- do.call(wm_graph_transport,
    c(a, list(representation = "conditional_weight_chart")))
  expect_identical(implicit, explicit)
  for (field in names(before)) expect_identical(implicit[[field]], before[[field]])
  expect_null(implicit$representation)
  score <- do.call(wm_graph_transport,
    c(a, list(representation = "score_measurable_weight")))
  expect_equal(score$graph_drift, implicit$graph_drift, tolerance = 1e-12)
  expect_equal(score$incoming_load_gradient, implicit$incoming_load_gradient, tolerance = 1e-12)
})

test_that("score-weight fitted specifications bind original coordinates and full parameter order", {
  a <- wm_score_weight_fixture()
  spec <- list(mode = "estimate", representation = "score_measurable_weight",
    raw_scores = a$raw_scores, tangents = a$tangents, cutoff = a$cutoff)
  parameters <- dimnames(a$tangents)[[3L]]
  validate <- function(s, scores = a$raw_scores) {
    WeightedMatching:::.wm_fi_transport_spec(s, scores, 12L, parameters, "transport0")
  }
  good <- validate(spec)
  expect_identical(good$representation, "score_measurable_weight")
  expect_identical(good$map_binding$method, "identity")
  expect_identical(good$tangents, a$tangents)
  bad <- spec; bad$raw_to_matching <- function(x) x
  expect_error(validate(bad), "requires identity matching coordinates", fixed = TRUE)
  bad <- spec; bad$representation <- "unknown"
  expect_error(validate(bad), "representation must", fixed = TRUE)
  bad <- spec; colnames(bad$raw_scores) <- rev(colnames(a$raw_scores))
  dimnames(bad$tangents)[[2L]] <- colnames(bad$raw_scores)
  expect_error(validate(bad), "coordinate names/order must match fit scores", fixed = TRUE)
  bad <- spec; bad$raw_scores[1, 1] <- bad$raw_scores[1, 1] + .01
  expect_error(validate(bad), "raw-to-matching data binding", fixed = TRUE)
  bad <- spec; dimnames(bad$tangents)[[3L]] <- rev(parameters)
  expect_error(validate(bad), "parameter names/order", fixed = TRUE)
  bad <- spec; rownames(bad$raw_scores) <- as.character(12:1)
  expect_error(validate(bad), "original row indices", fixed = TRUE)
  bad <- spec; bad$W <- a$W
  expect_error(validate(bad), "unsupported fields", fixed = TRUE)
  unnamed <- spec; colnames(unnamed$raw_scores) <- NULL
  dimnames(unnamed$tangents)[2L] <- list(NULL)
  expect_identical(validate(unnamed, unname(a$raw_scores))$representation,
    "score_measurable_weight")

  old <- spec; old$representation <- NULL
  explicit <- old; explicit$representation <- "conditional_weight_chart"
  expect_identical(validate(old), validate(explicit))
})

test_that("common inference and bootstrap retain the score-weight PATE/PATT bindings", {
  a <- wm_score_weight_fixture()
  n <- length(a$Z)
  x <- a$raw_scores[, 1L]
  u <- a$raw_scores[, 2L]
  v <- c(.7, .1, .9, .3, .6, .2, .5, .8, .1, .9, .4, .6)
  S1 <- cbind(x = x, v = v)
  Y <- a$residual + .4 * x + .2 * u + a$Z * (.7 + .3 * x)
  mean0 <- .2 + .5 * x - .1 * u
  mean1 <- -.3 + .7 * x + .2 * v
  influence <- cbind(ps = x - .4, pg = u - .5, center = v - .6,
                     scale = x * u - .2)
  D0 <- cbind(ps = .1 + x, pg = u^2, center = -.3 + v, scale = x * v)
  D1 <- cbind(ps = u, pg = .2 + x^2, center = -.4 + u, scale = v^2)
  counts <- matrix(1, n, 2L)
  counts[1:2, 1L] <- c(2, 0)
  counts[3:4, 2L] <- c(0, 2)
  parameters <- colnames(influence)

  # These finite maps/derivatives are supplied arrays, not fitted nuisance
  # models. Scope declarations below exercise API structure only and do not
  # certify population geometry, score measurability or bootstrap coverage.
  for (estimand in c("PATE", "PATT")) {
    pate <- identical(estimand, "PATE")
    fit <- wm_match(Y, a$Z, a$W, a$raw_scores,
      scores1 = if (pate) S1 else NULL, M = a$M, estimand = estimand,
      mean0 = mean0, mean1 = if (pate) mean1 else NULL, variance = FALSE)
    before_fit <- fit
    B0 <- a$tangents
    dimnames(B0)[2L] <- list(colnames(fit$graph$scores0))
    spec0 <- list(mode = "estimate", representation = "score_measurable_weight",
      raw_scores = fit$graph$scores0, tangents = B0, cutoff = a$cutoff,
      bandwidth = a$bandwidth, density_floor = a$density_floor)
    args <- list(fit = fit, nuisance_influence = influence,
      mean_derivative0 = D0, transport0 = spec0,
      covariance_scope = if (pate) "distinct_rarity" else "patt")
    r0 <- a
    r0$raw_scores <- fit$graph$scores0
    r0$tangents <- B0
    r0$residual <- Y - mean0
    j0 <- wm_score_weight_numeric_reference(r0)$drift
    gamma <- mean(if (pate) a$W else a$Z * a$W)
    graph <- -j0 / gamma
    if (pate) {
      B1 <- B0
      B1[, 2L, ] <- B1[, 2L, ] + outer(v, c(.2, -.3, .4, -.1))
      dimnames(B1)[2L] <- list(colnames(fit$graph$scores1))
      spec1 <- spec0
      spec1$raw_scores <- fit$graph$scores1
      spec1$tangents <- B1
      args$mean_derivative1 <- D1
      args$transport1 <- spec1
      r1 <- a
      r1$raw_scores <- fit$graph$scores1
      r1$tangents <- B1
      r1$residual <- Y - mean1
      r1$donor_arm <- 1L
      graph <- graph + wm_score_weight_numeric_reference(r1)$drift / gamma
    }
    out <- do.call(wm_fitted_inference, args)

    # Independent query-by-query fixed-edge derivative with original W;
    # normalization is the original target total, not an arm-specific n.
    smooth <- numeric(4L)
    queries <- if (pate) seq_len(n) else which(a$Z == 1L)
    for (i in queries) {
      donors <- fit$graph$neighbors[[i]]
      fraction <- a$W[donors] / sum(a$W[donors])
      D <- if (a$Z[i] == 1L) D0 else D1
      smooth <- smooth + a$W[i] * (1 - 2 * a$Z[i]) *
        (D[i, ] - colSums(D[donors, , drop = FALSE] * fraction)) / (n * gamma)
    }
    expect_identical(fit, before_fit)
    expect_identical(out$fit, before_fit)
    expect_identical(out$fit$weights, a$W)
    expect_identical(out$fit$graph, before_fit$graph)
    expect_identical(out$estimate, before_fit$estimate)
    expect_identical(out$status, paste0("conditional_", args$covariance_scope))
    expect_true(out$available)
    expect_identical(out$parameter_names, parameters)
    expect_identical(dim(out$nuisance_influence), c(n, 4L))
    expect_identical(dimnames(out$Sigma), list(parameters, parameters))
    expect_length(out$C, 4L)
    expect_identical(out$raw_target_weight_mean, gamma)
    expect_equal(unname(out$graph_sensitivity), unname(graph), tolerance = 2e-8)
    expect_equal(unname(out$smooth_sensitivity), unname(smooth), tolerance = 1e-12)
    expect_equal(unname(out$total_sensitivity), unname(smooth + graph), tolerance = 2e-8)
    expect_identical(names(out$graph_transport),
      if (pate) c("potential0", "potential1") else "potential0")
    expect_true(all(vapply(out$graph_transport, function(tr) {
      identical(tr$mode, "estimate") &&
        identical(tr$representation, "score_measurable_weight") &&
        identical(tr$parameter_names, parameters) &&
        identical(tr$map_binding$method, "identity")
    }, logical(1))))

    # Two deterministic count columns exercise the existing structural
    # validator/dispatch without RNG, nuisance refits or rematching.
    before_inference <- out
    boot <- wm_bootstrap(out, B = 2L, counts = counts, interval = "none")
    rows <- out$augmented_rows - mean(out$augmented_rows)
    expected_roots <- as.vector(crossprod(rows, counts - 1)) / sqrt(n)
    expect_identical(out, before_inference)
    expect_identical(boot$source_inference, before_inference)
    expect_identical(boot$source_inference$fit, before_fit)
    expect_identical(boot$estimate, before_fit$estimate)
    expect_identical(boot$B, 2L)
    expect_identical(boot$draw_input, "supplied_full_n_multinomial_counts")
    expect_identical(boot$supplied_counts, counts)
    expect_equal(boot$root_n_draws, expected_roots, tolerance = 1e-12)
    expect_equal(boot$draws, before_fit$estimate + expected_roots / sqrt(n),
      tolerance = 1e-12)
    expect_equal(boot$conditional_root_n_variance, out$root_n_variance, tolerance = 1e-12)
    expect_false(boot$supplied_draw_law_verified)
    expect_false(boot$population_assumptions_verified)
    expect_false(boot$original_refit_bootstrap)
    expect_null(boot$conf.int)
  }
})
