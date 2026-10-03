# One actual finite nested fit per target checks the new observable law inputs.
# This is regression/graph/covariance verification, not a coverage experiment.
wm_qc_fixture <- function(target) {
  set.seed(20391)
  n <- 320L
  x <- cbind(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1), x3 = runif(n, -1, 1), x4 = runif(n, -1, 1))
  small <- cbind(intercept = 1, x1 = x[, 1])
  large <- cbind(added = x[, 3], small, extra = x[, 4]) # shared columns deliberately reordered
  g0 <- cbind(intercept = 1, x2 = x[, 2])
  g1 <- cbind(g0, x3 = x[, 3])
  z <- rbinom(n, 1, plogis(-.2 + .5 * x[, 1]))
  y <- 1 + .7 * x[, 2] + z * (.8 - .3 * x[, 2] + .5 * x[, 3]) + runif(n, -.4, .4)
  w <- exp(.4 * z + .2 * x[, 1])
  wm_model_fit(y, z, w,
    ps_models = list(larger = list(design = large), smaller = list(design = small)),
    pg0_models = list(control = list(design = g0)),
    pg1_models = if (target == "PATE") list(treated = list(design = g1)) else NULL,
    M = 3L, estimand = target, correction = list(method = "quadratic_qr"))
}

wm_qc_manual <- function(object, h, small = 2L, large = 1L, graph = NULL) {
  s <- object$nuisance; n <- s$n; z <- s$inputs$Z; w <- s$inputs$weights / max(s$inputs$weights)
  D <- s$inputs$ps_models[[small]]$design
  large_design <- s$inputs$ps_models[[large]]$design
  added <- which(!colnames(large_design) %in% colnames(D))
  V <- large_design[, added, drop = FALSE]
  small_label <- s$layout$ps_names[small]; large_label <- s$layout$ps_names[large]
  baseline <- drop(D %*% s$parameter[s$blocks$ps[[small]]]) + s$inputs$ps_models[[small]]$offset
  p <- plogis(baseline)
  information <- crossprod(D, D * (w * p * (1 - p))) / n
  cross_information <- crossprod(D, V * (w * p * (1 - p))) / n
  information_map <- solve(information, cross_information)
  delta <- drop(V %*% h - D %*% information_map %*% h) / sqrt(n)
  raw <- s$raw_scores
  raw[, small_label] <- p
  raw[, large_label] <- plogis(baseline + delta)
  centered <- sweep(raw, 2L, colMeans(raw), "-")
  score <- sweep(centered, 2L, sqrt(colMeans(centered^2)), "/")
  phi <- -t(solve(s$jacobian, t(s$estimating_equations)))
  phi <- sweep(phi, 2L, colMeans(phi), "-")
  H <- phi[, s$blocks$ps[[large]][added], drop = FALSE]
  Omega <- crossprod(H) / n
  arms <- if (object$estimand == "PATE") c("0", "1") else "0"
  mu <- lapply(arms, function(a) {
    spec <- s$inputs$pg_models[[a]][[1L]]
    drop(spec$design %*% s$parameter[s$blocks$pg[[a]][[1L]]]) + spec$offset
  }); names(mu) <- arms
  gamma <- mean(w * if (object$estimand == "PATT") z else 1)
  row <- if (object$estimand == "PATE") w * (mu[["1"]] - mu[["0"]] - object$estimate) / gamma else
    z * w * (s$inputs$Y - mu[["0"]] - object$estimate) / gamma
  slope <- rep(0, ncol(phi))
  all_neighbors <- vector("list", n); K <- matrix(0, n, 2L)
  bases <- path_scores <- list()
  for (arm in arms) {
    zz <- as.integer(arm); donors <- which(z == zz); queries <- which(z != zz)
    scores <- score[, s$layout$scores[[arm]], drop = FALSE]
    path_scores[[arm]] <- scores
    pg <- score[, s$layout$pg_names[[arm]][1L]]
    v <- sqrt(n) * (score[, large_label] - score[, small_label]) / sqrt(sum(h^2))
    t <- score[, small_label]
    basis <- cbind(1, t, pg, v, t^2, t * pg, t * v, pg^2, pg * v, v^2)
    design <- s$inputs$pg_models[[arm]][[1L]]$design
    Lq <- numeric(10L); Ld <- numeric(ncol(design))
    for (i in queries) {
      dist <- rowSums(sweep(scores[donors, , drop = FALSE], 2L, scores[i, ], "-")^2)
      selected <- if (is.null(graph)) donors[order(dist, donors)[seq_len(object$M)]] else
        graph$neighbors[[i]]
      all_neighbors[[i]] <- selected
      share <- w[selected] / sum(w[selected])
      K[selected, zz + 1L] <- K[selected, zz + 1L] + w[i] * share
      Lq <- Lq + w[i] * (basis[i, ] - colSums(basis[selected, , drop = FALSE] * share))
      Ld <- Ld + w[i] * (design[i, ] - colSums(design[selected, , drop = FALSE] * share))
    }
    Lq <- Lq * (2 * zz - 1) / (n * gamma)
    Ld <- Ld * (2 * zz - 1) / (n * gamma)
    Gram <- crossprod(basis[donors, , drop = FALSE], basis[donors, , drop = FALSE] * w[donors]) / n
    beta <- solve(Gram, Lq)
    projection <- drop(basis %*% beta)
    Dcross <- crossprod(basis[donors, , drop = FALSE], design[donors, , drop = FALSE] * w[donors]) / n
    ix <- s$blocks$pg[[arm]][[1L]]
    slope[ix] <- Ld - drop(crossprod(Dcross, beta))
    residual <- s$inputs$Y - mu[[arm]]
    if (object$estimand == "PATE") row[donors] <- row[donors] +
      (2 * zz - 1) * (w[donors] + K[donors, zz + 1L]) * residual[donors] / gamma else
      row[donors] <- row[donors] - K[donors, zz + 1L] * residual[donors] / gamma
    row[donors] <- row[donors] + w[donors] * residual[donors] * projection[donors]
    bases[[arm]] <- basis
  }
  total <- row - mean(row) + drop(phi %*% slope)
  covariance <- drop(crossprod(H, total) / n)
  regression <- solve(Omega, covariance)
  list(neighbors = all_neighbors, incoming = K, rows = row, complete = total,
    slope = slope, Omega = Omega, covariance = covariance,
    path_scores = path_scores, information_map = information_map,
    mean = sum(h * regression), variance = mean((total - drop(H %*% regression))^2))
}

for (target in c("PATE", "PATT")) test_that(paste("nested quadratic inputs retain point and full root", target), {
    object <- wm_qc_fixture(target); original <- object
    control <- list(method = "nested_contrast", small_ps = "smaller", large_ps = "larger")
    scope <- if (target == "PATE") "full_x" else "patt"
    out <- wm_model_inference(object, scope, quadratic_control = control)
    expect_identical(object, original)
    for (f in c("fit", "nuisance", "estimate", "raw_estimate", "correction", "solver")) expect_identical(out[[f]], original[[f]])
    inputs <- out$inference$contrast_inputs
    expect_identical(out$inference$status, "model_quadratic_contrast_inputs_ready")
    expect_true(out$inference$numerically_available)
    expect_false(out$inference$sampling_inference_available)
    expect_false(out$inference$available)
    expect_true(all(is.na(out$conf.int)))
    expect_identical(inputs$contrast_parameters, c("ps1:added", "ps1:extra"))
    expect_identical(inputs$dimension, 2L)
    expect_identical(inputs$full_parameters, object$nuisance$parameter_names)
    expect_identical(out$inference_inputs$quadratic_contrast, inputs)
    expect_identical(inputs$matching_dimensions, object$nuisance$layout$matching_dimensions)
    h <- c(.65, -.35)
    expected <- wm_qc_manual(object, h)
    evaluated <- inputs$evaluate(h)
    expect_true(evaluated$available)
    expect_identical(evaluated$graph$neighbors, expected$neighbors)
    expect_equal(evaluated$incoming, expected$incoming, tolerance = 1e-12, ignore_attr = TRUE)
    expect_equal(evaluated$base_projection_rows, expected$rows, tolerance = 1e-10)
    expect_equal(evaluated$complete_rows, expected$complete, tolerance = 1e-10)
    expect_equal(unname(evaluated$sensitivity), expected$slope, tolerance = 1e-10)
    expect_equal(evaluated$conditional_mean, expected$mean, tolerance = 1e-10)
    expect_equal(evaluated$conditional_variance, expected$variance, tolerance = 1e-10)
    expect_gte(evaluated$conditional_variance, 0)
    expect_equal(evaluated$linear_root_variance -
      drop(crossprod(evaluated$contrast_cross_covariance,
        solve(inputs$contrast_covariance, evaluated$contrast_cross_covariance))),
      evaluated$conditional_variance, tolerance = 1e-10)
    expect_equal(mean(evaluated$complete_rows), 0, tolerance = 1e-12)
    expect_identical(names(evaluated$projections), if (target == "PATE") c("0", "1") else "0")
    expect_true(all(vapply(evaluated$projections, function(x) ncol(x$basis) == 10L, logical(1L))))
    # After retaining this actual fit, preparation/evaluation uses no fit solver.
    local_mocked_bindings(.wm_wdsm_fit_stack = function(...) stop("UNEXPECTED_FIT"),
      .wm_model_stack = function(...) stop("UNEXPECTED_FIT"),
      .wm_ns_ols = function(...) stop("UNEXPECTED_OUTCOME_FIT"), .package = "wdsmatch")
    again <- wm_model_inference(object, scope, quadratic_control = control)
    expect_true(again$inference$contrast_inputs$evaluate(-h)$available)
    # Reuse this one fitted object for the common all-draw API; no extra fit/DGP.
    retained <- again
    rng_before <- .Random.seed
    source <- again$inference
    expect_true(source$replication_available)
    expect_true(source$contrast_inputs$replication_available)
    sample <- wm_bootstrap(source, B = 4L, seed = 820, chunk_size = 1L)
    expect_identical(.Random.seed, rng_before)
    expect_identical(again, retained)
    expect_true(sample$available)
    expect_identical(sample$interval, "basic")
    expect_true(all(sample$draw_diagnostics$status == "complete"))
    expect_identical(sample$source_inference, source)
    expect_identical(sample$estimate, object$estimate)
    expect_false(sample$numerical_error_rate_verified)
    expect_false(sample$population_assumptions_verified)
    expect_false(sample$original_refit_bootstrap)
    expect_true(is.na(sample$conditional_variance))
    expect_true(is.na(sample$conditional_root_n_variance))
    # Independent finite-path algebra at every generated contrast, including any
    # outside the former radial cutoff. No empirical recentering of raw roots.
    expected_m <- expected_v <- numeric(sample$B)
    for (j in seq_len(sample$B)) {
      truth <- wm_qc_manual(object, sample$contrast_draws[, j])
      expected_m[j] <- truth$mean; expected_v[j] <- truth$variance
    }
    expected_roots <- expected_m + sqrt(expected_v) * sample$scalar_innovations
    expect_equal(sample$draw_diagnostics$conditional_mean, expected_m, tolerance = 1e-8)
    expect_equal(sample$draw_diagnostics$conditional_variance, expected_v, tolerance = 1e-8)
    expect_equal(sample$root_n_draws, expected_roots, tolerance = 1e-8)
    expect_equal(sample$draws, object$estimate + expected_roots / sqrt(object$n), tolerance = 1e-8)
    expected_normals <- (function() {
      previous <- .Random.seed
      on.exit(assign(".Random.seed", previous, .GlobalEnv), add = TRUE)
      set.seed(820)
      matrix(stats::rnorm((inputs$dimension + 1L) * 4L), inputs$dimension + 1L, 4L)
    })()
    expect_identical(unname(sample$standard_contrast_normals),
      expected_normals[seq_len(inputs$dimension), , drop = FALSE])
    expect_identical(sample$scalar_innovations, expected_normals[inputs$dimension + 1L, ])
    expect_identical(.Random.seed, rng_before)
    qrH <- inputs$contrast_factorization
    expect_equal(unname(sample$contrast_draws[qrH$pivot, , drop = FALSE]),
      unname(crossprod(qrH$R, sample$standard_contrast_normals)), tolerance = 1e-14)
    expect_equal(sample$variance, stats::var(sample$root_n_draws) / object$n, tolerance = 1e-14)
    expect_equal(unname(sample$conf.int), object$estimate -
      as.numeric(stats::quantile(sample$root_n_draws, c(.975, .025))) / sqrt(object$n), tolerance = 1e-14)
    expect_identical(sample$variance_divisor, "B-1")
    other_chunk <- wm_bootstrap(source, B = 4L, seed = 820, chunk_size = 7L)
    for (name in c("root_n_draws", "draws", "contrast_draws", "standard_contrast_normals",
      "scalar_innovations", "draw_diagnostics", "variance", "conf.int"))
      expect_identical(other_chunk[[name]], sample[[name]])
    expect_identical(.Random.seed, rng_before)
    one <- wm_bootstrap(source, B = 1L, seed = 820, interval = "none")
    expect_identical(one$root_n_draws, sample$root_n_draws[1L])
    expect_identical(one$status, "insufficient_replicates")
    expect_true(one$draws_available)
    expect_false(one$available)
    expect_true(all(is.na(c(one$variance, one$se, one$monte_carlo_variance))))
    expect_null(one$conf.int)
    absent_seed <- function() {
      previous <- .Random.seed
      on.exit(assign(".Random.seed", previous, .GlobalEnv), add = TRUE)
      rm(".Random.seed", envir = .GlobalEnv)
      answer <- wm_bootstrap(source, B = 1L, seed = 820, interval = "none")
      expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      answer
    }
    expect_identical(absent_seed()$root_n_draws, one$root_n_draws)
    expect_identical(.Random.seed, rng_before)
    # Explicitly unsupported routes and mutable inputs fail before RNG.
    expect_error(wm_bootstrap(source, B = 2L, seed = 820, interval = "normal"), "arg")
    expect_error(wm_bootstrap(source, B = 2L, seed = 820, counts = matrix(1, object$n, 2)), "rejects counts")
    expect_error(wm_bootstrap(source, B = 2L, seed = 820, method = "fixed_reuse"), "does not support")
    expect_error(wm_bootstrap(source, B = 2L, seed = 820, refit = function(...) stop("UNEXPECTED_REFIT")), "does not support")
    bad <- source; bad$fit$weights[1L] <- 2 * bad$fit$weights[1L]
    expect_error(wm_bootstrap(bad, B = 2L, seed = 820), "original quadratic fit binding")
    bad <- source; bad$contrast_inputs$context_binding$phi[1L, 1L] <- 1 + bad$contrast_inputs$context_binding$phi[1L, 1L]
    expect_error(wm_bootstrap(bad, B = 2L, seed = 820), "context binding changed")
    bad <- source; bad$contrast_inputs$contrast_factorization$R[1L, 1L] <- 2 * bad$contrast_inputs$contrast_factorization$R[1L, 1L]
    expect_error(wm_bootstrap(bad, B = 2L, seed = 820), "QR binding changed")
    bad <- source; bad$contrast_inputs$full_score_covariance[1L, 1L] <- 1 + bad$contrast_inputs$full_score_covariance[1L, 1L]
    expect_error(wm_bootstrap(bad, B = 2L, seed = 820), "full covariance")
    limited <- wm_model_inference(object, scope, quadratic_control =
      c(control, list(maximum_draw_entries = 1)))
    expect_error(wm_bootstrap(limited$inference, B = 2L, seed = 820), "budget exceeded before RNG")
    expect_error(wm_model_inference(object, scope, quadratic_control =
      c(control, list(completion = "zero_guard_cutoff"))), "requires method")
    expect_identical(.Random.seed, rng_before)
    # A failed slot, even with obsolete completion hints, is never zero-filled.
    # The next slot is evaluated and uses the original ordered innovations.
    fail_two <- function() {
      calls <- 0L; actual <- wdsmatch:::.wm_qc_evaluate
      local_mocked_bindings(.wm_qc_evaluate = function(context, contrast) {
        calls <<- calls + 1L
        if (calls == 2L) return(list(available = FALSE, failure_stage = "zero_contrast",
          unavailable_reason = "Forced numerical failure", completion_allowed = TRUE,
          completion_reason = "exact_zero_contrast"))
        if (calls == 3L) stop("Forced arithmetic error")
        actual(context, contrast)
      }, .package = "wdsmatch")
      answer <- wm_bootstrap(source, B = 4L, seed = 820, chunk_size = 2L)
      expect_identical(calls, 4L)
      answer
    }
    failed <- fail_two()
    expect_false(failed$available)
    expect_false(failed$sampling_inference_available)
    expect_identical(failed$status, "numerical_inference_unavailable")
    expect_true(all(is.na(failed$root_n_draws[2:3])))
    expect_identical(failed$root_n_draws[c(1, 4)], sample$root_n_draws[c(1, 4)])
    expect_true(all(is.na(c(failed$variance, failed$se, failed$root_n_variance,
      failed$monte_carlo_variance, failed$monte_carlo_root_n_variance, failed$conf.int))))
    for (name in c("contrast_draws", "standard_contrast_normals", "scalar_innovations"))
      expect_identical(failed[[name]], sample[[name]])
    expect_identical(failed$estimate, source$estimate)
    expect_identical(.Random.seed, rng_before)
    expect_false(inputs$evaluate(c(0, 0))$available)
    expect_identical(inputs$evaluate(c(0, 0))$failure_stage, "zero_contrast")
    # The previous radial cutoff is not an operational restriction.
    far <- inputs$evaluate(c(sqrt(log(object$n + 1)) + 1, 0))
    expect_true(far$available)
    expect_false(far$sampling_inference_available)
    expect_false(far$numerical_error_rate_verified)
    expect_error(wm_model_inference(object, scope,
      quadratic_control = c(control, list(guard_exponent = 2))), "no longer supported")
    expect_error(wm_model_inference(object, scope,
      quadratic_control = c(control, list(guard_exponent = NULL))), "no longer supported")
    expect_error(inputs$evaluate(stats::setNames(h, c("wrong1", "wrong2"))), "complete added-parameter order")
    bad <- object; bad$nuisance$inputs$ps_models[[1L]]$offset[1L] <- 1
    expect_error(wm_model_inference(bad, scope, quadratic_control = control), "same fixed known offset")
    bad <- object; bad$nuisance$inputs$ps_models[[1L]]$weighting <- "unit"
    expect_error(wm_model_inference(bad, scope, quadratic_control = control), "probability weights")
    bad <- object; bad$nuisance$nuisance_influence[1L, 1L] <- bad$nuisance$nuisance_influence[1L, 1L] + 1
    expect_error(wm_model_inference(bad, scope, quadratic_control = control), "equation/influence binding")
    expect_identical(wm_model_inference(object, scope, quadratic_control = control,
      max_influence_elements = 1)$inference$failure_stage, "quadratic_inputs")
    expect_error(wm_model_inference(object, scope, quadratic_control = c(control, list(ridge = 1))), "requires method")
})

# Literal retained arrays for arithmetic-boundary checks only. No statistical
# model or population claim is associated with this fixture; no fit is run.
wm_qc_numeric_fixture <- function() {
  x <- as.matrix(expand.grid(x1 = seq(-1, 1, length.out = 5L),
    x2 = seq(-.9, .9, length.out = 5L), x3 = seq(-.8, .8, length.out = 5L)))
  n <- nrow(x); z <- rep(0:1, length.out = n)
  D <- cbind(intercept = 1, x1 = x[, 1])
  L <- cbind(D, added = x[, 3]); G <- cbind(intercept = 1, x2 = x[, 2])
  theta <- c(.1, .3, .1, .3, .15, 1, .8)
  raw <- cbind(ps1 = plogis(drop(D %*% theta[1:2])),
    ps2 = plogis(drop(L %*% theta[3:5])), pg0 = drop(G %*% theta[6:7]))
  center <- colMeans(raw); variance <- colMeans(sweep(raw, 2L, center)^2)
  theta <- c(theta, center, variance)
  parameters <- c("ps1:intercept", "ps1:x1", "ps2:intercept", "ps2:x1",
    "ps2:added", "pg0:intercept", "pg0:x2", paste0("center:", colnames(raw)),
    paste0("variance:", colnames(raw)))
  names(theta) <- parameters; p <- length(parameters)
  phi <- matrix(0, n, p, dimnames = list(NULL, parameters))
  phi[, 5L] <- x[, 3]
  J <- -diag(p); dimnames(J) <- list(parameters, parameters)
  w <- exp(.2 * z + .1 * x[, 1]); y <- raw[, "pg0"] + .2 * z
  spec <- function(design) list(design = design, offset = numeric(n), weighting = "probability")
  inputs <- list(Y = y, Z = z, weights = w,
    ps_models = list(small = spec(D), large = spec(L)),
    pg_models = list("0" = list(control = spec(G))))
  layout <- list(matching_dimensions = c("0" = 3L), ps_names = c("ps1", "ps2"),
    pg_names = list("0" = "pg0"), scores = list("0" = c("ps1", "ps2", "pg0")))
  stack <- structure(list(n = n, score_only = TRUE, parameter_names = parameters,
    parameter = theta, inputs = inputs, layout = layout, multiplicity = rep(1, n),
    empirical_probability = rep(1 / n, n), estimand = "PATT",
    nuisance_influence = phi, estimating_equations = phi, jacobian = J,
    weight_derivative = phi * 0, raw_scores = raw,
    blocks = list(ps = list(1:2, 3:5), pg = list("0" = list(6:7)),
      center = 8:10, variance = 11:13)), class = ".wm_model_nuisance_stack")
  score <- sweep(sweep(raw, 2L, center), 2L, sqrt(variance), "/")
  list(n = n, M = 3L, estimand = "PATT", estimate = 0,
    nuisance = stack, fit = list(n = n, M = 3L, estimand = "PATT", estimate = 0,
    weights = w, data = list(Y = y, Z = z), graph = list(scores0 = unname(score)),
    predictions = list(mean0 = raw[, "pg0"])),
    correction_fit = list(method = "quadratic_qr", arms = list("0" = list(mean = raw[, "pg0"]))))
}

test_that("quadratic numerical covariance gates retain representable full moments", {
  object <- wm_qc_numeric_fixture(); n <- object$fit$n
  control <- list(method = "nested_contrast", small_ps = "small", large_ps = "large")
  signal <- rep(c(-1, 1), length.out = n)
  object$nuisance$nuisance_influence[, 1L] <- 1e154 * signal
  object$nuisance$estimating_equations <- object$nuisance$nuisance_influence
  centered <- signal - mean(signal)
  # The old sum-before-division overflows, despite a representable covariance.
  expect_true(is.infinite(drop(crossprod(1e154 * centered)) / n))
  local_mocked_bindings(.wm_wdsm_fit_stack = function(...) stop("UNEXPECTED_FIT"),
    wm_match = function(...) stop("UNEXPECTED_GRAPH"), .package = "wdsmatch")
  answer <- wdsmatch:::.wm_qc_prepare(object, control, 5e7)
  expect_true(answer$available)
  expect_true(all(is.finite(answer$full_score_covariance)))
  expect_equal(answer$full_score_covariance[1L, 1L] / 1e308,
    mean(centered^2), tolerance = 1e-13)
  object$nuisance$nuisance_influence[, 1L] <- 1e155 * signal
  object$nuisance$estimating_equations <- object$nuisance$nuisance_influence
  failed <- wdsmatch:::.wm_qc_prepare(object, control, 5e7)
  expect_false(failed$available)
  expect_identical(failed$failure_stage, "score_covariance")
})

test_that("quadratic numerical QR failures retain all columns and return unavailable", {
  factor <- wdsmatch:::.wm_qc_qr
  expect_equal(crossprod(factor(diag(c(1, 2)))$R), diag(c(4, 1)))
  expect_error(factor(diag(c(1, Inf))), "nonfinite")
  expect_error(factor(cbind(1:4, 1:4)), "unresolved")
  # Inject a numerical QR failure into a copy of the helper, without
  # modifying any package binding or running another statistical fixture.
  environment(factor) <- list2env(list(qr = function(...) stop("LAPACK failure")),
    parent = environment(factor))
  expect_error(factor(diag(2)), "LAPACK failure")
  object <- wm_qc_numeric_fixture()
  control <- list(method = "nested_contrast", small_ps = "small", large_ps = "large")
  prepared <- wdsmatch:::.wm_qc_prepare(object, control, 5e7)
  context <- environment(prepared$evaluate)$retained
  context$V <- context$V * 1e150
  context$information_map <- context$information_map * 1e150
  failed <- wdsmatch:::.wm_qc_evaluate(context, .6e-150)
  expect_false(failed$available)
  expect_identical(failed$failure_stage, "projection_qr")
  expect_identical(failed$contrast, .6e-150)
})

test_that("quadratic numerical return gates include the unconditional row moment", {
  object <- wm_qc_numeric_fixture()
  control <- list(method = "nested_contrast", small_ps = "small", large_ps = "large")
  prepared <- wdsmatch:::.wm_qc_prepare(object, control, 5e7)
  context <- environment(prepared$evaluate)$retained
  # Isolate a large row component that is removed by conditional projection.
  # Finite conditional residual moments do not certify its full second moment.
  H <- numeric(context$n)
  H[context$Z == 1L] <- rep(c(-1, 1), length.out = sum(context$Z == 1L))
  H <- matrix(H, ncol = 1L)
  context$mu[["0"]] <- numeric(context$n)
  context$derivatives[["0"]][] <- 0
  context$phi[] <- 0
  context$H <- H; context$Omega <- crossprod(H, H / context$n)
  context$contrast_qr <- wdsmatch:::.wm_qc_qr(H / sqrt(context$n))
  context$Y <- as.vector(H) * 1e153 / context$w
  accepted <- wdsmatch:::.wm_qc_evaluate(context, .6)
  expect_true(accepted$available)
  expect_true(is.finite(accepted$linear_root_variance))
  expect_true(is.finite(accepted$conditional_variance))
  expect_lt(accepted$conditional_variance / accepted$linear_root_variance, 1e-25)
  context$Y <- as.vector(H) * 1e155 / context$w
  failed <- wdsmatch:::.wm_qc_evaluate(context, .6)
  expect_false(failed$available)
  expect_identical(failed$failure_stage, "arithmetic")
})

test_that("untruncated QR inputs agree with independent finite rows at nonzero residuals", {
  object <- wm_qc_numeric_fixture(); original <- object
  # Modify only existing literal outcome records, without a fit, RNG or DGP.
  # This fixture's PG derivative is in the quadratic span, so its zero b is
  # an algebraic reduction; nonzero-b coverage stays in the existing fit tests.
  residual <- .07 * sin(seq_len(object$n)) + .03 * cos(2 * seq_len(object$n))
  object$nuisance$inputs$Y <- object$nuisance$inputs$Y + residual
  object$fit$data$Y <- object$nuisance$inputs$Y
  retained <- object
  control <- list(method = "nested_contrast", small_ps = "small", large_ps = "large")
  local_mocked_bindings(.wm_wdsm_fit_stack = function(...) stop("UNEXPECTED_FIT"),
    .wm_model_stack = function(...) stop("UNEXPECTED_FIT"),
    .wm_ns_ols = function(...) stop("UNEXPECTED_FIT"), .package = "wdsmatch")
  prepared <- wdsmatch:::.wm_qc_prepare(object, control, 5e7)
  value <- prepared$evaluate(.6)
  expect_true(value$available)
  # This literal grid has ties. Verify its path and information map separately,
  # then condition the independent projection algebra on the actual graph;
  # independently rounded near-tie distances need not select identical donors.
  expected <- wm_qc_manual(object, .6, small = 1L, large = 2L, graph = value$graph)
  expect_identical(object, retained)
  expect_identical(value$estimate, original$fit$estimate)
  expect_identical(dim(value$graph$scores0), dim(expected$path_scores[["0"]]))
  expect_lt(max(abs(value$graph$scores0 - expected$path_scores[["0"]])), 1e-13)
  expect_equal(prepared$contrast_information_map, expected$information_map,
    tolerance = 1e-13, ignore_attr = TRUE)
  expect_equal(value$incoming, expected$incoming, tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(value$base_projection_rows, expected$rows, tolerance = 1e-10)
  expect_equal(value$complete_rows, expected$complete, tolerance = 1e-10)
  expect_equal(unname(value$sensitivity), expected$slope, tolerance = 1e-10)
  expect_equal(unname(value$contrast_cross_covariance), unname(expected$covariance),
    tolerance = 1e-10)
  expect_equal(value$conditional_mean, expected$mean, tolerance = 1e-10)
  expect_equal(value$conditional_variance, expected$variance, tolerance = 1e-10)
  donor <- which(object$fit$data$Z == 0L)
  expect_gt(sum(abs(value$projections[["0"]]$residual[donor])), 0)
  expect_gt(sum(abs(value$projections[["0"]]$prediction[donor])), 0)
  expect_identical(ncol(value$projections[["0"]]$basis), 10L)
  expect_false(value$sampling_inference_available)
  expect_false(value$numerical_error_rate_verified)
  far <- prepared$evaluate(sqrt(log(object$n + 1)) + 1)
  expect_true(far$available)
  expect_false(far$sampling_inference_available)
  expect_error(wdsmatch:::.wm_qc_prepare(object,
    c(control, list(guard_exponent = NULL)), 5e7), "no longer supported")
  # A representable small contrast covariance is not compared with n^-2.
  small_covariance <- object
  small_covariance$nuisance$nuisance_influence[, 5L] <-
    small_covariance$nuisance$nuisance_influence[, 5L] * 1e-10
  small_covariance$nuisance$estimating_equations <- small_covariance$nuisance$nuisance_influence
  small_inputs <- wdsmatch:::.wm_qc_prepare(small_covariance, control, 5e7)
  expect_true(small_inputs$available)
  expect_gt(small_inputs$contrast_covariance[1L, 1L], 0)
  expect_lt(small_inputs$contrast_covariance[1L, 1L], object$n^-2)
})

test_that("tiny logistic contrasts use their finite high-precision value", {
  skip_if_not_installed("Rmpfr")
  object <- wm_qc_numeric_fixture()
  control <- list(method = "nested_contrast", small_ps = "small", large_ps = "large")
  prepared <- wdsmatch:::.wm_qc_prepare(object, control, 5e7)
  context <- environment(prepared$evaluate)$retained
  a <- drop(context$D %*% context$alpha) + context$offset
  direction <- drop(context$V - context$D %*% context$information_map)
  # The reference evaluates the exact finite expit difference at 256 bits.
  # It is not the derivative at zero and does not assert a general error bound.
  for (sign in c(-1, 1)) {
    step <- 1e-25 / sqrt(context$n)
    value <- wdsmatch:::.wm_qc_logistic_contrast(a, sign * direction, step)
    aa <- Rmpfr::mpfr(a, 256L)
    dd <- Rmpfr::mpfr(sign * direction, 256L)
    tt <- Rmpfr::mpfr(step, 256L)
    exact <- Rmpfr::asNumeric((1 / (1 + exp(-(aa + tt * dd))) -
      1 / (1 + exp(-aa))) / tt)
    expect_true(all(is.finite(value$divided_probability)))
    expect_equal(value$divided_probability, exact, tolerance = 2e-13)
    expect_equal(value$large_probability, plogis(a + step * sign * direction))
  }
  expect_error(wdsmatch:::.wm_qc_logistic_contrast(a, direction,
    .Machine$double.xmin * .Machine$double.eps), "unrepresentable")
})
