# Private independent Gaussian-Wick arithmetic reference.
# No coefficient/log-determinant recurrence, numerical integration, RNG or fits.
# shift_factor has 4 rows and any q >= 0 columns; its law is
# (Delta_0 H, Delta_1 H) = shift_factor %*% Z, Z ~ N_q(0,I).
.cp_wick_geometry_reference <- function(A0, A1, lambda,
                                        shift_factor = matrix(0, 4L, 0L)) {
  stopifnot(is.matrix(A0), is.matrix(A1),
            identical(dim(A0), c(2L, 2L)),
            identical(dim(A1), c(2L, 2L)),
            is.numeric(lambda), length(lambda) == 2L,
            all(is.finite(lambda)), all(lambda > 0),
            all(is.finite(A0)), all(is.finite(A1)),
            is.matrix(shift_factor), nrow(shift_factor) == 4L,
            all(is.finite(shift_factor)))
  E <- rbind(sqrt(lambda[1L]) * A0, sqrt(lambda[2L]) * A1)
  S <- crossprod(E)
  Pi <- E %*% solve(S, t(E))
  P <- diag(4L) - Pi
  R <- sweep(shift_factor, 1L, rep(sqrt(lambda), each = 2L), "*")
  q <- ncol(R)
  logdetS <- determinant(S, logarithm = TRUE)
  stopifnot(logdetS$sign == 1)
  if (q == 0L) {
    logdetF <- 0
    C <- Pi / 2
  } else {
    F <- diag(q) + 2 * crossprod(R, P %*% R)
    logdet <- determinant(F, logarithm = TRUE)
    stopifnot(logdet$sign == 1)
    logdetF <- as.numeric(logdet$modulus)
    PR <- P %*% R
    C <- Pi / 2 + PR %*% solve(F, t(PR))
  }
  # Exact covariance is symmetric; this only symmetrizes floating arithmetic.
  C <- (C + t(C)) / 2
  normalizer <- exp(log(pi) - (as.numeric(logdetS$modulus) + logdetF) / 2)
  stopifnot(all(is.finite(C)), is.finite(normalizer), normalizer > 0)
  list(normalizer = normalizer, covariance = C)
}

.cp_wick_moment_reference <- function(C, alpha) {
  stopifnot(is.matrix(C), identical(dim(C), c(4L, 4L)),
            all(is.finite(C)), length(alpha) == 4L,
            all(is.finite(alpha)), all(alpha >= 0),
            all(alpha == floor(alpha)))
  memo <- new.env(parent = emptyenv(), hash = TRUE)
  assign("0,0,0,0", 1, envir = memo)
  moment <- function(a) {
    key <- paste(a, collapse = ",")
    if (exists(key, envir = memo, inherits = FALSE))
      return(get(key, envir = memo, inherits = FALSE))
    if (sum(a) %% 2L == 1L) {
      ans <- 0
    } else {
      first <- which(a > 0L)[1L]
      beta <- a
      beta[first] <- beta[first] - 1L
      ans <- 0
      for (j in seq_len(4L)) {
        if (beta[j] > 0L) {
          remaining <- beta
          remaining[j] <- remaining[j] - 1L
          ans <- ans + beta[j] * C[first, j] * moment(remaining)
        }
      }
    }
    assign(key, ans, envir = memo)
    ans
  }
  moment(alpha)
}

.cp_m2_trace_reference <- function(normalizer, C) {
  a <- sum(diag(C)[1:2])
  b <- sum(diag(C)[3:4])
  ab <- a * b + 2 * sum(C[1:2, 3:4, drop = FALSE]^2)
  normalizer * (1 + a + b + ab)
}

.cp_wick_average_reference <- function(M, A0, A1, lambda,
                                       shift_factor = matrix(0, 4L, 0L)) {
  stopifnot(is.numeric(M), length(M) == 1L, is.finite(M),
            M >= 1, M == floor(M), M <= .Machine$integer.max)
  M <- as.integer(M)
  geometry <- .cp_wick_geometry_reference(A0, A1, lambda, shift_factor)
  C <- geometry$covariance
  polynomial_mean <- 0
  for (k in seq.int(0L, M - 1L)) {
    for (l in seq.int(0L, M - 1L)) {
      for (a in seq.int(0L, k)) {
        for (b in seq.int(0L, l)) {
          alpha <- 2 * c(a, k - a, b, l - b)
          polynomial_mean <- polynomial_mean +
            choose(k, a) * choose(l, b) /
            (factorial(k) * factorial(l)) *
            .cp_wick_moment_reference(C, alpha)
        }
      }
    }
  }
  value <- geometry$normalizer * polynomial_mean
  stopifnot(is.finite(value))
  list(value = value, normalizer = geometry$normalizer,
       covariance = C,
       M1_closed = geometry$normalizer,
       M2_closed = .cp_m2_trace_reference(geometry$normalizer, C))
}

# Exact finite-B identity check under SAME valid operative inputs.
# regression must obey Sigma %*% regression = C and q=C' regression.
# support_residual is disclosed rather than inferred from arbitrary arguments.
.cp_rb_raw_identity_reference <- function(vraw_draws, exact_vraw_mean,
                                          b, C, Sigma, regression) {
  stopifnot(length(vraw_draws) > 0L, all(is.finite(vraw_draws)),
            length(exact_vraw_mean) == 1L, is.finite(exact_vraw_mean),
            length(b) == length(C), length(C) == length(regression),
            is.matrix(Sigma), nrow(Sigma) == length(b),
            ncol(Sigma) == length(b))
  q <- sum(C * regression)
  beta <- b + regression
  linear_variance <- as.numeric(crossprod(beta, Sigma %*% beta))
  cross_and_fit <- 2 * sum(b * C) +
    as.numeric(crossprod(b, Sigma %*% b))
  adjustment <- pmax(q - vraw_draws, 0)
  rb <- mean(pmax(vraw_draws - q, 0)) + linear_variance
  aligned_raw <- exact_vraw_mean + cross_and_fit
  centering <- mean(vraw_draws) - exact_vraw_mean
  rhs <- centering + mean(adjustment)
  list(rb_variance = rb, aligned_raw_variance = aligned_raw,
       difference = rb - aligned_raw, identity_rhs = rhs,
       arithmetic_residual = rb - aligned_raw - rhs,
       sample_centering = centering, mean_cone_adjustment = mean(adjustment),
       support_residual = as.numeric(Sigma %*% regression - C))
}

dv_fixture <- function(M = 1L, unit = FALSE) {
  angle <- (0:7) * pi / 4
  S <- rbind(.45 * cbind(cos(angle), sin(angle)),
    c(-.06, .03), c(.07, -.03), c(0, -.04), c(.02, .05))
  n <- nrow(S); Z <- c(rep(0:1, 4), 1, 0, 1, 0)
  W <- if (unit) rep(1, n) else 1 + seq_len(n) / 20
  influence <- cbind(a = 3 * (seq_len(n) %% 2 - .5), b = 2 * S[, 1])
  B0 <- B1 <- array(0, c(n, 2, 2), dimnames = list(NULL, NULL, c("a", "b")))
  B0[, 1, 1] <- S[, 1]; B0[, 2, 2] <- Z / 3
  B1[, 1, 2] <- S[, 2]; B1[, 2, 1] <- (1 - Z) / 4
  dat <- list(Y = (-1)^seq_len(n) + S[, 1], Z = Z, weights = W,
    scores0 = S, scores1 = S, M = M, estimand = "PATE",
    mean0 = .1 + S[, 1], mean1 = .3 + S[, 2], variance = FALSE)
  control <- list(tangents0 = B0, tangents1 = B1, auxiliary_rows = 1:8,
    bandwidth = 1, density_bounds = c(.01, 10), marginal_mass_floor = .001,
    projection_floor = .2, tangent_cap = 10, residual_cap = 10,
    support = "smooth_stack_threshold", quadrature_tolerance = 1e-8,
    kernel_error_tolerance = 1e-4, maximum_laplace_terms = 200000000L)
  list(data = dat, influence = influence, D0 = cbind(a = S[, 2], b = S[, 1]),
    D1 = cbind(a = S[, 1], b = S[, 2]), control = control)
}
dv_fit_inference <- function(a, fit = NULL) {
  if (is.null(fit)) fit <- do.call(wm_match, a$data)
  zero <- list(mode = "zero", basis = "fixed_map")
  wm_fitted_inference(fit, a$influence, a$D0, a$D1,
    transport0 = zero, transport1 = zero, covariance_scope = "critical_planar",
    critical_control = a$control)
}
dv_A_exhaustive <- function(w, donor, probability, M) {
  if (M == 1) return(1)
  grid <- expand.grid(rep(list(seq_along(donor)), M - 1), KEEP.OUT.ATTRS = FALSE)
  answer <- 0
  for (r in seq_len(nrow(grid))) {
    index <- as.integer(unlist(grid[r, , drop = FALSE], use.names = FALSE))
    answer <- answer + prod(probability[index]) * w / (w + sum(donor[index]))
  }
  answer
}
dv_row_reference <- function(a, neighbors) {
  dat <- a$data; W <- dat$weights; Z <- dat$Z; n <- length(W)
  epsilon <- dat$Y - ifelse(Z == 1, dat$mean1, dat$mean0)
  load <- contrast <- numeric(n); derivative <- numeric(2)
  for (i in seq_len(n)) {
    donor <- neighbors[[i]]; fraction <- W[donor] / sum(W[donor])
    load[donor] <- load[donor] + W[i] * fraction
    mu <- if (Z[i] == 1) dat$mean0 else dat$mean1
    D <- if (Z[i] == 1) a$D0 else a$D1
    contrast[i] <- (2 * Z[i] - 1) * (dat$Y[i] - mu[i] - sum(fraction * epsilon[donor]))
    derivative <- derivative + W[i] * (1 - 2 * Z[i]) *
      (D[i, ] - colSums(D[donor, , drop = FALSE] * fraction))
  }
  estimate <- sum(W * contrast) / sum(W)
  base <- (W * (dat$mean1 - dat$mean0 - estimate) +
    (2 * Z - 1) * (W + load) * epsilon) / mean(W)
  base <- base - mean(base)
  influence <- sweep(a$influence, 2, colMeans(a$influence), "-")
  list(estimate = estimate, base = base, diagonal = mean(base^2),
    C = as.vector(crossprod(influence, base)) / n,
    Sigma = crossprod(influence) / n, b = unname(derivative / sum(W)))
}

dv_numeric_control <- function() list(
  projection_floor = .01, maximum_pair_entries = 2000000L,
  maximum_analytic_matrix_products = 5000000L,
  maximum_analytic_coefficient_terms = 200000000L,
  minimum_analytic_reciprocal_condition = 1e-12,
  analytic_precision_tolerance = 1e-10)

dv_root <- function(Sigma) {
  e <- eigen(Sigma, symmetric = TRUE)
  e$vectors %*% diag(sqrt(pmax(e$values, 0)), nrow = nrow(Sigma)) %*% t(e$vectors)
}

# Complete outer reference from raw inputs, not cached candidate pairs/laws/E.
dv_average_R_reference <- function(a, Sigma) {
  dat <- a$data; W <- dat$weights / max(dat$weights)
  T <- cbind(dat$scores0, dat$scores1); aux <- a$control$auxiliary_rows
  ev <- setdiff(seq_along(W), aux); h <- a$control$bandwidth; p <- nrow(Sigma)
  epsilon <- dat$Y - ifelse(dat$Z == 1L, dat$mean1, dat$mean0)
  factor <- dv_root(Sigma); result <- 0
  for (i in ev[dat$Z[ev] == 1L]) {
    v <- sweep(T[aux, , drop = FALSE], 2L, T[i, ], "-") / h
    weight <- pmax(0, 1 - sqrt(rowSums(v^2))); mass <- sum(weight)
    if (mass == 0) next
    center <- colSums(v * weight) / mass
    centered <- sweep(v, 2L, center, "-")
    Q <- crossprod(centered, centered * weight) / mass
    eigenvalues <- eigen(Q, symmetric = TRUE)$values
    E <- eigen(Q, symmetric = TRUE)$vectors[, 1:2, drop = FALSE]
    projections <- c(min(svd(E[1:2, , drop = FALSE], nu = 0, nv = 0)$d),
      min(svd(E[3:4, , drop = FALSE], nu = 0, nv = 0)$d))
    if (eigenvalues[2L] < .075 || any(projections < a$control$projection_floor)) next
    density <- numeric(2L); laws <- list()
    for (z in 0:1) {
      S <- if (z == 0L) dat$scores0 else dat$scores1
      donors <- aux[dat$Z[aux] == z]
      k <- pmax(0, 1 - sqrt(rowSums(sweep(S[donors, , drop = FALSE], 2L,
        S[i, ], "-")^2)) / h) / (pi / 3)
      mass <- sum(k); raw <- mass / (length(aux) * h^2)
      density[z + 1L] <- min(a$control$density_bounds[2L],
        max(a$control$density_bounds[1L], raw))
      positive <- k > 0
      laws[[z + 1L]] <- if (mass < a$control$marginal_mass_floor * length(aux) * h^2 ||
          !any(positive)) list(donor = W[1L], probability = 1) else
        list(donor = W[donors[positive]], probability = k[positive] / mass)
    }
    for (j in ev[dat$Z[ev] == 0L]) {
      k <- max(0, 1 - sqrt(sum((T[j, ] - T[i, ])^2)) / h)
      if (!k) next
      A1 <- dv_A_exhaustive(W[i], laws[[2L]]$donor, laws[[2L]]$probability, dat$M)
      A0 <- dv_A_exhaustive(W[j], laws[[1L]]$donor, laws[[1L]]$probability, dat$M)
      D0 <- matrix(a$control$tangents0[j, , ], 2L, p) -
        matrix(a$control$tangents0[i, , ], 2L, p)
      D1 <- matrix(a$control$tangents1[j, , ], 2L, p) -
        matrix(a$control$tangents1[i, , ], 2L, p)
      G <- .cp_wick_average_reference(dat$M, E[1:2, , drop = FALSE],
        E[3:4, , drop = FALSE], pi * density, rbind(D0, D1) %*% factor)$value
      result <- result + k * W[i] * epsilon[i] * A1 * W[j] * epsilon[j] * A0 * G
    }
  }
  result / (length(ev) * (length(ev) - 1L) * h^2 * (pi / 3))
}

test_that("fixed four-shift Gaussian average matches independent M1/M2/M3 Wick route", {
  E <- qr.Q(qr(matrix(c(1, .2, .6, .1, .3, 1.1, -.2, .9), 4L, 2L)))
  D0 <- matrix(c(.3, -.4, .1, -.2, .5, .2), 2L, 3L)
  D1 <- matrix(c(-.2, .3, .4, .1, -.5, .6), 2L, 3L)
  L <- rbind(c(.9, .1), c(-.2, .6), c(.3, -.4))
  rotation <- matrix(c(cos(.4), sin(.4), -sin(.4), cos(.4)), 2L, 2L)
  for (M in 1:3) for (Sigma in list(tcrossprod(L), diag(c(1, .4, .2)), matrix(0, 3L, 3L))) {
    reference <- .cp_wick_average_reference(M, E[1:2, ], E[3:4, ],
      pi * c(.3, .8), rbind(D0, D1) %*% dv_root(Sigma))
    got <- wdsmatch:::.wm_cp_G_average(E, c(.3, .8), D0, D1, Sigma, M,
      dv_numeric_control())
    expect_equal(got$value, reference$value, tolerance = 1e-10)
    if (M == 1L) expect_equal(got$value, reference$M1_closed, tolerance = 1e-11)
    if (M == 2L) expect_equal(got$value, reference$M2_closed, tolerance = 1e-10)
    expect_equal(got$Omega, rbind(D0, D1) %*% Sigma %*% t(rbind(D0, D1)), tolerance = 1e-13)
    expect_equal(got$omega_square_root %*% got$omega_square_root, got$Omega, tolerance = 1e-11)
    expect_equal(wdsmatch:::.wm_cp_G_average(E %*% rotation, c(.3, .8),
      D0, D1, Sigma, M, dv_numeric_control())$value, reference$value, tolerance = 1e-10)
    expect_identical(got$factor_columns, 4L); expect_identical(got$SPD_dimension, 6L)
    expect_identical(got$exact_arithmetic_geometric_error, 0)
    expect_false(got$roundoff_included)
  }
  # The reference uses original factor columns, including redundant zero columns.
  ref2 <- .cp_wick_average_reference(3L, E[1:2, ], E[3:4, ], pi * c(.3, .8),
    rbind(D0, D1) %*% L)$value
  ref5 <- .cp_wick_average_reference(3L, E[1:2, ], E[3:4, ], pi * c(.3, .8),
    cbind(rbind(D0, D1) %*% L, matrix(0, 4L, 3L)))$value
  expect_equal(ref2, ref5, tolerance = 1e-11)
  Sigma <- tcrossprod(L); Delta <- rbind(D0, D1); Omega <- Delta %*% Sigma %*% t(Delta)
  independent_arms <- matrix(0, 4L, 4L)
  independent_arms[1:2, 1:2] <- Omega[1:2, 1:2]
  independent_arms[3:4, 3:4] <- Omega[3:4, 3:4]
  wrong <- .cp_wick_average_reference(2L, E[1:2, ], E[3:4, ], pi * c(.3, .8),
    dv_root(independent_arms))$value
  correct <- .cp_wick_average_reference(2L, E[1:2, ], E[3:4, ], pi * c(.3, .8),
    Delta %*% L)$value
  expect_gt(abs(correct - wrong), 1e-5)
  # Only odd total degree vanishes; correlated individual odd exponents remain.
  G <- matrix(c(1, .4, 0, 0, .4, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1), 4L)
  expect_equal(.cp_wick_moment_reference(G, c(1, 1, 0, 0)), .4)
})

test_that("raw-only complete variance preserves all donor/root factors and full blocks", {
  for (M in 1:3) for (unit in c(FALSE, TRUE)) {
    a <- dv_fixture(M, unit); a$control$support <- NULL
    a$control$direct_variance <- "raw_only"
    # No unused Hermite rule/order/geometric workload is required by raw-only.
    a$control$maximum_hermite_order <- 1L
    a$control$maximum_geometric_evaluations <- 1L
    fit <- do.call(wm_match, a$data); frozen <- fit
    set.seed(7201); constructor_seed <- .Random.seed
    out <- dv_fit_inference(a, fit); raw <- out$direct_variance$raw_empirical
    expect_identical(.Random.seed, constructor_seed)
    row <- dv_row_reference(a, fit$graph$neighbors)
    R <- dv_average_R_reference(a, row$Sigma)
    expected <- row$diagonal - 2 * R / fit$gamma^2 +
      2 * sum(row$b * row$C) + sum(row$b * as.vector(row$Sigma %*% row$b))
    expect_true(raw$numerically_available)
    expect_lte(abs(raw$average_reciprocal - R), raw$reciprocal_quadrature_error_bound + 1e-10)
    expect_lte(abs(raw$raw_root_n_variance - expected), raw$root_n_variance_quadrature_error_bound + 1e-10)
    expect_equal(raw$raw_variance, raw$raw_root_n_variance / out$n)
    expect_equal(unname(raw$C), row$C, tolerance = 1e-12)
    expect_equal(raw$Sigma, row$Sigma, tolerance = 1e-12)
    expect_equal(unname(raw$total_sensitivity), row$b, tolerance = 1e-12)
    expect_identical(out$fit, frozen); expect_identical(fit, frozen)
    expect_null(out$critical$support); expect_null(out$critical$kernel$rule)
    expect_true(is.na(out$critical$kernel$diagnostics$hermite_order))
    expect_false(out$conditional_inference_available); expect_false(out$available)
    expect_true(is.na(out$variance)); expect_true(all(is.na(out$conf.int)))
    expect_match(raw$covariance_convention, "centered")
    expect_false(raw$Gaussian_interval_authorized)
    expect_true(is.na(raw$exact_operative_cone_root_n_variance))
    scaled <- a; scaled$data$weights <- 29 * scaled$data$weights
    other <- dv_fit_inference(scaled)$direct_variance$raw_empirical
    expect_equal(other$raw_root_n_variance, raw$raw_root_n_variance, tolerance = 1e-11)
    if (unit) expect_equal(raw$root_n_variance_quadrature_error_bound, 0)
    set.seed(2703); state <- .Random.seed
    expect_error(wm_bootstrap(out, B = 3L, seed = 17L), "does not authorize")
    expect_identical(.Random.seed, state)
  }
})

test_that("direct raw finite convention is centered and no support rate is imposed", {
  a <- dv_fixture(); a$influence[, 2L] <- 1e-5 * a$influence[, 2L]
  a$influence <- sweep(a$influence, 2L, c(2, -3), "+")
  a$control$support <- NULL; a$control$direct_variance <- "raw_only"
  out <- dv_fit_inference(a); raw <- out$direct_variance$raw_empirical
  centered <- sweep(a$influence, 2L, colMeans(a$influence), "-")
  expected <- crossprod(centered) / nrow(centered)
  expect_equal(raw$Sigma, expected, tolerance = 1e-13)
  expect_gt(max(abs(raw$Sigma - crossprod(a$influence) / nrow(a$influence))), 1)
  expect_null(out$critical$support)
  expect_true(raw$numerically_available)
  covariance <- wdsmatch:::.wm_cp_direct_psd(diag(c(1, 1e-15)), "tiny-positive", 1e-10)
  expect_equal(covariance$root[2L, 2L], sqrt(1e-15), tolerance = 1e-15)
  expect_identical(covariance$positive_eigenvalues_discarded, 0L)
  expect_null(covariance$support_threshold)
  # A full p=5 joint stack is retained; pair integration still uses four shifts.
  a <- dv_fixture(); n <- nrow(a$influence)
  a$influence <- cbind(a$influence, c = seq_len(n) / n,
    d = a$data$scores0[, 2L], e = a$data$Z)
  parameters <- colnames(a$influence)
  a$D0 <- cbind(a$D0, c = a$data$scores0[, 1L], d = a$data$Z, e = 1)
  a$D1 <- cbind(a$D1, c = a$data$scores0[, 2L], d = 1 - a$data$Z, e = 2)
  for (name in c("tangents0", "tangents1")) {
    old <- a$control[[name]]
    new <- array(0, c(n, 2L, 5L), dimnames = list(NULL, NULL, parameters))
    new[, , 1:2] <- old
    new[, 1L, 3L] <- a$data$scores0[, 1L]
    new[, 2L, 4L] <- a$data$scores0[, 2L]
    new[, 1L, 5L] <- a$data$Z / 4
    a$control[[name]] <- new
  }
  a$control$support <- NULL; a$control$direct_variance <- "raw_only"
  out <- dv_fit_inference(a); raw <- out$direct_variance$raw_empirical
  Sigma <- crossprod(sweep(a$influence, 2L, colMeans(a$influence), "-")) / n
  expect_identical(dim(raw$Sigma), c(5L, 5L))
  expect_identical(names(raw$total_sensitivity), parameters)
  R <- dv_average_R_reference(a, Sigma)
  expect_lte(abs(raw$average_reciprocal - R), raw$reciprocal_quadrature_error_bound + 1e-10)
})

test_that("mixture opt-in separates aligned raw variance from finite-B cone summaries", {
  a <- dv_fixture(3L); a$influence[, 2L] <- 10 * a$influence[, 2L]
  fit <- do.call(wm_match, a$data); ordinary <- dv_fit_inference(a, fit)
  a$control$direct_variance <- "with_mixture"
  out <- dv_fit_inference(a, fit); raw <- out$direct_variance$raw_empirical
  aligned <- out$direct_variance$support_aligned; support <- out$critical$support
  expect_true(raw$numerically_available); expect_true(aligned$numerically_available)
  expect_equal(aligned$Sigma, support$Sigma); expect_equal(aligned$C, support$C)
  reference_R <- dv_average_R_reference(a, support$Sigma)
  expect_lte(abs(aligned$average_reciprocal - reference_R),
    aligned$reciprocal_quadrature_error_bound + 1e-10)
  expected <- out$diagonal_root_n_variance - 2 * reference_R / fit$gamma^2 +
    2 * sum(out$total_sensitivity * support$C) +
    sum(out$total_sensitivity * as.vector(support$Sigma %*% out$total_sensitivity))
  expect_lte(abs(aligned$raw_root_n_variance - expected),
    aligned$root_n_variance_quadrature_error_bound + 1e-10)
  set.seed(993); before <- .Random.seed
  old_draws <- wm_bootstrap(ordinary, B = 7L, seed = 441L)
  new_draws <- wm_bootstrap(out, B = 7L, seed = 441L)
  old_draws$source_inference <- new_draws$source_inference <- NULL
  expect_identical(old_draws, new_draws); expect_identical(.Random.seed, before)
  identity <- .cp_rb_raw_identity_reference(new_draws$draw_diagnostics$raw_V0,
    out$diagonal_root_n_variance - 2 * aligned$average_reciprocal / fit$gamma^2,
    out$total_sensitivity, support$C, support$Sigma, support$regression)
  expect_equal(identity$rb_variance, new_draws$rao_blackwell_root_n_variance_estimate, tolerance = 1e-12)
  expect_equal(identity$aligned_raw_variance, aligned$raw_root_n_variance, tolerance = 1e-12)
  expect_equal(identity$arithmetic_residual, 0, tolerance = 1e-11)
  expect_equal(identity$support_residual, c(0, 0), tolerance = 1e-11)
  expect_true(is.na(new_draws$conditional_root_n_variance))
  expect_true(is.na(aligned$exact_operative_cone_root_n_variance))
})

test_that("signed raw variance and analytic failures are retained without rescue", {
  a <- dv_fixture(); a$control$support <- NULL; a$control$direct_variance <- "raw_only"
  out <- dv_fit_inference(a); kernel <- out$critical$kernel
  # Internal scalar arithmetic only: this does not modify a public fit or
  # assert that an artificial finite negative block satisfies population laws.
  calculation <- out; calculation$diagonal_root_n_variance <- -100
  signed <- wdsmatch:::.wm_cp_direct_one(calculation, kernel, out$raw_Sigma, out$raw_C, "literal")
  expect_true(signed$numerically_available); expect_true(signed$negative_raw_root_n_variance)
  expect_lt(signed$raw_root_n_variance, 0); expect_false(signed$Gaussian_interval_authorized)
  budget <- a; budget$control$maximum_analytic_matrix_products <- 1L
  rejected <- dv_fit_inference(budget)
  expect_identical(rejected$direct_variance$status, "direct_work_budget_exceeded")
  expect_false(rejected$direct_variance$available)
  expect_identical(rejected$fit$estimate, out$fit$estimate)
  strict <- a; strict$control$minimum_analytic_reciprocal_condition <- .99
  failed <- dv_fit_inference(strict)$direct_variance$raw_empirical
  expect_false(failed$numerically_available)
  expect_true(all(failed$pair_diagnostics$status == "numerical_evaluation_failed"))
  expect_match(paste(failed$pair_diagnostics$error, collapse = " "), "conditioning")
  expect_true(is.na(failed$raw_root_n_variance))
  expect_error(wdsmatch:::.wm_cp_direct_psd(diag(c(1, -.1)), "invalid", 1e-10), "PSD precision")
  rounded <- wdsmatch:::.wm_cp_direct_psd(diag(c(1, -1e-14)), "roundoff", 1e-10)
  expect_gt(max(rounded$negative_roundoff_correction), 0)
  expect_false(rounded$roundoff_included)
  wrong <- a; wrong$control$support <- "smooth_stack_threshold"
  expect_error(dv_fit_inference(wrong), "raw_only omits support")
  wrong <- a; wrong$control$direct_variance <- "unsupported"
  expect_error(dv_fit_inference(wrong), "direct_variance must")
  # High M cannot silently trim coefficient terms to meet a declared budget.
  expect_error(wdsmatch:::.wm_cp_G_average(matrix(c(1, 0, 1, 0, 0, 1, 0, 1), 4L, 2L),
    c(1, 1), matrix(0, 2L, 1L), matrix(0, 2L, 1L), matrix(0, 1L, 1L),
    1000L, dv_numeric_control()), "work budget")
})
