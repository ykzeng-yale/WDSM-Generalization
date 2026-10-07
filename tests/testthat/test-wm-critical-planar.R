# Finite arithmetic and conditional replication checks, not population/coverage
# certification. References do not call the candidate GH or whole-R evaluator.

cp_poly_add <- function(a, b) {
  for (key in names(b)) a[key] <- if (key %in% names(a)) a[key] + b[key] else b[key]
  a
}
cp_poly_multiply <- function(a, b) {
  out <- numeric()
  for (ka in names(a)) for (kb in names(b)) {
    power <- as.integer(strsplit(ka, ",", fixed = TRUE)[[1]]) +
      as.integer(strsplit(kb, ",", fixed = TRUE)[[1]])
    key <- paste(power, collapse = ",")
    value <- a[ka] * b[kb]
    out[key] <- if (key %in% names(out)) out[key] + value else value
  }
  out
}
cp_quadratic <- function(A, c, lambda) {
  H <- lambda * crossprod(A); b <- lambda * as.vector(crossprod(A, c))
  stats::setNames(c(lambda * sum(c^2), 2 * b, H[1, 1], 2 * H[1, 2], H[2, 2]),
    c("0,0", "1,0", "0,1", "2,0", "1,1", "0,2"))
}
cp_poisson_polynomial <- function(t, M) {
  sum <- power <- c("0,0" = 1)
  if (M > 1L) for (k in seq_len(M - 1L)) {
    power <- cp_poly_multiply(power, t)
    sum <- cp_poly_add(sum, power / factorial(k))
  }
  sum
}
cp_centered_moment <- function(a, b, V) {
  if (a < 0 || b < 0) return(0)
  if (a + b == 0) return(1)
  if ((a + b) %% 2) return(0)
  if (a) return((a - 1) * V[1, 1] * cp_centered_moment(a - 2, b, V) +
    b * V[1, 2] * cp_centered_moment(a - 1, b - 1, V))
  (b - 1) * V[2, 2] * cp_centered_moment(0, b - 2, V)
}
cp_normal_moment <- function(a, b, mu, V) {
  out <- 0
  for (k in 0:a) for (l in 0:b) out <- out + choose(a, k) * choose(b, l) *
    mu[1]^(a - k) * mu[2]^(b - l) * cp_centered_moment(k, l, V)
  out
}
cp_G_moment_reference <- function(E, density, c0, c1, M) {
  A0 <- E[1:2, , drop = FALSE]; A1 <- E[3:4, , drop = FALSE]
  lambda <- pi * density
  H <- lambda[1] * crossprod(A0) + lambda[2] * crossprod(A1)
  b <- lambda[1] * as.vector(crossprod(A0, c0)) + lambda[2] * as.vector(crossprod(A1, c1))
  mu <- -as.vector(solve(H, b)); V <- solve(2 * H)
  kappa <- lambda[1] * sum((A0 %*% mu + c0)^2) +
    lambda[2] * sum((A1 %*% mu + c1)^2)
  P <- cp_poly_multiply(cp_poisson_polynomial(cp_quadratic(A0, c0, lambda[1]), M),
    cp_poisson_polynomial(cp_quadratic(A1, c1, lambda[2]), M))
  expectation <- 0
  for (key in names(P)) {
    powers <- as.integer(strsplit(key, ",", fixed = TRUE)[[1]])
    expectation <- expectation + P[key] * cp_normal_moment(powers[1], powers[2], mu, V)
  }
  unname(pi * exp(-kappa) / sqrt(det(H)) * expectation)
}

test_that("finite Hermite geometry agrees with independent Gaussian moments", {
  E <- qr.Q(qr(matrix(c(1, .2, .6, .1, .3, 1.1, -.2, .9), 4, 2)))
  rotation <- matrix(c(cos(.4), sin(.4), -sin(.4), cos(.4)), 2, 2)
  for (M in 1:3) {
    rule <- WeightedMatching:::.wm_cp_hermite(M)
    expect_equal(sum(rule$weights), 1, tolerance = 1e-13)
    for (k in 0:(2 * rule$order - 1)) {
      exact <- if (k %% 2) 0 else gamma((k + 1) / 2) / sqrt(pi)
      expect_equal(sum(rule$weights * rule$nodes^k), exact, tolerance = 1e-11)
    }
    for (u in c(0, .2, -1)) {
      c0 <- u * c(.3, -.2); c1 <- u * c(-.1, .4)
      ref <- cp_G_moment_reference(E, c(.3, .8), c0, c1, M)
      got <- WeightedMatching:::.wm_cp_G(E, c(.3, .8), c0, c1, M, rule)
      expect_equal(got$value, ref, tolerance = 1e-10)
      expect_equal(WeightedMatching:::.wm_cp_G(E %*% rotation, c(.3, .8), c0, c1, M, rule)$value,
        ref, tolerance = 1e-10)
      expect_identical(got$quadrature_error_exact_arithmetic, 0)
      expect_false(got$roundoff_included)
    }
  }
})

cp_fixture <- function(M = 1L, unit = FALSE) {
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
cp_fit_inference <- function(a, fit = NULL) {
  if (is.null(fit)) fit <- do.call(wm_match, a$data)
  zero <- list(mode = "zero", basis = "fixed_map")
  wm_fitted_inference(fit, a$influence, a$D0, a$D1,
    transport0 = zero, transport1 = zero, covariance_scope = "critical_planar",
    critical_control = a$control)
}
cp_A_exhaustive <- function(w, donor, probability, M) {
  if (M == 1) return(1)
  grid <- expand.grid(rep(list(seq_along(donor)), M - 1), KEEP.OUT.ATTRS = FALSE)
  answer <- 0
  for (r in seq_len(nrow(grid))) {
    index <- as.integer(unlist(grid[r, , drop = FALSE], use.names = FALSE))
    answer <- answer + prod(probability[index]) * w / (w + sum(donor[index]))
  }
  answer
}
cp_R_reference <- function(a, u) {
  dat <- a$data; W <- dat$weights / max(dat$weights); S0 <- dat$scores0; S1 <- dat$scores1
  T <- cbind(S0, S1); aux <- a$control$auxiliary_rows; ev <- setdiff(seq_along(W), aux)
  h <- a$control$bandwidth; residual <- dat$Y - ifelse(dat$Z == 1, dat$mean1, dat$mean0)
  R <- 0
  for (i in ev[dat$Z[ev] == 1]) {
    v <- sweep(T[aux, , drop = FALSE], 2, T[i, ], "-") / h
    weight <- pmax(0, 1 - sqrt(rowSums(v^2)))
    mass <- sum(weight); center <- numeric(4); Q <- matrix(0, 4, 4)
    for (k in seq_along(aux)) center <- center + weight[k] * v[k, ] / mass
    for (k in seq_along(aux)) Q <- Q + weight[k] *
      tcrossprod(v[k, ] - center) / mass
    E <- eigen(Q, symmetric = TRUE)$vectors[, 1:2, drop = FALSE]
    eigenvalues <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values
    projection <- c(min(svd(E[1:2, ], nu = 0, nv = 0)$d),
      min(svd(E[3:4, ], nu = 0, nv = 0)$d))
    if (eigenvalues[2] < .075 || any(projection < a$control$projection_floor)) next
    laws <- list(); density <- numeric(2)
    for (z in 0:1) {
      S <- if (z == 0) S0 else S1; donor <- aux[dat$Z[aux] == z]
      k <- pmax(0, 1 - sqrt(rowSums(sweep(S[donor, , drop = FALSE], 2, S[i, ], "-")^2)) / h) /
        (pi / 3)
      density[z + 1] <- sum(k) / (length(aux) * h^2)
      laws[[z + 1]] <- list(donor = W[donor], probability = k / sum(k))
    }
    for (j in ev[dat$Z[ev] == 0]) {
      k <- max(0, 1 - sqrt(sum((T[j, ] - T[i, ])^2)) / h)
      A1 <- cp_A_exhaustive(W[i], laws[[2]]$donor, laws[[2]]$probability, dat$M)
      A0 <- cp_A_exhaustive(W[j], laws[[1]]$donor, laws[[1]]$probability, dat$M)
      c0 <- (matrix(a$control$tangents0[j, , ], 2, 2) -
        matrix(a$control$tangents0[i, , ], 2, 2)) %*% u
      c1 <- (matrix(a$control$tangents1[j, , ], 2, 2) -
        matrix(a$control$tangents1[i, , ], 2, 2)) %*% u
      R <- R + k * W[i] * residual[i] * A1 * W[j] * residual[j] * A0 *
        cp_G_moment_reference(E, density, c0, c1, dat$M)
    }
  }
  R / (length(ev) * (length(ev) - 1) * h^2 * (pi / 3))
}

# Direct donor imputation/load/derivative reference; stored contribution rows,
# load matrices, sensitivities and covariance helpers are not used.
cp_row_reference <- function(a, neighbors) {
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

test_that("whole reciprocal function matches exhaustive finite donor/pair sums", {
  for (M in 1:3) for (unit in c(FALSE, TRUE)) {
    a <- cp_fixture(M, unit); fit <- do.call(wm_match, a$data); previous <- fit
    out <- cp_fit_inference(a, fit)
    expect_true(out$available)
    expect_identical(fit, previous); expect_identical(out$fit, fit)
    expect_equal(out$estimate, fit$estimate)
    row_ref <- cp_row_reference(a, fit$graph$neighbors)
    expect_equal(out$estimate, row_ref$estimate, tolerance = 1e-12)
    expect_equal(out$base_rows, row_ref$base, tolerance = 1e-12)
    expect_equal(out$diagonal_root_n_variance, row_ref$diagonal, tolerance = 1e-12)
    expect_equal(unname(out$raw_C), row_ref$C, tolerance = 1e-12)
    expect_equal(out$raw_Sigma, row_ref$Sigma, tolerance = 1e-12)
    expect_equal(unname(out$total_sensitivity), row_ref$b, tolerance = 1e-12)
    expect_false(out$assumptions_verified); expect_false(out$application_verified)
    expect_true(is.na(out$root_n_variance)); expect_true(all(is.na(out$conf.int)))
    for (u in list(c(0, 0), c(.3, -.2), c(-1, .5))) {
      got <- WeightedMatching:::.wm_cp_evaluate(out$critical$kernel, u)
      ref <- cp_R_reference(a, u)
      expect_lte(abs(got$R - ref), got$quadrature_error_bound + 1e-10)
      expect_false(got$roundoff_included)
    }
    if (unit) expect_equal(out$critical$kernel$global_quadrature_error_bound, 0)
    scaled <- a; scaled$data$weights <- 29 * scaled$data$weights
    other <- cp_fit_inference(scaled)
    expect_equal(other$critical$origin, out$critical$origin, tolerance = 1e-12)
    expect_equal(other$C, out$C, tolerance = 1e-12)
    expect_equal(other$Sigma, out$Sigma, tolerance = 1e-12)
    expect_equal(other$total_sensitivity, out$total_sensitivity, tolerance = 1e-12)
  }
})

test_that("support and cone retain raw inputs including singular and zero rank", {
  a <- cp_fixture(); a$influence[, 2] <- a$influence[, 1]
  out <- cp_fit_inference(a); s <- out$critical$support
  expect_equal(s$threshold, out$n^(-1 / 3))
  expect_equal(s$rank, 1)
  expect_equal(s$Sigma %*% s$Sigma_plus, s$projector, tolerance = 1e-12)
  expect_equal(drop(s$projector %*% s$C), s$C, tolerance = 1e-12)
  expect_identical(dim(out$nuisance_influence), c(12L, 2L))
  # Internal function arithmetic at a declared covariance boundary, not an
  # altered public fit or a claim that this artificial input is population-valid.
  R0 <- cp_R_reference(a, c(0, 0))
  for (schur in c(-1, 0, 1)) {
    arithmetic <- out
    arithmetic$diagonal_root_n_variance <- s$q + schur + 2 * R0 / out$fit$gamma^2
    z <- WeightedMatching:::.wm_cp_schur(arithmetic, c(0, 0))
    expect_equal(z$raw_D, schur, tolerance = 1e-12)
    expect_equal(z$D, max(schur, 0), tolerance = 1e-12)
    expect_equal(z$cone_adjustment, max(schur, 0) - schur, tolerance = 1e-12)
  }
  for (u in list(c(0, 0), c(.2, -.3))) {
    z <- WeightedMatching:::.wm_cp_schur(out, u)
    expect_equal(z$raw_D, z$raw_V0 - s$q)
    expect_equal(z$D, max(z$raw_D, 0))
    expect_equal(z$cone_adjustment, z$D - z$raw_D)
    t <- .4; v <- c(.2, -.1)
    block <- rbind(c(z$V0, s$C), cbind(s$C, s$Sigma))
    expect_equal(as.numeric(crossprod(c(t, v), block %*% c(t, v))),
      t^2 * z$D + as.numeric(crossprod(v + t * s$regression,
        s$Sigma %*% (v + t * s$regression))), tolerance = 1e-11)
  }
  zero <- a; zero$influence[,] <- 0; zero <- cp_fit_inference(zero)
  expect_equal(zero$critical$support$rank, 0)
  expect_equal(zero$critical$support$q, 0)
  expect_equal(unname(zero$Sigma), matrix(0, 2, 2))
})

test_that("full-rank support and mixture agree with direct two-axis covariance", {
  a <- cp_fixture(); a$influence[, 2] <- 10 * a$influence[, 2]
  fit <- do.call(wm_match, a$data); out <- cp_fit_inference(a, fit)
  ref <- cp_row_reference(a, fit$graph$neighbors)
  Sigma <- unname(ref$Sigma); C <- ref$C
  expect_gt(min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values),
    fit$n^(-1 / 3))
  support <- out$critical$support
  expect_equal(support$rank, 2)
  expect_equal(unname(support$Sigma), Sigma, tolerance = 1e-12)
  expect_equal(unname(support$Sigma_plus), solve(Sigma), tolerance = 1e-12)
  expect_equal(unname(support$projector), diag(2), tolerance = 1e-12)
  expect_equal(unname(support$C), C, tolerance = 1e-12)
  regression <- as.vector(solve(Sigma, C))
  q <- sum(C * regression)
  expect_equal(unname(support$regression), regression, tolerance = 1e-12)
  expect_equal(support$q, q, tolerance = 1e-12)
  # Closed two-by-two principal square root, independent of spectral helper.
  determinant_root <- sqrt(det(Sigma))
  root <- (Sigma + determinant_root * diag(2)) /
    sqrt(sum(diag(Sigma)) + 2 * determinant_root)
  expect_equal(unname(support$square_root), root, tolerance = 1e-12)
  set.seed(824); previous <- .Random.seed
  set.seed(43); normals <- matrix(rnorm(3 * 4), 3, 4); .Random.seed <- previous
  H <- root %*% normals[1:2, , drop = FALSE]
  D <- expected <- numeric(4)
  for (b in 1:4) {
    R <- cp_R_reference(a, H[, b])
    D[b] <- max(ref$diagonal - 2 * R / fit$gamma^2 - q, 0)
    expected[b] <- sum((ref$b + regression) * H[, b]) +
      sqrt(D[b]) * normals[3, b]
  }
  got <- wm_bootstrap(out, B = 4, seed = 43, chunk_size = 2)
  expect_identical(.Random.seed, previous)
  expect_true(got$available)
  expect_equal(unname(got$nuisance_draws), H, tolerance = 1e-12)
  expect_equal(got$draw_diagnostics$D, D, tolerance = 1e-10)
  expect_equal(got$root_n_draws, expected, tolerance = 1e-10)
  expect_equal(got$rao_blackwell_root_n_variance_estimate,
    mean(D) + as.numeric(crossprod(ref$b + regression,
      Sigma %*% (ref$b + regression))), tolerance = 1e-10)
  expect_equal(got$rao_blackwell_root_n_variance_mcse, sd(D) / 2,
    tolerance = 1e-10)
})

test_that("common bootstrap produces complete nuisance Schur-mixture roots", {
  a <- cp_fixture(); out <- cp_fit_inference(a)
  set.seed(814); previous <- .Random.seed
  set.seed(26); normals <- matrix(rnorm(3 * 9), 3, 9); .Random.seed <- previous
  got <- wm_bootstrap(out, B = 9, seed = 26, chunk_size = 2)
  expect_identical(.Random.seed, previous)
  expect_identical(got$interval, "basic")
  s <- out$critical$support; H <- s$square_root %*% normals[1:2, , drop = FALSE]
  expected <- numeric(9); D <- numeric(9)
  for (b in 1:9) {
    R <- cp_R_reference(a, H[, b])
    raw_D <- out$diagonal_root_n_variance - 2 * R / out$fit$gamma^2 - s$q
    D[b] <- max(raw_D, 0)
    expected[b] <- sum((out$total_sensitivity + s$regression) * H[, b]) +
      sqrt(D[b]) * normals[3, b]
  }
  expect_equal(got$nuisance_draws, H, tolerance = 1e-12)
  expect_equal(got$root_n_draws, expected, tolerance = 1e-10)
  expect_equal(got$monte_carlo_root_n_variance, var(got$root_n_draws))
  expect_equal(got$rao_blackwell_root_n_variance_estimate,
    mean(D) + sum(crossprod(s$square_root, out$total_sensitivity + s$regression)^2),
    tolerance = 1e-10)
  expect_equal(got$rao_blackwell_root_n_variance_mcse, sd(D) / 3, tolerance = 1e-10)
  expect_true(is.na(got$conditional_root_n_variance))
  expect_equal(unname(got$conf.int), out$estimate -
    unname(quantile(got$root_n_draws, c(.975, .025))) / sqrt(out$n))
  expect_equal(wm_bootstrap(out, B = 9, seed = 26, chunk_size = 20)$root_n_draws,
    got$root_n_draws, tolerance = 0)
  single <- wm_bootstrap(out, B = 1, seed = 26, interval = "none")
  expect_true(is.na(single$monte_carlo_variance))
  expect_true(is.na(single$rao_blackwell_variance_mcse))
  expect_null(single$conf.int)
  expect_error(wm_bootstrap(out, interval = "normal"), "arg")
  expect_error(wm_bootstrap(out, counts = matrix(1, 12, 2)), "rejects counts")
  expect_error(wm_bootstrap(out, refit = function(...) stop("never")), "fixed_reuse")
  expect_identical(got$source_inference, out)
  expect_false(got$original_refit_bootstrap)
})

test_that("critical controls fail honestly and construction uses no RNG", {
  a <- cp_fixture(); fit <- do.call(wm_match, a$data)
  set.seed(528); previous <- .Random.seed; out <- cp_fit_inference(a, fit)
  expect_identical(.Random.seed, previous)
  small <- a; small$control$maximum_pairs <- 1
  fail <- cp_fit_inference(small, fit)
  expect_false(fail$available); expect_match(fail$unavailable_reason, "budget")
  expect_identical(fail$fit, fit)
  expect_error(wm_bootstrap(fail), "successful")
  changed <- out; changed$critical_control$tangents0[1, 1, 1] <- 1
  expect_error(wm_bootstrap(changed), "binding")
  changed <- out; changed$critical$kernel$pairs[[1]]$coefficient <- 100
  expect_error(wm_bootstrap(changed), "cache binding")
  changed <- out; changed$critical$support$regression[1] <-
    changed$critical$support$regression[1] + 1
  expect_error(wm_bootstrap(changed), "Critical support regression")
  changed <- a; dimnames(changed$control$tangents0)[[3]] <- c("b", "a")
  expect_error(cp_fit_inference(changed), "parameter order")
  changed <- a; changed$control$support <- "generic_L2"
  expect_error(cp_fit_inference(changed), "covariance-rate")
  expect_error(wm_fitted_inference(fit, a$influence, covariance_scope = "full_x",
    critical_control = a$control), "only supported")
  zero <- list(mode = "zero", basis = "fixed_map")
  expect_error(wm_fitted_inference(fit, a$influence,
    weight_derivative = matrix(1, 12, 2, dimnames = list(NULL, c("a", "b"))),
    covariance_scope = "critical_planar", critical_control = a$control,
    transport0 = zero, transport1 = zero), "known")
  budget <- out; budget$critical_control$maximum_draw_entries <- 1L
  budget$critical$kernel$control$maximum_draw_entries <- 1L
  expect_error(wm_bootstrap(budget, B = 2, seed = 1), "before RNG")
  expect_identical(.Random.seed, previous)
})

test_that("auxiliary guards and all-B evaluation failures remain explicit", {
  a <- cp_fixture(); original <- do.call(wm_match, a$data)
  high <- a; high$control$projection_floor <- .99
  rejected <- cp_fit_inference(high, original)
  expect_equal(rejected$critical$kernel$diagnostics$rejected_planes, 2)
  expect_equal(rejected$critical$kernel$diagnostics$retained_pairs, 0)
  expect_identical(rejected$fit, original)
  high <- a; high$control$marginal_mass_floor <- 100
  guarded <- cp_fit_inference(high, original)
  expect_equal(guarded$critical$kernel$diagnostics$marginal0_guard_count, 2)
  expect_equal(guarded$critical$kernel$diagnostics$marginal1_guard_count, 2)
  clipped <- a; clipped$control$residual_cap <- .2; clipped$control$tangent_cap <- .01
  clipped <- cp_fit_inference(clipped, original)
  expect_gt(clipped$critical$kernel$diagnostics$residual_cap_entries, 0)
  expect_gt(clipped$critical$kernel$diagnostics$tangent0_cap_entries, 0)
  expect_identical(clipped$fit, original)
  out <- cp_fit_inference(a, original)
  # Deliberate runtime error injection checks retention/availability only;
  # expected numerical roots elsewhere use independent unmodified formulas.
  ns <- asNamespace("WeightedMatching"); name <- ".wm_cp_schur"
  previous <- get(name, envir = ns, inherits = FALSE)
  locked <- bindingIsLocked(name, ns)
  on.exit({
    if (bindingIsLocked(name, ns)) unlockBinding(name, ns)
    assign(name, previous, ns)
    if (locked) lockBinding(name, ns)
  }, add = TRUE)
  if (locked) unlockBinding(name, ns)
  assign(name, function(...) stop("fixture evaluation failure"), ns)
  if (locked) lockBinding(name, ns)
  failed <- wm_bootstrap(out, B = 3, seed = 3)
  expect_equal(length(failed$root_n_draws), 3)
  expect_equal(ncol(failed$nuisance_draws), 3)
  expect_true(all(is.na(failed$root_n_draws)))
  expect_true(all(failed$draw_diagnostics$status == "numerical_evaluation_failed"))
  expect_false(failed$available)
  expect_true(all(is.na(failed$conf.int)))
  expect_true(is.na(failed$monte_carlo_root_n_variance))
  expect_true(is.na(failed$rao_blackwell_root_n_variance_estimate))
})

test_that("critical seed absence and unsupported dimensions remain guarded", {
  a <- cp_fixture(); out <- cp_fit_inference(a)
  saved_kind <- RNGkind()
  had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  saved_seed <- if (had_seed) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    do.call(RNGkind, as.list(saved_kind))
    if (had_seed) assign(".Random.seed", saved_seed, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind(normal.kind = "Inversion")
  if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
    rm(".Random.seed", envir = .GlobalEnv)
  invisible(wm_bootstrap(out, B = 1, seed = 28, interval = "none"))
  expect_false(exists(".Random.seed", .GlobalEnv, inherits = FALSE))
  RNGkind(normal.kind = "Box-Muller")
  expect_error(wm_bootstrap(out, B = 1, seed = 28), "Box-Muller")
  patt <- a; patt$data$estimand <- "PATT"
  patt$data[c("scores1", "mean1")] <- NULL
  patt_fit <- do.call(wm_match, patt$data)
  expect_error(cp_fit_inference(a, patt_fit), "PATT")
  scalar <- a; scalar$data$scores0 <- scalar$data$scores0[, 1, drop = FALSE]
  scalar$data$scores1 <- scalar$data$scores1[, 1, drop = FALSE]
  expect_error(cp_fit_inference(a, do.call(wm_match, scalar$data)), "exactly two")
})
