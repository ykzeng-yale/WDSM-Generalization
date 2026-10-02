# Fixed arrays and preflight checks only. The separately retained common-model
# validation fits the one new LP configuration; no legacy model is refitted here.
wm_lp_fixture <- function() {
  s <- as.matrix(expand.grid(x = seq(-1, 1, length.out = 9L),
                            u = seq(-1, 1, length.out = 9L)))
  n <- nrow(s)
  control <- wdsmatch:::.wm_lp_controls(list(method = "local_polynomial",
    degree = 2L, bandwidth_exponent = .05, guard_exponent = 6),
    c("0" = 2L), n, "PATT")[["0"]]
  list(scores = s, Z = rep(0, n), W = 1 + .15 * s[, 1] + .1 * s[, 2]^2,
    evaluation = rbind(c(0, 0), c(.12, -.18), c(.9, .1), c(-.95, -.7)),
    control = control,
    exponents = wdsmatch:::.wm_fit_polynomial_exponents(2L, 2L, 6L))
}

# Independent complete quadratic ordering and scalar intercept. Only evaluation
# s moves; the donor function, responses, W, h and training coordinates stay fixed.
wm_lp_scalar_reference <- function(a, Y, s) {
  delta <- sweep(a$scores, 2L, s, "-") / a$control$bandwidth
  radius <- rowSums(delta^2)
  K <- (4 / pi) * pmax(0, 1 - radius)^3
  r <- cbind(1, delta[, 1], delta[, 2], delta[, 1]^2,
             delta[, 1] * delta[, 2], delta[, 2]^2)
  w <- a$W * K / (length(Y) * a$control$bandwidth^2)
  G <- crossprod(r, r * w)
  moment <- crossprod(r, Y * w)
  list(value = as.numeric(solve(G, moment))[1L],
       minimum = min(eigen(G, symmetric = TRUE, only.values = TRUE)$values))
}

test_that("local-polynomial intercept and actual gradient reproduce a complete polynomial", {
  a <- wm_lp_fixture()
  polynomial <- function(s) 1 + 2 * s[, 1] - .5 * s[, 2] + .3 * s[, 1]^2 -
    .4 * s[, 1] * s[, 2] + .7 * s[, 2]^2
  Y <- polynomial(a$scores)
  out <- wdsmatch:::.wm_lp_arm(Y, a$Z, a$W, a$scores, 0L,
    a$control, a$exponents, a$evaluation)
  expect_identical(dim(a$exponents), c(6L, 2L))
  expect_identical(a$exponents[1L, ], c(score1 = 0L, score2 = 0L))
  expect_identical(nrow(unique(a$exponents)), 6L)
  expect_true(all(rowSums(a$exponents) <= 2L))
  expect_false(any(out$guard_failure))
  expect_false(any(out$derivative_boundary))
  expect_equal(out$mean, polynomial(a$evaluation), tolerance = 1e-9)
  exact <- cbind(x = 2 + .6 * a$evaluation[, 1] - .4 * a$evaluation[, 2],
                u = -.5 - .4 * a$evaluation[, 1] + 1.4 * a$evaluation[, 2])
  expect_equal(out$gradient, exact, tolerance = 1e-8)
  ref <- lapply(seq_len(nrow(a$evaluation)), function(i)
    wm_lp_scalar_reference(a, Y, a$evaluation[i, ]))
  expect_equal(out$minimum_gram_eigenvalue,
    vapply(ref, function(x) x$minimum, numeric(1L)), tolerance = 1e-10)
  expect_identical(out$donor_rows, seq_len(nrow(a$scores)))
  expect_identical(out$n, nrow(a$scores))
  expect_identical(out$ridge, 0)
  expect_false(out$assumptions_verified)
})

test_that("actual intercept derivatives include changing local weights and basis", {
  a <- wm_lp_fixture()
  Y <- sin(1.3 * a$scores[, 1]) + .3 * cos(2 * a$scores[, 2]) +
    .02 * sin(seq_len(nrow(a$scores)))
  out <- wdsmatch:::.wm_lp_arm(Y, a$Z, a$W, a$scores, 0L,
    a$control, a$exponents, a$evaluation)
  expect_true(all(out$gradient_available))
  fd <- matrix(0, nrow(a$evaluation), 2L)
  step <- 1e-5
  for (i in seq_len(nrow(a$evaluation))) for (k in 1:2) {
    plus <- minus <- a$evaluation[i, ]
    plus[k] <- plus[k] + step
    minus[k] <- minus[k] - step
    fd[i, k] <- (wm_lp_scalar_reference(a, Y, plus)$value -
                  wm_lp_scalar_reference(a, Y, minus)$value) / (2 * step)
  }
  expect_equal(unname(out$gradient), fd, tolerance = 1e-6)

  # A five-parameter reference derivative moves only its evaluation coordinate.
  B <- matrix(c(1, 0, 0, 1, .3, -.7, -.4, .2, .6, .5), 2L, 5L)
  contracted <- as.vector(crossprod(B, out$gradient[2L, ]))
  chain <- vapply(seq_len(ncol(B)), function(k) {
    s <- a$evaluation[2L, ]
    (wm_lp_scalar_reference(a, Y, s + step * B[, k])$value -
       wm_lp_scalar_reference(a, Y, s - step * B[, k])$value) / (2 * step)
  }, numeric(1L))
  expect_length(contracted, 5L)
  expect_equal(contracted, chain, tolerance = 1e-6)
})

test_that("local Gram uses full observed n rather than donor arm size", {
  a <- wm_lp_fixture()
  Y <- sin(a$scores[, 1]) + a$scores[, 2]^2
  original <- wdsmatch:::.wm_lp_arm(Y, a$Z, a$W, a$scores, 0L,
    a$control, a$exponents, a$evaluation)
  # Add fixed opposite-arm rows only. Keep donor data, evaluation points,
  # kernel amplitude, h and guard exactly fixed to isolate normalization.
  query_scores <- rbind(c(-.4, .2), c(.7, -.3), c(.1, .5))
  full_scores <- rbind(a$scores, query_scores)
  full_Y <- c(Y, 2, -1, .4)
  full_Z <- c(a$Z, 1, 1, 1)
  full_W <- c(a$W, .8, 1.5, 2)
  expanded <- wdsmatch:::.wm_lp_arm(full_Y, full_Z, full_W, full_scores, 0L,
    a$control, a$exponents, a$evaluation)
  expect_true(all(original$gradient_available))
  expect_true(all(expanded$gradient_available))
  expect_identical(expanded$donor_rows, original$donor_rows)
  expect_identical(expanded$n, length(full_Y))
  expect_equal(expanded$mean, original$mean, tolerance = 1e-9)
  expect_equal(expanded$gradient, original$gradient, tolerance = 1e-8)
  expect_equal(expanded$minimum_gram_eigenvalue,
    original$minimum_gram_eigenvalue * length(Y) / length(full_Y), tolerance = 1e-10)
})

test_that("finite-factor reference overflow is separate from a missing spatial gradient", {
  parameters <- c("score", "scale")
  first <- matrix(c(1e308, 3, 1, 1, 4, 2), 3L, 2L,
    dimnames = list(NULL, parameters))
  second <- matrix(c(1, 2, 1, 2, 3, 1), 3L, 2L,
    dimnames = list(NULL, parameters))
  gradient <- rbind(c(2, 1), c(.5, -.25), c(NA_real_, NA_real_))
  out <- wdsmatch:::.wm_lp_reference_derivative(list(first, second),
    gradient, c(TRUE, TRUE, FALSE), parameters)
  expect_true(all(is.finite(first)))
  expect_true(all(is.finite(gradient[1:2, ])))
  expect_identical(out$reference_derivative_failure, c(TRUE, FALSE, FALSE))
  expect_identical(out$reference_derivative_available, c(FALSE, TRUE, FALSE))
  expect_true(all(is.na(out$reference_mean_derivative[1L, ])))
  expect_equal(unname(out$reference_mean_derivative[2L, ]), c(1, 1.25), tolerance = 1e-12)
  expect_true(all(is.na(out$reference_mean_derivative[3L, ])))
  expect_identical(out$score_tangents[, 1L, ], first)
  expect_identical(dimnames(out$score_tangents)[[3L]], parameters)
})

test_that("Gram completion preserves raw-weight guard units and reports missing derivatives", {
  a <- wm_lp_fixture()
  Y <- sin(a$scores[, 1]) + a$scores[, 2]
  original <- wdsmatch:::.wm_lp_arm(Y, a$Z, a$W, a$scores, 0L,
    a$control, a$exponents, a$evaluation)
  scaled <- wdsmatch:::.wm_lp_arm(Y, a$Z, 7 * a$W, a$scores, 0L,
    a$control, a$exponents, a$evaluation)
  expect_equal(scaled$mean, original$mean, tolerance = 1e-9)
  expect_equal(scaled$gradient, original$gradient, tolerance = 1e-8)
  expect_equal(scaled$minimum_gram_eigenvalue,
    7 * original$minimum_gram_eigenvalue, tolerance = 1e-10)
  # Deliberately forced completion is a numerical unit check, not a new study.
  control <- a$control
  control$gram_guard <- 100
  failed <- wdsmatch:::.wm_lp_arm(Y, a$Z, a$W, a$scores, 0L,
    control, a$exponents, a$evaluation)
  expect_identical(failed$mean, numeric(nrow(a$evaluation)))
  expect_true(all(failed$guard_failure))
  expect_true(all(is.na(failed$gradient)))
  expect_false(any(failed$gradient_available))
  outside <- wdsmatch:::.wm_lp_arm(Y, a$Z, a$W, a$scores, 0L,
    a$control, a$exponents, matrix(c(10, 10), 1L))
  expect_identical(outside$mean, 0)
  expect_identical(outside$minimum_gram_eigenvalue, 0)
  expect_true(outside$guard_failure)
  expect_true(all(is.na(outside$gradient)))
})

test_that("LP controls and complete allocation fail before nuisance fitting", {
  base <- list(method = "local_polynomial", degree = 2L,
               bandwidth_exponent = .1)
  good <- wdsmatch:::.wm_lp_controls(base, c("0" = 3L), 500L, "PATT")[["0"]]
  expect_identical(good$gram_guard, 500^(-2))
  expect_identical(good$ridge, 0)
  expect_false(good$population_assumptions_verified)
  bad <- base; bad$degree <- 1L
  expect_error(wdsmatch:::.wm_lp_controls(bad, c("0" = 3L), 500L, "PATT"),
    "sufficient correction-transfer window", fixed = TRUE)
  bad <- base; bad$ridge <- .001
  expect_error(wdsmatch:::.wm_lp_controls(bad, c("0" = 3L), 500L, "PATT"),
    "no ridge", fixed = TRUE)
  bad <- base; bad$guard_exponent <- 1e300
  expect_error(wdsmatch:::.wm_lp_controls(bad, c("0" = 3L), 500L, "PATT"),
    "positive and representable", fixed = TRUE)
  design <- cbind(intercept = 1, x = seq_len(20L) / 20)
  ps <- rep(list(list(design = design, offset = rep(0, 20L), weighting = "probability")), 5L)
  pg <- list("0" = list(), "1" = list())
  layout <- wdsmatch:::.wm_model_layout(ps, pg, 0:1)
  controls <- wdsmatch:::.wm_lp_controls(list(method = "local_polynomial",
    degree = 10L, bandwidth_exponent = .04), layout$matching_dimensions, 20L, "PATE")
  expect_error(wdsmatch:::.wm_lp_allocation(20L, ps, pg, layout, controls, 5e7),
    "allocation exceeds max_stack_elements", fixed = TRUE)
})
