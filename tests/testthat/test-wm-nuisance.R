test_that("complete weighted polynomials recover known response surfaces", {
  grid <- expand.grid(x = seq(0, 1, length.out = 4),
                      v = seq(0, 1, length.out = 4), z = 0:1)
  s0 <- as.matrix(grid[c("x", "v")])
  s1 <- matrix(grid$x, ncol = 1)
  m0 <- 1 + 2 * grid$x - grid$v + grid$x^2 + 3 * grid$x * grid$v
  m1 <- 5 - grid$x
  y <- ifelse(grid$z == 0, m0, m1)
  w <- 1 + grid$x + 2 * grid$v
  fit <- wm_fit(y, grid$z, w, s0, scores1 = s1, M = 2,
                degree = c(2, 1))
  expect_s3_class(fit, "wm_fit")
  expect_s3_class(fit, "wm_match")
  expect_equal(fit$nuisance$predictions$mean0, m0, tolerance = 1e-10)
  expect_equal(fit$nuisance$predictions$mean1, m1, tolerance = 1e-10)
  expect_equal(fit$estimate, sum(w * (m1 - m0)) / sum(w),
               tolerance = 1e-10)
  control <- fit$nuisance$models[[1]]
  expect_equal(control$basis$n_columns, 6)
  expect_equal(control$mean$rank, 6)
  expect_true(all(rowSums(control$basis$exponents) <= 2))
  expect_equal(control$training_rows, which(grid$z == 0))
  scaled <- wm_fit(y, grid$z, 7 * w, s0, scores1 = s1, M = 2,
                   degree = c(2, 1))
  expect_equal(scaled$estimate, fit$estimate, tolerance = 1e-10)
  expect_equal(scaled$nuisance$predictions$mean0,
               fit$nuisance$predictions$mean0, tolerance = 1e-10)
})

test_that("outcome fits use W while conditional weight fits use unit weights", {
  z <- rep(0:1, each = 4)
  x <- rep(seq(0, 1, length.out = 4), 2)
  y <- c(1, 2, 4, 8, 2, 4, 7, 9)
  w <- c(1, 2, 4, 8, 2, 3, 5, 9)
  fit <- wm_fit(y, z, w, x, M = 1, degree = 0,
                method = "stabilized", rho_bounds = c(0.1, 20))
  for (arm in 0:1) {
    rows <- which(z == arm)
    expect_equal(fit$nuisance$predictions[[paste0("mean", arm)]],
                 rep(sum(w[rows] * y[rows]) / sum(w[rows]), length(y)))
    expect_equal(fit$nuisance$predictions[[paste0("rho", arm)]],
                 rep(mean(w[rows]), length(y)))
  }
  expect_true(fit$nuisance$models[[1]]$mean$weighted)
  expect_false(fit$nuisance$models[[1]]$rho$weighted)
  expect_equal(fit$nuisance$rho_clipped_count, 0)
})

test_that("PATT never fits an unnecessary treated regression", {
  x <- c(0, 0.25, 0.5, 1, 0.6)
  z <- c(0, 0, 0, 0, 1)
  y <- c(2 + 3 * x[1:4], 10)
  fit <- wm_fit(y, z, c(1, 2, 3, 2, 4), x, M = 2, estimand = "PATT")
  expect_equal(fit$estimate, 10 - (2 + 3 * 0.6), tolerance = 1e-12)
  expect_length(fit$nuisance$models, 1)
  expect_equal(fit$nuisance$models[[1]]$arm, 0)
  expect_null(fit$nuisance$predictions$mean1)
  expect_null(fit$nuisance$predictions$rho1)
  expect_error(wm_fit(y, z, rep(1, 5), x, scores1 = x,
                      estimand = "PATT"), "omit scores1")
})

test_that("tensor B-splines reproduce endpoint values and bilinear surfaces", {
  grid <- expand.grid(x = seq(0, 1, length.out = 4),
                      v = seq(0, 1, length.out = 4), z = 0:1)
  s <- as.matrix(grid[c("x", "v")])
  m0 <- 1 + grid$x + 2 * grid$v + 3 * grid$x * grid$v
  m1 <- 5 - 2 * grid$x + grid$v
  y <- ifelse(grid$z == 0, m0, m1)
  support <- rbind(c(0, 1), c(0, 1))
  fit <- wm_fit(y, grid$z, 1 + grid$x, s, M = 1,
                regression = "spline", spline_order = 2,
                mesh_exponent = 0.15, moment_order = Inf, support0 = support)
  expect_equal(fit$nuisance$predictions$mean0, m0, tolerance = 1e-10)
  expect_equal(fit$nuisance$predictions$mean1, m1, tolerance = 1e-10)
  spec <- fit$nuisance$models[[1]]$basis
  basis <- wdsmatch:::.wm_fit_basis(s, spec)
  expect_equal(rowSums(basis), rep(1, nrow(s)), tolerance = 1e-12)
  expect_equal(spec$knots[[1]], c(0, 0, 0.5, 1, 1))
  expect_equal(unname(spec$support), support)
  corners <- which(grid$x %in% c(0, 1) & grid$v %in% c(0, 1))
  expect_equal(fit$nuisance$predictions$mean0[corners], m0[corners],
               tolerance = 1e-10)
  expect_equal(spec$tensor_index[1, ], c(1L, 1L))
  expect_equal(spec$tensor_index[spec$n_columns, ], c(3L, 3L))
})

wm_nuisance_fold_data <- function() {
  dat <- expand.grid(x = seq(0, 1, length.out = 12), z = 0:1, fold = 1:2)
  dat$w <- 1 + dat$x + 0.2 * dat$z
  dat$y <- 2 + 3 * dat$x + dat$z + 0.1 * (dat$fold == 2) * sin(6 * dat$x)
  dat
}

wm_nuisance_fold_fit <- function(dat) {
  wm_fit(dat$y, dat$z, dat$w, dat$x, M = 2, method = "stabilized",
         regression = "spline", spline_order = 2, mesh_exponent = 0.2,
         moment_order = Inf, support0 = c(0, 1), rho_bounds = c(0.1, 10),
         fold_id = dat$fold)
}

test_that("held-out spline fits exclude the entire evaluation fold", {
  dat <- wm_nuisance_fold_data()
  fit <- wm_nuisance_fold_fit(dat)
  changed <- dat
  held_out <- dat$fold == 1
  changed$y[held_out] <- changed$y[held_out] + 20 * dat$x[held_out]^2
  changed$w[held_out] <- 2 * changed$w[held_out]
  refit <- wm_nuisance_fold_fit(changed)
  for (model in fit$nuisance$models) {
    expect_length(intersect(model$training_rows, model$evaluation_rows), 0)
    expect_true(all(dat$fold[model$training_rows] != as.numeric(model$fold)))
  }
  for (prediction in c("mean0", "mean1", "rho0", "rho1")) {
    expect_equal(fit$nuisance$predictions[[prediction]][held_out],
                 refit$nuisance$predictions[[prediction]][held_out],
                 tolerance = 1e-12)
  }
  expect_false(isTRUE(all.equal(fit$nuisance$predictions$rho0[!held_out],
                               refit$nuisance$predictions$rho0[!held_out])))
  expect_true(all(dat$fold[fit$graph$edges$query] == dat$fold[fit$graph$edges$donor]))
})

test_that("strata restrict nuisance training as well as matching", {
  dat <- expand.grid(x = seq(0, 1, length.out = 5),
                     z = 0:1, stratum = c("A", "B"))
  m0 <- ifelse(dat$stratum == "A", 1 + dat$x, 10 - dat$x)
  m1 <- ifelse(dat$stratum == "A", 7 + 2 * dat$x, -5 + 3 * dat$x)
  dat$y <- ifelse(dat$z == 0, m0, m1)
  fit <- wm_fit(dat$y, dat$z, 1 + dat$x, dat$x, M = 2,
                strata = dat$stratum)
  expect_equal(fit$nuisance$predictions$mean0, m0, tolerance = 1e-12)
  expect_equal(fit$nuisance$predictions$mean1, m1, tolerance = 1e-12)
  for (model in fit$nuisance$models) {
    expect_true(all(as.character(dat$stratum[model$training_rows]) == model$stratum))
  }
  expect_true(all(dat$stratum[fit$graph$edges$query] == dat$stratum[fit$graph$edges$donor]))
})

test_that("fold and stratum identities preserve close numbers and dotted labels", {
  dat <- expand.grid(x = seq(0, 1, length.out = 4), z = 0:1,
                     fold_index = 1:2, stratum_index = 1:2)
  numeric_fold <- ifelse(dat$fold_index == 1, 1, 1 + .Machine$double.eps)
  numeric_stratum <- ifelse(dat$stratum_index == 1, 2,
                            2 + 2 * .Machine$double.eps)
  dat$y <- 1 + dat$x + 4 * dat$z + 10 * dat$stratum_index
  fit <- wm_fit(dat$y, dat$z, rep(1, nrow(dat)), dat$x, M = 1,
                fold_id = numeric_fold, strata = numeric_stratum)
  expect_length(fit$nuisance$models, 8)
  expect_equal(fit$nuisance$fold_id, numeric_fold)
  expect_equal(fit$nuisance$strata, numeric_stratum)
  for (model in fit$nuisance$models) {
    expect_true(all(numeric_fold[model$training_rows] != model$fold))
    expect_true(all(numeric_stratum[model$training_rows] == model$stratum))
  }
  dotted_fold <- c("a.b", "a")[dat$fold_index]
  dotted_stratum <- c("c", "b.c")[dat$stratum_index]
  dotted <- wm_fit(dat$y, dat$z, rep(1, nrow(dat)), dat$x, M = 1,
                   fold_id = dotted_fold, strata = dotted_stratum)
  expect_length(dotted$nuisance$models, 8)
  expect_true(all(dotted_fold[dotted$graph$edges$query] ==
                    dotted_fold[dotted$graph$edges$donor]))
  expect_true(all(dotted_stratum[dotted$graph$edges$query] ==
                    dotted_stratum[dotted$graph$edges$donor]))
})

test_that("normalizer clipping is recorded without changing the fit", {
  # Deliberately tight bounds exercise clipping; this is an implementation
  # diagnostic, not a configuration satisfying the true-range assumption.
  z <- rep(0:1, each = 4)
  fit <- wm_fit(seq_len(8), z, rep(c(1, 2, 3, 4), 2),
                rep(seq(0, 1, length.out = 4), 2), M = 1,
                method = "stabilized", degree = 0, rho_bounds = c(0.5, 2))
  expect_equal(fit$nuisance$predictions$rho0, rep(2, 8))
  expect_equal(fit$nuisance$predictions$rho1, rep(2, 8))
  expect_equal(fit$nuisance$rho_clipped_count, 16)
  expect_equal(unname(fit$nuisance$models[[1]]$rho$coefficients), 2.5)
  expect_equal(fit$nuisance$models[[1]]$rho$clipped_rows, seq_len(8))
})

test_that("unsupported basis and sample configurations fail explicitly", {
  x <- rep(seq(0, 1, length.out = 4), 2)
  z <- rep(0:1, each = 4)
  y <- 1 + x + z
  expect_error(wm_fit(y, z, rep(1, 8), cbind(x, x), M = 1),
               "rank-deficient")
  expect_error(wm_fit(y, z, rep(1, 8), rep(0, 8), M = 1),
               "zero training basis")
  expect_error(wm_fit(y, z, rep(1, 8), x, M = 5),
               "insufficient evaluation donors")
  expect_error(wm_fit(y, z, rep(1, 8), x, M = 1, degree = 2, max_basis = 2),
               "max_basis")
  expect_error(wm_fit(y, z, rep(1, 8), x, M = 1, max_basis_elements = 2),
               "max_basis_elements")
  # Both evaluation folds have donors; one training arm has fewer rows than
  # the full quadratic basis, which must not trigger silent term deletion.
  fold <- c(1, 1, 1, 2, 1, 1, 2, 2)
  expect_error(wm_fit(y, z, rep(1, 8), x, M = 1, degree = 2,
                      fold_id = fold), "insufficient training donors")
  expect_error(wm_fit(y, z, rep(1, 8), x, M = 1,
                      fold_id = rep(1, 8)), "at least two")
})

test_that("spline contracts require declared support, rates and folds", {
  dat <- wm_nuisance_fold_data()
  args <- list(Y = dat$y, Z = dat$z, weights = dat$w, scores0 = dat$x,
               M = 1, regression = "spline", spline_order = 2,
               mesh_exponent = 0.2, moment_order = Inf, support0 = c(0, 1))
  bad <- args
  bad$support0 <- NULL
  expect_error(do.call(wm_fit, bad), "support0")
  bad <- args
  bad$support0 <- c(0.1, 1)
  expect_error(do.call(wm_fit, bad), "out-of-support")
  bad <- args
  bad$moment_order <- NULL
  expect_error(do.call(wm_fit, bad), "explicit")
  bad <- args
  bad$moment_order <- 2
  expect_error(do.call(wm_fit, bad), "moment_order")
  bad <- args
  bad$mesh_exponent <- 1 / 3
  expect_error(do.call(wm_fit, bad), "strict sufficient rate interval")
  bad <- args
  bad$method <- "stabilized"
  bad$rho_bounds <- c(0.1, 10)
  expect_error(do.call(wm_fit, bad), "supplied fixed fold_id")
  bad$fold_id <- dat$fold
  bad$rho_bounds <- NULL
  expect_error(do.call(wm_fit, bad), "rho_bounds")
})
