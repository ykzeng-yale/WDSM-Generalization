test_that("PS fitting weights are separated from effect and correction weights", {
  set.seed(29317)
  n <- 320L
  x <- cbind(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
  D <- cbind(intercept = 1, x)
  G <- cbind(D, product = x[, 1] * x[, 2])
  Z <- rbinom(n, 1, plogis(-0.2 + 0.5 * x[, 1] - 0.3 * x[, 2]))
  Y <- 1 + x[, 1] + 0.4 * x[, 2] + Z * (0.8 + x[, 2]) + rnorm(n, sd = 0.4)
  W <- exp(0.5 * Z + 0.4 * x[, 1] - 0.2 * x[, 2])
  args <- list(Y = Y, Z = Z, weights = W, ps_design = D,
               pg0_design = G, pg1_design = G, M = 3L)
  fits <- lapply(c("probability", "unit"), function(mode)
    do.call(wm_wdsm_fit, c(args, list(ps_weighting = mode, inference = "full_x"))))
  names(fits) <- c("probability", "unit")
  default <- do.call(wm_wdsm_fit, args)
  expect_identical(default$estimate, fits$probability$estimate)
  expect_identical(default$nuisance$parameter, fits$probability$nuisance$parameter)
  expect_gt(max(abs(fits$probability$nuisance$raw_scores[, "ps"] -
                    fits$unit$nuisance$raw_scores[, "ps"])), 0.01)
  for (mode in names(fits)) {
    a <- fits[[mode]]; st <- a$nuisance
    expected_weights <- if (mode == "probability") W else rep(1, n)
    independent <- glm.fit(D, Z, weights = expected_weights,
      family = quasibinomial(), start = rep(0, ncol(D)),
      control = glm.control(epsilon = 1e-13, maxit = 200))
    expect_equal(unname(st$parameter[st$blocks$ps]), unname(independent$coefficients), tolerance = 1e-9)
    expect_identical(a$ps_weighting, mode)
    expect_identical(st$inputs$weights, W)
    expect_equal(a$fit$analysis_weights * a$fit$weight_scale, W)
    # Every donor fraction retains its probability weights in either PS mode.
    for (i in c(1L, 17L, 233L)) {
      donors <- a$fit$graph$neighbors[[i]]
      edge <- a$fit$graph$edges[a$fit$graph$edges$query == i, , drop = FALSE]
      expect_equal(edge$share, W[donors] / sum(W[donors]))
    }
    # Check the whole generated-score Jacobian by differentiating its equations.
    evaluator <- WeightedMatching:::.wm_wdsm_nuisance_evaluate
    h <- 1e-5
    numeric_J <- vapply(seq_along(st$parameter), function(j) {
      plus <- minus <- st$parameter
      # Pooled PS variances can be much smaller than one. Scale the numerical
      # difference to the parameter instead of perturbing them by a fixed 1e-5.
      step <- h * max(abs(st$parameter[j]), 1e-4)
      plus[j] <- plus[j] + step; minus[j] <- minus[j] - step
      (evaluator(st, plus)$mean_equation - evaluator(st, minus)$mean_equation) / (2 * step)
    }, numeric(length(st$parameter)))
    expect_equal(unname(st$jacobian), numeric_J, tolerance = 1e-7)
    # A smooth empirical-mass perturbation independently checks the stacked IF.
    v <- sin(seq_len(n)); v <- v - mean(v)
    stack_args <- args[setdiff(names(args), "M")]
    stack_args$ps_weighting <- mode
    fit_stack <- WeightedMatching:::.wm_wdsm_nuisance_stack
    plus <- do.call(fit_stack, c(stack_args, list(multiplicity = 1 + h * v)))
    minus <- do.call(fit_stack, c(stack_args, list(multiplicity = 1 - h * v)))
    derivative <- (plus$parameter - minus$parameter) / (2 * h)
    expected <- colMeans(st$nuisance_influence * v)
    expect_equal(derivative, expected, tolerance = 2e-6)
    expect_equal(plus$scaled_ps_weights, if (mode == "probability") W / max(W) else rep(1, n))
    expect_false(a$assumptions_verified)
  }
  expect_error(do.call(wm_wdsm_fit, c(args, list(ps_weighting = "invalid"))), "arg")
  # With unit probability weights, both choices must give exactly the same stack.
  args$weights <- rep(1, n)
  u <- do.call(wm_wdsm_fit, c(args, list(ps_weighting = "unit")))
  p <- do.call(wm_wdsm_fit, c(args, list(ps_weighting = "probability")))
  expect_identical(u$nuisance$parameter, p$nuisance$parameter)
  expect_identical(u$estimate, p$estimate)
})
