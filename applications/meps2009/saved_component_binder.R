# SPDX-License-Identifier: GPL-3.0-only
# MEPS saved-only input binding. No model fitter, graph constructor or RNG.
# The common package evaluates inference; this file only assembles the retained
# model equations/derivatives using existing pure evaluators and saved inputs.
meps_cc_value <- function(x, label) {
  if (!is.list(x) || !isTRUE(x$ok) || !is.list(x$value))
    stop("Unavailable original component: ", label)
  x$value
}

meps_cc_bind <- function(a, components, fit, family, wm, saved_stack) {
  stopifnot(family %in% c("PS", "DSM", "full_X"), is.environment(wm),
    inherits(fit, "wm_match"), identical(fit$method, "self_normalized"),
    fit$M == 3L, fit$estimand %in% c("PATE", "PATT"), fit$n == a$n,
    identical(as.numeric(fit$data$Y), as.numeric(a$Y)),
    identical(as.numeric(fit$data$Z), as.numeric(a$Z)),
    identical(as.numeric(fit$weights), as.numeric(a$W)),
    identical(fit$graph$tie_rule, "Exact distance, then original row index"))
  pate <- fit$estimand == "PATE"; used <- if (pate) 0:1 else 0L
  n <- a$n; X <- a$X; D <- a$design; pD <- ncol(D)
  stopifnot(nrow(D) == n, nrow(X) == n, pD == ncol(X) + 1L,
    identical(colnames(D), c("(Intercept)", colnames(X))),
    identical(as.numeric(D[, 1L]), rep(1, n)),
    identical(as.numeric(D[, -1L, drop = FALSE]), as.numeric(X)))
  agree <- wm$.wm_reciprocal_agree
  bind <- function(actual, expected, label) {
    stopifnot(length(actual) == length(expected),
      all(is.finite(actual)), all(is.finite(expected)))
    agree(as.numeric(actual), as.numeric(expected), label)
  }
  ps <- if (family != "full_X") meps_cc_value(components$PS, "PS") else NULL
  if (!is.null(ps)) {
    stopifnot(identical(names(ps$coefficients), colnames(D)),
      length(ps$probability) == n, all(ps$probability > 0 & ps$probability < 1))
    bind(stats::plogis(as.vector(D %*% ps$coefficients)), ps$probability,
         "Saved no-offset propensity coefficient/probability binding")
  }
  if (family != "DSM") {
    labels <- paste0("mean", used)
    models <- stats::setNames(lapply(used, function(z) {
      model <- meps_cc_value(components$linear[[as.character(z)]], paste0("linear", z))
      stopifnot(length(model$coefficients) == pD)
      bind(as.vector(D %*% model$coefficients), model$mean, "Saved raw linear prediction")
      bind(model$mean, fit$predictions[[paste0("mean", z)]], "Original matching prediction")
      list(design = D, coefficients = model$coefficients)
    }), labels)
    data <- list(Y = a$Y, Z = a$Z, weights = a$W)
    if (family == "PS") {
      stack <- saved_stack(data, fit, models, family = "PS", ps = list(
        design = D, probability = ps$probability, weighting = "probability",
        coefficients = ps$coefficients), wm = wm)
    } else {
      # Exact complete degree-one affine coordinate conversion, not a refit.
      center <- colMeans(X); variance <- colMeans(sweep(X, 2L, center, "-")^2)
      stopifnot(all(is.finite(variance)), all(variance > 0))
      standardized <- sweep(sweep(X, 2L, center, "-"), 2L, sqrt(variance), "/")
      bind(standardized, components$full_X_scores$value, "Saved full-X coordinate binding")
      basis_names <- c("intercept", colnames(X))
      exponents <- rbind(rep(0, ncol(X)), diag(ncol(X)))
      dimnames(exponents) <- list(basis_names, colnames(X))
      std_design <- cbind(intercept = 1, standardized)
      for (label in labels) {
        alpha <- models[[label]]$coefficients
        delta <- c(alpha[1L] + sum(center * alpha[-1L]), sqrt(variance) * alpha[-1L])
        names(delta) <- basis_names
        bind(as.vector(std_design %*% delta), fit$predictions[[label]],
             "Degree-one affine prediction binding")
        models[[label]] <- list(design = std_design, coefficients = delta)
      }
      stack <- saved_stack(data, fit, models, family = "full_X", raw_X = X,
                           polynomial_exponents = exponents, wm = wm)
    }
  } else {
    # Build only the input record of the existing pure evaluator. Never call
    # .wm_wdsm_nuisance_stack: that constructor would perform new model fits.
    inputs <- list(Y = a$Y, Z = a$Z, weights = a$W, ps_design = D,
                   ps_offset = rep(0, n))
    parameters <- character(); k <- list()
    block <- function(prefix, labels) {
      at <- length(parameters) + seq_along(labels)
      parameters <<- c(parameters, paste0(prefix, ":", labels)); at
    }
    k$ps <- block("ps", colnames(D))
    for (z in used) {
      label <- paste0("pg", z)
      inputs[[paste0(label, "_design")]] <- D
      inputs[[paste0(label, "_offset")]] <- rep(0, n)
      k[[label]] <- block(label, colnames(D))
    }
    raw_names <- c("ps", paste0("pg", used))
    k$center <- stats::setNames(block("center", raw_names), raw_names)
    k$variance <- stats::setNames(block("variance", raw_names), raw_names)
    basis_names <- c("intercept", "ps", "pg", "ps_sq", "ps_pg", "pg_sq")
    for (z in used) k[[paste0("bc", z)]] <- block(paste0("bc", z), basis_names)
    theta <- stats::setNames(numeric(length(parameters)), parameters)
    theta[k$ps] <- ps$coefficients
    raw <- matrix(0, n, length(raw_names), dimnames = list(NULL, raw_names))
    raw[, "ps"] <- ps$probability
    for (z in used) {
      label <- paste0("pg", z); pg <- meps_cc_value(components$PG[[as.character(z)]], label)
      stopifnot(length(pg$coefficients) == pD)
      theta[k[[label]]] <- pg$coefficients
      raw[, label] <- as.vector(D %*% pg$coefficients)
      bind(raw[, label], pg$mean, "Saved unweighted PG prediction")
      bc <- meps_cc_value(components$DSM[[as.character(z)]], paste0("DSM", z))
      stopifnot(length(bc$coefficients) == 6L)
      theta[k[[paste0("bc", z)]]] <- bc$coefficients
    }
    theta[k$center] <- colMeans(raw)
    theta[k$variance] <- colMeans(sweep(raw, 2L, theta[k$center], "-")^2)
    stopifnot(all(theta[k$variance] > 0))
    w <- a$W / max(a$W)
    shell <- structure(list(parameter = theta, parameter_names = parameters,
      blocks = k, used_arms = used, estimand = fit$estimand, inputs = inputs,
      multiplicity = rep(1, n), empirical_probability = rep(1 / n, n),
      weight_scale = max(a$W), scaled_weights = w,
      ps_weighting = "probability", scaled_ps_weights = w),
      class = c(".wm_wdsm_nuisance_stack", "list"))
    evaluated <- wm$.wm_wdsm_nuisance_evaluate(shell)
    for (z in used) {
      arm <- as.character(z); score_label <- paste0("scores", z); mean_label <- paste0("mean", z)
      bind(evaluated[[score_label]], components$DSM[[arm]]$value$scores,
           "Saved DSM generated-coordinate binding")
      bind(evaluated[[score_label]], fit$graph[[score_label]], "Original DSM graph coordinates")
      bind(evaluated[[mean_label]], components$DSM[[arm]]$value$mean,
           "Saved DSM complete-quadratic prediction binding")
      bind(evaluated[[mean_label]], fit$predictions[[mean_label]], "Original DSM point predictions")
    }
    J <- evaluated$jacobian; psi <- evaluated$estimating_equations
    influence <- -t(solve(J, t(psi)))
    colnames(influence) <- parameters
    centered <- sweep(influence, 2L, colMeans(influence), "-")
    stack <- c(evaluated, list(parameter_names = parameters, blocks = k,
      nuisance_influence = influence, nuisance_covariance = crossprod(centered) / n,
      inference_arguments = list(fit = fit, nuisance_influence = influence,
        mean_derivative0 = evaluated$mean0_derivative,
        mean_derivative1 = if (pate) evaluated$mean1_derivative else NULL,
        covariance_scope = if (pate) "full_x" else "patt")))
    if (!pate) stack$inference_arguments$transport0 <- list(mode = "zero", basis = "current_centering")
  }
  parameters <- stack$parameter_names; p <- length(parameters)
  expected <- switch(family, PS = pD + 2L + length(used) * pD,
    full_X = 2L * ncol(X) + length(used) * pD,
    DSM = pD + length(used) * pD + 2L * (1L + length(used)) + 6L * length(used))
  stopifnot(p == expected, length(unique(parameters)) == p,
    identical(colnames(stack$estimating_equations), parameters),
    identical(dimnames(stack$jacobian), list(parameters, parameters)),
    identical(colnames(stack$nuisance_influence), parameters),
    identical(stack$inference_arguments$fit, fit))
  for (x in list(stack$estimating_equations, stack$jacobian, stack$nuisance_influence,
                 stack$nuisance_covariance, stack$inference_arguments$mean_derivative0,
                 stack$inference_arguments$mean_derivative1))
    if (!is.null(x)) stopifnot(all(is.finite(x)))
  for (label in c("mean_derivative0", if (pate) "mean_derivative1")) {
    derivative <- stack$inference_arguments[[label]]
    stopifnot(identical(dim(derivative), c(as.integer(n), as.integer(p))),
              identical(colnames(derivative), parameters))
  }
  stack$inference_arguments$weight_derivative <- matrix(0, n, p,
    dimnames = list(NULL, parameters))
  stack$diagnostics_saved_only <- list(parameter_dimension = p,
    parameter_names = parameters, blocks = stack$blocks,
    max_mean_equation = max(abs(colMeans(stack$estimating_equations))),
    mean_equation = stats::setNames(colMeans(stack$estimating_equations), parameters),
    max_influence_mean = max(abs(colMeans(stack$nuisance_influence))),
    jacobian_reciprocal_condition = rcond(stack$jacobian),
    conditional_population_premises_verified = FALSE, original_root_localization_verified = FALSE,
    original_predictions_and_scores_bound = TRUE, original_point_graph_retained = TRUE)
  stack
}
