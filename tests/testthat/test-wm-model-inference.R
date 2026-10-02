# Fixed literal arrays only. The constructed retained-object layout exercises
# dispatch contracts, not a fitted population model, regular root or coverage.
# No nuisance model/local polynomial is fitted in this file.
wm_model_handoff_fixture <- function(estimand = "PATE", lp = FALSE) {
  x <- c(.02, .13, .29, .48, .69, .91, .07, .21, .36, .56, .77, .96)
  u <- c(.2, .8, .4, .1, .9, .5, .6, .3, .7, .4, .2, .8)
  z <- rep(0:1, each = 6L)
  y <- c(1, -2, 4, 0, 3, 2, 4, 1, -1, 5, 2, 6)
  w <- 1 + seq_along(x) / 10
  s0 <- cbind(x, u); s1 <- cbind(x^2, 1-u)
  D0 <- cbind(score = x, center = u); D1 <- cbind(score = u, center = x^2)
  IF <- cbind(score = x - .4, center = u - .5)
  parameters <- colnames(IF); n <- length(x)
  J <- -diag(2L); dimnames(J) <- list(parameters, parameters)
  fit <- wm_match(y, z, w, s0, if (estimand == "PATE") s1 else NULL,
    M = 3L, estimand = estimand, mean0 = x,
    mean1 = if (estimand == "PATE") 1 + x^2 else NULL, variance = FALSE)
  used <- if (estimand == "PATE") c("0", "1") else "0"
  deriv <- list()
  for (a in used) deriv[[a]] <- list(D0, D1)
  stack <- structure(list(parameter = stats::setNames(c(.1, .2), parameters),
    parameter_names = parameters, estimand = estimand, multiplicity = rep(1, n),
    empirical_probability = rep(1 / n, n), estimating_equations = IF, jacobian = J,
    inputs = list(Y = y, Z = z, weights = w), nuisance_influence = IF,
    weight_derivative = matrix(0, n, 2L, dimnames = list(NULL, parameters)),
    scores0 = s0, scores1 = if (estimand == "PATE") s1 else NULL,
    mean0 = fit$predictions$mean0, mean1 = fit$predictions$mean1,
    mean0_derivative = D0, mean1_derivative = if (estimand == "PATE") D1 else NULL,
    score_derivatives = deriv), class = c(".wm_model_nuisance_stack", "list"))
  result <- structure(list(estimate = fit$estimate, raw_estimate = fit$raw_estimate,
    correction = fit$correction, n = fit$n, M = fit$M, estimand = estimand,
    fit = fit, nuisance = stack, assembly = list(original_marker = "retained"),
    inference = list(method = "none"), solver = list(original_marker = "retained"),
    model_recipe = list(original_marker = "retained"), assumptions_verified = FALSE,
    bootstrap_requested = FALSE), class = c("wm_model_fit", "list"))
  if (lp) {
    result$nuisance$score_only <- TRUE
    arms <- list()
    for (a in used) arms[[a]] <- list(mean = fit$predictions[[paste0("mean", a)]],
      reference_mean_derivative = if (a == "0") D0 else D1,
      guard_failure = rep(FALSE, n), derivative_boundary = rep(FALSE, n),
      reference_derivative_failure = rep(FALSE, n),
      reference_derivative_available = rep(TRUE, n))
    result$correction_fit <- list(method = "local_polynomial", arms = arms,
      reference_derivatives_available = TRUE)
    result$nuisance$mean0_derivative <- result$nuisance$mean1_derivative <- NULL
  }
  result
}

test_that("post-fit inference is the exact common result and preserves point/root", {
  for (estimand in c("PATE", "PATT")) for (lp in c(FALSE, TRUE)) {
    object <- wm_model_handoff_fixture(estimand, lp); before <- object
    scope <- if (estimand == "PATE") "full_x" else "patt"
    transport <- if (estimand == "PATT") list(mode = "zero", basis = "current_centering") else NULL
    args <- wdsmatch:::.wm_mi_bind(object)$arguments
    direct <- do.call(wm_fitted_inference, c(args,
      list(covariance_scope = scope, transport0 = transport, conf.level = .9)))
    out <- wm_model_inference(object, scope, transport0 = transport, conf.level = .9)
    expect_identical(object, before)
    expect_identical(out$inference, direct)
    for (name in c("fit", "nuisance", "estimate", "raw_estimate", "correction",
                   "solver", "model_recipe", "n", "M", "estimand")) {
      expect_identical(out[[name]], before[[name]])
    }
    expect_identical(out$conf.int, direct$conf.int)
    expect_null(out$assembly)
    expect_identical(out$inference_handoff$original_assembly, before$assembly)
    expect_identical(out$inference_handoff$original_inference, before$inference)
    expect_false(out$inference_handoff$models_refitted)
    expect_false(out$assumptions_verified)
    again <- wm_model_inference(out, scope, transport0 = transport, conf.level = .9)
    expect_identical(again$inference_handoff$original_assembly, before$assembly)
    expect_identical(again$inference_handoff$original_inference, before$inference)
    expect_identical(again$inference, direct)
  }
})

test_that("unavailable retained inputs stop before low-level inference without refitting", {
  local_mocked_bindings(wm_fitted_inference = function(...) stop("LOW_LEVEL_WAS_CALLED"), .package = "wdsmatch")
  a <- wm_model_handoff_fixture()
  a$nuisance$nuisance_influence <- NULL
  a$nuisance$estimating_equations <- NULL
  out <- wm_model_inference(a, "full_x")
  expect_false(out$inference_handoff$dispatched)
  expect_identical(out$fit, a$fit)
  expect_match(out$inference$unavailable_reason, "neither a complete IF nor both", fixed = TRUE)
  a <- wm_model_handoff_fixture(); a$nuisance$mean0_derivative <- NULL
  out <- wm_model_inference(a, "full_x")
  expect_false(out$inference$available)
  expect_match(out$inference$unavailable_reason, "Missing retained mean_derivative0", fixed = TRUE)
  a <- wm_model_handoff_fixture(lp = TRUE)
  a$correction_fit$reference_derivatives_available <- FALSE
  out <- wm_model_inference(a, "full_x")
  expect_identical(out$fit, a$fit)
  expect_false(out$inference_handoff$dispatched)
  expect_identical(out$inference$failure_stage, "correction")
  a$correction_fit$reference_derivatives_available <- TRUE
  a$correction_fit$arms[["0"]]$guard_failure[1] <- TRUE
  out <- wm_model_inference(a, "full_x")
  expect_false(out$inference_handoff$dispatched)
  a <- wm_model_handoff_fixture()
  a$fit$graph$tie_rule <- "Source-compatible random squared-distance boundary"
  out <- wm_model_inference(a, "full_x")
  expect_identical(out$fit, a$fit)
  expect_identical(out$inference$failure_stage, "tie_rule")
})

test_that("fitted sentinel binds only actual score-coordinate derivatives", {
  a <- wm_model_handoff_fixture(lp = TRUE); parameters <- a$nuisance$parameter_names
  spec <- list(mode = "estimate", representation = "score_measurable_weight",
    tangents = "fitted", cutoff = rep(1, a$n))
  resolved <- wdsmatch:::.wm_mi_transport(spec, a, "0", parameters)
  expect_identical(spec$tangents, "fitted")
  expect_identical(resolved$raw_scores, a$fit$graph$scores0)
  expect_identical(dim(resolved$tangents), c(12L, 2L, 2L))
  expect_identical(dimnames(resolved$tangents)[[3L]], parameters)
  expect_equal(unname(resolved$tangents[, 1L, ]), unname(a$nuisance$score_derivatives[["0"]][[1L]]))
  expect_identical(resolved$cutoff, spec$cutoff)
  expect_null(resolved$bandwidth)
  bad <- spec; bad$representation <- "conditional_weight_chart"
  expect_error(wdsmatch:::.wm_mi_transport(bad, a, "0", parameters), "not a raw-chart derivative", fixed = TRUE)
  bad <- spec; bad$raw_scores <- a$fit$graph$scores0
  expect_error(wdsmatch:::.wm_mi_transport(bad, a, "0", parameters), "Omit raw_scores", fixed = TRUE)
  bad <- spec; bad$raw_to_matching <- function(s) s
  expect_error(wdsmatch:::.wm_mi_transport(bad, a, "0", parameters), "no raw_to_matching", fixed = TRUE)
})

test_that("a caller raw-chart coordinate callback runs exactly once", {
  a <- wm_model_handoff_fixture(); calls <- 0L
  mapping <- function(x) { calls <<- calls + 1L; x }
  B <- wdsmatch:::.wm_mi_tangents(a,"0",a$nuisance$parameter_names)
  spec <- list(mode="estimate", representation="conditional_weight_chart",
    raw_scores=a$fit$graph$scores0, tangents=B, cutoff=rep(1,a$n),
    raw_to_matching=mapping)
  # Isolate dispatch/validation: no transport numerical experiment is needed.
  local_mocked_bindings(wm_graph_transport=function(...) {
    stop("literal unavailable transport for callback-count test")
  }, .package="wdsmatch")
  answer <- wm_model_inference(a,"distinct_rarity",transport0=spec,
    transport1=list(mode="zero",basis="fixed_map"))
  expect_identical(calls,1L)
  expect_identical(answer$fit,a$fit)
  expect_identical(answer$inference$status,"transport_unavailable")
  expect_match(answer$inference$unavailable_reason,"literal unavailable transport",fixed=TRUE)
  expect_true(answer$inference_handoff$dispatched)
})

test_that("API and protected-binding mistakes remain errors", {
  a <- wm_model_handoff_fixture()
  expect_error(wm_model_inference(a), "explicitly", fixed = TRUE)
  expect_error(wm_model_inference(a, "scalar_psm"), "not a generic", fixed = TRUE)
  expect_error(wm_model_inference(a, "patt"), "respective covariance", fixed = TRUE)
  expect_error(wm_model_inference(a, "full_x", transport0 = list(mode="zero",basis="fixed_map")), "omit transport", fixed=TRUE)
  bad <- a; bad$nuisance$inputs$weights[1] <- 99
  expect_error(wm_model_inference(bad, "full_x"), "Nuisance supplied W", fixed = TRUE)
  bad <- a; bad$nuisance$weight_derivative[1,1] <- 1
  expect_error(wm_model_inference(bad, "full_x"), "exact-zero", fixed = TRUE)
  bad <- a; colnames(bad$nuisance$nuisance_influence) <- rev(a$nuisance$parameter_names)
  expect_error(wm_model_inference(bad, "full_x"), "complete n-by-p", fixed = TRUE)
  bad <- a; bad$nuisance$multiplicity[1] <- 2
  expect_error(wm_model_inference(bad, "full_x"), "unit multiplicities", fixed = TRUE)
  patt <- wm_model_handoff_fixture("PATT")
  expect_error(wm_model_inference(patt, "patt", transport1=list(mode="zero",basis="fixed_map")), "only transport0", fixed = TRUE)
})

test_that("only recognized dispatch arithmetic errors become unavailable", {
  a <- wm_model_handoff_fixture()
  local_mocked_bindings(wm_fitted_inference = function(...) stop("Fitted-variance arithmetic exceeded numerical range."), .package="wdsmatch")
  out <- wm_model_inference(a, "full_x")
  expect_identical(out$fit, a$fit)
  expect_false(out$inference$available)
  expect_true(out$inference_handoff$dispatched)
  expect_true(all(is.na(out$conf.int)))
  local_mocked_bindings(wm_fitted_inference = function(...) stop("protected binding changed"), .package="wdsmatch")
  expect_error(wm_model_inference(a, "full_x"), "protected binding changed", fixed = TRUE)
})

test_that("retained equations complete the same full influence once without a fit", {
  a <- wm_model_handoff_fixture(); full <- a$nuisance$nuisance_influence
  a$nuisance$nuisance_influence <- NULL; before <- a
  out <- wm_model_inference(a, "full_x")
  expect_identical(a, before)
  expect_equal(out$nuisance$nuisance_influence, full, tolerance = 1e-14)
  centered <- sweep(full, 2L, colMeans(full), "-")
  expect_equal(out$nuisance$nuisance_covariance, crossprod(centered) / nrow(full), tolerance=1e-14)
  for (key in names(before$nuisance)) expect_identical(out$nuisance[[key]], before$nuisance[[key]])
  expect_false(out$nuisance$influence_completion$models_refitted)
  expect_true(out$nuisance$influence_completion$planned_elements > 0)
  again <- wm_model_inference(out, "full_x")
  expect_identical(again$nuisance, out$nuisance)
  expect_identical(again$inference, out$inference)
  capped <- wm_model_inference(a, "full_x", max_influence_elements=1)
  expect_false(capped$inference_handoff$dispatched)
  expect_null(capped$nuisance$nuisance_influence)
  expect_match(capped$inference$unavailable_reason,"max_influence_elements",fixed=TRUE)
  singular <- a; singular$nuisance$jacobian[,] <- 0
  unavailable <- wm_model_inference(singular, "full_x")
  expect_false(unavailable$inference_handoff$dispatched)
  expect_identical(unavailable$fit, singular$fit)
  expect_match(unavailable$inference$unavailable_reason, "solve failed", fixed=TRUE)
  wrong <- a; rownames(wrong$nuisance$jacobian) <- rev(colnames(full))
  expect_error(wm_model_inference(wrong,"full_x"),"complete named p-by-p",fixed=TRUE)
  wrong <- a; wrong$nuisance$empirical_probability[1] <- .2
  expect_error(wm_model_inference(wrong,"full_x"),"uniform empirical",fixed=TRUE)
})

test_that("critical fitted sentinel preserves exact common result and direct-variance conventions", {
  a <- wm_model_handoff_fixture(lp=TRUE)
  control <- list(tangents0="fitted", tangents1="fitted", auxiliary_rows=1:6,
    bandwidth=.5, density_bounds=c(.01,10), marginal_mass_floor=.001,
    projection_floor=.01, tangent_cap=10, residual_cap=10, direct_variance="raw_only")
  zero <- list(mode="zero",basis="current_centering")
  marker <- structure(list(estimate=a$estimate, root_n_variance=-.1, variance=-.1/a$n,
    se=NA_real_, conf.int=c(lower=NA_real_,upper=NA_real_), available=FALSE,
    status="literal_raw_direct_result", direct_variance=list(mode="raw_only")),
    class=c("wm_fitted_inference","list"))
  captured <- NULL
  local_mocked_bindings(wm_fitted_inference=function(...) {
    captured <<- list(...); marker
  }, .package="wdsmatch")
  out <- wm_model_inference(a,"critical_planar",zero,zero,control)
  expect_identical(out$inference,marker)
  expect_identical(out$root_n_variance,-.1)
  expect_true(all(is.na(out$conf.int)))
  expect_identical(captured$fit,a$fit)
  expect_identical(captured$nuisance_influence,a$nuisance$nuisance_influence)
  expect_true(all(captured$weight_derivative==0))
  expect_identical(captured$critical_control$tangents0,
    wdsmatch:::.wm_mi_tangents(a,"0",a$nuisance$parameter_names))
  expect_identical(captured$critical_control$tangents1,
    wdsmatch:::.wm_mi_tangents(a,"1",a$nuisance$parameter_names))
  expect_identical(captured$critical_control$direct_variance,"raw_only")
  expect_null(captured$critical_control$support)
  expect_identical(control$tangents0,"fitted")
})
