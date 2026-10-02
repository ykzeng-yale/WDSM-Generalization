# Private source-specification comparison engine; review receipts pin versions.
# Does not alter WM theory or claim design-based validity for WM intervals.
# Dependencies: survey, MatchIt, data.table, sandwich, and source-loaded WM R code.

lenis_capture <- function(fun) {
  warnings <- character()
  began <- proc.time()[[3L]]
  value <- tryCatch(withCallingHandlers(fun(), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = identity)
  list(ok = !inherits(value, "error"),
       value = if (inherits(value, "error")) NULL else value,
       error = if (inherits(value, "error")) conditionMessage(value) else "",
       condition_class = if (inherits(value, "error")) class(value) else character(),
       solver_attempts = if (inherits(value, "error")) value$attempts else NULL,
       partial_results = if (inherits(value, "error")) value$partial_results else NULL,
       warnings = warnings, elapsed = proc.time()[[3L]] - began)
}

lenis_response_data <- function(selected, mechanism) {
  stopifnot(mechanism %in% c("No", "MAR", "MARX", "MART"),
            all(c("id", "Prob", "z", "y", "ones", "Cluster", "Strata",
                  paste0("x", 1:7), paste0("r", 1:4)) %in% names(selected)),
            !anyDuplicated(selected$id))
  d <- selected
  d$s.wt <- 1/d$Prob
  response_fit <- NULL
  probability <- rep(1, nrow(d))
  response_name <- c(No = "r1", MAR = "r2", MARX = "r3", MART = "r4")[[mechanism]]
  if (mechanism != "No") {
    covariates <- c(paste0("x", 1:6),
      if (mechanism == "MARX") "x7" else if (mechanism == "MART") "z")
    formula <- stats::reformulate(covariates, response_name)
    design <- survey::svydesign(ids = ~1, data = d)
    response_fit <- survey::svyglm(formula, design, family = "binomial")
    probability <- as.numeric(stats::predict(response_fit, type = "response"))
    stopifnot(all(is.finite(probability)), all(probability > 0 & probability <= 1))
    d$ps.r <- probability
    d$s.wt <- d$s.wt/probability
  }
  # The original diagnostic takes this SD before respondent filtering.
  sd_weight <- stats::sd(d$s.wt)
  rows <- which(d[[response_name]] == 1)
  list(data = d[rows, , drop = FALSE], selected_rows = rows,
       selected_ids = selected$id, response_name = response_name,
       response_probability = probability, response_fit = response_fit,
       sd_weight_before_filtering = sd_weight, mechanism = mechanism)
}

lenis_ps_match <- function(d, ps_method) {
  formula <- stats::reformulate(c(paste0("x", 1:6),
    if (ps_method == "C") "s.wt"), "z")
  design <- if (ps_method == "W")
    survey::svydesign(ids = ~Cluster, strata = ~Strata, weights = ~s.wt, data = d) else
    survey::svydesign(ids = ~1, data = d)
  fit <- survey::svyglm(formula, design, family = "quasibinomial")
  d$ps.w <- as.numeric(stats::predict(fit, type = "response"))
  match <- MatchIt::matchit(formula, data = d, method = "nearest",
                           ratio = 1, distance = d$ps.w, replace = FALSE)
  stopifnot(is.matrix(match$match.matrix), ncol(match$match.matrix) == 1L)
  # Historical transfer code expects data.frame(match.matrix)$X1.
  if (is.null(colnames(match$match.matrix))) colnames(match$match.matrix) <- "1"
  list(data = d, ps_fit = fit, match = match,
       diagnostics = list(converged = isTRUE(fit$converged),
                          boundary = isTRUE(fit$boundary), rank = fit$rank))
}

lenis_matched_data <- function(ps, transfer) {
  d <- MatchIt::match_data(ps$match, data = ps$data)
  if (transfer == "WT") {
    d$id <- rownames(d)
    matched <- data.frame(ps$match$match.matrix)
    matched$id <- rownames(matched)
    d <- merge(d, matched, by = "id", all.x = TRUE)
    d$X1[d$z == 0] <- d$id[d$z == 0]
    # Preserve the original WT merge/order/lead convention.
    dt <- data.table::data.table(d)
    dt <- dt[order(dt$X1, dt$z), ]
    dt[, s.wt2 := data.table::shift(s.wt, 1, type = "lead"), by = X1]
    d <- as.data.frame(dt)
    d$s.wt[d$z == 0] <- d$s.wt2[d$z == 0]
  }
  stopifnot(all(is.finite(d$s.wt)), all(d$s.wt > 0),
            !anyDuplicated(as.character(d$id)))
  d
}

lenis_native_case <- function(response, target) {
  d <- response$data
  ps <- lapply(c("U", "W", "C"), function(k) lenis_capture(function() lenis_ps_match(d, k)))
  names(ps) <- c("U", "W", "C")
  settings <- data.frame(ps = c("U", "U", "U", "W", "W", "C", "C"),
    om = c("U", rep("W", 6)), transfer = c("OW", "WT", "OW", "WT", "OW", "WT", "OW"))
  results <- vector("list", nrow(settings))
  for (j in seq_len(nrow(settings))) {
    s <- settings[j, ]; id <- paste(s$ps, "PS", s$om, "OM", s$transfer, sep = "_")
    results[[j]] <- lenis_capture(function() {
      if (!ps[[s$ps]]$ok) stop(ps[[s$ps]]$error)
      md <- lenis_matched_data(ps[[s$ps]]$value, s$transfer)
      design <- if (s$om == "U") survey::svydesign(ids = ~1, data = md) else
        survey::svydesign(ids = ~Cluster, strata = ~Strata, weights = ~s.wt, data = md)
      formulas <- list(A = stats::reformulate(c(paste0("x", 1:6), "z"), "y"),
                       U = y ~ z)
      fits <- lapply(formulas, function(f) survey::svyglm(f, design))
      method <- paste(paste0(s$ps, ".PS"), paste0(s$om, ".OM"),
                      s$transfer, response$mechanism, sep = "|")
      balance_weights <- if (s$ps == "U" && s$om == "U") md$ones else md$s.wt
      smd <- vapply(paste0("x", 1:6), function(v)
        diff(tapply(md[[v]] * balance_weights, md$z, sum)/
               tapply(balance_weights, md$z, sum)), numeric(1))
      smd <- c(smd, diff(tapply(md$s.wt * md$ones, md$z, sum)/
                        tapply(md$ones, md$z, sum))/response$sd_weight_before_filtering)
      names(smd) <- paste("SMD", c(paste0("x", 1:6), "s.wt"), method, sep = ".")
      values <- smd
      for (adjustment in names(fits)) {
        fit <- fits[[adjustment]]
        est <- unname(stats::coef(fit)["z"])
        se <- unname(sqrt(diag(stats::vcov(fit)))["z"])
        triple <- c(est, se, as.numeric(est - 1.96 * se < target &
                                        est + 1.96 * se > target))
        names(triple) <- paste(c("ATT", "SE", "Cov"), adjustment, method, sep = "|")
        values <- c(values, triple)
      }
      result <- list(values = values, matched_data = md, fits = fits, specification = s,
           diagnostics = list(ps = ps[[s$ps]]$value$diagnostics,
             outcome_converged = vapply(fits, function(f) isTRUE(f$converged), logical(1)),
             outcome_boundary = vapply(fits, function(f) isTRUE(f$boundary), logical(1))))
      lenis_validate_native(result)
    })
    names(results)[j] <- id
  }
  list(results = results, ps = ps, settings = settings)
}

lenis_validate_native <- function(result) {
  if (any(!is.finite(result$values))) {
    # Preserve the source arithmetic and diagnostic fit state; never substitute
    # a different native estimator merely to obtain a finite result.
    stop(structure(list(message = "Native method returned nonfinite results.",
      call = NULL, partial_results = result),
      class = c("lenis_nonfinite_native", "error", "condition")))
  }
  result
}

# Use the declared weighted logistic-score root for the corrected quick
# comparator. Supplied weights are rescaled only inside the equivalent score
# equation; matching and outcome analysis retain the original supplied weights.
# Explicit starting values and controls reproduce the accepted correction.
lenis_quick_weighted_ps <- function(d) {
  diagnostics <- NULL
  result <- lenis_capture(function() {
    formula <- stats::reformulate(c(paste0("x", 1:6), "s.wt"), "z")
    X <- stats::model.matrix(formula, d); y <- d$z; w <- d$s.wt
    stopifnot(nrow(X) == nrow(d), all(is.finite(X)), all(y %in% 0:1),
      length(unique(y)) == 2L, all(is.finite(w)), all(w > 0))
    fit <- stats::glm(formula, data = d, weights = w/mean(w),
      family = stats::quasibinomial("logit"), start = rep(0, ncol(X)),
      control = stats::glm.control(epsilon = 1e-12, maxit = 100L),
      model = FALSE, x = FALSE, y = FALSE)
    beta <- stats::coef(fit)
    diagnostics <<- list(coefficients = beta, converged = fit$converged,
      boundary = fit$boundary, rank = fit$rank, iter = fit$iter,
      reported_deviance = fit$deviance,
      start = rep(0, ncol(X)), solver_weight_divisor = mean(w),
      epsilon = 1e-12, maxit = 100L, formula = deparse(formula))
    stopifnot(fit$rank == ncol(X), all(is.finite(beta)))
    eta <- as.vector(X %*% beta); p <- stats::plogis(eta)
    # Stable objective uses eta; it does not reuse clipped GLM deviance.
    loss <- pmax(eta, 0) - y * eta + log1p(exp(-abs(eta)))
    raw_weight_loglik <- -sum(w * loss)
    stable_loss_per_raw_weight <- -raw_weight_loglik/sum(w)
    column_scale <- sqrt(colSums(X^2 * w)/sum(w))
    stopifnot(all(is.finite(column_scale)), all(column_scale > 0))
    D <- sweep(X, 2L, column_scale, "/")
    score <- as.vector(crossprod(X, w * (y - p)))/sum(w)
    normalized_score <- max(abs(score)/pmax(colSums(abs(X) * w)/sum(w),
      .Machine$double.eps))
    information <- crossprod(D, D * (w * p * (1 - p)))/sum(w)
    information_eigenvalues <- eigen(information, symmetric = TRUE,
      only.values = TRUE)$values
    positive_information <- all(is.finite(information_eigenvalues)) &&
      min(information_eigenvalues) > 0
    newton <- if (positive_information) tryCatch(as.vector(D %*%
      solve(information, crossprod(D, w * (y - p))/sum(w))),
      error = function(e) rep(NA_real_, nrow(d))) else rep(NA_real_, nrow(d))
    boundary <- p <= .Machine$double.eps | p >= 1 - .Machine$double.eps
    diagnostics <<- c(diagnostics, list(linear_predictors = eta,
      probabilities = p, raw_weight_loglik = raw_weight_loglik,
      stable_loss_per_raw_weight = stable_loss_per_raw_weight,
      score_per_raw_weight = setNames(score, colnames(X)),
      normalized_score_max = normalized_score,
      score_normalization = "max |sum w Xj (Z-p)| / max(sum w |Xj|, eps sum w)",
      column_scale = column_scale, scaled_information = information,
      information_eigenvalues = information_eigenvalues,
      max_abs_newton_linear_predictor_step = max(abs(newton)),
      probability_range = range(p), boundary_count = sum(boundary),
      boundary_fraction = mean(boundary), distinct_probability_count = length(unique(p))))
    failures <- c(if (!isTRUE(fit$converged)) "GLM did not converge",
      if (!all(is.finite(c(eta, p, loss, score)))) "Nonfinite PS arithmetic",
      if (!all(is.finite(c(raw_weight_loglik, stable_loss_per_raw_weight, fit$deviance))))
        "Nonfinite aggregate objective or reported deviance",
      if (!positive_information) "Numerical information is not positive definite",
      if (!is.finite(normalized_score) || normalized_score > 1e-10)
        "Normalized weighted-score residual exceeds 1e-10",
      if (!all(is.finite(newton)) || max(abs(newton)) > 1e-8)
        "Newton-predicted linear-predictor change exceeds 1e-8",
      if (all(boundary)) "All fitted probabilities are at machine boundaries")
    if (length(failures)) stop(paste(failures, collapse = "; "))
    list(probability = p, formula = formula)
  })
  result$diagnostics <- diagnostics
  result
}

lenis_quick_case <- function(response, weighted_ps = TRUE) {
  d <- response$data
  formula <- stats::reformulate(c(paste0("x", 1:6), "s.wt"), "z")
  if (weighted_ps) {
    ps <- lenis_quick_weighted_ps(d)
    for (message in ps$warnings) warning(message, call. = FALSE)
    if (!ps$ok) stop(ps$error, call. = FALSE)
    match <- MatchIt::matchit(formula, data = d, method = "quick",
      distance = ps$value$probability, estimand = "ATT", s.weights = ~s.wt)
  } else {
    match <- MatchIt::matchit(formula, data = d, method = "quick",
      distance = "glm", estimand = "ATT")
    match <- MatchIt::add_s.weights(match, ~s.wt, data = d)
  }
  md <- MatchIt::match_data(match, data = d)
  expected <- match$weights[rownames(md)] * match$s.weights[rownames(md)]
  stopifnot(max(abs(md$weights - expected)) < 1e-10 * max(1, abs(expected)))
  formula_y <- y ~ z * (x1 + x2 + x3 + x4 + x5 + x6)
  fit <- stats::lm(formula_y, data = md, weights = weights)
  treated <- md[md$z == 1, , drop = FALSE]
  p1 <- p0 <- treated
  p1$z <- 1; p0$z <- 0
  design1 <- stats::model.matrix(stats::delete.response(stats::terms(fit)), p1)
  design0 <- stats::model.matrix(stats::delete.response(stats::terms(fit)), p0)
  contrast <- colSums((design1 - design0) * (treated$s.wt/sum(treated$s.wt)))
  covariance <- sandwich::vcovCL(fit, cluster = md$subclass, type = "HC1")
  estimate <- sum(contrast * stats::coef(fit))
  variance <- as.numeric(crossprod(contrast, covariance %*% contrast))
  stopifnot(all(is.finite(c(estimate, variance))), variance >= 0)
  list(estimate = estimate, variance = variance, se = sqrt(variance),
       conf.int = estimate + c(-1, 1) * stats::qnorm(.975) * sqrt(variance),
       contrast = contrast, match = match, matched_data = md,
       fit = fit, weighted_ps = weighted_ps,
       propensity = if (weighted_ps) ps else NULL,
       method = "Adapted official MatchIt generalized-full ATT procedure")
}

lenis_full_x_prediction <- function(d, multiplicity = rep(1, nrow(d))) {
  X <- cbind(intercept = 1, as.matrix(d[paste0("x", 1:6)]))
  w <- multiplicity * (d$s.wt/max(d$s.wt)) * (1 - d$z)
  fit <- stats::lm.wfit(X, d$y, w = w, tol = 1e-10)
  if (fit$rank != ncol(X) || any(!is.finite(fit$coefficients))) {
    stop("Full-X weighted control regression has no finite full-rank fit.")
  }
  list(mean0 = as.vector(X %*% fit$coefficients),
       coefficients = fit$coefficients, design = X)
}

# PATT outcomes may vary on treated rows only. Scores, fitted control means
# and donor graphs then remain unchanged in the Lenis six-multiplier design.
lenis_patt_outcome_batch <- function(fit, Y, counts = NULL, refit_means = NULL) {
  stopifnot(inherits(fit, "wm_match"), fit$estimand == "PATT",
            fit$method == "self_normalized", is.matrix(Y), nrow(Y) == fit$n,
            all(is.finite(Y)))
  z <- fit$data$Z; controls <- z == 0
  stopifnot(all(Y[controls, , drop = FALSE] == fit$data$Y[controls]))
  w <- fit$analysis_weights; K <- fit$loads$incoming[, 1L]
  a <- z * w - (1 - z) * K
  point <- as.vector(crossprod(a, Y - fit$predictions$mean0)/sum(z * w))
  raw <- as.vector(crossprod(a, Y)/sum(z * w))
  if (is.null(counts)) return(list(estimate = point, raw_estimate = raw, status = "point_only"))
  stopifnot(is.matrix(counts), nrow(counts) == fit$n, ncol(counts) >= 2,
            all(is.finite(counts)), all(counts >= 0 & counts == floor(counts)),
            all(colSums(counts) == fit$n),
            is.matrix(refit_means), identical(dim(refit_means), dim(counts)))
  denominator <- as.vector(crossprod(counts, z * w))
  failed <- which(denominator <= 0 |
    colSums(counts[controls, , drop = FALSE]) <= 0 |
    colSums(!is.finite(refit_means)) > 0)
  valid <- setdiff(seq_len(ncol(counts)), failed)
  draws <- matrix(NA_real_, ncol(counts), ncol(Y),
                  dimnames = list(colnames(counts), colnames(Y)))
  if (length(valid)) {
    numerator <- crossprod(counts[, valid, drop = FALSE], Y * a) -
      colSums(counts[, valid, drop = FALSE] * refit_means[, valid, drop = FALSE] * a)
    draws[valid, ] <- numerator/denominator[valid]
  }
  variance <- if (length(failed)) setNames(rep(NA_real_, ncol(Y)), colnames(Y)) else
    colMeans(sweep(draws, 2L, colMeans(draws), "-")^2)
  stopifnot(all(is.finite(c(point, raw, draws[valid, , drop = FALSE]))))
  list(estimate = point, raw_estimate = raw, draws = draws, variance = variance,
       se = sqrt(variance), lower = point - stats::qnorm(.975) * sqrt(variance),
       upper = point + stats::qnorm(.975) * sqrt(variance),
       variance_divisor = "B", B = ncol(counts), n = fit$n,
       status = if (length(failed)) "failed_interval" else "ok",
       failed_columns = failed, completed_columns = valid,
       inference_scope = "Fixed-reuse prediction refit; supplied response weights frozen; empirical external-design comparison")
}

lenis_refit_predictions <- function(d, wm, counts) {
  n <- nrow(d)
  stopifnot(is.matrix(counts), nrow(counts) == n, ncol(counts) >= 2,
            all(is.finite(counts)), all(counts >= 0 & counts == floor(counts)),
            all(colSums(counts) == n))
  D <- cbind(intercept = 1, as.matrix(d[paste0("x", 1:6)]))
  means <- list(full_x = matrix(NA_real_, n, ncol(counts)),
                DSM = matrix(NA_real_, n, ncol(counts)))
  diagnostics <- vector("list", ncol(counts))
  arguments <- list(Y = d$y, Z = d$z, weights = d$s.wt, ps_design = D,
    pg0_design = D, estimand = "PATT", ps_weighting = "probability")
  for (b in seq_len(ncol(counts))) {
    fx <- lenis_capture(function() {
      value <- lenis_full_x_prediction(d, counts[, b])
      stopifnot(length(value$mean0) == n, all(is.finite(value$mean0)))
      value
    })
    if (fx$ok) means$full_x[, b] <- fx$value$mean0
    fx$value <- NULL
    arguments$multiplicity <- counts[, b]
    ds <- lenis_capture(function() {
      value <- wm$.wm_wdsm_fit_stack(arguments,
        wm$.wm_wdsm_prediction_stack, wm$.wm_wdsm_controls(list()))
      if (length(value$stack$mean0) != n || any(!is.finite(value$stack$mean0))) {
        stop(structure(list(message = "DSM refit returned nonfinite predictions.",
          call = NULL, attempts = value$attempts),
          class = c("lenis_nonfinite_refit", "error", "condition")))
      }
      value
    })
    if (ds$ok) {
      means$DSM[, b] <- ds$value$stack$mean0
      ds$solver_attempts <- ds$value$attempts
      ds$refined <- ds$value$refined
    }
    ds$value <- NULL
    diagnostics[[b]] <- list(column = b,
      count_sha256 = digest::digest(counts[, b], algo = "sha256"),
      full_x = fx, DSM = ds)
  }
  list(means = means, diagnostics = diagnostics, B = ncol(counts),
       counts_sha256 = digest::digest(counts, algo = "sha256"))
}

lenis_wm_case <- function(response, wm, M = c(1L, 3L, 5L), counts = NULL) {
  d <- response$data; n <- nrow(d)
  D <- cbind(intercept = 1, as.matrix(d[paste0("x", 1:6)]))
  full_x <- lenis_full_x_prediction(d)
  base <- wm$wm_wdsm_fit(d$y, d$z, d$s.wt, D, D,
    M = max(M), estimand = "PATT", inference = "none", ps_weighting = "probability")
  # Reuse exactly the wrapper's fitted PS. An affine scalar standardization
  # leaves ordering unchanged; full X uses pooled variance with divisor n.
  X <- sweep(D[, -1, drop = FALSE], 2L, colMeans(D[, -1, drop = FALSE]), "-")
  X <- sweep(X, 2L, sqrt(colMeans(X^2)), "/")
  scores <- list(PS = base$nuisance$scores0[, 1L, drop = FALSE],
                 DSM = base$nuisance$scores0, X6 = X)
  means <- list(PS = full_x$mean0, DSM = base$nuisance$mean0, X6 = full_x$mean0)
  fits <- list()
  for (family in names(scores)) for (donors in M) {
    key <- paste(family, donors, sep = "_M")
    fits[[key]] <- if (family == "DSM" && donors == max(M)) base$fit else
      wm$wm_match(d$y, d$z, d$s.wt, scores[[family]], M = donors,
                  estimand = "PATT", mean0 = means[[family]], variance = FALSE)
  }
  refits <- list(full_x = NULL, DSM = NULL); diagnostics <- list()
  if (!is.null(counts)) {
    result <- lenis_refit_predictions(d, wm, counts)
    refits <- result$means
    diagnostics <- result$diagnostics
  }
  list(fits = fits, base = base, full_x = full_x, refit_means = refits,
       refit_diagnostics = diagnostics, response = response$mechanism,
       correction = c(PS = "weighted linear full-X", DSM = "source weighted quadratic double-score",
                      X6 = "weighted linear full-X"),
       scope = "External Lenis comparison; no new WM theorem or estimated-response-weight correction asserted")
}
