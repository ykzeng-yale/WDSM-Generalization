# SPDX-License-Identifier: GPL-3.0-only
# Weighted Matching application statistics; extracted from the validated
# local application adapter. Function bodies are unchanged.
# WDSM/wdsmatch author credits: Yukang Zeng, Guangyu Tong, Jiaqi Tong,
# Haidong Lu, Bhramar Mukherjee, and Fan Li.
# Sourcing defines functions only; callers supply data and all count columns.

app_capture <- function(action) {
  warnings <- character()
  value <- tryCatch(withCallingHandlers(action(), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = identity)
  if (inherits(value, "error")) return(list(ok = FALSE, error = conditionMessage(value),
    warnings = warnings, attempts = value$attempts, failed_columns = value$failed_columns, partial = value$partial))
  list(ok = TRUE, value = value, warnings = warnings)
}

app_count_columns <- function(a, counts) {
  stopifnot(is.matrix(counts), nrow(counts) == a$n, ncol(counts) == 200L)
  lapply(seq_len(200L), function(b) {
    m <- counts[, b]
    finite_integer <- all(is.finite(m) & m >= 0 & m == floor(m))
    total <- if (finite_integer) sum(m) else NA_real_
    arms <- if (finite_integer) vapply(0:1, function(z) sum(m[a$Z == z]), numeric(1)) else c(NA_real_, NA_real_)
    okay <- finite_integer && total == a$n && all(arms > 0)
    list(column = b, ok = okay, total = total, arm0 = arms[1L], arm1 = arms[2L],
      column_sha256 = digest::digest(m, algo = "sha256"),
      error = if (okay) "" else "Invalid multinomial count/arm mass; retained, never redrawn")
  })
}

app_ps <- function(a, modules, m = rep(1, a$n)) app_capture(function() {
  if (any(!is.finite(m) | m < 0) || sum(m) <= 0 ||
      any(vapply(0:1, function(z) sum(m[a$Z == z]) <= 0, logical(1))))
    stop("Invalid count/arm mass; no PS fit attempted")
  fit <- modules$source$wdsm_case_fit_propensity(stats::reformulate(a$variables, "A"), a$data,
    m * if (a$ps_weighting == "unit") 1 else a$W/mean(a$W))
  list(probability = as.numeric(stats::predict(fit, a$data, type = "response")),
    coefficients = stats::coef(fit), diagnostics = attr(fit, "wdsm_ps_diagnostics"))
})

app_components <- function(a, modules, m = rep(1, a$n), ps = NULL) {
  valid <- length(m) == a$n && all(is.finite(m) & m >= 0) && sum(m) > 0 &&
    all(vapply(0:1, function(z) sum(m[a$Z == z]) > 0, logical(1)))
  if (!valid) {
    failure <- list(ok = FALSE, error = "Invalid count/arm mass; no fit attempted", warnings = character())
    arms <- stats::setNames(list(failure, failure), c("0", "1"))
    return(list(PS = failure, PS_score = failure, PG = arms, DSM = arms,
      linear = arms, full_X_scores = failure, count_validation_failed = TRUE))
  }
  w <- a$W/max(a$W); cw <- a$n * m/sum(m)
  if (is.null(ps)) ps <- app_ps(a, modules, m)
  pg <- linear <- ds <- vector("list", 2L); names(pg) <- names(linear) <- names(ds) <- c("0", "1")
  standardize <- function(S) {
    centered <- sweep(S, 2L, colSums(S * (m/sum(m))), "-")
    scale <- sqrt(colSums(centered^2 * (m/sum(m))))
    if (any(!is.finite(scale) | scale <= 0)) stop("Degenerate pooled score variance.")
    sweep(centered, 2L, scale, "/")
  }
  for (z in 0:1) {
    k <- as.character(z)
    pg[[k]] <- app_capture(function() {
      beta <- modules$common$.wm_ns_ols(a$design, a$Y, cw * (a$Z == z), paste0("pg", z), 1e-10)
      list(mean = as.vector(a$design %*% beta), coefficients = beta)
    })
    linear[[k]] <- app_capture(function() {
      beta <- modules$common$.wm_ns_ols(a$design, a$Y, cw * w * (a$Z == z), paste0("full_X", z), 1e-10)
      list(mean = as.vector(a$design %*% beta), coefficients = beta)
    })
    ds[[k]] <- app_capture(function() {
      if (!ps$ok) stop(ps$error)
      if (!pg[[k]]$ok) stop(pg[[k]]$error)
      scores <- standardize(cbind(ps = ps$value$probability, pg = pg[[k]]$value$mean))
      basis <- modules$common$.wm_ns_basis(scores)
      beta <- modules$common$.wm_ns_ols(basis, a$Y, cw * w * (a$Z == z), paste0("bc", z), 1e-10)
      list(scores = scores, mean = as.vector(basis %*% beta), coefficients = beta)
    })
  }
  ps_score <- app_capture(function() {if (!ps$ok) stop(ps$error); standardize(matrix(ps$value$probability, ncol = 1L))})
  list(PS = ps, PS_score = ps_score, PG = pg, DSM = ds, linear = linear,
    full_X_scores = if (all(m == 1)) app_capture(function() standardize(a$X)) else NULL,
    scope = "Source PS equation; current shared OLS/basis helpers; empirical point/refit components only")
}

app_predictions <- function(components, family, estimand) {
  values <- if (family == "DSM") components$DSM else components$linear
  if (!values[["0"]]$ok) stop(values[["0"]]$error)
  result <- list(mean0 = values[["0"]]$value$mean)
  if (estimand == "PATE") {
    if (!values[["1"]]$ok) stop(values[["1"]]$error)
    result$mean1 <- values[["1"]]$value$mean
  }
  result
}

app_row <- function(a, family, M, estimand, point, replication = NULL) {
  okay <- point$ok && !is.null(replication) && replication$ok
  data.frame(study = a$study, outcome = a$outcome, baseline = a$baseline,
    method = paste0("WM_", family), estimand = estimand, M = M, n = a$n,
    n_treated = sum(a$Z), B = 200L, estimate = if (point$ok) point$value$estimate else NA_real_,
    raw_estimate = if (point$ok) point$value$raw_estimate else NA_real_,
    SE = if (okay) replication$value$se else NA_real_,
    lower = if (okay) replication$value$conf.int[1L] else NA_real_,
    upper = if (okay) replication$value$conf.int[2L] else NA_real_,
    status = if (!point$ok) "failed_point" else if (!okay) "failed_interval" else "empirical_complete",
    error = if (!point$ok) point$error else if (!okay) replication$error else "",
    units = a$units, variance_divisor = "B", inference_scope = "Empirical fixed-reuse; supplied W frozen; sampling assumptions unverified")
}

app_match <- function(a, modules, components, family, M, estimand) {
  predictions <- app_predictions(components, family, estimand)
  if (family == "PS" && !components$PS_score$ok) stop(components$PS_score$error)
  if (family == "full_X" && !components$full_X_scores$ok) stop(components$full_X_scores$error)
  scores0 <- switch(family, PS = components$PS_score$value,
    DSM = components$DSM[["0"]]$value$scores, full_X = components$full_X_scores$value)
  scores1 <- if (estimand == "PATE") switch(family, PS = scores0,
    DSM = components$DSM[["1"]]$value$scores, full_X = scores0) else NULL
  do.call(modules$common$wm_match, c(list(Y = a$Y, Z = a$Z, weights = a$W,
    scores0 = scores0, scores1 = scores1, M = M, estimand = estimand,
    variance = FALSE, tie_rule = "source_random", tie_seed = 20260917L + M), predictions))
}

app_replicate <- function(point, modules, counts, refits, family, estimand) {
  kind <- if (family == "DSM") "DSM" else "linear"
  cache <- refits$means[[kind]]
  failed <- which(colSums(!is.finite(cache$mean0)) > 0 |
    (estimand == "PATE" & colSums(!is.finite(cache$mean1)) > 0))
  validation <- app_count_columns(list(n = nrow(counts), Z = point$data$Z), counts)
  failed <- sort(unique(c(failed, which(!vapply(validation, `[[`, logical(1), "ok")))))
  if (length(failed)) stop(structure(list(message = paste("Requested count/prediction columns failed:",
    paste(failed, collapse = ",")), call = NULL, failed_columns = failed),
    class = c("application_column_error", "error", "condition")))
  b <- 0L
  callback <- function(m) {
    b <<- b + 1L; stopifnot(identical(m, counts[, b]))
    p <- list(mean0 = cache$mean0[, b])
    if (estimand == "PATE") p$mean1 <- cache$mean1[, b]
    p
  }
  value <- modules$common$wm_bootstrap_refit(point, counts, callback, conf.level = .95)
  stopifnot(b == 200L, value$variance_divisor == "B")
  value
}

