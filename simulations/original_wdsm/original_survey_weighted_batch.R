# Original-study integration draft. Independent review/pilot required before use.
# No population generation, sampling, count generation or method substitution.
ows_capture <- function(fun) {
  warnings <- character()
  started <- proc.time()[[3L]]
  x <- tryCatch(withCallingHandlers(fun(), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = identity)
  list(ok = !inherits(x, "error"), value = if (!inherits(x, "error")) x else NULL,
    error = if (inherits(x, "error")) conditionMessage(x) else "",
    solver_attempts = if (inherits(x, "error")) x$attempts else NULL,
    warnings = warnings, elapsed = proc.time()[[3L]] - started)
}

ows_standardize <- function(x) {
  x <- as.matrix(x); centered <- sweep(x, 2L, colMeans(x), "-")
  scale <- sqrt(colMeans(centered^2))
  stopifnot(all(is.finite(scale)), all(scale > 0))
  sweep(centered, 2L, scale, "/")
}

ows_quadratic <- function(x) {
  x <- as.matrix(x); out <- cbind(intercept = 1, x)
  for (j in seq_len(ncol(x))) for (k in j:ncol(x)) {
    out <- cbind(out, x[, j] * x[, k])
    colnames(out)[ncol(out)] <- paste0(colnames(x)[j], ":", colnames(x)[k])
  }
  stopifnot(ncol(out) == 1L + ncol(x) + ncol(x) * (ncol(x) + 1L) / 2L)
  out
}

ows_full_x_predictions <- function(data, design, multiplicity, estimand) {
  n <- nrow(data)
  stopifnot(nrow(design) == n, length(multiplicity) == n,
    all(is.finite(c(design, multiplicity))), all(multiplicity >= 0))
  output <- list(mean0 = NULL, mean1 = NULL, coefficients = list())
  for (arm in if (estimand == "PATE") 0:1 else 0L) {
    rows <- which(data$A == arm)
    weights <- data$survey_weight[rows] * multiplicity[rows]
    stopifnot(all(is.finite(weights)), sum(weights) > 0)
    fit <- stats::lm.wfit(design[rows, , drop = FALSE], data$Y[rows], weights)
    stopifnot(fit$rank == ncol(design), all(is.finite(fit$coefficients)))
    prediction <- as.vector(design %*% fit$coefficients)
    stopifnot(length(prediction) == n, all(is.finite(prediction)))
    output[[paste0("mean", arm)]] <- prediction
    output$coefficients[[paste0("arm", arm)]] <- fit$coefficients
  }
  output
}

# Algebra independent of the original source's replicate routine. The modern
# fit supplies its actual normalized weights and M-specific graph loads.
ows_wm_replicate <- function(fit, multiplicity, mean0, mean1 = NULL) {
  Y <- fit$data$Y; Z <- as.numeric(fit$data$Z); w <- fit$analysis_weights
  K <- fit$loads$incoming[cbind(seq_along(Y), Z + 1L)]
  stopifnot(length(multiplicity) == length(Y), length(mean0) == length(Y),
    all(is.finite(c(Y, Z, w, K, multiplicity, mean0))), all(multiplicity >= 0))
  if (fit$estimand == "PATE") {
    stopifnot(length(mean1) == length(Y), all(is.finite(mean1)))
    denominator <- sum(multiplicity * w)
    numerator <- sum(multiplicity * w * (mean1 - mean0)) +
      sum(multiplicity * (2 * Z - 1) * (w + K) * (Y - ifelse(Z == 1, mean1, mean0)))
  } else {
    denominator <- sum(multiplicity * Z * w)
    numerator <- sum(multiplicity * Z * w * (Y - mean0)) -
      sum(multiplicity * (1 - Z) * K * (Y - mean0))
  }
  stopifnot(is.finite(denominator), denominator > 0, is.finite(numerator))
  numerator / denominator
}

ows_count_summary <- function(point, draws, errors, alpha = .05) {
  stopifnot(length(draws) >= 2L, length(errors) == length(draws))
  okay <- isTRUE(point$ok) && all(is.finite(draws)) && all(errors == "")
  estimate <- if (isTRUE(point$ok)) point$value$estimate else NA_real_
  variance <- if (okay) mean((draws - mean(draws))^2) else NA_real_
  interval <- if (okay) estimate + c(-1, 1) * stats::qnorm(1 - alpha / 2) * sqrt(variance) else
    c(lower = NA_real_, upper = NA_real_)
  arithmetic_ok <- okay && is.finite(estimate) && is.finite(variance) && variance >= 0 &&
    all(is.finite(interval))
  list(estimate = estimate, draws = draws, errors = errors, B = length(draws),
    variance = if (arithmetic_ok) variance else NA_real_, variance_divisor = "B",
    conf.int = if (arithmetic_ok) interval else c(lower = NA_real_, upper = NA_real_),
    interval_error = if (!isTRUE(point$ok)) "Original point unavailable" else
      if (any(errors != "")) "One or more prescribed count fits failed" else
      if (any(!is.finite(draws))) "One or more prescribed draws nonfinite" else
      if (!arithmetic_ok) "Nonfinite interval arithmetic" else "",
    point_status = if (isTRUE(point$ok)) "ok" else "failed",
    interval_status = if (arithmetic_ok) "ok" else "failed",
    scope = "Source-design empirical fixed-reuse refit interval; all prescribed draws retained")
}

ows_weighted_batch <- function(data, specification, estimand, sampling_design,
                               replicate, counts, wm, original, M = c(1L, 3L, 5L)) {
  stopifnot(specification %in% c("CorCor", "CorMis", "MisCor", "MisMis"),
    estimand %in% c("PATE", "PATT"), sampling_design %in% c("retrospective", "prospective"),
    is.matrix(counts), nrow(counts) == nrow(data), ncol(counts) >= 2L,
    all(is.finite(counts)), all(counts >= 0 & counts == floor(counts)),
    all(colSums(counts) == nrow(data)), identical(as.integer(M), c(1L, 3L, 5L)))
  n <- nrow(data); B <- ncol(counts); Y <- data$Y; Z <- data$A; w <- data$survey_weight
  ps_terms <- c(paste0("X", 1:6), if (specification %in% c("CorCor", "CorMis")) "X1:X2")
  pg_terms <- c(paste0("X", 1:6), if (specification %in% c("CorCor", "MisCor")) "X1:X2")
  ps_formula <- stats::reformulate(ps_terms, "A")
  pg_formula <- stats::reformulate(pg_terms, "Y")
  Dps <- stats::model.matrix(stats::reformulate(ps_terms), data)
  Dpg <- stats::model.matrix(stats::reformulate(pg_terms), data)
  X6 <- ows_standardize(data[paste0("X", 1:6)]); D6 <- ows_quadratic(X6)
  arguments <- list(Y = Y, Z = Z, weights = w, ps_design = Dps, pg0_design = Dpg,
    pg1_design = if (estimand == "PATE") Dpg else NULL, estimand = estimand,
    ps_weighting = if (sampling_design == "retrospective") "probability" else "unit")
  controls <- wm$.wm_wdsm_controls(list())
  fit_modern <- function(m) wm$.wm_wdsm_fit_stack(c(arguments, list(multiplicity = m)),
    wm$.wm_wdsm_prediction_stack, controls)
  fit_original <- function(m) original$wdsm_final_nuisance(data, ps_formula, pg_formula,
    m, estimand, sampling_design)
  bases <- list(original = ows_capture(function() fit_original(rep(1, n))),
    DSM = ows_capture(function() fit_modern(rep(1, n))),
    PS_means = ows_capture(function() ows_full_x_predictions(data, Dpg, rep(1, n), estimand)),
    X6_means = ows_capture(function() ows_full_x_predictions(data, D6, rep(1, n), estimand)),
    PS = ows_capture(function() {
      fitted <- original$wdsm_final_fit_ps(ps_formula, data,
        if (sampling_design == "retrospective") w else rep(1, n))
      list(scores = ows_standardize(cbind(ps = fitted$probability)),
           probability = fitted$probability, diagnostics = fitted$diagnostics)
    }))
  points <- list()
  require_value <- function(capture) {
    if (!isTRUE(capture$ok)) stop(capture$error, call. = FALSE)
    capture$value
  }
  for (donors in M) {
    points[[paste0("WDSM_M", donors)]] <- ows_capture(function() {
      q <- require_value(bases$original)
      graph <- original$wdsm_make_matches(Z, q$D0, q$D1, donors, estimand,
        3000000L + 100L * replicate + donors, "random", 64 * .Machine$double.eps)
      K <- original$wdsm_final_reuse(Z, w, graph$matches_0, graph$matches_1, estimand)
      estimate <- original$wdsm_final_replicate(Y, Z, w, K, rep(1, n), q$q0, q$q1, estimand)
      raw <- original$wdsm_final_replicate(Y, Z, w, K, rep(1, n), numeric(n),
        if (estimand == "PATE") numeric(n) else NULL, estimand)
      list(estimate = estimate, raw_estimate = raw, K = K, graph = graph, M = donors)
    })
    for (family in c("PS", "DSM", "X6")) {
      key <- paste0("WM_", family, "_M", donors)
      points[[key]] <- ows_capture(function() {
        if (family == "DSM") {
          q <- require_value(bases$DSM)$stack; score0 <- q$scores0; score1 <- q$scores1
        } else {
          q <- require_value(bases[[paste0(family, "_means")]])
          score0 <- if (family == "PS") require_value(bases$PS)$scores else X6
          score1 <- if (estimand == "PATE") score0 else NULL
        }
        fit <- wm$wm_match(Y, Z, w, score0, scores1 = score1, M = donors, estimand = estimand,
                          mean0 = q$mean0, mean1 = q$mean1, variance = FALSE)
        reconstructed <- ows_wm_replicate(fit, rep(1, n), q$mean0, q$mean1)
        stopifnot(abs(reconstructed - fit$estimate) <= 1e-10 * (1 + abs(fit$estimate)))
        fit
      })
    }
  }
  draws <- matrix(NA_real_, B, length(points), dimnames = list(NULL, names(points)))
  errors <- matrix("Not evaluated", B, length(points), dimnames = dimnames(draws))
  diagnostics <- vector("list", B)
  for (b in seq_len(B)) {
    m <- counts[, b]
    refits <- list(original = ows_capture(function() fit_original(m)),
      DSM = ows_capture(function() fit_modern(m)),
      PS = ows_capture(function() ows_full_x_predictions(data, Dpg, m, estimand)),
      X6 = ows_capture(function() ows_full_x_predictions(data, D6, m, estimand)))
    for (key in names(points)) {
      result <- ows_capture(function() {
        point <- require_value(points[[key]])
        if (startsWith(key, "WDSM_")) {
          q <- require_value(refits$original)
          original$wdsm_final_replicate(Y, Z, w, point$K, m, q$q0, q$q1, estimand)
        } else {
          family <- sub("_M[0-9]+$", "", sub("^WM_", "", key))
          q <- require_value(refits[[family]])
          if (family == "DSM") q <- q$stack
          ows_wm_replicate(point, m, q$mean0, q$mean1)
        }
      })
      if (result$ok && length(result$value) == 1L && is.finite(result$value)) {
        draws[b, key] <- result$value; errors[b, key] <- ""
      } else errors[b, key] <- if (result$ok) "Nonfinite scalar refit draw" else result$error
    }
    # Keep compact fit evidence; do not duplicate all n-by-B prediction matrices.
    compact <- lapply(names(refits), function(family) {
      result <- refits[[family]]
      if (result$ok) {
        q <- result$value
        if (family == "DSM") {
          result$solver_attempts <- q$attempts
          result$parameter <- q$stack$parameter; q <- q$stack
        }
        result$prediction_sha256 <- digest::digest(if (family == "original")
          list(mean0 = q$q0, mean1 = q$q1) else list(mean0 = q$mean0, mean1 = q$mean1), algo = "sha256")
        result$coefficients <- q$coefficients
        if (family == "original") result$ps_diagnostics <- q$ps_diagnostics
      }
      result$value <- NULL
      result
    })
    names(compact) <- names(refits)
    diagnostics[[b]] <- list(column = b, count_sha256 = digest::digest(m, algo = "sha256"), fits = compact)
  }
  summaries <- lapply(names(points), function(key) ows_count_summary(points[[key]], draws[, key], errors[, key]))
  names(summaries) <- names(points)
  list(points = points, bases = bases, summaries = summaries, refit_diagnostics = diagnostics,
    input_data_sha256 = digest::digest(data, algo = "sha256"),
    counts_sha256 = digest::digest(counts, algo = "sha256"), B = B, n = n,
    estimand = estimand, specification = specification, sampling_design = sampling_design,
    declared_corrections = c(PS = "scenario full-X weighted regression", DSM = "source own-map quadratic",
      X6 = "complete quadratic in standardized X6"),
    review_status = "DRAFT_REQUIRES_INDEPENDENT_INTEGRATION_REVIEW")
}
