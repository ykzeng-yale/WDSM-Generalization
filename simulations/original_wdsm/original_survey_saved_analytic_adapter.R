# Companion analytic rows from ORIGINAL SAVED PS/X6 points and nuisance fits.
# Source simulations/wm_saved_full_x_stack.R explicitly before this file, or
# pass that helper as saved_stack. It is not an exported package function.
# No fit, matching, count generation, bootstrap or RNG call occurs here.
# This adapter does not load or replace the frozen original DSM adapter/results.
ows_saved_analytic_batch <- function(data, weighted, wm,
                                     saved_stack = wm_saved_full_x_stack,
                                     alpha = .05, M = NULL, authenticated_input = NULL) {
  if (!is.data.frame(data) || !all(c("Y", "A", "survey_weight", paste0("X", 1:6)) %in% names(data)) ||
      !is.list(weighted) || !identical(weighted$n, nrow(data)) ||
      !weighted$estimand %in% c("PATE", "PATT") ||
      !weighted$sampling_design %in% c("retrospective", "prospective") ||
      !weighted$specification %in% c("CorCor", "CorMis", "MisCor", "MisMis") ||
      !is.environment(wm) || !is.function(wm$wm_fitted_inference) ||
      !is.function(saved_stack) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("Supply original data, weighted batch, common package environment and saved-stack helper.", call. = FALSE)
  }
  if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required.", call. = FALSE)
  decoded_hash <- digest::digest(data, algo = "sha256")
  if (is.null(authenticated_input)) {
    if (!identical(weighted$input_data_sha256, decoded_hash)) {
      stop("Original data do not match the recorded weighted-batch data hash.", call. = FALSE)
    }
  } else {
    # The caller must have authenticated the original enclosing plan/shard/sample
    # and file/receipt hashes with the frozen lossless reader. Version-2 lossless
    # restoration guarantees strict object identity, not canonical serialization.
    # The fresh decoded hash binds this in-memory object only; it never replaces
    # the historical context/weighted hash or the original artifact identity.
    bound <- authenticated_input
    is_sha <- function(x) is.character(x) && length(x) == 1L && !is.na(x) &&
      grepl("^[0-9a-f]{64}$", x)
    required <- c("data", "context", "input_reference", "receipt", "provenance", "decoded_data_sha256")
    if (!is.list(bound) || !setequal(names(bound), required) || anyDuplicated(names(bound)) ||
        !is.data.frame(bound$data) || !identical(data, bound$data, attrib.as.set = FALSE) ||
        !identical(decoded_hash, bound$decoded_data_sha256) || !is.list(bound$context) ||
        !is.list(bound$input_reference) || !is.list(bound$receipt) || !is.list(bound$provenance)) {
      stop("Authenticated input binding does not identify the complete supplied data.", call. = FALSE)
    }
    ctx <- bound$context; ref <- bound$input_reference; receipt <- bound$receipt
    provenance <- bound$provenance
    chain_fields <- c("original_plan_sha256", "shard_completion_sha256", "sample_completion_sha256")
    if (!is_sha(ctx$input_data_sha256) || !identical(weighted$input_data_sha256, ctx$input_data_sha256) ||
        !identical(ctx$sample_n, nrow(data)) || !identical(ctx$sampling_design, weighted$sampling_design) ||
        !identical(ref$directory, "input") || !is_sha(ref$file_sha256) || !is_sha(ref$receipt_sha256) ||
        !identical(receipt$status, "ARTIFACT_COMMITTED") ||
        !identical(receipt$format, "ows_lossless_atomic_arrays_v2") ||
        !identical(receipt$file, "artifact.rds") || !identical(receipt$file_sha256, ref$file_sha256) ||
        !isTRUE(receipt$exact_roundtrip) || !isTRUE(receipt$attribute_order_preserved) ||
        !identical(provenance$schema, "original_ows_authenticated_input_v1") ||
        !all(vapply(chain_fields, function(k) is_sha(provenance[[k]]), logical(1))) ||
        !identical(provenance$decoder_sha256, "22ac24acdf0012d6735e30c07f4080aa401e7c01d40822083c2c2c06123ac218") ||
        !identical(provenance$reader_sha256, "8304cf8a333acc51f5fd3351d5bd1c362844c8703ffa3bd8b38af6ba7e6fb9c9")) {
      stop("Authenticated input provenance or historical weighted/context identity is inconsistent.", call. = FALSE)
    }
  }
  capture <- function(fun) {
    warnings <- character()
    value <- tryCatch(withCallingHandlers(fun(), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = identity)
    list(ok = !inherits(value, "error"),
      value = if (!inherits(value, "error")) value else NULL,
      error = if (inherits(value, "error")) conditionMessage(value) else "",
      warnings = warnings)
  }
  require_value <- function(x, label) {
    if (!is.list(x) || !isTRUE(x$ok) || is.null(x$value)) {
      reason <- if (is.list(x) && is.character(x$error) && length(x$error) == 1L) x$error else ""
      stop(label, " unavailable", if (nzchar(reason)) paste0(": ", reason), call. = FALSE)
    }
    x$value
  }
  recorded <- names(weighted$points)
  if (is.null(recorded) || anyDuplicated(recorded)) stop("Recorded point names must be unique.", call. = FALSE)
  keys <- recorded[grepl("^WM_(PS|X6)_M[1-9][0-9]*$", recorded)]
  recorded_M <- sort(unique(as.integer(sub(".*_M", "", keys))))
  if (!length(recorded_M) || anyNA(recorded_M)) stop("No valid recorded PS/X6 matching orders.", call. = FALSE)
  if (is.null(M)) M <- recorded_M else {
    if (!is.numeric(M) || is.complex(M) || !is.null(dim(M)) || !length(M) ||
        any(!is.finite(M)) || any(M < 1 | M != floor(M)) || anyDuplicated(M) ||
        any(!M %in% recorded_M)) stop("Requested M must contain unique positive recorded matching orders.", call. = FALSE)
    M <- as.integer(M)
  }
  pate <- identical(weighted$estimand, "PATE")
  labels <- if (pate) c("mean0", "mean1") else "mean0"
  source_data <- list(Y = data$Y, Z = data$A, weights = data$survey_weight)
  ps_terms <- c(paste0("X", 1:6), if (weighted$specification %in% c("CorCor", "CorMis")) "X1:X2")
  pg_terms <- c(paste0("X", 1:6), if (weighted$specification %in% c("CorCor", "MisCor")) "X1:X2")
  Dps <- stats::model.matrix(stats::reformulate(ps_terms), data)
  Dpg <- stats::model.matrix(stats::reformulate(pg_terms), data)
  raw_X <- as.matrix(data[paste0("X", 1:6)])
  quadratic <- function(x) {
    out <- cbind(intercept = 1, x)
    for (j in seq_len(ncol(x))) for (k in j:ncol(x)) {
      out <- cbind(out, x[, j] * x[, k])
      colnames(out)[ncol(out)] <- paste0(colnames(x)[j], ":", colnames(x)[k])
    }
    out
  }
  bind_fit <- function(fit, reference, donors) {
    if (!inherits(fit, "wm_match") || !identical(as.integer(fit$M), donors) ||
        !identical(fit$estimand, weighted$estimand)) stop("Recorded point has wrong M/estimand/class.", call. = FALSE)
    fields <- c("data", "weights", "weight_scale", "analysis_weights", "predictions", "estimand", "method", "n")
    if (!identical(fit[fields], reference[fields]) ||
        !identical(fit$graph$scores0, reference$graph$scores0) ||
        !identical(fit$graph$scores1, reference$graph$scores1)) {
      stop("Recorded M-specific fit does not share the original family data/maps/predictions/weights.", call. = FALSE)
    }
    invisible(TRUE)
  }
  scope <- paste("Companion empirical analytic comparison using the original corrected point and graph.",
    "Correct full-X mean, identification, complete regular influence, applicable current-graph",
    "sampling/moment and numerical-root conditions are substantive premises, not verified here.",
    "Source sampling designs and misspecified means are not certified by this calculation.",
    "No cluster/stratum correction; supplied probability weights held fixed.",
    "No automatic agreement with original fixed-reuse refit variance or finite-sample coverage claim.")
  rows <- assemblies <- nuisance <- list()
  for (family in c("PS", "X6")) {
    point_keys <- paste0("WM_", family, "_M", M)
    successful <- point_keys[vapply(point_keys, function(key) isTRUE(weighted$points[[key]]$ok), logical(1))]
    reference <- if (length(successful)) weighted$points[[successful[1L]]]$value else NULL
    solved <- capture(function() {
      if (is.null(reference)) stop("No successful recorded family point.", call. = FALSE)
      means <- require_value(weighted$bases[[paste0(family, "_means")]], "Recorded mean models")
      if (family == "PS") {
        ps <- require_value(weighted$bases$PS, "Recorded propensity scores")
        design <- Dpg
      } else {
        standardized <- reference$graph$scores0
        if (ncol(standardized) != ncol(raw_X)) stop("Original X6 dimension mismatch.", call. = FALSE)
        colnames(standardized) <- colnames(raw_X)
        design <- quadratic(standardized)
      }
      models <- stats::setNames(lapply(labels, function(label) list(design = design,
        coefficients = means$coefficients[[if (label == "mean0") "arm0" else "arm1"]])), labels)
      if (family == "PS") saved_stack(source_data, reference, models, family = "PS",
        ps = list(design = Dps, probability = ps$probability,
          weighting = if (weighted$sampling_design == "retrospective") "probability" else "unit"), wm = wm)
      else saved_stack(source_data, reference, models, family = "full_X", raw_X = raw_X, wm = wm)
    })
    nuisance[[family]] <- solved
    for (donors in M) {
      point_key <- paste0("WM_", family, "_M", donors)
      key <- paste0("WM_", family, "_analytic_M", donors)
      point <- weighted$points[[point_key]]
      result <- capture(function() {
        fit <- require_value(point, "Recorded matching point")
        stack <- require_value(solved, "Complete saved nuisance stack")
        bind_fit(fit, reference, donors)
        arguments <- stack$inference_arguments
        arguments$fit <- fit
        arguments$conf.level <- 1 - alpha
        do.call(wm$wm_fitted_inference, arguments)
      })
      assemblies[[key]] <- result
      inf <- if (result$ok) result$value else NULL
      point_ok <- isTRUE(point$ok) && is.list(point$value) &&
        length(point$value$estimate) == 1L && is.finite(point$value$estimate)
      aligned <- !is.null(inf) && point_ok && identical(inf$fit, point$value) &&
        identical(inf$estimate, point$value$estimate)
      okay <- aligned && isTRUE(inf$available) && length(inf$variance) == 1L &&
        is.finite(inf$variance) && inf$variance > 0 && length(inf$conf.int) == 2L && all(is.finite(inf$conf.int))
      reason <- if (!result$ok) result$error else if (!aligned) "Analytic result did not retain the recorded point and fit" else
        if (!okay) paste("Analytic variance or interval unavailable:", inf$unavailable_reason) else ""
      refit <- weighted$summaries[[point_key]]
      refit_variance <- if (is.list(refit) && length(refit$variance) == 1L && is.finite(refit$variance)) refit$variance else NA_real_
      rows[[key]] <- list(estimate = if (point_ok) point$value$estimate else NA_real_,
        point_status = if (point_ok) "ok" else "failed",
        interval_status = if (okay) "ok" else "failed", interval_error = reason,
        variance = if (okay) inf$variance else NA_real_,
        conf.int = if (okay) inf$conf.int else c(lower = NA_real_, upper = NA_real_),
        variance_method = "complete_saved_stack_full_X_adjusted_analytic",
        point_source = point_key, aligned_with_recorded_point = aligned,
        saved_refit_variance = refit_variance,
        saved_refit_interval_status = if (is.list(refit)) refit$interval_status else "missing",
        saved_refit_interval_error = if (is.list(refit)) refit$interval_error else "Recorded refit summary missing",
        saved_refit_B = if (is.list(refit)) refit$B else NA_integer_,
        saved_refit_variance_divisor = if (is.list(refit)) refit$variance_divisor else NA_character_,
        analytic_to_saved_refit_variance = if (okay && identical(refit$interval_status, "ok") &&
          is.finite(refit_variance) && refit_variance > 0 && is.finite(inf$variance / refit_variance)) inf$variance / refit_variance else NA_real_,
        assumptions_verified = FALSE, original_refit_limit_agreement_declared = FALSE, scope = scope)
    }
  }
  list(rows = rows, assemblies = assemblies, nuisance = nuisance, M = M,
    n = weighted$n, estimand = weighted$estimand, specification = weighted$specification,
    sampling_design = weighted$sampling_design, input_data_sha256 = weighted$input_data_sha256,
    scope = scope, existing_results_modified = FALSE)
}
