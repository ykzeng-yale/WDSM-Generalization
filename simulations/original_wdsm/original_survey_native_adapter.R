# In-memory observation of the independently sourced original comparator.
# No package namespace or cached original source is modified.
# Source original_survey_weighted_batch.R for ows_capture first.
ows_replace_namespace_calls <- function(expression, replacements) {
  if (missing(expression)) return(quote(expr = ))
  if (!is.call(expression)) return(expression)
  head <- expression[[1L]]
  key <- if (is.call(head) && identical(head[[1L]], as.name("::")))
    paste(as.character(head[[2L]]), as.character(head[[3L]]), sep = "::") else ""
  parts <- lapply(as.list(expression), ows_replace_namespace_calls, replacements = replacements)
  if (key %in% names(replacements)) parts[[1L]] <- as.name(replacements[[key]])
  result <- as.call(parts)
  attributes(result) <- attributes(expression)
  result
}

ows_native_comparator <- function(data, method, estimand, sampling_design,
                                  specification, B, original, repo, alpha = .05) {
  stopifnot(method %in% c("PSM", "PGM", "DSM", "SWPSM"), estimand %in% c("PATE", "PATT"),
    sampling_design %in% c("retrospective", "prospective"),
    specification %in% c("CorCor", "CorMis", "MisCor", "MisMis"), B >= 2L, B == as.integer(B))
  state <- new.env(parent = emptyenv())
  state$matching <- list(); state$bootstrap <- NULL; state$point <- NA_real_
  state$source_ledger <- NULL; state$full_match <- NULL; state$matched_data <- NULL
  state$contrast <- NULL; state$source_warnings <- list()
  sha_object <- function(x) digest::digest(x, algo = "sha256")
  rng_state <- function() if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  before <- rng_state()

  record_match <- function(...) {
    args <- list(...)
    result <- do.call(Matching::Match, args)
    # Keep actual input geometry and ordered donor edges, including the two
    # directional Match calls used by native ATE. Do not reinterpret its metric.
    state$matching[[length(state$matching) + 1L]] <- list(
      score = args$X, treatment = args$Tr, outcome = args$Y,
      settings = args[setdiff(names(args), c("X", "Tr", "Y"))],
      index_treated = result$index.treated, index_control = result$index.control,
      weights = result$weights, estimate = result$est,
      matched_outcomes = result$mdata,
      score_sha256 = sha_object(args$X))
    result
  }

  record_boot <- function(data, statistic, R, ...) {
    frame <- parent.frame()
    point_name <- switch(method, DSM = "est.ds", PSM = "est.ps", PGM = "est.pg")
    stopifnot(exists(point_name, envir = frame, inherits = FALSE))
    state$point <- as.numeric(get(point_name, envir = frame, inherits = FALSE))
    stopifnot(length(state$point) == 1L, is.finite(state$point), R == B)
    calls <- list(); original_statistic <- statistic
    wrapped_statistic <- function(data, indices, ...) {
      id <- length(calls)  # boot first evaluates t0, then each prescribed draw.
      result <- ows_capture(function() original_statistic(data, indices, ...))
      scalar_ok <- result$ok && length(result$value) == 1L && is.finite(result$value)
      if (result$ok && !scalar_ok) {result$ok <- FALSE; result$error <- "Nonfinite native scalar statistic"}
      returned <- if (scalar_ok) as.numeric(result$value) else NA_real_
      result$value <- NULL
      result$call_index <- id; result$estimate <- returned
      result$indices_sha256 <- sha_object(indices)
      result$counts_sha256 <- sha_object(tabulate(indices, nbins = nrow(data)))
      calls[[id + 1L]] <<- result
      # Continue the already prescribed index draws after a model error. The
      # unchanged source variance guard subsequently refuses any failed draw.
      returned
    }
    bootstrap <- ows_capture(function() boot::boot(data = data, statistic = wrapped_statistic, R = R, ...))
    state$bootstrap <- list(input_data = data, input_data_sha256 = sha_object(data),
      requested = as.integer(R), calls = calls,
      status = bootstrap[c("ok", "error", "warnings", "elapsed")])
    if (!bootstrap$ok) stop(bootstrap$error, call. = FALSE)
    result <- bootstrap$value
    stopifnot(length(calls) == B + 1L, is.matrix(result$t), nrow(result$t) == B, ncol(result$t) == 1L)
    state$bootstrap$seed <- result$seed
    state$bootstrap$t0 <- result$t0; state$bootstrap$draws <- as.numeric(result$t[, 1L])
    state$bootstrap$sim <- result$sim; state$bootstrap$stype <- result$stype
    if (any(!vapply(calls, function(x) isTRUE(x$ok), logical(1))))
      stop("Native bootstrap original-sample statistic or prescribed draw failed", call. = FALSE)
    result
  }

  record_matchit <- function(...) {
    args <- list(...); result <- do.call(MatchIt::matchit, args)
    state$full_match <- list(weights = result$weights, sampling_weights = result$s.weights,
      subclass = result$subclass, distance = result$distance, discarded = result$discarded,
      estimand = result$estimand, method = "full",
      ps_coefficients = if (!is.null(result$model)) stats::coef(result$model) else NULL)
    result
  }
  record_match_data <- function(...) {
    result <- MatchIt::match.data(...)
    state$matched_data <- result[, intersect(c("Y", "A", paste0("X", 1:6), "sw", "weights", "subclass"),
                                             names(result)), drop = FALSE]
    result
  }
  record_contrast <- function(model, newdata, variables, vcov, wts, conf_level) {
    stopifnot(identical(variables, "A"), identical(wts, "sw"),
      all(newdata$A %in% 0:1), all(is.finite(newdata$sw)), all(newdata$sw > 0))
    # Preserve the finite point if a later subclass covariance computation fails.
    # This is the same interacted-linear-model contrast, evaluated directly.
    d1 <- d0 <- newdata; d1$A <- 1; d0$A <- 0
    terms <- stats::delete.response(stats::terms(model))
    X1 <- stats::model.matrix(terms, d1, contrasts.arg = model$contrasts, xlev = model$xlevels)
    X0 <- stats::model.matrix(terms, d0, contrasts.arg = model$contrasts, xlev = model$xlevels)
    gradient <- colSums((X1 - X0) * (newdata$sw / sum(newdata$sw)))
    state$point <- sum(gradient * stats::coef(model))
    stopifnot(is.finite(state$point))
    state$contrast <- list(point = state$point, gradient = gradient,
      coefficients = stats::coef(model), target_rows = rownames(newdata),
      target_sampling_weights = newdata$sw, covariance = "subclass sandwich")
    answer <- marginaleffects::avg_comparisons(model, newdata = newdata, variables = variables,
      vcov = vcov, wts = wts, conf_level = conf_level)
    stopifnot(abs(as.numeric(answer$estimate) - state$point) <= 1e-9 * (1 + abs(state$point)))
    answer
  }

  adapter <- new.env(parent = original)
  # The source handler muffles its own warnings. Observe its message conversion
  # so those warnings also survive an error before the source returns a result.
  adapter$conditionMessage <- function(condition) {
    message <- base::conditionMessage(condition)
    if (inherits(condition, "warning")) state$source_warnings[[length(state$source_warnings)+1L]] <-
      list(message=message, class=class(condition))
    message
  }
  adapter$.ows_record_matchit <- record_matchit
  adapter$.ows_record_match_data <- record_match_data
  adapter$.ows_record_contrast <- record_contrast
  adapter$.comparator_final_vendor <- function(repo, B, ledger) {
    state$source_ledger <- ledger
    vendor <- original$.comparator_final_vendor(repo, B, ledger)
    vendor$.ows_record_match <- record_match
    vendor$.ows_record_boot <- record_boot
    for (name in c("dsmatchATE", "dsmatchATT")) {
      fun <- vendor[[name]]
      body(fun) <- ows_replace_namespace_calls(body(fun), c("Matching::Match" = ".ows_record_match",
                                                          "boot::boot" = ".ows_record_boot"))
      vendor[[name]] <- fun
    }
    vendor
  }
  compare <- original$compare_final
  environment(compare) <- adapter
  body(compare) <- ows_replace_namespace_calls(body(compare), c("MatchIt::matchit" = ".ows_record_matchit",
    "MatchIt::match.data" = ".ows_record_match_data", "marginaleffects::avg_comparisons" = ".ows_record_contrast"))
  result <- ows_capture(function() compare(data, method = method, estimand = estimand,
    sampling_design = sampling_design, M = 1L, B = B, alpha = alpha,
    use_interaction_ps = specification %in% c("CorCor", "CorMis"),
    use_interaction_pg = specification %in% c("CorCor", "MisCor"), repo = repo))
  point <- if (result$ok) result$value$estimate else state$point
  point_ok <- is.numeric(point) && length(point) == 1L && is.finite(point)
  native_order <- if (method=="SWPSM") seq_len(nrow(data)) else sort(data$A,index.return=TRUE)$ix
  warning_inventory <- unique(c(result$warnings,
    vapply(state$source_warnings, function(w) w$message, character(1)),
    state$bootstrap$status$warnings,
    unlist(lapply(state$bootstrap$calls, function(x) x$warnings),use.names=FALSE)))
  list(result = result, estimate = if (point_ok) point else NA_real_,
    input_data_sha256 = sha_object(data),
    point_status = if (point_ok) "ok" else "failed",
    interval_status = if (result$ok) "ok" else "failed",
    matching = state$matching, bootstrap = state$bootstrap,
    full_match = state$full_match, matched_data = state$matched_data, contrast = state$contrast,
    source_ledger = if (is.environment(state$source_ledger)) as.list(state$source_ledger) else NULL,
    source_warnings=state$source_warnings, warning_inventory=warning_inventory,
    native_row_to_input_row=native_order,
    native_row_to_id=if ("id" %in% names(data)) data$id[native_order] else NULL,
    rng_before = before, rng_after = rng_state(),
    source_procedure_unchanged_on_success = TRUE,
    failure_extension = "Retain source point and all pre-generated boot index draws; any failed prescribed draw invalidates its interval",
    instrumentation = "Private source-function namespace dispatch only; package namespaces and source files unchanged")
}
