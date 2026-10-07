# SPDX-License-Identifier: GPL-3.0-only
# Portable MEPS 2009 pipeline. Sourcing defines functions only.
# The original fixed-reuse/refit and complete-contribution operators stay separate.

meps_decoded_sha <- function(path) {
  connection <- gzfile(path, "rb")
  on.exit(close(connection), add = TRUE)
  chunks <- list()
  repeat {
    chunk <- readBin(connection, "raw", n = 1048576L)
    if (!length(chunk)) break
    chunks[[length(chunks) + 1L]] <- chunk
  }
  digest::digest(do.call(c, chunks), algo = "sha256", serialize = FALSE)
}

meps_runtime <- function(release_dir) {
  release_dir <- normalizePath(release_dir, mustWork = TRUE)
  for (package in c("digest", "jsonlite", "WeightedMatching"))
    if (!requireNamespace(package, quietly = TRUE)) stop("Missing dependency: ", package)
  helper <- new.env(parent = globalenv())
  common <- file.path(release_dir, "applications", "common")
  helper_files <- file.path(common, c("case_propensity_solver.R", "application_statistics.R",
    "application_modules.R", "application_pipeline.R", "weighted_balance.R"))
  for (path in helper_files) sys.source(path, envir = helper)
  modules <- helper$wm_application_modules()
  ns <- modules$common
  required <- c("wm_fitted_inference", "wm_bootstrap", ".wm_wdsm_nuisance_evaluate",
                ".wm_reciprocal_agree")
  for (name in required) if (!exists(name, ns, inherits = FALSE))
    stop("Installed WeightedMatching lacks the complete-contribution API: ", name)
  # Bind the normal installed namespace to this public release. No namespace
  # modification, alternate library selection or private runtime overlay occurs.
  source_signatures <- list()
  for (path in sort(list.files(file.path(release_dir, "R"), "[.]R$", full.names = TRUE))) {
    for (node in parse(path, keep.source = FALSE)) {
      if (!is.call(node) || !identical(node[[1L]], as.name("<-")) ||
          !is.symbol(node[[2L]]) || !is.call(node[[3L]]) ||
          !identical(node[[3L]][[1L]], as.name("function"))) next
      name <- as.character(node[[2L]])
      if (name %in% names(source_signatures)) stop("Duplicate public source function: ", name)
      expected <- helper$wm_application_function_identity(eval(node[[3L]], baseenv()))
      if (!exists(name, ns, inherits = FALSE) ||
          !identical(helper$wm_application_function_identity(get(name, ns)), expected))
        stop("Installed WeightedMatching differs from this public release: ", name,
             ". Install the current release before resuming; no fallback is used.")
      source_signatures[[name]] <- expected
    }
  }
  if (!length(source_signatures)) stop("No package R source found in release_dir")
  saved_helpers <- new.env(parent = globalenv())
  stack_path <- file.path(release_dir, "simulations", "wm_saved_full_x_stack.R")
  sys.source(stack_path, envir = saved_helpers)
  meps_dir <- file.path(release_dir, "applications", "meps2009")
  sys.source(file.path(meps_dir, "saved_component_binder.R"), envir = saved_helpers)
  sys.source(file.path(meps_dir, "analysis_spec.R"), envir = saved_helpers)
  sys.source(file.path(meps_dir, "saved_graph_balance.R"), envir = saved_helpers)
  code_paths <- c(helper_files, stack_path,
    file.path(meps_dir, c("meps_analysis.R", "saved_component_binder.R", "analysis_spec.R",
      "saved_graph_balance.R")))
  list(helper = helper, modules = modules, saved_helpers = saved_helpers,
    spec = saved_helpers$meps_analysis_spec(),
    identity = list(R_version = R.version.string, R_platform = R.version$platform,
      package = helper$wm_application_package_identity(ns),
      source_function_sha256 = source_signatures,
      application_source_sha256 = stats::setNames(lapply(code_paths, meps_sha),
        c(basename(helper_files), "wm_saved_full_x_stack.R", "meps_analysis.R",
          "saved_component_binder.R", "analysis_spec.R", "saved_graph_balance.R")),
      digest_version = as.character(utils::packageVersion("digest")),
      jsonlite_version = as.character(utils::packageVersion("jsonlite"))))
}

meps_input <- function(prepared_dir, spec, pair) {
  entry <- spec$cohorts[[pair]]
  if (is.null(entry)) stop("Unknown frozen cohort: ", pair)
  input <- file.path(prepared_dir, entry$design_file)
  if (!identical(meps_decoded_sha(input), entry$decoded_sha256))
    stop("Prepared CSV content differs from the frozen cohort")
  schema <- jsonlite::read_json(file.path(prepared_dir, "common_design_schema.json"))[[entry$schema_key]]
  variables <- unlist(schema$columns, use.names = FALSE)
  stopifnot(identical(variables, unlist(spec$design$columns, use.names = FALSE)))
  connection <- gzfile(input, "rt")
  on.exit(close(connection), add = TRUE)
  records <- utils::read.csv(connection, colClasses = "character", check.names = FALSE)
  metadata <- unlist(spec$metadata_columns, use.names = FALSE)
  stopifnot(identical(names(records), c(metadata, variables)), length(variables) == 35L)
  for (name in setdiff(names(records), c("DUPERSID", "race_group"))) {
    converted <- suppressWarnings(as.numeric(records[[name]]))
    if (anyNA(converted)) stop("Non-numeric frozen field: ", name)
    records[[name]] <- converted
  }
  X <- as.matrix(records[, variables, drop = FALSE]); storage.mode(X) <- "double"
  Z <- as.integer(records$Z); Y <- records$Y; W <- records$W
  stopifnot(nrow(records) == entry$n, sum(Z == 1L) == entry$n_white,
    sum(Z == 0L) == entry$n_comparator, sum(Y == 0) == entry$zero_expenditure_n,
    !anyDuplicated(records$DUPERSID), !anyDuplicated(records$source_row),
    identical(order(records$source_row), seq_len(nrow(records))),
    all(records$Z %in% 0:1), all(Z %in% 0:1), identical(Z == 1L, records$race_group == "White"),
    all(records$race_group[Z == 0L] == entry$comparator),
    all(is.finite(X)), all(is.finite(Y)), all(Y >= 0), all(is.finite(W)), all(W > 0),
    identical(Y, records$TOTEXP09), identical(W, records$SAQWT09F))
  data <- data.frame(A = Z, as.data.frame(X, check.names = FALSE), check.names = FALSE)
  list(data = data, X = X, design = cbind("(Intercept)" = 1, X), variables = variables,
    Y = Y, Z = Z, W = W, n = nrow(X), id = records$DUPERSID,
    source_row = records$source_row, metadata = records[, metadata, drop = FALSE],
    study = "MEPS2009", outcome = "TOTEXP09", units = "2009 USD",
    baseline = pair, ps_weighting = "probability",
    sample_sha256 = digest::digest(list(id = records$DUPERSID, source_row = records$source_row,
      X = X, Y = Y, Z = Z, W = W), algo = "sha256"))
}

meps_sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)

meps_atomic_rds <- function(value, path) {
  if (file.exists(path)) stop("Refuse to overwrite checkpoint: ", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  if (file.exists(temporary)) stop("Inspect existing temporary checkpoint: ", temporary)
  saveRDS(value, temporary, compress = "gzip", version = 3L)
  if (!file.rename(temporary, path)) stop("Could not commit checkpoint: ", path)
  invisible(path)
}

meps_atomic_csv <- function(value, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  if (file.exists(temporary)) stop("Inspect existing temporary CSV: ", temporary)
  utils::write.csv(value, temporary, row.names = FALSE, na = "NA")
  if (file.exists(path)) {
    if (!identical(meps_sha(path), meps_sha(temporary)))
      stop("Existing CSV differs from its saved summary: ", path)
    unlink(temporary)
  } else if (!file.rename(temporary, path)) stop("Could not commit CSV: ", path)
  invisible(path)
}

meps_checkpoint <- function(path, key, action = NULL) {
  if (file.exists(path)) {
    saved <- readRDS(path)
    if (!identical(saved$key, key)) stop("Checkpoint contract differs: ", path)
    return(saved$value)
  }
  if (is.null(action)) stop("Required checkpoint is absent: ", path)
  started <- paste0(path, ".started.rds")
  if (file.exists(started)) stop("Interrupted action requires review before retry: ", started)
  meps_atomic_rds(list(key = key, pid = Sys.getpid(),
    started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)), started)
  value <- action()
  meps_atomic_rds(list(key = key, value = value,
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)), path)
  value
}

meps_match <- function(a, runtime, components, family, estimand) {
  predictions <- runtime$helper$app_predictions(components, family, estimand)
  if (family == "PS" && !components$PS_score$ok) stop(components$PS_score$error)
  if (family == "full_X" && !components$full_X_scores$ok) stop(components$full_X_scores$error)
  scores0 <- switch(family, PS = components$PS_score$value,
    DSM = components$DSM[["0"]]$value$scores, full_X = components$full_X_scores$value)
  scores1 <- if (estimand == "PATE") switch(family, PS = scores0,
    DSM = components$DSM[["1"]]$value$scores, full_X = scores0) else NULL
  do.call(runtime$modules$common$wm_match, c(list(Y = a$Y, Z = a$Z, weights = a$W,
    scores0 = scores0, scores1 = scores1, M = 3L, estimand = estimand,
    method = "self_normalized", variance = FALSE, tie_rule = "row_order"), predictions))
}

meps_psw <- function(a, probability, estimand, m = rep(1, a$n)) {
  if (length(probability) != a$n || any(!is.finite(probability)) ||
      any(probability <= 0 | probability >= 1) || length(m) != a$n ||
      any(!is.finite(m) | m < 0) || sum(m) != a$n) stop("Invalid PSW inputs")
  if (estimand == "PATE") {
    weight1 <- m * a$W * a$Z / probability
    weight0 <- m * a$W * (1 - a$Z) / (1 - probability)
  } else if (estimand == "PATT") {
    weight1 <- m * a$W * a$Z
    weight0 <- m * a$W * (1 - a$Z) * probability / (1 - probability)
  } else stop("Unknown estimand")
  denominator <- c(White = sum(weight1), Comparator = sum(weight0))
  means <- c(White = sum(weight1 * a$Y), Comparator = sum(weight0 * a$Y)) / denominator
  estimate <- unname(means[1L] - means[2L])
  if (any(!is.finite(c(denominator, means, estimate))) || any(denominator <= 0))
    stop("Nonfinite PSW mean or invalid Hájek denominator")
  list(estimate = estimate, means = means, denominator = denominator,
    method = "Supplied-weight PSW Hájek", estimand = estimand)
}

meps_counts <- function(a, helper, seed) {
  old_kind <- RNGkind(); had_seed <- exists(".Random.seed", globalenv(), inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", globalenv(), inherits = FALSE)
  on.exit({do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, globalenv())
    else if (exists(".Random.seed", globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv())}, add = TRUE)
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); set.seed(as.integer(seed))
  before <- get(".Random.seed", globalenv(), inherits = FALSE)
  counts <- stats::rmultinom(200L, size = a$n, prob = rep(1 / a$n, a$n))
  validation <- helper$app_count_columns(a, counts)
  list(matrix = counts, seed = seed, RNGkind = RNGkind(), seed_before = before,
    seed_after = get(".Random.seed", globalenv(), inherits = FALSE), validation = validation,
    matrix_sha256 = digest::digest(counts, algo = "sha256"),
    all_columns_valid = all(vapply(validation, `[[`, logical(1), "ok")))
}

meps_component_ledger <- function(components, column, pair) {
  blocks <- c(list(PS = components$PS, PS_score = components$PS_score),
    stats::setNames(components$PG, c("PG0", "PG1")),
    stats::setNames(components$linear, c("linear0", "linear1")),
    stats::setNames(components$DSM, c("DSM0", "DSM1")))
  do.call(rbind, lapply(names(blocks), function(name) {
    value <- blocks[[name]]
    data.frame(pair = pair, column = column, component = name, ok = isTRUE(value$ok),
      error = if (isTRUE(value$ok)) "" else value$error,
      warnings = paste(value$warnings, collapse = " | "), stringsAsFactors = FALSE)
  }))
}

meps_check_counts <- function(a, runtime, counts) {
  stopifnot(is.list(counts), is.matrix(counts$matrix),
    identical(dim(counts$matrix), c(as.integer(a$n), 200L)),
    identical(digest::digest(counts$matrix, algo = "sha256"), counts$matrix_sha256))
  checked <- runtime$helper$app_count_columns(a, counts$matrix)
  stopifnot(identical(checked, counts$validation))
  checked
}

# Pure saved-component constructor: no fitted model, graph or count is created.
# The returned capture can be checkpointed before count replication starts.
meps_contribution_inference <- function(a, runtime, components, point, family, estimand) {
  runtime$helper$app_capture(function() {
    if (!isTRUE(point$ok)) stop(point$error)
    fit <- point$value; wm <- runtime$modules$common
    stopifnot(identical(fit$estimand, estimand))
    before <- digest::digest(fit, algo = "sha256")
    stack <- runtime$saved_helpers$meps_cc_bind(a, components, fit, family, wm,
      runtime$saved_helpers$wm_saved_full_x_stack)
    inf <- do.call(wm$wm_fitted_inference,
      c(stack$inference_arguments, list(conf.level = .95)))
    stopifnot(identical(inf$fit, fit), identical(inf$estimate, fit$estimate))
    if (!isTRUE(inf$available))
      stop("Complete fitted inference unavailable: ", inf$unavailable_reason)
    wm$.wm_reciprocal_agree(inf$Sigma, stack$nuisance_covariance,
      "Complete saved-stack covariance binding")
    matrices <- list(psi = stack$estimating_equations, J = stack$jacobian,
      IF = stack$nuisance_influence, D0 = stack$inference_arguments$mean_derivative0,
      D1 = stack$inference_arguments$mean_derivative1,
      DW = stack$inference_arguments$weight_derivative)
    input_digests <- lapply(matrices, function(x) if (is.null(x)) NULL else list(
      dimensions = dim(x), sha256 = digest::digest(x, algo = "sha256")))
    stopifnot(identical(digest::digest(fit, algo = "sha256"), before))
    list(inference = inf, diagnostics = stack$diagnostics_saved_only,
      input_digests = input_digests, original_fit_object_sha256 = before)
  })
}

# Same-point complete-contribution operator, distinct from wm_bootstrap_refit.
# Every requested count column is retained; any invalid column leaves inference
# unavailable. The conditional row variance is not the finite-B sample variance.
meps_contribution_replication <- function(a, runtime, point, counts, inference,
                                         family, estimand, original_row) {
  checked <- meps_check_counts(a, runtime, counts)
  stopifnot(nrow(original_row) == 1L,
    identical(as.character(original_row$method), paste0("WM_", family)),
    identical(as.character(original_row$estimand), estimand))
  wm <- runtime$modules$common
  point_ok <- isTRUE(point$ok)
  fit <- if (point_ok) point$value else NULL
  fit_digest <- if (point_ok) digest::digest(fit, algo = "sha256") else NULL
  if (point_ok) wm$.wm_reciprocal_agree(as.numeric(original_row$estimate), fit$estimate,
    "Original summary and saved point")
  had_seed <- exists(".Random.seed", globalenv(), inherits = FALSE)
  if (had_seed) seed_before <- get(".Random.seed", globalenv(), inherits = FALSE)
  failed_columns <- which(!vapply(checked, `[[`, logical(1), "ok"))
  result <- runtime$helper$app_capture(function() {
    if (!point_ok) stop(point$error)
    if (!isTRUE(inference$ok)) stop(inference$error)
    if (length(failed_columns)) stop(structure(list(message = paste(
      "Complete-contribution count columns failed:", paste(failed_columns, collapse = ",")),
      call = NULL, failed_columns = failed_columns),
      class = c("application_column_error", "error", "condition")))
    inf <- inference$value$inference
    stopifnot(identical(inf$fit, fit), identical(inf$estimate, fit$estimate),
      identical(inference$value$original_fit_object_sha256, fit_digest))
    replication <- wm$wm_bootstrap(inf, counts = counts$matrix,
      interval = "normal", conf.level = .95)
    stopifnot(identical(replication$estimate, fit$estimate), replication$B == 200L,
      identical(replication$source_inference, inf), !isTRUE(replication$original_refit_bootstrap),
      length(replication$root_n_draws) == 200L, all(is.finite(replication$root_n_draws)))
    wm$.wm_reciprocal_agree(replication$conditional_variance, inf$variance,
      "Complete row variance")
    list(diagnostics = inference$value$diagnostics,
      input_digests = inference$value$input_digests,
      sensitivity = inf$total_sensitivity, smooth_sensitivity = inf$smooth_sensitivity,
      graph_sensitivity = inf$graph_sensitivity,
      covariance = list(V0 = inf$V0, C = inf$C, Sigma = inf$Sigma,
        cross_term = inf$cross_term, nuisance_variance = inf$nuisance_variance,
        root_n_variance = inf$root_n_variance, formula_error = inf$covariance_formula_error),
      variance = replication$conditional_variance, se = sqrt(replication$conditional_variance),
      conf.int = replication$conf.int, root_n_draws = replication$root_n_draws,
      monte_carlo_variance = replication$monte_carlo_variance,
      monte_carlo_root_n_variance = replication$monte_carlo_root_n_variance,
      common_contract = inf$contract, replication_contract = replication$replication_contract,
      counts_sha256 = counts$matrix_sha256, B = replication$B,
      original_fit_object_sha256 = fit_digest, source_estimate = fit$estimate,
      all_B_retained = TRUE, point_graph_prediction_weight_unchanged = TRUE,
      full_IF_columns_retained = TRUE, numerical_available = TRUE,
      population_assumptions_verified = FALSE,
      original_refit_available = identical(as.character(original_row$status), "empirical_complete"))
  })
  okay <- isTRUE(result$ok)
  result$row <- data.frame(pair = a$baseline, method = paste0("WM_", family),
    estimand = estimand, M = 3L, n = a$n, n_treated = sum(a$Z), B = 200L,
    estimate = if (point_ok) fit$estimate else NA_real_,
    SE = if (okay) result$value$se else NA_real_,
    lower = if (okay) result$value$conf.int[1L] else NA_real_,
    upper = if (okay) result$value$conf.int[2L] else NA_real_,
    status = if (okay) "conditional_contribution_available" else "unavailable",
    error = if (okay) "" else result$error,
    original_refit_status = original_row$status, original_refit_error = original_row$error,
    conditional_variance = if (okay) result$value$variance else NA_real_,
    monte_carlo_variance_B_minus_1 = if (okay) result$value$monte_carlo_variance else NA_real_,
    counts_sha256 = counts$matrix_sha256, units = a$units, target_scope = original_row$target_scope,
    inference_scope = paste("Same-point complete-contribution normal interval; original full-n counts;",
      "correct full-X working predictions and applicable actual-graph/root/row-moment conditions",
      "are substantive unverified premises; no Gaussian expenditure or complex-design claim"),
    scientific_acceptance = FALSE, stringsAsFactors = FALSE)
  result$source_binding <- list(pair = a$baseline, sample_sha256 = a$sample_sha256,
    original_fit_object_sha256 = fit_digest, counts_sha256 = counts$matrix_sha256)
  result$count_validation <- checked
  result$inference_warnings <- inference$warnings
  stopifnot(identical(had_seed, exists(".Random.seed", globalenv(), inherits = FALSE)),
    identical(digest::digest(counts$matrix, algo = "sha256"), counts$matrix_sha256))
  if (point_ok) stopifnot(identical(digest::digest(fit, algo = "sha256"), fit_digest))
  if (had_seed) stopifnot(identical(seed_before, get(".Random.seed", globalenv(), inherits = FALSE)))
  result
}

meps_summarize <- function(a, runtime, output, contract, counts, points) {
  capture <- runtime$helper$app_capture
  empty <- matrix(NA_real_, a$n, 200L)
  refits <- list(means = list(linear = list(mean0 = empty, mean1 = empty),
                             DSM = list(mean0 = empty, mean1 = empty)))
  psw_draws <- matrix(NA_real_, 200L, 2L, dimnames = list(NULL, c("PATE", "PATT")))
  psw_columns <- vector("list", 200L); ledger <- vector("list", 200L)
  for (b in seq_len(200L)) {
    key <- list(contract = contract, column = b,
      count_sha256 = counts$validation[[b]]$column_sha256)
    saved <- meps_checkpoint(file.path(output, "refits", sprintf("column_%03d.rds", b)), key)
    components <- saved$components
    ledger[[b]] <- meps_component_ledger(components, b, a$baseline)
    for (kind in c("linear", "DSM")) for (z in 0:1) {
      value <- components[[kind]][[as.character(z)]]
      if (isTRUE(value$ok)) refits$means[[kind]][[paste0("mean", z)]][, b] <- value$value$mean
    }
    psw_columns[[b]] <- saved$PSW
    for (estimand in c("PATE", "PATT")) if (isTRUE(saved$PSW[[estimand]]$ok))
      psw_draws[b, estimand] <- saved$PSW[[estimand]]$value$estimate
  }
  rows <- results <- list()
  for (family in c("PS", "DSM", "full_X")) for (estimand in c("PATE", "PATT")) {
    name <- paste(family, estimand, sep = "_"); point <- points[[name]]
    result <- meps_checkpoint(file.path(output, "replication", paste0(name, ".rds")),
      list(contract = contract, counts = counts$matrix_sha256, method = name), function() {
        replication <- capture(function() {
          if (!isTRUE(point$ok)) stop(point$error)
          # Existing helper uses the same n-row count vector for the cached
          # predictions and the exported original fixed-reuse WM calculation.
          runtime$helper$app_replicate(point$value, runtime$modules, counts$matrix,
            refits, family, estimand)
        })
        list(replication = replication,
          row = runtime$helper$app_row(a, family, 3L, estimand, point, replication))
      })
    results[[name]] <- result; rows[[name]] <- result$row
  }
  for (estimand in c("PATE", "PATT")) {
    name <- paste0("PSW_", estimand); point <- points[[name]]
    draws <- psw_draws[, estimand]
    failed <- which(!is.finite(draws))
    replication <- capture(function() {
      if (!isTRUE(point$ok)) stop(point$error)
      if (length(failed)) stop(structure(list(message = paste("PSW count refits failed:",
        paste(failed, collapse = ",")), call = NULL, failed_columns = failed),
        class = c("application_column_error", "error", "condition")))
      variance <- mean((draws - mean(draws))^2); se <- sqrt(variance)
      interval <- point$value$estimate + c(-1, 1) * stats::qnorm(.975) * se
      if (any(!is.finite(c(variance, se, interval)))) stop("Nonfinite PSW uncertainty")
      list(draws = draws, variance = variance, se = se, conf.int = interval,
        B = 200L, variance_divisor = "B", predictions = "same-count refitted PS")
    })
    okay <- isTRUE(point$ok) && isTRUE(replication$ok)
    rows[[name]] <- data.frame(study = a$study, outcome = a$outcome, baseline = a$baseline,
      method = "PSW_Hajek", estimand = estimand, M = NA_integer_, n = a$n,
      n_treated = sum(a$Z), B = 200L,
      estimate = if (point$ok) point$value$estimate else NA_real_, raw_estimate = NA_real_,
      SE = if (okay) replication$value$se else NA_real_,
      lower = if (okay) replication$value$conf.int[1L] else NA_real_,
      upper = if (okay) replication$value$conf.int[2L] else NA_real_,
      status = if (!point$ok) "failed_point" else if (!okay) "failed_interval" else "empirical_complete",
      error = if (!point$ok) point$error else if (!okay) replication$error else "",
      units = a$units, variance_divisor = "B", inference_scope = paste(
        "Smooth same-count refitted-PS Hajek interval; supplied W fixed;",
        "requires its own iid/root/overlap conditions; no WM fixed-reuse theorem claim"))
    results[[name]] <- list(point = point, replication = replication,
      draws = draws, failed_columns = failed, refit_results = lapply(psw_columns, `[[`, estimand))
  }
  table <- do.call(rbind, rows); rownames(table) <- NULL
  wm <- startsWith(table$method, "WM_")
  table$inference_scope[wm] <- paste("Original same-count fixed-reuse/refit nominal interval;",
    "supplied W/observed graphs/original K fixed; applicable sampling/refit agreement",
    "and working-model assumptions are not established by this application")
  table$target_scope <- ifelse(table$estimand == "PATE",
    "Combined-pair supplied-weight standardized complete-record disparity",
    "White-reference supplied-weight standardized complete-record disparity")
  table$contrast <- "White minus comparator"; table$causal_intervention_on_race <- FALSE
  table$complex_design_inference <- FALSE; table$exact_author_reproduction <- FALSE
  table$tie_rule <- ifelse(wm, "row_order", NA_character_)
  table$matching_columns <- ifelse(table$method == "WM_PS", 1L,
    ifelse(table$method == "WM_DSM", 2L, ifelse(table$method == "WM_full_X", 35L, NA_integer_)))
  list(table = table, results = results, component_failures = do.call(rbind, ledger),
    count_validation = counts$validation, counts_sha256 = counts$matrix_sha256,
    all_requested_points_available = all(vapply(points, function(x) isTRUE(x$ok), logical(1))),
    all_requested_intervals_available = all(table$status == "empirical_complete"),
    scientific_acceptance = FALSE,
    scope = "Observed execution output only; no finite-sample nominal coverage or complex-design certification")
}

meps_method_names <- function() {
  as.vector(outer(c("PS", "DSM", "full_X"), c("PATE", "PATT"), paste, sep = "_"))
}

meps_run <- function(prepared_dir, output_dir, pair, release_dir,
                     phase = c("preflight", "points", "balance", "counts", "refits",
                               "original-summary", "contribution", "summary"),
                     columns = seq_len(200L), methods = meps_method_names()) {
  phase <- match.arg(phase)
  if (length(pair) != 1L || !pair %in% c("white_asian", "white_hispanic"))
    stop("pair must be white_asian or white_hispanic")
  if (!is.numeric(columns) || !length(columns) || anyNA(columns) ||
      any(columns != as.integer(columns)) || any(!columns %in% seq_len(200L)) ||
      anyDuplicated(columns)) stop("columns must be distinct integers in 1:200")
  if (!is.character(methods) || !length(methods) || anyNA(methods) ||
      any(!methods %in% meps_method_names()) || anyDuplicated(methods))
    stop("methods must be distinct names among PS/DSM/full_X crossed with PATE/PATT")
  prepared_dir <- normalizePath(prepared_dir, mustWork = TRUE)
  runtime <- meps_runtime(release_dir)
  spec <- runtime$spec
  a <- meps_input(prepared_dir, spec, pair)
  output <- file.path(output_dir, pair)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  output <- normalizePath(output, mustWork = TRUE)
  lock <- file.path(output, ".phase.lock")
  if (!dir.create(lock, showWarnings = FALSE))
    stop("Another phase is active, or a stale lock requires process inspection: ", lock)
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  saveRDS(list(pid = Sys.getpid(), phase = phase,
    started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)), file.path(lock, "owner.rds"))
  contract <- digest::digest(list(specification = spec, runtime = runtime$identity,
    cohort = pair, sample_sha256 = a$sample_sha256), algo = "sha256")
  context <- meps_checkpoint(file.path(output, "00_context.rds"), contract, function()
    list(contract_sha256 = contract, runtime = runtime$identity, specification = spec,
      n = a$n, n_white = sum(a$Z), original_row_ids = a$id, source_row = a$source_row,
      sample_sha256 = a$sample_sha256, predictor_columns = a$variables,
      units = a$units, target = spec$target, session = utils::sessionInfo(),
      no_fits_graphs_counts_at_preflight = TRUE))
  if (phase == "preflight") return(invisible(context))
  capture <- runtime$helper$app_capture
  components_path <- file.path(output, "point_components.rds")
  if (phase == "points") {
    components <- meps_checkpoint(components_path, contract,
      function() runtime$helper$app_components(a, runtime$modules))
    for (family in c("PS", "DSM", "full_X")) for (estimand in c("PATE", "PATT")) {
      name <- paste(family, estimand, sep = "_")
      meps_checkpoint(file.path(output, "points", paste0(name, ".rds")),
        list(contract = contract, method = name), function()
          capture(function() meps_match(a, runtime, components, family, estimand)))
    }
    for (estimand in c("PATE", "PATT")) {
      name <- paste0("PSW_", estimand)
      meps_checkpoint(file.path(output, "points", paste0(name, ".rds")),
        list(contract = contract, method = name), function() capture(function() {
          if (!components$PS$ok) stop(components$PS$error)
          meps_psw(a, components$PS$value$probability, estimand)
        }))
    }
    return(invisible(list(status = "POINT_ACTIONS_RECORDED", contract = contract)))
  }
  if (phase == "balance") {
    points <- stats::setNames(lapply(meps_method_names(), function(name)
      meps_checkpoint(file.path(output, "points", paste0(name, ".rds")),
        list(contract = contract, method = name))), meps_method_names())
    schema <- jsonlite::read_json(file.path(prepared_dir, "common_design_schema.json"))[[spec$cohorts[[pair]]$schema_key]]
    diagnostic <- meps_checkpoint(file.path(output, "BALANCE.rds"), contract,
      function() runtime$saved_helpers$meps_saved_balance(a, schema, points, runtime$helper))
    meps_atomic_csv(diagnostic$statuses, file.path(output, "requested_graph_status.csv"))
    if (!is.null(diagnostic$covariate_balance))
      meps_atomic_csv(diagnostic$covariate_balance, file.path(output, "covariate_balance.csv"))
    if (!is.null(diagnostic$balance_summary))
      meps_atomic_csv(diagnostic$balance_summary, file.path(output, "balance_summary.csv"))
    return(invisible(diagnostic))
  }
  counts_path <- file.path(output, "shared_counts.rds")
  if (phase == "counts") return(invisible(meps_checkpoint(counts_path, contract,
    function() meps_counts(a, runtime$helper, spec$cohorts[[pair]]$count_seed))))
  counts <- meps_checkpoint(counts_path, contract)
  meps_check_counts(a, runtime, counts)
  if (phase == "refits") {
    for (b in as.integer(columns)) {
      key <- list(contract = contract, column = b,
        count_sha256 = counts$validation[[b]]$column_sha256)
      meps_checkpoint(file.path(output, "refits", sprintf("column_%03d.rds", b)), key,
        function() {
          m <- counts$matrix[, b]
          components <- runtime$helper$app_components(a, runtime$modules, m)
          PSW <- stats::setNames(lapply(c("PATE", "PATT"), function(estimand)
            capture(function() {
              if (!components$PS$ok) stop(components$PS$error)
              meps_psw(a, components$PS$value$probability, estimand, m)
            })), c("PATE", "PATT"))
          list(components = components, PSW = PSW, count_validation = counts$validation[[b]])
        })
    }
    return(invisible(list(status = "REQUESTED_COUNT_ACTIONS_RECORDED", columns = columns,
      contract = contract)))
  }
  method_names <- c(meps_method_names(), "PSW_PATE", "PSW_PATT")
  points <- stats::setNames(lapply(method_names, function(name)
    meps_checkpoint(file.path(output, "points", paste0(name, ".rds")),
      list(contract = contract, method = name))), method_names)
  original_path <- file.path(output, "ORIGINAL_SUMMARY.rds")
  if (phase == "original-summary") {
    original <- meps_checkpoint(original_path, contract,
      function() meps_summarize(a, runtime, output, contract, counts, points))
    meps_atomic_csv(original$table, file.path(output, "original_table.csv"))
    meps_atomic_csv(original$component_failures, file.path(output, "component_failures.csv"))
    return(invisible(original))
  }
  original <- meps_checkpoint(original_path, contract)
  contribution_key <- function(name) list(contract = contract,
    counts = counts$matrix_sha256, method = name, operator = "complete_contribution")
  if (phase == "contribution") {
    components <- meps_checkpoint(components_path, contract)
    for (name in methods) {
      estimand <- if (endsWith(name, "_PATE")) "PATE" else "PATT"
      family <- sub(paste0("_", estimand, "$"), "", name)
      folder <- file.path(output, "contributions", name)
      key <- contribution_key(name)
      # Completed method records (including failures) skip construction entirely.
      if (file.exists(file.path(folder, "record.rds"))) {
        meps_checkpoint(file.path(folder, "record.rds"), key)
        next
      }
      inference <- meps_checkpoint(file.path(folder, "inference.rds"),
        list(contract = contract, method = name, operator = "complete_inference"),
        function() meps_contribution_inference(a, runtime, components,
          points[[name]], family, estimand))
      oldrow <- original$table[original$table$method == paste0("WM_", family) &
        original$table$estimand == estimand, , drop = FALSE]
      meps_checkpoint(file.path(folder, "record.rds"), key,
        function() meps_contribution_replication(a, runtime, points[[name]], counts,
          inference, family, estimand, oldrow))
      rm(inference); gc(verbose = FALSE)
    }
    return(invisible(list(status = "REQUESTED_CONTRIBUTIONS_RECORDED", methods = methods,
      contract = contract)))
  }
  final <- meps_checkpoint(file.path(output, "FINAL_SUMMARY.rds"), contract, function() {
    contributions <- stats::setNames(lapply(meps_method_names(), function(name)
      meps_checkpoint(file.path(output, "contributions", name, "record.rds"),
        contribution_key(name))), meps_method_names())
    contribution_table <- do.call(rbind, lapply(contributions, `[[`, "row"))
    rownames(contribution_table) <- NULL
    list(original_table = original$table, contribution_table = contribution_table,
      original_results = original$results, contributions = contributions,
      component_failures = original$component_failures,
      count_validation = counts$validation, counts_sha256 = counts$matrix_sha256,
      inference_operators_combined_or_relabelled = FALSE, scientific_acceptance = FALSE)
  })
  meps_atomic_csv(final$original_table, file.path(output, "original_table.csv"))
  meps_atomic_csv(final$contribution_table, file.path(output, "complete_contribution_table.csv"))
  meps_atomic_csv(final$component_failures, file.path(output, "component_failures.csv"))
  invisible(final)
}
