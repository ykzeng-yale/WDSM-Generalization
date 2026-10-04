# Archived saved-array calculation source; execution status and scope are documented in results/README.md.
# Read authenticated saved arrays; never select donors or fit a statistical root.
# lns_bind() performs all source/data/point checks before lns_augment() evaluates
# any J, psi, IF or complete-contribution variance. No bootstrap is called.

.lns_require <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}
.lns_sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
.lns_hash <- function(x) digest::digest(x, algo = "sha256")
.lns_fields <- function(x, fields, label, exact = FALSE) {
  .lns_require(is.list(x) && !is.null(names(x)) && !anyDuplicated(names(x)) &&
    all(fields %in% names(x)) && (!exact || setequal(names(x), fields)),
    paste(label, "has missing, duplicate or unexpected fields"))
}
.lns_finite <- function(x, n, label) {
  .lns_require(is.numeric(x) && !is.complex(x) && length(x) == n &&
    all(is.finite(x)), paste(label, "has invalid numeric values/length"))
}
.lns_agree <- function(x, y, label, tolerance = 1e-10) {
  .lns_require(is.numeric(x) && is.numeric(y) && length(x) == length(y) &&
    all(is.finite(x)) && all(is.finite(y)) &&
    all(abs(x - y) <= tolerance * pmax(1, abs(x), abs(y))),
    paste(label, "does not reproduce the authenticated value"))
}
.lns_capture <- function(fun) {
  warnings <- character()
  value <- tryCatch(withCallingHandlers(fun(), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = identity)
  list(ok = !inherits(value, "error"), value = if (!inherits(value, "error")) value else NULL,
    error = if (inherits(value, "error")) conditionMessage(value) else "", warnings = warnings)
}
.lns_rng <- function() list(kind = RNGkind(), present = exists(".Random.seed", .GlobalEnv,
  inherits = FALSE), seed = get0(".Random.seed", .GlobalEnv, inherits = FALSE))

# contract is the reviewed JSON object, rooted in the existing workspace.
# Hash only consumed compact inputs, small provenance files and one population
# cache. Never read/hash the large compressed shard archive.
lns_read <- function(contract, root) {
  .lns_require(identical(contract$response, "No") && contract$M == 3L &&
    identical(contract$estimand, "PATT"), "Unsupported comparison contract")
  rng <- .lns_rng()
  pins <- contract$input_pins
  .lns_require(is.list(pins) && length(pins) > 0L && !is.null(names(pins)), "Input pins missing")
  paths <- vapply(names(pins), function(k) normalizePath(file.path(root, k), mustWork = TRUE), "")
  for (k in names(pins)) .lns_require(identical(.lns_sha(paths[[k]]), pins[[k]]), paste("Hash changed:", k))
  for (k in names(contract$calculation_pins)) .lns_require(
    identical(.lns_sha(file.path(root, k)), contract$calculation_pins[[k]]),
    paste("Calculation source changed:", k))
  path <- function(key) paths[[contract$files[[key]]]]
  json <- function(key) jsonlite::read_json(path(key), simplifyVector = FALSE)
  collection <- json("collection_receipt"); verification <- json("archive_verification")
  inventory <- json("archive_inventory"); manifest <- json("shard_manifest")
  receipt <- json("shard_receipt")
  for (key in c("archive_verification", "archive_inventory", "shard_manifest", "shard_receipt")) {
    .lns_require(identical(collection$consumed_input_pins[[path(key)]], .lns_sha(path(key))),
      paste("Accepted collection does not bind", key))
  }
  .lns_require(identical(verification$inventory_sha256, .lns_sha(path("archive_inventory"))) &&
    identical(verification$archive_sha256, verification$remote$sha256) &&
    identical(verification$remote$manifest_sha256, collection$manifest_sha256),
    "Original archive receipt chain is inconsistent")
  .lns_require(receipt$status %in% c("SHARD_COMPLETE", "SHARD_COMPLETE_WITH_FAILURES_RETAINED") &&
    isTRUE(receipt$source_pins_unchanged) && contract$population %in% unlist(manifest$job$population_scenarios) &&
    contract$replicate %in% unlist(manifest$job$replicates), "Original shard completion/job mismatch")
  inv_paths <- vapply(inventory$files, `[[`, "", "path")
  .lns_require(!anyDuplicated(inv_paths), "Archive inventory contains duplicate paths")
  for (key in c("selected", "response", "compact", "outcomes", "records")) {
    member <- contract$archive_members[[key]]
    at <- match(member, inv_paths)
    .lns_require(!is.na(at), paste("Missing archive member", member))
    entry <- inventory$files[[at]]
    .lns_require(identical(entry$sha256, .lns_sha(path(key))) &&
      entry$size == file.info(path(key))$size, paste("Original archive member changed:", key))
  }
  cache_name <- paste0("scenario_", contract$population, ".rds")
  .lns_require(identical(manifest$cache_receipt$cache_hashes[[cache_name]], pins[[contract$files$cache]]),
    "The population cache differs from the executed source cache")
  for (key in names(contract$frozen_sources)) {
    rel <- contract$frozen_sources[[key]]$historical_path
    .lns_require(identical(manifest$source_pins[[rel]],
      pins[[contract$frozen_sources[[key]]$path]]), paste("Historical source identity mismatch:", key))
  }
  # Evaluate only this pure function definition from the authenticated writer.
  expr <- parse(path("batch_source"))
  take <- vapply(expr, function(e) is.call(e) && identical(e[[1L]], as.name("<-")) &&
    identical(e[[2L]], as.name("lenis_selected_from_cache")), logical(1))
  .lns_require(sum(take) == 1L, "Pure saved-selection helper definition is missing or duplicated")
  holder <- new.env(parent = baseenv()); eval(expr[[which(take)]], holder)
  selected <- readRDS(path("selected")); cache <- readRDS(path("cache"))
  input <- holder$lenis_selected_from_cache(cache, selected)
  rm(cache)
  schema <- c(population_scenario = "integer", multiplier = "integer", response = "character",
    replicate = "integer", scenario = "character", estimand = "character", selected_n = "integer",
    respondent_n = "integer", method = "character", target = "numeric", estimate = "numeric",
    raw_estimate = "numeric", variance = "numeric", lower = "numeric", upper = "numeric",
    source_coverage = "numeric", status = "character", point_status = "character",
    interval_status = "character", warning_count = "integer", error = "character", inference_scope = "character")
  header <- names(read.csv(path("records"), nrows = 0L, check.names = FALSE, colClasses = "character"))
  .lns_require(identical(header, names(schema)), "Original record schema changed")
  records <- read.csv(path("records"), check.names = FALSE, colClasses = unname(schema),
    stringsAsFactors = FALSE, na.strings = "NA", row.names = NULL, fill = FALSE)
  result <- list(input = input, selected = selected, response = readRDS(path("response")),
    compact = readRDS(path("compact")), outcomes = readRDS(path("outcomes")), records = records,
    contract = contract, input_pins = pins, original_manifest = manifest,
    source_binding = "Accepted original collection, verified inventory membership and file hashes")
  .lns_require(identical(rng, .lns_rng()), "Saved-only reading changed RNG state")
  result
}

# This shell reproduces wm_match's deterministic PATT decomposition from the
# recorded edges. Selected-edge distance/order validation is allowed; no donor
# search, candidate-distance matrix, ranking or nearest-neighbor constructor runs.
.lns_fit <- function(Y, Z, W, scores, mean0, graph, estimate, raw_estimate, wm) {
  n <- length(Y); w <- W / max(W); e <- graph$edges
  queries <- which(Z == 1L)
  .lns_require(is.data.frame(e) && identical(names(e),
    c("query", "donor", "arm", "distance", "share", "outcome_share")) &&
    nrow(e) == 3L * length(queries) && all(e$query == rep(queries, each = 3L)) &&
    all(e$donor %in% which(Z == 0L)) && all(e$arm == 0L), "Saved M3 PATT edge structure is invalid")
  neighbors <- vector("list", n)
  for (position in seq_along(queries)) {
    block <- (as.double(position) - 1) * 3L + seq_len(3L)
    neighbors[[queries[position]]] <- e$donor[block]
  }
  sum_at <- wm$.wm_sum_at
  .lns_require(is.function(sum_at), "Common summation helper missing")
  K <- sum_at(e$donor, w[e$query] * e$share, n)
  .lns_agree(as.numeric(graph$analysis_weights), w, "Original common weight scaling")
  .lns_require(is.matrix(graph$loads) && identical(dim(graph$loads), c(n, 2L)), "Incoming-load shape mismatch")
  .lns_agree(as.numeric(graph$loads), as.numeric(cbind(K, 0)), "Original incoming loads")
  imputed <- raw <- matrix(NA_real_, n, 2L, dimnames = list(NULL, c("Y0", "Y1")))
  imputed[cbind(seq_len(n), Z + 1L)] <- Y; raw[cbind(seq_len(n), Z + 1L)] <- Y
  missing <- sum_at(e$query, e$outcome_share * (Y[e$donor] - mean0[e$donor]), n)
  raw_missing <- sum_at(e$query, e$outcome_share * Y[e$donor], n)
  imputed[queries, 1L] <- mean0[queries] + missing[queries]
  raw[queries, 1L] <- raw_missing[queries]
  outer <- Z * w; denominator <- sum(outer); gamma <- denominator / n
  numerator <- sum(outer * (Y - imputed[, 1L]))
  direct_raw <- sum(outer * (Y - raw[, 1L])) / denominator
  coefficient <- Z * w - (1L - Z) * K
  uncentered <- coefficient * (Y - mean0)
  .lns_agree(numerator / denominator, estimate, "Saved corrected point (query form)")
  .lns_agree(sum(uncentered) / denominator, estimate, "Saved corrected point (load form)")
  .lns_agree(direct_raw, raw_estimate, "Saved raw point (query form)")
  .lns_agree(sum(coefficient * Y) / denominator, raw_estimate, "Saved raw point (load form)")
  actual_numerator <- uncentered - outer * estimate; row <- actual_numerator / gamma
  fit <- structure(list(estimate = estimate, raw_estimate = raw_estimate,
    correction = estimate - raw_estimate, n = n, M = 3L, estimand = "PATT", method = "self_normalized",
    root_n_variance = NA_real_, variance = NA_real_, se = NA_real_, numerator = numerator,
    denominator = denominator, gamma = gamma, weights = W, analysis_weights = graph$analysis_weights,
    weight_scale = max(W), data = list(Y = Y, Z = as.integer(Z)), imputed = imputed, raw_imputed = raw,
    predictions = list(mean0 = mean0, mean1 = NULL, rho0 = NULL, rho1 = NULL),
    loads = list(incoming = graph$loads, mean_incoming = numeric(n)),
    contributions = list(row = row, edge = numeric(), actual = row, nuisance = numeric(n)),
    numerator_identity = list(direct = numerator, reconstructed = sum(uncentered),
      difference = numerator - sum(uncentered), centered_sum = sum(actual_numerator),
      row_edge_sum = sum(actual_numerator)),
    graph = list(edges = e, neighbors = neighbors,
      unordered = data.frame(i = integer(), j = integer(), i_to_j = logical(), j_to_i = logical()),
      cell = rep(1L, n), fold_id = rep(1L, n), strata = rep(1L, n), scores0 = scores, scores1 = NULL,
      transform = "Supplied coordinates; no automatic scaling", tie_rule = "Exact distance, then original row index"),
    info = list(inference_status = "point_estimate_only", corrected = TRUE, nuisance_correction = FALSE,
      contract = "Rehydrated authenticated original point; no population assumptions verified")),
    class = c("wm_match", "list"))
  wm$.wm_reciprocal_validated_state(fit)
  fit
}

# Binding stage: every original identity is checked before any nuisance matrix
# or inference calculation. Retain all six saved points, not only multiplier 1.
lns_bind <- function(saved, wm) {
  rng <- .lns_rng(); d <- saved$input$selected; Y <- saved$input$Y
  response <- saved$response; compact <- saved$compact; n <- nrow(d)
  .lns_require(is.data.frame(d) && all(c("id", "Prob", "s.wt", "r1", "y", "z",
    paste0("x", 1:6)) %in% names(d)), "Original selected-data columns missing")
  .lns_require(is.environment(wm) && is.function(wm$.wm_reciprocal_validated_state), "Common validator missing")
  .lns_fields(response, c("ok", "selected_ids", "respondent_rows", "respondent_ids", "supplied_weights",
    "response_probability", "targets", "outcomes_sha256"), "Saved response")
  .lns_require(n == 5000L && isTRUE(response$ok) && is.null(response$response_fit) &&
    identical(response$selected_ids, saved$selected$ids) && identical(response$respondent_ids, d$id) &&
    identical(response$respondent_rows, seq_len(n)) && !anyDuplicated(d$id) &&
    all(response$response_probability == 1) && length(response$response_probability) == n &&
    identical(as.numeric(response$supplied_weights), as.numeric(1 / saved$selected$Prob)) &&
    identical(as.numeric(d$s.wt), as.numeric(response$supplied_weights)) &&
    identical(as.numeric(d$Prob), as.numeric(saved$selected$Prob)) && all(d$r1 == 1),
    "No-response IDs, inclusion probabilities or original supplied weights disagree")
  .lns_require(is.matrix(Y) && identical(dim(Y), c(n, 6L)) && all(is.finite(Y)) &&
    identical(.lns_hash(Y), response$outcomes_sha256) &&
    all(Y[d$z == 0, , drop = FALSE] == d$y[d$z == 0]) &&
    identical(as.numeric(response$targets), as.numeric(saved$input$targets)),
    "Original six-outcome matrix, controls or finite-population targets disagree")
  .lns_finite(response$targets, 6L, "Original population targets")
  .lns_fields(compact, c("correction", "nuisance_parameter", "nuisance_solver", "full_x_coefficients",
    "graphs", "scores", "base_predictions"), "Compact WM")
  .lns_require(identical(unname(compact$correction), c("weighted linear full-X",
    "source weighted quadratic double-score", "weighted linear full-X")) &&
    identical(names(compact$correction), c("PS", "DSM", "X6")), "Original correction recipe changed")
  theta_names <- c(paste0("ps:", c("intercept", paste0("x", 1:6))),
    paste0("pg0:", c("intercept", paste0("x", 1:6))), "center:ps", "center:pg0",
    "variance:ps", "variance:pg0", paste0("bc0:", c("intercept", "ps", "pg", "ps_sq", "ps_pg", "pg_sq")))
  theta <- compact$nuisance_parameter
  .lns_finite(theta, 24L, "Saved DSM root")
  .lns_require(identical(names(theta), theta_names), "Saved complete 24-parameter order changed")
  D <- cbind(intercept = 1, as.matrix(d[paste0("x", 1:6)])); X <- D[, -1L, drop = FALSE]
  beta <- compact$full_x_coefficients
  .lns_finite(beta, 7L, "Saved full-X linear coefficients")
  .lns_require(identical(names(beta), colnames(D)), "Saved linear coefficient order changed")
  ps <- stats::plogis(as.vector(D %*% theta[1:7])); pg <- as.vector(D %*% theta[8:14])
  .lns_require(all(ps > 0 & ps < 1) && all(theta[17:18] > 0), "Saved root has nonregular finite scores")
  raw_scores <- cbind(ps, pg)
  .lns_agree(as.numeric(theta[15:16]), colMeans(raw_scores), "Saved pooled centers")
  centered <- sweep(raw_scores, 2L, theta[15:16], "-")
  .lns_agree(as.numeric(theta[17:18]), colMeans(centered^2), "Saved pooled variances")
  S <- sweep(centered, 2L, sqrt(theta[17:18]), "/")
  Xc <- sweep(X, 2L, colMeans(X), "-"); Xv <- colMeans(Xc^2)
  .lns_require(all(Xv > 0 & is.finite(Xv)), "Original X6 coordinate variance is not positive")
  XS <- sweep(Xc, 2L, sqrt(Xv), "/")
  Q <- cbind(1, S[, 1L], S[, 2L], S[, 1L]^2, S[, 1L] * S[, 2L], S[, 2L]^2)
  .lns_agree(as.vector(D %*% beta), compact$base_predictions$full_x, "Saved raw linear control predictions")
  .lns_agree(as.vector(Q %*% theta[19:24]), compact$base_predictions$DSM, "Saved quadratic DSM predictions")
  .lns_require(is.list(compact$scores) && length(compact$scores) == 3L && is.null(names(compact$scores)),
    "Expected original writer's unnamed PS/DSM/X6 score list")
  records <- saved$records
  .lns_require(nrow(records) == 150L && !anyDuplicated(records[c("method", "multiplier")]) &&
    all(records$population_scenario == saved$contract$population) &&
    all(records$replicate == saved$contract$replicate) && all(records$response == "No") &&
    all(records$estimand == "PATT") && all(records$scenario == "source_specification") &&
    all(records$selected_n == n & records$respondent_n == n), "Original case record identities disagree")
  fits <- list(); reference <- list(); families <- c("PS", "DSM", "X6")
  reconstructed_scores <- list(S[, 1L, drop = FALSE], S, XS)
  for (k in seq_along(families)) {
    family <- families[k]; key <- paste0(family, "_M3"); graph <- compact$graphs[[key]]
    .lns_fields(graph, c("M", "estimate", "raw_estimate", "edges", "scores_sha256", "mean0_sha256",
      "analysis_weights", "loads"), key, exact = TRUE)
    score <- compact$scores[[k]]; mu <- compact$base_predictions[[if (family == "DSM") "DSM" else "full_x"]]
    .lns_require(graph$M == 3L && is.matrix(score) && identical(dim(score), c(n, c(1L, 2L, 6L)[k])) &&
      identical(.lns_hash(score), graph$scores_sha256) && identical(.lns_hash(mu), graph$mean0_sha256),
      paste(key, "original score/prediction hashes or dimensions disagree"))
    .lns_agree(as.numeric(score), as.numeric(reconstructed_scores[[k]]), paste(key, "actual score map"))
    batch <- saved$outcomes[[key]]
    .lns_fields(batch, c("estimate", "raw_estimate", "draws", "variance", "se", "lower", "upper",
      "variance_divisor", "B", "n", "status", "failed_columns", "completed_columns"), key)
    .lns_finite(batch$estimate, 6L, paste(key, "six points"))
    .lns_finite(batch$raw_estimate, 6L, paste(key, "six raw points"))
    for (field in c("variance", "se", "lower", "upper")) .lns_require(
      is.numeric(batch[[field]]) && length(batch[[field]]) == 6L && is.null(dim(batch[[field]])) &&
        all(is.na(batch[[field]]) | is.finite(batch[[field]])), paste("Saved", field, "shape/range changed"))
    .lns_require(batch$B == 200L && batch$n == n && identical(batch$variance_divisor, "B") &&
      is.matrix(batch$draws) && identical(dim(batch$draws), c(200L, 6L)) &&
      identical(colnames(batch$draws), as.character(1:6)), "Saved refit metadata/draw dimensions changed")
    failed <- batch$failed_columns; good <- setdiff(seq_len(200L), failed)
    .lns_require(is.numeric(failed) && all(failed %in% seq_len(200L)) && !anyDuplicated(failed) &&
      identical(as.integer(batch$completed_columns), as.integer(good)) &&
      all(is.finite(batch$draws[good, , drop = FALSE])) &&
      (!length(failed) || all(is.na(batch$draws[failed, , drop = FALSE]))), "Saved refit failure slots changed")
    if (length(failed)) {
      .lns_require(identical(batch$status, "failed_interval") &&
        all(is.na(c(batch$variance, batch$se, batch$lower, batch$upper))), "Failed refit advertised finite summary")
    } else {
      .lns_require(identical(batch$status, "ok"), "Refit success status disagrees with all saved draws")
      v <- colMeans(sweep(batch$draws, 2L, colMeans(batch$draws), "-")^2)
      .lns_agree(batch$variance, v, "Original refit divisor-B variance")
      .lns_agree(batch$se, sqrt(v), "Saved refit SE")
      .lns_agree(batch$lower, batch$estimate - stats::qnorm(.975) * sqrt(v), "Saved refit lower interval")
      .lns_agree(batch$upper, batch$estimate + stats::qnorm(.975) * sqrt(v), "Saved refit upper interval")
    }
    .lns_agree(graph$estimate, batch$estimate[1L], "Original multiplier-1 corrected point")
    .lns_agree(graph$raw_estimate, batch$raw_estimate[1L], "Original multiplier-1 raw point")
    rows <- records[records$method == paste0("WM_", key), , drop = FALSE]
    .lns_require(nrow(rows) == 6L && setequal(rows$multiplier, 1:6), "Missing original method–multiplier rows")
    rows <- rows[match(1:6, rows$multiplier), , drop = FALSE]
    for (field in c("estimate", "raw_estimate", "variance", "lower", "upper")) {
      .lns_require(identical(is.na(rows[[field]]), is.na(as.numeric(batch[[field]]))), paste(field, "missingness changed"))
      keep <- !is.na(rows[[field]])
      .lns_agree(rows[[field]][keep], as.numeric(batch[[field]])[keep], paste("CSV", field), 32 * .Machine$double.eps)
    }
    .lns_agree(rows$target, as.numeric(response$targets), "CSV finite-population target", 32 * .Machine$double.eps)
    .lns_require(all(rows$point_status == "ok") && all(rows$interval_status == batch$status) &&
      all(rows$status == if (batch$status == "ok") "ok" else "failed"), "Recorded point/interval status changed")
    fits[[family]] <- lapply(1:6, function(q) .lns_fit(as.numeric(Y[, q]), as.integer(d$z),
      as.numeric(d$s.wt), score, mu, graph, batch$estimate[q], batch$raw_estimate[q], wm))
    reference[[family]] <- rows
  }
  result <- list(fits = fits, records = reference, data = d, Y = Y, compact = compact,
    ps = ps, design = D, raw_X = X, input_pins = saved$input_pins, contract = saved$contract,
    saved_outcomes = saved$outcomes, parameter_names = theta_names,
    binding_stage = "all_original_arrays_and_18_points_checked_before_nuisance_evaluation")
  .lns_require(identical(rng, .lns_rng()), "Binding changed RNG state")
  result$mutation_binding <- .lns_hash(result)
  class(result) <- c("lenis_no_saved_binding", "list")
  result
}

# Pure J/psi/derivative evaluation at the saved 24-vector. This function does not
# call .wm_wdsm_nuisance_stack (the fitting constructor).
.lns_dsm_stack <- function(bound, wm) {
  d <- bound$data; n <- nrow(d); D <- bound$design; theta <- bound$compact$nuisance_parameter
  k <- list(ps = 1:7, pg0 = 8:14, center = c(ps = 15L, pg0 = 16L),
    variance = c(ps = 17L, pg0 = 18L), bc0 = 19:24)
  w <- d$s.wt / max(d$s.wt)
  shell <- structure(list(parameter = theta, parameter_names = names(theta), blocks = k,
    used_arms = 0L, estimand = "PATT", inputs = list(Y = d$y, Z = d$z, weights = d$s.wt,
      ps_design = D, pg0_design = D, ps_offset = numeric(n), pg0_offset = numeric(n)),
    multiplicity = rep(1, n), empirical_probability = rep(1 / n, n),
    weight_scale = max(d$s.wt), scaled_weights = w, ps_weighting = "probability", scaled_ps_weights = w),
    class = c(".wm_wdsm_nuisance_stack", "list"))
  value <- wm$.wm_wdsm_nuisance_evaluate(shell)
  .lns_agree(as.numeric(value$scores0), as.numeric(bound$fits$DSM[[1L]]$graph$scores0), "Evaluated saved DSM score")
  .lns_agree(value$mean0, bound$fits$DSM[[1L]]$predictions$mean0, "Evaluated saved DSM prediction")
  influence <- -t(solve(value$jacobian, t(value$estimating_equations)))
  colnames(influence) <- names(theta)
  .lns_require(all(is.finite(influence)), "Complete DSM influence calculation is nonfinite")
  centered <- sweep(influence, 2L, colMeans(influence), "-")
  list(parameter_names = names(theta), parameter = theta, jacobian = value$jacobian,
    estimating_equations = value$estimating_equations, nuisance_influence = influence,
    nuisance_covariance = crossprod(centered) / n, mean0_derivative = value$mean0_derivative,
    inference_arguments = list(nuisance_influence = influence, mean_derivative0 = value$mean0_derivative,
      covariance_scope = "patt", transport0 = list(mode = "zero", basis = "current_centering")))
}

# Exactly three family assemblies are shared across the six treated-outcome
# multipliers. Each multiplier still uses its own observed Y, point and target.
lns_augment <- function(bound, wm, saved_stack) {
  .lns_require(inherits(bound, "lenis_no_saved_binding"), "Run strict binding before augmentation")
  plain <- unclass(bound); hash <- plain$mutation_binding; plain$mutation_binding <- NULL
  .lns_require(identical(.lns_hash(plain), hash), "Bound original arrays were mutated")
  .lns_require(is.function(saved_stack) && is.function(wm$wm_fitted_inference), "Common helper/API missing")
  rng <- .lns_rng(); d <- bound$data; D <- bound$design; X <- bound$raw_X
  source_data <- list(Y = as.numeric(d$y), Z = as.numeric(d$z), weights = as.numeric(d$s.wt))
  coefficients <- bound$compact$full_x_coefficients
  scope <- paste("Saved original Lenis No-response finite-population dependent-design comparison.",
    "Known W supplied once; correct full-X reference mean/current centering is a working inference premise.",
    "No iid nominal-calibration, cluster/stratum, S12 raw-PS, transport or nested-law claim.",
    "No refit, matching or bootstrap count/draw generation; no automatic variance agreement.")
  stacks <- assemblies <- rows <- list()
  for (family in c("PS", "DSM", "X6")) {
    fit <- bound$fits[[family]][[1L]]
    assembled <- .lns_capture(function() {
      if (family == "DSM") return(.lns_dsm_stack(bound, wm))
      if (family == "PS") return(saved_stack(source_data, fit,
        list(mean0 = list(design = D, coefficients = coefficients)), family = "PS",
        ps = list(design = D, probability = bound$ps, weighting = "probability",
          coefficients = stats::setNames(as.numeric(bound$compact$nuisance_parameter[1:7]), colnames(D))), wm = wm))
      # Original X6 correction is LINEAR in raw X. Express the same coefficients
      # in its pooled standardized degree-one basis, without fitting or deleting terms.
      center <- colMeans(X); scale <- sqrt(colMeans(sweep(X, 2L, center, "-")^2))
      standard <- fit$graph$scores0; colnames(standard) <- colnames(X)
      design <- cbind(intercept = 1, standard)
      beta_std <- c(coefficients[1L] + sum(center * coefficients[-1L]), scale * coefficients[-1L])
      names(beta_std) <- colnames(design)
      exponents <- rbind(rep(0, ncol(X)), diag(ncol(X)))
      dimnames(exponents) <- list(colnames(design), colnames(X))
      saved_stack(source_data, fit, list(mean0 = list(design = design, coefficients = beta_std)),
        family = "full_X", raw_X = X, polynomial_exponents = exponents, wm = wm)
    })
    stacks[[family]] <- assembled
    expected_p <- c(PS = 16L, DSM = 24L, X6 = 19L)[[family]]
    if (assembled$ok && length(assembled$value$parameter_names) != expected_p) {
      assembled <- list(ok = FALSE, value = NULL, error = "Complete nuisance dimension mismatch", warnings = character())
      stacks[[family]] <- assembled
    }
    for (q in 1:6) {
      key <- paste0("WM_", family, "_analytic_M3_q", q); fit <- bound$fits[[family]][[q]]
      result <- .lns_capture(function() {
        if (!assembled$ok) stop(assembled$error, call. = FALSE)
        args <- assembled$value$inference_arguments; args$fit <- fit; args$conf.level <- .95
        # Known individual weights: no fitted-weight column or weight update.
        args$weight_derivative <- matrix(0, fit$n, expected_p,
          dimnames = list(NULL, assembled$value$parameter_names))
        inf <- do.call(wm$wm_fitted_inference, args)
        .lns_require(identical(inf$fit, fit) && identical(inf$estimate, fit$estimate), "Common inference changed original point/graph")
        inf
      })
      assemblies[[key]] <- result
      inf <- result$value; ref <- bound$records[[family]][q, , drop = FALSE]
      okay <- result$ok && isTRUE(inf$available) && length(inf$variance) == 1L && is.finite(inf$variance) && inf$variance > 0 &&
        length(inf$conf.int) == 2L && all(is.finite(inf$conf.int))
      reason <- if (!result$ok) result$error else if (!okay) paste("Analytic inference unavailable", inf$unavailable_reason) else ""
      rows[[key]] <- data.frame(population_scenario = ref$population_scenario, replicate = ref$replicate,
        multiplier = q, response = "No", estimand = "PATT", family = family, M = 3L, n = fit$n,
        parameter_count = expected_p, target = ref$target, estimate = fit$estimate, raw_estimate = fit$raw_estimate,
        point_status = "ok", analytic_status = if (okay) "ok" else "failed", error = reason,
        analytic_variance = if (okay) inf$variance else NA_real_,
        analytic_lower = if (okay) inf$conf.int[1L] else NA_real_, analytic_upper = if (okay) inf$conf.int[2L] else NA_real_,
        saved_refit_variance = ref$variance, saved_refit_lower = ref$lower, saved_refit_upper = ref$upper,
        saved_refit_status = ref$interval_status, saved_refit_error = ref$error,
        saved_refit_B = 200L, saved_refit_divisor = "B", assumptions_verified = FALSE,
        source_method = ref$method, scope = scope, stringsAsFactors = FALSE)
    }
  }
  .lns_require(identical(rng, .lns_rng()), "Saved analytic augmentation changed RNG state")
  list(records = do.call(rbind, rows), stacks = stacks, assemblies = assemblies,
    input_pins = bound$input_pins, scope = scope, original_results_modified = FALSE)
}
