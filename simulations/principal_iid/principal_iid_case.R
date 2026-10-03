# A single synthetic sample, using the existing common WM implementation.
# Definitions only; call through reproduce.R for fresh-output recording.
.principal_pack_count <- function(x) {
  if (isTRUE(x$ok)) {
    x$value$mean0 <- NULL
    x$value$mean1 <- NULL
  }
  x
}

principal_iid_case <- function(model, targets, overlap, n, replicate_id, out,
                               wm, original_helpers, stack_helpers, model_env,
                               progress = function(phase) invisible(NULL)) {
  if (length(n) != 1L || !is.numeric(n) || !is.finite(n) || !n %in% c(1000L, 5000L) ||
      length(overlap) != 1L || !is.character(overlap) || is.na(overlap) ||
      !overlap %in% c("GoodOverlap", "PoorOverlap") ||
      length(replicate_id) != 1L || !is.numeric(replicate_id) ||
      !is.finite(replicate_id) || replicate_id < 1 || replicate_id > 1000 ||
      replicate_id != floor(replicate_id)) stop("Use one declared overlap, n=1000/5000 and replicate=1:1000.")
  n <- as.integer(n); replicate_id <- as.integer(replicate_id)
  group_index <- if (overlap == "GoodOverlap") match(n, c(1000L, 5000L)) else
    2L + match(n, c(1000L, 5000L))
  case_id <- sprintf("g%d_%s_n%d_rep%04d", group_index, overlap, n, replicate_id)
  stream_key <- group_index * 100000L + replicate_id
  data_seed <- 104000000L + stream_key; count_seed <- 204000000L + stream_key
  target_rows <- targets[targets$overlap == overlap, , drop = FALSE]
  if (nrow(target_rows) != 2L || anyDuplicated(target_rows$estimand) ||
      !setequal(target_rows$estimand, c("PATE", "PATT")) ||
      !all(target_rows$admitted) || any(!is.finite(target_rows$target)))
    stop("Both targets must be admitted by the unchanged truth calculation.")
  sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
  objsha <- function(x) digest::digest(x, algo = "sha256")
  write_json <- function(x, path) jsonlite::write_json(x, path, auto_unbox = TRUE,
    pretty = TRUE, null = "null", na = "null", digits = NA)
  fatal <- function(message) stop(structure(list(message = message, call = NULL),
    class = c("principal_integrity_error", "error", "condition")))
  guard <- function(label) progress(label)
  capture <- function(fun) {
    warnings <- character(); started <- proc.time()[[3L]]
    value <- tryCatch(withCallingHandlers(fun(), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) {
      if (inherits(e, "principal_integrity_error")) stop(e)
      e
    })
    list(ok = !inherits(value, "error"), value = if (!inherits(value, "error")) value else NULL,
         error = if (inherits(value, "error")) conditionMessage(value) else "",
         attempts = if (inherits(value, "error")) value$attempts else NULL,
         warnings = warnings, elapsed_seconds = proc.time()[[3L]] - started)
  }
  require_value <- function(x) { if (!isTRUE(x$ok)) stop(x$error, call. = FALSE); x$value }
  check_close <- function(x, y, label) {
    if (length(x) != length(y) || any(!is.finite(c(x, y))) ||
        max(abs(x - y)) > 1e-10 * max(1, abs(x), abs(y))) fatal(paste(label, "disagrees"))
  }
  save_new <- function(x, relative) {
    path <- file.path(out, relative)
    if (file.exists(path)) fatal(paste("Refusing to overwrite", relative))
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(x, path, compress = "gzip", version = 3)
    list(path = relative, sha256 = sha(path), bytes = unname(file.info(path)$size))
  }
  records <- list(); files <- list()
  sampled <- model_env$principal_iid_sample(model, overlap, n, seed = data_seed)
  dat <- sampled$data
  if (nrow(dat) != n || !all(c(paste0("X", 1:6), "Z", "W", "Y", "mu0", "mu1", "pi", "e") %in% names(dat)) ||
      any(!is.finite(as.matrix(dat[c(paste0("X", 1:6), "Z", "W", "Y")])))) fatal("Generated sample schema invalid.")
  if (!identical(as.numeric(dat$W), as.numeric(1/dat$pi))) fatal("Supplied W is not the known reciprocal selection probability.")
  files$data <- save_new(sampled, "sample.rds")
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(count_seed)
  C <- stats::rmultinom(200L, n, rep(1/n, n))
  rownames(C) <- as.character(seq_len(n)); colnames(C) <- sprintf("count%03d", seq_len(200L))
  rng_after_counts <- list(kind = RNGkind(), seed = get(".Random.seed", envir = .GlobalEnv))
  if (!identical(dim(C), c(n, 200L)) || any(colSums(C) != n)) fatal("Count generation failed.")
  files$counts <- save_new(list(counts = C, seed = count_seed, rng_kind = RNGkind()), "counts.rds")
  data_hash <- objsha(dat); counts_hash <- objsha(C)
  source_data <- list(Y = dat$Y, Z = as.numeric(dat$Z), weights = dat$W)
  work <- dat; work$A <- dat$Z; work$survey_weight <- dat$W
  formula_ps <- A ~ X1 + X2 + X3 + X4 + X5 + X6 + X1:X2
  D <- stats::model.matrix(~ X1 + X2 + X3 + X4 + X5 + X6 + X1:X2, work)
  raw_X <- as.matrix(dat[paste0("X", 1:6)])
  SX <- original_helpers$ows_standardize(raw_X); QX <- original_helpers$ows_quadratic(SX)
  controls <- wm$.wm_wdsm_controls(list())
  # Literal legacy PS fitting convention, retaining the actual coefficients.
  # One fit supplies both estimands; it is not a scalar-specific inference path.
  ps <- capture(function() {
    fit <- stats::glm(formula_ps, data = work, weights = work$survey_weight,
      family = stats::quasibinomial(link = "logit"), start = rep(0, ncol(D)),
      na.action = stats::na.fail, control = stats::glm.control(maxit = 100L, epsilon = 1e-10))
    e <- as.numeric(stats::predict(fit, newdata = work, type = "response")); beta <- stats::coef(fit)
    if (!isTRUE(fit$converged) || fit$rank != ncol(D) || any(!is.finite(beta)) ||
        any(!is.finite(e)) || any(e <= .Machine$double.eps | e >= 1-.Machine$double.eps))
      stop("PS fit did not select a finite nonsaturated full-rank root.")
    eq <- drop(crossprod(D, work$survey_weight * (work$A-e)))
    scale <- pmax(colSums(work$survey_weight * abs(D)), .Machine$double.eps * sum(work$survey_weight))
    residual <- max(abs(eq)/scale)
    if (!is.finite(residual) || residual > 1e-7) stop("PS normalized score exceeds original 1e-7 threshold.")
    J <- -crossprod(D, D * (work$survey_weight * e * (1-e))) / n
    step <- solve(J, eq/n)
    list(probability = e, coefficients = beta,
      scores = original_helpers$ows_standardize(cbind(ps = e)),
      diagnostics = list(converged = fit$converged, iterations = fit$iter, rank = fit$rank,
        normalized_score = residual, normalized_score_tolerance = 1e-7,
        newton_step = step, newton_relative_residual = max(abs(step)/(1+abs(beta))),
        exact_statistical_root_certified = FALSE))
  })
  files$ps <- save_new(ps, "original/shared_PS_root.rds")
  for (estimand in c("PATE", "PATT")) for (family in c("PS", "DSM", "X6")) {
    key <- paste(family, estimand, sep = "_"); guard(paste0(key, ":original"))
    pate <- identical(estimand, "PATE")
    expected_p <- switch(family, PS = if (pate) 26L else 18L,
      DSM = if (pate) 42L else 26L, X6 = if (pate) 68L else 40L)
    args_dsm <- list(Y = dat$Y, Z = as.numeric(dat$Z), weights = dat$W,
      ps_design = D, pg0_design = D, pg1_design = if (pate) D else NULL,
      estimand = estimand, ps_weighting = "probability")
    # Separate point construction from inference: inference failures retain point/stack.
    original_models <- capture(function() {
      if (family == "DSM") {
        solved <- wm$.wm_wdsm_fit_stack(args_dsm, wm$.wm_wdsm_nuisance_stack, controls)
        return(list(stack = solved$stack, original_solver = solved$attempts,
                    coefficients = solved$stack$parameter))
      }
      design <- if (family == "PS") D else QX
      q <- original_helpers$ows_full_x_predictions(work, design, rep(1, n), estimand)
      pscore <- if (family == "PS") require_value(ps) else NULL
      labels <- if (pate) c("mean0", "mean1") else "mean0"
      models <- stats::setNames(lapply(labels, function(label) list(design = design,
        coefficients = q$coefficients[[if (label == "mean0") "arm0" else "arm1"]])), labels)
      list(coefficients = q$coefficients, mean_models = models, predictions = q,
           pscore = pscore, stack = NULL)
    })
    original <- capture(function() {
      base <- require_value(original_models)
      if (family == "DSM") {
        stack <- base$stack
        point <- wm$wm_match(dat$Y, dat$Z, dat$W, stack$scores0,
          scores1 = stack$scores1, M = 3L, estimand = estimand,
          mean0 = stack$mean0, mean1 = stack$mean1, variance = FALSE)
        inf_args <- list(fit = point, nuisance_influence = stack$nuisance_influence,
          mean_derivative0 = stack$mean0_derivative, mean_derivative1 = stack$mean1_derivative,
          weight_derivative = stack$weight_derivative,
          covariance_scope = if (pate) "full_x" else "patt")
        if (!pate) inf_args$transport0 <- list(mode = "zero", basis = "current_centering")
        base$point <- point; base$inference_arguments <- inf_args
        return(base)
      }
      q <- base$predictions
      scores <- if (family == "PS") base$pscore$scores else SX
      point <- wm$wm_match(dat$Y, dat$Z, dat$W, scores,
        scores1 = if (pate) scores else NULL, M = 3L, estimand = estimand,
        mean0 = q$mean0, mean1 = q$mean1, variance = FALSE)
      # Assembly is captured separately below so the original point is never discarded.
      base$point <- point; base$inference_arguments <- NULL
      base
    })
    assembled <- capture(function() {
      base <- require_value(original)
      if (family == "DSM") return(base$stack)
      if (family == "PS") stack_helpers$wm_saved_full_x_stack(source_data, base$point,
        base$mean_models, family = "PS", ps = list(design = D,
          probability = base$pscore$probability, coefficients = base$pscore$coefficients,
          weighting = "probability"), wm = wm)
      else stack_helpers$wm_saved_full_x_stack(source_data, base$point,
        base$mean_models, family = "full_X", raw_X = raw_X, wm = wm)
    })
    inf <- capture(function() {
      base <- require_value(original); stack <- require_value(assembled)
      ia <- if (family == "DSM") base$inference_arguments else stack$inference_arguments
      if (ncol(ia$nuisance_influence) != expected_p ||
          !identical(colnames(ia$nuisance_influence), stack$parameter_names)) fatal("Complete p/order mismatch.")
      result <- do.call(wm$wm_fitted_inference, ia)
      if (!identical(result$fit, base$point) || !identical(result$estimate, base$point$estimate))
        fatal("Inference changed original point/graph.")
      ws <- result$smooth_derivative$weight_sensitivity
      if (length(ws) != expected_p || any(!is.finite(ws)) || any(ws != 0))
        fatal("Known-W sensitivity is missing or not zero.")
      result
    })
    # One saved original point, full stack and inference; no new calculation.
    retained_inf <- inf
    if (isTRUE(retained_inf$ok)) retained_inf$value$fit <- NULL
    checkpoint <- list(point = if (original$ok) original$value$point else NULL,
      model_status = original_models[c("ok", "error", "warnings", "attempts")],
      point_status = original[c("ok", "error", "warnings")],
      stack = assembled, inference = retained_inf, family = family,
      estimand = estimand, p = expected_p, sample_sha256 = data_hash,
      counts_sha256 = counts_hash, inference_fit_is_saved_point = TRUE)
    files[[paste0(key, "_original")]] <- save_new(checkpoint, paste0("original/", key, ".rds"))
    rm(checkpoint)
    guard(paste0(key, ":contribution"))
    contribution <- capture(function() {
      result <- wm$wm_bootstrap(require_value(inf), B = 200L, counts = C,
        method = "contribution", interval = "normal")
      check_close(result$conditional_variance, inf$value$variance, "Contribution/analytic variance")
      check_close(result$conf.int, inf$value$conf.int, "Contribution/analytic interval")
      if (!identical(result$estimate, original$value$point$estimate)) fatal("Contribution point changed.")
      result$source_inference <- NULL; result$supplied_counts <- NULL
      result$source_inference_artifact <- files[[paste0(key, "_original")]]
      result$counts_artifact <- files$counts
      result
    })
    cache <- vector("list", 200L); cache_files <- vector("list", 4L)
    compact_slots <- vector("list", 50L)
    dir.create(file.path(out, "refit_cache", key), recursive = TRUE)
    manual_draws <- rep(NA_real_, 200L)
    draw_errors <- rep("Not evaluated", 200L)
    for (b in seq_len(200L)) {
      guard(paste0(key, ":count", b, ":before"))
      cache[[b]] <- capture(function() {
        if (family == "DSM") {
          solved <- wm$.wm_wdsm_fit_stack(c(args_dsm, list(multiplicity = C[, b])),
            wm$.wm_wdsm_prediction_stack, controls)
          q <- solved$stack
          list(mean0 = q$mean0, mean1 = q$mean1, coefficients = q$parameter,
               diagnostics = q$diagnostics, attempts = solved$attempts)
        } else original_helpers$ows_full_x_predictions(work,
          if (family == "PS") D else QX, C[, b], estimand)
      })
      cache[[b]]$column <- b; cache[[b]]$count_sha256 <- objsha(C[, b])
      draw <- capture(function() {
        q <- require_value(cache[[b]]); point <- require_value(original)$point
        original_helpers$ows_wm_replicate(point, C[, b], q$mean0, q$mean1)
      })
      if (draw$ok) { manual_draws[b] <- draw$value; draw_errors[b] <- "" }
      else draw_errors[b] <- draw$error
      compact_slots[[(b - 1L) %% 50L + 1L]] <- list(
        prediction_fit = .principal_pack_count(cache[[b]]), manual_scalar = draw)
      if (b %% 50L == 0L) {
        chunk <- b %/% 50L
        cache_files[[chunk]] <- save_new(list(schema = "principal_iid_count_chunk_v1",
          columns = seq.int(b - 49L, b), slots = compact_slots,
          original_artifact = files[[paste0(key, "_original")]],
          sample_sha256 = data_hash, counts_sha256 = counts_hash),
          sprintf("refit_cache/%s/columns%03d-%03d.rds", key, b - 49L, b))
        files[[sprintf("%s_count_chunk%02d", key, chunk)]] <- cache_files[[chunk]]
        compact_slots <- vector("list", 50L)
      }
      guard(paste0(key, ":count", b, ":after"))
    }
    files[[paste0(key, "_cache")]] <- save_new(list(schema = "principal_iid_count_cache_v1", chunk_files = cache_files,
      manual_draws = manual_draws, draw_errors = draw_errors,
      sample_sha256 = data_hash, counts_sha256 = counts_hash,
      no_survivor_variance = TRUE), paste0("refit_cache/", key, ".rds"))
    invoked <- 0L
    callback <- function(m) {
      invoked <<- invoked + 1L
      if (invoked > 200L || !identical(as.numeric(m), as.numeric(C[, invoked])))
        fatal("Public refit callback count/order mismatch.")
      q <- require_value(cache[[invoked]])
      list(mean0 = q$mean0, mean1 = q$mean1)
    }
    refit <- capture(function() wm$wm_bootstrap(require_value(original)$point, B = 200L,
      counts = C, method = "fixed_reuse", refit = callback, interval = "normal"))
    if (refit$ok) {
      if (invoked != 200L || any(draw_errors != "") || any(!is.finite(manual_draws)))
        fatal("Public refit success disagrees with all-slot status.")
      check_close(refit$value$draws, manual_draws, "Public/manual fixed-reuse draws")
      if (!identical(refit$value$variance_divisor, "B") ||
          !identical(refit$value$estimate, original$value$point$estimate)) fatal("Refit divisor/point mismatch.")
    }
    public_scalar_status <- rep("not_attempted", 200L)
    if (refit$ok) public_scalar_status[] <- "completed"
    else if (invoked > 0L) {
      if (invoked > 1L) public_scalar_status[seq_len(invoked - 1L)] <- "completed_before_next_callback"
      public_scalar_status[invoked] <- "not_returned_public_error"
    }
    files[[paste0(key, "_replication")]] <- save_new(list(contribution = contribution,
      original_refit = refit, callback_columns_consumed = invoked,
      public_scalar_status = public_scalar_status,
      manual_draws_are_separate_algebra_checks = TRUE,
      all_200_prediction_fit_slots_preserved = TRUE,
      cache_artifact = files[[paste0(key, "_cache")]]), paste0("replication/", key, ".rds"))
    target <- target_rows$target[match(estimand, target_rows$estimand)]
    records[[key]] <- list(key = key, case_id = case_id, overlap = overlap, estimand = estimand, family = family,
      n = n, replicate = replicate_id, group_index = group_index,
      M = 3L, p = expected_p, B = 200L, target = target,
      original_artifact = files[[paste0(key, "_original")]],
      point_status = if (original$ok) "ok" else "failed", point_error = original$error,
      estimate = if (original$ok) original$value$point$estimate else NA_real_,
      raw_estimate = if (original$ok) original$value$point$raw_estimate else NA_real_,
      analytic_available = inf$ok && isTRUE(inf$value$available), analytic_error = inf$error,
      analytic_variance = if (inf$ok) inf$value$variance else NA_real_,
      analytic_interval = if (inf$ok) inf$value$conf.int else c(NA_real_, NA_real_),
      contribution_ok = contribution$ok, contribution_error = contribution$error,
      refit_available = refit$ok, refit_error = refit$error,
      refit_variance = if (refit$ok) refit$value$variance else NA_real_,
      refit_interval = if (refit$ok) refit$value$conf.int else c(NA_real_, NA_real_),
      count_fit_failures = which(!vapply(cache, function(x) isTRUE(x$ok), logical(1L))),
      count_draw_failures = which(draw_errors != ""), sample_sha256 = data_hash,
      counts_sha256 = counts_hash, population_premises_verified_by_arrays = FALSE)
    write_json(records[[key]], file.path(out, paste0(key, ".json")))
    if (!identical(data_hash, objsha(dat)) || !identical(counts_hash, objsha(C))) fatal("Data/count mutation.")
    rm(original_models, original, assembled, inf, contribution, cache, refit); gc(verbose = FALSE)
    guard(paste0(key, ":completed"))
  }
  if (length(records) != 6L) fatal("Six planned methods were not all recorded.")
  post_count_RNG_unchanged <- identical(rng_after_counts,
    list(kind = RNGkind(), seed = get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)))
  if (!post_count_RNG_unchanged) fatal("Post-count scientific calls changed RNG state.")
  status <- if (all(vapply(records, function(x) isTRUE(x$analytic_available) &&
    isTRUE(x$contribution_ok) && isTRUE(x$refit_available), logical(1L))))
    "COMPLETE_CASE" else "COMPLETE_WITH_STATISTICAL_FAILURES"
  list(schema = "principal_iid_public_case_v1", status = status, case_id = case_id,
    overlap = overlap, n = n, replicate = replicate_id, group_index = group_index,
    data_seed = data_seed, count_seed = count_seed, M = 3L, B = 200L,
    rows = records, files = files, post_count_RNG_unchanged = post_count_RNG_unchanged,
    population_premises_verified_by_arrays = FALSE)
}
