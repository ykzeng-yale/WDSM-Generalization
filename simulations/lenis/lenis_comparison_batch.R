# Six literal Lenis outcome multipliers share all covariates, treatment,
# response indicators, control outcomes, and matching stages. Source identities
# are verified by cache_lenis_populations.R; this is PATT-specific batching.

lenis_find_assignment <- function(x, name) {
  if (missing(x)) return(NULL)
  if (is.call(x) && as.character(x[[1L]])[1L] %in% c("<-", "=") &&
      identical(x[[2L]], as.name(name))) return(x)
  if (is.call(x) || is.expression(x)) for (y in as.list(x)) {
    found <- lenis_find_assignment(y, name)
    if (!is.null(found)) return(found)
  }
  NULL
}

lenis_source_sample <- function(cache, seed,
    source = "research/source-cache/lenis-2019/biosts-17039-File013.R") {
  # Only design fields enter the literal source sampling expressions. Keeping
  # the frame narrow avoids copying all potential outcomes through mstage.
  frame <- cache$common[c("Cluster", "Strata")]
  frame$ones <- 1; frame$id <- seq_len(nrow(frame))
  env <- new.env(parent = globalenv()); env$data <- frame
  src <- parse(source)
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rounding"))
  for (name in c("size1", "size2", "size3", "s", "sample")) {
    e <- lenis_find_assignment(src, name); stopifnot(!is.null(e))
    if (name == "s") set.seed(seed)
    eval(e, env)
  }
  list(ids = env$sample$id, Prob = env$sample$Prob,
    design_fields = env$sample[setdiff(names(env$sample), names(frame))],
    rng_after = .Random.seed, seed = seed,
    frame_sha256 = digest::digest(frame, algo = "sha256"))
}

lenis_selected_from_cache <- function(cache, sample) {
  d <- cache$common[sample$ids, , drop = FALSE]
  y1 <- cache$y1[sample$ids, , drop = FALSE]
  Y <- y1 * d$z + d$y0 * (1 - d$z)
  d$y1 <- y1[, 1L]; d$y <- Y[, 1L]; d$Multiplier <- 1L
  d <- d[cache$original_columns]
  d$ones <- 1; d$id <- sample$ids
  for (k in names(sample$design_fields)) d[[k]] <- sample$design_fields[[k]]
  d$s.wt <- 1/d$Prob
  stopifnot(!anyDuplicated(d$id), nrow(d) == 5000L,
    all(Y[d$z == 0, , drop = FALSE] == d$y[d$z == 0]))
  list(selected = d, Y = Y, y1 = y1, targets = cache$truth$PATT)
}

lenis_numeric_fit <- function(fit) {
  if (is.null(fit)) return(NULL)
  list(coefficients = stats::coef(fit), covariance = stats::vcov(fit),
    rank = fit$rank, converged = fit$converged, boundary = fit$boundary)
}

lenis_compact_status <- function(x) {
  result <- x[intersect(c("ok", "error", "condition_class", "solver_attempts",
    "warnings", "elapsed", "refined"), names(x))]
  if (!is.null(x$partial_results)) result$partial_results <-
    x$partial_results[intersect(c("values", "diagnostics", "specification"),
                                names(x$partial_results))]
  result
}

lenis_native_outcome_batch <- function(native, response, Y, targets) {
  d <- response$data
  stopifnot(is.matrix(Y), nrow(Y) == nrow(d), ncol(Y) == length(targets),
            identical(as.numeric(Y[, 1L]), as.numeric(d$y)))
  output <- vector("list", length(native$results)); names(output) <- names(native$results)
  for (j in seq_along(native$results)) {
    first <- native$results[[j]]; s <- native$settings[j, ]
    output[[j]] <- vector("list", ncol(Y))
    for (q in seq_len(ncol(Y))) {
      value <- if (q == 1L) first else if (!first$ok) lenis_capture(function() {
        # A failed first outcome does not imply failure of different treated
        # outcomes. Rare fallback evaluates the prescribed later outcome; it
        # neither retries q=1 nor changes data, counts, or estimator settings.
        changed <- response; changed$data$y <- Y[, q]
        later_case <- lenis_native_case(changed, targets[q])
        later <- later_case$results[[j]]
        for (message in later$warnings) warning(message, call. = FALSE)
        if (!later$ok) stop(structure(list(message = later$error, call = NULL,
          partial_results = later$partial_results),
          class = c(later$condition_class, "error", "condition")))
        value <- later$value
        value$fallback_matching <- lenis_compact_native(later_case)
        value
      }) else lenis_capture(function() {
        md <- first$value$matched_data
        rows <- if (s$transfer == "WT") match(as.character(md$id), rownames(d)) else
          match(as.character(md$id), as.character(d$id))
        stopifnot(!anyNA(rows), identical(as.numeric(md$y), as.numeric(Y[rows, 1L])))
        md$y <- Y[rows, q]
        design <- if (s$om == "U") survey::svydesign(ids = ~1, data = md) else
          survey::svydesign(ids = ~Cluster, strata = ~Strata, weights = ~s.wt, data = md)
        fits <- lapply(list(A = stats::reformulate(c(paste0("x", 1:6), "z"), "y"),
                           U = y ~ z), function(f) survey::svyglm(f, design))
        values <- first$value$values[startsWith(names(first$value$values), "SMD.")]
        method <- paste(paste0(s$ps, ".PS"), paste0(s$om, ".OM"),
                        s$transfer, response$mechanism, sep = "|")
        for (adjustment in names(fits)) {
          fit <- fits[[adjustment]]; estimate <- unname(stats::coef(fit)["z"])
          se <- unname(sqrt(diag(stats::vcov(fit)))["z"])
          triple <- c(estimate, se,
            as.numeric(estimate - 1.96 * se < targets[q] & estimate + 1.96 * se > targets[q]))
          names(triple) <- paste(c("ATT", "SE", "Cov"), adjustment, method, sep = "|")
          values <- c(values, triple)
        }
        lenis_validate_native(list(values = values, fits = fits, specification = s,
          diagnostics = list(ps = first$value$diagnostics$ps,
            outcome_converged = vapply(fits, function(f) isTRUE(f$converged), logical(1)),
            outcome_boundary = vapply(fits, function(f) isTRUE(f$boundary), logical(1)))))
      })
      compact <- lenis_compact_status(value)
      if (value$ok) {
        compact$values <- value$value$values
        compact$fits <- lapply(value$value$fits, lenis_numeric_fit)
        compact$diagnostics <- value$value$diagnostics
        if (!is.null(value$value$fallback_matching))
          compact$fallback_matching <- value$value$fallback_matching
      }
      output[[j]][[q]] <- compact
    }
  }
  output
}

lenis_quick_outcome_batch <- function(quick, response, Y, weighted_ps) {
  stopifnot(is.logical(weighted_ps), length(weighted_ps) == 1L, !is.na(weighted_ps))
  output <- vector("list", ncol(Y))
  for (q in seq_len(ncol(Y))) {
    fallback <- q > 1L && !quick$ok
    value <- if (q == 1L) quick else if (fallback) lenis_capture(function() {
      changed <- response; changed$data$y <- Y[, q]
      lenis_quick_case(changed, weighted_ps = weighted_ps)
    }) else lenis_capture(function() {
      base <- quick$value; md <- base$matched_data
      rows <- match(as.character(md$id), as.character(response$data$id))
      stopifnot(!anyNA(rows), identical(as.numeric(md$y), as.numeric(Y[rows, 1L])))
      md$y <- Y[rows, q]
      fit <- stats::lm(y ~ z * (x1 + x2 + x3 + x4 + x5 + x6), data = md, weights = weights)
      covariance <- sandwich::vcovCL(fit, cluster = md$subclass, type = "HC1")
      estimate <- sum(base$contrast * stats::coef(fit))
      variance <- as.numeric(crossprod(base$contrast, covariance %*% base$contrast))
      stopifnot(all(is.finite(c(estimate, variance))), variance >= 0)
      list(estimate = estimate, variance = variance, se = sqrt(variance),
        conf.int = estimate + c(-1, 1) * stats::qnorm(.975) * sqrt(variance), fit = fit,
        sandwich_covariance = covariance)
    })
    compact <- lenis_compact_status(value)
    if (value$ok) {
      compact$values <- value$value[c("estimate", "variance", "se", "conf.int")]
      compact$fit <- lenis_numeric_fit(value$value$fit)
      compact$sandwich_covariance <- if (q == 1L || fallback)
        sandwich::vcovCL(value$value$fit, cluster = value$value$matched_data$subclass,
                         type = "HC1") else value$value$sandwich_covariance
      if (fallback) compact$fallback_matching <- lenis_compact_quick(value)
    }
    output[[q]] <- compact
  }
  output
}

lenis_compact_native <- function(native) {
  list(settings = native$settings,
    first_values = lapply(native$results, function(x) if (x$ok) x$value$values else NULL),
    ps = lapply(native$ps, function(x) {
      result <- lenis_compact_status(x)
      if (x$ok) {
        result$fit <- lenis_numeric_fit(x$value$ps_fit)
        result$distance <- x$value$match$distance
        result$match_matrix <- x$value$match$match.matrix
        result$weights <- x$value$match$weights
        result$row_names <- rownames(x$value$data)
        result$population_ids <- x$value$data$id
      }
      result
    }),
    matching = lapply(native$results, function(x) {
      result <- lenis_compact_status(x)
      if (x$ok) result$data <- x$value$matched_data[
        intersect(c("id", "z", "s.wt", "weights", "X1", "Cluster", "Strata"),
                  names(x$value$matched_data))]
      result
    }))
}

lenis_compact_quick <- function(quick) {
  result <- lenis_compact_status(quick)
  if (quick$ok) {
    q <- quick$value
    result$distance <- q$match$distance
    result$matching_weights <- q$match$weights
    result$sampling_weights <- q$match$s.weights
    result$subclass <- q$match$subclass
    result$target_contrast <- q$contrast
    result$first_values <- q[c("estimate", "variance", "se", "conf.int")]
    result$population_ids <- q$matched_data$id
    result$analysis_weights <- q$matched_data$weights
    result$target_weights <- q$matched_data$s.wt[q$matched_data$z == 1]
  }
  result
}

lenis_compact_wm <- function(case) {
  hash <- function(x) digest::digest(x, algo = "sha256")
  list(correction = case$correction, scope = case$scope,
    nuisance_parameter = case$base$nuisance$parameter,
    nuisance_solver = case$base$solver,
    full_x_coefficients = case$full_x$coefficients,
    graphs = lapply(case$fits, function(f) list(M = f$M,
      estimate = f$estimate, raw_estimate = f$raw_estimate, edges = f$graph$edges,
      scores_sha256 = hash(f$graph$scores0), mean0_sha256 = hash(f$predictions$mean0),
      analysis_weights = f$analysis_weights, loads = f$loads$incoming)),
    scores = lapply(c("PS_M1", "DSM_M1", "X6_M1"), function(k) case$fits[[k]]$graph$scores0),
    base_predictions = list(DSM = case$base$nuisance$mean0, full_x = case$full_x$mean0),
    refit_prediction_sha256 = lapply(case$refit_means, hash),
    refit_diagnostics = case$refit_diagnostics)
}

lenis_assert_plain <- function(x) {
  if (is.environment(x) || is.function(x) || is.language(x) || inherits(x, "formula"))
    stop("Compact artifact contains an executable object or environment")
  if (is.list(x)) invisible(lapply(x, lenis_assert_plain))
  invisible(TRUE)
}
