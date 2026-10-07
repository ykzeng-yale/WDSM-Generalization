# Sourceable first-stage validation DGP and paired-run helpers.
# Source package R files (or load WeightedMatching) before calling the run helper.

.wm_fs_integer <- function(x, name, lower = 1L, upper = .Machine$integer.max) {
  if (!is.numeric(x) || is.complex(x) || length(x) != 1L ||
      !is.finite(x) || x != floor(x) || x < lower || x > upper)
    stop(name, " must be an integer between ", lower, " and ", upper)
  as.integer(x)
}

wm_first_stage_pilot <- function() {
  data.frame(scenario = c("scores_same_d2_n200", "scores_independent_d1_n200",
                          "weights_d3_n200"),
    branch = c("same_sample", "independent_training", "estimated_weights"),
    n = 200L, m = 200L, d = c(2L, 1L, 3L), M = 3L,
    seed = c(929101L, 929102L, 929103L), replications = 40L,
    stringsAsFactors = FALSE)
}

wm_first_stage_config <- function(config) {
  fields <- c("scenario", "branch", "n", "m", "d", "M", "seed", "replications")
  if (!is.data.frame(config) || nrow(config) < 1L || nrow(config) > 100L ||
      !identical(sort(names(config)), sort(fields)) ||
      anyNA(config$scenario) || !is.character(config$scenario) ||
      any(!nzchar(config$scenario)) || anyDuplicated(config$scenario) ||
      anyNA(config$branch) ||
      any(!config$branch %in% c("same_sample", "independent_training", "estimated_weights", "gaussian_same_sample")))
    stop("invalid first-stage configuration or duplicate scenario")
  for (row in seq_len(nrow(config))) {
    for (name in c("n", "m", "d", "M", "seed", "replications")) {
      config[[name]][row] <- .wm_fs_integer(config[[name]][row], name,
        lower = if (name == "seed") 0L else if (name %in% c("n", "m")) 2L else 1L,
        upper = if (name == "replications") 100000L else .Machine$integer.max)
    }
    if (as.double(config$n[row]) * config$d[row] > 1e7 ||
        as.double(config$m[row]) * config$d[row] > 1e7)
      stop("first-stage allocation guard exceeded")
    if (config$branch[row] != "independent_training" && config$m[row] != config$n[row])
      stop("same-sample and estimated-weight configurations require m = n")
    if (config$branch[row] == "same_sample" && config$d[row] < 2L)
      stop("same-sample changing scalar scores are outside the declared theorem")
  }
  if (anyDuplicated(config$seed)) stop("first-stage scenarios require distinct base seeds")
  config
}

wm_first_stage_stream <- function(base_seed, replication) {
  base_seed <- .wm_fs_integer(base_seed, "base_seed", 0L)
  replication <- .wm_fs_integer(replication, "replication", 1L, 100000L)
  kind <- RNGkind()
  if (identical(kind[2L], "Box-Muller")) stop("use a noncached normal RNG before indexing streams")
  had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  saved <- if (had_seed) get(".Random.seed", .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    do.call(RNGkind, as.list(kind))
    if (had_seed) assign(".Random.seed", saved, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  })
  set.seed(base_seed, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion",
           sample.kind = "Rejection")
  stream <- .Random.seed
  if (replication > 1L) for (i in seq_len(replication - 1L))
    stream <- parallel::nextRNGStream(stream)
  stream
}

wm_first_stage_generate <- function(n, d,
    branch = c("same_sample", "independent_training", "estimated_weights", "gaussian_same_sample")) {
  branch <- match.arg(branch)
  n <- .wm_fs_integer(n, "n", 2L)
  d <- .wm_fs_integer(d, "d")
  if (as.double(n) * d > 1e7) stop("first-stage allocation guard exceeded")
  S <- matrix(stats::runif(n * d), n, d)
  epsilon <- if (branch == "gaussian_same_sample") stats::rnorm(n, sd = .25) else
    0.25 * (2 * stats::rbinom(n, 1L, 0.5) - 1)
  if (branch == "estimated_weights") {
    t <- S[, 1L] - 0.5
    T <- cbind(intercept = 1, t = t)
    q <- stats::plogis(3 * t)
    Z <- stats::rbinom(n, 1L, q)
    weights <- ifelse(Z == 1L, 0.5 / q, 0.5 / (1 - q))
    mean0 <- 1 + S[, 1L] + rowSums(S^2) / (2 * d)
    mean1 <- mean0 + 1 + S[, 1L]
    D <- U <- NULL
    target <- 1.5
  } else {
    U <- 2 * stats::rbinom(n, 1L, 0.5) - 1
    q <- (1 + 0.6 * U) / 2
    Z <- stats::rbinom(n, 1L, q)
    weights <- if (branch == "gaussian_same_sample") 1 + Z else rep(1, n)
    mean0 <- S[, 1L]
    mean1 <- mean0 + 1
    D <- U * S[, 1L] * (1 - S[, 1L])
    T <- NULL
    target <- 1
  }
  list(Y = ifelse(Z == 1L, mean1, mean0) + epsilon, Z = Z, scores = S,
       weights = weights, mean0 = mean0, mean1 = mean1, epsilon = epsilon,
       U = U, D = D, design = T, propensity = q, target = target,
       sigma2 = 1 / 16, r = if (branch == "estimated_weights") NA_real_ else 0.6,
       branch = branch)
}

wm_first_stage_score_fit <- function(dat) {
  offset_residual <- dat$Y - dat$scores[, 1L] - dat$Z
  gram <- mean(dat$D^2)
  if (!is.finite(gram) || gram <= 0) stop("nonpositive first-stage prediction Gram")
  theta <- mean(dat$D * offset_residual) / gram
  if (!is.finite(theta) || abs(theta) >= 0.8)
    stop("prediction parameter failed the declared abs(theta) < 0.8 interior rule")
  influence <- matrix(dat$D * (offset_residual - theta * dat$D) / gram,
                       ncol = 1L, dimnames = list(NULL, "theta"))
  if (any(!is.finite(influence))) stop("nonfinite prediction influence")
  list(parameters = c(theta = theta), nuisance_influence = influence,
       gram = gram, parameter_bound = 0.9, accepted_bound = 0.8,
       estimating_equation = mean(dat$D * (offset_residual - theta * dat$D)))
}

.wm_fs_capture <- function(fun) {
  warnings <- character()
  started <- proc.time()[["elapsed"]]
  value <- tryCatch(withCallingHandlers(fun(), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = identity)
  ok <- !inherits(value, "error")
  list(ok = ok, value = if (ok) value else NULL,
       error = if (ok) "" else conditionMessage(value),
       warning = paste(warnings, collapse = " | "),
       elapsed = proc.time()[["elapsed"]] - started)
}

.wm_fs_failed <- function(message) {
  list(ok = FALSE, value = NULL, error = message, warning = "", elapsed = 0)
}

.wm_fs_match <- function(dat, estimand, M, scores = dat$scores,
                          weights = dat$weights, mean0 = dat$mean0,
                          mean1 = dat$mean1) {
  wm_match(dat$Y, dat$Z, weights, scores, M = M, estimand = estimand,
           mean0 = mean0, mean1 = if (estimand == "PATE") mean1 else NULL)
}

# OLS computed by moments and by QR can differ at floating-point precision.
# Close neighbors amplify this harmless discrepancy if distances are compared
# relative to the distance itself. Keep graph identities/fractions exact and
# coordinates at the original tolerance, and validate distances on their
# coordinate scale. This diagnostic never changes a score, edge, or estimate.
.wm_fs_gaussian_same_graph <- function(a, b) {
  if (!identical(attributes(a), attributes(b)) ||
      !identical(attributes(a$edges), attributes(b$edges))) return(FALSE)
  other <- setdiff(names(a), c("scores0", "scores1", "edges"))
  edge_other <- setdiff(names(a$edges), "distance")
  if (!identical(a[other], b[other]) ||
      !identical(a$edges[edge_other], b$edges[edge_other])) return(FALSE)
  if (!identical(attributes(a$edges$distance), attributes(b$edges$distance))) return(FALSE)
  for (name in c("scores0", "scores1")) {
    for (value in list(a[[name]], b[[name]])) {
      if (!is.null(value) && (!is.matrix(value) || !is.numeric(value) ||
          is.complex(value) || any(!is.finite(value)))) return(FALSE)
    }
    if (!identical(attributes(a[[name]]), attributes(b[[name]])) ||
        !isTRUE(all.equal(a[[name]], b[[name]], tolerance = 1e-13))) return(FALSE)
  }
  norm_rows <- function(delta) {
    largest <- abs(delta[, 1L])
    if (ncol(delta) > 1L) for (column in 2:ncol(delta))
      largest <- pmax(largest, abs(delta[, column]))
    largest * sqrt(rowSums((delta / ifelse(largest > 0, largest, 1))^2))
  }
  for (arm in unique(a$edges$arm)) {
    at <- which(a$edges$arm == arm)
    name <- if (arm == 0L) "scores0" else "scores1"
    sa <- a[[name]]; sb <- b[[name]]
    if (!is.matrix(sa) || !is.matrix(sb) || !identical(dim(sa), dim(sb))) return(FALSE)
    query <- a$edges$query[at]; donor <- a$edges$donor[at]
    aq <- sa[query, , drop = FALSE]; aj <- sa[donor, , drop = FALSE]
    bq <- sb[query, , drop = FALSE]; bj <- sb[donor, , drop = FALSE]
    da <- a$edges$distance[at]; db <- b$edges$distance[at]
    # Recompute each stored edge using the same scaled Euclidean arithmetic.
    # Even a small corruption of a recorded distance must be rejected.
    if (!identical(unname(da), unname(norm_rows(aj - aq))) ||
        !identical(unname(db), unname(norm_rows(bj - bq)))) return(FALSE)
    coordinate_change <- norm_rows(aq - bq) + norm_rows(aj - bj)
    # Arithmetic allowance for this bounded simulation design, not a claim
    # certifying arbitrary floating-point inputs or extreme dimensions.
    roundoff <- 8 * ncol(sa) * .Machine$double.eps *
      (norm_rows(aq) + norm_rows(aj) + norm_rows(bq) + norm_rows(bj))
    if (any(!is.finite(coordinate_change + roundoff)) ||
        any(abs(da - db) > coordinate_change + roundoff)) return(FALSE)
  }
  TRUE
}

wm_first_stage_run_rep <- function(config, replication, stream = NULL,
                                   retain_fits = FALSE) {
  cfg <- wm_first_stage_config(config)
  if (nrow(cfg) != 1L) stop("select exactly one first-stage scenario")
  replication <- .wm_fs_integer(replication, "replication", 1L, cfg$replications)
  if (!is.logical(retain_fits) || length(retain_fits) != 1L || is.na(retain_fits))
    stop("retain_fits must be TRUE or FALSE")
  if (is.null(stream)) stream <- wm_first_stage_stream(cfg$seed, replication)
  if (!is.integer(stream) || length(stream) != 7L || anyNA(stream) || stream[1L] != 10407L)
    stop("stream must be a L'Ecuyer-CMRG/Inversion/Rejection state from wm_first_stage_stream")
  kind <- RNGkind()
  if (identical(kind[2L], "Box-Muller")) stop("use a noncached normal RNG before running")
  had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  saved <- if (had_seed) get(".Random.seed", .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    do.call(RNGkind, as.list(kind))
    if (had_seed) assign(".Random.seed", saved, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  })
  assign(".Random.seed", stream, .GlobalEnv)
  training_stream <- if (cfg$branch == "independent_training")
    parallel::nextRNGSubStream(stream) else NULL
  generated <- .wm_fs_capture(function() {
    evaluation <- wm_first_stage_generate(cfg$n, cfg$d, cfg$branch)
    training <- if (is.null(training_stream)) evaluation else {
      assign(".Random.seed", training_stream, .GlobalEnv)
      wm_first_stage_generate(cfg$m, cfg$d, cfg$branch)
    }
    list(evaluation = evaluation, training = training)
  })
  stage <- if (!generated$ok) .wm_fs_failed(paste("Data generation:", generated$error)) else
    .wm_fs_capture(function() {
      dat <- generated$value$training
      if (cfg$branch == "estimated_weights")
        wm_logistic_weights(dat$Z, dat$design, pi = 0.5, parameter_bound = 5) else
          wm_first_stage_score_fit(dat)
    })
  records <- list()
  fitted_objects <- list()
  for (estimand in c("PATE", "PATT")) {
    dat <- if (generated$ok) generated$value$evaluation else NULL
    known <- if (!generated$ok) .wm_fs_failed(generated$error) else
      .wm_fs_capture(function() .wm_fs_match(dat, estimand, cfg$M))
    naive <- if (!stage$ok) .wm_fs_failed(paste("First stage:", stage$error)) else
      .wm_fs_capture(function() {
        if (cfg$branch == "estimated_weights")
          .wm_fs_match(dat, estimand, cfg$M, weights = stage$value$weights) else {
            theta <- unname(stage$value$parameters[1L])
            current_scores <- dat$scores
            current_scores[, 1L] <- current_scores[, 1L] + theta * dat$D
            .wm_fs_match(dat, estimand, cfg$M, scores = current_scores,
                          mean0 = dat$mean0 + theta * dat$D,
                          mean1 = dat$mean1 + theta * dat$D)
          }
      })
    adjusted <- if (!naive$ok) .wm_fs_failed(paste("Fitted matching:", naive$error)) else
      .wm_fs_capture(function() {
        fit <- naive$value
        if (cfg$branch == "estimated_weights") {
          result <- wm_weight_adjust(fit, stage$value$weight_derivative,
                                       stage$value$nuisance_influence)
        } else {
          D <- matrix(dat$D, ncol = 1L, dimnames = list(NULL, "theta"))
          if (cfg$branch == "gaussian_same_sample") {
            result <- wm_gaussian_match(dat$Y, dat$Z, dat$weights, D,
              score_map = function(theta) {
                if (abs(theta[1L]) >= .8) stop("Gaussian score parameter outside accepted interior")
                current <- dat$scores; current[, 1L] <- current[, 1L] + theta[1L] * dat$D
                current
              }, offset0 = dat$mean0, offset1 = dat$mean1, M = cfg$M, estimand = estimand)
          } else if (cfg$branch == "same_sample")
            result <- wm_prediction_adjust(fit, stage$value$nuisance_influence,
              mean_derivative0 = D, mean_derivative1 = if (estimand == "PATE") D else NULL,
              changing_scores = rep(TRUE, if (estimand == "PATE") 2L else 1L)) else
                result <- wm_training_adjust(fit, stage$value$nuisance_influence,
                  mean_derivative0 = D, mean_derivative1 = if (estimand == "PATE") D else NULL)
        }
        same_graph <- if (cfg$branch == "gaussian_same_sample")
          .wm_fs_gaussian_same_graph(result$graph, fit$graph) else
          identical(result$graph, fit$graph)
        stopifnot(isTRUE(all.equal(result$estimate, fit$estimate, tolerance = 1e-13)), same_graph)
        result
      })
    comparisons <- list(known_first_stage = known, fitted_naive = naive,
                         fitted_adjusted = adjusted)
    for (comparison in names(comparisons)) {
      result <- comparisons[[comparison]]
      fit <- result$value
      ok <- result$ok
      estimate <- if (ok) fit$estimate else NA_real_
      se <- if (ok) fit$se else NA_real_
      sensitivity <- if (ok && comparison == "fitted_adjusted") {
        if (cfg$branch == "estimated_weights") fit$weight_adjustment$sensitivity else
          if (cfg$branch == "gaussian_same_sample") fit$gaussian_adjustment$sensitivity else
          if (cfg$branch == "same_sample") fit$prediction_adjustment$sensitivity else
            fit$training_adjustment$sensitivity
      } else numeric()
      parameters <- if (stage$ok) as.numeric(stage$value$parameters) else numeric()
      records[[length(records) + 1L]] <- data.frame(
        schema_version = 1L, scenario = cfg$scenario, replication = replication,
        branch = cfg$branch, n = cfg$n, m = cfg$m, d = cfg$d, M = cfg$M,
        base_seed = cfg$seed, evaluation_rng_state = paste(stream, collapse = ":"),
        training_rng_state = paste(training_stream, collapse = ":"),
        estimand = estimand, comparison = comparison,
        pairing_id = paste(cfg$scenario, replication, estimand, sep = "/"),
        target = if (cfg$branch == "estimated_weights") 1.5 else 1,
        status = if (ok) "ok" else "error",
        inference_status = if (comparison == "fitted_naive")
          "first_stage_omission_diagnostic" else "within_declared_parametric_scope",
        estimate = estimate, root_n_variance = if (ok) fit$root_n_variance else NA_real_,
        sampling_variance = if (ok) fit$variance else NA_real_,
        lower = estimate - stats::qnorm(0.975) * se,
        upper = estimate + stats::qnorm(0.975) * se,
        first_stage_status = if (stage$ok) "ok" else "error",
        parameter1 = if (length(parameters)) parameters[1L] else NA_real_,
        parameter2 = if (length(parameters) > 1L) parameters[2L] else NA_real_,
        sensitivity1 = if (length(sensitivity)) unname(sensitivity[1L]) else NA_real_,
        sensitivity2 = if (length(sensitivity) > 1L) unname(sensitivity[2L]) else NA_real_,
        sigma2 = 1 / 16, r = if (cfg$branch == "estimated_weights") NA_real_ else 0.6,
        evaluation_training_ratio = if (cfg$branch == "independent_training") cfg$n / cfg$m else NA_real_,
        elapsed_seconds = result$elapsed, first_stage_seconds = stage$elapsed,
        generation_seconds = generated$elapsed, warning = result$warning, error = result$error,
        first_stage_warning = stage$warning, first_stage_error = stage$error,
        stringsAsFactors = FALSE)
      if (retain_fits) fitted_objects[paste(estimand, comparison, sep = "/")] <- list(fit)
    }
  }
  list(records = do.call(rbind, records), fits = if (retain_fits) fitted_objects else NULL,
       diagnostics = list(first_stage_ok = stage$ok,
         parameters = if (stage$ok) stage$value$parameters else NULL,
         score_gram = if (stage$ok) stage$value$gram else NULL,
         estimating_equation = if (stage$ok) stage$value$estimating_equation else NULL,
         evaluation_rng_state = stream, training_rng_state = training_stream,
         parameter_truth = if (cfg$branch == "estimated_weights") c(0, 3) else 0))
}
