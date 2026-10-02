# Standalone calibration engine. Source the four reviewed reference files first.
# No fitting or random generation occurs on source(). Not a package export.

.wm_cal_attempt <- function(data, model, estimand, multiplicity, controls,
                            fitter = wm_wdsm_nuisance_stack_reference) {
  began <- proc.time()[[3L]]
  warnings <- character()
  value <- tryCatch(withCallingHandlers(fitter(
    Y = data$Y, Z = data$Z, weights = data$weights,
    ps_design = if (model == "CorCor") data$ps_correct else data$ps_main,
    pg0_design = data$pg_correct,
    pg1_design = if (estimand == "PATE") data$pg_correct else NULL,
    estimand = estimand, multiplicity = multiplicity,
    glm_maxit = controls$maxit, glm_epsilon = controls$epsilon,
    rank_tolerance = 1e-10), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = identity)
  elapsed <- proc.time()[[3L]] - began
  if (inherits(value, "error")) {
    message <- conditionMessage(value)
    # This exact error is raised only after all structural stack checks pass.
    accuracy <- identical(message, "Weighted propensity root fails the normalized score check.")
    return(list(stack = NULL, record = list(status = "failed", error = message,
      accuracy_only = accuracy, controls = controls, elapsed_seconds = elapsed,
      warnings = warnings)))
  }
  step <- tryCatch(as.vector(solve(value$jacobian, value$mean_equation)), error = identity)
  inverse_ok <- !inherits(step, "error") && all(is.finite(step)) &&
    is.finite(value$diagnostics$jacobian_reciprocal_condition) &&
    value$diagnostics$jacobian_reciprocal_condition > 0
  residual <- if (inverse_ok) max(abs(step)/(1 + abs(value$parameter))) else Inf
  threshold <- min(1e-9, 1/length(data$Y))
  score <- value$diagnostics$normalized_ps_score
  okay <- inverse_ok && is.finite(score) && score <= 1e-10 && residual <= threshold
  record <- list(status = if (okay) "completed" else "failed",
    error = if (okay) "" else if (!inverse_ok) "Nonfinite inverse-Jacobian diagnostic" else
      "Root accuracy criterion failed", accuracy_only = !okay && inverse_ok && is.finite(score),
    controls = controls, elapsed_seconds = elapsed, warnings = warnings,
    parameter = value$parameter, mean_equation = value$mean_equation,
    newton_step = if (inverse_ok) step else NULL,
    newton_residual = residual, newton_threshold = threshold,
    diagnostics = value$diagnostics)
  list(stack = if (okay) value else NULL, record = record)
}

# Exactly one accuracy-only refinement. Both attempt records precede use of the fit.
wm_cal_fit <- function(data, model, estimand, multiplicity = rep(1, length(data$Y)),
                       save_attempt = function(i, x) invisible(NULL),
                       fitter = wm_wdsm_nuisance_stack_reference) {
  stopifnot(model %in% c("CorCor", "MisCor"), estimand %in% c("PATE", "PATT"))
  attempts <- list()
  for (i in 1:2) {
    controls <- if (i == 1L) list(maxit = 100L, epsilon = 1e-12) else
      list(maxit = 200L, epsilon = 1e-14)
    result <- .wm_cal_attempt(data, model, estimand, multiplicity, controls, fitter)
    attempts[[i]] <- result$record
    save_attempt(i, result$record)
    if (!is.null(result$stack)) return(list(status = "completed", error = "",
      stack = result$stack, attempts = attempts, refined = i == 2L))
    if (!isTRUE(result$record$accuracy_only)) break
  }
  list(status = "failed", error = attempts[[length(attempts)]]$error,
       stack = NULL, attempts = attempts, refined = length(attempts) == 2L)
}

.wm_cal_contract <- function() list(mode = "declared_full_x",
  justification = paste("Predeclared source_good r=.25 iid selected law; known W=1/pi_Z;",
    "correct source PG and full quadratic correction. Proofs: bounded_wdsm_design_reference.md,",
    "bounded_wdsm_explicit_radius_geometry.md and fitted_wdsm_pipeline_bridge.md.",
    "These are analytic declarations, never inferred from a fitted dataset."),
  assumptions = setNames(rep(TRUE, 7L), c("iid_bounded_sampling", "known_probability_weights",
    "target_identification_overlap", "correct_full_x_means", "current_score_geometry",
    "regular_joint_nuisance_influence", "positive_limit_variance")),
  conditional_refit_expansion = TRUE)

# Source ratio algebra; fixed original w and K, with all count-fitted predictions.
wm_cal_original_draw <- function(fit, multiplicity, mean0, mean1 = NULL) {
  n <- fit$n; Z <- fit$data$Z; Y <- fit$data$Y; w <- fit$analysis_weights
  if (length(multiplicity) != n || any(!is.finite(multiplicity)) ||
      any(multiplicity < 0) || any(multiplicity != floor(multiplicity)) ||
      sum(multiplicity) != n || any(vapply(0:1, function(z)
        sum(multiplicity[Z == z]) <= 0, logical(1L)))) stop("Invalid or missing-arm count column.")
  if (length(mean0) != n || any(!is.finite(mean0))) stop("Invalid control predictions.")
  K <- fit$loads$incoming[cbind(seq_len(n), Z + 1L)]
  if (fit$estimand == "PATE") {
    if (length(mean1) != n || any(!is.finite(mean1))) stop("Invalid treated predictions.")
    numerator <- sum(multiplicity*w*(mean1-mean0)) +
      sum(multiplicity*(2*Z-1)*(w+K)*(Y-ifelse(Z == 1, mean1, mean0)))
    denominator <- sum(multiplicity*w)
  } else {
    numerator <- sum(multiplicity*Z*w*(Y-mean0)) -
      sum(multiplicity*(1-Z)*K*(Y-mean0))
    denominator <- sum(multiplicity*Z*w)
  }
  value <- numerator/denominator
  if (!is.finite(value) || denominator <= 0) stop("Original ratio arithmetic failed.")
  value
}

.wm_cal_point <- function(branch, model, estimand, M, truth, n, fit = NULL,
                          inference = NULL, error = "") {
  point_ok <- !is.null(fit) && is.finite(fit$estimate)
  V <- if (!is.null(inference)) inference$root_n_variance else NA_real_
  if (is.null(inference) && point_ok && branch != "actual_fitted") {
    rows <- fit$contributions$actual
    V <- mean((rows-mean(rows))^2)
  }
  interval_ok <- point_ok && is.finite(V) && V > 0 &&
    (is.null(inference) || isTRUE(inference$sampling_inference_available))
  unavailable_reason <- if (interval_ok) "" else if (!point_ok) "Original point unavailable" else
    if (!is.null(inference) && !isTRUE(inference$sampling_inference_available))
      "Sampling inference contract unavailable" else if (!is.finite(V))
        "Root-n variance unavailable or nonfinite" else "Root-n variance nonpositive"
  estimate <- if (point_ok) fit$estimate else NA_real_
  se <- if (interval_ok) sqrt(V/n) else NA_real_
  V0 <- if (!is.null(inference)) inference$V0 else if (point_ok && branch != "actual_fitted") V else NA_real_
  naive_se <- if (is.finite(V0) && V0 > 0) sqrt(V0/n) else NA_real_
  z <- stats::qnorm(.975)
  data.frame(branch = branch, model = model, estimand = estimand, M = M, n = n,
    truth = truth, point_status = if (point_ok) "completed" else "failed",
    interval_status = if (interval_ok) "completed" else "unavailable",
    interval_unavailable_reason = unavailable_reason, error = error,
    estimate = estimate, root_n_variance = V, se = se,
    lower = estimate-z*se, upper = estimate+z*se,
    covered = if (interval_ok) abs(estimate-truth) <= z*se else NA,
    V0 = V0, cross_term = if (!is.null(inference)) inference$cross_term else NA_real_,
    nuisance_variance = if (!is.null(inference)) inference$nuisance_variance else NA_real_,
    naive_se = if (branch == "actual_fitted") naive_se else NA_real_,
    naive_lower = if (branch == "actual_fitted") estimate-z*naive_se else NA_real_,
    naive_upper = if (branch == "actual_fitted") estimate+z*naive_se else NA_real_,
    stringsAsFactors = FALSE)
}

.wm_cal_oracle <- function(data, estimand, M, scores0, scores1) {
  wdsmatch::wm_match(Y = data$Y, Z = data$Z, weights = data$weights,
    scores0 = scores0, scores1 = if (estimand == "PATE") scores1 else NULL,
    M = M, estimand = estimand, method = "self_normalized", mean0 = data$mu0,
    mean1 = if (estimand == "PATE") data$mu1 else NULL, variance = FALSE)
}

wm_cal_known_scores <- function(data) {
  k <- .wm_bd_constants("source_good"); r <- .25
  center <- k$a0+k$delta*k$a7
  half <- k$delta*r*sqrt(1+k$a7^2)
  ps <- (data$population_propensity-stats::plogis(center))/
    (stats::plogis(center+half)-stats::plogis(center-half))
  list(scores0 = cbind(ps, (data$mu0-.75)/(2*r*sqrt(.3^2+.75^2))),
       scores1 = cbind(ps, (data$mu1-2.05)/(2*r*sqrt(.5^2+1.05^2))))
}

# Inner count summaries preserve scale and the all-prescribed-draw gate.
wm_cal_count_summary <- function(draws, original_estimate, analytic_V, n) {
  B <- nrow(draws)
  complete <- B > 1L && all(draws$status == "completed") &&
    all(is.finite(draws$estimate))
  out <- list(B = B, successful_draws = sum(draws$status == "completed"),
    interval_status = "unavailable", interval_unavailable_reason =
      if (B <= 1L) "Insufficient prescribed count draws" else if (!complete)
        "At least one prescribed draw failed or is nonfinite" else if (!is.finite(original_estimate))
          "Original point unavailable or nonfinite" else "Source count variance not yet evaluated",
    source_point_variance_B = NA_real_,
    root_variance_unbiased_success_law = NA_real_, original_refit_se = NA_real_,
    lower = NA_real_, upper = NA_real_, one_step_exact_root_variance = analytic_V,
    one_step_root_variance_B = NA_real_, one_step_root_variance_Bminus1 = NA_real_,
    root_difference_mean = NA_real_, root_difference_rms = NA_real_,
    root_difference_mean_mcse = NA_real_, root_difference_rms_mcse = NA_real_,
    refit_root_variance_mcse = NA_real_, one_step_root_variance_mcse = NA_real_,
    paired_root_variance_difference = NA_real_, paired_variance_difference_mcse = NA_real_,
    scale_note = "Refit source variance is point scale/divisor B; one-step roots are already root-n.",
    failure_law_note = paste("B/(B-1) targets successful-draw conditional law on all-B-success admission;",
      "no unconditional failed-draw completion or successful-subset replacement is made."))
  # One-step remains available even when a count nuisance fit fails.
  if (B > 1L && all(is.finite(draws$one_step_root))) {
    one <- draws$one_step_root
    out$one_step_root_variance_B <- mean((one-mean(one))^2)
    out$one_step_root_variance_Bminus1 <- stats::var(one)
    out$one_step_root_variance_mcse <- (B/(B-1))*stats::sd((one-mean(one))^2)/sqrt(B)
  }
  if (!complete || !is.finite(original_estimate)) return(out)
  root <- sqrt(n)*(draws$estimate-original_estimate)
  q <- (root-mean(root))^2
  vpoint <- mean((draws$estimate-mean(draws$estimate))^2)
  out$source_point_variance_B <- vpoint
  out$root_variance_unbiased_success_law <- n*B/(B-1)*vpoint
  out$original_refit_se <- sqrt(vpoint)
  if (is.finite(vpoint) && vpoint > 0) {
    out$interval_status <- "completed"
    out$interval_unavailable_reason <- ""
    out$lower <- original_estimate-stats::qnorm(.975)*sqrt(vpoint)
    out$upper <- original_estimate+stats::qnorm(.975)*sqrt(vpoint)
  } else out$interval_unavailable_reason <- if (!is.finite(vpoint))
    "Source count variance nonfinite" else "Source count variance nonpositive"
  out$refit_root_variance_mcse <- (B/(B-1))*stats::sd(q)/sqrt(B)
  if (!all(is.finite(draws$one_step_root))) return(out)
  one <- draws$one_step_root; qone <- (one-mean(one))^2; delta <- root-one
  out$root_difference_mean <- mean(delta)
  out$root_difference_rms <- sqrt(mean(delta^2))
  out$root_difference_mean_mcse <- stats::sd(delta)/sqrt(B)
  out$root_difference_rms_mcse <- if (out$root_difference_rms > 0)
    stats::sd(delta^2)/(2*sqrt(B)*out$root_difference_rms) else 0
  out$one_step_root_variance_mcse <- (B/(B-1))*stats::sd(qone)/sqrt(B)
  out$paired_root_variance_difference <- stats::var(root)-stats::var(one)
  out$paired_variance_difference_mcse <- (B/(B-1))*stats::sd(q-qone)/sqrt(B)
  out
}

# All saves are supplied by the exclusive output writer. Count stacks are
# transient: save compact attempts and predictions' resulting draws, not IF arrays.
wm_cal_analyze_dataset <- function(data, counts, save_object, save_table,
                                    event = function(...) invisible(NULL)) {
  n <- length(data$Y); B <- ncol(counts)
  if (!identical(data$index, "source_good") || data$radius != .25 ||
      !is.matrix(counts) || nrow(counts) != n || any(!is.finite(counts)) ||
      any(counts < 0) || any(counts != floor(counts)) ||
      (B > 0L && any(colSums(counts) != n))) stop("Invalid predeclared dataset/counts.")
  point_rows <- list(); add <- function(x) { point_rows[[length(point_rows)+1L]] <<- x }
  # Independent controls execute even if every fitted model fails.
  known <- wm_cal_known_scores(data)
  for (estimand in c("PATE", "PATT")) for (M in c(1L,3L)) {
    id <- paste("known", estimand, M, sep = "-"); event("known_oracle", id)
    fit <- tryCatch(.wm_cal_oracle(data, estimand, M, known$scores0, known$scores1), error = identity)
    save_object(paste0(id,".rds"), fit)
    row <- .wm_cal_point("known_map_oracle", "shared", estimand, M, data$truth[[estimand]], n,
      fit = if (!inherits(fit,"error")) fit else NULL,
      error = if (inherits(fit,"error")) conditionMessage(fit) else "")
    save_table(paste0(id,"-point.csv"), row); add(row)
  }
  for (model in c("CorCor", "MisCor")) for (estimand in c("PATE", "PATT")) {
    id <- paste(model,estimand,sep="-"); event("base_fit", id)
    base <- wm_cal_fit(data, model, estimand, save_attempt = function(i,x)
      save_object(paste0(id,"-base-attempt",i,".rds"),x))
    if (!is.null(base$stack)) save_object(paste0(id,"-base-stack.rds"),base$stack)
    objects <- list(); point_fits <- list(); one_steps <- list(); draw_tables <- list()
    for (M in c(1L,3L)) {
      tag <- paste(id,M,sep="-"); event("fitted_graph",tag)
      obj <- if (is.null(base$stack)) simpleError(base$error) else tryCatch(
        wm_wdsm_fitted_pipeline_reference(base$stack,M,.wm_cal_contract()),error=identity)
      objects[[as.character(M)]] <- if (!inherits(obj,"error")) obj else NULL
      actual_fit <- if (!inherits(obj,"error")) obj$fit else NULL
      # Preserve a valid point if downstream gradient/variance validation fails.
      if (is.null(actual_fit) && !is.null(base$stack)) actual_fit <- tryCatch(
        wdsmatch::wm_match(Y=data$Y,Z=data$Z,weights=data$weights,
          scores0=base$stack$scores0,scores1=base$stack$scores1,M=M,estimand=estimand,
          method="self_normalized",mean0=base$stack$mean0,mean1=base$stack$mean1,
          variance=FALSE), error=identity)
      if (inherits(actual_fit,"error")) {
        save_object(paste0(tag,"-point-fallback-error.rds"),actual_fit)
        actual_fit <- NULL
      }
      point_fits[[as.character(M)]] <- actual_fit
      stored <- obj
      if (!inherits(stored,"error")) stored$stack <- NULL # separately pinned base stack
      save_object(paste0(tag,"-actual.rds"),list(result=stored,point_fit=actual_fit))
      row <- .wm_cal_point("actual_fitted",model,estimand,M,data$truth[[estimand]],n,
        fit=actual_fit,inference=if (!inherits(obj,"error")) obj$inference else NULL,
        error=if (inherits(obj,"error")) conditionMessage(obj) else "")
      save_table(paste0(tag,"-actual-point.csv"),row); add(row)
      oracle <- if (is.null(base$stack)) simpleError(base$error) else tryCatch(
        .wm_cal_oracle(data,estimand,M,base$stack$scores0,base$stack$scores1),error=identity)
      if (!inherits(oracle,"error") && !is.null(actual_fit) &&
          !identical(oracle$graph$neighbors,actual_fit$graph$neighbors)) stop("Oracle fitted graph changed.")
      save_object(paste0(tag,"-oracle-fitted.rds"),oracle)
      row <- .wm_cal_point("oracle_means_fitted_map",model,estimand,M,data$truth[[estimand]],n,
        fit=if (!inherits(oracle,"error")) oracle else NULL,
        error=if (inherits(oracle,"error")) conditionMessage(oracle) else "")
      save_table(paste0(tag,"-oracle-fitted-point.csv"),row); add(row)
      one <- if (B > 0L && !inherits(obj,"error")) tryCatch(
        wm_wdsm_fitted_count_reference(obj,counts),error=identity) else NULL
      if (inherits(one,"error")) {
        save_object(paste0(tag,"-one-step-error.rds"),one)
        one <- NULL
      }
      one_steps[[as.character(M)]] <- one
      if (!is.null(one)) save_object(paste0(tag,"-one-step.rds"),one)
      draw_tables[[as.character(M)]] <- data.frame(draw=seq_len(B),status=rep("pending",B),
        error=rep("",B),estimate=rep(NA_real_,B),root_n=rep(NA_real_,B),
        one_step_root=if (!is.null(one)) one$root_n_draws else rep(NA_real_,B),
        refit_attempts=integer(B),refined=rep(FALSE,B),elapsed_seconds=rep(NA_real_,B))
    }
    if (B > 0L) for (b in seq_len(B)) {
      event("count_fit",id,b)
      refit <- if (is.null(base$stack)) list(status="failed",error=paste("Base fit failed:",base$error),
        stack=NULL,attempts=list(),refined=FALSE) else wm_cal_fit(data,model,estimand,counts[,b],
          fitter=wm_cal_prediction_stack, save_attempt=function(i,x) save_object(sprintf("%s-count%03d-attempt%d.rds",id,b,i),x))
      for (M in c(1L,3L)) {
        key <- as.character(M); point_fit <- point_fits[[key]]; table <- draw_tables[[key]]
        value <- if (is.null(point_fit)) simpleError("Original fitted point/graph unavailable") else
          if (is.null(refit$stack)) simpleError(refit$error) else tryCatch(
            wm_cal_original_draw(point_fit,counts[,b],refit$stack$mean0,refit$stack$mean1),error=identity)
        table$status[b] <- if (inherits(value,"error")) "failed" else "completed"
        table$error[b] <- if (inherits(value,"error")) conditionMessage(value) else ""
        table$estimate[b] <- if (!inherits(value,"error")) value else NA_real_
        table$root_n[b] <- if (!inherits(value,"error")) sqrt(n)*(value-point_fit$estimate) else NA_real_
        table$refit_attempts[b] <- length(refit$attempts); table$refined[b] <- refit$refined
        table$elapsed_seconds[b] <- sum(vapply(refit$attempts,function(x)x$elapsed_seconds,numeric(1L)))
        save_table(sprintf("%s-%d-count%03d.csv",id,M,b),table[b,,drop=FALSE])
        draw_tables[[key]] <- table
      }
      rm(refit) # no count-specific sandwich is serialized
    }
    for (M in c(1L,3L)) {
      tag <- paste(id,M,sep="-"); obj <- objects[[as.character(M)]]
      point_fit <- point_fits[[as.character(M)]]
      table <- draw_tables[[as.character(M)]]
      save_table(paste0(tag,"-draws.csv"),table)
      if (B > 0L) save_object(paste0(tag,"-count-summary.rds"),wm_cal_count_summary(table,
        if (!is.null(point_fit)) point_fit$estimate else NA_real_,
        if (!is.null(obj)) obj$inference$root_n_variance else NA_real_,n))
    }
    rm(base,objects,point_fits,one_steps)
  }
  records <- do.call(rbind,point_rows)
  stopifnot(nrow(records)==20L)
  save_table("point_records.csv",records)
  event("dataset_completed","",NA_integer_)
  invisible(records)
}
