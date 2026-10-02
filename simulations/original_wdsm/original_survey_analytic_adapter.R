# Empirical full-X analytic rows for the original WDSM comparison.
# This does not certify full-X correctness, iid sampling, or a bootstrap theorem.
# The supplied weighted batch is the authoritative point-estimator calculation.
ows_analytic_batch <- function(data, weighted, wm, alpha = .05) {
  stopifnot(identical(weighted$n, nrow(data)), weighted$estimand %in% c("PATE", "PATT"),
    weighted$sampling_design %in% c("retrospective", "prospective"),
    weighted$specification %in% c("CorCor", "CorMis", "MisCor", "MisMis"),
    length(alpha) == 1L, is.finite(alpha), alpha > 0, alpha < 1)
  specification <- weighted$specification; estimand <- weighted$estimand
  ps_terms <- c(paste0("X", 1:6), if (specification %in% c("CorCor", "CorMis")) "X1:X2")
  pg_terms <- c(paste0("X", 1:6), if (specification %in% c("CorCor", "MisCor")) "X1:X2")
  ps <- stats::model.matrix(stats::reformulate(ps_terms), data)
  pg <- stats::model.matrix(stats::reformulate(pg_terms), data)
  # Derivatives and influence rows are computed once, then reused across M.
  solved <- ows_capture(function() {
    reference <- weighted$bases$DSM
    if (!isTRUE(reference$ok)) stop("Recorded DSM prediction stack unavailable",call.=FALSE)
    input <- reference$value$stack$inputs
    ps_weighting <- if (weighted$sampling_design=="retrospective") "probability" else "unit"
    if (!identical(input$Y,as.numeric(data$Y)) || !identical(input$Z,as.numeric(data$A)) ||
        !identical(input$weights,as.numeric(data$survey_weight)) ||
        !identical(input$ps_design,ps) || !identical(input$pg0_design,pg) ||
        !identical(input$pg1_design,if(estimand=="PATE")pg else NULL) ||
        !identical(reference$value$stack$ps_weighting,ps_weighting))
      stop("Analytic data/design does not match recorded DSM prediction inputs",call.=FALSE)
    wm$.wm_wdsm_fit_stack(list(Y=data$Y, Z=data$A,
      weights=data$survey_weight, ps_design=ps, pg0_design=pg,
      pg1_design=if (estimand=="PATE") pg else NULL, estimand=estimand,
      ps_weighting=ps_weighting),
      wm$.wm_wdsm_nuisance_stack, wm$.wm_wdsm_controls(list()))
  })
  rows <- list(); assemblies <- list()
  for (M in c(1L, 3L, 5L)) {
    key <- paste0("WM_DSM_analytic_M", M)
    point <- weighted$points[[paste0("WM_DSM_M", M)]]
    stopifnot(is.list(point), is.logical(point$ok), length(point$ok)==1L, !is.na(point$ok))
    result <- ows_capture(function() {
      if (!point$ok) stop("Matching point unavailable", call.=FALSE)
      if (!solved$ok) stop(solved$error, call.=FALSE)
      wm$.wm_wdsm_fitted_pipeline(solved$value$stack, M=M, conf.level=1-alpha)
    })
    assemblies[[key]] <- result
    aligned <- result$ok &&
      identical(result$value$fit$graph$neighbors, point$value$graph$neighbors) &&
      abs(result$value$estimate-point$value$estimate) <= 1e-10*(1+abs(point$value$estimate))
    inf <- if (result$ok) result$value$inference else NULL
    finite_interval <- !is.null(inf) && isTRUE(inf$available) &&
      length(inf$variance)==1L && is.finite(inf$variance) && inf$variance>0 &&
      length(inf$conf.int)==2L && all(is.finite(inf$conf.int))
    okay <- aligned && finite_interval
    reason <- if (!result$ok) result$error else if (!aligned)
      "Analytic assembly point or donor graph differs from recorded WM DSM point" else if (!finite_interval)
      "Analytic variance or interval unavailable/nonfinite/nonpositive" else ""
    rows[[key]] <- list(estimate=if (point$ok) point$value$estimate else NA_real_,
      point_status=if (point$ok) "ok" else "failed",
      interval_status=if (okay) "ok" else "failed", interval_error=reason,
      variance=if (okay) inf$variance else NA_real_,
      conf.int=if (okay) inf$conf.int else c(lower=NA_real_,upper=NA_real_),
      variance_method="fitted_stack_full_X_analytic_approximation",
      point_source=paste0("WM_DSM_M",M), aligned_with_recorded_point=aligned,
      assumptions_verified=FALSE, original_refit_limit_agreement_declared=FALSE,
      scope=paste("Separate empirical analytic approximation; full-X correction and iid premises",
        "are unverified under these source designs/specifications. No cluster/stratum correction.",
        "Supplied probability weights held fixed."))
  }
  list(rows=rows, assemblies=assemblies, nuisance=solved,
    review_status="DRAFT_REQUIRES_INDEPENDENT_INTEGRATION_REVIEW")
}
