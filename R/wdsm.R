wdsm_validate_formula <- function(formula, response, X) {
  if (is.null(formula)) {
    terms <- if (ncol(X)) paste(sprintf("`%s`", names(X)), collapse = " + ") else "1"
    formula <- stats::as.formula(paste(response, "~", terms))
  }
  if (!inherits(formula, "formula") || length(formula) != 3L ||
      !identical(formula[[2L]], as.name(response)))
    stop("Supply a model formula with response ", response, call. = FALSE)
  # Expand '.' using X only. In the fitting frame both Y and Z exist, so
  # expanding there would leak the outcome into a propensity model.
  formula_data <- cbind(data.frame(rep(0, nrow(X))), X)
  names(formula_data)[1L] <- response
  formula <- stats::formula(stats::terms(formula, data = formula_data))
  if (!all(all.vars(formula[[3L]]) %in% names(X)))
    stop("Model predictors must be columns of X", call. = FALSE)
  formula
}

wdsm_run <- function(Y, X, Z, weights, M, ps, pg, model.ps, model.pg,
                     sampling, use.bias.correction, varest, boots, alpha,
                     estimand, call, bootstrap_counts = NULL,
                     tie_seed = 20260917L, tie_tolerance = 64 * .Machine$double.eps) {
  scalar_integer <- function(x, minimum) is.numeric(x) && length(x) == 1L &&
    is.finite(x) && x >= minimum && x <= .Machine$integer.max && x == floor(x)
  scalar_logical <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)
  if (is.null(weights)) stop("'weights' are required for WDSM.", call. = FALSE)
  if (!is.numeric(Y) || !is.null(dim(Y)) || any(!is.finite(Y)))
    stop("'Y' must be a finite numeric vector.", call. = FALSE)
  n <- length(Y)
  if (!(is.numeric(Z) || is.logical(Z)) || !is.null(dim(Z)) ||
      length(Z) != n || anyNA(Z) || !all(Z %in% c(0, 1)))
    stop("'Z' must be a binary (0/1) vector with the same length as 'Y'.", call. = FALSE)
  Z <- as.numeric(Z)
  if (!is.numeric(weights) || !is.null(dim(weights)) || length(weights) != n ||
      any(!is.finite(weights) | weights <= 0) || !is.finite(sum(weights)))
    stop("'weights' must be finite, positive, and have the same length as 'Y'.", call. = FALSE)
  if (!(is.data.frame(X) || is.matrix(X)))
    stop("'X' must be a numeric matrix or data frame.", call. = FALSE)
  X <- as.data.frame(X)
  if (nrow(X) != n || !all(vapply(X, is.numeric, logical(1L))) ||
      any(!is.finite(as.matrix(X))))
    stop("'X' must have one finite numeric row per outcome.", call. = FALSE)
  if (anyDuplicated(names(X)) || any(!nzchar(names(X))) ||
      any(names(X) %in% c("Y", "Z")) || any(grepl("`", names(X), fixed = TRUE)))
    stop("X column names must be unique, nonempty, and exclude Y, Z, and backticks.", call. = FALSE)
  if (!scalar_integer(M, 1L)) stop("'M' must be a positive integer.", call. = FALSE)
  if (!scalar_logical(varest) || !scalar_logical(use.bias.correction))
    stop("'varest' and 'use.bias.correction' must be TRUE or FALSE.", call. = FALSE)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1)
    stop("'alpha' must be strictly between zero and one.", call. = FALSE)
  if (varest && !scalar_integer(boots, 2L))
    stop("'boots' must be an integer of at least two when varest = TRUE.", call. = FALSE)
  if (n < 2L || sum(Z == 1) < 1L || sum(Z == 0) < M ||
      (estimand == "PATE" && sum(Z == 1) < M))
    stop("Insufficient treatment-arm observations or opposite-arm donors for M.", call. = FALSE)
  if (!is.null(ps) && (!is.numeric(ps) || !is.null(dim(ps)) || length(ps) != n ||
      any(!is.finite(ps) | ps <= 0 | ps >= 1)))
    stop("'ps' must contain one finite probability strictly between zero and one per outcome.", call. = FALSE)
  if (!is.null(pg)) {
    if (is.data.frame(pg)) pg <- as.matrix(pg)
    allowed <- if (estimand == "PATE") 2L else 1:2
    if (!is.matrix(pg) || !is.numeric(pg) || nrow(pg) != n ||
        !ncol(pg) %in% allowed || any(!is.finite(pg)))
      stop("'pg' must be a finite numeric matrix: two columns for PATE, one or two for PATT.", call. = FALSE)
  }
  # An explicitly supplied score takes precedence; no unrelated model is fit.
  if (is.null(ps)) model.ps <- wdsm_validate_formula(model.ps, "Z", X)
  if (is.null(pg)) model.pg <- wdsm_validate_formula(model.pg, "Y", X)
  M <- as.integer(M)
  if (varest) boots <- as.integer(boots)
  scores <- estimate_scores(Y, X, Z, weights, ps, pg, model.ps, model.pg,
                            sampling, estimand, use.bias.correction = use.bias.correction)
  pt <- wdsm_point(Y, Z, weights, scores, M, estimand, tie_seed, tie_tolerance)
  pt$X_internal <- X
  result <- list(estimate = pt$estimate, se = NA_real_, ci = c(NA_real_, NA_real_),
    boot.estimates = NULL, M = M, n = n, n.treated = sum(Z == 1),
    n.control = sum(Z == 0), call = call, variance = NA_real_, alpha = alpha,
    interval.type = if (varest) "Wald" else NA_character_,
    variance.divisor = if (varest) "B" else NA_character_, n.boot = 0L,
    estimand = estimand, sampling = sampling,
    settings = list(bias.correction = use.bias.correction,
      propensity.coordinate = "probability",
      score.standardization = "pooled unweighted; multinomial multiplicity in replicates",
      sieve = if (use.bias.correction) "complete quadratic in each arm's own double score" else "none",
      supplied.ps.fixed = !is.null(ps), supplied.pg.fixed = !is.null(pg),
      replication = "original matches and raw weighted reuse held fixed"),
    diagnostics = list(ties = pt$tie_diagnostics, identity.error = pt$identity_error,
                       point.ps = scores$ps_diagnostics))
  if (varest) {
    boot <- wdsm_bootstrap(pt, boots, alpha, model.ps, model.pg, sampling,
                           ps, pg, use.bias.correction, bootstrap_counts)
    result$variance <- boot$variance
    result$se <- boot$se
    result$ci <- boot$ci
    result$boot.estimates <- boot$boot_estimates
    result$n.boot <- boots
    result$diagnostics$bootstrap.ps <- boot$ps_diagnostics
  }
  class(result) <- "wdsmatch"
  result
}
