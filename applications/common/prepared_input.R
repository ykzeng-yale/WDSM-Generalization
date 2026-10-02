# SPDX-License-Identifier: GPL-3.0-only
# Explicit data-only adapters. Sourcing performs no data import or analysis.
wm_application_input <- function(data, id, Y, Z, W, X, design, variables,
                                 study, outcome, baseline, units, ps_weighting,
                                 provenance = list()) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required")
  n <- length(Y)
  vector_ok <- function(v) is.numeric(v) && !is.complex(v) && is.null(dim(v)) &&
    length(v) == n && all(is.finite(v))
  if (n < 2L || !vector_ok(Y) || !vector_ok(Z) || !vector_ok(W) ||
      !all(Z %in% 0:1) || length(unique(Z)) != 2L || any(W <= 0))
    stop("Y/Z/W must be aligned finite vectors, both binary arms, and positive supplied W")
  if (!(is.character(id) || is.numeric(id)) || !is.null(dim(id)) ||
      length(id) != n || anyNA(id) || anyDuplicated(id) ||
      (is.character(id) && any(!nzchar(id))) ||
      (is.numeric(id) && any(!is.finite(id)))) stop("Explicit unique source-order IDs required")
  matrix_ok <- function(v) is.matrix(v) && is.numeric(v) && !is.complex(v) &&
    nrow(v) == n && ncol(v) > 0L && all(is.finite(v)) &&
    !is.null(colnames(v)) && !anyDuplicated(colnames(v)) && all(nzchar(colnames(v)))
  if (!matrix_ok(X) || !matrix_ok(design)) stop("Finite named X and separate nuisance design required")
  if (!is.data.frame(data) || nrow(data) != n || !is.character(variables) ||
      !length(variables) || anyNA(variables) || anyDuplicated(variables) ||
      any(!nzchar(variables)) || !all(c("A", variables) %in% names(data)) ||
      !identical(as.numeric(data$A), as.numeric(Z)) ||
      !all(vapply(data[variables], function(v) vector_ok(v), logical(1))))
    stop("Explicit prepared PS formula data and variables must agree with Z and row order")
  labels <- list(study, outcome, baseline, units)
  if (!all(vapply(labels, function(v) is.character(v) && length(v) == 1L &&
      !is.na(v) && nzchar(v), logical(1))) ||
      !is.character(ps_weighting) || length(ps_weighting) != 1L ||
      !ps_weighting %in% c("unit", "probability")) stop("Explicit labels/PS weighting required")
  list(data = data, X = X, design = design, variables = variables, Y = Y,
    Z = Z, W = W, n = n, study = study, outcome = outcome, baseline = baseline,
    units = units, ps_weighting = ps_weighting, id = id, provenance = provenance,
    sample_sha256 = digest::digest(list(id = id, Y = Y, Z = Z, W = W, X = X), algo = "sha256"),
    dataset_sha256 = digest::digest(list(id = id, Z = Z, W = W, X = X), algo = "sha256"),
    nuisance_sha256 = digest::digest(list(
      schema = "wm_application_nuisance_input_v1",
      design_columns = enc2utf8(colnames(design)),
      design_values = matrix(as.double(design), nrow = n, ncol = ncol(design)),
      formula_response = "A", formula_intercept = TRUE, formula_variables = enc2utf8(variables),
      formula_data = setNames(lapply(data[c("A", variables)], function(v) as.double(v)),
        enc2utf8(c("A", variables))),
      ps_weighting = ps_weighting), algo = "sha256"))
}

wm_application_nhanes_input <- function(data, baseline, provenance = list()) {
  if (!is.character(baseline) || length(baseline) != 1L || is.na(baseline) ||
      !baseline %in% c("manuscript_primary", "canonical_source_3790"))
    stop("Declare the separately prepared NHANES primary or canonical sample")
  vars <- c("ridageyr", "riagendr", "ridreth1", "dmdeduc2", "indfmpir", "bmxbmi", "smq020")
  required <- c("seqn", "A", "Y", "survey_weight", vars)
  if (!is.data.frame(data) || !all(required %in% names(data))) stop("Missing prepared NHANES columns")
  expected <- if (baseline == "manuscript_primary") c(3783L, 1001L) else c(3790L, 1002L)
  if (nrow(data) != expected[1L] || sum(data$A) != expected[2L]) stop("Wrong declared NHANES cohort")
  X <- as.matrix(data[vars]); storage.mode(X) <- "double"
  wm_application_input(data, data$seqn, data$Y, data$A, data$survey_weight,
    X, cbind(intercept = 1, X), vars, "NHANES", "HbA1c", baseline,
    "HbA1c percentage points", "probability", provenance)
}

wm_application_nsduh_input <- function(data, outcome, provenance = list()) {
  if (!is.character(outcome) || length(outcome) != 1L || is.na(outcome) ||
      !outcome %in% c("AMIPY", "SMIPY")) stop("Declare AMIPY or SMIPY")
  vars <- c("male", "white", "black", "hispanic", "income_low", "income_medium",
    "income_high", "large_metro", "small_metro", "early_alcohol", "ever_cigarette",
    "age_21_23", "age_24_25", "age_26_29")
  required <- c("QUESTID2", "A", "ANALWT2_C", "AMIPY", "SMIPY", vars)
  if (!is.data.frame(data) || !all(required %in% names(data)) || nrow(data) != 17471L)
    stop("Explicit shared two-outcome prepared NSDUH17471 data required")
  if (sum(data$A) != 1877L || !all(data$AMIPY %in% 0:1) ||
      !all(data$SMIPY %in% 0:1) || !all(as.matrix(data[vars]) %in% 0:1))
    stop("Prepared NSDUH arm, outcome or declared binary-coordinate contract differs")
  data$Y <- data[[outcome]]; data$survey_weight <- data$ANALWT2_C
  X <- as.matrix(data[vars]); storage.mode(X) <- "double"
  wm_application_input(data, data$QUESTID2, data$Y, data$A, data$survey_weight,
    X, cbind(intercept = 1, X), vars, "NSDUH", outcome, "manuscript_primary",
    "risk difference", "unit", provenance)
}

wm_application_ecls_input <- function(selected, provenance = list()) {
  if (!is.list(selected) ||
      !identical(selected$case_id, "ECLS_P2SPECND_Kyear_C6R4MSCL_corrected_observed_v1") ||
      length(selected$diagnostics$defects) != 0L || !is.character(selected$CHILDID) ||
      !is.character(selected$S2_ID) || nrow(selected$X) != 7220L || sum(selected$Z) != 725L ||
      ncol(selected$X) != 38L || ncol(selected$nuisance_design) != 37L ||
      !identical(selected$diagnostics$reference_omissions, c("P1ELHS", "TWOPAR")) ||
      !isTRUE(selected$diagnostics$expectation_identity) || !isTRUE(selected$diagnostics$family_identity) ||
      !isTRUE(selected$diagnostics$identical_span) || selected$diagnostics$rank_nuisance_global != 37L ||
      !all(selected$diagnostics$ranks_nuisance_arms == 37L)) stop("Corrected ECLS7220 construction contract differs")
  X <- selected$X; D <- selected$nuisance_design; nvars <- selected$nuisance_features
  if (!identical(colnames(D), c("intercept", nvars)) || !identical(colnames(X), selected$diagnostics$feature_order))
    stop("ECLS38 matching order and separate nuisance37 order differ")
  d <- as.data.frame(X, check.names = FALSE)
  rownames(d) <- as.character(selected$source_row_index)
  d$CHILDID <- selected$CHILDID; d$S2_ID <- selected$S2_ID
  d$TREAT <- as.integer(selected$Z); d$C6R4MSCL <- selected$Y; d$C1_6FC0 <- selected$W
  d$w <- selected$W; d$A <- d$TREAT; d$Y <- d$C6R4MSCL; d$survey_weight <- d$C1_6FC0
  wm_application_input(d, selected$CHILDID, selected$Y, d$A, selected$W,
    X, D, nvars, "ECLS-K", "C6R4MSCL", "corrected_selected_observed_7220",
    "fifth-grade math IRT scale score points", "probability", provenance)
}
