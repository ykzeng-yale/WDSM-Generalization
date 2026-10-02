# SPDX-License-Identifier: GPL-3.0-only
# Data-only portable wrapper of the corrected ECLS construction.
# Algorithm blocks are retained from the accepted constructor; no work at source time.
wm_ecls_construct_selected <- function(decoded, missing_reasons, companion, schema,
                                       provenance = list()) {
  d <- decoded; reasons <- missing_reasons; comp <- companion; pins <- provenance
  CASE_ID <- "ECLS_P2SPECND_Kyear_C6R4MSCL_corrected_observed_v1"
  if (!is.list(schema) || !identical(schema$Case_ID, CASE_ID)) stop("Corrected ECLS scientific schema required")
  if (!is.data.frame(d) || ncol(d) != 39L || nrow(d) != 21409L ||
      !all(c("CHILDID", "S1_ID", "S2_ID") %in% names(d)) ||
      !all(vapply(d[c("CHILDID", "S1_ID", "S2_ID")], is.character, logical(1))) ||
      !all(vapply(d[setdiff(names(d), c("CHILDID", "S1_ID", "S2_ID"))], is.numeric, logical(1))) ||
      anyNA(d$CHILDID) || any(!nzchar(d$CHILDID)) || anyDuplicated(d$CHILDID) ||
      !is.data.frame(reasons) || !identical(names(d), names(reasons)) ||
      !is.data.frame(comp) || !identical(d$CHILDID, reasons$CHILDID) ||
      !identical(d$CHILDID, comp$CHILDID)) stop("Decoded39 types and all-row ID/order binding differ")
FEATURES <- c("FEMALE", "WHITE", "WKSESL", "C1R4RSCL", "C1R4MSCL", "S2KPUPRI",
  "P1ELHS", "P1EHS", "P1ESC", "P1EC", "P1EMS", "P1EPHD", "P1FIRKDG",
  "P1AGEENT", "T1LEARN", "P1HSEVER", "FKCHGSCH", "S2KMINOR", "P1FSTAMP",
  "SGLPAR", "TWOPAR", "P1NUMSIB", "P1HMAFB", "WKCAREPK", "PRETERM_GT2WK",
  "BIRTH_TOTAL_OZ", "C1FMOTOR", "C1GMOTOR", "P1HSCALE", "P1SADLON",
  "P1IMPULS", "P1ATTENI", "P1SOLVE", "P1PRONOU", "DISABILITY_REVISED",
  "AVG_C1READ_CONTEXT", "AVG_C1MATH_CONTEXT", "AVG_SES_CONTEXT")
RAW_REQUIRED <- c("GENDER", "WKWHITE", "WKSESL", "C1R4RSCL", "C1R4MSCL",
  "S2KPUPRI", "P1EXPECT", "P1FIRKDG", "P1AGEENT", "T1LEARN", "P1HSEVER",
  "FKCHGSCH", "S2KMINOR", "P1FSTAMP", "P1HFAMIL", "P1NUMSIB", "P1HMAFB",
  "WKCAREPK", "P1PREMAT", "P1WEIGHP", "P1WEIGHO", "C1FMOTOR", "C1GMOTOR",
  "P1HSCALE", "P1SADLON", "P1IMPULS", "P1ATTENI", "P1SOLVE", "P1PRONOU")
BINARY <- c("FEMALE", "WHITE", "S2KPUPRI", "P1ELHS", "P1EHS", "P1ESC",
  "P1EC", "P1EMS", "P1EPHD", "P1FIRKDG", "P1HSEVER", "FKCHGSCH",
  "P1FSTAMP", "SGLPAR", "TWOPAR", "WKCAREPK", "PRETERM_GT2WK", "DISABILITY_REVISED")
ORDERED <- c("S2KMINOR", "P1HSCALE", "P1ATTENI", "P1SOLVE", "P1PRONOU")

qr_rank <- function(M) {
  if (!nrow(M) || !ncol(M)) return(0L)
  if (any(!is.finite(M))) return(NA_integer_)
  qr(M, tol = 1e-10, LAPACK = FALSE)$rank
}
dependent_columns <- function(M) {
  if (!nrow(M)) return(colnames(M))
  q <- qr(M, tol = 1e-10, LAPACK = FALSE)
  if (q$rank == ncol(M)) character() else colnames(M)[q$pivot[seq.int(q$rank + 1L, ncol(M))]]
}
    n_all <- nrow(d)
    school_ids <- unique(d$S2_ID[!is.na(d$S2_ID) & nzchar(d$S2_ID)])
    context <- data.frame(S2_ID = school_ids)
    context_fields <- c(AVG_C1READ_CONTEXT = "C1R4RSCL", AVG_C1MATH_CONTEXT = "C1R4MSCL", AVG_SES_CONTEXT = "WKSESL")
    for (feature in names(context_fields)) {
      value <- d[[context_fields[[feature]]]]
      valid <- !is.na(d$S2_ID) & nzchar(d$S2_ID) & is.finite(value)
      counts <- tabulate(match(d$S2_ID[valid], school_ids), nbins = length(school_ids))
      sums <- rep(NA_real_, length(school_ids))
      if (any(valid)) {
        grouped <- rowsum(matrix(value[valid], ncol = 1L), d$S2_ID[valid], reorder = FALSE)
        sums[match(rownames(grouped), school_ids)] <- grouped[, 1L]
      }
      context[[feature]] <- ifelse(counts > 0L, sums / counts, NA_real_)
      context[[paste0(feature, "_N")]] <- counts
    }
    X <- matrix(NA_real_, nrow = n_all, ncol = length(FEATURES), dimnames = list(NULL, FEATURES))
    binary <- function(value, yes, valid) ifelse(is.finite(value) & value %in% valid, as.numeric(value %in% yes), NA_real_)
    X[, "FEMALE"] <- binary(d$GENDER, 2, 1:2)
    X[, "WHITE"] <- binary(d$WKWHITE, 1, 1:2)
    for (name in c("S2KPUPRI", "P1FIRKDG", "P1HSEVER", "P1FSTAMP", "WKCAREPK")) X[, name] <- binary(d[[name]], 1, 1:2)
    for (i in seq_len(6L)) X[, c("P1ELHS", "P1EHS", "P1ESC", "P1EC", "P1EMS", "P1EPHD")[i]] <- binary(d$P1EXPECT, i, 1:6)
    X[, "SGLPAR"] <- binary(d$P1HFAMIL, 3:4, 1:5)
    X[, "TWOPAR"] <- binary(d$P1HFAMIL, 1:2, 1:5)
    X[, "PRETERM_GT2WK"] <- binary(d$P1PREMAT, 1, 1:2)
    X[, "DISABILITY_REVISED"] <- binary(comp$RP1DISAB_official_revised, 1, 1:2)
    X[, "BIRTH_TOTAL_OZ"] <- comp$birth_weight_total_ounces
    direct <- setdiff(FEATURES, c(BINARY, "BIRTH_TOTAL_OZ", names(context_fields)))
    for (name in direct) X[, name] <- d[[name]]
    X[, "FKCHGSCH"] <- d$FKCHGSCH
    group_index <- match(d$S2_ID, context$S2_ID)
    for (name in names(context_fields)) X[, name] <- context[[name]][group_index]
    Z <- binary(d$P2SPECND, 1, 1:2); Y <- d$C6R4MSCL; W <- d$C1_6FC0
    exclusion <- rep("", n_all)
    append_reason <- function(mask, label) {
      idx <- which(mask)
      if (!length(idx)) return(invisible(NULL))
      text <- if (length(label) == 1L) rep(label, length(idx)) else label[idx]
      exclusion[idx] <<- ifelse(nzchar(exclusion[idx]), paste(exclusion[idx], text, sep = ";"), text)
      invisible(NULL)
    }
    reason_of <- function(name) ifelse(is.na(reasons[[name]]) | !nzchar(reasons[[name]]), "missing_unclassified", reasons[[name]])
    for (name in RAW_REQUIRED) append_reason(!is.finite(d[[name]]), paste0(name, ":", reason_of(name)))
    append_reason(!is.finite(Z), paste0("P2SPECND:", reason_of("P2SPECND")))
    append_reason(!is.finite(Y), paste0("C6R4MSCL:", reason_of("C6R4MSCL")))
    append_reason(!is.finite(W), paste0("C1_6FC0:", reason_of("C1_6FC0")))
    append_reason(is.finite(W) & W <= 0, "C1_6FC0:nonpositive_supplied_weight")
    append_reason(is.na(d$S2_ID) | !nzchar(d$S2_ID), paste0("S2_ID:", reason_of("S2_ID")))
    append_reason(!is.finite(comp$birth_weight_total_ounces), paste0("BIRTH_TOTAL_OZ:", comp$birth_weight_status))
    append_reason(!is.finite(X[, "DISABILITY_REVISED"]), paste0("DISABILITY_REVISED:", comp$RP1DISAB_missing_reason))
    for (name in FEATURES) append_reason(!is.finite(X[, name]), paste0("feature:", name, ":unavailable"))
    selected <- !nzchar(exclusion)
    if (!identical(selected, is.finite(Z) & is.finite(Y) & is.finite(W) & W > 0 &
        !is.na(d$S2_ID) & nzchar(d$S2_ID) & rowSums(!is.finite(X)) == 0L)) stop("Eligibility/reason mismatch")
    idx <- which(selected); Xs <- X[idx, , drop = FALSE]; zs <- as.integer(Z[idx]); ws <- W[idx]
    n <- length(idx)
    means <- if (n) colMeans(Xs) else setNames(rep(NA_real_, 38L), FEATURES)
    centered <- sweep(Xs, 2L, means, "-")
    rms <- if (n) sqrt(colMeans(centered^2)) else setNames(rep(NA_real_, 38L), FEATURES)
    zero_scale <- !is.finite(rms) | rms <= 0
    scaled <- matrix(NA_real_, n, 38L, dimnames = list(NULL, FEATURES))
    if (n && any(!zero_scale)) scaled[, !zero_scale] <- sweep(centered[, !zero_scale, drop = FALSE], 2L, rms[!zero_scale], "/")
    expectation_identity <- n > 0L && all(rowSums(Xs[, c("P1ELHS", "P1EHS", "P1ESC", "P1EC", "P1EMS", "P1EPHD"), drop = FALSE]) == 1)
    family_identity <- n > 0L && all(Xs[, "SGLPAR"] + Xs[, "TWOPAR"] == 1)
    nuisance_features <- setdiff(FEATURES, c("P1ELHS", if (family_identity) "TWOPAR" else character()))
    design <- cbind(intercept = rep(1, n), Xs[, nuisance_features, drop = FALSE])
    full_design <- cbind(intercept = rep(1, n), Xs)
    rank_full <- qr_rank(full_design); rank_reduced <- qr_rank(design)
    rank_union <- qr_rank(cbind(full_design, design))
    ranks_arm <- setNames(vapply(0:1, function(z) qr_rank(design[zs == z, , drop = FALSE]), integer(1)), c("0", "1"))
    nuisance_support <- do.call(rbind, lapply(c("all", "0", "1"), function(g) {
      G <- design[if (g == "all") rep(TRUE, n) else zs == as.integer(g), , drop = FALSE]
      data.frame(column = colnames(G), group = g, n = nrow(G),
        nonzero_n = colSums(G != 0),
        distinct = vapply(seq_len(ncol(G)), function(j) length(unique(G[, j])), integer(1)),
        min = vapply(seq_len(ncol(G)), function(j) if (nrow(G)) min(G[, j]) else NA_real_, numeric(1)),
        max = vapply(seq_len(ncol(G)), function(j) if (nrow(G)) max(G[, j]) else NA_real_, numeric(1)), row.names = NULL)
    }))
    identical_span <- n > 0L && expectation_identity && rank_full == rank_reduced && rank_union == rank_full
    defects <- character()
    if (!n || any(tabulate(zs + 1L, nbins = 2L) == 0L)) defects <- c(defects, "empty_selected_cohort_or_arm")
    if (any(zero_scale)) defects <- c(defects, "full_X_zero_or_undefined_RMS")
    if (!expectation_identity || !identical_span) defects <- c(defects, "reference_span_not_verified")
    if (is.na(rank_reduced) || rank_reduced != ncol(design)) defects <- c(defects, "global_nuisance_rank_defect")
    if (any(is.na(ranks_arm) | ranks_arm != ncol(design))) defects <- c(defects, "arm_nuisance_rank_defect")
    status <- if (length(defects)) "DATA_CONSTRUCTED_EXPLICIT_DIAGNOSTIC_DEFECT_REVIEW_REQUIRED" else "DATA_CONSTRUCTED_PRE_EXECUTION_SCIENTIFIC_REVIEW_REQUIRED"
    support <- list()
    add_support <- function(name, value, levels, source) {
      for (g in c("all", "0", "1")) {
        keep <- if (g == "all") rep(TRUE, n) else zs == as.integer(g)
        support[[length(support) + 1L]] <<- data.frame(variable = name, source = source, group = g,
          category = as.character(levels), n = vapply(levels, function(level) sum(value[keep] == level), integer(1)))
      }
    }
    for (name in BINARY) add_support(name, Xs[, name], 0:1, "matching_coordinate")
    for (name in ORDERED) add_support(name, Xs[, name], if (name %in% c("S2KMINOR", "P1HSCALE")) 1:5 else 1:4, "matching_ordinal_code")
    add_support("P1EXPECT_raw", d$P1EXPECT[idx], 1:6, "expectation_source_category")
    add_support("P1HFAMIL_raw", d$P1HFAMIL[idx], 1:5, "family_source_category")
    feature_summary <- data.frame(feature = FEATURES,
      type = ifelse(FEATURES %in% BINARY, "binary0/1", ifelse(FEATURES %in% ORDERED, "ordered_numeric_code", "numeric_score_count_or_mean")),
      mean_unit = means, rms_unit = rms, zero_or_undefined_RMS = zero_scale,
      min = vapply(seq_len(38L), function(j) if (n) min(Xs[, j]) else NA_real_, numeric(1)),
      max = vapply(seq_len(38L), function(j) if (n) max(Xs[, j]) else NA_real_, numeric(1)),
      distinct = vapply(seq_len(38L), function(j) length(unique(Xs[, j])), integer(1)))
    q_w <- function(value) if (length(value)) as.numeric(quantile(value, c(0, .25, .5, .75, 1), names = FALSE)) else rep(NA_real_, 5L)
    weight_summary <- lapply(c("all", "0", "1"), function(g) {
      value <- ws[if (g == "all") rep(TRUE, n) else zs == as.integer(g)]
      list(group = g, n = length(value), sum = sum(value), min_q25_median_q75_max = q_w(value))
    })
    diagnostics <- list(status = status, defects = defects, n_all = n_all, n_selected = n,
      n_excluded = n_all - n, arm_n = setNames(tabulate(zs + 1L, nbins = 2L), c("0", "1")),
      input_source_order_selected = TRUE, feature_count = 38L, feature_order = FEATURES,
      matching_no_columns_dropped = TRUE, reference_omissions = setdiff(FEATURES, nuisance_features),
      nuisance_columns = colnames(design), nuisance_ncol = ncol(design), full_intercept_ncol = 39L,
      expectation_identity = expectation_identity, family_identity = family_identity,
      rank_full_global = rank_full, rank_nuisance_global = rank_reduced, rank_union = rank_union,
      ranks_nuisance_arms = ranks_arm, identical_span = identical_span, rank_tol = 1e-10,
      dependent_columns_by_declared_QR_order = list(global = dependent_columns(design),
        arm0 = dependent_columns(design[zs == 0L, , drop = FALSE]),
        arm1 = dependent_columns(design[zs == 1L, , drop = FALSE])),
      dependent_column_diagnostic_is_not_a_drop_rule = TRUE,
      zero_scale_features = FEATURES[zero_scale], feature_summary = feature_summary,
      weights = weight_summary, exclusion_patterns = sort(table(exclusion[!selected]), decreasing = TRUE),
      raw_missing_counts = setNames(vapply(RAW_REQUIRED, function(name) sum(!is.finite(d[[name]])), integer(1)), RAW_REQUIRED),
      selected_prematurity_status_counts = table(comp$prematurity_status[idx], useNA = "ifany"),
      selected_prematurity_flag_counts = table(comp$prematurity_consistency_flags[idx], useNA = "ifany"),
      school_context_definition = schema$Context_Cohort,
      selected_school_count = length(unique(d$S2_ID[idx])), context_school_count = nrow(context),
      rank_matrices_covariates_only_no_Y_or_W_columns = TRUE,
      rank_diagnostics_condition_on_selected_cohort_and_Z_arm = TRUE,
      selection_uses_outcome_availability_not_magnitude = TRUE,
      no_models_matching_RNG_effects = TRUE, fixed_supplied_weight_unchanged = TRUE,
      not_fit_admission = TRUE, estimator_identification_consistency_or_inference_certified = FALSE)
    payload <- list(schema_version = 1L, case_id = CASE_ID, status = status,
      notation = c(observed = "Q", target = "P_S"), notation_map = schema$Notation_Map,
      functional = schema$Primary_Functional, inputs = pins, schema = schema,
      source_row_index = idx, CHILDID = d$CHILDID[idx], S2_ID = d$S2_ID[idx],
      Y = Y[idx], Z = zs, W = ws, X = Xs, X_standardized_unit_RMS = scaled,
      nuisance_design = design, nuisance_features = nuisance_features,
      context_contributor_n = context[group_index[idx], paste0(names(context_fields), "_N"), drop = FALSE],
      diagnostics = diagnostics, no_fit_admission = TRUE)
    list(selected = payload,
      selection = data.frame(source_row_index = seq_len(n_all), CHILDID = d$CHILDID,
        selected = selected, exclusion_reasons = exclusion),
      school_context = context, feature_standardization = feature_summary,
      category_support = do.call(rbind, support), nuisance_column_support = nuisance_support,
      diagnostics = diagnostics, provenance = provenance,
      scientific_acceptance = FALSE)
}

wm_ecls_read_views <- function(decoded_csv, missing_reason_csv, companion_csv, schema_dcf,
                               provenance) {
  if (!is.list(provenance) || !length(provenance)) stop("Explicit decoder/source provenance required")
  if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required for view input pins")
  paths <- list(decoded = decoded_csv, reasons = missing_reason_csv, companions = companion_csv)
  for (name in names(paths)) {
    pin <- provenance$outputs[[name]]
    if (!is.list(pin) || !is.character(pin$sha256) || length(pin$sha256) != 1L ||
        !identical(digest::digest(file = paths[[name]], algo = "sha256", serialize = FALSE), pin$sha256))
      stop("Decoder output SHA differs: ", name)
  }
  header <- names(read.csv(decoded_csv, nrows = 0L, check.names = FALSE))
  types <- setNames(rep("numeric", length(header)), header)
  types[c("CHILDID", "S1_ID", "S2_ID")] <- "character"
  d <- read.csv(decoded_csv, colClasses = types, na.strings = "", check.names = FALSE)
  reasons <- read.csv(missing_reason_csv, colClasses = "character", na.strings = "", check.names = FALSE)
  c_header <- names(read.csv(companion_csv, nrows = 0L, check.names = FALSE))
  c_types <- setNames(rep("numeric", length(c_header)), c_header)
  c_types[c("CHILDID", "birth_weight_status", "prematurity_status", "prematurity_consistency_flags", "RP1DISAB_missing_reason")] <- "character"
  comp <- read.csv(companion_csv, colClasses = c_types, na.strings = "", check.names = FALSE)
  sc <- read.dcf(schema_dcf)
  if (nrow(sc) != 1L) stop("One ECLS scientific schema record required")
  wm_ecls_construct_selected(d, reasons, comp, as.list(sc[1L, ]), provenance)
}
