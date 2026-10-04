# SPDX-License-Identifier: GPL-3.0-only
# Frozen scientific recipe; no filesystem or execution approval contract.
meps_analysis_spec <- function() {
  list(
    `purpose` = "MEPS 2009 HC-129 supplied-weight standardized expenditure disparities.",
    `target` = "Supplied-weight standardized complete-record population for each White-comparator pair. PATE averages the adjusted White-minus-comparator disparity over the combined pair; PATT uses the White reference distribution. No causal intervention on race, all-national-adult target, or exact original-author pipeline reproduction is claimed.",
    `metadata_columns` = list("source_row", "DUPERSID", "PANEL", "race_group", "Z", "Y", "W", "TOTEXP09", "SAQWT09F", "PERWT09F", "HISPANX", "RACEX", "RACETHNX", "SAQELIG", "SFFLAG42", "VARSTR", "VARPSU"),
    `cohorts` = list(
      `white_asian` = list(
        `schema_key` = "White-Asian",
        `n` = 11729,
        `n_white` = 10215,
        `n_comparator` = 1514,
        `comparator` = "Asian",
        `zero_expenditure_n` = 1419,
        `count_seed` = 202610031,
        `design_file` = "white_asian_numeric_design.csv.gz",
        `accepted_gzip_sha256` = "30e77177a671f9c28a90adf436746cdc6a117720a1c39e0559485aeca9378cd0",
        `decoded_sha256` = "472f1a82c2d2d2c7ffeaff0d2fbc9916f2c58386eab8543d4c2c6cf093868d90"
      ),
      `white_hispanic` = list(
        `schema_key` = "White-Hispanic",
        `n` = 15808,
        `n_white` = 10215,
        `n_comparator` = 5593,
        `comparator` = "Hispanic",
        `zero_expenditure_n` = 2945,
        `count_seed` = 202610032,
        `design_file` = "white_hispanic_numeric_design.csv.gz",
        `accepted_gzip_sha256` = "e3d10ca29a3f4f3a42d0929da2b099ed1f25b1ed8ef13a187d67665f16bce257",
        `decoded_sha256` = "e08cda805f4a2a8baa3a070b2791527a33e2eecd5f9360ec06f1eca1b4bd3150"
      )
    ),
    `outcome` = list(
      `field` = "TOTEXP09",
      `alias` = "Y",
      `units` = "2009 USD",
      `zero_expenditure` = "retained",
      `transformation` = "none; no logarithm, clipping, winsorization, or outcome-based exclusion"
    ),
    `group` = list(
      `Z1` = "White",
      `Z0` = "Asian or Hispanic according to pair",
      `contrast` = "White minus comparator"
    ),
    `known_weight` = list(
      `field` = "SAQWT09F",
      `alias` = "W",
      `times_applied` = 1,
      `PERWT09F` = "provenance only; never multiplied into W",
      `status` = "supplied and frozen; no estimated/new nonresponse weight"
    ),
    `design` = list(
      `source_fields` = 19,
      `numeric_main_effect_predictors` = 35,
      `regression_columns_with_intercept` = 36,
      `columns` = list("AGE42X", "SEX__2", "MARRY42X__2", "MARRY42X__3", "MARRY42X__4", "MARRY42X__5", "MARRY42X__7", "MARRY42X__8", "MARRY42X__9", "MARRY42X__10", "BMINDX53", "PCS42", "MCS42", "RTHLTH42__2", "RTHLTH42__3", "RTHLTH42__4", "RTHLTH42__5", "DIABDX__2", "HIBPDX__2", "ASTHDX__2", "MIDX__2", "STRKDX__2", "POVCAT09__2", "POVCAT09__3", "POVCAT09__4", "POVCAT09__5", "EDUCYR", "INSCOV09__2", "INSCOV09__3", "REGION42__2", "REGION42__3", "REGION42__4", "ACTLIM53__2", "SOCLIM53__2", "COGLIM53__2"),
      `categorical_encoding` = "Frozen reference indicators; category codes are not continuous measurements",
      `continuous_encoding` = "Unchanged Age42, BMI53, PCS42, MCS42, education",
      `interactions` = "none in PS, PG, PS correction, or full-X correction",
      `rank_failure` = "Retain failure; no dropping/releveling/ridge/model replacement based on fitted rank or effects"
    ),
    `matching` = list(
      `families` = list("PS", "DSM", "full_X"),
      `estimands` = list("PATE", "PATT"),
      `M` = 3,
      `replacement` = TRUE,
      `method` = "self_normalized",
      `tie_rule` = "row_order",
      `tie_tolerance` = "unchanged wm_match default: 64 * .Machine$double.eps",
      `metric` = "Euclidean after pooled unweighted unit-RMS column standardization; no imposed exact categorical matching",
      `scores` = list(
        `PS` = "Weighted-logistic fitted probability, standardized in pooled observed rows",
        `DSM` = "Distinct arm-specific (weighted-logistic probability, unweighted-OLS PG_z), pooled-standardized separately by donor arm map",
        `full_X` = "All 35 frozen numeric predictor columns, pooled-standardized; no intercept in matching coordinates"
      ),
      `dimension_scope` = "35 is the ambient coded coordinate count, not a claim that 35 continuous density coordinates or the paper geometry conditions hold."
    ),
    `nuisance` = list(
      `PS` = "Same 36-column supplied-W Bernoulli/logistic root; accepted zero-start solver, maxit100, no clipping/retry/new model",
      `PG` = "Unweighted arm-specific OLS on the same 36-column design",
      `PS_and_full_X_correction` = "Supplied-W arm-specific OLS main-effect predictions on the same 36-column design",
      `DSM_correction` = "Supplied-W arm-specific complete six-term quadratic in its own two pooled-standardized scores; unchanged .wm_ns_basis order",
      `count_refits` = "Same original n rows and same count vector for PS, both PG, pooled score means/RMS, and correction regressions. W retained. Predictions always returned for all original rows."
    ),
    `PSW` = list(
      `PATE` = list(
        `White` = "m W Z/e",
        `Comparator` = "m W (1-Z)/(1-e)"
      ),
      `PATT` = list(
        `White` = "m W Z",
        `Comparator` = "m W (1-Z)e/(1-e)"
      ),
      `estimator` = "Difference of two separately normalized Hajek outcome means",
      `point_probability` = "Reuses the exact same weighted-logistic fit as WM",
      `replication` = "Same full-n count matrix and each count-refitted PS; separate smooth PSW procedure, not WM fixed-reuse inference",
      `variance_divisor` = "B",
      `interval` = "Original PSW point plus/minus qnorm(.975) times empirical replicate SD"
    ),
    `replication` = list(
      `B` = 200,
      `count_distribution` = "Multinomial(n,rep(1/n,n)), independently predeclared seed per pair",
      `RNG_kind` = list("L'Ecuyer-CMRG", "Inversion", "Rejection"),
      `shared_within_pair` = "One exact n-by-200 count matrix shared by all six WM results and both PSW results",
      `original_WM_rule` = "Exported wm_bootstrap_refit through unchanged app_replicate cached-prediction callback. Original graph, supplied W and original K are frozen; per-count target denominator uses counts times original W. No ordinary rematching.",
      `variance_divisor` = "B",
      `interval` = "Original WM point plus/minus qnorm(.975) times empirical fixed-reuse SD",
      `failure` = "Keep every count column and every component prediction/error. No filtering, redraw, successful-only variance, silent model fallback, or invalidation by an unused component. app_replicate prechecks all required WM columns: any required failure leaves the full WM interval unavailable and no WM scalar draw vector is evaluated for that method. The original point and all counts/component records remain. PSW separately retains every per-column scalar draw or failed-slot marker.",
      `scope` = "Nominal working-model iid/known-W calculation only. WM sampling/refit-agreement and score/geometry/centering conditions, and PSW root/overlap conditions, are not established by fitting this application. B200 has Monte Carlo variation."
    ),
    `complex_design` = list(
      `VARSTR_VARPSU` = "Preserved in private inputs; not used to introduce a cluster/stratum variance theorem",
      `reported_inference` = "No complex-design validity claim"
    ),
    `complete_contribution` = list(
      `operator` = "Complete fitted row contribution evaluated at the same saved original point, graph, supplied W and fitted components; wm_fitted_inference followed by wm_bootstrap with original full-n counts.",
      `nuisance_stack` = "Every used PS, PG, pooled score center/variance, and correction coefficient; known-weight derivative zero; unused treatment-arm blocks omitted only for PATT.",
      `variance` = "Exact finite-row conditional variance from centered complete contributions divided by n; B-minus-1 empirical count-draw variance retained separately as Monte Carlo output.",
      `interval` = "Original point plus/minus qnorm(.975) times square root of conditional variance.",
      `distinct_from_original_refit` = TRUE,
      `assumptions` = "Correct full-X working predictions and applicable actual-graph/root/row-moment conditions remain substantive unverified premises; no Gaussian expenditure or complex-design claim."
    ),
    `checkpoint_policy` = "Read completed immutable checkpoints, including failures. Retain all 200 count/refit columns. A started action without its completed checkpoint stops for inspection; no silent refitting or redraw."
  )
}
