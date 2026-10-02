# Reference components for fitted weighted matching inference

This standalone reference implements the covariance and one-step replication algebra for the reviewed fitted-map branches. It is not an exported package API, an automatic nuisance estimator, or a calibration result. The original corrected self-normalized point estimate is preserved.

## Interfaces and units

- wm_graph_gradient_reference differentiates the exact query-level estimator while holding donor lists fixed. Supply n-by-p derivatives of the original raw individual weights and fitted mean predictions. It includes both donor normalization and the target denominator. It returns separate weight and prediction components and their sum. It does not estimate graph transport.
- wm_fitted_variance_reference takes the joint nuisance influence matrix and separately supplied smooth and graph sensitivities. Both sensitivities use point-estimator units per nuisance parameter; do not add an additional gamma division. Parameter columns must use one order. All nuisance cross covariances and the covariance with the complete matching row are retained.
- wm_fitted_count_reference takes the fitted variance result and an n-by-B count matrix. Each column must sum to n. It implements full-slope one-step multinomial draws; it does not reproduce the original nuisance-refit ratio.

The exact fitted sensitivity is the sum of the smooth derivative and graph transport. The root-n variance is V0 + 2 b'C + b'Sigma b. Estimator variance divides by total n for both PATE and PATT. The PATT denominator is already incorporated in the complete rows and sensitivity units.

Weight derivatives are normalized by the fit's fixed positive common scale. The derivative of a common weight scale cancels from every donor fraction and target ratio. Tests explicitly verify this cancellation and invariance to multiplication of raw weights by 1e-100 or 1e100.

## Theorem boundary

The reference requires all used matching dimensions to be at least two. For PATE, callers explicitly choose distinct_rarity, full_x, or common_field; PATT uses patt. The full_x branch requires zero graph transport. A string declaration is not proof of identification, centering, current-map geometry, nuisance regularity or consistent transport estimation. The result records assumptions_verified = FALSE.

These branches have zero or asymptotically vanishing reciprocal covariance. The reference therefore uses complete diagonal rows and no reciprocal subtraction. It calls the existing reciprocal inference validator solely for strict reconstruction of the finite fit; that function's fixed-map interval and variance are not used.

Generic scalar outcome-fitted maps, nonzero reciprocal branches, growing nuisance dimension, arbitrary learners, and dependent survey designs are outside this reference. Separate scalar theory or later interfaces must be explicitly qualified. Code accepting numerical matrices does not make their influence or derivative interpretation valid.

## Replication contract

With centered augmented rows v_i, one-step roots are sum_i (N_i - 1) v_i / sqrt(n), and draws add the root divided by sqrt(n) to the original estimate. Under Multinomial(n; 1/n,...,1/n) counts their exact conditional root-n variance is mean(v_i^2). No Monte Carlo is needed to compute it.

The reported Monte Carlo variance instead summarizes the supplied B draws, with divisor B. For B = 1 that descriptive quantity is zero, while the analytic conditional variance remains available. Empty-arm count draws remain defined because the one-step scheme uses the original denominator. This differs from the original refit ratio's possible fit or denominator failures.

Counts are validated algebraically; the function cannot verify that a supplied count matrix was randomly generated from the stated multinomial law. Sampling validity additionally requires consistent full slopes/influences, the fitted matching theorem and empirical row Lindeberg conditions.

## Deterministic validation

From the code-release directory, run Rscript validation/check_fitted_wm_components.R with an optional output JSON path.

The checks compare analytic derivatives with an independent frozen-query finite difference for PATE/PATT and M = 1,3,5, using distinct two- and three-dimensional maps. They compare weight derivatives with the existing independently checked weight API; check scaling, parameter alignment, covariance blocks and guards; and enumerate every one of the 4^4 ordered multinomial assignments for n = 4. This is formula verification, not a simulation of bias, efficiency or coverage.
