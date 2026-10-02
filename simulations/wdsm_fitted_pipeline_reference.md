# Qualified fitted WDSM pipeline reference

This standalone bridge connects `wdsm_nuisance_stack_reference.R`, the package's
original self-normalized matching estimator, and `fitted_wm_components.R`.
It implements the known-weight, full-X-correct branch for PATE and PATT.
It is not a package export and does not infer population assumptions from data,
model names, fitted residuals, or a `CorCor`/`MisCor` label.

Source the two reference dependencies and this file. Supply the original
unit-multiplicity output of `wm_wdsm_nuisance_stack_reference()`:

```r
source("simulations/wdsm_nuisance_stack_reference.R")
source("simulations/fitted_wm_components.R")
source("simulations/wdsm_fitted_pipeline_reference.R")

out <- wm_wdsm_fitted_pipeline_reference(stack, M = 3L,
  contract = list(mode = "algebra_only",
    justification = "Deterministic implementation check; no population claim."))
replication <- wm_wdsm_fitted_count_reference(out, counts)
```

`algebra_only` computes the point, slopes, covariance blocks, augmented rows,
and supplied-count draws, but withholds sampling standard errors and intervals.
It is the appropriate mode for finite-data parity tests. A positive empirical
variance alone cannot establish a positive limiting variance or model validity.

## Declared theorem contract

For inference, pass `mode = "declared_full_x"`, a nonempty `justification`
identifying the applicable proof/model conditions, the following named logical
vector `assumptions` with every entry explicitly `TRUE`, and a separate logical
`conditional_refit_expansion`. These entries are caller declarations, never
tests performed by the bridge. Every returned object says
`assumptions_verified = FALSE`.

| Declaration | Required content |
|---|---|
| `iid_bounded_sampling` | Independent identically distributed observed rows; bounded outcomes/designs and weights bounded above and away from zero, as in the proved bounded theorem. This excludes treating clusters as independent individual rows. |
| `known_probability_weights` | The supplied positive weights are the fixed known target-law weights, possibly treatment dependent. There is no estimated-weight parameter or omitted weight derivative. |
| `target_identification_overlap` | The target law, potential-outcome identification, positive treatment overlap and the relevant target normalizer hold. |
| `correct_full_x_means` | The population quadratic correction is the full-X conditional outcome mean in both arms for PATE, or in the control arm for PATT. Correct PG plus its inclusion in the correction basis is one sufficient model construction. A correct target PS alone does not assert this. |
| `current_score_geometry` | The two-dimensional current maps satisfy the proved fitted-map support, density, bounded derivative/displacement and local-index conditions, including the source-form standardization coordinates. Baseline covariance rank by itself is insufficient. Matching uses fixed M and exact Euclidean distances with deterministic row-index ties. |
| `regular_joint_nuisance_influence` | All PS/PG/pooled-moment/correction roots are regular, nonsingular and statistically accurate, with the full joint influence and fitted-map transfer conditions. Numerical parameter error is negligible at root-n scale. |
| `positive_limit_variance` | The full adjusted limiting root-n variance is strictly positive. |

`conditional_refit_expansion = TRUE` additionally declares the conditional
regular expansion of the **same-count full-stack** original nuisance refit.
It permits the qualified statement that the original fixed-reuse refit and
sampling law have the same limiting variance in this branch. With `FALSE`,
the bridge makes no original-refit validity declaration. One-step draws do
not actually fit any bootstrap nuisance model.

The constructive bounded iid CorCor/MisCor application is proved in the
project's `actual_wdsm_fitted_pipeline_conditions.md`, Section 8, and
`fitted_wdsm_pipeline_bridge.md`. The historical Gaussian/cluster experiments
and generic CorMis/MisMis models are not certified by this bridge. PATT
requires only the control PG/correction; no unused treated model is fitted.

## Source fitting and donor conventions

The supplied stack preserves weighted logistic PS, ordinary unweighted arm
PG regression, unweighted pooled centering/scaling with divisor n, and weighted
full quadratic bias correction. Its shared PS moment blocks eliminate duplicate
columns while preserving all joint covariances. Known offsets and explicitly
supplied design matrices follow the stack API. There is no extra standardization
inside matching.

`wm_match(..., method = "self_normalized", variance = FALSE)` uses these
actual fitted maps and predictions. Each query selects M nearest opposite-arm
donors and uses donor fractions W_j divided by their selected weight sum.
This is the predecessor's mathematical donor rule with `tie_rule = "row_order"`
and exact distances. It is not a claim of equality to the predecessor's default
random near-tie resolution, or of bitwise equivalence between stable Euclidean
norm arithmetic and direct squared-distance arithmetic at floating-point ties.
The focused check compares exact donor identities against the pinned source.

One fixed common normalization `w = W / max(W)` cancels from the point, rows
and gradient. `original_reuse_normalized` uses these same units; multiply by
`original_weight_scale` to obtain the predecessor's raw-weight K when that
product is representable. The package point and PATT denominator retain their
original ratio normalization; both sampling variances are V/n, where n is the
total observed sample size.

The stack's tighter PS solver tolerance and normalized equation weights were
already independently compared with the pinned predecessor. This bridge does
not change those controls or promise identical numerical iterations. Nonunit
stack multiplicities are rejected here: a reweighted empirical nuisance fit
is not automatically an original iid-row sampling fit.

## Complete slope and variance

Write d_zi for the full parameter derivative of the fitted quadratic mean,
including the generated PS, PG, pooled centers and variances. Let K_i be the
incoming normalized weighted reuse. The bridge independently computes the
original fixed-K refit slope

\[
 \widehat a_A = \frac{\sum_i\{[(1-Z_i)K_i-Z_iw_i]d_{0i}
                  +[(1-Z_i)w_i-Z_iK_i]d_{1i}\}}{\sum_iw_i},
 \qquad
 \widehat a_T = \frac{\sum_i[(1-Z_i)K_i-Z_iw_i]d_{0i}}
                     {\sum_iZ_iw_i}.
\]

It compares this with the exact query-level frozen-graph gradient. Known
weights give zero weight derivative. Full-X baseline centering gives zero
population graph-transport drift, so the supplied branch uses
`b_hat = a_hat`. The finite equality of these two algebraic gradients is
checked; the population sampling interpretation is a separate theorem claim.
It is not an assertion of exact conditional centering on an outcome-fitted
random graph.

The nuisance influence is the complete stacked sandwich
`ell_i = -J_hat^{-1} psi_i`. Let t_i be the actual normalized centered
contribution from the fitted matching graph. Center both t and ell and form
`u_i = t_i + ell_i' b_hat`. Then

\[
 \widehat V=\frac1n\sum_i u_i^2
 =\widehat V_0+2\widehat b^\top\widehat C
       +\widehat b^\top\widehat\Sigma\widehat b.
\]

All cross-block covariances remain in the calculation. No Fisher-information
replacement or universal variance-reduction claim is used. The full-X fitted
theorem uses the diagonal limiting variance, with zero reciprocal component;
the bridge does not subtract a finite noisy reciprocal product from these rows.

For an n-by-B integer count matrix with each column summing to n, the one-step
root is `sum_i (N_ib - 1) u_i / sqrt(n)` and its draw is the point plus that
root divided by `sqrt(n)`. Under multinomial(n, 1/n) counts its exact conditional
root variance is `mean(u^2)`. The code checks count values but does not certify
their generating law. Empty-arm count columns remain defined, since no count
ratio or refit is evaluated. These finite one-step draws need not equal the
original nonlinear refit draws, even where their asymptotic laws agree.
