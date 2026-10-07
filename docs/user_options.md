# User options and predecessor matching interfaces

`wdsmatch` implements a specified fixed-count weighted estimator. It is not a
drop-in replacement for every matching design in another package. Start with
`wm_fit()` for supplied matching coordinates, `wm_wdsm_fit()` for fitted
arm-specific double scores, or `wm_model_fit()` for finite fitted model lists.
Use `wm_match()` when supplying mean predictions yourself. The installed R help
gives every argument, default, return component and inference condition.

## Main choices

| Choice | Current WM interface and default | Interpretation |
| --- | --- | --- |
| Rows and weights | `Y`, binary `Z`, strictly positive finite `weights` | Same row order throughout; missing values must be handled before fitting. W is supplied and known. Use `rep(1, length(Y))` for unweighted reductions. |
| Target | `estimand="PATE"` or `"PATT"`; default PATE | PATE uses all-row target weights, PATT the treated-weight denominator. These are population targets under the stated weighted identification assumptions. |
| Donor count | `M=3L` | Exactly M distinct donors per query, with reuse across queries. No variable ratio or matching without replacement. M is fixed in the asymptotic framework. |
| Matching variables | `scores0`, optional `scores1` for PATE | One or more coordinates per arm; maps may differ. PATT uses only the control-direction map and rejects unused arm1 inputs. `wm_match()` requires matrices, including a one-column matrix for PSM; `wm_fit()` also accepts a scalar-score vector. |
| Distance and scaling | Euclidean distance on supplied coordinates | No implicit scaling in supplied-map APIs. A fixed coordinate transformation can express another fixed quadratic metric. Estimating a transformation requires its applicable fitted-map argument. |
| Fitted DSM | `wm_wdsm_fit(ps_design, pg0_design, pg1_design)` | Uses separate `(PS, PG0)` and `(PS, PG1)` maps. Designs include explicitly supplied intercepts. Pooled score centering/scaling and complete quadratic correction are retained. |
| Fitted model lists | `wm_model_fit(ps_models, pg0_models, pg1_models)` | Each descriptor supplies `design`, optional known `offset`, and for PS optional `weighting`. Matching dimension counts model outputs, not regression coefficients. PS-only model lists still receive correction by default; they are not raw PSM. |
| PS fitting convention | `ps_weighting="probability"` or `"unit"` in WDSM; per-PS descriptor in model lists | Changes PS fitting only. Supplied W remains in donor fractions, target weighting and weighted correction. PG fits use ordinary arm OLS. |
| Donor rule | `method="self_normalized"` by default | Uses `W_j/sum(W_selected)` in each donor set. `"stabilized"` is a separately documented residual-imputation estimator with conditional-weight normalizers; it is not a weight-normalization switch. |
| Bias correction | Supplied `mean0`/`mean1`, or `wm_fit(regression="polynomial", degree=1L)` | Correction fits the donor arm's own matching variables. For a raw supplied-map point use `wm_match(..., variance=FALSE)` and omit means. Raw and corrected estimates are returned separately. |
| Spline correction | `regression="spline"`, declared order, mesh exponent, moments and support | Requires explicit rate/support choices and the documented training/evaluation scheme. No automatic bandwidth/model tuning or term deletion. |
| Fitted correction alternatives | `wm_model_fit(correction=...)` | Default complete quadratic; opt-in `local_polynomial` or `quadratic_qr` has its own inputs and inference restrictions. Other learners cannot be passed as an unchecked model flag. |
| Exact restrictions and folds | `strata`, `fold_id` in supplied-map APIs | Restricts donor graphs and corresponding nuisance fits. These matching cells do not implement dependent survey-stratum or cluster inference. |
| Ties | `tie_rule="row_order"` by default | Exact distance then original row index; no RNG. `"source_random"` requires `tie_seed`, retains its own squared-distance tolerance policy, and supports point/fixed-reuse empirical algebra only. Sampling inference is unavailable for that policy. |
| Inference request | Supplied-map `variance=TRUE`; fitted wrappers `inference="none"` by default | Supplied-map variance requires relevant means and its stated residual/nuisance conditions. `"full_x"` on fitted wrappers is conditional on correct full-X predictions and complete fitted-score conditions. Obtaining an SE does not verify those assumptions. |
| Replication | `wm_bootstrap(method="contribution", B=999L)` or `method="fixed_reuse"` | Contribution draws use the selected retained law. Fixed reuse requires full-n count columns and preserves original W and donors. Default fitted wrappers refit the retained complete quadratic prediction pipeline; nested matching fits hold means fixed unless `refit` is supplied. |
| Confidence intervals | `conf.level=.95`; bootstrap `interval="normal"`, `"basic"` or `"none"` where supported | The selected inference route governs permitted intervals. Fixed reuse supports normal/none and divisor B. Nested quadratic mixtures use their documented basic interval and reject normal. Analytic conditional variance and finite-B empirical variance are distinct outputs. |
| Reproducibility | `seed`, `tie_seed`, supplied `counts`; retained calls/recipes | Explicit seeds restore the caller's represented RNG state; documented Box–Muller restrictions apply. Counts retain original row order and must satisfy the route's checks. Failed prescribed draws are not deleted or redrawn. |
| Numerical controls | `controls`, `qr_tol`, basis/stack allocation limits | Solver/rank/allocation controls do not change the target, drop requested terms or prove asymptotic numerical availability. Unsupported options error rather than being silently forwarded. |

Known weights may depend on `(Z,X)` or only X. Products of sampling/response
weights require the corresponding recovery assumptions and are used once.
Estimated-weight uncertainty is outside the current framework; the historical
`wm_logistic_weights()`/`wm_weight_adjust()` compatibility helpers are not part
of this paper's supplied-weight inference. No cluster-robust SE option is implied.

## What the predecessor arguments mean

| Predecessor | Related options | Consequence for comparison |
| --- | --- | --- |
| `Matching::Match()` | `M`, `estimand`, `X`, observation `weights`, `Weight`/`Weight.matrix`, `replace`, `ties`, `caliper`, `exact`, `BiasAdjust`, `Var.calc` | Its observation-weighted path selects a donor-weight mass rather than simply M distinct records. `Weight` controls distance separately. Bias adjustment and Abadie–Imbens variance options are not WM bootstrap aliases. Unit weights recover ordinary donor counts, but distance, ties and correction still need alignment. |
| `MatchIt::matchit(method="nearest")` | `ratio`, `replace`, `reuse.max`, `distance`/`link`, `mahvars`, `exact`, `antiexact`, `caliper`, `discard`, `reestimate`, `m.order`, `s.weights`, `normalize` | Provides a matched-design object, with nearest-neighbor ATT/ATC directions. Defaults differ from WM. Sampling weights can enter model fitting and later combine with matching weights; they do not automatically implement WM donor-set normalization. Calipers/discarding can alter the analyzed target. |
| `dsmatchATE()` / `dsmatchATT()` (1.7.1) | `Y/X/A`, `method`, supplied PS/PG, model selectors/design matrices, `varest/boots/mc/ncpus` | ATE has no M argument or forwarding ellipsis. ATT can forward M on several Matching-backed paths; its DSM path forwards only without a caliper. Its retained correction/reuse calculations do not generally normalize M>1. Native comparisons therefore use M=1. ATE's DSM uses a common three-coordinate map; ATT's is two-coordinate. Fitted PS uses the linear predictor. These differ from WDSM's arm-specific probability/PG maps. |
| Original `wdsmatchATE()` / `wdsmatchATT()` (0.2.1) | `Y/X/Z/weights`, optional `model.ps/model.pg` formulas or supplied scores, `M=5`, `sampling`, `use.bias.correction`, `varest/boots/alpha`, `tie.seed/tie.tolerance` | The legacy APIs and defaults remain available. Bootstrap uses the ambient `set.seed()` stream; there is no bootstrap `seed` formal. Score, correction, tie, PS-weighting and fixed-reuse choices must be retained for exact WDSM reductions. General WM defaults to M=3 and adds explicit maps/model lists and qualified contribution inference. |

These are intentionally unsupported in the general WM interface: optimal/full/
genetic/cardinality matching, ATC/ATO targets, matching without replacement,
variable M, calipers or data-driven overlap trimming, and unrestricted imputation
or machine-learning nuisance learners. Adding a similarly named argument would
require a specified estimator, target and inference result. Use the predecessor
package when its different matching design is the intended analysis.

With W=1, verify finite point equality only after aligning the rows, target,
coordinates and scale, metric, replacement, exact donor count, ties, mean
correction and normalization. Variance equality additionally requires the same
inference procedure. Equal point estimates do not make an analytic variance,
contribution replication and prediction-refit variance identical in finite data.

## Inspect the result

```r
print(fit)                         # compact estimate and inference status
s <- summary(fit)
s$statistics                      # retained estimate, correction, SE, variances
s$donor_usage                     # graph counts, not a balance assessment
fit$graph$edges                    # wm_match/wm_fit: original row identifiers
fitted_wrapper$fit$graph$edges     # fitted wrappers retain the nested graph
replication$conf.int               # interval from the chosen replication route
```

`summary()` does not calculate a new interval or repair unavailable uncertainty.
Balance must be assessed against the appropriate weighted target; donor-use counts
alone do not establish balance or identification. The application workflows
include target-specific saved-graph balance code. Record package version, source
commit, model/scaling choices and inference operator in an analysis report.

Primary references: [Matching manual](https://cran.r-project.org/web/packages/Matching/Matching.pdf),
[MatchIt reference](https://kosukeimai.github.io/MatchIt/reference/matchit.html),
[nearest-neighbor options](https://kosukeimai.github.io/MatchIt/reference/method_nearest.html),
[sampling-weight workflow](https://kosukeimai.github.io/MatchIt/articles/sampling-weights.html),
and [wdsmatch manual](https://cran.r-project.org/web/packages/wdsmatch/wdsmatch.pdf).
The [dsmatch comparison](https://github.com/Yunshu7/dsmatch/tree/006a939ad84664eacf80f933ad3e4fc8fc12511d)
is bound to the retained 1.7.1 source; it is not a claim
about every future version. Installed help and [the method guide](theory_method_guide.md)
give the full WM contracts.
