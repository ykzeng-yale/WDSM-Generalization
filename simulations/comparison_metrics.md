# Performance metrics for the next comparator study

Source `simulations/comparison_metrics.R` and call
`wm_comparison_metrics(records, plan, reference_scenario = ..., reference_method = ...)`
with the study's predeclared reference. The defaults retain the original WDSM
`CorCor`/`PSM_M1` convention. Those two-score PS/PG specification labels are not
a required scenario taxonomy for arbitrary-dimensional WM. This helper does not
create a missing reference, change old simulation outputs, or establish
inference validity.

The plan has one row per `estimand`, `design`, `overlap`, `n`, `scenario`, and
`method`, with an integer `requested >= 2` and logical `interval_required`.
Use `group_cols` to add other DGP identifiers (e.g. dimension or population
version). Every factor defining a different simulation population must be in
the grouping columns; scenarios within a group change working models, not
the target population. Keep a distinct method identifier for each estimator
and inference procedure. Each planned cell requests replicate IDs
`1:requested`.

Records contain those keys, `replicate`, `target`, `estimate`, `lower`, `upper`,
`variance`, and `status` (`ok` or `failed`). Retain failed records and their
targets. `variance` is the variance of the treatment-effect estimator, not the
variance of its root-n limit; convert the latter by the actual dataset sample
size before supplying it. `NA` is allowed when variance is not provided. A
point-only method must have `interval_required=FALSE` declared in the plan.

The return value contains `records_complete`, `publication_ready`, a cell audit,
and `metrics`. These flags concern the reporting gate only; they do not certify
the simulation design, inference theorem, or manuscript readiness.
The audit reports requested/received/completed/failed/invalid/missing records
and usable interval counts. If any required record or interval fails, or the
reference cell is absent, `metrics` is NULL. Undefined relative efficiency or
nonfinite derived metrics also block the table, even if every raw record is
complete. This implements the original
study's complete-success publication gate; it does not silently publish
successful-only tables. Raw failed records remain available for a separately
labeled operational/failure analysis.

When the gate passes:

- Signed relative bias is `100*(mean(estimate)-target)/abs(target)`. Absolute
  relative bias is the absolute value of that aggregate quantity. Both are
  undefined for target zero, with an explicit status.
- Relative efficiency is empirical variance of the predeclared reference
  scenario/method divided by empirical variance of the method. The original
  WDSM study uses **CorCor PSM_M1**; other WM studies supply their own explicit
  reference. That same reference is used for every analysis specification
  within the group. Both variances use divisor R-1. A zero/nonfinite variance
  yields an explicit audit status and blocks the publication table. This
  compares variance, not bias or MSE.
- Variance calibration is mean reported estimator variance divided by
  empirical variance. It has its own name and report count; it is not RE.
- Coverage uses the actual endpoints and causal target. Its interval count,
  binomial MCSE and exact pointwise 95% binomial limits are returned. For a
  declared point-only method with some optional intervals, coverage is
  explicitly labeled as conditional on available optional intervals.
- Bias, bias MCSE, empirical SD, RMSE, mean reported SE and mean interval
  length are also reported. Replications must be independent for the MCSE
  interpretation. Different methods may be paired within a replication.

No relative-efficiency MCSE is currently supplied; a paired ratio MCSE requires
the actual cross-method and cross-scenario replication linkage. This helper
does not infer such linkage merely from equal row counts or IDs.

Run the hand-calculated reporting checks from the code checkout:

```sh
Rscript validation/check_comparison_metrics.R
```

These checks use small deterministic records, not additional simulation
experiments. Historical production without a PSM reference cannot be promoted
to the new comparison study by relabeling variance-calibration ratios.
