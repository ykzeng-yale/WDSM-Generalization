# Weighted Matching application statistical interface

These modules provide eight unchanged application statistical functions,
one unchanged checked propensity solver, and a thin installed-package provider.
They do not supply a complete raw-data-to-results application pipeline. No
participant data, count matrix, fitted object, scientific report, private cache
or execution controller is included.

## Load and dependencies

Distribute these three R files together within the GPL-3 `wdsmatch` code
distribution and its `COPYING` file. From the repository root, source the
three application files into the same environment:

```r
source("applications/common/case_propensity_solver.R")
source("applications/common/application_statistics.R")
source("applications/common/application_modules.R")
modules <- wm_application_modules()
```

Sourcing only defines functions. The last explicit call checks installed
`wdsmatch` and `digest` namespaces; it does not install packages, select a
private library, fit a model or generate counts. Base R and `stats` supply the
remaining direct dependencies. The provider requires these package functions:

- Exported `wm_match`, `wm_bootstrap_refit` and `wm_wdsm_fit`.
- Lower-level, unexported `.wm_ns_ols` and `.wm_ns_basis` in the current
  validated package. Their use is explicit: they are not promised to remain
  stable across future package changes. The provider checks their availability
  and required formal arguments.

`wm_wdsm_fit` is the package's DSM wrapper; these modules do not redefine it.
It supports distinct arm maps through its nuisance systems. Its point wrapper
uses `inference="none"` for `tie_rule="source_random"`; a `wm_wdsm_fit` result's
matching object is its `fit` component. The extracted applications use the
lower-level functions to share components across PS, DSM and full-X graphs.
This provider is an API-availability check, not a proof that an arbitrary
installed version has the identical validated numerical implementation. The
release must bind the adopted package source separately.

## Exact input contract

The statistical input `a` is an explicitly prepared list. All observation-level
objects retain the same selected rows and order. This module does not choose a
cohort, recode missingness, infer survey variables or validate raw eligibility.

| Field | Requirement and use |
|---|---|
| `n` | Number of selected observations, consistent with every row/vector. |
| `Y` | Finite numeric outcome vector of length `n`. |
| `Z` | 0/1 treatment vector of length `n`, containing both arms. |
| `W` | Finite strictly positive supplied weight vector of length `n`; kept fixed in all refits and replicates. |
| `X` | Finite numeric `n`-row matrix of the full matching coordinates, in the declared order. |
| `design` | Finite numeric `n`-row nuisance-regression design, including the intended intercept; its declared span is retained and rank failures are not repaired. |
| `variables` | Character vector naming the propensity formula's columns in `data`. |
| `data` | Prepared data.frame with the same rows, `A` equal to `Z`, and every column named by `variables`. |
| `ps_weighting` | Explicitly `"unit"` or `"probability"`; callers must validate this choice. |
| `study`, `outcome`, `baseline`, `units` | Nonempty scalar labels used by `app_row`; they do not select algorithms. |

`variables`, `design` and `X` need not have identical dimensions. In the corrected
ECLS application, all 38 matching coordinates remain in `X`; the equivalent
nuisance span uses 36 named columns plus an intercept, omitting two exact
aliases. A separate cohort loader must justify this span and preserve ordering.

Supply an explicit numeric/integer **`n` by 200** count matrix `counts`. Every
entry is finite, nonnegative and integral; every column sums to `n` and assigns
positive multiplicity to both arms. `app_count_columns(a, counts)` returns all
200 validation records and serialized-vector SHA-256 values using `digest`.
The fixed-200 contract is deliberate: the package's general
`wm_bootstrap_refit` permits other B, but these exact application helpers do not.

No `app_counts` or `app_rng` function is exported. The source application streams
are distinct: the NHANES primary case uses 3783 rows and seed 20261011; its
canonical source case uses 3790 rows and the separately recovered Linux source
manifest. Corrected ECLS uses 7220 rows and seed 20261024. Both primary recipes
declare Mersenne-Twister/Inversion/Rejection, B=200, and full-n uniform
multinomial counts. Those facts do not construct a count matrix or certify
cross-platform equality. A release driver must require and bind the correct
case's actual matrix, row order, provenance and seed-state records. It must not
reuse a matrix of another size, silently generate one or redraw failed columns.

## Functions and retained statistical behavior

| Function | Operation |
|---|---|
| `app_capture(action)` | Captures success, warnings and failure metadata, including failed columns and partial results. |
| `app_count_columns(a, counts)` | Checks all 200 supplied columns and records totals, arm masses, hashes and errors. |
| `app_ps(a, modules, m)` | Fits the checked zero-start logistic score using multiplicities `m`, with unit PS weights or `m * W / mean(original W)`. |
| `app_components(a, modules, m, ps)` | Computes prognostic and correction components, reusing a supplied captured PS fit when available. |
| `app_predictions(components, family, estimand)` | Selects required correction predictions; PATT only needs control predictions. |
| `app_match(a, modules, components, family, M, estimand)` | Calls the common `wm_match` with supplied W, fixed M and explicit `source_random` tie seed `20260917 + M`. |
| `app_replicate(point, modules, counts, refits, family, estimand)` | Uses every supplied count column and cached prediction column through common `wm_bootstrap_refit`; no rematching or weight fitting. |
| `app_row(a, family, M, estimand, point, replication)` | Retains corrected/raw points, interval availability, errors and labels in one comparison row. |

Use `family` exactly `"PS"`, `"DSM"` or `"full_X"`, and `estimand` exactly
`"PATE"` or `"PATT"`. The caller must validate these labels. The PS graph is
one-dimensional. DSM uses separate `(PS, PG_z)` maps for the two arms; full-X
uses all declared `X` columns. Score coordinates use pooled multiplicity-RMS
standardization, not W-weighted standardization. Matching is with replacement
under the common donor-normalized estimator. Point components use all-one
multiplicities; full-X point scores are computed only for that all-one case.

Prognostic LS uses `n * m / sum(m)` and no W. Linear full-X/PS correction and
DSM quadratic correction use those multiplicities times `W / max(original W)`.
DSM's six-column quadratic basis is the current `.wm_ns_basis`; linear
correction uses the supplied `design`. The checked PS solver retains zero-start,
convergence, full-rank, finite-probability and normalized score-residual guards.
Supplied nonconstant W is never re-estimated.

`refits` is an explicit input to `app_replicate`, with this structure:

```r
refits$means$DSM$mean0       # n-by-200 numeric matrix
refits$means$DSM$mean1       # n-by-200 numeric matrix
refits$means$linear$mean0    # n-by-200 numeric matrix
refits$means$linear$mean1    # n-by-200 numeric matrix
```

Column b must correspond exactly to `app_components(a, modules, counts[, b])`
or the same equations with an explicitly cached PS fit for that column. DSM
uses the DSM pair; PS/full-X use the linear pair. Every failed component is
retained with its diagnostics and nonfinite prediction column. Both `mean0`
and `mean1` must be present as aligned `n`-by-200 numeric matrices for the chosen
cache kind, including PATT: the unchanged helper calls `colSums()` on both matrices.
PATT uses only `mean0` values; its `mean1` values need not be finite and may be
retained as nonfinite placeholders. PATE also requires finite `mean1` values.
The driver must retain the complete diagnostic record, rather than treating
a prediction matrix as a certification of success.

Replication keeps the original graph, W and incoming reuse masses fixed and
uses full-n count columns for both PATE and PATT. Its empirical variance divisor
is B=200. Any invalid count or required nonfinite prediction column fails the
interval with the full failed-column list. No survivor-only variance, redraw,
ridge, pseudoinverse, coefficient deletion or changed correction is supplied.

For an explicit already-prepared input/count/refit set, the minimal calls are:

```r
components <- app_components(a, modules)
point <- app_capture(function()
  app_match(a, modules, components, family, M, estimand))
replication <- app_capture(function() {
  if (!point$ok) stop(point$error)
  app_replicate(point$value, modules, counts, refits, family, estimand)
})
row <- app_row(a, family, M, estimand, point, replication)
```

These calls fit and calculate only when explicitly invoked by the caller.
Pass a captured replication result, including a captured failure; this example
does not use `app_row`'s nullable default as a replacement for failure metadata.
The example assumes `a`, `counts`, `refits`, `family`, `M` and `estimand` have
already been supplied and validated; it is not a runnable raw-data driver.

## Provenance, checks and remaining pipeline components

The eight statistical functions originate from the local validated
`application_engine.R` snapshot (SHA-256
`102d0db2ced4535ba898d44ea5346d2533e174cda4263b39e3fbd4d12797cdbb`).
The solver originates from the supplied revised WDSM `WDSM_adaptive.R` snapshot
(SHA-256 `757ed59932b57adaf22cbe9500324a9299b83ddbb191a48343e37b5b50090f67`).
These are **local** source identities. The historical projects
`SW_DSM` and `wdsmatch` are provenance context, not a claim that these revised
bytes equal an unchanged public commit. GPL-3 and the WDSM/wdsmatch author
credits are retained. The nine extracted assignments pass identical parsed-AST
comparison with the source, with source-location metadata omitted; their source
span bytes, including final line feeds, are identical. All ten top-level
expressions in the three R files are function definitions. Verification
only parsed code and inspected dependencies: it never sourced those files,
invoked functions, fitted models, generated counts or changed an RNG stream.

This extraction reuses the previously reviewed application arithmetic. It does
not establish causal identification, the population centering/model/geometry
conditions, national transport, complex survey-design inference, or nominal
coverage. In particular, `source_random` has its explicitly qualified empirical
fixed-reuse inference scope; it does not inherit a continuous-score sampling
theorem merely because replication returns an interval.

The separate [common application workflow](../APPLICATION_PIPELINE.md) now
supplies selected39 ECLS decoding/construction, prepared NHANES/NSDUH adapters,
explicit supplied-count/refit assembly, complete failure reporting and saved
graph balance integration. Its complete-bundle reuse path was compared with
all72 adopted WM rows and all12600 available B200 draws, retaining all nine
unavailable ECLS PATE intervals and their69 failed columns. Full agency raw
NHANES/NSDUH cohort preparation, count-stream construction and predecessor
comparators remain explicit external/source requirements. Private ownership,
resource admission, caches and process workflow are not source dependencies;
internal verification scripts and extraction records are not distributed.
