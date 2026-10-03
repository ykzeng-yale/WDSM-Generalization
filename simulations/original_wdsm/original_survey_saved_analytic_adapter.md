The companion adapter adds analytic PS/X6 rows to an already saved original-study
weighted batch. It does not replace or rerun the original DSM analytic adapter,
the 19-method report, fitting, matching, or fixed-reuse count bootstrap.

```r
source("code-release/simulations/wm_saved_full_x_stack.R")
source("code-release/simulations/original_wdsm/original_survey_saved_analytic_adapter.R")
companion <- ows_saved_analytic_batch(
  data, weighted, asNamespace("wdsmatch"),
  saved_stack = wm_saved_full_x_stack, alpha = .05, M = 3L
)
```

Use the exact original `data` and `weighted` objects; the default path checks
the recorded data hash. For artifacts restored by the original lossless reader,
the study runner can pass `authenticated_input`: the complete decoded data and
context, original input reference and receipt, and the authenticated source chain.
The adapter then checks strict data/attribute identity and the unchanged
historical context/batch hash. This requires the runner's original file, receipt,
payload and block authentication; the binding is not itself a substitute for
those checks. A fresh serialization hash identifies the decoded in-memory
representation only and never replaces the historical hash. The original
case/input references and pairing checks remain required.
The helper is a sourced simulation function, not a package export.
One complete saved nuisance stack is evaluated per family and reused across
the requested recorded M values. Omit `M` to use every recorded order. Each inference call receives that M's original point and
graph. Complete same-row nuisance covariance, the finite prediction slope, and
the full-sample PATT normalization are retained. Supplied weights are known.

PS retains its original raw-covariate main-effect/interaction correction as
declared by the saved scenario. This is the common correct-full-X prediction
branch for a scalar matching variable, not the uncorrected `scalar_psm` route.
X6 retains its complete quadratic correction through an exactly equivalent raw
polynomial parameterization. No matching distance is changed. The adapter does
not use saved DSM success as a prerequisite for either family.

`rows` contains `WM_PS_analytic_M*` and `WM_X6_analytic_M*`, original point
sources, variance/interval statuses, saved refit variances and their recorded
B/divisor, and an analytic/refit variance ratio. A failed saved refit status
is retained even when its recorded variance is numeric; its ratio is unavailable.
Missing or failed points,
nuisance stacks and inference calculations remain failed rows. `assemblies`
retains the common inference results; `nuisance` retains the two complete stacks.
Save these as companion results without overwriting the frozen study outputs.

These numerical results do not establish correct outcome models, the applicable
sampling/geometry/moment conditions, complex-design variance validity, or
agreement with the source fixed-reuse bootstrap. Such assumptions remain
explicit and unverified. In particular, a successful misspecified-mean row is
an empirical comparison, not a theorem admission. No new coverage claim follows
from this saved-array calculation.
