# Balance from an existing weighted matching graph

`weighted_balance.R` contains the unchanged diagnostic functions used on the
reviewed NHANES, ECLS-K and NSDUH saved matching graphs. It requires base R and an
existing `wm_match` object from `wdsmatch`. Sourcing it defines functions;
it does not fit models, construct neighbors, draw counts or change the RNG.
This module belongs to the application code and adds no package API.

## Inputs and use

Supply the original selected rows in their original order, numeric/integer
`Z` coded 0/1, the same positive supplied `W`, and a diagnostic data.frame
with finite numeric columns. The saved object
must use `method="self_normalized"` and the same rows, W and PATE/PATT label.
For a `wm_wdsm_fit` result, the matching object is its `fit` component.
No participant data, fitted objects or study-specific private paths are
included with this module.

```r
source("applications/common/weighted_balance.R")
# point: an existing wm_match object; data: its rows' diagnostic covariates.
diagnostic <- wb_diagnostics(data, schema)
reference <- wb_reference(diagnostic$X, Z, W, diagnostic$binary,
                          point$estimand)
balance <- wb_saved_graph(point, diagnostic$X, Z, W, reference)
balance$table
```

`schema` is a list with one entry per diagnostic. Every entry provides
`source_variable`, unique `diagnostic`, `label` and `kind`. Supported kinds:

- `binary_numeric`: an existing 0/1 column.
- `continuous`: a finite numeric column. This name specifies the numeric
  diagnostic convention; an ordinal numeric-code column can use it without
  asserting a continuous matching distribution.
- `binary_indicator`: one declared category, with `category` and all valid
  `allowed_codes`. Supply every declared ordinal level as its own entry if
  evaluating category proportions, including levels with zero frequency.

The diagnostic matrix is separate from the matching coordinates. Adding a
category indicator here does not alter the saved graph. Missing or unexpected
codes fail explicitly; this helper does not define analytic eligibility,
impute missing values, merge sparse categories or filter unfavorable results.

## Definition of the reported balance

Before matching, means use W within each arm. Each denominator stays fixed:
PATE uses the pooled pre-match weighted arm SD, and PATT uses the pre-match
weighted treated SD. For numeric columns the weighted covariance uses
normalized weights with correction `1 - sum(p^2)`; binary indicators use
Bernoulli variance `mu_z * (1 - mu_z)`, where `mu_z` is the pre-match
W-weighted proportion in arm z. These conventions should be aligned before
comparing another program's SMDs. A zero or undefined denominator produces
an unavailable SMD with its retained reason, not a reported zero imbalance.

The stored incoming reuse masses are converted back from the package's
scaled-W units. After matching, PATE arm distributions use original W plus
their incoming mass; PATT uses treated W and control incoming mass. Outputs
retain weighted means, fixed SDs, signed/absolute SMDs, target-mean deviations,
ESS and mass residuals. Different targets have different fixed denominators;
these SMDs are not a cross-target ranking.

The summaries describe raw matching distributions. Outcome bias correction
does not change their graph weights. Improvement is not guaranteed, and all
worsened or unavailable diagnostics should remain in a complete report.
This module does not estimate an effect variance or establish identification,
conditional centering, overlap, model truth, complex-design inference or
confidence-interval coverage. Full-distribution balance does not follow from
these mean/proportion diagnostics.

## Provenance and scope

The source bytes are identical to the function module adopted for 18 NHANES
and 18 corrected ECLS-K graphs and then36 NSDUH graphs. Those reviews checked saved weighted moments
and graph-mass arithmetic independently; they did not certify a causal
interpretation or survey-design standard errors. The common package retains
its separate inference contracts. Study-specific cohort preparation,
comparators and full reproduction entrypoints are not supplied by this shared
helper; their source/data requirements must be documented with those exports.
