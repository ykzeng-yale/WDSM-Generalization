# MEPS 2009 preparation, analysis and inference

These two Python CLIs reproduce the declared HC-129 preparation before any
model fitting, matching, effect estimation or inference. They target macOS/Linux
with Python 3.8+ and its standard library; preserve LF line endings in the
distributed JSON specification. The current reproduction was run on macOS.
Windows newline conversion has not been validated. Run without `-O`, since source and format assertions
are part of validation. No participant records or generated outputs are
included in this directory. Keep all generated files in a private location.

## Official inputs

Obtain the ASCII data archive `h129dat.zip` and SAS programming statements
`h129su.txt` from the [AHRQ HC-129 download page](https://meps.ahrq.gov/mepsweb/data_stats/download_data_files_detail.jsp?cboPufNumber=HC-129).
The [official SAS statements](https://meps.ahrq.gov/mepsweb/data_stats/download_data/pufs/h129/h129su.txt)
supply field positions, implied decimals, labels and codes; SAS is not required.
The [codebook](https://meps.ahrq.gov/mepsweb/data_stats/download_data_files_codebook.jsp?PUFId=H129)
documents the source variables.

The decoder requires these exact source identities and stops on a mismatch:

| Input | SHA-256 |
| --- | --- |
| `h129dat.zip` | `a1615e9d0b271e0a605dcdfac84091fc61c6e9b6cc577df9d960bc86fe2a1788` |
| `h129su.txt` | `84a9d48d443619cc85c79ed1b8b02faddaa4eb633bf803df669b6814e146e20e` |

## Run the two stages

From the release repository root, replace the input and output paths below
with local paths. Each output directory must not already exist; both commands
refuse to overwrite a previous preparation. Neither command relies on the
current working directory to find input data.

```sh
python3 applications/meps2009/decode_source.py \
  --raw-zip /path/to/private/raw/h129dat.zip \
  --sas-program /path/to/private/raw/h129su.txt \
  --output-dir /path/to/private/meps2009/source

python3 applications/meps2009/prepare_complete_records.py \
  --source-selected /path/to/private/meps2009/source/source_selected.csv.gz \
  --dictionary /path/to/private/meps2009/source/source_selected_dictionary.json \
  --output-dir /path/to/private/meps2009/complete-record19
```

The first command decodes selected fields in source order, retaining original
missing codes, participant identifiers and row numbers. It writes the source
CSV, a dictionary with official formats and per-variable codebook links, two
SAQ-adult pair-domain CSVs and `SOURCE_DOMAIN_RECEIPT.json`. CSV gzip containers
omit filename and timestamp headers for reproducibility.

The second command reads the adjacent `FROZEN_SPEC.json`. Its `source.path`
and `dictionary.path` are descriptive basenames: actual inputs always come
from the two required CLI arguments. The spec hash, dictionary file hash and
SHA-256 of the exact decompressed selected CSV bytes are checked before output
creation. Different gzip headers, compression levels or zlib versions may
encode the same CSV differently; their compressed-byte hashes need not agree.
The accepted gzip hash remains in `source.sha256` as provenance, while
`source.decoded_sha256` is the required content identity:
`2fe7c07f7b764c9cc36b404fb5f6db8b110df4d74111fdb0cc4a12b844a42a00`.
Receipts record both the actual input gzip hash and the decoded-content hash.
The official raw archive and SAS statement hashes remain strict. Keep the
script and unchanged spec together.

## Fixed scientific definitions

The declared pairs are White–Asian and White–Hispanic. The source rules are
`HISPANX == 2` with `RACEX == 1` for White, `HISPANX == 2` with `RACEX == 4`
for Asian, and `HISPANX == 1` for Hispanic. Filters are applied in their frozen
order: pair membership, `AGE42X >= 18`, `SAQWT09F > 0`, `TOTEXP09 >= 0`, then
complete records in the 19 declared fields.

`Z = 1` denotes White; `Y = TOTEXP09` is unchanged, including zero. The supplied
individual weight is `W = SAQWT09F`, used once. It is not multiplied by
`PERWT09F`, reestimated, normalized, trimmed or used to duplicate records.
`PERWT09F`, source flags, `VARSTR` and `VARPSU` remain provenance fields. This
preparation does not implement cluster or stratum inference.

The covariates, in order, are `AGE42X`, `SEX`, `MARRY42X`, `BMINDX53`, `PCS42`,
`MCS42`, `RTHLTH42`, `DIABDX`, `HIBPDX`, `ASTHDX`, `MIDX`, `STRKDX`, `POVCAT09`,
`EDUCYR`, `INSCOV09`, `REGION42`, `ACTLIM53`, `SOCLIM53`, and `COGLIM53`.
Official missing/inapplicable codes cause complete-record exclusion;
`MARRY42X = 6` is also inapplicable. Unknown or out-of-codebook values cause a
validation error. No imputation, clipping, winsorization or outcome-driven
choice of preprocessing rules is performed.

Continuous fields remain unchanged. Categorical fields retain source codes
and labels; the numeric design uses one-hot columns excluding the smallest
observed valid reference level. Zero-count levels and constant fields are
documented without unidentified columns. One schema and cohort per pair are
shared by all later methods. The exported design has no intercept or scaling;
matching distance and model fitting belong to the subsequent analysis.

## Outputs and interpretation

Complete-record preparation writes, for each pair, `*_complete_records.csv.gz`
and `*_numeric_design.csv.gz`, plus:

- `common_design_schema.json` and a copy of `FROZEN_SPEC.json`;
- `cohort_flow.csv`, `row_disposition_ledger.csv.gz`,
  `covariate_exclusion_ledger.csv.gz`, and `per_field_exclusions.csv`;
- `weighted_baseline_summary.csv` with weighted covariate summaries;
- `PREPROCESSING_RECEIPT.json` with file hashes and run metadata.

All of these outputs belong in private storage, including the record ledgers.
The target is the supplied-weight standardized complete-record population for
each declared pair. It is an adjusted observed-outcome disparity, not a causal
intervention on race or an estimate for all national adults. These declared
rules do not claim exact reproduction of unpublished author preprocessing.
The Python workflow reproduces input preparation only; it does not certify
later models, matching assumptions, bootstrap validity or analysis results.

## Analysis and inference phases

`run_analysis.R` runs the declared MEPS analysis from the complete-record
directory above. Install this release of `WeightedMatching` and the `digest` and
`jsonlite` packages first. It uses the normal R library search path; no library
is installed or selected by the runner. Public source/namespace function
identities must agree, including the adopted complete-contribution covariance
and centering checks. The runner does not alter the namespace or replace a
failed validation with another implementation.

Each invocation runs one requested phase for one pair. From the release root:

```sh
Rscript --vanilla applications/meps2009/run_analysis.R \
  --prepared-dir /path/to/private/meps2009/complete-record19 \
  --output-dir /path/to/private/meps2009/analysis \
  --pair white_asian --phase preflight
```

Use `white_hispanic` for the second pair. Keep the same input/output paths and
pair when proceeding through these phases:

| Phase | Calculation | Required completed phase |
| --- | --- | --- |
| `preflight` | Check source content, row order, recipe and runtime; save context | Preparation |
| `points` | Fit shared original components once; save six WM points/graphs and two PSW points | Preflight |
| `balance` | Read saved graphs and report all 49 balance fields | Points |
| `counts` | Generate the pair's one shared full-row matrix of 200 count columns | Preflight |
| `refits` | Refit shared components and PSW means for requested count columns | Counts |
| `original-summary` | Assemble original fixed-reuse/refit and separate PSW results, retaining failures | Points and all 200 refit checkpoints |
| `contribution` | Construct and replicate the complete fitted row contribution at each saved WM point | Original summary |
| `summary` | Export both distinct inference tables and the failure ledger | All six contribution records |

`balance` is independent of count and inference phases. `preflight` is also
created automatically by the first phase if absent. The source recipe is
`analysis_spec.R`; it records the predictors, source content identities,
methods, estimands and inference definitions without private execution paths.

All WM points use fixed `M = 3`, replacement, donor-normalized weights and the
original row-order tie rule. PS is the supplied-weight logistic probability;
DSM uses its distinct arm-specific unweighted-OLS prognostic score alongside
that probability. Full-X matching uses the 35 frozen numeric columns.
Matching coordinates use pooled unweighted centering and unit-RMS scaling.
PS/full-X corrections are supplied-weight arm-specific main-effect OLS;
DSM corrections use the full quadratic basis in each arm's own two scores.
No cutoff, trimming, rank-based predictor deletion or fallback model is added.

Counts have fixed `B = 200`, multinomial size equal to the original row count,
and `L'Ecuyer-CMRG/Inversion/Rejection` RNG configuration. Seeds are `202610031`
for White–Asian and `202610032` for White–Hispanic. Every method within a pair
uses that exact matrix. Refits apply the same count vector to all components
on the original rows while retaining supplied W. Each count column keeps its
predictions, component errors/warnings and PSW result or failure marker.

To bound work, `refits` accepts `--columns 1:20` or `--columns 1,2,3`;
`contribution` accepts, for example, `--methods DSM_PATE,DSM_PATT`. These select
checkpoint tasks, not a successful subset for inference. The original summary
requires all 200 refit records, and final summary requires all six contribution
records. Required failures keep the corresponding interval unavailable: there
is no filtering, redraw or successful-only variance.

## Distinct inference outputs

`original_table.csv` retains the original fixed-reuse/refit WM calculation and
the separate smooth refitted-PS Hájek comparator. Original WM point graphs,
weights and reuse loads remain fixed, with the count-weighted target
denominator. If any required count/refit column fails, its WM interval remains
unavailable and no partial WM draw vector substitutes for it. PSW retains all
200 scalar draws or failed slots and uses variance divisor B when available.

`complete_contribution_table.csv` contains a different operator: the same
original WM point, fitted components and graph feed `wm_fitted_inference`, then
`wm_bootstrap` with the original 200 counts. The saved stack retains every used
PS, PG, pooled center/variance and correction coefficient, including cross
covariances; the supplied-weight derivative is zero. For PATT only, unused
arm-1 nuisance blocks are omitted. The unchanged public
`simulations/wm_saved_full_x_stack.R` and MEPS saved-component binder supply the
PS/full-X and DSM constructions. No new fit or graph is constructed.

This table uses the exact conditional complete-row variance for its normal
interval, and reports the B-minus-1 Monte Carlo variance separately. It retains
the original-refit status/error; availability of a complete-contribution
interval does not turn an unavailable original-refit interval into an available
one. Correct working predictions and the applicable root, graph and row-moment
conditions remain substantive premises, not facts established by finite output.
Neither table claims nominal coverage from one observed dataset or complex
survey-design inference.

## Checkpoints, resumption and saved-only APIs

All outputs are written below `<output-dir>/<pair>/` and contain private
participant data or derived arrays. They must not be published. Completed RDS
checkpoints, including recorded failures, are immutable and reused. Each shared
refit column has its own checkpoint. Contribution construction is saved in
`contributions/<method>/inference.rds` before replication is attempted;
`record.rds` stores the corresponding result or failure. Summary phases only
read the required saved results and do not rerun model fitting.

A cohort lock prevents simultaneous phases. An unfinished `*.started.rds` or
stale lock stops execution rather than automatically replaying work. Inspect
the process and partial files before resolving an interrupted action. The
runner never deletes unfinished checkpoints or silently changes their contract.
Changes to the cohort, recipe, helper code or installed implementation require
a separate output directory; they cannot be mixed into a prior run.

For saved-output validation, sourcing the R files only defines functions. The
following calls operate on explicitly supplied saved objects:

```r
source("applications/meps2009/meps_analysis.R")
runtime <- meps_runtime("/path/to/release")
a <- meps_input("/path/to/private/prepared", runtime$spec, "white_asian")

# point is the captured list(ok=..., value=...) from its original checkpoint.
# components, counts and original_row are the corresponding saved records.
inference <- meps_contribution_inference(
  a, runtime, components, point, "DSM", "PATE")
result <- meps_contribution_replication(
  a, runtime, point, counts, inference, "DSM", "PATE", original_row)
```

These calls preserve the original point, graph and counts. If complete inference
was already checkpointed, pass that captured object directly to the replication
function rather than reconstructing it. They do not import or reinterpret
arbitrary legacy checkpoint keys; any external saved-object adoption must bind
its original rows, recipe, components and point explicitly.

Saved-graph balance uses the unchanged `applications/common/weighted_balance.R`
with `runtime$saved_helpers$meps_saved_balance(a, schema, points, runtime$helper)`.
`points` contains all six named captured WM points. The helper reconstructs
14 reference-category indicators for diagnostics only, giving five continuous
fields and 44 categorical indicators. The original 35-column matching metric
is unchanged. All 49 diagnostic fields are retained, including undefined SMDs.
PATE diagnostics use arm weights W+K; PATT uses White W and comparator K0.
SMDs keep their original weighted pre-match reference SD. The files
`covariate_balance.csv`, `balance_summary.csv` and `requested_graph_status.csv`
describe the saved matching distributions; they do not alter effect estimates
or establish population inference.

## Validation scope

The portable preparation reproduced the saved source and complete-record
outputs. Using the preserved model fits, graphs and counts, the public analysis
functions reproduced all 12 complete-contribution results and full covariance
terms within numerical tolerance; all 200 root draws per result had zero
numerical differences. Both original eight-row summaries, all 800 PSW draws, 3,200
component-status rows and 588 balance rows also agreed with the retained
analysis. The CLI preflight ran for both cohorts. The package was installed
from this release, and all 199 source function identities agreed.

These are preparation and saved-input code-port checks, not a fresh rerun of
the fitting/count-generation phases or evidence of nominal statistical
coverage. Original fitting/count functions were preserved and independently
source-reviewed. The original 51/200 White–Asian and 1/200 White–Hispanic
refit failures remain in the original summaries; contribution inference is
reported separately. Participant-level validation inputs and outputs remain
private.
