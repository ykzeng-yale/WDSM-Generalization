# Portable common Weighted Matching application workflow

This code slice supplies the common WM application path, with external data and
count inputs. It adds deterministic ECLS selected39 decoding/construction,
prepared NHANES/NSDUH adapters, shared nuisance-refit assembly, 18 WM rows and
descriptive saved-graph reporting. It does not download data, generate counts,
reconstruct source RNG streams, reproduce third-party comparator tables or
provide complex survey-design inference. No participant data or fitted object
is distributed. Follow the GPL-3 distribution's `COPYING` and author credits.

## Load the functions

Install the adopted `wdsmatch` package source. From this code repository root:

```r
source("applications/common/case_propensity_solver.R")
source("applications/common/application_statistics.R")
source("applications/common/application_modules.R")
source("applications/common/weighted_balance.R")
source("applications/common/prepared_input.R")
source("applications/common/application_pipeline.R")
source("applications/common/application_reporting.R")
modules <- wm_application_modules()
```

Sourcing defines functions; it does not fit, match or draw counts. The explicit
provider call checks installed `wdsmatch`/`digest` APIs. It also uses the current
package's unexported `.wm_ns_ols`/`.wm_ns_basis`; those are explicit version
dependencies, not a stable future public API. Loading `balance_schema.json`
requires `jsonlite`. Base R/stats supply other direct statistical dependencies.
API presence does not prove identity to the adopted package source. Keep a
source/version manifest for the installed package and these application files.
The default driver authenticates the loaded solver/provider definitions against
the pinned sibling `case_propensity_solver.R`/`application_modules.R` before
calling the provider. It checks the eight loaded statistical definitions and
their lexical helper resolution against the pinned `application_statistics.R`,
and requires the actual `wdsmatch` namespace. Supplied alternate module functions
cannot use the default recipe.

## Explicit prepared inputs and counts

`wm_application_input` binds `data`, exact source-order `id`, finite `Y`, binary
`Z` with both arms, positive supplied known `W`, raw matching matrix `X`, separate
nuisance `design`, propensity formula `variables`, explicit study/outcome/baseline
labels, units and PS weighting. No row selection or recoding occurs in this
generic adapter. Every observation-level object has the same rows/order.
`data$A` equals Z; all formula variables are finite numeric columns. X and design
retain their declared coordinates/spans, including the intended intercept.
The adapter records serialized-object sample/dataset SHA-256 values with digest;
these differ from file SHA-256 and depend on object representations.
It also records `nuisance_sha256` from the explicitly ordered design columns and
double-valued matrix, formula response/intercept/variables and ordered formula
data, and PS weighting. This canonical nuisance input omits row labels but
retains source row order; it supplements the unchanged sample/dataset hashes.

The named adapters in `prepared_input.R` implement the exact accepted case
contracts in the adjacent case documentation. They retain the source probability
PS fit for NHANES and ECLS and unit PS fit for NSDUH. W remains supplied and frozen
under both choices. General new prepared cases can use `wm_application_input`
with an explicitly justified specification.

Supply a numeric/integer n-by-200 matrix, an **identical** `count_ids` vector and
nonempty `count_provenance` list. Every column is checked for integral nonnegative
finite entries, total n and positive count mass in both arms. Invalid columns
remain records/failure columns; they are never redrawn or removed. Malformed
matrix shape or wrong IDs fails before analysis. A recommended provenance list
includes case ID, source file SHA, matrix object SHA, original seed/RNGkind and
start/end-state source records when available, and the declared uniform
Multinomial(n;1/n,...,1/n) protocol. Providing labels is not proof that a supplied
matrix follows that protocol; authenticate the actual source/count artifact.

There is no count generator or seed-replay helper in this slice. NHANES primary
3783, canonical3790, corrected ECLS7220 and NSDUH17471 have distinct case/count
contracts. NSDUH outcomes share the same dataset, IDs and count matrix; their
outcome-dependent nuisances remain different.

## Run the same 18 WM rows

Given explicitly prepared `a`, `counts`, `count_ids`, `count_provenance`:

```r
target_scope <- c(PATE = "declared all-record supplied-W observed target",
                  PATT = "declared recipient supplied-W observed target")
inference_scope <- paste("Nominal fixed-graph/refit algorithm interval; known W;",
  "causal, national, complex-design and consistency assumptions require separate justification")
result <- wm_application_run(a, counts, count_ids, count_provenance,
  modules = modules,
  component_source = "applications/common/application_statistics.R",
  target_scope = target_scope,
  inference_scope = inference_scope)
schema <- wm_application_read_balance_schema("applications/ecls/balance_schema.json")
balance <- wm_application_balance(result, schema)  # choose this case's schema
wm_application_save(result, "my-new-local-result-directory", balance)
```

The explicitly invoked driver fits shared components once at the point and once
per supplied count column, then calls the existing functions for PS, distinct-arm
DSM and full-X, M=1/3/5, PATE/PATT. It uses self-normalized donor weights,
replacement matching and `source_random` ties with seed20260917+M. Point maps use
pooled unit-RMS scaling. DSM has separate (PS,PG0)/(PS,PG1) maps, and six-column
quadratic correction; PS/full-X have the source W-weighted linear correction.
PG is fitted without W, while correction LS uses multiplicities times W/maxW.
These equations are those of `application_statistics.R`; no new estimator is
defined by this driver. User source tie matching may internally use its explicit
tie RNG, restoring the caller stream according to the package contract. Count
generation and ordinary rematching bootstrap are never invoked.

`component_callback=function(a, modules, m) ...` can supply the same captured
component structure, including PS, PS_score, PG, DSM and linear arm components,
and full_X_scores for the all-one point. The default calls `app_components`.
Custom callbacks declare a different nuisance specification unless they retain
the same equations; they require separate scientific validation. Every run must
supply `component_source`. The default verifies that source file's SHA against
the accepted `application_statistics.R` bytes. For a custom callback, also supply
`component_recipe=list(recipe_id="my_recipe_v1", source_sha256=SOURCE_FILE_SHA,
metadata=list(...))` and the actual implementation file as `component_source`.
The file's SHA is recomputed before and after execution, without sourcing it.
Source paths and callback closure environments are never serialized into this
binding. Recipe metadata uses named lists of plain finite numeric, logical or
character vectors (or NULL), with list keys sorted and numeric values normalized
to doubles; no attributes, objects, missing values or absolute filesystem paths
are allowed. Ordered vectors retain their declared order.

The default recipe also binds the actual R version/platform and installed
`wdsmatch` implementation: DESCRIPTION/NAMESPACE, every installed R/lazy-load
and `libs` payload file, and canonical formals/body signatures for loaded
namespace functions. The loaded solver/provider and eight application helpers
are checked against source-only parsed function definitions. Signatures contain
syntax hashes, without srcref attributes, native pointer addresses, environment
serialization or installation paths. Both bytes and loaded signatures are
recomputed before/after a run, including fully reused runs. A changed installed
or loaded implementation therefore changes the exact recipe binding or fails
the accepted default source check.

Custom callbacks require explicitly supplied `modules`,
`application_statistics_source="applications/common/application_statistics.R"`
and nonempty `component_dependencies`, for example
`list(cache_adapter=list(source_file="my-cache-adapter.R", source_sha256=FILE_SHA))`.
These dependency files are verified before/after and only named SHA values are
persisted. The common module must be an actual package namespace or a named list
of the five required functions (`.wm_ns_ols`, `.wm_ns_basis`, `wm_match`,
`wm_bootstrap_refit`, `wm_wdsm_fit`); its actual function syntax, and namespace
payload when applicable, enter the recipe. The source module must explicitly
supply `wdsm_case_fit_propensity`. The canonical common statistical helper file
remains fixed even when the nuisance callback changes.

The source descriptor authenticates supplied file bytes, not the semantic
identity of an arbitrary in-memory callback to that file. The caller must load
the callback from that accepted source, declare all changed model/cache options
and relevant source/package dependencies in metadata, and authenticate those
dependencies against prior accepted records. The installed package/module
source manifest is still required. Function syntax/payload identity does not
identify arbitrary mutable state in custom closure environments. Declare every
relevant custom model/cache/runtime option and dependency, and authenticate the
callback/modules against the original accepted specification; this is a caller
responsibility rather than an inferred guarantee. Shared refit
assembly captures callback errors, retains all200 column diagnostics and aligned
mean0/mean1 placeholders, and never silently changes an arm model's rank/span.
An externally captured PS point/refit cache can be passed to `app_components`
through such an explicit callback, with its column/source binding retained by
the caller; the driver does not infer cache keys from a private filesystem.

All PATE/PATT calculations share the original fixed graph, supplied W and reuse
loads, through `app_replicate`/`wm_bootstrap_refit`. B=200 is requested even if
columns fail; empirical variance uses divisor B, not B-1. Both arm matrices must
retain n-by-200 shape, including PATT's unused mean1 matrix. A PATE interval fails
when either required arm prediction/count column fails; PATT needs control
predictions and valid counts. The full failed-column list, captured diagnostics,
raw/corrected points and unavailable interval remain available. No survivor-only
variance, changed basis, pseudoinverse, ridge, redraw or failure exclusion is
performed. A returned interval does not establish nominal coverage in real data.

## Reuse accepted outputs without duplicate fits or matching

Pass existing `components`, `refits` and all18 `precomputed_points` captured
records to `wm_application_run`, with `reuse_binding` containing:

```r
list(sample_sha256=a$sample_sha256, dataset_sha256=a$dataset_sha256,
     counts_sha256=digest::digest(counts, algo="sha256"),
     nuisance_sha256=a$nuisance_sha256,
     component_recipe_sha256=wm_application_recipe_binding(
       NULL, NULL, "applications/common/application_statistics.R",
       modules=modules)$sha256,
     components_sha256=digest::digest(components, algo="sha256"),
     refits_sha256=digest::digest(refits, algo="sha256"),
     points_sha256=digest::digest(precomputed_points, algo="sha256"))
```

Authenticate this binding against accepted source records; the caller cannot
turn an arbitrary object into scientific acceptance by hashing it. The driver
recomputes the nuisance-input hash and implementation-source/recipe binding on
every run; any reused category must match both exact values. A custom callback
uses its declared recipe/source in `wm_application_recipe_binding` instead of
the default example above. Changing the design, PS formula variables/data,
PS weighting, callback source, authenticated runtime implementation or declared recipe cannot silently mix newly
computed point components with old refits. The driver also checks
data/order/W/target/M/tie/maps/predictions against the same components,
full200 count/prediction shapes and object hashes before reusing points. Exactly
all18 captured point keys are required, so a missing supplied point does not
trigger rematching. Captured failed points remain unavailable. With all three
supplied, no nuisance callback or matching call executes; fixed-graph replicate
algebra/reporting can be compared to saved accepted results without new fits,
count draws or neighbor search. Supplying only some component/refit categories
explicitly allows the remaining required category to be computed only under the
same accepted nuisance-input and source/recipe binding. Historical reuse
bindings without the two new fields fail before a callback, fit or matching
operation. Adding those fields requires documentary authentication of the
original nuisance specification; it is not a numerical refit or validation.

Source graph/load arithmetic and the provenance of these supplied objects still
need prior acceptance. Structural/hash checks are not an independent proof that
a cached graph is the nearest-neighbor solution. The package refit calculation
also checks that original graph loads/predictions reconstruct the point.

## Descriptive reporting and output scope

`wm_application_balance` uses the saved point objects and prepared data only.
Schemas describe diagnostic columns; they do not change matching X. NHANES uses
17 diagnostics, corrected ECLS has all38 coordinates plus22 declared ordinal
category indicators, and the NSDUH schema declares its14 binary coordinates.
Absent/sparse categories stay declared. Fixed pre-W PATE pooled-arm or PATT
treated SDs use the unchanged balance helper: numeric reliability variance and
Bernoulli mu(1-mu) for indicators. Post graphs use W+K for PATE and treated W/
control K0 for PATT. Zero/undefined denominator labels, NA SMDs, raw differences,
category proportions, masses, ESS and increases in imbalance remain visible.
No standardization by postmatching SD or outcome correction is substituted.

The writer creates a fresh local directory, refuses overwrite/resume, saves
components/refits/counts, every captured point/interval row and optional balance,
and records partial failures. The prepared input retains its nuisance hash; the
completion object retains the path-free nuisance recipe and reuse binding for
later authentication. It is simple local IO, not an execution controller.
Output files contain participant-level data and must remain with the analyst;
only these source files/documentation are part of the code distribution.

These application results do not certify causal identification, conditioning or
centering assumptions, overlap, score-model correctness, target transport,
complex-design variance or coverage. Source_random/mixed support retains its
qualified empirical inference scope. All empirical/scientific interpretation
must distinguish the case's actual exposure timing and target.


## Validated reproduction scope

Complete supplied point/component/refit/count bundles were replayed for NHANES
primary3783, corrected ECLS7220, and both NSDUH17471 outcomes through this driver
with the authenticated default helper/package recipe. All72 rows and all12600
available B200 draws retain their original statuses. Native row estimates,
draws, variances, SEs, intervals and root-n draws agree exactly; the maximum
finite CSV row difference is5.11e-15. Nine ECLS PATE intervals remain unavailable,
with the same69 failed columns each; no variance from surviving columns is used.
The replay generates no fits, neighbors or counts and does not invoke nuisance
callbacks. It validates complete-bundle reproduction. Separately, the first
default new-fit workflow was run once on each of these four prepared
study/outcome cases, using the authenticated current default recipe and supplied
accepted counts without callbacks or cached components/refits/points. All72
point rows and14400 prescribed draw slots were accounted for:63 rows retain
all12600 available B200 draws, while nine ECLS PATE intervals remain unavailable
with the same69 failed columns each. No variance from the131 surviving columns
is substituted. Shared components, complete200-column arm-refit
matrices/diagnostics, point maps/graphs/loads, variance/SE/CI/root-n outputs and
1890 descriptive balance rows were compared. Native numerical gaps are0;
maximum finite CSV rounding difference is5.10702591327572e-15. Keys, statuses,
nonfinite categories and numerical shape/order/value comparisons retain their
fixed exact or scaled1e-10 rules. This is recursive numeric/shape/status parity,
not complete storage-class, all-attribute or serialized-byte identity. It
validates these default prepared-input workflows, not arbitrary custom recipes,
all possible inputs, upstream agency-data construction, causal identification
or sampling-inference assumptions. The default score guard applies the package's
same double/no-dimname canonicalization before strict value/dimension/order
comparison; numerical perturbations remain rejected.
