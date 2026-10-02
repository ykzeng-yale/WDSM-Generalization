# Bounded fitted WDSM calibration runner

This standalone runner implements the reviewed bounded source-form calibration
plan. It is not a package export, and creating its configuration does not launch
an experiment. The runtime pilot and production allocations are distinct:

| Configuration | Independent datasets | Original-refit count columns |
|---|---:|---:|
| `config/bounded_wdsm_runtime_pilot.json` | 4 at each n=500/2000 | B=39 on all 8 datasets |
| `config/bounded_wdsm_calibration.json` | 1000 at each n=500/2000 | B=399 on IDs 1,5,...,997 at each n |

The pilot is for runtime, resources, numerical roots and failure diagnostics.
Its few datasets must not be used to assess bias, coverage or variance calibration.
The scientific design is source_good, r=.25, CorCor/MisCor, PATE/PATT and M=1/3,
with the exact bounded-population truths. No source_poor, additional dimension,
or model grid is silently added.

## Replaying after the WDSM package interface was added

The completed September 29 study used
`config/bounded_wdsm_calibration_sharded.json` and its pinned 24-file package
inventory. That historical configuration remains unchanged. Exact historical
execution also requires the matching archived installed lazy-load artifacts;
reinstalling the old sources on another host does not establish those hashes.
The current checkout instead provides a separately identified forward replay.

Its static source admission is
`config/bounded_wdsm_current_source_admission_20261001.json`. It lists the
35 admitted R sources, DESCRIPTION/NAMESPACE and the seven unchanged scientific
modules. The preparer pins this admission's own SHA256 and rejects altered,
missing or extra sources. Inventory admission does not extend the statistical
scope of the bounded study or endorse historical estimated-weight helpers.
The preceding admissions remain available as historical records. The October 1
admission refreshes five reviewed R comment/returned-contract text hashes.
Estimator arithmetic, function signatures, guards, count handling, dispatch and
the seven scientific modules are unchanged; its package check skipped tests
and examples and reuses the earlier numerical validations. The current
admission includes the qualified scalar/full-X inference interface; the bounded
study's two-dimensional analysis grid, seeds and scientific modules are unchanged.

Install `jsonlite` and `digest` if needed, then install this checkout into a
separate R library. For example, from the checkout root, using an existing
empty library directory and an existing output parent:

```sh
R CMD INSTALL --library=/existing/path/wm-library .
R_LIBS=/existing/path/wm-library Rscript simulations/prepare_bounded_wdsm_replay.R /existing/path/new-replay.json
```

This prepares a new configuration and does not run simulations. It verifies the
complete admitted source inventory and all installed function names, bodies
and formals, then resolves that installation's runtime hashes. It retains the
original scientific configuration, seed allocation and historical configuration
hash. The new file identifies current-source forward execution; it is not
evidence of historical runtime identity, numerical parity or calibration.
Existing output or partial files are rejected. Later source changes require a
separately reviewed admission rather than dynamic acceptance of the checkout.

The unchanged runner can check admission and reconstruct only the stream/grid
allocation without generating observations, counts or fits:

```sh
R_LIBS=/existing/path/wm-library Rscript simulations/run_bounded_wdsm_calibration.R --config /existing/path/new-replay.json --output-root /existing/new-preflight-root --shard-index 1 --shard-count 20 --preflight-only
```

The output directory must be new and empty. Preflight performs deterministic
target quadrature and allocates RNG streams; it is not a simulation or a
statistical validation. The 20-shard allocation matches the completed study.
Executing the full study is a separate resource-bounded operation. Neither
preparation nor preflight starts it. Runtime hashes can differ across local
installations; unchanged source bodies and designs do not promise bitwise
floating-point agreement across all R versions or platforms.

## Invocation and admission

Run the script with an existing supervisor-owned output directory:

```
Rscript simulations/run_bounded_wdsm_calibration.R --config simulations/config/bounded_wdsm_runtime_pilot.json --output-root /absolute/new-run-root
```

The caller must supply the reviewed local resource supervisor: one worker,
60 minutes, aggregate memory 6 GiB, output 1 GiB for the pilot. The script itself
does not enforce aggregate process-tree resource bounds. It creates `records/`
exclusively and refuses to resume or overwrite that tree. The script path locates
the checkout, so the working directory does not have to be the checkout.

Appending `--preflight-only` verifies source/package pins, prepares the exact
requested grid, streams and targets, and returns **before any dataset, count
column or fitted model is generated**. Preflight still creates a new records
tree; its directory cannot subsequently be reused for execution.

The JSON freezes the helper/runner sources, every package R source, installed
package lazy-load artifacts, package version, reviewed-plan hash and full design
allocation. Before any generation, the runner compares all installed package
function bodies and formals with the pinned checkout definitions. R version and
dependency versions are recorded. A different installation may require a newly
reviewed configuration even when its version string matches. Source/config/package
hashes are checked again after the requested grid is processed.

## Randomness and failure ownership

L'Ecuyer-CMRG streams are allocated in sample-size/replication-ID order, with
one independent stream for data and one for counts per dataset. Pilot and
calibration have distinct fixed master seeds. All stream states and the requested
grid are saved before the first dataset. Each dataset and its full count matrix
are saved, with pre/post stream states, before fitting. Counts are shared across
both models, estimands and M. No failed dataset or count draw is replaced.

The helper first processes the independent known-map oracle. Within each model
and estimand it fits the full nuisance stack once and reuses each count refit's
predictions across M=1/3. PATT never fits an unused treated prognostic model.
Original count ratios retain original W/K and do not rematch. The separate
one-step count calculation may fail without erasing valid original-refit draws.
A failed inference assembly may retain a valid original matching point; fallback
point errors are saved separately. Model failures do not erase independent
known-map oracle results.

Every fit records the complete mean equation, Jacobian-derived Newton residual,
PS score residual, parameter vector, rank/conditioning diagnostics, warnings,
elapsed time and solver controls. Exactly one zero-start precision refinement
is allowed when only root accuracy fails: 100/1e-12 becomes 200/1e-14. Structural
failures do not receive a numerical retry. The count-specific constructor is
copied from the pinned full constructor with only unused influence-matrix and
sandwich calculations removed; predictions and complete equation/Jacobian
checks are unchanged and independently compared on the existing n=256 fixture.

## Records and units

Top-level records include `config.json`, `provenance.json`, `streams.rds`,
`requested_grid.csv`, `target.rds` and a final `completion.json`. Each dataset
directory contains saved rows/counts, compact attempt records, fitted objects,
20 point-method records, each requested original-refit draw, one-step results,
count summaries and timestamped stage events. Terminal errors preserve partial
records and add explicit requested-row/draw accounting. A supervisor interruption
may leave a dataset without completion; its missing records remain identifiable
from the prewritten requested grid and allocation.

Analytic `root_n_variance` uses the total observed n for PATE and PATT.
Original-refit `source_point_variance_B` uses point-estimate units and divisor B.
Its source SE/interval retains that divisor. The separately labeled
`root_variance_unbiased_success_law` equals n B/(B-1) times that variance and
requires every prescribed draw to succeed; it is not an unconditional failed-draw
completion. One-step draws already use root-n units. Their exact conditional
variance is an algebraic identity, not independent evidence of sampling validity.
Within-count MCSEs preserve paired covariance for refit/one-step differences.

Execution completion means the requested grid was processed; individual fit or
draw failures may remain. Neither execution completion nor preflight success
certifies the theoretical assumptions, coverage, moment convergence or empirical
calibration. Outer-dataset aggregation and its successful-intersection/MCSE audit
remain a separate step before interpreting production results.


The later explicit `ps_weighting` wrapper addition is admitted by the current
forward-replay preparer after independent default-path equivalence review.
Historical saved configurations retain their original pins and are not
rewritten. Use a new prepared forward configuration for the current checkout;
older invocation examples below or above document their historical runtime.
No additional bounded calibration is needed to validate the source-faithful
prospective option, which has separate original-survey reduction evidence.
