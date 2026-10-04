# Reproducing the aggregate results

This page covers the historical `primary/`, `supplement/`, `first_stage/` and
`bounded_wdsm/` archives. Current principal iid, original WDSM and Lenis compact
record commands and their execution scope are in the
[results index and reproduction matrix](README.md).

These commands reproduce preserved historical archives. The optional
estimated-weight cells in `first_stage/` are outside the current supplied-weight
WM framework; reproducing their arithmetic is not validation of that framework
and those cells are excluded from the current paper. No new estimated-weight
study is needed for the active research goal.

Run the commands below from the repository root, where `DESCRIPTION`, `R/`, `simulations/` and `results/` reside. R is required; the deterministic aggregate check uses only base/recommended R packages and does not require package installation. Results and simulation scripts belong to the source checkout and are excluded from the installed package.

## Saved aggregates, benchmarks and figures

```sh
Rscript --vanilla results/reproduce_benchmarks.R
Rscript --vanilla results/reproduce_first_stage.R
Rscript --vanilla results/reproduce_bounded_wdsm.R
```

This checks the public artifact MD5/size manifest, reconstructs the saved beta vectors and covariance matrices from CSV, checks their moment identities, and repeats both benchmark joins. It compares every column of all 576 primary and 54 supplemental rows with the saved tables at a tolerance allowing CSV rounding. It performs no estimator fitting or random-number generation. SHA256 fingerprints are also provided in `results/artifact_hashes.csv` for independent verification. The retained R 4.4.2 check passed; its output is in `results/provenance/aggregate_check_v2.log`.

The second command additionally reconstructs all 234 first-stage benchmark rows, their shared numerical-geometry covariance, and aggregate recovery accounting. It selects the declared first-stage geometry without replacing the primary/supplemental inputs. Those first two commands use saved aggregate data; reproducing their sampling summaries requires regenerating the synthetic replications.

The third command reconstructs the bounded fitted-WDSM outer sampling summaries directly from 60,000 preserved synthetic estimator records and 4,000 count-set diagnostics. It checks all 72 summary cells, 32 paired comparisons, 16 count-comparison cells and complete memberships, without generating data, refitting models or changing RNG state. It does not rederive individual count-draw moments or repeat the raw archive audit. See [bounded_wdsm/README.md](bounded_wdsm/README.md) for the qualified results, paired-subset requirement and missing outer-monitor receipt. Its separate artifact manifest verifies the exported files and pinned helper.

Optionally save the reconstructed geometry object and repeated joins to a **new** directory:

```sh
Rscript --vanilla results/reproduce_benchmarks.R reproduced-aggregates
Rscript --vanilla simulations/plot_summary.R results/primary/summary_joined.csv reproduced-primary-figures
Rscript --vanilla simulations/plot_extension_summary.R results/supplement/summary_joined.csv reproduced-supplement-figures
Rscript --vanilla results/reproduce_first_stage.R . reproduced-first-stage-aggregates
Rscript --vanilla simulations/plot_extension_summary.R results/first_stage/summary_joined.csv reproduced-first-stage-figures
```

The plotting commands consume saved aggregates and write PDF/PNG outputs. On an R installation without a PNG graphics device, the extension plotter accepts `--no-png`. Review rendered pages before using them in a publication. These commands do not overwrite existing output directories.

The reconstructed `geometry.rds` contains the saved fields required by the benchmark interfaces, with an explicit CSV-reconstruction marker. It is not claimed to be byte-identical to the original complete integration object. Its estimates, covariance and precision flags are retained at the CSV precision used in validation.

## Package and selected synthetic replications

Install the package with `R CMD INSTALL .`. For a full source check, first install all suggested dependencies listed in `DESCRIPTION` (currently `testthat`, `knitr` and `rmarkdown`), then run `R CMD check --no-manual .`; checking the PDF manual separately requires a TeX toolchain. The package's declared R dependencies are listed in `DESCRIPTION`.

To exercise the exact declared DGP/seed indexing on the current source without running a full production study, use new output names:

```sh
Rscript --vanilla simulations/run.R results/primary/config.csv strong_d1_m1_n200 1 2 fresh-run/primary.csv .
Rscript --vanilla simulations/supplement_run.R results/supplement/config.csv spline_d1_n800 1 2 fresh-run/supplement.csv .
Rscript --vanilla simulations/first_stage_run.R results/first_stage/config.csv same_sample_d2_m3_n200 1 2 fresh-run/first-stage.csv .
```

The final `.` explicitly sources this checkout's R files. Omitting it selects an installed `wdsmatch`, which may be a different version. These commands generate two datasets in each named scenario, retaining every requested estimator record and failure. Their outputs are smoke checks, not estimates of coverage precision. The scenario seed and replication ID determine the simulation stream independently of task chunking; paired estimators share the generated dataset. A complete first-stage regeneration uses IDs 1–1,000 in all 39 rows of `results/first_stage/config.csv`; the current Gaussian assertion includes the documented numerical repair.

For a full regeneration, use all IDs `1:1000` in every row of `results/primary/config.csv` and `1:200` in every row of `results/supplement/config.csv`, partitioned into disjoint bounded chunks. Source hashes and configuration must remain fixed across chunks. The supplied `simulations/production.R` and `simulations/extension_production.R` planners also preserve task manifests and refuse existing outputs. They require a newly frozen source manifest through `WDSM_RUN_MANIFEST`; **do not relabel a new run with a historical snapshot manifest**. The extension supplemental planner first imports newly generated IDs 1–5 from `simulations/supplement_pilot.R`, then schedules IDs 6–200. Its CSV configuration is supplied here as `results/supplement/config.csv`. Scheduler account, partition, wall time and concurrency are environment-specific and must be set explicitly; no script submits jobs automatically. Historical primary timeouts show that a fixed 30-minute allowance was insufficient for some original task sizes.

Create a manifest for a **new current-checkout run** before preparing its task plan. This optional production setup step uses Python 3's standard library; Python is not needed for the saved-aggregate R check. The example records relative source/configuration paths with SHA256 fingerprints and refuses to overwrite an existing manifest:

~~~sh
python3 - <<'PYMANIFEST'
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json

manifest_path = Path("fresh-source-manifest.json")
if manifest_path.exists():
    raise SystemExit("Use a new manifest name; preserve the existing run.")
required = [
    Path("DESCRIPTION"), Path("NAMESPACE"),
    Path("results/primary/config.csv"),
    Path("results/supplement/config.csv"),
    Path("results/first_stage/config.csv"),
    Path("results/geometry/config.csv"),
]
if any(not path.is_file() for path in required):
    raise SystemExit("Run this example from the source-checkout root.")
files = set(required)
for directory in ("R", "simulations"):
    files.update(path for path in Path(directory).rglob("*")
                 if path.is_file() and path.suffix in (".R", ".csv", ".md"))
record = {
    "schema_version": 1,
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "origin": "New run of the current source checkout",
    "hash_algorithm": "SHA256",
    "files": {path.as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
              for path in sorted(files, key=lambda path: path.as_posix())},
}
with manifest_path.open("x", encoding="utf-8") as stream:
    json.dump(record, stream, indent=2)
    stream.write("\n")
print("Created", manifest_path, "with", len(record["files"]), "file fingerprints.")
PYMANIFEST
export WDSM_RUN_MANIFEST="$PWD/fresh-source-manifest.json"
~~~

Keep this manifest and the source checkout unchanged throughout that run. The production runners record the manifest's digest and separately enforce their frozen source/configuration maps; they do not parse this JSON. Before starting or resuming tasks, verify the manifest against the checkout:

~~~sh
python3 - <<'PYMANIFEST'
from pathlib import Path
import hashlib
import json
import os

record = json.loads(Path(os.environ["WDSM_RUN_MANIFEST"]).read_text())
for name, expected in record["files"].items():
    path = Path(name)
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise SystemExit("Source/configuration differs: " + name)
print("Manifest agrees with the source checkout.")
PYMANIFEST
~~~

Then use the planner commands in [simulations/EXTENSION_PRODUCTION.md](../simulations/EXTENSION_PRODUCTION.md), retaining the same environment variable for each task. This setup creates no simulation, submits no job and does not change the historical results.

After complete regeneration, use `simulations/summarize.R` for complete primary chunks, `simulations/supplement_summarize.R` for the complete eight-scenario supplemental grid, and `simulations/first_stage_summarize.R` for first-stage chunks; their CLI headers give the arguments. They validate record grids and metadata, so partial smoke outputs cannot be passed off as completed studies. To join newly summarized outputs to the exported geometry, use `simulations/join_benchmarks.R` or `simulations/supplement_join.R` with `reproduced-aggregates/geometry.rds`.

For the first-stage join, extract the reconstructed geometry list from its covariance wrapper, then join a newly generated summary:

```sh
Rscript --vanilla -e 'x <- readRDS("reproduced-first-stage-aggregates/first_stage_geometry.rds"); stopifnot(!file.exists("reproduced-first-stage-aggregates/geometry-only.rds")); saveRDS(x$geometry, "reproduced-first-stage-aggregates/geometry-only.rds")'
Rscript --vanilla simulations/first_stage_join.R fresh-run/first-stage-summary.csv fresh-run/first-stage-joined.csv reproduced-first-stage-aggregates/geometry-only.rds
```

To repeat the numerical integrations themselves, the independent seeds and draw counts are fixed in `results/geometry/config.csv`:

```sh
Rscript --vanilla simulations/geometry_boundary_constants.R results/geometry/config.csv fresh-geometry .
```

This is a separate computation, not necessary for the saved-aggregate check. The saved geometry was generated with R 4.4.2 and `Mersenne-Twister/Inversion/Rejection`; simulations use `L'Ecuyer-CMRG/Inversion/Rejection`. Cross-platform floating-point arithmetic can change final digits. The numerical precision qualifications in [RESULTS.md](RESULTS.md) apply even when random streams reproduce exactly.

## Source-version boundary

`results/provenance/*_snapshot.json` records historical relative-file SHA256 mappings. Those fingerprints document which source versions produced the results; the older source archives and raw per-replication records for `primary/`, `supplement/` and `first_stage/` are not included in those aggregate-only exports. This exclusion does not describe the compact estimator records in the bounded, principal iid, original WDSM or Lenis archives. Later documentation, provenance handling or new interfaces can differ in the current checkout. The deterministic check establishes compatibility of current benchmark helpers with the saved aggregate calculations; it does not claim a new current-source run occurred or certify every historical source byte is present. A bit-for-bit historical execution requires its matching source archive and runtime. Independent reviewers checked the full retained records before this aggregate export was prepared.
