# Complete original WDSM source-setting comparison

These unchanged compact synthetic reports retain the full608-cell plan with
1,000 requested replications per cell:4,000 selected samples,32,000 sample/case
evaluations,608,000 estimator records and96,000 source-reduction rows. They use
the source design in [simulations/original_wdsm](../../simulations/original_wdsm/README.md),
including both selection designs, overlap settings, source nuisance recipes,
PATE/PATT and M=1/3/5. CorCor/CorMis/MisCor/MisMis describe the original two-score
model specifications; they are not general WM correctness or theorem labels.

`metrics.csv` contains every cell and its original reporting quantities;
`audit.csv` retains requested/received, point and interval denominators;
`point_success_diagnostics.csv` keeps availability/status diagnostics.
`method_failures.csv` and `reduction_mismatches.csv` retain the complete failure
schemas (zero data rows in this study). `main_comparison.csv` is the predefined
64-cell CorCor M=3 comparison, and `all_cells.csv` retains all608 display-labelled
cells. These are aggregate estimator reports, not participant samples, count
matrices, fitted objects, manuscript tables or private execution records.

Reporting follows the existing signed-relative-bias denominator abs(truth),
same-design/overlap/estimand CorCor PSM_M1 empirical-variance reference and paired
relative-efficiency MCSE. Interval availability, full-requested operational
coverage and available-interval conditional coverage stay distinct. Coverage
MCSE and exact pointwise intervals use the original requested denominator.
Known-map analytic variance, fitted variance and empirical prediction-refit
variance are separate quantities. Native PSM/PGM/DSM use their source B-1/1.96
intervals, WM/WDSM use full B=200/qnorm intervals, and source SWPSM retains its
subclass construction. The 96 analytic WM_DSM rows use an empirical
fitted-stack/full-X approximation for two-dimensional double-score matching.
Its population, fitted-model, regular-root and graph premises remain unverified;
dependent source selection and mixed support do not certify the
independent-observation theory or complex-design variance.

All source WDSM/WM-DSM paired graph, point, draw, variance and interval
comparisons meet the unchanged scaled1e-10 criterion. Native versus frozen
saved-report values have maximum numeric difference 2.66454e-14; native versus
public saved-report values have maximum difference 2.19380e-13. Both comparisons
pass their declared floating-point tolerance, without byte-identical R
serialization. These are
implementation/reporting checks, not mathematical or causal certification.

Finite-sample departures remain visible. CorCor WDSM/WM-DSM M=3 coverage ranges
from92.2% to94.8%; its relative efficiency ranges1.58–10.97 against the declared
reference. The prospective good-overlap PATE cell covers92.2%, with MCSE0.85
percentage points and variance calibration0.89446. Full-X M=3 coverage is
93.2–94.8%. Several unweighted comparators retain substantial bias despite
high coverage. Pointwise Monte Carlo diagnostics flag266 undercoverage cells,
36 overcoverage cells and405 finite-sample bias signals; these are neither
familywise tests nor correctness decisions. No blanket method superiority or
forced93–97% coverage criterion is imposed.

The [record archive](records/manifest.json) contains all 608,000 per-replication
synthetic estimator records in eight gzipped CSV parts (20,697,223 compressed
bytes in total). Each decoded part is at most 32 MiB. The 23-column schema retains
the original estimates, variances, intervals, status/error fields, interval
scope and actual dataset pairing IDs. Numbers and row order are preserved as
original CSV bytes, with the same header repeated in each part. The manifest
pins every compressed and decoded part and the reconstructed source CSV.
Participant samples, fitted models, matching graphs and bootstrap count arrays
are excluded. Historical aggregate hashes and source links remain in
`provenance.json` and `artifact_hashes.csv`; the new archive has its own manifest.

From the repository root, with R packages `jsonlite` and `digest` installed:

```sh
Rscript simulations/reaggregate_compact_archive.R original_wdsm /tmp/wm-original-summary
```

The output directory must be new, with an existing parent outside the checkout.
The wrapper authenticates the full original CSV stream, reads explicit column
types, and calls the unchanged `ows_report` and `wm_paired_comparison_metrics`.
It writes `metrics.csv`, `audit.csv`, `point_success_diagnostics.csv`,
`failure_records.csv`, `report.rds` and `aggregation_receipt.json`. It uses only
the public compact records, with no model fits or private execution files.
Recomputed RB, empirical variance, RE, paired MCSE, coverage and denominator
audits can be compared with the existing summaries. Original CSV decimal
serialization permits small roundoff differences from native RDS-based
summaries; the previous comparison used a declared 1e-12 numeric tolerance,
with status, missingness and denominators exact. Display labels in
`all_cells.csv` add no statistical calculations. The 64 main rows select CorCor
and methods WDSM_M3, WM_PS_M3, WM_DSM_M3, WM_X6_M3, PSM_M1, PGM_M1, DSM_M1 and
SWPSM_full across the eight design/overlap/estimand combinations.

The public archive command was executed once on all 608,000 records. All 608
metric, audit and diagnostic rows and the 49 common metric fields in
`all_cells.csv` and the 64-row `main_comparison.csv` agreed with the accepted
outputs. The largest absolute numeric difference was 2.203e-13; the audit CSV
was byte-identical. Comparison allowed absolute differences no greater than
max(1e-12, 1e-10 times the larger absolute value), with text, missingness and
denominators preserved. This validates summary reconstruction from the public
archive, separately from the producer comparison below.

Reproduction from new synthetic samples uses the separate generation/run
commands and external source dependencies. The public producer was executed
for replicate 1 in all four design/overlap settings: 608 method records and 32
source cases. Of these, 456 historical input/count-paired records met the
comparison tolerance. The 152 PoorOverlap retrospective records retained
point/raw/target compatibility; 56 intervals independent of the supplied WM
count matrix also matched. The other 96 historical shared-count intervals are
unpaired descriptive comparisons because the last two count rows differed.
The strict all 608 comparison failure is retained. All 96 within-run WDSM/WM-DSM
reductions passed. This is qualified software integration evidence, not a
fresh 1,000-replication study or complete historical graph/draw identity; see
the [producer README](../../simulations/original_wdsm/README.md).
