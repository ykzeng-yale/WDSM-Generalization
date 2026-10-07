# Original WDSM comparison source

This directory exposes the original population generators, verbatim sampler and
existing statistical comparison functions. Their transferred bytes and hashes
are in `provenance.json`. The small entry points parameterize user paths; they do
not change matching, model recipes, correction, counts, comparator inference or
failure retention. The first direct generation and replicate1 run checks are
described below; source-transfer, parse and preflight checks are separate evidence.

The complete design is in `design.json`: Good/Poor overlap, retrospective and
prospective sampling, PATE/PATT, all four original PS/PG specifications,
M=1/3/5, B=200 and 1,000 replicates. Each selected sample runs the original WDSM
and WM PS/DSM/X6 methods at all three M values, separate WM DSM analytic rows,
native PSM/PGM/DSM at M=1 and source full-matching SWPSM. Correct PS/PG includes
X1:X2; the source misspecified recipes omit it. The actual native coordinates,
metric and inference are retained. WM PS uses the scenario full-X correction,
and WM X6 uses a complete quadratic correction. Failed prescribed draws are
retained and prevent a complete refit interval; no replacement samples/counts
are drawn. See the design for shared-count coupling and the separate divisors.

Install the package from this checkout into a separate existing R library.
Additional simulation dependencies are digest, jsonlite, sampling, Matching,
boot, MatchIt, optmatch, marginaleffects and sandwich. The original runtime's
reference versions are in `dependency_versions.json`.
The namespace guard requires the package version declared in this checkout
and all required study interfaces. A matching version alone is insufficient.
Supply the source checkout of [SW_DSM](https://github.com/ykzeng-yale/SW_DSM/tree/e96dcc82baa9f8f8e50c2c5d56dadaa756cf6731).
`upstream_source_pins.json` checks the 14 required statistical source files;
that dependency tree is not redistributed here. Its use follows its upstream
terms. The authored code in this repository uses the root GPL-3 license.

From the repository root, with existing input/output-parent directories:

```sh
R CMD INSTALL --library=/existing/wm-library .
Rscript simulations/original_wdsm/generate_populations.R preflight /existing/new-populations
Rscript simulations/original_wdsm/generate_populations.R generate /existing/new-populations
R_LIBS=/existing/wm-library Rscript simulations/original_wdsm/run.R preflight /existing/SW_DSM /existing/new-populations/GoodOverlap_population.rds /existing/new-output GoodOverlap retrospective 1,2
R_LIBS=/existing/wm-library Rscript simulations/original_wdsm/run.R run /existing/SW_DSM /existing/new-populations/GoodOverlap_population.rds /existing/new-output GoodOverlap retrospective 1,2
```

Generation is an explicit million-row-per-overlap workload. Preflight does not
evaluate generators, load population contents, draw counts or fit models and
does not create output. `run` is an explicit sequential scientific workload,
with the specified replicate IDs and full B=200 calculation. Each output path
must be new, with a writable existing parent, outside the package and supplied
source/population/record input directories. The same IDs/settings reproduce the source-defined experiment;
cross-platform bitwise historical-file identity is not asserted. Existing
populations may instead be supplied as explicit RDS paths with the original
column schema. Neither populations
nor raw sample/count archives are included in this export.

`report.R NEW_OUTPUT_DIR RECORDS_CSV [RECORDS_CSV ...]` applies the unchanged
original reporting contract to supplied per-sample estimator records. It keeps
the complete 608-cell/1,000-replicate plan and exposes missing/failed records.
Its CorCor PSM_M1 empirical-variance reference remains within the same
design/overlap/estimand; paired relative-efficiency/MCSE and actual source
interval coverage are kept distinct. Case outputs also retain exact graph,
point, draw and interval reduction checks from the transferred implementation.

The reporting entry point was exercised on 38 existing pilot records from one
replicate in two original setting cells. Its output exactly matched the
unchanged reporting function, with the complete 608-cell/1,000-request plan
retained as incomplete; duplicate inputs were rejected. This checks saved-record
reporting and missing denominators, not the full formal grid. The installed
package's 127 function bodies and formals matched current package source without
calling estimators. RNG state and kind were unchanged. The later complete formal study is now collected and independently assessed:
the unchanged public reporting entry agrees with the native full608-cell report
within floating-point tolerance, with1,000 requested replicates in every cell.
See [the compact full-study reporting export](../../results/original_wdsm/README.md).
The public generator was separately executed for both million-row overlap
populations; its RDS payloads are byte-identical to the accepted
source-reconstructed populations. Direct `run` checks evaluated replicate1 in
each of the four overlap/design settings, all eight source PS/PG/target cases,
M=1/3/5 and B=200. Of the608 keyed records,456 passed historical full-input/count
pairing and numeric/status comparisons. In PoorOverlap retrospective sampling,
all152 point records are compatible within the fixed scaled1e-10 rule, with
four supplied weights differing at machine precision. Another56 intervals are
compatible independently of the supplied WM counts:24 WM DSM analytic
approximations,24 native own-bootstrap intervals and8 SWPSM subclass sandwich
intervals. The remaining96 historical shared-count intervals are unpaired:
the200 count columns reverse the last two rows, and their actual interval
differences are retained rather than labeled reproduction. Their maximum
absolute and scaled descriptive gaps are0.0004934918079100292 and
0.00024013170437741632. The strict all608 comparison failed and remains part of
the validation history; qualified closure does not replace that failed gate.
All96 current within-run WDSM/WM DSM point, draw, variance and interval
reductions passed on their shared inputs/counts. The32 historical raw graph/draw
bundles were not locally compared; whole historical data/pairing-object byte
identity and an independently emitted original-run namespace-path manifest are
not asserted. The platform/compiler cause of the count reversal is not proved.
These first-workflow integration checks use one replicate per setting and do
not establish repeated-sampling bias, efficiency, coverage or analytic premises.

These are empirical source-design comparisons. Dependent cluster/stratum
designs and the separate analytic approximations do not establish main WM iid
theory, design-based variance, original-refit validity or scientific acceptance.
The complete formal report retains its actual finite-sample bias, variance
calibration and coverage departures; numerical source reduction and reporting
agreement do not certify all-cell nominal inference or the main theory's
independent-observation premises.
