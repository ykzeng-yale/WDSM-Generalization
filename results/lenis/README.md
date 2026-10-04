# Complete Lenis-inspired synthetic comparison

All 1,800 method-setting summaries and audit rows are retained, with 1,000
requested replications per cell. The plan crosses three source populations,
six effect multipliers, four response settings and 25 methods. Multipliers share
source samples; 72,000 setting/replicate evaluations are not 72,000 independent
samples. This extends the supplied benchmark and does not claim exact numerical
replication of the published paper. Plots are in `figures/`; see the
[code and source scope](../../simulations/lenis/README.md).

The [record archive](records/manifest.json) contains all 1,800,000 final corrected
estimator records in 19 gzipped CSV parts (88,681,452 compressed bytes in total).
Each decoded part is at most 32 MiB. These are the final saved records after the
72,000 weighted-quick records were corrected; the other 1,728,000 records were
retained. All 23 columns and their original CSV bytes/order are preserved,
including actual dataset pairing IDs, source coverage, warnings, error/status
fields and inference scope. The manifest pins compressed/decoded parts and
the reconstructed original CSV stream. No selected/respondent participant
records, fitted objects, graphs, count matrices or private run manifests are
included. The separate [no-response inference archive](saved_inference/) is
a narrower comparison and does not substitute for these complete records.

From the repository root, with R packages `jsonlite` and `digest` installed:

```sh
Rscript simulations/reaggregate_compact_archive.R lenis /tmp/wm-lenis-summary
```

The output directory must be new, with an existing parent outside the checkout.
This archive-only wrapper authenticates the reconstructed source CSV, uses
explicit column types, and calls the unchanged `lenis_reporting_plan` and
`wm_paired_comparison_metrics`. It needs no source simulation objects or private
execution contracts. Outputs are `metrics.csv`, `audit.csv`,
`point_success_diagnostics.csv`, `failure_records.csv`, `report.rds` and
`aggregation_receipt.json`. The existing metrics and audit CSVs derive from
these corrected records. Reaggregation is a numerical summary calculation,
not a new simulation, fit or inference validation.

The public archive command was executed once on these complete records. Its
1,800-row metrics and audit files are byte-identical to the accepted files.
Its point-success diagnostics also match the corrected source byte for byte.
All 48 failure records remain in the generated failure table. This validates
the saved-record reporting path and does not imply a new complete producer run.

The 48 failed native method records remain in the archive: 24 have the source
error `Stratum (7) has only one PSU at stage 1` and 24 have the corresponding
stratum 8 error. Their unavailable point/interval values are not replaced.
Unconditional point metrics remain withheld in the 24 affected cells;
successful-point summaries are emitted only as separately labeled diagnostics.
Coverage uses all 1,000 requested replications, counting unavailable intervals
as noncoverage; available-interval coverage is separately labeled. Native
Lenis intervals retain their strict endpoint rule; adapted/WM intervals use
the declared inclusive rule. Relative bias is 100 times signed bias divided by
the absolute finite-population PATT. RE divides the empirical variance of
Lenis_U_PS_U_OM_OW_U by the method's empirical variance; paired MCSE uses the
retained dataset identities.

Supplied fitted response weights remain frozen in this external empirical
comparison. The comparison does not add estimated-weight or cluster/stratum
variance theory to WM. Public generation was executed for all 18 population/effect
combinations; the composed one-replicate check and repaired PS path are described
in the producer README. There was no new complete 1,800-cell postpatch producer
rerun. Public compact-record reaggregation requires no such rerun.
