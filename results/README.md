# Published synthetic simulation results

This directory contains the completed synthetic studies, including failures and adverse coverage. Participant data and manuscripts are excluded. The existing historical result directories remain available.

| Study | Code | Results |
|---|---|---|
| Principal iid Weighted Matching | [principal_iid](../simulations/principal_iid/) | [principal_iid](principal_iid/): 4,000 samples, 24,000 method records, 24 point and 48 inference summaries |
| Original WDSM survey benchmark | [original_wdsm](../simulations/original_wdsm/) | [original_wdsm](original_wdsm/): 608,000 compact method records, 608 summary cells; [saved inference comparison](original_wdsm/saved_inference/): 16,000 paired records, 32 inference summaries |
| Lenis comparison | [lenis](../simulations/lenis/) | [lenis](lenis/): 1,800,000 compact method records, 72 settings, 25 methods, 1,800 summary rows and figures |
| Lenis no-response saved inference comparison | [saved-array adapter](../simulations/lenis/lenis_no_saved_adapter.R) | [saved_inference](lenis/saved_inference/): 54,000 method records, 108 inference and 54 paired summaries |

Gzipped CSVs contain compact per-replicate synthetic results; CSVs preserve the accepted numerical values. SHA-256 hashes and row counts are in [publication_results_manifest.json](publication_results_manifest.json). Figures are synthetic simulation figures, not manuscript PDFs.

The [principal archive reaggregation command](principal_iid/README.md) regenerates
all eight principal summary tables directly from these public records and fixed
target metadata, with no private execution files or new model fits.
The [original WDSM](original_wdsm/README.md) and [Lenis](lenis/README.md) record
archives have separate manifests and use the shared
[compact archive reaggregation command](../simulations/reaggregate_compact_archive.R)
to recover their complete performance metrics, paired MCSE and failure audits.
These commands need R packages `jsonlite` and `digest`; they do not fit models.

| Reproduction path | Saved-result reconstruction checked | Producer execution scope |
|---|---|---|
| [Principal iid](principal_iid/README.md) | 24,000 records; all eight summary tables | Saved-record validation does not imply a fresh 4,000-sample producer run |
| [Original WDSM](original_wdsm/README.md) | 608,000 records; 608 metric/audit/diagnostic rows and main-comparison metric fields | Four-setting replicate-1 integration; qualified historical count comparison, not a new full study |
| [Lenis](lenis/README.md) | 1,800,000 corrected records; 1,800 metrics/audit rows byte-identical, all 48 failures retained | 18 population/effect combinations generated; composed one-replicate integration and PS repair, no complete postpatch producer rerun |

The linked READMEs give the exact commands, dependencies, comparison tolerances
and retained failures. [REPRODUCE.md](REPRODUCE.md) covers older historical
primary, supplemental, first-stage and bounded archives.

The principal iid study concerns the paper's iid theory. The original survey benchmark and Lenis studies are external empirical comparisons; they do not establish cluster/stratum inference or estimated-weight theory. Lenis response-weight fits are held fixed in that external comparison. Analytic and fixed-reuse refit variances are distinct procedures and need not agree. Coverage uses all requested replicates; unavailable intervals remain failures rather than deleted observations. Low coverage is retained.

The code and saved results are published together. This is not a claim that every public entrypoint has been freshly rerun end to end. Historical source-bound saved-array calculations also consumed large local simulation objects and authentication receipts; these private execution artifacts are not published. The Lenis saved adapter's statistical calculations are unchanged; its comment header and raw-PS section-reference label have been updated. Its authentication layer expects caller-supplied input contracts and saved arrays. The original saved-inference summary script accepts a companion-record CSV and optional original point-summary CSV. Existing production generators and simulation engines are in the code directories above.
