# Published synthetic simulation results

This directory contains the completed synthetic studies, including failures and adverse coverage. Participant data and manuscripts are excluded. The existing historical result directories remain available.

| Study | Code | Results |
|---|---|---|
| Principal iid Weighted Matching | [principal_iid](../simulations/principal_iid/) | [principal_iid](principal_iid/): 4,000 samples, 24,000 method records, 24 point and 48 inference summaries |
| Original WDSM survey benchmark | [original_wdsm](../simulations/original_wdsm/) | [original_wdsm](original_wdsm/) and [saved inference comparison](original_wdsm/saved_inference/): 16,000 paired records, 32 inference summaries |
| Lenis comparison | [lenis](../simulations/lenis/) | [lenis](lenis/): 72 settings, 25 methods, 1,800 summary rows and figures |
| Lenis no-response saved inference comparison | [saved-array adapter](../simulations/lenis/lenis_no_saved_adapter.R) | [saved_inference](lenis/saved_inference/): 54,000 method records, 108 inference and 54 paired summaries |

Gzipped CSVs contain compact per-replicate synthetic results; CSVs preserve the accepted numerical values. SHA-256 hashes and row counts are in [publication_results_manifest.json](publication_results_manifest.json). Figures are synthetic simulation figures, not manuscript PDFs.

The [principal archive reaggregation command](principal_iid/README.md) regenerates
all eight principal summary tables directly from these public records and fixed
target metadata, with no private execution files or new model fits.

The principal iid study concerns the paper's iid theory. The original survey benchmark and Lenis studies are external empirical comparisons; they do not establish cluster/stratum inference or estimated-weight theory. Lenis response-weight fits are held fixed in that external comparison. Analytic and fixed-reuse refit variances are distinct procedures and need not agree. Coverage uses all requested replicates; unavailable intervals remain failures rather than deleted observations. Low coverage is retained.

The code and saved results are published together. This is not a claim that every public entrypoint has been freshly rerun end to end. Historical source-bound saved-array calculations also consumed large local simulation objects and authentication receipts; these private execution artifacts are not published. The Lenis saved adapter's statistical functions are archived unchanged, with only its outdated comment header replaced. Its authentication layer expects caller-supplied input contracts and saved arrays. The original saved-inference summary script accepts a companion-record CSV and optional original point-summary CSV. Existing production generators and simulation engines are in the code directories above.
