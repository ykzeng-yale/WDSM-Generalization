# Simulation workflows

The current paper's principal comparison and designated predecessor benchmarks
have separate reproducible entry points:

| Experiment | Generate data and fit methods | Public record replay |
| --- | --- | --- |
| Principal iid PS / distinct-arm DSM / full-X | [principal_iid](principal_iid/README.md), using its frozen producer checkout | [Principal results](../results/principal_iid/README.md) |
| Original WDSM survey designs | [original_wdsm](original_wdsm/README.md) | [WDSM results](../results/original_wdsm/README.md) |
| Lenis survey/response-weight benchmark | [lenis](lenis/README.md) | [Lenis results](../results/lenis/README.md) |

Saved-record replay verifies reporting arithmetic. Fresh producer execution
requires the stated source versions, external comparator code/dependencies and
computing resources. The principal producer's source hashes deliberately remain
frozen; it is run in a dedicated checkout rather than against a changed package.
Do not reuse a completed output directory or replace failed datasets/draws.

Other scripts in this directory support earlier pilots, component/reference
calculations or historical studies indexed in [results/REPRODUCE.md](../results/REPRODUCE.md).
They are retained for traceability and are not additional principal studies of
the current paper. No cluster/stratum or estimated-weight inference follows
from running a weighted benchmark.

Current unpinned WM entry points use the separate `WeightedMatching` package.
Original WDSM comparison functions remain sourced from their pinned upstream
code; they are not replaced by a renamed package. Recorded results and source
receipts retain the package names and versions that actually produced them.

The principal iid producer and older bounded-calibration runners require their
frozen source checkouts and historical package artifacts. In particular,
`run_bounded_wdsm_calibration.R`, `prepare_bounded_wdsm_replay.R`,
`fitted_wm_components.R`, `wdsm_fitted_pipeline_reference.R` and
`bounded_wdsm_calibration.R` retain their original `wdsmatch` namespace and
checksum admissions. Installing `WeightedMatching` does not satisfy those
historical admissions. Do not change their pins or regenerate completed studies
to accommodate the package separation. For the principal study, use its
[documented frozen checkout](principal_iid/README.md); use the exact source and
installed-artifact receipts named by each bounded-calibration config for that
historical workflow. Current saved-fit callers of
[`wm_saved_full_x_stack.R`](wm_saved_full_x_stack.md) pass
`wm = asNamespace("WeightedMatching")` explicitly, because its frozen default
remains historical.
