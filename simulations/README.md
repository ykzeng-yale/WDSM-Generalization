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
