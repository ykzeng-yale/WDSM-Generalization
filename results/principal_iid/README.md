Completed principal iid study. See [results index](../README.md) for counts, code, reporting definitions and reproducibility scope. Formal plots are in `figures/`; compressed per-replicate results are in `bound_records.csv.gz`.

Reproduce the reported summary arithmetic from this public archive, using
Python 3 and the standard library:

```sh
python3 -I results/reproduce_principal_iid.py /your/new/summary-directory
```

Run from the release root, or use an absolute script path. The output directory
must be new and its parent must exist. The script verifies published file
hashes, all 4,000 sample keys and 24,000 method records, complete nuisance
dimensions, paired samples/counts, availability and count-failure positions.
It calls the unchanged `simulations/principal_iid/formal_metrics.py`, writes
all eight summary tables, and compares every field with the published tables.
The receipt records numerical discrepancies; the executed archive reproduction
matched all reported numbers exactly.

Fixed target values and their arithmetic sensitivity bounds come from the
published target metadata. This command does not recompute targets, observations,
model fits or inference arrays. For those operations, use the separate
[truth and case entrypoints](../../simulations/principal_iid/README.md).
No private controller, authentication file or original RDS archive is needed
to reproduce these summary calculations.
