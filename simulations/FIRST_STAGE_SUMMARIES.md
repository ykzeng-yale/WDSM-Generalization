# First-stage summary and benchmark join

These postprocessing scripts use base R and do not refit, rematch, simulate, or integrate geometry. They are separate from the frozen first-stage generation/runner sources.

```
Rscript simulations/first_stage_summarize.R summary.csv batch1.csv batch2.csv
Rscript simulations/first_stage_join.R summary.csv joined.csv [geometry.rds ...]
Rscript validation/check_first_stage_summaries.R
```

The summary requires every batch's adjacent `.metadata.rds`, a completion marker, the exact declared six-record grid and all planned replication IDs for each scenario. It rejects duplicate records/files, inconsistent scenario/configuration declarations, mismatched code or design hashes, broken paired RNG/first-stage identities, or a changed adjusted-versus-naive point estimate. It records input and metadata hashes plus original manifests in its own companion metadata. Hash equality certifies source-version consistency across the supplied batches; it does not re-execute historical code. The CLI refuses existing output or metadata.

Sourceable interfaces are `wm_first_stage_read_batches(inputs, require_complete=TRUE)`, `wm_first_stage_validate_records(records)`, `wm_first_stage_summarize(records)`, and `wm_first_stage_join(summary, geometry=list(), geometry_relative_goal=.01)`. The read function's explicit `require_complete=FALSE` permits a deliberately partial source-level diagnostic; the CLI uses complete scenarios only.

Each replication yields six dependent records: two estimands by three comparisons. A comparison's requested/success/failed counts and all MCSE calculations count independent replications, never six times that count. Numerical summaries condition on successful fits and retain failures. The summary reports bias and MCSE, empirical and mean fitted root-n variance and their MCSE, a paired delta-method MCSE for their ratio, coverage and binomial MCSE, exact pointwise 95% Clopper–Pearson coverage limits, and mean interval width. Exact binomial limits remain informative when all 40 pilot intervals cover; a zero plug-in MCSE does not establish certain coverage. These intervals are not adjusted for examining multiple scenarios.

Adjusted-minus-naive reported root-n variance is computed on joint successes after verifying point identity. This paired summary is repeated on the three comparison rows for a given scenario/estimand; those copies are one result, not independent observations. First-stage fitting failures may leave known-first-stage fits intact. The empirical variance and ratio MCSE are asymptotic approximations, especially imprecise at pilot replication counts.

The join accepts named `M%d_d%d` geometry lists or individual `wm_geometry`/`wm_geometry_alpha` objects. Dimension one uses the exact alpha; other dimensions preserve missing, failed or unresolved integration states. It keeps all summary fields and failure counts. Actual sampling variance and the deliberately naive reported-variance limit are separate columns. Geometry uncertainty is propagated by the affine alpha coefficient and remains separate from simulation MCSE and deterministic quadrature diagnostics. Combined discrepancy MCSE assumes independent geometric simulation. Descriptive z scores are omitted if the geometric precision goal is unresolved, or if a numerical estimate violates the known Jensen lower bound. They are finite-simulation discrepancies, not tests of an asymptotic theorem. Shared geometric draws and paired data induce dependence across table rows.

The expected adjusted-minus-naive reported variance change is available even without geometry because alpha cancels. Quadrature's reported absolute errors are retained as diagnostics; they are not a certified propagated error bound for the matrix-inverse variance reduction.
