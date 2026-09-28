# Joining and plotting the prespecified benchmarks

Run the deterministic checks before using the comparison stage:

```sh
Rscript validation/check_benchmarks.R
Rscript validation/check_join_benchmarks.R
```

From the package source directory, join a completed simulation summary with
independently generated geometric constants:

```sh
Rscript simulations/join_benchmarks.R results/summary.csv results/compared.csv geometry/geometry.rds
```

The final argument can be omitted for one-dimensional results, or replaced by
one or more individual `M{M}_d{d}.rds` files. Duplicate geometric keys are rejected;
do not provide both an aggregate and its constituent files. Both output paths
must be new. Input/code hashes and the independent-stream contract are saved in
the output's `.metadata.rds` sidecar. The geometry RDS contains covariance of the
estimated beta vector, already divided by the integration draw count. It is not
the covariance of a single integrand draw.

The join preserves all input rows, columns, and requested/completed/failed counts.
It adds the exact probability limit and target bias in every supported primary
scenario, plus corrected-estimator variance benchmarks for oracle and polynomial
rows. Raw estimators have no variance benchmark. Original-rule weak-centering
rows retain their known probability limits but have no established variance
benchmark; their reported intervals and contribution variances remain labeled
assumption-failure diagnostics.

One-dimensional geometry is computed from the exact closed form with zero
geometric Monte Carlo error. In higher dimensions the full beta covariance is
propagated into each benchmark. Missing geometry, failed runs, unresolved
zero-hit components, numerical zero empirical component variances, and a missing
covariance matrix are separately retained or flagged. A displayed numerical zero
from an unresolved integration is not treated as known zero error: its usable
`benchmark_geometry_mcse` is `NA`, while `geometry_reported_mcse` retains the
original numerical value.

`empirical_root_n_variance` is evaluation n times the empirical sampling variance;
its MCSE receives the same scaling. `mean_root_n_variance` and its MCSE are kept
from the original summary. Differences from the explicit benchmark use
`sqrt(simulation_MCSE^2 + geometry_MCSE^2)`, requiring independent simulation and
integration streams. Standardized discrepancies are provided only when geometry
precision is resolved and its empirical MCSE is at most 1% of the positive
benchmark. The sourceable function permits a different explicit
`geometry_relative_goal`, which must be recorded as an analysis revision.

These are approximate Monte Carlo diagnostics, not hypothesis tests or a
coverage certification. Finite-n variances can differ from an asymptotic limit.
The geometric MCSE is itself empirical and is not a rigorous numerical error
bound. The shared geometric benchmark also induces dependence among different
comparison rows; do not treat their standardized discrepancies as independent.

Generate scientific figures from the joined aggregate output with:

```sh
Rscript simulations/plot_summary.R results/compared.csv results/figures
```

The output directory must be new. The script writes a multipage vector PDF and
300-dpi PNG pages, grouped by design, estimand, dimension and donor count. Each
page shows reported-to-empirical variance ratios, coverage and target bias against
sample size. Error bars are approximate estimate ± 1.96 Monte Carlo SE; ratios
use the paired delta-method SE from the summarizer. Raw results appear only in
the bias panel. Original weak-centering intervals are explicitly diagnostic, and
the original-rule probability-limit bias is drawn separately. Failure counts are
record counts across the methods shown, not independent replication counts.

Inspect the rendered figures before publication. This script creates no new
simulated observations and does not change the analysis results.
