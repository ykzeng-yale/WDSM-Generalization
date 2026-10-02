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

Hashes and source links are in `provenance.json` and `artifact_hashes.csv`.
The code reporting entry accepts analyst-supplied estimator-record CSVs and
reconstructs the full plan without fits. This compact export does not include
those608,000 per-replication records, so its CSVs alone cannot recompute empirical
variances or paired MCSE from individual replications. Reproduction from new
synthetic samples uses the explicit generation/run commands and their external
source dependencies; the public direct-run branch has not been executed here.
