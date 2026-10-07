# Software validation

Version 0.3.1 adds compact WM print/summary methods, an explicit
[user-options guide](user_options.md), and installed tests for the three scalar
interfaces previously covered only by standalone validation scripts. Existing
donor-normalized estimation, graph selection, analytic variance and replication
routines are unchanged. The original synthetic scalar fixture is reused byte
for byte; publication simulations and participant analyses were not rerun.

The final source-bound package build/check receipt is recorded in
`software_validation.json`. Unit tests, independent reference comparisons and
public saved-result replay are separate checks. Linux/macOS/Windows CI is
defined in `.github/workflows/R-CMD-check.yaml`; consult the observed Actions
status rather than treating workflow configuration as a completed test.

The current local macOS/R 4.4.2 source archive passes `R CMD check --no-manual`
with zero errors, warnings and notes: 5,165 passed expectations, no failures,
warnings or skips. The 119 package-input identities and source/archive binding
are recorded in the receipt. Help usage/defaults, examples and the executed,
installed and rebuilt HTML vignette passed. The additional 141 expectations
cover retained summary statistics/intervals and scalar option/input contracts.
All 27 exports have installed help and test-call bindings; this is not exhaustive
coverage of every option combination. The earlier 0.3.0 receipt remains bound to
its own tag and source; its reference comparisons/replay were not newly rerun.

Independent reference checks cover PATE/PATT, d=1/2/3/6, M=1/3/7, nonconstant and
unit weights, distinct/unequal arm maps, direct weighted imputation and full
quadratic WLS. External Matching/MatchIt point comparisons align their actual
distance scaling, replacement, donor count, ties and target before comparison.
Weighted WDSM source comparisons include known/fitted scores and both original
PS-fitting conventions. These establish finite-array point reductions, not an
asymptotic sampling theorem.

New integration tests compare automatic complete same-count fitted-wrapper
refits against an independently written PS/PG/pooling/quadratic prediction
recipe and the original PATE/PATT counted ratio. Contribution replication uses
the retained complete covariance. The tests retain failed support/root requests
and reject inconsistent count labels/options. No failed draw is deleted or
silently replaced.

The principal public record replay reconstructs summary arithmetic without
fitting or drawing data. Fresh production uses the frozen source checkout
specified in its README. Later API/diagnostic changes do not update historical
study hashes or turn old results into evidence for a new algorithm.

These checks do not certify causal identification, all geometry/centering/root
premises, eventual numerical availability, dependent-survey or estimated-weight
inference, original refit-law agreement in a new application, or nominal
coverage. Detailed method contracts are in R help and the theory/method guide.
Historical check records in SOURCE_PROVENANCE.json retain their original pins.
