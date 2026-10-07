# Software validation

## Current package

The current source is the independent `WeightedMatching` 0.1.0 package, with
25 public `wm_*` statistical interfaces. The original `wdsmatch` 0.2.1 remains
separately installed; its `wdsmatchATE()`/`wdsmatchATT()` functions and legacy
S3 methods are not current `WeightedMatching` exports.

The local macOS/R 4.4.2 `WeightedMatching` 0.1.0 archive passed
`R CMD check --no-manual` with zero errors, warnings or notes: 5,172 test
expectations passed, with no failures, warnings or skips. The source-bound
receipt pins 117 package inputs and the checked archive, and records successful
help/example checks and vignette execution, installation and rebuild.
The package name, version, archive and input identities are recorded in
[software_validation.json](software_validation.json). A check of the
earlier `wdsmatch`-named development archive is not a check of the new package's
installation, namespace separation or installed help. Unit tests, independent
reference comparisons and public saved-result replay remain separate checks.

Identity regression checks retain all 209 common function bodies and formals,
12 exact supplied-map/seeded-contribution cases and two fitted-WDSM point cases
against the earlier development archive. They establish finite software parity,
not a new simulation performance result. The original CRAN `wdsmatch` 0.2.1 was
restored separately; both package load orders were checked without export
collision or replacement of its legacy S3 methods, and installing the new
package left the original installed files unchanged.

Linux/macOS/Windows CI is defined in `.github/workflows/R-CMD-check.yaml`;
consult observed Actions results rather than treating its configuration as a
completed test. Use `help(package="WeightedMatching")` for current interfaces.

## Historical check before package separation

The following describes the earlier `wdsmatch`-named 0.3.1 development source,
before the independent package was established. Its check counts and source
identities retain that historical scope; they do not certify the current
`WeightedMatching` archive.

Version 0.3.1 added compact WM print/summary methods, an explicit
[user-options guide](user_options.md), and installed tests for the three scalar
interfaces previously covered only by standalone validation scripts. Existing
donor-normalized estimation, graph selection, analytic variance and replication
routines are unchanged. The original synthetic scalar fixture is reused byte
for byte; publication simulations and participant analyses were not rerun.

The earlier local macOS/R 4.4.2 source archive passed `R CMD check --no-manual`
with zero errors, warnings and notes: 5,165 passed expectations, no failures,
warnings or skips. The 119 package-input identities and source/archive binding
were recorded in its historical receipt. Help usage/defaults, examples and the executed,
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
