# Weighted Matching

[![R CMD check](https://github.com/ykzeng-yale/WDSM-Generalization/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ykzeng-yale/WDSM-Generalization/actions/workflows/R-CMD-check.yaml)

R software for fixed-count nearest-neighbor matching with supplied positive
probability weights. The package is named **`wdsmatch`** and retains the original
`wdsmatchATE()` and `wdsmatchATT()` interfaces. The general WM interfaces support
PATE and PATT, a chosen positive integer `M`, supplied matching maps of any fixed
dimension, and different maps for the two treatment arms.

Weights can depend on treatment and covariates, `W=w_Z(X)`, with `w(X)` as a
special case. They define the target law `dP/dQ=W/E_Q(W)` for iid observed rows.
Identification of a causal effect and the assumptions for the selected inference
method must be justified separately. Estimated-weight uncertainty and dependent
cluster/stratum inference are outside the current framework.

## Install

This repository contains version 0.3.0 of `wdsmatch`. The GitHub source version
is separate from the previously published CRAN version. Installing it
updates an existing package of that name. R 3.6.0 or later is required; package
dependencies and optional dependencies are listed in [DESCRIPTION](DESCRIPTION).

```sh
git clone https://github.com/ykzeng-yale/WDSM-Generalization.git
cd WDSM-Generalization
R CMD INSTALL .
```

```r
library(wdsmatch)
```

The simulation, application and validation workflows use additional dependencies
listed in their own instructions. They belong to the source checkout and are
excluded from the installed R package.

## Quick start: a supplied matching map

This synthetic example illustrates matching on two covariates with a quadratic
outcome correction. Its donor-normalized weights and matching graph are retained
in the returned object.

```r
library(wdsmatch)
set.seed(928)
n <- 400
X <- matrix(runif(2 * n), n, 2)
Z <- rbinom(n, 1, plogis(-0.2 + 0.5 * X[, 1]))
W <- ifelse(Z == 1, 1 + X[, 1], 1 + 0.5 * X[, 2])
Y <- 1 + X[, 1] + X[, 2]^2 + Z + runif(n, -0.5, 0.5)

fit <- wm_fit(Y, Z, W, scores0 = X, M = 3,
              estimand = "PATE", regression = "polynomial", degree = 2)
fit$estimate
fit$se
wm_bootstrap(fit, B = 999, seed = 17)$conf.int
```

`wm_fit()` uses coordinates in the supplied scale; it does not standardize them.
Its ordinary row variance and contribution replication require the documented
residual-centering and feasible-correction conditions. This example is an
interface illustration, not a simulation coverage assessment.

Read the [step-by-step tutorial source](vignettes/weighted-matching.Rmd).
Direct `R CMD INSTALL .` does not build its HTML output. To include the installed
tutorial, first install `knitr` and `rmarkdown`, make Pandoc available, then run:

```sh
R CMD build .
R CMD INSTALL wdsmatch_0.3.0.tar.gz
```

After this built-archive installation, run
`vignette("weighted-matching", package="wdsmatch")`.

## Choose an estimator and inference route

| Task | Interface | Instructions |
| --- | --- | --- |
| Supplied scalar, distinct-arm or full-X coordinates and supplied means | `wm_match()` | `?wm_match` |
| Supplied coordinates with polynomial or tensor-spline correction | `wm_fit()` | `?wm_fit` |
| Fitted arm-specific double scores | `wm_wdsm_fit()` | `?wm_wdsm_fit` |
| Finite lists of fitted propensity/prognostic models | `wm_model_fit()` | `?wm_model_fit` |
| Original fixed-reuse count replication with optional prediction refits | `wm_bootstrap(..., method="fixed_reuse")`, `wm_bootstrap_refit()` | `?wm_bootstrap_refit` |
| Complete fitted-score contribution inference | `wm_fitted_inference()`, then `wm_bootstrap()` | `?wm_fitted_inference` |
| Raw fitted scalar PSM | `wm_scalar_logistic_match()` or the `scalar_psm` route | [Theory and method guide](docs/theory_method_guide.md) |

Matching selects exactly `M` distinct opposite-arm donors per query, permitting
reuse across queries. Under the original self-normalized rule, the donor fraction
is `W_j/sum(W_selected)`. For PATE, `scores0` and `scores1` can differ; PATT uses
only the control matching direction. Matching distance, scaling, tie handling
and correction models must align before comparing package outputs.

Fitted-score wrappers default to point-only estimation. `inference="full_x"`
requires correct full-X predictions and the applicable complete nuisance,
current-graph and numerical-root conditions. A correct PS by itself does not
establish that branch. `wm_model_fit()` fits score lists; it is not a generic
raw-covariate fitting wrapper. Raw scalar PSM and corrected full-X PS matching
are different statistics.

Contribution replication and original fixed-reuse refitting are different
operators. `wm_bootstrap(fitted_wrapper, method="contribution")` binds the
complete retained inference after an explicit full-X or post-fit inference
request. `method="fixed_reuse"` on a default quadratic fitted wrapper refits
its complete retained PS/PG/scaling/correction pipeline under each supplied
full-row count vector, while keeping the original graph and W. On a nested
`wm_match` object, `refit=NULL` instead holds predictions fixed. Other correction
methods require an explicit compatible callback. Weights remain fixed. See the
[theory and method guide](docs/theory_method_guide.md) for complete inputs,
covariance scopes, branch conditions and fixed-B versus growing-B statements.

The logistic fitters use a fixed numerical acceptance interval for fitted
probabilities. Under an unbounded Gaussian propensity index, eventual numerical
availability does not follow from exact-root asymptotic theory. Failed fits or
requested replicates remain failures: probabilities are not clipped and failed
draws are not silently replaced. Full details are in the guide and R help.

## Reproduce the published synthetic results

The [results index](results/README.md) provides the study designs, compact
per-replicate records, complete tables, failure accounting and hash manifests.

| Study | Simulation code | Saved results |
| --- | --- | --- |
| Principal iid PS / distinct-arm DSM / full-X comparison | [Recipe and producer commands](simulations/principal_iid/README.md) | [Tables and replay](results/principal_iid/README.md) |
| Original WDSM survey benchmark | [Source-design workflow](simulations/original_wdsm/README.md) | [Results and replay](results/original_wdsm/README.md) |
| Lenis survey/response-weight benchmark | [Source-design workflow](simulations/lenis/README.md) | [Results and replay](results/lenis/README.md) |

For example, regenerate all eight principal summary tables from the public
records using Python 3 and its standard library:

```sh
python3 -I results/reproduce_principal_iid.py /your/new/summary-directory
```

The output parent must exist and the output directory must be new. This verifies
saved records and reconstructs summary arithmetic; it does not rerun model fits.
The producer commands are separate. The principal producer has a frozen
53-file source inventory: use its documented frozen checkout
`f4acca0689b25625b0494c5daf965bc38e4a59bc` for a fresh declared-study run.
The current package source contains later diagnostic/documentation changes and
must not be relabelled as that frozen source. Saved-record replay works from the
current checkout and needs no historical model fits. Coverage uses all requested replicates and
retains unavailable intervals and adverse results. The original survey and Lenis
benchmarks do not establish dependent-design or estimated-weight inference.
[Historical archives](results/REPRODUCE.md) remain available separately.

## Real-data workflows

- [ECLS-K](applications/ecls/README.md): official data acquisition, selected-field
  decoding and corrected-cohort construction; then the
  [common analysis driver](applications/APPLICATION_PIPELINE.md) with explicitly
  supplied counts and source-order IDs.
- [MEPS 2009](applications/meps2009/README.md): official input preparation and
  phased point estimation, balance, counts, refits and contribution reporting.

These applications estimate supplied-weight standardized observed-outcome
contrasts in their declared selected cohorts. Their instructions state the
identification and inference limits. Participant data, fitted objects and count
archives remain private. Historical prepared NHANES/NSDUH helpers remain in
`applications/` for compatibility.

## Software checks and provenance

Build and check the source package after installing its suggested dependencies:

```sh
R CMD build .
R CMD check --no-manual wdsmatch_0.3.0.tar.gz
```

The executable HTML vignette requires Pandoc to build. The reference manual
requires a TeX toolchain when checked separately.
[Validation index](validation/README.md) covers estimator identities, aligned unit-weight
reductions, fitted covariance, count-refit arithmetic and failure handling.
Individual script dependencies and source/installed-package scopes are stated
in their headers. [SOURCE_PROVENANCE.json](SOURCE_PROVENANCE.json) records the
upstream source and version-bound checks. Software tests, saved-result replay
and repeated-sampling performance provide different evidence. See the
[current source-bound checks](docs/software_validation.md). GitHub Actions runs
package checks on Linux, macOS and Windows; its observed status is available in
the [Actions tab](https://github.com/ykzeng-yale/WDSM-Generalization/actions).

## License and reporting issues

The project-authored software is distributed under [GPL-3](COPYING). Author
credits are retained in [DESCRIPTION](DESCRIPTION). External source code and
application data must be obtained under their own terms as described in the
workflow instructions. Manuscripts and private research records are not part of
this repository.

For a reproducible bug report, use the [issue tracker](https://github.com/ykzeng-yale/WDSM-Generalization/issues)
and include `sessionInfo()`, the package commit and a minimal synthetic example.
Do not upload participant records or credentials.
