# Software validation entry points

The installed-package unit tests run during `R CMD check`. The standalone
scripts in this directory check reference calculations, source-study adapters
and archived reporting separately; read each script's header for dependencies
and source/installed-package mode.

| Check | Main entry points |
| --- | --- |
| Original donor-normalized estimator and target arithmetic | `check_core_reference.R`, `wm_core_reference.R` |
| Properly aligned all-ones PSM/full-X/DSM and WDSM reductions | `check_unweighted_reductions.R`, `check_legacy_parity.R` |
| Complete fitted inference and contribution dispatch | `check_fitted_wm_components.R`, `check_general_fitted_bootstrap.R`, `check_wm_bootstrap_dispatch.R` |
| Original fixed-reuse arithmetic | `tests/testthat/test-wm-bootstrap-refit.R`, `tests/testthat/test-wm-bootstrap-fit.R` |
| Comparison and saved reporting definitions | `check_comparison_metrics.R`, `check_paired_comparison_metrics.R`; [result replay](../results/README.md) |

Tests compare finite calculations, option validation and failure handling. They
do not establish causal identification, geometric/root assumptions, asymptotic
numerical availability or nominal coverage in a new population. Diagnostic and
historical scripts are not required steps for ordinary package users. Current
release checks are recorded in [software validation](../docs/software_validation.md).
