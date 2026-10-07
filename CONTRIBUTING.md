# Contributing

This checkout develops the independent `WeightedMatching` package, version
0.1.0. The original `wdsmatch` 0.2.1 is a separately installed predecessor.
Current public statistical interfaces use `wm_*`; legacy WDSM code retained
internally for source comparisons is not a public replacement for that package.

Use a minimal synthetic example when reporting an error. Include the package
version, Git commit, `sessionInfo()`, target, matching coordinates/scaling, M,
correction and inference method. Do not upload participant records, credentials
or unpublished manuscripts.

Install the suggested dependencies in DESCRIPTION, then run:

```sh
R CMD build .
R CMD check --no-manual WeightedMatching_0.1.0.tar.gz
```

Building the executable HTML vignette requires Pandoc. Checking the PDF reference
manual separately requires TeX. The package tests use synthetic fixtures; source
studies, result replay and external-data workflows have separate dependencies.

The installed reference in `man/*.Rd` and `NAMESPACE` is maintained alongside
the R source. Source roxygen comments are not complete enough to regenerate all
manual pages and exports: do not run an unreviewed `roxygen2::roxygenise()` over
this checkout. For an API change, update its Rd usage, argument/default and return
contracts, namespace/S3 registration, tutorial and tests together. `R CMD check`
checks current usage/documentation consistency; existing help does not establish
exhaustive behavioral coverage. The small scalar fixture under
`tests/testthat/fixtures/` is an exact copy of the documented synthetic validation
fixture, so those interface checks also run after source-package installation.

For estimator/inference changes, add a check against independently written
equations or a properly aligned comparator. Retain the complete nuisance
covariance, original graph and supplied weights where the method requires them.
Test PATE/PATT and failed-draw handling for affected routes. A unit test or an
available interval does not verify the population assumptions of a theorem.

Do not overwrite frozen simulation arrays, failure denominators or source
manifests for a package change. Reproduction of the original principal producer
uses the public frozen checkout specified in its README. A new algorithm or
study needs a separate declared source/version and results. Documentation fixes
do not require a new Monte Carlo campaign.
