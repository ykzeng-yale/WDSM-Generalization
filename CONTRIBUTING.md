# Contributing

Use a minimal synthetic example when reporting an error. Include the package
version, Git commit, `sessionInfo()`, target, matching coordinates/scaling, M,
correction and inference method. Do not upload participant records, credentials
or unpublished manuscripts.

Install the suggested dependencies in DESCRIPTION, then run:

```sh
R CMD build .
R CMD check --no-manual wdsmatch_0.3.0.tar.gz
```

Building the executable HTML vignette requires Pandoc. Checking the PDF reference
manual separately requires TeX. The package tests use synthetic fixtures; source
studies, result replay and external-data workflows have separate dependencies.

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
