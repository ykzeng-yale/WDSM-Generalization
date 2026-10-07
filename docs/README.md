# Documentation

Current documentation is for the independent `WeightedMatching` 0.1.0 package.
The original `wdsmatch` 0.2.1 remains a separately installed predecessor.

- [User options and predecessor interfaces](user_options.md): donor count,
  distance/scaling, weights, targets, correction, ties and inference, including
  differences from Matching, MatchIt, dsmatch and the original WDSM API.
- [Getting started](../vignettes/weighted-matching.Rmd): executable package
  tutorial for supplied coordinates, fitted DSM, PATE/PATT and the two
  replication operators. Build/install the package archive with the vignette
  prerequisites in the root README, then use
  `vignette("weighted-matching", package="WeightedMatching")`. A direct raw-source
  `R CMD INSTALL .` does not build the HTML tutorial.
- [Theory and method guide](theory_method_guide.md): estimator, nuisance inputs,
  covariance scopes, assumptions and numerical limits.
- [Detailed API reference and earlier development notes](advanced_reference.md):
  the previous long README from the old `wdsmatch`-named development checkout,
  retained for method detail and historical provenance. Its old installation
  commands and legacy export claims do not describe `WeightedMatching`.
  The current root README and R help govern current interfaces.
- [Software validation](software_validation.md): release checks and their scope.

`help(package="WeightedMatching")` lists the installed function reference. Use
`citation("WeightedMatching")` for package metadata, and record the actual Git commit.
The current source release is not a CRAN release of this version or a public
distribution of the manuscript.
