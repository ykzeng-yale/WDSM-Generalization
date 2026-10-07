# Documentation

- [User options and predecessor interfaces](user_options.md): donor count,
  distance/scaling, weights, targets, correction, ties and inference, including
  differences from Matching, MatchIt, dsmatch and the original WDSM API.
- [Getting started](../vignettes/weighted-matching.Rmd): executable package
  tutorial for supplied coordinates, fitted DSM, PATE/PATT and the two
  replication operators. Build/install the package archive with the vignette
  prerequisites in the root README, then use
  `vignette("weighted-matching", package="wdsmatch")`. A direct raw-source
  `R CMD INSTALL .` does not build the HTML tutorial.
- [Theory and method guide](theory_method_guide.md): estimator, nuisance inputs,
  covariance scopes, assumptions and numerical limits.
- [Detailed API reference and earlier development notes](advanced_reference.md):
  the previous long README, retained for specialized branches and historical
  provenance. The current root README and R help govern updated interfaces.
- [Software validation](software_validation.md): release checks and their scope.

`help(package="wdsmatch")` lists the installed function reference. Use
`citation("wdsmatch")` for package metadata, and record the actual Git commit.
The current source release is not a CRAN release of this version or a public
distribution of the manuscript.
