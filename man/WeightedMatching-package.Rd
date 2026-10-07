\name{WeightedMatching-package}
\alias{WeightedMatching-package}
\alias{WeightedMatching}
\docType{package}
\title{Weighted Matching for Causal Inference}
\description{
Fixed-count nearest-neighbor matching with supplied positive probability weights
for population average and treated-population effects. Matching maps may have
any fixed positive dimension and may differ between treatment arms.
}
\details{
Use \code{\link{wm_match}} for supplied matching maps,
\code{\link{wm_fit}} for fitted outcome correction,
\code{\link{wm_wdsm_fit}} for arm-specific double scores, and
\code{\link{wm_model_fit}} for model-list matching maps.
Inference procedures and their assumptions are documented with each function.
The weights are treated as supplied; estimated-weight uncertainty and general
clustered or stratified survey inference are outside this framework.

This is a separate package from \pkg{wdsmatch}. It does not replace that
package's \code{wdsmatchATE} or \code{wdsmatchATT} interfaces.
}
\seealso{\code{\link{wm_bootstrap}}, \code{\link{wm_model_inference}}}
\keyword{package}
