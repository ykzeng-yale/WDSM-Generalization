\name{wm_geometry}
\alias{wm_geometry}
\title{Universal Euclidean Matching Overlap Constants}
\usage{
wm_geometry(M, d, draws = 4096L, seed = NULL, chunk_size = 256L)
}
\arguments{
  \item{M}{Positive integer fixed number of matches.}
  \item{d}{Positive integer Euclidean matching dimension.}
  \item{draws}{Number of coupled angular Monte Carlo draws per component.
    Must be at least two when d exceeds one. Ignored for the exact d=1 branch.}
  \item{seed}{Optional nonnegative integer seed. An explicit seed restores
    the caller's RNG state. With seed = NULL, numerical integration advances
    the current stream. In the numerical branch, explicit seeds are unsupported with the Box--Muller
    normal generator because its hidden cache cannot be restored.}
  \item{chunk_size}{Requested maximum number of coupled draws per batch.
    Automatically reduced to respect the allocation cap.}
}
\description{
Returns the universal overlap vector beta for fixed-M Euclidean matching and
its sum alpha. The one-dimensional formula is exact. In larger dimensions,
the function estimates the proved finite angular integral and reports empirical
Monte Carlo standard errors, including dependence among its component estimates.
It does not estimate a data-specific matching variance.
}
\details{
The component with index ell counts common nonroot donors in two M-nearest
sets that both contain the root. Thus ell runs from zero through M-1.
For d=1,
\deqn{\beta_{M,1,\ell}=2(\ell+1)\quad(\ell<M-1),
       \qquad \beta_{M,1,M-1}=3M/2,}
and \eqn{\alpha(M,1)=M^2+M/2}. This branch makes no random draws.

For every dimension, write \eqn{a=M-1-\ell}, \eqn{k=2M-\ell} and \eqn{D=dk}.
A point on the unit sphere in D coordinates consists of two d-dimensional
queries y,t, ell shared donors, a donors private to y and a donors private to t.
The root is the origin. Each query's radius is the largest distance to the
root and its own specified donors. Let L be the volume of the union of these
two radius balls. Let I indicate that every private donor lies strictly
outside the other query's ball. The exact integral definition is
\deqn{\beta_{M,d,\ell}=
 \frac{\Gamma(k)}{d\,\ell!\,(a!)^2}
 \frac{2\pi^{D/2}}{\Gamma(D/2)}
 E\{I(\Theta)L(\Theta)^{-k}\},}
where Theta is uniform on the unit sphere. The sphere-area and factorial
factors are included in every numerical draw. Ball-union volume is evaluated
by adding the two volumes and subtracting their intersection, using the two
spherical caps and regularized incomplete beta probabilities. Log volumes
and log normalization factors avoid unnecessary underflow or overflow.

For each replicate, one vector of 2Md independent standard normals is drawn.
For each component, its first D entries are normalized to its own unit sphere.
This coupling gives the correct marginal uniform-sphere law in every component.
The reported covariance matrix is the estimated covariance of the component
sample means. Alpha's Monte Carlo standard error is calculated from the
per-draw sum, retaining all within-draw covariances; adding component variances
alone would be incorrect for this coupling.

The default run is bounded. The numerical branch requires M at most 128 and
at most 50 million component coordinates, calculated as
\eqn{\mathrm{draws}\,d(3M^2+M)/2}. A batch allocates at most one million
Gaussian coordinates. The exact branch permits at most one million output
components. Numerically unrepresentable integrands or sample moments are
rejected. Larger or more efficient integrations require a separate procedure.

Monte Carlo standard errors are empirical diagnostics, not rigorous error
bounds or evidence that rare influential configurations have been sampled.
A component with no nonzero samples is marked unresolved; its observed zero
MCSE is uninformative. The exact inequality alpha >= M squared is reported
as a diagnostic. Numerical estimates are never clipped to satisfy it, and
previously published numerical tables are not used as exact lookup values.
}
\value{
An object of class \code{wm_geometry}, containing:
\describe{
  \item{M, d, overlap}{Match count, dimension, and nonroot-overlap indices.}
  \item{beta, alpha}{Named component vector and its sum.}
  \item{beta_mcse, alpha_mcse}{Empirical Monte Carlo standard errors; zero for
    the exact one-dimensional formula.}
  \item{covariance}{Estimated covariance matrix of the beta estimates for
    numerical integration. NULL for the exact branch.}
  \item{method, precision_status, component_precision}{Explicit method and
    precision labels. Numerical estimates are never labeled exact.}
  \item{draws, requested_draws, seed}{Actual and requested draws and seed.
    Actual draws are zero for d=1.}
  \item{chunk_size, requested_chunk_size}{Actual and requested batch sizes.}
  \item{definition}{The exact angular-integral definition.}
  \item{normalization}{For the numerical branch, each component's point count,
    sphere dimension, log sphere area and complete log prefactor.}
  \item{diagnostics}{Nonzero counts and maximum observed integrand values for
    numerical integration, the Jensen lower bound, whether the estimate is
    below it, RNG usage, coupling and resource-cap metadata.}
}
}
\examples{
wm_geometry(M = 3, d = 1)$beta
\dontrun{
numerical <- wm_geometry(M = 2, d = 2, draws = 4096, seed = 2026)
cbind(estimate = numerical$beta, MCSE = numerical$beta_mcse)
}
}
