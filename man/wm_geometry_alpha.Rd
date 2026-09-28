\name{wm_geometry_alpha}
\alias{wm_geometry_alpha}
\title{Reduced Integration of the Universal Catchment Second Moment}
\usage{
wm_geometry_alpha(M, d, draws = 4096L, seed = NULL, chunk_size = 256L)
}
\arguments{
  \item{M}{Positive integer fixed donor count.}
  \item{d}{Positive integer Euclidean matching dimension.}
  \item{draws}{Positive integer requested draw count. Numerical integration
  requires at least two draws. The exact one-dimensional branch uses none.}
  \item{seed}{Optional nonnegative integer seed. An explicit seed restores
  the caller's RNG state, including its prior absence. The cached
  Box--Muller normal generator is rejected when an explicit seed is used.}
  \item{chunk_size}{Positive integer maximum batch size. The effective
  value is reduced by the allocation guard.}
}
\description{
Computes \eqn{\alpha_d(M)}, the second moment of a distinguished donor's
order-\eqn{M} Poisson catchment volume, by analytically integrating all
nonroot donor positions in the joint root-inclusion probability.
Returns alpha alone. The separate \code{\link{wm_geometry}} function
returns the complete-neighbor overlap beta vector.
}
\details{
For query positions \eqn{x,y}, form balls with respective radii
\eqn{\|x\|,\|y\|}, so their interiors contain donors closer than the
distinguished donor at zero. Let \eqn{I,E_1,E_2,U} be intersection,
two exclusive, and union volumes, and put
\eqn{p=(I,E_1,E_2)/U}. With \eqn{(x,y)} uniform on the unit sphere
in \eqn{\mathbb R^{2d}},
\deqn{\alpha_d(M)=
 E\left[\frac{|S^{2d-1}|}{dU^2}
 \sum_{\ell=0}^{M-1}\sum_{a,b=0}^{M-1-\ell}
 \frac{(\ell+a+b+1)!}{\ell!a!b!}p_0^\ell p_1^a p_2^b\right].}
The finite count sum includes every pair of counts at most \eqn{M-1}.
Restricting both counts to exactly \eqn{M-1} is incorrect when \eqn{M>1}.
The intersection count here is not the overlap of the two complete
nearest-neighbor sets and must not be used as an overlap-beta index.

The integrand is positive and at most
\deqn{B_{d,M}=
 \frac{2^{d+1}\Gamma(d/2+1)^2}{d\Gamma(d)}M(2M-1).}
There are no exclusion-based zero-hit components. The returned
empirical Monte Carlo SE is accompanied by the conservative upper
bound \eqn{B_{d,M}/(2\sqrt{\mathrm{draws}})} on the true Monte Carlo SE
and a Hoeffding 95 percent half-width
\eqn{B_{d,M}\sqrt{\log(40)/(2\,\mathrm{draws})}}.
These analytic bounds concern the mathematical integrand and do not
include floating-point error. A zero empirical variance is labeled
unresolved, not treated as proof of exactness. Estimates are not clipped
to the Jensen lower bound \eqn{M^2}.

For \eqn{d=1}, the exact expression \eqn{M^2+M/2} is returned, with no
RNG use. Otherwise independent Gaussian vectors normalized onto the
\eqn{2d}-dimensional sphere generate the draws. Each replicate's
coordinates are generated consecutively, so changing batch size
preserves the sampled integrand values up to accumulation round-off.
Runtime is guarded by \eqn{M,d\le128}, at most 50 million count terms
and at most 50 million generated coordinates. The count sum has
\eqn{M(M+1)(2M+1)/6} terms per draw.

This routine supplies stable direct alpha references for formulas
depending only on alpha and for checking the sum of a separately
estimated beta vector. It does not replace the overlap components
needed for general marked nonlinear donor-weight formulas.
}
\value{
A \code{wm_geometry_alpha} list containing \code{alpha},
\code{alpha_mcse}, method and precision labels, dimensions, donor
count, draw/seed/batch metadata and the explicit alpha-only scope.
Numerical diagnostics include integrand and Monte Carlo SE bounds,
a Hoeffding half-width, observed integrand range, work counts and the
unmodified estimate's Jensen-bound diagnostic. No beta vector or
beta-component covariance matrix is returned.
}
\seealso{\code{\link{wm_geometry}}}
\examples{
wm_geometry_alpha(M = 3, d = 1)$alpha
}
