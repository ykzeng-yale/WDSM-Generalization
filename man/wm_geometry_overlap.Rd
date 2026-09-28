\name{wm_geometry_overlap}
\alias{wm_geometry_overlap}
\title{Experimental Boundary-Rank Integration of Matching Overlap Constants}
\usage{
wm_geometry_overlap(M, d, draws = 4096L, seed = NULL, chunk_size = 256L)
}
\arguments{
  \item{M}{Positive integer number of distinct nearest donors per query.}
  \item{d}{Positive integer Euclidean matching dimension.}
  \item{draws}{Requested Monte Carlo draws. Numerical integration needs at
  least two. One-dimensional public calls use the exact formula instead.}
  \item{seed}{Optional nonnegative integer seed. An explicit seed restores
  the previous global random state, including its absence. The Box-Muller
  normal generator is rejected with an explicit seed because its cache
  cannot be restored by restoring the R seed vector alone.}
  \item{chunk_size}{Requested batch size, reduced when necessary to respect
  the one-million-element working-array cap.}
}
\description{
An experimental, separate numerical method for the complete
order-\eqn{M} neighbor-set overlap constants. The original
\code{\link{wm_geometry}} method is unchanged. A successful numerical
call is a Monte Carlo estimate, not an exact constant or a precision guarantee.
}
\details{
The overlap index \eqn{\ell=0,\ldots,M-1} counts shared donors other than
the Palm root. The integral partitions the farthest selected donors
in the two query neighborhoods into four families: both maxima at the
root; one at the root and the other at a nonroot donor; two distinct
nonroot boundary donors; or the same shared nonroot donor at both boundaries.
Symmetry gives multiplicity two to the second family.

All remaining shared and private donor positions are integrated into
intersection and exclusive-ball volumes. Surface additions are positive;
they are not ordinary derivatives of the private-volume exclusions.
The shared boundary donor is retained explicitly, avoiding a singular
simultaneous-sphere derivative. The four sphere dimensions are
\eqn{2d}, \eqn{2d+1}, \eqn{2d+2}, and \eqn{3d}.
With \eqn{a=M-1-\ell} and \eqn{k=2M-\ell}, each term has factor
\deqn{\frac{\Gamma(k)\,|S^{D-1}|}{d\,\ell!\,(a!)^2},\qquad
 |S^{D-1}|=\frac{2\pi^{D/2}}{\Gamma(D/2)},}
with the second family's multiplicity included. Structural zero terms
are omitted before evaluating powers; in particular a zero exponent
does not multiply a negative infinite log volume.

Every family is sampled by normalizing a prefix of the same Gaussian
vector. Consequently all overlap and boundary terms are correlated.
Their full covariance is retained and propagated into the returned
\code{covariance}, component Monte Carlo standard errors, and the
standard error of their sum \eqn{\alpha}.
The \code{one_nonroot} term estimate already includes its factor two.

All these integrands are bounded for each fixed \eqn{M,d}, but the bound
can be large. An empirical Monte Carlo standard error is not a rigorous
error bound. Active zero-hit or zero empirical-variance terms are
explicitly unresolved. Constants are never projected onto
\eqn{\alpha\ge M^2}, and numerical values are not selected or tuned to
match a desired reference.

Numerical calls require \eqn{M\le64}, \eqn{d\le128}, at most 50 million
total sphere coordinates and 100 million covariance products.
The covariance-product count is \eqn{\mathtt{draws}(4M)^2}.
Only active sphere families count toward geometric coordinate work.
These limits are checked before drawing. One-dimensional exact calls
use the allocation restrictions of \code{wm_geometry}.

For \eqn{d=1}, the public function returns the exact vector:
\eqn{\beta_\ell=2(\ell+1)} for \eqn{\ell<M-1} and
\eqn{\beta_{M-1}=3M/2}. Its sum is \eqn{M^2+M/2}.
Standalone validation can exercise the private numerical integrator
in one dimension independently of this public reduction.
}
\value{
A \code{wm_geometry} object, compatible with the original geometry
result's \code{M}, \code{d}, \code{overlap}, \code{beta}, \code{alpha},
\code{beta_mcse}, \code{alpha_mcse}, \code{covariance}, seed/draw,
method and precision fields.

Numerical results have method \code{boundary_rank_monte_carlo}.
Additional \code{term_estimates} and \code{term_mcse} matrices have
overlap rows and four boundary-family columns.
\code{term_covariance} contains covariance of the estimated term means,
ordered by family and then overlap, including cross-family covariance.
\code{covariance} is covariance of the estimated beta means.
\code{term_precision}, nonzero counts, maximum sampled values,
normalization factors, work counts and coupling information are retained.
The returned validation status identifies this as an experimental
method requiring independent numerical validation.

Exact one-dimensional results have the original exact-method fields,
zero Monte Carlo standard errors, no sampling and \code{NULL} covariance.
}
\seealso{\code{\link{wm_geometry}}, \code{\link{wm_geometry_alpha}}}
\examples{
wm_geometry_overlap(M = 3, d = 1)
wm_geometry_overlap(M = 1, d = 2, draws = 128, seed = 17)
}
