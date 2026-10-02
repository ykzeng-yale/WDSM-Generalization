\name{wm_scalar_root_certificate}
\alias{wm_scalar_root_certificate}
\title{Certify a Stored Scalar Logistic Root and Matching Donor Sets}
\usage{
wm_scalar_root_certificate(object, design, center = NULL,
                          radius = 1e-4, precision = 128L)
}
\arguments{
\item{object}{An uncorrected, unrestricted raw result from
\code{wm_scalar_logistic_match()}, with its original data, fitted probabilities
and donor graph retained. The result must have status \code{point_computed}.}
\item{design}{The complete numeric design matrix used for propensity fitting,
in the original row and column order, including the leading intercept.}
\item{center}{Optional finite numeric coefficient vector. By default QR applied
to the stored probability logits proposes a center. This is an arbitrary
proposal, not recovery or certification of the solver's original coefficients.}
\item{radius}{Positive Euclidean radius in the proposed preconditioned
coefficient coordinates. Failure at this radius is inconclusive.}
\item{precision}{Integer binary precision from 64 to 4096 for outward rational
enclosures. Increasing precision may tighten a certificate; it does not change
the stored fit or the sufficient mathematical conditions.}
}
\description{
Performs a separate, optional finite-data check of the weighted logistic
score root, the actual stored nearest-neighbor donor sets, and the numerical
point contrast. Uses the optional \pkg{gmp} and \pkg{Rmpfr} packages.
It neither fits a model nor changes or filters the supplied result.
}
\details{
The exact problem interprets the stored binary64 design, treatment, positive
weights and outcomes as exact real inputs. A common exact weight rescaling
preserves that problem. Floating linear algebra proposes a coordinate matrix
and center; it is not trusted as a certificate.

In coordinates \eqn{\beta=\beta_0+C t}, the routine encloses the gradient at
zero and a symmetric Hessian lower bound throughout the specified ball.
Exact rational row calculations, bounds on logistic derivatives and a
Gershgorin inequality give a lower curvature bound \eqn{c}. A gradient-norm
upper bound \eqn{b<c r} proves existence and global uniqueness of a finite
score root, with \eqn{\|t_*\|\le b/c}. Coefficient bounds refer to the proposed
center, which need not be the unretained original solver coefficient.

Exact GMP integer floor/ceiling operations construct outward dyadic
intervals. MPFR exponential results are padded using their precision with
explicit exponent-range checks. Conversion of dyadic endpoints to MPFR is
checked by exact reverse conversion. Rational endpoint strings in each
bound record are authoritative; double endpoints are converted outward and
may be infinite when they are not representable as finite doubles. The
\code{numeric_display_finite} field records this distinction.

The root enclosure bounds the discrepancy of each stored probability from
its exact-root value. Selected donor sets are checked against all eligible
opposite-arm donors using exact stored-probability distances. A strict gap
larger than four times the largest probability-error bound certifies set
membership. Boundary ties or an insufficient gap are inconclusive. If all
eligible donors are selected, membership is automatic once the indices and
root are certified. This checks membership, not donor ordering or bitwise
identity of differently ordered floating sums.

Separate interval evaluation of the raw weighted contrast bounds the
stored point estimate's arithmetic error. Only if every required donor set
is certified does this also bound its discrepancy from the exact-root
matching statistic. PATE checks both directions and PATT checks treated
queries with control donors. This purely finite-data check has no
three-column restriction; the separate statistical variance theorem retains
its own design and population restrictions.

A passing dataset does not prove high-probability certification success,
an asymptotic solver-accuracy rate, accurate variance-moment arithmetic,
population assumptions, or a bootstrap theorem. A failing sufficient
inequality is not proof that the estimator or exact root is invalid. Keep
all scientific results, including failed or inconclusive checks. The
original object's \code{numerical_root_bridge} field is not modified;
inspect this separate certificate's component statuses.
}
\value{
A list of class \code{wm_scalar_root_certificate}. The status is
\code{not_certified}, \code{root_contained_graph_unresolved}, or
\code{root_and_donor_sets_certified}. It includes bound input values,
center provenance, arithmetic runtime and precision, separate root and
graph flags, rigorous bound records and unresolved query indices when
available. Bounds contain outward numeric endpoints and exact rational
endpoint strings. \code{point_root_error_bound} is present only when the
root and all donor sets are certified. Malformed or mismatched objects
raise an error; a failed certificate calculation returns a reason with
the stored inputs and leaves the supplied object unchanged.
}
\seealso{\code{\link{wm_scalar_logistic_match}}, \code{\link{wm_match}}}
