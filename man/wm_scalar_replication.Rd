\name{wm_scalar_replication}
\alias{wm_scalar_replication}
\title{Current-Graph Contribution Replication for Qualified Scalar Matching}
\usage{
wm_scalar_replication(object, design = NULL, counts = NULL,
  multipliers = NULL, conf.level = 0.95)
}
\arguments{
\item{object}{A raw \code{wm_scalar_logistic_match} result with the stored
fitting-design binding and complete scalar inference predictions and blocks.}
\item{design}{Optional original design matrix. If supplied, its numeric values
and column order must equal the stored binding. No rows are reordered.}
\item{counts}{Optional finite numeric n-by-B nonnegative integer matrix, B at
least one, with every column summing to n.}
\item{multipliers}{Optional finite numeric n-by-B matrix, B at least one.
Supply counts or multipliers, not both.}
\item{conf.level}{Confidence level strictly between zero and one.}
}
\description{
Constructs complete matching contributions on the stored scalar graph,
centers them and the ordinary logistic influences, and adds the signed full
matching drift. Reports their empirical variance and optional supplied-input
replication algebra while retaining the original raw point and fit.}
\details{
The default uses no draws. No propensity or outcome model is fitted, no
neighbor is searched, and no random number is generated. The own-arm mean
predictions are used only for uncertainty. Their contribution rows need not
sum to zero around the raw estimate and are explicitly centered. PATT keeps
its complete direct treated-outcome row and uses only the control prediction
in that row. All normalizations, coefficient influences and replication roots
use the full observed sample size n for both estimands.

If Lc denotes the centered normalized matching row, phic the centered
logistic influence, and b the stored full signed matching drift, the augmented
row is A = Lc + b' phic. The root variance is mean(A^2), using divisor n.
The point variance is this quantity divided by n. The complete centered
V0, C and Sigma row blocks are returned. They are distinct finite-sample
estimates from the explicit scalar moment blocks stored in the input;
the two variances need not agree exactly. Neither estimate replaces the
other in the input object.

Count roots equal sum((counts-1)*A)/sqrt(n); multiplier roots equal
sum(multipliers*A)/sqrt(n). Point replicates add root/sqrt(n) to the unchanged
point. Independent standard Gaussian multipliers have exactly the reported
conditional root variance. Full-n multinomial counts with uniform cell
probabilities have the same conditional variance. Input validation does not
verify a supplied matrix's distribution or independence. A count column
concentrated in one treatment arm is allowed: this linear procedure refits
nothing and uses no replicate Hajek denominator. Arm-conditioned or
treated-only draws have a different law. If B is at least two, the separately
named empirical draw variance uses divisor B; it is only a Monte Carlo
diagnostic, not the default reported variance.

Sampling validity is restricted to the bounded ordinary-logistic scalar
known-score-weight model and its proved current-graph row moments and
maximum-row bounds. The original fitted-probability distance rule, positive
smooth known weights, population identification, accepted conditional-score
chart, outcome model, mean-estimation rates and exact-root or numerical-bridge
conditions remain required. This does not address arbitrary within-score
weights, estimated weights, general outcome-adaptive scalar fitting,
dependent sampling or the original nuisance-refit bootstrap. The inherited
numerical-root status and unverified population assumptions are retained;
empirical centering is not a numerical certificate.

The stored design binding is required. Legacy objects without it return an
unavailable result even when a design is supplied: J/B agreement alone does
not identify the original row-wise design. Missing scalar prediction or
nuisance blocks similarly produce an unavailable result without fitting
anything. Malformed or inconsistent inputs raise an error. Runtime arithmetic
failures preserve the raw fit and a reason. A finite nonpositive original
explicit variance does not prevent calculation when its required blocks exist.
Nonpositive or unrepresentable row-based interval arithmetic withholds that
interval without flooring or altering the original point.
}
\value{
A \code{wm_scalar_replication} list containing the preserved fit and explicit
inference record, binding/scope information, centered and uncentered rows,
full drift and empirical row blocks, analytic root and point variances,
interval status, an explicit-versus-row variance comparison, and any supplied
draws and calculated roots. \code{row_variance_finite} checks finite
nonnegative arithmetic only; it does not certify population assumptions.
}
\seealso{\code{\link{wm_scalar_logistic_match}},
\code{\link{wm_scalar_root_certificate}}, \code{\link{wm_bootstrap_refit}}}
