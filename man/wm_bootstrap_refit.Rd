\name{wm_bootstrap_refit}
\alias{wm_bootstrap_refit}
\title{Fixed-Reuse Multinomial Replication}
\usage{
wm_bootstrap_refit(object, counts, refit = NULL, conf.level = 0.95)
}
\arguments{
\item{object}{A self-normalized \code{wm_match} or \code{wm_fit} result. Variance estimation need not have been enabled.}
\item{counts}{An n-by-B matrix of nonnegative integer multinomial counts, with each column summing to n and B at least two. Each arm must have positive total multiplicity in every column. Optional row names must be literal original indices 1,...,n in order. For example, use \code{rmultinom(B, n, rep(1/n, n))}. Supplying counts permits exact paired comparisons.}
\item{refit}{NULL to hold predictions fixed, or a function accepting one length-n multiplicity vector and returning a list with finite length-n \code{mean0} and, for PATE, \code{mean1}. The callback must implement the intended nuisance-fitting rule, including any score refitting, standardization and bias-correction fit, and predict at the original n rows. PATT uses only \code{mean0}. Raw matching requires NULL.}
\item{conf.level}{Confidence level strictly between zero and one.}
}
\description{
Implements the original WDSM fixed-reuse replication calculation with arbitrary matching coordinates, including distinct arm-specific maps. Matching graphs, original individual weights and weighted donor reuse remain fixed. Nuisance predictions may be refitted through a callback. The common \code{wm_bootstrap} interface delegates to this calculation with \code{method = "fixed_reuse"}; its \code{method = "contribution"} uses fixed complete contributions instead.
}
\details{
Objects created with \code{tie_rule = "source_random"} retain that policy's
unsupported sampling-inference scope. Their returned variance and interval
are empirical fixed-reuse calculations, not a discrete-tie or positive-boundary
sampling theorem; \code{sampling_inference_available} is \code{FALSE}.

For PATE, a replicate is the sum of \eqn{m_i W_i(q_{1i}^*-q_{0i}^*) + m_i(2Z_i-1)(W_i+K_i)(Y_i-q_{Z_i i}^*)}, divided by \eqn{\sum_i m_i W_i}. For PATT, its numerator is \eqn{\sum_i m_i Z_i W_i(Y_i-q_{0i}^*)-\sum_i m_i(1-Z_i)K_i(Y_i-q_{0i}^*)}, and its denominator is \eqn{\sum_i m_i Z_i W_i}. Here K is the original weighted reuse, including the normalized donor fraction, with no extra division by M. Calculations retain the common weight normalization from the fitted object.

For the original retrospective WDSM protocol, the callback refits PS using multiplicity times individual weight, arm-specific PG using multiplicity, and the complete quadratic bias correction using multiplicity times individual weight on refitted standardized double scores. Externally supplied scores stay fixed if that was the original analysis. The callback must not rebuild the original matching graph or reuse coefficients.

The returned Monte Carlo variance uses divisor B, and the normal Wald interval is centered on the original estimate. All replicates must succeed; failed fits raise an error and are not discarded. The function itself draws no random numbers; any callback randomness is controlled by the caller.

For a retained fitted-score wrapper, use
\code{wm_bootstrap(fitted, method="fixed_reuse", counts=counts)}.
Its default quadratic adapter performs the full same-count prediction refit;
this low-level function continues to require an ordinary \code{wm_match} or
\code{wm_fit}, and NULL here continues to hold its predictions fixed.

This is an exact replication algorithm, not an automatic inference theorem. Valid inference additionally requires the relevant nuisance-refit expansion and agreement between sampling and replication variance. Empirical equivalence to another implementation does not establish those conditions. Weights remain fixed; estimated-weight uncertainty is not refitted. Influence corrections already attached to the object are not added to these refitted replicates. Stabilized, potential-mean and split-block procedures use different representations and are not accepted.

For the supported correct-full-X scalar sampling alternatives, agreement
with the sampling limit additionally requires the actual callback to refit
the complete regular estimating stack using the same full-n multiplicity
vector in every block, including shared/generated score, prediction,
standardization and correction parameters. The selected original and
nearby conditional statistical roots must satisfy the regular-root and
actual graph-row requirements. All same-row matching/nuisance and nuisance
cross-covariances remain. The actual finite-M frozen-edge prediction slope
then equals the sampling slope. Returning only mean predictions is compatible
with an internal full refit, but neither that return type nor a convergence
flag verifies it. A compatible prediction-only callback may instead return the
corresponding same-count raw weighted-least-squares roots when its fixed raw
prediction equations are independent of inactive score coordinates and the
complete-system and actual-row conditions hold. Its effect expansion retains
the complete influence and covariance, without returning inactive increments.
Complete polynomial regressions preserve their count-refitted predictions under
invertible affine standardization, whether the basis is held fixed or recomputed.
For full-covariate matching in a nonsingular affine chart, continuous prediction
derivatives depending only on those matched covariates give zero prediction and
sampling slopes under the applicable current-graph and full-design sampling laws;
the marked weight variance remains. Other partial callbacks require their actual
influence and covariance. A separate full-X prediction regression under scalar
propensity matching need not have zero slope.

A numerical conditional refit needs \eqn{o_{p^*}(n^{-1/2})} parameter error
relative to the specified nearby statistical root and the same localized
stack. For the outcome-fitted scalar alternative, the original graph still
needs its exact graph root or a separately justified numerical graph
transfer; an original \eqn{o_p(n^{-1/2})} parameter error alone is insufficient.
This function checks finite prediction/ratio algebra, not those population
or numerical-root conditions.

Conditional weak agreement does not make divisor-B variance consistent for
fixed B or by itself validate its finite-B normal interval. A separate empirical
variance result uses deterministic \eqn{B_n\to\infty} with
\eqn{B_n\le n^b} eventually for a fixed finite b, conditionally independent
full-n multinomial columns, the same complete iid regular-root system and
finite-moment actual-row/refit conditions. Prediction derivative and Hessian
envelopes must have a finite moment of some order q greater than two. For
callbacks selecting the corresponding exact nearby count roots,
\deqn{\frac{n}{B_n}\sum_{j=1}^{B_n}
 (\widehat\tau^*_{j,\mathrm{refit}}-
  \overline{\widehat\tau^*_{\mathrm{refit}}})^2
 \to_p V_{\mathrm{rep}}.}
When the proved sampling and replication variances agree and their common limit
is positive, this empirical variance and the sampling CLT justify normal
intervals. Numerical callbacks require simultaneous selection, success and
compatible effect-error control over all prescribed draws; single-draw numerical
accuracy alone is insufficient. This convergence-in-probability result does not
assert convergence of the nonlinear conditional second moment. Arbitrary
prescribed growing B outside this result retains the separate conditional
moment/uniform-integrability and failure/completion requirements. Failed columns must be
retained with their supplied counts rather than silently redrawn. This
implementation stops on failure; it does not provide a completion or certify
raw finite-sample sampling-variance convergence.
}
\value{
A list containing the original \code{estimate}, replicate \code{draws}, \code{root_n_draws}, Monte Carlo \code{variance} and \code{se}, centered normal \code{conf.int}, \code{variance_divisor}, \code{B}, \code{n}, \code{conf.level}, \code{method}, prediction-refitting status and an inference-contract description.
}
\seealso{\code{\link{wm_match}}, \code{\link{wm_bootstrap}}, \code{\link{wdsmatchATE}}}
