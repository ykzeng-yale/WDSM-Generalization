\name{wm_model_inference}
\alias{wm_model_inference}
\title{Qualified Inference for a Retained Finite-Model Matching Fit}
\usage{
wm_model_inference(object, covariance_scope,
  transport0 = NULL, transport1 = NULL,
  critical_control = NULL, conf.level = 0.95,
  max_influence_elements = 5e7)
}
\arguments{
\item{object}{A retained \code{wm_model_fit} result. The complete nuisance influence
matrix, or its complete named estimating equations and Jacobian, and the required
prediction/reference derivatives must already be retained. A point-only
quadratic fit that omitted its influence matrix is not refitted.}
\item{covariance_scope}{An explicit scope supported by \code{wm_fitted_inference}:
PATE \code{"distinct_rarity"}, \code{"full_x"}, \code{"common_field"},
\code{"regular_joint_d_gt2"} or \code{"critical_planar"}; PATT \code{"patt"}.
The separate \code{"scalar_psm"} method is not a generic model handoff.}
\item{transport0, transport1}{Explicit direction-specific transport specifications
for \code{wm_fitted_inference}. PATT omits \code{transport1}; \code{full_x} omits
both. No zero-transport declaration is inferred from fitted arrays. In explicit
\code{mode="estimate", representation="score_measurable_weight"} only,
\code{tangents="fitted"} binds the retained actual matching-coordinate derivatives
and graph coordinates. Omit \code{raw_scores} and \code{raw_to_matching} with this
sentinel. It is unsupported for \code{conditional_weight_chart}; supply that
chart's genuine coordinates, derivatives and map explicitly.}
\item{critical_control}{The explicit \code{critical_planar} controls accepted by
\code{wm_fitted_inference}. Either endpoint's \code{tangents0} or \code{tangents1}
may equal \code{"fitted"} to use the corresponding retained actual-coordinate
derivatives. All other controls remain caller supplied.}
\item{conf.level}{Confidence level strictly between zero and one, forwarded
unchanged. A selected branch need not provide a normal confidence interval.}
\item{max_influence_elements}{Positive conservative element cap before completing
a missing influence matrix from retained equations. The planned bound is
\eqn{8np+8p^2+4p}, for full nuisance dimension p. Inputs, solver/R overhead and
total resident memory are outside this cap. No parameter column is removed.}
}
\description{
Binds the original matching fit, full nuisance influence, complete
prediction/reference derivatives and exact-zero known-weight derivative to
\code{wm_fitted_inference}. It performs no model, local-regression, graph or
weight refitting. The covariance and transport premises remain caller declarations.
}
\details{
The quadratic route uses its complete retained prediction derivatives and
nuisance influence, including correction coefficients. The local-polynomial
route uses its complete score-only influence and the plug-in reference derivative
\eqn{B^T\nabla\widehat m}. This is not differentiation through refitting local
coefficients. These different parameter layouts must never be mixed.

All original-row and parameter bindings are preserved. Supplied weights remain
known and may depend on treatment and covariates; their estimation uncertainty
and cluster/stratum inference are unsupported. No tangent sentinel supplies a
population chart, score-measurability condition, cutoff, support, density bound,
bandwidth or covariance law. The dimension restrictions and separate scalar
conditions of \code{wm_fitted_inference} remain in force.

If the influence matrix was omitted by point-only fitting, the helper may compute
\code{-t(solve(J, t(psi)))} and its centered empirical covariance from the retained
complete named Jacobian J and estimating-equation rows psi. It requires original
unit multiplicities, uniform empirical probabilities, the full parameter order,
and a zero known-weight derivative. The arrays are retained for subsequent
handoffs. This is linear algebra at the stored root, with no coefficient fitting,
root iteration, prediction evaluation or graph rebuilding. An allocation limit,
missing source arrays, singular solve or nonfinite result yields unavailable
inference; no ridge, column deletion or retry is used.

The local-polynomial correction transfer and reference-derivative consistency
require their bounded-outcome, positive bounded-weight, smoothness, support and
rate premises. They do not inherit unrelated unbounded Gaussian/full-X
alternatives by selecting \code{full_x}. Row/slope laws, transport, covariance
kernel and complete joint sampling premises remain separately required.

Missing complete influence and its source equations/Jacobian, missing required
reference derivatives, any retained
local-polynomial guard/boundary/reference failure, and \code{source_random} ties
return an unavailable inference while retaining the original point. No rows,
models or basis terms are removed, and no missing derivative is set to zero.
Structural misuse, invalid scope/controls and source-binding errors are errors.
Explicitly recognized low-level arithmetic overflows preserve the point with an
unavailable inference; existing transport/critical unavailability is retained.

No bootstrap is run. If the returned low-level result and its selected branch
support replication, use \code{wm_bootstrap(answer$inference, ...)} subject to
that existing method's requirements. This helper adds no bootstrap route or
original-refit validity claim; critical direct-variance-only output does not
authorize replication or a normal interval.
}
\value{
The same model class, original matching fit, nuisance parameters and equations, point values,
recipe, call and solver records, with \code{inference} equal to the exact returned
\code{wm_fitted_inference} object when dispatched. The top-level variance,
standard-error and interval fields mirror that object. Pre-dispatch failures or
recognized arithmetic failures receive explicit unavailable status and NA
interval fields. Population assumptions are never automatically verified.

The active \code{assembly} is NULL to avoid retaining a stale full-X assembly
under a different covariance scope. The first original assembly and inference
are preserved in \code{inference_handoff}; repeated handoffs retain this first
provenance without nesting it. The handoff also records the request, complete
parameter order, whether low-level inference was called, and its correction scope.
When completed, the nuisance stack additionally retains the complete influence,
centered covariance and \code{influence_completion} source/allocation record.
}
\seealso{\code{wm_model_fit}, \code{wm_fitted_inference}, \code{wm_bootstrap}}
