\name{wm_model_inference}
\alias{wm_model_inference}
\title{Qualified Inference for a Retained Finite-Model Matching Fit}
\usage{
wm_model_inference(object, covariance_scope,
  transport0 = NULL, transport1 = NULL,
  critical_control = NULL, conf.level = 0.95,
  max_influence_elements = 5e7, quadratic_control = NULL)
}
\arguments{
\item{object}{A retained \code{wm_model_fit} or \code{wm_wdsm_fit} result. The complete nuisance influence
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
\item{quadratic_control}{NULL retains the existing handoff. For a
\code{quadratic_qr} point, an explicit list with
\code{method="nested_contrast"}, \code{small_ps} and \code{large_ps} names
prepares observable candidate-law inputs for subsequent
\code{wm_bootstrap(answer$inference, ...)}. Optional
\code{maximum_draw_entries = 2000000} and
\code{maximum_distance_evaluations = 100000000} are positive computational
allocation/work caps checked before RNG. The former bounds planned retained
numeric entries; it excludes R/solver overhead and is not an RSS bound. The latter
bounds all requested local-graph distance evaluations. Fixed caps do not establish
asymptotic success; callers must provide adequate resources for predeclared n and B.
They do not filter draws or impose a statistical cutoff. The obsolete \code{guard_exponent}
field is rejected, including an explicitly supplied NULL. It requires PATE
\code{full_x} or PATT \code{patt},
exactly two probability-weighted PS models with identical shared design columns
and fixed offsets, and one PG per used arm. Omit transport and critical inputs.
Preparation alone does not compute sampling variance, draws or an interval.}
}
\description{
Binds the original matching fit, full nuisance influence, complete
prediction/reference derivatives and exact-zero known-weight derivative to
\code{wm_fitted_inference}. It performs no model, local-regression, graph or
weight refitting. The covariance and transport premises remain caller declarations.
}
\details{
The default quadratic route uses its complete retained prediction derivatives and
nuisance influence, including correction coefficients. The local-polynomial
route uses its complete score-only influence and the plug-in reference derivative
\eqn{B^T\nabla\widehat m}. This is not differentiation through refitting local
coefficients. These different parameter layouts must never be mixed.

The \code{quadratic_qr} route retains the same complete quadratic WLS point
and a score-only root. Its correction-projection contribution and complete
joint law are not implemented by the regular handoff. Without explicit
\code{quadratic_control}, it returns
\code{model_correction_law_unavailable} without fitting, deleting
coordinates, supplying zero derivatives or dispatching an existing variance
formula. Numerical stability of the point fit does not verify sampling inference.

With explicit \code{quadratic_control}, the same API retains the original
point/root and prepares \code{inference$contrast_inputs}. Its
\code{evaluate(contrast)} function recomputes pooled standardization and
local-index graphs from the retained smaller fit and the declared added-parameter
perturbation. It fits no PS, PG or quadratic outcome model. Matching uses every
original coordinate and the same metric, M, donor-normalization rule and supplied W.
The reported point remains the original complete-quadratic fit.

The evaluator retains all ten quadratic terms, total-n donor-weighted projections,
signed PATE/PATT functionals, correction-projection rows, complete score influence
and full sensitivity. Its \code{conditional_mean} and
\code{conditional_variance} evaluate an algebraically exact representation of
the conditional Gaussian parameters, subject to floating-point error.
Full-column donor-weighted QR replaces normal equations. The basis uses
standardized smaller PS, standardized arm PG and the exact finite probability
divided difference evaluated with sign-stable exponential identities. It spans
the same complete quadratic space for each nonzero contrast; no limiting
directional approximation replaces the finite function. Matching continues
to use the original three standardized scores.

A full QR of the contrast influence, scaled by the square root of total n,
computes the conditional mean and the residual sum of squares without
inverting its covariance or subtracting variances. Singular complete score
covariance is permitted. This single-index conditional variance is not the
marginal sampling variance: the latter also contains variation in the conditional
mean over the Gaussian contrast law.

There is no radial cutoff or statistical Gram threshold. Full-column pivoting
does not delete a basis term or contrast direction. The numerical QR resolution
check requires the triangular factor's reciprocal condition estimate to exceed
\code{64 * .Machine$double.eps * max(nrow(A), ncol(A))}. This is a
machine-resolution diagnostic, not an asymptotic Gram guard or an error
certificate. An exact zero contrast, unresolved numerical rank, zero score
scale or nonrepresentable arithmetic returns unavailable evaluated inputs,
retaining the requested contrast. No failed index is silently replaced by
zero, redrawn or deleted. Returned QR reconstruction and orthogonality
diagnostics do not verify forward prediction error or simultaneous all-draw
accuracy; \code{numerical_error_rate_verified} remains FALSE. Numerical availability does not
establish correct candidate models, identification, bounded support/density,
sampling/input limits or a numerical prediction-error rate. Preparation provides
no interval, sampling variance or draws. Its \code{replication_available} flag
permits explicit \code{wm_bootstrap} dispatch to the scoped untruncated
Gaussian-mixture law. That call retains every requested contrast/innovation pair;
any failed slot makes all aggregate variance, standard-error and interval fields
unavailable. No runtime zero-completed bootstrap is supported.

This route requires the bounded outcome-correct nested weighted-logistic PS and
arm-specific linear PG sampling/input conditions. It does not admit arbitrary
model-union correctness or every singular multiple-score stack. The complete
nuisance covariance remains shared across matching and correction components.
The same-draw numerical RMS requirement and the original point's root-n numerical
error requirement remain separate unverified premises; QR diagnostics do not
certify either. See \code{wm_bootstrap} for growing-B variance and quantile
requirements and the distinction from full finite-data conditional variance.

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
that method's requirements. This helper makes no original-refit validity claim;
critical direct-variance-only output does not
authorize replication or a normal interval.
}
\value{
The same model class, original matching fit, nuisance parameters and equations, point values,
recipe, call and solver records, with \code{inference} equal to the exact returned
\code{wm_fitted_inference} object when dispatched. The top-level variance,
standard-error and interval fields mirror that object. Pre-dispatch failures or
recognized arithmetic failures receive explicit unavailable status and NA
interval fields. Population assumptions are never automatically verified.

For explicit nested quadratic preparation, status is
\code{model_quadratic_contrast_inputs_ready},
\code{conditional_law_inputs_available} and \code{numerically_available} are
TRUE, \code{replication_available} is TRUE, and
\code{sampling_inference_available} remains FALSE until a separate replication
call. The
\code{contrast_inputs} record contains the deterministic evaluator, full and
contrast parameter orders/covariances, information map, complete QR factors,
numerical-resolution diagnostics and allocation record. Evaluated projections
retain the complete basis, signed query/donor functional, triangular factors,
coefficients, predictions and projected derivative; no Gram inverse is formed.
Top-level sampling variance/interval fields remain NA. This is also
retained at \code{inference_inputs$quadratic_contrast}.

The active \code{assembly} is NULL to avoid retaining a stale full-X assembly
under a different covariance scope. The first original assembly and inference
are preserved in \code{inference_handoff}; repeated handoffs retain this first
provenance without nesting it. The handoff also records the request, complete
parameter order, whether low-level inference was called, and its correction scope.
When completed, the nuisance stack additionally retains the complete influence,
centered covariance and \code{influence_completion} source/allocation record.
}
\seealso{\code{wm_model_fit}, \code{wm_fitted_inference}, \code{wm_bootstrap}}
