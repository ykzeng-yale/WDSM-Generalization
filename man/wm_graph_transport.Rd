\name{wm_graph_transport}
\alias{wm_graph_transport}
\title{Estimate an Incoming Graph-Transport Component}
\usage{
wm_graph_transport(raw_scores, Z, W, residual, tangents, cutoff,
  donor_arm, M, bandwidth = NULL, density_floor = NULL,
  quadrature_tolerance = 1e-7, max_nodes = 131073L)
}
\arguments{
  \item{raw_scores}{Finite numeric \eqn{n} by \eqn{d} matrix of the
  fitted raw chart for one donor direction, evaluated on every original
  observation. Supply the chart used in the signed-transport model;
  the function does not recover it from standardized matching scores.}
  \item{Z}{Numeric length-\eqn{n} vector with values zero and one. Both
  arms must be represented. Rows agree with all other arguments.}
  \item{W}{Positive finite numeric length-\eqn{n} weights in their original
  raw units, treated as supplied known weights. Statistical consistency requires population weights bounded
  above and away from zero; sample positivity does not verify that premise.}
  \item{residual}{Finite length-\eqn{n} vector, equal on donor rows to
  \eqn{Y_i-f_z(X_i;\widehat\theta)}. Other rows are ignored in the final
  drift average but must be finite.}
  \item{tangents}{Finite numeric array of shape \eqn{n} by \eqn{d} by
  \eqn{p}. The third axis indexes the common nuisance parameter, not the
  matching-coordinate dimension. Its entries are the effective raw-chart
  tangents specified below. Optional third-axis names must be unique and
  nonempty. Unnamed parameters are interpreted by position.}
  \item{cutoff}{Finite length-\eqn{n} values in \eqn{[0,1]} of a fixed
  smooth interior cutoff evaluated at the supplied raw coordinates.
  It must equal one on the signed-density support. Arbitrary trimming or
  an inferred support cutoff does not have the stated interpretation.}
  \item{donor_arm}{Zero or one. All opposite-arm observations are queries.
  For PATT use only donor arm zero.}
  \item{M}{Positive integer fixed donor count. There must be at least
  \eqn{M} donor-arm observations and one opposite-arm observation.}
  \item{bandwidth}{Positive scalar bandwidth in raw-coordinate units for
  the product kernel. The default is \eqn{n^{-1/[2(d+4)]}}. No automatic
  rescaling or bandwidth selection is performed.}
  \item{density_floor}{Positive donor-subdensity threshold. The default
  is \eqn{0.01 n^{-1/4}}. A density at or below this threshold at any
  active donor raises an error; the function does not silently trim it.}
  \item{quadrature_tolerance}{Positive componentwise absolute tolerance
  for the tagged fraction and each spatial derivative. The propagated
  bound on the final drift is returned separately and may differ in size.
  A fixed tolerance is not an asymptotic error guarantee.}
  \item{max_nodes}{Maximum number of Simpson quadrature nodes per active
  donor, at least three. Exceeding the budget raises an error.}
}
\description{
Estimates the unnormalized incoming graph-transport vector for the
original self-normalized donor rule, retaining the random finite-M donor
denominator. The calculation uses a conditional empirical Laplace
transform and deterministic quadrature; it makes no random draws.
It produces a variance component input, not a variance or confidence interval.
}
\details{
For a raw chart \eqn{R_\theta=A_\theta V+B_\theta L}, with invertible
\eqn{A_\theta}, the effective tangent for parameter \eqn{k} is
\deqn{v_k(X)=\{\partial_k B_\theta-
  (\partial_k A_\theta)A_\theta^{-1}B_\theta\}L.}
Evaluate this expression at the fitted parameter. The within-chart
reparameterization is removed: \eqn{\partial_k R_\theta} alone is generally
not the required tangent. The actual matching map may be an invertible
transformation of the raw chart. The invariant first incoming moment can
be estimated in raw coordinates; this does not assert invariance of a
finite matching graph under that transformation.

Write \eqn{H_z(s)} for the donor raw-score subdensity and
\eqn{H_{1-z}^W(s)} for the weighted query subdensity. For the conditional
donor weight law at \eqn{s}, define
\deqn{A_M(w;s)=E\left\{\frac{w}{w+\sum_{r=1}^{M-1}W_r}
  \mathrel{\Big|} Z=z,R=s\right\},\qquad
  k_z(s,w)=M\frac{H_{1-z}^W(s)}{H_z(s)}A_M(w;s).}
The other \eqn{M-1} weights in this expectation are independent draws
from that conditional donor law. At \eqn{M=1}, \eqn{A_1=1}, but the
spatial density-ratio derivative can still be nonzero. Replacing the
finite donor sum by its expectation changes this formula in general.

The returned component is
\deqn{\widehat j_{z,k}=\frac1n\sum_{i:Z_i=z}
  \widehat\varepsilon_i\chi(\widehat R_i)
  \widehat v_{k,i}^{\mathsf T}
  \nabla_s\widehat k_z(\widehat R_i,W_i).}
The positive sign follows from integration by parts under compact
interior signed support. The kernel is the product of
\eqn{(35/32)(1-u^2)^3 1(|u|<1)}. Conditional donor weights are the
kernel-weighted empirical law; there is no kernel in the weight coordinate.
Discrete, continuous and mixed weights are computationally supported.
The tagged fraction and its spatial derivatives are integrated jointly
by a one-dimensional Laplace formula. M=1 and a constant local donor
weight law are exact branches. Otherwise a tail bound and composite
Simpson fourth-derivative bound are propagated through the load gradient
and drift average. These bounds exclude floating-point roundoff.

For PATE, compute each arm's component in its own raw chart and use
\eqn{(\widehat j_1-\widehat j_0)/\operatorname{mean}(W)}.
For PATT use \eqn{-\widehat j_0/\operatorname{mean}(ZW)}.
The numerator components and denominator must use the same weight scale.
There is no extra M or n factor, or second division by the target mean.
Then add the complete fixed-edge prediction derivative once, with the
supplied weights unchanged. Adding prediction sensitivity again would
double count it. Parameter order must agree across arms, prediction
derivatives and nuisance influence columns.

The statistical interpretation requires iid rows, fixed d, p and M,
the point theorem's identification and geometric assumptions, bounded
residuals and tangents, and own-(raw score,weight) residual centering.
The signed residual-tangent density must have compact support strictly
inside the positive donor-density region. The cutoff equals one there
and has compact support within a larger regular interior neighborhood.
Boundary transport is not estimated by arbitrary cutoff values.

An additional ordinary marked-density condition is necessary. On that
neighborhood, the joint arm/raw-score/weight measures must have densities
\eqn{h_a(s,w)} relative to spatial Lebesgue measure times a fixed finite
dominating weight measure. Spatial derivatives through order three
exist and are continuous, with integrable uniform bounds over weights;
the donor marginal subdensity is bounded away from zero. Finite weight
support is a special case. These premises do not follow merely from
positive bounded weights, including weights depending on treatment and
covariates. Singular conditional charts can fail them.

For uniformly root-n consistent generated coordinates, residuals and
effective tangents, a sufficient bandwidth choice is
\eqn{b_n=n^{-\alpha}}, \eqn{0<\alpha<1/[2(d+2)]}. The uniform gradient
error has the conservative order
\deqn{\Delta_n=b_n^2+
  \sqrt{\frac{\log n}{n b_n^{d+2}}}+
  \frac{\log n}{n b_n^{d+1}}+
  n^{-1/2}b_n^{-(d+2)}.}
With a vanishing propagated numerical bound, the drift is consistent.
A separate uniformly consistent mean estimate adds its evaluated error.
The supplied weights remain unchanged. The Laplace construction estimates
the graph term for fitted coordinates; it does not supply an estimated-weight
influence correction or a theorem accounting for weight estimation.

The component formula accepts any fixed positive d. Complete generic
same-sample fitted-map inference requires the separately established
branch with all used dimensions at least two, a valid covariance branch,
joint nuisance influence and consistent complete sensitivities. Scalar
fitted inference needs its separate model-specific theorem. In particular,
successful computation does not establish a generic fitted-PS result,
nonzero-reciprocal covariance treatment, cluster or stratified-design
inference, or validity of the original nuisance-refit bootstrap.

This function does not modify a matching fit. The existing
\code{\link{wm_prediction_adjust}} full-X branch has zero graph transport;
nonzero transport cannot be supplied to that interface by relabeling the
centering assumption. A fixed-map reciprocal interval is likewise not
a fitted-map interval. A separate complete-row variance assembly must
retain all matching/nuisance cross covariance under its proved branch.

The computation evaluates each active donor against all observations.
Quadrature temporary arrays are chunked, but total work can be quadratic
in sample size and grows with the node count. No large-sample runtime
guarantee is supplied.
}
\value{
An object of class \code{wm_graph_transport}. The numerical fields are:
\describe{
  \item{graph_drift, quadrature_error_bound}{Length-p unnormalized
  component and its componentwise quadrature-only absolute bound.}
  \item{evaluated_donor_rows}{Original row indices where donor-arm cutoff
  values are positive.}
  \item{incoming_load, incoming_load_gradient, gradient_error_bound}{First
  incoming weighted-load estimates and their spatial gradients/bounds,
  in the order of \code{evaluated_donor_rows}.}
  \item{donor_subdensity, quadrature_nodes}{Density and node diagnostics
  for each active donor.}
  \item{n, d, M, donor_arm, bandwidth, density_floor}{Used dimensions,
  donor count, direction and smoothing controls.}
  \item{status}{\code{component_computed}, or \code{no_active_donors} for
  the computational zero from an all-zero donor cutoff. The latter does
  not establish that population graph transport is zero.}
  \item{numerically_available, sampling_inference_available}{A successful
  return has TRUE and FALSE respectively. Numerical or input failures
  raise errors; no unavailable component is silently replaced by zero.}
  \item{assumptions_verified, roundoff_included}{Both FALSE.}
  \item{parameter_names, parameter_order}{Supplied third-axis names, if
  present, and named or positional order.}
  \item{weight_units, normalization, bound_scope, scope, contract}{Explicit
  units, assembly instructions and unverified scientific/numerical premises.}
}
}
\seealso{\code{\link{wm_match}}, \code{\link{wm_prediction_adjust}},
  \code{\link{wm_reciprocal_inference}}}
