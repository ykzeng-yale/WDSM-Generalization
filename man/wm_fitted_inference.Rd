\name{wm_fitted_inference}
\alias{wm_fitted_inference}
\title{Qualified Fitted-Map Inference with Complete Nuisance Sensitivity}
\usage{
wm_fitted_inference(fit, nuisance_influence = NULL,
  mean_derivative0 = NULL, mean_derivative1 = NULL,
  weight_derivative = NULL, transport0 = NULL, transport1 = NULL,
  covariance_scope, conf.level = 0.95, scalar_control = NULL)
}
\arguments{
  \item{fit}{One unadjusted, bias-corrected, self-normalized
  \code{\link{wm_match}} or \code{\link{wm_fit}} object. Retain original
  outcomes, treatment, weights, predictions and the full unrestricted
  graph. The original point estimate and object are preserved. Point-only
  fits with valid correction means are accepted; their fixed-map variance
  is not used. Generic branches require either all used dimensions at least
  two, or all used dimensions equal to one in the corrected full-X branch
  under at least one sampling-law alternative described below.
  Mixed scalar/higher-dimensional graphs are unsupported.
  For \code{covariance_scope = "scalar_psm"}, supply the raw completed
  \code{\link{wm_scalar_logistic_match}} object with its original fitting
  design binding, preferably obtained with \code{inference = FALSE}.}
  \item{nuisance_influence}{Finite numeric n-by-p matrix estimating one
  joint same-sample parameter influence. Column names are required,
  nonempty and unique. Include scores, predictions and
  empirical scaling in this one parameter order, retaining their covariance.
  Rows are in original observation order. Optional row names must be the
  literal original row indices 1,...,n in order. Required for every generic
  branch; omit it for \code{scalar_psm}, which derives its own influence.}
  \item{mean_derivative0, mean_derivative1}{Finite numeric n-by-p complete
  parameter derivatives of potential-outcome predictions, including the
  chain rule. Columns must have exactly the influence names/order; optional
  row names use original row indices. NULL declares an identically zero
  derivative. PATT requires \code{mean_derivative1 = NULL}.}
  \item{weight_derivative}{Compatibility argument for supplied known weights:
  use NULL or a finite n-by-p exact zero matrix in the influence parameter
  names/order. Nonzero derivatives are rejected before inference because
  estimated-weight uncertainty is outside this WM scope.}
  \item{transport0, transport1}{Explicit donor-direction lists described
  below. PATT uses only \code{transport0}. PATE \code{full_x} requires both
  arguments NULL because that branch supplies zero graph transport. Every
  other used direction requires an explicit estimation or zero declaration;
  missing transport is not interpreted as zero. Scalar PATT requires
  \code{list(mode = "zero", basis = "current_centering")} and correct
  full-X control predictions; scalar estimated or fixed-map declarations
  are not accepted.}
  \item{covariance_scope}{Required branch declaration: for PATE exactly
  \code{"distinct_rarity"}, \code{"full_x"}, \code{"common_field"} or
  \code{"regular_joint_d_gt2"};
  for PATT exactly \code{"patt"}. These are unverified population premises,
  not model-selection or goodness-of-fit labels. When every used matching
  dimension is one, only PATE \code{"full_x"} or PATT \code{"patt"} is
  supported, under at least one stated scalar/full-X sampling-law alternative.
  The separate
  \code{"scalar_psm"} branch accepts PATE or PATT and uses the same
  fitted logistic probability in each used direction.
  \code{"regular_joint_d_gt2"} requires PATE and equal used dimensions d>2
  under the bounded known-weight regular baseline joint-sheet conditions below.}
  \item{conf.level}{Finite confidence level strictly between zero and one.}
  \item{scalar_control}{NULL or a uniquely named list for \code{scalar_psm}
  only. Supported fields are \code{bandwidth_constant} (1),
  \code{density_constant} (0.01), \code{outcome_clip_constant} (1),
  \code{quadrature_error_constant} (0.01), \code{max_nodes} (131072),
  \code{maximum_matrix_entries} (1000000),
  \code{maximum_laplace_terms} (2000000000), optional
  \code{weight_lower_bound}, \code{design} and \code{root_certificate}.
  Numeric controls must be positive; the three computational caps must be
  representable integers, with \code{max_nodes >= 8}. A supplied design
  must be identical to the stored fitting design. A supplied
  \code{\link{wm_scalar_root_certificate}} must bind the exact same
  design, outcomes, treatment, weights, probability and graph.}
}
\description{
Assembles a conditional fitted-map Wald interval using the complete matching
rows, joint nuisance influence, the exact fixed-edge smooth derivative and
the signed incoming graph-transport component. The point estimator is
unchanged. The separate \code{scalar_psm} branch computes the finite-M integrated score drift
and complete contribution rows for raw weighted-logistic PSM. This interface
implements the stated iid fitted-map branches with supplied known weights;
it does not establish their assumptions for an application.
}
\details{
For \code{scalar_psm}, omit \code{nuisance_influence}, all supplied
derivatives and both transport specifications. The bound raw scalar object
must use known positive weights, an unrestricted common probability graph,
fixed M and a full-rank leading-intercept logistic design with p>=2.
The stored nearest donors and row-priority ties are validated without
changing the graph. No propensity or outcome model is refitted.

A compact nonnegative C2 triweight kernel estimates scalar subdensities,
weighted query subdensities, outcome means and their score derivatives.
The bandwidth, density threshold and training-outcome clip are respectively
\eqn{c_h n^{-1/12}}, \eqn{c_\delta n^{-1/256}} and
\eqn{c_Y n^{1/8}}. The default positive constants are declared engineering
choices, not coverage-selected tuning. The raw point and graph are never
trimmed. PATT estimates only the control outcome mean.

Conditional empirical Laplace transforms give
\deqn{A_z(w;s)=w\int_0^\infty e^{-tw}\phi_z(t;s)^{M-1}\,dt,}
\deqn{\partial_s A_z(w;s)=(M-1)w\int_0^\infty
 e^{-tw}\phi_z(t;s)^{M-2}\partial_s\phi_z(t;s)\,dt.}
The partial derivative holds the root weight w fixed. M=1 and constant
donor-weight laws use exact expressions. Other laws use chunked composite
eight-point Gauss-Legendre integration, a mixture-derivative error bound and
an exponential tail bound. These bounds are propagated through the complete
empirical drift; their sum in each parameter must not exceed
\eqn{c_q/\sqrt n}. The initial cutoff is
\eqn{4\log(n)/(M w_{\min})}, where weights are in the same normalized
units as all rows; empirical donor/root minima justify the empirical tail
bound. The optional raw \code{weight_lower_bound} only controls this
initial cutoff and cannot verify population boundedness. Cutoff enlargement
and node refinement obey the recorded finite caps. The bounds concern
ideal real-arithmetic integration and exclude floating-point roundoff,
kernel bias, sampling error and model error.

The scalar drift is the donor-direction residual formula, with the signed
PATE difference or negative control-only PATT normalization. There is no
additional fixed-edge prediction derivative. Feasible means are used only
in the complete inference rows, with retained values clipped at log(n).
The branch computes \eqn{H=\operatorname{mean}\{W e(1-e)XX^{\mathsf T}\}}
and influence \eqn{\ell_i=H^{-1}W_i(Z_i-e_i)X_i}, retaining its full
sandwich covariance. It reports the complete centered-row variance identity
below and an interval around the original point, with total n for both
estimands. Complete centered rows are returned for separately justified
contribution replication; no original refit-bootstrap agreement is asserted.

The scalar population premises require iid rows, known bounded positive
smooth W(Z,X), correct target logistic propensity and overlap, true
own-(score,weight) centering and target identities (control only for PATT),
bounded general covariates with an active continuous coordinate, the stated
globally C2 zero-extended mixed latent subdensities, smooth conditional means,
the weighted-root derivative-span condition and uniform conditional outcome
moments of order 2+eta. These are assumptions, not conclusions of observed
fit checks. The theorem concerns a regular exact logistic root; optional
finite-data root/donor/point certificates do not verify its asymptotic
numerical bridge or these population assumptions. Simultaneously estimated
weights and dependent survey-design inference are outside this branch.

Kernel/conditional-mark work is O(n^2); shared Laplace products are
prospectively capped by \code{maximum_laplace_terms}. Chunk sizes bound
temporary kernel matrices, not total process memory. Empty density trims,
integration-cap failures or other runtime nuisance failures return explicit
unavailability and preserve the raw point. Malformed or mismatched bound
inputs are rejected before nuisance computation. Finite caps do not prove
an asymptotic probability of numerical success.

For the generic corrected branches, including the restricted all-scalar
full-X branch, the following smooth and graph sensitivities apply.

Let \eqn{L_i} be the complete normalized matching row, centered by its
sample mean, and \eqn{U_i} the supplied joint influence centered by its
column means. Write \eqn{q} for the complete fixed-edge ratio derivative.
With supplied weights fixed, it includes the complete prediction derivative
once; donor fractions, outer weights and the target denominator have zero
parameter derivatives. The graph sensitivity is
\deqn{g_A=(j_1-j_0)/\operatorname{mean}(W),\qquad
       g_T=-j_0/\operatorname{mean}(ZW),\qquad b=q+g.}
Here \eqn{j_z} uses the original raw weight units returned by
\code{\link{wm_graph_transport}}. No extra M, n or target-normalizer
factor is introduced. PATT requires no treated-arm mean or transport.

The variance calculation is
\deqn{\widehat V_0=\operatorname{mean}(L_i^2),\quad
  \widehat C=\operatorname{mean}(U_iL_i),\quad
  \widehat\Sigma=\operatorname{mean}(U_iU_i^{\mathsf T}),}
\deqn{\widehat V=\widehat V_0+2b^{\mathsf T}\widehat C+
  b^{\mathsf T}\widehat\Sigma b
  =\operatorname{mean}\{(L_i+U_i^{\mathsf T}b)^2\}.}
All cross-covariances remain. In unnormalized notation
\eqn{C=\Gamma/\gamma}, not \eqn{\Gamma/\gamma^2}. The estimator variance
is \eqn{\widehat V/n} for both PATE and PATT. The positive limiting
variance condition is a population premise; a nonpositive assembled
variance does not produce an interval.

For PATE \code{regular_joint_d_gt2}, both matching dimensions are equal
and strictly greater than two. This bounded iid known-weight branch retains
baseline own-score-and-weight centering and target identities, fixed M,
the complete regular root-n score/prediction/centering/scale stack, bounded
smooth predictions and influence, empirical L2 influence consistency, and
the stated F1--F4 current projection/transport and compact-index laws.
The baseline joint-score law has a finite decomposition into compact
embedded C2 sheets with piecewise smooth Hausdorff-null boundaries,
dimension r>=d and full-rank score projections. Retain regular marginal
geometry and full donor subdensities continuous off Lebesgue-null sets,
bounded arm sheet densities, finite regular atlases with uniformly bounded
coordinate Jacobians and local inverses, and uniform component small-ball
bounds. Distinct d-sheets have Hausdorff-d-null intersections. Coincident
positive-measure sheets are merged with their arm measures and root-mark
mixtures only when a finite smooth decomposition is admissible; other
positive-measure coincidences are outside this branch. The full marginal
conditional donor-weight law, including every sheet and preimage, retains
its stated a.e. Wasserstein-1 continuity. Current sheets need not remain
coincident.
Finite arrays do not verify these population conditions. Estimated
transport retains its separate signed-chart/density and smoothing
premises; this scope does not inherit unbounded full-X alternatives.

The signed reciprocal products are computed from the actual fitted graph,
supplied weights and its own-arm fitted residuals, counting each pair once.
For normalized analysis weights w and gamma=mean(w), let R be their
reciprocal numerator average. The corrected calculation is
\deqn{\widehat V_0=\operatorname{mean}(L_i^2)-2R/\gamma^2,}
\deqn{\widehat V=\operatorname{mean}\{(L_i+U_i^{\mathsf T}b)^2\}
 -2R/\gamma^2=\widehat V_0+2b^{\mathsf T}\widehat C+
 b^{\mathsf T}\widehat\Sigma b.}
The subtraction is in root-n variance units and invariant to common weight
rescaling; there is no additional M factor. Signed R and raw corrected V0
are retained. A nonpositive finite-sample V0 does not gate a positive full
contrast V. Only final V positivity governs usable SE/intervals; otherwise
raw signed fields remain and inference is explicitly unavailable. No floor
or absolute-value correction is applied.

Through \code{wm_bootstrap(..., method = "contribution")}, this scope uses
explicit marginal Gaussian calibration with variance V rather than ordinary
row/count multipliers. Supplied counts are rejected even if the observed
R happens to be zero. The complete C and Sigma blocks remain in the inference
object; no joint nuisance draws, rematching or refitting are performed.
Positive limiting variance and the qualified sampling/plug-in premises
justify this marginal Gaussian law. Fixed B retains Monte Carlo error;
actual sampling-variance convergence still requires separate L2/UI conditions.

The generic all-scalar branch retains the original corrected point and
actual graph, with PATE \code{full_x} or PATT \code{patt}. At least one of
the following complete sampling-law alternatives, or the separate finite-moment
current-graph/Gaussian route below, must be justified.
They are population and fitting premises, not classifications verified
by the supplied arrays. Distinct directional scalar maps are allowed.

Under the design-fitted alternative, the scalar branch requires iid
rows, bounded outcomes and known positive bounded weights \eqn{W=w_Z(X)},
identified targets, correct full-X outcome predictions at the prediction
root, and bounded first two prediction-parameter derivatives. PATT requires
only the control prediction. One complete regular root-n joint influence
includes all used score, prediction and scale parameters; its limit is
bounded and its estimates are empirically L2 consistent. The graph parameter
and its influence component must be measurable in the complete observed
design \eqn{B=(X,Z)}. Prediction parameters may use the same outcomes.
Score-parameter paths through predictions must remain in their derivatives.

For every used direction, X must be reconstructible from a compact monotone
scalar chart \eqn{(U_z,V_z)}. The two arms share one reference measure for
that direction; their joint chart densities are bounded above and away from
zero, continuously differentiable in U with a common derivative bound.
The actual maps are uniformly C2 in U and the graph parameter on a fixed
neighborhood, equal to U at the population parameter, with a positive
bounded U derivative. The same direction map is applied to queries and
donors, with fixed M, replacement, scalar Euclidean distance and original
row ties. Directional charts may differ; no joint density or independence
of the two baseline coordinates is required. Bounded root marks may be Borel
and need no spatial continuity. These are population/design assumptions;
finite arrays cannot certify them.

Under the exact outcome-fitted-root alternative, the graph may use the same
outcomes. The graph parameter solves a closed fixed-dimensional estimating
equation, or determines the graph through an exactly equivalent
parameterization. On a fixed neighborhood its graph estimating function
has its baseline population zero, an invertible population Jacobian,
uniformly bounded first and second
parameter derivatives, and a finite \eqn{2+\delta} score moment for some
\eqn{0<\delta\le2}. The selected statistical root is the unique zero in a
fixed sufficiently small ball when it exists, with a fixed bounded
completion otherwise. This convention is not an arbitrary numerical solver.

The graph equation additionally requires accessible smooth score-block
patches. For graph-root dimension r, there are r row patches whose laws
dominate positive multiples of fixed normalized smooth bumps
\eqn{b_j\in C_c^\infty(I_j)}, with \eqn{\int_{I_j}b_j(t)\,dt=1},
in accessible continuous coordinates on bounded intervals \eqn{I_j}, times
fixed probability laws for the remaining coordinates. The allocation is
independent of the parameter. On the convex product \eqn{\prod_{j=1}^r I_j}
and the declared remaining-coordinate patches, their block-score map is
jointly twice differentiable with uniformly bounded score and needed
coordinate derivatives. The norm of its coordinate Jacobian minus one fixed
invertible matrix is uniformly less than half that matrix's smallest
singular value. These are substantive domination and injective-submersion
conditions; a continuous covariate or a nonsingular estimating Jacobian
alone does not establish them. They give nonsingular graph-root influence
covariance. The covariance of the complete prediction/correction influence
may remain singular and is not inverted by this API.

Each scalar direction uses a common fixed compact enclosing scalar interval
and finite-piecewise chart for its two arm densities and required bounded
joint marks. Endpoints and breakpoints move
Lipschitz continuously in the graph parameter. Within corresponding pieces,
densities are bounded and jointly Lipschitz in score and parameter, and
query/donor densities are mutually comparable. Maps are parameter-Lipschitz;
the conditional laws of the bounded joint weight, prediction-derivative,
variance and mixed covariance marks are uniformly Wasserstein-1 Lipschitz
within pieces. No cross-jump continuity or positive density lower bound is
required. These quantitative conditions are not implied by arbitrary Borel
weights or a qualitative scalar-chart law.

Supplied \eqn{W=w_Z(X)} remains known, positive and bounded. Baseline
predictions equal the bounded true full-X means, with bounded first two
parameter derivatives. One complete regular joint prediction
influence must include every shared or generated score, prediction and scale
coordinate, with empirically L2-consistent estimates. One admissible
construction is a complete fixed-dimensional estimating stack with a unique
interior population zero, finite score second moment, nonsingular square
population Jacobian, square-integrable first-derivative envelope and integrable
Hessian envelope on a fixed compact convex smooth neighborhood. Selected
nearby statistical roots, or exactly smooth generated outputs, give its
full influence. Other fitting influences require their own verified joint
expansion and empirical L2 consistency. Residuals and complete
influences have finite \eqn{2+\delta} moments, and their centered innovations
given \eqn{B=(X,Z)} have uniformly bounded conditional \eqn{2+\delta} moments.
The influence's conditional design mean need only have its finite ordinary
moment. PATT requires only the correct control prediction. The scalar
actual-root selection and actual square/cross/maximum laws follow under
these premises; they are not inferred from a finite covariance identity.

The original outcome-fitted scalar graph must use that exact root, an
exactly graph-equivalent parameterization, or a separately proved numerical
graph transfer. An \eqn{o_p(n^{-1/2})} original parameter error alone does
not establish that transfer. The raw logistic root certificate supplied by
\code{wm_scalar_root_certificate} is not a generic PG or coupled-root
certificate. Numerical availability here concerns the assembled arithmetic;
population and application verification remain false.

Under the applicable stated sampling-law alternative, correct full-X means and
the corresponding actual-graph law give zero graph drift and zero
reciprocal subtraction. PATE supplies zero transport under \code{full_x};
PATT requires the explicit zero/current-centering list. Scalar
\code{distinct_rarity}, \code{common_field}, mixed dimensions and estimated
graph transport remain rejected. Complete rows retain the finite-M
prediction slope and every same-row nuisance covariance block. Positive
limiting variance permits feasible studentization and complete-row
contribution inference. This result does not certify an arbitrary callback,
ordinary rematching bootstrap, growing-B moment conditions or raw sampling
variance convergence. Original same-count prediction refitting has the
separate full-root requirements in \code{wm_bootstrap_refit}; it is not
performed or verified by this function. The raw \code{scalar_psm} branch
retains its separate estimator and premises.

For a d>=2 zero transport declaration use
\code{list(mode = "zero", basis = "fixed_map")} to assert that the
direction's score map is parameter independent, or
\code{list(mode = "zero", basis = "current_centering")} to assert exact
current own-field centering of the baseline residual throughout a parameter
neighborhood. These premises give zero graph transport without estimating
a signed-density chart. No other fields are accepted in zero mode.

For an estimated direction supply a uniquely named list containing
\code{mode = "estimate"}, \code{raw_scores}, \code{tangents} and
\code{cutoff}. Optional fields are \code{raw_to_matching},
\code{bandwidth}, \code{density_floor}, \code{quadrature_tolerance}
and \code{max_nodes}. The last four controls follow
\code{\link{wm_graph_transport}}. Alternate outcome, treatment, weight,
residual, donor-arm and M fields are rejected: these come from the same
validated fit. The residual is the fit's outcome minus that direction's
prediction, and M is its actual fixed donor count.

Raw coordinates must be n-by-d with the same d as the direction's retained
matching coordinates. Tangents are n-by-d-by-p, with third-axis names equal
to the influence parameter names/order. Raw-coordinate column names and
tangent second-axis names must agree exactly, or both be absent. Optional
row labels of coordinates, tangents or cutoff use the original row indices.
For \eqn{R_\theta=A_\theta V+B_\theta L}, the effective raw tangent is
\eqn{(\partial_k B_\theta-(\partial_k A_\theta)A_\theta^{-1}B_\theta)L};
it is generally different from the unadjusted coordinate derivative.

NULL \code{raw_to_matching} uses identity. Otherwise it must be a function
evaluated once on the raw matrix. Its finite n-by-d result must match the
fit's retained scores in original row order within the existing package
engineering tolerance \eqn{10^{-10}\max(1,|x|,|x_{\rm fit}|)} componentwise.
This finite observed check rejects mismatched coordinates but does not prove
population invertibility, the neighborhood chart relation, or the derivative
interpretation. A function that ignores its input is not evidence for those
assumptions. Raw-to-actual coordinate transformations need not preserve the
finite graph; the transport construction uses the invariant first incoming
load moment in the specified raw chart.

The bounded d>=2 own-field alternative requires iid observations; bounded outcomes, residuals, weights
bounded above and away from zero; fixed d, p and M; the point theorem's target
identities and own-(baseline score,baseline weight) centering; bounded smooth
predictions equal to its own-score population functions at the true parameter;
and its regular baseline geometry, bounded current donor densities,
compact-index process modulus and current-row laws. The common parameter has
a root-n asymptotic linear expansion with bounded joint influence and
empirically L2 consistent influence estimates. Prediction/score derivatives
and the regularity of any fitted scale obey that same theorem. Finite arrays
cannot establish these premises. The fitted/baseline law is a distributional
transfer, not same-realization equality of their graph noises.


The correct-full-X alternative permits unbounded outcomes and complete nuisance
influences when every used direction has fixed dimension at least two. PATE
uses \code{covariance_scope = "full_x"}; PATT uses \code{"patt"}, the correct
control prediction and \code{transport0 = list(mode = "zero", basis =
"current_centering")}. Known weights remain positive, bounded functions
\eqn{W=w_Z(X)}. At the population parameter, the used predictions equal
\eqn{E(Y\mid X,Z=z)}; these prediction functions and their first two parameter
derivatives remain uniformly bounded on a neighborhood. PATT requires no
fitted treated prediction.

Each direction retains compact convex full-dimensional baseline donor
support, positive bounded baseline donor density, bounded query density
supported there, uniformly bounded first two parameter derivatives of its
actual Euclidean map, and a uniformly bounded deterministic current donor
subdensity on every root-n parameter ball. These are the F1/F2 and current-density
conditions of the finite-moment full-X theorem. Arm maps may differ; this
alternative needs no joint density of the two maps. Fixed M, row priority,
original donor-normalized fractions and target identification are unchanged.

Let \eqn{B=(X,Z)}, \eqn{\epsilon=Y-E(Y\mid B)}, \eqn{\ell} be the complete
mean-zero joint influence and \eqn{m(B)=E(\ell\mid B)}. For some
\eqn{0<\delta\le2}, require \eqn{E\|m(B)\|^{2+\delta}<\infty} and uniformly
bounded conditional \eqn{2+\delta} moments of \eqn{\epsilon} and
\eqn{\ell-m(B)} given B. An outcome moment alone does not establish this
complete-influence condition. Every fitted score, prediction, pooled centering
and scale belongs to the regular joint root-n expansion; influence estimates
are empirically L2 consistent. A sufficient regular estimating stack has an
isolated population root, nonsingular square Jacobian, square-integrable score
and first-derivative envelope, integrable second-derivative envelope, and the
specified nearby statistical root. Generated coordinates must still obey the
actual-map geometry above; finite solver diagnostics do not establish the
statistical-root or numerical localization conditions.

Under these premises the deterministic-current-graph martingale and
compact-index modulus arguments yield the actual sampling law and actual row square, cross and
maximum laws. The finite-M donor ratio stays inside the prediction-sensitivity
expectation; complete matching/nuisance and nuisance/nuisance covariances are
retained. The full nuisance covariance may be singular, while positive total
variance remains required. The bounded generic own-field and estimated-transport
routes are not extended by replacing their boundedness conditions with these
moments. No application premise, numerical-root rate or original full/partial
refit callback is verified by this interface.

A separate correct-full-X finite-moment current-graph route permits unbounded
supplied weights and predictions when the applicable sampling theorem supplies
the joint oracle/nuisance normal law, actual integrated row square,
matching/influence cross and maximum laws,
the prediction-gradient limit, prediction derivative/Hessian envelopes, and a
complete regular joint influence with empirically L2-consistent estimates.
For the noncompact design-fitted scalar propensity source, its common-design
likelihood-ratio moments, scalar density/marked local laws, VC/tie control and
cap-removal conditions remain substantive application premises. Correct
propensity specification alone does not supply those laws or correct full-X
predictions. This is the corrected full-X branch, not the raw \code{scalar_psm}
estimator or an uncorrected propensity-matching variance reduction.

The complete Gaussian sufficient route assumes iid full design rows B and,
given those rows, independent Gaussian innovations with variances bounded above
and away from zero. Used predictions equal the true full-X means. The complete
regular influence has form \eqn{\ell=m(B)+A(B)\epsilon}, with mean-zero m and
finite global second moments of W, the target heterogeneity row, m and A.
The actual graph, including original donor fractions, must be exactly measurable
given all design rows and the complete Gaussian noise projection. A
design-measurable graph is a special case; replacing the graph parameter by an
influence approximation inside donor rankings is insufficient.

On every fixed root-n parameter ball, the normalized conditional residual
variance and residual/influence cross-variance graph functionals must converge
uniformly to their marked limits; current row-square averages are tight and
the maximum current residual coefficient divided by square root n vanishes.
The actual graph parameter is root-n tight. Square-integrable prediction
derivative and Hessian envelopes, the actual prediction-slope limit and
empirical-L2 influence consistency complete this route. It retains the original
finite-M donor denominators and all fitted score, scale and prediction covariance
blocks, including singular complete nuisance covariance. It does not require a
uniform conditional moment of the complete influence. Positive limiting total
variance permits complete-row studentization and contribution replication.
Verified source double-score matching and full-covariate matching with complete
quadratic predictions retain their applicable score-chart/minor and interaction
conditions. Numerical
statistical-root/graph equivalence and actual refit callbacks require their
separate conditions; finite solver diagnostics do not establish them.

PATE \code{distinct_rarity} uses baseline maps \eqn{(A,B_0)} and
\eqn{(A,B_1)} after the allowed fixed regular coordinate changes, with a
bounded full-dimensional joint density and
\eqn{\dim(A)<\max(d_0,d_1)}; the actual-index reciprocal-rarity theorem
applies. Unequal fixed d are allowed when both are at least two.
\code{full_x} requires baseline residual centering given all covariates and
treatment, with baseline weights a function of these; it gives exact current
centering and zero graph transport without a joint-map density condition.
\code{common_field} requires the two current maps to generate the same field
for every parameter in a neighborhood, for example invertible transforms of
one common current map. This stated branch uses equal dimensions at least two;
merely observing similar score matrices is insufficient. PATT needs only its
control-direction conditions and no two-arm reciprocal condition. The
assembly uses the proved diagonal block in all these branches, without a
fixed-map reciprocal sample subtraction.

Estimated transport additionally requires the signed chart, compact interior
signed support, smooth cutoff equal to one on that support, ordinary marked
C3 density assumptions and bandwidth/rate conditions documented in
\code{\link{wm_graph_transport}}. These conditions are not implied by a
label such as survey or balancing weights. For statistical consistency the
propagated quadrature bound must vanish; the fixed default numerical tolerance
alone does not guarantee this. Returned quadrature bounds exclude sampling
error, kernel bias, nuisance error and floating-point roundoff.

Positive supplied weights define the target and remain fixed throughout
score, prediction and scale estimation. They may depend on treatment and
covariates. Each sampling alternative retains its stated weight bounds or
integrated moment conditions; the bounded own-field and estimated-transport
routes are not extended by the separate full-X/Gaussian alternatives.
This interface does not include weight-estimation parameters or
their uncertainty in the joint influence stack. NULL or an exact zero
\code{weight_derivative} matrix is required. The same known-weight restriction
applies to the separate \code{scalar_psm} branch, which requires that argument to be omitted.

Malformed data/specifications and unsupported branches raise errors. After
those guards, a numerical density/quadrature failure retains the point fit and
completed components with \code{status = "transport_unavailable"}, an explicit
reason, and unavailable variance/interval. There is no zero fallback. An
all-zero active cutoff is only the computational zero marked by the component;
it cannot establish that a nonempty signed support was handled correctly.

This interface does not fit a model, rematch at a perturbed parameter, bootstrap,
or mutate the supplied fit. Generic branches reject raw estimators and
unsupported scalar maps; \code{scalar_psm} requires its bound raw logistic object.
All branches reject unsupported previous adjustments, stabilized rules,
restricted/fold graphs and covariance specifications. It makes no generic learner,
unrestricted nonzero-reciprocal fitted-law or cluster/stratum design claim, and does not
certify an actual original-refit callback or raw sampling-variance convergence.
}
\value{
For \code{regular_joint_d_gt2}, additional fields retain
\code{reciprocal_numerator}, the once-oriented \code{reciprocal_pairs},
\code{diagonal_root_n_variance}, \code{augmented_row_root_n_variance},
\code{reciprocal_subtraction} in root-n variance units,
\code{unregularized_root_n_variance}, \code{unregularized_variance},
\code{analysis_target_weight_mean}, \code{weight_scale} and explicit
unfloored diagnostics. Its \code{V0} is the raw reciprocal-corrected oracle
block, including nonpositive values; unavailable usable variance fields are NA.

A separate \code{wm_fitted_inference} list preserving \code{estimate},
\code{fit}, n, M and estimand. On successful assembly it contains:
\describe{
  \item{smooth_sensitivity, graph_sensitivity, total_sensitivity}{Complete
  smooth slope, signed normalized graph slope and their sum in estimator units
  for generic corrected branches. The separate \code{scalar_psm} branch supplies \code{total_sensitivity} directly
  from its integrated residual formula.}
  \item{V0, C, Sigma, cross_term, nuisance_variance}{The complete covariance
  blocks and their contributions to root-n variance.}
  \item{root_n_variance, variance, se, conf.int}{Assembled root-n variance,
  estimator variance and conditional Wald interval.}
  \item{base_rows, augmented_rows, nuisance_influence}{Centered complete
  rows, their joint nuisance augmentation and centered influence.}
  \item{smooth_derivative, graph_transport}{Reviewed finite-edge derivative
  details and each used direction's component or explicit zero premise for
  generic corrected branches. Separate \code{scalar_psm} smoothing and donor-direction drift details are in
  \code{completed_components}.}
  \item{graph_sensitivity_quadrature_error_bound}{Componentwise propagated
  quadrature-only bound after target normalization.}
  \item{drift_error_bound, drift_tail_error_bound, controls, active,
  feasible_means, H, laplace_terms}{Separate \code{scalar_psm} total integration and tail
  bounds, recorded controls, common density trim, variance-only mean predictions,
  weighted-logistic derivative matrix and prospective Laplace work accounting.}
  \item{numerical_certificate, observed_logistic_check}{Separate \code{scalar_psm} exact
  optional certificate binding and observed score/nearest-donor checks.
  Neither verifies a population premise or asymptotic numerical-root rate.}
  \item{status, available, conditional_inference_available}{Positive variance
  gives \code{conditional_<branch>} and available TRUE, strictly under the
  unverified recorded premises. Otherwise the status/reason identifies the
  numerical failure or nonpositive variance; no interval is available.}
  \item{assumptions_verified, application_verified}{Both FALSE, including
  when a conditional interval is returned.}
  \item{contract, dimensions, parameter_names, covariance_scope}{Explicit
  model/theorem premises and dimensions/parameter order.}
}
For transport or scalar nuisance failure, covariance/row fields have not been assembled; the
returned fit, smooth derivative, completed components and failure stage/reason
remain available. Variance, standard error and interval are NA.
}
\seealso{\code{\link{wm_graph_transport}}, \code{\link{wm_match}},
  \code{\link{wm_reciprocal_inference}}, \code{\link{wm_scalar_logistic_match}},
  \code{\link{wm_scalar_root_certificate}}}
