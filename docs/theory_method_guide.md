# Weighted Matching: theory and method guide

The current framework treats supplied `W=w_Z(X)`, including `w(X)`, as known.
For iid observed rows with law Q, the target law is `dP/dQ=W/E_Q(W)`.
Identification of a causal PATE or PATT remains a separate requirement.
Each query selects fixed `M` opposite-arm donors, with replacement across
queries. The original donor share is `W_j/sum(W_selected)`; its target,
query weights and donor normalization remain in every original-rule API.

The matching dimension `d_z` and nuisance dimension `p` differ. `p` contains
every coefficient, center, scale and finite-dimensional prediction/correction
parameter in the declared root, in one named order. Auxiliary spline/local-
polynomial coefficients need separate value, gradient and current-graph/refit
transfer bounds; they do not enlarge this finite score root. All matching--nuisance
and nuisance cross-covariances must remain. PATT omits unused treated
prediction blocks. Root-n variance `V` is divided by total observed `n` for
both targets.

## Main method routes

| Route and public APIs | Correction and complete nuisance inputs | Variance or replication operator | Required conditions |
| --- | --- | --- | --- |
| Supplied maps: `wm_match`, `wm_fit` | Supply the intended scalar, distinct-arm or full-X coordinates and their scaling. Means are supplied or fitted by the documented polynomial/spline rule. This fixed-map interface does not create a fitted-score joint `p`; its feasible prediction remainder still needs justification. | The ordinary row variance/contribution bootstrap requires its strong graph-field centering. For distinct fixed maps, `wm_reciprocal_inference` retains the signed reciprocal-pair term; PATT has no reciprocal term. | Known maps and W, iid sampling, identification, branch-specific own/graph-field residual centering, donor geometry, moments and negligible feasible correction error. The fixed-map reciprocal implementation has its separate bounded contract. |
| Corrected full-X PS, DSM or X: generic `wm_match` plus `wm_fitted_inference`; fitted DSM through `wm_wdsm_fit`; finite score lists through `wm_model_fit` | Supply correct full-X conditional-mean predictions and complete joint influence/prediction derivatives. A quadratic in reduced scores is full-X correct only when the mean belongs to that space. `wm_wdsm_fit` retains PS, arm PG, pooled centers/variances and complete quadratic correction. | PATE `covariance_scope="full_x"`; PATT `"patt"` with the documented zero transport. Complete augmented rows estimate `V0+2*b'C+b'Sigma*b`; `wm_bootstrap` contribution replication uses those rows. | The applicable design-fitted, outcome-fitted or Gaussian actual-graph/root theorem, correct full-X predictions, current-graph row/slope laws and moments. Scalar maps retain separate chart and exact-root/numerical-transfer conditions. A fitted or correct PS alone does not establish this route. |
| Raw scalar PSM: `wm_match` with no point correction, then `wm_fitted_inference(..., covariance_scope="scalar_psm")`; bounded logistic wrapper `wm_scalar_logistic_match` | Means/smoothers estimate variance moments only. `p` counts every propensity-family coefficient. Configure the declared probability-family descriptor with its full derivative; this branch derives the supplied-weight Bernoulli influence internally. Omit generic influence/derivative/transport inputs; retain the finite-M weighted drift and full covariance. | Explicit moment variance; supported complete-row contribution replication through `wm_bootstrap` or the narrower `wm_scalar_replication` interface. The two finite-data variance estimates need not coincide. | Correct target PS family, known-weight compatibility/centering, specified scalar inverse chart, accessible derivatives, moments/smoothing rates, consistent exact root and numerical bridge. Generic smooth families require their actual descriptor and premises; the logistic wrapper has a narrower bounded model. This is distinct from corrected full-X PS matching. |
| Regular fitted transport: `wm_graph_transport`, `wm_fitted_inference`, or retained `wm_model_inference` | Complete `p`, prediction imbalance and arm-specific graph transport. Supply the genuine donor-weight chart, or the separate smooth score-measurable weight representation and matching-coordinate tangents. | Assemble the full drift, matching/nuisance covariance and sandwich term. The declared covariance scope determines any reciprocal subtraction. The transport component alone supplies no variance or interval. | Bounded iid known-W fitted-stack conditions, baseline and current geometry, centering/target identities, chart/derivative validity and consistent plug-in moments. `regular_joint_d_gt2` needs equal d>2 and its complete joint-component conditions. Transport modes and covariance scopes are caller premises. |
| Original fixed reuse: `wm_bootstrap(..., method="fixed_reuse")` or `wm_bootstrap_refit` | Keep original W, donors and incoming loads. A callback refits the intended predictions using the same full-n counts in every active block; complete `p` is that of the original system. A proved compatible prediction-only callback is separate. | Nonlinear count replicates retain the random PATE/PATT denominator and divisor B. Contribution adjustments are not appended to these replicates. | Nearby same-count roots, actual graph-row/refit expansions, simultaneous numerical success/error control and sampling/replication-law agreement. Returning predictions or obtaining convergence does not verify these conditions. `refit=NULL` holds predictions fixed and does not recover omitted fitted-score uncertainty. |

`wm_model_fit` constructs fitted PS/PG lists; it is not a raw-covariate fitting
wrapper. The public [principal-study recipe](../simulations/principal_iid/README.md)
shows corrected PS, distinct-arm DSM and full-X matching in one aligned
experiment. Its complete nuisance dimensions are 26/18, 42/26 and 68/40 for
PATE/PATT, respectively. Those counts describe that recipe, not every model.
The generic fitted-inference interface supports all-scalar used maps or used
maps all of dimension at least two under their respective contracts; mixed
scalar/vector maps remain point-only in that interface.

For fixed B, original divisor-B variance retains Monte Carlo variation; it is
not a consistent variance estimate merely because n grows. The separate
polynomial-growing-B result needs its complete regular finite-p, envelope and
all-draw conditions. The count-refitted local-polynomial result does not
automatically inherit that finite-p growing-B statement.

## Numerical and scope conditions

The WDSM, model-list and legacy propensity fitters require every fitted PS
to be strictly inside
`(.Machine$double.eps, 1-.Machine$double.eps)`. This fixed acceptance rule can
reject a finite regular root. For an unbounded Gaussian propensity index its
pass probability tends to zero along a consistent root sequence. Exact-root
statistical results therefore do not establish eventual availability of this
finite-data double-precision implementation. Root localization, numerical
parameter/effect error, original graph transfer and success over prescribed
refits require separate justification. No clipping, distance change, term
deletion or silent redraw is used to repair a failure.

The `critical_planar` reciprocal-function mixture/direct-variance calculations,
outcome-fitted scalar exact-root alternative, local-polynomial correction and
singular nested-quadratic contrast preparation have their own substantive
conditions and numerical/work contracts. They are explicit advanced branches,
not empirical defaults for arbitrary fitted scores. Selecting a branch does
not verify its assumptions, turn a non-Gaussian limit into a normal pivot, or
establish original-refit validity. Estimated weights, growing matching
dimension and dependent cluster/stratum inference remain outside the current
known-weight framework.

Read the R help for the selected API before interpreting an interval. Available
numerical output, implementation reductions and simulation coverage are
distinct evidence; none certifies all population assumptions.
