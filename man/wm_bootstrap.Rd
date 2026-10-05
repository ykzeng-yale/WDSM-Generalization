\name{wm_bootstrap}
\alias{wm_bootstrap}
\title{Bootstrap Replication for Weighted Matching}
\usage{
wm_bootstrap(object, B = 999L, seed = NULL, conf.level = 0.95,
  interval = c("normal", "basic", "none"), chunk_size = 65536L,
  counts = NULL, method = c("contribution", "fixed_reuse"), refit = NULL)
}
\arguments{
\item{object}{For contribution replication, a \code{wm_match} result with variance estimation enabled,
or a successful \code{\link{wm_fitted_inference}} result with available complete
augmented rows in \code{scalar_psm} or a qualified zero-reciprocal scope:
for d>=2, \code{distinct_rarity}, \code{full_x}, \code{common_field} or
\code{patt}; for all-scalar corrected full-X graphs, PATE
\code{full_x} or PATT \code{patt} with declared zero/current-centering,
under at least one stated scalar sampling-law alternative.
PATE \code{regular_joint_d_gt2} fitted inference instead supplies its
complete corrected marginal Gaussian variance; counts are unsupported in that scope.
PATE \code{critical_planar} instead supplies a complete nuisance/Schur-mixture
law through this same contribution interface, with counts unsupported.
A prepared \code{model_quadratic_contrast_inputs_ready} inference object instead
supplies the scoped nested quadratic Gaussian-mixture law described below.
Mixed scalar/higher-dimensional graphs and other generic scalar scopes
are unsupported.
For fixed-reuse replication, an ordinary self-normalized \code{wm_match}
result with supplied matching maps of any dimension. Retained
\code{wm_wdsm_fit} and \code{wm_model_fit} wrappers are also supported:
contribution replication requires their declared full-X or explicitly prepared
public inference, while fixed-reuse replication uses the retained original fit.}
\item{B}{Positive integer number of draws, at least two for fixed-reuse replication.
For supplied counts, B is
inferred from their number of columns when omitted. An explicit B must equal
the count column number.}
\item{seed}{Optional nonnegative integer seed. An explicit seed restores the
caller's \code{.Random.seed}, including its prior absence. The Box--Muller
normal generator does not support this option because its hidden cached state
cannot be restored. Use \code{seed = NULL} with that generator or select another
normal generator before calling this function.}
\item{conf.level}{Confidence level strictly between zero and one.}
\item{interval}{For contribution replication, \code{"normal"} uses the analytic
conditional contribution variance and \code{"basic"} uses empirical root-n draw
quantiles. For fixed-reuse replication, \code{"normal"} uses the original
B-divisor replicate variance; \code{"basic"} is unsupported.
\code{"none"} omits the interval under either method. For \code{critical_planar}
and the nested quadratic mixture, an omitted interval selects \code{"basic"};
explicitly requested \code{"normal"}
is rejected because this law is generally non-Gaussian.}
\item{chunk_size}{Positive integer maximum number of multiplier values
generated together for contribution replication. Ignored for fixed-reuse
replication. For regular_joint_d_gt2 this bounds scalar Gaussian draw generation.
For the nested quadratic law this bounds normal generation within each slot;
contrast normals precede its scalar innovation, independent of chunk size.
No contributions-by-B matrix is allocated.}
\item{counts}{Required for fixed-reuse replication; optional for fitted
contribution replication in zero-reciprocal/scalar_psm scopes only.
\code{regular_joint_d_gt2} and nested quadratic replication reject supplied counts. Entries must be finite nonnegative integers; every column must
sum to full observed n. Rows correspond to original observations; optional
row names must be literal indices 1,...,n in order. Do not supply seed with
frozen counts. Fixed-reuse replication requires positive multiplicity in both
arms in every column; complete-row contribution replication allows collapsed
columns. The array alone cannot verify its multinomial draw law.}
\item{method}{\code{"contribution"} (default) or \code{"fixed_reuse"}.}
\item{refit}{Optional function of one full-n multiplicity vector returning a
list of mean0 and, for PATE, mean1 predictions. Only supported with
\code{method = "fixed_reuse"}. For a \code{wm_match} or \code{wm_fit}, NULL
retains the original predictions. For a fitted \code{wm_wdsm_fit} or
\code{wm_model_fit} wrapper, NULL refits the retained default complete quadratic
prediction pipeline using each supplied count vector. Other correction methods
require an explicit compatible callback.
The callback performs any intended nuisance refitting without changing the
original graph, individual weights or weighted reuse.}
}
\description{
Provides contribution multipliers and fixed-reuse multinomial replication
through one interface. Both retain the original matching graph and individual
weights. The default contribution method applies independent Gaussian multipliers to normalized row contributions and,
where present, unordered-edge and independent-training contributions.
Graphs, predictions, weights,
and nuisance influence corrections remain fixed on these contribution routes.
Successful zero-reciprocal/scalar_psm fitted inference uses complete centered
augmented rows and also accepts supplied full-n fixed-contribution counts.
The regular_joint_d_gt2 fitted scope uses corrected-variance Gaussian replication.
The critical_planar scope uses full nuisance Gaussian draws and an independent
conditional Schur innovation, retaining the original fit and every B slot.
Nested quadratic preparation instead uses its contrast Gaussian law and
covariance-only local-index graphs; the original estimator graph and point are
retained.
}
\details{
Pass a fitted wrapper directly, for example
\code{wm_bootstrap(fitted, method = "fixed_reuse", counts = counts)}.
For its default quadratic correction and \code{refit = NULL}, the adapter
refits every declared PS and used-arm PG model, pooled centers and variances,
and complete weighted quadratic correction with the same full-n multiplicities.
PS fits retain their declared probability or unit fitting weights, PG fits use
multiplicity only, and corrections use multiplicity times supplied weight.
Known offsets, individual weights, original matching graph and incoming loads
remain fixed. Original solver controls and accuracy-only refinement are used;
any failed prescribed refit stops the whole calculation. No graph is rebuilt.
Local-polynomial and quadratic-QR corrections have no automatic count-refit
adapter; supply their justified compatible callback explicitly. Neither this
execution nor a converged refit verifies the sampling/refit-law or simultaneous
numerical-success premises.

For \code{wm_bootstrap(fitted, method = "contribution")}, an initial
\code{inference = "full_x"} request is bound through
\code{wm_model_inference} to the retained complete nuisance root, predictions
and graph. There is no model fit or matching rerun. A point-only wrapper requires
an explicit supported inference handoff before contribution replication;
uncertainty is not inferred from a point alone. An existing public
\code{wm_fitted_inference} component from \code{wm_model_inference} retains its
declared covariance/mixture scope and unchanged downstream option checks.
Its source \code{fit} must be exactly the wrapper's retained original fit;
replacing the inference with another fit's valid object is rejected before
replication, including the advanced mixture branches.
Pass the whole initial fitted wrapper, rather than its historical private
\code{inference} component. Contribution replication and original nonlinear
prediction refits remain different procedures and their finite-data variances
need not coincide. Returned wrapper metadata records the refit mode and source
class; population assumptions remain unverified.

For a nested quadratic preparation from \code{wm_model_inference}, the default
contribution dispatch uses its literal untruncated Gaussian-mixture law. If the
retained contrast influence satisfies
\eqn{(H/\sqrt n)\Pi = QR}, each slot draws independent standard normals zeta
and eta in that order, forms \eqn{h=\Pi R^T\zeta}, evaluates m(h) and s2(h),
and retains the raw root \eqn{L=m(h)+\sqrt{s2(h)}\eta}. The full nuisance
influence and its cross-covariances remain in these evaluated inputs. The
covariance-only local-index graph changes with h; the original estimator,
weights, correction and matching graph are retained without model refitting.

This branch accepts \code{"basic"} (the omitted-interval default) or
\code{"none"}; an explicit normal interval, counts, refit callback or
\code{method="fixed_reuse"} is rejected before RNG. No radial or statistical
Gram threshold, zero completion, column deletion, draw rejection, redraw or
survivor summary is used. Every ordered requested slot, including its normal
innovations and failures, is retained. Any failed slot makes all aggregate
variance, standard-error and interval fields unavailable. B=1 returns draws
without variance or an interval. Numerical evaluation failures retain the
requested contrast and cannot be interpreted as zero roots.

With all B slots successful and B>1, the reported Monte Carlo root variance is
\eqn{\sum_b(L_b-\bar L)^2/(B-1)}, including variation of m(h); the estimator-scale
variance divides this by full observed n. Reverse-basic intervals subtract the
upper/lower empirical quantiles of the \emph{uncentered} raw roots, divided by
\eqn{\sqrt n}, from the original point. The exact full finite-data conditional
variance is not evaluated, and its existence or convergence is not asserted:
\code{conditional_root_n_variance} and \code{conditional_variance} remain NA.

The scoped bounded outcome-correct nested weighted-logistic PS and arm-specific
linear PG sampling/input assumptions remain required. Limit-law quantiles use
predeclared deterministic \eqn{B_n\to\infty} and unique continuous limiting
quantiles. Consistency of the all-B empirical variance additionally requires
\eqn{B_n\le n^\kappa} for a fixed finite kappa. These are theorem premises, not
runtime statistical cutoffs. Same-draw numerical RMS and original-point root-n
numerical error premises are separate; \code{numerical_error_rate_verified} and
population-assumption flags remain FALSE. No finite-sample coverage or original
fixed-reuse bootstrap validity is claimed.

Original context and cached full QR bind before RNG. Computational output/work
caps come from \code{quadratic_control} at preparation; exceeding a cap rejects
the whole request before RNG. These caps neither filter draws nor establish
asymptotic success probability. Callers may supply adequate resources for their
predeclared n and B. Explicit seeds preserve prior RNG state including absence;
normal generation follows the same ordered zeta-then-eta sequence for every
chunk size.

For \code{critical_planar}, H*=Sigmahat^(1/2)Z* and an independent standard
Normal innovation generate the root (bhat+Sigmahat^+Chat)'H*+sqrt(Dhat(H*))xi*.
The full nuisance dimension and raw kernel/Schur/cone diagnostics are retained.
No participant refitting, rematching, count multipliers or failed-draw filtering
is performed. Any failed evaluation makes aggregate inference unavailable.
The exact conditional variance is an unevaluated nuisance-Gaussian expectation;
\code{conditional_root_n_variance} and \code{conditional_variance} are NA.
Root-draw sample variances use B-1. The separately labeled Rao--Blackwell Monte
Carlo estimate is mean(Dhat(H*))+betahat'Sigmahat betahat, with
MCSE sd(Dhat(H*))/sqrt(B). These are Monte Carlo estimates, not exact finite-B
conditional variances. B=1 leaves sample variance and MCSE NA.
An explicit seed preserves its prior value or absence. Within each draw,
p nuisance normals precede the scalar innovation; chunk size does not change
this ordering. Complete theorem/input and strict-CDF premises remain unverified.

A \code{critical_planar} object with
\code{critical_control$direct_variance = "raw_only"} is rejected before RNG:
a rate-free raw scalar variance does not authorize replication. The
\code{"with_mixture"} opt-in preserves the mixture outputs and draws;
direct diagnostics remain in \code{source_inference}. Its support-aligned
raw variance is not the exact operative cone variance, whose positive-part
expectation remains unevaluated.


With \code{method = "fixed_reuse"}, supplied counts are passed to the unchanged
\code{\link{wm_bootstrap_refit}} calculation. Matching maps may be scalar,
distinct-arm double scores or full covariates. Donor weights and original
weighted reuse loads remain fixed; an optional callback replaces the mean
predictions. PATE divides by the replicate all-row weight total and PATT by
the replicate treated-weight total. Replicate variance uses divisor B and
normal intervals use that variance around the original point. All prescribed
columns must succeed. The matrix alone does not verify its draw law, and
sampling validity requires the applicable nuisance-refit expansion and
sampling/replication variance agreement. Basic intervals are not implemented
for this method. Fitted-inference objects are not automatically extracted or
converted to prediction-refit replication. No random counts are generated.
In particular, raw estimated propensity-score matching with \code{refit = NULL}
does not add the fitted-score drift term or its covariance. This option exposes
the original fixed-K replicate algebra; the common interface alone does not
establish validity for estimated-score uncertainty. Complete prediction-refit
WDSM and adjusted contribution replication are distinct procedures.

Compatible prediction-only callbacks may use the same-count raw weighted
least-squares embedding under the complete-system and actual-row conditions
in \code{wm_bootstrap_refit}. For deterministic polynomially growing B, that
function's separate prediction-envelope and simultaneous selected-root/error/
success conditions justify consistency of the empirical nonlinear variance.
The fixed-B Monte Carlo law and the nonlinear conditional-second-moment question
remain separate.

The following formulas apply to \code{method = "contribution"}.

For a successful zero-reciprocal/scalar_psm fitted inference object, write \eqn{a_i} for its centered
complete augmented contribution. Gaussian draws are
\eqn{n^{-1/2}\sum_i a_i\xi_i}. Supplied count columns give
\deqn{T^*=n^{-1/2}\sum_i (N_i-1)a_i.}
Both PATE and PATT use the full observed n, the same source estimate and the fixed
inference rows; there is no replicate treated-weight denominator, nuisance
refit or matching rerun. The exact conditional variance under independent
standard Gaussian multipliers or multinomial(n;1/n,...,1/n) counts is
\eqn{n^{-1}\sum_i a_i^2}. It is computed from the contributions and remains
distinct from the empirical variance of any supplied columns. Supplied count
matrices do not establish independent multinomial draws. Basic intervals based
on arbitrary deterministic columns have no automatic sampling justification.
For d>=2, PATE uses the declared distinct_rarity, full_x or common_field scope;
PATT uses patt. All-scalar corrected graphs also use PATE full_x or PATT
patt under the inference object's stated design-fitted, exact
outcome-fitted-root or complete Gaussian alternative, with known positive supplied weights
and that alternative's weight bounds or integrated moments,
correct full-X predictions, complete joint influence and the corresponding
chart/moment/root conditions. Scalar PATT retains its explicit
zero/current-centering declaration. Population and original numerical
scalar-root graph transfer cannot be verified by count arrays.
The complete rows include the full joint sensitivity and all
matching/nuisance cross-covariances, with zero reciprocal subtraction. The
stored source fit, point, used dimensions, parameter axes, base rows and
covariance identities are checked. A row vector given a fitted-inference label
is insufficient. Private WDSM pipeline results are not automatically coerced.
Population assumptions remain those of the inference object; scalar_psm
also retains its numerical-root requirements. Its general propensity path
retains the raw \code{wm_match} source and supplied \code{scalar_model}
descriptor, rechecks the probability/derivative/known-weight binding and
complete influence and contribution rows, and does not fabricate a logistic
object or fit a model. Sampling transfer additionally
requires the applicable empirical row variance and Lindeberg conditions.
This contribution calculation does not validate the original refit bootstrap.
An unavailable or inconsistent inference object is rejected.

PATE \code{regular_joint_d_gt2} uses the same contribution-method dispatch
with direct roots \eqn{T^*=\sqrt{\widehat V}\xi}, where xi is independent
standard Gaussian and V is the complete fitted contrast variance including
the signed reciprocal subtraction and all nuisance terms. Stored source
point, graph, dimensions, pair products, parameter axes, transport/slope and
corrected covariance identities are reconstructed before drawing. This is
marginal variance-calibrated Gaussian replication, not augmented-row/count
replication; supplied counts are rejected. No fixed-map base interval flag
or hidden base floor gates it: final corrected V must be strictly positive.
Qualified bounded known-W equal-d>2 regular baseline-sheet sampling and
feasible plug-in premises remain unverified. This does not establish joint
nuisance draws, rematching or original nuisance-refit validity. The normal
interval uses analytic V/n and is independent of B. Finite-B empirical
quantiles/variances retain Monte Carlo error; B=1 has NA Monte Carlo variance.
Any growing B gives Monte Carlo variance consistency for these exact Gaussian
roots, given variance consistency, but not for nonlinear refits from weak
convergence alone. No actual finite-sample sampling-variance UI claim is made.
For this path chunk_size bounds normal draws generated together.

For the original wm_match Gaussian branch, the following decomposition applies.

For normalized row contributions \eqn{r_i}, unordered-edge contributions
\eqn{e_j}, and optional independent-training contributions \eqn{t_a}, a root-n draw is
\deqn{T^* = (\sum_i r_i \xi_i + \sum_j e_j \zeta_j +
              \sum_a t_a \omega_a)/\sqrt{n},}
with mutually independent standard-normal multipliers. Its conditional variance
is \eqn{(\sum_i r_i^2 + \sum_j e_j^2 + \sum_a t_a^2)/n}.
The value \eqn{n} is always the evaluation sample size. Training coefficients
must already contain the required sample-size scaling: with \eqn{m} independent
training rows, the training adjustment supplies \eqn{n/m} times the centered
scalar training influence. A missing or \code{NULL}
\code{contributions$training} is an empty vector. Components are not centered separately.
The fitted contributions already contain any supplied nuisance influence
correction. The procedure inherits the fitted estimator's inference assumptions;
it does not refit nuisance models or rematch observations.

The normal interval uses the exact conditional variance, not Monte Carlo
quantiles, so increasing B does not change it. Its sampling validity still
depends on the fitted estimator's asymptotic inference contract. Monte Carlo
variance estimates use the usual sample variance with denominator B minus one.
They are unavailable when B equals one.

With \code{seed = NULL}, draws advance the caller's RNG stream. Within each
replicate, multipliers are drawn in row, edge, then training order. An absent
training component preserves the existing draw sequence. Changing \code{chunk_size}
preserves that order, apart from floating-point summation differences.
Supplied count calculations use no RNG; chunk_size bounds their temporary
row-block arithmetic. Explicit B must agree with the count columns, and B=1
still retains the analytic conditional variance with unavailable empirical
draw variance.
}
\value{
For \code{method = "fixed_reuse"}, the unchanged \code{wm_bootstrap_refit}
list: \code{estimate}, \code{draws}, \code{root_n_draws}, \code{variance},
\code{se}, \code{conf.int}, \code{variance_divisor = "B"}, \code{B},
\code{n}, \code{conf.level}, \code{method}, \code{predictions} and
\code{inference_contract}. With \code{interval = "none"}, only
\code{conf.int} is removed from that list.

For \code{method = "contribution"},
A list with \code{root_n_draws}, multiplier or raw mixture draws on the root-n scale;
\code{draws}, equal to the fitted estimate plus those draws divided by
\code{sqrt(n)}; exact \code{conditional_root_n_variance} and
\code{conditional_variance} in the supported Gaussian/count branches, and
NA in \code{critical_planar} and nested quadratic mixtures; Monte Carlo sample variances
\code{monte_carlo_root_n_variance} and \code{monte_carlo_variance}; and
\code{conf.int}, \code{interval}, \code{conf.level}, \code{B}, \code{n},
\code{seed}, \code{estimate}, and \code{method}. The Monte Carlo variances are
\code{NA} for B equal to one. The interval is \code{NULL} for
\code{interval = "none"}.
For the nested quadratic law, results additionally retain \code{contrast_draws},
\code{standard_contrast_normals}, \code{scalar_innovations}, RNG-kind/seed
provenance, every slot's \code{draw_diagnostics}, and allocation controls.
\code{variance_divisor} is \code{"B-1"}; \code{variance}, \code{se} and the
basic interval are available only after every slot succeeds and B>1. Preparation
and original fit remain in \code{source_inference}. A failed slot is preserved
as unavailable with its failure stage; no roots or diagnostics are replaced.
For \code{critical_planar}, the result additionally retains full
\code{nuisance_draws}, every B-slot \code{draw_diagnostics} with raw signed
kernel/Schur values, cone adjustments, numerical bounds and failure statuses,
and \code{rao_blackwell_root_n_variance_estimate} with its
\code{rao_blackwell_root_n_variance_mcse} and estimator-scale counterparts.
These are finite-B Monte Carlo estimates, not exact conditional variances;
the MCSE is NA for B=1. Any failed draw makes aggregate variance/intervals
unavailable without filtering surviving draws. Basic quantiles or no interval
are supported; an explicitly requested Normal interval is rejected.
For supported fitted inference the result also preserves \code{source_inference},
\code{estimand}, \code{covariance_scope} and the explicit
\code{replication_contract}. Supplied-count results include
\code{supplied_counts} and \code{supplied_draw_law_verified = FALSE};
\code{population_assumptions_verified} and \code{original_refit_bootstrap}
are FALSE. The source inference and its original point are unchanged.
}
\seealso{\code{\link{wm_match}}, \code{\link{wm_fitted_inference}},
  \code{\link{wm_bootstrap_refit}}}
