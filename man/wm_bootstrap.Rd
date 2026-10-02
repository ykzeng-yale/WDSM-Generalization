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
Mixed scalar/higher-dimensional graphs and other generic scalar scopes
are unsupported.
For fixed-reuse replication, an ordinary self-normalized \code{wm_match}
result with supplied matching maps of any dimension.}
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
\code{"none"} omits the interval under either method.}
\item{chunk_size}{Positive integer maximum number of multiplier values
generated together for contribution replication. Ignored for fixed-reuse
replication. For regular_joint_d_gt2 this bounds scalar Gaussian draw generation.
No contributions-by-B matrix is allocated.}
\item{counts}{Required for fixed-reuse replication; optional for fitted
contribution replication in zero-reciprocal/scalar_psm scopes only.
\code{regular_joint_d_gt2} rejects supplied counts. Entries must be finite nonnegative integers; every column must
sum to full observed n. Rows correspond to original observations; optional
row names must be literal indices 1,...,n in order. Do not supply seed with
frozen counts. Fixed-reuse replication requires positive multiplicity in both
arms in every column; complete-row contribution replication allows collapsed
columns. The array alone cannot verify its multinomial draw law.}
\item{method}{\code{"contribution"} (default) or \code{"fixed_reuse"}.}
\item{refit}{Optional function of one full-n multiplicity vector returning a
list of mean0 and, for PATE, mean1 predictions. Only supported with
\code{method = "fixed_reuse"}. NULL retains the original predictions.
The callback performs any intended nuisance refitting without changing the
original graph, individual weights or weighted reuse.}
}
\description{
Provides contribution multipliers and fixed-reuse multinomial replication
through one interface. Both retain the original matching graph and individual
weights. The default contribution method applies independent Gaussian multipliers to normalized row contributions and,
where present, unordered-edge and independent-training contributions.
Graphs, predictions, weights,
and nuisance influence corrections remain fixed.
Successful zero-reciprocal/scalar_psm fitted inference uses complete centered
augmented rows and also accepts supplied full-n fixed-contribution counts.
The regular_joint_d_gt2 fitted scope uses corrected-variance Gaussian replication.
}
\details{
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
also retains its numerical-root requirements. Sampling transfer additionally
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
A list with \code{root_n_draws}, centered multiplier draws on the root-n scale;
\code{draws}, equal to the fitted estimate plus those draws divided by
\code{sqrt(n)}; exact \code{conditional_root_n_variance} and
\code{conditional_variance}; Monte Carlo sample variances
\code{monte_carlo_root_n_variance} and \code{monte_carlo_variance}; and
\code{conf.int}, \code{interval}, \code{conf.level}, \code{B}, \code{n},
\code{seed}, \code{estimate}, and \code{method}. The Monte Carlo variances are
\code{NA} for B equal to one. The interval is \code{NULL} for
\code{interval = "none"}.
For supported fitted inference the result also preserves \code{source_inference},
\code{estimand}, \code{covariance_scope} and the explicit
\code{replication_contract}. Supplied-count results include
\code{supplied_counts} and \code{supplied_draw_law_verified = FALSE};
\code{population_assumptions_verified} and \code{original_refit_bootstrap}
are FALSE. The source inference and its original point are unchanged.
}
\seealso{\code{\link{wm_match}}, \code{\link{wm_fitted_inference}},
  \code{\link{wm_bootstrap_refit}}}
