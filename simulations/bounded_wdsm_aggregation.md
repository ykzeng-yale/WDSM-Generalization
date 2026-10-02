# Read-only bounded-WDSM calibration aggregation

Source `summarize_bounded_wdsm_calibration.R`, call
`wm_cal_collect_records(records_roots)` with the available shard `records/`
directories, and pass its return value to `wm_cal_summarize()`. The collector
reconstructs and checks the entire predeclared scientific grid and common
source/configuration pins. Duplicate ownership is an error. Missing shards,
point records and prescribed count draws remain visible on the full requested
grid. `execution_status` labels unfinished retrieval as provisional. Terminal
dataset errors are distinct from absent completion records. Saved configuration
bytes must match their provenance hash, and recorded source/package pin objects
must agree with that configuration. The fixed scientific fields, full stream
allocation, per-shard ownership and per-dataset allocation/identity are checked.

A dataset marked completed must retain every final requested point and complete
draw table; missing records or draw IDs are integrity errors. A terminal failed
dataset must retain its full point/draw accounting tables. Unproduced entries
must be marked false and retain the exact terminal error. They become explicit
`terminal_unproduced` rows, distinct from missing/active records. A final draw
table can preserve a draw whose individual checkpoint was archived; if both are
present they must agree. No corruption is silently reclassified as a statistical
fit failure.

`execution_complete` requires all global dataset terminals and valid final shard
inventories. The separate `assessment_ready` flag requires calibration mode as
well; it means only that the record accounting is complete and internally
consistent, not that the scientific findings have passed independent review.

The collector accepts runtime-pilot records for schema/identity audits, but
`wm_cal_summarize()` rejects runtime-pilot mode. It must not be bypassed to
turn a runtime pilot into a performance experiment. Sourcing the helper or
running its validator performs no new fits, count draws or datasets.

## Estimands, memberships and failures

`summary` groups branch/model/estimand/M/n. Original fitted, oracle means on
fitted maps and known-map oracle points retain their own empirical sampling
variances. `fitted_naive` changes only the fitted point's variance to V0;
`original_refit` uses that same fitted point on the predeclared count subset.
There is no second one-step coverage method: its exact variance duplicates
the analytic fitted interval.

Every cell keeps three distinct counts: requested datasets, successful points,
and valid intervals. `membership` exposes all dataset IDs and eligibility.
`failures` retains original point and interval reasons. `draw_status` contains
every planned count ID, including absent/failed entries. `paired_membership`
lists the precise analytic/comparator intersection. No successful count subset
is substituted for a complete B-column set.

Point means, bias, RMSE and empirical SD use all successful points and explicitly
state when this is a conditional population. Mean SE, interval length and
conditional coverage use valid intervals. Success-and-cover uses the requested
denominator and assigns zero to unavailable intervals; it is a procedure
reliability measure, distinct from coverage conditional on success. Unprocessed
datasets remain missing and make the result provisional, rather than becoming
unlabelled final failures. Coverage intervals are exact pointwise binomial
intervals; MCSEs use the Bernoulli plug-in formula.

Signed relative bias is `100 * (mean(estimate)-truth)/truth`; it is undefined
at truth zero. Its MCSE divides the bias MCSE by `abs(truth)`. No relative
efficiency is reported because this design allocated no PSM reference.

## Exact scale conventions

For original refit estimates t_b, b=1,...,B, define

    v_B = mean((t_b - mean(t_b))^2).

The original source interval uses `sqrt(v_B)` on the point scale. The separately
labelled root-variance comparison uses `n * B/(B-1) * v_B`. Here n is the total
observed sample size for both PATE and PATT. One-step roots are already in
root-n units; their sample variance uses divisor B-1 with no additional n.
Saved original roots are checked against `sqrt(n)*(t_b-original_point)`.

The original interval and root-variance comparison require all prescribed
refits and the original point to be available. One-step count diagnostics may
remain available when an original refit fails. The B/(B-1) correction is a
divisor correction for the admitted successful-count law; it does not justify
an unconditional failed-draw completion.

## Outer and inner Monte Carlo uncertainty

All formulas below are empirical-influence/delta-method MCSE estimates, not
finite-sample exact confidence statements. For an outer common set of R valid
point/variance records, let

    q_i = (x_i - mean(x))^2,
    v = n * var(x),
    a_i = n * R/(R-1) * (q_i - mean(q)),
    c_i = V_i - mean(V).

The sampling root-variance MCSE is `sd(a)/sqrt(R)`. For the difference
`mean(V)-v` it is `sd(c-a)/sqrt(R)`; for the ratio `mean(V)/v` it is

    sd(c/v - mean(V)*a/v^2) / sqrt(R).

This retains the covariance between the variance estimate and squared centered
point errors. A common-set variance comparison is unavailable with fewer than
two points, and its ratio is unavailable if the empirical variance is zero.
The all-success point sampling summary is retained separately from this common
set. Original-refit versus analytic summaries recompute both sides on the same
joint valid dataset intersection.

Bias MCSE is `sd(x)/sqrt(R)`. RMSE MCSE is
`sd((x-truth)^2)/(2*sqrt(R)*RMSE)` when the RMSE is positive. Empirical-SD MCSE
uses the variance delta method. Mean SE and interval-length MCSEs use their
outer dataset SDs. Paired coverage contrasts use `sd(I_right-I_left)/sqrt(R)`
on the joint valid set, retaining covariance.

Within a dataset, let r_b be original root-n draws, o_b the paired one-step
roots, and Qr_b/Qo_b their squared deviations from their respective draw means.
The root-variance MCSE is `B/(B-1)*sd(Qr)/sqrt(B)`. For the paired variance
difference it is `B/(B-1)*sd(Qr-Qo)/sqrt(B)`. Root-difference mean and RMS use
`sd(r-o)/sqrt(B)` and `sd((r-o)^2)/(2*sqrt(B)*RMS)` respectively. These formulas
use empirical fourth moments and preserve the shared-count covariance.

The outer MCSE of an average refit variance or paired variance difference
already includes inner count randomness. Its separately reported inner
component is `sqrt(sum(inner_mcse_i^2))/R`; it must not be added again to the
outer MCSE. R datasets, not B*R counts, determine the outer sample size.
The per-dataset count table retains both uncertainty levels and the exact
one-step analytic variance identity.

Neither these diagnostics nor their agreement prove uniform integrability,
unconditional finite-sample variance convergence, or the correctness of a
theorem. They compare the actual declared estimator and variance procedures
within their explicitly recorded success populations.

## Focused verification

`validation/check_bounded_wdsm_aggregation.R` uses deterministic arrays to
check units, signed bias, failure denominators, partial/duplicate draw handling,
shared-count covariance, same-set comparisons and the pilot performance guard.
A fixture where reported variance is proportional to the squared centered
point error forces the variance-ratio and difference MCSE to vanish exactly;
an erroneous independent-error MCSE would fail it. The saved runtime pilot is
read solely for schema, target, root conversion and existing count-summary
identity checks. No pilot bias, coverage or calibration summary is generated.
Deterministic corruption fixtures additionally exercise configuration/provenance
mismatches, altered scientific allocations, missing completed outputs, wrong
dataset IDs, allocation mismatches, missing shard completion and false terminal
production flags. They operate on temporary copies of saved schema/CSV records.
