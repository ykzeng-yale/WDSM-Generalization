# Paired simulation metrics with failure accounting

`wm_paired_comparison_metrics()` reports the predeclared, paired Monte Carlo
comparison. It does not replace or change the older all-success
`wm_comparison_metrics()` helper. It fits no models and draws no random numbers.

Each plan row specifies the comparison group, working-model scenario, method,
requested replication count R, whether an interval is required, and a `strict`
or `inclusive` coverage-endpoint rule. The required record statuses distinguish
point availability from interval availability. Replicate IDs must be unique
members of1:R. Reference and candidate must prescribe the same IDs; the helper
never joins only their successful subsets. A comparison group has one fixed
population target. For Lenis, population scenario, effect multiplier and response
setting are separate grouping variables; six multipliers are not six Monte Carlo
replications. For original WDSM, the original design/overlap/estimand groups and
CorCor PSM_M1 reference remain applicable.

Each record also requires a `pairing_id` constructed from the actual data inputs
by the audited collector. Candidate and reference must agree for every common
replicate. Equal integer replicate labels alone do not establish data pairing.
Execution status is `ok`, `failed`, or `not_run`; a present unrun placeholder
remains unevaluated, not an observed failure. Missing rows and unrun placeholders
are counted separately. A literal `Not evaluated` error must not be marked failed.
Normalized point statuses are `ok`, `unavailable`, `failed`, or `failed_point`;
interval statuses are `ok`, `unavailable`, `failed`, `failed_interval`, or
`not_requested`. Unknown/nonterminal component statuses are rejected. Detailed
native errors remain separate diagnostics. Overall `failed` requires an actual
failed required component; overall `ok` requires all required components.
Point-only procedures may be `ok` with an unavailable/not-requested interval.
An unrun row has an unavailable point and no observed component failure.

## Point summaries

For all R valid estimates T_r and target tau, report mean(T), signed bias
mean(T)-tau, its absolute value, signed relative bias100*bias/abs(tau),
RMSE=sqrt(mean((T-tau)^2)), and empirical variance v with divisor R-1.
The signed relative bias is undefined at zero target. Bias MCSE is sqrt(v/R).
The RMSE delta MCSE is sd((T-tau)^2)/(2*sqrt(R)*RMSE), requiring positive RMSE.
It is not assigned a false delta-method value at RMSE=0.

Main point summaries are withheld if any prescribed point is missing, failed or
invalid. A separate `point_success_diagnostics` table explicitly labels the
available-point denominator and conditioning; these summaries cannot be
presented as the unconditional simulation performance. Finite points retained
after an interval failure still belong to the unconditional point comparison.

## Paired relative efficiency and differences

Relative efficiency is RE=v_ref/v_method. Both empirical variances must be
finite and positive for the reported efficiency/MCSE. Let

    q_m,r = R/(R-1) * (T_m,r - mean(T_m))^2,
    phi_r = (q_ref,r - RE*q_m,r) / v_method.

The empirical mean of q is v. Differentiating a variance functional gives its
centered squared residual as the first-order influence; differentiating the
variance ratio gives phi. Its Monte Carlo standard error is sd(phi)/sqrt(R).
The R/(R-1) factor matches the finite-sample variance convention. This is a
first-order delta estimate, not an exact finite-R uncertainty formula; its
usual justification needs independent replications, finite fourth moments and
positive limiting method variance. The reference's own RE is1 with MCSE0.
An underflowed zero ratio or overflowed/nonfinite derived uncertainty receives
an explicit computation status; it is never reported as a zero efficiency gain.
For the fixed Lenis population cache, these Monte Carlo statements are
conditional on those source populations.

The same replicate IDs give candidate-minus-reference differences of estimates,
squared errors and operational coverage indicators. Report each mean difference
and sd(per-replicate difference)/sqrt(R). These retain the covariance between
methods. Separate marginal standard errors cannot replace paired ones.

## Intervals, failures and warnings

An available interval must have an `ok` interval status, finite ordered endpoints,
a nonnegative finite reported estimator variance, and a valid point. The interval
procedure itself remains whatever the simulation declared; the helper does not
replace it or infer its theoretical validity. Native Lenis rows retain their
source strict coverage rule and1.96 endpoints; adapted quick/WM rows retain their
declared qnorm intervals and endpoint convention. Supplied `source_coverage`
values are checked against that convention.

Once all R records exist, availability is the fraction of valid intervals and
operational coverage is the fraction of all R replications that return a valid
covering interval. An observed interval failure contributes0 to operational
coverage. An unrun or missing replication is not a failure: until all records
exist the unconditional coverage/availability summaries are withheld. Conditional
coverage among valid intervals has its own denominator. Coverage MCSE uses the
Bernoulli plug-in standard error; the exact binomial interval is also reported
for operational coverage, retaining uncertainty when its empirical rate is0 or1.

Mean interval lengths and reported variances are explicitly summaries among
valid intervals. Their variance-calibration ratio is reported only when all
points and intervals are available. Warning-bearing finite native outputs remain
included; warning rates are reported, with original messages kept in case files.
Warnings alone are not silently promoted to failures or erased. Invalid claimed
successes and nonfinite derived arithmetic are flagged in the audit. Inventory
completeness, valid reporting inputs, and scientific interpretation are separate
questions; this function grants no publication or theorem certification.
