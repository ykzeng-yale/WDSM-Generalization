# Historical pre-specified validation design

**Scope amendment, 2026-09-30:** The current WM framework uses supplied known
probability weights and studies estimated scores, coordinate scales and outcome
corrections. The estimated-weight part of question 5 and Section C below is an
excluded historical extension, not a current theory or validation requirement.
The original plan and archived configurations/results are retained for provenance;
this amendment does not authorize a new run or change a frozen study.

This plan is written before estimator-performance results are inspected. The software checks are separate from these experiments. The statistical target is the positive-weight tilted iid law; the designs do not assert validity for dependent survey samples.

## Primary questions

1. Do the explicit marked fixed-M limits predict empirical sampling variance and fitted contribution variance, across matching dimensions 1, 2, 3 and 5?
2. Does bias correction recover root-n centering where the raw matching discrepancy matters?
3. Do the original and stabilized estimators separate as predicted when only weighted residual centering holds?
4. Do feasible regression fits and Gaussian contribution inference reproduce the oracle limits under their stated assumptions?
5. Do estimated-weight and supported estimated-coordinate influence corrections include the required covariance?

## A. Strong-centering marked design

Under the observed law Q, S is uniform on [0,1]^d, Z is Bernoulli(0.4), and U is Bernoulli(0.5), mutually independent. The original covariates include (S,U). Let W0 be 0.5 or 1.5 and W1 be 1 or 3 as U changes from 0 to 1. Conditional residual variances are V0=(0.25,1) and V1=(0.5,1.5). Independent equiprobable signed errors have these variances.

The potential means are m0(S)=1+S1+sum(Sj^2)/(2d) and m1(S)=m0(S)+1+S1. All observed outcomes are bounded. Both PATE and PATT equal 1.5 exactly. Weights and residual variance remain dependent within a score level, so the marked variance formula cannot replace E[V times a weight function] by separate means. The score support and density satisfy the fixed-map conditions in every tested dimension.

Compare original and stabilized rules, both estimands, raw points, oracle correction and correctly specified degree-two polynomial correction. Use fixed M in {1,3}; add M=5 as a sensitivity run if the benchmark precision and measured compute cost allow it. Oracle rho is (1,2). Stabilized feasible rho uses the correctly specified constant model (or the same complete polynomial basis, whose extra coefficients have zero population values), with declared bounds strictly containing the true conditional means.

## B. Weighted-centering design with known failure limits

Start with the target law P in which Z and U are independent fair Bernoulli variables and S is independent uniform on [0,1]^d. Let W0=1+U and W1=1+2U, and define Q by density proportional to 1/W. Then Q(Z=1)=8/17, Q(U=1|Z=0)=1/3, Q(U=1|Z=1)=1/4, and E_Q W=24/17. Potential outcomes are Y(0)=Y(1)=U. Hence PATE=PATT=0 and both weighted score means are 1/2. Strong residual centering conditional on W fails deliberately.

Let v_z(M) be the exact expectation of the weighted mean of U in M iid arm-z donors. The original PATE converges to {v1(M)-v0(M)}/2 and original PATT to 1/2-v0(M). For M=1 these are -1/24 and 1/6. Report errors both from the true causal target and from these derived probability limits. Do not call original-rule noncoverage an implementation failure here. The stabilized rules target zero; the common-score PATE has nonzero reciprocal covariance. Compare its correct row-plus-edge variance with the intentionally incomplete diagonal-only expression, clearly labeled as a diagnostic comparator.

## Explicit variance benchmarks

Compute all discrete donor-mark expectations by finite binomial enumeration, including shared and private donor counts in C_ell. For d=1 use exact beta_ell=2(ell+1) for ell<M-1 and beta_(M-1)=3M/2. For d>1 use the reviewed angular integral with its own independent seed, draw count and Monte Carlo standard error. Propagate the covariance of the beta estimates into the reported benchmark error; never treat a noisy geometric estimate as an exact target. The ordinary constant-weight and M=1 reductions provide additional checks. Estimated contribution variances do not depend on numerical beta evaluation.

## C. Feasible and first-stage extensions

After primary checks pass, add supplied-support spline fits satisfying the recorded strict dimension/order/moment rate bounds, with all requested basis terms retained. For stabilized splines, use fixed data-independent folds and match only within evaluation folds. Record rank failures and clipping counts; do not drop failing replications. Finite-stratum scenarios fit and match within each cell and retain between-cell composition variance.

Estimated-weight and fitted-coordinate experiments require a separate configuration stating the exact nuisance model, influence function, score dimension and applicable theorem. Include the full corrected-row covariance; compare its removal only as a labeled diagnostic. The original plan makes no same-sample scalar estimated-coordinate claim; a later separately qualified Gaussian scalar branch is specified in `first_stage_design.md`. Independent-training and split-potential-mean variants use their separate sample-size and variance scaling.

## Execution and analysis

A pilot uses n=200 and 40 replications per selected scenario to check correctness, runtime, failure capture and approximate effect size. It is not a coverage certification. Primary production sizes are n in {200,800,2000}, with 1,000 independent replications per selected scenario. At 95% coverage, the binomial Monte Carlo SE is about 0.0069 for 1,000 replications. A larger 4,000-replication run is allowed only for conclusions requiring a smaller error and after a resource estimate based on the pilot. Record any design revision before executing the revised configuration, with its reason.

Use the same generated data across competing methods within each replication and independent reproducible streams across replication IDs. Save one result record per requested method and replication, including failures. Seeds must be determined by the frozen configuration and replication ID, not by completion order, worker count, or the success of another method. No replicate is replaced because of an unfavorable estimate, noncoverage, donor deficiency, singularity or a numerical exception.

Report requested/completed/failed counts, bias with Monte Carlo SE, empirical SD, root mean squared error, mean estimated variance, the ratio of mean estimated variance to empirical variance, coverage and its binomial Monte Carlo SE, and mean interval width. Sampling-variance estimates and squared-error estimates need their own empirical Monte Carlo SEs. Use exact Gaussian conditional quantiles for primary intervals; drawing more Gaussian bootstrap replicates does not change these intervals. A separately seeded multiplier check verifies the conditional law.

Only code and synthetic, aggregated results may enter the public repository. Every run records a code/configuration hash, runtime versions, seed definition, resource allocation, elapsed time and all task-owned job IDs. Preserve complete failure logs locally. Full-scale jobs follow independent implementation/statistical review and a successful bounded pilot.
