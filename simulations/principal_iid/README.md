# Principal iid weighted-matching reproduction

These scripts reproduce the fitted-score comparison using the original WDSM
six-covariate Gaussian-mixture model, independently sampled from its selected
law. They contain no participant data and require no external data account.
They do not reproduce the separate frozen finite-population sampling experiment.

The common implementation is loaded from the adjacent package `R/` sources.
Point estimates use `wm_match`; inference uses `wm_fitted_inference` and
`wm_bootstrap`. The existing WDSM nuisance fitter and saved full-X stack helper
supply the complete fitted parameter. This directory adds the model, one-sample
recipe and reporting entrypoints, not another matching or variance algorithm.
The source inventory binds the package sources and the two existing simulation
helpers used here. Keep this directory in `simulations/principal_iid/` of the
same package checkout.

## Model, targets and fixed recipe

The 200 equally weighted source means are reconstructed from seed 20220930:
first a 10-by-7 matrix of N(0,0.35^2) offsets, then a 20-by-7 matrix of
N(0,0.15^2) offsets, in R column order. All 210 draws are retained before
selecting the first six columns. Observations have independent unit-variance
normal coordinates conditional on their component. The GoodOverlap and
PoorOverlap treatment models and selection coefficients retain the source
values in `principal_iid_model.R`.

The sampler uses exact exponential-normal tilting and rejection to draw iid
rows from Q, proportional to the source law times its selection probability.
It retains the first n accepted rows. Supplied weights are exactly 1/pi_Z(X),
with no clipping or weight fitting. Independent E0,E1 give
Y0=mu0(X)+E0 and Y1=Y0+delta(X)+E1. The source mixture averages define the targets;
there is no newly simulated finite-frame target.

PATE is analytic. PATT uses the unchanged three-dimensional Gaussian
Gauss–Hermite/Stein calculation at orders 32,64,96,128. Both last successive
changes must be at most 1e-7 before its order-128 value is used. This operational
convergence check is not a rigorous numerical-error enclosure. All order and
component results are saved. Both targets retain the predeclared +/-1e-6
arithmetic sensitivity display; PATE remains analytic, without quadrature
uncertainty.

Each sample uses M=3, row-order ties, pooled unweighted centering and unit RMS
scaling, PATE and PATT, and these three methods:

| Method | Matching coordinates | Correction | Full nuisance p (PATE/PATT) |
| --- | --- | --- | --- |
| PS | Weighted logistic propensity, one coordinate | Supplied-W WLS on intercept, X1–X6, X1:X2 | 26 / 18 |
| DSM | Distinct-arm propensity/prognostic pairs | Supplied-W complete six-term own-score quadratic | 42 / 26 |
| X6 | Six standardized covariates | Supplied-W complete 28-term quadratic in standardized X | 68 / 40 |

The propensity model includes X1:X2. DSM prognostic regressions are unweighted
arm OLS on the same eight-column design. PATT does not fit an unused treated
prognostic/correction model. All fitted scale, prediction and shared score
blocks remain in the complete covariance. These PS estimates are corrected
full-X predictions, not raw scalar propensity matching.

One n-by-200 multinomial count matrix is shared by all six method/target
combinations. Complete-contribution inference uses the common full-p fitted
rows. The original fixed-reuse callback keeps W, donors and incoming loads,
refits the declared predictions under those same counts and recomputes the
counted target denominator. PS/X6 callbacks refit their full-X WLS predictions;
DSM refits its complete prediction stack. The source nonlinear variance has
divisor B=200. Contribution Monte Carlo draws are retained separately; they are
not a second sampling-variance method or a substitute for the original refit.

The original numerical solver controls and accuracy checks are unchanged.
A numerical root check does not certify an exact statistical root or the
population assumptions. No failed slot is dropped, zero-filled or redrawn.
All 200 prediction-fit attempts and their coefficients/errors are retained;
the public refit API may stop at its first callback error, in which case its
later scalar positions are explicitly unattempted and its full interval is
unavailable. Independently retained manual scalars do not become a
survivor-only variance. Sample, graph, point and inference failures remain
separate.

## Commands

Use R 4.4.2 with `digest` and `jsonlite` already installed, and Python 3 with its
standard library. The scripts install nothing. Set numerical threads before R:

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
mkdir -p reproduction/cases
Rscript --vanilla simulations/principal_iid/reproduce.R truth reproduction/truth
Rscript --vanilla simulations/principal_iid/reproduce.R case reproduction/truth GoodOverlap 1000 1 reproduction/cases/g1_GoodOverlap_n1000_rep0001
```

Each output directory must be new. The second command computes one declared
sample, not a coverage experiment. It saves the synthetic sample and counts,
original points and complete stacks/inference, four compact 50-slot coefficient
chunks per method, replication results, six small method JSON files and a
terminal result with source/input hashes. Interrupted outputs remain incomplete;
the scripts have no automatic retry or job scheduler.

The complete comparison contains 1,000 independent replications in each of
GoodOverlap/PoorOverlap crossed with n=1000/5000. The group order is Good1000=1,
Good5000=2,Poor1000=3,Poor5000=4. For group g and replication r, the data seed is
104000000+100000*g+r and the count seed is 204000000+100000*g+r. Invoke the same
case command for the remaining keys, choosing execution resources separately.
Do not rerun a completed key or treat a partial directory as a completed case.
Truth calculation and each case run in one process; the user chooses an
appropriate external memory/time bound. No resource monitor is embedded in
the statistical recipe.

After all 4,000 cases are complete:

```sh
python3 simulations/principal_iid/summarize.py reproduction/truth reproduction/cases reproduction/summary
```

The summary refuses missing/duplicate keys and nonterminal or integrity-failed
cases. It checks all 24,000 method rows, full-p labels, fixed targets, paired
sample/count hashes and the six JSON records per case. It uses small JSON/CSV
outputs; it does not decode or rehash the large RDS arrays. Creation-time RDS
hashes remain in the case receipts.

Signed relative bias is 100*(mean estimate-target)/abs(target), empirical
variance uses R-1, and empirical relative efficiency uses PS within the same
overlap/n/estimand/M. Coverage is inclusive and its operational denominator is
all 1,000 requested samples; failed intervals contribute zero to operational
coverage. If points fail, full-estimator point performance is withheld and
available-only diagnostics are labeled separately. If intervals fail, the
full mean-variance/calibration summary is withheld rather than silently using
survivors. Paired analytic/refit coverage differences and variance-ratio MCSEs
retain the pairing. `formal_metrics.py` contains these reporting definitions.

Outputs include a 24-row method table, 48-row analytic/refit table and figure
data, paired summaries, failure records, and target sensitivity. Their
interpretation remains repeated-sampling evidence under this synthetic law,
not a proof of theorem assumptions or a guarantee of 93–97% coverage in every
finite setting. Platform BLAS and floating-point details may affect numerical
fits and failures; preserve the recorded warnings, R version and source hashes.
