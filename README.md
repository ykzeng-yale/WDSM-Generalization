# Individually weighted matching

Development extension of the `wdsmatch` R package. This repository contains package source, reproducible synthetic simulations and validation code. The manuscript and research source documents are maintained separately and are not included here.

The new functions match each query to exactly `M` opposite-arm donors using supplied Euclidean coordinates. Weights enter the target population, the query contributions and, for the original rule, each donor-set normalization. The original `wdsmatchATE()` and `wdsmatchATT()` interfaces remain available.

For iid observations with law Q and fixed weights W=w_Z(X), the target law is dP/dQ=W/E_Q(W). PATE averages the causal effect under P; PATT averages it among treated individuals under P. Weights may depend on treatment and covariates, with covariate-only weights as a special case. Recovering an external population using sampling or response weights requires the corresponding identification assumptions. Balancing weights also define a target through this law; their name alone does not establish identification or the required residual centering.

## Installation

From this source checkout, run `R CMD INSTALL .` in a terminal. In R, load the package with `library(wdsmatch)`. The package name remains `wdsmatch`; installing this development version updates an existing installation of that name.

## Fixed matching maps

```r
library(wdsmatch)
set.seed(928)
n <- 400
X <- matrix(runif(2 * n), n, 2)
Z <- rbinom(n, 1, plogis(-0.2 + 0.5 * X[, 1]))
W <- ifelse(Z == 1, 1 + X[, 1], 1 + 0.5 * X[, 2])
Y <- 1 + X[, 1] + X[, 2]^2 + Z + runif(n, -0.5, 0.5)

fit <- wm_fit(Y, Z, W, scores0 = X, M = 3,
              estimand = "PATE", regression = "polynomial", degree = 2)
fit$estimate
fit$se
wm_bootstrap(fit, B = 999, seed = 17)$conf.int
```

`wm_match()` accepts supplied mean predictions; `wm_fit()` fits complete weighted polynomials or tensor B-splines. Neither changes the supplied coordinate scaling. Bias correction imputes a query's predicted mean plus matched donor residuals. PATT needs only the control mean. Graphs, donor fractions, incoming loads, contribution rows and numerator identities are returned for inspection.

The default `self_normalized` rule requires the stated strong residual centering conditional on weights and the matching graph's score field. Weighted score centering alone is insufficient. The `stabilized` rule uses an estimated conditional-weight mean and has a different estimator and inference formula. Full-sample stabilized PATE currently requires a common score map and the bounded-outcome theorem; its variance includes both row and unordered-edge contributions. `wm_potential()` and `wm_split_pate()` provide the separately justified one-direction and independent-block constructions.

## Inference and nuisance inputs

- `wm_bootstrap()` applies Gaussian multipliers to fixed contribution rows and, where required, unordered edges and independent-training components. It does not rematch or refit. The default normal interval uses the exact conditional Gaussian variance and does not change with the number of simulation draws.
- `wm_prediction_adjust()` includes qualified parametric prediction/score uncertainty and its covariance with matching contributions. `wm_training_adjust()` uses a separate contribution component for an independent training sample. Their theorem contracts are explicit in the documentation.
- `wm_gaussian_match()` fits a correct full-covariate linear outcome model with conditionally Gaussian errors, evaluates supplied score maps at the actual OLS coefficient, and includes prediction uncertainty. Weights must be fixed functions of treatment and those full covariates, including any outcome-informative weight marks. Its qualified contract permits scalar scores; it does not establish generic non-Gaussian scalar-score inference.
- `wm_logistic_weights()` fits a regular bounded-domain logistic weight model and returns the weight derivatives and nuisance influence vectors. `wm_weight_adjust()` includes their covariance with the matching contributions on a fixed graph.
- Spline fits require a declared fixed support, moment order, spline order and mesh exponent satisfying the documented strict rate window. Stabilized spline fits train outside fixed evaluation folds and match wholly inside each evaluation fold. Singular bases and insufficient donors produce explicit errors.
- `strata` restricts matches and nuisance fits. It does not specify a complex-survey variance estimator. The theory concerns iid observations with positive bounded individual weights, appropriate overlap and continuous matching geometry within each stratum.

The software cannot establish causal identification, smoothness or an influence-function expansion from the observed data. Supplied learned coordinates or weights do not automatically include their estimation uncertainty. See the function documentation for each supported contract and limitation. Generic scalar fitted-score inference, dependent survey designs, growing dimension, variable donor counts and unrestricted discrete ties are outside the new general inference claims.

## Explicit variance constants

```r
wm_geometry(M = 3, d = 1)  # exact one-dimensional formula
wm_geometry_overlap(M = 3, d = 2, draws = 20000, seed = 92) # complete beta vector
wm_geometry_alpha(M = 3, d = 2, draws = 20000, seed = 93)   # alpha only
```

For general dimensions, the boundary-rank integral returns component Monte Carlo errors, their covariance, and precision diagnostics. The original `wm_geometry()` angular integral remains available as a separate implementation, but can have large numerical error when M and d increase. `wm_geometry_alpha()` uses a separate bounded-integrand formula and does not supply the marked-overlap vector needed with variable donor weights. Zero-hit components are unresolved. Numerical estimates are not projected onto the Jensen lower bound. Contribution-based standard errors do not require evaluating these constants numerically.

## Validation and simulations

The matching, regression, nuisance-adjustment, geometry and bootstrap snapshot `c2817610a917` passed `R CMD check --no-manual`, including 2,347 test expectations with no failures, warnings or skips, and seven standalone pipeline checks. Subsequent documentation corrections passed independent agent review and parsing of all 18 R help files; they do not change estimator arithmetic. Collector and plotting changes passed separate regression checks. The optimized matching engine also passed 147 complete comparisons against its frozen reference and four parity cases against the legacy public interfaces. The PDF reference manual was not compiled in that check. Numerical geometry precision and statistical simulation performance are assessed separately from software tests. Reviews in this development workflow were performed by separate AI agents; statistical reaggregation used independently implemented checks.

The pre-specified primary plan is in [simulations/DESIGN.md](simulations/DESIGN.md); first-stage and supplemental plans are in [simulations/first_stage_design.md](simulations/first_stage_design.md) and [simulations/supplement_design.md](simulations/supplement_design.md). They include exact targets, marked-weight variance benchmarks, assumption-failure examples, independent random streams, explicit failure records and Monte Carlo errors. The primary study completed 48,000 datasets and 576,000 comparison records; the supplemental study completed 1,600 datasets and 10,800 records; the first-stage study completed 39,000 datasets and 234,000 records. All three received independent raw-record, summary and benchmark checks. The first-stage study retains 28 genuine fit-failure records from seven datasets and explicitly documents recovery from a Gaussian graph-comparison assertion defect. Its adjusted coverage at n=2,000 ranges from 0.931 to 0.964 across the declared cells; uniform finite-sample calibration is not claimed. Qualified findings and deterministic aggregate reconstruction are in [results/RESULTS.md](results/RESULTS.md) and [results/REPRODUCE.md](results/REPRODUCE.md). Pilot results are not a coverage certification.

## Provenance

The starting package source is the author's GPL-3 `wdsmatch` version 0.2.1, pinned in [SOURCE_PROVENANCE.json](SOURCE_PROVENANCE.json). Existing author credits and interfaces are retained. The supplied example dataset is synthetic and was recovered exactly from the upstream Git revision; its generator is included in `data-raw/`. No real study data or private source cache is included.

The software is distributed under GPL-3. [COPYING](COPYING) contains the unmodified [GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.txt).
