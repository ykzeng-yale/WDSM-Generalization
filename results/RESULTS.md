# Validated aggregate results

This export contains the completed primary, supplemental and first-stage synthetic studies, their exact configurations and numerical geometry inputs. Independent automated audits reproduced the retained records, sampling summaries and variance benchmarks before this aggregate export. These checks do not certify proofs or uniform finite-sample calibration.

| Study | Frozen source snapshot | Independent datasets | Paired estimator records | Aggregate rows |
|---|---|---:|---:|---:|
| Primary | `a8f9cf66188d` | 48,000 | 576,000 | 576 |
| Supplemental | `7fa0a62d1f9f` | 1,600 | 10,800 | 54 |
| First-stage | `7fa0a62d1f9f`, with separately recorded guard correction | 39,000 | 234,000 | 234 |

Every primary scenario contains replication IDs 1–1,000; every supplemental scenario contains IDs 1–200. Both studies have zero estimator failures. The primary execution retained completed prefixes after infrastructure timeouts and generated only the 168 missing dataset IDs in 25 recovery batches. The supplemental execution retained its eight initial five-replication batches. Configuration, seed and recovery tables preserve those decisions.

## Primary findings and limitations

The grid crosses `n = 200, 800, 2000`, matching dimension `d = 1, 2, 3, 5`, `M = 1, 3`, and strong/weak centering designs. The joined table distinguishes 288 corrected-estimator variance comparisons, 192 raw point-only rows and 96 original-rule weak-centering diagnostics. **The original rule's weak-centering variance theorem is not asserted.** Its known wrong probability limits remain reported; a seemingly reasonable variance estimate cannot remove that target bias. For example, original PATT at `M=1,d=1,n=2000` has bias 0.16642, near the exact limit 1/6, and coverage 0.001.

Among the valid corrected comparisons, mean reported variance divided by its asymptotic benchmark ranges from 0.7991–1.0458 at `n=200`, 0.9530–1.0275 at `n=800`, and 0.9876–1.0220 at `n=2000`. Corresponding coverage ranges are 0.918–0.964, 0.920–0.966 and 0.940–0.964. These are descriptive ranges across paired cells, not simultaneous calibration guarantees. At `d=5,M=3,n=200`, strong-design PATE polynomial fits have reported/empirical variance ratios about 0.808 and coverage 0.918–0.919. The original fit exhibits this deficit without normalizer clipping, so clipping alone does not explain it.

Some departures from the asymptotic variance remain many Monte Carlo standard errors from zero even at `n=2000`. These are finite-sample departures, not proof that the limiting variance is wrong or that every difference is simulation noise. The experiments do not isolate their mechanism. Clipping remains recorded: 1,865 valid primary estimator records at `n=200`, three at `n=800`, and none at `n=2000`; paired records are not distinct datasets.

## Supplemental findings and limitations

Eight prespecified scenarios examine nonlinear order-three spline fits in dimensions one and two, finite random strata, and independent potential-mean blocks with different score dimensions. There are 200 datasets per scenario. The strata targets are exactly `393/281` for PATE and `41/27` for PATT; their variances include random stratum composition.

Two-dimensional stabilized spline clipping is a material finite-sample feature. Counts of fits with any clipping are:

| n | PATE / 200 | PATT / 200 |
|---:|---:|---:|
| 800 | 193 | 126 |
| 2000 | 82 | 23 |

For these feasible fits at `n=800`, empirical/limit and reported/limit variance ratios are 1.3024/1.3624 for PATE and 1.1588/1.2440 for PATT. At `n=2000` they are 1.0238/1.0718 and 0.9783/1.0655. Thus reported variances remain about 6.6–7.2% above the limit at the larger sample size. These runs neither justify describing clipping as inactive nor identify how much of the inflation comes from each fitting step.

The half-split contrast covers in 198/200 replicates: 0.990, with exact pointwise 95% interval `[0.96435, 0.99879]`. Its component sample correlation is 0.25114; the empirical contrast variance is 0.80562 times its benchmark, whereas the mean reported variance is 0.99252 times it. Independent data blocks can have nonzero realized sample covariance. No concrete implementation route for population cross-block dependence was found. The recorded covariance must remain in the empirical contrast variance. The one-third split has correlation −0.00494, empirical/limit ratio 1.02552 and coverage 0.940 (`[0.89754, 0.96862]`). These different prespecified outcomes are both retained. Component rows use their own `analysis_n`; the contrast uses total `n`.

## First-stage findings and limitations

The 39 prespecified scenarios have 1,000 datasets each and compare known-first-stage, fitted-naive and fitted-adjusted inference for PATE and PATT on paired data. Same-sample score fitting uses dimensions 2, 3 and 5; independent training and fitted individual weights use dimensions 1, 2, 3 and 5. The separately qualified Gaussian same-sample branch uses scalar scores and M=1 or 3. Other branches use M=3; sample sizes are 200, 800 and 2,000. Score-calibration scenarios use declared known offsets, including the constant effect, to isolate the first-stage correction; they are analytically checkable validation constructions.

At n=2,000, coverage ranges across each branch's estimands and matching settings are:

| Branch | Fitted naive | Fitted adjusted | Adjusted mean reported variance / its limit |
|---|---:|---:|---:|
| Same-sample scores, d=2,3,5 | 0.977–0.985 | 0.951–0.959 | 0.9985–1.0278 |
| Independent score training | 0.900–0.926 | 0.931–0.955 | 1.0001–1.0281 |
| Estimated individual weights | 0.951–0.977 | 0.944–0.964 | 0.9945–1.0091 |
| Qualified Gaussian scalar scores | 0.966–0.974 | 0.934–0.956 | 0.9962–0.9983 |

These descriptive ranges are not simultaneous uncertainty intervals. The paired comparisons retain the sign of the first-stage effect: covariance can reduce variance for same-sample fitting, whereas independent training adds a separate variance component. The naive fitted estimator's sampling variance and its reported-variance limit are different quantities; both appear in the joined table.

There are material finite-sample departures. For independent training with d=5, M=3, n=2,000, adjusted PATE coverage is 0.931, with exact pointwise 95% interval [0.91348, 0.94592]. Empirical variance is 1.10422 times its asymptotic benchmark, while mean reported variance is 0.92123 times empirical variance. Gaussian scalar PATT with M=3 and n=2,000 has coverage 0.934 [0.91679, 0.94859]. These outcomes remain reported; the experiment does not establish their mechanism or justify calling every departure simulation noise. Cells share paired estimators and numerical geometry, and the many pointwise intervals are not a joint test.

Seven estimated-weight datasets at n=200 reached the declared logistic parameter boundary: three in d=2, three in d=3 and one in d=5. Their 28 fitted-comparison failures remain, together with the valid known-first-stage rows. There are 233,972 successful records in total. Sampling summaries condition on successful fits and expose the requested, successful and failed denominators; no failed dataset was replaced.

Recovery preserved 222,813 original rows exactly, replaced only 327 eligible Gaussian assertion-error rows, and completed 1,810 missing dataset IDs (10,860 records) using unchanged configurations and seeds. The Gaussian harness had compared machine-scale distance differences from equivalent OLS calculations too strictly. Its corrected assertion retains exact donor identity and fractions, and bounds distance differences by the actual coordinate perturbation. Estimator formulas were unchanged. Full six-row replay blocks checked successful siblings; only the declared failed rows entered the derived results. The original 355 failure records remain in the private audit history, including all 28 genuine fit failures. A narrowly declared parameter-pairing allowance applies only to 11 cross-platform Gaussian replay keys; it checks machine-scale coefficient differences without altering any record. Separate recovery/export script corrections affected execution and metadata handling only. Aggregate recovery counts and both scientific source maps are in `first_stage/`.

## Numerical and Monte Carlo uncertainty

Geometry snapshot `4485ffbfce20` uses exact one-dimensional constants and independent boundary-rank Monte Carlo in higher dimensions. The saved covariance is the covariance **of the estimated beta vector**, already divided by the integration draw count. Keep its off-diagonal entries when propagating benchmark uncertainty.

The largest propagated relative geometry MCSE is 0.92664% for the primary variance benchmarks and 0.14772% for the supplemental benchmarks. This satisfies the declared 1% *benchmark* goal as an empirical diagnostic; it does not establish 1% accuracy of every beta component. At `M=3,d=5`, the largest beta-component relative MCSE is about 3.375%, and alpha relative MCSE is about 2.562%. Empirical integration errors are not rigorous bounds. Numerical constants are not projected onto a theoretical lower bound.

The first-stage assembly uses a separately seeded 200,000-draw alpha-only integral for M=3,d=5, selected before inspecting estimator performance. Its alpha estimate is 9.25670473 with empirical MCSE 0.01135639. The maximum propagated relative geometry MCSE across first-stage benchmarks is 0.56579%. The saved selection table identifies the actual geometry used; it does not replace the primary or supplemental beta-vector inputs. The scalar alpha object supplies no marked beta vector.

Simulation and geometry streams are independent. Their MCSEs are stored separately and combined only for the indicated descriptive discrepancies. Shared geometry induces correlation across benchmark rows. Empirical-variance MCSEs are large-replication approximations, variance-ratio MCSEs are paired delta-method approximations, and coverage intervals are pointwise. Supplemental column names abbreviate these qualifiers. No across-cell independence or simultaneous coverage claim follows from the tables.

## Files and reproduction

- `primary/` and `supplement/`: unchanged aggregate summaries and joins; exact seeds/configurations and task plans. `primary/recovery_audit.csv` replaces execution paths with artifact basenames.
- `first_stage/`: independently audited summaries and joins, configuration/task plan, aggregate recovery counts, actual geometry selection and sanitized provenance. `reproduce_first_stage.R` reconstructs its 234 benchmark rows and shared geometry covariance.
- `geometry/`: eight geometry cases, 16 beta components, full beta covariance, and compact boundary-term diagnostics/covariance. One-dimensional zero covariance is exact, not an imputation for missing numerical precision.
- `provenance.json`, `provenance/` and `artifact_hashes.csv`: source snapshot fingerprints, transformation notes and public-file hashes. Snapshot IDs identify execution freezes, not Git commits or the current source tree.
- [REPRODUCE.md](REPRODUCE.md): repeat the aggregate joins and figures from a fresh source checkout, or run explicitly selected synthetic replications. The aggregate check uses no random draws.

No individual-subject data, per-replication estimator records, manuscripts, proofs, private execution metadata or internal reviews are included. The aggregate-only export reproduces benchmark joins and figures; reproducing sampling summaries requires regenerating the declared synthetic replications.
