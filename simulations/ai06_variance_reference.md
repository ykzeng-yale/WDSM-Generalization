# Abadie–Imbens variance comparator

Source `ai06_variance_reference.R` for a standalone, unit-weight, common-known-map reference. It uses base R and independent direct neighbor ordering. It does not call the WM matching engine. This helper is comparison infrastructure, not an exported package inference method.

The implemented formulas are Abadie and Imbens (2006), *Econometrica* 74, 235–267, Section 4: equation (14), Theorem 6 and Theorem 7. The user-supplied detailed proof version numbers these equation (19), Theorem 7 and Theorem 8. Its treatment notation is W; this helper uses Z and all target weights equal one.

- M counts opposite-arm imputation donors. J independently counts other same-arm neighbors used to estimate conditional outcome variance.
- The latter estimate is J/(J+1) times the squared difference between the root outcome and its J-neighbor average. For J>1 it is not the sample variance of all J+1 outcomes.
- `conditional_variance` estimates the estimator-scale conditional residual variance. It excludes population effect heterogeneity.
- `marginal_variance` estimates the estimator-scale population-total variance using imputed contrast squares and the source noise correction. It does not add the full residual variance a second time.
- All `root_n` values multiply estimator variance by the total n. ATT's original source scale uses the treated count n1; those values are also returned explicitly as `root_target_count` quantities. The reported SEs are on the point-estimator scale.

The point output is raw matching. A variance formula does not remove its matching bias. Population confidence intervals require negligible root-n bias or separately justified bias correction. The source variance estimator may consistently estimate the limiting variance used by a valid corrected estimator, but that does not make arbitrary correction algorithms valid. This helper deliberately produces no automatic confidence interval.

Use one common supplied score matrix with a fixed Euclidean metric, fixed M and replacement. No automatic standardization occurs. Discrete ties use original-row priority for reproducibility; the continuous-density asymptotic theorem does not thereby extend to unrestricted discrete ties. Each arm must contain at least J+1 rows, and each used opposite-arm donor pool must contain at least M rows.

The helper does not accept probability weights, distinct arm maps, fitted-score adjustments, variable ratios, calipers or no-replacement matching. The target-law and weighted-mark issues in the WM theory require their own estimators. Supplying estimated scores here computes a known-map formula; it does not establish that formula's validity for estimated-score matching.

For native `Matching::Match` comparisons, align coordinates/metric, M, replacement, tie handling, bias correction and `sample`. At the inspected Matching 4.10-15 source, `Var.calc=J` uses the local sample variance including the root. With J=1 and no boundary ties it equals the displayed AI06 squared-difference estimator. J>1 generally differs at finite n. `sample=TRUE` selects the conditional residual variance; `sample=FALSE` selects its population variance calculation. Native variance comparisons require `ties=TRUE` and the standard version: `ties=FALSE` forces the fast version and returns `se=NULL`, even when `Var.calc` is supplied. The deterministic comparison data have no boundary ties, so `ties=TRUE` retains exactly M donors. `se.standard` is a separate statistic and must not be silently substituted for `se`.

This helper and its deterministic external checks do not constitute a Monte Carlo calibration study, a new variance theorem, or proof of agreement with an arbitrary bootstrap. Later comparisons must keep the known/fitted-score regime, target, scale, matching bias and replication law aligned.
