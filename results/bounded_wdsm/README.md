# Bounded fitted WDSM calibration

This completed synthetic study evaluates the source-form WDSM fitting pipeline
in a bounded iid selection model. The original PS, separate arm PG, pooled
standardization, weighted quadratic correction and finite-donor ratios are
retained. The covariate law and sampling experiment are changed explicitly
from the historical clustered/fixed-population survey design to satisfy the
reviewed bounded full-X branch. It is not a replication of that dependent
survey experiment or a universal weighted-matching validation.

The fixed design is `source_good`, radius 0.25, known W=1/pi_Z(X),
CorCor/MisCor, PATE/PATT and M=1/3. There are 1,000 independent datasets at each
of n=500 and n=2000. The PATE and PATT target effects are 1.3 and approximately
1.300182961805330, respectively, at both sample sizes; the target integrator preserves
the PATT integration bound separately from floating-point error. Models share
datasets. Original fixed-reuse refits use only predetermined IDs 1,5,...,997
at each n, with 399 shared multinomial count columns per dataset.

All 2,000 datasets, 40,000 original requests, 60,000 expanded method records,
4,000 count sets and 1,596,000 prescribed draws were accounted for. There were
no recorded point, interval, dataset or draw failures. An independent automated
audit passed 218 direct statistical calculations plus ownership/accounting
checks. These are separate AI-agent assessments, not human peer review.

## Findings and limitations

Actual fitted analytic variance ratios are 0.9167–0.9822 at n=500, with
coverage 0.933–0.948. At n=2000 they are 0.9821–1.0347, with coverage
0.948–0.961. The ratio compares mean reported root variance with that
estimator's empirical root variance on the same valid datasets. Smaller-n
underestimation is retained. Several larger-n signed biases are two to three
pointwise MCSEs above zero; exact finite-sample unbiasedness is not established.

Compare refit and analytic results on the same 250 datasets, not their
different 250- and 1,000-dataset headline summaries. In MisCor/PATT/M=1/n=500,
the paired refit and analytic variance ratios are 0.8552 and 0.8569. The latter
is 0.9822 on all 1,000 datasets. Refit minus analytic reported root variance
on the shared subset is -0.0244 (outer MCSE 0.1036); coverage difference is
+0.008 (0.008). That subset variation is not evidence of a refit-specific
failure. `paired_comparisons.csv` retains every such paired comparison.

Refit and one-step count variances differ detectably in some finite-sample
cells. Their difference shrinks in several larger-n comparisons but does not
prove a rate or exact equality. Outer MCSE already contains inner-count
randomness; the separately labelled inner component is not added again.
Source intervals use divisor B, while root-variance comparisons use
n*B/(B-1) times that source point variance. PATT also uses total observed n.
One-step exact variance equals the analytic definition algebraically and
is not independent sampling-calibration evidence. No PSM comparator was
allocated, so PSM-relative efficiency is unavailable.

## Reconstructing the summaries

From the repository root:

```sh
Rscript --vanilla results/reproduce_bounded_wdsm.R
```

This verifies the MD5/size manifest, reads the compressed preserved estimator
records, and reconstructs all 72 summary cells, 32 paired comparisons,
16 count-comparison cells and their memberships with the pinned summary
helper. Every saved column is compared, allowing CSV rounding. No subject
data, fits or count draws are generated; RNG state is unchanged. Gzip files
contain the original CSV bytes. SHA256 values are also provided for independent
verification.

The compact export includes count-set moments and draw-completeness counts,
not individual count draws. It cannot independently rederive raw draw moments
or the original archive audit. The historical configuration and scientific
source pins are in `simulations/config/bounded_wdsm_calibration_sharded.json`.
For a new execution on the extended package, use the separately labelled
forward-replay preparation described in
`simulations/bounded_wdsm_calibration.md`; do not relabel a new run as the
historical execution.

The initial production outer monitor has no completion receipt. Its exit
status, disappearance cause and aggregate peak resources remain unverified.
A separately reviewed reconciliation used terminal controller evidence,
all twenty successful child supervisors and verified raw archives. Final
aggregation completed under a new supervisor. No scientific job was rerun to
repair that missing receipt. This limitation is preserved in `provenance.json`.
