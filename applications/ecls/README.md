# Corrected ECLS-K application: raw selected fields to common WM

This is the corrected parent-reported kindergarten-school-year services case,
with spring-fifth-grade C6R4MSCL mathematics IRT scale outcome. It differs from
the literal7259 File016 replay and from the spring2002 school-administrative
exposure in Keller/Tipton and Morgan. Its primary/secondary targets are
recipient/all-selected-record supplied-W standardized **observed** conditional
outcome contrasts. Temporal overlap, unidentified treatment initiation,
selection and national/causal transport remain substantive limitations.

## Obtain the original files locally

Obtain ICPSR28023-v1 DS0001's `28023-0001-Data.tsv` and codebook from
[ICPSR28023-v1](https://www.icpsr.umich.edu/web/ICPSR/studies/28023/versions/V1).
The supported official TSV identity is MD5
`ebb102e230b0ebd64b506d152ad8972a`,21409 rows,18949 unique column names. Its raw
source36 tokens/order were aligned to the NCES K–8 public-use source. This does
not establish identity to the original authors' unavailable `data.1.RData`
import/version. Preserve raw files; use the source codebook's supplied missing
codes rather than treating all negative values as missing.

Separately obtain the official disability correction
[NCES2007031 report](https://nces.ed.gov/pubs2007/2007031.pdf) and
[correction ZIP](https://nces.ed.gov/pubs2007/data/2007031.zip).
The supported ZIP SHA-256 is
`0004e7d4330e45eb2d2ef0b047cafe018992ac0d2d19d62c820949030bf9e636`.
Its `RP1DISAB.dat` member is255167 bytes, CRC32`847b993c`, SHA-256
`e192801ca944703b08214a035aa491892cded6dd1b78615cd5e3eadf7317308d`.
The actual tab-delimited header is `CHILDID,rp1disab`; the report labels the
revised variable RP1DISABL. It has21260 unique child IDs. The loader accepts this
ZIP or its unchanged extracted `RP1DISAB.dat`, reads only that small member and
joins exact character IDs. It preserves149 unmatched raw child IDs as missing;
it never uses row position, numeric ID casting or old P1DISABL fallback.

No registration credentials, downloaded participant files or extracts are part
of this code distribution. A different source/version fails the known identity
checks and needs explicit source review rather than a silent recoding change.

## Extract and decode without a dense full dataset

With Python3 and standard-library dependencies, use your own local paths:

```bash
python3 applications/ecls/decode_selected_fields.py extract \
  --tsv /your/local/28023-0001-Data.tsv \
  --selected-csv /your/local/new-selected39.csv \
  --receipt /your/local/new-extraction-receipt.json
python3 applications/ecls/decode_selected_fields.py decode \
  --selected-csv /your/local/new-selected39.csv \
  --tokens-sha256 SHA_FROM_EXTRACTION_RECEIPT \
  --revised-disability /your/local/2007031.zip \
  --output /your/local/new-decoded-directory
```

The first step streams the full TSV once, hashing original bytes while writing
only39 selected tokens; it retains source order and every row. The second step
creates a decoded39 CSV, per-cell missing-reason CSV, companion CSV and receipt.
It retains documented skip/refusal/don't-know/not-ascertained versus system
missing, invalid categories and observed values. WKSESL's valid negative values,
including -1, remain numeric. The policy CSV's replication notes describe old
source choices; the corrected construction/schema below determine eligibility.

Birth weight is16*pounds+ounces with valid ounce16 retained; any missing
component leaves total missing. Prematurity uses the P1PREMAT >2weeks lead;
quantity/unit payloads and inconsistencies are descriptive. No lead means
not>2weeks early, not exact0days. Unknown lead or unit is not assigned0.
C1_6FC0 values are unchanged: this supplied longitudinal weight already includes
joint sampling/nonresponse/mover adjustments/trimming/raking. No additional
nonresponse multiplier or estimated-weight model is added.

## Construct the corrected selected data

Install `digest` and `jsonlite` before this step. Run the following block from
the repository root; both the ECLS constructor and common prepared-input
adapter must be loaded.

```r
source("applications/ecls/construct_selected_ecls.R")
source("applications/common/prepared_input.R")
decoder_receipt <- jsonlite::fromJSON("/your/local/new-decoded-directory/semantic_decoding_receipt.json",
                                    simplifyVector=FALSE)
constructed <- wm_ecls_read_views(
  "/your/local/new-decoded-directory/codebook_decoded_selected39.csv",
  "/your/local/new-decoded-directory/selected39_missing_reasons.csv",
  "/your/local/new-decoded-directory/scientific_companion_derivations.csv",
  "applications/ecls/selected_data_schema.dcf", decoder_receipt)
a <- wm_application_ecls_input(constructed$selected,
       provenance=list(decoder=decoder_receipt, data_policy="corrected observed ECLS case"))
```

The constructor reproduces the accepted source statistical block. It defines
separate reading/math/SES school means on every baseline-valid child before
final eligibility, grouped by spring-K S2_ID, unweighted and including the index
child. Field-specific contributor counts/singletons are retained. It requires
actual predictors/derived features, valid Z/Y/W>0 and S2_ID. It does not impose
all39 complete cases, impute unknown expectations, require unused S1_ID/
P1HPARNT/old disability/follow-up quantities or filter on sample size/overlap.
Unknown P1EXPECT remains missing, not six all-zero attainment indicators.

All38 raw matching coordinates stay in X. The nuisance span uses an intercept
and36 columns, omitting expectation reference P1ELHS and TWOPAR only under the
exact family identity. It retains globally/arm-specific QR diagnostics and
explicit defects, without rank repairs. Only the accepted7220/725-recipient
same-span case passes `wm_application_ecls_input`; a different selection remains
its own data construction requiring review, not an automatically relabelled
canonical case. The data adapter does not certify identification or inference.

Use the common driver with the correct supplied7220-by-200 count matrix and
matching IDs/provenance, then `ecls/balance_schema.json` for all38 coordinates
plus22 category indicators. In the accepted original execution,69 PATE refit
columns lacked treated support for the declared nuisance span:68 omitted the
sole treated school-change category, plus one omitted all7 treated expectation
reference-category records. All9 PATE intervals stayed unavailable; PATT's200
control refits remained available. The code preserves those failures rather
than reducing B or dropping columns. It does not promise the same count pattern
for a separately supplied stream, or nominal coverage in this real dataset.

The original participant/count archives are not distributed. For a new analysis,
declare a new stream and generate the full-row multinomial counts once after
checking the adapter's retained row order:

```r
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(20261005)
count_start <- .Random.seed
counts <- rmultinom(200, a$n, rep(1 / a$n, a$n))
count_ids <- a$id
count_provenance <- list(seed = 20261005, RNGkind = RNGkind(),
  R_version = R.version.string, start_state = count_start,
  end_state = .Random.seed, protocol = "Multinomial(n;1/n,...,1/n)",
  matrix_sha256 = digest::digest(counts, algo = "sha256"))
```

Bind those rows, counts and provenance in the common driver as described in
[`APPLICATION_PIPELINE.md`](../APPLICATION_PIPELINE.md). This creates a new
count stream; it does not reproduce the historical 69-column failure pattern.
Do not redraw columns that lose model support or compute intervals from only
successful columns.

The scale-score fields C1R4RSCL/C1R4MSCL/C6R4MSCL are not theta scores. The
[NCES2010-052 theta erratum](https://nces.ed.gov/pubs2010/2010052.pdf) does not
replace them; R4 is calibration notation, not observation-round timing.
