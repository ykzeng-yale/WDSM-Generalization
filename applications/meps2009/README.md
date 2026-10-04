# MEPS 2009 source and complete-record preparation

These two Python CLIs reproduce the declared HC-129 preparation before any
model fitting, matching, effect estimation or inference. They target macOS/Linux
with Python 3.8+ and its standard library; preserve LF line endings in the
distributed JSON specification. The current reproduction was run on macOS.
Windows newline conversion has not been validated. Run without `-O`, since source and format assertions
are part of validation. No participant records or generated outputs are
included in this directory. Keep all generated files in a private location.

## Official inputs

Obtain the ASCII data archive `h129dat.zip` and SAS programming statements
`h129su.txt` from the [AHRQ HC-129 download page](https://meps.ahrq.gov/mepsweb/data_stats/download_data_files_detail.jsp?cboPufNumber=HC-129).
The [official SAS statements](https://meps.ahrq.gov/mepsweb/data_stats/download_data/pufs/h129/h129su.txt)
supply field positions, implied decimals, labels and codes; SAS is not required.
The [codebook](https://meps.ahrq.gov/mepsweb/data_stats/download_data_files_codebook.jsp?PUFId=H129)
documents the source variables.

The decoder requires these exact source identities and stops on a mismatch:

| Input | SHA-256 |
| --- | --- |
| `h129dat.zip` | `a1615e9d0b271e0a605dcdfac84091fc61c6e9b6cc577df9d960bc86fe2a1788` |
| `h129su.txt` | `84a9d48d443619cc85c79ed1b8b02faddaa4eb633bf803df669b6814e146e20e` |

## Run the two stages

From the release repository root, replace the input and output paths below
with local paths. Each output directory must not already exist; both commands
refuse to overwrite a previous preparation. Neither command relies on the
current working directory to find input data.

```sh
python3 applications/meps2009/decode_source.py \
  --raw-zip /path/to/private/raw/h129dat.zip \
  --sas-program /path/to/private/raw/h129su.txt \
  --output-dir /path/to/private/meps2009/source

python3 applications/meps2009/prepare_complete_records.py \
  --source-selected /path/to/private/meps2009/source/source_selected.csv.gz \
  --dictionary /path/to/private/meps2009/source/source_selected_dictionary.json \
  --output-dir /path/to/private/meps2009/complete-record19
```

The first command decodes selected fields in source order, retaining original
missing codes, participant identifiers and row numbers. It writes the source
CSV, a dictionary with official formats and per-variable codebook links, two
SAQ-adult pair-domain CSVs and `SOURCE_DOMAIN_RECEIPT.json`. CSV gzip containers
omit filename and timestamp headers for reproducibility.

The second command reads the adjacent `FROZEN_SPEC.json`. Its `source.path`
and `dictionary.path` are descriptive basenames: actual inputs always come
from the two required CLI arguments. The spec hash, dictionary file hash and
SHA-256 of the exact decompressed selected CSV bytes are checked before output
creation. Different gzip headers, compression levels or zlib versions may
encode the same CSV differently; their compressed-byte hashes need not agree.
The accepted gzip hash remains in `source.sha256` as provenance, while
`source.decoded_sha256` is the required content identity:
`2fe7c07f7b764c9cc36b404fb5f6db8b110df4d74111fdb0cc4a12b844a42a00`.
Receipts record both the actual input gzip hash and the decoded-content hash.
The official raw archive and SAS statement hashes remain strict. Keep the
script and unchanged spec together.

## Fixed scientific definitions

The declared pairs are White–Asian and White–Hispanic. The source rules are
`HISPANX == 2` with `RACEX == 1` for White, `HISPANX == 2` with `RACEX == 4`
for Asian, and `HISPANX == 1` for Hispanic. Filters are applied in their frozen
order: pair membership, `AGE42X >= 18`, `SAQWT09F > 0`, `TOTEXP09 >= 0`, then
complete records in the 19 declared fields.

`Z = 1` denotes White; `Y = TOTEXP09` is unchanged, including zero. The supplied
individual weight is `W = SAQWT09F`, used once. It is not multiplied by
`PERWT09F`, reestimated, normalized, trimmed or used to duplicate records.
`PERWT09F`, source flags, `VARSTR` and `VARPSU` remain provenance fields. This
preparation does not implement cluster or stratum inference.

The covariates, in order, are `AGE42X`, `SEX`, `MARRY42X`, `BMINDX53`, `PCS42`,
`MCS42`, `RTHLTH42`, `DIABDX`, `HIBPDX`, `ASTHDX`, `MIDX`, `STRKDX`, `POVCAT09`,
`EDUCYR`, `INSCOV09`, `REGION42`, `ACTLIM53`, `SOCLIM53`, and `COGLIM53`.
Official missing/inapplicable codes cause complete-record exclusion;
`MARRY42X = 6` is also inapplicable. Unknown or out-of-codebook values cause a
validation error. No imputation, clipping, winsorization or outcome-driven
choice of preprocessing rules is performed.

Continuous fields remain unchanged. Categorical fields retain source codes
and labels; the numeric design uses one-hot columns excluding the smallest
observed valid reference level. Zero-count levels and constant fields are
documented without unidentified columns. One schema and cohort per pair are
shared by all later methods. The exported design has no intercept or scaling;
matching distance and model fitting belong to the subsequent analysis.

## Outputs and interpretation

Complete-record preparation writes, for each pair, `*_complete_records.csv.gz`
and `*_numeric_design.csv.gz`, plus:

- `common_design_schema.json` and a copy of `FROZEN_SPEC.json`;
- `cohort_flow.csv`, `row_disposition_ledger.csv.gz`,
  `covariate_exclusion_ledger.csv.gz`, and `per_field_exclusions.csv`;
- `weighted_baseline_summary.csv` with weighted covariate summaries;
- `PREPROCESSING_RECEIPT.json` with file hashes and run metadata.

All of these outputs belong in private storage, including the record ledgers.
The target is the supplied-weight standardized complete-record population for
each declared pair. It is an adjusted observed-outcome disparity, not a causal
intervention on race or an estimate for all national adults. These declared
rules do not claim exact reproduction of unpublished author preprocessing.
This workflow reproduces input preparation only; it does not certify later
models, matching assumptions, bootstrap validity or analysis results.
