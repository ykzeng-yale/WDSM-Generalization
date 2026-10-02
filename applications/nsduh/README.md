# NSDUH prepared-input contract

The caller obtains the NSDUH2023 public source release, its codebook and the
original study's age/exposure/outcome/weight construction, preserving QUESTID2
and source order. The supported prepared study has17471 records,1877 treated,
and the shared AMIPY/SMIPY complete-case sample. This adapter does not interpret
arbitrary agency raw sentinels, infer the release/year or recreate upstream
cohort eligibility. Record original source URLs/releases/checksums and the
prepared sample/covariate/weight provenance explicitly. Public data access starts
at [SAMHSA NSDUH data releases](https://www.samhsa.gov/data/data-we-collect/nsduh-national-survey-drug-use-and-health).
This pointer is not a claim that every release has the same source variables.

The closed source reconstruction used `NSDUH_2023.Rdata`, object `data`,56705
rows/2636 columns, MD5`a2736da9213ec34c0717475dcb96c5a0`, local file SHA-256
`41b298c38a34524e9d5ce2f9d7d05627481dda9ce767c23b24cb9eb7ff174a5d`.
These are local conversion identities, not a claim that an arbitrary SAMHSA
download has identical R serialization. Its source cohort rule is AGE3 in4:7;
early_mj_use=1 for observed IRMJAGE<15,0 for observed15<=age<900 or code991,
otherwise missing. Both outcome/known-weight availability are required. The
caller must obtain/build that prepared cohort from its authentic source;
unknown use-age codes are not controls. This export's prepared adapter does
not load that entire RData or regenerate the cohort.

`wm_application_nsduh_input` requires QUESTID2,A,ANALWT2_C,AMIPY,SMIPY and these14
prepared binary coordinates in order: male,white,black,hispanic,income_low,
income_medium,income_high,large_metro,small_metro,early_alcohol,ever_cigarette,
age_21_23,age_24_25,age_26_29. A is the source early_mj_use indicator. W is supplied
known ANALWT2_C, strictly positive, never re-estimated. Both outcomes are0/1.

The original prepared-input recipe defines male(IRSEX1), white/black/Hispanic
(NEWRACE2=1/2/7), income_low(IRFAMIN3<=2), medium(2<code<=4), high(4<code<=6),
metro(COUTYP4=1/2), cigarette(CIGEVER1), early_alcohol(IRALCAGE<15 and<900),
and ages AGE3=5/6/7. It assigns missing early_alcohol0 after that comparison.
That is an inherited source prepared-data choice, not a general scientifically
valid missing-data rule. Verify that upstream source codes/availability are
compatible before applying or supplying these indicators; this adapter never
silently converts raw sentinel codes to0 or assumes unavailable exposure data.

```r
a_ami <- wm_application_nsduh_input(prepared_data, "AMIPY", provenance=your_source_manifest)
a_smi <- wm_application_nsduh_input(prepared_data, "SMIPY", provenance=your_source_manifest)
stopifnot(identical(a_ami$dataset_sha256, a_smi$dataset_sha256),
          identical(a_ami$id, a_smi$id))
```

The common driver keeps source **unit** PS weights, intercept-plus14 nuisance
design, raw14-coordinate full-X metric and known W in WM effects/corrections.
Both outcomes use one authenticated17471-by-200 count matrix, bound to exact
IDs. The matrix must be supplied from its accepted source/count reconstruction;
this slice provides no Linux RNG replay or count generator. Outcome-dependent
nuisance predictions are separate. Explicit callbacks may reuse a previously
accepted PS cache while retaining its count-column binding; there is no implicit
private cache lookup.

`balance_schema.json` declares all14 binary matching coordinates. It is a
diagnostic declaration. The adopted saved graphs have a complete36-graph,
504-row raw-coordinate balance report with independently checked weighted
moments, target-specific fixed scales, incoming masses and ESS. All worsening
and undefined diagnostics are retained; this does not certify overlap or
inference. National/causal transport,
temporality/confounding, mixed support, residual centering and complex-design
inference require their own justification. Historical source point/count/variance
comparisons are optional predecessor reproduction and are not supplied by the
common18-row driver.
