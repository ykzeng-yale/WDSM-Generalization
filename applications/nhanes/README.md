# NHANES prepared-input contract

The caller obtains and merges the official NHANES2015–2016 and2017–2018 inputs
by exact SEQN, applies the adopted adult/exposure/outcome/covariate eligibility
and preserves selected source order. This slice accepts the prepared case; it
does not supply or claim completion of the upstream raw XPT download/merge.
Record those input versions, source checksums, eligibility reasons and selected
ID/order checksum in the `provenance` argument. Do not upload participant data.

The closed upstream reconstruction used the local merged
`nhanes_2007_2018.rda`, object `nhanes_2007_2018`,59842 rows/1048 columns,
MD5`affadc09f1c6390c5e546281cae5bd08`, local SHA-256
`3a0dd263afe41c30afc8a394bfab1a7b25910cb1cae5659209f0054f1c89c2bf`.
These are author-supplied local merge identities, not public-download byte
claims. Its relevant sleep observations are cycles9/10. Eligibility is age40–70,
observed sleep/HbA1c, DIQ010 not1 with known response, and required covariates/W;
DIQ010/SMQ020/DMDEDUC2 refused/unknown codes7/9 become missing. DIQ010 borderline
code3 remains. Source SLD0122/14 in cycle10 are valid bottom/top-coded values,
not nonresponse. The upstream clean3790 retains167 borderline-diabetes records.
The primary pregnancy exclusion is separately keyed to the accepted original
pregnancy-status table; this adapter does not infer it from an unverified raw
column. The full author merge and pregnancy preprocessing remain caller-supplied
and must be authenticated; a newly downloaded/merged XPT collection does not
gain exact reproduction status from having the expected row count.

Required data.frame columns are `seqn`, `A`, `Y`, `survey_weight`, and seven raw
numeric matching coordinates in this order: ridageyr,riagendr,ridreth1,dmdeduc2,
indfmpir,bmxbmi,smq020. Y is LBXGH/HbA1c in percentage points; A is the adopted
short-sleep indicator SLD012<7. W is unchanged known WTMEC2YR from the prepared
examination sample. A common four-year scaling WTMEC2YR/2 is proportional and
cancels from this self-normalized estimator; it is not a fitted weight.

The primary case3783/1001 treated removes the seven separately recorded pregnant
adults from the clean3790 case. The canonical source3790/1002 treated remains a
different sample with its own count stream. A row-count check alone does not
establish the actual pregnancy-ID/source/order rule; authenticate the prepared
input against accepted preprocessing. The primary adapter never re-filters a
3790 file or substitutes another pregnancy/missingness policy.

```r
a <- wm_application_nhanes_input(prepared_data, "manuscript_primary",
       provenance=list(source_files=your_pinned_source_manifest,
                       sample_file_sha256=your_prepared_file_sha))
```

This retains source probability-weighted PS fitting, raw seven-coordinate
matching and intercept-plus-seven nuisance design. The separate balance schema
uses17 diagnostics: three numeric variables and every declared gender/race/
education/smoking category indicator. It does not factor-expand the matching
metric or nuisance model.

Primary official documentation is
[SLQ2015–2016](https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/SLQ_I.htm),
[SLQ2017–2018](https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/SLQ_J.htm),
[GHB2015–2016](https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/GHB_I.htm),
[GHB2017–2018](https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/GHB_J.htm),
[DEMO2015–2016](https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/DEMO_I.htm), and
[DEMO2017–2018](https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/DEMO_J.htm).
The caller also needs the BMI/smoking source modules to produce their named
prepared covariates. Data/download permission is separate from code distribution.

Short sleep and HbA1c are observational and contemporaneous. The mixed categorical
metric, residual-centering assumptions and complex survey design do not gain
causal/continuous-score asymptotic or coverage certification from this execution.
The common fixed-graph/count interval is qualified empirical algebra. Historical
WDSM/DSM/PGM/PSM/SWPSM comparisons and subclass inference are a separate optional
source comparator path, not output of the common18-row driver.
