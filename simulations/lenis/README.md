This directory supplies the code-only Lenis source-specification comparison. It retains the statistical functions and numerical design, with reviewed path wrappers and separate direct-workflow validation described below. No populations, samples, counts, fits, private execution manifests, or external supplement files are distributed.

Cite Lenis, Nguyen, Dong and Stuart (2019), “It’s all about balance: propensity score matching in the context of complex survey data,” *Biostatistics* 20(1), 147–163, [doi:10.1093/biostatistics/kxx063](https://academic.oup.com/biostatistics/article/20/1/147/4780267). Obtain the supplement from the publisher under its applicable terms. Supply `biosts-17039-File013.R` and `biosts-17039-File015.R` in your own source directory; their required SHA256 values are in [external_sources.json](external_sources.json). Redistribution permission for those external files has not been established, so this export does not include or relicense them. Project-authored code is covered by the package's [GPL-3 license](../../COPYING).

Install this checkout of `WeightedMatching` and make these R packages available: `digest`, `jsonlite`, `sampling`, `survey`, `MatchIt`, `data.table`, `sandwich`, and `quickmatch`. Use an existing R library or configure your own `R_LIBS_USER`. The entry points do not install packages. They use the installed `asNamespace("WeightedMatching")`, including the existing WDSM prediction-stack helpers. Historical source-transfer/result receipts retain the original package identity. New executions record the actual `WeightedMatching` version and namespace path; loading packages alone does not establish numerical reproduction.

Run from any working directory. Replace the illustrative user paths below with your own. Output parents must already exist. Outputs must be new directories outside this package, supplied-source directories, population directories and input runs. `preflight` loads dependencies, parses pinned sources and checks paths/input hashes, without creating output directories or generating populations, samples, counts or fits. Every computational action has a separate explicit verb.

```sh
Rscript simulations/lenis/generate.R preflight /user/lenis-supplement /user/results/populations
Rscript simulations/lenis/generate.R generate /user/lenis-supplement /user/results/populations
Rscript simulations/lenis/run.R preflight /user/lenis-supplement /user/results/populations /user/job.json /user/results/run-001
Rscript simulations/lenis/run.R run /user/lenis-supplement /user/results/populations /user/job.json /user/results/run-001
```

These commands assume the package root is the working directory solely to shorten the script paths; absolute script paths work as well. Copy [job.example.json](job.example.json) to your own input location and choose explicit scenario and replicate IDs. It illustrates two replicates for the three original scenarios, with the formal numerical settings. The complete prescribed comparison uses IDs 1:1000, all scenarios 1:3, six multipliers, response mechanisms No/MAR/MARX/MART, M = 1/3/5 and B = 200. A small subset does not establish Monte Carlo performance. Each invocation evaluates its supplied IDs once; avoid overlapping IDs across inputs, and preserve interrupted outputs rather than overwriting them. These scripts do not schedule, resume, retry or reconcile a frozen study.

`generate` evaluates the literal File015 expressions through its scenario loop, under the existing Mersenne-Twister/Inversion/Rounding RNG convention and source population seed 1357. It retains all 18 populations in the existing three-cache representation and checks the exact common-column and control-outcome identities. `run` evaluates the literal File013 multistage sampling assignments through the unchanged `lenis_source_sample(..., source=)` helper. The historical relative default inside that byte-identical helper is retained only for source provenance; the public run wrapper always supplies the pinned user path, and direct helper callers must do the same. The full original external scripts are never sourced.

The target is finite-source-population PATT for each population, not sample ATT. Native comparators retain source survey regressions, weight transfer, no-replacement nearest matching and strict estimate ± 1.96 SE coverage. Adapted MatchIt quick matching retains its ATT contrast and subclass HC1 sandwich. The weighted-PS repair uses the checked normalized-weight, zero-start quasibinomial logit calculation; supplied final response weights remain fixed. Its separately adopted composed validation is described below. WM retains fixed-M replacement matching, supplied Euclidean coordinates, source corrections, shared multiplier graphs, divisor-B fixed-reuse prediction-refit intervals and qnorm(.975) endpoints. All prescribed count columns and failures are retained. Fitted response weights are held fixed in WM intervals: this is an empirical supplied-fixed external comparison, not a theorem for estimated response weights or dependent cluster/stratum inference. Original printed Table 1 population identity remains unresolved; this code makes no historical numerical-identity claim and does not search for seeds.

For reporting, create your own JSON input with `{"runs":["/user/results/run-001"],"requested":1000}`. Runs must come from this unchanged candidate; original private frozen runs remain the responsibility of their existing collector. A requested value below 1000 can describe a complete initial subset starting at replicate 1, but does not turn that subset into the complete formal study.

```sh
Rscript simulations/lenis/report.R preflight /user/report-input.json /user/results/report
Rscript simulations/lenis/report.R report /user/report-input.json /user/results/report
```

Reporting reuses the existing saved-record schema, provenance adapter and paired metrics. It preserves strict/native and inclusive/WM endpoints, signed relative bias, reference variance ratios, paired MCSE, failures and missing denominators. Partial reporting retains diagnostics and returns failure. Saved-report checks used a private temporary projection of actual original replicate IDs 1 and 2 (3,600 unchanged records) into this wrapper's metadata layout: the complete initial subset output matched the unchanged private reporting functions exactly. A request for 1,000 replicates retained incomplete diagnostics; duplicate paths/record identities were rejected; actual replicate 367 retained all 24 failed point/interval records. The temporary metadata is a wrapper-layout fixture, not original execution or generation evidence. The installed package's 127 function bodies and formals matched current source without calling estimators; RNG state/kind and original inputs remained unchanged.

The public generator was executed for all18 source populations. Its three cache
files, population truths, final RNG record, evaluated source and
multiplier-identity file are byte-identical to the accepted source-reconstructed
outputs. The first direct run evaluated replicate1 in all three scenarios and
four response mechanisms, retaining all six multipliers, M=1/3/5 and B=200.
Independent saved-output review found no scientific record differences in the
1728 unaffected records at the fixed scaled1e-10 rule. The earlier full
comparison did not close as PASS:72 weighted-PS quick records differed from the
accepted corrected reference because its PS solver was different, and those
mismatches remain in the validation history. The exact engine repair now has
adopted composed validation of12 PS fits,12 partitions and72 outcomes using
the already fitted supplied final weights:the first six actual outcomes were
reused without refitting, and the remaining66 were evaluated once. Complete
scientific PS/partition/outcome comparisons agree at the fixed scaled1e-10
rule; maximum absolute and scaled component gaps are5.32907051820075e-15 and
2.41912196384032e-15. The initial validator failed on compact subclass names
and absent startup RNG-seed initialization. The continuation still returned
failure for two all-NA CSV storage-type guards; both original R failure receipts
remain retained. A separate saved-only semantic reconciliation passed all504
numeric record cells:360 finite cells have gap0, and144 literal NA cells match
exactly. This is qualified composed value/status/schema validation, not a claim
that the original R receipts passed or that all RDS storage classes/bytes agree.
The first-case RNG trajectory is not retrospectively certified; the remaining
11 cases preserve RNG exactly after the disclosed accepted-source checkpoint.
There is no postpatch full1800 rerun claim. Precise elapsed-time fields are
separated from scientific comparisons. Count provenance uses saved matrix/column hashes and RNG records;
the old full matrix is not independently scanned. Corrected PS/partition/outcome
fields are mapped explicitly, without claiming identical serialized schemas.
One-replicate integration does not establish Monte Carlo performance,
cross-platform reproduction, printed Table1 identity, estimated-weight theory
or dependent-design validity. [SOURCE_TRANSFER.json](SOURCE_TRANSFER.json)
retains original source origins and the separate repair/validation scopes.
