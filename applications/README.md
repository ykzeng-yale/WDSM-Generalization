# Real-data applications

The current paper uses two applications:

- [ECLS-K](ecls/README.md): obtain the documented ICPSR/NCES sources, stream and
  decode selected fields, construct the corrected cohort, then use the
  [common analysis pipeline](APPLICATION_PIPELINE.md) with declared counts.
- [MEPS 2009](meps2009/README.md): official source preparation and phased point,
  balance, count-refit and contribution-inference commands.

Participant data, fitted objects and historical count archives are not
distributed. These are supplied-weight standardized observed-outcome contrasts
in the stated selected cohorts; their documentation retains identification,
support and inference limits. The application driver accepts newly declared
counts, but those counts do not reproduce a withheld historical stream.

[Common modules](common/README.md) supply prepared-input, analysis, reporting and
saved-graph balance functions. Historical NHANES and NSDUH prepared-input/schema
helpers remain in their directories for compatibility; they are not additional
applications of the current paper. Never infer dependent-survey design standard
errors from the supplied probability weights alone.
