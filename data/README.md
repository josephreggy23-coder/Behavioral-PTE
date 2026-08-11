# Data

`raw/behavioral_pte_source.xlsx` is the canonical supplied workbook. Treat it as immutable input:
do not clean, relabel, round, or restyle it in place. The pipeline writes all derived artifacts to
`../results/`.

`manifest.json` records the input checksum, byte size, required sheet names, row counts, columns,
and declared keys. A changed checksum means the source workbook changed and the results must be
regenerated and reviewed.

## Data status

The tables are internally well formed: seven flat sheets, no formulas, no missing cells in the used
ranges, no duplicate rows, and no duplicate declared keys. This is a structural validation only.
The repository does not contain enough primary documentation to authenticate how the observations
were acquired or labeled.

An author-supplied acquisition protocol is now documented in
[`docs/experimental-methods.md`](../docs/experimental-methods.md). Before using this dataset in a
competition or publication, link the primary records that substantiate it, including:

- the data owner and permitted reuse terms;
- dates and research site;
- organism source, collection dates, cohort disposition, and welfare/oversight records;
- exact injury calibration, tracking thresholds, outcome-label generation, qPCR QC, and PTZ records;
- raw video or instrument-export inventory and immutable checksums;
- randomization, blinding, exclusions, and outcome-labeling rules;
- applicable prior approvals and the public-safe approval identifiers.

Do not commit signed competition or animal-approval forms, personal contact information,
animal-facility credentials, or other private records. The workbook's embedded Office metadata
should also be reviewed with the data owner before redistribution.

## Workbook notes

- `cfos_cohort_features` contains 86 unique fish across 249 sessions, but the current analysis does
  not use this sheet to reconstruct the supplied qPCR risk bins.
- `cfos_pools` contains 72 unique fish in 18 pools. Fourteen `cf_*` fish are not pooled, and the
  selection rule is not documented.
- `converted` is a supplied label. Its operational definition is not recoverable from the workbook.
- The `fish_*`, `cf_*`, and `ptz_*` prefixes do not overlap; together they represent 253 unique IDs.

See [the data dictionary](../docs/data-dictionary.md) for field-level definitions and unresolved
units.
