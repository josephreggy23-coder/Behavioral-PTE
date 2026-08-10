# Data provenance and authenticity status

## Current status: unresolved

The repository establishes that the workbook is structurally consistent with the analysis, but it
does not establish who collected the observations, when or where they were collected, which protocol
was approved, or how raw measurements became the workbook's labels and summaries.

The following facts are verifiable from the file itself:

- it is an `.xlsx` workbook with seven flat tables;
- it contains 19,500 trial rows and the cohort records summarized in `data/manifest.json`;
- it was programmatically assembled with an Excel-writing library;
- the nonlinear pipeline nearly reproduces the rounded supplied curve-fit parameters;
- its SHA-256 checksum is
  `8c0b1d903ab058129653ca722683a4d66a99e1937a323150df209c545828599d`.

None of those facts proves that the observations are authentic or fabricated. The correct status is
**unverified** until primary records are linked.

## Evidence required to resolve provenance

| Area | Public-safe evidence to add | Current repository status |
|---|---|---|
| Ownership | Data owner, contributor roles, reuse permission | Missing |
| Chronology | Dated notebook/log and collection window | Missing |
| Source measurements | Inventory and checksums for videos/instrument exports | Missing |
| Blast exposure | Apparatus, calibration, dose units, sham procedure | Missing |
| Habituation | Stimulus waveform/intensity, interval, tracking and response threshold | Missing |
| Outcome | Observation window, burst definition, blinded scoring/adjudication | Missing |
| qPCR | Tissue timing, primers, efficiencies, raw Ct replicates, plate/QC records | Missing |
| PTZ | Exposure timing, scoring rule, censoring, disposition and SDS/risk record | Missing |
| Oversight | Research site and applicable approval identifiers | Missing |
| Independence | Student, mentor, and external-support contribution record | Missing |

Signed forms and private approval documents should be shown to the relevant SRC, not published in
Git. A public repository should include only the minimum non-sensitive identifiers needed to make the
record auditable.

## Claim rule while status is unresolved

Use: “In the supplied dataset, the model achieved ...”

Avoid: “I experimentally demonstrated ...,” “the first validated biomarker,” “real-world diagnostic,”
or any statement that implies the repository authenticates animal work it does not document.
