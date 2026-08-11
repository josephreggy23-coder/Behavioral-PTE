# Data provenance and authenticity status

## Current status: unresolved

The repository establishes that the workbook is structurally consistent with the analysis and now
includes an author-supplied experimental protocol. It still does not establish who collected the
observations, when or where they were collected, which protocol was approved, or how raw
measurements became the workbook's labels and summaries.

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
| Animals | Strain, husbandry, age and screening protocol | Supplied; source/site and welfare records missing |
| Pressure-wave exposure | Apparatus and sham procedure | Supplied; exact dose geometry, calibration and traces missing |
| Habituation | Visual stimulus, interval, rig and acquisition rate | Supplied; videos, tracking and response threshold missing |
| Outcome | Intended observation window and threshold concept | Supplied; workbook labels do not reproduce from the stated rule |
| qPCR | Tissue handling, extraction, target/reference and QC plan | Supplied; primers, efficiencies, raw Ct and plate/QC records missing |
| PTZ | Concentration, exposure window and censoring concept | Supplied; cohort assignment, vehicle, scoring and SDS records missing |
| Oversight | Research site and applicable approval identifiers | Missing |
| Independence | Student, mentor, and external-support contribution record | Missing |

The protocol and workbook conflict on sample size, timepoints, pressure variables, outcome-label
derivation, qPCR pooling, and PTZ cohort structure. See the
[reconciled experimental methods](experimental-methods.md) for the exact boundary between supplied
methods and data-visible facts.

Signed forms and private approval documents should be shown to the relevant SRC, not published in
Git. A public repository should include only the minimum non-sensitive identifiers needed to make the
record auditable.

## Claim rule while status is unresolved

Use: “In the supplied dataset, the model achieved ...”

Avoid: “I experimentally demonstrated ...,” “the first validated biomarker,” “real-world diagnostic,”
or any statement that implies the repository authenticates animal work it does not document.
