# File and naming conventions

| Content | Canonical format | Location | Rule |
|---|---|---|---|
| Supplied source data | `.xlsx` | `data/raw/` | Immutable; checksum in `data/manifest.json` |
| Derived tables | UTF-8 `.csv` | `results/tables/` | One header row, one variable per column, one observation per row |
| Statistics ledger | UTF-8 `.csv` | `results/all_statistics.csv` | Every reportable statistic with test, effect size, interval, and notes |
| Machine metadata | `.json` | `data/` or `results/` | Valid JSON; explicit schema/version fields |
| Human documentation | UTF-8 Markdown | repository root and `docs/` | Descriptive kebab-case filenames; relative links |
| Figures | 300-dpi `.png` | `results/figures/` | Generated, numbered, and never used as the only source of a value |
| Python | `.py` | `src/`, `tests/`, root runner | Deterministic seed passed explicitly; no hidden notebook state |

Use lowercase snake_case for data fields and generated file stems. Avoid ambiguous authenticity terms
such as `realdata`, spaces in filenames, or manually edited copies named `final_v2`. Generated files
must be reproducible from the canonical input and code.

The source workbook remains `.xlsx` because it arrived as a seven-sheet file. Analysis outputs are
CSV so they are open, reviewable, and easy to diff. Do not create duplicate “cleaned” workbooks
unless the transformation is scripted and the lineage is recorded.
