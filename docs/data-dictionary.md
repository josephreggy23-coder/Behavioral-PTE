# Data dictionary

The canonical input is `data/raw/behavioral_pte_source.xlsx`. Each worksheet is a flat table with one
header row. Units below are taken from column names where available; “not documented” means the
workbook and repository do not provide a defensible definition.

## Shared identifiers

| Field | Meaning |
|---|---|
| `fish_id` | Cohort-specific animal identifier. Prefixes: `fish_` longitudinal, `cf_` molecular, `ptz_` PTZ. |
| `group` | `sham`, `low_impact`, or `high_impact`. Physical dose/impulse definitions are not documented. |
| `clutch` | Biological/recording block: `clutch_A`, `clutch_B`, or `clutch_C`. |
| `timepoint_h` | Hours relative to injury; `-1` is pre-injury baseline. Longitudinal records also use 0.5, 1, 5, and 24 h. |
| `well` | Plate-well number. The supplied protocol describes individual 24-well rearing and 96-well recording, but the mapping is absent. |

## `habituation_trials` — 19,500 rows

One row per fish-session-trial; declared key `(fish_id, timepoint_h, trial)`.

| Field | Type / unit | Meaning |
|---|---|---|
| `trial` | integer, 1-30 | Repeated trial number. The supplied protocol identifies a visual dark flash: 1 s dark with a 15 s interval. |
| `block` | integer, 1-6 | Five-trial block. |
| `distance_mm` | numeric, mm | Distance attributed to the response on that trial. Tracking window is not documented. |
| `responded` | binary | Supplied response classification; threshold is not documented. |
| `baseline_locomotion` | numeric, unit not documented | Session-level locomotion value repeated across the 30 trial rows. |

## `fish_features` — 650 rows

One row per longitudinal fish-session; declared key `(fish_id, timepoint_h)`.

| Field | Type / unit | Meaning |
|---|---|---|
| `amplitude` | numeric, apparently mm | Supplied `A` from `A * exp(-(k-1)/tau) + C`. |
| `decay_constant` | numeric, trials | Supplied habituation decay constant `tau`. |
| `offset` | numeric, apparently mm | Supplied asymptotic offset `C`. |
| `baseline_locomotion` | numeric, unit not documented | Session-level locomotion summary. |
| `resp_prob_block1`, `resp_prob_block6` | proportion | Supplied response proportions for the first and sixth blocks. |
| `mean_distance` | numeric, mm | Mean trial distance for the session. |

## `outcomes_6dpf` — 133 rows

One row per longitudinal fish; declared key `fish_id`.

| Field | Type / unit | Meaning |
|---|---|---|
| `converted` | binary | Supplied 6-dpf seizure-like behavioral label. The protocol proposes a sham-95th-percentile rule, but that rule disagrees with 9–11 stored labels depending on percentile convention. |
| `burst_events_per_hour_6dpf` | events/hour | Supplied burst rate. Values overlap between outcome classes, so this is not itself the label rule. |
| `mean_burst_duration_s` | seconds | Mean supplied burst duration. |

## `cfos_cohort_features` — 249 rows

Same feature fields as `fish_features`, for 86 separate `cf_*` fish at up to three timepoints. The
current code loads this required sheet but does not use it to reproduce risk scores or qPCR pool
assignment.

## `cfos_pools` — 18 rows

One row per qPCR pool; declared key `pool_id`. The table forms nine nominal high/low pairs by
`(group, clutch)`.

| Field | Type / unit | Meaning |
|---|---|---|
| `risk_pool` | category | Supplied `high_risk` or `low_risk` label; scoring/cutoff procedure is missing. |
| `n_larvae_in_pool` | count | Four for every listed pool. |
| `target` | text | `fosab`. |
| `reference_gene` | text | `rpl13a`. |
| `delta_delta_ct` | numeric | Supplied delta-delta Ct summary; raw Ct replicates are absent. |
| `cfos_fold_change` | ratio | Supplied fold change, consistent with `2^(-delta_delta_ct)`. |
| `pooled_fish_ids` | semicolon-delimited text | Four source IDs. This should eventually be normalized to one pool-fish row per record. |

The supplied protocol describes silica-column RNA extraction, DNase treatment, 500-ng reverse
transcription input, SYBR Green chemistry, three technical replicates, and 2^(−ΔΔCt)
quantification. Raw Ct replicates and QC files are not included.

## `ptz_challenge` — 34 rows

One row per separate PTZ-cohort fish; declared key `fish_id`.

| Field | Type / unit | Meaning |
|---|---|---|
| `ptz_mM` | millimolar | PTZ concentration; 2.5 for all records. |
| `latency_s` | seconds | Supplied latency during a reported 30-minute exposure; 1,800 s is the specified right-censoring time. |
| `seized` | binary | Supplied seizure-like response label; scoring rule is missing. |
| `baseline_locomotion` | numeric, unit not documented | Baseline locomotion summary. |

## `session_log` — 15 rows

One row per clutch-timepoint session; declared key `session_id`.

| Field | Type / unit | Meaning |
|---|---|---|
| `n_fish_recorded` | count | Fish recorded in the session. |
| `setup_min`, `acquisition_min`, `analysis_min` | minutes | Supplied session-duration components. |
| `operator_hands_on_min` | minutes | Supplied operator time. |
| `consumables_cost_usd` | USD | Supplied session-level consumables cost. Receipts/cost basis are absent. |
| `fish_lost_this_session` | count | Supplied loss count. “Lost” and disposition are not defined. |
