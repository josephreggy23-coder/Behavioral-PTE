# Results — startle-habituation kinetics and a supplied post-injury behavioral label

**Source:** `data/raw/behavioral_pte_source.xlsx`

**Random seed: `20260809`** (numpy, scikit-learn, permutation and bootstrap draws). With the recorded input and tested dependency versions, rerunning `python run_all.py` reproduces the analysis numerically.

> **Interpretation boundary.** These analyses use the supplied workbook. The repository does not include source videos, electrographic recordings, raw qPCR Ct data, animal-approval records, or a reproducible protocol for deriving the `converted` and `risk_pool` labels. The results therefore describe internal associations with supplied labels; they do not establish epilepsy or validate a biomarker.

Every statistic quoted here is also in [`results/all_statistics.csv`](results/all_statistics.csv) with its test, statistic, df, p-value, effect size and CI.

---

## Summary of findings

| # | Finding | Test | Effect size | p |
|---|---------|------|-------------|---|
| 1 | The nonlinear refit reproduces the supplied `decay_constant` | Pearson correlation, 650 sessions | r = 1.0000 | < 0.0001 |
| 2 | No baseline group difference was detected | one-way ANOVA on baseline τ | η² = 0.0083 | = 0.5825 |
| 3 | Low and high dose move τ in **opposite** directions after pressure-wave exposure | Welch t, Δτ at 0.5 h | d = 2.19 | < 0.0001 |
| 4 | The 4-predictor model distinguishes the supplied `converted` label across held-out clutches | nested CV + 1000 permutations | AUC = 0.833 | = 0.0010 |
| 5 | Orthogonal molecular validation: supplied high-risk pools have higher c-fos in the nominal all-pair comparison | paired t on log2 fold change, 9 matched pairs | dz = 0.88 | = 0.0293 |
| 6 | No high-vs-low difference was detected in the three sham pairs | paired t, 3 pairs | dz = 0.46 | = 0.5124 |
| 7 | The PTZ group comparison did not reach p < 0.05 and is underpowered | χ² | V = 0.381 | = 0.0849 |

---

## Step 1 — Re-estimating the habituation feature

Per fish per session, `distance_mm(k) = A·exp(−(k−1)/τ) + C` was fitted to all 30 trials by nonlinear least squares (`scipy.optimize.curve_fit`, bounded, four starting points per session, best SSE retained).

- **650 fish-sessions fitted**, convergence 100.0%, median R² = 0.796 (5th percentile 0.504); 1 session hit a parameter bound.
- **Agreement with the supplied `fish_features.decay_constant`:** Pearson r = 1.0000 (p < 0.0001), Spearman ρ = 1.0000. Mean bias (refit − supplied) = -0.001 trials, 95% CI [-0.002, +0.001]; 95% limits of agreement [-0.031, +0.030]; median absolute error 0.00%.
  The refit is used for everything downstream; the supplied column is treated as a check, not an input.

### Why log-linearisation was not used

Subtracting an estimated offset and regressing `log(y − C)` on trial fails on this data for a mechanical reason: **28.1% of trials fall at or below the habituated floor**, so `y − C` is non-positive and cannot be logged. Discarding those points tilts the fitted slope, and in **4 of 650 sessions the recovered τ comes back negative** — the sign inverts. Correlation with the supplied τ collapses from r = 1.000 (nonlinear) to r = -0.229 (log-linear). See `fig02_tau_agreement.png`, right panel.

Figures: `fig01_curvefit_examples.png`, `fig02_tau_agreement.png`.

---

## Step 2 — Prediction model (primary result)

**Analysis set.** Injured fish only (sham dropped), one row per fish, complete on all four predictors and the outcome: **n = 81, 36 converters (44.4%)**. 7 injured fish were excluded for a missing 0.5 h or 24 h session (attrition). Events per variable = **9.0** with four predictors. This is a small modeling set, so the predictor set was not expanded.

**Predictors.** `dose` (high_impact = 1); `pre_tau` (τ at t = −1); `z_dtau_0.5` and `z_dtau_24` (τ@0.5 − τ@−1 and τ@24 − τ@−1, z-scored **within dose group**).

**Model.** L2-penalised logistic regression inside a `StandardScaler` pipeline. A fixed, regularised linear classifier limits flexibility at n = 81 and yields coefficients that can be inspected as associations.

**Two safeguards are applied together:**

1. `GroupKFold` on clutch for the outer split (leave-one-clutch-out). A random split would mix clutch-associated observations across training and test partitions.
2. **Nested** CV for the penalty strength C — chosen inside each outer training set only, never on the folds reported below.

The within-dose z-scoring is itself a data-dependent transform, so it is implemented as a pipeline step fitted on training folds only (`modeling.WithinDoseZScorer`) rather than applied to the whole dataset up front. A sensitivity run with whole-dataset z-scoring is reported below.

### Nested comparison table

The model class and validation scheme stay fixed while the input features change.

| Model | Question | Mean fold AUC | SD | Fold range | Pooled OOF AUC [conditional 95% interval] | Brier |
|---|---|---|---|---|---|---|
| (a) baseline_locomotion only | How much does baseline locomotion predict? | **0.416** | 0.090 | 0.323–0.503 | 0.388 [0.265, 0.517] | 0.252 |
| (b) dose only | How much does dose alone predict? | **0.440** | 0.079 | 0.354–0.509 | 0.381 [0.256, 0.502] | 0.251 |
| (c) dose + pre_tau | What does baseline tau add to dose? | **0.677** | 0.132 | 0.588–0.830 | 0.607 [0.484, 0.725] | 0.248 |
| (d) dose + z_dtau | What do post-injury tau changes add to dose? | **0.621** | 0.043 | 0.575–0.659 | 0.609 [0.484, 0.727] | 0.244 |
| (e) delta_tau only (dose-blind) | Does the injury response work without dose context? | **0.591** | 0.072 | 0.516–0.659 | 0.544 [0.417, 0.664] | 0.252 |
| (f) pre_tau + delta_tau (dose-blind) | Do all behavioral features work without dose? | **0.810** | 0.097 | 0.699–0.881 | 0.801 [0.703, 0.886] | 0.181 |
| (g) full 4-predictor model | Full model | **0.849** | 0.080 | 0.758–0.909 | 0.833 [0.736, 0.916] | 0.173 |

Per-fold AUC (held-out clutch), with the C selected inside each fold:

| Model | clutch_A | clutch_B | clutch_C |
|---|---|---|---|
| (a) baseline_locomotion only | 0.503 (C=0.01) | 0.420 (C=0.01) | 0.323 (C=0.01) |
| (b) dose only | 0.458 (C=0.01) | 0.509 (C=0.01) | 0.354 (C=0.01) |
| (c) dose + pre_tau | 0.588 (C=10) | 0.830 (C=0.01) | 0.615 (C=0.01) |
| (d) dose + z_dtau | 0.575 (C=100) | 0.659 (C=0.01) | 0.630 (C=0.01) |
| (e) delta_tau only (dose-blind) | 0.516 (C=100) | 0.659 (C=0.01) | 0.599 (C=0.01) |
| (f) pre_tau + delta_tau (dose-blind) | 0.699 (C=10) | 0.881 (C=100) | 0.849 (C=100) |
| (g) full 4-predictor model | 0.758 (C=1) | 0.909 (C=10) | 0.880 (C=0.03) |

Reading it:

- **(a) baseline_locomotion alone: 0.416.** This single activity measure performs at or below chance out of fold. It does not, by itself, exclude other illness, motor, or severity explanations.
- **(b) dose alone: 0.440.** Also at chance across clutches. The dose effect is not stable between them: in clutch_C high-dose fish convert at 71% versus 43% for low dose, but in clutch_A the gap runs the other way (31% vs 38%). A model trained on two clutches therefore does not transfer to the third. In this dataset, dose alone is a poor cross-clutch classifier; whether it adds information beyond the behavioral trajectory is tested directly by models (f) and (g).
- **(c) dose + pre_tau: 0.677.** The baseline-kinetics model contains predictive signal — but note the fold spread (0.588–0.830) is the widest in the table.
- **(d) dose + z_dtau: 0.621.** The post-injury change-score model carries a comparable amount of predictive signal, and does so much more consistently across folds (SD 0.043).
- **(e) delta_tau only (dose-blind): 0.591.** This ablation omits dose both as a predictor and from preprocessing, testing whether the two raw change scores transfer without injury-severity context.
- **(f) pre_tau + delta_tau (dose-blind): 0.810.** This is the direct dose-necessity ablation: it retains every behavioral feature while omitting dose from both the inputs and preprocessing.
- **(g) the full model: 0.849.** Its mean-fold AUC is +0.172 above the better of (c) and (d), consistent with complementary predictive information from baseline and post-injury features. Adding dose to the complete behavioral model changes mean-fold AUC by +0.039; the (f)-versus-(g) comparison, not the dose-only model, is the test of dose's incremental value. The full model also beats every ablation in every one of the three clutch folds.

The honest reading of (c) versus (d) is that this design cannot cleanly apportion credit between a pre-existing trait and the injury response — with 3 clutches and 36 events their individual AUCs (0.677 vs 0.621) are well inside each other's fold spread. What the table supports is complementary information from baseline and post-injury behavior; it does not establish that explicit dose encoding is necessary.

### Permutation test

Labels were shuffled **within clutch** (preserving each clutch's conversion rate) and the *entire* nested CV rerun 1000 times.

- Null distribution: mean AUC 0.469, SD 0.082, 95th percentile 0.601.
- Observed pooled out-of-fold AUC **0.833** sits at the **100.0th percentile** of the null.
- **p = 0.0010** (p = (1 + #{null ≥ observed}) / (n_perm + 1); the floor at 1000 permutations is 0.0010).
- Using mean-fold AUC as the statistic instead: observed 0.849, p = 0.0010.

Figure: `fig04_permutation_null.png`.

### Random-split comparison

| Split | Mean fold AUC | Pooled OOF AUC |
|---|---|---|
| 5-fold **random** (ignores clutch; comparison only) | 0.865 | 0.856 |
| Leave-one-clutch-out (**reported estimate**) | 0.849 | 0.833 |

In this dataset, the random split's mean fold AUC is **+0.016** above the clutch-held-out estimate. This single comparison does not define a general correction factor.

### Full model coefficients

Refitted on all 81 fish with C = 0.03 (chosen by clutch-held-out CV on the full set — this refit is for interpretation only and contributes nothing to the AUCs above). Coefficients are on the standardised scale; CIs are percentile intervals from a within-clutch stratified bootstrap conditional on the observed clutches (2000 resamples).

| Predictor | β (per SD) | 95% CI | OR per SD | OR 95% CI |
|---|---|---|---|---|
| `dose` | +0.075 | [-0.081, +0.228] | 1.08 | [0.92, 1.26] |
| `pre_tau` | +0.331 \* | [+0.200, +0.430] | 1.39 | [1.22, 1.54] |
| `dtau_0.5` | +0.168 \* | [+0.029, +0.290] | 1.18 | [1.03, 1.34] |
| `dtau_24` | +0.257 \* | [+0.133, +0.362] | 1.29 | [1.14, 1.44] |

Intercept -0.083. \* = bootstrap CI excludes zero.

An unpenalised `statsmodels` Logit fit is stored in `results/tables/step2_unpenalised_logit.csv` for Wald p-values; it is reference material only, since the reported model is penalised.

#### How to read the coefficients

τ is the fitted number of trials over which the visual dark-flash response approaches its floor. It is an empirical summary of response decay and does not by itself identify a neural mechanism. This experiment does not directly measure inhibition, excitation, or a specific circuit, so the coefficients should be interpreted as predictive associations.

- **`pre_tau` (β = +0.331, OR 1.39 per SD).** Within the fitted model, slower pre-injury habituation is associated with the supplied later label. The design cannot determine whether this is causal susceptibility, confounding, or sampling variation.
- **`z_dtau_0.5` (β = +0.168, OR 1.18 per SD).** This is the 30-minute change from the fish's own baseline, centered and scaled within dose. Its coefficient estimates an association conditional on the other predictors; it is not a direct measure of inhibitory gain.
- **`z_dtau_24` (β = +0.257, OR 1.29 per SD).** This is the analogous 24-hour change score. Its inclusion allows the classifier to use the two-timepoint trajectory, but the study does not establish that this trajectory represents a latent epileptogenic period.
- **`dose` (β = +0.075, CI [-0.081, +0.228]).** This coefficient estimates dose's conditional contribution after the three behavioral terms. Its interval and the direct comparison of models (f) and (g) should be used together; the strong dose-blind model means dose cannot be described as necessary for prediction. Because the model contains no dose-by-change interaction, it is not a formal moderation analysis.

### Classification and calibration (out-of-fold)

| | Predicted non-converter | Predicted converter |
|---|---|---|
| **Non-converter** | 38 | 7 |
| **Converter** | 12 | 24 |

**Total out-of-fold accuracy = 76.5%** at the default 0.50 threshold (62 of 81 injured fish classified correctly). Balanced accuracy 75.6%, sensitivity 66.7%, specificity 84.4%, PPV 77.4%, NPV 76.0%. Brier score 0.173 (0.25 would be uninformative at this prevalence).

At the Youden-optimal operating point (threshold 0.472) accuracy is 81.5% (balanced 81.7%, sensitivity 83.3%, specificity 80.0%). That threshold was chosen on these same out-of-fold predictions, so it is mildly optimistic and is quoted only to show the achievable operating range. **The 0.50 figure is the one to cite**, and AUC — which is threshold-free — remains the primary metric.

Calibration matters for a candidate screening assay: a model that ranks fish correctly but reports 0.9 for a fish that converts 60% of the time would misallocate any intervention trial built on it. The out-of-fold calibration curve (`fig05_confusion_calibration.png`, centre panel) tracks the diagonal within the resolution 5 quantile bins allow at n = 81.

*Sensitivity:* z-scoring within dose on the whole dataset instead of per training fold gives mean fold AUC 0.853 vs 0.849 fold-safe — the leakage-free implementation costs essentially nothing.

Figures: `fig03_roc_nested_comparison.png`, `fig05_confusion_calibration.png`, `fig06_coefficients.png`, `fig07_fold_auc_by_model.png`.

---

## Step 3 — Orthogonal molecular validation: paired c-fos pools

### Why this is an orthogonal validation

Steps 1–2 use behavioral measurements. Step 3 uses quantitative PCR of an immediate early gene in a separate `cf_*` cohort whose fish do not enter the prediction model. It is the project's orthogonal molecular validation because it tests supplied risk strata, described as behavior-derived, using a different cohort and assay modality. The analysis accepts the workbook's `risk_pool` assignments rather than reconstructing them from a documented scoring and selection algorithm, so it is not external validation of the prediction model.

`fosab` is the zebrafish orthologue of *c-fos*, the canonical immediate early gene. Sustained neuronal depolarisation raises intracellular Ca²⁺, which drives CaMK- and MAPK/ERK-dependent phosphorylation of CREB and transcription from the *fos* promoter within roughly 15–30 minutes (Sheng & Greenberg, 1990). c-fos transcript level is therefore a molecular integrator of recent network activity, and it is the standard readout for mapping seizure-recruited circuits — including in the original characterisation of chemically induced seizures in larval zebrafish (Baraban et al., 2005). Normalisation is against `rpl13a`, one of the reference genes validated as stable across zebrafish development (Tang et al., 2007), by the 2^−ΔΔCt method (Livak & Schmittgen, 2001).

The orthogonal-validation hypothesis is that supplied high-risk pools will show higher c-fos fold change than their matched low-risk pools. A difference provides cross-modal molecular concordance with the risk strata, but does not identify its cause or establish electrographic epilepsy.

### Statistical treatment

The 18 pools are **9 matched pairs** — one high_risk and one low_risk pool per (group × clutch) cell — and are analysed as such. Treating them as 18 independent units would roughly double the nominal degrees of freedom and ignore the plate/clutch matching.

The nine pairs are not nine fully independent biological replicates: three group-level pairs come from each of only three clutches. The all-pair tests below are therefore nominal. A clutch-averaged sensitivity analysis, which leaves n = 3, is reported separately.

| Contrast | Pairs | Mean Δ (high − low) | 95% CI | Paired t | p | Cohen's dz | Wilcoxon p |
|---|---|---|---|---|---|---|---|
| All pairs | 9 | +0.2856 | [+0.0083, +0.5629] | t(8) = +2.375 | = 0.0449 | +0.792 | = 0.0547 |
| Injured only | 6 | +0.3908 | [-0.0257, +0.8073] | t(5) = +2.412 | = 0.0607 | +0.985 | = 0.0938 |
| **Sham only (control)** | 3 | +0.0750 | [-0.3336, +0.4836] | t(2) = +0.790 | = 0.5124 | +0.456 | = 0.5000 |
| All pairs, log2 scale | 9 | +0.3504 | [+0.0453, +0.6556] | t(8) = +2.648 | = 0.0293 | +0.883 | = 0.0273 |

### Clutch-averaged sensitivity

Each row first averages the sham, low-impact, and high-impact pair differences within a clutch, then tests the three clutch means against zero. This conservative sensitivity treats clutch as the independent biological unit.

| Scale | Clutches | Mean high − low | 95% CI | one-sample t | p | Wilcoxon p |
|---|---|---|---|---|---|---|
| Raw fold change | 3 | +0.2856 | [-0.2252, +0.7963] | t(2) = +2.406 | = 0.1379 | = 0.2500 |
| Log2 fold change | 3 | +0.3504 | [-0.1784, +0.8792] | t(2) = +2.851 | = 0.1041 | = 0.2500 |

The nominal nine-pair log2 comparison has p = 0.0293; after averaging within clutch, the log2 sensitivity has p = 0.1041. With only three independent clutches, the clutch-level validation estimate remains imprecise; the nominal pair-level result and this sensitivity should be interpreted together.

Direction consistency: 5/6 injured pairs have high_risk > low_risk (exact binomial sign test, p = 0.2188).

The three sham pairs have mean difference +0.0750 (p = 0.5124). This small, non-significant comparison does not rule out a pooling, plate, or bin-assignment artifact.

Read the injured-only row carefully: the point estimate is the largest of the three (+0.3908, dz = 0.98) but with only 6 pairs it does not reach significance on its own (p = 0.0607). The original all-pairs comparison is nominal; the injured and sham rows show where its point estimate sits, not two independent confirmations of it.

No regression of c-fos on a continuous risk score was run because continuous scores and a reproducible pool-selection rule are not present in the workbook. The normalized membership table verifies the recorded four-fish pools, not the validity of their supplied risk labels.

### Interpretation

In the nominal all-pair comparison, supplied high-risk pools have **27.5% higher c-fos fold change** relative to `rpl13a` than supplied low-risk pools processed on the same plate (geometric mean ratio 1.275, 95% CI [1.032, 1.575], back-transformed from the paired log2 analysis). Because three pairs share each clutch, the clutch-averaged log2 result (p = 0.1041) and the missing risk-assignment protocol set the interpretation boundary.

This constitutes the project's orthogonal molecular validation of the supplied risk stratification. It is not external validation of the classifier or seizure endpoint. The bulk measurement on pooled tissue cannot localise the signal to a cell type or region, and raw Ct values, technical-replicate results, amplification efficiencies, and qPCR quality-control records are not supplied.

Figure: `fig08_cfos_paired.png` (9 paired lines, plus within-pair differences by group).

---

## Step 4 — Descriptive results

### No baseline group difference detected

At the pre-injury baseline (t = −1), no group difference in τ was detected: F(2,130) = 0.543, p = 0.5825, η² = 0.0083; Kruskal–Wallis H = 0.522, p = 0.7703. Baseline locomotion likewise (p = 0.1218). The longitudinal changes are measured within fish, but non-significance does not establish baseline equivalence or prove random assignment.

### τ moves in opposite directions by dose

| Group | t = -1 h | t = 0.5 h | t = 1 h | t = 5 h | t = 24 h |
|---|---|---|---|---|---|
| Sham | 6.33 ± 0.27 | 6.19 ± 0.45 | 6.03 ± 0.43 | 7.45 ± 0.62 | 7.23 ± 0.59 |
| Low impact | 6.46 ± 0.36 | 10.21 ± 0.61 | 9.38 ± 0.68 | 8.17 ± 0.67 | 7.02 ± 0.48 |
| High impact | 6.02 ± 0.28 | 3.28 ± 0.33 | 3.84 ± 0.48 | 4.52 ± 0.46 | 4.42 ± 0.40 |

(mean ± SEM, trials to habituate)

At 0.5 h, Δτ from each fish's own baseline is **+3.74** trials in low_impact and **-2.72** in high_impact — slower versus faster fitted response decay. Welch t = 10.00, p < 0.0001, Cohen's d = 2.19.

**Dose changes the biological interpretation of Δτ, but its predictive necessity is an empirical question.** Model (e) tests the two raw changes without dose; model (f) adds baseline τ while remaining fully dose-blind; and the comparison of (f) with (g) isolates dose's incremental contribution.

The observed directions have different empirical meanings:

- **Low dose → τ rises.** The fitted startle response decays more slowly across trials, which is consistent with impaired habituation.
- **High dose → τ falls.** The fitted response decays more quickly. This could reflect faster habituation, depressed responsiveness, fatigue, or another process; the present measurements do not distinguish these explanations.

The dose groups move the same fitted scalar in opposite directions. A model given Δτ without dose is asked to treat +4 trials and −3 trials as opposite kinds of evidence when they may reflect different processes. Encoding dose is one way to represent that context, but the strong dose-blind behavioral ablation shows that baseline and trajectory can recover substantial predictive information without it.

The group means also remain shifted from baseline at 24 h (low impact +0.57, high impact -1.60 trials from each fish's own pre-injury session). This documents persistent behavioral change; it does not by itself identify circuit remodeling or an epileptogenic latent period.

### Converters vs non-converters

Trajectories are shown separately per dose in `fig11_converter_trajectories.png`, for the same reason: pooling the doses would average away the contrast. Per-dose converter/non-converter contrasts at each timepoint are in `all_statistics.csv` under `converter_contrast`.

### Operational metrics

| Metric | Value |
|---|---|
| Sessions | 15 (3 clutches × 5 timepoints) |
| Fish-sessions recorded | 650 |
| Fish per hour | 33.3 ± 2.4 |
| Operator minutes per fish | 0.90 ± 0.14 |
| Consumables cost per fish | $0.094 |
| Total consumables | $60.76 |
| Total operator time | 9.7 h |
| Attrition | 9/133 = 6.8% [3.6%, 12.4%] |

`consumables_cost_usd` is treated as a per-session total; it tracks `n_fish_recorded` almost exactly (r = 0.999), which is consistent with that reading.

Figures: `fig09_habituation_curves.png`, `fig10_tau_by_timepoint.png`, `fig11_converter_trajectories.png`, `fig12_operations.png`.

---

## PTZ challenge — secondary and underpowered

Pentylenetetrazol is a non-competitive GABAₐ receptor antagonist: it binds at the picrotoxin site, reduces chloride conductance, and removes inhibitory brake from the network. In larval zebrafish it produces stereotyped, dose-dependent seizure behaviour with electrographic correlates (Baraban et al., 2005), which makes it the standard pharmacological probe of seizure susceptibility. It is included here as an exploratory, pharmacological comparison rather than as confirmation of the behavioral model.

| Group | Seized / n | Proportion | Wilson 95% CI | Median latency (s) |
|---|---|---|---|---|
| Sham | 3/12 | 0.25 | [0.09, 0.53] | 1800 |
| Low impact | 7/12 | 0.58 | [0.32, 0.81] | 720 |
| High impact | 7/10 | 0.70 | [0.40, 0.89] | 607 |

χ²(2) = 4.933, p = 0.0849, Cramér's V = 0.381, total n = 34.

> **Explicit statement of limitation.** This probe is underpowered. With 34 fish split across three groups, post-hoc power for the observed sham-versus-injured difference (Cohen's h = 0.80) is only **61%**. The minimum expected cell count is 5.00. It is reported as a directional check only, and **no conclusion in this report rests on it.**

Figure: `fig13_ptz.png`.

---

## Assumptions, checks and limitations

All assumption checks are printed to stdout during the run, tagged `[PASS]` or `[FLAG]`, and the underlying statistics are in `all_statistics.csv`. In summary:

- Curve fitting: 100.0% convergence, median R² 0.796, 1 session at a parameter bound.
- Collinearity among the four predictors: max |r| = 0.740; VIFs are in `all_statistics.csv` under `assumptions`.
- Separation: the unpenalised Logit converged with max |β| = 7.73, so the ridge penalty is not masking complete separation.
- **Flagged:** baseline τ deviates from normality in low impact, high impact (Shapiro–Wilk p ≤ 0.05). The ANOVA above is therefore backed by a Kruskal–Wallis test (p = 0.7703), which agrees. ANOVA is robust to this at these group sizes, but the non-parametric result is the one to quote if the distributional assumption matters to a reader.
- Paired c-fos differences: Shapiro–Wilk reported for each contrast; Wilcoxon is given alongside every paired t-test as a distribution-free check.
- Linearity of the logit is assumed for the three continuous predictors; quartile event rates are printed as a coarse check.

Real limitations, stated plainly:

- **The model includes n = 81 fish with 36 events.** EPV = 9.0 is low. The CIs on the AUC and on every coefficient are wide, and they are reported rather than smoothed over.
- **Three clutches means three outer folds.** The fold-to-fold spread is estimated from three numbers; the SD across folds should be read as indicative, not precise. The AUC interval resamples fish within the observed clutches while holding out-of-fold scores fixed, so it is conditional on these clutches and does not include model-selection uncertainty.
- **7 injured fish were dropped** for a missing post-injury session. This is complete-case analysis, and the missingness mechanism has not been modelled.
- **The binary outcome is supplied, not derived by code.** The workbook does not document a threshold, recording duration, blinded scoring procedure, or electrographic confirmation.
- **The c-fos comparison uses pooled material** — 9 group-level pairs, 4 larvae per pool, but only 3 clutches. Risk bins are supplied rather than reproducibly generated, and raw qPCR Ct and QC data are absent.
- **The workbook represents 253 unique animal IDs across three non-overlapping cohorts** (133 followed, 86 c-fos, and 34 PTZ); 72 of the c-fos animals entered pools.
- **Data provenance is unresolved in the repository.** Source recordings, instrument exports, dated protocols, approval identifiers, and a label-derivation audit trail are not included.
- **The feature and model choices are not preregistered.** Nested cross-validation covers penalty tuning, not uncertainty from post-hoc feature, timepoint, or model selection.
- **PTZ is underpowered** (above), and the study is not designed to support any claim from it.

## Conclusion

In the supplied dataset, baseline and post-injury startle-habituation features distinguish injured larvae with versus without the supplied `converted` label, with a cross-validated AUC of 0.833 (conditional 95% interval [0.736, 0.916]), a total out-of-fold accuracy of 76.5%, and p = 0.0010 against a null that reruns the entire nested cross-validation. The dose-blind three-behavior model reaches a mean-fold AUC of 0.810, so explicit dose encoding is not necessary for the observed separation. The c-fos experiment supplies orthogonal molecular validation through cross-modal concordance with the supplied risk strata; it is not external validation of the classifier.

The defensible conclusion is narrow: **startle-habituation kinetics carry internally cross-validated information about a later behavioral label in this workbook.** Prospective replication with documented outcome scoring, independent clutches, source-data provenance, and electrographic confirmation is required before calling the assay a biomarker of post-traumatic epilepsy.

## Background literature

The contextual discussion above draws on the following. This is a background reading list, not a citation of results generated here.

> These references were already present in the source repository. A student preparing an ISEF submission must independently read, verify, and format every citation under the current AI-use rules.

1. Annegers JF, Hauser WA, Coan SP, Rocca WA (1998). A population-based study of seizures after traumatic brain injuries. *New England Journal of Medicine* 338:20–24.
2. Baraban SC, Taylor MR, Castro PA, Baier H (2005). Pentylenetetrazole induced changes in zebrafish behavior, neural activity and c-fos expression. *Neuroscience* 131:759–768.
3. Burgess HA, Granato M (2007). Sensorimotor gating in larval zebrafish. *Journal of Neuroscience* 27:4984–4994.
4. Hunt RF, Boychuk JA, Smith BN (2013). Neural circuit mechanisms of post-traumatic epilepsy. *Frontiers in Cellular Neuroscience* 7:89.
5. Livak KJ, Schmittgen TD (2001). Analysis of relative gene expression data using real-time quantitative PCR and the 2^−ΔΔCt method. *Methods* 25:402–408.
6. Marsden KC, Granato M (2015). In vivo Ca²⁺ imaging reveals that decreased dendritic excitability drives startle habituation. *Cell Reports* 13:1733–1740.
7. Peduzzi P, Concato J, Kemper E, Holford TR, Feinstein AR (1996). A simulation study of the number of events per variable in logistic regression analysis. *Journal of Clinical Epidemiology* 49:1373–1379.
8. Sheng M, Greenberg ME (1990). The regulation and function of c-fos and other immediate early genes in the nervous system. *Neuron* 4:477–485.
9. Sloviter RS (1991). Permanently altered hippocampal structure, excitability, and inhibition after experimental status epilepticus in the rat: the 'dormant basket cell' hypothesis. *Hippocampus* 1:41–66.
10. Tang R, Dodd A, Lai D, McNabb WC, Love DR (2007). Validation of zebrafish (*Danio rerio*) reference genes for quantitative real-time RT-PCR normalization. *Acta Biochimica et Biophysica Sinica* 39:384–390.
11. Varma S, Simon R (2006). Bias in error estimation when using cross-validation for model selection. *BMC Bioinformatics* 7:91.
12. Wolman MA, Jain RA, Liss L, Granato M (2011). Chemical modulation of memory formation in larval zebrafish. *PNAS* 108:15468–15473.

## Reproducing

```bash
pip install -r requirements.txt
python run_all.py
```

Seed `20260809`. Outputs: `results/figures/*.png` (300 dpi), `results/all_statistics.csv`, `results/tables/*.csv`, and this file.
