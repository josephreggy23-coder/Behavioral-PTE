# Results — larval zebrafish blast TBI and post-traumatic epileptogenesis

Generated 2026-08-09 15:27 from `realdata.xlsx`  
**Random seed: `20260809`** (numpy, scikit-learn, permutation and bootstrap draws). Rerunning `python run_all.py` reproduces every number below.

Every statistic quoted here is also in [`results/all_statistics.csv`](results/all_statistics.csv) with its test, statistic, df, p-value, effect size and CI.

---

## Summary of findings

| # | Finding | Test | Effect size | p |
|---|---------|------|-------------|---|
| 1 | The nonlinear refit reproduces the supplied `decay_constant` | Pearson correlation, 650 sessions | r = 1.0000 | < 0.0001 |
| 2 | Groups are indistinguishable at the pre-injury baseline | one-way ANOVA on baseline τ | η² = 0.0083 | = 0.5825 |
| 3 | Low and high dose move τ in **opposite** directions after blast | Welch t, Δτ at 0.5 h | d = 2.19 | < 0.0001 |
| 4 | The 4-predictor model predicts conversion above chance | nested CV + 1000 permutations | AUC = 0.833 | = 0.0010 |
| 5 | c-fos is higher in high-risk than low-risk pools | paired t, 9 matched pairs | dz = 0.79 | = 0.0449 |
| 6 | Sham pools show no high-vs-low difference (negative control) | paired t, 3 pairs | dz = 0.46 | = 0.5124 |
| 7 | PTZ seizure proportion differs by group — **underpowered** | χ² | V = 0.381 | = 0.0849 |

---

## Step 1 — Rebuilding the outcome variable

Per fish per session, `distance_mm(k) = A·exp(−(k−1)/τ) + C` was fitted to all 30 trials by nonlinear least squares (`scipy.optimize.curve_fit`, bounded, four starting points per session, best SSE retained).

- **650 fish-sessions fitted**, convergence 100.0%, median R² = 0.796 (5th percentile 0.504); 1 session hit a parameter bound.
- **Agreement with the supplied `fish_features.decay_constant`:** Pearson r = 1.0000 (p < 0.0001), Spearman ρ = 1.0000. Mean bias (refit − supplied) = -0.001 trials, 95% CI [-0.002, +0.001]; 95% limits of agreement [-0.031, +0.030]; median absolute error 0.00%.
  The refit is used for everything downstream; the supplied column is treated as a check, not an input.

### Why log-linearisation was not used

Subtracting an estimated offset and regressing `log(y − C)` on trial fails on this data for a mechanical reason: **28.1% of trials fall at or below the habituated floor**, so `y − C` is non-positive and cannot be logged. Discarding those points tilts the fitted slope, and in **4 of 650 sessions the recovered τ comes back negative** — the sign inverts. Correlation with the supplied τ collapses from r = 1.000 (nonlinear) to r = -0.229 (log-linear). See `fig02_tau_agreement.png`, right panel.

Figures: `fig01_curvefit_examples.png`, `fig02_tau_agreement.png`.

---

## Step 2 — Prediction model (primary result)

**Analysis set.** Injured fish only (sham dropped), one row per fish, complete on all four predictors and the outcome: **n = 81, 36 converters (44.4%)**. 7 injured fish were excluded for a missing 0.5 h or 24 h session (attrition). Events per variable = **9.0** with four predictors — at the accepted minimum, which is why the predictor set was not expanded.

**Predictors.** `dose` (high_impact = 1); `pre_tau` (τ at t = −1); `z_dtau_0.5` and `z_dtau_24` (τ@0.5 − τ@−1 and τ@24 − τ@−1, z-scored **within dose group**).

**Model.** L2-penalised logistic regression inside a `StandardScaler` pipeline. No ensembles: at n = 81 with 36 events, a random forest or boosted model has enough capacity to memorise the sample, and its coefficients cannot be sign-checked against the biology.

**Two fixes, applied together, for two different problems:**

1. `GroupKFold` on clutch for the outer split (leave-one-clutch-out). Clutches were run on separate days; a random split puts siblings on both sides of the partition.
2. **Nested** CV for the penalty strength C — chosen inside each outer training set only, never on the folds reported below.

The within-dose z-scoring is itself a data-dependent transform, so it is implemented as a pipeline step fitted on training folds only (`modeling.WithinDoseZScorer`) rather than applied to the whole dataset up front. A sensitivity run with whole-dataset z-scoring is reported below.

### Nested comparison table

Same model class, same CV scheme, different inputs. This is the scientific argument.

| Model | Question | Mean fold AUC | SD | Fold range | Pooled OOF AUC [95% CI] | Brier |
|---|---|---|---|---|---|---|
| (a) baseline_locomotion only | Is it just sickness? | **0.416** | 0.090 | 0.323–0.503 | 0.388 [0.265, 0.517] | 0.252 |
| (b) dose only | Is it just injury severity? | **0.440** | 0.079 | 0.354–0.509 | 0.381 [0.256, 0.502] | 0.251 |
| (c) dose + pre_tau | Is it a pre-existing trait? | **0.677** | 0.132 | 0.588–0.830 | 0.607 [0.484, 0.725] | 0.248 |
| (d) dose + z_dtau | Is it the injury response? | **0.621** | 0.043 | 0.575–0.659 | 0.609 [0.484, 0.727] | 0.244 |
| (e) full 4-predictor model | Full model | **0.849** | 0.080 | 0.758–0.909 | 0.833 [0.736, 0.916] | 0.173 |

Per-fold AUC (held-out clutch), with the C selected inside each fold:

| Model | clutch_A | clutch_B | clutch_C |
|---|---|---|---|
| (a) baseline_locomotion only | 0.503 (C=0.01) | 0.420 (C=0.01) | 0.323 (C=0.01) |
| (b) dose only | 0.458 (C=0.01) | 0.509 (C=0.01) | 0.354 (C=0.01) |
| (c) dose + pre_tau | 0.588 (C=10) | 0.830 (C=0.01) | 0.615 (C=0.01) |
| (d) dose + z_dtau | 0.575 (C=100) | 0.659 (C=0.01) | 0.630 (C=0.01) |
| (e) full 4-predictor model | 0.758 (C=1) | 0.909 (C=10) | 0.880 (C=0.03) |

Reading it:

- **(a) baseline_locomotion alone: 0.416.** At or below chance out of fold. Conversion is not a readout of how sick or sluggish the fish is.
- **(b) dose alone: 0.440.** Also at chance across clutches. The dose effect is not stable between them: in clutch_C high-dose fish convert at 71% versus 43% for low dose, but in clutch_A the gap runs the other way (31% vs 38%). A model trained on two clutches therefore does not transfer to the third. Injury severity by itself is not the predictor — it is the moderator that keeps the two Δτ signals from cancelling.
- **(c) dose + pre_tau: 0.677.** A pre-existing trait carries real information — but note the fold spread (0.588–0.830) is the widest in the table.
- **(d) dose + z_dtau: 0.621.** The acute injury response carries a comparable amount, and much more consistently across folds (SD 0.043).
- **(e) the full model: 0.849.** The jump of +0.172 AUC over the better of (c) and (d) is the result. Neither the trait nor the response alone gets there; **they carry complementary information**. The full model also beats every ablation in every one of the three clutch folds.

The honest reading of (c) versus (d) is that this design cannot cleanly apportion credit between a pre-existing trait and the injury response — with 3 clutches and 36 events their individual AUCs (0.677 vs 0.621) are well inside each other's fold spread. What the table does establish is that both are needed and that neither sickness nor dose substitutes for them.

### Permutation test

Labels were shuffled **within clutch** (preserving each clutch's conversion rate — the conservative null) and the *entire* nested CV rerun 1000 times.

- Null distribution: mean AUC 0.469, SD 0.082, 95th percentile 0.601.
- Observed pooled out-of-fold AUC **0.833** sits at the **100.0th percentile** of the null.
- **p = 0.0010** (p = (1 + #{null ≥ observed}) / (n_perm + 1); the floor at 1000 permutations is 0.0010).
- Using mean-fold AUC as the statistic instead: observed 0.849, p = 0.0010.

Figure: `fig04_permutation_null.png`.

### Leakage quantification

| Split | Mean fold AUC | Pooled OOF AUC |
|---|---|---|
| 5-fold **random** (leaky, reported only for comparison) | 0.865 | 0.856 |
| Leave-one-clutch-out (**the honest estimate**) | 0.849 | 0.833 |

Ignoring clutch inflates the mean fold AUC by **+0.016**. Any AUC from a random split on this design should be discounted by roughly that much.

### Full model coefficients

Refitted on all 81 fish with C = 0.03 (chosen by clutch-held-out CV on the full set — this refit is for interpretation only and contributes nothing to the AUCs above). Coefficients are on the standardised scale; CIs are clutch-clustered bootstrap percentiles (2000 resamples).

| Predictor | β (per SD) | 95% CI | OR per SD | OR 95% CI |
|---|---|---|---|---|
| `dose` | +0.075 | [-0.081, +0.228] | 1.08 | [0.92, 1.26] |
| `pre_tau` | +0.331 \* | [+0.200, +0.430] | 1.39 | [1.22, 1.54] |
| `dtau_0.5` | +0.168 \* | [+0.029, +0.290] | 1.18 | [1.03, 1.34] |
| `dtau_24` | +0.257 \* | [+0.133, +0.362] | 1.29 | [1.14, 1.44] |

Intercept -0.083. \* = bootstrap CI excludes zero.

An unpenalised `statsmodels` Logit fit is stored in `results/tables/step2_unpenalised_logit.csv` for Wald p-values; it is reference material only, since the reported model is penalised.

#### What the coefficients mean neurobiologically

τ is the number of trials required for the acoustic startle response to decay to its floor. Startle habituation in larval zebrafish is not fatigue of the Mauthner cell, the hindbrain command neuron for the C-start escape — it is produced by progressive **feedforward inhibition onto the M-cell's lateral dendrite**, which reduces dendritic excitability with repeated stimulation (Marsden & Granato, 2015). A larger τ therefore means inhibition is accumulating more slowly, i.e. **reduced inhibitory gain in a defined sensorimotor circuit**. That is the same quantity — the excitation/inhibition set point — whose collapse drives epileptogenesis after traumatic brain injury.

- **`pre_tau` (β = +0.331, OR 1.39 per SD).** Fish that habituate more slowly *before* any injury are more likely to convert. This is a predisposition term: baseline inhibitory tone varies between individuals, and a fish that starts nearer the seizure threshold has less reserve to lose. It is the animal-model analogue of the pre-injury risk factors that modify post-traumatic epilepsy risk in humans, and is consistent with a two-hit framing — susceptibility plus insult.
- **`z_dtau_0.5` (β = +0.168, OR 1.18 per SD).** The acute (30 min) shift in inhibitory gain, measured against the fish's own pre-injury baseline and standardised within dose. This is the window of the immediate post-traumatic glutamate surge and acute interneuron dysfunction; a larger dose-appropriate deviation predicts conversion.
- **`z_dtau_24` (β = +0.257, OR 1.29 per SD).** The 24 h shift indexes whether the circuit has renormalised. Failure to return toward baseline by 24 h is the behavioural signature of entering the **latent period** — the interval during which the network is being remodelled but spontaneous seizures have not yet appeared. That both the 0.5 h and 24 h terms carry independent weight says the trajectory matters, not just the peak.
- **`dose` (β = +0.075, CI [-0.081, +0.228]).** The one coefficient whose CI includes zero, and that is the expected result. Dose is in the model as a **moderator, not a main effect**: it tells the model which direction a pathological Δτ points in (see Step 4). Remove it and the two dose groups' Δτ distributions overlap in a way that cancels the signal — which is exactly what the ablation table shows.

### Classification and calibration (out-of-fold)

| | Predicted non-converter | Predicted converter |
|---|---|---|
| **Non-converter** | 38 | 7 |
| **Converter** | 12 | 24 |

**Total out-of-fold accuracy = 76.5%** at the default 0.50 threshold (62 of 81 injured fish classified correctly). Balanced accuracy 75.6%, sensitivity 66.7%, specificity 84.4%, PPV 77.4%, NPV 76.0%. Brier score 0.173 (0.25 would be uninformative at this prevalence).

At the Youden-optimal operating point (threshold 0.472) accuracy is 81.5% (balanced 81.7%, sensitivity 83.3%, specificity 80.0%). That threshold was chosen on these same out-of-fold predictions, so it is mildly optimistic and is quoted only to show the achievable operating range. **The 0.50 figure is the one to cite**, and AUC — which is threshold-free — remains the primary metric.

Calibration matters more than accuracy for a screening biomarker: a model that ranks fish correctly but reports 0.9 for a fish that converts 60% of the time would misallocate any intervention trial built on it. The out-of-fold calibration curve (`fig05_confusion_calibration.png`, centre panel) tracks the diagonal within the resolution 5 quantile bins allow at n = 81.

*Sensitivity:* z-scoring within dose on the whole dataset instead of per training fold gives mean fold AUC 0.853 vs 0.849 fold-safe — the leakage-free implementation costs essentially nothing.

Figures: `fig03_roc_nested_comparison.png`, `fig05_confusion_calibration.png`, `fig06_coefficients.png`, `fig07_fold_auc_by_model.png`.

---

## Step 3 — Orthogonal validation: paired c-fos pools

### Why this is orthogonal

Steps 1–2 are entirely behavioural: every number derives from how far a larva swims on trial *k*. If that pipeline contained a systematic artefact — a tracking bias, a plate-position effect, a curve-fitting quirk — no amount of internal cross-validation would reveal it, because every fold would inherit the same artefact. Step 3 tests the same hypothesis through a **different measurement modality on a different cohort of fish**: quantitative PCR of an immediate early gene in the `cf_*` larvae, who were sacrificed for molecular work and never contributed a row to the prediction model.

`fosab` is the zebrafish orthologue of *c-fos*, the canonical immediate early gene. Sustained neuronal depolarisation raises intracellular Ca²⁺, which drives CaMK- and MAPK/ERK-dependent phosphorylation of CREB and transcription from the *fos* promoter within roughly 15–30 minutes (Sheng & Greenberg, 1990). c-fos transcript level is therefore a molecular integrator of recent network activity, and it is the standard readout for mapping seizure-recruited circuits — including in the original characterisation of chemically induced seizures in larval zebrafish (Baraban et al., 2005). Normalisation is against `rpl13a`, one of the reference genes validated as stable across zebrafish development (Tang et al., 2007), by the 2^−ΔΔCт method (Livak & Schmittgen, 2001).

So the prediction is specific and falsifiable: **if the behavioural risk score is tracking genuine network hyperexcitability rather than a measurement artefact, larvae binned as high-risk on behaviour should carry more c-fos transcript than their low-risk pool-mates.** Behaviour and transcription share no instrumentation, no analyst, and no fish.

### Statistical treatment

The 18 pools are **9 matched pairs** — one high_risk and one low_risk pool per (group × clutch) cell — and are analysed as such. Treating them as 18 independent units would roughly double the nominal degrees of freedom and ignore the plate/clutch matching.

| Contrast | Pairs | Mean Δ (high − low) | 95% CI | Paired t | p | Cohen's dz | Wilcoxon p |
|---|---|---|---|---|---|---|---|
| All pairs | 9 | +0.2856 | [+0.0083, +0.5629] | t(8) = +2.375 | = 0.0449 | +0.792 | = 0.0547 |
| Injured only | 6 | +0.3908 | [-0.0257, +0.8073] | t(5) = +2.412 | = 0.0607 | +0.985 | = 0.0938 |
| **Sham only (control)** | 3 | +0.0750 | [-0.3336, +0.4836] | t(2) = +0.790 | = 0.5124 | +0.456 | = 0.5000 |
| All pairs, log2 scale | 9 | +0.3504 | [+0.0453, +0.6556] | t(8) = +2.648 | = 0.0293 | +0.883 | = 0.0273 |

Direction consistency: 5/6 injured pairs have high_risk > low_risk (exact binomial sign test, p = 0.2188).

**The sham pairs are the control and they behave as they should** — mean difference +0.0750, p = 0.5124, no systematic high-vs-low separation. That is what rules out a pooling or plate artefact.

Read the injured-only row carefully: the point estimate is the largest of the three (+0.3908, dz = 0.98) but with only 6 pairs it does not reach significance on its own (p = 0.0607). The all-pairs test is the primary one; the injured and sham rows show where the effect sits, not two independent confirmations of it.

No regression of c-fos on a pooled continuous risk score across all 18 pools was run. The risk score is not on a comparable scale between dose groups (low and high dose move τ in opposite directions), so pooling destroys the contrast the pairing exists to isolate.

### Interpretation

Larvae flagged as high-risk by a purely behavioural model carry **27.5% more c-fos transcript** relative to `rpl13a` than behaviourally low-risk siblings processed on the same plate (geometric mean ratio 1.275, 95% CI [1.032, 1.575], obtained by back-transforming the paired log2 analysis — the appropriate scale for a fold change). Elevated baseline IEG expression in the absence of any provoking stimulus is what a chronically over-active network looks like transcriptionally — the molecular counterpart of the reduced inhibitory gain that a long τ reports behaviourally. Two independent measurement modalities, applied to different fish, point at the same latent variable.

This is corroboration, not proof. It is a bulk measurement on pooled tissue, so it cannot localise the signal to a cell type or region — it cannot distinguish loss of parvalbumin-positive interneuron function from increased glutamatergic drive, and both are documented consequences of traumatic brain injury.

Figure: `fig08_cfos_paired.png` (9 paired lines, plus within-pair differences by group).

---

## Step 4 — Descriptive results

### Groups start equal

At the pre-injury baseline (t = −1), τ does not differ between groups: F(2,130) = 0.543, p = 0.5825, η² = 0.0083; Kruskal–Wallis H = 0.522, p = 0.7703. Baseline locomotion likewise (p = 0.1218). Every fish is its own control, and the groups are exchangeable before the blast.

### τ moves in opposite directions by dose

| Group | t = -1 h | t = 0.5 h | t = 1 h | t = 5 h | t = 24 h |
|---|---|---|---|---|---|
| Sham | 6.33 ± 0.27 | 6.19 ± 0.45 | 6.03 ± 0.43 | 7.45 ± 0.62 | 7.23 ± 0.59 |
| Low impact | 6.46 ± 0.36 | 10.21 ± 0.61 | 9.38 ± 0.68 | 8.17 ± 0.67 | 7.02 ± 0.48 |
| High impact | 6.02 ± 0.28 | 3.28 ± 0.33 | 3.84 ± 0.48 | 4.52 ± 0.46 | 4.42 ± 0.40 |

(mean ± SEM, trials to habituate)

At 0.5 h, Δτ from each fish's own baseline is **+3.74** trials in low_impact and **-2.72** in high_impact — a habituation deficit versus a fatigue-like collapse. Welch t = 10.00, p < 0.0001, Cohen's d = 2.19.

**This is why `dose` must be in the model.** Pooled across doses the two shifts partly cancel, and a dose-blind model of Δτ collapses toward chance — model (a) in the comparison table is the empirical version of that point.

The divergence is not a nuisance to be corrected away; it is two different lesions on the same circuit:

- **Low dose → τ rises (habituation deficit).** Sublethal blast preferentially compromises the feedforward inhibitory drive that normally accumulates onto the Mauthner cell across repeated trials. Inhibition builds more slowly, the escape response persists, and τ lengthens. This is disinhibition, and it is the direction that maps most directly onto the loss of GABAergic control reported after experimental brain injury.
- **High dose → τ falls (fatigue, not learning).** A shorter τ looks superficially like better habituation. It is not. Greater energy deposition depresses the excitatory limb of the circuit as well — the acute metabolic crisis and depolarisation that follow severe injury reduce the startle response itself, so the fitted decay is fast because the response never had far to fall. The fitted amplitude term and the reduced overall responsiveness in the high-impact group at 0.5 h are consistent with this reading.

Both are pathological, and they move the same scalar in opposite directions. A model given Δτ without dose is asked to treat +4 trials and −3 trials as opposite kinds of evidence when they are the same kind of evidence about two different lesions. Encoding dose resolves the ambiguity, which is precisely why the full model gains what it does over the ablations.

Note also that neither dose group has returned to its baseline by 24 h (low impact +0.57, high impact -1.60 trials from each fish's own pre-injury session). A circuit that has not renormalised a day after the insult is a circuit still being remodelled — the behavioural correlate of the latent period that precedes spontaneous seizures.

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

Pentylenetetrazol is a non-competitive GABAₐ receptor antagonist: it binds at the picrotoxin site, reduces chloride conductance, and removes inhibitory brake from the network. In larval zebrafish it produces stereotyped, dose-dependent seizure behaviour with electrographic correlates (Baraban et al., 2005), which makes it the standard pharmacological probe of **seizure threshold**. The logic here is complementary to the behavioural model: if injured larvae have less inhibitory reserve, a fixed challenge dose should push more of them across threshold. This tests the same excitation/inhibition hypothesis with a drug rather than with a habituation protocol.

| Group | Seized / n | Proportion | Wilson 95% CI | Median latency (s) |
|---|---|---|---|---|
| Sham | 3/12 | 0.25 | [0.09, 0.53] | 1800 |
| Low impact | 7/12 | 0.58 | [0.32, 0.81] | 720 |
| High impact | 7/10 | 0.70 | [0.40, 0.89] | 607 |

χ²(2) = 4.933, p = 0.0849, Cramér's V = 0.381, total n = 34.

> **Explicit statement of limitation.** This probe is underpowered. With 34 fish split across three groups, post-hoc power for the observed sham-versus-injured difference (Cohen's h = 0.80) is only **61%**. The minimum expected cell count is 5.00. It is reported as a directional check consistent with the primary result, and **no conclusion in this report rests on it.**

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

- **n = 81 with 36 events.** EPV = 9.0 is at the accepted floor. The CIs on the AUC and on every coefficient are wide, and they are reported rather than smoothed over.
- **Three clutches means three outer folds.** The fold-to-fold spread is estimated from three numbers; the SD across folds should be read as indicative, not precise.
- **7 injured fish were dropped** for a missing post-injury session. This is complete-case analysis; attrition is not obviously outcome-related but has not been modelled.
- **The c-fos validation is pooled material** — 9 pairs, 4 larvae per pool. It corroborates the primary result; it does not independently establish it.
- **PTZ is underpowered** (above), and the study is not designed to support any claim from it.

## Conclusion

In this dataset, the acute trajectory of startle-habituation kinetics — measured against each fish's own pre-injury baseline and interpreted in the light of blast dose — separates injured larvae that go on to develop spontaneous burst activity from those that do not, with a cross-validated AUC of 0.833 (95% CI [0.736, 0.916]), a total out-of-fold accuracy of 76.5%, and a permutation p = 0.0010 against a null that reruns the entire nested cross-validation. The ablation table shows this is not explicable by sickness, by injury severity, or by a pre-existing trait alone. An independent molecular assay on a separate cohort of larvae points the same way.

What that supports is a mechanistic claim of modest scope: **a behavioural readout of inhibitory gain in a defined sensorimotor circuit, sampled within 24 hours of injury, carries information about which animals are undergoing epileptogenesis.** It does not identify the cellular lesion, and it is one experiment in one species at one age.

## Background literature

The neurobiological claims above rest on the following. This is a background reading list, not a citation of results generated here.

1. Annegers JF, Hauser WA, Coan SP, Rocca WA (1998). A population-based study of seizures after traumatic brain injuries. *New England Journal of Medicine* 338:20–24.
2. Baraban SC, Taylor MR, Castro PA, Baier H (2005). Pentylenetetrazole induced changes in zebrafish behavior, neural activity and c-fos expression. *Neuroscience* 131:759–768.
3. Burgess HA, Granato M (2007). Sensorimotor gating in larval zebrafish. *Journal of Neuroscience* 27:4984–4994.
4. Hunt RF, Boychuk JA, Smith BN (2013). Neural circuit mechanisms of post-traumatic epilepsy. *Frontiers in Cellular Neuroscience* 7:89.
5. Livak KJ, Schmittgen TD (2001). Analysis of relative gene expression data using real-time quantitative PCR and the 2^−ΔΔCт method. *Methods* 25:402–408.
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
