# Results: larval zebrafish blast TBI and post-traumatic epileptogenesis

Generated 2026-08-10 02:37 from `tbidataset.xlsx`  
**Random seed: `20260809`** (numpy, scikit-learn, permutation and bootstrap draws). Rerunning `python run_all.py` reproduces every number below.

Every statistic quoted here is also in [`results/all_statistics.csv`](results/all_statistics.csv), with its test, statistic, degrees of freedom, p-value, effect size and confidence interval.

---

## Summary of findings

| # | Finding | Test | Effect size | p |
|---|---------|------|-------------|---|
| 1 | The nonlinear refit reproduces the supplied `decay_constant` | Pearson correlation, 390 sessions | r = 1.0000 | < 0.0001 |
| 2 | Groups are indistinguishable at the pre-injury baseline | one-way ANOVA on baseline τ | η² = 0.0316 | = 0.1319 |
| 3 | Low and high dose move τ in **opposite** directions after blast | Welch t, Δτ at 0.5 h | d = 1.97 | < 0.0001 |
| 4 | Conversion is predictable from pre-injury τ and dose | nested CV with model selection, 1000 permutations | AUC = 0.787 | = 0.0010 |
| 5 | The acute Δτ terms do **not** add predictive value here | ablation, same CV | ΔAUC = -0.051 | n/a |
| 6 | Converter pools carry more c-fos than matched non-converter pools | paired t, 6 matched cells | dz = 1.52 | = 0.0136 |
| 7 | c-fos does not rise with injury alone (specificity control) | Welch t, non-converter pools | d = 0.39 | = 0.6072 |
| 8 | Injured larvae seize more readily under PTZ | χ², n = 54 | V = 0.515 | = 0.0008 |
| 9 | PTZ seizure and conversion coincide in the same larvae | Fisher exact | OR = 17.0 | = 0.0019 |

> **Headline.** Conversion is predictable, but **not from the variable this study was designed around.** The pre-injury habituation constant plus blast dose carries the signal; the acute change in habituation does not add to it in this dataset. Section 2.2 sets out what that does and does not license.

---

## Step 1. Rebuilding the outcome variable

Per fish per session, `distance_mm(k) = A·exp(-(k-1)/τ) + C` was fitted to all 30 trials by nonlinear least squares (`scipy.optimize.curve_fit`, bounded, four starting points per session, best sum of squared errors retained).

- **390 fish-sessions fitted**, convergence 100.0%, median R² = 0.836 (5th percentile 0.659); 0 sessions hit a parameter bound.
- **Agreement with the supplied `fish_features.decay_constant`:** Pearson r = 1.0000 (p < 0.0001), Spearman ρ = 1.0000. Mean bias (refit minus supplied) = +0.000 trials, 95% CI [+0.000, +0.000]; 95% limits of agreement [-0.001, +0.001].
  The refit is what every downstream analysis uses; the supplied column is a check, not an input.

### Why log-linearisation was not used

Subtracting an estimated offset and regressing `log(y - C)` on trial number requires y > C. Once responses reach the habituated floor they scatter below it: **26.2% of trials are unusable**. Discarding them tilts the fitted slope so far that correlation with the reference τ collapses from r = 1.000 (nonlinear) to r = -0.252 (log-linear), i.e. it inverts.
This diagnostic runs on every execution rather than being asserted.

Figures: `fig01_curvefit_examples.png`, `fig02_tau_agreement.png`.

---

## Step 2. Prediction model

**Analysis set.** Injured larvae only (sham dropped), one row per fish, complete on all four predictors and the outcome: **n = 82, 34 converters (41.5%)**. 6 injured larvae were excluded for a missing 0.5 h or 24 h session.

> **Events per variable = 8.5, below the conventional floor of 9 to 10.** The predictor set is held at four by design rather than trimmed to fit the rule, because each term answers a distinct pre-specified question. The cost is wider intervals on every coefficient, and they are reported rather than hidden. Section 2.3 shows the model that the data actually support, which uses two predictors and comfortably clears the rule.

**Predictors.** `dose` (high_impact = 1); `pre_tau` (τ at t = -1 h); `z_dtau_0.5` and `z_dtau_24` (τ@0.5 - τ@-1 and τ@24 - τ@-1, z-scored **within dose group**).

**Model.** L2-penalised logistic regression in a `StandardScaler` pipeline. No ensembles: at n = 82 with 34 events, a random forest or boosted model has enough capacity to memorise the sample, and its output cannot be sign-checked against the biology.

**Validation.** `GroupKFold` on clutch for the outer split (leave-one-clutch-out), never random. The penalty strength C is tuned by **nested** cross-validation inside each outer training set, never on the folds reported. Within-dose z-scoring is a pipeline step fitted on training folds only, so no test-fold information reaches the scaling constants.

### 2.1 Ablation table

Same model class, same cross-validation, different inputs.

| Model | Question | Mean fold AUC | SD | Fold range | Pooled OOF AUC [95% CI] | Brier |
|---|---|---|---|---|---|---|
| (a) baseline_locomotion only | Is it just sickness? | **0.392** | 0.066 | 0.322-0.453 | 0.418 [0.292, 0.541] | 0.249 |
| (b) dose only | Is it just injury severity? | **0.623** | 0.044 | 0.592-0.674 | 0.604 [0.477, 0.728] | 0.245 |
| (c) dose + pre_tau | Is it a pre-existing trait? | **0.788** | 0.061 | 0.722-0.841 | 0.773 [0.664, 0.869] | 0.239 |
| (d) dose + z_dtau | Is it the injury response? | **0.627** | 0.097 | 0.547-0.735 | 0.620 [0.491, 0.747] | 0.238 |
| (e) full 4-predictor model | Full model | **0.738** | 0.083 | 0.667-0.829 | 0.756 [0.647, 0.856] | 0.198 |

Per-fold AUC (held-out clutch), with the C selected inside each fold:

| Model | clutch_A | clutch_B | clutch_C |
|---|---|---|---|
| (a) baseline_locomotion only | 0.322 (C=0.01) | 0.453 (C=0.01) | 0.400 (C=0.01) |
| (b) dose only | 0.592 (C=0.01) | 0.604 (C=0.01) | 0.674 (C=0.01) |
| (c) dose + pre_tau | 0.722 (C=0.01) | 0.802 (C=0.01) | 0.841 (C=0.01) |
| (d) dose + z_dtau | 0.600 (C=0.3) | 0.547 (C=10) | 0.735 (C=3) |
| (e) full 4-predictor model | 0.717 (C=1) | 0.667 (C=0.3) | 0.829 (C=100) |

Reading it:

- **(a) baseline locomotion alone: 0.392.** At or below chance out of fold. Conversion is not a readout of how sick or sluggish the larva is.
- **(b) dose alone: 0.623.** Above chance. Blast dose is a genuine main effect in this cohort, which was not true in the pilot data.
- **(c) dose + pre_tau: 0.788.** The strongest row in the table, and the most stable across folds (SD 0.061).
- **(d) dose + z_dtau: 0.627.** Barely above dose alone. The acute injury response adds little.
- **(e) the pre-specified four-predictor model: 0.738.** Adding the two Δτ terms to (c) **costs -0.051 AUC**.

### 2.2 The result this study did not expect

The design was built on the hypothesis that the *acute change* in habituation kinetics carries the predictive signal. In this dataset it does not. The two Δτ terms:

- do not separate converters from non-converters within either dose group (low impact p = 0.4748, high impact p = 0.6306);
- add nothing to dose on their own (model (d), 0.627, versus dose alone, 0.623);
- and actively degrade the model when added to `pre_tau` (model (e), 0.738, versus model (c), 0.788).

What does carry signal is `pre_tau`, the **pre-injury** habituation constant, together with dose. That is a different scientific claim, and a weaker one for a biomarker: a variable measured before the injury cannot be a readout of the injury response. It points to **susceptibility** rather than to acute pathophysiology. Section 4.2 discusses what that means mechanistically, and the limitations section says plainly why the original framing is not supported here.

The group-level Δτ effect is not in doubt: low and high dose still move τ in opposite directions with a very large effect (Cohen's d = 1.97, p < 0.0001). Injury clearly perturbs the circuit. What fails is the step from that group difference to **individual** prediction: within a dose group, the size of a larva's acute Δτ does not tell you whether that larva will convert.

### 2.3 Honest model selection

Quoting the best row of an ablation table is selection on the test folds: the winner is partly chosen by the noise in the folds it is scored on. To get an unbiased number, the choice of predictor set is treated as one more hyper-parameter and tuned on the **inner** folds only.

| Held-out clutch | Selected by inner CV | C | Outer-fold AUC |
|---|---|---|---|
| clutch_B | `c_dose_pretau` | 0.01 | 0.802 |
| clutch_C | `c_dose_pretau` | 0.01 | 0.841 |
| clutch_A | `e_full` | 1 | 0.717 |

- **Honestly selected performance: mean fold AUC 0.787** (SD 0.064), pooled out-of-fold AUC 0.746 [0.629, 0.848].
- Naively quoting the best table row gives 0.788. The optimism from picking it by eye is +0.002 AUC, which is small here only because the same model wins in most folds.
- **Total out-of-fold accuracy 70.7%** at threshold 0.50 (58 of 82 correct), balanced accuracy 69.0%, sensitivity 58.8%, specificity 79.2%.
- Permutation test on the **entire selection procedure** (selection re-run inside every one of the 1000 shuffles): null mean 0.487, p = 0.0010.

This is the number to quote for best achievable accuracy on this design, because it is the performance of a procedure that could be applied to a new clutch without knowing the answer first.

### 2.4 Pre-specified model: full detail

Reported because it was pre-specified, not because it is the best.

- Mean fold AUC 0.738 (SD 0.083), pooled out-of-fold AUC 0.756 [0.647, 0.856].
- Permutation p = 0.0010 against a null that reruns the whole nested CV 1000 times (null mean 0.493, SD 0.079).
- Total out-of-fold accuracy 70.7% at threshold 0.50, balanced accuracy 69.0%, sensitivity 58.8%, specificity 79.2%, Brier 0.198.

Confusion matrix, out-of-fold, threshold 0.50:

| | Predicted non-converter | Predicted converter |
|---|---|---|
| **Non-converter** | 38 | 10 |
| **Converter** | 14 | 20 |

At the Youden-optimal threshold (0.357) accuracy is 69.5% with sensitivity 82.4%. That threshold was chosen on these same out-of-fold predictions, so it is mildly optimistic and is shown only to indicate the achievable operating range.

**Coefficients** (refit on all data at C = 100, standardised scale, clutch-clustered bootstrap CIs, 2000 resamples):

| Predictor | β (per SD) | 95% CI | OR per SD | OR 95% CI |
|---|---|---|---|---|
| `dose` | +0.936 \* | [+0.405, +1.851] | 2.55 | [1.50, 6.37] |
| `pre_tau` | +1.300 \* | [+0.708, +2.594] | 3.67 | [2.03, 13.39] |
| `dtau_0.5` | +0.100 | [-0.484, +0.637] | 1.11 | [0.62, 1.89] |
| `dtau_24` | +0.733 \* | [+0.021, +1.850] | 2.08 | [1.02, 6.36] |

Intercept -0.448. \* = bootstrap CI excludes zero.

#### What the coefficients mean neurobiologically

τ is the number of trials the acoustic startle response takes to decay to its floor. Startle habituation in larval zebrafish is not fatigue of the Mauthner cell, the hindbrain command neuron for the C-start escape. It is produced by progressive **feedforward inhibition onto the M-cell lateral dendrite**, which reduces dendritic excitability with repeated stimulation (Marsden & Granato, 2015). A larger τ therefore means inhibition accumulates more slowly, i.e. **reduced inhibitory gain in a defined sensorimotor circuit**. That is the same excitation/inhibition set point whose collapse drives epileptogenesis after brain injury.

- **`pre_tau` (β = +1.300, OR 3.67 per SD), the dominant term.** Larvae that habituate more slowly *before* any injury are more likely to convert. Baseline inhibitory tone varies between individuals, and an animal that starts nearer the seizure threshold has less reserve to lose. This is a susceptibility term, the animal analogue of pre-injury risk factors that modify post-traumatic epilepsy risk in humans, and it fits a two-hit framing of predisposition plus insult.
- **`dose` (β = +0.936).** A real main effect here: conversion rises monotonically with blast dose.
- **The Δτ terms.** Both intervals are wide and the 0.5 h term straddles zero. Given that the ablation shows they cost AUC, the honest reading is that they carry no individual-level information in this cohort, not that they carry a small amount.

### 2.5 Leakage quantification

| Split | Mean fold AUC | Pooled OOF AUC |
|---|---|---|
| 5-fold **random** (leaky, for comparison only) | 0.828 | 0.744 |
| Leave-one-clutch-out (**the honest estimate**) | 0.738 | 0.756 |

Ignoring clutch inflates the mean fold AUC by **+0.091**. Clutches are sibling groups run on separate days; a random split puts siblings on both sides of the partition and the model learns clutch identity. Any AUC from a random split on this design should be discounted by roughly this much.

Figures: `fig03_roc_nested_comparison.png`, `fig04_permutation_null.png`, `fig05_confusion_calibration.png`, `fig06_coefficients.png`, `fig07_fold_auc_by_model.png`.

---

## Step 3. Molecular validation: paired c-fos pools

### 3.1 What this validates, and what it does not

In this dataset the qPCR pools are labelled by **realised outcome** (`converter` / `non_converter`), not by a predicted risk score, and they are drawn from the same `zf_*` larvae that train the model. Two consequences, both stated rather than glossed:

1. This is **not** an independent test of the prediction model. It is a test of the **outcome variable the model predicts**: does a larva scored as converted carry a molecular signature of elevated network activity, or is `converted` merely a behavioural scoring threshold?
2. Model and assay share animals, so sample independence is not claimed.

What makes it worth doing is that the modality is genuinely different. Everything in Step 2 derives from how far a larva swam. If that pipeline carried a systematic artefact, no amount of internal cross-validation would reveal it, because every fold would inherit the artefact. Transcript abundance shares no instrumentation with swim tracking.

`fosab` is the zebrafish orthologue of *c-fos*, the canonical immediate early gene. Sustained depolarisation raises intracellular Ca²⁺, driving CaMK- and MAPK/ERK-dependent phosphorylation of CREB and transcription from the *fos* promoter within roughly 15 to 30 minutes (Sheng & Greenberg, 1990). c-fos transcript is a molecular integrator of recent network activity and the standard readout for mapping seizure-recruited circuits, including in the original characterisation of chemically induced seizures in larval zebrafish (Baraban et al., 2005). Normalisation is against `rpl13a`, a validated zebrafish reference gene (Tang et al., 2007), by the 2^-ΔΔCт method (Livak & Schmittgen, 2001).

### 3.2 Why six pairs and not nine

There are 27 pools, three larvae each. They are not independent units: they are matched within (group × clutch) cells. The counts are unbalanced by biology, because conversion is rare in sham:

| Group | Converter pools per clutch | Non-converter pools per clutch |
|---|---|---|
| sham | 0 | 2 |
| low_impact | 1 | 2 |
| high_impact | 2 | 2 |

A balanced nine-pair design is therefore impossible: **sham cells contain no converter pool at all**, because too few sham larvae converted to fill one. The pairing the data support is **6 matched cells**, the injured group × clutch combinations, with replicate pools within a cell averaged before pairing. Sham is used instead as an unpaired reference level, which answers a complementary question in section 3.4.

### 3.3 Paired result

| Contrast | Pairs | Mean Δ (converter - non-converter) | 95% CI | Paired t | p | Cohen's dz | Wilcoxon p |
|---|---|---|---|---|---|---|---|
| Fold-change scale (primary) | 6 | +0.9549 | [+0.2962, +1.6136] | t(5) = +3.727 | = 0.0136 | +1.521 | = 0.0312 |
| log₂ scale | 6 | +0.7951 | [+0.3070, +1.2831] | t(5) = +4.188 | = 0.0086 | +1.710 | = 0.0312 |

Converter pools carry **73.5% more c-fos transcript** relative to `rpl13a` than matched non-converter pools (geometric mean ratio 1.735, 95% CI [1.237, 2.434], back-transformed from the paired log₂ analysis, which is the appropriate scale for a fold change).

Direction is consistent in **6 of 6** pairs (exact binomial sign test, p = 0.0312).

### 3.4 Specificity control: does c-fos track conversion, or just injury?

Since sham cannot be paired, the control question is turned around. If injury alone raised c-fos, then **non-converter** pools from injured groups should sit above non-converter pools from sham. They do not:

- Injured non-converter pools 1.256 versus sham non-converter pools 1.181; difference +0.074, 95% CI [-0.242, +0.391], Welch t = +0.556, p = 0.6072, d = 0.39.

**A null result here is the desired one**, and it is what rules out the trivial explanation. The c-fos elevation tracks *which larvae converted*, not *which larvae were hit*. That is the specific claim, and it is the one that makes the assay informative about epileptogenesis rather than about injury exposure.

Elevated immediate early gene expression in the absence of any provoking stimulus is what a chronically over-active network looks like transcriptionally. This is corroboration, not proof: it is a bulk measurement on pooled tissue, so it cannot localise the signal to a cell type or region, and it cannot distinguish loss of parvalbumin-positive interneuron function from increased glutamatergic drive. Both are documented consequences of brain injury.

No regression of c-fos on a pooled continuous risk score is run. The pools are outcome-labelled, and the risk score is not comparable between dose groups.

Figure: `fig08_cfos_paired.png`.

---

## Step 4. Descriptive results

### 4.1 Groups start equal

At the pre-injury baseline, τ does not differ between groups: F(2,126) = 2.059, p = 0.1319, η² = 0.0316; Kruskal-Wallis H = 3.749, p = 0.1534. Baseline locomotion likewise (p = 0.0866). Randomisation held, and every larva is its own control thereafter.

### 4.2 τ moves in opposite directions by dose

| Group | t = -1 h | t = 0.5 h | t = 24 h |
|---|---|---|---|
| Sham | 6.74 ± 0.29 | 6.74 ± 0.47 | 6.48 ± 0.41 |
| Low impact | 6.44 ± 0.31 | 9.80 ± 0.46 | 7.36 ± 0.47 |
| High impact | 5.94 ± 0.23 | 3.96 ± 0.33 | 5.29 ± 0.48 |

(mean ± SEM, trials to habituate)

At 0.5 h, Δτ from each larva's own baseline is **+3.29** trials in low impact and **-1.94** in high impact. Welch t = 9.01, p < 0.0001, Cohen's d = 1.97.

The divergence is two different lesions on the same circuit, not a nuisance to correct away:

- **Low dose, τ rises (habituation deficit).** Sublethal blast preferentially compromises the feedforward inhibition that normally accumulates onto the Mauthner cell across trials. Inhibition builds more slowly, the escape response persists, τ lengthens. This is disinhibition, and it maps onto the loss of GABAergic control reported after experimental brain injury.
- **High dose, τ falls (fatigue, not learning).** A shorter τ looks superficially like better habituation. It is not. Greater energy deposition depresses the excitatory limb as well: the acute metabolic crisis and depolarisation that follow severe injury reduce the startle response itself, so the fitted decay is fast because the response never had far to fall.

Both are pathological and they move the same scalar in opposite directions, which is why any model using Δτ must also encode dose. **This group-level effect is robust. What this dataset shows is that it does not translate into individual prediction** (section 2.2).

### 4.3 Operational metrics

| Metric | Value |
|---|---|
| Sessions | 12 (9 habituation, 3 outcome) |
| Fish-sessions recorded | 522 |
| Fish per hour | 29.1 ± 6.2 |
| Operator minutes per fish | 0.99 ± 0.21 |
| Consumables cost per fish | $0.093 |
| Total consumables | $48.70 |
| Total operator time | 8.6 h |
| Attrition | 4/129 = 3.1% [1.2%, 7.7%] |

Figures: `fig09_habituation_curves.png`, `fig10_tau_by_timepoint.png`, `fig11_converter_trajectories.png`, `fig12_operations.png`.

---

## Step 5. PTZ seizure threshold

Pentylenetetrazol is a non-competitive GABAₐ receptor antagonist: it binds at the picrotoxin site, reduces chloride conductance and removes inhibitory brake from the network. In larval zebrafish it produces stereotyped, dose-dependent seizure behaviour with electrographic correlates (Baraban et al., 2005), which makes it the standard pharmacological probe of **seizure threshold**. If injured larvae have less inhibitory reserve, a fixed challenge dose should push more of them across it.

| Group | Seized / n | Proportion | Wilson 95% CI | Median latency (s) |
|---|---|---|---|---|
| Sham | 5/17 | 0.29 | [0.13, 0.53] | 1800 |
| Low impact | 11/15 | 0.73 | [0.48, 0.89] | 902 |
| High impact | 19/22 | 0.86 | [0.67, 0.95] | 589 |

χ²(2) = 14.300, p = 0.0008, Cramér's V = 0.515, total n = 54.

**On power.** At this n the sham-versus-injured contrast is adequately powered (96%, Cohen's h = 1.10), so unlike the pilot cohort this probe is not merely directional. It remains a **secondary** outcome for a different reason: it is a group-level comparison, and the primary claim of this project is about predicting individual animals, which a group difference does not address. The low-versus-high dose contrast is still not powered, and the two injured groups are not distinguishable here (p = 0.4081).

### 5.1 Is the challenge a confound?

PTZ is a proconvulsant given to a subset of the **same** larvae that supply the conversion outcome, so the obvious worry is that the drug caused the outcome. It did not: conversion rates are effectively identical between challenged and unchallenged larvae within every group.

| Group | Challenged | Not challenged | Fisher p |
|---|---|---|---|
| sham | 1/17 | 1/27 | = 1.0000 |
| low_impact | 5/15 | 8/29 | = 0.7367 |
| high_impact | 12/22 | 12/22 | = 1.0000 |

### 5.2 Seizure threshold and conversion coincide

Because PTZ and conversion are measured on the same larvae, they can be crossed directly. This was not pre-specified and is exploratory.

- Larvae that seized under PTZ were far more likely to have converted: Fisher exact odds ratio **17.0**, p = 0.0019, n = 54.
- Converters reached their first seizure sooner: median converter 405 s, non-converter 1775 s, Mann-Whitney U = 144, p = 0.0008.

Two assays that share no measurement apparatus, a pharmacological threshold test and a spontaneous-activity recording, agree at the level of the individual animal. That is the strongest evidence in this report that `converted` denotes a real hyperexcitable state rather than a scoring artefact. It says nothing about whether that state is *predictable* in advance, which is Step 2's job.

Figure: `fig13_ptz.png`.

---

## Assumptions, checks and limitations

Assumption checks print during every run tagged `[PASS]` or `[FLAG]`, and the underlying statistics are in `all_statistics.csv`.

- Curve fitting: 100.0% convergence, median R² 0.836, 0 sessions at a parameter bound.
- Collinearity: max |r| among predictors = 0.706; all VIFs are in `all_statistics.csv` under `assumptions`.
- Separation: the unpenalised Logit converged with max |β| = 6.22, so the ridge penalty is not masking complete separation.
- **Flagged:** baseline τ departs from normality in low impact, high impact (Shapiro-Wilk p ≤ 0.05). The ANOVA is therefore backed by Kruskal-Wallis (p = 0.1534), which agrees.
- Paired c-fos differences: Shapiro-Wilk reported for each contrast, Wilcoxon alongside every paired t-test.

Real limitations, stated plainly:

- **The study's central hypothesis is not supported.** The acute Δτ terms, which the design was built around, add nothing to individual prediction and slightly degrade it. What predicts conversion is a pre-injury trait plus dose. Reporting this is the point of running the ablation rather than fitting one model and describing it.
- **n = 82 with 34 events, EPV = 8.5**, below the conventional floor. Coefficient intervals are wide. The two-predictor model that the data support clears the rule comfortably (EPV = 17.0).
- **Three clutches means three outer folds.** Fold-to-fold spread is estimated from three numbers and should be read as indicative.
- **6 injured larvae were dropped** for a missing session. This is a complete-case analysis; attrition has not been modelled.
- **Conversion is a behavioural proxy.** Spontaneous burst activity is not electrographically confirmed epilepsy, though the PTZ concordance in section 5.2 supports it.
- **The c-fos and PTZ analyses share animals with the model cohort.** They validate the outcome label, not the model's generalisation.
- **`pre_tau` as a predictor cannot support the biomarker framing.** A variable measured before injury is a susceptibility marker, not an acute readout, and it could not be used to triage patients after an injury that has already happened.
- **One species, one age, one injury model.** Generalisation to mammalian TBI is a hypothesis.

## Conclusion

Conversion to spontaneous burst activity is predictable in this cohort at a cross-validated AUC of 0.787 and total out-of-fold accuracy of 70.7%, with the predictor set chosen inside the cross-validation and a permutation p = 0.0010 against a null that reruns the entire procedure. Baseline locomotion is ruled out as an explanation, and clutch-aware splitting shows a random split would have inflated the estimate by +0.091 AUC.

The signal lies in the **pre-injury** habituation constant together with blast dose, not in the acute change in habituation that this study set out to test. Injury does perturb the circuit strongly and in dose-dependent opposite directions, but that group-level effect does not carry individual-level predictive information here.

Two assays that share no instrumentation with the behavioural pipeline agree that `converted` denotes a genuinely hyperexcitable state: converter pools carry 73.5% more c-fos transcript than matched non-converter pools, with injured non-converters indistinguishable from sham, and larvae that seize under PTZ are 17 times more likely to have converted. The outcome variable is real. Whether it can be predicted from the acute injury response remains open, and on this evidence the answer is no.

## Background literature

Reading behind the neurobiological claims above. These support the framing; they are not citations of results generated here.

1. Annegers JF, Hauser WA, Coan SP, Rocca WA (1998). A population-based study of seizures after traumatic brain injuries. *New England Journal of Medicine* 338:20-24.
2. Baraban SC, Taylor MR, Castro PA, Baier H (2005). Pentylenetetrazole induced changes in zebrafish behavior, neural activity and c-fos expression. *Neuroscience* 131:759-768.
3. Burgess HA, Granato M (2007). Sensorimotor gating in larval zebrafish. *Journal of Neuroscience* 27:4984-4994.
4. Hunt RF, Boychuk JA, Smith BN (2013). Neural circuit mechanisms of post-traumatic epilepsy. *Frontiers in Cellular Neuroscience* 7:89.
5. Livak KJ, Schmittgen TD (2001). Analysis of relative gene expression data using real-time quantitative PCR and the 2^-ΔΔCт method. *Methods* 25:402-408.
6. Marsden KC, Granato M (2015). In vivo Ca²⁺ imaging reveals that decreased dendritic excitability drives startle habituation. *Cell Reports* 13:1733-1740.
7. Peduzzi P, Concato J, Kemper E, Holford TR, Feinstein AR (1996). A simulation study of the number of events per variable in logistic regression analysis. *Journal of Clinical Epidemiology* 49:1373-1379.
8. Sheng M, Greenberg ME (1990). The regulation and function of c-fos and other immediate early genes in the nervous system. *Neuron* 4:477-485.
9. Sloviter RS (1991). Permanently altered hippocampal structure, excitability, and inhibition after experimental status epilepticus in the rat: the 'dormant basket cell' hypothesis. *Hippocampus* 1:41-66.
10. Tang R, Dodd A, Lai D, McNabb WC, Love DR (2007). Validation of zebrafish (*Danio rerio*) reference genes for quantitative real-time RT-PCR normalization. *Acta Biochimica et Biophysica Sinica* 39:384-390.
11. Varma S, Simon R (2006). Bias in error estimation when using cross-validation for model selection. *BMC Bioinformatics* 7:91.
12. Wolman MA, Jain RA, Liss L, Granato M (2011). Chemical modulation of memory formation in larval zebrafish. *PNAS* 108:15468-15473.

## Reproducing

```bash
pip install -r requirements.txt
python run_all.py
```

Seed `20260809`. Outputs: `results/figures/*.png` at 300 dpi, `results/all_statistics.csv`, `results/tables/*.csv`, and this file.
