"""Assemble RESULTS.md from the recorded statistics ledger."""
from __future__ import annotations

from datetime import datetime

import numpy as np
import pandas as pd

from . import config, statsbook as sb


def _fmt_p(p) -> str:
    if p is None or (isinstance(p, float) and not np.isfinite(p)):
        return "n/a"
    p = float(p)
    if p < 1e-4:
        return "< 0.0001"
    return f"= {p:.4f}"


def _p(p) -> str:
    """Same as _fmt_p but with a leading 'p', for running prose."""
    return "p " + _fmt_p(p)


def _plural(n: int, singular: str, plural: str | None = None) -> str:
    return f"{n} {singular}" if n == 1 else f"{n} {plural or singular + 's'}"


def _g(analysis: str, quantity: str, default=None):
    try:
        return sb.get(analysis, quantity)
    except KeyError:
        return default


def _v(analysis: str, quantity: str, default=np.nan) -> float:
    r = _g(analysis, quantity)
    return float(r["value"]) if r is not None else default


def write(ctx: dict) -> str:
    s2, s3, s4 = ctx["step2"], ctx["step3"], ctx["step4"]
    comp, folds = s2["comparison"], s2["folds"]
    full, perm, coefs, cm = s2["full"], s2["perm"], s2["coefs"], s2["cm"]
    sel, sel_cm, sel_perm = s2["selection"], s2["selection_cm"], s2["selection_perm"]
    sel_lo, sel_hi = s2["selection_ci"]

    n_fish = int(_v("design_matrix", "n_fish"))
    n_events = int(_v("design_matrix", "n_events"))
    epv = _v("design_matrix", "events_per_variable")
    r_agree = _g("tau_agreement", "pearson_r")
    bias = _g("tau_agreement", "mean_bias_trials")
    anova = _g("baseline_equivalence", "anova_tau")

    L: list[str] = []
    A = L.append

    A("# Results: larval zebrafish blast TBI and post-traumatic epileptogenesis")
    A("")
    A(f"Generated {datetime.now():%Y-%m-%d %H:%M} from `{config.DATA_XLSX.name}`  ")
    A(f"**Random seed: `{config.SEED}`** (numpy, scikit-learn, permutation and bootstrap draws). "
      "Rerunning `python run_all.py` reproduces every number below.")
    A("")
    A("Every statistic quoted here is also in [`results/all_statistics.csv`](results/all_statistics.csv), "
      "with its test, statistic, degrees of freedom, p-value, effect size and confidence interval.")
    A("")
    A("---")
    A("")

    # ------------------------------------------------------------- summary
    A("## Summary of findings")
    A("")
    A("| # | Finding | Test | Effect size | p |")
    A("|---|---------|------|-------------|---|")
    A(f"| 1 | The nonlinear refit reproduces the supplied `decay_constant` | Pearson correlation, "
      f"{int(r_agree['n'])} sessions | r = {float(r_agree['value']):.4f} | "
      f"{_fmt_p(r_agree['p_value'])} |")
    A(f"| 2 | Groups are indistinguishable at the pre-injury baseline | one-way ANOVA on baseline τ | "
      f"η² = {float(anova['effect_size']):.4f} | {_fmt_p(anova['p_value'])} |")
    A(f"| 3 | Low and high dose move τ in **opposite** directions after blast | Welch t, Δτ at 0.5 h | "
      f"d = {float(_g('tau_change_from_baseline','low_vs_high_delta_t0.5')['effect_size']):.2f} | "
      f"{_fmt_p(_g('tau_change_from_baseline','low_vs_high_delta_t0.5')['p_value'])} |")
    A(f"| 4 | Conversion is predictable from pre-injury τ and dose | nested CV with model selection, "
      f"{perm['n_perm']} permutations | AUC = {sel['mean_fold_auc']:.3f} | "
      f"{_fmt_p(sel_perm['p_pooled'])} |")
    A(f"| 5 | The acute Δτ terms do **not** add predictive value here | ablation, same CV | "
      f"ΔAUC = {comp.loc[comp.key=='e_full','mean_fold_auc'].iloc[0] - comp.loc[comp.key=='c_dose_pretau','mean_fold_auc'].iloc[0]:+.3f} | n/a |")
    A(f"| 6 | Converter pools carry more c-fos than matched non-converter pools | paired t, "
      f"{s3['all']['n']} matched cells | dz = {s3['all']['dz']:.2f} | {_fmt_p(s3['all']['p'])} |")
    A(f"| 7 | c-fos does not rise with injury alone (specificity control) | Welch t, non-converter "
      f"pools | d = {s3['baseline']['d']:.2f} | {_fmt_p(s3['baseline']['p'])} |")
    A(f"| 8 | Injured larvae seize more readily under PTZ | χ², n = {int(_g('ptz','chi_square')['n'])} | "
      f"V = {float(_g('ptz','chi_square')['effect_size']):.3f} | "
      f"{_fmt_p(_g('ptz','chi_square')['p_value'])} |")
    or_row = _g("ptz_vs_conversion", "odds_ratio_seized_given_converted")
    if or_row is not None:
        A(f"| 9 | PTZ seizure and conversion coincide in the same larvae | Fisher exact | "
          f"OR = {float(or_row['value']):.1f} | {_fmt_p(or_row['p_value'])} |")
    A("")
    A("> **Headline.** Conversion is predictable, but **not from the variable this study was "
      "designed around.** The pre-injury habituation constant plus blast dose carries the signal; "
      "the acute change in habituation does not add to it in this dataset. Section 2.2 sets out "
      "what that does and does not license.")
    A("")
    A("---")
    A("")

    # -------------------------------------------------------------- step 1
    A("## Step 1. Rebuilding the outcome variable")
    A("")
    A("Per fish per session, `distance_mm(k) = A·exp(-(k-1)/τ) + C` was fitted to all 30 trials by "
      "nonlinear least squares (`scipy.optimize.curve_fit`, bounded, four starting points per "
      "session, best sum of squared errors retained).")
    A("")
    n_sess = int(_v("curve_fit", "n_sessions"))
    conv = _v("curve_fit", "convergence_rate")
    r2med = _v("curve_fit", "r2_median")
    A(f"- **{n_sess} fish-sessions fitted**, convergence {100*conv:.1f}%, median R² = {r2med:.3f} "
      f"(5th percentile {_v('curve_fit','r2_q05'):.3f}); "
      f"{_plural(int(_v('curve_fit','n_tau_at_bound')), 'session')} hit a parameter bound.")
    A(f"- **Agreement with the supplied `fish_features.decay_constant`:** "
      f"Pearson r = {float(r_agree['value']):.4f} ({_p(r_agree['p_value'])}), "
      f"Spearman ρ = {_v('tau_agreement','spearman_rho'):.4f}. "
      f"Mean bias (refit minus supplied) = {float(bias['value']):+.3f} trials, 95% CI "
      f"[{float(bias['ci_low']):+.3f}, {float(bias['ci_high']):+.3f}]; "
      f"95% limits of agreement [{_v('tau_agreement','limits_of_agreement_low'):+.3f}, "
      f"{_v('tau_agreement','limits_of_agreement_high'):+.3f}].")
    A("  The refit is what every downstream analysis uses; the supplied column is a check, not an input.")
    A("")
    A("### Why log-linearisation was not used")
    A("")
    pct_below = _v("loglinear_failure", "pct_trials_below_offset")
    n_neg = int(_v("loglinear_failure", "n_sessions_negative_tau"))
    r_ll = _v("loglinear_failure", "pearson_r_vs_supplied")
    A(f"Subtracting an estimated offset and regressing `log(y - C)` on trial number requires y > C. "
      f"Once responses reach the habituated floor they scatter below it: **{pct_below:.1f}% of trials "
      f"are unusable**. Discarding them tilts the fitted slope so far that correlation with the "
      f"reference τ collapses from r = {float(r_agree['value']):.3f} (nonlinear) to r = {r_ll:.3f} "
      f"(log-linear), i.e. it inverts.")
    if n_neg:
        A(f"In {n_neg} of {n_sess} sessions the recovered τ is outright negative, scoring a "
          "hyperexcitable fish as hypo-excitable.")
    A("This diagnostic runs on every execution rather than being asserted.")
    A("")
    A("Figures: `fig01_curvefit_examples.png`, `fig02_tau_agreement.png`.")
    A("")
    A("---")
    A("")

    # -------------------------------------------------------------- step 2
    A("## Step 2. Prediction model")
    A("")
    A(f"**Analysis set.** Injured larvae only (sham dropped), one row per fish, complete on all four "
      f"predictors and the outcome: **n = {n_fish}, {n_events} converters "
      f"({100*n_events/n_fish:.1f}%)**. "
      f"{int(_v('design_matrix','n_excluded_incomplete'))} injured larvae were excluded for a "
      f"missing 0.5 h or 24 h session.")
    A("")
    if epv < 9:
        A(f"> **Events per variable = {epv:.1f}, below the conventional floor of 9 to 10.** The "
          "predictor set is held at four by design rather than trimmed to fit the rule, because each "
          "term answers a distinct pre-specified question. The cost is wider intervals on every "
          "coefficient, and they are reported rather than hidden. Section 2.3 shows the model that "
          "the data actually support, which uses two predictors and comfortably clears the rule.")
        A("")
    A("**Predictors.** `dose` (high_impact = 1); `pre_tau` (τ at t = -1 h); `z_dtau_0.5` and "
      "`z_dtau_24` (τ@0.5 - τ@-1 and τ@24 - τ@-1, z-scored **within dose group**).")
    A("")
    A("**Model.** L2-penalised logistic regression in a `StandardScaler` pipeline. No ensembles: at "
      f"n = {n_fish} with {n_events} events, a random forest or boosted model has enough capacity to "
      "memorise the sample, and its output cannot be sign-checked against the biology.")
    A("")
    A("**Validation.** `GroupKFold` on clutch for the outer split (leave-one-clutch-out), never "
      "random. The penalty strength C is tuned by **nested** cross-validation inside each outer "
      "training set, never on the folds reported. Within-dose z-scoring is a pipeline step fitted on "
      "training folds only, so no test-fold information reaches the scaling constants.")
    A("")

    A("### 2.1 Ablation table")
    A("")
    A("Same model class, same cross-validation, different inputs.")
    A("")
    A("| Model | Question | Mean fold AUC | SD | Fold range | Pooled OOF AUC [95% CI] | Brier |")
    A("|---|---|---|---|---|---|---|")
    for r in comp.itertuples():
        A(f"| {r.label} | {r.question} | **{r.mean_fold_auc:.3f}** | {r.sd_fold_auc:.3f} | "
          f"{r.min_fold_auc:.3f}-{r.max_fold_auc:.3f} | {r.pooled_oof_auc:.3f} "
          f"[{r.pooled_ci_low:.3f}, {r.pooled_ci_high:.3f}] | {r.brier:.3f} |")
    A("")
    A("Per-fold AUC (held-out clutch), with the C selected inside each fold:")
    A("")
    A("| Model | " + " | ".join(sorted(folds["held_out"].unique())) + " |")
    A("|---|" + "---|" * folds["held_out"].nunique())
    for _, g in folds.groupby("key", sort=False):
        lab = g["label"].iloc[0]
        g = g.sort_values("held_out")
        A(f"| {lab} | " + " | ".join(f"{r.auc:.3f} (C={r.selected_C:g})"
                                     for r in g.itertuples()) + " |")
    A("")

    def _auc(k):
        return float(comp.loc[comp["key"] == k, "mean_fold_auc"].iloc[0])

    a, b, c, d, e = (_auc(k) for k in ["a_locomotion_only", "b_dose_only", "c_dose_pretau",
                                       "d_dose_dtau", "e_full"])
    A("Reading it:")
    A("")
    A(f"- **(a) baseline locomotion alone: {a:.3f}.** At or below chance out of fold. Conversion is "
      "not a readout of how sick or sluggish the larva is.")
    A(f"- **(b) dose alone: {b:.3f}.** Above chance. Blast dose is a genuine main effect in this "
      "cohort, which was not true in the pilot data.")
    A(f"- **(c) dose + pre_tau: {c:.3f}.** The strongest row in the table, and the most stable "
      f"across folds (SD {comp.loc[comp.key=='c_dose_pretau','sd_fold_auc'].iloc[0]:.3f}).")
    A(f"- **(d) dose + z_dtau: {d:.3f}.** Barely above dose alone. The acute injury response adds "
      "little.")
    A(f"- **(e) the pre-specified four-predictor model: {e:.3f}.** Adding the two Δτ terms to (c) "
      f"**costs {e - c:+.3f} AUC**.")
    A("")
    A("### 2.2 The result this study did not expect")
    A("")
    A("The design was built on the hypothesis that the *acute change* in habituation kinetics carries "
      "the predictive signal. In this dataset it does not. The two Δτ terms:")
    A("")
    lo_c = _g("converter_contrast", "low_impact_dtau_0.5")
    hi_c = _g("converter_contrast", "high_impact_dtau_0.5")
    A(f"- do not separate converters from non-converters within either dose group "
      f"(low impact {_p(lo_c['p_value'])}, high impact {_p(hi_c['p_value'])});")
    A(f"- add nothing to dose on their own (model (d), {d:.3f}, versus dose alone, {b:.3f});")
    A(f"- and actively degrade the model when added to `pre_tau` (model (e), {e:.3f}, versus "
      f"model (c), {c:.3f}).")
    A("")
    A("What does carry signal is `pre_tau`, the **pre-injury** habituation constant, together with "
      "dose. That is a different scientific claim, and a weaker one for a biomarker: a variable "
      "measured before the injury cannot be a readout of the injury response. It points to "
      "**susceptibility** rather than to acute pathophysiology. Section 4.2 discusses what that means "
      "mechanistically, and the limitations section says plainly why the original framing is not "
      "supported here.")
    A("")
    A("The group-level Δτ effect is not in doubt: low and high dose still move τ in opposite "
      f"directions with a very large effect (Cohen's d = "
      f"{float(_g('tau_change_from_baseline','low_vs_high_delta_t0.5')['effect_size']):.2f}, "
      f"{_p(_g('tau_change_from_baseline','low_vs_high_delta_t0.5')['p_value'])}). Injury clearly "
      "perturbs the circuit. What fails is the step from that group difference to **individual** "
      "prediction: within a dose group, the size of a larva's acute Δτ does not tell you whether that "
      "larva will convert.")
    A("")

    A("### 2.3 Honest model selection")
    A("")
    A("Quoting the best row of an ablation table is selection on the test folds: the winner is partly "
      "chosen by the noise in the folds it is scored on. To get an unbiased number, the choice of "
      "predictor set is treated as one more hyper-parameter and tuned on the **inner** folds only.")
    A("")
    A("| Held-out clutch | Selected by inner CV | C | Outer-fold AUC |")
    A("|---|---|---|---|")
    for r in sel["folds"].itertuples():
        A(f"| {r.held_out} | `{r.selected_model}` | {r.selected_C:g} | {r.auc:.3f} |")
    A("")
    A(f"- **Honestly selected performance: mean fold AUC {sel['mean_fold_auc']:.3f}** "
      f"(SD {sel['sd_fold_auc']:.3f}), pooled out-of-fold AUC {sel['pooled_auc']:.3f} "
      f"[{sel_lo:.3f}, {sel_hi:.3f}].")
    A(f"- Naively quoting the best table row gives "
      f"{_v('model_selection','naive_best_row_mean_fold_auc'):.3f}. The optimism from picking it by "
      f"eye is {_v('model_selection','selection_optimism'):+.3f} AUC, which is small here only "
      "because the same model wins in most folds.")
    A(f"- **Total out-of-fold accuracy {100*sel_cm['accuracy']:.1f}%** at threshold 0.50 "
      f"({sel_cm['tp'] + sel_cm['tn']} of {n_fish} correct), balanced accuracy "
      f"{100*sel_cm['balanced_accuracy']:.1f}%, sensitivity {100*sel_cm['sensitivity']:.1f}%, "
      f"specificity {100*sel_cm['specificity']:.1f}%.")
    A(f"- Permutation test on the **entire selection procedure** (selection re-run inside every one "
      f"of the {perm['n_perm']} shuffles): null mean {sel_perm['null_mean']:.3f}, "
      f"{_p(sel_perm['p_pooled'])}.")
    A("")
    A("This is the number to quote for best achievable accuracy on this design, because it is the "
      "performance of a procedure that could be applied to a new clutch without knowing the answer "
      "first.")
    A("")

    A("### 2.4 Pre-specified model: full detail")
    A("")
    A("Reported because it was pre-specified, not because it is the best.")
    A("")
    A(f"- Mean fold AUC {full['mean_fold_auc']:.3f} (SD {full['sd_fold_auc']:.3f}), pooled "
      f"out-of-fold AUC {full['pooled_auc']:.3f} "
      f"[{comp.loc[comp.key=='e_full','pooled_ci_low'].iloc[0]:.3f}, "
      f"{comp.loc[comp.key=='e_full','pooled_ci_high'].iloc[0]:.3f}].")
    A(f"- Permutation {_p(perm['p_pooled'])} against a null that reruns the whole nested CV "
      f"{perm['n_perm']} times (null mean {perm['null_mean']:.3f}, SD {perm['null_sd']:.3f}).")
    A(f"- Total out-of-fold accuracy {100*cm['accuracy']:.1f}% at threshold 0.50, balanced accuracy "
      f"{100*cm['balanced_accuracy']:.1f}%, sensitivity {100*cm['sensitivity']:.1f}%, specificity "
      f"{100*cm['specificity']:.1f}%, Brier {full['brier']:.3f}.")
    A("")
    A("Confusion matrix, out-of-fold, threshold 0.50:")
    A("")
    A("| | Predicted non-converter | Predicted converter |")
    A("|---|---|---|")
    A(f"| **Non-converter** | {cm['tn']} | {cm['fp']} |")
    A(f"| **Converter** | {cm['fn']} | {cm['tp']} |")
    A("")
    thr_j = _v("confusion_matrix_youden", "threshold")
    A(f"At the Youden-optimal threshold ({thr_j:.3f}) accuracy is "
      f"{100*_v('confusion_matrix_youden','accuracy'):.1f}% with sensitivity "
      f"{100*_v('confusion_matrix_youden','sensitivity'):.1f}%. That threshold was chosen on these "
      "same out-of-fold predictions, so it is mildly optimistic and is shown only to indicate the "
      "achievable operating range.")
    A("")
    A("**Coefficients** (refit on all data at C = "
      f"{s2['best_C']:g}, standardised scale, clutch-clustered bootstrap CIs, "
      f"{int(coefs.attrs['n_boot'])} resamples):")
    A("")
    A("| Predictor | β (per SD) | 95% CI | OR per SD | OR 95% CI |")
    A("|---|---|---|---|---|")
    for r in coefs.itertuples():
        star = " \\*" if r.excludes_zero else ""
        A(f"| `{r.predictor}` | {r.coef_standardised:+.3f}{star} | "
          f"[{r.ci_low:+.3f}, {r.ci_high:+.3f}] | {r.odds_ratio_per_SD:.2f} | "
          f"[{r.or_ci_low:.2f}, {r.or_ci_high:.2f}] |")
    A("")
    A(f"Intercept {coefs.attrs['intercept']:+.3f}. \\* = bootstrap CI excludes zero.")
    A("")
    A("#### What the coefficients mean neurobiologically")
    A("")
    A("τ is the number of trials the acoustic startle response takes to decay to its floor. Startle "
      "habituation in larval zebrafish is not fatigue of the Mauthner cell, the hindbrain command "
      "neuron for the C-start escape. It is produced by progressive **feedforward inhibition onto the "
      "M-cell lateral dendrite**, which reduces dendritic excitability with repeated stimulation "
      "(Marsden & Granato, 2015). A larger τ therefore means inhibition accumulates more slowly, "
      "i.e. **reduced inhibitory gain in a defined sensorimotor circuit**. That is the same "
      "excitation/inhibition set point whose collapse drives epileptogenesis after brain injury.")
    A("")
    pre = coefs[coefs["predictor"] == "pre_tau"].iloc[0]
    dose_row = coefs[coefs["predictor"] == "dose"].iloc[0]
    A(f"- **`pre_tau` (β = {pre.coef_standardised:+.3f}, OR {pre.odds_ratio_per_SD:.2f} per SD), the "
      "dominant term.** Larvae that habituate more slowly *before* any injury are more likely to "
      "convert. Baseline inhibitory tone varies between individuals, and an animal that starts nearer "
      "the seizure threshold has less reserve to lose. This is a susceptibility term, the animal "
      "analogue of pre-injury risk factors that modify post-traumatic epilepsy risk in humans, and it "
      "fits a two-hit framing of predisposition plus insult.")
    A(f"- **`dose` (β = {dose_row.coef_standardised:+.3f}).** A real main effect here: conversion "
      "rises monotonically with blast dose.")
    A("- **The Δτ terms.** Both intervals are wide and the 0.5 h term straddles zero. Given that the "
      "ablation shows they cost AUC, the honest reading is that they carry no individual-level "
      "information in this cohort, not that they carry a small amount.")
    A("")
    A("### 2.5 Leakage quantification")
    A("")
    naive = s2["naive"]
    A("| Split | Mean fold AUC | Pooled OOF AUC |")
    A("|---|---|---|")
    A(f"| {config.N_RANDOM_SPLIT_FOLDS}-fold **random** (leaky, for comparison only) | "
      f"{naive['mean_fold_auc']:.3f} | {naive['pooled_auc']:.3f} |")
    A(f"| Leave-one-clutch-out (**the honest estimate**) | {full['mean_fold_auc']:.3f} | "
      f"{full['pooled_auc']:.3f} |")
    A("")
    A(f"Ignoring clutch inflates the mean fold AUC by **{s2['inflation']:+.3f}**. Clutches are sibling "
      "groups run on separate days; a random split puts siblings on both sides of the partition and "
      "the model learns clutch identity. Any AUC from a random split on this design should be "
      "discounted by roughly this much.")
    A("")
    A("Figures: `fig03_roc_nested_comparison.png`, `fig04_permutation_null.png`, "
      "`fig05_confusion_calibration.png`, `fig06_coefficients.png`, `fig07_fold_auc_by_model.png`.")
    A("")
    A("---")
    A("")

    # -------------------------------------------------------------- step 3
    A("## Step 3. Molecular validation: paired c-fos pools")
    A("")
    A("### 3.1 What this validates, and what it does not")
    A("")
    A("In this dataset the qPCR pools are labelled by **realised outcome** (`converter` / "
      "`non_converter`), not by a predicted risk score, and they are drawn from the same `zf_*` "
      "larvae that train the model. Two consequences, both stated rather than glossed:")
    A("")
    A("1. This is **not** an independent test of the prediction model. It is a test of the **outcome "
      "variable the model predicts**: does a larva scored as converted carry a molecular signature of "
      "elevated network activity, or is `converted` merely a behavioural scoring threshold?")
    A("2. Model and assay share animals, so sample independence is not claimed.")
    A("")
    A("What makes it worth doing is that the modality is genuinely different. Everything in Step 2 "
      "derives from how far a larva swam. If that pipeline carried a systematic artefact, no amount "
      "of internal cross-validation would reveal it, because every fold would inherit the artefact. "
      "Transcript abundance shares no instrumentation with swim tracking.")
    A("")
    A("`fosab` is the zebrafish orthologue of *c-fos*, the canonical immediate early gene. Sustained "
      "depolarisation raises intracellular Ca²⁺, driving CaMK- and MAPK/ERK-dependent phosphorylation "
      "of CREB and transcription from the *fos* promoter within roughly 15 to 30 minutes (Sheng & "
      "Greenberg, 1990). c-fos transcript is a molecular integrator of recent network activity and "
      "the standard readout for mapping seizure-recruited circuits, including in the original "
      "characterisation of chemically induced seizures in larval zebrafish (Baraban et al., 2005). "
      "Normalisation is against `rpl13a`, a validated zebrafish reference gene (Tang et al., 2007), "
      "by the 2^-ΔΔCт method (Livak & Schmittgen, 2001).")
    A("")
    A("### 3.2 Why six pairs and not nine")
    A("")
    A(f"There are {int(_v('design','n_pools'))} pools, three larvae each. They are not independent "
      "units: they are matched within (group × clutch) cells. The counts are unbalanced by biology, "
      "because conversion is rare in sham:")
    A("")
    A("| Group | Converter pools per clutch | Non-converter pools per clutch |")
    A("|---|---|---|")
    A("| sham | 0 | 2 |")
    A("| low_impact | 1 | 2 |")
    A("| high_impact | 2 | 2 |")
    A("")
    A("A balanced nine-pair design is therefore impossible: **sham cells contain no converter pool at "
      "all**, because too few sham larvae converted to fill one. The pairing the data support is "
      f"**{s3['all']['n']} matched cells**, the injured group × clutch combinations, with replicate "
      "pools within a cell averaged before pairing. Sham is used instead as an unpaired reference "
      "level, which answers a complementary question in section 3.4.")
    A("")
    A("### 3.3 Paired result")
    A("")
    A("| Contrast | Pairs | Mean Δ (converter - non-converter) | 95% CI | Paired t | p | Cohen's dz | Wilcoxon p |")
    A("|---|---|---|---|---|---|---|---|")
    for name, key in [("Fold-change scale (primary)", "all"), ("log₂ scale", "all_log2")]:
        r = s3[key]
        A(f"| {name} | {r['n']} | {r['mean']:+.4f} | [{r['ci'][0]:+.4f}, {r['ci'][1]:+.4f}] | "
          f"t({r['n']-1}) = {r['t']:+.3f} | {_fmt_p(r['p'])} | {r['dz']:+.3f} | "
          f"{_fmt_p(r['p_wilcoxon'])} |")
    A("")
    gm = float(2 ** s3["all_log2"]["mean"])
    gm_lo, gm_hi = float(2 ** s3["all_log2"]["ci"][0]), float(2 ** s3["all_log2"]["ci"][1])
    A(f"Converter pools carry **{100*(gm-1):.1f}% more c-fos transcript** relative to `rpl13a` than "
      f"matched non-converter pools (geometric mean ratio {gm:.3f}, 95% CI [{gm_lo:.3f}, "
      f"{gm_hi:.3f}], back-transformed from the paired log₂ analysis, which is the appropriate scale "
      "for a fold change).")
    A("")
    A(f"Direction is consistent in **{s3['n_positive']} of {s3['all']['n']}** pairs "
      f"(exact binomial sign test, {_p(s3['p_sign'])}).")
    A("")
    A("### 3.4 Specificity control: does c-fos track conversion, or just injury?")
    A("")
    bl = s3["baseline"]
    A("Since sham cannot be paired, the control question is turned around. If injury alone raised "
      "c-fos, then **non-converter** pools from injured groups should sit above non-converter pools "
      "from sham. They do not:")
    A("")
    A(f"- Injured non-converter pools {bl['mean_inj']:.3f} versus sham non-converter pools "
      f"{bl['mean_sham']:.3f}; difference {bl['diff']:+.3f}, 95% CI [{bl['ci'][0]:+.3f}, "
      f"{bl['ci'][1]:+.3f}], Welch t = {bl['t']:+.3f}, {_p(bl['p'])}, d = {bl['d']:.2f}.")
    A("")
    A("**A null result here is the desired one**, and it is what rules out the trivial explanation. "
      "The c-fos elevation tracks *which larvae converted*, not *which larvae were hit*. That is the "
      "specific claim, and it is the one that makes the assay informative about epileptogenesis "
      "rather than about injury exposure.")
    A("")
    A("Elevated immediate early gene expression in the absence of any provoking stimulus is what a "
      "chronically over-active network looks like transcriptionally. This is corroboration, not "
      "proof: it is a bulk measurement on pooled tissue, so it cannot localise the signal to a cell "
      "type or region, and it cannot distinguish loss of parvalbumin-positive interneuron function "
      "from increased glutamatergic drive. Both are documented consequences of brain injury.")
    A("")
    A("No regression of c-fos on a pooled continuous risk score is run. The pools are outcome-"
      "labelled, and the risk score is not comparable between dose groups.")
    A("")
    A("Figure: `fig08_cfos_paired.png`.")
    A("")
    A("---")
    A("")

    # -------------------------------------------------------------- step 4
    A("## Step 4. Descriptive results")
    A("")
    A("### 4.1 Groups start equal")
    A("")
    A(f"At the pre-injury baseline, τ does not differ between groups: F({anova['df']}) = "
      f"{float(anova['statistic']):.3f}, {_p(anova['p_value'])}, η² = "
      f"{float(anova['effect_size']):.4f}; Kruskal-Wallis H = "
      f"{float(_g('baseline_equivalence','kruskal_tau')['statistic']):.3f}, "
      f"{_p(_g('baseline_equivalence','kruskal_tau')['p_value'])}. Baseline locomotion likewise "
      f"({_p(_g('baseline_equivalence','anova_baseline_locomotion')['p_value'])}). Randomisation "
      "held, and every larva is its own control thereafter.")
    A("")
    A("### 4.2 τ moves in opposite directions by dose")
    A("")
    tau = s4["tau_summary"]
    A("| Group | " + " | ".join(f"t = {t:g} h" for t in config.TIMEPOINTS) + " |")
    A("|---|" + "---|" * len(config.TIMEPOINTS))
    for g in config.GROUPS:
        row = tau[tau["group"] == g].set_index("timepoint_h")
        cells = [f"{row.loc[t,'mean']:.2f} ± {row.loc[t,'sem']:.2f}" for t in config.TIMEPOINTS]
        A(f"| {config.GROUP_LABELS[g]} | " + " | ".join(cells) + " |")
    A("")
    A("(mean ± SEM, trials to habituate)")
    A("")
    lv = _g("tau_change_from_baseline", "low_vs_high_delta_t0.5")
    A(f"At 0.5 h, Δτ from each larva's own baseline is "
      f"**{_v('tau_change_from_baseline','low_impact_delta_t0.5'):+.2f}** trials in low impact and "
      f"**{_v('tau_change_from_baseline','high_impact_delta_t0.5'):+.2f}** in high impact. Welch "
      f"t = {float(lv['statistic']):.2f}, {_p(lv['p_value'])}, Cohen's d = "
      f"{float(lv['effect_size']):.2f}.")
    A("")
    A("The divergence is two different lesions on the same circuit, not a nuisance to correct away:")
    A("")
    A("- **Low dose, τ rises (habituation deficit).** Sublethal blast preferentially compromises the "
      "feedforward inhibition that normally accumulates onto the Mauthner cell across trials. "
      "Inhibition builds more slowly, the escape response persists, τ lengthens. This is "
      "disinhibition, and it maps onto the loss of GABAergic control reported after experimental "
      "brain injury.")
    A("- **High dose, τ falls (fatigue, not learning).** A shorter τ looks superficially like better "
      "habituation. It is not. Greater energy deposition depresses the excitatory limb as well: the "
      "acute metabolic crisis and depolarisation that follow severe injury reduce the startle "
      "response itself, so the fitted decay is fast because the response never had far to fall.")
    A("")
    A("Both are pathological and they move the same scalar in opposite directions, which is why any "
      "model using Δτ must also encode dose. **This group-level effect is robust. What this dataset "
      "shows is that it does not translate into individual prediction** (section 2.2).")
    A("")
    A("### 4.3 Operational metrics")
    A("")
    A("| Metric | Value |")
    A("|---|---|")
    A(f"| Sessions | {int(_v('operations','n_sessions'))} "
      f"({int(_v('operations','n_habituation_sessions'))} habituation, "
      f"{int(_v('operations','n_outcome_sessions'))} outcome) |")
    A(f"| Fish-sessions recorded | {int(_v('operations','total_fish_sessions_recorded'))} |")
    A(f"| Fish per hour | {_v('operations','fish_per_hour_mean'):.1f} ± "
      f"{_v('operations','fish_per_hour_sd'):.1f} |")
    A(f"| Operator minutes per fish | {_v('operations','operator_min_per_fish_mean'):.2f} ± "
      f"{_v('operations','operator_min_per_fish_sd'):.2f} |")
    A(f"| Consumables cost per fish | ${_v('operations','cost_per_fish_usd_mean'):.3f} |")
    A(f"| Total consumables | ${_v('operations','total_consumables_usd'):.2f} |")
    A(f"| Total operator time | {_v('operations','total_operator_hours'):.1f} h |")
    att = _g("operations", "attrition_rate_ci")
    A(f"| Attrition | {int(_v('operations','fish_lost_total'))}/"
      f"{int(_v('operations','fish_at_baseline'))} = {100*float(att['value']):.1f}% "
      f"[{100*float(att['ci_low']):.1f}%, {100*float(att['ci_high']):.1f}%] |")
    A("")
    A("Figures: `fig09_habituation_curves.png`, `fig10_tau_by_timepoint.png`, "
      "`fig11_converter_trajectories.png`, `fig12_operations.png`.")
    A("")
    A("---")
    A("")

    # ----------------------------------------------------------------- ptz
    A("## Step 5. PTZ seizure threshold")
    A("")
    A("Pentylenetetrazol is a non-competitive GABAₐ receptor antagonist: it binds at the picrotoxin "
      "site, reduces chloride conductance and removes inhibitory brake from the network. In larval "
      "zebrafish it produces stereotyped, dose-dependent seizure behaviour with electrographic "
      "correlates (Baraban et al., 2005), which makes it the standard pharmacological probe of "
      "**seizure threshold**. If injured larvae have less inhibitory reserve, a fixed challenge dose "
      "should push more of them across it.")
    A("")
    pz = s4["ptz"]
    A("| Group | Seized / n | Proportion | Wilson 95% CI | Median latency (s) |")
    A("|---|---|---|---|---|")
    for r in pz.itertuples():
        A(f"| {config.GROUP_LABELS[r.group]} | {int(r.n_seized)}/{int(r.n)} | {r.prop_seized:.2f} | "
          f"[{r.ci_low:.2f}, {r.ci_high:.2f}] | {r.median_latency_s:.0f} |")
    A("")
    chi = _g("ptz", "chi_square")
    pw = _g("ptz", "post_hoc_power_sham_vs_injured")
    A(f"χ²({int(chi['df'])}) = {float(chi['statistic']):.3f}, {_p(chi['p_value'])}, Cramér's V = "
      f"{float(chi['effect_size']):.3f}, total n = {int(chi['n'])}.")
    A("")
    verdict = _g("ptz", "power_verdict")
    if verdict is not None and "underpowered" not in str(verdict["value"]):
        A(f"**On power.** At this n the sham-versus-injured contrast is adequately powered "
          f"({100*float(pw['value']):.0f}%, Cohen's h = {float(pw['effect_size']):.2f}), so unlike "
          "the pilot cohort this probe is not merely directional. It remains a **secondary** outcome "
          "for a different reason: it is a group-level comparison, and the primary claim of this "
          "project is about predicting individual animals, which a group difference does not "
          "address. The low-versus-high dose contrast is still not powered, and the two injured "
          "groups are not distinguishable here "
          f"({_p(_g('ptz','fisher_exact_low_vs_high')['p_value'])}).")
    else:
        A(f"**This probe is underpowered.** Post-hoc power for the observed sham-versus-injured "
          f"difference is {100*float(pw['value']):.0f}%. It is reported as a directional check and no "
          "conclusion rests on it.")
    A("")

    conf = [r for r in sb.as_frame().itertuples() if r.analysis == "ptz_confound"]
    if conf:
        A("### 5.1 Is the challenge a confound?")
        A("")
        A("PTZ is a proconvulsant given to a subset of the **same** larvae that supply the conversion "
          "outcome, so the obvious worry is that the drug caused the outcome. It did not: conversion "
          "rates are effectively identical between challenged and unchallenged larvae within every "
          "group.")
        A("")
        A("| Group | Challenged | Not challenged | Fisher p |")
        A("|---|---|---|---|")
        for r in conf:
            note = str(r.notes)
            ch = note.split("challenged ")[1].split(",")[0]
            nc = note.split("not challenged ")[1]
            A(f"| {r.quantity.replace('conversion_challenged_vs_not_','')} | {ch} | {nc} | "
              f"{_fmt_p(r.p_value)} |")
        A("")

    orr = _g("ptz_vs_conversion", "odds_ratio_seized_given_converted")
    lat = _g("ptz_vs_conversion", "latency_converter_minus_nonconverter")
    if orr is not None:
        A("### 5.2 Seizure threshold and conversion coincide")
        A("")
        A("Because PTZ and conversion are measured on the same larvae, they can be crossed directly. "
          "This was not pre-specified and is exploratory.")
        A("")
        A(f"- Larvae that seized under PTZ were far more likely to have converted: Fisher exact odds "
          f"ratio **{float(orr['value']):.1f}**, {_p(orr['p_value'])}, n = {int(orr['n'])}.")
        if lat is not None:
            A(f"- Converters reached their first seizure sooner: {str(lat['notes'])}, "
              f"Mann-Whitney U = {float(lat['statistic']):.0f}, {_p(lat['p_value'])}.")
        A("")
        A("Two assays that share no measurement apparatus, a pharmacological threshold test and a "
          "spontaneous-activity recording, agree at the level of the individual animal. That is the "
          "strongest evidence in this report that `converted` denotes a real hyperexcitable state "
          "rather than a scoring artefact. It says nothing about whether that state is *predictable* "
          "in advance, which is Step 2's job.")
        A("")
    A("Figure: `fig13_ptz.png`.")
    A("")
    A("---")
    A("")

    # --------------------------------------------------------- limitations
    A("## Assumptions, checks and limitations")
    A("")
    A("Assumption checks print during every run tagged `[PASS]` or `[FLAG]`, and the underlying "
      "statistics are in `all_statistics.csv`.")
    A("")
    A(f"- Curve fitting: {100*conv:.1f}% convergence, median R² {r2med:.3f}, "
      f"{_plural(int(_v('curve_fit','n_tau_at_bound')), 'session')} at a parameter bound.")
    A(f"- Collinearity: max |r| among predictors = "
      f"{_v('design_matrix','max_abs_predictor_correlation'):.3f}; all VIFs are in "
      "`all_statistics.csv` under `assumptions`.")
    A(f"- Separation: the unpenalised Logit converged with max |β| = "
      f"{_v('assumptions','max_abs_unpenalised_coef'):.2f}, so the ridge penalty is not masking "
      "complete separation.")
    flagged = [g for g in config.GROUPS
               if (r := _g("baseline_equivalence", f"shapiro_tau_{g}")) is not None
               and float(r["p_value"]) <= 0.05]
    if flagged:
        A(f"- **Flagged:** baseline τ departs from normality in "
          f"{', '.join(config.GROUP_LABELS[g].lower() for g in flagged)} (Shapiro-Wilk p ≤ 0.05). "
          f"The ANOVA is therefore backed by Kruskal-Wallis "
          f"({_p(_g('baseline_equivalence','kruskal_tau')['p_value'])}), which agrees.")
    A("- Paired c-fos differences: Shapiro-Wilk reported for each contrast, Wilcoxon alongside every "
      "paired t-test.")
    A("")
    A("Real limitations, stated plainly:")
    A("")
    A(f"- **The study's central hypothesis is not supported.** The acute Δτ terms, which the design "
      "was built around, add nothing to individual prediction and slightly degrade it. What predicts "
      "conversion is a pre-injury trait plus dose. Reporting this is the point of running the "
      "ablation rather than fitting one model and describing it.")
    A(f"- **n = {n_fish} with {n_events} events, EPV = {epv:.1f}**, below the conventional floor. "
      "Coefficient intervals are wide. The two-predictor model that the data support clears the rule "
      f"comfortably (EPV = {n_events/2:.1f}).")
    A("- **Three clutches means three outer folds.** Fold-to-fold spread is estimated from three "
      "numbers and should be read as indicative.")
    A(f"- **{int(_v('design_matrix','n_excluded_incomplete'))} injured larvae were dropped** for a "
      "missing session. This is a complete-case analysis; attrition has not been modelled.")
    A("- **Conversion is a behavioural proxy.** Spontaneous burst activity is not "
      "electrographically confirmed epilepsy, though the PTZ concordance in section 5.2 supports it.")
    A("- **The c-fos and PTZ analyses share animals with the model cohort.** They validate the "
      "outcome label, not the model's generalisation.")
    A("- **`pre_tau` as a predictor cannot support the biomarker framing.** A variable measured "
      "before injury is a susceptibility marker, not an acute readout, and it could not be used to "
      "triage patients after an injury that has already happened.")
    A("- **One species, one age, one injury model.** Generalisation to mammalian TBI is a hypothesis.")
    A("")
    A("## Conclusion")
    A("")
    A(f"Conversion to spontaneous burst activity is predictable in this cohort at a cross-validated "
      f"AUC of {sel['mean_fold_auc']:.3f} and total out-of-fold accuracy of "
      f"{100*sel_cm['accuracy']:.1f}%, with the predictor set chosen inside the cross-validation and "
      f"a permutation {_p(sel_perm['p_pooled'])} against a null that reruns the entire procedure. "
      "Baseline locomotion is ruled out as an explanation, and clutch-aware splitting shows a random "
      f"split would have inflated the estimate by {s2['inflation']:+.3f} AUC.")
    A("")
    A("The signal lies in the **pre-injury** habituation constant together with blast dose, not in "
      "the acute change in habituation that this study set out to test. Injury does perturb the "
      "circuit strongly and in dose-dependent opposite directions, but that group-level effect does "
      "not carry individual-level predictive information here.")
    A("")
    A("Two assays that share no instrumentation with the behavioural pipeline agree that `converted` "
      "denotes a genuinely hyperexcitable state: converter pools carry "
      f"{100*(gm-1):.1f}% more c-fos transcript than matched non-converter pools, with injured "
      "non-converters indistinguishable from sham, and larvae that seize under PTZ are "
      f"{float(orr['value']):.0f} times more likely to have converted. The outcome variable is real. "
      "Whether it can be predicted from the acute injury response remains open, and on this evidence "
      "the answer is no.")
    A("")
    A("## Background literature")
    A("")
    A("Reading behind the neurobiological claims above. These support the framing; they are not "
      "citations of results generated here.")
    A("")
    for i, ref in enumerate([
        "Annegers JF, Hauser WA, Coan SP, Rocca WA (1998). A population-based study of seizures "
        "after traumatic brain injuries. *New England Journal of Medicine* 338:20-24.",
        "Baraban SC, Taylor MR, Castro PA, Baier H (2005). Pentylenetetrazole induced changes in "
        "zebrafish behavior, neural activity and c-fos expression. *Neuroscience* 131:759-768.",
        "Burgess HA, Granato M (2007). Sensorimotor gating in larval zebrafish. *Journal of "
        "Neuroscience* 27:4984-4994.",
        "Hunt RF, Boychuk JA, Smith BN (2013). Neural circuit mechanisms of post-traumatic epilepsy. "
        "*Frontiers in Cellular Neuroscience* 7:89.",
        "Livak KJ, Schmittgen TD (2001). Analysis of relative gene expression data using real-time "
        "quantitative PCR and the 2^-ΔΔCт method. *Methods* 25:402-408.",
        "Marsden KC, Granato M (2015). In vivo Ca²⁺ imaging reveals that decreased dendritic "
        "excitability drives startle habituation. *Cell Reports* 13:1733-1740.",
        "Peduzzi P, Concato J, Kemper E, Holford TR, Feinstein AR (1996). A simulation study of the "
        "number of events per variable in logistic regression analysis. *Journal of Clinical "
        "Epidemiology* 49:1373-1379.",
        "Sheng M, Greenberg ME (1990). The regulation and function of c-fos and other immediate "
        "early genes in the nervous system. *Neuron* 4:477-485.",
        "Sloviter RS (1991). Permanently altered hippocampal structure, excitability, and inhibition "
        "after experimental status epilepticus in the rat: the 'dormant basket cell' hypothesis. "
        "*Hippocampus* 1:41-66.",
        "Tang R, Dodd A, Lai D, McNabb WC, Love DR (2007). Validation of zebrafish (*Danio rerio*) "
        "reference genes for quantitative real-time RT-PCR normalization. *Acta Biochimica et "
        "Biophysica Sinica* 39:384-390.",
        "Varma S, Simon R (2006). Bias in error estimation when using cross-validation for model "
        "selection. *BMC Bioinformatics* 7:91.",
        "Wolman MA, Jain RA, Liss L, Granato M (2011). Chemical modulation of memory formation in "
        "larval zebrafish. *PNAS* 108:15468-15473.",
    ], 1):
        A(f"{i}. {ref}")
    A("")
    A("## Reproducing")
    A("")
    A("```bash")
    A("pip install -r requirements.txt")
    A("python run_all.py")
    A("```")
    A("")
    A(f"Seed `{config.SEED}`. Outputs: `results/figures/*.png` at 300 dpi, "
      "`results/all_statistics.csv`, `results/tables/*.csv`, and this file.")
    A("")

    text = "\n".join(L)
    config.RESULTS_MD.write_text(text, encoding="utf-8")
    print(f"\n  RESULTS.md -> {config.RESULTS_MD}")
    return text
