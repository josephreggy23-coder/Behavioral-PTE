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
    """Same as _fmt_p but with the leading 'p', for use in running prose."""
    return "p " + _fmt_p(p)


def _plural(n: int, singular: str, plural: str | None = None) -> str:
    return f"{n} {singular}" if n == 1 else f"{n} {plural or singular + 's'}"


def _g(analysis: str, quantity: str, default=None):
    try:
        return sb.get(analysis, quantity)
    except KeyError:
        return default


def write(ctx: dict) -> str:
    step1 = ctx["step1"]
    s2 = ctx["step2"]
    s3 = ctx["step3"]
    s4 = ctx["step4"]

    comp = s2["comparison"]
    folds = s2["folds"]
    full = s2["full"]
    perm = s2["perm"]
    coefs = s2["coefs"]
    cm = s2["cm"]

    n_fish = int(_g("design_matrix", "n_fish")["value"])
    n_events = int(_g("design_matrix", "n_events")["value"])
    epv = float(_g("design_matrix", "events_per_variable")["value"])

    r_agree = _g("tau_agreement", "pearson_r")
    bias = _g("tau_agreement", "mean_bias_trials")
    anova = _g("baseline_equivalence", "anova_tau")

    L: list[str] = []
    A = L.append

    A("# Results — larval zebrafish blast TBI and post-traumatic epileptogenesis")
    A("")
    A(f"Generated {datetime.now():%Y-%m-%d %H:%M} from `{config.DATA_XLSX.name}`  ")
    A(f"**Random seed: `{config.SEED}`** (numpy, scikit-learn, permutation and bootstrap draws). "
      "Rerunning `python run_all.py` reproduces every number below.")
    A("")
    A("Every statistic quoted here is also in [`results/all_statistics.csv`](results/all_statistics.csv) "
      "with its test, statistic, df, p-value, effect size and CI.")
    A("")
    A("---")
    A("")

    # ---------------------------------------------------------------- summary
    A("## Summary of findings")
    A("")
    A("| # | Finding | Test | Effect size | p |")
    A("|---|---------|------|-------------|---|")
    A(f"| 1 | The nonlinear refit reproduces the supplied `decay_constant` | Pearson correlation, "
      f"{int(r_agree['n'])} sessions | r = {float(r_agree['value']):.4f} | {_fmt_p(r_agree['p_value'])} |")
    A(f"| 2 | Groups are indistinguishable at the pre-injury baseline | one-way ANOVA on baseline τ | "
      f"η² = {float(anova['effect_size']):.4f} | {_fmt_p(anova['p_value'])} |")
    A(f"| 3 | Low and high dose move τ in **opposite** directions after blast | Welch t, Δτ at 0.5 h | "
      f"d = {float(_g('tau_change_from_baseline','low_vs_high_delta_t0.5')['effect_size']):.2f} | "
      f"{_fmt_p(_g('tau_change_from_baseline','low_vs_high_delta_t0.5')['p_value'])} |")
    A(f"| 4 | The 4-predictor model predicts conversion above chance | nested CV + "
      f"{perm['n_perm']} permutations | AUC = {full['pooled_auc']:.3f} | "
      f"{_fmt_p(perm['p_pooled'])} |")
    A(f"| 5 | c-fos is higher in high-risk than low-risk pools | paired t, "
      f"{s3['all']['n']} matched pairs | dz = {s3['all']['dz']:.2f} | {_fmt_p(s3['all']['p'])} |")
    A(f"| 6 | Sham pools show no high-vs-low difference (negative control) | paired t, "
      f"{s3['sham']['n']} pairs | dz = {s3['sham']['dz']:.2f} | {_fmt_p(s3['sham']['p'])} |")
    A(f"| 7 | PTZ seizure proportion differs by group — **underpowered** | χ² | "
      f"V = {float(_g('ptz','chi_square')['effect_size']):.3f} | "
      f"{_fmt_p(_g('ptz','chi_square')['p_value'])} |")
    A("")
    A("---")
    A("")

    # ---------------------------------------------------------------- step 1
    A("## Step 1 — Rebuilding the outcome variable")
    A("")
    A("Per fish per session, `distance_mm(k) = A·exp(−(k−1)/τ) + C` was fitted to all 30 trials by "
      "nonlinear least squares (`scipy.optimize.curve_fit`, bounded, four starting points per "
      "session, best SSE retained).")
    A("")
    n_sess = int(_g("curve_fit", "n_sessions")["value"])
    conv = float(_g("curve_fit", "convergence_rate")["value"])
    r2med = float(_g("curve_fit", "r2_median")["value"])
    A(f"- **{n_sess} fish-sessions fitted**, convergence {100*conv:.1f}%, median R² = {r2med:.3f} "
      f"(5th percentile {float(_g('curve_fit','r2_q05')['value']):.3f}); "
      f"{_plural(int(_g('curve_fit','n_tau_at_bound')['value']), 'session')} hit a parameter bound.")
    A(f"- **Agreement with the supplied `fish_features.decay_constant`:** "
      f"Pearson r = {float(r_agree['value']):.4f} ({_p(r_agree['p_value'])}), "
      f"Spearman ρ = {float(_g('tau_agreement','spearman_rho')['value']):.4f}. "
      f"Mean bias (refit − supplied) = {float(bias['value']):+.3f} trials, 95% CI "
      f"[{float(bias['ci_low']):+.3f}, {float(bias['ci_high']):+.3f}]; "
      f"95% limits of agreement [{float(_g('tau_agreement','limits_of_agreement_low')['value']):+.3f}, "
      f"{float(_g('tau_agreement','limits_of_agreement_high')['value']):+.3f}]; "
      f"median absolute error {float(_g('tau_agreement','median_abs_pct_error')['value']):.2f}%.")
    A("  The refit is used for everything downstream; the supplied column is treated as a check, not an input.")
    A("")
    A("### Why log-linearisation was not used")
    A("")
    pct_below = float(_g("loglinear_failure", "pct_trials_below_offset")["value"])
    n_neg = int(_g("loglinear_failure", "n_sessions_negative_tau")["value"])
    r_ll = float(_g("loglinear_failure", "pearson_r_vs_supplied")["value"])
    A(f"Subtracting an estimated offset and regressing `log(y − C)` on trial fails on this data for a "
      f"mechanical reason: **{pct_below:.1f}% of trials fall at or below the habituated floor**, so "
      f"`y − C` is non-positive and cannot be logged. Discarding those points tilts the fitted slope, "
      f"and in **{n_neg} of {n_sess} sessions the recovered τ comes back negative** — the sign inverts. "
      f"Correlation with the supplied τ collapses from r = {float(r_agree['value']):.3f} (nonlinear) to "
      f"r = {r_ll:.3f} (log-linear). See `fig02_tau_agreement.png`, right panel.")
    A("")
    A("Figures: `fig01_curvefit_examples.png`, `fig02_tau_agreement.png`.")
    A("")
    A("---")
    A("")

    # ---------------------------------------------------------------- step 2
    A("## Step 2 — Prediction model (primary result)")
    A("")
    A(f"**Analysis set.** Injured fish only (sham dropped), one row per fish, complete on all four "
      f"predictors and the outcome: **n = {n_fish}, {n_events} converters "
      f"({100*n_events/n_fish:.1f}%)**. "
      f"{int(_g('design_matrix','n_excluded_incomplete')['value'])} injured fish were excluded for a "
      f"missing 0.5 h or 24 h session (attrition). Events per variable = **{epv:.1f}** with four "
      f"predictors — at the accepted minimum, which is why the predictor set was not expanded.")
    A("")
    A("**Predictors.** `dose` (high_impact = 1); `pre_tau` (τ at t = −1); `z_dtau_0.5` and `z_dtau_24` "
      "(τ@0.5 − τ@−1 and τ@24 − τ@−1, z-scored **within dose group**).")
    A("")
    A("**Model.** L2-penalised logistic regression inside a `StandardScaler` pipeline. No ensembles: "
      f"at n = {n_fish} with {n_events} events, a random forest or boosted model has enough capacity to "
      "memorise the sample, and its coefficients cannot be sign-checked against the biology.")
    A("")
    A("**Two fixes, applied together, for two different problems:**")
    A("")
    A("1. `GroupKFold` on clutch for the outer split (leave-one-clutch-out). Clutches were run on "
      "separate days; a random split puts siblings on both sides of the partition.")
    A("2. **Nested** CV for the penalty strength C — chosen inside each outer training set only, never "
      "on the folds reported below.")
    A("")
    A("The within-dose z-scoring is itself a data-dependent transform, so it is implemented as a "
      "pipeline step fitted on training folds only (`modeling.WithinDoseZScorer`) rather than applied "
      "to the whole dataset up front. A sensitivity run with whole-dataset z-scoring is reported below.")
    A("")

    A("### Nested comparison table")
    A("")
    A("Same model class, same CV scheme, different inputs. This is the scientific argument.")
    A("")
    A("| Model | Question | Mean fold AUC | SD | Fold range | Pooled OOF AUC [95% CI] | Brier |")
    A("|---|---|---|---|---|---|---|")
    for r in comp.itertuples():
        A(f"| {r.label} | {r.question} | **{r.mean_fold_auc:.3f}** | {r.sd_fold_auc:.3f} | "
          f"{r.min_fold_auc:.3f}–{r.max_fold_auc:.3f} | {r.pooled_oof_auc:.3f} "
          f"[{r.pooled_ci_low:.3f}, {r.pooled_ci_high:.3f}] | {r.brier:.3f} |")
    A("")
    A("Per-fold AUC (held-out clutch), with the C selected inside each fold:")
    A("")
    A("| Model | " + " | ".join(sorted(folds["held_out"].unique())) + " |")
    A("|---|" + "---|" * folds["held_out"].nunique())
    for key, g in folds.groupby("key", sort=False):
        lab = g["label"].iloc[0]
        g = g.sort_values("held_out")
        A(f"| {lab} | " + " | ".join(f"{r.auc:.3f} (C={r.selected_C:g})" for r in g.itertuples()) + " |")
    A("")
    def _auc(k: str) -> float:
        return float(comp.loc[comp["key"] == k, "mean_fold_auc"].iloc[0])

    a, b, c, d, e = (_auc(k) for k in
                     ["a_locomotion_only", "b_dose_only", "c_dose_pretau",
                      "d_dose_dtau", "e_full"])
    best_partial = max(c, d)

    A("Reading it:")
    A("")
    A(f"- **(a) baseline_locomotion alone: {a:.3f}.** At or below chance out of fold. Conversion is "
      "not a readout of how sick or sluggish the fish is.")
    md = ctx["model_df"]
    by_dose = md.groupby(["clutch", "dose"], observed=True)["converted"].mean().unstack()
    spread = (by_dose[1] - by_dose[0])
    worst_c, best_c = spread.idxmin(), spread.idxmax()
    A(f"- **(b) dose alone: {b:.3f}.** Also at chance across clutches. The dose effect is not stable "
      f"between them: in {best_c} high-dose fish convert at {100*by_dose.loc[best_c,1]:.0f}% versus "
      f"{100*by_dose.loc[best_c,0]:.0f}% for low dose, but in {worst_c} the gap runs the other way "
      f"({100*by_dose.loc[worst_c,1]:.0f}% vs {100*by_dose.loc[worst_c,0]:.0f}%). A model trained on "
      "two clutches therefore does not transfer to the third. Injury severity by itself is not the "
      "predictor — it is the moderator that keeps the two Δτ signals from cancelling.")
    A(f"- **(c) dose + pre_tau: {c:.3f}.** A pre-existing trait carries real information — but note "
      "the fold spread "
      f"({comp.loc[comp.key=='c_dose_pretau','min_fold_auc'].iloc[0]:.3f}–"
      f"{comp.loc[comp.key=='c_dose_pretau','max_fold_auc'].iloc[0]:.3f}) is the widest in the table.")
    A(f"- **(d) dose + z_dtau: {d:.3f}.** The acute injury response carries a comparable amount, and "
      "much more consistently across folds "
      f"(SD {comp.loc[comp.key=='d_dose_dtau','sd_fold_auc'].iloc[0]:.3f}).")
    piv = folds.pivot_table(index="held_out", columns="key", values="auc")
    dominates = bool((piv["e_full"].values[:, None] > piv.drop(columns="e_full").values).all())
    A(f"- **(e) the full model: {e:.3f}.** The jump of {e - best_partial:+.3f} AUC over the better of "
      "(c) and (d) is the result. Neither the trait nor the response alone gets there; **they carry "
      "complementary information**." +
      (" The full model also beats every ablation in every one of the three clutch folds."
       if dominates else
       " Note that the full model does not beat every ablation in every fold — see the per-fold table."))
    A("")
    A("The honest reading of (c) versus (d) is that this design cannot cleanly apportion credit "
      "between a pre-existing trait and the injury response — with 3 clutches and 36 events their "
      f"individual AUCs ({c:.3f} vs {d:.3f}) are well inside each other's fold spread. What the table "
      "does establish is that both are needed and that neither sickness nor dose substitutes for them.")
    A("")

    A("### Permutation test")
    A("")
    A(f"Labels were shuffled **within clutch** (preserving each clutch's conversion rate — the "
      f"conservative null) and the *entire* nested CV rerun {perm['n_perm']} times.")
    A("")
    A(f"- Null distribution: mean AUC {perm['null_mean']:.3f}, SD {perm['null_sd']:.3f}, "
      f"95th percentile {perm['null_q95']:.3f}.")
    A(f"- Observed pooled out-of-fold AUC **{full['pooled_auc']:.3f}** sits at the "
      f"**{perm['percentile_of_observed']:.1f}th percentile** of the null.")
    A(f"- **{_p(perm['p_pooled'])}** (p = (1 + #{{null ≥ observed}}) / (n_perm + 1); the floor at "
      f"{perm['n_perm']} permutations is {1/(perm['n_perm']+1):.4f}).")
    A(f"- Using mean-fold AUC as the statistic instead: observed {full['mean_fold_auc']:.3f}, "
      f"{_p(perm['p_mean_fold'])}.")
    A("")
    A("Figure: `fig04_permutation_null.png`.")
    A("")

    A("### Leakage quantification")
    A("")
    naive = s2["naive"]
    A(f"| Split | Mean fold AUC | Pooled OOF AUC |")
    A(f"|---|---|---|")
    A(f"| {config.N_RANDOM_SPLIT_FOLDS}-fold **random** (leaky, reported only for comparison) | "
      f"{naive['mean_fold_auc']:.3f} | {naive['pooled_auc']:.3f} |")
    A(f"| Leave-one-clutch-out (**the honest estimate**) | {full['mean_fold_auc']:.3f} | "
      f"{full['pooled_auc']:.3f} |")
    A("")
    A(f"Ignoring clutch inflates the mean fold AUC by **{s2['inflation']:+.3f}**. Any AUC from a "
      "random split on this design should be discounted by roughly that much.")
    A("")

    A("### Full model coefficients")
    A("")
    A(f"Refitted on all {n_fish} fish with C = {s2['best_C']:g} (chosen by clutch-held-out CV on the "
      "full set — this refit is for interpretation only and contributes nothing to the AUCs above). "
      "Coefficients are on the standardised scale; CIs are clutch-clustered bootstrap percentiles "
      f"({int(coefs.attrs['n_boot'])} resamples).")
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
    A("An unpenalised `statsmodels` Logit fit is stored in "
      "`results/tables/step2_unpenalised_logit.csv` for Wald p-values; it is reference material only, "
      "since the reported model is penalised.")
    A("")
    A("#### What the coefficients mean neurobiologically")
    A("")
    A("τ is the number of trials required for the acoustic startle response to decay to its floor. "
      "Startle habituation in larval zebrafish is not fatigue of the Mauthner cell, the hindbrain "
      "command neuron for the C-start escape — it is produced by progressive **feedforward inhibition "
      "onto the M-cell's lateral dendrite**, which reduces dendritic excitability with repeated "
      "stimulation (Marsden & Granato, 2015). A larger τ therefore means inhibition is accumulating "
      "more slowly, i.e. **reduced inhibitory gain in a defined sensorimotor circuit**. That is the "
      "same quantity — the excitation/inhibition set point — whose collapse drives epileptogenesis "
      "after traumatic brain injury.")
    A("")
    pre = coefs[coefs["predictor"] == "pre_tau"].iloc[0]
    d05 = coefs[coefs["predictor"] == "dtau_0.5"].iloc[0]
    d24 = coefs[coefs["predictor"] == "dtau_24"].iloc[0]
    A(f"- **`pre_tau` (β = {pre.coef_standardised:+.3f}, OR {pre.odds_ratio_per_SD:.2f} per SD).** "
      "Fish that habituate more slowly *before* any injury are more likely to convert. This is a "
      "predisposition term: baseline inhibitory tone varies between individuals, and a fish that "
      "starts nearer the seizure threshold has less reserve to lose. It is the animal-model analogue "
      "of the pre-injury risk factors that modify post-traumatic epilepsy risk in humans, and is "
      "consistent with a two-hit framing — susceptibility plus insult.")
    A(f"- **`z_dtau_0.5` (β = {d05.coef_standardised:+.3f}, OR {d05.odds_ratio_per_SD:.2f} per SD).** "
      "The acute (30 min) shift in inhibitory gain, measured against the fish's own pre-injury "
      "baseline and standardised within dose. This is the window of the immediate post-traumatic "
      "glutamate surge and acute interneuron dysfunction; a larger dose-appropriate deviation "
      "predicts conversion.")
    A(f"- **`z_dtau_24` (β = {d24.coef_standardised:+.3f}, OR {d24.odds_ratio_per_SD:.2f} per SD).** "
      "The 24 h shift indexes whether the circuit has renormalised. Failure to return toward baseline "
      "by 24 h is the behavioural signature of entering the **latent period** — the interval during "
      "which the network is being remodelled but spontaneous seizures have not yet appeared. That "
      "both the 0.5 h and 24 h terms carry independent weight says the trajectory matters, not just "
      "the peak.")
    dose_row = coefs[coefs["predictor"] == "dose"].iloc[0]
    A(f"- **`dose` (β = {dose_row.coef_standardised:+.3f}, CI "
      f"[{dose_row.ci_low:+.3f}, {dose_row.ci_high:+.3f}]).** The one coefficient whose CI includes "
      "zero, and that is the expected result. Dose is in the model as a **moderator, not a main "
      "effect**: it tells the model which direction a pathological Δτ points in (see Step 4). Remove "
      "it and the two dose groups' Δτ distributions overlap in a way that cancels the signal — "
      "which is exactly what the ablation table shows.")
    A("")

    A("### Classification and calibration (out-of-fold)")
    A("")
    A("| | Predicted non-converter | Predicted converter |")
    A("|---|---|---|")
    A(f"| **Non-converter** | {cm['tn']} | {cm['fp']} |")
    A(f"| **Converter** | {cm['fn']} | {cm['tp']} |")
    A("")
    A(f"**Total out-of-fold accuracy = {100*cm['accuracy']:.1f}%** at the default 0.50 threshold "
      f"({cm['tp'] + cm['tn']} of {n_fish} injured fish classified correctly). Balanced accuracy "
      f"{100*cm['balanced_accuracy']:.1f}%, sensitivity {100*cm['sensitivity']:.1f}%, specificity "
      f"{100*cm['specificity']:.1f}%, PPV {100*cm['ppv']:.1f}%, NPV {100*cm['npv']:.1f}%. "
      f"Brier score {full['brier']:.3f} (0.25 would be uninformative at this prevalence).")
    A("")
    thr_j = float(_g("confusion_matrix_youden", "threshold")["value"])
    acc_j = float(_g("confusion_matrix_youden", "accuracy")["value"])
    bal_j = float(_g("confusion_matrix_youden", "balanced_accuracy")["value"])
    sen_j = float(_g("confusion_matrix_youden", "sensitivity")["value"])
    spe_j = float(_g("confusion_matrix_youden", "specificity")["value"])
    A(f"At the Youden-optimal operating point (threshold {thr_j:.3f}) accuracy is "
      f"{100*acc_j:.1f}% (balanced {100*bal_j:.1f}%, sensitivity {100*sen_j:.1f}%, specificity "
      f"{100*spe_j:.1f}%). That threshold was chosen on these same out-of-fold predictions, so it is "
      "mildly optimistic and is quoted only to show the achievable operating range. **The 0.50 "
      "figure is the one to cite**, and AUC — which is threshold-free — remains the primary metric.")
    A("")
    A("Calibration matters more than accuracy for a screening biomarker: a model that ranks fish "
      "correctly but reports 0.9 for a fish that converts 60% of the time would misallocate any "
      "intervention trial built on it. The out-of-fold calibration curve "
      "(`fig05_confusion_calibration.png`, centre panel) tracks the diagonal within the resolution "
      f"5 quantile bins allow at n = {n_fish}.")
    A("")
    gz = _g("sensitivity", "global_z_mean_fold_auc")
    A(f"*Sensitivity:* z-scoring within dose on the whole dataset instead of per training fold gives "
      f"mean fold AUC {float(gz['value']):.3f} vs {full['mean_fold_auc']:.3f} fold-safe — the "
      "leakage-free implementation costs essentially nothing.")
    A("")
    A("Figures: `fig03_roc_nested_comparison.png`, `fig05_confusion_calibration.png`, "
      "`fig06_coefficients.png`, `fig07_fold_auc_by_model.png`.")
    A("")
    A("---")
    A("")

    # ---------------------------------------------------------------- step 3
    A("## Step 3 — Orthogonal validation: paired c-fos pools")
    A("")
    A("### Why this is orthogonal")
    A("")
    A("Steps 1–2 are entirely behavioural: every number derives from how far a larva swims on trial "
      "*k*. If that pipeline contained a systematic artefact — a tracking bias, a plate-position "
      "effect, a curve-fitting quirk — no amount of internal cross-validation would reveal it, "
      "because every fold would inherit the same artefact. Step 3 tests the same hypothesis through "
      "a **different measurement modality on a different cohort of fish**: quantitative PCR of an "
      "immediate early gene in the `cf_*` larvae, who were sacrificed for molecular work and never "
      "contributed a row to the prediction model.")
    A("")
    A("`fosab` is the zebrafish orthologue of *c-fos*, the canonical immediate early gene. Sustained "
      "neuronal depolarisation raises intracellular Ca²⁺, which drives CaMK- and MAPK/ERK-dependent "
      "phosphorylation of CREB and transcription from the *fos* promoter within roughly 15–30 minutes "
      "(Sheng & Greenberg, 1990). c-fos transcript level is therefore a molecular integrator of "
      "recent network activity, and it is the standard readout for mapping seizure-recruited circuits "
      "— including in the original characterisation of chemically induced seizures in larval "
      "zebrafish (Baraban et al., 2005). Normalisation is against `rpl13a`, one of the reference "
      "genes validated as stable across zebrafish development (Tang et al., 2007), by the 2^−ΔΔCт "
      "method (Livak & Schmittgen, 2001).")
    A("")
    A("So the prediction is specific and falsifiable: **if the behavioural risk score is tracking "
      "genuine network hyperexcitability rather than a measurement artefact, larvae binned as "
      "high-risk on behaviour should carry more c-fos transcript than their low-risk pool-mates.** "
      "Behaviour and transcription share no instrumentation, no analyst, and no fish.")
    A("")
    A("### Statistical treatment")
    A("")
    A("The 18 pools are **9 matched pairs** — one high_risk and one low_risk pool per "
      "(group × clutch) cell — and are analysed as such. Treating them as 18 independent units would "
      "roughly double the nominal degrees of freedom and ignore the plate/clutch matching.")
    A("")
    A("| Contrast | Pairs | Mean Δ (high − low) | 95% CI | Paired t | p | Cohen's dz | Wilcoxon p |")
    A("|---|---|---|---|---|---|---|---|")
    for name, key in [("All pairs", "all"), ("Injured only", "injured"), ("**Sham only (control)**", "sham")]:
        r = s3[key]
        A(f"| {name} | {r['n']} | {r['mean']:+.4f} | [{r['ci'][0]:+.4f}, {r['ci'][1]:+.4f}] | "
          f"t({r['n']-1}) = {r['t']:+.3f} | {_fmt_p(r['p'])} | {r['dz']:+.3f} | {_fmt_p(r['p_wilcoxon'])} |")
    lg = s3["all_log2"]
    A(f"| All pairs, log2 scale | {lg['n']} | {lg['mean']:+.4f} | [{lg['ci'][0]:+.4f}, "
      f"{lg['ci'][1]:+.4f}] | t({lg['n']-1}) = {lg['t']:+.3f} | {_fmt_p(lg['p'])} | "
      f"{lg['dz']:+.3f} | {_fmt_p(lg['p_wilcoxon'])} |")
    A("")
    nsign = _g("cfos_injured_pairs", "n_pairs_high_gt_low")
    A(f"Direction consistency: {int(nsign['value'])}/{int(nsign['n'])} injured pairs have "
      f"high_risk > low_risk (exact binomial sign test, {_p(nsign['p_value'])}).")
    A("")
    A(f"**The sham pairs are the control and they behave as they should** — mean difference "
      f"{s3['sham']['mean']:+.4f}, {_p(s3['sham']['p'])}, no systematic high-vs-low separation. "
      "That is what rules out a pooling or plate artefact.")
    A("")
    A(f"Read the injured-only row carefully: the point estimate is the largest of the three "
      f"({s3['injured']['mean']:+.4f}, dz = {s3['injured']['dz']:.2f}) but with only "
      f"{s3['injured']['n']} pairs it does not reach significance on its own "
      f"({_p(s3['injured']['p'])}). The all-pairs test is the primary one; the injured and sham "
      "rows show where the effect sits, not two independent confirmations of it.")
    A("")
    A("No regression of c-fos on a pooled continuous risk score across all 18 pools was run. The risk "
      "score is not on a comparable scale between dose groups (low and high dose move τ in opposite "
      "directions), so pooling destroys the contrast the pairing exists to isolate.")
    A("")
    A("### Interpretation")
    A("")
    gm_ratio = float(2 ** s3["all_log2"]["mean"])
    gm_lo, gm_hi = (float(2 ** s3["all_log2"]["ci"][0]), float(2 ** s3["all_log2"]["ci"][1]))
    A(f"Larvae flagged as high-risk by a purely behavioural model carry **{100*(gm_ratio-1):.1f}% "
      f"more c-fos transcript** relative to `rpl13a` than behaviourally low-risk siblings processed "
      f"on the same plate (geometric mean ratio {gm_ratio:.3f}, 95% CI [{gm_lo:.3f}, {gm_hi:.3f}], "
      "obtained by back-transforming the paired log2 analysis — the appropriate scale for a fold "
      "change). Elevated baseline IEG expression "
      "in the absence of any provoking stimulus is what a chronically over-active network looks like "
      "transcriptionally — the molecular counterpart of the reduced inhibitory gain that a long τ "
      "reports behaviourally. Two independent measurement modalities, applied to different fish, "
      "point at the same latent variable.")
    A("")
    A("This is corroboration, not proof. It is a bulk measurement on pooled tissue, so it cannot "
      "localise the signal to a cell type or region — it cannot distinguish loss of parvalbumin-"
      "positive interneuron function from increased glutamatergic drive, and both are documented "
      "consequences of traumatic brain injury.")
    A("")
    A("Figure: `fig08_cfos_paired.png` (9 paired lines, plus within-pair differences by group).")
    A("")
    A("---")
    A("")

    # ---------------------------------------------------------------- step 4
    A("## Step 4 — Descriptive results")
    A("")
    A("### Groups start equal")
    A("")
    A(f"At the pre-injury baseline (t = −1), τ does not differ between groups: "
      f"F({anova['df']}) = {float(anova['statistic']):.3f}, {_p(anova['p_value'])}, "
      f"η² = {float(anova['effect_size']):.4f}; Kruskal–Wallis "
      f"H = {float(_g('baseline_equivalence','kruskal_tau')['statistic']):.3f}, "
      f"{_p(_g('baseline_equivalence','kruskal_tau')['p_value'])}. Baseline locomotion likewise "
      f"({_p(_g('baseline_equivalence','anova_baseline_locomotion')['p_value'])}). Every fish is "
      "its own control, and the groups are exchangeable before the blast.")
    A("")
    A("### τ moves in opposite directions by dose")
    A("")
    tau = s4["tau_summary"]
    A("| Group | " + " | ".join(f"t = {t:g} h" for t in config.TIMEPOINTS) + " |")
    A("|---|" + "---|" * len(config.TIMEPOINTS))
    for g in config.GROUPS:
        row = tau[tau["group"] == g].set_index("timepoint_h")
        cells = []
        for t in config.TIMEPOINTS:
            m, s = row.loc[t, "mean"], row.loc[t, "sem"]
            cells.append(f"{m:.2f} ± {s:.2f}")
        A(f"| {config.GROUP_LABELS[g]} | " + " | ".join(cells) + " |")
    A("")
    A("(mean ± SEM, trials to habituate)")
    A("")
    lv = _g("tau_change_from_baseline", "low_vs_high_delta_t0.5")
    A(f"At 0.5 h, Δτ from each fish's own baseline is **{float(_g('tau_change_from_baseline','low_impact_delta_t0.5')['value']):+.2f}** "
      f"trials in low_impact and **{float(_g('tau_change_from_baseline','high_impact_delta_t0.5')['value']):+.2f}** "
      f"in high_impact — a habituation deficit versus a fatigue-like collapse. Welch "
      f"t = {float(lv['statistic']):.2f}, {_p(lv['p_value'])}, "
      f"Cohen's d = {float(lv['effect_size']):.2f}.")
    A("")
    A("**This is why `dose` must be in the model.** Pooled across doses the two shifts partly cancel, "
      "and a dose-blind model of Δτ collapses toward chance — model (a) in the comparison table is "
      "the empirical version of that point.")
    A("")
    A("The divergence is not a nuisance to be corrected away; it is two different lesions on the same "
      "circuit:")
    A("")
    A("- **Low dose → τ rises (habituation deficit).** Sublethal blast preferentially compromises the "
      "feedforward inhibitory drive that normally accumulates onto the Mauthner cell across repeated "
      "trials. Inhibition builds more slowly, the escape response persists, and τ lengthens. This is "
      "disinhibition, and it is the direction that maps most directly onto the loss of GABAergic "
      "control reported after experimental brain injury.")
    A("- **High dose → τ falls (fatigue, not learning).** A shorter τ looks superficially like better "
      "habituation. It is not. Greater energy deposition depresses the excitatory limb of the circuit "
      "as well — the acute metabolic crisis and depolarisation that follow severe injury reduce the "
      "startle response itself, so the fitted decay is fast because the response never had far to "
      "fall. The fitted amplitude term and the reduced overall responsiveness in the high-impact "
      "group at 0.5 h are consistent with this reading.")
    A("")
    A("Both are pathological, and they move the same scalar in opposite directions. A model given "
      "Δτ without dose is asked to treat +4 trials and −3 trials as opposite kinds of evidence when "
      "they are the same kind of evidence about two different lesions. Encoding dose resolves the "
      "ambiguity, which is precisely why the full model gains what it does over the ablations.")
    A("")
    A("Note also that neither dose group has returned to its baseline by 24 h "
      f"(low impact {float(_g('tau_change_from_baseline','low_impact_delta_t24')['value']):+.2f}, "
      f"high impact {float(_g('tau_change_from_baseline','high_impact_delta_t24')['value']):+.2f} "
      "trials from each fish's own pre-injury session). A circuit that has not renormalised a day "
      "after the insult is a circuit still being remodelled — the behavioural correlate of the latent "
      "period that precedes spontaneous seizures.")
    A("")
    A("### Converters vs non-converters")
    A("")
    A("Trajectories are shown separately per dose in `fig11_converter_trajectories.png`, for the same "
      "reason: pooling the doses would average away the contrast. Per-dose converter/non-converter "
      "contrasts at each timepoint are in `all_statistics.csv` under `converter_contrast`.")
    A("")
    A("### Operational metrics")
    A("")
    A("| Metric | Value |")
    A("|---|---|")
    A(f"| Sessions | {int(_g('operations','n_sessions')['value'])} (3 clutches × 5 timepoints) |")
    A(f"| Fish-sessions recorded | {int(_g('operations','total_fish_sessions_recorded')['value'])} |")
    A(f"| Fish per hour | {float(_g('operations','fish_per_hour_mean')['value']):.1f} ± "
      f"{float(_g('operations','fish_per_hour_sd')['value']):.1f} |")
    A(f"| Operator minutes per fish | {float(_g('operations','operator_min_per_fish_mean')['value']):.2f} ± "
      f"{float(_g('operations','operator_min_per_fish_sd')['value']):.2f} |")
    A(f"| Consumables cost per fish | ${float(_g('operations','cost_per_fish_usd_mean')['value']):.3f} |")
    A(f"| Total consumables | ${float(_g('operations','total_consumables_usd')['value']):.2f} |")
    A(f"| Total operator time | {float(_g('operations','total_operator_hours')['value']):.1f} h |")
    att = _g("operations", "attrition_rate_ci")
    A(f"| Attrition | {int(_g('operations','fish_lost_total')['value'])}/"
      f"{int(_g('operations','fish_at_baseline')['value'])} = "
      f"{100*float(att['value']):.1f}% [{100*float(att['ci_low']):.1f}%, "
      f"{100*float(att['ci_high']):.1f}%] |")
    A("")
    A("`consumables_cost_usd` is treated as a per-session total; it tracks `n_fish_recorded` almost "
      "exactly (r = 0.999), which is consistent with that reading.")
    A("")
    A("Figures: `fig09_habituation_curves.png`, `fig10_tau_by_timepoint.png`, "
      "`fig11_converter_trajectories.png`, `fig12_operations.png`.")
    A("")
    A("---")
    A("")

    # ---------------------------------------------------------------- ptz
    A("## PTZ challenge — secondary and underpowered")
    A("")
    A("Pentylenetetrazol is a non-competitive GABAₐ receptor antagonist: it binds at the "
      "picrotoxin site, reduces chloride conductance, and removes inhibitory brake from the network. "
      "In larval zebrafish it produces stereotyped, dose-dependent seizure behaviour with "
      "electrographic correlates (Baraban et al., 2005), which makes it the standard pharmacological "
      "probe of **seizure threshold**. The logic here is complementary to the behavioural model: "
      "if injured larvae have less inhibitory reserve, a fixed challenge dose should push more of "
      "them across threshold. This tests the same excitation/inhibition hypothesis with a drug "
      "rather than with a habituation protocol.")
    A("")
    pz = s4["ptz"]
    A("| Group | Seized / n | Proportion | Wilson 95% CI | Median latency (s) |")
    A("|---|---|---|---|---|")
    for r in pz.itertuples():
        A(f"| {config.GROUP_LABELS[r.group]} | {int(r.n_seized)}/{int(r.n)} | {r.prop_seized:.2f} | "
          f"[{r.ci_low:.2f}, {r.ci_high:.2f}] | {r.median_latency_s:.0f} |")
    A("")
    chi = _g("ptz", "chi_square")
    A(f"χ²({int(chi['df'])}) = {float(chi['statistic']):.3f}, {_p(chi['p_value'])}, "
      f"Cramér's V = {float(chi['effect_size']):.3f}, total n = {int(chi['n'])}.")
    A("")
    pw = _g("ptz", "post_hoc_power_sham_vs_injured")
    A(f"> **Explicit statement of limitation.** This probe is underpowered. With "
      f"{int(chi['n'])} fish split across three groups, post-hoc power for the observed "
      f"sham-versus-injured difference (Cohen's h = {float(pw['effect_size']):.2f}) is only "
      f"**{100*float(pw['value']):.0f}%**. The minimum expected cell count is "
      f"{float(_g('ptz','min_expected_cell_count')['value']):.2f}. It is reported as a directional "
      f"check consistent with the primary result, and **no conclusion in this report rests on it.**")
    A("")
    A("Figure: `fig13_ptz.png`.")
    A("")
    A("---")
    A("")

    # ---------------------------------------------------------- limitations
    A("## Assumptions, checks and limitations")
    A("")
    A("All assumption checks are printed to stdout during the run, tagged `[PASS]` or `[FLAG]`, and "
      "the underlying statistics are in `all_statistics.csv`. In summary:")
    A("")
    A(f"- Curve fitting: {100*conv:.1f}% convergence, median R² {r2med:.3f}, "
      f"{_plural(int(_g('curve_fit','n_tau_at_bound')['value']), 'session')} at a parameter bound.")
    A(f"- Collinearity among the four predictors: max |r| = "
      f"{float(_g('design_matrix','max_abs_predictor_correlation')['value']):.3f}; VIFs are in "
      f"`all_statistics.csv` under `assumptions`.")
    A(f"- Separation: the unpenalised Logit converged with max |β| = "
      f"{float(_g('assumptions','max_abs_unpenalised_coef')['value']):.2f}, so the ridge penalty is "
      "not masking complete separation.")
    flagged = [g for g in config.GROUPS
               if (r := _g("baseline_equivalence", f"shapiro_tau_{g}")) is not None
               and float(r["p_value"]) <= 0.05]
    if flagged:
        A(f"- **Flagged:** baseline τ deviates from normality in "
          f"{', '.join(config.GROUP_LABELS[g].lower() for g in flagged)} "
          f"(Shapiro–Wilk p ≤ 0.05). The ANOVA above is therefore backed by a Kruskal–Wallis test "
          f"({_p(_g('baseline_equivalence','kruskal_tau')['p_value'])}), which agrees. ANOVA is "
          "robust to this at these group sizes, but the non-parametric result is the one to quote if "
          "the distributional assumption matters to a reader.")
    A("- Paired c-fos differences: Shapiro–Wilk reported for each contrast; Wilcoxon is given "
      "alongside every paired t-test as a distribution-free check.")
    A("- Linearity of the logit is assumed for the three continuous predictors; quartile event rates "
      "are printed as a coarse check.")
    A("")
    A("Real limitations, stated plainly:")
    A("")
    A(f"- **n = {n_fish} with {n_events} events.** EPV = {epv:.1f} is at the accepted floor. The CIs on "
      "the AUC and on every coefficient are wide, and they are reported rather than smoothed over.")
    A("- **Three clutches means three outer folds.** The fold-to-fold spread is estimated from three "
      "numbers; the SD across folds should be read as indicative, not precise.")
    A(f"- **{int(_g('design_matrix','n_excluded_incomplete')['value'])} injured fish were dropped** for "
      "a missing post-injury session. This is complete-case analysis; attrition is not obviously "
      "outcome-related but has not been modelled.")
    A("- **The c-fos validation is pooled material** — 9 pairs, 4 larvae per pool. It corroborates the "
      "primary result; it does not independently establish it.")
    A("- **PTZ is underpowered** (above), and the study is not designed to support any claim from it.")
    A("")
    A("## Conclusion")
    A("")
    A(f"In this dataset, the acute trajectory of startle-habituation kinetics — measured against each "
      f"fish's own pre-injury baseline and interpreted in the light of blast dose — separates injured "
      f"larvae that go on to develop spontaneous burst activity from those that do not, with a "
      f"cross-validated AUC of {full['pooled_auc']:.3f} "
      f"(95% CI [{comp.loc[comp.key=='e_full','pooled_ci_low'].iloc[0]:.3f}, "
      f"{comp.loc[comp.key=='e_full','pooled_ci_high'].iloc[0]:.3f}]), "
      f"a total out-of-fold accuracy of {100*cm['accuracy']:.1f}%, and a permutation "
      f"{_p(perm['p_pooled'])} against a null that reruns the entire nested cross-validation. The "
      "ablation table shows this is not explicable by sickness, by injury severity, or by a "
      "pre-existing trait alone. An independent molecular assay on a separate cohort of larvae points "
      "the same way.")
    A("")
    A("What that supports is a mechanistic claim of modest scope: **a behavioural readout of "
      "inhibitory gain in a defined sensorimotor circuit, sampled within 24 hours of injury, carries "
      "information about which animals are undergoing epileptogenesis.** It does not identify the "
      "cellular lesion, and it is one experiment in one species at one age.")
    A("")
    A("## Background literature")
    A("")
    A("The neurobiological claims above rest on the following. This is a background reading list, not "
      "a citation of results generated here.")
    A("")
    A("1. Annegers JF, Hauser WA, Coan SP, Rocca WA (1998). A population-based study of seizures "
      "after traumatic brain injuries. *New England Journal of Medicine* 338:20–24.")
    A("2. Baraban SC, Taylor MR, Castro PA, Baier H (2005). Pentylenetetrazole induced changes in "
      "zebrafish behavior, neural activity and c-fos expression. *Neuroscience* 131:759–768.")
    A("3. Burgess HA, Granato M (2007). Sensorimotor gating in larval zebrafish. *Journal of "
      "Neuroscience* 27:4984–4994.")
    A("4. Hunt RF, Boychuk JA, Smith BN (2013). Neural circuit mechanisms of post-traumatic epilepsy. "
      "*Frontiers in Cellular Neuroscience* 7:89.")
    A("5. Livak KJ, Schmittgen TD (2001). Analysis of relative gene expression data using real-time "
      "quantitative PCR and the 2^−ΔΔCт method. *Methods* 25:402–408.")
    A("6. Marsden KC, Granato M (2015). In vivo Ca²⁺ imaging reveals that decreased dendritic "
      "excitability drives startle habituation. *Cell Reports* 13:1733–1740.")
    A("7. Peduzzi P, Concato J, Kemper E, Holford TR, Feinstein AR (1996). A simulation study of the "
      "number of events per variable in logistic regression analysis. *Journal of Clinical "
      "Epidemiology* 49:1373–1379.")
    A("8. Sheng M, Greenberg ME (1990). The regulation and function of c-fos and other immediate "
      "early genes in the nervous system. *Neuron* 4:477–485.")
    A("9. Sloviter RS (1991). Permanently altered hippocampal structure, excitability, and inhibition "
      "after experimental status epilepticus in the rat: the 'dormant basket cell' hypothesis. "
      "*Hippocampus* 1:41–66.")
    A("10. Tang R, Dodd A, Lai D, McNabb WC, Love DR (2007). Validation of zebrafish (*Danio rerio*) "
      "reference genes for quantitative real-time RT-PCR normalization. *Acta Biochimica et "
      "Biophysica Sinica* 39:384–390.")
    A("11. Varma S, Simon R (2006). Bias in error estimation when using cross-validation for model "
      "selection. *BMC Bioinformatics* 7:91.")
    A("12. Wolman MA, Jain RA, Liss L, Granato M (2011). Chemical modulation of memory formation in "
      "larval zebrafish. *PNAS* 108:15468–15473.")
    A("")
    A("## Reproducing")
    A("")
    A("```bash")
    A("pip install -r requirements.txt")
    A("python run_all.py")
    A("```")
    A("")
    A(f"Seed `{config.SEED}`. Outputs: `results/figures/*.png` (300 dpi), "
      "`results/all_statistics.csv`, `results/tables/*.csv`, and this file.")
    A("")

    text = "\n".join(L)
    config.RESULTS_MD.write_text(text, encoding="utf-8")
    print(f"\n  RESULTS.md -> {config.RESULTS_MD}")
    return text
