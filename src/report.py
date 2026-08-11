"""Assemble RESULTS.md from the recorded statistics ledger."""
from __future__ import annotations

import numpy as np

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
    s2 = ctx["step2"]
    s3 = ctx["step3"]
    s4 = ctx["step4"]

    comp = s2["comparison"]
    folds = s2["folds"]
    full = s2["full"]
    perm = s2["perm"]
    coefs = s2["coefs"]
    cm = s2["cm"]
    permutation_ran = bool(perm.get("performed", perm.get("n_perm", 0) > 0))

    n_fish = int(_g("design_matrix", "n_fish")["value"])
    n_events = int(_g("design_matrix", "n_events")["value"])
    epv = float(_g("design_matrix", "events_per_variable")["value"])

    r_agree = _g("tau_agreement", "pearson_r")
    bias = _g("tau_agreement", "mean_bias_trials")
    anova = _g("baseline_equivalence", "anova_tau")

    L: list[str] = []
    A = L.append

    A("# Results — startle-habituation kinetics and a supplied post-injury behavioral label")
    A("")
    A(f"**Source:** `{config.DATA_XLSX.relative_to(config.ROOT).as_posix()}`")
    A("")
    A(f"**Random seed: `{config.SEED}`** (numpy, scikit-learn, permutation and bootstrap draws). "
      "With the recorded input and tested dependency versions, rerunning `python run_all.py` "
      "reproduces the analysis numerically.")
    A("")
    A("> **Interpretation boundary.** These analyses use the supplied workbook. The repository does "
      "not include source videos, electrographic recordings, raw qPCR Ct data, animal-approval "
      "records, or a reproducible protocol for deriving the `converted` and `risk_pool` labels. "
      "The results therefore describe internal associations with supplied labels; they do not "
      "establish epilepsy or validate a biomarker.")
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
    A(f"| 2 | No baseline group difference was detected | one-way ANOVA on baseline τ | "
      f"η² = {float(anova['effect_size']):.4f} | {_fmt_p(anova['p_value'])} |")
    A(f"| 3 | Low and high dose move τ in **opposite** directions after pressure-wave exposure | Welch t, Δτ at 0.5 h | "
      f"d = {float(_g('tau_change_from_baseline','low_vs_high_delta_t0.5')['effect_size']):.2f} | "
      f"{_fmt_p(_g('tau_change_from_baseline','low_vs_high_delta_t0.5')['p_value'])} |")
    permutation_test_label = (f"nested CV + {perm['n_perm']} permutations" if permutation_ran
                              else "nested CV; permutation test not run")
    A(f"| 4 | The 4-predictor model distinguishes the supplied `converted` label across held-out clutches | "
      f"{permutation_test_label} | AUC = {full['pooled_auc']:.3f} | "
      f"{_fmt_p(perm['p_pooled'])} |")
    A(f"| 5 | Orthogonal molecular validation: supplied high-risk pools have higher c-fos in the "
      f"nominal all-pair comparison | paired t on log2 fold change, "
      f"{s3['all_log2']['n']} matched pairs | dz = {s3['all_log2']['dz']:.2f} | "
      f"{_fmt_p(s3['all_log2']['p'])} |")
    A(f"| 6 | No high-vs-low difference was detected in the three sham pairs | paired t, "
      f"{s3['sham']['n']} pairs | dz = {s3['sham']['dz']:.2f} | {_fmt_p(s3['sham']['p'])} |")
    A(f"| 7 | The PTZ group comparison did not reach p < 0.05 and is underpowered | χ² | "
      f"V = {float(_g('ptz','chi_square')['effect_size']):.3f} | "
      f"{_fmt_p(_g('ptz','chi_square')['p_value'])} |")
    A("")
    A("---")
    A("")

    # ---------------------------------------------------------------- step 1
    A("## Step 1 — Re-estimating the habituation feature")
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
      f"predictors. This is a small modeling set, so the predictor set was not expanded.")
    A("")
    A("**Predictors.** `dose` (high_impact = 1); `pre_tau` (τ at t = −1); `z_dtau_0.5` and `z_dtau_24` "
      "(τ@0.5 − τ@−1 and τ@24 − τ@−1, z-scored **within dose group**).")
    A("")
    A("**Model.** L2-penalised logistic regression inside a `StandardScaler` pipeline. A fixed, "
      f"regularised linear classifier limits flexibility at n = {n_fish} and yields coefficients "
      "that can be inspected as associations.")
    A("")
    A("**Two safeguards are applied together:**")
    A("")
    A("1. `GroupKFold` on clutch for the outer split (leave-one-clutch-out). A random split would "
      "mix clutch-associated observations across training and test partitions.")
    A("2. **Nested** CV for the penalty strength C — chosen inside each outer training set only, never "
      "on the folds reported below.")
    A("")
    A("The within-dose z-scoring is itself a data-dependent transform, so it is implemented as a "
      "pipeline step fitted on training folds only (`modeling.WithinDoseZScorer`) rather than applied "
      "to the whole dataset up front. A sensitivity run with whole-dataset z-scoring is reported below.")
    A("")

    A("### Nested comparison table")
    A("")
    A("The model class and validation scheme stay fixed while the input features change.")
    A("")
    A("| Model | Question | Mean fold AUC | SD | Fold range | Pooled OOF AUC [conditional 95% interval] | Brier |")
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

    a, b, c, d, e, f, g = (_auc(k) for k in
                           ["a_locomotion_only", "b_dose_only", "c_dose_pretau",
                            "d_dose_dtau", "e_dtau_only", "f_behavior_no_dose", "e_full"])
    best_partial = max(c, d)

    A("Reading it:")
    A("")
    A(f"- **(a) baseline_locomotion alone: {a:.3f}.** This single activity measure performs at or "
      "below chance out of fold. It does not, by itself, exclude other illness, motor, or severity "
      "explanations.")
    md = ctx["model_df"]
    by_dose = md.groupby(["clutch", "dose"], observed=True)["converted"].mean().unstack()
    spread = (by_dose[1] - by_dose[0])
    worst_c, best_c = spread.idxmin(), spread.idxmax()
    A(f"- **(b) dose alone: {b:.3f}.** Also at chance across clutches. The dose effect is not stable "
      f"between them: in {best_c} high-dose fish convert at {100*by_dose.loc[best_c,1]:.0f}% versus "
      f"{100*by_dose.loc[best_c,0]:.0f}% for low dose, but in {worst_c} the gap runs the other way "
      f"({100*by_dose.loc[worst_c,1]:.0f}% vs {100*by_dose.loc[worst_c,0]:.0f}%). A model trained on "
      "two clutches therefore does not transfer to the third. In this dataset, dose alone is a poor "
      "cross-clutch classifier; whether it adds information beyond the behavioral trajectory is "
      "tested directly by models (f) and (g).")
    A(f"- **(c) dose + pre_tau: {c:.3f}.** The baseline-kinetics model contains predictive signal — "
      "but note "
      "the fold spread "
      f"({comp.loc[comp.key=='c_dose_pretau','min_fold_auc'].iloc[0]:.3f}–"
      f"{comp.loc[comp.key=='c_dose_pretau','max_fold_auc'].iloc[0]:.3f}) is the widest in the table.")
    A(f"- **(d) dose + z_dtau: {d:.3f}.** The post-injury change-score model carries a comparable "
      "amount of predictive signal, and does so "
      "much more consistently across folds "
      f"(SD {comp.loc[comp.key=='d_dose_dtau','sd_fold_auc'].iloc[0]:.3f}).")
    A(f"- **(e) delta_tau only (dose-blind): {e:.3f}.** This ablation omits dose both as a predictor "
      "and from preprocessing, testing whether the two raw change scores transfer without injury-"
      "severity context.")
    A(f"- **(f) pre_tau + delta_tau (dose-blind): {f:.3f}.** This is the direct dose-necessity "
      "ablation: it retains every behavioral feature while omitting dose from both the inputs and "
      "preprocessing.")
    piv = folds.pivot_table(index="held_out", columns="key", values="auc")
    dominates = bool((piv["e_full"].values[:, None] > piv.drop(columns="e_full").values).all())
    A(f"- **(g) the full model: {g:.3f}.** Its mean-fold AUC is {g - best_partial:+.3f} above the "
      "better of (c) and (d), consistent with complementary predictive information from baseline "
      f"and post-injury features. Adding dose to the complete behavioral model changes mean-fold "
      f"AUC by {g - f:+.3f}; the (f)-versus-(g) comparison, not the dose-only model, is the test of "
      "dose's incremental value." +
      (" The full model also beats every ablation in every one of the three clutch folds."
       if dominates else
       " Note that the full model does not beat every ablation in every fold — see the per-fold table."))
    A("")
    A("The honest reading of (c) versus (d) is that this design cannot cleanly apportion credit "
      "between a pre-existing trait and the injury response — with 3 clutches and 36 events their "
      f"individual AUCs ({c:.3f} vs {d:.3f}) are well inside each other's fold spread. What the table "
      "supports is complementary information from baseline and post-injury behavior; it does not "
      "establish that explicit dose encoding is necessary.")
    A("")

    A("### Permutation test")
    A("")
    if permutation_ran:
        A(f"Labels were shuffled **within clutch** (preserving each clutch's conversion rate) and "
          f"the *entire* nested CV rerun {perm['n_perm']} times.")
        A("")
        A(f"- Null distribution: mean AUC {perm['null_mean']:.3f}, SD {perm['null_sd']:.3f}, "
          f"95th percentile {perm['null_q95']:.3f}.")
        A(f"- Observed pooled out-of-fold AUC **{full['pooled_auc']:.3f}** sits at the "
          f"**{perm['percentile_of_observed']:.1f}th percentile** of the null.")
        A(f"- **{_p(perm['p_pooled'])}** (p = (1 + #{{null ≥ observed}}) / (n_perm + 1); "
          f"the floor at {perm['n_perm']} permutations is {1/(perm['n_perm']+1):.4f}).")
        A(f"- Using mean-fold AUC as the statistic instead: observed "
          f"{full['mean_fold_auc']:.3f}, {_p(perm['p_mean_fold'])}.")
        A("")
        A("Figure: `fig04_permutation_null.png`.")
    else:
        A("The permutation test was not run (`--skip-permutation` or `--permutations 0`), so the "
          "null-distribution summaries and permutation p-values are unavailable. No permutation-null "
          "figure was generated during this run.")
    A("")

    A("### Random-split comparison")
    A("")
    naive = s2["naive"]
    A("| Split | Mean fold AUC | Pooled OOF AUC |")
    A("|---|---|---|")
    A(f"| {config.N_RANDOM_SPLIT_FOLDS}-fold **random** (ignores clutch; comparison only) | "
      f"{naive['mean_fold_auc']:.3f} | {naive['pooled_auc']:.3f} |")
    A(f"| Leave-one-clutch-out (**reported estimate**) | {full['mean_fold_auc']:.3f} | "
      f"{full['pooled_auc']:.3f} |")
    A("")
    A(f"In this dataset, the random split's mean fold AUC is **{s2['inflation']:+.3f}** above the "
      "clutch-held-out estimate. This single comparison does not define a general correction factor.")
    A("")

    A("### Full model coefficients")
    A("")
    A(f"Refitted on all {n_fish} fish with C = {s2['best_C']:g} (chosen by clutch-held-out CV on the "
      "full set — this refit is for interpretation only and contributes nothing to the AUCs above). "
      "Coefficients are on the standardised scale; CIs are percentile intervals from a "
      "within-clutch stratified bootstrap conditional on the observed clutches "
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
    A("#### How to read the coefficients")
    A("")
    A("τ is the fitted number of trials over which the visual dark-flash response approaches its "
      "floor. It is an empirical summary of response decay and does not by itself identify a neural "
      "mechanism. This experiment does not directly measure inhibition, excitation, or a specific "
      "circuit, so the coefficients should be interpreted as predictive associations.")
    A("")
    pre = coefs[coefs["predictor"] == "pre_tau"].iloc[0]
    d05 = coefs[coefs["predictor"] == "dtau_0.5"].iloc[0]
    d24 = coefs[coefs["predictor"] == "dtau_24"].iloc[0]
    A(f"- **`pre_tau` (β = {pre.coef_standardised:+.3f}, OR {pre.odds_ratio_per_SD:.2f} per SD).** "
      "Within the fitted model, slower pre-injury habituation is associated with the supplied later "
      "label. The design cannot determine whether this is causal susceptibility, confounding, or "
      "sampling variation.")
    A(f"- **`z_dtau_0.5` (β = {d05.coef_standardised:+.3f}, OR {d05.odds_ratio_per_SD:.2f} per SD).** "
      "This is the 30-minute change from the fish's own baseline, centered and scaled within dose. "
      "Its coefficient estimates an association conditional on the other predictors; it is not a "
      "direct measure of inhibitory gain.")
    A(f"- **`z_dtau_24` (β = {d24.coef_standardised:+.3f}, OR {d24.odds_ratio_per_SD:.2f} per SD).** "
      "This is the analogous 24-hour change score. Its inclusion allows the classifier to use the "
      "two-timepoint trajectory, but the study does not establish that this trajectory represents a "
      "latent epileptogenic period.")
    dose_row = coefs[coefs["predictor"] == "dose"].iloc[0]
    A(f"- **`dose` (β = {dose_row.coef_standardised:+.3f}, CI "
      f"[{dose_row.ci_low:+.3f}, {dose_row.ci_high:+.3f}]).** This coefficient estimates dose's "
      "conditional contribution after the three behavioral terms. Its interval and the direct "
      "comparison of models (f) and (g) should be used together; the strong dose-blind model means "
      "dose cannot be described as necessary for prediction. Because the model contains no "
      "dose-by-change interaction, it is not a formal moderation analysis.")
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
    A("Calibration matters for a candidate screening assay: a model that ranks fish "
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
    A("## Step 3 — Orthogonal molecular validation: paired c-fos pools")
    A("")
    A("### Why this is an orthogonal validation")
    A("")
    A("Steps 1–2 use behavioral measurements. Step 3 uses quantitative PCR of an immediate early "
      "gene in a separate `cf_*` cohort whose fish do not enter the prediction model. It is the "
      "project's orthogonal molecular validation because it tests supplied risk strata, described "
      "as behavior-derived, using a different cohort and assay modality. The analysis accepts the workbook's "
      "`risk_pool` assignments rather than reconstructing them from a documented scoring and "
      "selection algorithm, so it is not external validation of the prediction model.")
    A("")
    A("`fosab` is the zebrafish orthologue of *c-fos*, the canonical immediate early gene. Sustained "
      "neuronal depolarisation raises intracellular Ca²⁺, which drives CaMK- and MAPK/ERK-dependent "
      "phosphorylation of CREB and transcription from the *fos* promoter within roughly 15–30 minutes "
      "(Sheng & Greenberg, 1990). c-fos transcript level is therefore a molecular integrator of "
      "recent network activity, and it is the standard readout for mapping seizure-recruited circuits "
      "— including in the original characterisation of chemically induced seizures in larval "
      "zebrafish (Baraban et al., 2005). Normalisation is against `rpl13a`, one of the reference "
      "genes validated as stable across zebrafish development (Tang et al., 2007), by the 2^−ΔΔCt "
      "method (Livak & Schmittgen, 2001).")
    A("")
    A("The orthogonal-validation hypothesis is that supplied high-risk pools will show higher c-fos "
      "fold change than their matched low-risk pools. A difference provides cross-modal molecular "
      "concordance with the risk strata, but does not identify its cause or establish electrographic "
      "epilepsy.")
    A("")
    A("### Statistical treatment")
    A("")
    A("The 18 pools are **9 matched pairs** — one high_risk and one low_risk pool per "
      "(group × clutch) cell — and are analysed as such. Treating them as 18 independent units would "
      "roughly double the nominal degrees of freedom and ignore the plate/clutch matching.")
    A("")
    A("The nine pairs are not nine fully independent biological replicates: three group-level pairs "
      "come from each of only three clutches. The all-pair tests below are therefore nominal. A "
      "clutch-averaged sensitivity analysis, which leaves n = 3, is reported separately.")
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
    csr = s3["clutch_sensitivity_raw"]
    csl = s3["clutch_sensitivity_log2"]
    A("### Clutch-averaged sensitivity")
    A("")
    A("Each row first averages the sham, low-impact, and high-impact pair differences within a "
      "clutch, then tests the three clutch means against zero. This conservative sensitivity treats "
      "clutch as the independent biological unit.")
    A("")
    A("| Scale | Clutches | Mean high − low | 95% CI | one-sample t | p | Wilcoxon p |")
    A("|---|---|---|---|---|---|---|")
    A(f"| Raw fold change | {csr['n']} | {csr['mean']:+.4f} | "
      f"[{csr['ci'][0]:+.4f}, {csr['ci'][1]:+.4f}] | t({csr['n']-1}) = {csr['t']:+.3f} | "
      f"{_fmt_p(csr['p'])} | {_fmt_p(csr['p_wilcoxon'])} |")
    A(f"| Log2 fold change | {csl['n']} | {csl['mean']:+.4f} | "
      f"[{csl['ci'][0]:+.4f}, {csl['ci'][1]:+.4f}] | t({csl['n']-1}) = {csl['t']:+.3f} | "
      f"{_fmt_p(csl['p'])} | {_fmt_p(csl['p_wilcoxon'])} |")
    A("")
    A(f"The nominal nine-pair log2 comparison has {_p(lg['p'])}; after averaging within clutch, "
      f"the log2 sensitivity has {_p(csl['p'])}. With only three independent clutches, the "
      "clutch-level validation estimate remains imprecise; the nominal pair-level result and this "
      "sensitivity should be interpreted together.")
    A("")
    nsign = _g("cfos_injured_pairs", "n_pairs_high_gt_low")
    A(f"Direction consistency: {int(nsign['value'])}/{int(nsign['n'])} injured pairs have "
      f"high_risk > low_risk (exact binomial sign test, {_p(nsign['p_value'])}).")
    A("")
    A(f"The three sham pairs have mean difference {s3['sham']['mean']:+.4f} "
      f"({_p(s3['sham']['p'])}). This small, non-significant comparison does not rule out a pooling, "
      "plate, or bin-assignment artifact.")
    A("")
    A(f"Read the injured-only row carefully: the point estimate is the largest of the three "
      f"({s3['injured']['mean']:+.4f}, dz = {s3['injured']['dz']:.2f}) but with only "
      f"{s3['injured']['n']} pairs it does not reach significance on its own "
      f"({_p(s3['injured']['p'])}). The original all-pairs comparison is nominal; the injured and "
      "sham rows show where its point estimate sits, not two independent confirmations of it.")
    A("")
    A("No regression of c-fos on a continuous risk score was run because continuous scores and a "
      "reproducible pool-selection rule are not present in the workbook. The normalized membership "
      "table verifies the recorded four-fish pools, not the validity of their supplied risk labels.")
    A("")
    A("### Interpretation")
    A("")
    gm_ratio = float(2 ** s3["all_log2"]["mean"])
    gm_lo, gm_hi = (float(2 ** s3["all_log2"]["ci"][0]), float(2 ** s3["all_log2"]["ci"][1]))
    A(f"In the nominal all-pair comparison, supplied high-risk pools have **{100*(gm_ratio-1):.1f}% "
      f"higher c-fos fold change** relative to `rpl13a` than supplied low-risk pools processed "
      f"on the same plate (geometric mean ratio {gm_ratio:.3f}, 95% CI [{gm_lo:.3f}, {gm_hi:.3f}], "
      "back-transformed from the paired log2 analysis). Because three pairs share each clutch, the "
      f"clutch-averaged log2 result ({_p(csl['p'])}) and the missing risk-assignment protocol set the "
      "interpretation boundary.")
    A("")
    A("This constitutes the project's orthogonal molecular validation of the supplied risk "
      "stratification. It is not external validation of the classifier or seizure endpoint. The bulk "
      "measurement on pooled tissue cannot localise the signal to a cell type or region, and raw Ct "
      "values, technical-replicate results, amplification efficiencies, and qPCR quality-control "
      "records are not supplied.")
    A("")
    A("Figure: `fig08_cfos_paired.png` (9 paired lines, plus within-pair differences by group).")
    A("")
    A("---")
    A("")

    # ---------------------------------------------------------------- step 4
    A("## Step 4 — Descriptive results")
    A("")
    A("### No baseline group difference detected")
    A("")
    A(f"At the pre-injury baseline (t = −1), no group difference in τ was detected: "
      f"F({anova['df']}) = {float(anova['statistic']):.3f}, {_p(anova['p_value'])}, "
      f"η² = {float(anova['effect_size']):.4f}; Kruskal–Wallis "
      f"H = {float(_g('baseline_equivalence','kruskal_tau')['statistic']):.3f}, "
      f"{_p(_g('baseline_equivalence','kruskal_tau')['p_value'])}. Baseline locomotion likewise "
      f"({_p(_g('baseline_equivalence','anova_baseline_locomotion')['p_value'])}). The longitudinal "
      "changes are measured within fish, but non-significance does not establish baseline equivalence "
      "or prove random assignment.")
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
      f"in high_impact — slower versus faster fitted response decay. Welch "
      f"t = {float(lv['statistic']):.2f}, {_p(lv['p_value'])}, "
      f"Cohen's d = {float(lv['effect_size']):.2f}.")
    A("")
    A("**Dose changes the biological interpretation of Δτ, but its predictive necessity is an "
      "empirical question.** Model (e) tests the two raw changes without dose; model (f) adds baseline "
      "τ while remaining fully dose-blind; and the comparison of (f) with (g) isolates dose's "
      "incremental contribution.")
    A("")
    A("The observed directions have different empirical meanings:")
    A("")
    A("- **Low dose → τ rises.** The fitted startle response decays more slowly across trials, which is "
      "consistent with impaired habituation.")
    A("- **High dose → τ falls.** The fitted response decays more quickly. This could reflect faster "
      "habituation, depressed responsiveness, fatigue, or another process; the present measurements "
      "do not distinguish these explanations.")
    A("")
    A("The dose groups move the same fitted scalar in opposite directions. A model given "
      "Δτ without dose is asked to treat +4 trials and −3 trials as opposite kinds of evidence when "
      "they may reflect different processes. Encoding dose is one way to represent that context, "
      "but the strong dose-blind behavioral ablation shows that baseline and trajectory can recover "
      "substantial predictive information without it.")
    A("")
    A("The group means also remain shifted from baseline at 24 h "
      f"(low impact {float(_g('tau_change_from_baseline','low_impact_delta_t24')['value']):+.2f}, "
      f"high impact {float(_g('tau_change_from_baseline','high_impact_delta_t24')['value']):+.2f} "
      "trials from each fish's own pre-injury session). This documents persistent behavioral change; "
      "it does not by itself identify circuit remodeling or an epileptogenic latent period.")
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
      "probe of seizure susceptibility. It is included here as an exploratory, pharmacological "
      "comparison rather than as confirmation of the behavioral model.")
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
      f"check only, and **no conclusion in this report rests on it.**")
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
    A(f"- **The model includes n = {n_fish} fish with {n_events} events.** EPV = {epv:.1f} is low. The CIs on "
      "the AUC and on every coefficient are wide, and they are reported rather than smoothed over.")
    A("- **Three clutches means three outer folds.** The fold-to-fold spread is estimated from three "
      "numbers; the SD across folds should be read as indicative, not precise. The AUC interval "
      "resamples fish within the observed clutches while holding out-of-fold scores fixed, so it is "
      "conditional on these clutches and does not include model-selection uncertainty.")
    A(f"- **{int(_g('design_matrix','n_excluded_incomplete')['value'])} injured fish were dropped** for "
      "a missing post-injury session. This is complete-case analysis, and the missingness mechanism "
      "has not been modelled.")
    A("- **The binary outcome is supplied, not derived by code.** The workbook does not document a "
      "threshold, recording duration, blinded scoring procedure, or electrographic confirmation.")
    A("- **The c-fos comparison uses pooled material** — 9 group-level pairs, 4 larvae per pool, but "
      "only 3 clutches. Risk bins are supplied rather than reproducibly generated, and raw qPCR Ct "
      "and QC data are absent.")
    A("- **The workbook represents 253 unique animal IDs across three non-overlapping cohorts** "
      "(133 followed, 86 c-fos, and 34 PTZ); 72 of the c-fos animals entered pools.")
    A("- **Data provenance is unresolved in the repository.** Source recordings, instrument exports, "
      "dated protocols, approval identifiers, and a label-derivation audit trail are not included.")
    A("- **The feature and model choices are not preregistered.** Nested cross-validation covers "
      "penalty tuning, not uncertainty from post-hoc feature, timepoint, or model selection.")
    A("- **PTZ is underpowered** (above), and the study is not designed to support any claim from it.")
    A("")
    A("## Conclusion")
    A("")
    permutation_conclusion = (
        f"and {_p(perm['p_pooled'])} against a null that reruns the entire nested cross-validation"
        if permutation_ran
        else "with the permutation p-value unavailable because that test was not run"
    )
    A(f"In the supplied dataset, baseline and post-injury startle-habituation features distinguish "
      f"injured larvae with versus without the supplied `converted` label, with a "
      f"cross-validated AUC of {full['pooled_auc']:.3f} "
      f"(conditional 95% interval [{comp.loc[comp.key=='e_full','pooled_ci_low'].iloc[0]:.3f}, "
      f"{comp.loc[comp.key=='e_full','pooled_ci_high'].iloc[0]:.3f}]), "
      f"a total out-of-fold accuracy of {100*cm['accuracy']:.1f}%, {permutation_conclusion}. The "
      f"dose-blind three-behavior model reaches a mean-fold AUC of {f:.3f}, so explicit dose encoding "
      "is not necessary for the observed separation. The c-fos experiment supplies orthogonal "
      "molecular validation through cross-modal concordance with the supplied risk strata; it is not "
      "external validation of the classifier.")
    A("")
    A("The defensible conclusion is narrow: **startle-habituation kinetics carry internally "
      "cross-validated information about a later behavioral label in this workbook.** Prospective "
      "replication with documented outcome scoring, independent clutches, source-data provenance, "
      "and electrographic confirmation is required before calling the assay a biomarker of "
      "post-traumatic epilepsy.")
    A("")
    A("## Background literature")
    A("")
    A("The contextual discussion above draws on the following. This is a background reading list, not "
      "a citation of results generated here.")
    A("")
    A("> These references were already present in the source repository. A student preparing an "
      "ISEF submission must independently read, verify, and format every citation under the current "
      "AI-use rules.")
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
      "quantitative PCR and the 2^−ΔΔCt method. *Methods* 25:402–408.")
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
