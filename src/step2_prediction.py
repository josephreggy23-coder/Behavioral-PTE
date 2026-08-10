"""STEP 2 driver: nested model comparison, permutation test, and split sensitivity."""
from __future__ import annotations

import numpy as np
import pandas as pd
import statsmodels.api as smapi
from statsmodels.stats.outliers_influence import variance_inflation_factor

from . import config, features, modeling, plotting, statsbook as sb


def _assumption_checks(model_df: pd.DataFrame) -> None:
    sb.banner("STEP 2b -- assumption checks for the prediction model")

    X = model_df[features.PREDICTORS].copy()
    y = model_df["converted"].to_numpy(int)

    # variance inflation
    Xc = smapi.add_constant(X.astype(float))
    vifs = {c: float(variance_inflation_factor(Xc.to_numpy(), i))
            for i, c in enumerate(Xc.columns) if c != "const"}
    for c, v in vifs.items():
        sb.record("2", "assumptions", f"vif_{c}", v, n=len(X))
        print(f"  VIF {c:<10} = {v:.2f}")
    sb.check("all VIF < 5", max(vifs.values()) < 5, f"max VIF = {max(vifs.values()):.2f}")

    # complete/quasi separation: an unpenalised fit that fails to converge or
    # produces huge coefficients would signal separation
    with np.errstate(all="ignore"):
        try:
            res = smapi.Logit(y, smapi.add_constant(X.astype(float))).fit(disp=0)
            maxabs = float(np.max(np.abs(res.params)))
            sb.record("2", "assumptions", "max_abs_unpenalised_coef", maxabs, n=len(X),
                      notes="separation diagnostic; |coef| > 10 suggests quasi-separation")
            sb.check("no (quasi-)separation, max |unpenalised coef| < 10", maxabs < 10,
                     f"max |b| = {maxabs:.2f}")
            unpen = pd.DataFrame(
                {"predictor": res.params.index, "coef": res.params.values,
                 "se": res.bse.values, "z": res.tvalues.values, "p": res.pvalues.values,
                 "ci_low": res.conf_int()[0].values, "ci_high": res.conf_int()[1].values}
            )
            unpen.to_csv(config.TABLES / "step2_unpenalised_logit.csv", index=False)
            for r in unpen.itertuples():
                sb.record("2", "unpenalised_logit", f"coef_{r.predictor}", float(r.coef),
                          ci_low=float(r.ci_low), ci_high=float(r.ci_high), n=len(X),
                          test="Wald z (statsmodels Logit, unpenalised)",
                          statistic=float(r.z), p_value=float(r.p),
                          effect_size_name="odds ratio (raw units)",
                          effect_size=float(np.exp(r.coef)),
                          notes="reference only; the reported model is ridge-penalised")
        except Exception as exc:  # pragma: no cover
            print(f"  [FLAG] unpenalised Logit failed: {exc}")

    # linearity of the logit for continuous predictors (Box-Tidwell style check
    # via a quartile-binned event rate is more readable at this n)
    for col in ["pre_tau", "dtau_0.5", "dtau_24"]:
        q = pd.qcut(model_df[col], 4, labels=False, duplicates="drop")
        rates = model_df.groupby(q, observed=True)["converted"].mean()
        mono = bool(rates.is_monotonic_increasing or rates.is_monotonic_decreasing)
        print(f"  event rate by quartile of {col:<9}: "
              + " ".join(f"{v:.2f}" for v in rates.values)
              + ("  (monotone)" if mono else "  (non-monotone)"))
        sb.record("2", "assumptions", f"quartile_event_rates_{col}",
                  ";".join(f"{v:.3f}" for v in rates.values), n=len(X),
                  notes="monotone" if mono else "non-monotone (linear logit is an approximation)")

    # class balance per outer fold
    per = model_df.groupby("clutch", observed=True)["converted"].agg(["size", "sum"])
    print("  outer folds (leave-one-clutch-out):")
    for c, r in per.iterrows():
        print(f"    {c}: n = {int(r['size'])}, events = {int(r['sum'])}")


def run(model_df: pd.DataFrame, *, n_perm: int = config.N_PERMUTATIONS) -> dict:
    _assumption_checks(model_df)

    y = model_df["converted"].to_numpy(int)
    groups = model_df["clutch"].to_numpy()
    X = model_df.copy()

    # ---------------------------------------------------------------- table
    sb.banner("STEP 2c -- nested comparison table (leave-one-clutch-out, nested C tuning)")
    rows, folds_long, roc_data = [], [], []
    results: dict[str, dict] = {}

    for key, spec in features.MODEL_SPECS.items():
        res = modeling.nested_cv(X, y, groups, spec["features"], seed=config.SEED)
        results[key] = res
        lo, hi = modeling.bootstrap_auc_ci(y, res["oof"], groups, seed=config.SEED)
        rows.append(
            {
                "key": key,
                "label": spec["label"],
                "question": spec["question"],
                "n_predictors": len(spec["features"]),
                "predictors": " + ".join(spec["features"]),
                "mean_fold_auc": res["mean_fold_auc"],
                "sd_fold_auc": res["sd_fold_auc"],
                "min_fold_auc": res["min_fold_auc"],
                "max_fold_auc": res["max_fold_auc"],
                "pooled_oof_auc": res["pooled_auc"],
                "pooled_ci_low": lo,
                "pooled_ci_high": hi,
                "brier": res["brier"],
                "selected_C": ";".join(f"{c:g}" for c in res["folds"]["selected_C"]),
            }
        )
        f = res["folds"].copy()
        f["key"] = key
        f["label"] = spec["label"]
        folds_long.append(f)
        fpr, tpr = modeling.roc_points(y, res["oof"])
        roc_data.append({"key": key, "label": spec["label"], "fpr": fpr, "tpr": tpr,
                         "auc": res["pooled_auc"]})

        sb.record("2", "nested_comparison", f"{key}_mean_fold_auc", res["mean_fold_auc"],
                  n=len(y), test="nested CV (GroupKFold on clutch)",
                  notes=f"{spec['label']}; folds: " +
                        ", ".join(f"{r.held_out}={r.auc:.3f}" for r in res["folds"].itertuples()))
        sb.record("2", "nested_comparison", f"{key}_sd_fold_auc", res["sd_fold_auc"], n=len(y))
        sb.record("2", "nested_comparison", f"{key}_pooled_oof_auc", res["pooled_auc"],
                  ci_low=lo, ci_high=hi, n=len(y),
                  notes=("within-clutch stratified bootstrap conditional on the observed "
                         "clutches; percentile CI"))
        sb.record("2", "nested_comparison", f"{key}_brier", res["brier"], n=len(y))

    comparison = pd.DataFrame(rows)
    folds_long = pd.concat(folds_long, ignore_index=True)
    comparison.to_csv(config.TABLES / "step2_nested_comparison.csv", index=False)
    folds_long.to_csv(config.TABLES / "step2_fold_detail.csv", index=False)

    print(f"\n  {'model':<32}{'mean fold AUC':>15}{'SD':>7}{'range':>16}"
          f"{'pooled OOF AUC [conditional 95% interval]':>45}")
    for r in comparison.itertuples():
        print(f"  {r.label:<32}{r.mean_fold_auc:>15.3f}{r.sd_fold_auc:>7.3f}"
              f"{f'{r.min_fold_auc:.3f}-{r.max_fold_auc:.3f}':>16}"
              f"{f'{r.pooled_oof_auc:.3f} [{r.pooled_ci_low:.3f}, {r.pooled_ci_high:.3f}]':>28}")
    print("\n  per-fold detail (held-out clutch):")
    for r in folds_long.itertuples():
        print(f"    {r.label:<32} {r.held_out}  n={r.n_test:2d} events={r.n_events_test:2d} "
              f"C={r.selected_C:<6g} AUC={r.auc:.3f}")

    full = results["e_full"]

    # ------------------------------------------ cluster-ignoring split comparison
    sb.banner("STEP 2d -- random split vs clutch-held-out split")
    naive = modeling.nested_cv(X, y, groups, features.MODEL_SPECS["e_full"]["features"],
                               outer_splits=config.N_RANDOM_SPLIT_FOLDS,
                               seed=config.SEED, random_outer=True)
    inflation = naive["mean_fold_auc"] - full["mean_fold_auc"]
    sb.record("2", "leakage", "random_split_mean_fold_auc", naive["mean_fold_auc"],
              n=len(y), test=f"{config.N_RANDOM_SPLIT_FOLDS}-fold random StratifiedKFold",
              notes="ignores clutch; reported only to quantify inflation")
    sb.record("2", "leakage", "random_split_pooled_auc", naive["pooled_auc"], n=len(y))
    sb.record("2", "leakage", "clutch_split_mean_fold_auc", full["mean_fold_auc"], n=len(y),
              test="GroupKFold on clutch")
    sb.record("2", "leakage", "inflation_auc", inflation, n=len(y),
              notes="random-split AUC minus clutch-split AUC")
    print(f"  random split (ignores clutch): mean fold AUC = {naive['mean_fold_auc']:.3f} "
          f"(pooled {naive['pooled_auc']:.3f})")
    print(f"  clutch-held-out split:         mean fold AUC = {full['mean_fold_auc']:.3f} "
          f"(pooled {full['pooled_auc']:.3f})")
    print(f"  inflation from ignoring clutch = {inflation:+.3f} AUC")

    # ---------------------------------------------------- permutation test
    if n_perm < 0:
        raise ValueError("n_perm must be non-negative")
    if n_perm == 0:
        sb.banner("STEP 2e -- permutation test skipped")
        perm = modeling.skipped_permutation_result()
        sb.record("2", "permutation", "status", "skipped", n=0,
                  notes="disabled by --skip-permutation or --permutations 0; p-value unavailable")
        print("  no label shuffles run; permutation p-values are unavailable")
    else:
        sb.banner(f"STEP 2e -- permutation test ({n_perm} full nested-CV reruns)")
        perm = modeling.permutation_test(
            X, y, groups, features.MODEL_SPECS["e_full"]["features"],
            full["pooled_auc"], full["mean_fold_auc"], n_perm=n_perm, seed=config.SEED,
        )
        sb.record("2", "permutation", "observed_pooled_auc", full["pooled_auc"], n=len(y))
        sb.record("2", "permutation", "null_mean_auc", perm["null_mean"], n=n_perm)
        sb.record("2", "permutation", "null_sd_auc", perm["null_sd"], n=n_perm)
        sb.record("2", "permutation", "null_95th_percentile", perm["null_q95"], n=n_perm)
        sb.record("2", "permutation", "p_value_pooled_auc", perm["p_pooled"], n=len(y),
                  test=f"permutation test, {n_perm} label shuffles within clutch",
                  statistic=full["pooled_auc"], p_value=perm["p_pooled"],
                  notes="p = (1 + #{null >= observed}) / (n_perm + 1)")
        sb.record("2", "permutation", "p_value_mean_fold_auc", perm["p_mean_fold"], n=len(y),
                  test=f"permutation test, {n_perm} shuffles", statistic=full["mean_fold_auc"],
                  p_value=perm["p_mean_fold"])
        sb.record("2", "permutation", "percentile_of_observed_in_null",
                  perm["percentile_of_observed"], n=n_perm)
        print(f"  null: mean {perm['null_mean']:.3f}, SD {perm['null_sd']:.3f}, "
              f"95th pct {perm['null_q95']:.3f}")
        print(f"  observed pooled AUC {full['pooled_auc']:.3f} sits at the "
              f"{perm['percentile_of_observed']:.1f}th percentile of the null; "
              f"p = {perm['p_pooled']:.4f}")
        print(f"  (mean-fold-AUC statistic: observed {full['mean_fold_auc']:.3f}, "
              f"p = {perm['p_mean_fold']:.4f})")
        sb.check("null centred near 0.5 (permutation machinery is unbiased)",
                 abs(perm["null_mean"] - 0.5) < 0.05,
                 f"null mean = {perm['null_mean']:.3f}")

    # --------------------------------------------- final model & diagnostics
    sb.banner("STEP 2f -- full model: coefficients, confusion matrix, calibration")
    final_pipe, best_C = modeling.fit_final_model(
        X, y, groups, features.MODEL_SPECS["e_full"]["features"], seed=config.SEED
    )
    coefs = modeling.coefficient_table(
        final_pipe, X, y, groups, features.MODEL_SPECS["e_full"]["features"], best_C,
        seed=config.SEED,
    )
    coefs.to_csv(config.TABLES / "step2_full_model_coefficients.csv", index=False)
    sb.record("2", "full_model", "selected_C_full_data", best_C, n=len(y),
              notes="chosen by clutch-held-out CV on the full set; the reported AUC does NOT use this")
    print(f"  selected C (full data refit) = {best_C:g}; "
          f"intercept = {coefs.attrs['intercept']:+.3f}")
    for r in coefs.itertuples():
        sb.record("2", "full_model", f"coef_{r.predictor}", r.coef_standardised,
                  ci_low=r.ci_low, ci_high=r.ci_high, n=len(y),
                  test=("ridge logistic (within-clutch stratified bootstrap conditional on the "
                        "observed clutches)"),
                  effect_size_name="odds ratio per SD", effect_size=r.odds_ratio_per_SD,
                  notes=f"OR 95% CI [{r.or_ci_low:.3f}, {r.or_ci_high:.3f}]; "
                        f"{'excludes' if r.excludes_zero else 'includes'} zero")
        print(f"    {r.predictor:<10} b = {r.coef_standardised:+.3f} "
              f"[{r.ci_low:+.3f}, {r.ci_high:+.3f}]  OR/SD = {r.odds_ratio_per_SD:.2f} "
              f"[{r.or_ci_low:.2f}, {r.or_ci_high:.2f}]")

    cm = modeling.confusion_at(y, full["oof"], 0.5)
    for k in ["accuracy", "balanced_accuracy", "sensitivity", "specificity", "ppv", "npv"]:
        sb.record("2", "confusion_matrix", k, cm[k], n=len(y),
                  notes="out-of-fold predictions, threshold 0.50")
    sb.record("2", "confusion_matrix", "tn_fp_fn_tp",
              f"{cm['tn']};{cm['fp']};{cm['fn']};{cm['tp']}", n=len(y))
    sb.record("2", "calibration", "brier_score", full["brier"], n=len(y),
              notes="out-of-fold; 0.25 = uninformative at 50% prevalence")
    print(f"  confusion (thr 0.50): TN {cm['tn']} FP {cm['fp']} FN {cm['fn']} TP {cm['tp']}")
    print(f"  accuracy {cm['accuracy']:.3f}, balanced accuracy {cm['balanced_accuracy']:.3f}, "
          f"sens {cm['sensitivity']:.3f}, spec {cm['specificity']:.3f}")
    print(f"  Brier score {full['brier']:.3f}")

    # secondary operating point: Youden-optimal threshold
    thr_j = modeling.youden_threshold(y, full["oof"])
    cm_j = modeling.confusion_at(y, full["oof"], thr_j)
    sb.record("2", "confusion_matrix_youden", "threshold", thr_j, n=len(y),
              notes="Youden-optimal operating point, selected on the out-of-fold predictions; "
                    "mildly optimistic, reported as an operating-point illustration only")
    for k in ["accuracy", "balanced_accuracy", "sensitivity", "specificity", "ppv", "npv"]:
        sb.record("2", "confusion_matrix_youden", k, cm_j[k], n=len(y),
                  notes=f"threshold {thr_j:.3f} (Youden)")
    sb.record("2", "confusion_matrix_youden", "tn_fp_fn_tp",
              f"{cm_j['tn']};{cm_j['fp']};{cm_j['fn']};{cm_j['tp']}", n=len(y))
    print(f"  Youden-optimal threshold {thr_j:.3f}: accuracy {cm_j['accuracy']:.3f}, "
          f"balanced accuracy {cm_j['balanced_accuracy']:.3f}, sens {cm_j['sensitivity']:.3f}, "
          f"spec {cm_j['specificity']:.3f}")

    # ---------------------------------------- sensitivity: global z-scoring
    gz = features.add_global_within_dose_z(model_df)
    gz_feats = ["dose", "pre_tau", "z_dtau_0.5", "z_dtau_24"]
    gz_res = modeling.nested_cv(gz, y, groups, gz_feats, seed=config.SEED)
    sb.record("2", "sensitivity", "global_z_mean_fold_auc", gz_res["mean_fold_auc"], n=len(y),
              notes="within-dose z fitted on the WHOLE dataset (mildly leaky); "
                    "reported to show the fold-safe scaler costs nothing")
    print(f"  sensitivity -- whole-dataset within-dose z: mean fold AUC "
          f"{gz_res['mean_fold_auc']:.3f} vs fold-safe {full['mean_fold_auc']:.3f}")

    # ----------------------------------------------------------- figures
    sb.banner("STEP 2g -- figures")
    plotting.fig_roc_comparison(roc_data)
    if perm["performed"]:
        plotting.fig_permutation_null(
            perm["null_pooled"], full["pooled_auc"], perm["p_pooled"], n_perm
        )
    plotting.fig_confusion_and_calibration(y, full["oof"], cm, full["brier"])
    plotting.fig_coefficients(coefs)
    plotting.fig_fold_auc(comparison, folds_long)

    oof_df = model_df[["fish_id", "group", "clutch", "converted"]].copy()
    oof_df["oof_probability"] = full["oof"]
    oof_df.to_csv(config.TABLES / "step2_out_of_fold_predictions.csv", index=False)

    return {
        "comparison": comparison,
        "folds": folds_long,
        "full": full,
        "perm": perm,
        "naive": naive,
        "coefs": coefs,
        "cm": cm,
        "best_C": best_C,
        "inflation": inflation,
    }
