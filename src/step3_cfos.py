"""STEP 3 -- molecular validation of the conversion phenotype (paired c-fos pools).

WHAT CHANGED FROM THE EARLIER DESIGN, AND WHY IT MATTERS
--------------------------------------------------------
In this dataset the qPCR pools are labelled ``converter`` / ``non_converter``,
i.e. by the larva's *realised outcome*, not by a *predicted* risk score. Two
consequences follow and both are stated in the report rather than glossed:

1. This is no longer an independent validation of the prediction model. It is a
   validation of the **outcome variable the model predicts**: it asks whether
   larvae scored as converted carry a molecular signature of elevated network
   activity, or whether "converted" is just a behavioural scoring threshold.
   That is a weaker claim about the model and a stronger claim about the label.
2. The pools are drawn from the same ``zf_*`` cohort that trains the model, so
   the two analyses share animals. Sample independence is not claimed.

The pool counts are also unbalanced by construction, because conversion is rare
in sham:

    sham         0 converter pools, 2 non_converter pools per clutch
    low_impact   1 converter pool,  2 non_converter pools per clutch
    high_impact  2 converter pools, 2 non_converter pools per clutch

So a balanced 9-pair design is impossible. The pairing that the data actually
support is **6 matched pairs**: the (group x clutch) cells of the two injured
groups, each contributing one converter value and one non-converter value, with
technical replicate pools within a cell averaged first. Sham cells carry no
converter pool at all and therefore cannot be paired; they are used instead as
an unpaired reference level for the non-converter baseline, which tests the
complementary question -- does c-fos track conversion, or merely injury?
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from . import config, io_data, plotting, statsbook as sb

VALUE_COL = "cfos_fold_change"
CONVERTER, NON_CONVERTER = "converter", "non_converter"


def build_pairs(pools: pd.DataFrame | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Collapse replicate pools within each (group x clutch, pool_type) cell.

    Returns (cell_table, pair_table). ``pair_table`` holds only the cells where
    both pool types are present, which is the paired analysis set.
    """
    if pools is None:
        pools = io_data.cfos_pools()

    cells = (
        pools.groupby(["group", "clutch", "pool_type"], observed=True)
        .agg(
            fold_change=(VALUE_COL, "mean"),
            ddct=("delta_ddct", "mean"),
            n_pools=(VALUE_COL, "size"),
            n_larvae=("n_larvae_in_pool", "sum"),
        )
        .reset_index()
    )

    wide = cells.pivot_table(
        index=["group", "clutch"], columns="pool_type", values="fold_change", observed=True
    ).reset_index()
    counts = cells.pivot_table(
        index=["group", "clutch"], columns="pool_type", values="n_pools", observed=True
    ).reset_index().rename(columns={CONVERTER: "n_conv_pools", NON_CONVERTER: "n_nonconv_pools"})
    wide = wide.merge(counts, on=["group", "clutch"], how="left")

    have_both = wide[CONVERTER].notna() & wide[NON_CONVERTER].notna()
    pairs = wide[have_both].copy()
    pairs["diff"] = pairs[CONVERTER] - pairs[NON_CONVERTER]
    pairs["log2_conv"] = np.log2(pairs[CONVERTER])
    pairs["log2_nonconv"] = np.log2(pairs[NON_CONVERTER])
    pairs["log2_diff"] = pairs["log2_conv"] - pairs["log2_nonconv"]
    pairs["ratio"] = pairs[CONVERTER] / pairs[NON_CONVERTER]

    return cells, pairs.sort_values(["group", "clutch"]).reset_index(drop=True)


def _paired_test(d: pd.Series, label: str, scope: str, value_name: str) -> dict:
    n = len(d)
    mean = float(d.mean())
    sd = float(d.std(ddof=1))
    t, p = stats.ttest_1samp(d, 0.0)  # paired t == one-sample t on the differences
    ci = stats.t.interval(0.95, n - 1, loc=mean, scale=sd / np.sqrt(n))
    dz = mean / sd if sd > 0 else np.nan

    try:
        w, pw = stats.wilcoxon(d, zero_method="wilcox", alternative="two-sided")
    except ValueError:
        w, pw = np.nan, np.nan
    sh_w, sh_p = stats.shapiro(d) if n >= 3 else (np.nan, np.nan)

    sb.record("3", scope, "n_pairs", n, n=n, notes=label)
    sb.record("3", scope, f"mean_{value_name}_difference", mean, ci_low=float(ci[0]),
              ci_high=float(ci[1]), n=n, test="paired t-test (one-sample on differences)",
              statistic=float(t), df=n - 1, p_value=float(p),
              effect_size_name="Cohen's dz", effect_size=float(dz), notes=label)
    sb.record("3", scope, f"wilcoxon_{value_name}", mean, n=n,
              test="Wilcoxon signed-rank (robustness check)", statistic=float(w),
              p_value=float(pw), notes=label)
    sb.record("3", scope, f"shapiro_{value_name}_differences", float(sh_w), n=n,
              test="Shapiro-Wilk on within-pair differences", statistic=float(sh_w),
              p_value=float(sh_p), notes="normality assumption for the paired t-test")

    print(f"  {label} (n = {n} pairs)")
    print(f"    mean difference = {mean:+.4f}  95% CI [{ci[0]:+.4f}, {ci[1]:+.4f}]")
    print(f"    paired t({n-1}) = {t:+.3f}, p = {p:.4g}, Cohen's dz = {dz:+.3f}")
    print(f"    Wilcoxon W = {w}, p = {pw:.4g}")
    sb.check(f"{label}: differences normal (Shapiro p > 0.05)",
             bool(np.isfinite(sh_p) and sh_p > 0.05), f"W = {sh_w:.3f}, p = {sh_p:.3g}")
    return {"n": n, "mean": mean, "ci": (float(ci[0]), float(ci[1])), "t": float(t),
            "p": float(p), "dz": float(dz), "w": float(w), "p_wilcoxon": float(pw),
            "shapiro_p": float(sh_p)}


def _non_converter_baseline(cells: pd.DataFrame) -> dict:
    """Does c-fos track conversion, or merely injury?

    Compares non-converter pools of injured groups against non-converter pools of
    sham. If injury alone raised c-fos, injured non-converters would sit above
    sham non-converters. If the signal is specific to conversion, they will not.
    """
    nc = cells[cells["pool_type"] == NON_CONVERTER]
    sham = nc[nc["group"] == "sham"]["fold_change"]
    inj = nc[nc["group"].isin(config.INJURED)]["fold_change"]

    t, p = stats.ttest_ind(inj, sham, equal_var=False)
    psd = np.sqrt((inj.var(ddof=1) + sham.var(ddof=1)) / 2)
    d = (inj.mean() - sham.mean()) / psd if psd > 0 else np.nan
    diff = float(inj.mean() - sham.mean())
    se = np.sqrt(inj.var(ddof=1) / len(inj) + sham.var(ddof=1) / len(sham))
    dof = len(inj) + len(sham) - 2
    ci = stats.t.interval(0.95, dof, loc=diff, scale=se)

    sb.record("3", "cfos_noncoverter_baseline", "injured_minus_sham_noncoverters", diff,
              ci_low=float(ci[0]), ci_high=float(ci[1]), n=len(inj) + len(sham),
              test="Welch t-test, non-converter pools only (injured vs sham)",
              statistic=float(t), p_value=float(p),
              effect_size_name="Cohen's d", effect_size=float(d),
              notes="a NULL result here is the desired outcome: it shows c-fos tracks "
                    "conversion rather than injury exposure")
    print(f"  Non-converter baseline (injury without conversion): injured "
          f"{inj.mean():.3f} vs sham {sham.mean():.3f}, Welch t = {t:+.3f}, p = {p:.4g}, "
          f"d = {d:+.3f}")
    sb.check("c-fos does not rise with injury alone (non-converters: injured == sham)",
             bool(p > 0.05), f"p = {p:.4f}")
    return {"n_inj": len(inj), "n_sham": len(sham), "mean_inj": float(inj.mean()),
            "mean_sham": float(sham.mean()), "diff": diff,
            "ci": (float(ci[0]), float(ci[1])), "t": float(t), "p": float(p), "d": float(d)}


def run(pools: pd.DataFrame | None = None) -> dict:
    sb.banner("STEP 3 -- c-fos pools: converter vs non-converter, matched within group x clutch")
    raw = io_data.cfos_pools() if pools is None else pools
    cells, pairs = build_pairs(raw)
    raw.to_csv(config.TABLES / "step3_cfos_pools_raw.csv", index=False)
    cells.to_csv(config.TABLES / "step3_cfos_cells.csv", index=False)
    pairs.to_csv(config.TABLES / "step3_cfos_pairs.csv", index=False)

    n_pools = len(raw)
    print(f"  {n_pools} pools, {int(raw['n_larvae_in_pool'].sum())} larvae "
          f"({raw['n_larvae_in_pool'].unique().tolist()} per pool)")
    print(f"  pool counts by cell:")
    tab = raw.pivot_table(index=["group", "clutch"], columns="pool_type",
                          values="pool_id", aggfunc="size", observed=True).fillna(0).astype(int)
    print(tab.to_string())
    print(f"  -> {len(pairs)} cells contain BOTH pool types and can be paired; "
          f"{len(cells[cells['pool_type'] == CONVERTER])} converter cells exist in total")

    sb.record("3", "design", "n_pools", n_pools,
              notes="labelled by realised outcome (converter / non_converter)")
    sb.record("3", "design", "n_paired_cells", len(pairs),
              notes="group x clutch cells containing both pool types; replicate pools averaged")
    sb.record("3", "design", "n_sham_converter_pools", 0 if raw[
        (raw["group"] == "sham") & (raw["pool_type"] == CONVERTER)].empty else int(
        len(raw[(raw["group"] == "sham") & (raw["pool_type"] == CONVERTER)])),
              notes="sham conversion is too rare to form a pool, so sham cannot be paired")

    sb.check("every paired cell has both pool types", bool(len(pairs) == 6), f"{len(pairs)} pairs")

    out = {"pairs": pairs, "cells": cells, "raw": raw}
    out["all"] = _paired_test(pairs["diff"], f"All {len(pairs)} injured pairs", "cfos_all_pairs",
                              "fold_change")
    out["all_log2"] = _paired_test(pairs["log2_diff"], f"All {len(pairs)} injured pairs, log2 scale",
                                   "cfos_all_pairs_log2", "log2_fold_change")

    for grp in config.INJURED:
        g = pairs[pairs["group"] == grp]
        sb.record("3", f"cfos_{grp}", "mean_fold_change_difference", float(g["diff"].mean()),
                  n=len(g), notes=f"{grp}: descriptive only, {len(g)} pairs is too few to test")
        print(f"  {config.GROUP_LABELS[grp]}: mean converter - non-converter = "
              f"{g['diff'].mean():+.4f} ({len(g)} pairs, descriptive only)")

    out["baseline"] = _non_converter_baseline(cells)

    n_pos = int((pairs["diff"] > 0).sum())
    p_sign = float(stats.binomtest(n_pos, len(pairs), 0.5).pvalue)
    sb.record("3", "cfos_all_pairs", "n_pairs_converter_gt_nonconverter", n_pos, n=len(pairs),
              test="exact binomial sign test", p_value=p_sign,
              notes="direction consistency across paired cells")
    print(f"  direction: {n_pos}/{len(pairs)} pairs have converter > non-converter "
          f"(sign test p = {p_sign:.4g})")
    out["n_positive"] = n_pos
    out["p_sign"] = p_sign

    print("  NOTE: no pooled continuous risk-score regression is run; pools are labelled by "
          "realised outcome, and the risk score is not comparable across dose groups.")
    sb.record("3", "design", "pooled_risk_regression_run", "no",
              notes="deliberately omitted: pools are outcome-labelled and the risk score is "
                    "not comparable between dose groups")
    sb.record("3", "design", "shares_animals_with_model_cohort", "yes",
              notes="single zf_* cohort; this validates the OUTCOME LABEL, not the model's "
                    "predictions on unseen animals")

    plotting.fig_cfos_paired(pairs, cells)
    return out
