"""STEP 3 -- orthogonal molecular validation on the paired c-fos pools.

The 18 pools are NOT 18 independent observations.  They are 9 matched pairs:
one high_risk and one low_risk pool per (group x clutch) cell, processed on the
same plate with the same reference gene.  The correct unit of analysis is the
within-pair difference, so the primary test is a paired t-test across the 9
pairs, with Wilcoxon signed-rank as a distribution-free robustness check.

A pooled regression of c-fos on a continuous risk score across all 18 pools is
deliberately NOT run: the risk score is not on a comparable scale between dose
groups (low and high dose move tau in opposite directions), so pooling collapses
the very contrast the pairing was designed to isolate.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from . import config, io_data, plotting, statsbook as sb

VALUE_COL = "cfos_fold_change"


def build_pairs(pools: pd.DataFrame | None = None) -> pd.DataFrame:
    if pools is None:
        pools = io_data.cfos_pools()
    wide = pools.pivot_table(
        index=["group", "clutch"], columns="risk_pool", values=VALUE_COL, observed=True
    ).reset_index()
    ddct = pools.pivot_table(
        index=["group", "clutch"], columns="risk_pool", values="delta_ddct", observed=True
    ).reset_index().rename(columns={"high_risk": "ddct_high", "low_risk": "ddct_low"})
    wide = wide.merge(ddct, on=["group", "clutch"])
    wide["diff"] = wide["high_risk"] - wide["low_risk"]
    wide["log2_high"] = np.log2(wide["high_risk"])
    wide["log2_low"] = np.log2(wide["low_risk"])
    wide["log2_diff"] = wide["log2_high"] - wide["log2_low"]
    wide["ratio"] = wide["high_risk"] / wide["low_risk"]
    return wide.sort_values(["group", "clutch"]).reset_index(drop=True)


def _paired_test(d: pd.Series, label: str, scope: str, value_name: str) -> dict:
    n = len(d)
    mean = float(d.mean())
    sd = float(d.std(ddof=1))
    # a paired t-test is a one-sample t-test on the within-pair differences
    t, p = stats.ttest_1samp(d, 0.0)
    ci = stats.t.interval(0.95, n - 1, loc=mean, scale=sd / np.sqrt(n))
    dz = mean / sd if sd > 0 else np.nan

    if n >= 3:
        try:
            w, pw = stats.wilcoxon(d, zero_method="wilcox", alternative="two-sided")
        except ValueError:
            w, pw = np.nan, np.nan
    else:
        w, pw = np.nan, np.nan

    sh_w, sh_p = (stats.shapiro(d) if n >= 3 else (np.nan, np.nan))

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


def run(pools: pd.DataFrame | None = None) -> dict:
    sb.banner("STEP 3 -- paired c-fos pools (9 matched high_risk / low_risk pairs)")
    pairs = build_pairs(pools)
    pairs.to_csv(config.TABLES / "step3_cfos_pairs.csv", index=False)

    n_cells = pairs.groupby(["group", "clutch"], observed=True).size()
    sb.check("exactly one pair per group x clutch cell",
             bool((n_cells == 1).all()) and len(pairs) == 9, f"{len(pairs)} pairs")
    sb.record("3", "design", "n_pools", 9 * 2, notes="analysed as 9 pairs, NOT 18 independent units")
    print(f"  pools: {len(pairs)*2} -> {len(pairs)} pairs (group x clutch)")
    print(f"  larvae per pool: {io_data.cfos_pools()['n_larvae_in_pool'].unique().tolist()}")

    out = {"pairs": pairs}
    out["all"] = _paired_test(pairs["diff"], "All 9 pairs", "cfos_all_pairs", "fold_change")
    out["all_log2"] = _paired_test(pairs["log2_diff"], "All 9 pairs, log2 scale",
                                   "cfos_all_pairs_log2", "log2_fold_change")

    inj = pairs[pairs["group"].isin(config.INJURED)]
    out["injured"] = _paired_test(inj["diff"], "Injured pairs only (low + high impact)",
                                  "cfos_injured_pairs", "fold_change")

    sham = pairs[pairs["group"] == "sham"]
    out["sham"] = _paired_test(sham["diff"], "Sham pairs only (negative control)",
                               "cfos_sham_pairs", "fold_change")

    for grp in config.INJURED:
        g = pairs[pairs["group"] == grp]
        sb.record("3", f"cfos_{grp}", "mean_fold_change_difference", float(g["diff"].mean()),
                  n=len(g), notes=f"{grp}: descriptive only, 3 pairs is too few to test")
        print(f"  {config.GROUP_LABELS[grp]}: mean high-low = {g['diff'].mean():+.4f} "
              f"(3 pairs, descriptive only)")

    # how many pairs run in the expected direction
    n_pos = int((inj["diff"] > 0).sum())
    p_sign = float(stats.binomtest(n_pos, len(inj), 0.5).pvalue)
    sb.record("3", "cfos_injured_pairs", "n_pairs_high_gt_low", n_pos, n=len(inj),
              test="exact binomial sign test", p_value=p_sign,
              notes="direction consistency across injured pairs")
    print(f"  direction: {n_pos}/{len(inj)} injured pairs have high_risk > low_risk "
          f"(sign test p = {p_sign:.4g})")

    print("  NOTE: no pooled continuous risk-score regression across the 18 pools is run "
          "(risk score is not comparable across dose groups).")
    sb.record("3", "design", "pooled_risk_regression_run", "no",
              notes="deliberately omitted: risk score not comparable between dose groups")

    plotting.fig_cfos_paired(pairs)
    return out
