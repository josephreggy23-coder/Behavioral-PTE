"""STEP 3 -- orthogonal molecular validation using the paired c-fos pools.

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
POOL_MEMBERSHIP_COLUMNS = ["pool_id", "risk_pool", "group", "clutch", "fish_id"]
CLUTCH_SENSITIVITY_COLUMNS = (
    ["clutch", "n_group_pairs"]
    + [f"{group}_raw_high_minus_low" for group in config.GROUPS]
    + ["mean_raw_high_minus_low"]
    + [f"{group}_log2_high_minus_low" for group in config.GROUPS]
    + ["mean_log2_high_minus_low"]
)


def build_pool_membership(pools: pd.DataFrame | None = None) -> pd.DataFrame:
    """Normalize recorded pool membership; this does not reconstruct risk scoring."""
    if pools is None:
        pools = io_data.cfos_pools()

    required = ["pool_id", "risk_pool", "group", "clutch", "pooled_fish_ids"]
    missing = [col for col in required if col not in pools.columns]
    if missing:
        raise ValueError(f"c-fos pools are missing membership columns: {missing}")

    rows: list[dict[str, object]] = []
    fish_to_pool: dict[str, object] = {}
    for row in pools[required].itertuples(index=False):
        raw_ids = "" if pd.isna(row.pooled_fish_ids) else str(row.pooled_fish_ids)
        fish_ids = [fish_id.strip() for fish_id in raw_ids.split(";") if fish_id.strip()]
        if len(fish_ids) != 4 or len(set(fish_ids)) != 4:
            raise ValueError(
                f"Pool {row.pool_id!r} must contain exactly 4 unique fish IDs; "
                f"found {fish_ids!r}"
            )

        reused = sorted(fish_id for fish_id in fish_ids if fish_id in fish_to_pool)
        if reused:
            previous = {fish_id: fish_to_pool[fish_id] for fish_id in reused}
            raise ValueError(
                f"Fish IDs are reused in pool {row.pool_id!r}; previous pools: {previous}"
            )

        for fish_id in fish_ids:
            fish_to_pool[fish_id] = row.pool_id
            rows.append(
                {
                    "pool_id": row.pool_id,
                    "risk_pool": row.risk_pool,
                    "group": row.group,
                    "clutch": row.clutch,
                    "fish_id": fish_id,
                }
            )

    return pd.DataFrame(rows, columns=POOL_MEMBERSHIP_COLUMNS)


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


def build_clutch_sensitivity_table(pairs: pd.DataFrame) -> pd.DataFrame:
    """Average the three group-level paired differences within each clutch."""
    required = ["group", "clutch", "diff", "log2_diff"]
    missing = [col for col in required if col not in pairs.columns]
    if missing:
        raise ValueError(f"c-fos pair table is missing clutch-sensitivity columns: {missing}")

    expected_groups = set(config.GROUPS)
    rows: list[dict[str, object]] = []
    for clutch, block in pairs.groupby("clutch", observed=True, sort=False):
        observed_groups = set(block["group"].dropna().tolist())
        if len(block) != len(config.GROUPS) or observed_groups != expected_groups:
            raise ValueError(
                f"Clutch {clutch!r} must have one paired difference for every group; "
                f"found groups {sorted(observed_groups)!r}"
            )

        by_group = block.set_index("group")
        row: dict[str, object] = {"clutch": clutch, "n_group_pairs": len(block)}
        for group in config.GROUPS:
            row[f"{group}_raw_high_minus_low"] = float(by_group.loc[group, "diff"])
            row[f"{group}_log2_high_minus_low"] = float(by_group.loc[group, "log2_diff"])
        row["mean_raw_high_minus_low"] = float(block["diff"].mean())
        row["mean_log2_high_minus_low"] = float(block["log2_diff"].mean())
        rows.append(row)

    table = pd.DataFrame(rows, columns=CLUTCH_SENSITIVITY_COLUMNS)
    expected_clutches = set(config.CLUTCHES)
    observed_clutches = set(table["clutch"].tolist())
    if len(table) != len(config.CLUTCHES) or observed_clutches != expected_clutches:
        raise ValueError(
            "Clutch sensitivity requires exactly the configured clutches; "
            f"found {sorted(observed_clutches)!r}"
        )
    clutch_order = {clutch: i for i, clutch in enumerate(config.CLUTCHES)}
    return table.sort_values("clutch", key=lambda s: s.map(clutch_order)).reset_index(drop=True)


def _paired_test(
    d: pd.Series,
    label: str,
    scope: str,
    value_name: str,
    *,
    n_quantity: str = "n_pairs",
    unit_label: str = "pairs",
) -> dict:
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
    test_name = (
        "paired t-test (one-sample on within-pair differences)"
        if unit_label == "pairs"
        else "one-sample t-test on clutch-mean differences"
    )
    console_test_name = "paired t" if unit_label == "pairs" else "one-sample t"

    sb.record("3", scope, n_quantity, n, n=n, notes=label)
    sb.record("3", scope, f"mean_{value_name}_difference", mean, ci_low=float(ci[0]),
              ci_high=float(ci[1]), n=n, test=test_name,
              statistic=float(t), df=n - 1, p_value=float(p),
              effect_size_name="Cohen's dz", effect_size=float(dz), notes=label)
    sb.record("3", scope, f"wilcoxon_{value_name}", mean, n=n,
              test="Wilcoxon signed-rank (robustness check)", statistic=float(w),
              p_value=float(pw), notes=label)
    sb.record("3", scope, f"shapiro_{value_name}_differences", float(sh_w), n=n,
              test="Shapiro-Wilk on within-pair differences", statistic=float(sh_w),
              p_value=float(sh_p), notes="normality assumption for the paired t-test")

    print(f"  {label} (n = {n} {unit_label})")
    print(f"    mean difference = {mean:+.4f}  95% CI [{ci[0]:+.4f}, {ci[1]:+.4f}]")
    print(f"    {console_test_name}({n-1}) = {t:+.3f}, p = {p:.4g}, Cohen's dz = {dz:+.3f}")
    print(f"    Wilcoxon W = {w}, p = {pw:.4g}")
    sb.check(f"{label}: differences normal (Shapiro p > 0.05)",
             bool(np.isfinite(sh_p) and sh_p > 0.05), f"W = {sh_w:.3f}, p = {sh_p:.3g}")
    return {"n": n, "mean": mean, "ci": (float(ci[0]), float(ci[1])), "t": float(t),
            "p": float(p), "dz": float(dz), "w": float(w), "p_wilcoxon": float(pw),
            "shapiro_p": float(sh_p)}


def run(pools: pd.DataFrame | None = None) -> dict:
    sb.banner("STEP 3 -- orthogonal molecular validation: paired c-fos pools")
    source_pools = io_data.cfos_pools() if pools is None else pools.copy()
    membership = build_pool_membership(source_pools)
    membership.to_csv(config.TABLES / "step3_cfos_pool_membership.csv", index=False)

    pairs = build_pairs(source_pools)
    pairs.to_csv(config.TABLES / "step3_cfos_pairs.csv", index=False)

    n_cells = pairs.groupby(["group", "clutch"], observed=True).size()
    sb.check("exactly one pair per group x clutch cell",
             bool((n_cells == 1).all()) and len(pairs) == 9, f"{len(pairs)} pairs")
    sb.record("3", "design", "n_pools", 9 * 2, notes="analysed as 9 pairs, NOT 18 independent units")
    sb.record("3", "design", "n_pool_membership_rows", len(membership), n=len(source_pools),
              notes="4 unique recorded fish IDs per pool; no fish ID reused across pools")
    print(f"  pools: {len(pairs)*2} -> {len(pairs)} pairs (group x clutch)")
    print(f"  larvae per pool: {source_pools['n_larvae_in_pool'].unique().tolist()}")

    clutch_sensitivity = build_clutch_sensitivity_table(pairs)
    clutch_sensitivity.to_csv(
        config.TABLES / "step3_cfos_clutch_aggregated_sensitivity.csv", index=False
    )
    sb.record("3", "design", "n_clutch_sensitivity_units", len(clutch_sensitivity),
              n=len(clutch_sensitivity),
              notes="one mean of the 3 group-level paired differences per observed clutch")

    out = {
        "pairs": pairs,
        "pool_membership": membership,
        "clutch_sensitivity_table": clutch_sensitivity,
    }
    out["all"] = _paired_test(pairs["diff"], "All 9 pairs", "cfos_all_pairs", "fold_change")
    out["all_log2"] = _paired_test(pairs["log2_diff"], "All 9 pairs, log2 scale",
                                   "cfos_all_pairs_log2", "log2_fold_change")
    out["clutch_sensitivity_raw"] = _paired_test(
        clutch_sensitivity["mean_raw_high_minus_low"],
        "Clutch-aggregated sensitivity, raw fold-change scale",
        "cfos_clutch_aggregated",
        "fold_change",
        n_quantity="n_clutches",
        unit_label="clutches",
    )
    out["clutch_sensitivity_log2"] = _paired_test(
        clutch_sensitivity["mean_log2_high_minus_low"],
        "Clutch-aggregated sensitivity, log2 scale",
        "cfos_clutch_aggregated_log2",
        "log2_fold_change",
        n_quantity="n_clutches",
        unit_label="clutches",
    )

    inj = pairs[pairs["group"].isin(config.INJURED)]
    out["injured"] = _paired_test(inj["diff"], "Injured pairs only (low + high impact)",
                                  "cfos_injured_pairs", "fold_change")

    sham = pairs[pairs["group"] == "sham"]
    out["sham"] = _paired_test(sham["diff"], "Sham pairs only (descriptive comparison)",
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
