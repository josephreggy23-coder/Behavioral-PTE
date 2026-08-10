"""STEP 4 -- descriptive analyses, operational metrics, and the PTZ probe."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.proportion import proportion_confint

from . import config, io_data, plotting, statsbook as sb


# --------------------------------------------------------------------------
# Baseline equivalence and tau dynamics
# --------------------------------------------------------------------------
def baseline_equivalence(fits: pd.DataFrame) -> None:
    sb.banner("STEP 4a -- do the groups start equal at the pre-injury baseline?")
    base = fits[fits["timepoint_h"] == config.BASELINE_TP]
    samples = [base.loc[base["group"] == g, "decay_constant_fit"].dropna() for g in config.GROUPS]

    lev_s, lev_p = stats.levene(*samples, center="median")
    sb.record("4", "baseline_equivalence", "levene_tau", float(lev_s), n=len(base),
              test="Levene (median-centred) on baseline tau", statistic=float(lev_s),
              p_value=float(lev_p))
    sb.check("equal variances at baseline (Levene p > 0.05)", lev_p > 0.05,
             f"W = {lev_s:.3f}, p = {lev_p:.3f}")

    for g, s in zip(config.GROUPS, samples):
        w, p = stats.shapiro(s)
        sb.record("4", "baseline_equivalence", f"shapiro_tau_{g}", float(w), n=len(s),
                  test="Shapiro-Wilk", statistic=float(w), p_value=float(p))
        sb.check(f"baseline tau approx normal in {g}", p > 0.05, f"W = {w:.3f}, p = {p:.3f}")

    f, p = stats.f_oneway(*samples)
    k, n = len(samples), len(base)
    eta2 = (f * (k - 1)) / (f * (k - 1) + (n - k))
    sb.record("4", "baseline_equivalence", "anova_tau", float(f), n=n,
              test="one-way ANOVA on baseline tau", statistic=float(f),
              df=f"{k-1},{n-k}", p_value=float(p),
              effect_size_name="eta squared", effect_size=float(eta2),
              notes="a NON-significant result is the desired outcome here")
    h, ph = stats.kruskal(*samples)
    sb.record("4", "baseline_equivalence", "kruskal_tau", float(h), n=n,
              test="Kruskal-Wallis (robustness)", statistic=float(h), df=k - 1, p_value=float(ph))
    print(f"  baseline tau by group: " +
          ", ".join(f"{g} {s.mean():.2f}+/-{s.std(ddof=1):.2f}" for g, s in zip(config.GROUPS, samples)))
    print(f"  ANOVA F({k-1},{n-k}) = {f:.3f}, p = {p:.4f}, eta^2 = {eta2:.4f}; "
          f"Kruskal H = {h:.3f}, p = {ph:.4f}")
    sb.check("groups are equivalent at baseline (ANOVA p > 0.05)", p > 0.05, f"p = {p:.4f}")

    # baseline locomotion too
    loc = [base.loc[base["group"] == g, "baseline_locomotion"].dropna() for g in config.GROUPS]
    fl, pl = stats.f_oneway(*loc)
    sb.record("4", "baseline_equivalence", "anova_baseline_locomotion", float(fl), n=len(base),
              test="one-way ANOVA on pre-injury baseline_locomotion", statistic=float(fl),
              df=f"{k-1},{len(base)-k}", p_value=float(pl))
    sb.check("baseline locomotion equivalent across groups", pl > 0.05, f"p = {pl:.4f}")


def tau_dynamics(fits: pd.DataFrame) -> pd.DataFrame:
    sb.banner("STEP 4b -- tau by timepoint and group (the opposite-direction signature)")
    summ = (
        fits.groupby(["group", "timepoint_h"], observed=True)["decay_constant_fit"]
        .agg(n="size", mean="mean", sd=lambda s: s.std(ddof=1), sem="sem")
        .reset_index()
    )
    summ.to_csv(config.TABLES / "step4_tau_by_group_timepoint.csv", index=False)
    for r in summ.itertuples():
        sb.record("4", "tau_dynamics", f"{r.group}_t{r.timepoint_h:g}_mean", r.mean,
                  ci_low=r.mean - 1.96 * r.sem, ci_high=r.mean + 1.96 * r.sem, n=r.n,
                  notes=f"SD = {r.sd:.3f}")

    print(f"  {'group':<13}" + "".join(f"{f't={t:g}h':>12}" for t in config.TIMEPOINTS))
    for g in config.GROUPS:
        row = summ[summ["group"] == g].set_index("timepoint_h")["mean"].reindex(config.TIMEPOINTS)
        print(f"  {g:<13}" + "".join(f"{v:>12.2f}" for v in row.values))

    # within-fish change from baseline, per group, at each post-injury timepoint
    wide = fits.pivot_table(index=["fish_id", "group"], columns="timepoint_h",
                            values="decay_constant_fit", observed=True)
    for g in config.GROUPS:
        w = wide.xs(g, level="group")
        for tp in [t for t in config.TIMEPOINTS if t != config.BASELINE_TP]:
            d = (w[tp] - w[config.BASELINE_TP]).dropna()
            t, p = stats.ttest_1samp(d, 0.0)
            ci = stats.t.interval(0.95, len(d) - 1, loc=d.mean(),
                                  scale=d.std(ddof=1) / np.sqrt(len(d)))
            dz = d.mean() / d.std(ddof=1)
            sb.record("4", "tau_change_from_baseline", f"{g}_delta_t{tp:g}", float(d.mean()),
                      ci_low=float(ci[0]), ci_high=float(ci[1]), n=len(d),
                      test="paired t-test vs own pre-injury baseline", statistic=float(t),
                      df=len(d) - 1, p_value=float(p),
                      effect_size_name="Cohen's dz", effect_size=float(dz))
    d05 = {g: (wide.xs(g, level="group")[0.5] - wide.xs(g, level="group")[config.BASELINE_TP]).dropna()
           for g in config.INJURED}
    t_between, p_between = stats.ttest_ind(d05["low_impact"], d05["high_impact"], equal_var=False)
    pooled_sd = np.sqrt((d05["low_impact"].var(ddof=1) + d05["high_impact"].var(ddof=1)) / 2)
    dcohen = (d05["low_impact"].mean() - d05["high_impact"].mean()) / pooled_sd
    sb.record("4", "tau_change_from_baseline", "low_vs_high_delta_t0.5",
              float(d05["low_impact"].mean() - d05["high_impact"].mean()),
              n=len(d05["low_impact"]) + len(d05["high_impact"]),
              test="Welch t-test, low vs high impact delta-tau at 0.5 h",
              statistic=float(t_between), p_value=float(p_between),
              effect_size_name="Cohen's d", effect_size=float(dcohen),
              notes="the two doses move tau in opposite directions; this is why dose must be a predictor")
    print(f"  delta-tau at 0.5 h: low {d05['low_impact'].mean():+.2f} vs "
          f"high {d05['high_impact'].mean():+.2f} trials, Welch t = {t_between:.2f}, "
          f"p = {p_between:.3g}, d = {dcohen:.2f}")

    # what a dose-blind model would see: pooled |correlation| of delta-tau with outcome
    return summ


def converter_contrast(fits: pd.DataFrame) -> None:
    sb.banner("STEP 4c -- converters vs non-converters within injured fish")
    oc = io_data.outcomes()[["fish_id", "converted"]]
    d = fits.merge(oc, on="fish_id").query("group in @config.INJURED")
    wide = d.pivot_table(index=["fish_id", "group", "converted"], columns="timepoint_h",
                         values="decay_constant_fit", observed=True).reset_index()
    for g in config.INJURED:
        sub = wide[wide["group"] == g]
        for tp in [config.BASELINE_TP, 0.5, 24.0]:
            a = sub.loc[sub["converted"] == 1, tp].dropna()
            b = sub.loc[sub["converted"] == 0, tp].dropna()
            t, p = stats.ttest_ind(a, b, equal_var=False)
            psd = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
            dd = (a.mean() - b.mean()) / psd
            sb.record("4", "converter_contrast", f"{g}_tau_t{tp:g}", float(a.mean() - b.mean()),
                      n=len(a) + len(b), test="Welch t-test (converter - non-converter)",
                      statistic=float(t), p_value=float(p),
                      effect_size_name="Cohen's d", effect_size=float(dd),
                      notes=f"converters n = {len(a)}, non-converters n = {len(b)}")
        # delta tau at 0.5 h
        dl = (sub[0.5] - sub[config.BASELINE_TP])
        a, b = dl[sub["converted"] == 1].dropna(), dl[sub["converted"] == 0].dropna()
        t, p = stats.ttest_ind(a, b, equal_var=False)
        psd = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
        sb.record("4", "converter_contrast", f"{g}_dtau_0.5", float(a.mean() - b.mean()),
                  n=len(a) + len(b), test="Welch t-test on delta-tau at 0.5 h",
                  statistic=float(t), p_value=float(p),
                  effect_size_name="Cohen's d", effect_size=float((a.mean() - b.mean()) / psd))
        print(f"  {g}: delta-tau(0.5h) converters {a.mean():+.2f} vs non {b.mean():+.2f}, "
              f"t = {t:.2f}, p = {p:.3g}")
    plotting.fig_converter_trajectories(fits, io_data.outcomes())


# --------------------------------------------------------------------------
# Operational metrics
# --------------------------------------------------------------------------
def operations() -> pd.DataFrame:
    sb.banner("STEP 4d -- operational metrics from session_log")
    sl = io_data.session_log().copy()
    sl["total_min"] = sl["setup_min"] + sl["acquisition_min"] + sl["analysis_min"]
    sl["fish_per_hour"] = sl["n_fish_recorded"] / (sl["total_min"] / 60.0)
    sl["operator_min_per_fish"] = sl["operator_hands_on_min"] / sl["n_fish_recorded"]
    sl["cost_per_fish_usd"] = sl["consumables_cost_usd"] / sl["n_fish_recorded"]
    sl.to_csv(config.TABLES / "step4_session_metrics.csv", index=False)

    hab = sl[sl["session_kind"] == "habituation"]
    total_lost = int(sl["fish_lost_this_session"].sum())
    base = hab[hab["timepoint_h"] == config.BASELINE_TP]
    n_start = int(base["n_fish_recorded"].sum()) if len(base) else int(sl["n_fish_recorded"].max())
    attrition = total_lost / n_start if n_start else np.nan

    metrics = {
        "n_sessions": len(sl),
        "n_habituation_sessions": len(hab),
        "n_outcome_sessions": int((sl["session_kind"] == "outcome").sum()),
        "total_fish_sessions_recorded": int(sl["n_fish_recorded"].sum()),
        "fish_per_hour_mean": float(sl["fish_per_hour"].mean()),
        "fish_per_hour_sd": float(sl["fish_per_hour"].std(ddof=1)),
        "fish_per_hour_habituation": float(hab["fish_per_hour"].mean()),
        "operator_min_per_fish_mean": float(sl["operator_min_per_fish"].mean()),
        "operator_min_per_fish_sd": float(sl["operator_min_per_fish"].std(ddof=1)),
        "cost_per_fish_usd_mean": float(sl["cost_per_fish_usd"].mean()),
        "cost_per_fish_usd_sd": float(sl["cost_per_fish_usd"].std(ddof=1)),
        "total_consumables_usd": float(sl["consumables_cost_usd"].sum()),
        "total_operator_hours": float(sl["operator_hands_on_min"].sum() / 60.0),
        "fish_lost_total": total_lost,
        "fish_at_baseline": n_start,
        "attrition_rate": float(attrition),
    }
    for k, v in metrics.items():
        sb.record("4", "operations", k, v, n=len(sl),
                  notes="consumables_cost_usd treated as a per-session total (it scales with "
                        "n_fish_recorded)" if "cost" in k else "")
        print(f"  {k:<32} {v:.3f}" if isinstance(v, float) else f"  {k:<32} {v}")

    # the outcome recording is a longer, more expensive session than habituation
    for kind in ("habituation", "outcome"):
        sub = sl[sl["session_kind"] == kind]
        if len(sub):
            sb.record("4", "operations", f"{kind}_operator_min_per_fish",
                      float(sub["operator_min_per_fish"].mean()), n=len(sub),
                      notes=f"{kind} sessions only")
            print(f"  {kind:<12} operator min/fish = {sub['operator_min_per_fish'].mean():.2f}, "
                  f"acquisition {sub['acquisition_min'].mean():.0f} min")

    ci = proportion_confint(total_lost, n_start, method="wilson")
    sb.record("4", "operations", "attrition_rate_ci", float(attrition), ci_low=float(ci[0]),
              ci_high=float(ci[1]), n=n_start, test="Wilson 95% CI for a proportion")
    print(f"  attrition {total_lost}/{n_start} = {100*attrition:.1f}% "
          f"[{100*ci[0]:.1f}%, {100*ci[1]:.1f}%]")

    plotting.fig_operations(sl)
    return sl


# --------------------------------------------------------------------------
# PTZ probe -- secondary, underpowered
# --------------------------------------------------------------------------
def ptz_confound_check() -> None:
    """PTZ is now given to a subset of the SAME larvae that supply the outcome.

    That raises an obvious question: does receiving a proconvulsant change the
    conversion outcome? If it did, the PTZ subset would convert at a different
    rate and the primary analysis would be contaminated. This checks it.
    """
    sb.banner("STEP 4e -- is the PTZ challenge a confound for the conversion outcome?")
    oc = io_data.outcomes()
    if "ptz_challenged" not in oc.columns:
        print("  no ptz_challenged flag in outcomes; skipping")
        return
    print(oc.groupby(["group", "ptz_challenged"], observed=True)["converted"]
          .agg(["size", "sum"]).to_string())
    for g in config.GROUPS:
        d = oc[oc["group"] == g]
        a = d[d["ptz_challenged"] == 1]["converted"]
        b = d[d["ptz_challenged"] == 0]["converted"]
        table = np.array([[int(a.sum()), len(a) - int(a.sum())],
                          [int(b.sum()), len(b) - int(b.sum())]])
        _, p = stats.fisher_exact(table)
        sb.record("4", "ptz_confound", f"conversion_challenged_vs_not_{g}",
                  float(a.mean() - b.mean()), n=len(d),
                  test="Fisher exact, converted by ptz_challenged", p_value=float(p),
                  notes=f"challenged {int(a.sum())}/{len(a)}, not challenged "
                        f"{int(b.sum())}/{len(b)}")
        print(f"  {g:<12} challenged {a.mean():.3f} vs not {b.mean():.3f}, Fisher p = {p:.4f}")
        sb.check(f"PTZ challenge does not shift conversion in {g}", bool(p > 0.05),
                 f"p = {p:.4f}")


def ptz_vs_conversion() -> None:
    """PTZ and conversion are measured on the same animals, so they can be crossed."""
    oc = io_data.outcomes().set_index("fish_id")
    pz = io_data.ptz_challenge()
    d = pz.join(oc["converted"], on="fish_id").dropna(subset=["converted"])
    if d.empty:
        return
    table = pd.crosstab(d["converted"].astype(int), d["seized"])
    print("  seized x converted (PTZ subset):")
    print(table.to_string())
    if table.shape == (2, 2):
        odds, p = stats.fisher_exact(table.to_numpy())
        sb.record("4", "ptz_vs_conversion", "odds_ratio_seized_given_converted", float(odds),
                  n=len(d), test="Fisher exact, seized x converted within the PTZ subset",
                  p_value=float(p), effect_size_name="odds ratio", effect_size=float(odds),
                  notes="both measured on the same larvae; exploratory, not pre-specified")
        print(f"  Fisher exact OR = {odds:.2f}, p = {p:.4f}")
    a = d[d["converted"] == 1]["latency_s"]
    b = d[d["converted"] == 0]["latency_s"]
    if len(a) > 1 and len(b) > 1:
        u, pu = stats.mannwhitneyu(a, b, alternative="two-sided")
        sb.record("4", "ptz_vs_conversion", "latency_converter_minus_nonconverter",
                  float(a.median() - b.median()), n=len(d),
                  test="Mann-Whitney U on PTZ latency by conversion status",
                  statistic=float(u), p_value=float(pu),
                  notes=f"median converter {a.median():.0f} s, non-converter {b.median():.0f} s")
        print(f"  PTZ latency: converters median {a.median():.0f} s vs non-converters "
              f"{b.median():.0f} s, Mann-Whitney U = {u:.0f}, p = {pu:.4f}")


def ptz() -> pd.DataFrame:
    sb.banner("STEP 4f -- PTZ threshold probe (SECONDARY)")
    pz = io_data.ptz_challenge()
    rows = []
    for g in config.GROUPS:
        d = pz[pz["group"] == g]
        n, s = len(d), int(d["seized"].sum())
        lo, hi = proportion_confint(s, n, method="wilson")
        rows.append({"group": g, "n": n, "n_seized": s, "prop_seized": s / n,
                     "ci_low": float(lo), "ci_high": float(hi),
                     "median_latency_s": float(d["latency_s"].median()),
                     "ptz_mM": float(d["ptz_mM"].iloc[0])})
        sb.record("4", "ptz", f"prop_seized_{g}", s / n, ci_low=float(lo), ci_high=float(hi),
                  n=n, test="Wilson 95% CI", notes=f"{s}/{n} seized")
    summary = pd.DataFrame(rows)
    summary.to_csv(config.TABLES / "step4_ptz_summary.csv", index=False)
    print(summary.to_string(index=False))

    table = np.array([[r["n_seized"], r["n"] - r["n_seized"]] for r in rows])
    chi2, p, dof, expected = stats.chi2_contingency(table, correction=False)
    n_total = int(table.sum())
    cramers_v = float(np.sqrt(chi2 / (n_total * (min(table.shape) - 1))))
    sb.record("4", "ptz", "chi_square", float(chi2), n=n_total,
              test="Pearson chi-square, 3 groups x seized/not", statistic=float(chi2),
              df=int(dof), p_value=float(p),
              effect_size_name="Cramer's V", effect_size=cramers_v,
              notes="SECONDARY and UNDERPOWERED -- see power estimate below")
    print(f"  chi2({dof}) = {chi2:.3f}, p = {p:.4f}, Cramer's V = {cramers_v:.3f}")

    min_exp = float(expected.min())
    sb.record("4", "ptz", "min_expected_cell_count", min_exp, n=n_total,
              notes="chi-square is unreliable below 5; Fisher reported alongside")
    sb.check("all expected cell counts >= 5 (chi-square valid)", min_exp >= 5,
             f"minimum expected = {min_exp:.2f}")
    try:
        _, p_fisher = stats.fisher_exact(table[[1, 2], :]) if table.shape[0] > 2 else (np.nan, np.nan)
        sb.record("4", "ptz", "fisher_exact_low_vs_high", np.nan, n=int(table[[1, 2], :].sum()),
                  test="Fisher exact, low vs high impact", p_value=float(p_fisher))
        print(f"  Fisher exact (low vs high impact): p = {p_fisher:.4f}")
    except Exception:
        pass

    # crude post-hoc power for the observed sham-vs-injured contrast
    p_sham = summary.loc[summary["group"] == "sham", "prop_seized"].iloc[0]
    p_inj = float(pz[pz["group"].isin(config.INJURED)]["seized"].mean())
    n_sham = int(summary.loc[summary["group"] == "sham", "n"].iloc[0])
    n_inj = int(pz["group"].isin(config.INJURED).sum())
    h = 2 * np.arcsin(np.sqrt(p_inj)) - 2 * np.arcsin(np.sqrt(p_sham))
    n_eff = 2 / (1 / n_sham + 1 / n_inj)
    z = abs(h) * np.sqrt(n_eff / 2) - stats.norm.ppf(0.975)
    power = float(stats.norm.cdf(z))
    sb.record("4", "ptz", "post_hoc_power_sham_vs_injured", power, n=n_sham + n_inj,
              effect_size_name="Cohen's h", effect_size=float(h),
              notes="two-proportion z, alpha = 0.05, two-sided; reported to make the "
                    "underpowering explicit, not to reinterpret the p-value")
    print(f"  post-hoc power for sham ({p_sham:.2f}) vs injured ({p_inj:.2f}) at this n: "
          f"{power:.2f} (Cohen's h = {h:.2f})")

    # The power statement must follow the data, not a prior expectation.
    if power >= 0.80:
        print(f"  STATEMENT: at n = {n_sham + n_inj} the sham-versus-injured contrast is "
              f"adequately powered ({100*power:.0f}%). It remains a SECONDARY outcome: the "
              "primary claim is about individual prediction, which this group-level test does "
              "not address.")
        sb.record("4", "ptz", "power_verdict", "adequately powered for sham vs injured",
                  n=n_sham + n_inj, notes="secondary outcome regardless of power")
    else:
        print("  STATEMENT: the PTZ probe is underpowered at this n and is reported as a "
              "directional check only; no conclusion rests on it.")
        sb.record("4", "ptz", "power_verdict", "underpowered", n=n_sham + n_inj)

    # the low-vs-high dose contrast is a separate question and is not powered
    n_low = int(summary.loc[summary["group"] == "low_impact", "n"].iloc[0])
    n_high = int(summary.loc[summary["group"] == "high_impact", "n"].iloc[0])
    print(f"  NOTE: the low-versus-high dose contrast (n = {n_low} vs {n_high}) remains "
          "underpowered; the two injured groups are not distinguishable on this probe.")
    sb.record("4", "ptz", "low_vs_high_powered", "no", n=n_low + n_high,
              notes="dose-discrimination within injured groups is not powered at this n")

    lat = {g: pz.loc[pz['group'] == g, 'latency_s'].to_numpy() for g in config.GROUPS}
    plotting.fig_ptz(summary, lat)
    return summary


def run(fits: pd.DataFrame) -> dict:
    baseline_equivalence(fits)
    summ = tau_dynamics(fits)
    converter_contrast(fits)
    plotting.fig_habituation_curves(io_data.habituation_trials())
    plotting.fig_tau_trajectories(fits)
    sl = operations()
    ptz_confound_check()
    pz = ptz()
    ptz_vs_conversion()
    return {"tau_summary": summ, "sessions": sl, "ptz": pz}
