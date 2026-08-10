"""All figures. Every figure is written as a 300 dpi PNG into results/figures."""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.calibration import calibration_curve

from . import config

config.apply_style()


def _save(fig, name: str) -> str:
    path = config.FIGURES / name
    fig.savefig(path, dpi=config.DPI)
    plt.close(fig)
    print(f"  figure -> {path.relative_to(config.ROOT)}")
    return str(path)


# ==========================================================================
# STEP 1
# ==========================================================================
def fig_curvefit_examples(trials: pd.DataFrame, fits: pd.DataFrame, model_fn) -> str:
    """Example fitted curves: one fish per group at baseline and at 24 h."""
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 5.6), sharex=True, sharey=True)
    for j, grp in enumerate(config.GROUPS):
        fid = fits.loc[fits["group"] == grp, "fish_id"].iloc[0]
        for i, tp in enumerate([config.BASELINE_TP, 24.0]):
            ax = axes[i, j]
            d = trials[(trials["fish_id"] == fid) & (trials["timepoint_h"] == tp)].sort_values("trial")
            row = fits[(fits["fish_id"] == fid) & (fits["timepoint_h"] == tp)]
            ax.plot(d["trial"], d["distance_mm"], "o", ms=3, alpha=0.6,
                    color=config.GROUP_COLORS[grp], label="trials")
            if len(row):
                r = row.iloc[0]
                k = np.linspace(1, config.N_TRIALS, 200)
                ax.plot(k, model_fn(k, r["amplitude_fit"], r["decay_constant_fit"], r["offset_fit"]),
                        "-", color="k", lw=1.6, label="NLS fit")
                ax.axhline(r["offset_fit"], ls=":", lw=1, color="0.5")
                ax.set_title(f"{config.GROUP_LABELS[grp]} | t = {tp:g} h\n"
                             f"$\\tau$ = {r['decay_constant_fit']:.2f}, $R^2$ = {r['r2']:.2f}",
                             fontsize=9)
            if i == 1:
                ax.set_xlabel("Trial")
            if j == 0:
                ax.set_ylabel("Distance (mm)")
    axes[0, 0].legend(fontsize=7, loc="upper right")
    fig.suptitle("Step 1: nonlinear fits of $A\\,e^{-(k-1)/\\tau}+C$ (dotted line = fitted offset $C$)",
                 fontsize=10)
    fig.tight_layout()
    return _save(fig, "fig01_curvefit_examples.png")


def fig_tau_agreement(merged: pd.DataFrame) -> str:
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))

    ax = axes[0]
    for grp in config.GROUPS:
        d = merged[merged["group"] == grp]
        ax.scatter(d["decay_constant_supplied"], d["decay_constant_fit"], s=8, alpha=0.55,
                   color=config.GROUP_COLORS[grp], label=config.GROUP_LABELS[grp])
    lim = [0, float(np.nanmax(merged[["decay_constant_fit", "decay_constant_supplied"]].to_numpy())) * 1.05]
    ax.plot(lim, lim, "k--", lw=1)
    ax.set(xlabel="Supplied decay_constant", ylabel="Refit $\\tau$ (NLS)", xlim=lim, ylim=lim,
           title="Refit vs supplied")
    ax.legend(fontsize=7)

    ax = axes[1]
    mean = (merged["decay_constant_fit"] + merged["decay_constant_supplied"]) / 2
    diff = merged["tau_diff"]
    bias, sd = diff.mean(), diff.std(ddof=1)
    ax.scatter(mean, diff, s=8, alpha=0.5, color="#333333")
    ax.axhline(bias, color="C3", lw=1.3, label=f"bias {bias:+.3f}")
    ax.axhline(bias + 1.96 * sd, color="C3", ls="--", lw=1)
    ax.axhline(bias - 1.96 * sd, color="C3", ls="--", lw=1, label="95% LoA")
    ax.set(xlabel="Mean of the two estimates", ylabel="Refit - supplied",
           title="Bland-Altman")
    ax.legend(fontsize=7)

    ax = axes[2]
    ll = merged["tau_loglin"].replace([np.inf, -np.inf], np.nan)
    ax.scatter(merged["decay_constant_supplied"], ll, s=8, alpha=0.5, color="#8c564b")
    ax.axhline(0, color="k", lw=1)
    n_neg = int((ll < 0).sum())
    ax.set(xlabel="Supplied decay_constant", ylabel="Log-linearised $\\tau$",
           title=f"Log-linearisation inverts sign\n({n_neg} sessions with $\\tau$ < 0)")
    fig.tight_layout()
    return _save(fig, "fig02_tau_agreement.png")


# ==========================================================================
# STEP 2
# ==========================================================================
def fig_roc_comparison(roc_data: list[dict]) -> str:
    """All five nested-CV models on one axes."""
    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    cmap = plt.get_cmap("viridis")
    for i, d in enumerate(roc_data):
        col = "k" if d["key"] == "e_full" else cmap(i / max(len(roc_data) - 1, 1) * 0.85)
        lw = 2.2 if d["key"] == "e_full" else 1.4
        ax.plot(d["fpr"], d["tpr"], color=col, lw=lw,
                label=f"{d['label']}  AUC = {d['auc']:.3f}")
    ax.plot([0, 1], [0, 1], "--", color="0.6", lw=1, label="Chance")
    ax.set(xlabel="False positive rate", ylabel="True positive rate",
           title="Nested comparison: out-of-fold ROC\n(leave-one-clutch-out, C tuned inside each fold)")
    ax.legend(fontsize=7.5, loc="lower right")
    fig.tight_layout()
    return _save(fig, "fig03_roc_nested_comparison.png")


def fig_permutation_null(null: np.ndarray, observed: float, p_value: float, n_perm: int) -> str:
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.hist(null, bins=40, color="#9ecae1", edgecolor="white", label=f"Null ({n_perm} permutations)")
    ax.axvline(observed, color="#d62728", lw=2.2,
               label=f"Observed AUC = {observed:.3f}\npermutation p = {p_value:.4f}")
    ax.axvline(float(np.quantile(null, 0.95)), color="k", ls="--", lw=1.2,
               label=f"Null 95th pct = {np.quantile(null, 0.95):.3f}")
    ax.set(xlabel="Pooled out-of-fold AUC", ylabel="Permutations",
           title="Permutation null: labels shuffled within clutch,\nentire nested CV rerun each time")
    ax.legend(fontsize=7.5)
    fig.tight_layout()
    return _save(fig, "fig04_permutation_null.png")


def fig_confusion_and_calibration(y: np.ndarray, oof: np.ndarray, cm: dict, brier: float) -> str:
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.0))

    ax = axes[0]
    m = cm["matrix"]
    ax.imshow(m, cmap="Blues")
    for (i, j), v in np.ndenumerate(m):
        ax.text(j, i, str(v), ha="center", va="center",
                color="white" if v > m.max() / 2 else "black", fontsize=13)
    ax.set(xticks=[0, 1], yticks=[0, 1],
           xticklabels=["pred non-conv.", "pred converter"],
           yticklabels=["non-converter", "converter"],
           title=(f"Confusion matrix (out-of-fold, thr = {cm['threshold']:.2f})\n"
                  f"sens {cm['sensitivity']:.2f} | spec {cm['specificity']:.2f} | "
                  f"bal.acc {cm['balanced_accuracy']:.2f}"))
    ax.grid(False)

    ax = axes[1]
    frac_pos, mean_pred = calibration_curve(y, oof, n_bins=5, strategy="quantile")
    ax.plot([0, 1], [0, 1], "--", color="0.6", lw=1, label="Perfect calibration")
    ax.plot(mean_pred, frac_pos, "o-", color="#d62728", label=f"Model (Brier = {brier:.3f})")
    ax.set(xlabel="Mean predicted probability", ylabel="Observed conversion fraction",
           title="Calibration (5 quantile bins, out-of-fold)", xlim=(0, 1), ylim=(0, 1))
    ax.legend(fontsize=7.5)

    ax = axes[2]
    for cls, lbl in [(0, "Non-converters"), (1, "Converters")]:
        ax.hist(oof[y == cls], bins=np.linspace(0, 1, 16), alpha=0.6,
                color=config.OUTCOME_COLORS[cls], label=lbl)
    ax.axvline(cm["threshold"], color="k", ls="--", lw=1)
    ax.set(xlabel="Out-of-fold predicted probability", ylabel="Fish",
           title="Score separation")
    ax.legend(fontsize=7.5)
    fig.tight_layout()
    return _save(fig, "fig05_confusion_calibration.png")


def fig_coefficients(coefs: pd.DataFrame) -> str:
    fig, ax = plt.subplots(figsize=(6.0, 3.4))
    y = np.arange(len(coefs))[::-1]
    ax.errorbar(coefs["coef_standardised"], y,
                xerr=[coefs["coef_standardised"] - coefs["ci_low"],
                      coefs["ci_high"] - coefs["coef_standardised"]],
                fmt="o", color="#333333", capsize=3, lw=1.3)
    ax.axvline(0, color="0.5", lw=1, ls="--")
    ax.set(yticks=y, yticklabels=coefs["predictor"],
           xlabel="Ridge logistic coefficient (per SD, standardised scale)",
           title=("Full model coefficients\n(within-clutch stratified bootstrap 95% CI; "
                  "conditional on observed clutches)"))
    fig.tight_layout()
    return _save(fig, "fig06_coefficients.png")


def fig_fold_auc(comparison: pd.DataFrame, folds_long: pd.DataFrame) -> str:
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    keys = comparison["key"].tolist()
    xpos = np.arange(len(keys))
    ax.bar(xpos, comparison["mean_fold_auc"], color="#c6dbef", edgecolor="#3182bd", zorder=1)
    for i, k in enumerate(keys):
        f = folds_long[folds_long["key"] == k]
        ax.scatter(np.full(len(f), i) + np.linspace(-0.12, 0.12, len(f)), f["auc"],
                   color="k", s=22, zorder=3)
    ax.axhline(0.5, color="0.5", ls="--", lw=1)
    ax.set(xticks=xpos, ylabel="AUC", ylim=(0.0, 1.0),
           title="Nested CV AUC by model (bars = mean of folds, dots = individual clutch folds)")
    ax.set_xticklabels(
        [label.replace(" ", "\n", 1) for label in comparison["label"]], fontsize=7.5
    )
    fig.tight_layout()
    return _save(fig, "fig07_fold_auc_by_model.png")


# ==========================================================================
# STEP 3
# ==========================================================================
def fig_cfos_paired(pairs: pd.DataFrame) -> str:
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.2))

    ax = axes[0]
    for i, r in enumerate(pairs.itertuples()):
        col = config.GROUP_COLORS[r.group]
        ax.plot([0, 1], [r.low_risk, r.high_risk], "-o", color=col, ms=5, alpha=0.85, lw=1.5)
        ax.annotate(f"{r.clutch[-1]}", (1.02, r.high_risk), fontsize=6.5, color=col,
                    va="center")
    ax.set(xticks=[0, 1], xticklabels=["low_risk pool", "high_risk pool"],
           ylabel="c-fos (fosab) fold change vs rpl13a",
           title=f"Paired pools: {len(pairs)} group x clutch cells")
    ax.axhline(1.0, color="0.6", ls=":", lw=1)
    ax.set_xlim(-0.25, 1.25)
    handles = [Line2D([], [], color=config.GROUP_COLORS[g], marker="o",
                      label=config.GROUP_LABELS[g]) for g in config.GROUPS]
    ax.legend(handles=handles, fontsize=7.5)

    ax = axes[1]
    order = [g for g in config.GROUPS]
    for i, g in enumerate(order):
        d = pairs[pairs["group"] == g]["diff"]
        ax.scatter(np.full(len(d), i) + np.linspace(-0.07, 0.07, len(d)), d,
                   color=config.GROUP_COLORS[g], s=45, zorder=3)
        ax.hlines(d.mean(), i - 0.2, i + 0.2, color="k", lw=2, zorder=4)
    ax.axhline(0, color="0.5", ls="--", lw=1)
    ax.set(xticks=range(len(order)), xticklabels=[config.GROUP_LABELS[g] for g in order],
           ylabel="high_risk - low_risk fold change",
           title="Within-pair difference by group\n(black bar = group mean, n = 3 pairs each)")
    fig.tight_layout()
    return _save(fig, "fig08_cfos_paired.png")


# ==========================================================================
# STEP 4
# ==========================================================================
def fig_habituation_curves(trials: pd.DataFrame) -> str:
    tps = config.TIMEPOINTS
    fig, axes = plt.subplots(1, len(tps), figsize=(3.0 * len(tps), 3.3), sharey=True)
    for ax, tp in zip(axes, tps):
        d = trials[trials["timepoint_h"] == tp]
        for grp in config.GROUPS:
            g = d[d["group"] == grp].groupby("trial")["distance_mm"]
            m, se = g.mean(), g.sem()
            ax.plot(m.index, m.values, color=config.GROUP_COLORS[grp],
                    label=config.GROUP_LABELS[grp])
            ax.fill_between(m.index, m - se, m + se, color=config.GROUP_COLORS[grp], alpha=0.22,
                            lw=0)
        title = "t = -1 h (PRE-INJURY)" if tp == config.BASELINE_TP else f"t = {tp:g} h"
        ax.set(xlabel="Trial", title=title)
    axes[0].set_ylabel("Distance (mm)")
    axes[0].legend(fontsize=7.5)
    fig.suptitle("Habituation curves by group and timepoint -- groups are superimposed at baseline",
                 fontsize=10)
    fig.tight_layout()
    return _save(fig, "fig09_habituation_curves.png")


def fig_tau_trajectories(fits: pd.DataFrame) -> str:
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    x = np.arange(len(config.TIMEPOINTS))
    for k, grp in enumerate(config.GROUPS):
        d = fits[fits["group"] == grp]
        g = d.groupby("timepoint_h")["decay_constant_fit"]
        m = g.mean().reindex(config.TIMEPOINTS)
        se = g.sem().reindex(config.TIMEPOINTS)
        ax.errorbar(x + (k - 1) * 0.05, m.values, yerr=se.values, fmt="-o", capsize=3,
                    color=config.GROUP_COLORS[grp], label=config.GROUP_LABELS[grp], ms=5)
    ax.axvline(0.5, color="0.7", ls=":", lw=1.2)
    ax.text(0.52, ax.get_ylim()[1], " blast", fontsize=7, va="top", color="0.4")
    ax.set(xticks=x, xticklabels=[f"{t:g}" for t in config.TIMEPOINTS],
           xlabel="Timepoint (h; -1 = pre-injury)", ylabel="Decay constant $\\tau$ (trials)",
           title="$\\tau$ by timepoint and group (mean $\\pm$ SEM)\n"
                 "Low and high dose move in OPPOSITE directions")
    ax.legend(fontsize=8)
    fig.tight_layout()
    return _save(fig, "fig10_tau_by_timepoint.png")


def fig_converter_trajectories(fits: pd.DataFrame, outcomes: pd.DataFrame) -> str:
    d = fits.merge(outcomes[["fish_id", "converted"]], on="fish_id", how="inner")
    d = d[d["group"].isin(config.INJURED)]
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.2), sharey=True)
    x = np.arange(len(config.TIMEPOINTS))
    for ax, grp in zip(axes, config.INJURED):
        dd = d[d["group"] == grp]
        for cls, lbl in [(0, "Non-converter"), (1, "Converter")]:
            s = dd[dd["converted"] == cls]
            for fid, f in s.groupby("fish_id"):
                f = f.set_index("timepoint_h").reindex(config.TIMEPOINTS)
                ax.plot(x, f["decay_constant_fit"].values, color=config.OUTCOME_COLORS[cls],
                        alpha=0.16, lw=0.8)
            g = s.groupby("timepoint_h")["decay_constant_fit"]
            m = g.mean().reindex(config.TIMEPOINTS)
            se = g.sem().reindex(config.TIMEPOINTS)
            ax.errorbar(x, m.values, yerr=se.values, fmt="-o", capsize=3, lw=2.4, ms=6,
                        color=config.OUTCOME_COLORS[cls],
                        label=f"{lbl} (n = {s['fish_id'].nunique()})")
        ax.set(xticks=x, xticklabels=[f"{t:g}" for t in config.TIMEPOINTS],
               xlabel="Timepoint (h)", title=config.GROUP_LABELS[grp])
        ax.legend(fontsize=7.5)
    axes[0].set_ylabel("Decay constant $\\tau$ (trials)")
    fig.suptitle("Converter vs non-converter $\\tau$ trajectories, within injured fish", fontsize=10)
    fig.tight_layout()
    return _save(fig, "fig11_converter_trajectories.png")


def fig_operations(sl: pd.DataFrame) -> str:
    fig, axes = plt.subplots(1, 4, figsize=(14.0, 3.4))
    x = np.arange(len(sl))
    labels = [f"{r.clutch[-1]}\n{r.timepoint_h:g}h" for r in sl.itertuples()]
    # neutral clutch palette -- the group colours are reserved for dose groups
    clutch_palette = dict(zip(config.CLUTCHES, ["#8da0cb", "#66c2a5", "#fc8d62"]))
    colors = [clutch_palette.get(c, "#999999") for c in sl["clutch"]]

    for ax, col, ylab, title in [
        (axes[0], "fish_per_hour", "Fish / hour", "Throughput"),
        (axes[1], "operator_min_per_fish", "Operator min / fish", "Hands-on time"),
        (axes[2], "cost_per_fish_usd", "USD / fish", "Consumables cost"),
        (axes[3], "fish_lost_this_session", "Fish lost", "Attrition"),
    ]:
        ax.bar(x, sl[col], color=colors)
        ax.axhline(sl[col].mean(), color="k", ls="--", lw=1,
                   label=f"mean {sl[col].mean():.2f}")
        ax.set(xticks=x, ylabel=ylab, title=title)
        ax.set_xticklabels(labels, fontsize=6)
        ax.legend(fontsize=7)
    fig.suptitle("Operational metrics per session (clutch / timepoint)", fontsize=10)
    fig.tight_layout()
    return _save(fig, "fig12_operations.png")


def fig_ptz(summary: pd.DataFrame) -> str:
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.8))
    ax = axes[0]
    x = np.arange(len(summary))
    ax.bar(x, summary["prop_seized"],
           color=[config.GROUP_COLORS[g] for g in summary["group"]])
    for i, r in enumerate(summary.itertuples()):
        ax.errorbar(i, r.prop_seized, yerr=[[r.prop_seized - r.ci_low], [r.ci_high - r.prop_seized]],
                    color="k", capsize=4, lw=1.2)
        ax.text(i, 0.03, f"{int(r.n_seized)}/{int(r.n)}", ha="center", fontsize=8, color="white")
    ax.set(xticks=x, xticklabels=[config.GROUP_LABELS[g] for g in summary["group"]],
           ylabel="Proportion seizing", ylim=(0, 1.05),
           title="PTZ challenge (SECONDARY, UNDERPOWERED)\nWilson 95% CI")
    ax = axes[1]
    ax.axis("off")
    ax.text(0.0, 0.95,
             "Underpowered at this sample size:\n"
             f"n = {int(summary['n'].sum())} across 3 groups\n"
             "Interpret as a directional check only;\n"
             "no confirmatory claim rests on it.",
            fontsize=9, va="top")
    fig.tight_layout()
    return _save(fig, "fig13_ptz.png")
