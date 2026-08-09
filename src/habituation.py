"""STEP 1 -- rebuild the outcome variable from trial-level data.

Model, per fish per session, over trials k = 1..30:

    distance_mm(k) = A * exp(-(k - 1) / tau) + C

fitted by NONLINEAR least squares (scipy.optimize.curve_fit).

Why not log-linearise
---------------------
The usual shortcut is to subtract an estimated offset, take logs and regress
log(y - C) on k.  That fails here for a mechanical reason: once responses reach
the habituated floor, individual trials scatter *below* the offset, so y - C is
negative for a sizeable fraction of late trials.  Those points are either
dropped (biasing the tail upward, inflating tau) or produce NaNs; if the offset
is under-estimated to avoid that, the residual floor tilts the regression slope
positive and the recovered tau comes back NEGATIVE.  ``loglinear_diagnostic``
below quantifies exactly how often that happens in this dataset, and it is
printed as an assumption check rather than asserted.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from scipy.optimize import OptimizeWarning, curve_fit
from scipy import stats

from . import config, io_data, statsbook as sb


# --------------------------------------------------------------------------
# Model
# --------------------------------------------------------------------------
def habituation_model(k, amplitude, tau, offset):
    """A * exp(-(k-1)/tau) + C, with k the 1-based trial index."""
    return amplitude * np.exp(-(k - 1.0) / tau) + offset


def _initial_guesses(k: np.ndarray, y: np.ndarray) -> list[tuple[float, float, float]]:
    """Several starting points; curve_fit on 30 noisy points can find local minima."""
    early = float(np.mean(y[: max(3, len(y) // 6)]))
    late = float(np.mean(y[-max(3, len(y) // 6) :]))
    amp0 = max(early - late, 1e-3)
    off0 = max(late, 1e-6)
    return [
        (amp0, 5.0, off0),
        (amp0, 2.0, off0),
        (amp0, 12.0, off0),
        (max(float(y.max() - y.min()), 1e-3), 8.0, max(float(y.min()), 1e-6)),
    ]


def fit_one_session(trials: np.ndarray, distance: np.ndarray) -> dict:
    """Fit one fish-session. Returns parameters plus fit diagnostics."""
    k = np.asarray(trials, dtype=float)
    y = np.asarray(distance, dtype=float)
    ok = np.isfinite(k) & np.isfinite(y)
    k, y = k[ok], y[ok]

    out = {
        "amplitude_fit": np.nan,
        "decay_constant_fit": np.nan,
        "offset_fit": np.nan,
        "tau_se": np.nan,
        "r2": np.nan,
        "rmse": np.nan,
        "n_points": int(len(y)),
        "converged": False,
        "tau_at_bound": False,
    }
    if len(y) < 5:
        return out

    lower = [config.AMP_BOUNDS[0], config.TAU_BOUNDS[0], config.OFF_BOUNDS[0]]
    upper = [config.AMP_BOUNDS[1], config.TAU_BOUNDS[1], config.OFF_BOUNDS[1]]

    best = None
    for p0 in _initial_guesses(k, y):
        p0 = [float(np.clip(v, lo, hi)) for v, lo, hi in zip(p0, lower, upper)]
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", OptimizeWarning)
                popt, pcov = curve_fit(
                    habituation_model,
                    k,
                    y,
                    p0=p0,
                    bounds=(lower, upper),
                    maxfev=config.CURVEFIT_MAXFEV,
                )
        except (RuntimeError, ValueError):
            continue
        resid = y - habituation_model(k, *popt)
        sse = float(np.sum(resid**2))
        if best is None or sse < best[0]:
            best = (sse, popt, pcov)

    if best is None:
        return out

    sse, popt, pcov = best
    sst = float(np.sum((y - y.mean()) ** 2))
    amp, tau, off = (float(v) for v in popt)
    tau_se = float(np.sqrt(pcov[1, 1])) if np.all(np.isfinite(pcov)) else np.nan

    out.update(
        amplitude_fit=amp,
        decay_constant_fit=tau,
        offset_fit=off,
        tau_se=tau_se,
        r2=1.0 - sse / sst if sst > 0 else np.nan,
        rmse=float(np.sqrt(sse / len(y))),
        converged=True,
        tau_at_bound=bool(
            np.isclose(tau, config.TAU_BOUNDS[0], rtol=1e-3)
            or np.isclose(tau, config.TAU_BOUNDS[1], rtol=1e-3)
        ),
    )
    return out


def loglinear_tau(trials: np.ndarray, distance: np.ndarray) -> dict:
    """The log-linearised estimator, kept ONLY to demonstrate that it breaks.

    Offset is estimated as the mean of the last block, then log(y - C) is
    regressed on (k-1) and tau taken as -1/slope.  Points with y <= C cannot be
    logged and are dropped -- we count them.
    """
    k = np.asarray(trials, dtype=float)
    y = np.asarray(distance, dtype=float)
    c_hat = float(np.mean(y[-5:]))
    resid = y - c_hat
    usable = resid > 0
    n_dropped = int((~usable).sum())
    if usable.sum() < 3:
        return {"tau_loglin": np.nan, "n_dropped": n_dropped, "frac_dropped": n_dropped / len(y)}
    slope, *_ = stats.linregress((k[usable] - 1.0), np.log(resid[usable]))
    tau = -1.0 / slope if slope != 0 else np.nan
    return {"tau_loglin": float(tau), "n_dropped": n_dropped, "frac_dropped": n_dropped / len(y)}


# --------------------------------------------------------------------------
# Batch refit
# --------------------------------------------------------------------------
def refit_all_sessions(trials: pd.DataFrame | None = None) -> pd.DataFrame:
    """Refit every fish-session in habituation_trials. One row per session."""
    if trials is None:
        trials = io_data.habituation_trials()

    rows = []
    keys = ["fish_id", "group", "clutch", "timepoint_h"]
    for key, g in trials.groupby(keys, observed=True, sort=False):
        g = g.sort_values("trial")
        rec = dict(zip(keys, key))
        rec.update(fit_one_session(g["trial"].to_numpy(), g["distance_mm"].to_numpy()))
        rec.update(loglinear_tau(g["trial"].to_numpy(), g["distance_mm"].to_numpy()))
        rec["baseline_locomotion"] = float(g["baseline_locomotion"].iloc[0])
        rec["mean_distance_obs"] = float(g["distance_mm"].mean())
        rec["resp_prob_block1_obs"] = float(g.loc[g["block"] == 1, "responded"].mean())
        rec["resp_prob_block6_obs"] = float(
            g.loc[g["block"] == config.N_BLOCKS, "responded"].mean()
        )
        rows.append(rec)

    fits = pd.DataFrame(rows)
    fits["group"] = pd.Categorical(fits["group"], categories=config.GROUPS, ordered=True)
    return fits.sort_values(["fish_id", "timepoint_h"]).reset_index(drop=True)


# --------------------------------------------------------------------------
# Agreement with the supplied feature table
# --------------------------------------------------------------------------
def compare_with_supplied(fits: pd.DataFrame) -> pd.DataFrame:
    """Merge refit tau against fish_features.decay_constant and quantify agreement."""
    supplied = io_data.fish_features()[
        ["fish_id", "timepoint_h", "amplitude", "decay_constant", "offset"]
    ].rename(
        columns={
            "amplitude": "amplitude_supplied",
            "decay_constant": "decay_constant_supplied",
            "offset": "offset_supplied",
        }
    )
    merged = fits.merge(supplied, on=["fish_id", "timepoint_h"], how="left", validate="1:1")
    merged["tau_diff"] = merged["decay_constant_fit"] - merged["decay_constant_supplied"]
    return merged


def report_step1(merged: pd.DataFrame) -> pd.DataFrame:
    """Print assumption checks and record every Step 1 statistic."""
    sb.banner("STEP 1 -- nonlinear refit of the habituation decay constant")

    n = len(merged)
    n_conv = int(merged["converged"].sum())
    sb.record("1", "curve_fit", "n_sessions", n, n=n, notes="fish_* cohort, all timepoints")
    sb.record("1", "curve_fit", "n_converged", n_conv, n=n)
    sb.record(
        "1", "curve_fit", "convergence_rate", n_conv / n, n=n,
        notes="scipy.optimize.curve_fit, bounded, 4 starting points per session",
    )
    print(f"  sessions fitted: {n}   converged: {n_conv} ({100*n_conv/n:.1f}%)")

    r2 = merged["r2"].dropna()
    for q, v in [("r2_median", r2.median()), ("r2_q05", r2.quantile(0.05)), ("r2_min", r2.min())]:
        sb.record("1", "curve_fit", q, float(v), n=len(r2))
    print(f"  R^2 median {r2.median():.3f}  (5th pct {r2.quantile(0.05):.3f}, min {r2.min():.3f})")

    n_bound = int(merged["tau_at_bound"].sum())
    sb.record("1", "curve_fit", "n_tau_at_bound", n_bound, n=n)

    sb.check("all sessions converged", n_conv == n, f"{n - n_conv} failures")
    sb.check("tau not pinned at a parameter bound", n_bound == 0,
             f"{n_bound}/{n} sessions hit a bound")
    sb.check("fit quality R^2 > 0.5 for >=95% of sessions", float((r2 > 0.5).mean()) >= 0.95,
             f"{100*float((r2 > 0.5).mean()):.1f}% above 0.5")

    # ---- agreement with supplied features -------------------------------
    ok = merged[["decay_constant_fit", "decay_constant_supplied"]].dropna()
    r_p, p_p = stats.pearsonr(ok["decay_constant_fit"], ok["decay_constant_supplied"])
    r_s, p_s = stats.spearmanr(ok["decay_constant_fit"], ok["decay_constant_supplied"])
    diff = merged["tau_diff"].dropna()
    bias, sd = float(diff.mean()), float(diff.std(ddof=1))
    loa = (bias - 1.96 * sd, bias + 1.96 * sd)
    ci = stats.t.interval(0.95, len(diff) - 1, loc=bias, scale=sd / np.sqrt(len(diff)))

    sb.record("1", "tau_agreement", "pearson_r", float(r_p), n=len(ok),
              test="Pearson correlation", statistic=float(r_p), p_value=float(p_p),
              notes="refit tau vs supplied fish_features.decay_constant")
    sb.record("1", "tau_agreement", "spearman_rho", float(r_s), n=len(ok),
              test="Spearman correlation", statistic=float(r_s), p_value=float(p_s))
    sb.record("1", "tau_agreement", "mean_bias_trials", bias, ci_low=float(ci[0]),
              ci_high=float(ci[1]), n=len(diff), test="one-sample t (paired difference)",
              notes="refit minus supplied")
    sb.record("1", "tau_agreement", "limits_of_agreement_low", loa[0], n=len(diff),
              notes="Bland-Altman, bias -/+ 1.96 SD")
    sb.record("1", "tau_agreement", "limits_of_agreement_high", loa[1], n=len(diff))
    sb.record("1", "tau_agreement", "sd_of_difference", sd, n=len(diff))
    med_abs_pct = float((diff.abs() / merged["decay_constant_supplied"]).median() * 100)
    sb.record("1", "tau_agreement", "median_abs_pct_error", med_abs_pct, n=len(diff))

    print(f"  agreement with supplied tau: r = {r_p:.4f} (p = {p_p:.3g}), rho = {r_s:.4f}")
    print(f"  bias (refit - supplied) = {bias:+.3f} trials, 95% CI [{ci[0]:+.3f}, {ci[1]:+.3f}]")
    print(f"  limits of agreement [{loa[0]:+.3f}, {loa[1]:+.3f}]; median |error| {med_abs_pct:.2f}%")
    sb.check("refit reproduces supplied tau (r > 0.95)", r_p > 0.95, f"r = {r_p:.4f}")
    sb.check("no systematic bias (|bias| < 0.25 trials)", abs(bias) < 0.25, f"bias = {bias:+.3f}")

    # ---- the log-linearisation failure ----------------------------------
    ll = merged["tau_loglin"]
    n_neg = int((ll < 0).sum())
    n_nan = int(ll.isna().sum())
    frac_dropped = float(merged["frac_dropped"].mean())
    valid = merged.dropna(subset=["tau_loglin"])
    r_ll, p_ll = stats.pearsonr(valid["tau_loglin"], valid["decay_constant_supplied"])

    sb.record("1", "loglinear_failure", "pct_trials_below_offset", 100 * frac_dropped, n=n,
              notes="points unusable after offset subtraction (y <= C)")
    sb.record("1", "loglinear_failure", "n_sessions_negative_tau", n_neg, n=n,
              notes="sign-inverted tau from log-linearisation")
    sb.record("1", "loglinear_failure", "n_sessions_undefined_tau", n_nan, n=n)
    sb.record("1", "loglinear_failure", "pearson_r_vs_supplied", float(r_ll), n=len(valid),
              test="Pearson correlation", statistic=float(r_ll), p_value=float(p_ll),
              notes="log-linearised tau vs supplied tau; compare with the NLS r")
    print(f"  log-linearisation: {100*frac_dropped:.1f}% of trials fall below the offset; "
          f"{n_neg} sessions return a NEGATIVE tau, {n_nan} undefined; r vs supplied = {r_ll:.3f}")
    sb.check("log-linearisation is unusable (as predicted)", (n_neg + n_nan) > 0 or r_ll < r_p,
             f"NLS r = {r_p:.3f} vs log-linear r = {r_ll:.3f}")

    merged.to_csv(config.TABLES / "step1_session_fits.csv", index=False)
    return merged
