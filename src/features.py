"""Build the one-row-per-fish design matrix for STEP 2.

Predictors (exactly four, by protocol):
    dose         high_impact = 1, low_impact = 0
    pre_tau      decay_constant at t = -1 (the pre-injury baseline session)
    dtau_0.5     tau@0.5 - tau@-1   -> z-scored WITHIN dose group
    dtau_24      tau@24  - tau@-1   -> z-scored WITHIN dose group

The within-dose z-scoring is deliberately NOT applied here.  It is a
data-dependent transform, so it lives inside the modelling pipeline
(``modeling.WithinDoseZScorer``) where it is fitted on training folds only.
Applying it to the whole dataset up front would leak test-fold information into
the scaling constants.  A globally z-scored copy is emitted alongside so the
leakage-free and leaky variants can be compared as a sensitivity check.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import config, io_data, statsbook as sb

PREDICTORS = ["dose", "pre_tau", "dtau_0.5", "dtau_24"]

MODEL_SPECS: dict[str, dict] = {
    "a_locomotion_only": {
        "label": "(a) baseline_locomotion only",
        "features": ["baseline_locomotion_pre"],
        "question": "Is it just sickness?",
    },
    "b_dose_only": {
        "label": "(b) dose only",
        "features": ["dose"],
        "question": "Is it just injury severity?",
    },
    "c_dose_pretau": {
        "label": "(c) dose + pre_tau",
        "features": ["dose", "pre_tau"],
        "question": "Is it a pre-existing trait?",
    },
    "d_dose_dtau": {
        "label": "(d) dose + z_dtau",
        "features": ["dose", "dtau_0.5", "dtau_24"],
        "question": "Is it the injury response?",
    },
    "e_full": {
        "label": "(e) full 4-predictor model",
        "features": ["dose", "pre_tau", "dtau_0.5", "dtau_24"],
        "question": "Full model",
    },
}

# columns that must be z-scored within dose group wherever they appear
WITHIN_DOSE_Z_COLS = ["dtau_0.5", "dtau_24"]


def build_fish_table(fits: pd.DataFrame, tau_col: str = "decay_constant_fit") -> pd.DataFrame:
    """Wide, one row per followed fish, all groups (sham dropped later)."""
    wide = fits.pivot_table(index="fish_id", columns="timepoint_h", values=tau_col)
    wide.columns = [f"tau_{c:g}" for c in wide.columns]

    loco = fits.pivot_table(index="fish_id", columns="timepoint_h", values="baseline_locomotion")
    loco.columns = [f"loco_{c:g}" for c in loco.columns]

    meta = (
        fits.groupby("fish_id", observed=True)[["group", "clutch"]].first()
    )
    tbl = meta.join(wide).join(loco)

    tbl["pre_tau"] = tbl[f"tau_{config.BASELINE_TP:g}"]
    tbl["baseline_locomotion_pre"] = tbl[f"loco_{config.BASELINE_TP:g}"]
    for tp in config.DELTA_TPS:
        tbl[f"dtau_{tp:g}"] = tbl[f"tau_{tp:g}"] - tbl["pre_tau"]

    keep = ["converted", "burst_events_per_hour"]
    if "ptz_challenged" in io_data.outcomes().columns:
        keep.append("ptz_challenged")
    oc = io_data.outcomes().set_index("fish_id")[keep]
    tbl = tbl.join(oc, how="left")
    tbl["dose"] = (tbl["group"] == "high_impact").astype(int)
    return tbl.reset_index()


def injured_modeling_set(tbl: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Injured fish only, complete on the four predictors and the outcome.

    Returns (analysis_set, dropped_rows).
    """
    inj = tbl[tbl["group"].isin(config.INJURED)].copy()
    needed = PREDICTORS + ["baseline_locomotion_pre", "converted"]
    complete = inj.dropna(subset=needed)
    dropped = inj[~inj["fish_id"].isin(complete["fish_id"])]
    return complete.reset_index(drop=True), dropped.reset_index(drop=True)


def add_global_within_dose_z(df: pd.DataFrame) -> pd.DataFrame:
    """Whole-dataset within-dose z-scores (leaky; sensitivity analysis only)."""
    out = df.copy()
    for col in WITHIN_DOSE_Z_COLS:
        out[f"z_{col}"] = out.groupby("dose")[col].transform(
            lambda s: (s - s.mean()) / s.std(ddof=0)
        )
    return out


def report_features(tbl: pd.DataFrame, model_df: pd.DataFrame, dropped: pd.DataFrame) -> None:
    sb.banner("STEP 2a -- design matrix")

    n, ev = len(model_df), int(model_df["converted"].sum())
    epv = ev / len(PREDICTORS)
    print(f"  injured fish with complete data: n = {n}, events = {ev} "
          f"({100*ev/n:.1f}%), non-events = {n-ev}")
    print(f"  events per variable = {ev}/{len(PREDICTORS)} = {epv:.1f}")
    if len(dropped):
        miss = ", ".join(f"{r.fish_id} ({r.group})" for r in dropped.itertuples())
        print(f"  excluded for missing sessions ({len(dropped)}): {miss}")

    sb.record("2", "design_matrix", "n_fish", n, n=n, notes="injured only, complete cases")
    sb.record("2", "design_matrix", "n_events", ev, n=n)
    sb.record("2", "design_matrix", "event_rate", ev / n, n=n)
    sb.record("2", "design_matrix", "events_per_variable", epv, n=n,
              notes="4 predictors; >=10 is the conventional target, >=9 the accepted minimum")
    sb.record("2", "design_matrix", "n_excluded_incomplete", len(dropped),
              notes="missing a 0.5 h or 24 h session (attrition)")

    sb.check("events per variable >= 10 (conventional target)", epv >= 10, f"EPV = {epv:.1f}")
    sb.check("events per variable >= 9 (accepted minimum)", epv >= 9, f"EPV = {epv:.1f}")
    if epv < 9:
        print(f"  NOTE: EPV = {epv:.1f} is below the conventional floor of 9-10 events per "
              f"variable. The predictor set is held at {len(PREDICTORS)} by design rather than "
              "trimmed, because each term answers a distinct pre-registered question; the "
              "consequence is wider coefficient intervals, which are reported as such.")
    sb.check("outcome not degenerate (10-90% events)", 0.10 < ev / n < 0.90,
             f"event rate {ev/n:.2f}")

    by = model_df.groupby(["group", "clutch"], observed=True)["converted"].agg(["size", "sum"])
    print("  events by group x clutch (outer-fold units):")
    for (g, c), row in by.iterrows():
        print(f"    {g:<12} {c}: n={int(row['size']):2d}  converted={int(row['sum']):2d}")
        sb.record("2", "design_matrix", f"n_{g}_{c}", int(row["size"]),
                  notes=f"converted = {int(row['sum'])}")
    per_clutch = model_df.groupby("clutch", observed=True)["converted"].agg(["size", "sum"])
    ok = bool((per_clutch["sum"] > 0).all() and (per_clutch["sum"] < per_clutch["size"]).all())
    sb.check("every clutch has both classes (GroupKFold folds are scorable)", ok)

    # collinearity among the four predictors (pre-scaling, raw deltas)
    corr = model_df[PREDICTORS].corr()
    corr.to_csv(config.TABLES / "step2_predictor_correlations.csv")
    worst = float(np.abs(corr.to_numpy()[np.triu_indices(len(PREDICTORS), 1)]).max())
    sb.record("2", "design_matrix", "max_abs_predictor_correlation", worst, n=n)
    sb.check("no severe collinearity (max |r| < 0.8)", worst < 0.8, f"max |r| = {worst:.3f}")

    model_df.to_csv(config.TABLES / "step2_model_matrix.csv", index=False)
    tbl.to_csv(config.TABLES / "fish_wide_table.csv", index=False)
