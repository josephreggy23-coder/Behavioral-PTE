"""Loading and light validation of the workbook.

Workbook: ``tbidataset.xlsx``. Six sheets, one cohort of ``zf_*`` larvae that
carry every measurement: habituation sessions at t = -1, 0.5 and 24 h, a
conversion outcome at 6.8 dpf, membership in a qPCR pool, and for a subset a PTZ
challenge.

Column aliases applied on load so downstream code has one name to use:
    outcomes.burst_events_per_hour   <- also accepted as burst_events_per_hour_6dpf
    cfos_pools.delta_ddct            <- delta_delta_ct
    session_log.session              <- parsed into session_kind / timepoint_h
"""
from __future__ import annotations

import re
from functools import lru_cache

import numpy as np
import pandas as pd

from . import config

SHEETS = [
    "habituation_trials",
    "fish_features",
    "outcomes",
    "cfos_pools",
    "ptz_challenge",
    "session_log",
]


def _parse_session(label: str) -> tuple[str, float]:
    """'habituation_t-1.0' -> ('habituation', -1.0); 'outcome_6.8dpf' -> ('outcome', nan)."""
    m = re.match(r"habituation_t(-?\d+(?:\.\d+)?)", str(label))
    if m:
        return "habituation", float(m.group(1))
    return ("outcome" if "outcome" in str(label) else "other"), float("nan")


@lru_cache(maxsize=1)
def _book() -> dict[str, pd.DataFrame]:
    if not config.DATA_XLSX.exists():
        raise FileNotFoundError(
            f"Dataset not found at {config.DATA_XLSX}. "
            "Place the zebrafish TBI workbook there (see README)."
        )
    xl = pd.ExcelFile(config.DATA_XLSX)
    missing = [s for s in SHEETS if s not in xl.sheet_names]
    if missing:
        raise ValueError(f"Workbook is missing required sheets: {missing}")
    book = {s: xl.parse(s) for s in SHEETS}

    oc = book["outcomes"]
    if "burst_events_per_hour" not in oc.columns:
        for alt in ("burst_events_per_hour_6dpf", "burst_events_per_h"):
            if alt in oc.columns:
                oc["burst_events_per_hour"] = oc[alt]
                break
    book["outcomes"] = oc

    pools = book["cfos_pools"]
    if "delta_ddct" not in pools.columns and "delta_delta_ct" in pools.columns:
        pools["delta_ddct"] = pools["delta_delta_ct"]
    if "pool_type" not in pools.columns and "risk_pool" in pools.columns:
        pools["pool_type"] = pools["risk_pool"]
    book["cfos_pools"] = pools

    sl = book["session_log"]
    if "session" in sl.columns:
        parsed = sl["session"].map(_parse_session)
        sl["session_kind"] = [p[0] for p in parsed]
        sl["timepoint_h"] = [p[1] for p in parsed]
    elif "timepoint_h" in sl.columns:
        sl["session_kind"] = "habituation"
    book["session_log"] = sl

    for df in book.values():
        if "group" in df.columns:
            df["group"] = pd.Categorical(df["group"], categories=config.GROUPS, ordered=True)
    return book


def sheet(name: str) -> pd.DataFrame:
    return _book()[name].copy()


def habituation_trials() -> pd.DataFrame:
    return sheet("habituation_trials")


def fish_features() -> pd.DataFrame:
    return sheet("fish_features")


def outcomes() -> pd.DataFrame:
    return sheet("outcomes")


def cfos_pools() -> pd.DataFrame:
    return sheet("cfos_pools")


def ptz_challenge() -> pd.DataFrame:
    return sheet("ptz_challenge")


def session_log() -> pd.DataFrame:
    return sheet("session_log")


def describe_design() -> pd.DataFrame:
    """Per group x clutch counts of larvae and conversion events."""
    oc = outcomes()
    tab = (
        oc.groupby(["group", "clutch"], observed=True)
        .agg(n_fish=("fish_id", "size"), n_converted=("converted", "sum"))
        .reset_index()
    )
    tab["conversion_rate"] = tab["n_converted"] / tab["n_fish"]
    return tab


def study_metadata() -> dict:
    """Ages and intervals, which are now recorded per fish in the outcomes sheet."""
    oc = outcomes()
    out = {}
    for col in ("age_at_injury_dpf", "age_at_outcome_dpf", "days_post_injury"):
        if col in oc.columns:
            vals = pd.unique(oc[col].dropna())
            out[col] = float(vals[0]) if len(vals) == 1 else sorted(float(v) for v in vals)
    out["n_ptz_challenged"] = (
        int(oc["ptz_challenged"].sum()) if "ptz_challenged" in oc.columns else np.nan
    )
    return out
