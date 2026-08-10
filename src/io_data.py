"""Loading and light validation of the workbook.

The workbook shipped for this experiment is
``data/raw/behavioral_pte_source.xlsx``. Sheet names are as specified in the
protocol; one column is named ``delta_delta_ct`` rather than ``delta_ddct`` and
is aliased here so downstream code can use either.
"""
from __future__ import annotations

from functools import lru_cache

import pandas as pd

from . import config

SHEETS = [
    "habituation_trials",
    "fish_features",
    "outcomes_6dpf",
    "cfos_cohort_features",
    "cfos_pools",
    "ptz_challenge",
    "session_log",
]


@lru_cache(maxsize=1)
def _book() -> dict[str, pd.DataFrame]:
    if not config.DATA_XLSX.exists():
        raise FileNotFoundError(
            f"Dataset not found at {config.DATA_XLSX}. "
            "Restore the source workbook as data/raw/behavioral_pte_source.xlsx "
            "relative to the repository root."
        )
    xl = pd.ExcelFile(config.DATA_XLSX)
    missing = [s for s in SHEETS if s not in xl.sheet_names]
    if missing:
        raise ValueError(f"Workbook is missing required sheets: {missing}")
    book = {s: xl.parse(s) for s in SHEETS}

    pools = book["cfos_pools"]
    if "delta_ddct" not in pools.columns and "delta_delta_ct" in pools.columns:
        pools["delta_ddct"] = pools["delta_delta_ct"]
    book["cfos_pools"] = pools

    for name, df in book.items():
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
    return sheet("outcomes_6dpf")


def cfos_cohort_features() -> pd.DataFrame:
    return sheet("cfos_cohort_features")


def cfos_pools() -> pd.DataFrame:
    return sheet("cfos_pools")


def ptz_challenge() -> pd.DataFrame:
    return sheet("ptz_challenge")


def session_log() -> pd.DataFrame:
    return sheet("session_log")


def describe_design() -> pd.DataFrame:
    """Per group x clutch counts of followed fish and conversion events."""
    oc = outcomes()
    tab = (
        oc.groupby(["group", "clutch"], observed=True)
        .agg(n_fish=("fish_id", "size"), n_converted=("converted", "sum"))
        .reset_index()
    )
    tab["conversion_rate"] = tab["n_converted"] / tab["n_fish"]
    return tab
