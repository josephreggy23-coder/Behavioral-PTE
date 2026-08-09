"""A single append-only ledger for every statistic the pipeline reports.

Any number that appears in RESULTS.md must first pass through ``record``.
The ledger is written once at the end of the run to results/all_statistics.csv.
"""
from __future__ import annotations

from typing import Any

import pandas as pd

from . import config

_ROWS: list[dict[str, Any]] = []

COLUMNS = [
    "step",
    "analysis",
    "quantity",
    "value",
    "ci_low",
    "ci_high",
    "n",
    "test",
    "statistic",
    "df",
    "p_value",
    "effect_size_name",
    "effect_size",
    "notes",
]


def record(
    step: str,
    analysis: str,
    quantity: str,
    value: Any = None,
    *,
    ci_low: Any = None,
    ci_high: Any = None,
    n: Any = None,
    test: str = "",
    statistic: Any = None,
    df: Any = None,
    p_value: Any = None,
    effect_size_name: str = "",
    effect_size: Any = None,
    notes: str = "",
) -> None:
    _ROWS.append(
        {
            "step": step,
            "analysis": analysis,
            "quantity": quantity,
            "value": value,
            "ci_low": ci_low,
            "ci_high": ci_high,
            "n": n,
            "test": test,
            "statistic": statistic,
            "df": df,
            "p_value": p_value,
            "effect_size_name": effect_size_name,
            "effect_size": effect_size,
            "notes": notes,
        }
    )


def as_frame() -> pd.DataFrame:
    return pd.DataFrame(_ROWS, columns=COLUMNS)


def get(analysis: str, quantity: str) -> Any:
    """Look a recorded value back up (used when writing RESULTS.md)."""
    for row in reversed(_ROWS):
        if row["analysis"] == analysis and row["quantity"] == quantity:
            return row
    raise KeyError(f"{analysis}/{quantity} not recorded")


def flush() -> pd.DataFrame:
    df = as_frame()
    df.to_csv(config.STATS_CSV, index=False)
    return df


def banner(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def check(label: str, passed: bool, detail: str = "") -> bool:
    """Print an assumption check in a uniform, greppable format."""
    tag = "PASS" if passed else "FLAG"
    print(f"  [{tag}] {label}" + (f" -- {detail}" if detail else ""))
    return passed
