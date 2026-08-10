"""Global configuration: paths, seed, constants, plot style.

Every stochastic component in this pipeline draws from SEED (or a documented
offset of it) so that a rerun reproduces every number in RESULTS.md exactly.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------
# Reproducibility
# --------------------------------------------------------------------------
SEED = 20260809

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[1]
DATA_XLSX = ROOT / "tbidataset.xlsx"

RESULTS = ROOT / "results"
FIGURES = RESULTS / "figures"
TABLES = RESULTS / "tables"

for _d in (RESULTS, FIGURES, TABLES):
    _d.mkdir(parents=True, exist_ok=True)

STATS_CSV = RESULTS / "all_statistics.csv"
RESULTS_MD = ROOT / "RESULTS.md"

# --------------------------------------------------------------------------
# Experimental design constants
# --------------------------------------------------------------------------
GROUPS = ["sham", "low_impact", "high_impact"]
INJURED = ["low_impact", "high_impact"]
CLUTCHES = ["clutch_A", "clutch_B", "clutch_C"]
TIMEPOINTS = [-1.0, 0.5, 24.0]
BASELINE_TP = -1.0  # pre-injury session: every fish is its own control
N_TRIALS = 30
N_BLOCKS = 6

# Timepoints entering the prediction model
DELTA_TPS = [0.5, 24.0]

# --------------------------------------------------------------------------
# Curve-fit settings (STEP 1)
# --------------------------------------------------------------------------
# distance_mm(k) = A * exp(-(k - 1) / tau) + C   for trial k = 1..30
TAU_BOUNDS = (0.1, 100.0)
AMP_BOUNDS = (0.0, 200.0)
OFF_BOUNDS = (0.0, 200.0)
CURVEFIT_MAXFEV = 20000

# --------------------------------------------------------------------------
# Model settings (STEP 2)
# --------------------------------------------------------------------------
C_GRID = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 100.0]
N_PERMUTATIONS = 1000
N_BOOTSTRAP = 2000
N_OUTER_FOLDS = 3  # == number of clutches -> leave-one-clutch-out
N_INNER_FOLDS = 2  # 2 clutches remain inside each outer training set
N_RANDOM_SPLIT_FOLDS = 5  # for the leakage-quantifying naive comparison

# --------------------------------------------------------------------------
# Plot style
# --------------------------------------------------------------------------
DPI = 300

GROUP_COLORS = {
    "sham": "#6e6e6e",
    "low_impact": "#1f77b4",
    "high_impact": "#d62728",
}
GROUP_LABELS = {
    "sham": "Sham",
    "low_impact": "Low impact",
    "high_impact": "High impact",
}
OUTCOME_COLORS = {0: "#4c9f70", 1: "#b5179e"}

# c-fos pools are now defined by realised outcome, not by predicted risk
POOL_TYPES = ["non_converter", "converter"]
POOL_LABELS = {"non_converter": "Non-converter pool", "converter": "Converter pool"}
POOL_COLORS = {"non_converter": "#4c9f70", "converter": "#b5179e"}


def apply_style() -> None:
    """Consistent, publication-ish matplotlib defaults."""
    plt.rcParams.update(
        {
            "figure.dpi": 110,
            "savefig.dpi": DPI,
            "savefig.bbox": "tight",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "grid.linewidth": 0.6,
            "legend.frameon": False,
            "lines.linewidth": 1.6,
            "figure.autolayout": False,
        }
    )
