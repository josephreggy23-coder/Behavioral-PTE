"""STEP 2 -- L2-penalised logistic regression with nested, clutch-aware CV.

Design decisions, all deliberate:

* **Ridge logistic regression only.**  n = 81 with 36 events, so a fixed,
  regularised linear classifier limits model flexibility and yields inspectable
  coefficients.  Those coefficients are associations, not mechanistic effects.
* **GroupKFold on clutch for the outer split.**  A random split mixes
  clutch-associated observations across training and test partitions.  With 3
  clutches this is leave-one-clutch-out.
* **Nested CV for C.**  Selecting the penalty on the same folds we report
  optimistically biases AUC at this sample size.  C is chosen inside each outer
  training set only.
* **Within-dose z-scoring inside the pipeline.**  The z-scores are fitted on the
  training fold and applied to the held-out fold, so no test-fold information
  reaches the scaling constants.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    brier_score_loss,
    confusion_matrix,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import GroupKFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from . import config, features

DOSE_COL = "dose"


# --------------------------------------------------------------------------
# Within-dose z-scoring as a fold-safe transformer
# --------------------------------------------------------------------------
class WithinDoseZScorer(BaseEstimator, TransformerMixin):
    """Z-score selected columns within each level of the dose column.

    Fitted on the training fold only.  If the dose column is absent from the
    feature set, or no target columns are present, this is a no-op.
    """

    def __init__(self, cols: tuple[str, ...] = (), dose_col: str = DOSE_COL):
        self.cols = cols
        self.dose_col = dose_col

    def fit(self, X: pd.DataFrame, y=None):
        self.cols_ = [c for c in self.cols if c in X.columns]
        self.active_ = bool(self.cols_) and self.dose_col in X.columns
        self.stats_: dict = {}
        if self.active_:
            for level, block in X.groupby(self.dose_col):
                for c in self.cols_:
                    mu = float(block[c].mean())
                    sd = float(block[c].std(ddof=0))
                    self.stats_[(level, c)] = (mu, sd if sd > 1e-12 else 1.0)
            self.global_ = {
                c: (float(X[c].mean()), max(float(X[c].std(ddof=0)), 1e-12)) for c in self.cols_
            }
        return self

    def transform(self, X: pd.DataFrame) -> pd.DataFrame:
        if not self.active_:
            return X.copy()
        out = X.copy()
        for c in self.cols_:
            vals = np.empty(len(out), dtype=float)
            for i, (level, v) in enumerate(zip(out[self.dose_col].to_numpy(), out[c].to_numpy())):
                mu, sd = self.stats_.get((level, c), self.global_[c])
                vals[i] = (v - mu) / sd
            out[c] = vals
        return out

    def get_feature_names_out(self, input_features=None):
        return np.asarray(input_features)


def make_pipeline(feature_cols: list[str], C: float = 1.0, seed: int = config.SEED) -> Pipeline:
    z_cols = tuple(c for c in features.WITHIN_DOSE_Z_COLS if c in feature_cols)
    return Pipeline(
        [
            ("withindose_z", WithinDoseZScorer(cols=z_cols)),
            ("scaler", StandardScaler()),
            (
                # L2 (ridge) penalty. This is scikit-learn's default penalty and
                # is left implicit on purpose: the explicit `penalty="l2"`
                # argument is deprecated from scikit-learn 1.8. C is the inverse
                # penalty strength and is what the inner CV loop tunes.
                "logreg",
                LogisticRegression(
                    C=C,
                    solver="liblinear",
                    max_iter=5000,
                    random_state=seed,
                ),
            ),
        ]
    )


# --------------------------------------------------------------------------
# Nested cross-validation
# --------------------------------------------------------------------------
def _safe_auc(y_true, score) -> float:
    y_true = np.asarray(y_true)
    if len(np.unique(y_true)) < 2:
        return np.nan
    return float(roc_auc_score(y_true, score))


def _select_C(X, y, groups, feature_cols, inner_splits, C_grid, seed):
    """Inner loop: pick C by mean AUC over clutch-held-out inner folds."""
    inner = GroupKFold(n_splits=inner_splits)
    splits = list(inner.split(X, y, groups))
    scores = []
    for C in C_grid:
        fold_scores = []
        for tr, va in splits:
            pipe = make_pipeline(feature_cols, C=C, seed=seed)
            pipe.fit(X.iloc[tr], y[tr])
            p = pipe.predict_proba(X.iloc[va])[:, 1]
            fold_scores.append(_safe_auc(y[va], p))
        scores.append(np.nanmean(fold_scores) if np.any(np.isfinite(fold_scores)) else np.nan)
    scores = np.asarray(scores, dtype=float)
    if not np.any(np.isfinite(scores)):
        return float(np.median(C_grid))
    # ties -> strongest penalty (smallest C), the conservative choice
    best = int(np.nanargmax(scores))
    best_val = scores[best]
    for i, s in enumerate(scores):
        if np.isfinite(s) and np.isclose(s, best_val, atol=1e-12):
            best = i
            break
    return float(C_grid[best])


def nested_cv(
    X: pd.DataFrame,
    y: np.ndarray,
    groups: np.ndarray,
    feature_cols: list[str],
    *,
    outer_splits: int = config.N_OUTER_FOLDS,
    inner_splits: int = config.N_INNER_FOLDS,
    C_grid: list[float] | None = None,
    seed: int = config.SEED,
    random_outer: bool = False,
) -> dict:
    """Run nested CV and return out-of-fold predictions plus per-fold detail."""
    C_grid = C_grid or config.C_GRID
    Xf = X[feature_cols]

    if random_outer:
        outer = StratifiedKFold(n_splits=outer_splits, shuffle=True, random_state=seed)
        splits = list(outer.split(Xf, y))
    else:
        outer = GroupKFold(n_splits=outer_splits)
        splits = list(outer.split(Xf, y, groups))

    oof = np.full(len(y), np.nan)
    fold_rows = []
    for k, (tr, te) in enumerate(splits):
        inner_groups = groups[tr]
        if random_outer or len(np.unique(inner_groups)) < inner_splits:
            # random-split control: keep the inner loop random too
            inner = StratifiedKFold(n_splits=inner_splits, shuffle=True, random_state=seed + k)
            inner_splits_list = list(inner.split(Xf.iloc[tr], y[tr]))
            scores = []
            for C in C_grid:
                fs = []
                for itr, iva in inner_splits_list:
                    pipe = make_pipeline(feature_cols, C=C, seed=seed)
                    pipe.fit(Xf.iloc[tr].iloc[itr], y[tr][itr])
                    fs.append(_safe_auc(y[tr][iva], pipe.predict_proba(Xf.iloc[tr].iloc[iva])[:, 1]))
                scores.append(np.nanmean(fs))
            best_C = float(C_grid[int(np.nanargmax(scores))])
        else:
            best_C = _select_C(
                Xf.iloc[tr], y[tr], inner_groups, feature_cols, inner_splits, C_grid, seed
            )

        pipe = make_pipeline(feature_cols, C=best_C, seed=seed)
        pipe.fit(Xf.iloc[tr], y[tr])
        p = pipe.predict_proba(Xf.iloc[te])[:, 1]
        oof[te] = p
        fold_rows.append(
            {
                "fold": k,
                "held_out": "random" if random_outer else str(np.unique(groups[te])[0]),
                "n_test": int(len(te)),
                "n_events_test": int(y[te].sum()),
                "selected_C": best_C,
                "auc": _safe_auc(y[te], p),
            }
        )

    folds = pd.DataFrame(fold_rows)
    return {
        "oof": oof,
        "folds": folds,
        "mean_fold_auc": float(np.nanmean(folds["auc"])),
        "sd_fold_auc": float(np.nanstd(folds["auc"], ddof=1)) if len(folds) > 1 else np.nan,
        "min_fold_auc": float(np.nanmin(folds["auc"])),
        "max_fold_auc": float(np.nanmax(folds["auc"])),
        "pooled_auc": _safe_auc(y, oof),
        "brier": float(brier_score_loss(y, oof)),
    }


# --------------------------------------------------------------------------
# Permutation test
# --------------------------------------------------------------------------
def _shuffle_within_groups(y: np.ndarray, groups: np.ndarray, rng: np.random.Generator):
    out = y.copy()
    for g in np.unique(groups):
        idx = np.where(groups == g)[0]
        out[idx] = rng.permutation(y[idx])
    return out


def _one_permutation(X, y, groups, feature_cols, seed_i, kwargs):
    rng = np.random.default_rng(seed_i)
    y_perm = _shuffle_within_groups(y, groups, rng)
    res = nested_cv(X, y_perm, groups, feature_cols, **kwargs)
    return res["pooled_auc"], res["mean_fold_auc"]


def skipped_permutation_result() -> dict:
    """Return the report context for a permutation test that was not run."""
    return {
        "null_pooled": np.array([], dtype=float),
        "null_mean_fold": np.array([], dtype=float),
        "p_pooled": None,
        "p_mean_fold": None,
        "null_mean": None,
        "null_sd": None,
        "null_q95": None,
        "percentile_of_observed": None,
        "n_perm": 0,
        "performed": False,
    }


def permutation_test(
    X: pd.DataFrame,
    y: np.ndarray,
    groups: np.ndarray,
    feature_cols: list[str],
    observed_pooled: float,
    observed_mean_fold: float,
    *,
    n_perm: int = config.N_PERMUTATIONS,
    seed: int = config.SEED,
    n_jobs: int = -1,
    **kwargs,
) -> dict:
    """Shuffle labels WITHIN clutch and rerun the entire nested CV n_perm times.

    Within-clutch shuffling preserves each clutch's conversion rate, so the null
    is 'no fish-level information beyond clutch prevalence' rather than 'no
    information at all'.
    """
    if n_perm < 0:
        raise ValueError("n_perm must be non-negative")
    if n_perm == 0:
        return skipped_permutation_result()

    out = Parallel(n_jobs=n_jobs, verbose=0)(
        delayed(_one_permutation)(X, y, groups, feature_cols, seed + 1 + i, kwargs)
        for i in range(n_perm)
    )
    null_pooled = np.array([o[0] for o in out], dtype=float)
    null_mean = np.array([o[1] for o in out], dtype=float)
    p_pooled = (1 + int(np.sum(null_pooled >= observed_pooled))) / (n_perm + 1)
    p_mean = (1 + int(np.sum(null_mean >= observed_mean_fold))) / (n_perm + 1)
    return {
        "null_pooled": null_pooled,
        "null_mean_fold": null_mean,
        "p_pooled": p_pooled,
        "p_mean_fold": p_mean,
        "null_mean": float(np.nanmean(null_pooled)),
        "null_sd": float(np.nanstd(null_pooled, ddof=1)),
        "null_q95": float(np.nanquantile(null_pooled, 0.95)),
        "percentile_of_observed": float(100 * np.mean(null_pooled < observed_pooled)),
        "n_perm": n_perm,
        "performed": True,
    }


# --------------------------------------------------------------------------
# Bootstrap CI for a pooled out-of-fold AUC, stratified within clutch and
# conditional on the observed clutches
# --------------------------------------------------------------------------
def bootstrap_auc_ci(
    y: np.ndarray, scores: np.ndarray, groups: np.ndarray,
    n_boot: int = config.N_BOOTSTRAP, seed: int = config.SEED
) -> tuple[float, float]:
    rng = np.random.default_rng(seed + 991)
    idx_by_group = {g: np.where(groups == g) [0] for g in np.unique(groups)}
    vals = []
    for _ in range(n_boot):
        take = np.concatenate([rng.choice(ix, size=len(ix), replace=True)
                               for ix in idx_by_group.values()])
        a = _safe_auc(y[take], scores[take])
        if np.isfinite(a):
            vals.append(a)
    if not vals:
        return (np.nan, np.nan)
    return (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)))


# --------------------------------------------------------------------------
# Coefficients of the full model
# --------------------------------------------------------------------------
def fit_final_model(X: pd.DataFrame, y: np.ndarray, groups: np.ndarray,
                    feature_cols: list[str], seed: int = config.SEED):
    """Refit on all data with C chosen by clutch-held-out CV over the full set."""
    best_C = _select_C(X[feature_cols], y, groups, feature_cols,
                       config.N_OUTER_FOLDS, config.C_GRID, seed)
    pipe = make_pipeline(feature_cols, C=best_C, seed=seed).fit(X[feature_cols], y)
    return pipe, best_C


def coefficient_table(pipe: Pipeline, X: pd.DataFrame, y: np.ndarray, groups: np.ndarray,
                      feature_cols: list[str], best_C: float,
                      n_boot: int = config.N_BOOTSTRAP, seed: int = config.SEED) -> pd.DataFrame:
    """Coefficient CIs from a within-clutch stratified bootstrap conditional on the observed clutches."""
    coefs = pipe.named_steps["logreg"].coef_.ravel()
    rng = np.random.default_rng(seed + 77)
    idx_by_group = {g: np.where(groups == g)[0] for g in np.unique(groups)}
    boot = []
    for _ in range(n_boot):
        take = np.concatenate([rng.choice(ix, size=len(ix), replace=True)
                               for ix in idx_by_group.values()])
        yb = y[take]
        if len(np.unique(yb)) < 2:
            continue
        try:
            p = make_pipeline(feature_cols, C=best_C, seed=seed).fit(X[feature_cols].iloc[take], yb)
        except Exception:
            continue
        boot.append(p.named_steps["logreg"].coef_.ravel())
    boot = np.asarray(boot) if boot else np.empty((0, len(coefs)))

    rows = []
    for j, name in enumerate(feature_cols):
        lo, hi = (np.percentile(boot[:, j], [2.5, 97.5]) if len(boot) else (np.nan, np.nan))
        rows.append(
            {
                "predictor": name,
                "coef_standardised": float(coefs[j]),
                "odds_ratio_per_SD": float(np.exp(coefs[j])),
                "ci_low": float(lo),
                "ci_high": float(hi),
                "or_ci_low": float(np.exp(lo)),
                "or_ci_high": float(np.exp(hi)),
                "excludes_zero": bool(np.isfinite(lo) and (lo > 0 or hi < 0)),
            }
        )
    tab = pd.DataFrame(rows)
    tab.attrs["intercept"] = float(pipe.named_steps["logreg"].intercept_[0])
    tab.attrs["n_boot"] = int(len(boot))
    return tab


def confusion_at(y: np.ndarray, scores: np.ndarray, threshold: float = 0.5) -> dict:
    pred = (scores >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()
    sens = tp / (tp + fn) if (tp + fn) else np.nan
    spec = tn / (tn + fp) if (tn + fp) else np.nan
    ppv = tp / (tp + fp) if (tp + fp) else np.nan
    npv = tn / (tn + fn) if (tn + fn) else np.nan
    return {
        "matrix": np.array([[tn, fp], [fn, tp]]),
        "tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp),
        "accuracy": float((tp + tn) / len(y)),
        "sensitivity": float(sens), "specificity": float(spec),
        "ppv": float(ppv), "npv": float(npv),
        "balanced_accuracy": float((sens + spec) / 2),
        "threshold": threshold,
    }


def roc_points(y: np.ndarray, scores: np.ndarray):
    fpr, tpr, _ = roc_curve(y, scores)
    return fpr, tpr


def youden_threshold(y: np.ndarray, scores: np.ndarray) -> float:
    """Threshold maximising Youden's J = sensitivity + specificity - 1.

    Reported alongside the 0.50 threshold because accuracy is threshold-dependent
    and 0.50 is an arbitrary operating point. This threshold is chosen on the
    out-of-fold predictions, so the accuracy it yields is mildly optimistic --
    it is quoted as an operating-point illustration, never as the headline.
    """
    fpr, tpr, thr = roc_curve(y, scores)
    j = tpr - fpr
    return float(thr[int(np.argmax(j))])
