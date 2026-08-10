from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import run_all
from src import features, modeling, step2_prediction, step3_cfos


def test_skip_permutation_cli_passes_zero_to_step2(monkeypatch, tmp_path):
    captured: dict[str, int] = {}
    design = pd.DataFrame(columns=["group", "clutch", "n_fish", "n_converted"])
    sentinel = object()

    monkeypatch.setattr(run_all.config, "TABLES", tmp_path)
    monkeypatch.setattr(run_all.io_data, "describe_design", lambda: design)
    monkeypatch.setattr(run_all.io_data, "habituation_trials", lambda: sentinel)
    monkeypatch.setattr(run_all.habituation, "refit_all_sessions", lambda _: sentinel)
    monkeypatch.setattr(run_all.habituation, "compare_with_supplied", lambda _: sentinel)
    monkeypatch.setattr(run_all.habituation, "report_step1", lambda _: None)
    monkeypatch.setattr(run_all.plotting, "fig_curvefit_examples", lambda *args: None)
    monkeypatch.setattr(run_all.plotting, "fig_tau_agreement", lambda *args: None)
    monkeypatch.setattr(run_all.features, "build_fish_table", lambda _: sentinel)
    monkeypatch.setattr(
        run_all.features, "injured_modeling_set", lambda _: (sentinel, sentinel)
    )
    monkeypatch.setattr(run_all.features, "report_features", lambda *args: None)

    def fake_step2(model_df, *, n_perm):
        assert model_df is sentinel
        captured["n_perm"] = n_perm
        return {}

    monkeypatch.setattr(run_all.step2_prediction, "run", fake_step2)
    monkeypatch.setattr(run_all.step3_cfos, "run", lambda: {})
    monkeypatch.setattr(run_all.step4_descriptive, "run", lambda _: {})
    monkeypatch.setattr(run_all.sb, "flush", lambda: pd.DataFrame())
    monkeypatch.setattr(run_all.report, "write", lambda _: "")

    assert run_all.main(["--skip-permutation"]) == 0
    assert captured["n_perm"] == 0


def test_zero_permutation_request_short_circuits_parallel(monkeypatch):
    def fail_parallel(*args, **kwargs):  # pragma: no cover - called only on regression
        raise AssertionError("Parallel permutation work must not start for n_perm=0")

    monkeypatch.setattr(modeling, "Parallel", fail_parallel)
    result = modeling.permutation_test(
        pd.DataFrame(),
        np.array([], dtype=int),
        np.array([], dtype=object),
        [],
        observed_pooled=0.75,
        observed_mean_fold=0.70,
        n_perm=0,
    )

    assert result["performed"] is False
    assert result["n_perm"] == 0
    assert result["p_pooled"] is None
    assert result["p_mean_fold"] is None
    assert result["null_pooled"].size == 0


def test_step2_skip_neither_runs_nor_plots_permutation(monkeypatch, tmp_path):
    expected_seed = 8675309
    observed_seeds = {"nested": [], "bootstrap": [], "final": [], "coefficient": []}
    model_df = pd.DataFrame(
        {
            "fish_id": [f"fish_{i}" for i in range(6)],
            "group": ["low_impact", "low_impact", "high_impact"] * 2,
            "clutch": ["clutch_A"] * 2 + ["clutch_B"] * 2 + ["clutch_C"] * 2,
            "converted": [0, 1, 0, 1, 0, 1],
            "dose": [0, 0, 1, 1, 0, 1],
            "pre_tau": [1.0, 1.2, 0.9, 1.4, 1.1, 1.3],
            "dtau_0.5": [-0.2, 0.2, -0.1, 0.3, -0.3, 0.4],
            "dtau_24": [-0.1, 0.1, -0.2, 0.2, -0.4, 0.3],
            "baseline_locomotion_pre": [10, 11, 12, 13, 14, 15],
        }
    )
    y = model_df["converted"].to_numpy()
    oof = np.where(y == 1, 0.7, 0.3)
    folds = pd.DataFrame(
        {
            "fold": [0, 1, 2],
            "held_out": ["clutch_A", "clutch_B", "clutch_C"],
            "n_test": [2, 2, 2],
            "n_events_test": [1, 1, 1],
            "selected_C": [1.0, 1.0, 1.0],
            "auc": [0.75, 0.75, 0.75],
        }
    )
    nested_result = {
        "oof": oof,
        "folds": folds,
        "mean_fold_auc": 0.75,
        "sd_fold_auc": 0.0,
        "min_fold_auc": 0.75,
        "max_fold_auc": 0.75,
        "pooled_auc": 0.75,
        "brier": 0.09,
    }
    coefs = pd.DataFrame(
        {
            "predictor": features.MODEL_SPECS["e_full"]["features"],
            "coef_standardised": [0.1, 0.2, 0.3, 0.4],
            "odds_ratio_per_SD": np.exp([0.1, 0.2, 0.3, 0.4]),
            "ci_low": [-0.2, -0.1, 0.0, 0.1],
            "ci_high": [0.4, 0.5, 0.6, 0.7],
            "or_ci_low": np.exp([-0.2, -0.1, 0.0, 0.1]),
            "or_ci_high": np.exp([0.4, 0.5, 0.6, 0.7]),
            "excludes_zero": [False, False, False, True],
        }
    )
    coefs.attrs.update(intercept=0.0, n_boot=10)
    cm = {
        "matrix": np.array([[3, 0], [0, 3]]),
        "threshold": 0.5,
        "tn": 3,
        "fp": 0,
        "fn": 0,
        "tp": 3,
        "accuracy": 1.0,
        "balanced_accuracy": 1.0,
        "sensitivity": 1.0,
        "specificity": 1.0,
        "ppv": 1.0,
        "npv": 1.0,
    }

    def fail(*args, **kwargs):  # pragma: no cover - called only on regression
        raise AssertionError("permutation test/plot must not run")

    monkeypatch.setattr(step2_prediction.config, "TABLES", tmp_path)
    monkeypatch.setattr(step2_prediction.config, "SEED", expected_seed)
    monkeypatch.setattr(step2_prediction, "_assumption_checks", lambda _: None)

    def fake_nested(*args, seed, **kwargs):
        observed_seeds["nested"].append(seed)
        return nested_result

    def fake_bootstrap(*args, seed, **kwargs):
        observed_seeds["bootstrap"].append(seed)
        return (0.6, 0.9)

    def fake_final(*args, seed, **kwargs):
        observed_seeds["final"].append(seed)
        return (object(), 1.0)

    def fake_coefficient(*args, seed, **kwargs):
        observed_seeds["coefficient"].append(seed)
        return coefs

    monkeypatch.setattr(modeling, "nested_cv", fake_nested)
    monkeypatch.setattr(modeling, "bootstrap_auc_ci", fake_bootstrap)
    monkeypatch.setattr(modeling, "roc_points", lambda *args: (np.array([0, 1]), np.array([0, 1])))
    monkeypatch.setattr(modeling, "permutation_test", fail)
    monkeypatch.setattr(modeling, "fit_final_model", fake_final)
    monkeypatch.setattr(modeling, "coefficient_table", fake_coefficient)
    monkeypatch.setattr(modeling, "confusion_at", lambda *args: cm)
    monkeypatch.setattr(modeling, "youden_threshold", lambda *args: 0.5)
    monkeypatch.setattr(step2_prediction.sb, "banner", lambda *args: None)
    monkeypatch.setattr(step2_prediction.sb, "record", lambda *args, **kwargs: None)
    monkeypatch.setattr(step2_prediction.sb, "check", lambda *args, **kwargs: True)
    monkeypatch.setattr(step2_prediction.plotting, "fig_roc_comparison", lambda *args: None)
    monkeypatch.setattr(step2_prediction.plotting, "fig_permutation_null", fail)
    monkeypatch.setattr(
        step2_prediction.plotting, "fig_confusion_and_calibration", lambda *args: None
    )
    monkeypatch.setattr(step2_prediction.plotting, "fig_coefficients", lambda *args: None)
    monkeypatch.setattr(step2_prediction.plotting, "fig_fold_auc", lambda *args: None)

    result = step2_prediction.run(model_df, n_perm=0)

    assert result["perm"]["performed"] is False
    assert result["perm"]["p_pooled"] is None
    assert observed_seeds["nested"] and set(observed_seeds["nested"]) == {expected_seed}
    assert observed_seeds["bootstrap"] and set(observed_seeds["bootstrap"]) == {expected_seed}
    assert observed_seeds["final"] == [expected_seed]
    assert observed_seeds["coefficient"] == [expected_seed]


def test_delta_tau_only_ablation_is_genuinely_dose_blind():
    spec = features.MODEL_SPECS["e_dtau_only"]

    assert spec["features"] == ["dtau_0.5", "dtau_24"]
    assert "dose" not in spec["features"]
    assert features.MODEL_SPECS["f_behavior_no_dose"]["features"] == [
        "pre_tau",
        "dtau_0.5",
        "dtau_24",
    ]
    assert "dose" not in features.MODEL_SPECS["f_behavior_no_dose"]["features"]
    assert features.MODEL_SPECS["e_full"]["features"] == [
        "dose",
        "pre_tau",
        "dtau_0.5",
        "dtau_24",
    ]

    transformer = modeling.make_pipeline(spec["features"]).named_steps["withindose_z"]
    transformer.fit(pd.DataFrame({"dtau_0.5": [1.0, 2.0], "dtau_24": [2.0, 4.0]}))
    assert transformer.active_ is False


def test_cfos_pool_membership_is_normalized_and_validated():
    pools = pd.DataFrame(
        {
            "pool_id": ["pool_1", "pool_2"],
            "risk_pool": ["high_risk", "low_risk"],
            "group": ["low_impact", "low_impact"],
            "clutch": ["clutch_A", "clutch_A"],
            "pooled_fish_ids": ["f1; f2;f3;f4", "f5;f6;f7;f8"],
        }
    )

    membership = step3_cfos.build_pool_membership(pools)

    assert membership.columns.tolist() == step3_cfos.POOL_MEMBERSHIP_COLUMNS
    assert membership["fish_id"].tolist() == [f"f{i}" for i in range(1, 9)]
    assert membership.groupby("pool_id")["fish_id"].nunique().to_dict() == {
        "pool_1": 4,
        "pool_2": 4,
    }

    duplicate_within_pool = pools.iloc[[0]].copy()
    duplicate_within_pool.loc[:, "pooled_fish_ids"] = "f1;f1;f2;f3"
    with pytest.raises(ValueError, match="exactly 4 unique"):
        step3_cfos.build_pool_membership(duplicate_within_pool)

    reused_between_pools = pools.copy()
    reused_between_pools.loc[1, "pooled_fish_ids"] = "f4;f6;f7;f8"
    with pytest.raises(ValueError, match="reused"):
        step3_cfos.build_pool_membership(reused_between_pools)


def test_cfos_clutch_dependency_sensitivity_is_saved_and_returned(monkeypatch, tmp_path):
    raw_diffs = {
        "clutch_A": {"sham": 0.05, "low_impact": 0.20, "high_impact": 0.35},
        "clutch_B": {"sham": 0.10, "low_impact": 0.30, "high_impact": 0.50},
        "clutch_C": {"sham": 0.15, "low_impact": 0.40, "high_impact": 0.65},
    }
    pool_rows = []
    fish_number = 1
    for clutch in step3_cfos.config.CLUTCHES:
        for group in step3_cfos.config.GROUPS:
            low_value = 1.0 + 0.1 * fish_number
            for risk_pool, value in (
                ("low_risk", low_value),
                ("high_risk", low_value + raw_diffs[clutch][group]),
            ):
                fish_ids = [f"cfos_{i}" for i in range(fish_number, fish_number + 4)]
                fish_number += 4
                pool_rows.append(
                    {
                        "pool_id": f"{group}_{clutch}_{risk_pool}",
                        "risk_pool": risk_pool,
                        "group": group,
                        "clutch": clutch,
                        "pooled_fish_ids": ";".join(fish_ids),
                        "n_larvae_in_pool": 4,
                        "cfos_fold_change": value,
                        "delta_ddct": -np.log2(value),
                    }
                )
    pools = pd.DataFrame(pool_rows)
    records = []

    monkeypatch.setattr(step3_cfos.config, "TABLES", tmp_path)
    monkeypatch.setattr(step3_cfos.sb, "banner", lambda *args: None)
    monkeypatch.setattr(step3_cfos.sb, "check", lambda *args, **kwargs: True)
    monkeypatch.setattr(
        step3_cfos.sb,
        "record",
        lambda *args, **kwargs: records.append((args, kwargs)),
    )
    monkeypatch.setattr(step3_cfos.plotting, "fig_cfos_paired", lambda *args: None)

    result = step3_cfos.run(pools)

    table = result["clutch_sensitivity_table"]
    assert table["n_group_pairs"].tolist() == [3, 3, 3]
    assert table["mean_raw_high_minus_low"].to_numpy() == pytest.approx([0.20, 0.30, 0.40])
    assert result["clutch_sensitivity_raw"]["n"] == 3
    assert result["clutch_sensitivity_log2"]["n"] == 3
    assert "all" in result and result["all"]["n"] == 9

    saved = pd.read_csv(tmp_path / "step3_cfos_clutch_aggregated_sensitivity.csv")
    pd.testing.assert_frame_equal(saved, table, check_dtype=False)
    scopes = {args[1] for args, _ in records if len(args) >= 3}
    assert "cfos_clutch_aggregated" in scopes
    assert "cfos_clutch_aggregated_log2" in scopes
