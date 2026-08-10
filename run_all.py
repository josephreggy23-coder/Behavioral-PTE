"""Runner for the larval zebrafish blast-TBI behavioral-analysis pipeline.

    python run_all.py                 # full run
    python run_all.py --permutations 200   # faster smoke run
    python run_all.py --skip-permutation   # skip the permutation test entirely

Every output lands in results/ (figures, tables, all_statistics.csv) plus
RESULTS.md at the repository root.
"""
from __future__ import annotations

import argparse
import random
import sys
import time

import numpy as np

from src import (
    config,
    features,
    habituation,
    io_data,
    plotting,
    report,
    statsbook as sb,
    step2_prediction,
    step3_cfos,
    step4_descriptive,
)


def set_seeds(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    print(f"Random seed set to {seed} (python random, numpy; scikit-learn estimators "
          f"receive random_state={seed}).")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--permutations", type=int, default=config.N_PERMUTATIONS,
                    help=f"permutation-test iterations (default {config.N_PERMUTATIONS})")
    ap.add_argument("--skip-permutation", action="store_true",
                    help="skip the permutation test (sets iterations to 0)")
    ap.add_argument("--seed", type=int, default=config.SEED)
    args = ap.parse_args(argv)

    if args.seed != config.SEED:
        config.SEED = args.seed
    n_perm = 0 if args.skip_permutation else args.permutations
    if n_perm < 0:
        ap.error("--permutations must be non-negative")

    t0 = time.time()
    sb.banner("larval zebrafish blast TBI -- behavioral outcome analysis pipeline")
    set_seeds(config.SEED)
    print(f"Dataset: {config.DATA_XLSX}")

    design = io_data.describe_design()
    print("\nDesign (followed fish_* cohort):")
    print(design.to_string(index=False))
    design.to_csv(config.TABLES / "design_summary.csv", index=False)

    # ---------------------------------------------------------------- step 1
    trials = io_data.habituation_trials()
    fits = habituation.refit_all_sessions(trials)
    merged = habituation.compare_with_supplied(fits)
    habituation.report_step1(merged)
    plotting.fig_curvefit_examples(trials, merged, habituation.habituation_model)
    plotting.fig_tau_agreement(merged)

    # ---------------------------------------------------------------- step 2
    tbl = features.build_fish_table(merged)
    model_df, dropped = features.injured_modeling_set(tbl)
    features.report_features(tbl, model_df, dropped)
    if n_perm == 0:
        print("\n  [note] permutation test disabled (0 iterations)")
    s2 = step2_prediction.run(model_df, n_perm=n_perm)

    # ---------------------------------------------------------------- step 3
    s3 = step3_cfos.run()

    # ---------------------------------------------------------------- step 4
    s4 = step4_descriptive.run(merged)

    # --------------------------------------------------------------- outputs
    sb.banner("Writing outputs")
    stats = sb.flush()
    print(f"  {len(stats)} statistics -> {config.STATS_CSV.relative_to(config.ROOT)}")
    report.write({"step1": merged, "step2": s2, "step3": s3, "step4": s4,
                  "model_df": model_df})

    print(f"\nDone in {time.time() - t0:.1f} s. Seed = {config.SEED}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
