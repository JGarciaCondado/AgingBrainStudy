"""Training metrics for the brain-age model of a cohort.

Reports:
  - Number of training samples (cognitively normal participants, cn == 1).
  - Slope and intercept of the linear age-bias correction.

The age-bias correction used by ageml is
    corrected = predicted + (age - (slope * age + intercept))
so for any subject
    predicted - corrected + age = slope * age + intercept
which is exactly linear for not_cn subjects (the model fitted on the CN).
We recover slope and intercept by regressing that quantity on age.

Usage (from the repository root):
    python scripts/training_metrics.py a4
    python scripts/training_metrics.py adni
    python scripts/training_metrics.py habs
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

COHORTS = ["a4", "adni", "habs"]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("cohort", choices=COHORTS, help="Cohort name")
    parser.add_argument("--root", default=".",
                        help="Repository root (default: current directory)")
    args = parser.parse_args()

    clinical_path = os.path.join(args.root, "data", "final", args.cohort,
                                 "processed", "clinical.csv")
    pred_path = os.path.join(args.root, "analysis", "final", args.cohort,
                             "ageml", "model_age", "predicted_age.csv")
    for p in (clinical_path, pred_path):
        if not os.path.exists(p):
            sys.exit(f"File not found: {p}")

    clinical = pd.read_csv(clinical_path, index_col="ID")
    pred = pd.read_csv(pred_path, index_col="ID")

    # Column is "cn" in A4/HABS and "CN" in ADNI
    cn_col = next((c for c in clinical.columns if c.lower() == "cn"), None)
    if cn_col is None:
        sys.exit(f"No cn/CN column in {clinical_path}")

    # --- Number of training samples ---
    cn_ids = clinical.index[clinical[cn_col] == 1]
    n_cn = len(cn_ids)

    # --- Age-bias slope and intercept (from not_cn subjects) ---
    # HABS has no not_cn column, so use everyone who is not cn
    not_cn_ids = clinical.index[clinical[cn_col] == 0]
    sub = pred.loc[pred.index.isin(not_cn_ids)].dropna(
        subset=["age", "predicted_age_all", "corrected_age_all"])
    if len(sub) < 2:
        sys.exit("Not enough not_cn subjects in predicted_age.csv to fit the correction.")

    age = sub["age"].to_numpy()
    y = (sub["predicted_age_all"] - sub["corrected_age_all"] + sub["age"]).to_numpy()
    slope, intercept = np.polyfit(age, y, 1)

    # Sanity check: reconstruct corrected age and report max error
    recon = sub["predicted_age_all"] + (age - (slope * age + intercept))
    max_err = float(np.max(np.abs(recon - sub["corrected_age_all"])))

    print(f"Cohort: {args.cohort}")
    print(f"Training samples (cn in clinical.csv): {n_cn}")
    print(f"Age-bias correction (fit on {len(sub)} not_cn subjects):")
    print(f"  predicted_age = {slope:.6f} * age + {intercept:.6f}")
    print(f"  slope:     {slope:.6f}")
    print(f"  intercept: {intercept:.6f}")
    print(f"  max reconstruction error: {max_err:.2e} years")


if __name__ == "__main__":
    main()
