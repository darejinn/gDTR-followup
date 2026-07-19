"""EXP2 §2.0 — verify the paper's App C AUROC = 0.844 reproduces from the cached features.

CLAUDE.md rule: reproduce paper number before running any new experiment on the same data.

Inputs:
  /NHNHOME/WORKSPACE/0526040123_A/darejinn/tdig/results_cached/phase3_ensemble/variants_features_full.csv
Outputs:
  exp2_variant_downstream/results/00_verify_paper.json
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

TDIG = Path("/NHNHOME/WORKSPACE/0526040123_A/darejinn/tdig")
FEATURES_CSV = TDIG / "results_cached/phase3_ensemble/variants_features_full.csv"
OUT_JSON = TDIG / "exp2_variant_downstream/results/00_verify_paper.json"

SEED = 42


def main() -> None:
    t0 = time.time()
    df = pd.read_csv(FEATURES_CSV)
    print(f"loaded features: shape={df.shape}")

    # Paper §App C uses SNVs (single-base substitutions) across 15 cancer genes,
    # dropping frameshifts. Cohort size reported: 8,008 P/LP vs B/LB SNVs.
    is_snv = (df["ref"].str.len() == 1) & (df["alt"].str.len() == 1)
    print(f"SNV filter: {is_snv.sum()} / {len(df)} rows")
    d = df.loc[is_snv].copy()

    # Category: P_LP (positive) vs B_LB (negative). Drop any other category.
    d = d[d["category"].isin(["P_LP", "B_LB"])].copy()
    print(f"category filter: {d.shape}")
    print("class counts:", d["category"].value_counts().to_dict())

    d["y"] = (d["category"] == "P_LP").astype(int)

    # 32-d ΔD_cos features
    cos_cols = [f"dD_cos_{i}" for i in range(32)]
    d = d.dropna(subset=cos_cols + ["y"]).copy()
    print(f"after dropna on cos_cols: {d.shape}")

    X = d[cos_cols].to_numpy(dtype=np.float32)
    y = d["y"].to_numpy(dtype=np.int32)
    print(f"X shape: {X.shape}, y positive rate: {y.mean():.4f}, n_pos={y.sum()}, n_neg={(y==0).sum()}")

    # Stratified 10-fold CV, seed 42 (paper convention).
    skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=SEED)
    fold_aurocs = []
    for i, (tr, te) in enumerate(skf.split(X, y)):
        pipe = Pipeline([
            ("scaler", StandardScaler()),
            ("clf", LogisticRegression(max_iter=1000, random_state=SEED, C=1.0)),
        ])
        pipe.fit(X[tr], y[tr])
        proba = pipe.predict_proba(X[te])[:, 1]
        auc = roc_auc_score(y[te], proba)
        fold_aurocs.append(auc)
        print(f"  fold {i}: AUROC = {auc:.4f}")

    mean_auroc = float(np.mean(fold_aurocs))
    std_auroc = float(np.std(fold_aurocs, ddof=1))
    n_boot = 1000
    rng = np.random.default_rng(SEED)
    boot = np.asarray([rng.choice(fold_aurocs, size=len(fold_aurocs), replace=True).mean() for _ in range(n_boot)])
    ci_lo, ci_hi = float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))

    # Also try best single-layer for reference (paper reports 0.729 at L=30 for cos)
    print()
    print("Best single-layer ΔD_cos AUROC:")
    single_layer_aurocs = {}
    for l, col in enumerate(cos_cols):
        aucs = []
        for tr, te in skf.split(X, y):
            pipe = Pipeline([
                ("scaler", StandardScaler()),
                ("clf", LogisticRegression(max_iter=1000, random_state=SEED)),
            ])
            pipe.fit(X[tr, l:l+1], y[tr])
            proba = pipe.predict_proba(X[te, l:l+1])[:, 1]
            aucs.append(roc_auc_score(y[te], proba))
        single_layer_aurocs[l] = float(np.mean(aucs))
    best_layer = max(single_layer_aurocs, key=single_layer_aurocs.get)
    best_single_auroc = single_layer_aurocs[best_layer]
    print(f"  best layer = {best_layer}, AUROC = {best_single_auroc:.4f}   (paper: 0.729 @ L=30)")

    out = {
        "cohort": {"n": int(len(d)), "n_pos": int(y.sum()), "n_neg": int((y == 0).sum())},
        "auroc_32d_cos": {
            "mean": mean_auroc,
            "std": std_auroc,
            "ci95": [ci_lo, ci_hi],
            "per_fold": [float(a) for a in fold_aurocs],
            "n_folds": 10,
            "seed": SEED,
            "paper_target": 0.844,
        },
        "best_single_layer_cos": {
            "layer": int(best_layer),
            "auroc": best_single_auroc,
            "paper_target": 0.729,
            "paper_target_layer": 30,
        },
        "per_layer_aurocs_cos": {str(l): a for l, a in single_layer_aurocs.items()},
        "runtime_sec": time.time() - t0,
    }
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(out, indent=2))
    print()
    print(f"32-d ΔD_cos AUROC = {mean_auroc:.4f} ± {std_auroc:.4f} (95% CI [{ci_lo:.4f}, {ci_hi:.4f}])")
    print(f"Paper §App C:       0.844 [0.831, 0.857]")
    print(f"Best single-layer:  {best_single_auroc:.4f} at L={best_layer}   (paper: 0.729 @ L=30)")
    print(f"Runtime: {time.time() - t0:.1f} s")
    print(f"Wrote: {OUT_JSON}")

    # Gate: paper reproduction OK if within 0.02 of paper value
    if abs(mean_auroc - 0.844) <= 0.02:
        print("\n✓ Verification passed — safe to proceed with EXP2 hypotheses.")
    else:
        print("\n✗ Verification FAILED — investigate before running EXP2.")


if __name__ == "__main__":
    main()
