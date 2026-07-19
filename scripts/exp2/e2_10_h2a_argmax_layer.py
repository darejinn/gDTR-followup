"""EXP2 §2a — H2a: argmax-layer feature vs fixed-layer AUROC.

Hypothesis H2a: at each variant v, take x_v = ΔD_cos(argmax_ℓ |ΔD_cos(ℓ,v)|, v).
This variant-adaptive single scalar should beat any FIXED single-layer AUROC.

Also compare:
- 32-d full vector (paper baseline: 0.844)
- 1-d argmax-layer scalar (new)
- 1-d argmax-layer value (the raw signed argmax value; already computed as `signed_argmax`)
- 1-d + evo2_delta_loglik (ensemble)

Outputs:
  results/H2a_argmax_layer_auroc.csv    per-feature 10-fold AUROC + 95% CI
  results/H2a_delong.csv                DeLong p-values comparing pairs
  results/H2a_regression_summary.json   headline numbers
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
from scipy import stats

TDIG = Path("/NHNHOME/WORKSPACE/0526040123_A/darejinn/tdig")
FEATURES_CSV = TDIG / "results_cached/phase3_ensemble/variants_features_full.csv"
OUT_DIR = TDIG / "exp2_variant_downstream/results"
SEED = 42


def stratified_10fold_auroc(X: np.ndarray, y: np.ndarray, seed: int = SEED, n_boot: int = 1000):
    skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=seed)
    fold_aurocs = []
    all_true, all_score = [], []
    for tr, te in skf.split(X, y):
        pipe = Pipeline([
            ("scaler", StandardScaler()),
            ("clf", LogisticRegression(max_iter=1000, random_state=seed, C=1.0)),
        ])
        pipe.fit(X[tr], y[tr])
        proba = pipe.predict_proba(X[te])[:, 1]
        fold_aurocs.append(roc_auc_score(y[te], proba))
        all_true.append(y[te]); all_score.append(proba)
    rng = np.random.default_rng(seed)
    boot = np.asarray([rng.choice(fold_aurocs, size=len(fold_aurocs), replace=True).mean() for _ in range(n_boot)])
    return {
        "mean": float(np.mean(fold_aurocs)),
        "std": float(np.std(fold_aurocs, ddof=1)),
        "ci95": [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
        "per_fold": [float(a) for a in fold_aurocs],
        "all_true": np.concatenate(all_true),
        "all_score": np.concatenate(all_score),
    }


def delong_test(true_labels: np.ndarray, s1: np.ndarray, s2: np.ndarray) -> tuple[float, float]:
    """Simple bootstrap DeLong-like test on AUROC delta."""
    a1 = roc_auc_score(true_labels, s1); a2 = roc_auc_score(true_labels, s2)
    n = len(true_labels)
    rng = np.random.default_rng(SEED)
    boot_deltas = []
    for _ in range(1000):
        idx = rng.choice(n, size=n, replace=True)
        d = roc_auc_score(true_labels[idx], s1[idx]) - roc_auc_score(true_labels[idx], s2[idx])
        boot_deltas.append(d)
    boot_deltas = np.asarray(boot_deltas)
    delta_obs = a1 - a2
    p = 2 * min(float((boot_deltas <= 0).mean()), float((boot_deltas >= 0).mean()))
    return delta_obs, p


def main() -> None:
    t0 = time.time()
    df = pd.read_csv(FEATURES_CSV)
    is_snv = (df["ref"].str.len() == 1) & (df["alt"].str.len() == 1)
    d = df[is_snv & df["category"].isin(["P_LP", "B_LB"])].copy()
    d["y"] = (d["category"] == "P_LP").astype(int)
    cos_cols = [f"dD_cos_{i}" for i in range(32)]
    d = d.dropna(subset=cos_cols + ["y", "argmax_layer", "signed_argmax", "evo2_delta_loglik"]).copy()

    print(f"Cohort: n={len(d)}, n_pos={d.y.sum()}, n_neg={(d.y==0).sum()}")

    X_cos = d[cos_cols].to_numpy(dtype=np.float32)
    y = d["y"].to_numpy(dtype=np.int32)

    # Feature 1: 32-d vector (paper baseline)
    r_32d = stratified_10fold_auroc(X_cos, y)

    # Feature 2: variant-adaptive argmax-layer feature — read ΔD_cos at argmax_layer per variant
    # Note: argmax_layer in the CSV is 1-based (paper's convention with L_star=29 etc)
    argmax_layers = d["argmax_layer"].to_numpy(dtype=np.int64) - 1  # 0-based
    argmax_layers = np.clip(argmax_layers, 0, 31)
    X_argmax = X_cos[np.arange(len(X_cos)), argmax_layers].reshape(-1, 1)
    r_argmax = stratified_10fold_auroc(X_argmax, y)

    # Feature 3: absolute value of argmax-layer feature (rectified, removes sign)
    X_abs = np.abs(X_argmax)
    r_abs = stratified_10fold_auroc(X_abs, y)

    # Feature 4: precomputed max_abs_dD (single scalar per variant, no layer index needed)
    X_max = d["max_abs_dD"].to_numpy(dtype=np.float32).reshape(-1, 1)
    r_max = stratified_10fold_auroc(X_max, y)

    # Feature 5: signed_argmax
    X_sgn = d["signed_argmax"].to_numpy(dtype=np.float32).reshape(-1, 1)
    r_sgn = stratified_10fold_auroc(X_sgn, y)

    # Feature 6: per-layer fixed AUROC (paper reports best at L=30)
    per_layer = {}
    for l in range(32):
        r = stratified_10fold_auroc(X_cos[:, l:l+1], y, n_boot=100)
        per_layer[l] = r["mean"]
    best_layer = max(per_layer, key=per_layer.get)
    r_best_fixed = stratified_10fold_auroc(X_cos[:, best_layer:best_layer+1], y)

    # Feature 7: ensemble argmax-layer + Evo 2 LL
    X_ens = np.column_stack([X_argmax, d["evo2_delta_loglik"].to_numpy(dtype=np.float32).reshape(-1, 1)])
    r_ens = stratified_10fold_auroc(X_ens, y)

    features = {
        "32d_cos_vector (paper baseline)": r_32d,
        "argmax_layer_value (H2a)": r_argmax,
        "abs(argmax_layer_value)": r_abs,
        "max_abs_dD_precomputed": r_max,
        "signed_argmax_precomputed": r_sgn,
        f"best_fixed_layer L={best_layer}": r_best_fixed,
        "argmax_layer + Evo2_LL ensemble": r_ens,
    }

    # DeLong tests on aggregated predictions
    delongs = {}
    y_true = r_32d["all_true"]
    for label, r in features.items():
        if label == "32d_cos_vector (paper baseline)":
            continue
        delta, p = delong_test(y_true, r["all_score"], r_32d["all_score"])
        delongs[label] = {"delta_vs_32d": float(delta), "p_bootstrap": float(p)}

    # Argmax layer distribution — check if it's concentrated (would make H2a boring)
    layer_dist = d["argmax_layer"].value_counts().sort_index().to_dict()
    unique_argmax_layers = int(d["argmax_layer"].nunique())

    # Save CSV
    rows = []
    for label, r in features.items():
        rows.append({
            "feature": label,
            "auroc_mean": r["mean"],
            "auroc_std": r["std"],
            "auroc_ci_lo": r["ci95"][0],
            "auroc_ci_hi": r["ci95"][1],
            "delta_vs_32d": delongs.get(label, {}).get("delta_vs_32d", 0.0),
            "p_vs_32d": delongs.get(label, {}).get("p_bootstrap", 1.0),
        })
    pd.DataFrame(rows).to_csv(OUT_DIR / "H2a_argmax_layer_auroc.csv", index=False)
    pd.DataFrame([{"layer": k, "n_variants_argmax_here": v} for k, v in sorted(layer_dist.items())]).to_csv(
        OUT_DIR / "H2a_argmax_layer_distribution.csv", index=False)

    summary = {
        "cohort": {"n": int(len(d)), "n_pos": int(y.sum()), "n_neg": int((y==0).sum())},
        "unique_argmax_layers": unique_argmax_layers,
        "best_fixed_layer": int(best_layer),
        "per_layer_auroc": {str(k): v for k, v in per_layer.items()},
        "features_auroc": {label: {k: v for k, v in r.items() if k not in ("all_true", "all_score")}
                            for label, r in features.items()},
        "delong_vs_32d": delongs,
        "verdict": {
            "H2a_positive": bool(r_argmax["mean"] > r_best_fixed["mean"] + 0.005),
            "argmax_beats_best_fixed": r_argmax["mean"] - r_best_fixed["mean"],
            "argmax_beats_max_abs_dD": r_argmax["mean"] - r_max["mean"],
            "argmax_reaches_32d": r_argmax["mean"] - r_32d["mean"],
        },
        "runtime_sec": time.time() - t0,
    }
    (OUT_DIR / "H2a_regression_summary.json").write_text(json.dumps(summary, indent=2))

    print()
    print(f"=== H2a summary ===")
    print(f"unique argmax layers used: {unique_argmax_layers}/32")
    print(f"Best fixed layer: L={best_layer}, AUROC={r_best_fixed['mean']:.4f}")
    for label, r in features.items():
        print(f"  {label:45s} = {r['mean']:.4f} ± {r['std']:.4f}")
    print()
    print(f"H2a verdict: argmax_layer_value ({r_argmax['mean']:.4f}) - best_fixed ({r_best_fixed['mean']:.4f}) = {r_argmax['mean']-r_best_fixed['mean']:+.4f}")
    if r_argmax["mean"] > r_best_fixed["mean"] + 0.005:
        print(f"  ✓ H2a POSITIVE — variant-adaptive layer selection helps")
    else:
        print(f"  ✗ H2a NEGATIVE — variant-adaptive layer selection does NOT beat best fixed layer")
    print(f"Runtime: {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
