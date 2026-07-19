"""EXP2 §4 — Integrated H2b: does adding Δoscil / Δn_enter / Δn_exit to the 32-d ΔD_cos
vector improve the 4-class variant subtype classification (macro-F1 baseline = 0.641)?

Reads:
  results_cached/phase3_ensemble/variants_features_full.csv  (32-d ΔD_cos + max_abs_dD + Evo2 LL)
  exp3_threshold_crossing/results/B_variants_oscil.parquet  (per-variant Δoscil, Δn_enter, Δn_exit,
                                                              Δbelow_frac, Δc_t, Δmin_D)
  ClinVar MC join (already worked out in e2_20_h2b_subtype.py)

Writes:
  exp2_variant_downstream/results/H2b_integrated_features.json
  exp2_variant_downstream/results/H2b_integrated_confusion_matrix.csv
  exp2_variant_downstream/figures/EXP2_H2b_integrated.png/pdf
"""
from __future__ import annotations
import gzip, json, time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import f1_score, confusion_matrix, roc_auc_score, balanced_accuracy_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

TDIG = Path("/NHNHOME/WORKSPACE/0526040123_A/darejinn/tdig")
FEATURES_CSV = TDIG / "results_cached/phase3_ensemble/variants_features_full.csv"
B_PARQUET   = TDIG / "exp3_threshold_crossing/results/B_variants_oscil.parquet"
CLINVAR_VCF = TDIG / "data_ref/clinvar/clinvar.vcf.gz"
OUT_JSON    = TDIG / "exp2_variant_downstream/results/H2b_integrated_features.json"
SEED = 42


# Reuse consequence mapping from e2_20_h2b_subtype.py
MC_TO_CLASS = {
    "SO:0001583|missense_variant": "missense",
    "SO:0001587|nonsense": "nonsense",
    "SO:0001587|stop_gained": "nonsense",
    "SO:0001589|frameshift_variant": "frameshift",
    "SO:0001819|synonymous_variant": "synonymous",
    "SO:0001629|splice_acceptor_variant": "canonical_splice",
    "SO:0001575|splice_donor_variant": "canonical_splice",
    "SO:0001627|intron_variant": "intron",
}


def parse_clinvar(vcf_path, chroms):
    rows = []
    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"): continue
            fs = line.rstrip("\n").split("\t")
            if len(fs) < 8 or fs[0] not in chroms: continue
            info = fs[7]
            mc = next((kv[3:] for kv in info.split(";") if kv.startswith("MC=")), None)
            rows.append({"chrom": fs[0], "pos": int(fs[1]), "ref": fs[3], "alt": fs[4], "MC": mc})
    return pd.DataFrame(rows)


def clinvar_class(mc):
    if not isinstance(mc, str) or not mc: return "unknown"
    priorities = ["frameshift", "nonsense", "canonical_splice", "missense", "synonymous", "intron", "other"]
    matched = {MC_TO_CLASS.get(kv, None) for kv in mc.split(",")} - {None}
    for p in priorities:
        if p in matched: return p
    return "unknown"


def run_cv_multiclass(X, y, seed=SEED):
    skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=seed)
    y_pred = cross_val_predict(
        Pipeline([("s", StandardScaler()), ("c", LogisticRegression(max_iter=2000, random_state=seed))]),
        X, y, cv=skf,
    )
    proba = cross_val_predict(
        Pipeline([("s", StandardScaler()), ("c", LogisticRegression(max_iter=2000, random_state=seed))]),
        X, y, cv=skf, method="predict_proba",
    )
    classes = sorted(np.unique(y))
    macro_f1 = f1_score(y, y_pred, average="macro")
    bal_acc = balanced_accuracy_score(y, y_pred)
    cm = confusion_matrix(y, y_pred, labels=classes)
    aucs = {c: float(roc_auc_score((y == c).astype(int), proba[:, i])) for i, c in enumerate(classes)}
    return {"macro_f1": float(macro_f1), "bal_acc": float(bal_acc),
            "confusion": cm.tolist(), "classes": classes, "ovr_auroc": aucs}


def main():
    t0 = time.time()
    df = pd.read_csv(FEATURES_CSV)
    df["chrom_norm"] = df["chrom"].astype(str).str.replace("^chr", "", regex=True)
    is_snv = (df["ref"].str.len() == 1) & (df["alt"].str.len() == 1)
    d = df[is_snv & df["category"].isin(["P_LP", "B_LB"])].copy()

    # Join ClinVar MC
    cv = parse_clinvar(CLINVAR_VCF, set(d["chrom_norm"].unique()))
    cv["chrom_norm"] = cv["chrom"].astype(str)
    d = d.merge(cv[["chrom_norm", "pos", "ref", "alt", "MC"]],
                on=["chrom_norm", "pos", "ref", "alt"], how="left")
    d["consequence"] = d["MC"].apply(clinvar_class)

    # Join EXP3 variant oscil features
    if not B_PARQUET.exists():
        raise SystemExit(f"H3b output not found at {B_PARQUET} — run e3_40_h3b first")
    b = pd.read_parquet(B_PARQUET)
    b["chrom_norm"] = b["chrom"].astype(str).str.replace("^chr", "", regex=True)
    print(f"H3b parquet: {b.shape}")
    d = d.merge(
        b[["chrom_norm", "pos", "ref", "alt",
           "oscil_ref", "oscil_alt", "d_oscil",
           "n_enter_ref", "n_enter_alt", "d_n_enter",
           "n_exit_ref", "n_exit_alt", "d_n_exit",
           "below_frac_ref", "below_frac_alt", "d_below_frac",
           "c_t_ref", "c_t_alt", "d_c_t",
           "min_D_ref", "min_D_alt", "d_min_D",
           "argmin_layer_ref", "argmin_layer_alt"]],
        on=["chrom_norm", "pos", "ref", "alt"], how="left")

    n_with_oscil = d["d_oscil"].notna().sum()
    print(f"variants with EXP3 features: {n_with_oscil}/{len(d)}")

    keep = ["missense", "nonsense", "synonymous", "canonical_splice"]
    dc = d[d["consequence"].isin(keep) & d["d_oscil"].notna()].copy()
    cos_cols = [f"dD_cos_{i}" for i in range(32)]
    dc = dc.dropna(subset=cos_cols + ["consequence"]).copy()
    y = dc["consequence"].to_numpy()
    print(f"final cohort: {dc.shape}, class counts: {dc['consequence'].value_counts().to_dict()}")

    # Feature blocks
    F_cos32 = dc[cos_cols].to_numpy(dtype=np.float32)
    F_exp3  = dc[["d_oscil", "d_n_enter", "d_n_exit", "d_below_frac", "d_c_t", "d_min_D",
                   "oscil_ref", "oscil_alt", "min_D_ref", "min_D_alt"]].to_numpy(dtype=np.float32)
    F_scalars = dc[["max_abs_dD", "signed_argmax", "argmax_layer", "dc_interp",
                     "evo2_delta_loglik"]].to_numpy(dtype=np.float32)

    configs = {
        "cos32_only (paper baseline)": F_cos32,
        "cos32 + EXP3_features": np.hstack([F_cos32, F_exp3]),
        "cos32 + scalars": np.hstack([F_cos32, F_scalars]),
        "cos32 + EXP3 + scalars (all)": np.hstack([F_cos32, F_exp3, F_scalars]),
        "EXP3_features only (H3b-derived)": F_exp3,
        "scalars only": F_scalars,
    }

    results = {}
    for name, X in configs.items():
        r = run_cv_multiclass(X, y)
        results[name] = {"n_features": int(X.shape[1]), **{k: v for k, v in r.items() if k != "confusion"}}
        print(f"  {name:40s} n_feat={X.shape[1]:>3d} macro-F1={r['macro_f1']:.4f} bal-acc={r['bal_acc']:.4f}")

    summary = {
        "cohort_shape": list(dc.shape),
        "class_counts": {c: int((y == c).sum()) for c in results[list(results.keys())[0]]["classes"]},
        "baseline_paper": 0.641,
        "results": results,
        "runtime_sec": time.time() - t0,
    }
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(summary, indent=2, default=str))
    print()
    baseline = results["cos32_only (paper baseline)"]["macro_f1"]
    integrated = results["cos32 + EXP3_features"]["macro_f1"]
    all_feats = results["cos32 + EXP3 + scalars (all)"]["macro_f1"]
    print(f"BASELINE (32-d ΔD_cos only):        macro-F1 = {baseline:.4f}")
    print(f"+ EXP3 features (10 new):            macro-F1 = {integrated:.4f}  (Δ = {integrated - baseline:+.4f})")
    print(f"+ EXP3 + scalars (all features):     macro-F1 = {all_feats:.4f}   (Δ = {all_feats - baseline:+.4f})")
    print(f"\nWrote {OUT_JSON}")


if __name__ == "__main__":
    main()
