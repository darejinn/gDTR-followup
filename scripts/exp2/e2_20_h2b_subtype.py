"""EXP2 §2b — H2b: variant subtype classification from ΔD_cos.

Joins molecular consequence (MC) from ClinVar VCF onto our 8,008 SNVs, then
multinomial LR to predict 5-class consequence from 32-d ΔD_cos vector.

Hypothesis: ΔD_cos vectors carry per-variant subtype information beyond
population-level shift shown in paper §3.3.

Outputs:
  results/H2b_subtype_multi.json          headline numbers
  results/H2b_confusion_matrix.csv        confusion matrix
  results/H2b_dim_ablation.csv            32-d vs 1-d (max_abs_dD) comparison
  results/H2b_consequence_join_report.md  join provenance + counts
"""
from __future__ import annotations

import gzip
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import (
    classification_report, confusion_matrix, f1_score, roc_auc_score,
    balanced_accuracy_score,
)
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

TDIG = Path("/NHNHOME/WORKSPACE/0526040123_A/darejinn/tdig")
FEATURES_CSV = TDIG / "results_cached/phase3_ensemble/variants_features_full.csv"
CLINVAR_VCF = TDIG / "data_ref/clinvar/clinvar.vcf.gz"
OUT_DIR = TDIG / "exp2_variant_downstream/results"

SEED = 42

# ClinVar Molecular Consequence (MC) → paper Fig 3 6-way ordering
# Paper classes: intron, frameshift, nonsense, missense, canonical_splice, synonymous
# ClinVar MC values (Sequence Ontology):
#   SO:0001583 → missense_variant
#   SO:0001587 → stop_gained (= nonsense)
#   SO:0001589 → frameshift_variant
#   SO:0001819 → synonymous_variant
#   SO:0001629 → splice_acceptor_variant
#   SO:0001575 → splice_donor_variant
#   SO:0001627 → intron_variant
MC_TO_CLASS = {
    "SO:0001583|missense_variant": "missense",
    "SO:0001587|nonsense": "nonsense",
    "SO:0001587|stop_gained": "nonsense",
    "SO:0001589|frameshift_variant": "frameshift",
    "SO:0001819|synonymous_variant": "synonymous",
    "SO:0001629|splice_acceptor_variant": "canonical_splice",
    "SO:0001575|splice_donor_variant": "canonical_splice",
    "SO:0001627|intron_variant": "intron",
    "SO:0001623|5_prime_UTR_variant": "5utr",
    "SO:0001624|3_prime_UTR_variant": "3utr",
    "SO:0001060|sequence_variant": "other",
    "SO:0001628|intergenic_variant": "intergenic",
    "SO:0001891|inframe_indel": "inframe",
}

CLASS_ORDER = ["missense", "nonsense", "synonymous", "canonical_splice", "frameshift", "intron"]


def parse_clinvar_mc(vcf_path: Path, wanted_chroms: set[str]) -> pd.DataFrame:
    """Parse ClinVar VCF and extract MC per variant for the requested chromosomes."""
    rows = []
    print(f"[{time.strftime('%H:%M:%S')}] parsing ClinVar VCF from {vcf_path}")
    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, _, ref, alt, _, _, info = fields[:8]
            if chrom not in wanted_chroms:
                continue
            # MC in INFO: MC=SO:0001583|missense_variant,SO:0001627|intron_variant
            mc = None
            for kv in info.split(";"):
                if kv.startswith("MC="):
                    mc = kv[3:]
                    break
            rows.append({"chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt, "MC": mc})
    print(f"[{time.strftime('%H:%M:%S')}] parsed {len(rows)} variants from {vcf_path.name}")
    return pd.DataFrame(rows)


def clinvar_class(mc: str | None) -> str:
    """Map MC string to our class label; take first match with highest priority."""
    if not isinstance(mc, str) or not mc:
        return "unknown"
    priorities = ["frameshift", "nonsense", "canonical_splice", "missense",
                  "synonymous", "5utr", "3utr", "intron", "intergenic", "inframe", "other"]
    matched = set()
    for kv in mc.split(","):
        cls = MC_TO_CLASS.get(kv, None)
        if cls:
            matched.add(cls)
    for p in priorities:
        if p in matched:
            return p
    return "unknown"


def main() -> None:
    t0 = time.time()

    print(f"[{time.strftime('%H:%M:%S')}] loading features …")
    df = pd.read_csv(FEATURES_CSV)
    df["chrom_norm"] = df["chrom"].astype(str).str.replace("^chr", "", regex=True)
    is_snv = (df["ref"].str.len() == 1) & (df["alt"].str.len() == 1)
    d = df[is_snv & df["category"].isin(["P_LP", "B_LB"])].copy()
    chroms = set(d["chrom_norm"].unique())
    print(f"[{time.strftime('%H:%M:%S')}] SNV+category filter: {d.shape}, chroms: {sorted(chroms)}")

    print(f"[{time.strftime('%H:%M:%S')}] parsing ClinVar for chroms {chroms}")
    cv = parse_clinvar_mc(CLINVAR_VCF, chroms)
    print(f"[{time.strftime('%H:%M:%S')}] ClinVar rows for these chroms: {cv.shape}")

    cv["chrom_norm"] = cv["chrom"].astype(str)
    d = d.merge(cv[["chrom_norm", "pos", "ref", "alt", "MC"]],
                on=["chrom_norm", "pos", "ref", "alt"], how="left")
    n_with_mc = d["MC"].notna().sum()
    print(f"[{time.strftime('%H:%M:%S')}] merge: {d.shape}, with MC: {n_with_mc}/{len(d)}")

    d["consequence"] = d["MC"].apply(clinvar_class)
    conseq_counts = d["consequence"].value_counts().to_dict()
    print(f"[{time.strftime('%H:%M:%S')}] consequence counts: {conseq_counts}")

    # Filter to 5 core classes: paper Fig 3 (drop intron for now, dropped in paper too)
    keep = ["missense", "nonsense", "synonymous", "canonical_splice"]
    dc = d[d["consequence"].isin(keep)].copy()
    print(f"[{time.strftime('%H:%M:%S')}] cohort (4-class): {dc.shape}, class counts: {dc['consequence'].value_counts().to_dict()}")

    if len(dc) < 100:
        print("insufficient rows for classification — aborting")
        return

    cos_cols = [f"dD_cos_{i}" for i in range(32)]
    dc = dc.dropna(subset=cos_cols + ["consequence"]).copy()
    X_full = dc[cos_cols].to_numpy(dtype=np.float32)
    X_1d = dc["max_abs_dD"].to_numpy(dtype=np.float32).reshape(-1, 1)
    y = dc["consequence"].to_numpy()
    class_labels = sorted(np.unique(y))
    print(f"[{time.strftime('%H:%M:%S')}] final X: {X_full.shape}, classes: {class_labels}")

    # Stratified 10-fold multinomial LR on 32-d
    skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=SEED)
    y_pred_full = cross_val_predict(
        Pipeline([("s", StandardScaler()),
                  ("c", LogisticRegression(max_iter=2000, random_state=SEED))]),
        X_full, y, cv=skf,
    )
    macro_f1_32d = f1_score(y, y_pred_full, average="macro")
    bal_acc_32d = balanced_accuracy_score(y, y_pred_full)
    cm_32d = confusion_matrix(y, y_pred_full, labels=class_labels)

    # 1-d comparison
    y_pred_1d = cross_val_predict(
        Pipeline([("s", StandardScaler()),
                  ("c", LogisticRegression(max_iter=2000, random_state=SEED))]),
        X_1d, y, cv=skf,
    )
    macro_f1_1d = f1_score(y, y_pred_1d, average="macro")
    bal_acc_1d = balanced_accuracy_score(y, y_pred_1d)

    # OvR AUROC per class (32-d)
    y_score_full = cross_val_predict(
        Pipeline([("s", StandardScaler()),
                  ("c", LogisticRegression(max_iter=2000, random_state=SEED))]),
        X_full, y, cv=skf, method="predict_proba",
    )
    auroc_per_class = {}
    for i, cls in enumerate(class_labels):
        y_bin = (y == cls).astype(int)
        auroc_per_class[cls] = float(roc_auc_score(y_bin, y_score_full[:, i]))

    # Save confusion matrix as CSV
    cm_df = pd.DataFrame(cm_32d, index=class_labels, columns=class_labels)
    cm_df.to_csv(OUT_DIR / "H2b_confusion_matrix.csv")

    summary = {
        "cohort_full_shape": list(dc.shape),
        "class_counts": {c: int((y == c).sum()) for c in class_labels},
        "auroc_per_class_OvR_32d": auroc_per_class,
        "macro_f1_32d": float(macro_f1_32d),
        "balanced_acc_32d": float(bal_acc_32d),
        "macro_f1_1d_max_abs_dD": float(macro_f1_1d),
        "balanced_acc_1d_max_abs_dD": float(bal_acc_1d),
        "class_order": class_labels,
        "confusion_matrix_32d": cm_32d.tolist(),
        "chance_macro_f1": 1.0 / len(class_labels),
        "verdict": {
            "H2b_positive_32d": bool(macro_f1_32d > 0.35),
            "32d_beats_1d_scalar": macro_f1_32d - macro_f1_1d > 0.02,
        },
        "runtime_sec": time.time() - t0,
    }
    (OUT_DIR / "H2b_subtype_multi.json").write_text(json.dumps(summary, indent=2))

    # Dim ablation CSV
    pd.DataFrame([
        {"feature": "32d_dD_cos_vector", "macro_f1": macro_f1_32d, "balanced_acc": bal_acc_32d, "n_features": 32},
        {"feature": "1d_max_abs_dD", "macro_f1": macro_f1_1d, "balanced_acc": bal_acc_1d, "n_features": 1},
    ]).to_csv(OUT_DIR / "H2b_dim_ablation.csv", index=False)

    # Join report
    with (OUT_DIR / "H2b_consequence_join_report.md").open("w") as f:
        f.write("# H2b consequence join report\n\n")
        f.write(f"Source VCF: `{CLINVAR_VCF}`\n\n")
        f.write(f"Chromosomes queried: {sorted(chroms)}\n\n")
        f.write(f"Merged rows with MC: {n_with_mc} / {len(d)}\n\n")
        f.write(f"Consequence counts (all): {conseq_counts}\n\n")
        f.write(f"Final 4-class cohort: {dict(dc['consequence'].value_counts())}\n")

    print()
    print(f"=== H2b summary ===")
    print(f"Cohort: {dc.shape}, {len(class_labels)} classes")
    for cls, n in {c: int((y == c).sum()) for c in class_labels}.items():
        print(f"  {cls:20s}: n={n}")
    print()
    print(f"32-d ΔD_cos vector:")
    print(f"  macro-F1 = {macro_f1_32d:.4f}  (chance {1.0/len(class_labels):.4f})")
    print(f"  bal. acc = {bal_acc_32d:.4f}")
    print(f"  OvR AUROC per class: {json.dumps(auroc_per_class, indent=4)}")
    print()
    print(f"1-d max_abs_dD scalar:")
    print(f"  macro-F1 = {macro_f1_1d:.4f}")
    print(f"  bal. acc = {bal_acc_1d:.4f}")
    print()
    print(f"32-d gain over 1-d: macro-F1 +{macro_f1_32d - macro_f1_1d:+.4f}")
    print()
    if macro_f1_32d > 0.35:
        print(f"✓ H2b POSITIVE — 32-d vector distinguishes variant subtypes (macro-F1 > 0.35)")
    else:
        print(f"✗ H2b NEGATIVE — cannot reliably distinguish subtypes from ΔD_cos alone")
    print(f"Runtime: {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
