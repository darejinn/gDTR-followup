"""EXP2 §5 — Variant scoring with WGS-wide oscil/c_t normalization.

Uses H3b variant re-forward outputs (B_variants_oscil.parquet) and normalizes
per-variant metrics against WGS-wide distributions computed from all chr's
per-position parquets.

The idea: variant Δoscil could look big (e.g. +3 oscillations added) but is that a
lot relative to the intron/coding_exon background? Normalize into z-scores or
percentiles based on the WGS distribution of the same feature.

Reads:
  exp3_threshold_crossing/results/B_variants_oscil.parquet    # per-variant ref/alt oscil, c_t, etc
  wgs/results/chr{N}/chr{N}_per_position_chunk*.parquet       # position-level features (sampled)
  wgs/data/labels/chr{N}_position_labels.npy                   # context labels
  data_ref/clinvar/clinvar.vcf.gz                              # for MC labels
Writes:
  exp2_variant_downstream/results/wgs_normalized_variant_features.csv
  exp2_variant_downstream/results/H2c_wgs_normalized_scoring.json
  exp2_variant_downstream/figures/EXP2_H2c_wgs_normalized.png
"""
from __future__ import annotations
import gzip, json, time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import f1_score, balanced_accuracy_score, roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TDIG = Path("/NHNHOME/WORKSPACE/0526040123_A/darejinn/tdig")
B_PARQUET   = TDIG / "exp3_threshold_crossing/results/B_variants_oscil.parquet"
CLINVAR_VCF = TDIG / "data_ref/clinvar/clinvar.vcf.gz"
FEATURES_CSV = TDIG / "results_cached/phase3_ensemble/variants_features_full.csv"
WGS = TDIG / "wgs/results"
LABEL_DIR = TDIG / "wgs/data/labels"
OUT_JSON = TDIG / "exp2_variant_downstream/results/H2c_wgs_normalized_scoring.json"
OUT_CSV = TDIG / "exp2_variant_downstream/results/wgs_normalized_variant_features.csv"
FIG_DIR = TDIG / "exp2_variant_downstream/figures"

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
SEED = 42
SAMPLE_PER_CHR = 200_000  # sample positions per chr for background distribution

LABEL_MAP = {0: "intergenic", 1: "intron", 2: "coding_exon", 3: "5utr", 4: "3utr",
             5: "splice_donor", 6: "splice_acceptor"}

# ClinVar MC → 4-class
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


def collect_wgs_background(sample_n=SAMPLE_PER_CHR) -> dict[str, dict[str, tuple[float, float]]]:
    """For each context, compute WGS-wide mean and SD of c_t, oscil, n_enter, n_exit, below_frac, min_D.
    Uses a random sample of positions per chromosome to keep memory manageable.
    Returns: {context: {feature: (mean, sd)}}."""
    print("=== Computing WGS-wide background distributions per context ===")
    per_ctx_samples = {ctx: {"c_t": [], "oscil": [], "n_enter": [], "n_exit": [],
                              "below_frac": [], "min_D": []} for ctx in LABEL_MAP.values()}
    rng = np.random.default_rng(SEED)
    for chrom in CHROMS:
        t0 = time.time()
        chunks = sorted((WGS / chrom).glob(f"{chrom}_per_position_chunk*.parquet"))
        lab = np.load(LABEL_DIR / f"{chrom}_position_labels.npy")
        # Sample uniformly across chunks
        n_per_chunk = max(1, sample_n // len(chunks))
        for cp in chunks:
            d = pd.read_parquet(cp)
            if len(d) > n_per_chunk:
                idx = rng.choice(len(d), size=n_per_chunk, replace=False)
                d = d.iloc[idx]
            d["label"] = lab[d["pos"].to_numpy()]
            d["context"] = d["label"].map(LABEL_MAP).fillna("unknown")
            for ctx in LABEL_MAP.values():
                sub = d[d["context"] == ctx]
                if len(sub) == 0: continue
                per_ctx_samples[ctx]["c_t"].extend(sub["c_t"].to_numpy())
                per_ctx_samples[ctx]["oscil"].extend(sub["oscil"].to_numpy())
                per_ctx_samples[ctx]["n_enter"].extend(sub["n_enter"].to_numpy())
                per_ctx_samples[ctx]["n_exit"].extend(sub["n_exit"].to_numpy())
                per_ctx_samples[ctx]["below_frac"].extend(sub["below_frac"].to_numpy())
                per_ctx_samples[ctx]["min_D"].extend(sub["min_D"].to_numpy())
        print(f"[{time.strftime('%H:%M:%S')}] {chrom} sampled in {time.time()-t0:.1f}s")

    background = {}
    for ctx, feats in per_ctx_samples.items():
        background[ctx] = {}
        for f, vals in feats.items():
            if not vals: continue
            arr = np.asarray(vals, dtype=np.float32)
            background[ctx][f] = (float(arr.mean()), float(arr.std(ddof=1) or 1e-9))
    return background


def main():
    t0 = time.time()

    # ---- 1. Load B_variants_oscil ----
    b = pd.read_parquet(B_PARQUET)
    b["chrom_norm"] = b["chrom"].astype(str).str.replace("^chr", "", regex=True)
    print(f"B_variants: {b.shape}")

    # ---- 2. ClinVar MC join ----
    cv = parse_clinvar(CLINVAR_VCF, set(b["chrom_norm"].unique()))
    cv["chrom_norm"] = cv["chrom"].astype(str)
    b = b.merge(cv[["chrom_norm", "pos", "ref", "alt", "MC"]],
                on=["chrom_norm", "pos", "ref", "alt"], how="left")
    b["consequence"] = b["MC"].apply(clinvar_class)
    conseq_counts = b["consequence"].value_counts().to_dict()
    print(f"consequence counts: {conseq_counts}")

    # ---- 3. WGS background ----
    background = collect_wgs_background()
    print(f"[{time.strftime('%H:%M:%S')}] background loaded, contexts: {list(background.keys())}")

    # Save background for reproducibility
    bg_json = {ctx: {f: {"mean": m, "sd": s} for f, (m, s) in feats.items()}
               for ctx, feats in background.items()}
    (TDIG / "wgs/results/genome_summary/wgs_variant_background.json").write_text(
        json.dumps(bg_json, indent=2, default=str))

    # ---- 4. Determine each variant's context from consequence ----
    # For scoring purposes we use intron background for intron variants, coding_exon for missense/nonsense/synonymous, splice_donor/acceptor for canonical_splice
    ctx_for_conseq = {
        "missense": "coding_exon", "nonsense": "coding_exon", "synonymous": "coding_exon",
        "canonical_splice": "splice_donor",  # approximate; could split
        "intron": "intron", "frameshift": "coding_exon", "other": "intron",
        "unknown": "intron",
    }
    b["ctx_for_norm"] = b["consequence"].map(ctx_for_conseq).fillna("intron")

    # ---- 5. WGS-normalized z-scores for ref, alt, delta ----
    for tag in ("ref", "alt"):
        for f in ("c_t", "oscil", "below_frac", "min_D"):
            col_val = b[f"{f}_{tag}"] if f != "below_frac" else b[f"below_frac_{tag}"]
            b[f"z_{f}_{tag}"] = [
                (val - background[ctx][f][0]) / background[ctx][f][1]
                for val, ctx in zip(col_val, b["ctx_for_norm"])
            ]
    # delta z-scores
    for f in ("c_t", "oscil", "below_frac", "min_D"):
        b[f"z_d_{f}"] = b[f"z_{f}_alt"] - b[f"z_{f}_ref"]

    # save
    b.to_csv(OUT_CSV, index=False)
    print(f"wrote {OUT_CSV} shape {b.shape}")

    # ---- 6. Merge with cos32 features from paper cohort ----
    fdf = pd.read_csv(FEATURES_CSV)
    fdf["chrom_norm"] = fdf["chrom"].astype(str).str.replace("^chr", "", regex=True)
    is_snv = (fdf["ref"].str.len() == 1) & (fdf["alt"].str.len() == 1)
    fdf = fdf[is_snv & fdf["category"].isin(["P_LP", "B_LB"])].copy()
    merged = fdf.merge(b, on=["chrom_norm", "pos", "ref", "alt"], how="inner")
    print(f"merged with cos32 cohort: {merged.shape}")

    # 4-class subtype filter
    keep = ["missense", "nonsense", "synonymous", "canonical_splice"]
    dc = merged[merged["consequence"].isin(keep)].copy()
    cos_cols = [f"dD_cos_{i}" for i in range(32)]
    dc = dc.dropna(subset=cos_cols + ["consequence"]).copy()
    y = dc["consequence"].to_numpy()
    print(f"final 4-class cohort: {dc.shape}, counts: {dc['consequence'].value_counts().to_dict()}")

    # ---- 7. Compare configurations ----
    X_cos = dc[cos_cols].to_numpy(dtype=np.float32)
    X_raw_exp3 = dc[["d_oscil", "d_n_enter", "d_n_exit", "d_below_frac", "d_c_t", "d_min_D",
                     "oscil_ref", "oscil_alt", "min_D_ref", "min_D_alt"]].to_numpy(dtype=np.float32)
    X_norm_exp3 = dc[["z_d_c_t", "z_d_oscil", "z_d_below_frac", "z_d_min_D",
                       "z_c_t_ref", "z_c_t_alt", "z_oscil_ref", "z_oscil_alt",
                       "z_min_D_ref", "z_min_D_alt"]].to_numpy(dtype=np.float32)
    X_scalars = dc[["max_abs_dD", "signed_argmax", "argmax_layer", "dc_interp",
                     "evo2_delta_loglik"]].to_numpy(dtype=np.float32)

    configs = {
        "cos32_only (paper baseline)": X_cos,
        "cos32 + raw EXP3": np.hstack([X_cos, X_raw_exp3]),
        "cos32 + WGS-normalized EXP3": np.hstack([X_cos, X_norm_exp3]),
        "cos32 + scalars": np.hstack([X_cos, X_scalars]),
        "cos32 + scalars + WGS-normalized": np.hstack([X_cos, X_scalars, X_norm_exp3]),
        "cos32 + scalars + raw EXP3 + WGS-normalized": np.hstack([X_cos, X_scalars, X_raw_exp3, X_norm_exp3]),
        "WGS-normalized only": X_norm_exp3,
        "raw EXP3 only": X_raw_exp3,
    }

    def run_multiclass(X, y):
        skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=SEED)
        pipe = Pipeline([("s", StandardScaler()),
                          ("c", LogisticRegression(max_iter=2000, random_state=SEED))])
        y_pred = cross_val_predict(pipe, X, y, cv=skf)
        proba = cross_val_predict(pipe, X, y, cv=skf, method="predict_proba")
        classes = sorted(np.unique(y))
        return {
            "macro_f1": float(f1_score(y, y_pred, average="macro")),
            "bal_acc": float(balanced_accuracy_score(y, y_pred)),
            "ovr_auroc": {c: float(roc_auc_score((y == c).astype(int), proba[:, i]))
                           for i, c in enumerate(classes)},
            "n_features": int(X.shape[1]),
        }

    print("\n=== Multi-class subtype classification ===")
    results = {}
    for name, X in configs.items():
        r = run_multiclass(X, y)
        results[name] = r
        print(f"  {name:50s} n_feat={r['n_features']:>3} macro-F1={r['macro_f1']:.4f} bal-acc={r['bal_acc']:.4f}")

    # ---- 8. Binary P/LP vs B/LB ----
    # After merge, category column has _x/_y suffix if it appeared in both frames.
    cat_col = "category_x" if "category_x" in merged.columns else "category"
    dc_bin = merged[merged[cat_col].isin(["P_LP", "B_LB"])].copy()
    dc_bin = dc_bin.dropna(subset=cos_cols).copy()
    yb = (dc_bin[cat_col] == "P_LP").astype(int).to_numpy()
    Xb_cos = dc_bin[cos_cols].to_numpy(dtype=np.float32)
    Xb_raw = dc_bin[["d_oscil", "d_n_enter", "d_n_exit", "d_below_frac", "d_c_t", "d_min_D"]].to_numpy(dtype=np.float32)
    Xb_norm = dc_bin[["z_d_c_t", "z_d_oscil", "z_d_below_frac", "z_d_min_D"]].to_numpy(dtype=np.float32)

    def run_binary(X, y):
        skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=SEED)
        pipe = Pipeline([("s", StandardScaler()),
                          ("c", LogisticRegression(max_iter=2000, random_state=SEED))])
        proba = cross_val_predict(pipe, X, y, cv=skf, method="predict_proba")[:, 1]
        return float(roc_auc_score(y, proba))

    binary = {}
    binary["cos32"] = run_binary(Xb_cos, yb)
    binary["cos32 + raw EXP3"] = run_binary(np.hstack([Xb_cos, Xb_raw]), yb)
    binary["cos32 + WGS-normalized"] = run_binary(np.hstack([Xb_cos, Xb_norm]), yb)
    binary["cos32 + both"] = run_binary(np.hstack([Xb_cos, Xb_raw, Xb_norm]), yb)
    binary["raw EXP3 only"] = run_binary(Xb_raw, yb)
    binary["WGS-normalized only"] = run_binary(Xb_norm, yb)
    print("\n=== Binary P/LP vs B/LB (AUROC) ===")
    for k, v in binary.items():
        print(f"  {k:35s} = {v:.4f}")

    # ---- 9. Save ----
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps({
        "cohort_shape": list(dc.shape),
        "wgs_background_summary_written_to": str(TDIG / "wgs/results/genome_summary/wgs_variant_background.json"),
        "multiclass": results,
        "binary_auroc": binary,
        "runtime_sec": time.time() - t0,
    }, indent=2, default=str))

    # ---- 10. Figure ----
    plt.rcParams.update({"font.family": "DejaVu Sans", "savefig.dpi": 300, "savefig.bbox": "tight",
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(15, 5))

    # (a) multi-class macro-F1
    ax = axes[0]
    order = ["raw EXP3 only", "WGS-normalized only",
             "cos32_only (paper baseline)",
             "cos32 + raw EXP3", "cos32 + WGS-normalized EXP3",
             "cos32 + scalars", "cos32 + scalars + WGS-normalized",
             "cos32 + scalars + raw EXP3 + WGS-normalized"]
    macro = [results[k]["macro_f1"] for k in order]
    n_feat = [results[k]["n_features"] for k in order]
    colors = ["#8bc34a", "#4dd0e1", "#0d47a1", "#1976d2", "#00acc1", "#7b1fa2", "#c2185b", "#c62828"]
    y_pos = np.arange(len(order))
    ax.barh(y_pos, macro, color=colors, edgecolor="black", lw=0.5)
    for yi, (m, n) in enumerate(zip(macro, n_feat)):
        ax.text(m + 0.005, yi, f"{m:.3f}   (n_feat={n})", va="center", fontsize=9, fontweight="bold")
    ax.axvline(0.25, ls="--", color="#555", label="chance (4-class)")
    ax.axvline(results["cos32_only (paper baseline)"]["macro_f1"], ls=":", color="#0d47a1",
               label=f"paper baseline {results['cos32_only (paper baseline)']['macro_f1']:.3f}")
    ax.set_yticks(y_pos)
    ax.set_yticklabels([o.replace("cos32", "cos32").replace(" + ", "\n+ ") for o in order], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0, 0.85)
    ax.set_xlabel("macro-F1 (10-fold, seed 42)")
    ax.set_title("(a) 4-class subtype classification with WGS-normalized features",
                 loc="left", fontweight="bold", fontfamily="Times New Roman", fontsize=12, pad=8)
    ax.legend(loc="lower right", fontsize=8)
    ax.grid(axis="x", color="#dddddd", lw=0.5)

    # (b) binary P/LP AUROC
    ax = axes[1]
    order_b = ["raw EXP3 only", "WGS-normalized only", "cos32",
               "cos32 + raw EXP3", "cos32 + WGS-normalized", "cos32 + both"]
    v_b = [binary[k] for k in order_b]
    y_pos = np.arange(len(order_b))
    ax.barh(y_pos, v_b, color=["#8bc34a", "#4dd0e1", "#0d47a1", "#1976d2", "#00acc1", "#c62828"],
            edgecolor="black", lw=0.5)
    for yi, v in enumerate(v_b):
        ax.text(v + 0.005, yi, f"{v:.4f}", va="center", fontsize=9, fontweight="bold")
    ax.axvline(0.5, ls="--", color="#555")
    ax.axvline(binary["cos32"], ls=":", color="#0d47a1", label=f"cos32 = {binary['cos32']:.3f}")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(order_b)
    ax.invert_yaxis()
    ax.set_xlim(0.4, 0.9)
    ax.set_xlabel("AUROC (P/LP vs B/LB, 10-fold)")
    ax.set_title("(b) Binary pathogenicity",
                 loc="left", fontweight="bold", fontfamily="Times New Roman", fontsize=12, pad=8)
    ax.legend(loc="lower right", fontsize=8)
    ax.grid(axis="x", color="#dddddd", lw=0.5)

    fig.suptitle("EXP2 §5 — Variant scoring with WGS-wide oscil/c_t normalization",
                 fontsize=13, fontweight="bold", fontfamily="Times New Roman", y=1.02)
    fig.tight_layout()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / "EXP2_H2c_wgs_normalized.png")
    fig.savefig(FIG_DIR / "EXP2_H2c_wgs_normalized.pdf")
    plt.close(fig)
    print(f"wrote {FIG_DIR/'EXP2_H2c_wgs_normalized.png'}")

    print(f"\nTotal runtime: {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
