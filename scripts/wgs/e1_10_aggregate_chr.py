"""EXP1 §Phase B — aggregate all chr{N} per-position parquet chunks + join with labels →
per-context summary, comparable to paper §3.1 Table.

Reads: wgs/results/{chrom}/{chrom}_per_position_chunk*.parquet
       wgs/data/labels/{chrom}_position_labels.npy
Writes: wgs/results/{chrom}/{chrom}_context_summary.csv
        wgs/results/{chrom}/{chrom}_context_summary.json
        + figure comparing per-context c_t and oscil to paper baseline
"""
from __future__ import annotations
import argparse, json, time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TDIG = Path("/NHNHOME/WORKSPACE/0526040123_A/darejinn/tdig")

LABEL_MAP = {0: "intergenic", 1: "intron", 2: "coding_exon", 3: "5utr", 4: "3utr",
             5: "splice_donor", 6: "splice_acceptor"}

# Paper §3.1 numbers (pooled chr17+chr22) for comparison
PAPER_C = {"intron": 27.72, "splice_donor": 25.55, "splice_acceptor": 25.96,
           "3utr": 27.74, "coding_exon": 28.40, "intergenic": 28.66, "5utr": 29.22}


def cohen_d_indep(x, y):
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2: return float("nan")
    vx, vy = float(np.var(x, ddof=1)), float(np.var(y, ddof=1))
    pooled = np.sqrt(((nx-1)*vx + (ny-1)*vy) / (nx+ny-2))
    return float((np.mean(x)-np.mean(y))/pooled) if pooled else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True)
    args = ap.parse_args()

    IN_DIR = TDIG / f"wgs/results/{args.chrom}"
    LAB_PATH = TDIG / f"wgs/data/labels/{args.chrom}_position_labels.npy"
    OUT_CSV = IN_DIR / f"{args.chrom}_context_summary.csv"
    OUT_JSON = IN_DIR / f"{args.chrom}_context_summary.json"
    FIG_PATH = IN_DIR / f"{args.chrom}_context_summary.png"

    chunks = sorted(IN_DIR.glob(f"{args.chrom}_per_position_chunk*.parquet"))
    print(f"[{time.strftime('%H:%M:%S')}] {args.chrom}: {len(chunks)} chunks")
    if not chunks:
        raise SystemExit("no chunks")

    print(f"[{time.strftime('%H:%M:%S')}] loading labels …")
    labels = np.load(LAB_PATH)
    print(f"[{time.strftime('%H:%M:%S')}]   labels shape: {labels.shape}, dist: {dict(zip(*np.unique(labels, return_counts=True)))}")

    # Stream aggregation to avoid holding all rows
    df_parts = []
    total_rows = 0
    t0 = time.time()
    for i, cp in enumerate(chunks):
        d = pd.read_parquet(cp, columns=["pos", "c_t", "oscil", "n_enter", "n_exit", "below_frac", "min_D", "argmin_layer"])
        # Attach label from position
        d["label"] = labels[d["pos"].to_numpy()]
        df_parts.append(d)
        total_rows += len(d)
        if (i + 1) % 20 == 0:
            print(f"[{time.strftime('%H:%M:%S')}]   loaded {i+1}/{len(chunks)} chunks, {total_rows:,} rows, {time.time()-t0:.1f}s")
    df = pd.concat(df_parts, ignore_index=True)
    print(f"[{time.strftime('%H:%M:%S')}] total rows: {len(df):,}, load time {time.time()-t0:.1f}s")
    del df_parts
    df["context"] = df["label"].map(LABEL_MAP).fillna("unknown")

    # Per-context summary
    summary_rows = []
    intron_ct = df[df["context"] == "intron"]["c_t"].to_numpy(dtype=float)
    intron_osc = df[df["context"] == "intron"]["oscil"].to_numpy(dtype=float)
    for ctx in ["intron", "splice_donor", "splice_acceptor", "coding_exon", "3utr", "5utr", "intergenic"]:
        sub = df[df["context"] == ctx]
        if len(sub) < 30:
            continue
        ct = sub["c_t"].to_numpy(dtype=float)
        osc = sub["oscil"].to_numpy(dtype=float)
        d_ct = cohen_d_indep(ct, intron_ct) if ctx != "intron" else 0.0
        d_osc = cohen_d_indep(osc, intron_osc) if ctx != "intron" else 0.0
        try:
            _, p_ct = stats.mannwhitneyu(ct, intron_ct, alternative="two-sided") if ctx != "intron" else (0.0, 1.0)
            _, p_osc = stats.mannwhitneyu(osc, intron_osc, alternative="two-sided") if ctx != "intron" else (0.0, 1.0)
        except ValueError:
            p_ct, p_osc = float("nan"), float("nan")
        summary_rows.append({
            "context": ctx,
            "n": int(len(sub)),
            "mean_c_t": float(ct.mean()),
            "median_c_t": float(np.median(ct)),
            "cohen_d_c_t_vs_intron": float(d_ct),
            "mwu_p_c_t": float(p_ct) if ctx != "intron" else None,
            "mean_oscil": float(osc.mean()),
            "cohen_d_oscil_vs_intron": float(d_osc),
            "mwu_p_oscil": float(p_osc) if ctx != "intron" else None,
            "mean_below_frac": float(sub["below_frac"].mean()),
            "mean_min_D": float(sub["min_D"].mean()),
            "paper_c_pooled_chr17_22": PAPER_C.get(ctx, None),
        })
    sdf = pd.DataFrame(summary_rows)
    sdf.to_csv(OUT_CSV, index=False)
    OUT_JSON.write_text(json.dumps(summary_rows, indent=2, default=str))
    print(f"\n=== {args.chrom} per-context summary (baseline = intron) ===")
    print(sdf.to_string(index=False))

    # Figure — 2 panels: (a) mean c_t vs paper, (b) Cohen's d on both axes
    plt.rcParams.update({"font.family": "DejaVu Sans", "savefig.dpi": 300, "savefig.bbox": "tight",
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))

    ax = axes[0]
    contexts = sdf["context"].tolist()
    x = np.arange(len(contexts))
    w = 0.36
    ax.bar(x - w/2, sdf["mean_c_t"], w, color="#0d47a1", edgecolor="black", lw=0.4, label=f"this run ({args.chrom})")
    paper_vals = [PAPER_C.get(c, np.nan) for c in contexts]
    ax.bar(x + w/2, paper_vals, w, color="#d62728", alpha=0.55, edgecolor="black", lw=0.4, label="paper §3.1 (chr17+chr22)")
    for xi, (v, n) in enumerate(zip(sdf["mean_c_t"], sdf["n"])):
        ax.text(xi - w/2, v + 0.1, f"{v:.2f}\n(n={n:,})", ha="center", fontsize=7, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([c.replace("_", "\n") for c in contexts], rotation=0)
    ax.axhline(27.72, ls="--", color="#555", lw=0.7, label="paper intron baseline (27.72)")
    ax.set_ylabel(r"mean settling depth $\bar c(t)$")
    ax.set_title(f"(a) {args.chrom} per-context $\\bar c(t)$ vs paper (chr17+chr22 pooled)",
                 loc="left", fontsize=12, fontweight="bold", fontfamily="Times New Roman", pad=8)
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(axis="y", color="#dddddd", lw=0.5)

    ax = axes[1]
    ax.barh(x - w/2, sdf["cohen_d_c_t_vs_intron"], w, color="#d62728", edgecolor="black", lw=0.4, label="$c(t)$ vs intron")
    ax.barh(x + w/2, sdf["cohen_d_oscil_vs_intron"], w, color="#1f77b4", edgecolor="black", lw=0.4, label="oscil vs intron")
    for xi, (d_ct, d_osc) in enumerate(zip(sdf["cohen_d_c_t_vs_intron"], sdf["cohen_d_oscil_vs_intron"])):
        ax.text(d_ct + (0.02 if d_ct >= 0 else -0.02), xi - w/2, f"{d_ct:+.3f}",
                va="center", ha="left" if d_ct >= 0 else "right", fontsize=8, color="#d62728", fontweight="bold")
        ax.text(d_osc + (0.02 if d_osc >= 0 else -0.02), xi + w/2, f"{d_osc:+.3f}",
                va="center", ha="left" if d_osc >= 0 else "right", fontsize=8, color="#1f77b4", fontweight="bold")
    ax.axvline(0, color="#555", lw=0.6)
    ax.set_yticks(x)
    ax.set_yticklabels(contexts)
    ax.invert_yaxis()
    ax.set_xlabel("Cohen's d")
    ax.set_title(f"(b) {args.chrom} Cohen's d: c(t) vs oscil (orthogonal axes)",
                 loc="left", fontsize=12, fontweight="bold", fontfamily="Times New Roman", pad=8)
    ax.legend(loc="lower right")
    ax.grid(axis="x", color="#dddddd", lw=0.5)

    fig.tight_layout()
    fig.savefig(FIG_PATH)
    fig.savefig(str(FIG_PATH).replace(".png", ".pdf"))
    plt.close(fig)
    print(f"wrote {FIG_PATH}")


if __name__ == "__main__":
    main()
