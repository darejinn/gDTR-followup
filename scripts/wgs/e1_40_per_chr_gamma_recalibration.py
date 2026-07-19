"""EXP1 §Phase D — Per-chromosome γ recalibration (paper transfer claim quantification).

Paper §2 uses γ_cos = 0.397 derived from chr22 penultimate layer (L=30) q70 quantile.
It claims this single value transfers directionally to chr17 (94.6% preserved).

To rigorously quantify the transfer claim on 24 chr, we'd need to recompute q70 on
each chr's penultimate layer. Since we don't have raw D_cos stored genome-wide, we
use a PROXY: infer what γ each chromosome would have chosen if the intron's fraction
below γ (below_frac ≈ 0.05) is the design target.

Simpler alternative (implemented here):
For each chromosome, find the γ' such that intron pseudo-below_frac reaches the same
target as chr22 with γ=0.397. This gives a per-chr γ' vs paper's frozen γ.

Actually since we stored `min_D` per position, we CAN compute per-chr q70 of min_D
across intron positions — this gives a MEANINGFUL calibration proxy (the value at
which 70% of intron positions have min_D ≤ γ).

Reads:
  wgs/results/chr{N}/chr{N}_per_position_chunk*.parquet (min_D column)
  wgs/data/labels/chr{N}_position_labels.npy
Writes:
  wgs/results/genome_summary/wgs_gamma_recalibration.csv
  wgs/results/genome_summary/wgs_gamma_recalibration.json
  wgs/results/genome_summary/wgs_gamma_recalibration.{png,pdf}
"""
from __future__ import annotations
import json, time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TDIG = Path("/NHNHOME/WORKSPACE/0526040123_A/darejinn/tdig")
WGS = TDIG / "wgs/results"
LABEL_DIR = TDIG / "wgs/data/labels"
OUT = WGS / "genome_summary"

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
PAPER_GAMMA = 0.397
TARGET_QUANTILE = 0.70  # q70 = paper's regional calibration


def main():
    t0 = time.time()
    rows = []

    for chrom in CHROMS:
        t_c = time.time()
        lab = np.load(LABEL_DIR / f"{chrom}_position_labels.npy")
        chunks = sorted((WGS / chrom).glob(f"{chrom}_per_position_chunk*.parquet"))

        # Collect intron min_D values (large: use sampling)
        intron_min_D = []
        for cp in chunks:
            d = pd.read_parquet(cp, columns=["pos", "min_D"])
            is_intron = lab[d["pos"].to_numpy()] == 1
            intron_min_D.extend(d.loc[is_intron, "min_D"].to_numpy())
        intron_min_D = np.asarray(intron_min_D, dtype=np.float32)

        # q70 of intron min_D
        gamma_recalib = float(np.quantile(intron_min_D, TARGET_QUANTILE))
        # median
        median_min_D = float(np.median(intron_min_D))
        mean_min_D = float(intron_min_D.mean())
        # frac already below paper γ = 0.397
        frac_below_paper_gamma = float((intron_min_D <= PAPER_GAMMA).mean())
        # what quantile does paper γ correspond to at this chr?
        gamma_quantile_at_paper = float((intron_min_D <= PAPER_GAMMA).mean())

        rows.append({
            "chrom": chrom,
            "n_intron_pos": int(len(intron_min_D)),
            "intron_mean_min_D": mean_min_D,
            "intron_median_min_D": median_min_D,
            "intron_q70_min_D_gamma_recalib": gamma_recalib,
            "delta_vs_paper_gamma_0.397": gamma_recalib - PAPER_GAMMA,
            "frac_intron_below_paper_gamma": frac_below_paper_gamma,
            "target_frac_q70": TARGET_QUANTILE,
        })
        print(f"[{time.strftime('%H:%M:%S')}] {chrom}: n={len(intron_min_D):,} "
              f"q70={gamma_recalib:.4f} (Δ vs paper={gamma_recalib-PAPER_GAMMA:+.4f}) "
              f"paper_γ_at_quantile={gamma_quantile_at_paper:.3f} "
              f"{time.time()-t_c:.1f}s")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "wgs_gamma_recalibration.csv", index=False)

    summary = {
        "paper_gamma_frozen": PAPER_GAMMA,
        "target_quantile": TARGET_QUANTILE,
        "n_chromosomes": len(df),
        "wgs_mean_recalib_gamma": float(df["intron_q70_min_D_gamma_recalib"].mean()),
        "wgs_sd_recalib_gamma": float(df["intron_q70_min_D_gamma_recalib"].std(ddof=1)),
        "wgs_min_recalib_gamma": float(df["intron_q70_min_D_gamma_recalib"].min()),
        "wgs_max_recalib_gamma": float(df["intron_q70_min_D_gamma_recalib"].max()),
        "wgs_mean_shift_vs_paper": float(df["delta_vs_paper_gamma_0.397"].mean()),
        "runtime_sec": time.time() - t0,
    }
    (OUT / "wgs_gamma_recalibration.json").write_text(json.dumps(summary, indent=2, default=str))

    print("\n=== WGS γ recalibration summary ===")
    for k, v in summary.items():
        if isinstance(v, float):
            print(f"  {k}: {v:.4f}")
        else:
            print(f"  {k}: {v}")

    # Figure
    plt.rcParams.update({"font.family": "DejaVu Sans", "savefig.dpi": 300, "savefig.bbox": "tight",
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    ax.bar(range(len(df)), df["intron_q70_min_D_gamma_recalib"], color="#455a64",
           edgecolor="black", lw=0.4)
    ax.axhline(PAPER_GAMMA, ls="--", color="#d62728", lw=1.0,
               label=f"paper γ (chr22) = {PAPER_GAMMA}")
    ax.axhline(float(df["intron_q70_min_D_gamma_recalib"].mean()), ls=":", color="#0d47a1", lw=1.0,
               label=f"WGS mean recalib γ = {df['intron_q70_min_D_gamma_recalib'].mean():.4f}")
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df["chrom"], rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("q70 of intron min_D per chr (proxy for γ recalib)")
    ax.set_title("(a) Per-chromosome γ recalibration proxy",
                 loc="left", fontweight="bold", fontfamily="Times New Roman", fontsize=12, pad=8)
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(axis="y", color="#dddddd", lw=0.5)

    ax = axes[1]
    ax.bar(range(len(df)), df["frac_intron_below_paper_gamma"], color="#c62828",
           edgecolor="black", lw=0.4)
    ax.axhline(TARGET_QUANTILE, ls="--", color="#0d47a1", lw=1.0,
               label=f"target q70 = {TARGET_QUANTILE}")
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df["chrom"], rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("frac intron positions with min_D ≤ paper γ=0.397")
    ax.set_title("(b) Paper γ vs per-chr distribution",
                 loc="left", fontweight="bold", fontfamily="Times New Roman", fontsize=12, pad=8)
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(axis="y", color="#dddddd", lw=0.5)

    fig.suptitle("γ recalibration: paper's frozen γ_cos = 0.397 vs per-chromosome q70",
                 fontsize=13, fontweight="bold", fontfamily="Times New Roman", y=1.02)
    fig.tight_layout()
    fig.savefig(OUT / "wgs_gamma_recalibration.png")
    fig.savefig(OUT / "wgs_gamma_recalibration.pdf")
    plt.close(fig)
    print(f"wrote {OUT/'wgs_gamma_recalibration.png'}")

    print(f"\nTotal runtime: {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
