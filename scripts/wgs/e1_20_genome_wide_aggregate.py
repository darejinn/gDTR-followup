"""EXP1 §Phase B — genome-wide aggregation across all 24 chromosomes.

Combines per-chromosome context summaries into WGS-scale statistics
comparable to paper §3.1 Table.

Inputs: wgs/results/chr{1..22,X,Y}/chr{N}_context_summary.csv
Outputs:
  wgs/results/genome_summary/wgs_context_summary.csv        # per-context WGS-scale d + p
  wgs/results/genome_summary/wgs_context_summary.json       # detailed
  wgs/results/genome_summary/wgs_per_chr_matrix.csv         # 24 chr × 7 context matrix
  wgs/results/genome_summary/wgs_context_summary.{png,pdf}  # figure vs paper §3.1
  wgs/results/genome_summary/wgs_calibration_robustness.csv # per-chr intron c̄ vs frozen γ=0.397
"""
from __future__ import annotations
import json, time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TDIG = Path("/NHNHOME/WORKSPACE/0526040123_A/darejinn/tdig")
WGS = TDIG / "wgs/results"
OUT = WGS / "genome_summary"
OUT.mkdir(parents=True, exist_ok=True)

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
CONTEXTS = ["intron", "splice_donor", "splice_acceptor", "coding_exon",
            "3utr", "5utr", "intergenic"]

# Paper §3.1 pooled chr17+chr22 numbers for comparison
PAPER = {
    "intron":          {"mean_c": 27.72, "d": 0.0},
    "splice_donor":    {"mean_c": 25.55, "d": -0.354},
    "splice_acceptor": {"mean_c": 25.96, "d": -0.340},
    "coding_exon":     {"mean_c": 28.40, "d": +0.08},
    "3utr":            {"mean_c": 27.74, "d": -0.02},
    "5utr":            {"mean_c": 29.22, "d": +0.20},
    "intergenic":      {"mean_c": 28.66, "d": +0.16},
}


def pooled_sd(n1, sd1, n2, sd2):
    return np.sqrt(((n1-1)*sd1**2 + (n2-1)*sd2**2) / (n1+n2-2))


def cohen_d_from_summary(mean1, sd1, n1, mean2, sd2, n2):
    return (mean1 - mean2) / pooled_sd(n1, sd1, n2, sd2)


def main():
    t0 = time.time()

    # ---- 1. Load per-chr summaries ----
    print("=== loading 24 per-chr context summaries ===")
    per_chr = {}
    for chrom in CHROMS:
        csv = WGS / chrom / f"{chrom}_context_summary.csv"
        if not csv.exists():
            print(f"  MISSING: {chrom}")
            continue
        per_chr[chrom] = pd.read_csv(csv)
    print(f"loaded {len(per_chr)} chromosomes")

    # ---- 2. Per-chr × context matrix of d_c_t and d_oscil ----
    print("=== building per-chr × context matrices ===")
    d_ct_rows, d_osc_rows, mean_ct_rows, n_rows = [], [], [], []
    for chrom in CHROMS:
        if chrom not in per_chr: continue
        df = per_chr[chrom]
        row_d_ct = {"chrom": chrom}
        row_d_osc = {"chrom": chrom}
        row_mean = {"chrom": chrom}
        row_n = {"chrom": chrom}
        for ctx in CONTEXTS:
            r = df[df["context"] == ctx]
            if len(r) == 0:
                row_d_ct[ctx] = np.nan; row_d_osc[ctx] = np.nan
                row_mean[ctx] = np.nan; row_n[ctx] = 0
                continue
            row_d_ct[ctx] = float(r["cohen_d_c_t_vs_intron"].iloc[0])
            row_d_osc[ctx] = float(r["cohen_d_oscil_vs_intron"].iloc[0])
            row_mean[ctx] = float(r["mean_c_t"].iloc[0])
            row_n[ctx] = int(r["n"].iloc[0])
        d_ct_rows.append(row_d_ct); d_osc_rows.append(row_d_osc)
        mean_ct_rows.append(row_mean); n_rows.append(row_n)

    matrix_d_ct = pd.DataFrame(d_ct_rows)
    matrix_d_osc = pd.DataFrame(d_osc_rows)
    matrix_mean = pd.DataFrame(mean_ct_rows)
    matrix_n = pd.DataFrame(n_rows)
    matrix_d_ct.to_csv(OUT / "wgs_per_chr_d_c_t.csv", index=False)
    matrix_d_osc.to_csv(OUT / "wgs_per_chr_d_oscil.csv", index=False)
    matrix_mean.to_csv(OUT / "wgs_per_chr_mean_c_t.csv", index=False)
    matrix_n.to_csv(OUT / "wgs_per_chr_n.csv", index=False)

    # ---- 3. WGS pooled per-context: sum n, recompute weighted mean & d ----
    # For pooling means from summary: n-weighted mean is exact.
    # For pooling d: we need to reconstruct via pooled variance.
    # Approach: total_mean = sum(n_i * mean_i) / sum(n_i)
    # We don't have per-chr SD for c_t directly but can compute from mean_below_frac etc.
    # Simpler: use per-chr n and mean, compute weighted mean; for SD, aggregate as if a
    # single population by pooling (n_i-1)*var_i approx. Get var from per-position parquets
    # would be more accurate but adds cost. For now use pooled-mean + per-chr d avg.
    print("=== WGS pooled per-context ===")
    wgs_summary = []
    for ctx in CONTEXTS:
        n_sum = int(matrix_n[ctx].sum())
        weighted_mean = float((matrix_mean[ctx] * matrix_n[ctx]).sum() / n_sum) if n_sum > 0 else np.nan
        # Per-chr d weighted average by chr n
        d_ct_wavg = float((matrix_d_ct[ctx] * matrix_n[ctx]).sum() / n_sum) if n_sum > 0 else np.nan
        d_osc_wavg = float((matrix_d_osc[ctx] * matrix_n[ctx]).sum() / n_sum) if n_sum > 0 else np.nan
        # Per-chr d min/max/std (variability across chromosomes)
        d_ct_std = float(np.nanstd(matrix_d_ct[ctx], ddof=1))
        d_osc_std = float(np.nanstd(matrix_d_osc[ctx], ddof=1))
        paper_d = PAPER[ctx]["d"] if ctx in PAPER else np.nan
        paper_c = PAPER[ctx]["mean_c"] if ctx in PAPER else np.nan
        wgs_summary.append({
            "context": ctx,
            "wgs_n_positions": n_sum,
            "wgs_mean_c_t": weighted_mean,
            "wgs_d_c_t_wavg": d_ct_wavg,
            "wgs_d_c_t_std_across_chr": d_ct_std,
            "wgs_d_oscil_wavg": d_osc_wavg,
            "wgs_d_oscil_std_across_chr": d_osc_std,
            "paper_pooled_c_t": paper_c,
            "paper_d_c_t": paper_d,
            "delta_wgs_vs_paper_mean_c_t": weighted_mean - paper_c if paper_c else np.nan,
            "delta_wgs_vs_paper_d_c_t": d_ct_wavg - paper_d if paper_d is not None else np.nan,
        })
    wgs_df = pd.DataFrame(wgs_summary)
    wgs_df.to_csv(OUT / "wgs_context_summary.csv", index=False)
    (OUT / "wgs_context_summary.json").write_text(json.dumps(wgs_summary, indent=2, default=str))

    print("\n=== WGS Genome-wide per-context summary ===")
    print(wgs_df.to_string(index=False))

    # ---- 4. Calibration robustness: per-chr intron c̄ vs frozen γ=0.397 ----
    intron_c_ts = matrix_mean["intron"].dropna().to_numpy()
    calib = pd.DataFrame({
        "chrom": matrix_mean["chrom"],
        "intron_mean_c_t": matrix_mean["intron"],
        "delta_vs_paper_intron_27.72": matrix_mean["intron"] - 27.72,
        "delta_vs_wgs_mean": matrix_mean["intron"] - float(intron_c_ts.mean()),
    })
    calib.to_csv(OUT / "wgs_calibration_robustness.csv", index=False)
    print("\n=== Calibration robustness (intron c̄ across chromosomes) ===")
    print(f"WGS-mean intron c̄:  {intron_c_ts.mean():.4f}")
    print(f"WGS-SD intron c̄:    {intron_c_ts.std(ddof=1):.4f}")
    print(f"Paper intron c̄:     27.72 (chr17+chr22 pool)")
    print(f"Shift vs paper:      {intron_c_ts.mean() - 27.72:+.4f} layers")
    print(f"Chr-to-chr range:    [{intron_c_ts.min():.4f}, {intron_c_ts.max():.4f}]")

    # ---- 5. Figures ----
    print("=== rendering figures ===")
    plt.rcParams.update({"font.family": "DejaVu Sans", "savefig.dpi": 300, "savefig.bbox": "tight",
                         "axes.spines.top": False, "axes.spines.right": False})

    fig, axes = plt.subplots(2, 2, figsize=(15, 10))

    # (a) WGS mean c_t vs paper per context
    ax = axes[0, 0]
    ctxs = wgs_df["context"].tolist()
    xw = np.arange(len(ctxs)); w = 0.36
    ax.bar(xw - w/2, wgs_df["wgs_mean_c_t"], w, color="#0d47a1",
           edgecolor="black", lw=0.4, label="WGS (24 chrs)")
    ax.bar(xw + w/2, wgs_df["paper_pooled_c_t"], w, color="#d62728", alpha=0.55,
           edgecolor="black", lw=0.4, label="paper §3.1 (chr17+22 pool)")
    for xi, (m, n) in enumerate(zip(wgs_df["wgs_mean_c_t"], wgs_df["wgs_n_positions"])):
        ax.text(xi - w/2, m + 0.1, f"{m:.2f}\n(n={n/1e6:.0f}M)", ha="center",
                fontsize=7, fontweight="bold")
    ax.set_xticks(xw)
    ax.set_xticklabels([c.replace("_", "\n") for c in ctxs])
    ax.axhline(27.72, ls="--", color="#555", lw=0.7, label="paper intron 27.72")
    ax.set_ylabel(r"mean settling depth $\bar c(t)$")
    ax.set_title("(a) WGS vs paper §3.1 — per-context mean $c(t)$",
                 loc="left", fontweight="bold", fontfamily="Times New Roman", fontsize=12, pad=8)
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(axis="y", color="#dddddd", lw=0.5)

    # (b) Cohen's d comparison
    ax = axes[0, 1]
    ax.barh(xw - w/2, wgs_df["wgs_d_c_t_wavg"], w,
            xerr=wgs_df["wgs_d_c_t_std_across_chr"],
            color="#d62728", edgecolor="black", lw=0.4, label="WGS d (c_t)")
    ax.barh(xw + w/2, wgs_df["wgs_d_oscil_wavg"], w,
            xerr=wgs_df["wgs_d_oscil_std_across_chr"],
            color="#1f77b4", edgecolor="black", lw=0.4, label="WGS d (oscil)")
    ax.scatter(wgs_df["paper_d_c_t"], xw - w/2, s=60, color="black",
               marker="D", label="paper §3.1 d (c_t)", zorder=5)
    ax.axvline(0, color="#555", lw=0.7)
    ax.set_yticks(xw)
    ax.set_yticklabels(ctxs)
    ax.invert_yaxis()
    ax.set_xlabel("Cohen's d vs intron")
    ax.set_title("(b) WGS effect sizes: c(t) + oscil, error bars = per-chr SD",
                 loc="left", fontweight="bold", fontfamily="Times New Roman", fontsize=12, pad=8)
    ax.legend(loc="lower right", fontsize=8)
    ax.grid(axis="x", color="#dddddd", lw=0.5)

    # (c) Per-chromosome intron c_t (calibration robustness)
    ax = axes[1, 0]
    valid = matrix_mean.dropna(subset=["intron"])
    ax.bar(range(len(valid)), valid["intron"], color="#455a64",
           edgecolor="black", lw=0.4)
    ax.axhline(27.72, ls="--", color="#d62728", lw=1.0, label="paper (chr17+22 pool)")
    ax.axhline(float(intron_c_ts.mean()), ls=":", color="#0d47a1", lw=1.0,
               label=f"WGS mean = {intron_c_ts.mean():.2f}")
    ax.set_xticks(range(len(valid)))
    ax.set_xticklabels(valid["chrom"], rotation=45, ha="right", fontsize=7)
    ax.set_ylabel(r"intron $\bar c(t)$")
    ax.set_title("(c) Calibration robustness: intron $\\bar c(t)$ per chromosome (γ_cos=0.397 frozen)",
                 loc="left", fontweight="bold", fontfamily="Times New Roman", fontsize=12, pad=8)
    ax.legend(loc="lower right", fontsize=8)
    ax.grid(axis="y", color="#dddddd", lw=0.5)

    # (d) Splice donor d per chromosome (H3a orthogonal-axes replication)
    ax = axes[1, 1]
    donor_ct = matrix_d_ct["splice_donor"].dropna()
    donor_osc = matrix_d_osc["splice_donor"].dropna()
    accept_ct = matrix_d_ct["splice_acceptor"].dropna()
    accept_osc = matrix_d_osc["splice_acceptor"].dropna()
    ax.scatter(donor_ct, donor_osc, s=90, color="#1f77b4",
               edgecolor="black", label="splice donor (per chr)", alpha=0.85)
    ax.scatter(accept_ct, accept_osc, s=90, color="#2ca02c",
               edgecolor="black", marker="s", label="splice acceptor (per chr)", alpha=0.85)
    ax.scatter([-0.354], [None], s=120, color="#1f77b4", marker="*",
               edgecolor="black", label="paper donor d_c_t = -0.354")
    ax.axhline(0, color="#555", lw=0.6)
    ax.axvline(0, color="#555", lw=0.6)
    ax.set_xlabel("Cohen's d — c(t) vs intron")
    ax.set_ylabel("Cohen's d — oscil vs intron")
    ax.set_title("(d) Per-chr splice sites in c(t) × oscil plane\n(H3a orthogonal-axes replication)",
                 loc="left", fontweight="bold", fontfamily="Times New Roman", fontsize=12, pad=8)
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(color="#dddddd", lw=0.5)

    fig.suptitle(f"EXP1 §Phase B — WGS genome-wide aggregation (24 chromosomes, ~2.7 billion positions)",
                 fontsize=14, fontweight="bold", fontfamily="Times New Roman", y=1.01)
    fig.tight_layout()
    fig.savefig(OUT / "wgs_context_summary.png")
    fig.savefig(OUT / "wgs_context_summary.pdf")
    plt.close(fig)
    print(f"wrote {OUT/'wgs_context_summary.png'}")

    # ---- 6. Text report ----
    report = f"""# WGS Genome-wide aggregation report — 24 chromosomes

**Date**: {time.strftime('%Y-%m-%d %H:%M:%S')}
**Total positions analyzed**: {matrix_n[CONTEXTS].sum().sum():,}
**Chromosomes**: {len(per_chr)}

## Per-context WGS-scale summary vs paper §3.1

| Context | n (WGS) | WGS mean c_t | Paper c_t | WGS d_c_t (mean±SD) | Paper d_c_t | WGS d_oscil |
|---|---|---|---|---|---|---|
"""
    for _, r in wgs_df.iterrows():
        report += (f"| {r['context']} | {int(r['wgs_n_positions']):,} | "
                   f"{r['wgs_mean_c_t']:.3f} | {r['paper_pooled_c_t']:.2f} | "
                   f"{r['wgs_d_c_t_wavg']:+.4f} ± {r['wgs_d_c_t_std_across_chr']:.4f} | "
                   f"{r['paper_d_c_t']:+.3f} | "
                   f"{r['wgs_d_oscil_wavg']:+.4f} ± {r['wgs_d_oscil_std_across_chr']:.4f} |\n")

    report += f"""

## Calibration robustness (γ_cos = 0.397 frozen from paper chr22)
- WGS-mean intron c̄: {intron_c_ts.mean():.4f}
- WGS-SD intron c̄:   {intron_c_ts.std(ddof=1):.4f}
- Paper intron c̄:    27.72 (chr17+chr22 pool)
- Shift vs paper:     {intron_c_ts.mean() - 27.72:+.4f} layers
- Chr-to-chr range:   [{intron_c_ts.min():.4f}, {intron_c_ts.max():.4f}]

## Directional replication of paper §3.1
"""
    for _, r in wgs_df.iterrows():
        wgs_sign = np.sign(r["wgs_d_c_t_wavg"])
        paper_sign = np.sign(r["paper_d_c_t"])
        match = "✓ same sign" if wgs_sign == paper_sign or paper_sign == 0 else "✗ opposite"
        report += f"- **{r['context']}**: WGS d_c_t = {r['wgs_d_c_t_wavg']:+.4f}, paper = {r['paper_d_c_t']:+.3f} → {match}\n"

    (OUT / "wgs_context_report.md").write_text(report)
    print(f"wrote {OUT/'wgs_context_report.md'}")
    print(f"\nTotal runtime: {time.time()-t0:.1f} s")


if __name__ == "__main__":
    main()
