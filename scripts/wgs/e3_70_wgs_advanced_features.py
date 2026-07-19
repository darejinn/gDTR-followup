"""EXP3 §Phase B — Advanced crossings applied at WGS scale (24 chr).

The WGS batch stored per-position: c_t, oscil, n_enter, n_exit, below_frac, min_D, argmin_layer.
It did NOT retain raw D_cos trajectory (would be ~230 GB genome-wide).

So at WGS scale we compute derived features that are available from the stored columns:
  - crosses_at_all:      first_below >= 0  (does D_cos ever fall below γ)
  - committed:           n_enter=1 AND n_exit=0 (clean single crossing that stays)
  - dipped:              n_enter=1 AND n_exit=1 (crosses in and out)
  - deliberating:        oscil >= 1
  - amplitude:           γ - min_D  (how far below γ, negative if never crosses)
  - late_argmin:         argmin_layer >= 22 (deep-layer minimum)
  - early_argmin:        argmin_layer < 10  (early-layer minimum)
  - broad_commitment:    below_frac >= 0.1  (10% or more layers below γ)

Per context, per chromosome, we compute enrichment stats vs intron baseline.

Aggregate across 24 chr → WGS-scale H3d validation.

Reads:  wgs/results/chr{N}/chr{N}_per_position_chunk*.parquet
        wgs/data/labels/chr{N}_position_labels.npy
Writes: wgs/results/genome_summary/wgs_h3d_context_features.csv
        wgs/results/genome_summary/wgs_h3d_summary.json
        wgs/results/genome_summary/wgs_h3d_features.{png,pdf}
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
LABEL_DIR = TDIG / "wgs/data/labels"
OUT = WGS / "genome_summary"

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
GAMMA = 0.397
LABEL_MAP = {0: "intergenic", 1: "intron", 2: "coding_exon", 3: "5utr", 4: "3utr",
             5: "splice_donor", 6: "splice_acceptor"}
CTX_ORDER = ["intron", "splice_donor", "splice_acceptor", "coding_exon", "3utr", "5utr", "intergenic"]


def load_and_derive_chr(chrom: str) -> pd.DataFrame:
    """Load per-chromosome chunks, join labels, compute derived features.
    Returns per-context aggregate DataFrame."""
    chunks = sorted((WGS / chrom).glob(f"{chrom}_per_position_chunk*.parquet"))
    lab = np.load(LABEL_DIR / f"{chrom}_position_labels.npy")

    total_rows = 0
    # Aggregators per context
    agg = {ctx: {"n": 0, "mean_amp": 0.0, "n_committed": 0, "n_dipped": 0,
                 "n_deliberating": 0, "n_crosses": 0, "n_late_argmin": 0,
                 "n_early_argmin": 0, "n_broad_commit": 0,
                 "mean_argmin": 0.0}
           for ctx in CTX_ORDER}

    for cp in chunks:
        d = pd.read_parquet(cp, columns=["pos", "c_t", "oscil", "n_enter", "n_exit",
                                          "below_frac", "min_D", "argmin_layer"])
        d["label"] = lab[d["pos"].to_numpy()]
        d["context"] = d["label"].map(LABEL_MAP).fillna("unknown")
        # Derived features (per row)
        d["amplitude"] = np.clip(GAMMA - d["min_D"], 0, None).astype(np.float32)
        d["committed"] = ((d["n_enter"] == 1) & (d["n_exit"] == 0)).astype(np.int8)
        d["dipped"] = ((d["n_enter"] == 1) & (d["n_exit"] == 1)).astype(np.int8)
        d["deliberating"] = (d["oscil"] >= 1).astype(np.int8)
        d["crosses"] = (d["min_D"] <= GAMMA).astype(np.int8)
        d["late_argmin"] = (d["argmin_layer"] >= 22).astype(np.int8)
        d["early_argmin"] = (d["argmin_layer"] < 10).astype(np.int8)
        d["broad_commit"] = (d["below_frac"] >= 0.1).astype(np.int8)

        for ctx in CTX_ORDER:
            sub = d[d["context"] == ctx]
            if len(sub) == 0: continue
            a = agg[ctx]
            n_before = a["n"]
            n_add = len(sub)
            a["n"] += n_add
            # running mean update
            a["mean_amp"] = (a["mean_amp"] * n_before + sub["amplitude"].sum()) / a["n"]
            a["mean_argmin"] = (a["mean_argmin"] * n_before + sub["argmin_layer"].sum()) / a["n"]
            # counts
            a["n_committed"] += int(sub["committed"].sum())
            a["n_dipped"] += int(sub["dipped"].sum())
            a["n_deliberating"] += int(sub["deliberating"].sum())
            a["n_crosses"] += int(sub["crosses"].sum())
            a["n_late_argmin"] += int(sub["late_argmin"].sum())
            a["n_early_argmin"] += int(sub["early_argmin"].sum())
            a["n_broad_commit"] += int(sub["broad_commit"].sum())

        total_rows += len(d)

    # Convert to per-context frac (rate) rows
    rows = []
    for ctx in CTX_ORDER:
        a = agg[ctx]
        if a["n"] == 0: continue
        rows.append({
            "chrom": chrom,
            "context": ctx,
            "n": a["n"],
            "mean_amplitude": a["mean_amp"],
            "mean_argmin_layer": a["mean_argmin"],
            "frac_committed": a["n_committed"] / a["n"],
            "frac_dipped": a["n_dipped"] / a["n"],
            "frac_deliberating": a["n_deliberating"] / a["n"],
            "frac_crosses": a["n_crosses"] / a["n"],
            "frac_late_argmin": a["n_late_argmin"] / a["n"],
            "frac_early_argmin": a["n_early_argmin"] / a["n"],
            "frac_broad_commit": a["n_broad_commit"] / a["n"],
        })
    return pd.DataFrame(rows), total_rows


def main():
    t0 = time.time()
    all_rows = []
    grand_total = 0
    for i, chrom in enumerate(CHROMS):
        t_c = time.time()
        df, n = load_and_derive_chr(chrom)
        grand_total += n
        all_rows.append(df)
        print(f"[{time.strftime('%H:%M:%S')}] {i+1}/{len(CHROMS)} {chrom}: {n:,} pos, {time.time()-t_c:.1f}s")

    per_chr_ctx = pd.concat(all_rows, ignore_index=True)
    per_chr_ctx.to_csv(OUT / "wgs_h3d_per_chr_per_context.csv", index=False)
    print(f"\nGrand total positions analyzed: {grand_total:,}")

    # WGS pooled per context — n-weighted average
    wgs_rows = []
    features = ["mean_amplitude", "mean_argmin_layer", "frac_committed", "frac_dipped",
                "frac_deliberating", "frac_crosses", "frac_late_argmin",
                "frac_early_argmin", "frac_broad_commit"]
    for ctx in CTX_ORDER:
        sub = per_chr_ctx[per_chr_ctx["context"] == ctx]
        if len(sub) == 0: continue
        n_sum = int(sub["n"].sum())
        row = {"context": ctx, "n": n_sum}
        for f in features:
            row[f] = float((sub[f] * sub["n"]).sum() / n_sum)
            row[f"{f}_std_across_chr"] = float(sub[f].std(ddof=1))
        wgs_rows.append(row)
    wgs = pd.DataFrame(wgs_rows)
    wgs.to_csv(OUT / "wgs_h3d_context_features.csv", index=False)
    print(f"\n=== WGS-scale H3d per-context features ===")
    print(wgs.round(4).to_string(index=False))

    # Save JSON
    (OUT / "wgs_h3d_summary.json").write_text(
        json.dumps({"grand_total_positions": grand_total,
                    "wgs_per_context": wgs.to_dict(orient="records"),
                    "runtime_sec": time.time() - t0}, indent=2, default=str))

    # Figure — bar chart of enrichment factor for each feature at splice sites vs intron
    plt.rcParams.update({"font.family": "DejaVu Sans", "savefig.dpi": 300, "savefig.bbox": "tight",
                         "axes.spines.top": False, "axes.spines.right": False})
    intron_row = wgs[wgs["context"] == "intron"].iloc[0]
    feats_plot = ["mean_amplitude", "mean_argmin_layer", "frac_committed",
                  "frac_deliberating", "frac_crosses", "frac_late_argmin",
                  "frac_broad_commit"]
    contexts_plot = ["splice_donor", "splice_acceptor", "coding_exon", "3utr", "5utr", "intergenic"]
    colors = {"splice_donor": "#1f77b4", "splice_acceptor": "#2ca02c", "coding_exon": "#78909c",
              "3utr": "#a1887f", "5utr": "#c5cae9", "intergenic": "#8e8e93"}

    fig, ax = plt.subplots(figsize=(13, 5.5))
    x = np.arange(len(feats_plot))
    w = 0.13
    for i, ctx in enumerate(contexts_plot):
        row = wgs[wgs["context"] == ctx]
        if len(row) == 0: continue
        row = row.iloc[0]
        # enrichment = ctx / intron (ratio)
        ratio = np.array([row[f] / max(intron_row[f], 1e-9) for f in feats_plot])
        ax.bar(x + (i - 2.5) * w, ratio, w, color=colors[ctx], edgecolor="black",
               lw=0.4, label=ctx.replace("_", " "))
    ax.axhline(1.0, color="#555", lw=0.6, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels([f.replace("mean_", "").replace("frac_", "").replace("_", "\n") for f in feats_plot], fontsize=9)
    ax.set_ylabel("ratio vs intron baseline (log scale)")
    ax.set_yscale("log")
    ax.set_title(f"WGS H3d — advanced crossing features per context (n_total = {grand_total/1e9:.2f} B positions, 24 chromosomes)",
                 fontweight="bold", fontfamily="Times New Roman", fontsize=12)
    ax.legend(loc="upper right", ncol=2, fontsize=8)
    ax.grid(axis="y", color="#dddddd", lw=0.5, which="both")

    fig.tight_layout()
    fig.savefig(OUT / "wgs_h3d_features.png")
    fig.savefig(OUT / "wgs_h3d_features.pdf")
    plt.close(fig)
    print(f"wrote {OUT/'wgs_h3d_features.png'}")

    print(f"\nTotal runtime: {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
