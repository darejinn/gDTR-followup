"""EXP3 §3.7 — Visualize H3a: c_t vs oscil per context (orthogonal axes)."""
from __future__ import annotations
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

TDIG = Path("/NHNHOME/WORKSPACE/0526040123_A/darejinn/tdig")
RES = TDIG / "exp3_threshold_crossing/results"
FIG = TDIG / "exp3_threshold_crossing/figures"

LABEL_MAP = {0: "intergenic", 1: "intron", 2: "coding_exon", 3: "5utr", 4: "3utr",
             5: "splice_donor", 6: "splice_acceptor"}
COLORS = {"intergenic": "#8e8e93", "intron": "#bdbdbd", "coding_exon": "#78909c",
          "3utr": "#a1887f", "5utr": "#c5cae9", "splice_donor": "#1f77b4",
          "splice_acceptor": "#2ca02c"}


def setup():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.titlesize": 10, "axes.labelsize": 9,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
        "axes.spines.top": False, "axes.spines.right": False,
        "savefig.dpi": 300, "savefig.bbox": "tight",
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def main():
    setup()
    df = pd.read_parquet(RES / "A_crossing_stats.parquet")
    df["context"] = df["label"].map(LABEL_MAP).fillna("unknown")
    h3a = json.loads((RES / "H3a_context_test.json").read_text())

    contexts = ["intergenic", "intron", "coding_exon", "splice_donor", "splice_acceptor"]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))

    # ==== Panel (a): scatter of mean c_t vs mean oscil per context (2D axis test) ====
    ax = axes[0]
    for ctx in contexts:
        sub = df[df["context"] == ctx]
        if len(sub) < 30: continue
        ax.scatter(sub["c_t"].mean(), sub["oscil"].mean(),
                   s=300, color=COLORS[ctx], edgecolor="black", lw=0.8, zorder=3, label=f"{ctx}\n(n={len(sub):,})")
        ax.annotate(ctx.replace("_", " "), (sub["c_t"].mean(), sub["oscil"].mean()),
                    xytext=(6, 6), textcoords="offset points", fontsize=8.5, fontweight="bold")
    ax.set_xlabel(r"mean settling depth $\bar c(t)$ — paper §3.1 axis")
    ax.set_ylabel(r"mean oscillation count — EXP3 new axis")
    ax.set_title("(a) Two orthogonal axes:\n$c(t)$ (paper) vs oscil (EXP3)",
                 loc="left", fontsize=12, fontweight="bold", fontfamily="Times New Roman", pad=8)
    ax.grid(color="#dddddd", lw=0.5)
    ax.axhline(df[df["context"]=="intron"]["oscil"].mean(), color="#d62728", ls="--", lw=0.7, alpha=0.5)
    ax.axvline(df[df["context"]=="intron"]["c_t"].mean(), color="#d62728", ls="--", lw=0.7, alpha=0.5, label="intron baseline")
    ax.legend(loc="upper left", fontsize=7.5)

    # ==== Panel (b): oscil distribution per context (histogram, log y) ====
    ax = axes[1]
    x = np.arange(6)  # oscil values 0..5
    width = 0.14
    for i, ctx in enumerate(contexts):
        sub = df[df["context"] == ctx]
        if len(sub) < 30: continue
        counts = np.zeros(6)
        for v, c in sub["oscil"].value_counts().items():
            if 0 <= int(v) < 6: counts[int(v)] = c
        counts = counts / counts.sum()  # normalize per context
        ax.bar(x + (i - 2) * width, counts, width, color=COLORS[ctx], edgecolor="black", lw=0.4,
               label=ctx.replace("_", " "))
    ax.set_xlabel("oscillation count (# extra crossings beyond a single dip)")
    ax.set_ylabel("fraction of positions")
    ax.set_yscale("log")
    ax.set_ylim(0.0005, 1.5)
    ax.set_xticks(x)
    ax.set_title("(b) Oscillation distribution by context",
                 loc="left", fontsize=12, fontweight="bold", fontfamily="Times New Roman", pad=8)
    ax.legend(loc="upper right", ncol=2, fontsize=7)
    ax.grid(axis="y", color="#dddddd", lw=0.5, which="both")

    # ==== Panel (c): Cohen's d for oscil and c_t per context ====
    ax = axes[2]
    contexts_test = [c for c in ["splice_donor", "splice_acceptor", "coding_exon", "intergenic"]
                     if c in h3a["contexts"]]
    d_osc = [h3a["contexts"][c]["oscil"]["cohen_d_vs_ref"] for c in contexts_test]
    d_ct = [h3a["contexts"][c]["c_t"]["cohen_d_vs_ref"] for c in contexts_test]
    y = np.arange(len(contexts_test))
    ax.barh(y - 0.2, d_ct, height=0.35, color="#d62728", edgecolor="black", lw=0.4, label="$c(t)$ (paper axis)")
    ax.barh(y + 0.2, d_osc, height=0.35, color="#1f77b4", edgecolor="black", lw=0.4, label="oscil (EXP3 axis)")
    for yi, (dc, do) in enumerate(zip(d_ct, d_osc)):
        ax.text(dc + (0.02 if dc >= 0 else -0.02), yi - 0.2, f"{dc:+.2f}",
                va="center", ha="left" if dc >= 0 else "right", fontsize=8, color="#d62728", fontweight="bold")
        ax.text(do + (0.02 if do >= 0 else -0.02), yi + 0.2, f"{do:+.2f}",
                va="center", ha="left" if do >= 0 else "right", fontsize=8, color="#1f77b4", fontweight="bold")
    ax.axvline(0, color="#555", lw=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels([c.replace("_", " ") for c in contexts_test])
    ax.set_xlim(-0.8, 0.8)
    ax.set_xlabel("Cohen's d vs intron baseline")
    ax.invert_yaxis()
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title(r"(c) Effect sizes: $c(t)$ vs oscil disagree in sign",
                 loc="left", fontsize=12, fontweight="bold", fontfamily="Times New Roman", pad=8)
    ax.grid(axis="x", color="#dddddd", lw=0.5)

    fig.suptitle("EXP3 — Threshold-crossing dynamics reveal orthogonal readout to $c(t)$",
                 fontsize=12, fontweight="bold", fontfamily="Times New Roman", y=1.02)

    FIG.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / "EXP3_H3a_orthogonal_axes.png")
    fig.savefig(FIG / "EXP3_H3a_orthogonal_axes.pdf")
    plt.close(fig)
    print(f"wrote {FIG / 'EXP3_H3a_orthogonal_axes.png'}")


if __name__ == "__main__":
    main()
