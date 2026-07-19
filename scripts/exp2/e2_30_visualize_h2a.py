"""EXP2 §3 — visualize H2a: layer feature comparison.

Reads results/H2a_argmax_layer_auroc.csv + H2a_regression_summary.json.
Produces: figures/EXP2_H2a_layer_auroc.png (+ pdf).
"""
from __future__ import annotations
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

TDIG = Path("/NHNHOME/WORKSPACE/0526040123_A/darejinn/tdig")
RES = TDIG / "exp2_variant_downstream/results"
FIG = TDIG / "exp2_variant_downstream/figures"


def setup():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.titlesize": 10, "axes.labelsize": 9,
        "xtick.labelsize": 8, "ytick.labelsize": 8,
        "legend.fontsize": 8, "axes.spines.top": False, "axes.spines.right": False,
        "savefig.dpi": 300, "savefig.bbox": "tight",
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def main():
    setup()
    aur = pd.read_csv(RES / "H2a_argmax_layer_auroc.csv")
    summ = json.loads((RES / "H2a_regression_summary.json").read_text())
    per_layer = summ["per_layer_auroc"]
    layers = sorted(int(k) for k in per_layer)
    aurocs = [per_layer[str(l)] for l in layers]
    best_l = summ["best_fixed_layer"]

    layer_dist = pd.read_csv(RES / "H2a_argmax_layer_distribution.csv")

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.2), gridspec_kw={"width_ratios": [1.4, 1.4, 1.0]})

    # Panel (a): per-layer AUROC + variant-adaptive comparison
    ax = axes[0]
    ax.plot(layers, aurocs, "o-", color="#1f77b4", lw=1.4, ms=4, label="Fixed single-layer $\\Delta D_{cos}(\\ell)$")
    ax.axhline(summ["features_auroc"]["32d_cos_vector (paper baseline)"]["mean"],
               color="#2ca02c", ls="--", lw=1.2, label=f"32-d vector (paper: 0.844)  = {summ['features_auroc']['32d_cos_vector (paper baseline)']['mean']:.4f}")
    ax.axhline(summ["features_auroc"]["max_abs_dD_precomputed"]["mean"],
               color="#d62728", ls="--", lw=1.2, label=f"$\\max_\\ell|\\Delta D_{{cos}}(\\ell)|$ (var-adaptive) = {summ['features_auroc']['max_abs_dD_precomputed']['mean']:.4f}")
    ax.axvline(best_l, color="#1f77b4", ls=":", lw=0.8)
    ax.scatter([best_l], [aurocs[layers.index(best_l)]], s=110, color="#1f77b4",
               edgecolor="black", zorder=5)
    ax.set_xlabel("layer $\\ell$"); ax.set_ylabel("AUROC (10-fold CV, seed 42)")
    ax.set_title("(a) Fixed-layer vs variant-adaptive $|\\Delta D_{cos}|$", loc="left",
                 fontsize=13, fontweight="bold", fontfamily="Times New Roman", pad=8)
    ax.set_ylim(0.45, 0.90)
    ax.legend(loc="lower right", framealpha=0.9)
    ax.grid(axis="y", color="#dddddd", lw=0.5)

    # Panel (b): feature bar comparison
    ax = axes[1]
    order = ["argmax_layer_value (H2a)",
             "signed_argmax_precomputed",
             "abs(argmax_layer_value)",
             f"best_fixed_layer L={best_l}",
             "max_abs_dD_precomputed",
             "argmax_layer + Evo2_LL ensemble",
             "32d_cos_vector (paper baseline)"]
    means = [summ["features_auroc"][k]["mean"] for k in order]
    stds = [summ["features_auroc"][k]["std"] for k in order]
    short_names = {
        "argmax_layer_value (H2a)": "signed $\\Delta D_{cos}$ at\nargmax layer",
        "signed_argmax_precomputed": "signed argmax\n(precomputed)",
        "abs(argmax_layer_value)": "$|\\Delta D_{cos}|$ at\nargmax layer",
        f"best_fixed_layer L={best_l}": f"fixed $\\ell$=L{best_l}\n(paper baseline)",
        "max_abs_dD_precomputed": "$\\max_\\ell |\\Delta D_{cos}(\\ell)|$\n(H2a positive)",
        "argmax_layer + Evo2_LL ensemble": "argmax_val +\nEvo2 LL",
        "32d_cos_vector (paper baseline)": "32-d vector\n(paper: 0.844)",
    }
    colors = ["#c62828", "#c62828", "#ef6c00", "#1976d2", "#2e7d32", "#6a1b9a", "#0d47a1"]
    y_pos = np.arange(len(order))
    bars = ax.barh(y_pos, means, xerr=stds, color=colors, edgecolor="black", lw=0.5, height=0.65)
    for yi, (m, s) in enumerate(zip(means, stds)):
        ax.text(m + 0.005, yi, f"{m:.3f}", va="center", fontsize=8, fontweight="bold")
    ax.set_yticks(y_pos)
    ax.set_yticklabels([short_names[o] for o in order])
    ax.set_xlim(0.45, 0.90)
    ax.invert_yaxis()
    ax.axvline(0.5, color="#888", lw=0.6)
    ax.set_xlabel("AUROC ± SD (10-fold CV)")
    ax.set_title("(b) Feature comparison", loc="left", fontsize=13, fontweight="bold",
                 fontfamily="Times New Roman", pad=8)
    ax.grid(axis="x", color="#dddddd", lw=0.5)

    # Panel (c): argmax layer distribution
    ax = axes[2]
    ax.bar(layer_dist["layer"], layer_dist["n_variants_argmax_here"], color="#455a64",
           edgecolor="black", lw=0.4)
    ax.set_xlabel("argmax layer")
    ax.set_ylabel("# variants")
    ax.set_title(f"(c) Variant argmax-layer distribution\n({summ['unique_argmax_layers']}/32 layers used)",
                 loc="left", fontsize=13, fontweight="bold", fontfamily="Times New Roman", pad=8)
    ax.grid(axis="y", color="#dddddd", lw=0.5)

    FIG.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / "EXP2_H2a_layer_auroc.png")
    fig.savefig(FIG / "EXP2_H2a_layer_auroc.pdf")
    plt.close(fig)
    print(f"wrote {FIG / 'EXP2_H2a_layer_auroc.png'}")
    print(f"wrote {FIG / 'EXP2_H2a_layer_auroc.pdf'}")


if __name__ == "__main__":
    main()
