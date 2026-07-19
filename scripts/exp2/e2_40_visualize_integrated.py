"""Visualize integrated H2b: feature-set contribution to 4-class subtype macro-F1."""
from __future__ import annotations
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

TDIG = Path("/NHNHOME/WORKSPACE/0526040123_A/darejinn/tdig")
JSON = TDIG / "exp2_variant_downstream/results/H2b_integrated_features.json"
FIG = TDIG / "exp2_variant_downstream/figures"


def main():
    plt.rcParams.update({"font.family": "DejaVu Sans", "savefig.dpi": 300, "savefig.bbox": "tight",
                         "axes.spines.top": False, "axes.spines.right": False})
    d = json.loads(JSON.read_text())
    r = d["results"]

    order = [
        "EXP3_features only (H3b-derived)",
        "scalars only",
        "cos32_only (paper baseline)",
        "cos32 + EXP3_features",
        "cos32 + scalars",
        "cos32 + EXP3 + scalars (all)",
    ]
    labels_short = {
        "EXP3_features only (H3b-derived)": "EXP3 only\n(Δoscil, Δn_enter, ...)",
        "scalars only": "scalars only\n(max_abs, argmax_layer, ΔLL, ...)",
        "cos32_only (paper baseline)": "32-d ΔD_cos\n(paper baseline)",
        "cos32 + EXP3_features": "cos32 + EXP3",
        "cos32 + scalars": "cos32 + scalars",
        "cos32 + EXP3 + scalars (all)": "cos32 + EXP3 + scalars\n(all combined)",
    }
    n_feats = [r[k]["n_features"] for k in order]
    macro_f1 = [r[k]["macro_f1"] for k in order]
    bal_acc = [r[k]["bal_acc"] for k in order]
    colors = ["#8bc34a", "#ffb74d", "#0d47a1", "#1976d2", "#7b1fa2", "#c62828"]

    baseline = r["cos32_only (paper baseline)"]["macro_f1"]
    all_combined = r["cos32 + EXP3 + scalars (all)"]["macro_f1"]

    fig, ax = plt.subplots(figsize=(10, 5))
    y = np.arange(len(order))
    width = 0.35
    bars_f1 = ax.barh(y - width / 2, macro_f1, width, color=colors, edgecolor="black", lw=0.5, label="macro-F1")
    bars_ba = ax.barh(y + width / 2, bal_acc, width, color=[c + "80" for c in colors], edgecolor="black", lw=0.5, label="balanced accuracy")
    for yi, (f, b, n) in enumerate(zip(macro_f1, bal_acc, n_feats)):
        ax.text(f + 0.005, yi - width / 2, f"{f:.3f}", va="center", fontsize=8.5, fontweight="bold")
        ax.text(b + 0.005, yi + width / 2, f"{b:.3f}  (n_feat={n})", va="center", fontsize=8)
    ax.axvline(0.25, ls="--", color="#555", lw=0.8, label="chance (4-class)")
    ax.axvline(baseline, ls=":", color="#0d47a1", lw=1.2, label=f"paper baseline = {baseline:.3f}")
    ax.set_yticks(y)
    ax.set_yticklabels([labels_short[k] for k in order])
    ax.invert_yaxis()
    ax.set_xlim(0, 0.85)
    ax.set_xlabel("score (10-fold stratified CV, seed 42, n=6,191)")
    ax.set_title(f"EXP2 §4 — Adding EXP3 features improves 4-class subtype classification\n"
                 f"Combined (cos32 + EXP3 + scalars) reaches macro-F1 = {all_combined:.3f} (Δ = {all_combined - baseline:+.3f} over paper 32-d baseline)",
                 loc="left", fontsize=11, fontfamily="Times New Roman", fontweight="bold", pad=8)
    ax.legend(loc="lower right", fontsize=8)
    ax.grid(axis="x", color="#dddddd", lw=0.5)

    FIG.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / "EXP2_H2b_integrated.png")
    fig.savefig(FIG / "EXP2_H2b_integrated.pdf")
    plt.close(fig)
    print(f"wrote {FIG / 'EXP2_H2b_integrated.png'}")


if __name__ == "__main__":
    main()
