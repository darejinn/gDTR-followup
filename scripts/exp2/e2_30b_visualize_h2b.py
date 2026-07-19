"""EXP2 §3b — visualize H2b: subtype confusion + OvR AUROC bar."""
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
    summ = json.loads((RES / "H2b_subtype_multi.json").read_text())
    cm = np.asarray(summ["confusion_matrix_32d"])
    classes = summ["class_order"]
    aucs = summ["auroc_per_class_OvR_32d"]

    # Row-normalize confusion matrix for readability
    cm_norm = cm.astype(float) / cm.sum(axis=1, keepdims=True)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), gridspec_kw={"width_ratios": [1.15, 1.0]})

    # Panel (a): confusion matrix
    ax = axes[0]
    im = ax.imshow(cm_norm, cmap="Blues", vmin=0, vmax=1)
    for i in range(len(classes)):
        for j in range(len(classes)):
            v = cm_norm[i, j]
            ax.text(j, i, f"{v:.2f}\n({cm[i,j]})", ha="center", va="center",
                    color="white" if v > 0.5 else "black", fontsize=9)
    ax.set_xticks(range(len(classes))); ax.set_yticks(range(len(classes)))
    ax.set_xticklabels([c.replace("_", "\n") for c in classes], rotation=0)
    ax.set_yticklabels([c.replace("_", "\n") for c in classes])
    ax.set_xlabel("predicted class")
    ax.set_ylabel("true class")
    ax.set_title(f"(a) Row-normalized confusion matrix\n(macro-F1 = {summ['macro_f1_32d']:.3f}, bal-acc = {summ['balanced_acc_32d']:.3f})",
                 loc="left", fontsize=12, fontweight="bold", fontfamily="Times New Roman", pad=8)
    plt.colorbar(im, ax=ax, fraction=0.045, label="row fraction")

    # Panel (b): OvR AUROC per class bar
    ax = axes[1]
    class_labels_sorted = sorted(aucs, key=aucs.get, reverse=True)
    aur_values = [aucs[c] for c in class_labels_sorted]
    n_values = [summ["class_counts"][c] for c in class_labels_sorted]
    colors = ["#1f77b4", "#2ca02c", "#d62728", "#ff7f0e"][:len(class_labels_sorted)]
    y_pos = np.arange(len(class_labels_sorted))
    ax.barh(y_pos, aur_values, color=colors, edgecolor="black", lw=0.5, height=0.62)
    for yi, (v, n) in enumerate(zip(aur_values, n_values)):
        ax.text(v + 0.01, yi, f"{v:.3f}   (n={n})", va="center", fontsize=9, fontweight="bold")
    ax.axvline(0.5, color="#888", lw=0.6, ls="--")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(class_labels_sorted)
    ax.invert_yaxis()
    ax.set_xlim(0.5, 1.0)
    ax.set_xlabel("OvR AUROC (10-fold CV, 32-d $\\Delta D_{cos}$)")
    ax.set_title("(b) One-vs-rest AUROC per class",
                 loc="left", fontsize=12, fontweight="bold", fontfamily="Times New Roman", pad=8)
    ax.grid(axis="x", color="#dddddd", lw=0.5)

    # Big finding text
    fig.suptitle(f"EXP2 H2b — variant subtype classification from $\\Delta D_{{cos}}$ trajectories\n"
                 f"32-d macro-F1 = {summ['macro_f1_32d']:.3f}  vs  1-d (max$|\\Delta D|$) macro-F1 = {summ['macro_f1_1d_max_abs_dD']:.3f}  vs  chance = {1/len(classes):.3f}",
                 fontsize=11, y=1.02, fontfamily="Times New Roman", fontweight="bold")

    FIG.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / "EXP2_H2b_subtype.png")
    fig.savefig(FIG / "EXP2_H2b_subtype.pdf")
    plt.close(fig)
    print(f"wrote {FIG / 'EXP2_H2b_subtype.png'}")


if __name__ == "__main__":
    main()
