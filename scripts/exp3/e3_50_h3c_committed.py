"""EXP3 §3.4b — H3c: committed vs deliberating stratification.

Reads:  results/A_crossing_stats.parquet
Writes: results/H3c_committed.json + figures/EXP3_H3c_committed_barchart.png

Hypothesis (from plan.md):
  Positions with n_enter=1 AND stable below-threshold suffix ("committed by layer ℓ")
  are enriched in canonical elements (splice sites, TSS-adjacent) vs positions with
  multiple enter events ("still deliberating").
"""
from __future__ import annotations
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

TDIG = Path("/NHNHOME/WORKSPACE/0526040123_A/darejinn/tdig")
PARQUET = TDIG / "exp3_threshold_crossing/results/A_crossing_stats.parquet"
OUT_JSON = TDIG / "exp3_threshold_crossing/results/H3c_committed.json"
FIG = TDIG / "exp3_threshold_crossing/figures"

LABEL_MAP = {0: "intergenic", 1: "intron", 2: "coding_exon", 3: "5utr", 4: "3utr",
             5: "splice_donor", 6: "splice_acceptor"}


def main():
    df = pd.read_parquet(PARQUET)
    df["context"] = df["label"].map(LABEL_MAP).fillna("unknown")

    # Analyze only positions that crossed threshold at all
    crossed = df[df["first_below"] >= 0].copy()
    print(f"total positions: {len(df):,}")
    print(f"positions with any crossing: {len(crossed):,} ({100*len(crossed)/len(df):.1f}%)")

    # 3 strata:
    #   committed:    n_enter=1 AND n_exit=0 (single clean crossing, never leaves)
    #   deliberating: n_enter+n_exit ≥ 2 (multi-crossings; specifically oscil ≥ 1)
    #   dips:         n_enter=1 AND n_exit=1 (dipped in then out; borderline)
    df["stratum"] = pd.cut(
        pd.Series(np.where(df["first_below"] < 0, -1,
                           np.where((df["n_enter"] == 1) & (df["n_exit"] == 0), 0,
                                    np.where((df["n_enter"] == 1) & (df["n_exit"] == 1), 1, 2))),
                  index=df.index),
        bins=[-1.5, -0.5, 0.5, 1.5, 2.5],
        labels=["never_crossed", "committed", "dipped", "deliberating"],
    )

    print("\nstrata counts per context:")
    ctb = pd.crosstab(df["context"], df["stratum"])
    print(ctb)

    # Fisher tests: splice_donor + splice_acceptor enriched for "deliberating"?
    contexts = ["intergenic", "intron", "coding_exon", "splice_donor", "splice_acceptor"]
    strata = ["committed", "dipped", "deliberating"]

    results = {"strata_counts": ctb.to_dict()}

    fig_data = []
    for ctx in contexts:
        sub = df[df["context"] == ctx]
        n_ctx = len(sub)
        for strat in strata:
            n_ctx_strat = int((sub["stratum"] == strat).sum())
            fig_data.append({"context": ctx, "stratum": strat, "frac": n_ctx_strat / n_ctx if n_ctx else 0, "n": n_ctx_strat})
    fig_df = pd.DataFrame(fig_data)

    # Fisher OR: for each context and each stratum, compute the OR vs intron
    intron = df[df["context"] == "intron"]
    n_intron = len(intron)
    or_results = {}
    for ctx in contexts:
        if ctx == "intron": continue
        sub = df[df["context"] == ctx]
        n_sub = len(sub)
        for strat in strata:
            a = int((sub["stratum"] == strat).sum())
            b = n_sub - a
            c = int((intron["stratum"] == strat).sum())
            d = n_intron - c
            try:
                oddsratio, pv = stats.fisher_exact([[a, b], [c, d]])
            except ValueError:
                oddsratio, pv = float("nan"), float("nan")
            or_results.setdefault(ctx, {})[strat] = {
                "n_ctx_stratum": a, "n_ctx_other": b,
                "n_intron_stratum": c, "n_intron_other": d,
                "odds_ratio": float(oddsratio), "p_value": float(pv),
                "frac_ctx": a / n_sub, "frac_intron": c / n_intron,
            }

    results["fisher_or_vs_intron"] = or_results

    # verdict
    verdicts = {}
    for ctx in ["splice_donor", "splice_acceptor"]:
        if ctx not in or_results: continue
        d = or_results[ctx].get("deliberating", {})
        verdicts[ctx + "_deliberating_enriched"] = bool(
            d.get("odds_ratio", 1.0) >= 2.0 and d.get("p_value", 1.0) < 1e-6)
    results["verdicts"] = verdicts

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(results, indent=2, default=str))

    # print summary
    print("\n=== H3c summary — Fisher OR vs intron for 'deliberating' stratum ===")
    for ctx in ["splice_donor", "splice_acceptor", "coding_exon", "intergenic"]:
        if ctx not in or_results: continue
        d = or_results[ctx].get("deliberating", {})
        print(f"  {ctx:20s}: OR={d.get('odds_ratio', float('nan')):>7.2f} p={d.get('p_value', float('nan')):>10.2e} "
              f"frac_ctx={d.get('frac_ctx', 0)*100:.2f}% frac_intron={d.get('frac_intron', 0)*100:.2f}%")

    # figure
    plt.rcParams.update({"font.family": "DejaVu Sans", "savefig.dpi": 300, "savefig.bbox": "tight"})
    fig, ax = plt.subplots(figsize=(9, 4))
    strat_colors = {"committed": "#2e7d32", "dipped": "#fbc02d", "deliberating": "#c62828"}
    contexts_show = contexts
    x = np.arange(len(contexts_show))
    width = 0.28
    for i, strat in enumerate(strata):
        vals = [fig_df[(fig_df["context"] == c) & (fig_df["stratum"] == strat)]["frac"].iloc[0] * 100
                for c in contexts_show]
        ax.bar(x + (i - 1) * width, vals, width, color=strat_colors[strat],
               edgecolor="black", lw=0.4, label=strat)
        for xi, v in zip(x, vals):
            ax.text(xi + (i - 1) * width, v + 0.05, f"{v:.2f}%",
                    ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels([c.replace("_", "\n") for c in contexts_show])
    ax.set_ylabel("% of context positions")
    ax.set_title("EXP3 H3c — Committed vs deliberating positions by context",
                 fontfamily="Times New Roman", fontweight="bold")
    ax.legend(loc="upper right")
    ax.grid(axis="y", color="#dddddd", lw=0.5)
    FIG.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / "EXP3_H3c_committed_barchart.png")
    fig.savefig(FIG / "EXP3_H3c_committed_barchart.pdf")
    plt.close(fig)
    print(f"wrote {FIG / 'EXP3_H3c_committed_barchart.png'}")


if __name__ == "__main__":
    main()
