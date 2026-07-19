"""EXP3 §3.6 — Advanced threshold-crossing pattern analysis.

Beyond n_enter / n_exit / oscil (already computed in e3_10), we now compute richer
features that characterize HOW the trajectory crosses γ, not just how many times:

  first_enter_layer:    layer index of the FIRST above→below transition (or -1 if never)
  last_exit_layer:      layer index of the LAST below→above transition (or -1 if never leaves)
  longest_below_streak: max # of consecutive layers where D_cos ≤ γ
  streak_start_layer:   layer where the longest streak begins
  amplitude_at_min:     γ - min(D_cos) (positive = how far below γ it dips)
  mid_layer_below_frac: fraction of layers 10..21 (middle third) below γ  (structural "commitment phase")
  early_below_frac:     fraction of layers 0..10 below γ
  late_below_frac:      fraction of layers 22..31 below γ
  d_at_gamma_slope:     mean slope of D_cos around the first crossing (rate of commitment)

Reads:  exp3_threshold_crossing/results/A_windows_chr22.npz
        wgs/results/chr1/chr1_per_position_chunk*.parquet (optional, second data source)
Writes: exp3_threshold_crossing/results/advanced_crossings_chr22.parquet
        exp3_threshold_crossing/results/H3d_advanced_context_test.json
        exp3_threshold_crossing/figures/EXP3_H3d_advanced_patterns.png
"""
from __future__ import annotations
import json, time, argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TDIG = Path("/NHNHOME/WORKSPACE/0526040123_A/darejinn/tdig")
NPZ = TDIG / "exp3_threshold_crossing/results/A_windows_chr22.npz"
OUT_PARQUET = TDIG / "exp3_threshold_crossing/results/advanced_crossings_chr22.parquet"
OUT_JSON = TDIG / "exp3_threshold_crossing/results/H3d_advanced_context_test.json"
FIG_DIR = TDIG / "exp3_threshold_crossing/figures"

GAMMA = 0.397
N_LAYERS = 32
EARLY_END = 10   # layers 0..9 = "early"
MID_END = 22     # layers 10..21 = "middle"
                 # layers 22..31 = "late"

LABEL_MAP = {0: "intergenic", 1: "intron", 2: "coding_exon", 3: "5utr", 4: "3utr",
             5: "splice_donor", 6: "splice_acceptor"}


def longest_below_streak(below: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Given below: [L, T] bool, return (streak_len [T], streak_start [T]) for the LONGEST
    consecutive True run per position."""
    L, T = below.shape
    best_len = np.zeros(T, dtype=np.int16)
    best_start = np.full(T, -1, dtype=np.int16)
    cur_len = np.zeros(T, dtype=np.int16)
    cur_start = np.full(T, -1, dtype=np.int16)
    for l in range(L):
        # if below at (l, t), extend or start new streak
        b = below[l]
        newly_starting = b & (cur_len == 0)
        cur_start = np.where(newly_starting, l, cur_start)
        cur_len = np.where(b, cur_len + 1, 0)
        # update best
        better = cur_len > best_len
        best_len = np.where(better, cur_len, best_len)
        best_start = np.where(better, cur_start, best_start)
    return best_len, best_start


def compute_advanced(D: np.ndarray, gamma: float = GAMMA) -> dict:
    """D: [N, L, T]. Returns dict of [N, T] arrays with advanced features."""
    N, L, T = D.shape
    below = D <= gamma                       # [N, L, T]
    # Flatten to [L, T*N] for streak computation? Simpler: iterate over N.
    first_enter = np.full((N, T), -1, dtype=np.int16)
    last_exit = np.full((N, T), -1, dtype=np.int16)
    longest_len = np.zeros((N, T), dtype=np.int16)
    longest_start = np.full((N, T), -1, dtype=np.int16)
    for n in range(N):
        below_n = below[n]                   # [L, T]
        # first_enter: for each t, first l where below_n[l, t] is True
        any_below = below_n.any(axis=0)      # [T]
        first_below_layer = below_n.argmax(axis=0)  # [T] first True index or 0 if all False
        first_enter[n] = np.where(any_below, first_below_layer, -1)
        # last_exit: reverse-argmax on ~below (i.e., last True in above)
        # more directly: last transition from below → above
        exit_transitions = below_n[:-1] & ~below_n[1:]   # [L-1, T] — TRUE only where exit happens
        any_exit = exit_transitions.any(axis=0)          # [T]
        last_exit_layer = L - 2 - exit_transitions[::-1].argmax(axis=0)  # 0-based layer BEFORE the exit
        last_exit[n] = np.where(any_exit, last_exit_layer + 1, -1)  # layer AT which we exited (0-based)
        ln, ls = longest_below_streak(below_n)
        longest_len[n], longest_start[n] = ln, ls

    # Layer-band below fractions
    early = below[:, :EARLY_END, :].mean(axis=1).astype(np.float32)   # [N, T]
    mid = below[:, EARLY_END:MID_END, :].mean(axis=1).astype(np.float32)
    late = below[:, MID_END:, :].mean(axis=1).astype(np.float32)

    # amplitude at min
    min_D = D.min(axis=1)                    # [N, T]
    amplitude = np.clip(gamma - min_D, 0, None).astype(np.float32)

    return {
        "first_enter_layer": first_enter,
        "last_exit_layer": last_exit,
        "longest_below_streak": longest_len,
        "streak_start_layer": longest_start,
        "amplitude_below_gamma": amplitude,
        "early_below_frac": early,
        "mid_below_frac": mid,
        "late_below_frac": late,
        "min_D": min_D.astype(np.float32),
    }


def cohen_d_indep(x, y):
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2: return float("nan")
    vx, vy = float(np.var(x, ddof=1)), float(np.var(y, ddof=1))
    pooled = np.sqrt(((nx-1)*vx + (ny-1)*vy) / (nx+ny-2))
    return float((np.mean(x)-np.mean(y))/pooled) if pooled else 0.0


def main():
    t0 = time.time()
    print(f"loading {NPZ} …")
    data = np.load(NPZ)
    D = data["D_cos"].astype(np.float32)
    positions = data["genomic_positions"].astype(np.int64)
    labels = data["labels"].astype(np.uint8)
    N, L, T = D.shape
    print(f"  D shape: {D.shape}")

    feats = compute_advanced(D, GAMMA)
    print(f"computed {len(feats)} features in {time.time()-t0:.2f}s")

    # Flatten to per-position parquet
    flat_win = np.broadcast_to(np.arange(N).reshape(-1, 1), (N, T)).ravel()
    flat_pos = positions.ravel()
    flat_lab = labels.ravel()
    df = pd.DataFrame({
        "win_idx": flat_win.astype(np.int32),
        "pos": flat_pos,
        "label": flat_lab,
        "first_enter_layer": feats["first_enter_layer"].ravel(),
        "last_exit_layer": feats["last_exit_layer"].ravel(),
        "longest_below_streak": feats["longest_below_streak"].ravel(),
        "streak_start_layer": feats["streak_start_layer"].ravel(),
        "amplitude_below_gamma": feats["amplitude_below_gamma"].ravel(),
        "early_below_frac": feats["early_below_frac"].ravel(),
        "mid_below_frac": feats["mid_below_frac"].ravel(),
        "late_below_frac": feats["late_below_frac"].ravel(),
        "min_D": feats["min_D"].ravel(),
    })
    df["context"] = df["label"].map(LABEL_MAP).fillna("unknown")

    OUT_PARQUET.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(OUT_PARQUET, index=False)
    print(f"wrote {OUT_PARQUET} ({df.shape})")

    # Per-context test — compare each advanced feature to intron baseline
    metrics = ["first_enter_layer", "last_exit_layer", "longest_below_streak",
               "streak_start_layer", "amplitude_below_gamma",
               "early_below_frac", "mid_below_frac", "late_below_frac", "min_D"]
    contexts = ["intron", "splice_donor", "splice_acceptor", "coding_exon", "intergenic"]
    intron = df[df["context"] == "intron"]
    results = {"n_intron": int(len(intron)), "contexts": {}}
    for ctx in contexts:
        if ctx == "intron": continue
        sub = df[df["context"] == ctx]
        if len(sub) < 30: continue
        entry = {"n": int(len(sub))}
        for m in metrics:
            # Filter out -1 sentinel for layer features
            if "layer" in m:
                x_full = sub[m].to_numpy(dtype=float)
                y_full = intron[m].to_numpy(dtype=float)
                x = x_full[x_full >= 0]
                y = y_full[y_full >= 0]
            else:
                x = sub[m].to_numpy(dtype=float)
                y = intron[m].to_numpy(dtype=float)
            if len(x) < 30:
                entry[m] = {"n_valid": int(len(x)), "note": "insufficient valid"}
                continue
            d = cohen_d_indep(x, y)
            try:
                _, p = stats.mannwhitneyu(x, y, alternative="two-sided")
            except ValueError:
                p = float("nan")
            entry[m] = {
                "n_valid": int(len(x)),
                "mean": float(x.mean()),
                "ref_mean": float(y.mean()),
                "cohen_d_vs_intron": d,
                "mwu_p": float(p),
            }
        results["contexts"][ctx] = entry

    OUT_JSON.write_text(json.dumps(results, indent=2, default=str))

    # Summary print
    print("\n=== Advanced pattern features: Cohen's d vs intron ===")
    header = f"{'ctx':16s} {'metric':22s} {'d':>8s} {'p':>10s} {'mean_ctx':>10s} {'mean_intron':>12s}"
    print(header); print("-" * len(header))
    for ctx, e in results["contexts"].items():
        for m in metrics:
            r = e.get(m, {})
            if "cohen_d_vs_intron" not in r: continue
            print(f"{ctx:16s} {m:22s} {r['cohen_d_vs_intron']:>+8.4f} {r['mwu_p']:>10.2e} "
                  f"{r['mean']:>10.4f} {r['ref_mean']:>12.4f}")
        print()

    # Figure — grouped bar chart of Cohen's d for splice sites on all advanced metrics
    plt.rcParams.update({"font.family": "DejaVu Sans", "savefig.dpi": 300, "savefig.bbox": "tight",
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(11, 5))
    metrics_show = ["first_enter_layer", "last_exit_layer", "longest_below_streak",
                    "amplitude_below_gamma", "early_below_frac", "mid_below_frac", "late_below_frac"]
    x = np.arange(len(metrics_show))
    w = 0.20
    plotted_ctx = ["splice_donor", "splice_acceptor", "coding_exon", "intergenic"]
    colors = {"splice_donor": "#1f77b4", "splice_acceptor": "#2ca02c",
              "coding_exon": "#78909c", "intergenic": "#8e8e93"}
    for i, ctx in enumerate(plotted_ctx):
        if ctx not in results["contexts"]: continue
        vals = [results["contexts"][ctx].get(m, {}).get("cohen_d_vs_intron", 0) for m in metrics_show]
        ax.bar(x + (i - 1.5) * w, vals, w, color=colors[ctx], edgecolor="black", lw=0.4, label=ctx)
    ax.axhline(0, color="#555", lw=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels([m.replace("_", "\n") for m in metrics_show], rotation=0, fontsize=8)
    ax.set_ylabel("Cohen's d vs intron")
    ax.set_title("EXP3 §3.6 — Advanced crossing pattern features by context (chr22 79 windows)",
                 fontweight="bold", fontfamily="Times New Roman", fontsize=12)
    ax.legend(loc="upper right", ncol=2, fontsize=9)
    ax.grid(axis="y", color="#dddddd", lw=0.5)

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "EXP3_H3d_advanced_patterns.png")
    fig.savefig(FIG_DIR / "EXP3_H3d_advanced_patterns.pdf")
    plt.close(fig)
    print(f"wrote {FIG_DIR / 'EXP3_H3d_advanced_patterns.png'}")

    print(f"\nTotal runtime: {time.time()-t0:.2f}s")


if __name__ == "__main__":
    main()
