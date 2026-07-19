"""EXP3 §3.4 — H3a: does oscil (or n_enter, n_exit) differ by GENCODE context?

Hypothesis (from user + plan.md):
  Splice donors / acceptors show CLEANER single-crossings (oscil ≈ 0) than intronic
  positions, revealing an orthogonal readout beyond monotone c(t).

  We test whether context labels {intron, splice_donor, splice_acceptor, coding_exon,
  intergenic} have significantly different oscil, n_enter, and below_frac distributions.

Reads:  results/A_crossing_stats.parquet
Writes: results/H3a_context_test.json
"""
from __future__ import annotations
import json, time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

TDIG = Path("/NHNHOME/WORKSPACE/0526040123_A/darejinn/tdig")
PARQUET = TDIG / "exp3_threshold_crossing/results/A_crossing_stats.parquet"
OUT_JSON = TDIG / "exp3_threshold_crossing/results/H3a_context_test.json"

# Same codebook as paper's tA_forward.py / data_ref/*_position_labels.npy
LABEL_MAP = {
    0: "intergenic",
    1: "intron",
    2: "coding_exon",
    3: "5utr",
    4: "3utr",
    5: "splice_donor",
    6: "splice_acceptor",
}


def cohen_d_indep(x: np.ndarray, y: np.ndarray) -> float:
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2: return float("nan")
    vx, vy = float(np.var(x, ddof=1)), float(np.var(y, ddof=1))
    pooled_sd = np.sqrt(((nx - 1) * vx + (ny - 1) * vy) / (nx + ny - 2))
    return float((np.mean(x) - np.mean(y)) / pooled_sd) if pooled_sd else 0.0


def main() -> None:
    t0 = time.time()
    df = pd.read_parquet(PARQUET)
    df["context"] = df["label"].map(LABEL_MAP).fillna("unknown")
    print(f"total positions: {len(df):,}")
    print(f"context counts: {df['context'].value_counts().to_dict()}")

    # Filter to contexts of interest (skip 5utr, 3utr for now, small n on chr22 100-window sample)
    keep = ["intron", "splice_donor", "splice_acceptor", "coding_exon", "intergenic"]

    metrics = ["oscil", "n_enter", "n_exit", "below_frac", "c_t", "min_D"]

    ref_context = "intron"  # baseline like paper §3.1
    ref = df[df["context"] == ref_context]
    print(f"\nreference '{ref_context}' n={len(ref):,}")

    results = {"reference": ref_context, "n_reference": int(len(ref)), "contexts": {}}
    for ctx in keep:
        if ctx == ref_context: continue
        sub = df[df["context"] == ctx]
        if len(sub) < 30:
            print(f"skip {ctx} (n={len(sub)} < 30)")
            continue
        entry = {"n": int(len(sub))}
        for m in metrics:
            x, y = sub[m].to_numpy(dtype=float), ref[m].to_numpy(dtype=float)
            d = cohen_d_indep(x, y)
            try:
                u, p = stats.mannwhitneyu(x, y, alternative="two-sided")
            except ValueError:
                u, p = np.nan, np.nan
            entry[m] = {
                "mean": float(x.mean()),
                "median": float(np.median(x)),
                "ref_mean": float(y.mean()),
                "cohen_d_vs_ref": d,
                "mwu_p": float(p) if p is not None else None,
                "delta_mean": float(x.mean() - y.mean()),
            }
        results["contexts"][ctx] = entry

    # Also compute per-context oscil distribution (histograms up to oscil=6)
    osc_hist = {}
    for ctx in keep:
        sub = df[df["context"] == ctx]
        vc = sub["oscil"].value_counts().sort_index().to_dict()
        osc_hist[ctx] = {int(k): int(v) for k, v in vc.items()}
    results["oscil_distribution_per_context"] = osc_hist

    # Verdict per H3a criteria
    verdicts = {}
    for ctx in ["splice_donor", "splice_acceptor"]:
        if ctx not in results["contexts"]: continue
        e = results["contexts"][ctx]
        d = e["oscil"]["cohen_d_vs_ref"]
        p = e["oscil"]["mwu_p"]
        verdicts[ctx] = {
            "H3a_positive": bool(abs(d) >= 0.15 and (p is not None and p < 1e-6)),
            "cohen_d": d,
            "p_value": p,
        }
    results["verdict_H3a"] = verdicts

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(results, indent=2, default=str))

    print("\n=== H3a summary vs intron baseline ===")
    print(f"{'context':20s} {'n':>8s} {'mean_oscil':>10s} {'d_oscil':>10s} {'p_oscil':>10s} {'mean_c_t':>10s} {'d_c_t':>10s}")
    ref_mean_osc = ref["oscil"].mean()
    ref_mean_ct = ref["c_t"].mean()
    print(f"{ref_context:20s} {len(ref):>8d} {ref_mean_osc:>10.4f} {'0':>10s} {'—':>10s} {ref_mean_ct:>10.4f} {'0':>10s}")
    for ctx, e in results["contexts"].items():
        print(f"{ctx:20s} {e['n']:>8d} {e['oscil']['mean']:>10.4f} {e['oscil']['cohen_d_vs_ref']:>+10.4f} {e['oscil']['mwu_p']:>10.2e} {e['c_t']['mean']:>10.4f} {e['c_t']['cohen_d_vs_ref']:>+10.4f}")

    print(f"\nH3a verdicts:")
    for ctx, v in verdicts.items():
        icon = "✓" if v["H3a_positive"] else "✗"
        print(f"  {icon} {ctx}: cohen_d={v['cohen_d']:+.4f}, p={v['p_value']:.2e}")

    print(f"\nWrote: {OUT_JSON}")
    print(f"Runtime: {time.time()-t0:.2f} s")


if __name__ == "__main__":
    main()
