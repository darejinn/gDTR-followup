"""EXP1 §Phase C — WGS-wide cCRE-ELS join (paper Fig 2 panel b replication).

Loads per-chromosome position parquets, masks by ENCODE SCREEN V3 cCRE-ELS BED
(both dELS and pELS subtypes), and computes per-chromosome and WGS-pooled c_t / oscil
means vs intron baseline.

Reads:
  wgs/results/chr{N}/chr{N}_per_position_chunk*.parquet
  wgs/data/labels/chr{N}_position_labels.npy
  data_ref/ccre_els_wgs.bed  (809k regions, all 24 chr)
Writes:
  wgs/results/genome_summary/wgs_ccre_els_context.csv
  wgs/results/genome_summary/wgs_ccre_els_summary.json
  wgs/results/genome_summary/wgs_ccre_els_figure.{png,pdf}
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
BED = TDIG / "data_ref/ccre_els_wgs.bed"
OUT = WGS / "genome_summary"

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

# Paper §3.1 panel b: cCRE-ELS (chr22 only) d = -0.118 (unified pooled intron baseline)
PAPER_CCRE_D = -0.118


def load_bed_mask(bed_path: Path, chrom: str, chrom_len: int, subtype_filter: set[str] | None = None) -> np.ndarray:
    mask = np.zeros(chrom_len, dtype=bool)
    n = 0
    with bed_path.open() as f:
        for ln in f:
            if not ln.strip() or ln.startswith("#") or ln.startswith("track"):
                continue
            fs = ln.split("\t")
            if len(fs) < 3 or fs[0] != chrom: continue
            # subtype column (10 or 11): "dELS,CTCF-bound" or "pELS" etc.
            if subtype_filter is not None and len(fs) > 10:
                sub = fs[10].strip()  # short subtype
                if not any(s in sub for s in subtype_filter):
                    continue
            try:
                a, b = int(fs[1]), int(fs[2])
            except ValueError: continue
            a = max(0, a); b = min(chrom_len, b)
            if 0 <= a < b <= chrom_len:
                mask[a:b] = True
                n += 1
    return mask, n


def main():
    t0 = time.time()
    rows = []
    for chrom in CHROMS:
        t_c = time.time()
        lab = np.load(LABEL_DIR / f"{chrom}_position_labels.npy")
        chrom_len = len(lab)

        # Load cCRE mask (ELS subtypes: dELS + pELS)
        ccre_mask, n_ccre_records = load_bed_mask(BED, chrom, chrom_len, subtype_filter={"dELS", "pELS"})
        ccre_bp = int(ccre_mask.sum())

        # Load position parquets and compute stats for ccre positions and intron baseline
        chunks = sorted((WGS / chrom).glob(f"{chrom}_per_position_chunk*.parquet"))
        ccre_c_t, ccre_osc = [], []
        intron_c_t, intron_osc = [], []
        for cp in chunks:
            d = pd.read_parquet(cp, columns=["pos", "c_t", "oscil"])
            pos = d["pos"].to_numpy()
            # ccre positions (in cCRE-ELS AND not classified as splice/coding by GENCODE mask)
            # But paper does not exclude — just intersect
            is_ccre = ccre_mask[pos]
            is_intron = lab[pos] == 1  # intron
            ccre_c_t.extend(d.loc[is_ccre, "c_t"].to_numpy())
            ccre_osc.extend(d.loc[is_ccre, "oscil"].to_numpy())
            intron_c_t.extend(d.loc[is_intron, "c_t"].to_numpy())
            intron_osc.extend(d.loc[is_intron, "oscil"].to_numpy())

        ccre_c_t = np.asarray(ccre_c_t, dtype=np.float32)
        ccre_osc = np.asarray(ccre_osc, dtype=np.float32)
        intron_c_t = np.asarray(intron_c_t, dtype=np.float32)
        intron_osc = np.asarray(intron_osc, dtype=np.float32)

        n_ccre = len(ccre_c_t); n_intron = len(intron_c_t)
        if n_ccre < 30 or n_intron < 30:
            print(f"[{chrom}] insufficient positions: ccre={n_ccre}, intron={n_intron}")
            continue

        m_c_ccre, m_c_int = float(ccre_c_t.mean()), float(intron_c_t.mean())
        m_o_ccre, m_o_int = float(ccre_osc.mean()), float(intron_osc.mean())
        v_c_ccre = float(np.var(ccre_c_t, ddof=1)); v_c_int = float(np.var(intron_c_t, ddof=1))
        v_o_ccre = float(np.var(ccre_osc, ddof=1)); v_o_int = float(np.var(intron_osc, ddof=1))
        pooled_c = np.sqrt(((n_ccre-1)*v_c_ccre + (n_intron-1)*v_c_int) / (n_ccre+n_intron-2))
        pooled_o = np.sqrt(((n_ccre-1)*v_o_ccre + (n_intron-1)*v_o_int) / (n_ccre+n_intron-2))
        d_c = (m_c_ccre - m_c_int) / pooled_c if pooled_c else 0.0
        d_o = (m_o_ccre - m_o_int) / pooled_o if pooled_o else 0.0

        rows.append({
            "chrom": chrom,
            "n_ccre_records": n_ccre_records,
            "n_ccre_bp_in_bed": ccre_bp,
            "n_ccre_pos_analyzed": n_ccre,
            "n_intron_pos_analyzed": n_intron,
            "mean_c_t_ccre": m_c_ccre,
            "mean_c_t_intron": m_c_int,
            "cohen_d_c_t_vs_intron": d_c,
            "mean_oscil_ccre": m_o_ccre,
            "mean_oscil_intron": m_o_int,
            "cohen_d_oscil_vs_intron": d_o,
        })
        print(f"[{time.strftime('%H:%M:%S')}] {chrom}: ccre={n_ccre:,} intron={n_intron:,} "
              f"d_c={d_c:+.4f} d_o={d_o:+.4f} {time.time()-t_c:.1f}s")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "wgs_ccre_els_context.csv", index=False)

    # WGS pooled — reconstruct from per-chr n-weighted means (for simplicity)
    n_ccre_total = int(df["n_ccre_pos_analyzed"].sum())
    n_intron_total = int(df["n_intron_pos_analyzed"].sum())
    m_c_ccre_pool = float((df["mean_c_t_ccre"] * df["n_ccre_pos_analyzed"]).sum() / n_ccre_total)
    m_c_int_pool = float((df["mean_c_t_intron"] * df["n_intron_pos_analyzed"]).sum() / n_intron_total)
    m_o_ccre_pool = float((df["mean_oscil_ccre"] * df["n_ccre_pos_analyzed"]).sum() / n_ccre_total)
    m_o_int_pool = float((df["mean_oscil_intron"] * df["n_intron_pos_analyzed"]).sum() / n_intron_total)
    d_c_wavg = float((df["cohen_d_c_t_vs_intron"] * df["n_ccre_pos_analyzed"]).sum() / n_ccre_total)
    d_o_wavg = float((df["cohen_d_oscil_vs_intron"] * df["n_ccre_pos_analyzed"]).sum() / n_ccre_total)

    summary = {
        "n_ccre_records_total": int(df["n_ccre_records"].sum()),
        "n_ccre_bp_total": int(df["n_ccre_bp_in_bed"].sum()),
        "n_ccre_pos_analyzed_total": n_ccre_total,
        "n_intron_pos_analyzed_total": n_intron_total,
        "wgs_mean_c_t_ccre": m_c_ccre_pool,
        "wgs_mean_c_t_intron_baseline": m_c_int_pool,
        "wgs_d_c_t_ccre_vs_intron_wavg": d_c_wavg,
        "wgs_d_c_t_std_across_chr": float(df["cohen_d_c_t_vs_intron"].std(ddof=1)),
        "wgs_mean_oscil_ccre": m_o_ccre_pool,
        "wgs_mean_oscil_intron_baseline": m_o_int_pool,
        "wgs_d_oscil_ccre_vs_intron_wavg": d_o_wavg,
        "wgs_d_oscil_std_across_chr": float(df["cohen_d_oscil_vs_intron"].std(ddof=1)),
        "paper_ccre_d_c_t_chr22_only": PAPER_CCRE_D,
        "delta_wgs_vs_paper_d_c_t": d_c_wavg - PAPER_CCRE_D,
        "runtime_sec": time.time() - t0,
    }
    (OUT / "wgs_ccre_els_summary.json").write_text(json.dumps(summary, indent=2, default=str))

    print("\n=== WGS cCRE-ELS summary ===")
    for k, v in summary.items():
        if isinstance(v, float):
            print(f"  {k}: {v:.4f}")
        else:
            print(f"  {k}: {v:,}" if isinstance(v, int) else f"  {k}: {v}")

    # Figure — per-chr d_c_t + WGS bar
    plt.rcParams.update({"font.family": "DejaVu Sans", "savefig.dpi": 300, "savefig.bbox": "tight",
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    ax.bar(range(len(df)), df["cohen_d_c_t_vs_intron"], color="#d62728",
           edgecolor="black", lw=0.4, label="d(c_t vs intron)")
    ax.axhline(PAPER_CCRE_D, ls="--", color="#0d47a1", lw=1.0,
               label=f"paper (chr22-only) = {PAPER_CCRE_D:+.3f}")
    ax.axhline(d_c_wavg, ls=":", color="#c62828", lw=1.0,
               label=f"WGS weighted avg = {d_c_wavg:+.4f}")
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df["chrom"], rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Cohen's d — c(t) vs intron baseline")
    ax.set_title("(a) Per-chromosome cCRE-ELS d(c_t)",
                 loc="left", fontweight="bold", fontfamily="Times New Roman", fontsize=12, pad=8)
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(axis="y", color="#dddddd", lw=0.5)

    ax = axes[1]
    ax.bar(range(len(df)), df["cohen_d_oscil_vs_intron"], color="#1f77b4",
           edgecolor="black", lw=0.4, label="d(oscil vs intron)")
    ax.axhline(0, color="#555", lw=0.6)
    ax.axhline(d_o_wavg, ls=":", color="#1976d2", lw=1.0,
               label=f"WGS weighted avg = {d_o_wavg:+.4f}")
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df["chrom"], rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Cohen's d — oscil vs intron baseline")
    ax.set_title("(b) Per-chromosome cCRE-ELS d(oscil) — H3a extension",
                 loc="left", fontweight="bold", fontfamily="Times New Roman", fontsize=12, pad=8)
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(axis="y", color="#dddddd", lw=0.5)

    fig.suptitle(f"cCRE-ELS WGS join — 24 chromosomes, {n_ccre_total:,} ELS positions vs {n_intron_total/1e6:.0f}M intron",
                 fontsize=13, fontweight="bold", fontfamily="Times New Roman", y=1.02)
    fig.tight_layout()
    fig.savefig(OUT / "wgs_ccre_els_figure.png")
    fig.savefig(OUT / "wgs_ccre_els_figure.pdf")
    plt.close(fig)
    print(f"wrote {OUT/'wgs_ccre_els_figure.png'}")

    print(f"\nTotal runtime: {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
