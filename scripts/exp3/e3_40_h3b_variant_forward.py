"""EXP3 §3.2 — H3b: variant re-forward. Compute Δoscil for 8,008 SNVs (ref+alt pairs).

Reads variants from results_cached/phase3_ensemble/variants_features_full.csv (paper's App C
cohort). For each SNV, we extract a 6 kb window centered on the variant, tokenize both
ref and alt sequences, forward through Evo 2 7B, capture raw D_cos [32 layers × 6000 pos],
and extract:
  - oscil_ref, oscil_alt at the variant position (central pos, index=3000)
  - n_enter_ref/alt, n_exit_ref/alt
  - Δoscil = oscil_alt - oscil_ref
  - Δn_enter, Δn_exit

Runs on 2× B200 with Evo 2 auto-shard (steady-state ~0.10s per 6000 bp window).
16,016 forwards × 0.10s ≈ 27 minutes.

Output: exp3_threshold_crossing/results/B_variants_oscil.parquet
"""
from __future__ import annotations
import os, time, argparse, logging
from pathlib import Path

os.environ.setdefault("HF_HOME", "/NHNHOME/WORKSPACE/0526040123_A/darejinn/tdig/env/hf_cache")

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F

TDIG = Path("/NHNHOME/WORKSPACE/0526040123_A/darejinn/tdig")
VARIANTS_CSV = TDIG / "results_cached/phase3_ensemble/variants_features_full.csv"
HG38_FA = TDIG / "wgs/data/reference/hg38.fa"
OUT_PARQUET = TDIG / "exp3_threshold_crossing/results/B_variants_oscil.parquet"
CACHE_INTERVAL = 500  # flush every N variants

GAMMA_COS = 0.397
WINDOW = 6000
CENTER = WINDOW // 2  # variant sits at index 3000 (0-based) → tokenizer pos = CENTER

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("h3b_fwd")


def cosine_lens(hidden_states, n_layers=32, bos_offset=0):
    ref_device = hidden_states["norm"].device
    h_norm = hidden_states["norm"].float()[:, bos_offset:, :].to(ref_device)
    h_norm_n = F.normalize(h_norm, p=2, dim=-1)
    D = torch.zeros((n_layers, h_norm.shape[1]), dtype=torch.float32, device=ref_device)
    for ell in range(n_layers):
        h_l = hidden_states[f"blocks.{ell}"].float()[:, bos_offset:, :].to(ref_device)
        h_l_n = F.normalize(h_l, p=2, dim=-1)
        D[ell] = (1.0 - (h_l_n * h_norm_n).sum(dim=-1)).clamp(min=0.0).mean(dim=0)
    return D


def crossings_at_pos(D_col: np.ndarray, gamma: float = GAMMA_COS) -> dict:
    """D_col: [L] float32 per-layer D_cos at one position. Returns per-position stats."""
    below = D_col <= gamma
    n_enter = int(below[0]) + int(np.sum(below[1:] & ~below[:-1]))
    n_exit  = int(np.sum(below[:-1] & ~below[1:]))
    oscil   = max(0, n_enter + n_exit - 1)
    below_frac = float(below.mean())
    argmin_layer = int(D_col.argmin())
    return {
        "n_enter": n_enter, "n_exit": n_exit, "oscil": oscil,
        "below_frac": below_frac, "min_D": float(D_col.min()),
        "argmin_layer": argmin_layer,
        "c_t": int(np.argmax(np.minimum.accumulate(D_col) <= gamma) + 1) if (D_col <= gamma).any() else len(D_col),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--start-idx", type=int, default=0)
    args = ap.parse_args()

    log.info("loading variant features …")
    df = pd.read_csv(VARIANTS_CSV,
                     usecols=["chrom", "pos", "ref", "alt", "gene", "category", "stars"])
    is_snv = (df["ref"].str.len() == 1) & (df["alt"].str.len() == 1)
    df = df[is_snv & df["category"].isin(["P_LP", "B_LB"])].reset_index(drop=True)
    log.info(f"cohort: {df.shape}, chroms: {sorted(df['chrom'].unique())}")
    if args.limit:
        df = df.iloc[args.start_idx:args.start_idx + args.limit]
        log.info(f"processing {args.start_idx}..{args.start_idx+len(df)}")

    log.info("loading hg38 …")
    from pyfaidx import Fasta
    fa = Fasta(str(HG38_FA), as_raw=True, sequence_always_upper=True)

    log.info("loading Evo 2 …")
    t_load = time.time()
    from evo2 import Evo2
    model = Evo2("evo2_7b")
    sh = model.model
    log.info(f"  loaded {time.time()-t_load:.1f}s")

    captured = {}
    handles = []
    def hook(name):
        def _hook(m, i, o):
            captured[name] = (o[0] if isinstance(o, tuple) else o).detach()
        return _hook
    for i in range(32):
        handles.append(sh.blocks[i].register_forward_hook(hook(f"blocks.{i}")))
    handles.append(sh.norm.register_forward_hook(hook("norm")))

    rows = []
    t0 = time.time()
    n_ok = 0
    n_skip = 0
    for i, r in df.iterrows():
        chrom = str(r["chrom"]).replace("chr", "")
        try:
            chrom_name = f"chr{chrom}"
            wstart = int(r["pos"]) - CENTER - 1  # 0-based BED
            wend = wstart + WINDOW
            if wstart < 0 or wend > len(fa[chrom_name]):
                n_skip += 1
                continue
            seq = str(fa[chrom_name][wstart:wend])
            if len(seq) != WINDOW or seq.count("N") > 500:
                n_skip += 1
                continue

            # Build ref/alt sequences by substituting central base
            ref_base = seq[CENTER]
            expected_ref = str(r["ref"]).upper()
            if ref_base != expected_ref:
                # sometimes coord off-by-one or strand — try adjacent
                if wstart + CENTER + 1 < len(fa[chrom_name]) and str(fa[chrom_name][wstart + CENTER + 1:wstart + CENTER + 2]) == expected_ref:
                    # shift +1
                    wstart += 1; wend += 1
                    seq = str(fa[chrom_name][wstart:wend])
                elif wstart + CENTER - 1 >= 0 and str(fa[chrom_name][wstart + CENTER - 1:wstart + CENTER]) == expected_ref:
                    wstart -= 1; wend -= 1
                    seq = str(fa[chrom_name][wstart:wend])
                else:
                    n_skip += 1
                    continue
                ref_base = seq[CENTER]
                if ref_base != expected_ref:
                    n_skip += 1
                    continue

            alt_seq = seq[:CENTER] + str(r["alt"]).upper() + seq[CENTER + 1:]

            stats = {}
            for tag, s in [("ref", seq), ("alt", alt_seq)]:
                tokens = model.tokenizer.tokenize(s)
                input_ids = torch.tensor(tokens, dtype=torch.int64).unsqueeze(0).cuda()
                with torch.no_grad():
                    _ = sh(input_ids)
                D = cosine_lens(captured, 32, 0)      # [32, 6000]
                D_at_var = D[:, CENTER].cpu().numpy().astype(np.float32)  # [32] at variant pos
                s = crossings_at_pos(D_at_var, GAMMA_COS)
                stats[tag] = s
                # cleanup
                for k in list(captured.keys()): del captured[k]
                del D
                torch.cuda.empty_cache()

            entry = {
                "chrom": r["chrom"], "pos": int(r["pos"]),
                "ref": r["ref"], "alt": r["alt"], "gene": r["gene"], "category": r["category"],
                "oscil_ref": stats["ref"]["oscil"], "oscil_alt": stats["alt"]["oscil"],
                "d_oscil": stats["alt"]["oscil"] - stats["ref"]["oscil"],
                "n_enter_ref": stats["ref"]["n_enter"], "n_enter_alt": stats["alt"]["n_enter"],
                "d_n_enter": stats["alt"]["n_enter"] - stats["ref"]["n_enter"],
                "n_exit_ref": stats["ref"]["n_exit"], "n_exit_alt": stats["alt"]["n_exit"],
                "d_n_exit": stats["alt"]["n_exit"] - stats["ref"]["n_exit"],
                "below_frac_ref": stats["ref"]["below_frac"], "below_frac_alt": stats["alt"]["below_frac"],
                "d_below_frac": stats["alt"]["below_frac"] - stats["ref"]["below_frac"],
                "c_t_ref": stats["ref"]["c_t"], "c_t_alt": stats["alt"]["c_t"],
                "d_c_t": stats["alt"]["c_t"] - stats["ref"]["c_t"],
                "min_D_ref": stats["ref"]["min_D"], "min_D_alt": stats["alt"]["min_D"],
                "d_min_D": stats["alt"]["min_D"] - stats["ref"]["min_D"],
                "argmin_layer_ref": stats["ref"]["argmin_layer"], "argmin_layer_alt": stats["alt"]["argmin_layer"],
            }
            rows.append(entry)
            n_ok += 1

            if (i + 1) % 100 == 0:
                elapsed = time.time() - t0
                rate = (i + 1) / elapsed
                eta = (len(df) - i - 1) / rate / 60
                log.info(f"  {i+1}/{len(df)} n_ok={n_ok} n_skip={n_skip} rate={rate:.2f}/s ETA={eta:.1f}min")

            if (i + 1) % CACHE_INTERVAL == 0 and rows:
                OUT_PARQUET.parent.mkdir(parents=True, exist_ok=True)
                pd.DataFrame(rows).to_parquet(OUT_PARQUET, index=False)

        except Exception as e:
            log.error(f"variant {r['chrom']}:{r['pos']} {r['ref']}>{r['alt']}: {type(e).__name__}: {e}")
            n_skip += 1

    for h in handles: h.remove()

    OUT_PARQUET.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_parquet(OUT_PARQUET, index=False)
    log.info(f"DONE n_ok={n_ok} n_skip={n_skip} total_time={time.time()-t0:.1f}s → {OUT_PARQUET}")


if __name__ == "__main__":
    main()
