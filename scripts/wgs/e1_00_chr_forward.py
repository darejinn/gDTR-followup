"""EXP1 — genome-scale chr forward pass, saves per-position c_t + oscil + n_enter/exit.

Adapts paper/scripts/tierA/tA_forward.py to arbitrary chromosome, with EXP3-style
crossing stats retained. Uses 2× B200 auto-shard (steady-state ~0.10s per 6 kb window).

For a full chr1 sweep at 6 kb stride 3 kb: ~82k windows × 0.10s ≈ 2.3 hours.

Outputs (chunked, resumable):
  wgs/results/{chr}/{chr}_per_position.parquet [central 3 kb positions × window]
    columns: window_idx, pos, c_t, oscil, n_enter, n_exit, below_frac, min_D, argmin_layer
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
HG38_FA = TDIG / "wgs/data/reference/hg38.fa"

GAMMA_COS = 0.397
WINDOW = 6000
STRIDE = 3000
CENTRAL_START = 1500
CENTRAL_END = 4500

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("wgs_fwd")


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


def crossings_all(D_np: np.ndarray, gamma: float = GAMMA_COS) -> dict:
    """D_np: [L, T]. Returns per-position stats as [T]-length dicts of arrays."""
    below = D_np <= gamma          # [L, T]
    # n_enter: layer 0 counts if below[0], plus each below-transition
    layer0 = below[0:1].astype(np.int8)
    en = below[1:] & ~below[:-1]
    n_enter = layer0.sum(axis=0) + en.sum(axis=0)
    ex = below[:-1] & ~below[1:]
    n_exit = ex.sum(axis=0)
    oscil = np.clip(n_enter + n_exit - 1, 0, None)
    below_frac = below.mean(axis=0).astype(np.float32)
    min_D = D_np.min(axis=0).astype(np.float32)
    argmin_layer = D_np.argmin(axis=0).astype(np.int8)
    # c_t = first ℓ where running-min ≤ γ (1-based, L if never)
    rmin = np.minimum.accumulate(D_np, axis=0)
    below_rmin = rmin <= gamma
    any_rmin = below_rmin.any(axis=0)
    first_rmin = below_rmin.argmax(axis=0).astype(np.int16) + 1
    c_t = np.where(any_rmin, first_rmin, D_np.shape[0]).astype(np.int16)
    return {
        "c_t": c_t, "n_enter": n_enter.astype(np.int8), "n_exit": n_exit.astype(np.int8),
        "oscil": oscil.astype(np.int8), "below_frac": below_frac, "min_D": min_D,
        "argmin_layer": argmin_layer,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True, help="e.g. chr1")
    ap.add_argument("--start-idx", type=int, default=0)
    ap.add_argument("--end-idx", type=int, default=None)
    ap.add_argument("--chunk-rows", type=int, default=500_000)
    args = ap.parse_args()

    OUT_DIR = TDIG / f"wgs/results/{args.chrom}"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    log.info(f"chrom={args.chrom} output_dir={OUT_DIR}")

    log.info("loading hg38 …")
    from pyfaidx import Fasta
    fa = Fasta(str(HG38_FA), as_raw=True, sequence_always_upper=True)
    seq = str(fa[args.chrom][:])
    chrom_len = len(seq)
    log.info(f"  {args.chrom} length: {chrom_len:,}")

    # Window schedule
    starts = np.arange(0, chrom_len - WINDOW + 1, STRIDE)
    log.info(f"  {len(starts):,} windows total (window={WINDOW}, stride={STRIDE})")
    if args.end_idx is None: args.end_idx = len(starts)
    starts_run = starts[args.start_idx:args.end_idx]
    log.info(f"  processing [{args.start_idx}, {args.end_idx}) → {len(starts_run):,} windows")

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

    chunks = []
    rows_win, rows_pos, rows_ct, rows_osc, rows_nen, rows_nex, rows_bf, rows_md, rows_al = [], [], [], [], [], [], [], [], []
    chunk_id = args.start_idx // args.chunk_rows if args.start_idx > 0 else 0
    total_rows = 0
    t0 = time.time()
    n_skip = 0
    n_ok = 0

    for i, wstart in enumerate(starts_run):
        try:
            wend = wstart + WINDOW
            seq_win = seq[wstart:wend]
            if len(seq_win) != WINDOW or seq_win.count("N") > WINDOW * 0.5:
                n_skip += 1
                continue

            tokens = model.tokenizer.tokenize(seq_win)
            input_ids = torch.tensor(tokens, dtype=torch.int64).unsqueeze(0).cuda()
            with torch.no_grad():
                _ = sh(input_ids)
            D = cosine_lens(captured, 32, 0)          # [32, 6000]
            D_central = D[:, CENTRAL_START:CENTRAL_END].cpu().numpy().astype(np.float32)  # [32, 3000]
            stats = crossings_all(D_central, GAMMA_COS)

            P = CENTRAL_END - CENTRAL_START
            win_idx = args.start_idx + i
            genomic_positions = (wstart + np.arange(CENTRAL_START, CENTRAL_END)).astype(np.int64)

            rows_win.append(np.full(P, win_idx, dtype=np.int32))
            rows_pos.append(genomic_positions)
            rows_ct.append(stats["c_t"])
            rows_osc.append(stats["oscil"])
            rows_nen.append(stats["n_enter"])
            rows_nex.append(stats["n_exit"])
            rows_bf.append(stats["below_frac"])
            rows_md.append(stats["min_D"])
            rows_al.append(stats["argmin_layer"])
            n_ok += 1
            total_rows += P

            for k in list(captured.keys()): del captured[k]
            del D
            torch.cuda.empty_cache()

            if (i + 1) % 200 == 0:
                elapsed = time.time() - t0
                rate = (i + 1) / elapsed
                eta = (len(starts_run) - i - 1) / rate / 60
                log.info(f"  win {i+1}/{len(starts_run)} n_ok={n_ok} n_skip={n_skip} "
                         f"rate={rate:.2f}/s ETA={eta:.1f}min")

            if total_rows >= args.chunk_rows:
                df = pd.DataFrame({
                    "window_idx": np.concatenate(rows_win),
                    "pos": np.concatenate(rows_pos),
                    "c_t": np.concatenate(rows_ct),
                    "oscil": np.concatenate(rows_osc),
                    "n_enter": np.concatenate(rows_nen),
                    "n_exit": np.concatenate(rows_nex),
                    "below_frac": np.concatenate(rows_bf),
                    "min_D": np.concatenate(rows_md),
                    "argmin_layer": np.concatenate(rows_al),
                })
                chunk_path = OUT_DIR / f"{args.chrom}_per_position_chunk{chunk_id:04d}.parquet"
                df.to_parquet(chunk_path, index=False)
                log.info(f"  flushed chunk {chunk_id}: {len(df):,} rows → {chunk_path}")
                chunk_id += 1
                rows_win, rows_pos, rows_ct, rows_osc, rows_nen, rows_nex, rows_bf, rows_md, rows_al = [], [], [], [], [], [], [], [], []
                total_rows = 0
        except Exception as e:
            log.error(f"win {wstart}: {type(e).__name__}: {e}")
            n_skip += 1

    for h in handles: h.remove()

    # Flush remaining
    if rows_win:
        df = pd.DataFrame({
            "window_idx": np.concatenate(rows_win),
            "pos": np.concatenate(rows_pos),
            "c_t": np.concatenate(rows_ct),
            "oscil": np.concatenate(rows_osc),
            "n_enter": np.concatenate(rows_nen),
            "n_exit": np.concatenate(rows_nex),
            "below_frac": np.concatenate(rows_bf),
            "min_D": np.concatenate(rows_md),
            "argmin_layer": np.concatenate(rows_al),
        })
        chunk_path = OUT_DIR / f"{args.chrom}_per_position_chunk{chunk_id:04d}.parquet"
        df.to_parquet(chunk_path, index=False)
        log.info(f"  final chunk {chunk_id}: {len(df):,} rows → {chunk_path}")

    log.info(f"DONE {args.chrom}: n_ok={n_ok}, n_skip={n_skip}, total_time={(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
