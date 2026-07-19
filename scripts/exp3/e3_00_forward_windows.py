"""EXP3 §3.1 — Forward pass over N random chr22 6-kb windows, saving RAW per-layer D_cos trajectory.

This is Option A from exp3/plan.md. Small subset (100 windows × 3 kb central = 300k positions)
to test threshold-crossing hypothesis H3a before committing to variant-scale (Option B) or genome-wide.

Outputs: results/A_windows_chr22.npz with:
  - windows_idx      [N] window index into chr22_windows.tsv
  - genomic_positions [N, 3000] 0-based genome positions
  - D_cos            [N, 32, 3000] float16 per-layer cosine distance to h_norm
  - labels           [N, 3000] uint8 GENCODE annotation code
  - runtime_sec

Also aggregates to a per-position table for the crossing analysis in exp3/scripts/e3_10.
"""
from __future__ import annotations
import os, sys, time, json, argparse, logging
from pathlib import Path

os.environ.setdefault("HF_HOME", "/NHNHOME/WORKSPACE/0526040123_A/darejinn/tdig/env/hf_cache")

import numpy as np
import torch
import torch.nn.functional as F

TDIG = Path("/NHNHOME/WORKSPACE/0526040123_A/darejinn/tdig")
REF_FA = TDIG / "wgs/data/reference/hg38.fa"
CHR22_LABELS = TDIG / "data_ref/chr22_position_labels.npy"
OUT_NPZ = TDIG / "exp3_threshold_crossing/results/A_windows_chr22.npz"
LOG_DIR = TDIG / "exp3_threshold_crossing/logs"

# Paper conventions — frozen
GAMMA_COS = 0.397
WINDOW = 6000
CENTRAL_START = 1500
CENTRAL_END = 4500
N_LAYERS = 32
BOS_OFFSET = 0

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("e3_forward")


def cosine_lens(hidden_states, n_layers=N_LAYERS, bos_offset=BOS_OFFSET):
    """D_cos(ell, t) = 1 - cos(h_ell(t), h_norm(t)). Returns float32 [n_layers, T_real].

    Evo 2 on 2× B200 auto-shards blocks across GPUs (0..~15 on cuda:0, ~16..31 on cuda:1).
    We move each block's hidden state to the norm-module's device before the cosine dot product.
    """
    # `norm` output lives on whichever device sh.norm was assigned to (typically cuda:1).
    ref_device = hidden_states["norm"].device
    h_norm = hidden_states["norm"].float()[:, bos_offset:, :].to(ref_device)
    h_norm_n = F.normalize(h_norm, p=2, dim=-1)
    B, T, H = h_norm.shape
    D = torch.zeros((n_layers, T), dtype=torch.float32, device=ref_device)
    for ell in range(n_layers):
        h_l = hidden_states[f"blocks.{ell}"].float()[:, bos_offset:, :].to(ref_device)
        h_l_n = F.normalize(h_l, p=2, dim=-1)
        cos = (h_l_n * h_norm_n).sum(dim=-1)
        D[ell] = (1.0 - cos).clamp(min=0.0).mean(dim=0)
    return D


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-windows", type=int, default=100)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--limit-positions", type=int, default=None,
                    help="(debug) analyze fewer positions per window")
    args = ap.parse_args()

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log.info(f"tdig root: {TDIG}")
    log.info(f"gamma_cos={GAMMA_COS} window={WINDOW} central=[{CENTRAL_START},{CENTRAL_END})")

    # Load chr22 labels for annotation coverage
    log.info("loading chr22 position labels ...")
    labels_full = np.load(CHR22_LABELS)
    log.info(f"  labels shape: {labels_full.shape}, unique: {dict(zip(*np.unique(labels_full, return_counts=True)))}")

    # Load reference
    log.info("loading hg38 chr22 ...")
    from pyfaidx import Fasta
    fa = Fasta(str(REF_FA), as_raw=True, sequence_always_upper=True)
    seq22 = str(fa["chr22"][:])
    log.info(f"  chr22 length: {len(seq22)}")
    assert labels_full.shape[0] == len(seq22), f"labels vs sequence length mismatch: {labels_full.shape[0]} vs {len(seq22)}"

    # Load Evo 2
    log.info("loading Evo 2 7B ...")
    t_load = time.time()
    from evo2 import Evo2
    model = Evo2("evo2_7b")
    sh = model.model
    log.info(f"  Evo 2 loaded in {time.time()-t_load:.1f}s")

    # Register hooks
    captured = {}
    handles = []
    def make_hook(name):
        def hook(module, inp, out):
            h = out[0] if isinstance(out, tuple) else out
            captured[name] = h.detach()
        return hook
    for i in range(N_LAYERS):
        handles.append(sh.blocks[i].register_forward_hook(make_hook(f"blocks.{i}")))
    handles.append(sh.norm.register_forward_hook(make_hook("norm")))

    # Sample N windows — start positions uniformly in valid range, seeded
    rng = np.random.default_rng(args.seed)
    max_start = len(seq22) - WINDOW
    starts = np.sort(rng.integers(0, max_start, size=args.n_windows))
    log.info(f"sampled {args.n_windows} windows, starts range [{starts[0]}, {starts[-1]}]")

    N = args.n_windows
    P = CENTRAL_END - CENTRAL_START  # 3000
    all_pos = np.zeros((N, P), dtype=np.int32)
    all_lab = np.zeros((N, P), dtype=np.uint8)
    all_Dcos = np.zeros((N, N_LAYERS, P), dtype=np.float16)

    n_ok = 0
    t_fwd_total = 0.0
    for i, wstart in enumerate(starts):
        wend = wstart + WINDOW
        if wend > len(seq22):
            continue
        seq_win = seq22[wstart:wend]
        # skip windows with lots of Ns
        if seq_win.count("N") > WINDOW * 0.1:
            log.info(f"  win {i} @ {wstart}: skipped ({seq_win.count('N')} Ns)")
            continue

        # Tokenize
        try:
            tokens = model.tokenizer.tokenize(seq_win)
            input_ids = torch.tensor(tokens, dtype=torch.int64).unsqueeze(0).cuda()

            t0 = time.time()
            with torch.no_grad():
                _ = sh(input_ids)
            D = cosine_lens(captured, N_LAYERS, BOS_OFFSET)  # [L, T] fp32 cuda
            t_fwd = time.time() - t0
            t_fwd_total += t_fwd

            # Slice central 3 kb positions [CENTRAL_START, CENTRAL_END)
            D_central = D[:, CENTRAL_START:CENTRAL_END].cpu().numpy().astype(np.float16)
            positions_in_window = np.arange(CENTRAL_START, CENTRAL_END)
            genomic_positions = (wstart + positions_in_window).astype(np.int32)
            label_central = labels_full[genomic_positions].astype(np.uint8)

            all_pos[n_ok] = genomic_positions
            all_lab[n_ok] = label_central
            all_Dcos[n_ok] = D_central
            n_ok += 1

            if (i + 1) % 10 == 0 or i == N - 1:
                log.info(f"  win {i+1}/{N} start={wstart} fwd={t_fwd:.2f}s avg={t_fwd_total/(i+1):.2f}s "
                         f"n_ok={n_ok} label_counts={dict(zip(*np.unique(label_central, return_counts=True)))}")

            # Clear GPU memory
            for k in list(captured.keys()):
                del captured[k]
            del D
            torch.cuda.empty_cache()

        except Exception as e:
            log.error(f"  win {i} @ {wstart} FAILED: {type(e).__name__}: {e}")

    for h in handles:
        h.remove()

    # Trim to actual n_ok
    all_pos = all_pos[:n_ok]
    all_lab = all_lab[:n_ok]
    all_Dcos = all_Dcos[:n_ok]

    OUT_NPZ.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT_NPZ,
                        windows_start=starts[:n_ok].astype(np.int64),
                        genomic_positions=all_pos,
                        labels=all_lab,
                        D_cos=all_Dcos,
                        gamma_cos=np.float32(GAMMA_COS))
    log.info(f"wrote {OUT_NPZ}: n_ok={n_ok} windows, D_cos shape={all_Dcos.shape}, total {all_Dcos.nbytes/1e6:.1f} MB")
    log.info(f"total forward time: {t_fwd_total:.1f}s (avg {t_fwd_total/max(1,n_ok):.2f}s per window)")

    manifest = {
        "n_windows_ok": int(n_ok),
        "gamma_cos": float(GAMMA_COS),
        "seed": int(args.seed),
        "avg_forward_sec_per_window": float(t_fwd_total / max(1, n_ok)),
        "produced_by": "exp3_threshold_crossing/scripts/e3_00_forward_windows.py",
        "produced_at": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    (OUT_NPZ.parent / "A_windows_chr22_manifest.json").write_text(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
