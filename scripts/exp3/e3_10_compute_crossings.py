"""EXP3 §3.3 — Compute enter/exit/oscil per position from raw D_cos trajectory.

Reads: results/A_windows_chr22.npz
Writes:
  results/A_crossing_stats.parquet   [n_positions, columns:
    win_idx, pos, label,
    c_t (running-min settling depth, same as paper c(t)),
    n_enter, n_exit, oscil, below_frac, min_D, argmin_layer,
    first_enter_layer, last_exit_layer
  ]

Definitions (paper's γ_cos = 0.397 threshold):
- n_enter(t) = # of layer transitions where D_cos(ℓ-1) > γ AND D_cos(ℓ) ≤ γ
- n_exit(t)  = # of layer transitions where D_cos(ℓ-1) ≤ γ AND D_cos(ℓ) > γ
- oscil(t)   = max(0, n_enter + n_exit - 1)  # 0 = clean single crossing; higher = oscillating
- below_frac(t) = mean_ℓ [D_cos(ℓ, t) ≤ γ]  # what fraction of layers are below threshold
- c_t(t) = first ℓ where running-min(D_cos)(ℓ, t) ≤ γ  (paper's settling depth)
"""
from __future__ import annotations
import time
from pathlib import Path

import numpy as np
import pandas as pd

TDIG = Path("/NHNHOME/WORKSPACE/0526040123_A/darejinn/tdig")
NPZ = TDIG / "exp3_threshold_crossing/results/A_windows_chr22.npz"
OUT_PARQUET = TDIG / "exp3_threshold_crossing/results/A_crossing_stats.parquet"

GAMMA = 0.397
N_LAYERS = 32


def main() -> None:
    t0 = time.time()
    data = np.load(NPZ)
    D_cos = data["D_cos"].astype(np.float32)          # [N, 32, T=3000]
    positions = data["genomic_positions"].astype(np.int64)  # [N, T]
    labels = data["labels"].astype(np.uint8)          # [N, T]
    windows_start = data["windows_start"]
    N, L, T = D_cos.shape
    print(f"loaded: N={N} windows, L={L} layers, T={T} positions per window")

    below = D_cos <= GAMMA                # [N, L, T] bool

    # Enter: D_cos crossed FROM above → below.
    # For layer ℓ (ℓ in 1..L-1): below[ℓ] and NOT below[ℓ-1].
    enter_mask = below[:, 1:, :] & (~below[:, :-1, :])  # [N, L-1, T]
    exit_mask  = below[:, :-1, :] & (~below[:, 1:, :])  # [N, L-1, T]
    # Layer 0 "enter" if D_cos(0) ≤ γ from the very start
    layer0_enter = below[:, 0:1, :]                     # [N, 1, T]

    n_enter = enter_mask.sum(axis=1) + layer0_enter.astype(np.int32).sum(axis=1)  # [N, T]
    n_exit  = exit_mask.sum(axis=1)                                                # [N, T]

    # oscil = # of full "enter → exit → re-enter" cycles roughly
    # = max(0, n_enter + n_exit - 1). A clean single crossing gives n_enter=1, n_exit=0 → oscil=0.
    # If crossed and re-crossed (dipped below then went back above then below again):
    #   n_enter=2, n_exit=1 → oscil=2 (2 "extra" crossings beyond the minimal single dip).
    oscil = np.clip(n_enter + n_exit - 1, 0, None)     # [N, T]

    # below_frac
    below_frac = below.mean(axis=1).astype(np.float32)  # [N, T]

    # min D_cos and its layer
    min_D = D_cos.min(axis=1)                          # [N, T]
    argmin_layer = D_cos.argmin(axis=1).astype(np.int8)  # [N, T]

    # first / last layer where below (0-based; -1 if never)
    below_int = below.astype(np.int8)                  # [N, L, T]
    any_below = below.any(axis=1)                      # [N, T]
    # first: argmax over below_int returns first True or 0 if all False; guard with any_below
    first_below = below_int.argmax(axis=1)             # [N, T]
    first_below = np.where(any_below, first_below, -1).astype(np.int16)
    # last: reverse-argmax
    last_below = L - 1 - below_int[:, ::-1, :].argmax(axis=1)
    last_below = np.where(any_below, last_below, -1).astype(np.int16)

    # Paper c(t): first ℓ where running-min(D_cos) ≤ γ (1-based, L if never)
    rmin = np.minimum.accumulate(D_cos, axis=1)        # [N, L, T]
    below_rmin = rmin <= GAMMA
    any_rmin = below_rmin.any(axis=1)
    first_rmin = below_rmin.argmax(axis=1).astype(np.int16) + 1
    c_t = np.where(any_rmin, first_rmin, L).astype(np.int16)

    # Flatten to per-position rows
    flat_win = np.broadcast_to(np.arange(N).reshape(-1, 1), (N, T)).ravel()
    flat_pos = positions.ravel()
    flat_lab = labels.ravel()

    df = pd.DataFrame({
        "win_idx": flat_win.astype(np.int32),
        "pos": flat_pos.astype(np.int64),
        "label": flat_lab.astype(np.uint8),
        "c_t": c_t.ravel(),
        "n_enter": n_enter.ravel().astype(np.int8),
        "n_exit":  n_exit.ravel().astype(np.int8),
        "oscil":   oscil.ravel().astype(np.int8),
        "below_frac": below_frac.ravel(),
        "min_D": min_D.ravel(),
        "argmin_layer": argmin_layer.ravel().astype(np.int8),
        "first_below": first_below.ravel(),
        "last_below":  last_below.ravel(),
    })
    print(f"per-position table: {df.shape}")

    OUT_PARQUET.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(OUT_PARQUET, index=False)
    print(f"wrote {OUT_PARQUET}")

    # Quick summary — just to eyeball
    print("\n=== overall summary ===")
    print(f"positions total: {len(df):,}")
    print(f"positions with any layer below γ: {(df['first_below'] >= 0).sum():,}  ({(df['first_below']>=0).mean()*100:.1f}%)")
    print(f"mean c_t: {df['c_t'].mean():.3f}")
    print(f"mean n_enter: {df['n_enter'].mean():.3f}")
    print(f"mean n_exit:  {df['n_exit'].mean():.3f}")
    print(f"mean oscil:   {df['oscil'].mean():.3f}")
    print(f"oscil distribution: {df['oscil'].value_counts().sort_index().to_dict()}")

    print(f"\nRuntime: {time.time()-t0:.2f} s")


if __name__ == "__main__":
    main()
