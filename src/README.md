# src/ — gDTR framework (frozen from paper commit)

These 14 modules are the **frozen gDTR framework** used by all follow-up analyses in this repo. They are byte-identical (except `__init__.py`) to the corresponding files in the original paper repository [`darejinn/gDTR/src/`](https://github.com/darejinn/gDTR/tree/main/src) at the time of paper acceptance.

Nothing in this directory should be modified. If you need to extend gDTR, add a new module to `scripts/` instead.

## Module summary

| File | Purpose |
|---|---|
| `gdtr.py` | Running-min settling depth. `settling_depth_discrete`, `settling_depth_interp`, `gdtr`, `deep_thinking_mask`. |
| `ur_gdtr_evo2.py` | UR-cosine lens for Evo 2. `cosine_lens(hidden_states)` → D_cos [L, T]. |
| `variant_delta.py` | Variant Δ metrics. `compute_delta_metrics(D_ref, D_alt)` → dict of ΔD, Δc, max_abs_dD, etc. |
| `calibration.py` | Regional q70 calibration. `compute_q70_per_region`. |
| `model_loader_evo2.py` | Evo 2 model loader with block/norm hooks. |
| `constants.py` | Framework-wide constants: `GAMMA_DEFAULT = 0.5` (paper uses 0.397 explicit), `L_DEFAULT`, `RHO_DEFAULT`. |
| `constants_evo2.py` | Evo 2-specific: `N_LAYERS = 32`, `MODEL_NAME`, `BOS_OFFSET = 0`. |
| `controls.py` | Motif edit + flank shuffle helpers (used by paper §3.2). |
| `block_type.py` | StripedHyena-2 block type labels (hcl / hcm / hcs / attn). |
| `tuned_lens.py` | 4096×4096 affine A_ℓ fit (paper App A.4). |
| `logit_lens_evo2.py` | JSD lens (paper legacy). |
| `stats.py` | Mann-Whitney U with effect size, rank-biserial r, etc. |
| `viz.py` | Shared plotting utilities. |
| `ur_gdtr.py` | UR-cosine base class (arch-agnostic). |

## Usage

All scripts in this repo import via:

```python
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from gdtr import settling_depth_discrete, running_min
from ur_gdtr_evo2 import cosine_lens
from variant_delta import compute_delta_metrics
from calibration import compute_q70_per_region
```

## Version pin

If regenerating from scratch, use torch 2.7.0+cu128 (B200) or torch 2.4.1+cu124 (H200 / paper reference). See `data/requirements_phase1.lock.txt` for the paper's exact frozen environment.
