# METHODS — Definitions and provenance for every metric

Every feature you see in `results/` is defined here. Whenever a result file references a column name (e.g. `oscil`, `d_oscil`, `mean_amplitude`), find its definition in this document.

## Notation

| Symbol | Meaning |
|---|---|
| $h_\ell(t)$ | residual-stream state at layer $\ell$, token position $t$ |
| $h_{\mathrm{norm}}(t)$ | post-final-RMSNorm state (Evo 2 output-ready frame) |
| $D_{\cos}(\ell, t)$ | $1 - \cos(h_\ell(t), h_{\mathrm{norm}}(t))$ — cosine distance to output frame |
| $\gamma_{\cos}$ | **0.397** (paper-locked, chr22 penultimate q70) |
| $L$ | 32 (Evo 2 7B block count) |
| $L^\star$ | 29 (paper's "canonical tap" — deepest pre-rotation layer) |

## Core paper feature — `c(t)`

**Definition** (paper §2, [`src/gdtr.py`](../src/gdtr.py)):

$$c(t) = \min\{\ell : \mathrm{run\text{-}min}\, D_{\cos}(\ell, t) \le \gamma_{\cos}\}$$

where $\mathrm{run\text{-}min}\, D_{\cos}(\ell, t) = \min_{k \le \ell} D_{\cos}(k, t)$.

`c(t) = 32` means the token never crossed γ.

**Where computed**: [`scripts/wgs/e1_00_chr_forward.py`](../scripts/wgs/e1_00_chr_forward.py) → column `c_t` in `chr{N}_per_position_chunk*.parquet`.

## Raw crossing features (EXP3, H3a base)

Computed from the [L, T] raw D_cos tensor per window:

### `n_enter`

Number of layer transitions where D_cos crosses from *above* γ to *below*.

$$n_{\mathrm{enter}}(t) = \mathbb{1}[D_{\cos}(0,t) \le \gamma] + \sum_{\ell=1}^{L-1} \mathbb{1}[D_{\cos}(\ell-1,t) > \gamma \wedge D_{\cos}(\ell,t) \le \gamma]$$

### `n_exit`

Number of layer transitions where D_cos crosses from below to above.

$$n_{\mathrm{exit}}(t) = \sum_{\ell=1}^{L-1} \mathbb{1}[D_{\cos}(\ell-1,t) \le \gamma \wedge D_{\cos}(\ell,t) > \gamma]$$

### `oscil` — the key new axis

Number of "extra" crossings beyond a single dip.

$$\mathrm{oscil}(t) = \max(0, n_{\mathrm{enter}} + n_{\mathrm{exit}} - 1)$$

- `oscil = 0` = at most one clean crossing (or none)
- `oscil = 1, 2` are rare (bimodal distribution — see H3a)
- `oscil ≥ 3` = actively re-crossing multiple times

### `below_frac`

Fraction of layers below γ.

$$\mathrm{below\_frac}(t) = \frac{1}{L} \sum_{\ell=0}^{L-1} \mathbb{1}[D_{\cos}(\ell,t) \le \gamma]$$

### `min_D`

$$\mathrm{min\_D}(t) = \min_\ell D_{\cos}(\ell,t)$$

### `argmin_layer`

$$\mathrm{argmin\_layer}(t) = \arg\min_\ell D_{\cos}(\ell,t)$$

**Where computed**: [`scripts/exp3/e3_10_compute_crossings.py`](../scripts/exp3/e3_10_compute_crossings.py) → `A_crossing_stats.parquet` for chr22; [`scripts/wgs/e1_00_chr_forward.py`](../scripts/wgs/e1_00_chr_forward.py) for all 24 chr per-position parquets.

## H3d advanced features (9 more)

Computed in [`scripts/exp3/e3_60_advanced_crossings.py`](../scripts/exp3/e3_60_advanced_crossings.py) from raw D_cos tensor:

### `first_enter_layer`

Layer of first above→below crossing (−1 if never).

### `last_exit_layer`

Layer at which token exited below-γ region (−1 if never leaves).

### `longest_below_streak`

Length of the longest consecutive run of layers below γ.

### `streak_start_layer`

Layer at which the longest streak begins.

### `amplitude_below_gamma`

$$\mathrm{amplitude}(t) = \max(0, \gamma - \min_\ell D_{\cos}(\ell,t))$$

How far below γ the trajectory dipped.

### `early_below_frac`, `mid_below_frac`, `late_below_frac`

Fraction of layers below γ in three layer bands:
- early: layers 0–9
- mid: layers 10–21
- late: layers 22–31

The `late_below_frac` is the most discriminating single feature at splice sites (d = +3.18 in H3d chr22 test).

## Committed / dipped / deliberating strata (H3c)

For positions that crossed γ at least once:

| Stratum | Definition |
|---|---|
| committed | `n_enter = 1 AND n_exit = 0` (clean single crossing that stays) |
| dipped | `n_enter = 1 AND n_exit = 1` (crosses in, then out) |
| deliberating | `oscil ≥ 1` (multi-crossings) |

## Variant features (EXP2 base)

The paper's 32-d `ΔD_cos` vector: for a variant at position $p$ with ref sequence and alt sequence,

$$\Delta D_{\cos}(\ell, p) = D_{\cos}^{\mathrm{alt}}(\ell, p) - D_{\cos}^{\mathrm{ref}}(\ell, p)$$

Additional scalars computed by paper:
- `max_abs_dD` = $\max_\ell |\Delta D_{\cos}(\ell, p)|$
- `argmax_layer` = $\arg\max_\ell |\Delta D_{\cos}(\ell, p)|$ (1-based)
- `signed_argmax` = $\Delta D_{\cos}(\mathrm{argmax\_layer}, p)$

## H3b variant re-forward features

For each variant, forward both ref and alt through Evo 2, extract per-layer D_cos, compute:

- `oscil_ref`, `oscil_alt`, `d_oscil = oscil_alt - oscil_ref`
- `n_enter_ref`, `n_enter_alt`, `d_n_enter`
- `n_exit_ref`, `n_exit_alt`, `d_n_exit`
- `below_frac_ref`, `below_frac_alt`, `d_below_frac`
- `c_t_ref`, `c_t_alt`, `d_c_t`
- `min_D_ref`, `min_D_alt`, `d_min_D`
- `argmin_layer_ref`, `argmin_layer_alt`

**Where computed**: [`scripts/exp3/e3_40_h3b_variant_forward.py`](../scripts/exp3/e3_40_h3b_variant_forward.py) → `B_variants_oscil.parquet` (regenerable).

## STEP 3 WGS-normalized variant features

For each variant with context $c$ (inferred from consequence: e.g. missense/nonsense/synonymous → `coding_exon` background, canonical_splice → `splice_donor` background, intron → `intron` background):

$$z_{f, \mathrm{ref}}(v) = \frac{f_{\mathrm{ref}}(v) - \mu_{c, f}^{\mathrm{WGS}}}{\sigma_{c, f}^{\mathrm{WGS}}}$$

where $\mu_{c, f}^{\mathrm{WGS}}$ and $\sigma_{c, f}^{\mathrm{WGS}}$ are computed on a **200 000-position random sample per chromosome** from the 24-chr per-position parquets, filtered by context.

Same for `z_{f, alt}`; then `z_d_f = z_{f, alt} - z_{f, ref}`.

Background stored at [`results/genome_summary/wgs_variant_background.json`](../results/genome_summary/wgs_variant_background.json).

## Cohen's d convention

$$d = \frac{\bar x_{\mathrm{ctx}} - \bar x_{\mathrm{intron}}}{\sigma_{\mathrm{pooled}}}$$

where $\sigma_{\mathrm{pooled}} = \sqrt{\frac{(n_1-1)\sigma_1^2 + (n_2-1)\sigma_2^2}{n_1 + n_2 - 2}}$.

Convention: **negative d = shallower than intron** (settles earlier); **positive d = deeper** (settles later).

## MC → variant subtype class map (paper Fig 3)

Priority-ordered:

```python
{
  "SO:0001583|missense_variant":     "missense",
  "SO:0001587|nonsense":              "nonsense",
  "SO:0001587|stop_gained":           "nonsense",
  "SO:0001589|frameshift_variant":    "frameshift",
  "SO:0001819|synonymous_variant":    "synonymous",
  "SO:0001629|splice_acceptor_variant": "canonical_splice",
  "SO:0001575|splice_donor_variant":  "canonical_splice",
  "SO:0001627|intron_variant":        "intron",
}
```

## Reproducibility parameters

Every stochastic step uses `SEED = 42`. Every LR pipeline uses `StandardScaler + LogisticRegression(max_iter=2000, random_state=42)`. Every CV uses `StratifiedKFold(n_splits=10, shuffle=True, random_state=42)`.

## Data provenance table

| Source | Version | File |
|---|---|---|
| Evo 2 7B weights | HF revision `bda0089f92582d5baabf0f22d9fc85f3588f6b58` | see `data/MODEL_REVISIONS.txt` |
| GRCh38 | UCSC hg38 | `hg38.fa` (3 GB, regenerable) |
| GENCODE annotation | v44 basic | `gencode.v44.annotation.gtf` (1.5 GB) |
| ClinVar variant classifications | 2026-04-18 archive (paper) / current weekly release (this repo) | `clinvar.vcf.gz` (184 MB) |
| ENCODE cCRE | SCREEN v3 (paper) / UCSC encodeCcreCombined.bb (this repo) | `encodeCcreCombined.bed` (145 MB) |
| GTEx | v8 cis-eQTL | not used here |
| GWAS | Catalog v1.0 | not used here |
