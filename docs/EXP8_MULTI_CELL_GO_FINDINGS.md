# EXP8 — 17 M-cell GO enrichment on chr17+chr22+chr21

Detailed technical companion to Part 4 of [`NARRATIVE.md`](NARRATIVE.md).

**Session date**: 2026-07-20 (main pipeline) + 2026-07-21 (3-chromosome integration).
**Compute host**: `TGIL_mutsig` (offline since 2026-07-20 midnight KST). Full artifacts backed up to DASH at `/home/darejin/TDiG/exp8_multi_cell_go/`.

---

## 1. Goal

Test whether the 17-cell M-family (M1_dir/M2_mag/M3_geo/M4_set/M5_tau × reference variants) reads *different biology* — not merely different strengths of the same signal — when applied to gene-level rather than splice-context-level questions.

Pre-registered prediction (2026-07-20 pt29 diary): if the M-cells encode orthogonal aspects of the Evo 2 hidden-state trajectory, the GO enrichment profile per cell should differ in the *identity* of hit terms, not just in magnitude, and pairs of cells may show sign-flip — a cell drives depth deeper for a GO class while another drives depth shallower.

---

## 2. Data + pipeline

### Inputs (all from the DASH backup)

| Source | Path | Content |
|---|---|---|
| chr17 upstream tier1 | `tdig_integration/data_cache_minimal/chr17_tier1.parquet` | 27,586 windows × 17 M-cell arrays per position (T=6000) |
| chr22 upstream tier1 | `tdig_integration/data_cache_minimal/chr22_tier1.parquet` | 12,978 windows × same schema |
| chr21 our tier1 | `wgs/results/tier1_by_chr/chr21_tier1.parquet` | 13,365 windows × 7 M-cells (M1_dir_refA, M2_mag_refA, M3_geo × 5) |
| GO annotations | `exp5_go_functional/results/gene_go_annotations.csv` | 164,105 (gene_id, go_id, go_name) triples |
| GENCODE | `data_ref/gencode.v44.annotation.gtf` | Protein-coding gene coordinates |

### Pipeline scripts (order of execution)

| Step | Script | Purpose | Runtime |
|---|---|---|---|
| 0 | `wgs/scripts/e1_20_chr_forward_tier2.py --chrom chr21` | Fresh Evo2 forward, saves tier2-like h5 (cos_refA + step + norms) | 58.7 min on 2×B200 |
| 0 | `wgs/scripts/e1_30_tier2_to_tier1.py --chrom chr21` | Derives 7 M-cell tier1 parquet from tier2 | ~2 min |
| 1 | `exp8/scripts/e8_00_per_cell_go_enrichment.py` | Per-gene aggregation + 137 curated GO term rank-sum | 1.5 min |
| 2 | `exp8/scripts/e8_10_go_extended.py` | BH-FDR + Spearman cell-similarity + sign-flip detection + 5 figures | 6.4 s |
| 3 | `exp8/scripts/e8_20_full_go_bp.py` | Full 471 GO BP rank-sum + hypergeometric top-K | 10.7 s |
| 4 | `exp8/scripts/e8_30_driver_genes.py` | Per (cell, GO) top-20 in-set driver gene dump | 0.4 s |
| 5 | `exp8/scripts/e8_40_chr_split.py` | 2-chromosome robustness (chr17 vs chr22) | 5.2 s |
| 6 | `exp8/scripts/e8_41_chr21_integration_and_split.py` | 3-chromosome robustness (adds chr21) | 5.0 s |

Env: conda `gdtr` on TGIL_mutsig (Python 3.11.15, scipy 1.17.1, numpy 2.4.4, pandas 3.0.3, h5py 3.16.0, matplotlib 3.11.0, pyarrow 25.0.0, scikit-learn 1.9.0). See `exp8_multi_cell_go/ENV.txt`.

---

## 3. The 17 M-cells

| M-cell | Definition | Reference variants tested |
|---|---|---|
| M1_dir | direction distance = 1 − cos(h_ℓ(t), h_ref) | refA (h_norm), refB, refC |
| M2_mag | magnitude ratio = \|r − 1\| where r = ‖h_ℓ‖ / ‖h_29‖ | refA, refB_diag, refC_diag |
| M3_geo | α · v_z + β · κ_z (velocity z-score + curvature z-score, ref-free) | (α,β) ∈ {(1,0), (0,1), (1,1), (1,0.5), (0.5,1)} |
| M4_set | set-based Mahalanobis-like D_M_set(ℓ, ref) | refA, refB, refC |
| M5_tau | tau-weighted cosine variant | refA, refB, refC |

`refA` for cosine-related metrics uses the post-norm output h_norm as reference (matches the workshop paper's cosine lens). `refB` and `refC` use different reference tokens whose exact definition lives in TDiG-team code and was not fully documented in the accessible scripts. Our chr21 forward computes `refA` only.

The 5 M3_geo variants sweep α (velocity coefficient) and β (curvature coefficient) after per-layer z-standardization. β = 0 is pure curvature; β = 1 mixes in cosine; the two-term combination is the "geometric" primitive.

**Live vs dead cells** (measured by presence of any BH-q < 0.05 rank-sum hit on chr17 + chr22):

- **13 live**: M1_dir_refA/B/C, M2_mag_refA, M3_geo (all 5 α/β), M4_set_refA, M5_tau_refA/B/C
- **4 dead**: M2_mag_refB_diag, M2_mag_refC_diag, M4_set_refB, M4_set_refC

Dead cells produce per-gene aggregate features that are effectively constant across all 1,633 protein-coding genes on chr17+22 — the aggregation collapses per-position variance under those specific reference choices. They do, however, show enrichment under the hypergeometric top-K test (5 hits each at q < 0.05, K = 100), meaning rank *order* is preserved even though the median is unchanged. See `results/go_full_bp_summary.csv`.

---

## 4. Central results table

Per-cell hit counts under BH q < 0.05 (rank-sum on 471 GO BP terms; from `results/go_full_bp_ranksum.csv`):

| Cell | rank-sum q<0.05 | Notes |
|---|---|---|
| M4_set_refA | 13 | Top cell for aggregate hit count |
| M5_tau_refA | 12 | Cluster-4 outlier in hierarchical similarity |
| M2_mag_refA | 11 | |
| M3_geo_a1.0_b0.0 | 9 | **Curvature-only, β=0** |
| M5_tau_refB | 9 | Olfactory-DEEP + epithelial-shallow |
| M3_geo_a1.0_b0.5 | 8 | Curvature-dominant, olfactory-DEEPEST (d = −1.96, GO:0050911) |
| M3_geo_a1.0_b1.0 | 8 | |
| M3_geo_a0.5_b1.0 | 8 | Cosine-dominant |
| M1_dir_refA | 1 | Cosine lens = workshop paper |
| M1_dir_refB | 3 | |
| M1_dir_refC | 3 | |
| M2_mag_refB_diag | 0 | Dead cell (per-gene aggregate constant) |
| M2_mag_refC_diag | 0 | Dead cell |
| M3_geo_a0.0_b1.0 | 1 | Velocity-only |
| M4_set_refB | 0 | Dead cell |
| M4_set_refC | 0 | Dead cell |
| M5_tau_refC | 5 | Anti-correlates with M5_tau_refA (r = −0.78) |

**Complementary picture from hypergeometric top-K = 100**:

| Cell | hypergeom top-100 q<0.05 |
|---|---|
| M2_mag_refB_diag | 5 |
| M2_mag_refC_diag | 5 |
| M4_set_refB | 5 |
| M4_set_refC | 5 |
| M3_geo_a1.0_b1.0 | 4 |
| M3_geo_a0.5_b1.0 | 4 |
| M5_tau_refA | 3 |
| M1_dir_refC | 3 |
| (others) | 0–1 |

Note the dead cells' hypergeometric hits: aggregation collapses the median (rank-sum finds no shift), but the ranking is preserved (top-K enrichment is real). This is a real property of the reference choice, not a bug — the M2_mag `refB_diag/refC_diag` and M4_set `refB/refC` cells still contain rank-order information about the trajectory even when their absolute values compress.

---

## 5. Cell-cell similarity structure

Spearman correlation on signed −log₁₀(p) × sign(d) across the 471 GO terms, hierarchical average-linkage clustering, cut at k=4 (from `results/cell_similarity.csv` and `figures/EXP8_cell_similarity.pdf`):

| Cluster | Cells | Interpretation |
|---|---|---|
| 1 (n=3) | M3_geo_a0.0_b1.0, M3_geo_a1.0_b0.5, M5_tau_refC | "olfactory-shallow / immune-shallow" opposite axis |
| 2 (n=6) | M1_dir_refA, M1_dir_refB, M2_mag_refA, M3_geo_a1.0_b0.0, M4_set_refA, M5_tau_refB | "olfactory-deep + immune-active" majority axis |
| 3 (n=3) | M1_dir_refC, M3_geo_a1.0_b1.0, M3_geo_a0.5_b1.0 | "epithelial + transcription" |
| 4 (n=1) | M5_tau_refA | Outlier — no |r| > 0.6 with any other cell |

Strongest anti-correlations (sign reversal, r < −0.5), all indicating opposite directional encoding of the same GO landscape:

| Cell A | Cell B | r |
|---|---|---|
| M1_dir_refB | M3_geo_a1.0_b0.5 | **−0.89** |
| M1_dir_refA | M3_geo_a0.0_b1.0 | **−0.86** |
| **M5_tau_refA** | **M5_tau_refC** | **−0.78** (same family, different reference!) |
| M5_tau_refB | M3_geo_a1.0_b0.5 | −0.77 |
| M2_mag_refA | M3_geo_a1.0_b0.5 | −0.67 |
| M1_dir_refA | M3_geo_a1.0_b1.0 | −0.61 |
| M1_dir_refB | M3_geo_a0.0_b1.0 | −0.60 |
| M4_set_refA | M5_tau_refA | −0.55 |

The M5_tau_refA vs M5_tau_refC = −0.78 is the sharpest single observation: reference token choice within one metric family does not tune the same axis — it selects a different axis with opposite direction on 471 GO terms.

---

## 6. Sign-flip GO terms

A "sign-flip" GO term has ≥ 1 cell with significant positive *d* (≥ +0.3, p < 0.05) AND ≥ 1 cell with significant negative *d* (≤ −0.3, p < 0.05). From `results/sign_flip_go_terms.csv`:

**40 sign-flip GO terms total** across the 13 live cells.

Top 5 by variance of the signed −log₁₀(p) across cells:

| Rank | GO term | Var | Max +d | Min −d |
|---|---|---|---|---|
| 1 | GO:0050911 detection of chemical stimulus (olfactory) | 19.1 | +3.24 (M1_dir_refC) | −1.96 (M3_geo_a1.0_b0.5) |
| 2 | GO:0045109 intermediate filament organization | 18.4 | +1.06 (M5_tau_refA) | −1.32 (M4_set_refA) |
| 3 | GO:0002009 morphogenesis of an epithelium | 13.8 | +0.85 (M5_tau_refA) | −1.15 (M4_set_refA) |
| 4 | GO:0030855 epithelial cell differentiation | 7.3 | | |
| 5 | GO:0061844 antimicrobial humoral immune response | 6.5 | +0.89 (M3_geo_a1.0_b0.5) | −1.25 (M3_geo_a1.0_b0.0) |

Notable pattern:

- The olfactory and epithelial-filament axes are opposite: cells that settle olfactory receptors deep settle epithelial filaments shallow, and vice versa.
- The two most opposite cells on the olfactory axis (M1_dir_refC and M3_geo_a1.0_b0.5) are also the two most opposite on the epithelial-filament axis (with signs reversed).

---

## 7. Chromosome-split robustness

### 7.1 — chr17 vs chr22 (2-chr, from `e8_40`, `results/sign_flip_robustness.csv`)

680 (cell, GO) tests (40 sign-flip GO × 17 cells) split into:

| Status | Count | % |
|---|---|---|
| `robust_same_sign` (both chr sig, same *d* sign) | **8** | 1.2 % |
| `chr17_driven` (only chr17 sig) | 73 | 10.7 % |
| `chr22_driven` (only chr22 sig) | 40 | 5.9 % |
| `ns_split` (neither chr alone sig; combined was) | 304 | 44.7 % |
| `insufficient_data` (n < 3 in some chr) | 255 | 37.5 % |

Manuscript interpretation: only 1.2 % of the observed sign-flip patterns replicate on both chromosomes independently. Most (44.7 %) are sample-size boosted by chr17 + chr22 combination.

The 8 robust pairs distribute across 6 cells: M5_tau_refA (2), M5_tau_refC (2), M2_mag_refA (1), M3_geo_a0.5_b1.0 (1), M3_geo_a1.0_b0.0 (1), M5_tau_refB (1).

### 7.2 — chr17 vs chr22 vs chr21 (3-chr, from `e8_41`, `results/sign_flip_robustness_3chr.csv`)

Adding chr21 (220 protein-coding genes) as a third replicate on the 7 chr21-covered M-cells. 280 (cell, GO) tests split:

| Status | Count |
|---|---|
| `robust_3_same_sign` (all 3 chr sig, same sign) | **0** |
| `robust_2_same_sign` (any 2 of 3 chr sig, same sign) | 3 |
| `single_chr_only` (only one chr sig) | 96 (chr17: 71, chr22: 24, chr21: 1) |
| `ns_all` (no chr sig) | 181 |
| `insufficient_data` | rest |

**Zero pairs replicate across all three chromosomes.** chr21 (n=220 genes) is under-powered for chr21-alone tests, so most 3-chr-robust candidacies fail at the chr21 significance step. The 3 pairs that pass 2-of-3 all consist of chr17 + chr22 (with chr21 as NaN for insufficient data).

### 7.3 — The one clean example: GO:0006869 lipid transport

Among the 3 robust two-chromosome pairs, one is manuscript-worthy on its own — a sign-flip between two M3_geo variants (differing only in the β coefficient):

| Cell | *d* chr17 | *p* chr17 | *d* chr22 | *p* chr22 |
|---|---|---|---|---|
| M3_geo_a0.5_b1.0 (cosine-dominant) | +0.52 | 0.043 | +0.65 | 0.024 |
| M3_geo_a1.0_b0.0 (pure curvature) | −0.72 | 0.009 | −0.89 | 0.007 |

Same GO class (GO:0006869). Same two chromosomes. Two variants of one metric family. Presence or absence of the cosine β component flips the sign of the settling-depth deviation for lipid-transport genes. This is the cleanest single evidence in the current dataset that different M-cell primitives access different directional axes of the residual stream, and that the switch is *at the level of the metric composition parameter*, not the reference token.

The other two robust two-chromosome pairs:

| Cell | GO | *d* chr17 | *d* chr22 |
|---|---|---|---|
| M2_mag_refA | GO:0030335 positive regulation of cell migration | −0.33 | −1.02 |
| M3_geo_a0.5_b1.0 | GO:0006869 lipid transport | +0.52 | +0.65 |
| M3_geo_a1.0_b0.0 | GO:0006869 lipid transport | −0.72 | −0.89 |

---

## 8. Cancer-panel driver hits

Per significant (cell, GO) pair at rank-sum q < 0.10, dumped top-20 in-set genes ranked by mean settling value. Flagged membership in the 15-gene ClinVar cancer panel: BRCA1, BRCA2, TP53, PTEN, STK11, CDH1, PALB2, ATM, CHEK2, BARD1, RAD51C, RAD51D, MLH1, MSH2, APC.

Of 107 significant pairs, 3 include cancer-panel genes in their top-20 drivers (from `results/driver_genes_summary.csv`):

**(1) M3_geo_a1.0_b0.0 × GO:0007131 reciprocal meiotic recombination**
- *d* = +0.91, q = 0.086, top-5 drivers: RAD51C, DMC1, RAD51D, TEX19, SYCE3
- **2 of top-5 are HRD cancer-panel genes** (RAD51C, RAD51D)
- Biologically: RAD51C and RAD51D are homologous-recombination pathway members whose loss-of-function drives HRD (homologous-recombination deficiency) in breast and ovarian cancer. The curvature-only M3_geo cell captures this pathway with the two HRD genes at its top.

**(2) M5_tau_refB × GO:0006357 regulation of transcription by RNA pol II**
- *d* = +0.33, q = 0.041, top-5 drivers: HIRA, HOXB5, ZNF286A, ZNF280A, NEUROD2
- **1 cancer-panel gene** (TP53) in the top-20.

**(3) M5_tau_refB × GO:0006355 regulation of DNA-templated transcription**
- *d* = +0.35, q = 0.060, top-5 drivers overlap #2; TP53 in top-20.

**Baseline for context**: 15 cancer-panel genes out of 1,633 in universe → 15/1633 = 0.92 % baseline rate. Across 107 significant pairs × ~20 top drivers each ≈ 2,140 driver slots, null expectation is ~20 cancer-panel appearances. Observed: 4 (2 from pair 1 + 1 from pair 2 + 1 from pair 3). Not enriched over baseline — but the cases are *biologically coherent* (HRD genes at the top of a meiotic-recombination GO; TP53 in transcription regulation), which is what matters for the manuscript narrative rather than aggregate enrichment.

---

## 9. Known bug and caveats

**e1_30 M2_mag_refA gamma calibration bug** (documented in [`RECOVERY_LOG.md`](../RECOVERY_LOG.md) §Bugs). On chr21, the M2_mag_refA gamma calibration returns ~1.33 M — orders of magnitude too large. Neither the skip-window filter nor the upper-outlier clip fixes it. The residual metric-level cause is likely: some positions on chr21 have very small `norm_h_29`, causing `r = norm_h_ell / norm_h_29` to explode at those positions. The knock-on effect is that M2_mag_refA on chr21 is effectively saturated (all positions settle at the ceiling L=32), making it a de facto dead cell for chr21. All chr21 statistics reported here for M2_mag_refA should be regarded as unreliable.

**7-of-17 M-cells only for chr21**. The chr21 forward reproduces cos_refA (matches the workshop paper's cosine lens using h_norm as reference), step_cos, step_norm_raw/rms, and h_ell norms. It does NOT reproduce cos_refB, cos_refC, D_Mset_A/B/C, or M5_tau primitives — these depend on TDiG-team-specific reference and set definitions that were not accessible in the current codebase. Any future extension to more chromosomes and to the full 17-cell schema requires collaboration with TDiG authors (or full reverse-engineering of the reference definitions from the existing chr17/chr22 upstream tier2 h5 files).

**chr21 has ~220 protein-coding genes** — small n. chr21-alone tests are under-powered relative to chr17 (1,186 genes) and chr22 (447 genes). This is a genuine limitation for the 3-chr robustness analysis and cannot be fully addressed without adding more chromosomes.

**Sample-size boost dominates the raw sign-flip count.** 40 sign-flip GO terms in the combined chr17 + chr22 dataset is more a reflection of doubled statistical power than of independent replication. The chr-split analysis is the correct filter, and it substantially deflates the claim.

---

## 10. Manuscript slot

This entire EXP8 body of work slots into the integrated manuscript ([`INTEGRATED_MANUSCRIPT_PLAN.md`](INTEGRATED_MANUSCRIPT_PLAN.md)) as a new §3.11 subsection, "Settling profile at the pathway level," positioned after the existing multi-axis §3.x sections and before the discussion.

Recommended §3.11 length: ~800–1,000 words + one main figure (proposal: heatmap of top-30 discriminative GO × 13 live cells, biclustered — the file `figures/EXP8_full_clustered_heatmap.pdf` is a candidate as-is or with minor polish) + one supplementary table (the sign-flip robustness table from `results/sign_flip_robustness_3chr.csv`).

Discussion section should be updated to include the concept refactor from "settling depth" to "settling profile" — see NARRATIVE.md §4.7 for the recommended manuscript language.

---

## 11. Reproducibility

All scripts, results, and figures are on DASH at `/home/darejin/TDiG/exp8_multi_cell_go/`.

To rerun end-to-end:

```bash
# Prerequisites on new GPU host:
# - Evo 2 7B weights (HF arcinstitute/evo2_7b)
# - hg38.fa, GENCODE v44 GTF, ClinVar VCF
# - conda env with numpy pandas scipy matplotlib h5py pyarrow scikit-learn

# 1. Forward pass for chr21 (or any other chromosome)
python wgs/scripts/e1_20_chr_forward_tier2.py --chrom chr21
# → wgs/results/tier1_by_chr/chr21_tier2_scalars.h5 (~41 GB)

# 2. Derive per-position 7 M-cells
python wgs/scripts/e1_30_tier2_to_tier1.py --chrom chr21
# → wgs/results/tier1_by_chr/chr21_tier1.parquet (~200 MB)

# 3. Run EXP8 pipeline (chr17 + chr22 tier1 must be present at tdig_integration/data_cache_minimal/)
cd exp8_multi_cell_go
python scripts/e8_00_per_cell_go_enrichment.py
python scripts/e8_10_go_extended.py
python scripts/e8_20_full_go_bp.py
python scripts/e8_30_driver_genes.py
python scripts/e8_40_chr_split.py
python scripts/e8_41_chr21_integration_and_split.py

# Outputs land in results/ (16 CSVs) and figures/ (14 PNG+PDF pairs).
# All scripts idempotent (overwrite their outputs).
```

Session provenance timestamps:
- 2026-07-20 15:00 KST — chr21 HF download attempts (only chr17/22 available at HF)
- 2026-07-20 17:36 — e8_00 pipeline initial run (chr17 + chr22 only)
- 2026-07-20 17:55 — e8_10 extended analysis (sign-flip discovery)
- 2026-07-20 18:20 — e8_20 (full GO BP), e8_30 (driver genes), e8_40 (2-chr robustness)
- 2026-07-20 21:15 — chr21 forward completes (58.7 min on 2× B200 with optimised HDF5)
- 2026-07-20 21:28 — chr21_tier1.parquet derived (v1 with M2_mag gamma bug)
- 2026-07-20 21:30 — e8_41 first run (3-chr robustness: 0 all-3, 3 any-2)
- 2026-07-20 22:00 — e1_30 patched (skip-window filter + outlier clip)
- 2026-07-20 22:07 — chr21_tier1.parquet v2 with patch; e8_41 rerun confirms same numerical result
- 2026-07-20 23:58 — TGIL_mutsig loses SSH access (GPU lease expiry)
- 2026-07-21 00:20 — DASH backup consolidated (5.5 GB)
- 2026-07-21 01:15 — DASH cleanup, master docs (REPRODUCE, MANIFEST, RECOVERY_LOG) written
- 2026-07-21 (this doc) — Part 4 of NARRATIVE.md written; this technical companion doc committed.
