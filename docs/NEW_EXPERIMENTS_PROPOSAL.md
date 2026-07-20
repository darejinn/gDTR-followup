# Proposed new experiments — metric interpretation + biological utility

Two new experiment families to add to the integrated manuscript, addressing the two most common reviewer concerns for interpretability-lens papers.

**Important correction on this revision**: an earlier draft proposed a "Why cosine?" ablation. But TDiG has already run a much more thorough metric comparison (5 metric families × 3 reference variants = 17 cells, chr22 with 1,000-window bootstrap CI). So this document has been rewritten around the *interpretation* of TDiG's existing metric-comparison result — a re-framing question — rather than adding a redundant ablation. The old framing is preserved in [`CHANGELOG.md`](CHANGELOG.md).

Reference documents:
- [`INTEGRATED_MANUSCRIPT_PLAN.md`](INTEGRATED_MANUSCRIPT_PLAN.md) — where these experiments land
- [`TDIG_INTEGRATION_STRATEGY.md`](TDIG_INTEGRATION_STRATEGY.md) — TDiG's contribution list (T1–T8)
- TDiG upstream: [`results/RESULTS_v3.md`](https://github.com/heoneyzi/TDiG/blob/main/results/RESULTS_v3.md), specifically §E bootstrap CI on 17 cells

---

## Experiment A — Interpreting TDiG's metric comparison for the workshop paper

### A.0 The already-existing result

TDiG ran the systematic metric comparison the workshop paper never did. On chr22 with 1,000-window bootstrap CI (`bootstrap_d_ci.csv`, RESULTS_v3.md §E), the splice_donor vs intron Cohen's *d* per metric cell is:

| Cell (metric × reference) | mean d | 95% CI | family |
|---|---|---|---|
| M3_geo_a0.5_b1.0 (curvature-weighted, ref-free) | **−0.801** | [−0.827, −0.775] | geometry |
| M3_geo_a0.0_b1.0 (curvature-only, ref-free) | −0.780 | [−0.810, −0.754] | geometry |
| M3_geo_a1.0_b1.0 (velocity + curvature) | −0.574 | [−0.603, −0.545] | geometry |
| M5_tau_refB (RMSNormed path efficiency) | −0.426 | [−0.447, −0.404] | path efficiency |
| M3_geo_a1.0_b0.5 (velocity-weighted) | −0.335 | [−0.357, −0.313] | geometry |
| **M1_dir_refA (raw cosine, paper baseline)** | **+0.301** | [+0.286, +0.318] | direction |
| M2_mag_refA (magnitude ratio) | +0.178 | [+0.156, +0.198] | magnitude |
| M5_tau_refA (path efficiency, raw ref) | −0.171 | [−0.188, −0.156] | path efficiency |
| M4_set_refA (Mahalanobis, whitened) | −0.158 | [−0.172, −0.147] | distribution |
| M1_dir_refB (both-normed cosine, DTR-style) | +0.148 | [+0.132, +0.164] | direction |
| M3_geo_a1.0_b0.0 (velocity-only, ref-free) | +0.147 | [+0.131, +0.163] | geometry |
| M1_dir_refC (workshop paper's default) | −0.069 | [−0.085, −0.051] | direction |
| M5_tau_refC (path efficiency vs h_norm) | +0.024 | [+0.014, +0.033] | path efficiency |

(TDiG's `M1_dir_refC` is what the workshop paper §2 calls `D_cos`. TDiG's chr22 bootstrap gives d = −0.069, small; but the workshop paper's own §3.1 reports d = −0.354 on the pooled chr17 + chr22 baseline — the difference is baseline definition + sample size, not framework contradiction. Discussed below.)

### A.1 The claim change this forces

- **Curvature-based M3_geo is the strongest single settler discriminator** at splice sites, not the workshop paper's cosine (M1_dir_refC).
- Cosine is not "wrong" — it captures a real signal (paper reproduced on chr17 + chr22 pool), but on the same chr22 sample the curvature axis captures a **~10× stronger** effect.
- Magnitude (M2_mag) is small in the *opposite direction* — matches TDiG's "residual accumulation artifact" acknowledgement.
- Mahalanobis (M4_set) is small — matches the observation that Σ_ref whitening collapses toward identity.

**The integrated manuscript should NOT frame `c(t)` as the primary settling metric.** It should frame it as *one axis* in a metric family where M3_geo turns out to be the strongest axis at splice sites, but c(t) remains valuable for downstream tasks where directionality matters (variant scoring — TDiG's own §G1 says M3_geo underperforms as a scorer despite winning on context discrimination). This "different axes win different tasks" is exactly the task-decomposition claim from [`TDIG_INTEGRATION_STRATEGY.md`](TDIG_INTEGRATION_STRATEGY.md) T8.

### A.2 New experiments (small, targeted — do not duplicate TDiG)

Because TDiG already did the systematic comparison, we do NOT re-run a full metric ablation. Instead we run three targeted analyses that TDiG has not yet done and that make the metric comparison actionable for the manuscript:

**A.2.a — Reconcile the sign discrepancy in cosine.**
TDiG's M1_dir_refC (workshop's cosine) gives d = −0.069 on chr22 bootstrap. Workshop paper §3.1 reports d = −0.354 on chr17 + chr22 pool. Two possible resolutions:
- baseline definition difference (per-window intron vs pooled genome intron)
- sample-size difference (TDiG 1,000 windows vs workshop 12,978 chr22 windows)

Test: recompute M1_dir_refC d on chr22 with (i) workshop paper's exact intron baseline and (ii) TDiG's per-window baseline, over the same 79-window sample where we have raw D_cos in [`../results/exp3/A_windows_chr22.npz`](../results/exp3/A_windows_chr22.npz). Report both values and explain the difference.

Cost: 5 min CPU. Deliverable: `results/exp4/cosine_baseline_reconciliation.json`.

**A.2.b — Add curvature axis to the WGS results.**
STEP 2's `argmin_layer` and `amplitude` cover part of TDiG's M3_geo territory (via crossing dynamics) but not directly velocity + curvature. If TDiG can share their per-position M3_geo on chr17 + chr22, we can:
- Cross-reference their curvature axis with our WGS `oscil` on the same positions
- Test: does M3_geo curvature correlate with our `oscil` at splice sites? Prediction: high curvature (M3_geo strong signal) coincides with high oscil (our finding), giving convergent evidence for "residual stream repeatedly reconsidering."

Cost: TDiG provides `chr22_tier1_settling_v2.parquet` (~883 MB) or its M3_geo column subset; ~10 min CPU to join with our chr22 crossing stats.

**A.2.c — Task-decomposition benchmark.**
Extend the T8 benchmark table from [`TDIG_INTEGRATION_STRATEGY.md`](TDIG_INTEGRATION_STRATEGY.md) into an executed table showing which axis wins which task:

| Task | Best-performing feature family | Runner-up |
|---|---|---|
| Splice-site context discrimination (Cohen d) | M3_geo curvature (d = −0.80) | M5_tau_refB (d = −0.43) |
| Variant pathogenicity binary AUROC | TDiG ΔH_norm 64-d (0.949) | Workshop cos32 (0.844) |
| Variant subtype 4-class macro-F1 | STEP 3 cos32 + WGS-norm H3 (0.815) | cos32 alone (0.641) |
| Consequence-layer split (per-layer AUROC) | TDiG ΔH per layer (0.855 at L=8) | (no comparable single-layer baseline) |
| Intron-outlier discovery (near-splice enrichment) | TDiG M5_tau_refB (2.70×) | our top-oscil intron subset (to compute) |

**The unified claim**: no single axis dominates. Different biological tasks recruit different geometric properties of the residual stream. Cosine wins for variant scoring; curvature wins for context discrimination; magnitude is an artifact; Mahalanobis whitening washes out signal. **§3.6 becomes a task-decomposition table, not a "cosine is best" table.**

Cost: aggregation only, no new forward passes. Existing per-task result files.

### A.3 What this does for the manuscript

- Removes a redundant experiment I initially proposed
- Elevates TDiG's already-done metric bootstrap to a manuscript centerpiece
- Reframes cosine as one member of a family, not as THE method
- Sets up §3.6 as a task-decomposition analysis rather than a scorer race

---

## Experiment B — Gene Ontology (GO) functional coherence

### B.0 Motivation

The workshop paper demonstrates gDTR's utility with a variant classifier — that is *per-position* utility. But biological interpretation happens at the *gene* level. A reviewer will naturally ask:

> "If gDTR features are biologically meaningful, do genes with similar gDTR profiles share similar Gene Ontology functions?"

This experiment answers that. If gDTR features aggregated per-gene are enriched for coherent GO categories, the framework is *functionally coherent* — not just position-level noise averaging out. It also opens a downstream use case: **gene function prediction from residual-stream trajectory alone**, no protein sequence needed.

### B.1 Hypotheses

**H_B1** (necessary condition): Per-gene aggregate gDTR features (mean c(t), mean oscil, mean amplitude, mean argmin_layer, aggregated over the gene's exons/introns) cluster genes into functionally coherent groups. Random shuffling of gene labels destroys this clustering.

**H_B2** (specific predictions from paper context findings):
- **GO:0000398 (mRNA splicing via spliceosome)** — genes annotated here should have shallower splice_donor mean c(t) and higher mean oscil than genome-wide splice_donor mean. Prediction magnitude: Δ ≥ 0.2 layers, Δ_oscil ≥ 0.05, both directions consistent with paper §3.1 + our H3a.
- **GO:0006412 (translation)** — genes here should show the coding_exon "late clean deep-commit" pattern from STEP 2: **mean argmin_layer closer to 27**, **frac_committed higher** than genome-wide coding_exon mean.
- **GO:0006351 (transcription, DNA-templated)** — exploratory; broader signal.

**H_B3** (downstream utility): Adding per-gene GO features to the variant subtype classifier improves macro-F1 modestly beyond STEP 3's 0.815.

### B.2 Data needed

- **WGS 24-chr per-position parquets** (already have — 37 GB on server)
- **GENCODE v44 gene body coordinates** (already have `gencode.v44.annotation.gtf`, 1.5 GB)
- **GO annotations per gene**: download from `mygene.info` API (~50 MB, batch query 20,000 genes)
- **Ensembl BioMart backup** if mygene.info incomplete

### B.3 Method

**Step 1 — per-gene feature aggregation** (~5 min CPU)

For each protein-coding gene in GENCODE v44 (~20,000):
- Read positions within gene body from the appropriate per-chr parquet
- Split by label (exon / intron / splice_donor / splice_acceptor)
- Compute per-region: `mean_c_t`, `mean_oscil`, `mean_amplitude`, `mean_argmin_layer`, `frac_committed`, `frac_deliberating`
- Combine into a per-gene 20-column feature vector

Output: `results/exp5/per_gene_gdtr_features.csv` (~20K rows × ~24 columns)

**Step 2 — GO term join** (~10 min including download)

- Query `mygene.info` for GO annotations of each Ensembl gene ID:
  - GO Biological Process (BP)
  - GO Molecular Function (MF)
  - GO Cellular Component (CC)
- Filter to "specific" terms (depth ≥ 4, or annotation term size 50-500 genes to avoid over-broad categories)

Output: `results/exp5/gene_go_annotations.csv`

**Step 3 — GO enrichment tests** (~30 min CPU)

Two complementary tests:

(a) **GSEA-style rank-sum**: For each GO term with ≥ 30 member genes, rank all genes by one of the gDTR features (e.g., splice_donor mean c(t) ascending). Test if GO-annotated genes are enriched at the top (Mann-Whitney U or GSEA). Bonferroni or Benjamini-Hochberg correct.

(b) **Predicted-term test**: Compute mean gDTR feature vector for the three predicted GO terms (0000398, 0006412, 0006351) and compare to genome-wide bootstrap mean. Report specific numbers matching or refuting H_B2.

Output: `results/exp5/go_enrichment_gsea.csv`, `results/exp5/go_predicted_test.json`

**Step 4 — Clustering + label preservation** (~15 min CPU)

- UMAP + HDBSCAN clustering on per-gene gDTR feature vectors
- Compute adjusted mutual information (AMI) between HDBSCAN clusters and top-100 GO terms
- Permutation control: shuffle gene labels 100× and recompute AMI → real AMI should be > 3σ above shuffled AMI

Output: `results/exp5/gene_clustering_ami.json`, `results/exp5/gene_umap_coords.csv`

**Step 5 — Downstream variant task** (~20 min CPU)

- For each variant in the 8,008-SNV cohort, add "parent gene's GO membership vector" as feature (multi-hot on top-100 enriched GO terms)
- Retrain STEP 3's 4-class subtype classifier with the added features
- Compare macro-F1 vs STEP 3 baseline (0.815)

Output: `results/exp5/subtype_task_with_go.json`

### B.4 Expected results

**H_B1**: HDBSCAN gives 20-50 clusters. AMI vs top-100 GO terms > 0.15 (real signal). Permuted labels give AMI < 0.05.

**H_B2**:
- mRNA splicing GO: splice_donor Δc(t) ≥ 0.2 layers below genome mean, Δoscil ≥ 0.05 above — expected direction consistent with paper §3.1 + H3a.
- Translation GO: argmin_layer shift toward 27 by ≥ 0.5 layers vs genome-wide coding_exon mean, frac_committed higher by ≥ 5 %. Direct extension of STEP 2 coding_exon finding.
- Transcription GO: exploratory, no specific prediction.

**H_B3**: modest but positive Δ macro-F1 (+0.01 to +0.03).

### B.5 Cost

- 5 min CPU per-gene aggregation
- 10 min GO annotation download
- 65 min CPU tests + clustering + downstream
- 20 min figures
- **~1.5 hr total, no GPU needed**

### B.6 Deliverables

- Scripts: `scripts/exp5/e5_00_per_gene_aggregate.py`, `e5_10_go_join_and_enrichment.py`, `e5_20_clustering.py`, `e5_30_variant_downstream.py`
- Results: `per_gene_gdtr_features.csv`, `gene_go_annotations.csv`, `go_enrichment_gsea.csv`, `go_predicted_test.json`, `gene_clustering_ami.json`, `subtype_task_with_go.json`
- Figures:
  - `EXP5_go_enrichment_volcano.png` — top GO terms (x = effect size, y = −log10 p)
  - `EXP5_gene_umap.png` — per-gene UMAP colored by top GO term
  - `EXP5_predicted_go_signatures.png` — 3 predicted GO categories × gDTR feature bar chart
- Manuscript sections:
  - **§3.11 "Gene-level functional coherence"** (new subsection)
  - **§3.6 Task E** (variant classifier with gene GO features, extension of Task B)

### B.7 What this adds to the manuscript

- Bridges **per-position gDTR features** to **per-gene biology** — a common reviewer concern for interpretability papers
- **Testable specific predictions** (splicing genes ↔ splice_donor shallow; translation genes ↔ argmin_layer 27) — makes gDTR falsifiable at a new level
- **Public tool value**: the per-gene feature CSV is a downloadable resource — 20,000 genes × 24 features
- **Novel downstream demonstration**: variant classifier improved by gene-level features derived from the same gDTR framework

---

## Combined execution plan

| Priority | Step | Cost | Blocks / blocked by |
|---|---|---|---|
| 1 | Experiment A.2.a — cosine baseline reconciliation | 5 min CPU | none |
| 2 | Experiment A.2.c — task-decomposition table (aggregate existing results) | 30 min | needs TDiG numbers imported |
| 3 | Experiment A.2.b — M3_geo × oscil join | 10 min CPU + TDiG parquet share | needs TDiG data access |
| 4 | Experiment B.1 — per-gene aggregation | 5 min CPU | uses existing WGS parquets |
| 5 | Experiment B.2 — GO join | 10 min | needs internet for mygene.info |
| 6 | Experiment B.3 — enrichment tests | 30 min CPU | after B.1 + B.2 |
| 7 | Experiment B.4 — clustering + AMI | 15 min CPU | after B.1 |
| 8 | Experiment B.5 — variant downstream w/ GO | 20 min CPU | after B.1 + B.3 |
| **Total** | | **~2 hours of compute** | mostly CPU |

Both experiments can be run in the next session as a single ~2-hour block.

---

## What each experiment defends against in review

| Reviewer concern | Defended by |
|---|---|
| "Why cosine and not one of the alternatives?" | Existing TDiG bootstrap CI + our reconciliation (A.2.a) + task-decomposition table (A.2.c) showing cosine wins the downstream task, curvature wins context discrimination — no single winner |
| "This is per-position toy analysis. What does it mean for biology?" | Experiment B GO enrichment shows gene-level functional coherence |
| "The variant classifier is on 15 cancer genes only. Generalise?" | Experiment B per-gene features are gene-agnostic (works for all ~20,000 genes) |
| "How would a biologist actually use gDTR?" | Experiment B provides a downloadable per-gene feature matrix + interpretable GO enrichment |

---

*Last updated: 2026-07-20.*
