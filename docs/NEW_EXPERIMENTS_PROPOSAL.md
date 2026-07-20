# Proposed new experiments — metric justification + biological utility

Two new experiment families to add to the integrated manuscript, addressing the two most common reviewer concerns for interpretability-lens papers:

- **"Why this metric and not another?"** — cosine distance is the paper's default. This section shows it is not arbitrary by benchmarking against L2, scaled-L2, Mahalanobis, angular, and correlation on the same positions, and by direct predictions the paper's Method makes about their relative behaviour.
- **"How would a biologist actually use this?"** — the paper's downstream demonstration is a variant classifier. This section adds a Gene Ontology (GO) enrichment analysis showing the gDTR signal is functionally coherent at the gene level, not just at the individual position level.

Both are designed to (a) run within the existing infrastructure, (b) produce results usable in §3 Results as standalone subsections, and (c) address specific limitations the workshop paper Appendix G already flags.

Reference documents:
- [`INTEGRATED_MANUSCRIPT_PLAN.md`](INTEGRATED_MANUSCRIPT_PLAN.md) — where these experiments land
- [`METHODS.md`](METHODS.md) — metric definitions
- [`RESULTS.md`](RESULTS.md) — existing results

---

## Experiment A — "Why cosine?" metric-choice justification

### A.0 Motivation

The workshop paper §2 defines `c(t)` using cosine distance
$$D_{\cos}(\ell, t) = 1 - \cos(h_\ell(t), h_{\mathrm{norm}}(t))$$
and Appendix G.1 argues cosine is chosen "on purpose" because Evo 2's residual norm grows monotonically across layers, so a raw distance conflates *directional commitment* with *residual accumulation* (an architectural artifact). But the paper does not directly *compare* cosine to alternatives.

Reviewers will ask: what if L2 works as well? What if Mahalanobis is better? What if the paper's whole story could be told with any distance metric?

This experiment answers: no. Different metrics give quantitatively different biological signals, and cosine is the one that separates splice sites cleanly without conflating direction with residual accumulation. It also makes TDiG's M4_set (Mahalanobis) an alternative worth reporting for completeness.

### A.1 Hypothesis

**H_A**: On the same chr22 positions, the per-context Cohen's d for splice_donor vs intron depends on which distance metric is used to define `c(t)`. Cosine and Mahalanobis (M4_set) give the strongest and cleanest signal; raw L2 conflates settling with residual accumulation and gives a weaker or opposite signal.

**H_A refutability**: if L2 and cosine give the same |d|, or if L2 is stronger, the paper's cosine choice is retrospectively unjustified.

### A.2 Method

**Sample**: 100 random chr22 6 kb windows (same seed 42 as H3a's original 100). Re-forward through Evo 2 7B with **raw h_ℓ retained** (~15 min on 2 × B200, ~1.5 GB storage per 100 windows × 32 layers × 3000 pos × 4096 dim × float16 = 78 MB per window × 100 = ~7.8 GB — reuse existing chr22 raw cache from TDiG if available, else regenerate).

Compute 6 distance metrics between each h_ℓ(t) and h_norm(t):

| Symbol | Formula | Property | Notes |
|---|---|---|---|
| D_cos (baseline) | 1 − cos(h_ℓ, h_norm) | scale-invariant | paper's default |
| D_L2 | ‖h_ℓ − h_norm‖₂ | raw | conflates direction + magnitude |
| D_L2_scaled | D_L2 / ‖h_norm‖ | scale-invariant | ratio, still sensitive to sign |
| D_Mahal (M4_set A) | √((h_ℓ − h_norm)ᵀ Σ⁻¹ (h_ℓ − h_norm)) | whitened | TDiG M4_set, requires Σ_ref |
| D_ang | arccos(1 − D_cos) / π | direction-only, [0,1] | monotone to cosine, sanity check |
| D_corr | 1 − corr(h_ℓ, h_norm) | mean-centered cosine | robust to constant offset |

For each metric, define `c_metric(t) = first ℓ where running-min(D_metric) ≤ γ_metric`, with γ_metric = q70 of D_metric at penultimate layer over chr22 intron positions (same protocol as paper §2 for cosine).

Compute per-context d for splice_donor vs intron for each metric.

**Additional analyses**:
- **Cross-metric Spearman ρ**: 6 × 6 matrix. Predict D_cos and D_ang should be perfectly correlated (they are monotone transforms of each other); D_L2 and D_L2_scaled should correlate strongly; D_Mahal should be closer to D_L2 than to D_cos in behavior but still capture direction.
- **Superposition sanity**: pairs of positions where cos(h_A, h_B) ≈ 1 but the two tokens have very different biological contexts. Fraction of these positions should be non-negligible — establishing that cosine ≠ semantic identity (paper §2 acknowledges this; we measure it).
- **Permutation control**: shuffle the biological labels, recompute d — should collapse to |d| < 0.05 for every metric. Establishes that the biological signal is not a metric artifact.

### A.3 Expected result (predictions)

Based on paper §2 + Appendix G + TDiG's M2 discussion:

| Metric | Predicted \|d\| at splice_donor vs intron | Rationale |
|---|---|---|
| D_cos | ≈ 0.30 (matches paper -0.354 at 100-window scale) | baseline |
| D_L2 | ≪ 0.30 or opposite sign | residual accumulation dominates |
| D_L2_scaled | 0.15 − 0.25 | direction signal recoverable but degraded |
| D_Mahal | 0.20 − 0.35 | competes with cosine, likely comparable |
| D_ang | ≈ 0.30 (identical to D_cos) | monotone transform, same information |
| D_corr | 0.20 − 0.30 | mean-centering doesn't matter much for h_norm ≈ mean-zero |

**Cross-metric ρ**: D_cos and D_ang = 1.00 exactly; D_L2 and D_L2_scaled ≥ 0.85; D_cos and D_L2 = 0.3 − 0.6 (large but not perfect — confirms they measure related but distinct things); D_cos and D_Mahal = 0.6 − 0.8.

**Permutation control**: |d| collapses to < 0.05 for every metric after label shuffle.

### A.4 Cost

- 15 min on 2 × B200 for re-forward with raw h_ℓ retained (or free if TDiG chr22 cache is accessible)
- 5 min CPU for metric computation
- 2 min for figure

### A.5 Deliverable

- Script: `scripts/exp4/e4_00_metric_ablation.py`
- Result: `results/exp4/metric_ablation_context_d.csv` (5 metrics × 7 contexts × d + p + n)
- Cross-metric ρ heatmap + per-context d bar chart: `results/exp4/EXP4_metric_ablation.{png,pdf}`
- Superposition sanity table: `results/exp4/superposition_check.json`
- Permutation control: `results/exp4/permutation_control.csv`
- Section in manuscript: **§3.1.1 "Why cosine? A metric-choice ablation"** (subsection of §3.1 baseline)

### A.6 What this claims to add to the manuscript

- **Retrospective justification** of paper's cosine choice (currently only motivated verbally in §2)
- **Direct comparison to TDiG's M4_set Mahalanobis** (integrates the two contributions on a common axis)
- **Reviewer defense**: pre-empts the "why not L2?" question with data

---

## Experiment B — Gene Ontology (GO) functional coherence

### B.0 Motivation

The workshop paper demonstrates gDTR's utility with a variant classifier (§App C). That's *per-position* utility — for a specific base substitution, gDTR features help predict pathogenicity. But biological interpretation happens at the *gene* level. A reviewer or biologist reading the paper will naturally ask:

> "If gDTR features are biologically meaningful, do genes with similar gDTR profiles share similar Gene Ontology functions?"

This experiment answers that question. If gDTR features aggregated per-gene are enriched for coherent GO categories, the framework is *functionally coherent*, not just position-level noise averaging out to a signal. It also opens a downstream use case: **gene function prediction from residual-stream trajectory alone**, no protein sequence needed.

### B.1 Hypotheses

**H_B1** (necessary condition): Per-gene aggregate gDTR features (mean c(t), mean oscil, mean amplitude across the gene's exons and introns) cluster genes into functionally coherent groups. Random shuffling of gene labels destroys this clustering.

**H_B2** (specific prediction): Genes annotated with certain GO categories — specifically GO:0006397 (**mRNA processing**), GO:0000398 (**mRNA splicing via spliceosome**), GO:0006412 (**translation**), GO:0006351 (**transcription, DNA-templated**) — show gDTR profile enrichment consistent with the paper's context-level findings. Splicing-related GO genes should show shallower `c(t)` at their splice sites and higher `oscil` (matching H3a); translation-related GO genes should show the coding-exon "late deep commit" pattern (matching STEP 2 argmin_layer 27.6).

**H_B3** (downstream utility): Adding per-gene gDTR features to a variant classifier improves subtype AUROC beyond per-position features alone.

**H_B refutability**: if random gene-label shuffling produces the same enrichment, or if predicted GO categories don't emerge, or if per-gene features add nothing to per-variant features, the functional coherence claim is undermined.

### B.2 Data needed

- **WGS 24-chr per-position parquets** (already have — 37 GB)
- **GENCODE v44 gene body coordinates** (already have gencode.v44.annotation.gtf, 1.5 GB)
- **GO annotations per gene**: download from Ensembl BioMart or GO Consortium (~50 MB, biomart or `mygene.info` API)
- **Optional: PhyloP** for conservation cross-reference (10 GB, only if H_B4 pursued)

### B.3 Method

**Step 1 — per-gene feature aggregation**
For each protein-coding gene in GENCODE v44:
- Find its position range on the appropriate chromosome
- Aggregate over all per-position parquet rows within the gene body:
  - `mean_c_t`, `mean_oscil`, `mean_amplitude`, `mean_argmin_layer`
  - Separately for exons vs introns (using position labels)
  - Per splice donor and splice acceptor within the gene (mean over 60,507 + 59,370 splice sites across chr1)
  - Fraction of positions that are "committed" (n_enter=1, n_exit=0)
  - Fraction "deliberating" (oscil ≥ 1)
- Output: gene × feature matrix, one row per gene, ~20,000 genes × 20 features

**Step 2 — GO term join**
- Query GO annotations via `mygene.info` (Python client) or Ensembl BioMart (bulk download):
  - GO Biological Process (BP)
  - GO Molecular Function (MF)
  - GO Cellular Component (CC)
- For each gene, get its list of associated GO terms
- Restrict to "specific" GO terms (depth ≥ 4 in DAG, or "leaf" annotations) to avoid over-broad categories like "protein binding"

**Step 3 — GO enrichment tests**

Two complementary tests:

(a) **GSEA-style enrichment**: For each GO term with ≥ 30 annotated genes, rank all genes by mean_c_t (ascending). Test whether GO-annotated genes are enriched at the top (shallower settling) using GSEA / rank-sum test. Bonferroni or FDR-correct across all GO terms tested.

(b) **Per-GO-term profile**: For each significantly enriched GO term, compute the mean gDTR feature vector across its member genes. Visualize as GO × feature heatmap.

**Step 4 — Clustering + label preservation test**
- UMAP + HDBSCAN clustering on per-gene gDTR feature vectors
- Compute clustering-adjusted mutual information (AMI) between HDBSCAN clusters and top-100 GO terms
- Permutation control: shuffle gene labels, recompute AMI — real AMI should be significantly above shuffled AMI

**Step 5 — Downstream variant task**
- For each variant in the 8,008-SNV cohort, add "parent gene's GO profile" as feature (e.g. one-hot encoding for top-10 enriched GO terms it belongs to)
- Compare 4-class subtype macro-F1: cos32 + WGS-norm + scalars (STEP 3 baseline 0.815) vs. + gene GO features
- Report ΔF1

### B.4 Expected results

**H_B1**: HDBSCAN on per-gene gDTR features should yield 20-50 clusters. AMI vs top-100 GO terms should be > 0.15 (real signal). Permuted labels should give AMI < 0.05.

**H_B2 predictions**:
- GO:0000398 (mRNA splicing) enrichment: genes annotated for this should have **splice_donor mean c_t** lower than the genome-wide splice_donor mean by 0.2-0.5 layers, and mean_oscil higher by 0.05-0.10. Confidence: moderate — splice-related genes may exemplify canonical splicing signal.
- GO:0006412 (translation) enrichment: gene body **mean_argmin_layer** should be closer to 27 (STEP 2 coding_exon signature). Confidence: high — direct extension of STEP 2 finding.
- GO:0006351 (transcription) enrichment: broader signal, could be either direction. Confidence: low; treated as exploratory.
- Enrichr / GProfiler independent replication for validation.

**H_B3**: adding per-gene GO features should give modest improvement (+0.01-0.03 macro-F1 on subtype). Small because variant-level features (Δoscil etc) already contain gene identity implicitly.

### B.5 Cost

- 5 min CPU to aggregate per-gene features from 37 GB WGS parquets
- 10 min for GO annotation download (`mygene.info` API for 20,000 genes)
- 30 min CPU for GSEA + enrichment tests + clustering + permutation controls
- 20 min for figures
- **~1.5 hr total, no GPU needed**

### B.6 Deliverable

- Script: `scripts/exp5/e5_00_per_gene_aggregate.py`
- Script: `scripts/exp5/e5_10_go_enrichment.py`
- Script: `scripts/exp5/e5_20_go_clustering.py`
- Script: `scripts/exp5/e5_30_go_variant_downstream.py`
- Results:
  - `results/exp5/per_gene_gdtr_features.csv` (~20K rows × 20 columns)
  - `results/exp5/go_term_enrichment.csv` (per-GO term rank-sum p, effect size, member count)
  - `results/exp5/go_gdtr_heatmap.csv` (GO × feature)
  - `results/exp5/gene_clustering_ami.json` (real vs permutation AMI)
  - `results/exp5/variant_downstream_with_go.json` (subtype F1 with vs without GO features)
- Figures:
  - `EXP5_go_term_enrichment_volcano.png` (top enriched GO terms)
  - `EXP5_gene_clustering_umap.png` (per-gene UMAP colored by top GO term)
  - `EXP5_go_gdtr_heatmap.png` (predicted category × feature signature)
- Sections in manuscript:
  - **§3.11 "Gene-level functional coherence"** (new subsection)
  - **§3.6 Task D** (extension of variant task with GO features)

### B.7 What this adds to the manuscript

- **Bridges per-position gDTR to per-gene biology** — closes a common reviewer concern for interpretability papers
- **Testable specific biological predictions** (splicing genes ↔ splice_donor shallow; translation genes ↔ argmin_layer 27) — makes gDTR falsifiable at a new level
- **Novel downstream demonstration**: variant classifier improved by gene-level features derived from the same gDTR framework
- **Public tool value**: `results/exp5/per_gene_gdtr_features.csv` is a resource — 20,000 genes × 20 features, downloadable, usable in any GO enrichment tool

---

## Combined execution plan

| Priority | Experiment | Cost | Blocks/blocked by |
|---|---|---|---|
| 1 | Experiment A "why cosine" | 15 min GPU + 20 min CPU | none (uses existing chr22 windows) |
| 2 | Experiment B — per-gene aggregation | 5 min CPU | uses existing WGS parquets |
| 3 | Experiment B — GO join + enrichment | 40 min CPU | requires GO download (10 min) |
| 4 | Experiment B — variant downstream w/ GO | 20 min CPU | requires B1 + B2 |
| **Total** | **~2 hours of compute** | | (mostly CPU) |

Both experiments produce results that map to specific manuscript sections and address specific reviewer concerns. They can be run in the next session as a single ~2-hour block.

---

## What each experiment defends against in review

| Reviewer concern | Defended by |
|---|---|
| "Why not L2? Why not Mahalanobis? Why not correlation?" | Experiment A metric ablation |
| "Cosine could be arbitrary; you should test alternatives" | Experiment A + integrates TDiG's M4_set |
| "This is per-position toy analysis. What does it mean for biology?" | Experiment B GO enrichment |
| "The variant classifier is on 15 cancer genes only. Generalise?" | Experiment B per-gene features are gene-agnostic |
| "How would a biologist actually use gDTR?" | Experiment B provides a downloadable per-gene feature matrix + interpretable GO enrichment |

---

*Last updated: 2026-07-20.*
