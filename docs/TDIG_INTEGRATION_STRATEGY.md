# TDiG integration strategy — how the follow-up strengthens TDiG's contributions

**Purpose**: For each substantive TDiG contribution, this document describes (a) what TDiG already established, (b) what limitations TDiG itself acknowledges, and (c) how the WGS + downstream + threshold-crossing work in this repo strengthens or extends it into the integrated journal manuscript.

The goal is a manuscript in which every TDiG headline gains additional evidentiary weight, not one where the two contributions sit side-by-side without interaction.

Reference documents:
- [`INTEGRATED_MANUSCRIPT_PLAN.md`](INTEGRATED_MANUSCRIPT_PLAN.md) — §-by-§ manuscript outline
- [`NARRATIVE.md`](NARRATIVE.md) — story arc
- [`RESULTS.md`](RESULTS.md) — this repo's numbers
- TDiG upstream: [`docs/PLAN.md`](https://github.com/heoneyzi/TDiG/blob/main/PLAN.md), [`docs/metric_definitions.md`](https://github.com/heoneyzi/TDiG/blob/main/docs/metric_definitions.md), [`results/RESULTS_v3.md`](https://github.com/heoneyzi/TDiG/blob/main/results/RESULTS_v3.md)

---

## Overview — TDiG's 8 substantive contributions and their strengthening path

| # | TDiG contribution | TDiG scope | TDiG acknowledged limitation | Strengthening from this repo |
|---|---|---|---|---|
| **T1** | M1–M5 metric family + 3 reference variants | chr22 + chr17, 17 cells total | Def 3 (h_norm) has no metric in family | `oscil` axis partially fills Def 3; adds a **6th axis family** to the taxonomy |
| **T2** | 98-cell context × metric heatmap | chr17 + chr22 | not WGS; 7 context labels only | WGS replication on cells common to both; extended to cCRE-ELS + 5'/3' UTR |
| **T3** | Chr17 replication (13/13 sign preserved, ρ=0.989) | chr22 → chr17 (2 chromosomes) | not genome-wide | WGS extension: 24 chromosomes, per-metric variance across chr |
| **T4** | Per-consequence variant layer split (intron/3'UTR L=8; missense/syn L=27) | 8,008 SNVs (paper cohort) | only ΔH_norm feature | Combined with `oscil`-based features; per-consequence stratified WGS-normalisation |
| **T5** | L29 phase transition (context 0.980→0.799; variant 0.85→0.79) | chr22 probing | descriptive; mechanism from SVD only | Cross-referenced with `argmin_layer = 27.6` late-commit finding on WGS coding_exon |
| **T6** | PCA (PC1 = 75.7 % variance; bidirectional settling geometric origin) | 4,096-d hidden state, chr22 | PC1 alone ≠ context axis (0.706 vs full 0.980) | WGS-scale generalisation test; PC1 vs `oscil` interaction |
| **T7** | Intron-outlier functional element discovery (M5_tau_refB 2.70× near-splice) | chr22 intron only | not validated against cCRE / eQTL | WGS + cCRE-ELS join provides validation cohort |
| **T8** | VUS reclassification AUROC 0.949 | 5-fold CV on 8,008 labelled → 2,902 VUS | plateaus at 0.949 with layer-resolved features | Adds WGS-normalised H3 features → test whether AUROC 0.949 → higher |

---

## T1 — M1–M5 metric family + 3 reference variants

### What TDiG established
- Five axes of settling: direction (M1), magnitude (M2, diagnostic), geometry (M3, ref-free), Mahalanobis distribution (M4_set), path efficiency (M5).
- Three reference variants (A no-norm / B both-norm / C current-gDTR) to isolate RMSNorm γ asymmetry.
- 17 settling cells total, each with q70 calibration + persistence window W=3.
- Explicit 3-settling-definition table (Def 1 stops evolving / Def 2 → h₂₉ / Def 3 → h_norm).
- Honest acknowledgement: **no metric in M1–M5 addresses Def 3** (convergence to h_norm) — cos(h₂₈, h_norm) ≈ −0.013 makes direction-based metrics blind, magnitude metrics unusable because ‖h_norm‖ is tiny.

### Strengthening from this repo

**(a) The `oscil` axis partially resolves the Def 3 gap.**
`oscil` is defined against γ (not against h_norm), so it inherits none of h_norm's noise. But it reveals structure that Def 3-agnostic metrics miss: **at splice sites, `oscil` moves in the *opposite direction* to `c(t)`**. If Def 3 settling were the sole mechanism, `oscil` and `c(t)` should co-move. Instead, `oscil` sees the "continued reconsideration after commitment" phase — which is exactly the missing piece of the Def 3 gap description.

Manuscript action: add a 6th row to TDiG's Table 2 (metric × reference variants):

| Family | Cells | Definitions targeted |
|---|---|---|
| M1_dir × {A,B,C} | 3 | Def 2 / 2 / 3 (C degenerate) |
| M2_mag × {A} | 1 | residual accumulation (diagnostic) |
| M3_geo × 5 α/β | 5 | Def 1 |
| M4_set × {A,B,C} | 3 | Def 2 |
| M5_tau × {A,B,C} | 3 | Def 2 (Option B lock) |
| **H3 crossing family (new)** | **9** | **partial Def 3** (γ-referenced, h_norm-agnostic) |
| **Total** | **24** | complementary axes on the residual stream |

**(b) Anti-symmetry test on the reference variant matrix.**
TDiG's Ref A/B/C isolates RMSNorm γ. `oscil` is γ-agnostic (uses only the sign of D_cos − γ). So we predict `oscil` should be *reference-variant-invariant* — the same value regardless of A/B/C. If TDiG computes `oscil` under all 3 references and finds equal values, that's a strong validation that `oscil` and the metric family are decorrelated. This becomes a §2.5 "reference-variant matrix ablation" table row.

**(c) Persistence window W ablation validated at scale.**
TDiG's design v2 chose W=3 as production, with W∈{1,3,5} in supplementary. Our WGS 24-chr `c(t)` and `oscil` were both computed at W=1 (running-min for c_t; single crossing for oscil). A joint ablation table (paper cell × W ∈ {1,3,5} × chr subset) becomes an Extended Data table jointly authored by both teams.

---

## T2 — 98-cell context × metric heatmap

### What TDiG established
- 14 settling cells × 7 biological contexts (intron / coding_exon / 5'UTR / 3'UTR / splice_donor / splice_acceptor / intergenic) = 98-cell heatmap.
- Reveals metric-context dissociations invisible to any single-metric collapse.
- M3_geo family wins 13/21 context pairs as best discriminator; M5_tau_refB wins 6 splice-related pairs.

### Strengthening from this repo

**(a) WGS extension of the cells common to both frameworks.**
For any cell that is reference-free or uses Ref C (h_norm), we already have WGS-scale computation: `c(t)` and `oscil` families cover this. Concretely:
- `c(t)` = TDiG's M1_dir_refC (up to sign convention) → chr17 replication ρ = 0.989 (TDiG) → **WGS 24-chr replication d = −0.283 ± 0.051 for splice donor**
- `oscil`, `amplitude`, `argmin_layer` = new axes → **WGS Cohen's d matrix (24 chr × 9 axes × 7 contexts = 1,512-cell table)** as Extended Data

**(b) Adding cCRE-ELS as the 8th context.**
TDiG's 98-cell heatmap has 7 contexts (intron / coding_exon / 5'/3'UTR / splice sites / intergenic). Our cCRE-ELS WGS join (809,429 records, STEP 4) provides an 8th biological context with:
- WGS *d*(c_t) = −0.132 ± 0.044 (matches paper chr22-only −0.118)
- WGS *d*(oscil) = +0.117 ± 0.038 (**not in paper, not in TDiG**)

Reruning TDiG's 14-cell analysis on cCRE-ELS gives a **14 × 8 = 112-cell heatmap** for the integrated paper.

**(c) Coding exon "late clean deep-commit" as the anomalous discovery.**
Our STEP 2 finds coding_exon has argmin_layer = 27.6 (deeper than splice donor 23.5) and frac_committed = 43 % (splice donor 21 %) — the strongest single-cell effect in the WGS scan. TDiG's coding_exon vs intron d = −0.94 on M3_geo_a0.0_b1.0 (curvature-only) is the mirror image on a different metric. Together they establish coding_exon as having a distinctive "late, decisive, curvature-driven" signature. **This becomes a §3.5 dedicated result subsection**, not just a bullet in each source paper.

---

## T3 — Chr17 replication (Spearman ρ = 0.989)

### What TDiG established
- 17-cell settling on chr17 (27,586 windows) using chr22-frozen γ_v2 thresholds.
- Spearman ρ(chr22_d, chr17_d) = 0.989 across 13 valid cells.
- Median magnitude retention = 97.2 % (vs original paper's single-cell 94.6 %).
- Sign preserved: 13/13.

### Strengthening from this repo

**(a) WGS-scale replication test.**
Chr22 → chr17 is a 2-chromosome test. Our WGS 24-chr run provides a *distribution* of magnitudes and directions. For cells computable from our stored features (`c(t)`, `oscil`, `amplitude`, etc.):
- 6 of 7 contexts (paper §3.1) retain direction across all 24 chromosomes with SD 0.02–0.07 per context
- 5'UTR is the one exception — **chromosome-independent sign flip** (paper +0.20 → WGS mean −0.011). This is a stronger claim than TDiG could make with 2 chromosomes.

Manuscript action: extend TDiG's Table 3 with a "WGS replication SD across 24 chr" column, per cell. Cells with SD < 0.05 are labelled "genome-invariant"; cells with SD 0.05–0.15 are "chromosome-variable"; cells with sign disagreement across chromosomes are flagged.

**(b) γ chromosome-invariance is much stronger than the workshop paper's Appendix A.2 grid.**
STEP 5 finds per-chromosome q70 recalibrated γ has SD = 0.15 % across 24 chromosomes. This is a stronger version of TDiG's chr17 transferability claim (single number matched) and of the workshop paper's A.2 flat-plateau claim (5×5 grid stable). **Combined, all three imply: γ = 0.397 is a hyperparameter with essentially zero cross-chromosome variance under the paper's calibration protocol.**

---

## T4 — Per-consequence variant layer split (intron/3'UTR L=8; missense/synonymous L=27)

### What TDiG established
- On 8,008 SNVs × 15 cancer genes, ΔH_norm_L1 AUROC per layer per consequence.
- Sequence-level information (intron, 3'UTR) → best AUROC at L=8.
- Protein-semantic information (missense, synonymous, 5'UTR) → best AUROC at L=27.
- Confirms workshop paper §3.3 at scale (paper showed group-median class ordering; TDiG shows per-variant AUROC).

### Strengthening from this repo

**(a) Our H2b + H3 features complete the per-consequence axis assignment.**
TDiG shows *when* each consequence class's ΔH peaks. Our H3d shows *which axis* (early_below vs late_below vs amplitude) each class dominates on. Combined stratification:

| Consequence | TDiG peak layer | Followup dominant axis | Together suggest |
|---|---|---|---|
| intron | L=8 (ΔH) | (baseline; no strong axis) | Early sequence-level flag |
| 3'UTR | L=8 (ΔH) | mild `d(oscil) = +0.086` | Early flag + mild reconsideration |
| synonymous | L=27 (ΔH) | high `late_below_frac` | Late clean commit (codon-frame) |
| missense | L=27 (ΔH) | `late_below_frac` + `amplitude` | Late deep commit (protein-level) |
| canonical_splice | L=27 (ΔH) + L=8 (splice donor position) | High `oscil` + late `first_enter` | Multi-axis: motif + integration + reconsideration |

This is a **new integrated result** neither team could produce alone.

**(b) WGS-normalisation improves the 4-class subtype classifier from 0.641 → 0.815 (Task B).**
TDiG's per-consequence analysis is descriptive (per-layer AUROC per class). Our STEP 3 makes it *predictive*: given a variant, classify its consequence subtype from Δ features alone. Combining:
- cos32 vector: macro-F1 = 0.641
- + WGS-normalised H3 features: 0.686
- + all combined: **0.815** (+17.4 %p over paper baseline)

The gain is *asymmetric* across consequences: canonical_splice benefits most from WGS-normalisation because its context background (splice_donor cells in the WGS distribution) is atypical relative to the intron reference the paper originally used. **This is a new §3.6 Task B result that TDiG's per-consequence analysis explains.**

---

## T5 — L29 phase transition

### What TDiG established
- Context probing AUROC on residual stream: **0.980 at L=27 → 0.799 at L=29** (drop of 0.181 in two layers).
- Variant AUROC: **0.85 → 0.79** at same layer.
- PC1-metric correlations collapse simultaneously.
- L29 is Evo 2's rotation layer where the residual stream is projected into h_norm's frame.
- SVD mechanism (§H3 in RESULTS_v3): L15→L16 R² = 0.997, cond# = 12.8, smooth rotation; presumably L28→L29 is anomalous.

### Strengthening from this repo

**(a) `argmin_layer = 27.6` for coding_exon is the mirror image of the L29 crash.**
TDiG shows *what breaks* at L29 (context probing). Our STEP 2 shows *what commits* just before L29 (coding_exon argmin at layer 27.6, splice donors at 23.5). The L29 phase transition is not a random crash — it's the point where the residual stream terminates its computation and rotates. Everything with biological meaning has committed by L28. **This unifies TDiG's phase-transition finding with our late-commit finding into one story.**

**(b) `first_enter_layer` distribution across contexts.**
Our WGS `first_enter_layer` mean for splice donors is 22.4 (chr22 sample); TDiG reports variant AUROC peaks at L=27–29 for protein-semantic classes. **These two numbers are compatible if the crossing → commit → rotation pipeline is: enter γ around L22, hold below γ until L27–28, rotate at L29–30, project to output.** This is a mechanistic story that neither dataset alone tells.

Manuscript action: dedicate §3.7 to L29 phase transition, showing (a) TDiG's context probing crash, (b) our late-commit distribution, (c) mechanistic proposal: L22 first-cross → L27 stable → L29 rotation. Include a Figure with layer-index x-axis and 3 stacked panels.

---

## T6 — PCA (PC1 = 75.7 % variance)

### What TDiG established
- On chr22 residual stream, layer-centered PCA gives PC1 = 75.7 % of variance.
- PC1 alone AUROC = 0.706 (max at L=30) vs full 4096-d AUROC = 0.980 at L=27 — so PC1 is NOT the context axis.
- PC1 |r| = 0.85 with M1_dir_refC at L=4 — PC1 encodes bidirectional settling in a single geometric direction.
- PC1+ encodes Def 1 (M3_geo); PC1− encodes Def 2 (M1_dir, M2, M5_tau_refB).

### Strengthening from this repo

**(a) `oscil` is expected to load on a *different* PC than PC1.**
PC1 = bidirectional settling axis (Def 1 vs Def 2). `oscil` measures *count of crossings*, which is orthogonal to settling direction. Prediction: `oscil` should have |r| < 0.3 with PC1 across layers, and > 0.5 with some higher PC (PC3 or PC4). **This is a testable falsifiable prediction that becomes a §3.8 subresult.**

**(b) WGS scale strengthens PCA generalisability.**
TDiG PCA is on chr22 residual stream. We don't store raw h_ell at WGS scale (would be 200 GB), but we can:
- Verify PC1 = 75.7 % variance replicates on chr17 (already have the h_ell there)
- Use our WGS-scale `c(t)` and `oscil` distributions to identify positions that lie on PC1 vs PC3, then compare their biological annotation enrichment

**(c) Sparse-autoencoder future work becomes concrete.**
TDiG acknowledges "superposition disentanglement" as open item (design v2 §9). Our downstream benchmark provides a *use case* for SAE features: if learned dictionary features carry more subtype-classification signal than raw PCA components, that's a concrete SAE motivation. **§4 Discussion adds this as a natural next step.**

---

## T7 — Intron-outlier functional element discovery (M5_tau_refB 2.70× near-splice enrichment)

### What TDiG established
- 31.5 million intron tokens on chr22; top 0.5 % of M5_tau_refB (i.e., 157,283 outliers) are enriched **2.70× vs random intron** for being within 200 bp of an annotated splice site.
- Candidate biology: splicing silencers/enhancers (ISS/ISE), branchpoint sites, cryptic splice donors, polypyrimidine tracts.
- **This is TDiG's flagship "previously-impossible discrimination" result.** Basic VEP/AlphaMissense don't distinguish within annotated "intron" regions.

### Strengthening from this repo

**(a) WGS 24-chr replication of the intron-outlier test.**
TDiG tested one chromosome (chr22). Our WGS parquets have per-position `c_t`, `oscil`, `amplitude`, `argmin_layer` for all 1.67 B intron positions across 24 chromosomes. Compute analogous outlier-enrichment tests for each of our axes:
- Top 0.5 % `c(t)` shallowest intron positions → near-splice enrichment per chromosome
- Top 0.5 % `oscil` highest intron positions → near-splice enrichment
- Top 0.5 % `amplitude` deepest intron positions → near-splice enrichment

Prediction: `oscil` outliers should also be near-splice enriched (splice donors have d(oscil) = +0.22 genome-wide → outlier tail hits real splice-adjacent positions). **A per-axis × per-chromosome outlier enrichment table is a strong Extended Data figure.**

**(b) Cross-validation against cCRE + external annotations.**
TDiG's future work notes "compare M5_tau_refB outlier positions to ENCODE cCRE, GTEx splice-QTL, SpliceAI-flagged". Our cCRE-ELS WGS join (STEP 4) already has these positions in the same coordinate system. **Direct join: for each TDiG-identified outlier on chr22, is it inside a cCRE-ELS region? Fraction should be > baseline.** This closes TDiG's own open item.

**(c) Novelty test with GENCODE `intergenic` positions.**
TDiG restricted the outlier analysis to annotated introns. Repeat on annotated intergenic positions (1.18 B WGS-wide). Enrichment for cCRE-ELS overlap, splice-adjacent, GTEx eQTL positions gives a "functional element density" score for the *unannotated* genome. **This is a natural §5 Discussion talking point about what the framework can discover.**

---

## T8 — VUS reclassification (AUROC 0.949)

### What TDiG established
- 5-fold CV on 8,008 labelled ClinVar SNVs (P/LP vs B/LB).
- Layer-resolved features (32-d log(ΔH_norm_L2) + 32-d Δcos = 64-d) with LR/GBM.
- **AUROC 0.949**, +9.4 %p over ΔH-at-L=8 baseline (0.855).
- Applied to 2,902 VUS: 1,303 Likely Pathogenic, 853 Uncertain, 746 Likely Benign.

### Strengthening from this repo

**(a) Add WGS-normalised H3 features → test whether 0.949 → higher.**
TDiG's 64-d feature is layer-resolved but H3-agnostic. Our STEP 3 shows for the 4-class subtype task, adding WGS-normalised H3 gives +17.4 %p. Question: does the same lift help pathogenicity (P/LP vs B/LB)? Our STEP 3 result on this shows **binary AUROC saturates at cos32 alone (0.844 → 0.849)** — but that's without TDiG's layer-resolved ΔH features. Combined feature set:

| Feature set | Binary AUROC | Task |
|---|---|---|
| ΔH L=8 (paper single-layer) | 0.856 | binary P/LP vs B/LB (TDiG) |
| ΔH 64-d resolved (TDiG G3) | 0.949 | binary (TDiG) |
| cos32 (paper §App C) | 0.844 | binary (paper + this repo verify) |
| cos32 + H3 + scalars + WGS-norm (this repo STEP 3) | 0.849 | binary |
| **ΔH 64-d + WGS-norm H3 (proposed)** | **?** | binary |
| **ΔH 64-d + WGS-norm H3 + cos32** | **?** | binary |

This combined benchmark table is a §3.6 headline. Even if the gain is smaller than the subtype task, showing 0.949 → 0.955+ demonstrates the integrated framework provides monotonic improvement.

**(b) Task-decomposition figure.**
Show that different feature families dominate different tasks:
- Binary pathogenicity: TDiG's ΔH 64-d wins (0.949) → magnitude-of-perturbation dominates
- 4-class subtype: our WGS-normalised H3 wins (+17.4 %p) → shape-of-trajectory dominates
- Consequence layer split: TDiG's per-layer ΔH wins → layer-index dominates

**One framework, three specialised readouts, each with a distinct dominant feature.** This is the manuscript's "unified but not homogeneous" claim.

**(c) VUS predictions as a supplementary resource.**
TDiG produces 2,902 VUS predictions. We recompute with the enhanced feature set → published as a supplementary CSV, curated for downstream users. Cross-check against ClinVar 2026-current for any VUS that has since been reclassified — if the two frameworks agree on VUS-→likely-P/LP, that's independent evidence of both.

---

## Priority order for merging TDiG artifacts

If we merge in stages (rather than all at once), the sequence that gives the manuscript the most immediate lift:

1. **T5 L29 phase transition** (§3.7) — pairs immediately with our late-commit finding; one figure for both.
2. **T4 per-consequence layer split** (§3.6 Task C) — pairs with our WGS-normalised subtype benchmark; one integrated Table.
3. **T2 98-cell heatmap** (§3.2) — flagship figure of the manuscript; needs both teams' cells side-by-side.
4. **T8 VUS reclassification** (§3.6 Task A) — extended feature set benchmark.
5. **T7 intron-outlier discovery** (§3.9) — moved from workshop paper §5 Q2 hint to headline result.
6. **T6 PCA** (§3.8) — theoretical bridge between metric axes and geometric structure.
7. **T3 chr17 replication** (§3.10 Calibration) — merged with our γ recalibration.
8. **T1 M1–M5 taxonomy** (§2 Methods) — the framework itself; needed but bookkeeping-heavy.

---

## Coordination checklist before importing TDiG artifacts

The following must be settled with the TDiG team (heoneyzi et al.) before any artifact is copied here:

- [ ] Co-authorship agreement (first-author sharing / middle-author list / senior)
- [ ] Access to TDiG chr22 raw cache (chr22_cache.h5, tier-3 raw hidden states, tier-2 scalars) for re-running M1–M5 with the same random seed used here
- [ ] Agreement on which TDiG results (from RESULTS_v3.md) go into the main text vs Extended Data
- [ ] Figure attribution convention
- [ ] Final metric naming (TDiG uses M1–M5; we could rename or preserve)
- [ ] Agreement on the manuscript's opening claim (whose framing "wins" — TDiG's "spread it out" or our "unified but not homogeneous")

Once these are settled, populate `tdig_integration/` with the imported artifacts and update this document with actual file paths.

---

*Last updated: 2026-07-20.*
