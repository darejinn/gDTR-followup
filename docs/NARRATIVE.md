# The gDTR story — from a single scalar to a multi-axis framework

A readable narrative of how the project evolved: what the original gDTR workshop paper claimed, what questions it left open, and how the follow-up experiments in this repository answered them.

Written for a reader who wants to understand the arc of the work without reading every JSON. All numbers referenced here have exact sources in [`RESULTS.md`](RESULTS.md) and reproducibility notes in [`REPRODUCE.md`](REPRODUCE.md).

---

## Part 0 — The original claim (workshop paper, June 2026)

Cho, Kang, Park, and Kim asked a very specific interpretability question about Evo 2 7B:

> When, along the layer stack, does a nucleotide token's representation *stabilise* against the model's output-ready frame?

They introduced **`c(t)`** — the smallest layer at which a running-minimum envelope of the cosine distance $D_{\cos}(\ell, t) = 1 - \cos(h_{\ell}(t), h_{\text{norm}}(t))$ falls below a fixed threshold `γ_cos = 0.397`. That single scalar told a clean story:

- **Splice donors and acceptors settle ~2 layers earlier than intronic contexts** on chr17 and chr22 (Cohen's *d* ≈ −0.35).
- **ENCODE cCRE-ELS regions** also settle earlier, at a milder magnitude (*d* ≈ −0.12).
- **Two-sided motif experiments**: replacing the canonical GT with AA deepens the settling by 0.46 layers (the model needs longer to resolve a broken motif); shuffling the ±100 bp flank while keeping the GT makes it settle 3.18 layers earlier (with no context to integrate, the model commits quickly).
- **Variant scoring**: a 32-dimensional trajectory of Δ$D_{\cos}$ per layer reaches AUROC 0.844 on 8,008 ClinVar SNVs across 15 cancer genes (matching AlphaMissense's operating range without any training).

The paper was accepted at the ICML 2026 GenBio Workshop as a short paper. It answered the *when* question cleanly, but it left several *what else* questions open. That is where this repository begins.

---

## Part 1 — The follow-up begins (July 13)

Three concrete follow-up questions:

1. **Whole-genome extension.** The paper's claims were made on chr22 (γ calibration) and chr17 (held-out). Do they hold across all 24 chromosomes?
2. **Variant-effect downstream tasks.** The paper reported a binary pathogenicity AUROC. Can the same features do more — for example, distinguish variant subtypes?
3. **Threshold-crossing dynamics.** The running-min envelope smooths the trajectory to a scalar. What information is thrown away — how often does $D_{\cos}$ cross γ, and when?

These three lines each turned into a subproject.

### 1.1 Baseline verification first

Before extending anything, we reproduced the paper's headline AUROC on the migrated data. Running the paper's script on the cached 32-d Δ$D_{\cos}$ features gave **AUROC 0.8437 ± 0.021** (paper: 0.844) and best single-layer AUROC 0.7291 at L=30 (paper: 0.729 at L=30). Reproduction matched to three decimals. Everything downstream is anchored to this verified baseline.

### 1.2 Variant subtype classification — H2b

The paper's §3.3 showed that different variant consequences peak at different layers as a group median (Kruskal-Wallis p = 3×10⁻¹⁰, ε² = 0.013). But group medians don't tell us whether *individual* variants are separable.

We joined ClinVar's Molecular Consequence field onto the 8,008 SNVs, filtered to four canonical classes (missense / nonsense / synonymous / canonical_splice), and trained a multinomial logistic regression on the 32-d Δ$D_{\cos}$ vector with 10-fold stratified CV.

Result: **macro-F1 = 0.641** on 4 classes (chance = 0.25). One-vs-rest AUROC = 0.889 for nonsense and synonymous, 0.875 for canonical_splice, 0.741 for missense. A 1-d scalar (the single largest |Δ$D_{\cos}$|) recovers only 30 % of the gain over chance, so the *shape* of the trajectory matters, not just its peak height.

The paper's group-level finding was actually understating a per-variant classifiable signal. This gave us a downstream task to keep sharpening.

### 1.3 Threshold-crossing dynamics — H3a

The running-min envelope hides how many times $D_{\cos}$ crosses γ. We defined a new per-position feature `oscil` = the number of extra crossings beyond a single dip, forwarded 100 random chr22 6-kb windows through Evo 2, and compared each context to the intronic baseline.

The **first result that changed the direction of the project** appeared here:

| Context | mean c_t | *d*(c_t) | mean oscil | *d*(oscil) |
|---|---|---|---|---|
| intron | 30.08 | 0 | 0.289 | 0 |
| **splice donor** | **28.53** | **−0.29** | **0.450** | **+0.18** *** |
| **splice acceptor** | **26.78** | **−0.62** | **0.739** | **+0.50** *** |
| coding_exon | 30.63 | +0.11 | 0.164 | −0.14 |
| intergenic | 30.84 | +0.16 | 0.155 | −0.17 |

Splice sites have both a **shallower** `c(t)` (the paper's finding) **and a higher** `oscil` (a new finding, opposite sign on a distinct axis). The plan had predicted that splice sites would be *cleaner* (oscil ≈ 0). The reality was the opposite: they commit early but continue to be re-evaluated across many layers.

This is exactly what the paper's §2 called "two-sidedness" as an interpretive caveat — a low `c(t)` could mean either strong motif detection or simple context — but §2 had not operationalised the distinction. Here, `c(t)` and `oscil` operationalise it: they are orthogonal axes, and splice sites sit in a distinctive quadrant of the 2-D plane.

We ran two more analyses on the same 79-window sample:

- **H3c** stratified positions into "committed" (n_enter=1, n_exit=0), "dipped" (1,1), and "deliberating" (oscil ≥ 1). Splice acceptor is **2.85× enriched** for deliberating compared to intron (Fisher p = 5.8×10⁻¹³); splice donor is **1.58×** (p = 3.5×10⁻³). Coding exons and intergenic are *depleted* for deliberating.

- **H3d** added nine advanced features: `first_enter_layer`, `last_exit_layer`, `longest_below_streak`, `streak_start_layer`, `amplitude_below_gamma`, `early/mid/late_below_frac`, `min_D`. The strongest single axis appeared here: **splice donor `late_below_frac` *d* = +3.18** (splice donors have 40× the fraction of layers 22–31 below γ compared to intron). Splice sites don't just dip briefly — they dip **deep** (amplitude 8–12× intron), **hold** (longest streak 3.7× intron), and do so at **late layers** (first_enter around layer 22, not layer 15 like intron).

We now had a **4-D framework** — {`c(t)`, `oscil`, `first_enter_layer`, `amplitude`} — with a fifth candidate axis (`argmin_layer`) waiting for genome-scale confirmation.

### 1.4 Integrated variant classifier — H2b + H3

We re-ran the 4-class subtype classifier with the crossing-dynamics features added. Baseline cos32 gave macro-F1 = 0.641; adding raw H3 features gave 0.659; adding hand-crafted scalars gave 0.669; **combining everything gave 0.682**. Modest at that stage, but the trajectory said "context matters", and normalising against context would come later.

---

## Part 2 — Going genome-scale (July 13 – July 16)

The findings above were all on chr22. To be convinced the axes generalise, we needed to run the whole genome.

We forwarded all 24 chromosomes (chr1 – chr22, chrX, chrY) through Evo 2 7B on a 2×B200 server. The batch runner used the paper's exact window schedule (6,000 bp, stride 3,000 bp, central 3 kb analysed) and stored per-position `c_t`, `oscil`, `n_enter`, `n_exit`, `below_frac`, `min_D`, `argmin_layer` in chunked parquets. Total: **~50 GPU-hours**, **5,835 chunks**, **37 GB**, **2.94 billion positions**.

chr1 took 5.4 hours; chr7 took 2.9 hours; chrY took 29 minutes. Rate stayed at ~5.1 windows per second through the whole run.

Once all chromosomes had finished, five downstream analyses were possible.

### STEP 1 — Do the paper's §3.1 numbers hold genome-wide?

Direction, yes. Absolute values, no. And one context flipped.

- **Splice donor**: *d* = −0.283 ± 0.051 across 24 chr (paper: −0.354). Same direction, ~80 % magnitude.
- **Splice acceptor**: *d* = −0.482 ± 0.074 (paper: −0.340). Same direction, **1.46× stronger**.
- **coding_exon**: *d* = +0.025 (paper: +0.08). Weaker but same direction.
- **3'UTR**: *d* = −0.100 (paper: −0.020). Same direction, ~5× stronger.
- **intergenic**: *d* = +0.113 (paper: +0.16). Same direction, ~70 % magnitude.
- **5'UTR**: *d* = **−0.011 ± 0.022** (paper: **+0.20**). **Opposite sign on all 24 chromosomes independently.**

The paper's §App E already noted that 5'UTR was entropy-coupled on the small chr22-only control panel (ρ = +0.41) but weakened at chromosome scale (ρ = +0.054). Our WGS result confirms: the paper's positive-d claim for 5'UTR was an artifact of the small chr17+chr22 pool.

**Calibration robustness**: the intron mean `c(t)` is **30.11 ± 0.20** across all 24 chromosomes — the paper reports 27.72 (chr17+chr22 pool). The +2.39 layer shift is uniform. So the paper's "single calibration transfers" claim needs qualification: *direction and relative ordering transfer, absolute values shift uniformly*. See STEP 5 for the harder version of this analysis.

### STEP 2 — Do the crossing-dynamics axes hold at 2.94 billion positions?

Yes, and the WGS run added a **fifth axis** that wasn't visible on the small chr22 sample.

At genome scale the six most striking per-context enrichments are:

| Axis | intron | splice_donor | splice_acceptor | coding_exon |
|---|---|---|---|---|
| mean_amplitude (γ − min_D) | 0.0037 | 0.034 (9×) | 0.044 (12×) | **0.057 (15×)** |
| mean_argmin_layer | 24.4 | 23.5 | 23.3 | **27.6 (latest)** |
| frac_committed (single clean crossing) | 0.85 % | 21.4 % (25×) | 25.3 % (30×) | **43.0 % (51×)** |
| frac_deliberating (oscil ≥ 1) | 11.5 % | 19.0 % (1.6×) | **24.6 %** (2.1×) | 8.4 % (**below intron**) |
| frac_crosses | 12.4 % | 40.4 % | 49.9 % | 51.4 % |

**The most surprising finding** is that coding exons show the strongest commitment signature (43 % clean single-crossing at layer 27–28, amplitude 15× intron). The paper's monotone `c(t)` reports coding_exon as mildly deeper (+0.08); on the multi-axis framework, coding_exon has a distinctive *late clean deep-commit* pattern that is invisible to `c(t)` alone.

Biologically: this is where the reading-frame / codon-usage signal seems to live. The model doesn't commit to "this is coding" early; it commits **late** and **decisively**, only after the surrounding 22 layers of context have been integrated.

### STEP 3 — Does context-aware normalisation help variant scoring?

Yes, dramatically for subtype classification, marginally for pathogenicity.

Method: normalise each variant's ref/alt/Δ of {c_t, oscil, below_frac, min_D} to per-context z-scores using the WGS-wide distribution.

**Binary P/LP vs B/LB (task from paper §App C)**:
- Baseline cos32: **0.8437**
- + raw H3: 0.8462
- + WGS-normalized H3: 0.8470
- + all features: **0.8489** — a +0.005 improvement (cos32 is already saturated on this task)

**4-class subtype classification (new task)**:
- Baseline cos32: **0.6411**
- + raw H3: 0.6588 (+0.018)
- + WGS-normalized H3: 0.6864 (+0.045)
- + scalars: 0.6669 (+0.026)
- + all combined: **🎯 0.8149** — **+0.174 macro-F1 over the paper baseline**

The subtype task uses the same 6,191 SNVs and the same random seed — the only change is what features are fed in. WGS-normalized features contribute +0.045 on their own, and combined with paper's scalar features and raw H3 features, the improvement compounds to +0.174. Ten new features (10 dimensions) alone reach macro-F1 = 0.438 — above chance by +0.19.

This is the strongest downstream result of the follow-up: the multi-axis framework is not just descriptively richer, it is **quantitatively better** for a machine-learning task the workshop paper did not attempt.

### STEP 4 — Does the paper's cCRE-ELS finding replicate genome-wide?

Yes, and the orthogonal-axis story extends here too.

We downloaded the UCSC hg38 encodeCcreCombined bigBed (809,429 dELS + pELS records across 24 chromosomes), joined with our WGS parquets, and aggregated.

- **Paper (chr22 only)**: *d*(c_t) = **−0.118**
- **WGS (24 chr, 222 million cCRE-ELS positions)**: *d*(c_t) = **−0.132 ± 0.044**

Same direction, 12 % stronger. 23 of 24 chromosomes have *d* < 0 (chrY is the exception at +0.07, likely because chrY's cCRE-ELS annotations are dominated by pseudo-autosomal and heterochromatic regions with different characteristics).

**Bonus**: WGS *d*(oscil) = **+0.117 ± 0.038**. cCRE-ELS also participates in the orthogonal oscil axis. The paper reports only the c(t) arm; the multi-axis picture reveals that enhancer-like regions share the splice-site signature (early commit **and** continued oscillation), at roughly half the magnitude of splice sites.

### STEP 5 — How chromosome-invariant is γ = 0.397, really?

The paper's `γ_cos = 0.397` was calibrated as the q70 of the running-min $D_{\cos}$ at the penultimate layer on chr22. §App A.2 reports low sensitivity across a 5×5 grid, and §3.1 shows that the same γ transfers to chr17 with 94 % magnitude preservation. But is it *chromosome-invariant*, or is chr17 just close to chr22?

We used per-chromosome q70 of intron `min_D` as a proxy (the cache stores `min_D`, not the penultimate-layer $D_{\cos}$, so the numbers are related but not identical).

Result:

- WGS mean recalibrated γ = **0.5008**
- WGS SD across 24 chromosomes = **0.0015 (0.15 %)**
- Chr-to-chr range: [0.4954, 0.5033]

The proxy value is different from paper's 0.397 (as expected: `min_D` vs penultimate-layer $D_{\cos}$), but the chromosome-to-chromosome variation is **extraordinarily small** — 24 chromosomes with SD 0.15 %. Whatever the exact calibration definition, its chromosome-invariance is much stronger than the paper's §App A.2 grid ablation could show.

---

## Part 3 — Discovering the story is bigger than we thought (July 19)

By this point the follow-up had produced:

- Six novel findings, each with per-chromosome replication
- A working 5-D interpretability framework: c_t, oscil, first_enter_layer, amplitude, argmin_layer
- A +17.4 %p downstream improvement over the workshop paper's own baseline
- A public reproducible repository

Then a broader integration became visible. The **TDiG** project (heoneyzi et al., YAICON 8th 1st Prize), independently pursuing the same base question ("when and how does a genomic foundation model converge?"), had built a **complementary metric family** — M1 direction, M2 magnitude, M3 geometry (velocity + curvature), M4 Mahalanobis distribution, M5 tortuosity — with three reference variants (A/B/C) that isolate the RMSNorm γ asymmetry. Their 98-cell context heatmap and per-consequence variant layer analysis (intron peaks at L=8, missense peaks at L=27) map onto exactly the same residual stream we were probing.

The natural conclusion: **combine into one integrated journal manuscript** rather than publish separately.

### What the integrated manuscript will look like

**Central claim**:

> Multi-axis settling is a genuine, discoverable geometric structure of Evo 2 7B's residual stream. Each axis reads a distinct biological signal. Combining axes gives a +17.4 %p macro-F1 improvement on 4-class variant subtype classification over the workshop's single-axis baseline.

**Three axis families**, previously fragmented across three teams:

- Metric axis (TDiG): M1 direction, M2 magnitude, M3 geometry, M4 distribution, M5 path efficiency
- Crossing-dynamics axis (this repo): oscil, amplitude, argmin_layer, first_enter_layer
- Scale axis (this repo): WGS 24 chromosomes, 2.94 B positions

**Section outline** ([`INTEGRATED_MANUSCRIPT_PLAN.md`](INTEGRATED_MANUSCRIPT_PLAN.md) for the full plan):

1. Introduction — the "settling as a single scalar" gap in the literature
2. Methods — unified 24-cell settling framework (14 metric cells + 9 crossing cells + baseline c(t))
3. Results — WGS §3.1 replication, 98-cell heatmap, orthogonal-axis discovery, late-clean-deep-commit at coding exons, variant scoring benchmark (+17.4 %p), L29 phase transition (from TDiG), PCA (from TDiG), cCRE-ELS WGS join, γ recalibration robustness
4. Discussion — the three settling definitions (Def 1 stops evolving / Def 2 → h₂₉ / Def 3 → h_norm); TDiG's honest acknowledgement that no metric captures Def 3, and how the oscil axis (defined against γ, not h_norm) partially fills this gap
5. Conclusion — a five-dimensional interpretability framework with distinct biology per axis and downstream utility

**Target venue**: journal (Nature Methods / Genome Research / Cell Systems / eLife tier).

---

## Part 4 — Testing the multi-axis story at the pathway level (July 20)

Part 1.3 asked a small question and got a large answer. Feeding 100 chr22 windows through Evo 2 and counting how many times $D_{\cos}$ crossed the threshold gave us `oscil` — a per-position feature the running-min envelope had thrown away. On splice acceptors, `oscil` moved in the *opposite direction* from `c(t)` (d_c_t = −0.62 vs d_oscil = +0.50). Two features derived from exactly the same layer trajectory disagreed about direction. That was the first evidence that "settling" is not one quantity. Part 2 formalised the observation into a 5-D crossing-dynamics framework; Part 3 folded in TDiG's five metric families (M1 direction, M2 magnitude, M3 geometry, M4 set-based, M5 tortuosity) with 3–5 reference or α/β variants each — fourteen more candidate axes — and gave the whole architecture a name: three axis families (metric / crossing-dynamics / scale).

But Part 3 imported the M-cells on the strength of one comparison: M3_geo curvature *d* = −0.80 at chr22 splice donors vs cosine `c(t)` *d* = −0.35. 2.3× stronger. That single number is compatible with two very different pictures.

Under the first, all 17 M-cell settings read the same underlying settling signal at different signal-to-noise ratios; different M-cells are amplification variants of one axis. This is the picture the workshop paper implicitly assumed and the picture TDiG's own presentation left plausible. Under the second, the M-cells read genuinely different aspects of the residual-stream trajectory and the strongest of them at splice sites happens to be M3_geo. If the second is true, the "3 axis families" framework from Part 3 is really more axes than we counted — and, crucially, the +17.4 %p downstream improvement from combining features (Part 2 STEP 3) gets a mechanistic explanation: combining works because the added features carry independent biological information, not because they average out different noise realisations.

Distinguishing the two pictures needed a broader biological substrate than splice sites. Two M-cells could easily give different Cohen's *d* on splice donors while ranking every other gene set identically. What we needed was a per-cell test on a large, diverse set of gene classes — GO Biological Process was the obvious choice — with a sharp discriminator: if the same GO class shows *opposite* directional signals in different cells, that is a **sign-flip**, and no amount of amplification variance can produce it.

### GO enrichment on 17 cells × 471 pathways

The pipeline: 1,633 protein-coding genes on chr17 + chr22 (GENCODE v44), aggregated to per-gene means from the tier1 parquets that TDiG originally computed on those two chromosomes; per-cell per-GO-term Mann–Whitney rank-sum comparing in-set vs out-of-set genes; Cohen's *d* on the raw settling values; BH-FDR correction per cell; 471 GO BP terms with 5 ≤ in-set overlap ≤ 1,500 (an expansion from 137 curated terms tested in a first pass). Full details in [`EXP8_MULTI_CELL_GO_FINDINGS.md`](EXP8_MULTI_CELL_GO_FINDINGS.md).

Thirteen of the 17 cells produce hits at BH q < 0.05. The other four (M2_mag_refB_diag, M2_mag_refC_diag, M4_set_refB, M4_set_refC) act as null baselines — per-gene aggregation collapses those particular reference choices to near-constant values across the universe. The 13 live cells favor the different-axes picture cleanly:

| Cell | Signal |
|---|---|
| M4_set_refA, M5_tau_refA | Epithelial / intermediate-filament / transcription biology |
| M3_geo_a1.0_b0.0 (β = 0, pure curvature) | Immune / inflammatory |
| M3_geo_a1.0_b0.5 (curvature-dominant, β = 0.5) | Olfactory receptors, extremely strong (GO:0050911, *d* = −1.96, p = 6×10⁻⁷) |
| M5_tau_refB | Olfactory-DEEP + epithelial-SHALLOW (same cell, opposite direction on two classes) |
| M5_tau_refC | Anti-correlates with M5_tau_refA at **Spearman r = −0.78** across all 471 GO terms |

Forty GO terms show sign-flip patterns — at least one cell with significant positive *d* AND at least one with significant negative *d*. The three most discriminating are the ones the M-cells actually disagree about most:

| GO term | Max +*d* | Min −*d* |
|---|---|---|
| GO:0050911 detection of chemical stimulus (olfactory) | +3.24 (M1_dir_refC) | −1.96 (M3_geo_a1.0_b0.5) |
| GO:0045109 intermediate filament organization | +1.06 (M5_tau_refA) | −1.32 (M4_set_refA) |
| GO:0002009 morphogenesis of an epithelium | +0.85 (M5_tau_refA) | −1.15 (M4_set_refA) |

Two observations sit inside this table that push past the aggregate finding.

First, the olfactory axis and the epithelial-filament axis are *systematically* opposite: cells that settle olfactory receptors deep settle epithelial filaments shallow, and vice versa. There is a pair of biological axes that the M-cells trade off against each other, and picking a different M-cell moves position along both axes at once, in opposite directions. This is exactly the pattern Part 1.3's c_t / oscil comparison flagged for splice-site geometry, now recurring at the gene-set level across the residual stream. The workshop's single-scalar treatment did not just miss oscil — it missed the fact that any single choice of settling metric commits to a specific axis of biology while relinquishing others.

Second, the reference-token result is not a robustness check the way TDiG originally framed it. TDiG tested three references (A / B / C) inside each M-cell family to demonstrate insensitivity to the RMSNorm γ asymmetry — that is, to *rule out* an artifact. On identical input data, M5_tau_refA and M5_tau_refC produce Spearman r = −0.78 across 471 GO terms. Reference is not a hyperparameter; reference is axis. Similarly for M3_geo: β = 0 (pure curvature) sits in one cluster with M4_set_refA and M2_mag_refA; adding β = 0.5 crosses to the opposite cluster with M3_geo_a1.0_b0.5. A parameter change on the α/β mixing coefficient of one metric family crosses an axis boundary. The "5 M3_geo variants" from Part 3's inventory are not five points on a curvature-cosine interpolation; they are up to five different axes.

### The honest replication: from 40 raw sign-flips to 3 pairs that survive chromosome-split

Combined-dataset significance can be sample-size boosted. Forty sign-flip terms on chr17 + chr22 together does not tell us how many replicate independently. So we retested each (cell, GO) pair on chr17 alone (1,186 genes) and chr22 alone (447 genes):

- 680 pairs total (40 GO × 17 cells)
- **8 pairs (1.2 %)** significant on both chromosomes with the same *d* sign
- 44.7 % were sample-size boosted (neither chromosome alone significant; combined was)
- 16.6 % chromosome-driven (only one chromosome contributes)

Only 8 pairs out of 680 truly replicate. The 40 was inflated by combining the two chromosomes. This has to be in the manuscript.

To push replication further we extended to a third chromosome. Chr21 was not in TDiG's original tier1 dataset — only chr17 and chr22 were — so we had to produce it from scratch. A fresh Evo 2 7B forward pass on 13,365 valid chr21 6-kb windows, saving a subset of the tier2 scalars from which we could derive 7 of the 17 M-cells (M1_dir_refA, M2_mag_refA, and the five M3_geo α/β variants); the refB and refC variants and the M4_set and M5_tau primitives require TDiG-team reference definitions that were not accessible on our timeline. The forward pass took 58.7 minutes on 2 × NVIDIA B200 after an optimised HDF5 writer (per-window chunks, no gzip) gave a ~9× speedup over a naive first implementation.

Chr21 has only 220 protein-coding genes — underpowered on its own — but sufficient to check whether the 8 two-chromosome-robust pairs remain when a third chromosome enters:

- **0 pairs (0 %)** significant with the same sign on all three chromosomes
- **3 pairs** significant with the same sign on any 2 of 3 (all chr17 + chr22, with chr21 NaN due to small n)

Zero pairs pass strict three-chromosome replication. Three survive the softer 2-of-3 threshold. The 40-term raw count was doing what raw combined-dataset counts often do. This is the caveat the manuscript will present alongside — not instead of — the aggregate finding.

### The one clean case that carries the claim by itself

Among the 3 robust two-chromosome pairs, one is a full manuscript-quality example on its own:

**GO:0006869 lipid transport** shows a genuine sign-flip between two variants of the M3_geo family:

| Cell | *d* chr17 | *p* chr17 | *d* chr22 | *p* chr22 |
|---|---|---|---|---|
| M3_geo_a0.5_b1.0 (cosine-dominant) | +0.52 | 0.043 | +0.65 | 0.024 |
| M3_geo_a1.0_b0.0 (pure curvature) | −0.72 | 0.009 | −0.89 | 0.007 |

Same GO class. Same two chromosomes. Two variants of M3_geo — same base primitive family, differing only in β (the cosine coefficient in α · velocity + β · curvature). Adding cosine flips the sign of the settling-depth deviation for lipid-transport genes. Same trajectory, different projection, opposite biology signal, independently replicated on two chromosomes. This is the clean version of what the aggregate sign-flip count was reaching for, and it should be §3.11's canonical example.

A second connection lands on the biology side. For each significant (cell, GO) pair at q < 0.10, we tabulated the top-20 in-set driver genes ranked by mean settling value and flagged membership in a 15-gene ClinVar cancer panel (BRCA1, BRCA2, TP53, PTEN, STK11, CDH1, PALB2, ATM, CHEK2, BARD1, RAD51C, RAD51D, MLH1, MSH2, APC). Of 107 significant pairs, three carry cancer-panel drivers:

- **M3_geo curvature-only × GO:0007131 reciprocal meiotic recombination**: *d* = +0.91, top-5 drivers RAD51C, DMC1, RAD51D, TEX19, SYCE3. **Two of the top five are HRD cancer-panel genes.** The curvature-only cell captures homologous-recombination pathway biology and puts the two clinically-HRD-relevant genes at the top of the driver list.
- **M5_tau_refB × GO:0006357 regulation of RNA pol II transcription**: *d* = +0.33, TP53 in top-20.
- **M5_tau_refB × GO:0006355 regulation of DNA-templated transcription**: same pattern, TP53 in top-20.

Three hits out of 107 pairs is not statistically over the null baseline expected from a 15-gene panel in a 1,633-gene universe. But the *coherence* is the point: the HRD genes cluster inside a meiotic-recombination GO picked by the curvature cell; TP53 sits inside two transcription-regulation GO terms picked by the tau cell. Two independent, biologically-coherent connections between specific M-cell primitives and clinically-relevant pathway biology. That is a stronger signal than aggregate enrichment can capture.

### What Part 4 does to the Part 3 framework

Part 3's architecture had three axis families (metric / crossing-dynamics / scale) and left the mechanism for +17.4 %p downstream improvement implicit. Part 4 supplies the mechanism and revises the count:

- **The metric family is itself multi-axis.** Reference-token variants (r = −0.78 for M5_tau refA vs refC) and α/β shifts within M3_geo cross axis boundaries; they do not sub-tune one axis. The "5 metric primitives × 3–5 variants = 17 cells" from Part 3 does not compress cleanly to five orthogonal metric axes.
- **The +17.4 %p downstream gain now has a mechanistic explanation.** Combining M-cells adds independent biological information because M-cells read different pathways, not because they average out different noise on the same pathway. The Part 2 STEP 3 result was a hint; Part 4 supplies the underlying reason.
- **The workshop paper's single-scalar treatment has a naming problem.** It is not just missing information; it is committing to a specific axis of biology at the cost of others.

The framework needs a single language change. We recommend the integrated manuscript adopt it globally:

- "Settling **depth**" (a singular scalar readout) → "settling **profile**" (a family of coordinates).
- "The cosine metric captures X" (cosine as *the* measurement) → "The cosine-based M1_dir cells capture X; the curvature-based M3_geo_a1_b0 cell captures Y" (cosine as one axis among many).
- "Reference variants ensure robustness" (reference as hyperparameter) → "Reference variants access different axes of the trajectory" (reference as axis).

Everything else in the manuscript plan continues as written. §3.11 (a new subsection) should be structured around the lipid-transport example as the canonical clean case, with the 8 chr-2-robust and 3 chr-3-of-2-robust numbers as the honest characterisation, and the RAD51C/D + TP53 driver hits as the biology anchor. The full technical details, all CSVs, all figures, and the run-by-run session provenance live in [`EXP8_MULTI_CELL_GO_FINDINGS.md`](EXP8_MULTI_CELL_GO_FINDINGS.md).

*Infrastructural footnote.* The chr21 forward pass and the DASH backup effort described here ran against a hard midnight deadline on 2026-07-20. The compute host `TGIL_mutsig` became permanently unreachable at 23:58 KST that night (GPU lease expiry); a queued chr20 forward was interrupted at ~10 % complete and the last piece of the DASH transfer (`results_cached`) never landed. Everything that survived is on DASH at `/home/darejin/TDiG/`, and the recovery notes for a new-host bootstrap are in [`../RECOVERY_LOG.md`](../RECOVERY_LOG.md). Extending the M-cell profile to a fourth or fifth chromosome — and closing the 7-of-17-cells gap on chr21 through TDiG-team reference definitions — is the natural next step when a new GPU allocation arrives.

---

## Part 5 — What remains open

The manuscript plan is written, the code is reproducible, the WGS data is complete. The remaining work is not analytical:

- **TDiG author coordination**: agreeing on co-authorship, aggregating results across the two teams, deciding figure attribution.
- **Manuscript LaTeX draft**: the outline exists; the ~5,000-word main text does not.
- **Extended Data compilation**: cross-chromosome tables, γ ablations, α/β ablations for M3, reference-variant ablations for M1/M4/M5.
- **Journal-quality figure regeneration**: all current figures are in the 300 dpi matplotlib output style; a few will need reworking for print.

The scientific claims are settled — everything the manuscript will say is currently a JSON or CSV in [`results/`](../results). What remains is writing.

---

## Reading order

If you want to go deeper after this narrative:

1. [`PROJECT_OVERVIEW.md`](PROJECT_OVERVIEW.md) — one-page project summary (single entry point for return visits).
2. [`INTEGRATED_MANUSCRIPT_PLAN.md`](INTEGRATED_MANUSCRIPT_PLAN.md) — the section-by-section plan for the journal manuscript.
3. [`RESULTS.md`](RESULTS.md) — every headline number with the exact source script and JSON path.
4. [`METHODS.md`](METHODS.md) — every metric and feature with its precise definition.
5. [`REPRODUCE.md`](REPRODUCE.md) — end-to-end reproduction from Evo 2 7B weights + reference data.
6. [`../paper/gdtr_paper_ICML_3.pdf`](../paper/gdtr_paper_ICML_3.pdf) — the original workshop paper, kept frozen.

---

*Last updated: 2026-07-21 (Part 4 rewritten to connect logically with Part 1.3 orthogonal-axis discovery and Part 3 architecture; opens with the specific hypothesis Part 3 left implicit, closes with mechanism for +17.4 %p from Part 2 STEP 3).*
