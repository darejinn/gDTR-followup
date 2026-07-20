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

## Part 4 — Extending the M-cell family to a third chromosome and the GO story (July 20)

The July-19 integration plan noted an outstanding item: *"M1–M5 metric family not yet extended to WGS (would require raw h_ell across 24 chromosomes — deferred as future work in §4.4 of the manuscript plan)."* This section reports one substantial step toward closing that gap and the biological pathway story it revealed.

### 4.1 — Motivation: does the M-cell family read *different biology*, not just *stronger signals*?

TDiG's headline claim about M1–M5 was that different metric primitives give quantitatively different Cohen's *d* at splice sites (M3_geo `curvature`, `d = -0.80`, vs cosine `c(t)` `d = -0.35` on chr22). That answers whether metrics differ in *strength*. It does not answer whether they differ in *what they read*.

Two hypotheses, both plausible in advance:

- **Same axis, different strengths.** All 17 M-cell settings (M1 direction × 3 refs, M2 magnitude × 3 refs, M3 geometry × 5 α/β configs, M4 set-based × 3 refs, M5 tau × 3 refs) read the same underlying "settling" signal at different amplitudes.
- **Distinct axes.** The 17 settings read partially or fully orthogonal biological signals — some settle olfactory genes deep, some settle epithelial genes deep, some settle immune genes deep, and the identity of which biology gets settled depends on which primitive is used.

We picked GO Biological Process as the discriminating substrate. If hypothesis 1 is right, all 17 cells should rank the same GO terms first (with different *d* magnitudes). If hypothesis 2 is right, top hits should differ by cell — and some GO terms may have *opposite* directional signals in different cells (a **sign-flip**).

### 4.2 — EXP8 pipeline: 17 M-cells × per-gene mean × GO rank-sum

Setup:

- Universe: 1,633 protein-coding genes on chr17 + chr22 (from GENCODE v44).
- Feature: per-gene mean of the per-position M-cell settling value (aggregated within gene bodies from the chr17 and chr22 tier1 parquets provided by the TDiG team).
- Test: per-cell, per-GO-term Mann–Whitney rank-sum comparing in-set vs out-of-set genes; Cohen's *d* on the raw settling values; BH-FDR correction per cell.
- Scope: 471 GO BP terms with 5 ≤ in-set overlap ≤ 1,500 (expanded from 137 curated terms).

Headline numbers, per cell (BH q < 0.05 hit counts):

| Cell | q<0.05 hits | Top axis |
|---|---|---|
| M4_set_refA | 13 | epithelial / transcription |
| M5_tau_refA | 12 | epithelial / ubiquitin |
| M2_mag_refA | 11 | immune / migration |
| M3_geo_a1.0_b0.0 (curvature-only) | 9 | **immune / inflammatory** |
| M5_tau_refB | 9 | **olfactory-DEEP + epithelial-shallow** |
| M3_geo_a1.0_b0.5 (curvature-dominant) | 8 | **olfactory-DEEP** |
| M3_geo_a1.0_b1.0 (equal weight) | 8 | virus response / transcription |
| M3_geo_a0.5_b1.0 (cosine-dominant) | 8 | immune / transcription |

Four cells are effectively *dead* (0 rank-sum hits): M2_mag_refB_diag, M2_mag_refC_diag, M4_set_refB, M4_set_refC. Per-gene aggregation of those primitives with `refB/refC` reference collapses to near-constant values across all 1,633 genes.

Immediately, the top-axis column looks like hypothesis 2: cells split cleanly into immune-dominant, epithelial-dominant, olfactory-dominant, and transcription-dominant clusters — not a single strength ranking.

### 4.3 — 40 sign-flip GO terms

We tabulated, per GO term, how many cells give significantly positive *d* (≥ +0.3) and how many give significantly negative *d* (≤ −0.3). Forty terms have at least one positive AND one negative significant cell — a **sign-flip**: the *same* GO class is encoded with opposite direction by different cells.

The strongest sign-flip GO terms by variance of the signed −log₁₀(p) across cells:

| GO term | Variance | Range | Interpretation |
|---|---|---|---|
| GO:0050911 detection of chemical stimulus (olfactory) | 19.1 | 15.5 | The olfactory-receptor family separates cells cleanly into "settles deep" (M1_dir, M5_tau_refB) vs "settles shallow" (M3_geo β>0) groups |
| GO:0045109 intermediate filament organization | 18.4 | 17.9 | The epithelial/keratin axis flips opposite to olfactory: cells that settle olfactory shallow settle filaments deep |
| GO:0002009 morphogenesis of an epithelium | 13.8 | 15.6 | Tracks intermediate filament (same +/− cell split) |

Cluster analysis on the 13 live cells (average-linkage hierarchical clustering on Spearman similarity of signed −log₁₀(p) across GO terms) gives four clusters, of which the most striking finding is that **`M5_tau_refA` and `M5_tau_refC` correlate at −0.78** — the same underlying M5_tau family with different reference tokens produces near-perfectly opposite biology maps. Reference token choice does not tune the same axis; it selects a different axis.

Similarly for M3_geo: β = 0 (curvature-only) sits in one cluster (majority, immune-dominant), while β = 1 (adding cosine) crosses to the opposite cluster. Small parameter shifts within one metric family cross an axis-boundary.

### 4.4 — Adding chr21 as a third replicate: 8 → 3 → 0 robust pairs

The 40 sign-flip terms are the raw pattern in the combined chr17 + chr22 dataset. The next question is whether they replicate on independent chromosomes.

**Chromosome-split analysis on the existing chr17 + chr22 data (`e8_40`):**

Each of the 40 sign-flip GO terms was retested on chr17 alone (1,186 genes) and chr22 alone (447 genes). For each (cell, GO) pair — 40 GO × 17 cells = 680 tests — status was classified as:

- `robust_same_sign` (chr17 and chr22 both significant with same *d* sign): **8 pairs (1.2 %)**
- `chr17_driven` (only chr17 significant): 73
- `chr22_driven`: 40
- `ns_split` (neither chr alone significant, but combined was): 304 (44.7 %)
- `insufficient_data` (< 3 in-set genes in one chr): 255

Most sign-flips are sample-size boosted, not chromosome-independent. This is a **required caveat** for the manuscript.

**Extending to chr21 (`e8_41`, this session's core work):**

We produced a chr21 tier1 parquet from scratch. This required a Evo 2 7B forward pass over ~13,365 valid 6-kb chr21 windows (chr21 has ~15,568 total; ~2,200 are heavy-N centromeric/telomeric regions), saving per-position tier2 scalars (cos_ref_A, step_cos, step_norm, norm_h_ell, norm_h_29), then deriving 7 of the original 17 M-cells (M1_dir_refA, M2_mag_refA, M3_geo × 5 α/β) via a downstream calibration step. The forward pass took 58.7 minutes on 2 × NVIDIA B200 with an optimised HDF5 writer (removing gzip and using per-window chunks gave a **~9× speedup** over an initial naive implementation).

The full 17-cell schema was not reproduced — the reference-token definitions used by TDiG for the refB and refC variants (and the exact set definition used for M4, and the tau weighting for M5) are held in TDiG-team-only code and were not accessible on the timeline available. Our chr21 covers 7 of the 17 cells with our own explicit refA definition (h_norm-based cosine, matching the workshop paper's cosine lens).

Reruns of the 40 sign-flip tests with chr21 added:

| Status | Count |
|---|---|
| `robust_3_same_sign` (all 3 chr significant, same sign) | **0** |
| `robust_2_same_sign` (any 2 of 3 chr significant, same sign) | 3 |
| `single_chr_only` | 96 |
| `ns_all` | 181 |
| `insufficient_data` (chr21 has only 220 protein-coding genes) | rest |

**No pair replicates across all three chromosomes.** The three pairs that replicate across two of three are all chr17 + chr22 pairs (chr21 has too few genes for chr21-alone significance).

### 4.5 — The one clean example: GO:0006869 lipid transport

Among the 3 robust two-chromosome pairs, one is manuscript-worthy on its own:

**GO:0006869 lipid transport** shows a genuine sign-flip between two M3_geo variants, replicating on both chr17 and chr22:

| Cell | *d* chr17 | *p* chr17 | *d* chr22 | *p* chr22 |
|---|---|---|---|---|
| M3_geo_a0.5_b1.0 (cosine-dominant) | **+0.52** | 0.043 | **+0.65** | 0.024 |
| M3_geo_a1.0_b0.0 (pure curvature, no cosine) | **−0.72** | 0.009 | **−0.89** | 0.007 |

Same GO class. Same two chromosomes. Two metric variants that share the same base primitive family (M3_geo, α · velocity_z + β · curvature_z) differ only in the presence or absence of a cosine (β) component. That single change flips the sign of the settling-depth deviation for lipid-transport genes, robustly across two independent chromosomes.

This is the cleanest single evidence that a small change in the metric primitive — adding a cosine component to a pure-curvature composition — accesses a different directional axis of the residual stream's biological encoding. Same trajectory, different projection, opposite biology signal.

### 4.6 — Cancer-panel driver hits

For each significant (cell, GO) pair at rank-sum q < 0.10, we dumped the top-20 in-set genes ranked by per-cell settling value. Out of 107 significant pairs, three include cancer-panel genes (a 15-gene ClinVar panel: BRCA1, BRCA2, TP53, PTEN, STK11, CDH1, PALB2, ATM, CHEK2, BARD1, RAD51C, RAD51D, MLH1, MSH2, APC) in their top-20 drivers:

- **M3_geo_a1.0_b0.0 × GO:0007131 (reciprocal meiotic recombination)**: *d* = +0.91, drivers include **RAD51C and RAD51D** (both HRD-related cancer-panel genes). The M3_geo curvature-only cell picks up homologous-recombination deficiency biology at the pathway level.
- **M5_tau_refB × GO:0006357 (regulation of transcription by RNA pol II)**: *d* = +0.33, driver includes **TP53**.
- **M5_tau_refB × GO:0006355 (regulation of DNA-templated transcription)**: *d* = +0.35, driver includes TP53 (same top-5 as above).

These are 3 hits in 107 significant pairs — with a 15-gene cancer panel against a 1,633-gene universe, the null expectation would be ~1 hit. Not statistically dispositive, but suggestive that specific M-cell settings pick up clinically-relevant biology beyond the group-level GO enrichment.

### 4.7 — The concept refactor: "settling depth" → "settling profile"

Given that:

1. Reference-token choice within one metric family flips the axis (M5_tau_refA vs refC: r = −0.78);
2. Adding a β component to M3_geo crosses axis-boundaries at small α/β shifts;
3. Forty GO terms show sign-flip patterns; one (lipid transport) replicates robustly across two chromosomes with opposite direction between two closely-related M3_geo variants,

the "settling depth" concept as a single scalar readout of the trajectory is inadequate. Each M-cell primitive probes a different aspect of the trajectory, and treating them as variations of one quantity misses the primary phenomenon.

Recommended manuscript reframing:

- "Settling **depth**" (singular, treats all M-cells as measures of the same thing) → "settling **profile**" (family, treats the M-cell set as coordinates on a multi-axis readout).
- "The cosine metric captures X" (framing cosine as *the* correct measurement) → "The cosine-based M1_dir_ref* cells capture X; the curvature-based M3_geo_a1_b0 cell captures Y" (framing each M-cell as one axis among many).
- Add an explicit subsection on sign-flip robustness — presenting the 40 raw pairs, the 8 chr-2-robust pairs, the 3 chr-3-robust pairs, and the one canonical clean example (lipid transport) as the honest characterization.

### 4.8 — Infrastructural note: TGIL_mutsig shutdown, July 20 midnight

The compute host used for the chr21 forward pass (`TGIL_mutsig`, 2 × B200) became permanently unavailable at midnight KST on 2026-07-20 (GPU lease expiry). This truncated the planned chr20 forward (~10 % complete when the host died) and the full BATCH-3 transfer to DASH (`results_cached` did not transfer; `data_ref` partial 60 %). Everything critical for reproducing the EXP8 findings above **is** on DASH at `/home/darejin/TDiG/` — see [`../RECOVERY_LOG.md`](../RECOVERY_LOG.md) for the incident narrative and next-session bootstrap. Extending the M-cell analysis to a fourth or fifth chromosome will require a new GPU host and re-running the `e1_20_chr_forward_tier2.py + e1_30_tier2_to_tier1.py` pipeline.

### 4.9 — What Part 4 adds to the manuscript

The integrated manuscript now has:

- One new empirical section on 17-M-cell GO enrichment at the pathway level (data on chr17 + chr22 + chr21, 7 of 17 M-cells for chr21)
- One conceptual refactor (depth → profile)
- One canonical clean sign-flip example (lipid transport, replicating on 2 chromosomes)
- One honest caveat (most sign-flip patterns do NOT replicate across chromosomes; only 3/280 pairs robust across 2 of 3 chromosomes, 0 across all 3)
- Two cancer-panel biology anchors (RAD51C/RAD51D via M3_geo curvature; TP53 via M5_tau)

Full technical details: [`EXP8_MULTI_CELL_GO_FINDINGS.md`](EXP8_MULTI_CELL_GO_FINDINGS.md). Reproducibility notes: [`../RECOVERY_LOG.md`](../RECOVERY_LOG.md) plus the on-DASH `exp8_multi_cell_go/README.md` and `exp8_multi_cell_go/NOTES_ON_REPRODUCIBILITY.md`.

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

*Last updated: 2026-07-21 (Part 4 added — 17-M-cell GO enrichment on chr17+22+21, sign-flip robustness, lipid transport clean example, concept refactor).*
