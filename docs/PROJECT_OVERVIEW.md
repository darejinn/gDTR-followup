# PROJECT OVERVIEW — gDTR integrated manuscript

**One-page summary of the entire project as of 2026-07-19.**

Use this as the single entry point when returning to the project. Every other document expands one section of this overview.

---

## 1. What this repository is

The working home of an **integrated journal manuscript** on gDTR (Genomic Deep-Thinking Ratio), currently in preparation. The manuscript will supersede a workshop paper that has been accepted but not yet published in main proceedings.

Repository URL: <https://github.com/darejinn/gDTR-followup>

**Not a standalone follow-up.** Not a fork. The base paper is being extended into a journal-length manuscript that combines three bodies of work into one coherent story.

## 2. The three bodies of work being merged

| Source | Contribution | Status |
|---|---|---|
| **gDTR workshop paper** (Cho, Kang, Park, Kim) — accepted at ICML 2026 GenBio Workshop | Base method: `c(t)` running-min settling depth, γ_cos = 0.397, chr17+chr22 replication, splice/enhancer shallowness, motif ↔ flank bidirectionality, ClinVar variant scoring 32-d AUROC 0.844 | Frozen PDF at `paper/gdtr_paper_ICML_3.pdf`. Will be superseded by the integrated manuscript. |
| **TDiG** (heoneyzi et al., YAICON 8th 1st Prize, [heoneyzi/TDiG](https://github.com/heoneyzi/TDiG)) | Extension: M1–M5 metric family (direction / magnitude / geometry / distribution / path-efficiency), 3 reference variants (A/B/C), 98-cell context heatmap, per-consequence variant analysis (intron L=8, missense L=27), PCA (PC1 = 75.7 % variance), L29 phase transition, γ ablation | Available at TDiG repo. Awaiting co-authorship and data-access coordination before merging. |
| **followup analyses** (this repo, July 2026) | Extension: WGS 24-chromosome replication (2.94 B positions), threshold-crossing dynamics (H3a–H3d), variant subtype benchmark macro-F1 **0.641 → 0.815** (+17.4 %p), cCRE-ELS genome-wide, γ recalibration SD 0.15 % | Complete. All scripts, results, and documentation in this repo. |

## 3. The central claim of the integrated manuscript

> Multi-axis settling is a genuine, discoverable geometric structure of Evo 2 7B's residual stream. Each axis reads a distinct biological signal. Combining axes gives a **+17.4 %p macro-F1 improvement** on 4-class variant subtype classification over the workshop's single-axis 32-d ΔD_cos baseline.

Three axis families in the framework:

- **Metric axis** (from TDiG): direction, magnitude, geometry, distribution, path-efficiency (M1–M5)
- **Crossing-dynamics axis** (from followup): oscil, amplitude, argmin_layer, first_enter_layer
- **Scale axis** (from followup): whole genome — 24 chromosomes, 2.94 B positions

## 4. Ten headline results (with attribution)

| # | Question | Result | Source |
|---|---|---|---|
| 1 | Do §3.1 direction claims hold at WGS scale? | 6 of 7 contexts ✓; **5'UTR flips sign uniformly** on 24 chr | followup STEP 1 |
| 2 | Is `c(t)` the only axis? | 5-D framework (c_t, oscil, first_enter, amplitude, argmin_layer) | followup H3d |
| 3 | c_t vs oscil — same signal or orthogonal? | **Orthogonal**: splice acceptor d_c_t = −0.62, d_oscil = +0.50 (opposite sign) | followup H3a |
| 4 | Does adding axes help downstream ML? | Subtype macro-F1 **0.641 → 0.815** (+17.4 %p) | followup STEP 3 |
| 5 | Does cCRE-ELS replicate genome-wide? | d_c_t = −0.132 ± 0.044 (paper chr22-only −0.118); **new d_oscil = +0.117** | followup STEP 4 |
| 6 | How chromosome-invariant is γ = 0.397? | q70 SD = **0.15 %** across 24 chr | followup STEP 5 |
| 7 | Best single variant feature? | ΔH_norm_L1 at L=8 — **AUROC 0.856** | TDiG §A1 (to be merged) |
| 8 | Does variant peak layer index disrupted information level? | intron/3'UTR peak L=8, missense/synonymous peak L=27 | TDiG §A2 (to be merged) |
| 9 | What happens at L29? | Context probing crashes 0.980→0.799; variant AUROC crashes 0.85→0.79 | TDiG §B1 (to be merged) |
| 10 | Low-dim structure behind M1–M5? | PC1 = 75.7 % of variance; single geometric axis encodes bidirectional settling | TDiG §C (to be merged) |

Rows 1–6 are executed and reproducible from this repo today. Rows 7–10 come from TDiG and are staged pending coordination.

## 5. Authorship (current understanding)

Not finalised. Working roster:

- **First authors** (co-first): Yoonjin Cho, [TDiG lead — heoneyzi]
- **Middle authors**: Jiheon Kang, Subin Park, [TDiG YAICON 8th team members]
- **Senior / corresponding**: Sangwoo Kim (Yonsei University College of Medicine)

Author order and contribution attribution to be confirmed by all parties before manuscript submission.

## 6. Target venue

Journal (Nature Methods / Genome Research / Cell Systems / eLife / Bioinformatics tier). Preference order captured in [`docs/INTEGRATED_MANUSCRIPT_PLAN.md`](INTEGRATED_MANUSCRIPT_PLAN.md) §4.

Expected length: ~5,000 words main text + unlimited Extended Data / Supplementary.

## 7. Repository layout (what to look at first)

```
gDTR-followup/
├── README.md                              # short introduction + headline results table
├── docs/
│   ├── PROJECT_OVERVIEW.md                # ← this file (single-page entry point)
│   ├── INTEGRATED_MANUSCRIPT_PLAN.md      # section-by-section outline of the journal manuscript
│   ├── RESULTS.md                         # every headline number with source script + JSON path
│   ├── METHODS.md                         # every feature and metric definition
│   ├── REPRODUCE.md                       # end-to-end reproduction recipe
│   └── CHANGELOG.md                       # dated log
├── src/                                   # frozen gDTR framework (17 modules from workshop paper)
├── scripts/                               # analysis pipelines
│   ├── wgs/                               # STEP 1/2/4/5 + forward + labels + batch
│   ├── exp2/                              # H2a / H2b / integrated / WGS-normalized
│   └── exp3/                              # H3a / H3b / H3c / H3d
├── tdig_integration/                      # staging area for TDiG artifacts (currently empty by design)
├── results/
│   ├── genome_summary/                    # STEP 1/2/4/5 outputs (25 files, 1.5 MB)
│   ├── exp2/                              # H2 JSON + CSV + figures
│   └── exp3/                              # H3 JSON + CSV + figures
├── paper/                                 # frozen workshop paper artifacts
└── data/                                  # data version pointers + requirements
```

Not committed (regenerable, see [`docs/REPRODUCE.md`](REPRODUCE.md)):
- 37 GB per-chromosome per-position parquets (5,835 chunks)
- 13 GB Evo 2 7B weights
- 4 GB per-chromosome position labels
- ClinVar VCF, GRCh38 fa, GENCODE v44 GTF, UCSC cCRE bigBed

## 8. Reproducibility

**Everything committed to this repo can be regenerated from Evo 2 7B weights + published reference data** by following [`docs/REPRODUCE.md`](REPRODUCE.md).

- WGS forward pass: **~50 GPU-hours** on 2 × B200 (or ~140 GPU-hours on 1 × H200)
- Downstream analyses (STEP 1 through 5, EXP2, EXP3 except forward passes): **~20 min on CPU** after WGS parquets exist

Locked parameters (do not modify):
- Evo 2 7B weights: `arcinstitute/evo2_7b`, SHA `bda0089f92582d5baabf0f22d9fc85f3588f6b58`, MD5 `359ef88ccac2a62644035578de8a7db4`
- γ_cos = 0.397, L* = 29, window / stride = 6,000 / 3,000 bp
- Random seed = 42 for every stochastic step
- Data versions: GENCODE v44, ClinVar 2026-04-18, GRCh38 UCSC, ENCODE SCREEN v3/v4

## 9. What is finished vs pending

### Finished
- Workshop paper submitted and accepted (frozen at `paper/`)
- WGS batch runner across all 24 chromosomes (completed 2026-07-16)
- 5-STEP downstream analysis (this repo)
- EXP2 (H2a, H2b, integrated H2b, WGS-normalized H2c)
- EXP3 (H3a, H3b, H3c, H3d)
- All documentation stack (README, RESULTS, METHODS, REPRODUCE, CHANGELOG)
- Integrated-manuscript plan document

### Pending
- **Co-author coordination with TDiG team** (heoneyzi et al., YAICON 8th) — must precede any TDiG-artifact import
- **TDiG data access** (chr22 raw tier-1/2/3 caches, ~120 GB, on TDiG's HF space)
- **Manuscript LaTeX draft** — the plan is written; the manuscript itself is not
- **Extended Data / Supplementary tables** compilation
- **Journal-quality figure regeneration**
- **Submission** to journal

### Deferred (future work, discussed in §4 of the manuscript plan)
- M1–M5 extension to WGS (requires raw h_ell at 24-chr scale — additional 200 GB storage + ~40 GPU-hours)
- Caduceus and other per-bp architecture replication
- Causal / SAE intervention grounding
- Non-cancer / whole-exome variant cohort

## 10. Where to look next depending on your goal

| I want to … | Read … |
|---|---|
| Understand what claim the paper makes | `docs/INTEGRATED_MANUSCRIPT_PLAN.md` §1–2 |
| See every headline number with source | `docs/RESULTS.md` |
| Understand every metric's precise definition | `docs/METHODS.md` |
| Reproduce a specific number | `docs/REPRODUCE.md` |
| Learn what changed and when | `docs/CHANGELOG.md` |
| Contribute a TDiG analysis | `tdig_integration/README.md` |
| Read the original workshop paper | `paper/gdtr_paper_ICML_3.pdf` |
| See the section-by-section manuscript outline | `docs/INTEGRATED_MANUSCRIPT_PLAN.md` §2–3 |

## 11. Contact

Corresponding author: **Sangwoo Kim** (`swkim@yuhs.ac`)
Yonsei University College of Medicine, Department of Biomedical Systems Informatics
Seoul, Korea

---

*Last updated: 2026-07-19*
