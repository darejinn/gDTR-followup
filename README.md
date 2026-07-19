# gDTR — Multi-Axis Settling in Genomic Foundation Models

**Repository home of the integrated gDTR journal manuscript** (in preparation).

This repository is the working home for the **extended, integrated version** of gDTR, currently being drafted for journal submission. It merges:

1. **gDTR original workshop paper** (Cho, Kang, Park, Kim — ICML 2026 GenBio Workshop, accepted). Base: layer-wise settling depth `c(t)` on Evo 2 7B.
2. **TDiG** (heoneyzi et al., YAICON 8th 1st Prize, [heoneyzi/TDiG](https://github.com/heoneyzi/TDiG)). Extension: M1–M5 multi-metric family, 3-reference variants, 98-cell context heatmap, per-consequence variant analysis, PCA interpretability, L29 phase transition.
3. **gDTR-followup analyses** (this repo's Jul 2026 sessions). Extension: WGS 24-chromosome replication (2.94 B positions), threshold-crossing dynamics (H3a–H3d), variant subtype benchmark (macro-F1 0.641 → 0.815), cCRE-ELS genome-wide, γ recalibration.

**Authorship (in progress)**: TDiG co-authors + original gDTR authors + additional contributors. Corresponding author: Sangwoo Kim (Yonsei).

**Manuscript status**: draft in preparation, target journal Nature Methods / Genome Research / Cell Systems tier.

**Original workshop paper (frozen)**: [`paper/gdtr_paper_ICML_3.pdf`](paper/gdtr_paper_ICML_3.pdf) and its LaTeX source, kept for reference. **The integrated manuscript will supersede this**.

---

## 1. The integrated framing

The workshop paper answered *when* a genomic residual-stream token stabilises (introduced `c(t)`) but treated "settling" as a single scalar. The integrated manuscript reframes this: **settling is a multi-axis phenomenon**, and each axis reads a distinct biological signature.

Three families of axes:

- **Metric axis** (from TDiG): direction / magnitude / geometry / distribution / path-efficiency (M1–M5)
- **Crossing-dynamics axis** (from followup): oscil, amplitude, argmin_layer, first_enter_layer (H3a–H3d)
- **Scale axis** (from followup): whole-genome (24 chromosomes, 2.94 B positions)

**Central claim**:

> Multi-axis settling is a genuine, discoverable geometric structure of Evo 2 7B's residual stream. Each axis reads a distinct biological signal. Combining axes gives a **+17.4 %p macro-F1 improvement** on 4-class variant subtype classification over the workshop's single-axis 32-d ΔD_cos baseline.

Full plan: [`docs/INTEGRATED_MANUSCRIPT_PLAN.md`](docs/INTEGRATED_MANUSCRIPT_PLAN.md).

## 2. Repository layout

```
gDTR-followup/
├── README.md                     # this file (integrated-manuscript home)
├── LICENSE                       # MIT
├── docs/
│   ├── INTEGRATED_MANUSCRIPT_PLAN.md   # section-by-section plan for the journal manuscript
│   ├── RESULTS.md                # every headline number with source script + JSON path
│   ├── METHODS.md                # every metric + feature definition
│   ├── REPRODUCE.md              # end-to-end reproduction recipe
│   └── CHANGELOG.md              # what was added when
├── src/                          # frozen gDTR framework (17 modules, byte-identical to paper)
├── scripts/                      # analysis pipelines
│   ├── wgs/                      # 8 files: STEP 1/2/4/5 + forward + labels + batch
│   ├── exp2/                     # 8 files: H2a/H2b/integrated/WGS-normalized
│   └── exp3/                     # 7 files: H3a/H3b/H3c/H3d
├── tdig_integration/             # (to be populated) TDiG scripts + results merged in
├── results/
│   ├── genome_summary/           # 25 files, 1.5 MB — 5-STEP outputs
│   ├── exp2/                     # H2 outputs
│   └── exp3/                     # H3 outputs
├── paper/
│   ├── gdtr_paper_ICML_3.pdf     # frozen workshop version
│   ├── gdtr_paper_ICML_3.tex     # frozen workshop LaTeX
│   ├── gdtr_paper.bib
│   └── (drafts/, figures/, tables/ — populated during manuscript work)
└── data/                         # data version pointers + requirements
```

## 3. Headline results (as of 2026-07-19)

| # | Question | Result | Source | Contribution |
|---|---|---|---|---|
| 1 | Do §3.1 direction claims hold at WGS scale? | **6 of 7** contexts ✓; 5'UTR flips sign uniformly | `results/genome_summary/wgs_context_summary.csv` | followup STEP 1 |
| 2 | Is `c(t)` the *only* axis, or are there more? | 5-D framework: c_t, oscil, first_enter, amplitude, argmin_layer | `results/exp3/H3d_advanced_context_test.json` | followup H3d |
| 3 | c_t vs oscil — same signal or orthogonal? | **Orthogonal**: splice acceptor d_c_t = −0.62, d_oscil = +0.50 | `results/exp3/H3a_context_test.json` | followup H3a |
| 4 | Does adding axes help downstream ML? | Subtype macro-F1 **0.641 → 0.815** (+17.4 %p) | `results/exp2/H2c_wgs_normalized_scoring.json` | followup STEP 3 |
| 5 | Does cCRE-ELS replicate genome-wide? | d_c_t = −0.132 ± 0.044 (paper chr22-only −0.118); **new d_oscil = +0.117** | `results/genome_summary/wgs_ccre_els_summary.json` | followup STEP 4 |
| 6 | How chromosome-invariant is γ = 0.397? | q70 SD = **0.15 %** across 24 chromosomes | `results/genome_summary/wgs_gamma_recalibration.json` | followup STEP 5 |
| 7 | What's the best single variant-effect feature? | ΔH_norm_L1 at L=8 — **AUROC 0.856** | TDiG §A1 (to be merged) | TDiG |
| 8 | Does variant peak layer index disrupted information level? | intron/3'UTR peak L=8; missense/synonymous peak L=27 | TDiG §A2 (to be merged) | TDiG |
| 9 | What happens at L29? | Context probing crashes 0.980→0.799; variant crashes 0.85→0.79 | TDiG §B1 (to be merged) | TDiG |
| 10 | Is there a low-dim structure behind M1–M5? | PC1 = 75.7 % variance; encodes bidirectional settling in a single direction | TDiG §C (to be merged) | TDiG |

**Full narrative**: [`docs/RESULTS.md`](docs/RESULTS.md)
**Integration status**: rows 7–10 marked "to be merged" require pulling TDiG artifacts (see [`tdig_integration/README.md`](tdig_integration/README.md) once created)

## 4. Locked parameters (paper reference; do not modify)

| Parameter | Value | Meaning |
|---|---|---|
| Evo 2 revision | `arcinstitute/evo2_7b`, SHA `bda0089f92582d5baabf0f22d9fc85f3588f6b58` | model weights |
| Weights MD5 | `359ef88ccac2a62644035578de8a7db4` | integrity check |
| `γ_cos` | 0.397 | frozen calibration (chr22 penultimate q70) |
| `L*` | 29 | canonical tap (Evo 2 pre-rotation) |
| Window / stride | 6,000 / 3,000 bp | forward pass config |
| Random seed | 42 | all CV / bootstraps |
| GENCODE | v44 | annotation |
| ClinVar | 2026-04-18 | variant set |
| GRCh38 | UCSC | reference |
| ENCODE | SCREEN v3 / v4 | regulatory annotations |

## 5. Citation

Once the integrated manuscript is submitted, this section will be updated. In the meantime, please cite the workshop paper as:

```bibtex
@inproceedings{cho2026gdtr,
  author = {Cho, Yoonjin and Kang, Jiheon and Park, Subin and Kim, Sangwoo},
  title = {GDTR: Layer-wise Settling Depth Reveals Biological Grammar in Genomic Foundation Models},
  booktitle = {Workshop on Generative and Agentic AI for Biology, ICML 2026},
  year = {2026}
}
```

TDiG's contributions will be jointly credited when the integrated preprint is posted.

## 6. License

MIT. See [`LICENSE`](LICENSE).

---
_Last updated: 2026-07-19_
