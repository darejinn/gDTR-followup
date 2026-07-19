# Integrated Manuscript Plan — gDTR (revised & extended)

**Target**: Journal (Nature Methods / Genome Research / Cell Systems tier)
**Status**: Planning document — the manuscript itself is being drafted from this plan.
**Prepared**: 2026-07-19

This document defines how three bodies of work merge into **one unified journal-length manuscript** on gDTR:

1. **gDTR original paper** (Cho, Kang, Park, Kim — ICML 2026 GenBio Workshop, 4 pages, accepted). Base: settling depth `c(t)`, γ = 0.397, splice/enhancer shallowness, motif↔flank bidirectionality, ClinVar variant scoring.
2. **TDiG** (heoneyzi et al., YAICON 8th 1st Prize). Extension: M1–M5 metric family, 3 reference variants, 98-cell context heatmap, per-consequence layer split, PCA interpretability, L29 phase transition.
3. **gDTR-followup** (this repo, July 2026 sessions). Extension: WGS 24-chromosome replication (2.94 B positions), threshold-crossing dynamics (H3a–H3d), variant scoring benchmark (macro-F1 0.641 → 0.815), cCRE-ELS genome-wide, γ recalibration.

Authorship (proposed):
- **First author(s)**: Yoonjin Cho + [TDiG lead] (co-first)
- **Middle authors**: Jiheon Kang, Subin Park, [TDiG team members from YAICON 8th]
- **Senior/corresponding**: Sangwoo Kim

---

## 1. Story arc (what changed from the workshop paper)

The workshop paper answered *when* a token stabilises (introduced `c(t)`) but treated "settling" as a single scalar. The integrated manuscript reframes this: settling is a **multi-axis phenomenon**, and each axis maps to a distinct biological signature.

The single-axis workshop framing is preserved as **§3 Baseline** and extended along three orthogonal directions:

- **Metric axis** (TDiG M1–M5): direction, magnitude, geometry, distribution, path-efficiency
- **Crossing-dynamics axis** (followup H3): oscil, amplitude, argmin_layer, first_enter_layer
- **Scale axis** (followup WGS): from chr17+chr22 to whole genome (24 chr, 2.94 B positions)

The claim of the integrated paper is:

> **Multi-axis settling is a genuine, discoverable geometric structure of Evo 2 7B's residual stream. Each axis reads a distinct biological signal. A downstream ML task (variant subtype classification) obtains a +17.4 %p macro-F1 improvement over the workshop baseline by combining axes.**

## 2. Section-by-section structure

### Abstract (200 words)
Extend the current 250-word abstract to reflect (a) multi-axis framing, (b) genome-scale replication, (c) downstream +17.4 %p, (d) TDiG collaboration.

### §1 Introduction
- Existing 4-axis interpretability gap (paper §1)
- Add: the settling problem specifically has **three distinct definitions** (Def 1 stops evolving / Def 2 converges to h_29 / Def 3 converges to h_norm). Existing literature conflates them.
- Preview: unified framework + WGS + downstream benchmark

### §2 Methods — Unified settling framework

#### §2.1 GDTR base (from workshop paper §2)
- `c(t)`, running-min envelope, γ_cos = 0.397, L* = 29 (Evo 2 idle L31)
- Two-sidedness argument (motif detection vs contextual determinacy)

#### §2.2 Three settling definitions (from TDiG design v2)
- Def 1 trajectory stops evolving
- Def 2 converges to h_29 (pre-final-block)
- Def 3 converges to h_norm (post-final-norm)
- Explicit table of which metric addresses which Def

#### §2.3 Metric family M1–M5 (from TDiG)
- M1 direction (`c_dir`) — 3 reference variants A/B/C
- M2 magnitude (`c_mag`) — diagnostic, residual-accumulation artifact
- M3 geometry (`c_geo`) — reference-free velocity + curvature, α/β ablation
- M4 distribution (`c_M`) — Mahalanobis with Σ_ref shrinkage
- M5 path efficiency (`c_τ`) — tortuosity

Persistence window W = 3 uniform. γ q70 calibration matrix, ablations at q50/q90 supplementary.

#### §2.4 Crossing-dynamics axes (from followup H3d)
- `oscil` — # of extra crossings beyond first dip
- `amplitude` — γ − min_D
- `first_enter_layer` — WHEN the first crossing happens
- `argmin_layer` — WHERE the trajectory bottoms out
- `late_below_frac`, `mid_below_frac`, `early_below_frac`
- `longest_below_streak`

**Key claim**: These axes are orthogonal to M1's running-min collapse. Splice acceptor oscil d = +0.502 while c_t d = −0.62 (opposite signs).

#### §2.5 Reference variants matrix (A/B/C, from TDiG)
- A: raw both sides — pure hidden-state comparison
- B: RMSNormed both — DTR-style, γ-symmetric
- C: h_ℓ raw vs h_norm — asymmetric (workshop paper's default)

Applied to every reference-dependent metric. Isolates RMSNorm γ asymmetry effect.

#### §2.6 Complete settling-cell inventory

| Family | Cells | Definitions targeted |
|---|---|---|
| M1_dir × {A,B,C} | 3 | Def 2/2/3 (C degenerate per §2.5 rationale) |
| M2_mag × {A} | 1 | residual accumulation (diagnostic) |
| M3_geo × 5 α/β | 5 | Def 1 |
| M4_set × {A,B,C} | 3 | Def 2 |
| M5_tau × {A,B,C} | 3 | Def 2 (via Option B lock) |
| **H3 crossing family** | 9 | axis (WHEN / HOW DEEP / HOW OFTEN) |
| **Total** | **24** | complementary axes on the residual stream |

### §3 Results

#### §3.1 Baseline reproduces at WGS scale (from followup STEP 1)
- 24 chromosomes, 2.94 B positions
- Splice donor d(c_t) = −0.283 ± 0.051 (paper: −0.354)
- Splice acceptor d(c_t) = −0.482 ± 0.074 (stronger than paper's −0.340)
- Intron c̄ shifts +2.39 layers uniformly (direction transfers, absolute does not)
- **New finding**: 5'UTR chromosome-independent sign flip vs paper (+0.20 → −0.011)

#### §3.2 98-cell context × metric heatmap (from TDiG headline)
14 metric cells × 7 biological contexts, split into:
- Panel a: c(t) family (workshop paper §3.1 + WGS extension)
- Panel b: M1–M5 family (TDiG contribution) at chr17 + chr22 (available today)
- Panel c: WGS pooled c(t) + oscil (this repo)

**Discovery axis dissociations**:
- M1_dir_refA and M3_geo track different contexts
- M2_mag captures residual accumulation, not settling — reported as artifact
- Splice sites cluster distinctively across all axes; coding_exon shows a *late clean deep-commit* pattern

#### §3.3 Bidirectionality on axes beyond c(t) (paper §3.2 extended)
- Original motif edit vs flank shuffle (paper Table 1)
- Extended: same edits on H3 crossing axes — do all axes see two-sidedness or just some?

#### §3.4 Orthogonal-axis discovery (from followup H3a)
- `c(t) × oscil` orthogonal axes at chr22 79 windows
- Genome-scale replication (this repo, STEP 2): full 24 chr, per-context d(oscil)
- Splice donors d(oscil) = +0.217 WGS (vs chr22 sample +0.180); splice_acceptor +0.451 WGS

#### §3.5 Late-clean-deep-commit at coding exons (from followup STEP 2)
- coding_exon `mean_argmin_layer = 27.6` (deeper than splice donors 23.5)
- coding_exon `frac_committed = 43 %` vs intron 0.85 % (51× enrichment)
- Interpretation: codon-usage / reading-frame commitments happen *late* and *decisively*, invisible to `c(t)` alone

#### §3.6 Variant scoring benchmark (paper §App C extended with H3 + WGS-normalized features)

Two tasks on 8,008 SNVs × 15 cancer genes:

Task A: binary P/LP vs B/LB (paper §App C)
- cos32 vector: AUROC 0.844 (paper reproduced to 3 decimals)
- + H3 + scalars: 0.849 (marginal — cos32 already saturated)

Task B: 4-class subtype (missense/nonsense/synonymous/canonical_splice)
- cos32 vector: macro-F1 **0.641**
- + H3 raw: **0.659**
- + H3 WGS-normalized: **0.686**
- + all (cos32 + scalars + raw H3 + WGS-normalized H3): **🎯 0.815** — **+17.4 %p over baseline**

Task C: **variant-consequence layer split** (from TDiG §A2)
- intron / 3'UTR peak at L=8 (sequence-level information)
- missense / synonymous / 5'UTR peak at L=27 (protein-semantic information)
- Per-consequence AUROC 0.809–0.927
- **Interpretation**: variant peak layer indexes the semantic level being disrupted

Task D: **best single feature** for pathogenicity
- ΔH_norm_L1 at L=8: AUROC = **0.856** (TDiG §A1)
- max_abs_dD (variant-adaptive magnitude): AUROC = **0.787** (followup H2a)

#### §3.7 L29 phase transition (from TDiG §B1)
- Context probing AUROC at L=29 crashes 0.980 → 0.799
- Variant AUROC crashes 0.85 → 0.79 at same layer
- PC1-metric correlations collapse
- **Interpretation**: L29 is Evo 2's architectural rotation layer where the residual stream is projected into h_norm's frame. All settling processes discontinuously terminate here.

#### §3.8 PCA interpretability (from TDiG §C)
- PC1 = 75.7 % of variance = "transcribed-region settling readiness axis"
- PC1+ encodes Def 1 (M3_geo); PC1− encodes Def 2 (M1_dir, M2, M5_tau_refB)
- Bidirectional settling has a geometric origin (single PC direction)

#### §3.9 cCRE-ELS WGS join (from followup STEP 4)
- Paper Fig 2 panel b replication at 809,429 dELS+pELS records × 24 chr
- WGS d(c_t) = −0.132 ± 0.044 (paper chr22-only: −0.118)
- **New**: WGS d(oscil) = +0.117 — cCRE-ELS also participates in orthogonal oscil axis

#### §3.10 Calibration robustness (from followup STEP 5)
- Per-chromosome q70 recalibration: SD 0.15 % across 24 chr
- Paper's "single calibration transfers" claim confirmed in a strong sense
- Discussion of paper §2 wording (q70 of penultimate running-min vs all-layer min_D)

### §4 Discussion & Limitations

#### 4.1 What this framework establishes
- Multi-axis settling is real and distinct from any single-metric collapse
- Each axis has a distinguishable biological correlate
- The framework is downstream-useful (variant subtype +17.4 %p)

#### 4.2 The three settling definitions gap (TDiG's key limitation)
- Def 3 (convergence to h_norm) has no metric in the current M1–M5 family
- **Partial resolution from followup**: oscil axis is defined against γ, not h_norm — but reveals structure that h_norm-agnostic metrics see

#### 4.3 Composition and off-manifold shuffle (paper §App G)
- Retained from workshop paper

#### 4.4 Cross-architecture / MLM transferability
- HyenaDNA replicates c(t) direction (workshop paper)
- NT-v2, DNABERT-2 tokenisation-limited (workshop paper)
- Full M1–M5 replication on Caduceus: deferred future work

#### 4.5 Causal / SAE grounding
- Present framework is geometric-descriptive; no causal intervention
- Deferred to companion / next paper

### §5 Conclusion
- Multi-axis settling framework for Evo 2 7B
- Distinct biology per axis; 5-D reading beats 1-D
- Public repro at [`darejinn/gDTR-followup`](https://github.com/darejinn/gDTR-followup)

### Extended Data / Supplementary

- Full γ ablation (q50/q70/q90 × 24 metric cells)
- M3 α/β ablation (5 cells × context)
- Reference-variant ablation (A/B/C × 3 metrics)
- Cross-architecture summary (paper Table 2 extended)
- Splice positional fine-profile (paper Fig A4 + WGS extension)
- Per-chromosome tables (all 24 × all metrics)
- Full downstream benchmark tables (all feature-set combinations)

## 3. Data and code artifacts required for the final manuscript

| Artifact | Status | Source |
|---|---|---|
| Workshop paper `c(t)` on chr17 + chr22 | ✅ available | `results_cached/phase2.4/` |
| 8,008 variant features (cos32) | ✅ available | `results_cached/phase3_ensemble/` |
| 8,008 variant H3 re-forward | ✅ available | `B_variants_oscil.parquet` (regenerable) |
| WGS 24-chr per-position | ✅ available | `wgs/results/chr{N}/` (37 GB) |
| WGS 24-chr context summaries | ✅ available | `results/genome_summary/` |
| cCRE-ELS WGS join | ✅ available | `results/genome_summary/wgs_ccre_els_*` |
| γ recalibration | ✅ available | `results/genome_summary/wgs_gamma_recalibration_*` |
| **TDiG M1–M5 at chr17+chr22** | ⚠ from TDiG repo | need to copy or re-run |
| **TDiG 98-cell heatmap** | ⚠ from TDiG repo | need to copy or re-run |
| **TDiG per-consequence variant analysis (§A2)** | ⚠ from TDiG repo | need to copy |
| **TDiG L29 phase transition analysis (§B1)** | ⚠ from TDiG repo | need to copy |
| **TDiG PCA analysis (§C)** | ⚠ from TDiG repo | need to copy |
| **TDiG γ ablation (§B2)** | ⚠ from TDiG repo | need to copy |
| M1–M5 extension to WGS | ✘ not run | future work (mentioned in §4) |

## 4. Journal target details

Preferred order:
1. **Nature Methods** (methodology focus, biology impact, code+data sharing culture)
2. **Genome Research** (genomics focus, permissive length, molecular biology audience)
3. **Cell Systems** (systems-biology framing works well for multi-axis)
4. **eLife** (rapid, open, no page limit)
5. **Bioinformatics** (methodology, shorter but respected)

Target length: 5,000 words main text + Extended Data (unlimited).

## 5. Next steps (execution roadmap)

### Immediate (this session)
- ✅ This plan document
- Reconfigure `darejinn/gDTR-followup` as manuscript home (this doc + code + results + paper draft)
- Update README to reflect integrated-manuscript vision
- Reference TDiG explicitly as sibling work merging into this repo

### Short-term (this week)
- Pull TDiG scripts + results into `gDTR-followup/tdig_integration/`
- Draft §2 Methods (unified settling framework)
- Draft §3 Results (with WGS + TDiG numbers side-by-side)

### Medium-term (2–3 weeks)
- Full manuscript draft
- Extended Data compilation
- Figure regeneration at journal-quality

### Author coordination checklist
- TDiG team (heoneyzi et al.): confirm co-authorship, get chr22 raw cache access, agree on figure attribution
- Original team (Cho, Kang, Park, Kim): confirm scope + author order
- Sangwoo Kim (senior): manuscript review + submission
