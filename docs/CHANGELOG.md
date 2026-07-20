# CHANGELOG

Append-only log of what was added to this repo and when.

## 2026-07-21 — EXP8 17-M-cell GO enrichment + 3-chromosome robustness + concept refactor

### New empirical content
- **`docs/NARRATIVE.md` Part 4 added** — extending the M-cell family to a third chromosome (chr21) and the GO story. Nine-subsection narrative covering motivation, EXP8 pipeline, 40 sign-flip GO terms, 2-chr and 3-chr robustness, the clean lipid transport example, cancer-panel driver hits, and the "settling depth → settling profile" concept refactor. Existing "What remains open" renumbered to Part 5.
- **`docs/EXP8_MULTI_CELL_GO_FINDINGS.md` (new, ~14 KB)** — technical companion. Pipeline description, all CSV outputs enumerated with schema, cluster analysis, sign-flip robustness tables, known bugs, reproduction commands, session timestamps.

### Headline additions (all data on DASH `/home/darejin/TDiG/exp8_multi_cell_go/`):
- **40 sign-flip GO terms** in combined chr17 + chr22 dataset (17 M-cells × 471 full BP GO terms)
- **8 pairs robust across chr17 + chr22 alone** (1.2 % of 680 tests) — sample-size-boost caveat
- **0 pairs robust across chr17 + chr22 + chr21** (of 280 tests); 3 robust across any 2 of 3
- **Clean example**: GO:0006869 lipid transport shows genuine sign-flip between M3_geo_a0.5_b1.0 (+) and M3_geo_a1.0_b0.0 (−), replicating on both chr17 and chr22
- **Cancer-panel driver hits**: RAD51C + RAD51D in M3_geo curvature × meiotic recombination; TP53 in M5_tau_refB × transcription regulation
- **Reference-token choice flips axis** within the same metric family: M5_tau_refA vs M5_tau_refC Spearman r = −0.78 on 471 GO terms

### Concept refactor
- **"Settling depth" (singular scalar) → "settling profile" (multi-axis family)**. Recommended manuscript language change; each M-cell reads a different axis, not a stronger version of the same axis. See NARRATIVE.md §4.7 and EXP8_MULTI_CELL_GO_FINDINGS.md §10.

### Infrastructure
- **New chr21 tier1 parquet** (13,365 windows × 7 M-cells) — first extension of the M-cell family beyond the TDiG team's chr17 + chr22 upstream data
- `wgs/scripts/e1_20_chr_forward_tier2.py` (new) — Evo 2 7B forward saving tier2-like scalars, optimised HDF5 (per-window chunks, no gzip) for ~9× speedup vs naive
- `wgs/scripts/e1_30_tier2_to_tier1.py` (new) — derives 7 of 17 M-cells from the tier2 h5 (M1_dir_refA, M2_mag_refA, M3_geo × 5 α/β variants). The other 10 M-cells (refB/refC + M4_set + M5_tau) require TDiG-team reference definitions.
- `exp8_multi_cell_go/scripts/e8_00 → e8_41` (six Python scripts totalling ~40 KB) — the EXP8 pipeline itself

### Infrastructural incident
- **TGIL_mutsig went offline at 2026-07-20 midnight KST** (GPU lease expiry). Truncated the planned chr20 forward (~10 % complete when host died) and the BATCH-3 DASH transfer (`results_cached` did not transfer). All critical artifacts backed up to DASH; see [`../RECOVERY_LOG.md`](../RECOVERY_LOG.md) for the incident narrative and next-session bootstrap.

### Known bug
- `e1_30` M2_mag_refA gamma calibration returns ~1.33 M on chr21 due to some intron positions having near-zero `norm_h_29`. Skip-window and outlier-clip filters didn't fully resolve; needs additional filter on `norm_h_29 < threshold` or a switch to median-based calibration for M2_mag primitives. Impact: M2_mag_refA on chr21 is de facto dead-cell. Other 6 M-cells on chr21 unaffected.

## 2026-07-19 — Reframed as integrated-manuscript repository

**Scope change**: this repository is no longer a "standalone follow-up analysis" but the **working home of the integrated gDTR journal manuscript**, currently in preparation.

### Reframing
- README rewritten to describe the integrated manuscript vision (workshop paper + TDiG + followup analyses → one unified journal submission).
- Added [`docs/INTEGRATED_MANUSCRIPT_PLAN.md`](INTEGRATED_MANUSCRIPT_PLAN.md) — section-by-section outline for the manuscript, with claim structure, results per section, and TDiG-integration checkpoints.
- Added [`tdig_integration/README.md`](../tdig_integration/README.md) — staging directory for TDiG artifacts + coordination checklist with TDiG team.
- Author roster is expanding to include TDiG co-authors (YAICON 8th 「띵디지놈」 team + original gDTR authors).
- Target venue: journal (Nature Methods / Genome Research / Cell Systems tier).
- Original workshop paper `paper/gdtr_paper_ICML_3.pdf` is kept as **frozen reference**; the integrated manuscript will supersede it.

### Followup analyses previously logged (still current)
- **Genome-wide aggregation (STEP 1)** — 2.94 B positions, direction transfer 6/7 contexts, 5'UTR chromosome-independent sign flip.
- **Advanced crossings WGS (STEP 2)** — 5-D framework: coding_exon late-clean-deep-commit at layer 27.6.
- **WGS-normalized variant scoring (STEP 3)** — 4-class subtype macro-F1 **0.6411 → 0.8149** (+17.4 %p).
- **cCRE-ELS WGS join (STEP 4)** — d(c_t) = −0.132 ± 0.044 + new d(oscil) = +0.117.
- **Per-chr γ recalibration (STEP 5)** — SD 0.15 % chromosome-invariance.
- EXP2 (H2a magnitude, H2b subtype, integrated H2b).
- EXP3 (H3a orthogonal axes, H3b variant re-forward, H3c committed enrichment, H3d 9-feature analysis).

### Documentation stack
- [`docs/RESULTS.md`](RESULTS.md) — every headline number with source script + JSON path.
- [`docs/METHODS.md`](METHODS.md) — every feature and metric definition.
- [`docs/REPRODUCE.md`](REPRODUCE.md) — end-to-end reproduction recipe.

### Not committed (regenerable)
- 37 GB per-chromosome per-position parquets (5,835 chunks)
- 4 GB per-chr position label npys (24 files)
- 13 GB Evo 2 7B weights (HF snapshot)
- 1.5 GB GRCh38 fa + GENCODE v44 GTF
- 184 MB ClinVar VCF + 145 MB UCSC cCRE bigBed

### Known outstanding items
- TDiG scripts + lightweight results not yet imported (pending coordination with TDiG team on authorship + data access — see [`tdig_integration/README.md`](../tdig_integration/README.md)).
- M1–M5 metric family not yet extended to WGS (would require raw h_ell across 24 chromosomes — deferred as future work in §4.4 of the manuscript plan).
- Manuscript draft itself not yet written; the [`INTEGRATED_MANUSCRIPT_PLAN.md`](INTEGRATED_MANUSCRIPT_PLAN.md) is the outline.

## 2026-07-19 (earlier) — Initial release: 5-STEP WGS downstream analysis
(See git commit `5af5f2b`.)
- First public commit of the followup analyses as an independent repo.
- Contained all scripts, results, documentation for followup work.
- Author scope at that point was original gDTR authors only.
- Superseded by the reframing above on the same day, once the decision was made to merge with TDiG for a joint journal manuscript.

## Original paper (2026-06-27, unchanged)
- Paper `paper/gdtr_paper_ICML_3.{pdf, tex}` frozen as accepted at ICML 2026 GenBio.
