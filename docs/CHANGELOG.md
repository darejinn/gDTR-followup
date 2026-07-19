# CHANGELOG

Append-only log of what was added to this repo and when.

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
