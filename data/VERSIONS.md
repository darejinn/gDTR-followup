# Version map — gDTR project (paper iterations + companion split)

This file is the canonical pointer for "which folder / file = which version".
Nothing here is deleted between revisions: every prior submission state
remains on disk for diff and audit purposes.

## TL;DR — what to read

- **Current ICML 2026 workshop short paper (Paper 1, gDTR mechanistic probe)**:
  [`ICML_0429_v3/gdtr_paper_ICML.pdf`](ICML_0429_v3/gdtr_paper_ICML.pdf).
  Source: [`ICML_0429_v3/gdtr_paper_ICML.tex`](ICML_0429_v3/gdtr_paper_ICML.tex).
  Build instructions: [`ICML_0429_v3/README.md`](ICML_0429_v3/README.md).
- **Companion paper (Paper 2, ‖Δh‖₂ scorer benchmark)**:
  [`Paper2_DeltaH.docx`](Paper2_DeltaH.docx) — kept as a separate manuscript
  per the 2026-04-29 split decision; **must not be merged back into Paper 1**.

## Paper 1 — gDTR (mechanistic probe)

The paper has been iterated through three "shape" rewrites and several
DOCX polish passes. All sources are kept.

| Stage | Folder / file | Date | Status |
|---|---|---|---|
| v1 (8-page Nature-style) | `ICML_0429_v1 2/` | 2026-04-29 | Original submission shape, retained for diff. Minor table-header polish committed without rebuild. |
| v3 (4-page ICML workshop, narrative restructure) | `ICML_0429_v3/` | 2026-04-28 → 2026-05-05 | **Canonical short-paper version.** Built PDF: `ICML_0429_v3/gdtr_paper_ICML.pdf`. Includes the 2026-05-05 v11 layout / audit polish. |
| DOCX writing path | `Paper1_gDTR.docx`, `Paper1_gDTR_0429.docx`, `Paper1_gDTR_0429_v5.docx`, `Paper1_gDTR_0429_v6.docx`, `Paper1_gDTR_0429_v7.docx` | 2026-04-28 → 2026-04-29 | DOCX iteration history. v7 is the latest DOCX. The LaTeX in `ICML_0429_v3/` is downstream of v7 plus the v11 polish. |

The v1 → v3 narrative reframing is logged in
[`ICML_0429_v3/V3_REWRITE_NOTES.md`](ICML_0429_v3/V3_REWRITE_NOTES.md);
the v3 → v11 audit / layout polish is logged in the "What changed in v11"
section of [`ICML_0429_v3/README.md`](ICML_0429_v3/README.md).

## Paper 2 — ΔH scorer benchmark (companion split)

[`Paper2_DeltaH.docx`](Paper2_DeltaH.docx) is the standalone DOCX skeleton
for the **separate** companion manuscript, decided on 2026-04-29.

> **Important.** Paper 2 (ΔH benchmark, AUROC 0.926 on 10,910 ClinVar
> variants) and Paper 1 (gDTR mechanistic probe) are intentionally split
> and must not be merged. Paper 1 keeps `‖Δh‖₂` only as a non-headline
> sanity check on cosine-trajectory information content
> (AUROC 0.844 in Appendix C); the headline ‖Δh‖₂ benchmark belongs to
> Paper 2.

## Reproducibility map

The figures and numbers in `ICML_0429_v3/gdtr_paper_ICML.pdf` are
backed by:

| Result file | Used for |
|---|---|
| `results/exp1_entropy/`, `results/exp1_entropy_meta.json` | §3.1 entropy-control paragraph (chr22 control panel; ρ=−0.079, d_raw=−0.452 → d_resid=−0.583). |
| `results/exp2_shuffled/`, `results/exp2_shuffled_meta.json` | §3.2 motif perturbation (real c̄=26.77; flank-shuffle d=+0.51, mut-GT d=−0.09). |
| `results/figures_v3/` | All v3-era figure regeneration outputs (PNG + PDF + meta JSON). |
| `results/phase4/per_model_summary.json` | Cross-architecture replication tables (Table 5–7) and Fig. 4 / Fig. 5. |
| `results/phase3_main/` | ClinVar variant features (cached, see .gitignore for the heavy CSV). |
| `scripts/make_v3_figures_remote.py` (within `ICML_0429_v3/scripts/`) | Master figure regeneration script (runs on the H200 server). |

Per-phase analysis docs live under `docs/findings/`. Top-level
`README.md` retains the full project narrative.

## Stray files at repo root

These exist in the working tree and are intentionally **not** in any
commit until the user labels them:

- `05042134_manuscript.pdf`, `document_pdf.pdf` — provenance unclear;
  candidates for deletion or for moving into a `vendored/` folder.
- `ICML_0429_v3.zip` — backup snapshot of the v3 folder; ignored by
  `.gitignore` (`ICML_*.zip`) since the unzipped folder is the
  canonical version.
