# tdig_integration/ — Merging TDiG (heoneyzi et al.) into the unified manuscript

This directory is the staging area for TDiG artifacts that will contribute to the integrated journal manuscript.

## Origin

TDiG source repository: [heoneyzi/TDiG](https://github.com/heoneyzi/TDiG)
- **Team**: YAICON 8th (「띵디지놈」 team), 1st Prize
- **Focus**: Multi-metric settling family (M1–M5), 3 reference variants, 98-cell context heatmap, per-consequence variant analysis, PCA interpretability, L29 phase transition
- **Documentation**: TDiG's [`docs/`](https://github.com/heoneyzi/TDiG/tree/main/docs) is exceptionally thorough (design decisions, metric definitions, reproduction guide, thesis)

## What lives here

**Currently empty** — this directory is staging for imported TDiG artifacts. The plan (per [`../docs/INTEGRATED_MANUSCRIPT_PLAN.md`](../docs/INTEGRATED_MANUSCRIPT_PLAN.md) §3) is to import:

| Artifact | Manuscript section | TDiG source location |
|---|---|---|
| M1–M5 metric family definitions | §2.3 Methods | [`docs/metric_definitions.md`](https://github.com/heoneyzi/TDiG/blob/main/docs/metric_definitions.md) |
| Design decisions log | manuscript revision history | [`docs/design_decisions.md`](https://github.com/heoneyzi/TDiG/blob/main/docs/design_decisions.md) |
| M1–M5 forward pass implementation | `scripts/tdig/` | [`scripts/15_chr22_forward.py`](https://github.com/heoneyzi/TDiG/blob/main/scripts/15_chr22_forward.py) |
| 98-cell context heatmap (chr22 + chr17) | §3.2 Results, headline figure | `results/context_separation/` |
| Per-layer probing AUROC (context) | §3.7 L29 phase transition | `results/analysis_BD/per_layer_auroc.csv` |
| γ ablation q50/q70/q90 | Extended Data | `results/gamma_ablation/` |
| Per-consequence variant AUROC | §3.6 Task C | `results/variant_analysis_scalars/` + `results/variant_per_consequence/` |
| PCA analysis (PC1 = 75.7 %) | §3.8 | `results/analysis_BD/metric_pca_corr.*` |
| 17-cell variant Δc analysis | §3.6 Task B extension | `results/variant_settling_cells/` (was running when TDiG last updated) |

## Import strategy (planned)

1. **Fork or clone** TDiG at a locked commit hash; record hash in [`docs/CHANGELOG.md`](../docs/CHANGELOG.md).
2. **Copy scripts** relevant to the manuscript into `scripts/tdig/` here.
3. **Copy lightweight results** (CSVs, JSONs, PNGs) into `results/tdig/`.
4. **Do NOT copy** the 883 MB tier-1 settling parquet or 108 GB tier-3 hidden states — these live on TDiG's HF space + minimal-archive split-file staging. Reference their download URL in [`../docs/REPRODUCE.md`](../docs/REPRODUCE.md).
5. **Reconcile shared code**: `src/` here contains a frozen copy of `src/ur_gdtr.py`, `src/tuned_lens.py`, etc. from the paper commit. TDiG uses these too; verify no divergence.
6. **Author attribution**: every imported script header must credit TDiG team (see [`docs/AUTHORSHIP.md`](../docs/AUTHORSHIP.md) once written).

## Coordination checklist with TDiG team

Before importing:

- [ ] Confirm co-authorship with TDiG team lead (heoneyzi)
- [ ] Agree on which TDiG results are ready for manuscript (some analyses were "running" in TDiG's last public commit)
- [ ] Get access to chr22 raw cache (`chr22_cache.h5`) for regeneration if needed
- [ ] Align on figure attribution (which figures TDiG generated vs adapted here)
- [ ] Confirm authorship order

Until this coordination is complete, this directory stays empty and the manuscript plan references TDiG results by URL only. The manuscript will not use TDiG-derived numbers as headline claims until authorship + data access are settled.

## Related

- Integrated manuscript plan: [`../docs/INTEGRATED_MANUSCRIPT_PLAN.md`](../docs/INTEGRATED_MANUSCRIPT_PLAN.md)
- Results narrative (this repo's contributions only): [`../docs/RESULTS.md`](../docs/RESULTS.md)
- TDiG own documentation: <https://github.com/heoneyzi/TDiG/tree/main/docs>
