# TDiG Recovery Log — TGIL_mutsig shutdown 2026-07-20

## Event

At approximately **2026-07-20 23:58 KST** the compute host `TGIL_mutsig`
(59.150.33.1:48203) became unreachable. All SSH attempts return
`Operation timed out` at the banner-exchange stage — the same pattern
previously observed on 2026-07-15.

Root cause: GPU lease expired at midnight boundary. Host was likely
taken offline as part of the routine lease-end procedure. **User has
confirmed that TGIL is completely shut down** (2026-07-21).

## Session activities interrupted

The following operations were in progress or planned:

### Interrupted mid-transfer
- **BATCH 3 (DASH sync)**: `data_ref` + `results_cached` transfer
  - `data_ref/` reached ~380 MB out of 623 MB (~60% complete)
  - `results_cached/` (729 MB) never transferred (0%)

### Interrupted mid-compute
- **chr20 extended forward** (`e1_20_chr_forward_tier2.py --chrom chr20`)
  - Started 21:55 KST, was at ~10-15% completion at TGIL death time
  - Partial `chr20_tier2_scalars.h5` (~6.4-10 GB) was on TGIL only
  - Not recoverable without host access

### Never started (deferred)
- **chr19 extended forward** — planned but skipped due to time budget
- Any subsequent per-chromosome forwards

## What is safe on DASH (as of 2026-07-21 00:20 KST)

**5.5 GB backed up** across the CRITICAL tier of the project. See
`REPRODUCE.md` §2 for the current directory structure and `TDIG_MANIFEST.md`
for the full artifact inventory.

Critical items confirmed present:
- Upstream chr17/22 tier1 parquets (from TDiG team)
- Our chr21 tier1 parquet + gamma calibration
- Full EXP8 analysis pipeline (5 scripts + 16 CSVs + 14 figures + docs)
- All experiment directories (exp2 through exp8)
- Paper drafts (139 MB)
- All key code

## What is NOT on DASH — recovery plan

| Item | Was on | Status | Recovery |
|---|---|---|---|
| chr20 partial tier2 h5 | TGIL | Lost | Rerun `e1_20 --chrom chr20` on new host (~1 h) |
| chr21_tier2_scalars.h5 (41 GB) | TGIL | Lost | Rerun `e1_20 --chrom chr21` on new host (~1 h) |
| wgs/results/chr1..chrY cosine forwards | TGIL | Lost | Rerun `e1_00 --chrom chrN` per chr (~2.3 h each) |
| tdig_integration/hf_cache | TGIL | Lost | Re-download HF `darejinn/TDiG-evo2-hidden-states` |
| env/ (Evo 2 weights) | TGIL | Lost | Re-download HF `arcinstitute/evo2_7b` |
| wgs/data/ (hg38, GENCODE, ClinVar) | TGIL | Lost | Re-download from UCSC / GENCODE / NCBI |
| data_ref/ remainder (~243 MB) | TGIL | Partial on DASH | Regenerate specific files, or recover from TGIL if resurrected |
| results_cached/ (729 MB) | TGIL | Lost | Regenerate via `code/scripts/` runs; may recover from git if code tracked outputs |

## Next-session bootstrap checklist

When resuming on a new GPU host:

```bash
# 1. Set up environment
git clone https://github.com/darejinn/gDTR-PoC          # phase1-5 script library
conda create -n gdtr python=3.11
conda activate gdtr
pip install numpy pandas scipy matplotlib h5py pyarrow scikit-learn torch
pip install evo2  # official evo2 package

# 2. Download public prerequisites
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
gunzip gencode.v44.annotation.gtf.gz
# Evo 2 weights auto-download via evo2 package

# 3. Pull working DASH backup
rsync -av darejin@dash.yuhs.ac:2010:/home/darejin/TDiG/ ./tdig_work/

# 4. Regenerate lost tier2 (if analyzing chr21+22 etc.)
python tdig_work/wgs/scripts/e1_20_chr_forward_tier2.py --chrom chr21
python tdig_work/wgs/scripts/e1_30_tier2_to_tier1.py --chrom chr21

# 5. Rerun EXP8 pipeline (should exactly match saved results)
cd tdig_work/exp8_multi_cell_go
python scripts/e8_00_per_cell_go_enrichment.py
# ...
```

## Open bugs to fix

1. **e1_30 M2_mag_refA gamma blowup** — gamma calibration returns 1.3M
   for chr21. Not fixed by skip-window filter or upper-outlier clip in
   pt37. Likely due to some intron positions having near-zero norm_h_29
   causing `r = norm_h_ell / norm_h_29` to explode.
   - Fix candidates:
     - Median instead of q70 for M2_mag calibration (more robust to
       right-tail outliers)
     - Explicit filter: skip positions where `norm_h_29 < 1e-3` (or
       some minimum threshold) before computing r
     - Consider whether M2_mag_refA is a legitimate settling depth at
       all — the "|r-1|" definition may not be the right primitive.

2. **e1_20 skip-window handling** — currently `done_mask[i] = 1` for
   both successful compute AND N-skip. Downstream e1_30 has to filter
   skips via `window_start > 0`. Better: use `done_mask == 1` for
   success only, `done_mask == 2` for skips.

## Reference

- Full session narrative: MemPalace diary entries `wing_claude-opus-4-7`
  for 2026-07-20 pt29 through pt38.
- Manuscript decisions: MemPalace drawer `wing_gdtr/decisions/` —
  "settle depth → settle profile" (pt32).
- Key findings: MemPalace drawer `wing_gdtr/findings/` for
  sign-flip analysis, cancer-panel driver hits, lipid transport
  robust finding.
