# REPRODUCE — End-to-end recipe

Step-by-step regeneration of every number in [`docs/RESULTS.md`](RESULTS.md), from raw model weights to final tables and figures.

Follows the same style as [heoneyzi/TDiG's `docs/reproduction.md`](https://github.com/heoneyzi/TDiG/blob/main/docs/reproduction.md).

---

## TL;DR

| Item | Value |
|---|---|
| Total wall clock | **~50 GPU-hours** for full 24-chr forward on 2 × B200; **~20 min CPU** for downstream analysis |
| Total storage on server | ~55 GB (37 GB parquets + 13 GB weights + 4 GB labels + ~1 GB refs) |
| Hardware tested | 2 × NVIDIA B200 (183 GB HBM each) at TGIL_mutsig cluster |
| Reference hardware | 1 × NVIDIA H200 (141 GB HBM) — paper's setup |
| Environment | conda env `gdtr`, python 3.11.15, torch 2.7.0+cu128, evo2 0.6.0, vortex 1.1.0, flash-attn 2.8.3 |
| Pipeline | 3 waves: (A) forward all chromosomes → (B) aggregate → (C) 5 downstream steps |

## 1. Hardware and environment

### 1.1 Server prerequisites

- ≥ 1 CUDA-compatible GPU with ≥ 141 GB HBM (H200 or B200 supported; B200 requires torch cu128 due to sm_100 Blackwell)
- ≥ 200 GB disk (for 37 GB parquets + 13 GB weights + 4 GB labels + 3 GB GRCh38 + 1.5 GB GENCODE + working scratch)
- python 3.11
- conda or mamba
- git, wget, tabix (optional)

### 1.2 conda env setup

```bash
mamba create -n gdtr python=3.11 -y
mamba activate gdtr

# torch — cu128 for B200 (Blackwell sm_100), cu126 for H200/A100
pip install torch==2.7.0 torchvision==0.22.0 --index-url https://download.pytorch.org/whl/cu128

# HuggingFace stack
pip install "transformers>=4.49,<4.52" huggingface_hub accelerate safetensors

# Scientific + genomics
pip install numpy==1.26.4 pandas==2.2.3 scipy scikit-learn matplotlib seaborn statsmodels pyarrow tqdm pyyaml
pip install pysam pyfaidx gffutils biopython pybedtools

# Evo 2 + vortex
pip install --no-build-isolation git+https://github.com/Zymrael/vortex.git
pip install --no-build-isolation git+https://github.com/ArcInstitute/evo2.git

# flash-attn (may need TMPDIR on same fs as pip cache)
TMPDIR=$HOME/tmp pip install flash-attn --no-build-isolation

# verify
python -c "import evo2, vortex, flash_attn, torch; \
  print(f'evo2 OK, torch {torch.__version__}, arch {torch.cuda.get_arch_list()}')"
```

### 1.3 Download Evo 2 7B weights

```bash
export HF_HOME=/YOUR/PATH/env/hf_cache
python -c "
from huggingface_hub import snapshot_download
snapshot_download(
  repo_id='arcinstitute/evo2_7b',
  revision='bda0089f92582d5baabf0f22d9fc85f3588f6b58',
  cache_dir='$HF_HOME/hub',
  allow_patterns=['*.json', '*.pt', '*.py', '*.md', '*.txt', 'tokenizer*'],
  max_workers=4,
)
"
```

Verify: `md5sum $HF_HOME/hub/models--arcinstitute--evo2_7b/snapshots/bda0089f.../evo2_7b.pt` → expect `359ef88ccac2a62644035578de8a7db4`.

### 1.4 Download reference data

```bash
cd data/reference_pointers/
# GRCh38
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip -k hg38.fa.gz

# GENCODE v44 basic annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
gunzip -k gencode.v44.annotation.gtf.gz

# ClinVar VCF (weekly, current release)
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi

# ENCODE cCRE (via UCSC bigBed)
wget https://hgdownload.soe.ucsc.edu/gbdb/hg38/encode3/ccre/encodeCcreCombined.bb
wget -O ~/bigBedToBed http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
chmod +x ~/bigBedToBed
~/bigBedToBed encodeCcreCombined.bb encodeCcreCombined.bed
grep -E "dELS|pELS" encodeCcreCombined.bed > ccre_els_wgs.bed
```

### 1.5 Verify environment

```bash
python -c "
import torch, evo2, vortex, flash_attn
print(f'torch {torch.__version__} cuda {torch.version.cuda}')
print(f'arch {torch.cuda.get_arch_list()}')
print(f'flash_attn {flash_attn.__version__}')
print(f'evo2 imports OK; vortex imports OK')
"
# quick GPU sanity — will load weights (~7 s) and run a 200-bp forward
python scripts/wgs/../env/test_evo2_b200.py  # skip if not on this repo layout
```

## 2. Wave A — WGS batch forward pass (**50 GPU-hours**)

For each chromosome:

```bash
python scripts/wgs/build_chr_position_labels_v2.py --chrom chr1
python scripts/wgs/e1_00_chr_forward.py --chrom chr1
python scripts/wgs/e1_10_aggregate_chr.py --chrom chr1
```

Or all 24 chromosomes in one shot (recommended, disconnect-safe):

```bash
tmux new -s wgs_batch
bash scripts/wgs/wgs_batch_runner.sh   # chr2..chrY sequential; chr1 run separately as above
# Runs for ~50 hours on 2 × B200. Check progress:
#   cat wgs/logs/batch_status.txt
```

Per-chromosome wall time (2 × B200 auto-shard):

| chr | wall | # windows | # positions |
|---|---|---|---|
| chr1 | 5.4 h | 82,984 | 230.5 M |
| chr2 | 4.4 h | 80,730 | 240.5 M |
| chr3 | 3.6 h | 66,097 | 198.1 M |
| ... | ~1.1 min / Mb | ... | ... |
| chrY | 29 min | 8,879 | 26.4 M |

**Output**: `wgs/results/chr{N}/chr{N}_per_position_chunk{0000..NNNN}.parquet` with columns `[window_idx, pos, c_t, oscil, n_enter, n_exit, below_frac, min_D, argmin_layer]`. Total 5,835 chunks, 37 GB, 2.94 B positions.

Also generates: `chr{N}_context_summary.csv` and `.json` per chromosome after aggregation.

## 3. Wave B — Genome-wide analysis (**~20 min CPU**)

### 3.1 STEP 1 — Genome-wide aggregation (< 2 s)

```bash
python scripts/wgs/e1_20_genome_wide_aggregate.py
```

Reads: `wgs/results/chr{N}/chr{N}_context_summary.csv` for all 24 chr.

Writes: `results/genome_summary/wgs_context_summary.{csv, json, png, pdf}`, `wgs_per_chr_*.csv`, `wgs_calibration_robustness.csv`, `wgs_context_report.md`.

Expected key result: splice_donor d(c_t) = **−0.283 ± 0.051**; intron shift **+2.39** vs paper.

### 3.2 STEP 2 — Advanced crossings genome-wide (~5.5 min)

```bash
python scripts/wgs/e3_70_wgs_advanced_features.py
```

Reads all 37 GB of per-position parquets.

Writes: `results/genome_summary/wgs_h3d_context_features.csv`, `wgs_h3d_per_chr_per_context.csv`, `wgs_h3d_summary.json`, `wgs_h3d_features.{png, pdf}`.

Expected key result: coding_exon `mean_argmin_layer = 27.6`, `frac_committed = 43 %`.

### 3.3 STEP 4 — cCRE-ELS WGS join (~6.8 min)

Requires `data_ref/ccre_els_wgs.bed` from step 1.4.

```bash
python scripts/wgs/e1_30_ccre_wgs_join.py
```

Reads 37 GB of parquets + 809 k cCRE-ELS records.

Writes: `results/genome_summary/wgs_ccre_els_{context.csv, summary.json, figure.png, figure.pdf}`.

Expected key result: WGS d(c_t) = **−0.132 ± 0.044**; WGS d(oscil) = **+0.117** (new).

### 3.4 STEP 5 — Per-chr γ recalibration (~2.2 min)

```bash
python scripts/wgs/e1_40_per_chr_gamma_recalibration.py
```

Reads all 37 GB of parquets (intron positions only).

Writes: `results/genome_summary/wgs_gamma_recalibration.{csv, json, png, pdf}`.

Expected key result: q70 = 0.5008 ± 0.0015 across 24 chr.

## 4. Wave C — Variant analysis (**~30 min GPU + ~2 min CPU**)

### 4.1 Verify paper AUROC (~30 s CPU)

Requires `results_cached/phase3_ensemble/variants_features_full.csv` from the original paper repo.

```bash
python scripts/exp2/e2_00_verify_paper.py
```

Expected: 32-d ΔD_cos AUROC = **0.8437** (paper 0.844).

### 4.2 H2a — variant-adaptive magnitude (~30 s CPU)

```bash
python scripts/exp2/e2_10_h2a_argmax_layer.py
python scripts/exp2/e2_30_visualize_h2a.py
```

Expected: `max_abs_dD` AUROC = **0.787** > best fixed L=30 (0.729).

### 4.3 H2b — subtype classification (~30 s CPU)

Requires ClinVar VCF (step 1.4).

```bash
python scripts/exp2/e2_20_h2b_subtype.py
python scripts/exp2/e2_30b_visualize_h2b.py
```

Expected: 4-class macro-F1 = **0.641**; OvR AUROC nonsense/synonymous = 0.889.

### 4.4 H3a — orthogonal axes on chr22 100 windows (~15 min GPU)

```bash
python scripts/exp3/e3_00_forward_windows.py --n-windows 100 --seed 42
python scripts/exp3/e3_10_compute_crossings.py
python scripts/exp3/e3_20_h3a_context.py
python scripts/exp3/e3_30_visualize.py
```

Expected: splice_donor d(oscil) = **+0.180**; splice_acceptor = **+0.502**.

### 4.5 H3b — variant re-forward (~67 min GPU on 2 × B200)

```bash
python scripts/exp3/e3_40_h3b_variant_forward.py
```

Expected: 8008 / 8008 variants, output `B_variants_oscil.parquet` (~260 KB, not committed to git).

### 4.6 H3c — committed enrichment (~1 s CPU)

```bash
python scripts/exp3/e3_50_h3c_committed.py
```

Expected: splice_acceptor deliberating OR = **2.85** vs intron.

### 4.7 H3d — advanced crossings on chr22 (~1 s CPU)

```bash
python scripts/exp3/e3_60_advanced_crossings.py
```

Expected: splice_donor `late_below_frac` d = **+3.18**.

### 4.8 Integrated H2b (~30 s CPU)

Requires 4.5 output.

```bash
python scripts/exp2/e2_40_integrated_h2b.py
python scripts/exp2/e2_40_visualize_integrated.py
```

Expected: cos32 + EXP3 + scalars macro-F1 = **0.6824**.

### 4.9 STEP 3 — WGS-normalized variant scoring (~1.4 min CPU)

Requires 3.1-3.3 outputs (WGS parquets + variant background) and 4.5 output.

```bash
python scripts/exp2/e2_50_variant_wgs_normalize.py
```

Expected: cos32 + scalars + raw + WGS-normalized macro-F1 = **0.8149**.

## 5. Directory layout expected during / after reproduction

```
tdig/
├── code/src/                                       # from this repo's src/
├── env/hf_cache/hub/models--arcinstitute--evo2_7b/ # 13 GB
├── data_ref/
│   ├── clinvar/clinvar.vcf.gz                       # 184 MB
│   ├── ccre_els_wgs.bed                             # ~20 MB filtered
│   └── encodeCcreCombined.bed                       # 145 MB
├── wgs/
│   ├── data/reference/hg38.fa                       # 3 GB
│   ├── data/reference/gencode.v44.annotation.gtf    # 1.5 GB
│   ├── data/labels/chr{N}_position_labels.npy       # ~4 GB total
│   ├── results/chr{N}/chr{N}_per_position_chunk*.parquet  # 37 GB total
│   ├── results/chr{N}/chr{N}_context_summary.csv    # ~1 MB total
│   ├── results/genome_summary/                      # 1.5 MB
│   ├── logs/                                         # forward-pass logs
│   └── scripts/                                      # from this repo
├── exp2_variant_downstream/{scripts,results,figures}
├── exp3_threshold_crossing/{scripts,results,figures}
└── results_cached/phase3_ensemble/variants_features_full.csv    # paper cache, 15 MB
```

## 6. Troubleshooting

### CUDA arch mismatch on B200 (Blackwell)

If `torch.cuda.get_arch_list()` does NOT include `sm_100`, you have the wrong torch build. Reinstall:

```bash
pip install --force-reinstall torch==2.7.0 --index-url https://download.pytorch.org/whl/cu128
```

Verify: should show `arch ['sm_75', 'sm_80', 'sm_86', 'sm_90', 'sm_100', 'sm_120', 'compute_120']`.

### flash-attn "cross-device link" error during build

```bash
mkdir -p ~/tmp
TMPDIR=~/tmp pip install flash-attn --no-build-isolation
```

### Cross-device tensor error in H3a forward (`cuda:0 vs cuda:1`)

Evo 2 auto-shards blocks 0–15 to cuda:0 and 16–31 to cuda:1. Ensure the `cosine_lens` function in [`scripts/exp3/e3_00_forward_windows.py`](../scripts/exp3/e3_00_forward_windows.py) does `.to(h_norm.device)` on each `h_l` before the cosine dot product. This is the v2 version in this repo.

### Chunk file resume

If a per-chr forward crashes mid-way, the chunks up to the last flush are already valid. Rerun `e1_00_chr_forward.py --chrom chrX --start-idx N` where N = 500,000 × (last_chunk_id + 1).

### bigBedToBed missing

Download from UCSC:

```bash
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
chmod +x bigBedToBed
```

## 7. Verifying reproduction

After running the full pipeline, compare `results/genome_summary/wgs_context_summary.csv` to this repo's committed version. Every d value should match to 3 decimals if the same random seed (42) and the same input BEDs / VCF are used. Small deviations (< 1 %) can arise from:

- Different ClinVar release (paper uses 2026-04-18 archive; this repo uses the current weekly release)
- Different ENCODE cCRE version (paper uses SCREEN v3; this repo uses UCSC encodeCcreCombined which is SCREEN v4)
- Different WGS chunk order due to variable Ns
