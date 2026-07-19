# WGS Genome-wide aggregation report — 24 chromosomes

**Date**: 2026-07-19 09:56:27
**Total positions analyzed**: 2,937,756,000
**Chromosomes**: 24

## Per-context WGS-scale summary vs paper §3.1

| Context | n (WGS) | WGS mean c_t | Paper c_t | WGS d_c_t (mean±SD) | Paper d_c_t | WGS d_oscil |
|---|---|---|---|---|---|---|
| intron | 1,669,197,918 | 30.112 | 27.72 | +0.0000 ± 0.0000 | +0.000 | +0.0000 ± 0.0000 |
| splice_donor | 624,245 | 28.681 | 25.55 | -0.2829 ± 0.0509 | -0.354 | +0.2119 ± 0.0516 |
| splice_acceptor | 610,208 | 27.656 | 25.96 | -0.4822 ± 0.0742 | -0.340 | +0.4384 ± 0.0831 |
| coding_exon | 26,451,317 | 30.284 | 28.40 | +0.0250 ± 0.0194 | +0.080 | -0.1015 ± 0.0176 |
| 3utr | 50,075,413 | 29.618 | 27.74 | -0.1003 ± 0.0377 | -0.020 | +0.0815 ± 0.0370 |
| 5utr | 10,058,124 | 30.104 | 29.22 | -0.0110 ± 0.0221 | +0.200 | +0.0004 ± 0.0217 |
| intergenic | 1,180,738,775 | 30.653 | 28.66 | +0.1127 ± 0.0252 | +0.160 | -0.1080 ± 0.0237 |


## Calibration robustness (γ_cos = 0.397 frozen from paper chr22)
- WGS-mean intron c̄: 30.1061
- WGS-SD intron c̄:   0.1980
- Paper intron c̄:    27.72 (chr17+chr22 pool)
- Shift vs paper:     +2.3861 layers
- Chr-to-chr range:   [29.6967, 30.7545]

## Directional replication of paper §3.1
- **intron**: WGS d_c_t = +0.0000, paper = +0.000 → ✓ same sign
- **splice_donor**: WGS d_c_t = -0.2829, paper = -0.354 → ✓ same sign
- **splice_acceptor**: WGS d_c_t = -0.4822, paper = -0.340 → ✓ same sign
- **coding_exon**: WGS d_c_t = +0.0250, paper = +0.080 → ✓ same sign
- **3utr**: WGS d_c_t = -0.1003, paper = -0.020 → ✓ same sign
- **5utr**: WGS d_c_t = -0.0110, paper = +0.200 → ✗ opposite
- **intergenic**: WGS d_c_t = +0.1127, paper = +0.160 → ✓ same sign
