# Result 2 — The scar program is larger and faster at P7/P8 than predicted by injury severity

**Question.** After myocardial infarction (MI), does the scar gene program rise more, and keep rising, in the
non-regenerative heart (P7/P8) than at P1, and does this excess remain after adjusting for how severe the injury is?

**Answer in the paper.** Yes: at P1 induction is modest and stable, at P8 it rises progressively; an independent
laboratory shows the same at day 3; the excess survives adjustment for a transcriptional injury-severity score in
both datasets (Figure 2).

## Data

- GSE123863 (Wang Z 2019): bulk RNA-seq of ventricles, P1 and P8, MI vs sham, days 1.5, 3 and 7, n = 3 per group.
- GSE198300 (Kankuri 2023): bulk RNA-seq of left ventricle (CD-1), P1 and P7, MI vs sham, days 3 (n = 4) and 21 (n = 3).

Download with `bash download_bulk_data.sh` (from `code/`).

## Method

Counts → CPM; genes with mean CPM < 1 removed; gene-set score = mean z-score of log2(CPM + 1) across all samples of a
dataset. Scar set: Postn, Col1a1, Col3a1, Lox, Fap, Cthrc1, Fn1. Injury-severity set: Nppa, Nppb, Ankrd1, Xirp2, Myh7.
Cell-means linear models (age × condition × day) with pooled residual variance; induction = MI − sham at the same
age and day; one-sided t tests.

## Steps (run from `code/`)

| Script | Output (in `output/`) | Reproduces |
|---|---|---|
| `result2_scar_kinetics/01_scar_kinetics_GSE123863.py` | `result2_scar_kinetics_GSE123863.json` | Figure 2A, kinetic interaction |
| `result2_scar_kinetics/02_scar_replication_GSE198300.py` | `result2_scar_GSE198300.json` | day-3 replication |
| `result2_scar_kinetics/03_injury_severity_adjustment.py` | `result2_scar_severity_adjustment.json` | Figure 2B, severity adjustment |

## Key numbers reproduced

| Manuscript statement | Value | File in `results/` (key) | Supplementary Table S1 ID |
|---|---|---|---|
| P1 induction +0.74, +1.00, +0.78 (days 1.5, 3, 7); P8 +0.93, +1.91, +2.51 | same | `result2_scar_kinetics_GSE123863.json` (`scar_induction`) | H12c |
| Kinetic interaction +1.54, SE 0.56, p = 0.005 | +1.536, 0.556, p = 0.0054 | same file (`scar_kinetics_contrast`) | H12c |
| Day 3, independent lab: P7 +1.66 vs P1 +0.57; difference +1.09, p = 4 × 10⁻⁵ | +1.66, +0.57, +1.093, p = 3.9e-5 | `result2_scar_GSE198300.json` (`day3`) | H13b |
| Injury-severity difference +0.45, p = 0.013 (first dataset) | +0.453, p = 0.013 | `result2_scar_severity_adjustment.json` (`GSE123863.severity_contrast`) | H14a |
| Adjusted excess +0.53 (p = 0.008) and +0.74 (p = 0.003); unadjusted +0.95 and +1.09 | +0.531 (0.0082), +0.740 (0.0025); +0.946, +1.093 | same file (`scar_adjusted`, `scar_unadjusted`) | H14a |
| Slope from P1 samples only: +0.72 (p = 0.0006) and +0.82 (p = 8 × 10⁻⁵) | +0.720 (6.4e-4), +0.819 (8.1e-5) | same file (`scar_adjusted_beta_P1_only`) | H14a |
| Slope of scar on severity among MI: P1 0.02, P8 ≈ 1.13 | 0.020, 1.135 | same file (`GSE123863.ancova_MI_only`) | H14a |

Step 02 also computes the pre-specified day-21 comparison in GSE198300 (`day21`), which was not supported (same
direction, p = 0.067).

## Verification

All three scripts were re-run on the public input files; each output was identical (all keys and values) to the
corresponding file in `results/`.
