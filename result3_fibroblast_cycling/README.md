# Result 3 — Injury recruits fibroblasts into the cell cycle more strongly at P7/P8

**Question.** Does the larger scar program at P7/P8 come from a stronger response per fibroblast, or from more
fibroblasts entering the response (cell cycle)?

**Answer in the paper.** Not per fibroblast: in sorted fibroblasts the scar response is not larger at P10 than at P2.
But the effect of injury on fibroblast cycling differs by age: positive age × injury interaction in all five
single-cell units, pooled log odds ratio +0.58; driven by low baseline cycling at P7/P8 that rises with injury
(Figure 3).

## Data

- GSE188844 (Misra 2021): bulk RNA-seq of sorted Tcf21-lineage fibroblasts, P2 and P10, apical resection vs sham,
  n = 3 (`download_bulk_data.sh`).
- GSE153480 (Wang Z 2020): scRNA-seq, P1 and P8, MI vs sham, days 1 and 3 (one library per condition).
- GSE247626 (Fan 2025): scRNA-seq of non-myocytes, P1 and P7, apical resection vs sham, days 1, 4 and 7.

## Steps (paths relative to `code/`)

| Step | What it does | Where it runs |
|---|---|---|
| `result3_fibroblast_cycling/01_sorted_fibroblasts_GSE188844.py` | scar response per fibroblast (AR − sham), P10 vs P2; cell-cycle response | local |
| `result3_fibroblast_cycling/notebooks/fibroblast_cell_cycle_scrna.ipynb` | QC (≥ 200 genes, < 20% mitochondrial reads), normalisation, Leiden clustering; fibroblast clusters by ≥ 2 of Pdgfra, Tcf21, Dpt, Gsn, Smoc2 in ≥ 50% of cells; S/G2M phase (Tirosh genes); per-unit interaction log odds ratio (0.5 continuity correction, Woolf variance); DerSimonian–Laird pooling | Kaggle, CPU |

## Key numbers reproduced

| Manuscript statement | Value | File in `results/` (key) | Supplementary Table S1 ID |
|---|---|---|---|
| Scar induction per fibroblast P10 +0.62 vs P2 +0.85; p = 0.67 | +0.616 vs +0.846, p = 0.67 | `result3_sorted_fibroblasts_GSE188844.json` (`scar_response`) | H14b |
| Interaction positive in every unit (+0.36 to +0.84); pooled +0.58, SE 0.099, p = 2 × 10⁻⁹, τ² = 0.017 | +0.36…+0.84; +0.580, 0.099, 2.3e-9, 0.017; 5/5 positive | `result3_fibroblast_cycling_scrna.json` (`fibroblast_cycling_meta_analysis.units`, `.pooled`) | H15c |
| P8 cycling 8.5% → 19.1% (day 1) and 10.9% → 26.0% (day 3) | same | same file (`units.GSE153480_d*.cyc_rate`) | H15c |
| Robust to Mki67/Top2a definition (+0.57, p = 0.0002) and to author annotations (+0.57) | +0.569, p = 0.00023; +0.565 | same file (`sensitivity_Mki67_Top2a`, `sensitivity_author_annotation_GSE247626`) | H15c |

Step 01 also computes the cell-cycle score in sorted fibroblasts (`cell_cycle_response`). It was pre-specified in the
opposite direction (P2 more than P10) and was not supported; the reversed pattern motivated the single-cell test.

## Verification

Step 01 was re-run on the public input file; its output was identical (all keys and values) to
`results/result3_sorted_fibroblasts_GSE188844.json`. The single-cell notebook was not re-run; it is the fibroblast
cell-cycle part of the Kaggle notebook that produced `results/result3_fibroblast_cycling_scrna.json`, with code for
analyses not reported in the paper removed and English comments and output keys.
