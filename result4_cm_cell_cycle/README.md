# Result 4 — The cardiomyocyte cell-cycle program is maintained through P7

**Question.** Within cardiomyocytes, do cell-cycle transcripts already fall between P1 and P7 (when the
regenerative window closes), or only later?

**Answer in the paper.** They do not decline between P1 and P7 (equivalent within pre-specified bounds) and collapse
only by P21, confirming the original report; scar and cell-cycle module scores both peak at P7 (descriptive;
Figure 4).

## Data

GSE155658 (Li Z 2022): bulk RNA-seq of isolated cardiomyocytes at P1, P7 and P21 from three mouse strains (C57, KM,
NIH; n = 6 per age, 2 per strain). Entrez IDs are mapped to symbols with the NCBI mouse gene_info table. Download
with `bash download_bulk_data.sh` (from `code/`). This dataset had not been opened before the test was fixed.

## Method

Counts → CPM → log2(CPM + 1); genes kept if maximum CPM ≥ 5; z-scores within strain, then pooled.
Cell-cycle panel: Mki67, Top2a, Ccnb1, Cdk1, Birc5, Aurkb, Plk1, Anln, Ect2. Equivalence of the mean log2 fold
change P7/P1 by two one-sided tests with bounds −0.5 and +1.0, plus a pre-specified decline of ≥ 2 log2 units from
P7 to P21.

## Step (run from `code/`)

`python result4_cm_cell_cycle/01_cardiomyocyte_cell_cycle_GSE155658.py` →
`output/result4_cardiomyocyte_cell_cycle_GSE155658.json` and `output/result4_module_trajectories.json`.

## Key numbers reproduced

| Manuscript statement | Value | File in `results/` (key) | Supplementary Table S1 ID |
|---|---|---|---|
| Mean log2FC of nine genes P7/P1 +0.65; equivalent, TOST p = 0.0003 | +0.651; p = 0.00031 | `result4_cardiomyocyte_cell_cycle_GSE155658.json` (`cell_cycle_panel.equivalence_test_TOST`) | H18c |
| Collapse by P21: mean −3.39 log2 units | −3.388 | same file (`cell_cycle_panel.mean_log2FC_p21_p7`) | H18c |
| Scar and cell-cycle module scores both peak at P7 (descriptive) | cell cycle +0.43 / +0.80 / −1.23; scar +0.40 / +0.82 / −1.22 | `result4_module_trajectories.json` (`module_trajectories`) | — |
| ECM–integrin module drawn in Figure 4B | rises P1 → P7, then flat | `result4_cardiomyocyte_cell_cycle_GSE155658.json` (`ecm_integrin_module`) | H18b |

## Verification

The script was re-run on the public input files; its outputs were identical (all keys and values) to the two files
in `results/`. The NCBI gene_info table changes over time; the original run used the version downloaded on
2026-10-14.
