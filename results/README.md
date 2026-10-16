# Analysis outputs

The JSON files below hold the values from which the numbers in the manuscript were taken; the PNG files are
Figures 1–4. `figures/make_figures.py` draws the figures from these files and prints the key numbers. Each script
writes a file with the same structure to `output/` when it is re-run.

| File | Produced by (paths relative to `code/`) | Used for |
|---|---|---|
| `result1_zebrafish_heart_spacing.json`, `result1_killifish_fin_spacing.json` | Result 1 notebooks + `02_ap1_spacing_test.py` | Figure 1A |
| `result1_zebrafish_heart_enhancers_summary.json`, `result1_killifish_fin_enhancers_summary.json`, `result1_zebrafish_fin_enhancers_summary.json` | `01_call_injury_responsive_enhancers.py` (in the notebooks) | numbers of enhancers |
| `result1_zebrafish_fin_spacing.json` | `03_spacing_sets_and_repeat_masking.py` | Figure 1A |
| `result1_repeat_masked_zebrafish_heart_killifish_fin.json` | `03_spacing_sets_and_repeat_masking.py` | repeat-removed test |
| `result1_local_composition_null.json` | `04_local_composition_null.py` | local composition null |
| `result1_jaspar_motif_panel.json` | `05_jaspar_motif_panel.py` | motif specificity |
| `result1_polypterus_units_summary.json`, `result1_polypterus_spacing.json` | `06_polypterus_prepare_units.py`, then `03_...` | Figure 1A |
| `result1_mouse_heart_p1_p8_spacing.json` | mouse notebooks + `07_mouse_p1_p8_summary.py` | Figure 1A |
| `result1_spacing_vs_induction.json` | `08_spacing_vs_induction.py` | spacing vs induction |
| `result1_chrombpnet_cooperativity.json` | `09_chrombpnet_cooperativity.py` (Kaggle GPU notebook) | Figure 1B |
| `result2_scar_kinetics_GSE123863.json` | `result2_scar_kinetics/01_...` | Figure 2A |
| `result2_scar_GSE198300.json` | `result2_scar_kinetics/02_...` | day-3 replication |
| `result2_scar_severity_adjustment.json` | `result2_scar_kinetics/03_...` | Figure 2B |
| `result3_sorted_fibroblasts_GSE188844.json` | `result3_fibroblast_cycling/01_...` | per-fibroblast response |
| `result3_fibroblast_cycling_scrna.json` | `result3_fibroblast_cycling/notebooks/fibroblast_cell_cycle_scrna.ipynb` | Figure 3 |
| `result4_cardiomyocyte_cell_cycle_GSE155658.json` | `result4_cm_cell_cycle/01_...` | Figure 4A |
| `result4_module_trajectories.json` | `result4_cm_cell_cycle/01_...` | Figure 4B |
| `Figure1.png` … `Figure4.png` | `figures/make_figures.py` | manuscript figures |

Common fields: `est`, `se`, `t`, `df`, `p` = estimate, standard error, t statistic, degrees of freedom and one-sided
p-value of a contrast; `frac_7_60_rre` / `frac_7_60_ctl` = fraction of AP-1 pairs spaced 7–60 bp in injury-responsive
and matched control enhancers; `perm_p_one_sided` = permutation p-value.
