# Result 1 — Paired AP-1 sites spaced 7–60 bp mark injury-responsive enhancers across species

**Question.** Among enhancers that switch on after injury (injury-responsive enhancers, RREs) and carry at least two
canonical AP-1 sites (TGA[C/G]TCA), are AP-1–AP-1 pairs more often spaced 7–60 bp apart than in control enhancers
matched for AP-1 count and GC content? Is this specific to AP-1, and does it differ between the regenerative P1 and
the non-regenerative P8 mouse heart?

**Answer in the paper.** Yes in all four fish datasets (excess of 16–19 percentage points), robust to repeat
removal, to a local composition null and to 218 other motifs; predicted by ChromBPNet in heart and fibroblast
models; equally present at P1 and P8 in the mouse heart (Figure 1).

## Data

| Set | Accession | Notebook / input |
|---|---|---|
| Zebrafish heart (confirmatory) | PRJNA305621 | `notebooks/zebrafish_heart_h3k27ac_kang2016.ipynb` |
| Killifish fin (confirmatory) | PRJNA559885 | `notebooks/killifish_fin_h3k27ac_wang2020.ipynb` |
| Zebrafish fin (replication) | PRJNA559885 | `notebooks/zebrafish_fin_h3k27ac_wang2020.ipynb` |
| Polypterus fin (replication) | GSE315999 (authors' differential ATAC peaks) | step 06 |
| Mouse heart P1 and P8, MI vs sham, day 3 | GSE123867 / PRJNA510092 | `notebooks/mouse_heart_p1_h3k27ac.ipynb`, `..._p8_...` |
| ChromBPNet models (ENCODE, human) | Hugging Face `kundajelab/encode-chrombpnet-*` | `notebooks/chrombpnet_ap1_spacing.ipynb` |

## Steps (in order; paths relative to `code/`)

| Step | What it does | Where it runs |
|---|---|---|
| notebooks | FASTQ → fastp → Bowtie 2 → MAPQ ≥ 30, duplicates removed → MACS3 broad peaks; then step 01 (and step 02 for the mouse sets) | Kaggle, CPU |
| `01_call_injury_responsive_enhancers.py` | RREs (log2FC ≥ 1 in both replicate pairs, signal ≥ 75th percentile, > 2 kb from TSS) and control pool (\|log2FC\| < 0.3, signal ≥ median) | inside the notebooks |
| `02_ap1_spacing_test.py` | primary test: matching on AP-1 count then GC, fraction of pairs at 7–60 bp, 5,000 label permutations (seed 7) | local (zebrafish heart, killifish fin) / notebook (mouse) |
| `03_spacing_sets_and_repeat_masking.py` | same test for several sets; repeat-masked version | local |
| `04_local_composition_null.py` | excess over a local 7-mer null; difference-in-differences, pooled over the two confirmatory sets | local |
| `05_jaspar_motif_panel.py` | spacing statistic for 879 JASPAR 2024 motifs (FIMO p < 1e-4); empirical p for AP-1 | local |
| `06_polypterus_prepare_units.py` | Polypterus RREs/controls from the authors' table; then step 03 | local |
| `07_mouse_p1_p8_summary.py` | combines the P1 and P8 notebook outputs | local |
| `08_spacing_vs_induction.py` | are spaced RREs more strongly induced? (stratified permutation, 10,000) | local |
| `09_chrombpnet_cooperativity.py` | in silico two-site spacing sweep, C(d) = Δ2(d) − 2Δ1 | Kaggle, GPU |

`run_local.sh` (run from `code/`) lists the local commands with the reference downloads. `lib/spacing_core.py`
and `lib/spacing_followup.py` hold the shared code (steps 02 and 03 are thin wrappers around them). Set names used
on the command line (`zebrafish_heart`, `killifish_fin`, `zebrafish_fin`, `polypterus`) become the keys of the
output files.

## Key numbers reproduced

| Manuscript statement | Value in `results/` | File (key) | Supplementary Table S1 ID |
|---|---|---|---|
| Zebrafish heart 29.4% vs 11.8%, 34 pairs, p = 0.043 | 0.294 vs 0.118, 34, p = 0.0428 | `result1_zebrafish_heart_spacing.json` (`ap1_pair_spacing_test`) | H2 |
| Killifish fin 28.7% vs 12.6%, 87 pairs, p = 0.008 | 0.287 vs 0.126, 87, p = 0.0080 | `result1_killifish_fin_spacing.json` (`ap1_pair_spacing_test`) | H2 |
| Zebrafish fin 28.0% vs 11.0%, p = 0.0002 | 0.280 vs 0.110, p = 0.0002 | `result1_zebrafish_fin_spacing.json` (`zebrafish_fin.all_ap1_sites`) | H2s |
| Polypterus 31.3% vs 11.9%, p = 0.008 | 0.313 vs 0.119, p = 0.0080 | `result1_polypterus_spacing.json` (`polypterus.all_ap1_sites`) | D1 |
| Repeat-removed: zebrafish heart 31% vs 10%, p = 0.031; killifish 33% vs 16%, p = 0.020 | 0.310 vs 0.103, p = 0.031; 0.329 vs 0.164, p = 0.020 | `result1_repeat_masked_zebrafish_heart_killifish_fin.json` (`*.repeat_masked`) | H2r |
| Local composition null: difference-in-differences 0.183, p = 0.001 | 0.183, p = 0.0010 | `result1_local_composition_null.json` (`pooled_confirmatory`) | H3a |
| 2 of 219 other motifs ≥ AP-1, empirical p = 0.014 | 2 / 219, p = 0.0136 | `result1_jaspar_motif_panel.json` (`motif_specificity_test`) | H3b |
| ChromBPNet heart LV p = 4 × 10⁻⁷, IMR-90 p = 1 × 10⁻¹⁰; shuffled controls not significant | 4.0e-7 and 1.5e-10; controls p ≥ 0.12 | `result1_chrombpnet_cooperativity.json` | H4'a |
| Mouse P1 26.5% vs 10.3%, p = 0.006; P8 28.6% vs 13.0%, p = 0.002 | 0.265 vs 0.103, p = 0.0058; 0.286 vs 0.130, p = 0.0022 | `result1_mouse_heart_p1_p8_spacing.json` | H5a-1, H5a-2 |
| Spaced RREs more strongly induced: +0.081, p = 0.013 | +0.081, p = 0.013 | `result1_spacing_vs_induction.json` (`pooled`) | H6a |

ChromBPNet: a third model (foreskin keratinocyte) was part of the plan. There the AP-1 test was also significant
(p = 6 × 10⁻¹²) but so were all three shuffled-motif controls, so the model was judged non-specific; its results are
in the same file. The paper reports the heart and IMR-90 models.

## Verification

Steps 02–08 were re-run with these scripts on the original inputs; every output was identical (all keys and values)
to the corresponding file in `results/` (JASPAR panel: about 18,700 values). The Kaggle notebooks (alignment, peak
calling, ChromBPNet) were not re-run. Two notes on how the original runs differed: the zebrafish heart and killifish
fin spacing tests were run locally on the notebook's `rre.bed` / `control_pool.bed` (the in-notebook run stopped on
an unrelated import), and for the zebrafish fin set the enhancer call with median-signal controls was run locally on
the notebook's BAM and peak files; the code is the same.
