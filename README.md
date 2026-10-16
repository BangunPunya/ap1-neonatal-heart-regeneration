# Shared AP-1 spacing grammar, amplified fibroblast response

Code and analysis outputs for the manuscript **"Shared AP-1 spacing grammar, amplified fibroblast response:
transcriptomic and epigenomic reanalysis of the closing regenerative window in the neonatal mouse heart"** by
Adityana Arta Bangun and Didik Priyandoko (Universitas Pendidikan Indonesia). The neonatal mouse heart regenerates
after injury at postnatal day 1 (P1) but forms a fibrotic scar by P7/P8. This secondary analysis of public
transcriptomic and epigenomic data asks whether the transition reflects a change in the regulatory logic of the
injury response or in its magnitude and cellular source. Each test was pre-specified before the test data were
opened; all pre-specified tests and outcomes are listed in Supplementary Table S1 of the manuscript.

## Installation

- Python 3.11 or newer, then `pip install -r requirements.txt` (numpy, pandas, scipy, matplotlib).
- Single-cell, ChIP-seq and ChromBPNet steps run as Kaggle notebooks (free account; Internet on; CPU, or GPU T4 for
  ChromBPNet). The notebooks install their own tools (fastp, Bowtie 2, samtools, bedtools, MACS3, Scanpy,
  TensorFlow).
- The local follow-up steps of Result 1 need Linux or WSL with UCSC `twoBitToFa`/`faToTwoBit` and MEME suite FIMO
  5.5.x.
- On Windows, set `PYTHONUTF8=1` before running the scripts.

## Data

All data are public and are not stored in this repository. The scripts download them into `code/data/` (ignored
by git); the bulk count tables for Results 2–4 come from `code/download_bulk_data.sh`.

| Accession | Source | Organism / tissue | Design | Assay | Used for |
|---|---|---|---|---|---|
| PRJNA305621 (GSE75894) | Kang 2016 | Zebrafish heart | Regenerating vs uninjured; n = 2 vs 2 | H3K27ac ChIP-seq | AP-1 spacing (confirmatory) |
| PRJNA559885 | Wang W 2020 | Killifish fin; zebrafish fin | 0 vs 1 dpa; n = 2 vs 2 | H3K27ac ChIP-seq | AP-1 spacing (confirmatory; replication) |
| GSE315999 | F Sousa 2026 | Polypterus fin | Intact vs 3 dpa; n = 3 vs 3 | ATAC-seq (author peaks) | AP-1 spacing (replication) |
| GSE123867 | Wang Z 2019 | Mouse heart | P1, P8; MI vs sham, day 3; n = 2 | H3K27ac ChIP-seq | AP-1 spacing at P1 and P8 |
| GSE123863 | Wang Z 2019 | Mouse ventricle | P1, P8; MI vs sham; days 1.5, 3, 7; n = 3 | Bulk RNA-seq | Scar kinetics; severity adjustment |
| GSE198300 | Kankuri 2023 | Mouse LV (CD-1) | P1, P7; MI vs sham; days 3, 21; n = 3–4 | Bulk RNA-seq | Replication; severity adjustment |
| GSE153480 | Wang Z 2020 | Mouse heart (ICR/CD1) | P1, P8; MI vs sham; days 1, 3; 1 library each | scRNA-seq | Fibroblast cycling (2 units) |
| GSE247626 | Fan 2025 | Mouse heart non-myocytes | P1, P7; resection vs sham; days 1, 4, 7 | scRNA-seq | Fibroblast cycling (3 units) |
| GSE188844 | Misra 2021 | Sorted Tcf21+ fibroblasts (C57BL/6J) | P2, P10; resection vs sham; n = 3 | Bulk RNA-seq | Per-fibroblast scar response |
| GSE155658 | Li Z 2022 | Isolated cardiomyocytes (3 strains) | P1, P7, P21; n = 6 per age | Bulk RNA-seq | Cardiomyocyte cell cycle |

## Usage

Run everything from `code/`:

```bash
cd code

# Figures 1-4 and the key numbers, from the saved outputs in results/
python figures/make_figures.py

# Results 2-4 from the public count tables (minutes on a laptop); outputs go to output/
bash download_bulk_data.sh
python result2_scar_kinetics/01_scar_kinetics_GSE123863.py
python result2_scar_kinetics/02_scar_replication_GSE198300.py
python result2_scar_kinetics/03_injury_severity_adjustment.py
python result3_fibroblast_cycling/01_sorted_fibroblasts_GSE188844.py
python result4_cm_cell_cycle/01_cardiomyocyte_cell_cycle_GSE155658.py

# Result 1: run the Kaggle notebooks in result1_ap1_spacing/notebooks/, extract their outputs into
# data/result1/<set>/, then (Linux/WSL)
bash result1_ap1_spacing/run_local.sh

# Result 3 single-cell part: Kaggle notebook result3_fibroblast_cycling/notebooks/fibroblast_cell_cycle_scrna.ipynb
```

Each result folder has a README with the steps in order and the numbers each step reproduces. Re-running a script
writes a file with the same structure as the corresponding file in `results/`.

## Results

- **Paired AP-1 sites spaced 7–60 bp mark injury-responsive enhancers across species** (Figure 1). Injury-responsive
  enhancers contained more AP-1 pairs spaced 7–60 bp than matched control enhancers: 29.4% vs 11.8% in regenerating
  zebrafish heart (p = 0.043) and 28.7% vs 12.6% in killifish fin (p = 0.008), replicated in zebrafish fin (28.0% vs
  11.0%, p = 0.0002) and Polypterus fin (31.3% vs 11.9%, p = 0.008). The excess survived repeat removal and a local
  composition null (difference-in-differences 0.183, p = 0.001), was specific to AP-1 (2 of 219 other motifs,
  p = 0.014) and was predicted by ChromBPNet heart and fibroblast models. It was equally strong at P1 (26.5% vs 10.3%,
  p = 0.006) and P8 (28.6% vs 13.0%, p = 0.002).
- **The scar program is larger and faster at P7/P8 than predicted by injury severity** (Figure 2). Scar induction at
  P1 was +0.74, +1.00 and +0.78 (days 1.5, 3, 7) versus +0.93, +1.91 and +2.51 at P8 (kinetic interaction +1.54,
  p = 0.005); day 3 in an independent laboratory +1.66 at P7 vs +0.57 at P1 (p = 4 × 10⁻⁵). After adjustment for a
  transcriptional injury-severity score the excess remained +0.53 (p = 0.008) and +0.74 (p = 0.003).
- **Injury recruits fibroblasts into the cell cycle more strongly at P7/P8** (Figure 3). The scar response per sorted
  fibroblast was not larger at P10 (+0.62 vs +0.85, p = 0.67), but the age × injury interaction for cycling
  fibroblasts was positive in all five single-cell units, pooled log odds ratio +0.58 (SE 0.099, p = 2 × 10⁻⁹), driven
  by low baseline cycling at P7/P8 (for example 8.5% → 19.1% one day after MI at P8).
- **The cardiomyocyte cell-cycle program is maintained through P7** (Figure 4). In isolated cardiomyocytes nine
  cell-cycle transcripts did not decline between P1 and P7 (mean log2 fold change +0.65; equivalent, TOST p = 0.0003)
  and collapsed only by P21 (−3.39).

![Figure 1](code/results/Figure1.png)
![Figure 2](code/results/Figure2.png)
![Figure 3](code/results/Figure3.png)
![Figure 4](code/results/Figure4.png)

## Citation

Bangun AA, Priyandoko D. Shared AP-1 spacing grammar, amplified fibroblast response: transcriptomic and epigenomic
reanalysis of the closing regenerative window in the neonatal mouse heart. Manuscript in preparation.

## License

MIT License. See [LICENSE](LICENSE).
