"""Result 1, step 06: RREs and control pool for Polypterus senegalus from the authors' differential ATAC-seq
peak table (GSE315999_Ps_ATAC_seq_diff_peaks.txt.gz).

RRE (primary): FDR < 0.05, Fold > 0, > 2 kb from TSS, Conc_3dpa >= 75th percentile (of all peaks).
Sensitivity (a): Conc_3dpa >= median; sensitivity (b): FDR < 0.05 and Fold >= 1.
Control: |Fold| < 0.3, FDR >= 0.05, Conc_3dpa and Conc_Uninjured >= median, > 2 kb from TSS.
Each 1 kb unit is taken directly from the genome (soft-mask preserved) and stored as a "pseudo-genome" with one
sequence per unit, so that step 03 can be used unchanged (BED: unit_name 0 1000).

Usage: python 06_polypterus_prepare_units.py --table GSE315999_Ps_ATAC_seq_diff_peaks.txt.gz
       --genome GCF_016835505.1_ASM1683550v1_genomic.fna.gz --out DIR
Writes rre.bed, rre_median.bed, rre_fold1.bed, control_pool.bed, units.2bit, units.sizes, prep_summary.json.
Requires UCSC faToTwoBit on PATH.
"""
import argparse
import gzip
import json
import os
import subprocess

import pandas as pd

COLS = ["chrom", "start", "end", "name", "lfc", "inj", "tss_dist"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    d = pd.read_csv(a.table, sep="\t", low_memory=False)
    q75, med3, med0 = d.Conc_3dpa.quantile(0.75), d.Conc_3dpa.median(), d.Conc_Uninjured.median()  # thresholds from all peaks
    d = d[d["Distance to TSS"].notna()].copy()
    d["tss_dist"] = d["Distance to TSS"].abs()
    distal = d.tss_dist > 2000
    up = (d.FDR < 0.05) & (d.Fold > 0) & distal
    sets = {
        "rre": up & (d.Conc_3dpa >= q75),                                     # primary
        "rre_median": up & (d.Conc_3dpa >= med3),                             # sensitivity (a)
        "rre_fold1": (d.FDR < 0.05) & (d.Fold >= 1) & distal,                 # sensitivity (b)
        "control_pool": (d.Fold.abs() < 0.3) & (d.FDR >= 0.05) & (d.Conc_3dpa >= med3)
                        & (d.Conc_Uninjured >= med0) & distal,
    }
    mid = (d.start + d.end) // 2
    d["u0"] = (mid - 500).clip(lower=1)  # GEO table uses 1-based coordinates
    d["uname"] = [f"u{p}" for p in d.Peak]
    need = d[pd.concat(sets.values(), axis=1).any(axis=1)]
    by_chrom = {c: g for c, g in need.groupby("seqnames")}

    seqs, chrom, buf = {}, None, []

    def flush():
        if chrom in by_chrom:
            s = "".join(buf)
            for r in by_chrom[chrom].itertuples():
                w = s[r.u0 - 1: r.u0 - 1 + 1000]
                if len(w) == 1000:
                    seqs[r.uname] = w
    with gzip.open(a.genome, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                chrom, buf = line[1:].split()[0], []
            else:
                buf.append(line.rstrip())
    flush()

    fa = os.path.join(a.out, "units.fa")
    with open(fa, "w") as fh:
        for n, s in seqs.items():
            fh.write(f">{n}\n{s}\n")
    subprocess.run(f"faToTwoBit {fa} {os.path.join(a.out, 'units.2bit')}", shell=True, check=True)
    os.remove(fa)
    with open(os.path.join(a.out, "units.sizes"), "w") as fh:
        fh.writelines(f"{n}\t1000\n" for n in seqs)
    summ = {}
    for nm, m in sets.items():
        t = d[m & d.uname.isin(seqs)]
        out = pd.DataFrame({"chrom": t.uname, "start": 0, "end": 1000, "name": t.uname,
                            "lfc": t.Fold, "inj": t.Conc_3dpa, "tss_dist": t.tss_dist})
        out[COLS].to_csv(os.path.join(a.out, f"{nm}.bed"), sep="\t", header=False, index=False)
        summ[nm] = int(len(out))
    summ["units_fetched"] = len(seqs)
    json.dump(summ, open(os.path.join(a.out, "prep_summary.json"), "w"), indent=1)
    print(json.dumps(summ))


if __name__ == "__main__":
    main()
