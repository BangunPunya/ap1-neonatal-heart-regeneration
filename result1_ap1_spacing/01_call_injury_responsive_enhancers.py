"""Result 1, step 01: injury-responsive enhancers (RREs) and the control pool from H3K27ac ChIP-seq.

Input : MACS3 broadPeak files per sample, filtered BAM per sample (uninjured R1 R2, injured R1 R2), TSS BED.
Consensus peaks = merged peaks of all samples; signal = CPM per kb per replicate.
RRE    : log2FC(injured/uninjured) >= 1 in BOTH replicate pairs (R1-R1, R2-R2), injured signal >= 75th percentile,
         > 2 kb from the nearest TSS.
Control: |mean log2FC| < 0.3 and signal >= the --ctl-quantile quantile at both time points, > 2 kb from TSS.
         The paper uses --ctl-quantile 0.5 (median) for every ChIP-seq set.
Output : rre.bed, control_pool.bed (chrom start end name lfc inj tss_dist), consensus.bed, summary.json.

Run inside the Kaggle notebooks in notebooks/ (zebrafish heart, killifish fin, mouse P1/P8) and on the
zebrafish fin BAM files (see README). Requires bedtools and samtools on PATH.
"""
import argparse
import io
import json
import os
import subprocess

import numpy as np
import pandas as pd


def sh(cmd):
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True).stdout


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--peaks", nargs="+", required=True)
    ap.add_argument("--bam0", nargs=2, required=True, help="uninjured BAM, replicate 1 and 2")
    ap.add_argument("--bam1", nargs=2, required=True, help="injured BAM, replicate 1 and 2")
    ap.add_argument("--tss", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ctl-quantile", dest="ctl_q", type=float, default=0.75,
                    help="control signal threshold (the paper uses 0.5)")
    ap.add_argument("--chrom-regex", dest="chrom_re", default=r"^chr[0-9]+\b",
                    help="chromosomes used (UCSC: ^chr[0-9]+; NCBI: ^NC_)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    cons = os.path.join(a.out, "consensus.bed")
    sh(f"cat {' '.join(a.peaks)} | cut -f1-3 | grep -E '{a.chrom_re}' | sort -k1,1 -k2,2n | bedtools merge -i - > {cons}")
    bams = a.bam0 + a.bam1
    cov = sh(f"bedtools multicov -bams {' '.join(bams)} -bed {cons}")
    df = pd.read_csv(io.StringIO(cov), sep="\t", header=None,
                     names=["chrom", "start", "end", "c0a", "c0b", "c1a", "c1b"])
    lib = {c: float(sh(f"samtools view -c {b}").strip()) for c, b in zip(["c0a", "c0b", "c1a", "c1b"], bams)}
    for c in lib:
        df[c] = df[c] / lib[c] * 1e6 / ((df.end - df.start) / 1000)  # CPM per kb
    ps = 0.1
    df["lfc_a"] = np.log2((df.c1a + ps) / (df.c0a + ps))
    df["lfc_b"] = np.log2((df.c1b + ps) / (df.c0b + ps))
    df["lfc"] = (df.lfc_a + df.lfc_b) / 2
    df["inj"] = (df.c1a + df.c1b) / 2
    df["uni"] = (df.c0a + df.c0b) / 2
    p75_1, p75_0 = df.inj.quantile(0.75), df.uni.quantile(0.75)
    q1, q0 = df.inj.quantile(a.ctl_q), df.uni.quantile(a.ctl_q)

    tss = pd.read_csv(a.tss, sep="\t", header=None, usecols=[0, 1], names=["chrom", "pos"])
    df["tss_dist"] = np.inf
    for c, g in df.groupby("chrom"):
        t = np.sort(tss.loc[tss.chrom == c, "pos"].values)
        if len(t) < 2:
            continue
        mid = ((g.start + g.end) // 2).values
        j = np.searchsorted(t, mid).clip(1, len(t) - 1)
        df.loc[g.index, "tss_dist"] = np.minimum(np.abs(mid - t[j - 1]), np.abs(mid - t[j]))
    distal = df.tss_dist > 2000
    rre = df[(df.lfc_a >= 1) & (df.lfc_b >= 1) & (df.inj >= p75_1) & distal].copy()
    ctl = df[(df.lfc.abs() < 0.3) & (df.inj >= q1) & (df.uni >= q0) & distal].copy()
    rre["name"] = [f"RRE_{i}" for i in range(len(rre))]
    ctl["name"] = [f"CTL_{i}" for i in range(len(ctl))]
    cols = ["chrom", "start", "end", "name", "lfc", "inj", "tss_dist"]
    rre[cols].to_csv(os.path.join(a.out, "rre.bed"), sep="\t", header=False, index=False)
    ctl[cols].to_csv(os.path.join(a.out, "control_pool.bed"), sep="\t", header=False, index=False)
    s = {"consensus_peaks": len(df), "rre": len(rre), "control_pool": len(ctl), "lib_sizes": lib,
         "p75_inj": float(p75_1), "p75_uni": float(p75_0), "ctl_quantile": a.ctl_q,
         "rep_concordance_spearman": float(df[["lfc_a", "lfc_b"]].corr(method="spearman").iloc[0, 1])}
    json.dump(s, open(os.path.join(a.out, "summary.json"), "w"), indent=2)
    print(json.dumps(s, indent=2))


if __name__ == "__main__":
    main()
