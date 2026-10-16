"""AP-1 spacing test for one enhancer set (the primary test of Result 1).

Units: 1 kb centred on each enhancer. Canonical AP-1 sites: regex TGA[CG]TCA (overlapping matches allowed).
Every injury-responsive enhancer (RRE) with >= 2 AP-1 sites is matched, without replacement, to a control
enhancer with exactly the same number of AP-1 sites and the closest GC content (RREs visited in random order,
seed 7). Statistic: fraction of AP-1-AP-1 pairs (start-to-start distance) spaced 7-60 bp, RRE minus control.
Null: permutation of group labels at the enhancer level (5,000 permutations, one-sided). Secondary: KS test,
sub-ranges 7-30 and 30-60 bp.

Usage (step 02): python spacing_core.py --rre RRE.bed --ctrl CONTROL_POOL.bed --genome-2bit GENOME.2bit
                 --genome-sizes GENOME.sizes --out spacing.json
Input BED columns (from 01_call_injury_responsive_enhancers.py): chrom start end name lfc inj tss_dist.
Output: {"ap1_pair_spacing_test": {...}}.
Requires UCSC twoBitToFa on PATH.
"""
import argparse
import json
import os
import re
import subprocess

import numpy as np
import pandas as pd
from scipy.stats import ks_2samp

AP1 = re.compile(r"(?=(TGA[CG]TCA))")


def fetch_ap1_stream(bed, twobit, tmp):
    """Memory-light: stream the FASTA and keep only AP-1 positions and GC content per unit."""
    bed[["chrom", "start", "end", "name"]].to_csv(tmp + ".bed", sep="\t", header=False, index=False)
    subprocess.run(f"twoBitToFa -bed={tmp}.bed {twobit} {tmp}.fa", shell=True, check=True)
    out, name, buf = {}, None, []

    def flush():
        if name is not None:
            sq = "".join(buf)
            n_acgt = sum(sq.count(b) for b in "ACGT")
            out[name] = (ap1(sq), (sq.count("G") + sq.count("C")) / max(1, n_acgt))

    with open(tmp + ".fa") as fh:
        for l in fh:
            l = l.strip()
            if l.startswith(">"):
                flush()
                name, buf = l[1:].split(":")[0], []
            else:
                buf.append(l.upper())
    flush()
    os.remove(tmp + ".fa")
    return out


def ap1(s):
    return [m.start() for m in AP1.finditer(s)]


def pairs(pos):
    return [b - a for i, a in enumerate(pos) for b in pos[i + 1:]]


def frac_close(ds, lo=7, hi=60):
    ds = np.asarray(ds)
    return float(((ds >= lo) & (ds <= hi)).mean()) if len(ds) else np.nan


def l2_test(rre_pos, ctl_pos, nperm, seed):
    """rre_pos/ctl_pos: list of lists of AP-1 positions per unit. Permutation of labels at the unit level."""
    units = rre_pos + ctl_pos
    lab = np.array([1] * len(rre_pos) + [0] * len(ctl_pos), bool)
    up = [pairs(p) for p in units]

    def stat(l):
        a = [d for u, x in zip(up, l) if x for d in u]
        b = [d for u, x in zip(up, l) if not x for d in u]
        return frac_close(a) - frac_close(b), a, b

    obs, a, b = stat(lab)
    rng = np.random.default_rng(seed)
    cnt = sum(stat(rng.permutation(lab))[0] >= obs for _ in range(nperm))
    return {"frac_7_60_rre": frac_close(a), "frac_7_60_ctl": frac_close(b), "diff": obs,
            "sens_7_30": [frac_close(a, 7, 30), frac_close(b, 7, 30)],
            "sens_30_60": [frac_close(a, 30, 60), frac_close(b, 30, 60)],
            "perm_p_one_sided": (cnt + 1) / (nperm + 1), "ks_p": float(ks_2samp(a, b).pvalue),
            "n_units": [len(rre_pos), len(ctl_pos)], "n_pairs": [len(a), len(b)]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rre", required=True)
    ap.add_argument("--ctrl", required=True)
    ap.add_argument("--genome-2bit", dest="genome_2bit", required=True)
    ap.add_argument("--genome-sizes", dest="genome_sizes", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--nperm", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    tmp = a.out + ".tmp"
    res = {}

    # controls matched on AP-1 count from the whole control pool
    cols = ["chrom", "start", "end", "name", "lfc", "inj", "tss_dist"]
    rre = pd.read_csv(a.rre, sep="\t", header=None, names=cols)
    pool = pd.read_csv(a.ctrl, sep="\t", header=None, names=cols)
    for t_ in (rre, pool):
        mid = (t_.start + t_.end) // 2
        t_["start"] = (mid - 500).clip(lower=0)
        t_["end"] = t_.start + 1000
    zsz = pd.read_csv(a.genome_sizes, sep="\t", header=None, names=["chrom", "len"]).set_index("chrom").len
    rre = rre[rre.end <= rre.chrom.map(zsz)].reset_index(drop=True)
    pool = pool[pool.end <= pool.chrom.map(zsz)].reset_index(drop=True)
    fr = fetch_ap1_stream(rre, a.genome_2bit, tmp + "_zr")
    fp = fetch_ap1_stream(pool, a.genome_2bit, tmp + "_zp")
    rre["pos"], rre["gc"] = rre.name.map(lambda n: fr[n][0]), rre.name.map(lambda n: fr[n][1])
    pool["pos"], pool["gc"] = pool.name.map(lambda n: fp[n][0]), pool.name.map(lambda n: fp[n][1])
    del fr, fp
    rre["k"], pool["k"] = rre.pos.map(len), pool.pos.map(len)
    R = rre[rre.k >= 2].reset_index(drop=True)
    P = pool[pool.k >= 2].reset_index(drop=True)
    used = np.zeros(len(P), bool)
    ctl_pos, rre_pos, unmatched = [], [], 0
    for r in R.sample(frac=1, random_state=a.seed).itertuples():
        cand = np.where((P.k.values == r.k) & ~used)[0]
        if len(cand) == 0:
            unmatched += 1
            continue
        j = cand[np.argmin(np.abs(P.gc.values[cand] - r.gc))]
        used[j] = True
        rre_pos.append(r.pos)
        ctl_pos.append(P.pos[j])
    l2 = l2_test(rre_pos, ctl_pos, a.nperm, a.seed)
    l2.update({"rre_units_ge2": int(len(R)), "pool_units_ge2": int(len(P)), "unmatched": unmatched,
               "matched_fraction": (len(R) - unmatched) / max(len(R), 1)})
    res["ap1_pair_spacing_test"] = l2
    json.dump(res, open(a.out, "w"), indent=2, default=float)
    print(json.dumps(res, indent=2, default=float))


if __name__ == "__main__":
    main()
