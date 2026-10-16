"""AP-1 spacing test across several enhancer sets, with and without AP-1 sites in repeat sequence (Result 1).

For each set the primary test of spacing_core.py is reproduced exactly ("all_ap1_sites"), and repeated after
removing AP-1 sites that touch soft-masked (repeat) sequence, with re-matching ("repeat_masked").
Also reports how many closely spaced RRE pairs touch repeat sequence ("close_rre_pairs_in_repeat": [n_repeat, n_close]).
Shared helpers (unit loading, sequence retrieval, site calling, matching) are imported by steps 04, 05 and 08.

Usage (step 03): python spacing_followup.py --set NAME RRE.bed CONTROL.bed GENOME.2bit GENOME.sizes [--set ...]
                 --out result.json
Requires UCSC twoBitToFa on PATH.
"""
import argparse
import json
import os
import re
import subprocess

import numpy as np
import pandas as pd

from spacing_core import frac_close, l2_test, pairs

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
AP1_HALF = "TGA"


def motif_regex(half):
    """Palindromic 7-mer half[CG]revcomp(half); for AP1_HALF = 'TGA' this is TGA[CG]TCA."""
    rc = "".join(COMP[b] for b in reversed(half))
    return re.compile(f"(?=({half}[CG]{rc}))")


def load_units(bed, sizes):
    cols = ["chrom", "start", "end", "name", "lfc", "inj", "tss_dist"]
    t = pd.read_csv(bed, sep="\t", header=None, names=cols)
    mid = (t.start + t.end) // 2
    t["start"] = (mid - 500).clip(lower=0)
    t["end"] = t.start + 1000
    zsz = pd.read_csv(sizes, sep="\t", header=None, names=["chrom", "len"]).set_index("chrom").len
    return t[t.end <= t.chrom.map(zsz)].reset_index(drop=True)


def fetch_case(t, twobit, tmp):
    """Sequence per unit with lower case (soft-mask) preserved."""
    t[["chrom", "start", "end", "name"]].to_csv(tmp + ".bed", sep="\t", header=False, index=False)
    subprocess.run(f"twoBitToFa -bed={tmp}.bed {twobit} {tmp}.fa", shell=True, check=True)
    d, n, buf = {}, None, []
    with open(tmp + ".fa") as fh:
        for l in fh:
            l = l.strip()
            if l.startswith(">"):
                if n is not None:
                    d[n] = "".join(buf)
                n, buf = l[1:].split(":")[0], []
            else:
                buf.append(l)
    d[n] = "".join(buf)
    os.remove(tmp + ".fa")
    os.remove(tmp + ".bed")
    return [d[x] for x in t.name]


def sites(seq, up, rx, mask):
    p = [m.start() for m in rx.finditer(up)]
    if mask:
        p = [x for x in p if seq[x:x + 7].isupper()]
    return p


def gc(up):
    n = sum(up.count(b) for b in "ACGT")
    return (up.count("G") + up.count("C")) / max(1, n)


def match(rre_pos, rre_gc, pool_pos, pool_gc, seed):
    """As in spacing_core.py: RREs (k >= 2) in random order, control with the same k then closest GC,
    without replacement."""
    R = pd.DataFrame({"pos": rre_pos, "gc": rre_gc})
    R["k"] = R.pos.map(len)
    R = R[R.k >= 2].reset_index(drop=True)
    P = pd.DataFrame({"pos": pool_pos, "gc": pool_gc})
    P["k"] = P.pos.map(len)
    P = P[P.k >= 2].reset_index(drop=True)
    used = np.zeros(len(P), bool)
    rp, cp, unmatched = [], [], 0
    for r in R.sample(frac=1, random_state=seed).itertuples():
        cand = np.where((P.k.values == r.k) & ~used)[0]
        if len(cand) == 0:
            unmatched += 1
            continue
        j = cand[np.argmin(np.abs(P.gc.values[cand] - r.gc))]
        used[j] = True
        rp.append(r.pos)
        cp.append(P.pos[j])
    return rp, cp, {"rre_units_ge2": int(len(R)), "pool_units_ge2": int(len(P)), "unmatched": unmatched,
                    "matched_fraction": (len(R) - unmatched) / max(len(R), 1)}


def quick_diff(rp, cp):
    a = [d for u in rp for d in pairs(u)]
    b = [d for u in cp for d in pairs(u)]
    return frac_close(a) - frac_close(b), len(a), len(b)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", nargs=5, action="append", required=True,
                    metavar=("NAME", "RRE", "CONTROL", "TWOBIT", "SIZES"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--nperm", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    res = {}
    for name, rbed, cbed, twobit, sizes in a.set:
        rre, pool = load_units(rbed, sizes), load_units(cbed, sizes)
        sr = fetch_case(rre, twobit, a.out + f".{name}_r")
        sp = fetch_case(pool, twobit, a.out + f".{name}_p")
        ur, upl = [s.upper() for s in sr], [s.upper() for s in sp]
        gr, gp = [gc(s) for s in ur], [gc(s) for s in upl]
        out = {}
        ap1 = motif_regex(AP1_HALF)
        # primary test (must equal spacing_core.py) and repeat-masked version
        for tag, mask in (("all_ap1_sites", False), ("repeat_masked", True)):
            rp, cp, info = match([sites(s, u, ap1, mask) for s, u in zip(sr, ur)], gr,
                                 [sites(s, u, ap1, mask) for s, u in zip(sp, upl)], gp, a.seed)
            out[tag] = {**l2_test(rp, cp, a.nperm, a.seed), **info}
        # closely spaced RRE pairs that touch repeat sequence
        n_close = n_rep = 0
        for s, u in zip(sr, ur):
            p = sites(s, u, ap1, False)
            for i, x in enumerate(p):
                for y in p[i + 1:]:
                    if 7 <= y - x <= 60:
                        n_close += 1
                        n_rep += not (s[x:x + 7].isupper() and s[y:y + 7].isupper())
        out["close_rre_pairs_in_repeat"] = [n_rep, n_close]
        res[name] = out
        print(name, json.dumps({k: out[k] for k in ("all_ap1_sites", "repeat_masked")}, default=float), flush=True)
        del sr, sp, ur, upl
    json.dump(res, open(a.out, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
