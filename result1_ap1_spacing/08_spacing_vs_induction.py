"""Result 1, step 08: are RREs with AP-1 pairs spaced 7-60 bp more strongly induced?

Within RREs with >= 2 AP-1 sites: percentile of the injury log2FC within its set, "spaced" (>= 1 pair at 7-60 bp)
vs "not spaced"; permutation (10,000) within strata of set x number of AP-1 sites (capped at 4).
The pooled test uses the primary sets; sets listed with --desc are reported per set only.

Usage: python 08_spacing_vs_induction.py --set NAME RRE.bed GENOME.2bit GENOME.sizes [--set ...] [--desc NAME ...]
       --out result.json
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from spacing_core import pairs  # noqa: E402
from spacing_followup import AP1_HALF, fetch_case, gc, load_units, motif_regex, sites  # noqa: E402


def units_table(name, bed, twobit, sizes, out, rx):
    t = load_units(bed, sizes)
    up = [s.upper() for s in fetch_case(t, twobit, out + f".{name}")]
    t["pct"] = t.lfc.rank(pct=True)
    pos = [sites(u, u, rx, False) for u in up]
    t["k"] = [len(p) for p in pos]
    t["n_close"] = [sum(7 <= d <= 60 for d in pairs(p)) for p in pos]
    t["gc"] = [gc(u) for u in up]
    t["set"] = name
    t = t[t.k >= 2].copy()
    t["spaced"] = t.n_close > 0
    t["stratum"] = name + ":" + t.k.clip(upper=4).astype(str)
    return t


def diff(t, lab):
    return t.pct.values[lab].mean() - t.pct.values[~lab].mean()


def perm_test(t, nperm, seed):
    rng = np.random.default_rng(seed)
    lab = t.spaced.values
    obs = diff(t, lab)
    idx = [np.where(t.stratum.values == s)[0] for s in t.stratum.unique()]
    ge = 0
    for _ in range(nperm):
        l = lab.copy()
        for i in idx:
            l[i] = rng.permutation(lab[i])
        ge += diff(t, l) >= obs
    return float(obs), float((ge + 1) / (nperm + 1))


def summarize(t):
    s, n = t[t.spaced], t[~t.spaced]
    rho = [spearmanr(g.pct, g.n_close).statistic for _, g in t.groupby("stratum") if g.n_close.nunique() > 1 and len(g) >= 5]
    return {"n_spaced": int(len(s)), "n_not": int(len(n)),
            "mean_pct": [float(s.pct.mean()), float(n.pct.mean())],
            "median_lfc": [float(s.lfc.median()), float(n.lfc.median())],
            "mean_gc": [float(s.gc.mean()), float(n.gc.mean())],
            "k_dist_spaced": s.k.clip(upper=4).value_counts().sort_index().to_dict(),
            "k_dist_not": n.k.clip(upper=4).value_counts().sort_index().to_dict(),
            "mean_within_stratum_spearman": float(np.nanmean(rho)) if rho else None}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", nargs=4, action="append", required=True, metavar=("NAME", "RRE", "TWOBIT", "SIZES"))
    ap.add_argument("--desc", nargs="*", default=[], help="sets reported descriptively only (not in the pooled test)")
    ap.add_argument("--nperm", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rx = motif_regex(AP1_HALF)
    T = {n: units_table(n, b, tb, sz, a.out, rx) for n, b, tb, sz in a.set}
    res = {"per_set": {}}
    for n, t in T.items():
        obs, p = perm_test(t, a.nperm, a.seed) if t.spaced.nunique() == 2 else (None, None)
        res["per_set"][n] = {**summarize(t), "diff_pct": obs, "perm_p": p, "primary": n not in a.desc}
        print(n, json.dumps(res["per_set"][n]), flush=True)
    P = pd.concat([t for n, t in T.items() if n not in a.desc], ignore_index=True)
    obs, p = perm_test(P, a.nperm, a.seed)
    res["pooled"] = {**summarize(P), "sets": [n for n in T if n not in a.desc], "diff_pct": obs, "perm_p": p,
                  "supported": p < 0.05, "low_power": int((~P.spaced).sum()) < 20}
    pd.concat(T.values()).drop(columns=["inj"]).to_csv(a.out.replace(".json", "_units.tsv"), sep="\t", index=False)
    json.dump(res, open(a.out, "w"), indent=1)
    print(json.dumps(res["pooled"], indent=1))


if __name__ == "__main__":
    main()
