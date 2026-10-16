"""Result 1, step 04: does close AP-1 spacing in RREs exceed what the local sequence composition predicts?

Local null per unit = the other 7-mers (+ strand) occurring exactly k times in the same unit (k = number of AP-1
sites), excluding 7-mers with internal period <= 3. Excess = observed fraction of pairs at 7-60 bp minus the local
null. Primary statistic = difference-in-differences (excess in RREs minus excess in matched controls), pooled over
the two confirmatory sets (zebrafish heart + killifish fin); permutation of labels within set.

Usage: python 04_local_composition_null.py --set NAME RRE CONTROL TWOBIT SIZES [...]
       --primary zebrafish_heart killifish_fin --out result.json
Output: {"sets": {name: per-set statistics}, "pooled_confirmatory": {...}, "per_set_permutation": {name: {...}}}.
"""
import argparse
import json
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import binomtest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from spacing_followup import AP1_HALF, fetch_case, gc, load_units, motif_regex, sites  # noqa: E402

AP1_SEQS = {"TGACTCA", "TGAGTCA"}


def low_period(s):
    return any(s[p:] == s[:-p] for p in (1, 2, 3))


def local_null(up, k):
    occ = defaultdict(list)
    for i in range(len(up) - 6):
        occ[up[i:i + 7]].append(i)
    fr = []
    for km, pos in occ.items():
        if len(pos) != k or "N" in km or km in AP1_SEQS or low_period(km):
            continue
        ds = [b - a for i, a in enumerate(pos) for b in pos[i + 1:]]
        fr.append(np.mean([(7 <= d <= 60) for d in ds]))
    return (float(np.mean(fr)), len(fr)) if fr else (np.nan, 0)


def match_idx(rre_pos, rre_gc, pool_pos, pool_gc, seed):
    """Same as spacing_followup.match, but returns unit indices."""
    R = pd.DataFrame({"pos": rre_pos, "gc": rre_gc, "i": range(len(rre_pos))})
    R["k"] = R.pos.map(len)
    R = R[R.k >= 2].reset_index(drop=True)
    P = pd.DataFrame({"pos": pool_pos, "gc": pool_gc, "i": range(len(pool_pos))})
    P["k"] = P.pos.map(len)
    P = P[P.k >= 2].reset_index(drop=True)
    used = np.zeros(len(P), bool)
    out = []
    for r in R.sample(frac=1, random_state=seed).itertuples():
        cand = np.where((P.k.values == r.k) & ~used)[0]
        if len(cand) == 0:
            continue
        j = cand[np.argmin(np.abs(P.gc.values[cand] - r.gc))]
        used[j] = True
        out.append((r.i, int(P.i[j])))
    return out


def unit_rec(pos, up):
    ds = [b - a for i, a in enumerate(pos) for b in pos[i + 1:]]
    nf, nk = local_null(up, len(pos))
    return {"w": len(ds), "close": sum(7 <= d <= 60 for d in ds), "null": nf, "n_null_kmers": nk}


def excess(units):
    u = [x for x in units if not np.isnan(x["null"])]
    w = np.array([x["w"] for x in u], float)
    if w.sum() == 0:
        return np.nan
    return float((np.array([x["close"] for x in u]).sum() - (w * np.array([x["null"] for x in u])).sum()) / w.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", nargs=5, action="append", required=True)
    ap.add_argument("--primary", nargs="+", default=["zebrafish_heart", "killifish_fin"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--nperm", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    rx = motif_regex(AP1_HALF)
    res, data = {"sets": {}}, {}
    for name, rbed, cbed, twobit, sizes in a.set:
        rre, pool = load_units(rbed, sizes), load_units(cbed, sizes)
        sr, sp = fetch_case(rre, twobit, a.out + f".{name}_r"), fetch_case(pool, twobit, a.out + f".{name}_p")
        ur, upl = [s.upper() for s in sr], [s.upper() for s in sp]
        pr, pp = [sites(s, u, rx, False) for s, u in zip(sr, ur)], [sites(s, u, rx, False) for s, u in zip(sp, upl)]
        pairs_idx = match_idx(pr, [gc(u) for u in ur], pp, [gc(u) for u in upl], a.seed)
        R = [unit_rec(pr[i], ur[i]) for i, _ in pairs_idx]
        C = [unit_rec(pp[j], upl[j]) for _, j in pairs_idx]
        data[name] = R + C, np.array([1] * len(R) + [0] * len(C), bool)
        # check: reproduces the primary spacing fractions (RRE / control)
        fr = sum(x["close"] for x in R) / max(1, sum(x["w"] for x in R))
        fc = sum(x["close"] for x in C) / max(1, sum(x["w"] for x in C))
        ex_r = [x["close"] / x["w"] - x["null"] for x in R if not np.isnan(x["null"])]
        st = binomtest(sum(e > 0 for e in ex_r), sum(e != 0 for e in ex_r), 0.5, alternative="greater") if ex_r else None
        res["sets"][name] = {"n_matched": len(R), "spacing_frac_rre": fr, "spacing_frac_ctl": fc,
                     "null_frac_rre": float(np.nanmean([x["null"] for x in R])),
                     "null_frac_ctl": float(np.nanmean([x["null"] for x in C])),
                     "excess_rre": excess(R), "excess_ctl": excess(C), "did": excess(R) - excess(C),
                     "units_without_null": [sum(np.isnan(x["null"]) for x in R), sum(np.isnan(x["null"]) for x in C)],
                     "sign_test_rre_excess": {"pos": int(sum(e > 0 for e in ex_r)), "neg": int(sum(e < 0 for e in ex_r)),
                                              "p_one_sided": float(st.pvalue) if st else None}}
        print(name, json.dumps(res["sets"][name], default=float), flush=True)

    rng = np.random.default_rng(a.seed)

    def did(names, perm=False):
        R_all, C_all = [], []
        for n in names:
            units, lab = data[n]
            if perm:
                lab = rng.permutation(lab)
            R_all += [u for u, x in zip(units, lab) if x]
            C_all += [u for u, x in zip(units, lab) if not x]
        return excess(R_all) - excess(C_all)

    res["per_set_permutation"] = {}
    for key, names in (("pooled", a.primary),) + tuple((n, [n]) for n in data):
        if not all(n in data for n in names):
            continue
        obs = did(names)
        cnt = sum(did(names, True) >= obs for _ in range(a.nperm))
        r = {"did": obs, "perm_p_one_sided": (cnt + 1) / (a.nperm + 1)}
        if key == "pooled":
            res["pooled_confirmatory"] = {"sets": list(a.primary), **r}
        else:
            res["per_set_permutation"][key] = r
        print(key, r, flush=True)
    res = {k: res[k] for k in ("sets", "pooled_confirmatory", "per_set_permutation") if k in res}
    json.dump(res, open(a.out, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
