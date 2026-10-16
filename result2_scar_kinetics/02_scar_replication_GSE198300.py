"""Result 2, step 2: scar-program induction at P1 and P7 in an independent laboratory (GSE198300, bulk RNA-seq).

Day 3 after MI (n = 4 per group; reported in the paper) and day 21 (n = 3 per group; pre-specified test, not
supported: same direction, p = 0.067).
Score = mean z-score of log2(CPM + 1) across all samples; genes with mean CPM < 1 removed.
Contrasts from a cell-means model with pooled residual variance, one-sided t tests.

Usage: python 02_scar_replication_GSE198300.py [--dir data/GSE198300] [--out output/result2_scar_GSE198300.json]
Reproduces results/result2_scar_GSE198300.json.
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "shared"))
from bulk_rnaseq import G198, cellmeans_contrast, load_198300  # noqa: E402

SCAR = "Postn Col1a1 Col3a1 Lox Fap Cthrc1 Fn1".split()


def zscore_set(lg, genes, floor_mask):
    g = [x for x in genes if x in lg.index and floor_mask.get(x, False)]
    z = lg.loc[g].sub(lg.loc[g].mean(1), axis=0).div(lg.loc[g].std(1).replace(0, np.nan), axis=0)
    return z.mean(0), g


def scar_198300(d):
    cnt, _ = load_198300(d)
    cpm = cnt / cnt.sum() * 1e6
    lg = np.log2(cpm + 1)
    floor = (cpm.mean(1) >= 1).to_dict()
    grp = pd.Series({c: G198[c] for c in cnt.columns})
    out = {"n_samples": int(cnt.shape[1]), "groups": grp.value_counts().to_dict()}
    sc, used = zscore_set(lg, SCAR, floor)
    out["scar_genes"] = used
    for day in ("3", "21"):
        w = {f"P7|{day}|MI": 1, f"P7|{day}|Sham": -1, f"P1|{day}|MI": -1, f"P1|{day}|Sham": 1}
        sub = grp[grp.str.contains(f"\\|{day}\\|")]
        out[f"day{day}"] = {"P7_minus_P1": cellmeans_contrast(sc[sub.index], sub, w)}
        for age in ("P1", "P7"):
            out[f"day{day}"][f"induction_{age}"] = cellmeans_contrast(sc[sub.index], sub,
                                                                      {f"{age}|{day}|MI": 1, f"{age}|{day}|Sham": -1})
    c = out["day21"]["P7_minus_P1"]
    out["day21_verdict"] = "supported" if (c["est"] > 0 and c["p"] < 0.05) else "not supported"
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=os.path.join(ROOT, "data", "GSE198300"))
    ap.add_argument("--out", default=os.path.join(ROOT, "output", "result2_scar_GSE198300.json"))
    a = ap.parse_args()
    res = scar_198300(a.dir)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump(res, open(a.out, "w"), indent=1)
    print(json.dumps({d: {k: {"est": round(x["est"], 3), "p": x["p"]} for k, x in res[d].items()}
                      for d in ("day3", "day21")}, indent=1))


if __name__ == "__main__":
    main()
