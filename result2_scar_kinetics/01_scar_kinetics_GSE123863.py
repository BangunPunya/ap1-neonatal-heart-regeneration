"""Result 2, step 1: scar-program kinetics after MI at P1 and P8 (GSE123863, bulk RNA-seq, n = 3 per group).

Per-sample gene-set score = mean z-score (across all 36 samples) of log2(CPM + 1); genes with mean CPM < 1 removed.
Cell-means model with 12 cells (age x day x condition) and pooled residual variance (df = 24); one-sided t tests.
Induction I(age, day) = mean(MI) - mean(Sham).
Kinetics contrast: [I(P8, 7) - I(P8, 1.5)] - [I(P1, 7) - I(P1, 1.5)] > 0.

Usage: python 01_scar_kinetics_GSE123863.py [--counts data/GSE123863/GSE123863_RNASeq_RawCountTable.txt.gz]
       [--out output/result2_scar_kinetics_GSE123863.json]
Reproduces results/result2_scar_kinetics_GSE123863.json.
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
import scipy.stats as ss

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SETS = {"SCAR": "Postn Col1a1 Col3a1 Lox Fap Cthrc1 Fn1"}
CELLS = [(a, d, c) for a in ("P1", "P8") for d in ("1.5", "3", "7") for c in ("MI", "Sham")]


def contrast(score, w):
    """w: dict (age, day, condition) -> weight; returns estimate, SE, t, df and one-sided p (> 0)."""
    m, n, ss_res = {}, {}, 0.0
    for k in CELLS:
        x = score[k]
        m[k], n[k] = x.mean(), len(x)
        ss_res += ((x - x.mean()) ** 2).sum()
    df = sum(n.values()) - len(CELLS)
    s2 = ss_res / df
    est = sum(w.get(k, 0) * m[k] for k in CELLS)
    se = np.sqrt(s2 * sum(w.get(k, 0) ** 2 / n[k] for k in CELLS))
    t = est / se
    return {"est": float(est), "se": float(se), "t": float(t), "df": int(df), "p": float(ss.t.sf(t, df))}


def I(a, d, s=1.0):
    return {(a, d, "MI"): s, (a, d, "Sham"): -s}


def add(*ws):
    out = {}
    for w in ws:
        for k, v in w.items():
            out[k] = out.get(k, 0) + v
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", default=os.path.join(ROOT, "data", "GSE123863", "GSE123863_RNASeq_RawCountTable.txt.gz"))
    ap.add_argument("--out", default=os.path.join(ROOT, "output", "result2_scar_kinetics_GSE123863.json"))
    a = ap.parse_args()
    raw = pd.read_csv(a.counts, sep="\t")
    cols = [c for c in raw.columns if c.startswith(("P1+", "P8+"))]
    cnt = raw.groupby("GeneName")[cols].sum()
    cpm = cnt / cnt.sum() * 1e6
    lg = np.log2(cpm + 1)
    z = lg.sub(lg.mean(1), axis=0).div(lg.std(1).replace(0, np.nan), axis=0)
    key = {c: (c.split("+")[0], c.split("+")[1].split("MI")[0].split("Sham")[0], "MI" if "MI_" in c else "Sham") for c in cols}
    res = {}
    sc = {}
    for s, genes in SETS.items():
        g = [x for x in genes.split() if x in cpm.index and cpm.loc[x].mean() >= 1]
        res["scar_gene_set"] = {"used": g, "dropped": [x for x in genes.split() if x not in g]}
        v = z.loc[g].mean(0)
        sc[s] = {k: v[[c for c in cols if key[c] == k]].values for k in CELLS}
        res["scar_induction"] = {f"{ag}|{d}": contrast(sc[s], I(ag, d)) for ag in ("P1", "P8") for d in ("1.5", "3", "7")}
    # kinetics contrast
    cc = contrast(sc["SCAR"], add(I("P8", "7"), I("P8", "1.5", -1), I("P1", "7", -1), I("P1", "1.5")))
    res["scar_kinetics_contrast"] = {"contrast": cc, "verdict": "supported" if (cc["est"] > 0 and cc["p"] < 0.05) else "not supported"}
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump(res, open(a.out, "w"), indent=1)
    print(json.dumps({"scar_induction": {k: round(v["est"], 3) for k, v in res["scar_induction"].items()},
                      "kinetics_contrast": cc}, indent=1))


if __name__ == "__main__":
    main()
