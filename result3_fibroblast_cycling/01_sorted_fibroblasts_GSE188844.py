"""Result 3, step 1: response per fibroblast in sorted Tcf21-lineage fibroblasts (GSE188844, bulk RNA-seq).

P2 vs P10; apical resection (AR) vs sham, 3 days after surgery; n = 3 per group (plus uninjured baseline).
Score = mean z-score of log2(CPM + 1) across all samples; genes with mean CPM < 1 removed. I = AR - sham.
  Scar response per fibroblast: I_SCAR(P10) - I_SCAR(P2) > 0   -> reported in the paper.
  Cell-cycle response: I_CYC(P2) - I_CYC(P10) > 0               -> not supported (reversed direction); this
  pattern motivated the single-cell test in step 2.

Usage: python 01_sorted_fibroblasts_GSE188844.py [--counts data/GSE188844/GSE188844_deSeq2_counts.txt.gz]
       [--out output/result3_sorted_fibroblasts_GSE188844.json]
Reproduces results/result3_sorted_fibroblasts_GSE188844.json.
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "shared"))
from bulk_rnaseq import cellmeans_contrast, resolve  # noqa: E402

SETS = {
    "SCAR": "Postn Col1a1 Col3a1 Lox Fap Cthrc1 Fn1",
    "CYC": "Mki67 Top2a Ccnb1 Cdk1 Birc5",
}


def scores(cnt):
    """cnt: genes x samples. Returns z-score set scores, the genes used and the gene-level z matrix."""
    cpm = cnt / cnt.sum() * 1e6
    lg = np.log2(cpm + 1)
    keep = cpm.mean(1) >= 1
    z = lg.sub(lg.mean(1), axis=0).div(lg.std(1).replace(0, np.nan), axis=0)
    sc, used = {}, {}
    for s, genes in SETS.items():
        g = [x for x in resolve(genes, cnt.index) if keep.get(x, False)]
        used[s] = g
        sc[s] = z.loc[g].mean(0)
    return sc, used, z


def label_188844(c):
    """Column names of the GEO count table; from the GSM descriptions, column '21' = P2 AR #1, '101' = P10 AR #1."""
    if c in ("21", "101"):
        return ("P2" if c == "21" else "P10") + "|AR"
    age = "P10" if c.startswith("10") else "P2"
    cond = {"a": "AR", "s": "sham", "u": "base"}[c[len(age) - 1]]
    return f"{age}|{cond}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", default=os.path.join(ROOT, "data", "GSE188844", "GSE188844_deSeq2_counts.txt.gz"))
    ap.add_argument("--out", default=os.path.join(ROOT, "output", "result3_sorted_fibroblasts_GSE188844.json"))
    a = ap.parse_args()
    t = pd.read_csv(a.counts, sep="\t", index_col=0)
    t = t.groupby(level=0).sum()
    integer = bool(np.allclose(t.values, np.round(t.values)))
    grp = pd.Series({c: label_188844(c) for c in t.columns})
    sc, used, z = scores(t)
    res = {}
    I = lambda age, s=1: {f"{age}|AR": s, f"{age}|sham": -s}  # noqa: E731
    wb = {**I("P10"), **I("P2", -1)}
    b = cellmeans_contrast(sc["SCAR"], grp, wb)
    keep = [c for c in t.columns if c not in ("21", "101")]
    b_s = cellmeans_contrast(sc["SCAR"][keep], grp[keep], wb)
    per = {x: cellmeans_contrast(z.loc[x], grp, wb)["est"] for x in used["SCAR"]}
    res["scar_response"] = {"contrast_P10_minus_P2": b,
                            "induction": {ag: cellmeans_contrast(sc["SCAR"], grp, I(ag))["est"] for ag in ("P2", "P10")},
                            "sensitivity_without_samples_21_101": b_s, "per_gene": per,
                            "verdict": "supported" if (b["est"] > 0 and b["p"] < 0.05) else "not supported"}
    wc = {**I("P2"), **I("P10", -1)}
    cc = cellmeans_contrast(sc["CYC"], grp, wc)
    res["cell_cycle_response"] = {"contrast_P2_minus_P10": cc,
                                  "induction": {ag: cellmeans_contrast(sc["CYC"], grp, I(ag))["est"] for ag in ("P2", "P10")},
                                  "sensitivity_without_samples_21_101": cellmeans_contrast(sc["CYC"][keep], grp[keep], wc),
                                  "verdict": "supported" if (cc["est"] > 0 and cc["p"] < 0.05) else "not supported"}
    res["dataset"] = {"n_per_group": grp.value_counts().to_dict(), "integer_counts": integer,
                      "genes": {"scar": used["SCAR"], "cell_cycle": used["CYC"]},
                      "group_means": {"scar": sc["SCAR"].groupby(grp).mean().to_dict(),
                                      "cell_cycle": sc["CYC"].groupby(grp).mean().to_dict()}}
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump(res, open(a.out, "w"), indent=1, default=float)
    print(json.dumps({"scar_induction": res["scar_response"]["induction"], "scar_contrast_P10_minus_P2": b,
                      "cell_cycle_induction": res["cell_cycle_response"]["induction"]}, indent=1, default=float))


if __name__ == "__main__":
    main()
