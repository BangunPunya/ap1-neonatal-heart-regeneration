"""Shared helpers for the bulk RNA-seq analyses (Results 2-4).

Loaders for the public count tables, the GSE198300 sample annotation, gene-symbol resolution and the
cell-means linear contrast used throughout (pooled residual variance, one-sided t test).
"""
import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.stats as ss

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUTPUT = ROOT / "output"

# GSE198300 (Kankuri 2023): GSM5944529-GSM5944556 in GEO order -> "age|day|condition"
G198 = {}
for i, lab in zip(range(29, 57), ["P1|3|Sham"] * 4 + ["P1|3|MI"] * 4 + ["P7|3|Sham"] * 4 + ["P7|3|MI"] * 4
                  + ["P1|21|Sham"] * 3 + ["P1|21|MI"] * 3 + ["P7|21|Sham"] * 3 + ["P7|21|MI"] * 3):
    G198[f"GSM59445{i:02d}"] = lab


def cellmeans_contrast(score, groups, w):
    """Cell-means model with pooled residual variance.

    score: Series per sample; groups: Series of group labels; w: dict label -> weight.
    Returns estimate, SE, t, df and the one-sided p-value for the alternative contrast > 0.
    """
    m = score.groupby(groups).mean()
    n = score.groupby(groups).size()
    res = score - groups.map(m)
    df = len(score) - len(m)
    s2 = (res ** 2).sum() / df
    est = sum(w[k] * m[k] for k in w)
    se = np.sqrt(s2 * sum(w[k] ** 2 / n[k] for k in w))
    t = est / se
    return {"est": float(est), "se": float(se), "t": float(t), "df": int(df), "p": float(ss.t.sf(t, df)),
            "means": {k: float(m[k]) for k in m.index}, "n": {k: int(n[k]) for k in n.index}}


def resolve(genes, index):
    """Space-separated gene list; 'A|B' means use A, or alias B if A is absent."""
    out = []
    for g in genes.split():
        hit = [a for a in g.split("|") if a in index]
        if hit:
            out.append(hit[0])
    return out


def load_123863(path):
    """GSE123863 raw count table -> (genes x samples, Series 'age|day|condition')."""
    raw = pd.read_csv(path, sep="\t")
    cols = [c for c in raw.columns if c.startswith(("P1+", "P8+"))]
    cnt = raw.groupby("GeneName")[cols].sum()
    lab = {c: f"{c.split('+')[0]}|{c.split('+')[1].split('MI')[0].split('Sham')[0]}|{'MI' if 'MI_' in c else 'Sham'}"
           for c in cols}
    return cnt, pd.Series(lab)


def load_198300(d):
    """GSE198300 per-sample count files (extracted from GSE198300_RAW.tar into directory d)."""
    cols = {}
    for f in sorted(glob.glob(os.path.join(d, "*.counts.txt.gz"))):
        gsm = os.path.basename(f).split("_")[0]
        s = pd.read_csv(f, sep="\t").groupby("external_gene_name")["counts"].sum()
        cols[gsm] = cols[gsm].add(s, fill_value=0) if gsm in cols else s
    cnt = pd.DataFrame(cols).fillna(0)
    return cnt, pd.Series({c: G198[c] for c in cnt.columns})
