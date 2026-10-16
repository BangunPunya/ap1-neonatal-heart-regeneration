"""Result 4: cell-cycle transcripts in isolated cardiomyocytes at P1, P7 and P21 (GSE155658, Li et al. 2022).

Three mouse strains (C57, KM, NIH; 2 samples per strain and age; strain treated as a blocking factor).
Raw counts -> CPM -> log2(CPM + 1); genes kept if max CPM across samples >= 5; z-scores computed within strain
before pooling; set score = mean z across genes. P56 (C57 only, n = 2) is excluded.

  Cell-cycle panel (primary): mean log2FC(P7/P1) of 9 genes within [-0.5, +1.0] (TOST)
  and mean log2FC(P21/P7) <= -2.0.
  ECM-integrin module: module score P7 > P1 with >= 5/8 genes log2FC >= 1, then flat
  (|P21 - P7| <= 0.6); Agrn predicted |log2FC| < 0.5.
  Module trajectories (descriptive, Figure 4B): mean module score at P1, P7, P21 for SCAR, ECM, CYC, DMG.

Entrez IDs are mapped to symbols with the NCBI mouse gene_info table (data/ncbi/Mus_musculus.gene_info.gz).

Usage: python 01_cardiomyocyte_cell_cycle_GSE155658.py [--dir data/GSE155658] [--gene-info ...]
       [--out output/result4_cardiomyocyte_cell_cycle_GSE155658.json]
Writes the output JSON and, next to it, result4_module_trajectories.json.
Reproduces results/result4_cardiomyocyte_cell_cycle_GSE155658.json and results/result4_module_trajectories.json.
"""
import argparse
import gzip
import json
import os
import sys

import numpy as np
import pandas as pd
import scipy.stats as ss

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "shared"))
from bulk_rnaseq import resolve  # noqa: E402

ECM = "Vtn Itga7 Itgb6 Mfap5 Aspn Pcolce2 Col4a4 Hmcn2"
SCAR = "Postn Col1a1 Col3a1 Lox Fap Cthrc1 Fn1"
DMG = "Nppa Nppb Ankrd1 Xirp2 Myh7"
CYC = "Mki67 Top2a Ccnb1 Cdk1 Birc5 Aurkb Plk1 Anln Ect2"
AGES = ["p1", "p7", "p21"]


def entrez2symbol(ids, gene_info):
    """Entrez -> symbol via the NCBI mouse gene_info table."""
    if not os.path.exists(gene_info):
        raise SystemExit(f"need {gene_info}; download the NCBI mouse gene_info table first (see download_bulk_data.sh)")
    m = {}
    with gzip.open(gene_info, "rt", errors="replace") as f:
        hdr = f.readline().lstrip("#").rstrip("\n").split("\t")
        gi, sy = hdr.index("GeneID"), hdr.index("Symbol")
        for line in f:
            p = line.rstrip("\n").split("\t")
            m[p[gi]] = p[sy]
    return {i: m.get(str(i).strip('"')) for i in ids}


def load_155658(d, gene_info):
    """Three mouse strains; returns counts (genes x samples) + labels (strain, age)."""
    specs = [("C57", os.path.join(d, "GSE155658_RNAseq_genecountsraw.txt.gz")),
             ("KM", os.path.join(d, "GSE155658_RNAseq_genecountsraw_KM.txt.gz")),
             ("NIH", os.path.join(d, "GSE155658_RNAseq_genecountsraw_NIH.txt.gz"))]
    frames, lab = [], {}
    for strain, path in specs:
        t = pd.read_csv(path, sep="\t", index_col=0)
        t.index = [str(i).strip('"') for i in t.index]
        ren = {}
        for c in t.columns:
            s = c.strip('"').lower()
            age = next((a for a in ("p56", "p21", "p7", "p1") if f"_{a}" in s or s.endswith(a)), None)
            if age is None:
                continue
            new = f"{strain}|{age}|{len(ren)}"
            ren[c] = new
            lab[new] = (strain, age)
        frames.append(t[list(ren)].rename(columns=ren))
    cnt = pd.concat(frames, axis=1).dropna()
    sym = entrez2symbol(cnt.index, gene_info)
    cnt = cnt.loc[[i for i in cnt.index if sym.get(i)]]
    n_unmapped = len(sym) - len(cnt)
    cnt.index = [sym[i] for i in cnt.index]
    return cnt.groupby(level=0).sum(), pd.DataFrame(lab, index=["strain", "age"]).T, n_unmapped


def zwithin(lg, strain):
    """Per-gene z-score WITHIN each strain, then pooled (strain as blocking factor)."""
    out = lg.copy().astype(float)
    for s in strain.unique():
        c = strain.index[strain == s]
        sub = lg[c]
        out[c] = sub.sub(sub.mean(1), axis=0).div(sub.std(1).replace(0, np.nan), axis=0)
    return out


def setscore(z, genes):
    g = resolve(genes, z.index)
    return (z.loc[g].mean(0) if g else pd.Series(np.nan, index=z.columns)), g


def welch(a, b, alt="greater"):
    t = ss.ttest_ind(a, b, equal_var=False, alternative=alt)
    return {"mean_a": float(np.mean(a)), "mean_b": float(np.mean(b)),
            "diff": float(np.mean(a) - np.mean(b)), "t": float(t.statistic), "p": float(t.pvalue)}


def lfc(lgm, age, a, b):
    """Mean log2 fold change (a - b) on log2 CPM."""
    return float(lgm[age == a].mean() - lgm[age == b].mean())


def tost(vals, lo, hi):
    """Two one-sided tests: is the mean inside [lo, hi]?"""
    m, se, n = np.mean(vals), ss.sem(vals), len(vals)
    p1 = ss.t.sf((m - lo) / se, n - 1)
    p2 = ss.t.cdf((m - hi) / se, n - 1)
    return {"mean": float(m), "p_above_lower_bound": float(p1), "p_below_upper_bound": float(p2),
            "equivalent": bool(max(p1, p2) < 0.05)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=os.path.join(ROOT, "data", "GSE155658"))
    ap.add_argument("--gene-info", dest="gene_info", default=os.path.join(ROOT, "data", "ncbi", "Mus_musculus.gene_info.gz"))
    ap.add_argument("--out", default=os.path.join(ROOT, "output", "result4_cardiomyocyte_cell_cycle_GSE155658.json"))
    a = ap.parse_args()
    res = {}

    cnt, meta, n_unmapped = load_155658(a.dir, a.gene_info)
    cpm = cnt / cnt.sum() * 1e6
    lg = np.log2(cpm + 1)
    keep = cpm.max(1) >= 5
    strain, age = meta["strain"], meta["age"]
    use = age.index[age.isin(AGES)]
    lgk = lg.loc[keep[keep].index, use]
    z = zwithin(lgk, strain[use])
    ag = age[use]
    res["dataset"] = {"n_genes_detected": int(keep.sum()), "n_entrez_unmapped": int(n_unmapped),
                      "samples_per_age": {k: int((ag == k).sum()) for k in AGES},
                      "n_P56_samples_C57_only": int((age == "p56").sum())}

    def agemean(s):
        return s.groupby(ag).mean()

    # ---------- ECM-integrin module ----------
    b = {}
    se, ge = setscore(z, ECM)
    b["genes"] = ge
    b["score_test_P7_vs_P1"] = welch(se[ag == "p7"], se[ag == "p1"])
    per = {}
    for g in ge:
        per[g] = {"log2FC_p7_p1": lfc(lgk.loc[g], ag, "p7", "p1"),
                  "log2FC_p21_p7": lfc(lgk.loc[g], ag, "p21", "p7"),
                  "cpm_p1": float(cpm.loc[g, use][ag == "p1"].mean()),
                  "cpm_p7": float(cpm.loc[g, use][ag == "p7"].mean())}
    b["per_gene"] = per
    n_up = sum(1 for v in per.values() if v["log2FC_p7_p1"] >= 1.0)
    b["n_genes_up_ge_1_log2"] = n_up
    b["rise_test_passed"] = bool(b["score_test_P7_vs_P1"]["p"] < 0.05 and n_up >= 5)
    sk = agemean(se)
    b["score_per_age"] = {k: float(sk[k]) for k in AGES if k in sk}
    b["flatness_P21_minus_P7"] = float(sk.get("p21", np.nan) - sk.get("p7", np.nan))
    b["flatness_test_passed"] = bool(abs(b["flatness_P21_minus_P7"]) <= 0.6)
    if "Agrn" in lgk.index:
        b["agrn_check"] = {"log2FC_p7_p1": lfc(lgk.loc["Agrn"], ag, "p7", "p1"),
                           "cpm_p1": float(cpm.loc["Agrn", use][ag == "p1"].mean())}
        b["agrn_check"]["prediction_met"] = bool(abs(b["agrn_check"]["log2FC_p7_p1"]) < 0.5)
    b["verdict"] = "supported" if b["rise_test_passed"] and b["flatness_test_passed"] else "not supported"
    res["ecm_integrin_module"] = b

    # ---------- cell-cycle panel (primary) ----------
    c = {}
    gc = resolve(CYC, lgk.index)
    c["genes"] = gc
    per = {g: lfc(lgk.loc[g], ag, "p7", "p1") for g in gc}
    c["per_gene_log2FC_p7_p1"] = per
    c["equivalence_test_TOST"] = tost(list(per.values()), -0.5, 1.0)
    per21 = {g: lfc(lgk.loc[g], ag, "p21", "p7") for g in gc}
    c["per_gene_log2FC_p21_p7"] = per21
    c["mean_log2FC_p21_p7"] = float(np.mean(list(per21.values())))
    c["decline_test_passed"] = bool(c["mean_log2FC_p21_p7"] <= -2.0)
    c["verdict"] = "supported" if c["equivalence_test_TOST"]["equivalent"] and c["decline_test_passed"] else "not supported"
    res["cell_cycle_panel"] = c

    # ---------- module trajectories (descriptive) ----------
    traj = {}
    for nm, gg in (("scar", SCAR), ("ecm_integrin", ECM), ("cell_cycle", CYC), ("injury_severity", DMG)):
        s, _ = setscore(z, gg)
        m = agemean(s)
        traj[nm] = {k: float(m[k]) for k in AGES if k in m}

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump(res, open(a.out, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    json.dump({"module_trajectories": traj}, open(os.path.join(os.path.dirname(a.out), "result4_module_trajectories.json"), "w",
                                       encoding="utf-8"), indent=1, ensure_ascii=False)
    print(json.dumps({"cell_cycle_TOST": c["equivalence_test_TOST"], "mean_log2FC_P21_P7": c["mean_log2FC_p21_p7"],
                      "trajectories": traj}, indent=1))


if __name__ == "__main__":
    main()
