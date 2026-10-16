"""Result 2, step 3: does the excess scar induction in older hearts survive adjustment for injury severity?

Datasets: GSE123863 (P1 vs P8, days 1.5/3/7) and GSE198300 day 3 (P1 vs P7).
Severity = transcriptional injury score DMG (cardiomyocyte stress genes Nppa, Nppb, Ankrd1, Xirp2, Myh7).
Adjusted score A = SCAR - beta * DMG, beta = OLS slope across all samples of the dataset (primary);
sensitivity: beta estimated from P1 samples only. Test: I_A(older) - I_A(P1) > 0 in both datasets,
averaged over days for GSE123863. Descriptive ANCOVA: SCAR ~ 1 + older + DMG + older:DMG.

Usage: python 03_injury_severity_adjustment.py [--counts123863 ...] [--dir198300 ...]
       [--out output/result2_severity_adjustment.json]
Reproduces results/result2_scar_severity_adjustment.json.
"""
import argparse
import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "shared"))
from bulk_rnaseq import cellmeans_contrast, load_123863, load_198300, resolve  # noqa: E402

SETS = {
    "SCAR": "Postn Col1a1 Col3a1 Lox Fap Cthrc1 Fn1",
    "DMG": "Nppa Nppb Ankrd1 Xirp2 Myh7",
}


def scores(cnt):
    """cnt: genes x samples. Returns z-score set scores and the genes used."""
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


def beta(y, x):
    return float(np.polyfit(x.values, y.values, 1)[0])


def ancova(sc, grp, old):
    """Descriptive: SCAR ~ 1 + older + DMG + older:DMG."""
    o = grp.str.startswith(old).astype(float).values
    d = sc["DMG"].values
    X = np.column_stack([np.ones_like(d), o, d, o * d])
    coef = np.linalg.lstsq(X, sc["SCAR"].values, rcond=None)[0]
    return dict(zip(["intercept", "old", "slope_DMG", "old_x_DMG"], map(float, coef)))


def adjust_dataset(cnt, grp, old, young, days, ctrl):
    sc, used, _ = scores(cnt)
    out = {"genes": {"scar": used["SCAR"], "injury_severity": used["DMG"]}}
    sub = grp.index[grp.str.contains("|".join(f"\\|{d}\\|" for d in days))]
    out["induction_severity"] = {f"{a}|{d}": cellmeans_contrast(sc["DMG"][sub], grp[sub], {f"{a}|{d}|MI": 1, f"{a}|{d}|{ctrl}": -1})["est"]
                    for a in (young, old) for d in days}
    out["induction_scar"] = {f"{a}|{d}": cellmeans_contrast(sc["SCAR"][sub], grp[sub], {f"{a}|{d}|MI": 1, f"{a}|{d}|{ctrl}": -1})["est"]
                     for a in (young, old) for d in days}

    def w():
        k = 1 / len(days)
        r = {}
        for d in days:
            r.update({f"{old}|{d}|MI": k, f"{old}|{d}|{ctrl}": -k, f"{young}|{d}|MI": -k, f"{young}|{d}|{ctrl}": k})
        return r
    out["severity_contrast"] = cellmeans_contrast(sc["DMG"][sub], grp[sub], w())
    out["scar_unadjusted"] = cellmeans_contrast(sc["SCAR"][sub], grp[sub], w())
    b_all = beta(sc["SCAR"], sc["DMG"])
    yi = grp.index[grp.str.startswith(young)]
    b_y = beta(sc["SCAR"][yi], sc["DMG"][yi])
    out["beta_all_samples"], out["beta_P1_only"] = b_all, b_y
    out["scar_adjusted"] = cellmeans_contrast((sc["SCAR"] - b_all * sc["DMG"])[sub], grp[sub], w())
    out["scar_adjusted_beta_P1_only"] = cellmeans_contrast((sc["SCAR"] - b_y * sc["DMG"])[sub], grp[sub], w())
    mi = grp.index[grp.str.endswith("|MI")]
    out["ancova_all_samples"] = ancova({k: sc[k] for k in ("SCAR", "DMG")}, grp, old)
    out["ancova_MI_only"] = ancova({k: sc[k][mi] for k in ("SCAR", "DMG")}, grp[mi], old)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts123863", default=os.path.join(ROOT, "data", "GSE123863", "GSE123863_RNASeq_RawCountTable.txt.gz"))
    ap.add_argument("--dir198300", default=os.path.join(ROOT, "data", "GSE198300"))
    ap.add_argument("--out", default=os.path.join(ROOT, "output", "result2_scar_severity_adjustment.json"))
    a = ap.parse_args()
    c1, g1 = load_123863(a.counts123863)
    c2, g2 = load_198300(a.dir198300)
    r1 = adjust_dataset(c1, g1, "P8", "P1", ["1.5", "3", "7"], "Sham")
    r2 = adjust_dataset(c2, g2, "P7", "P1", ["3"], "Sham")
    ok = [r["scar_adjusted"]["est"] > 0 and r["scar_adjusted"]["p"] < 0.05 for r in (r1, r2)]
    res = {"GSE123863": r1, "GSE198300_day3": r2,
           "verdict": "supported" if all(ok) else ("partially supported" if any(ok) else "not supported")}
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump(res, open(a.out, "w"), indent=1, default=float)
    print(json.dumps({ds: {k: {"est": round(r[k]["est"], 3), "p": r[k]["p"]}
                           for k in ("severity_contrast", "scar_unadjusted", "scar_adjusted", "scar_adjusted_beta_P1_only")}
                      for ds, r in (("GSE123863", r1), ("GSE198300_day3", r2))}, indent=1))


if __name__ == "__main__":
    main()
