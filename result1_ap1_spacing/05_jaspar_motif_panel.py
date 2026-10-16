"""Result 1, step 05: is the spacing excess specific to AP-1? The spacing statistic for every JASPAR 2024 CORE
vertebrate non-redundant motif.

FIMO p < 1e-4, both strands; hits on opposite strands whose centres differ by <= 1 bp are merged.
Matching as in the primary test (number of sites, then GC). AP-1-family motifs (consensus containing
TGA[CG]TCA) are excluded. A motif is eligible if >= 15 matched RRE units in both confirmatory sets
(sets named zebrafish_heart and killifish_fin). Empirical p = (1 + number of eligible motifs with mean difference
>= AP-1) / (1 + number eligible), where the AP-1 value is the regex-based mean difference of the two confirmatory
sets from step 02.

Usage: python 05_jaspar_motif_panel.py --meme JASPAR2024_CORE_vertebrates_non-redundant.meme
       --set NAME RRE CONTROL TWOBIT SIZES [...] --workdir DIR --out result.json
Requires MEME suite FIMO (5.5.x) and UCSC twoBitToFa on PATH.
"""
import argparse
import concurrent.futures as cf
import json
import os
import re
import subprocess
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from spacing_followup import fetch_case, gc, load_units, match, quick_diff  # noqa: E402

AP1_RE = re.compile("TGA[CG]TCA")
AP1_REGEX_MEAN_DIFF = 0.16869506423258956  # step 02: mean difference of the two confirmatory sets (regex AP-1)
CONF = ("zebrafish_heart", "killifish_fin")


def parse_meme(path):
    """-> header (text before the first MOTIF), list of (id, alt name, text block, consensus)."""
    txt = open(path).read()
    head, *blocks = re.split(r"(?m)^(?=MOTIF )", txt)
    out = []
    for b in blocks:
        f = b.split("\n")[0].split()
        rows = [list(map(float, l.split())) for l in b.split("\n")
                if re.match(r"^\s*[0-9.]+\s+[0-9.]+\s+[0-9.]+\s+[0-9.]+\s*$", l)]
        cons = "".join("ACGT"[int(np.argmax(r))] for r in rows)
        out.append((f[1], f[2] if len(f) > 2 else f[1], b, cons))
    return head, out


def run_fimo(args):
    meme, fa, out = args
    if not os.path.exists(out):
        subprocess.run(f"fimo --thresh 1e-4 --text --skip-matched-sequence --verbosity 1 {meme} {fa} > {out}.part"
                       f" && mv {out}.part {out}", shell=True, check=True)
    return out


def load_hits(paths):
    """-> {motif: {sequence: [positions]}} with the strand-merging rule."""
    res = defaultdict(dict)
    for p in paths:
        df = pd.read_csv(p, sep="\t", usecols=[0, 2, 3, 4, 5], names=["motif", "seq", "start", "stop", "strand"],
                         header=0, dtype={"motif": "category", "seq": "category", "strand": "category"})
        if df.empty:
            continue
        df["c"] = (df.start + df.stop) // 2
        df = df.sort_values(["motif", "seq", "c"], ignore_index=True)
        # drop hits on the opposite strand of the previous row whose centre differs by <= 1 bp
        same = (df.motif.values[1:] == df.motif.values[:-1]) & (df.seq.values[1:] == df.seq.values[:-1])
        dup = same & (df.strand.values[1:] != df.strand.values[:-1]) & (np.abs(np.diff(df.c.values)) <= 1)
        df = df[~np.concatenate([[False], dup])]
        for (m, s), v in df.groupby(["motif", "seq"], observed=True).c:
            res[m][s] = v.tolist()
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meme", required=True)
    ap.add_argument("--set", nargs=5, action="append", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--chunks", type=int, default=16)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    os.makedirs(a.workdir, exist_ok=True)
    head, motifs = parse_meme(a.meme)
    ap1fam = {m for m, _, _, c in motifs if AP1_RE.search(c)}
    alt = {m: al for m, al, _, _ in motifs}
    # split the motif file into chunks for parallel FIMO runs
    chunk_files = []
    for i in range(a.chunks):
        fn = os.path.join(a.workdir, f"motifs_{i:02d}.meme")
        with open(fn, "w") as fh:
            fh.write(head + "".join(b for _, _, b, _ in motifs[i::a.chunks]))
        chunk_files.append(fn)
    res = {"n_motifs": len(motifs), "ap1_family_excluded": sorted(f"{m} {alt[m]}" for m in ap1fam), "sets": {}}
    diffs = {}
    for name, rbed, cbed, twobit, sizes in a.set:
        rre, pool = load_units(rbed, sizes), load_units(cbed, sizes)
        ur = [s.upper() for s in fetch_case(rre, twobit, os.path.join(a.workdir, f"{name}_r"))]
        upl = [s.upper() for s in fetch_case(pool, twobit, os.path.join(a.workdir, f"{name}_p"))]
        fa = os.path.join(a.workdir, f"{name}.fa")
        with open(fa, "w") as fh:
            for i, s in enumerate(ur):
                fh.write(f">R{i}\n{s}\n")
            for j, s in enumerate(upl):
                fh.write(f">P{j}\n{s}\n")
        gr, gp = [gc(s) for s in ur], [gc(s) for s in upl]
        jobs = [(m, fa, os.path.join(a.workdir, f"{name}_fimo_{i:02d}.tsv")) for i, m in enumerate(chunk_files)]
        with cf.ThreadPoolExecutor(a.jobs) as ex:
            outs = list(ex.map(run_fimo, jobs))
        hits = load_hits(outs)
        d = {}
        for m, _, _, _ in motifs:
            h = hits.get(m, {})
            rp = [sorted(h.get(f"R{i}", [])) for i in range(len(ur))]
            pp = [sorted(h.get(f"P{j}", [])) for j in range(len(upl))]
            r_pos, c_pos, info = match(rp, gr, pp, gp, a.seed)
            df_, na, nb = quick_diff(r_pos, c_pos)
            d[m] = {"alt": alt[m], "diff": df_, "n_matched": len(r_pos), "n_pairs": [na, nb],
                    "matched_fraction": info["matched_fraction"], "rre_units_ge2": info["rre_units_ge2"]}
        diffs[name] = d
        res["sets"][name] = {"n_rre_units": len(ur), "n_pool_units": len(upl)}
        print(name, "done", flush=True)
        del hits

    if set(CONF) <= set(diffs):
        elig = [m for m, _, _, _ in motifs if m not in ap1fam
                and all(diffs[s][m]["n_matched"] >= 15 and not np.isnan(diffs[s][m]["diff"]) for s in CONF)]
        md = {m: float(np.mean([diffs[s][m]["diff"] for s in CONF])) for m in elig}
        top = sorted(md.items(), key=lambda kv: -kv[1])
        res["motif_specificity_test"] = {
            "ap1_regex_mean_diff": AP1_REGEX_MEAN_DIFF, "n_eligible": len(elig),
            "n_ge_ap1": int(sum(v >= AP1_REGEX_MEAN_DIFF for v in md.values())),
            "p_empirical": (1 + sum(v >= AP1_REGEX_MEAN_DIFF for v in md.values())) / (1 + len(elig)),
            "decidable": len(elig) >= 19,
            "null_median_diff": float(np.median(list(md.values()))) if md else None,
            "top20": [{"motif": m, "alt": alt[m], "mean_diff": v,
                       CONF[0]: diffs[CONF[0]][m]["diff"], CONF[1]: diffs[CONF[1]][m]["diff"],
                       "n": [diffs[CONF[0]][m]["n_matched"], diffs[CONF[1]][m]["n_matched"]]} for m, v in top[:20]]}
        # descriptive: AP-1 as JASPAR PWMs
        res["ap1_as_jaspar_pwm"] = {
            f"{m} {alt[m]}": {s: diffs[s][m] for s in diffs} for m in ap1fam if alt[m] in ("FOS::JUN", "FOSL1::JUN", "JUN")}
    res["per_set_motif_differences"] = diffs
    json.dump(res, open(a.out, "w"), indent=1, default=float)
    print(json.dumps({k: res[k] for k in ("motif_specificity_test",) if k in res}, indent=1, default=float))


if __name__ == "__main__":
    main()
