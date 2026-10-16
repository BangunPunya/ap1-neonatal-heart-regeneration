"""Result 1, step 07: combine the mouse heart P1 and P8 spacing tests (GSE123867, H3K27ac, MI vs sham, day 3).

Input: the output folders of the two mouse Kaggle notebooks (each with spacing.json and summary.json).
  P1 test: spacing excess in P1 RREs, one-sided p < 0.05.
  P8 (descriptive): is the pattern specific to the regenerative age? (P1 significant, P8 not, P1 > P8)

Usage: python 07_mouse_p1_p8_summary.py --p1 DIR_P1 --p8 DIR_P8 --out mouse_p1_p8.json
"""
import argparse
import json
import os


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--p1", required=True)
    ap.add_argument("--p8", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--spacing-file", dest="spacing_file", default="spacing.json",
                    help="name of the step-02 output inside each folder")
    a = ap.parse_args()
    res = {}
    for age, d in (("P1", a.p1), ("P8", a.p8)):
        f = json.load(open(os.path.join(d, a.spacing_file)))["ap1_pair_spacing_test"]
        s = json.load(open(os.path.join(d, "summary.json")))
        res[f"spacing_{age}"] = {**f, "rre": s["rre"], "control_pool": s["control_pool"],
                                 "rep_concordance": s["rep_concordance_spearman"]}
    p1 = res["spacing_P1"]
    res["p1_test"] = {"p": p1["perm_p_one_sided"], "supported": p1["perm_p_one_sided"] < 0.05,
                   "low_power": p1["rre_units_ge2"] < 20}
    p8 = res["spacing_P8"]
    res["regenerative_age_specificity"] = {"P1_diff": p1["diff"], "P8_diff": p8["diff"], "P8_p": p8["perm_p_one_sided"],
                               "regen_age_pattern": bool(p1["perm_p_one_sided"] < 0.05 and p8["perm_p_one_sided"] >= 0.05
                                                         and p1["diff"] > p8["diff"])}
    json.dump(res, open(a.out, "w"), indent=1)
    print(json.dumps({k: res[k] for k in ("p1_test", "regenerative_age_specificity")}, indent=1))


if __name__ == "__main__":
    main()
