"""Result 1, step 09: in silico spacing sweep of two AP-1 sites with pre-trained ENCODE ChromBPNet models
(bias-free "nobias" models, folds 0-4, averaged).

Backgrounds: 100 random hg38 windows (2,114 bp, GC 0.35-0.50, < 50% soft-masked, no N), dinucleotide-shuffled
until free of AP-1 sites (seed 7). One site (TGACTCA) is inserted at position 1,050, the second at 1,050 + d
(d = 7-100 bp). C(d) = [f(pair at d) - f(background)] - 2 [f(single) - f(background)], f = predicted log counts.
Test: one-sided Wilcoxon signed-rank on mean C(7-60) - mean C(61-100) over backgrounds.
Negative controls: three shuffled 7-mers without AP-1 (neg1-neg3) analysed identically.
Models: foreskin keratinocyte, IMR-90 fibroblast, heart left ventricle. The paper reports the heart and IMR-90
models; in the keratinocyte model the negative controls also passed, so it was judged non-specific (see README).

Usage: python 09_chrombpnet_cooperativity.py --hg38 hg38.2bit --sizes hg38.chrom.sizes --models MODEL_CACHE_DIR --out DIR
       (--predict-only with --folds splits the prediction work across GPUs; a final call without it computes the
       statistics). Writes DIR/cooperativity.json. Requires TensorFlow, huggingface_hub and UCSC twoBitToFa.
"""
import argparse
import json
import os
import re
import subprocess

import numpy as np
import pandas as pd
from scipy import stats

SEQLEN, POS1 = 2114, 1050
DISTS = list(range(7, 101))
AP1_RE = re.compile("TGA[CG]TCA")
INSERTS = {"AP1": "TGACTCA", "neg1": "TCAGTGA", "neg2": "ATCGACT", "neg3": "GTACTCA"}
MODELS = {
    "keratinocyte": "kundajelab/encode-chrombpnet-ATAC-foreskin-keratinocyte-ENCSR356KRQ-ENCSR280GJX",
    "IMR90_fibroblast": "kundajelab/encode-chrombpnet-ATAC-IMR-90-ENCSR200OML-ENCSR978WIX",
    "heart_LV": "kundajelab/encode-chrombpnet-ATAC-heart-left-ventricle-ENCSR925LGW-ENCSR923QPG",
}


def dinuc_shuffle(seq, rng):
    """Dinucleotide-preserving shuffle (Eulerian algorithm, as in deeplift)."""
    arr = np.frombuffer(seq.encode(), dtype=np.int8)
    chars, tokens = np.unique(arr, return_inverse=True)
    nxt = [np.where(tokens[:-1] == t)[0] + 1 for t in range(len(chars))]
    for t in range(len(chars)):
        if len(nxt[t]) > 1:
            idx = np.arange(len(nxt[t]))
            idx[:-1] = rng.permutation(len(idx) - 1)
            nxt[t] = nxt[t][idx]
    cnt = [0] * len(chars)
    ind, res = 0, np.empty_like(tokens)
    res[0] = tokens[0]
    for j in range(1, len(tokens)):
        t = tokens[ind]
        ind = nxt[t][cnt[t]]
        cnt[t] += 1
        res[j] = tokens[ind]
    return chars[res].tobytes().decode()


def make_backgrounds(hg38, sizes, n, seed, tmp, chrom_re=r"^chr[0-9]+$", gc_lo=0.35, gc_hi=0.50):
    rng = np.random.default_rng(seed)
    sz = pd.read_csv(sizes, sep="\t", header=None, names=["chrom", "len"])
    sz = sz[sz.chrom.str.match(chrom_re)]
    out = []
    while len(out) < n:
        ch = rng.choice(sz.chrom.values, 400, p=sz.len.values / sz.len.sum())
        st = [int(rng.integers(0, sz.set_index("chrom").len[c] - SEQLEN)) for c in ch]
        with open(tmp + ".bed", "w") as fh:
            for i, (c, s) in enumerate(zip(ch, st)):
                fh.write(f"{c}\t{s}\t{s + SEQLEN}\tr{i}\n")
        subprocess.run(f"twoBitToFa -bed={tmp}.bed {hg38} {tmp}.fa", shell=True, check=True)
        seqs, cur = [], []
        for l in open(tmp + ".fa"):
            if l.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                cur = []
            else:
                cur.append(l.strip())
        seqs.append("".join(cur))
        for s in seqs:
            if len(out) >= n:
                break
            up = s.upper()
            if len(s) != SEQLEN or "N" in up or sum(c.islower() for c in s) / SEQLEN >= 0.5:
                continue
            g = (up.count("G") + up.count("C")) / SEQLEN
            if not gc_lo <= g <= gc_hi:
                continue
            for _ in range(20):
                sh = dinuc_shuffle(up, rng)
                if not AP1_RE.search(sh):
                    out.append(sh)
                    break
    return out


def insert(bg, motif, positions):
    s = list(bg)
    for p in positions:
        s[p:p + len(motif)] = motif
    return "".join(s)


def onehot(seqs):
    lut = np.zeros((256, 4), np.float32)
    for i, c in enumerate("ACGT"):
        lut[ord(c), i] = 1
    return lut[np.frombuffer("".join(seqs).encode(), np.uint8)].reshape(len(seqs), SEQLEN, 4)


def build_batch(bg):
    """Sequences per background: [background] + per insert [single] + [pair at d for every d]."""
    seqs = [bg]
    for m in INSERTS.values():
        seqs.append(insert(bg, m, [POS1]))
        seqs += [insert(bg, m, [POS1, POS1 + d]) for d in DISTS]
    return seqs


def load_predictor(tf, dl, repo, fn, cache):
    """Try the .h5 model (Keras 2); if that fails (e.g. Keras 3) use the SavedModel (.tar). Returns f -> log counts."""
    try:
        model = tf.keras.models.load_model(dl(repo, fn, cache_dir=cache), compile=False)
        return lambda x: np.asarray(model.predict(x, batch_size=256, verbose=0)[1]).reshape(-1)
    except Exception as e:  # noqa: BLE001
        print("h5 failed, using SavedModel:", type(e).__name__, str(e)[:200], flush=True)
        import tarfile
        tar = dl(repo, fn.replace(".h5", ".tar"), cache_dir=cache)
        d = tar + "_x"
        if not os.path.isdir(d):
            tarfile.open(tar).extractall(d)
        sub = [os.path.join(d, x) for x in os.listdir(d)][0]
        sig = tf.saved_model.load(sub).signatures["serving_default"]
        key = None

        def f(x):
            nonlocal key
            outs = []
            for i in range(0, len(x), 256):
                o = sig(tf.constant(x[i:i + 256]))
                if key is None:  # the (N, 1) output = log counts
                    key = [k for k, v in o.items() if v.shape[-1] == 1][0]
                outs.append(o[key].numpy().reshape(-1))
            return np.concatenate(outs)
        return f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hg38", required=True)
    ap.add_argument("--sizes", required=True)
    ap.add_argument("--models", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--nbg", type=int, default=100)
    ap.add_argument("--folds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    ap.add_argument("--only", nargs="+", default=list(MODELS))
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--predict-only", dest="predict_only", action="store_true",
                    help="only save per-fold predictions (to split work across GPUs); statistics in the final call")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    for nm, m in INSERTS.items():
        assert nm == "AP1" or not AP1_RE.search(m), m
    bgf = os.path.join(a.out, "backgrounds.txt")
    if not os.path.exists(bgf):
        open(bgf, "w").write("\n".join(make_backgrounds(a.hg38, a.sizes, a.nbg, a.seed, os.path.join(a.out, "tmp"))))
    bgs = open(bgf).read().split("\n")[:a.nbg]
    print("backgrounds:", len(bgs), flush=True)

    import tensorflow as tf
    from huggingface_hub import hf_hub_download
    from tensorflow.keras.utils import get_custom_objects
    get_custom_objects().update({"tf": tf})

    ns = 1 + len(DISTS)
    res = {}
    for name in a.only:
        repo = MODELS[name]
        enc = repo.split("-ENCSR")[1].split("-")[0]
        F = []
        for fold in a.folds:
            pf = os.path.join(a.out, f"pred_{name}_fold{fold}.npy")
            if not os.path.exists(pf):
                fn = f"fold_{fold}/model.chrombpnet_nobias.fold_{fold}.ENCSR{enc}.h5"
                predict = load_predictor(tf, hf_hub_download, repo, fn, a.models)
                pred = np.zeros((len(bgs), 1 + len(INSERTS) * ns), np.float32)
                for i, bg in enumerate(bgs):
                    pred[i] = predict(onehot(build_batch(bg)))
                np.save(pf, pred)
                print(name, "fold", fold, "done", flush=True)
                tf.keras.backend.clear_session()
            F.append(np.load(pf))
        if a.predict_only:
            continue
        P = np.mean(F, axis=0)  # backgrounds x sequences
        base = P[:, 0]
        r = {}
        for k, nm in enumerate(INSERTS):
            off = 1 + k * ns
            d1 = P[:, off] - base
            C = (P[:, off + 1: off + ns] - base[:, None]) - 2 * d1[:, None]  # backgrounds x d
            dist = np.array(DISTS)
            near, far = C[:, (dist >= 7) & (dist <= 60)].mean(1), C[:, dist > 60].mean(1)
            w = stats.wilcoxon(near - far, alternative="greater")
            mc = C.mean(0)
            r[nm] = {"delta1_mean": float(d1.mean()), "valid": bool(d1.mean() > 0),
                     "C_near_7_60": float(near.mean()), "C_far_61_100": float(far.mean()),
                     "wilcoxon_p_near_vs_far": float(w.pvalue), "median_diff_near_vs_far": float(np.median(near - far)),
                     "C_7_30": float(C[:, (dist >= 7) & (dist <= 30)].mean()),
                     "C_31_60": float(C[:, (dist > 30) & (dist <= 60)].mean()),
                     "d_max_C": int(dist[np.argmax(mc)]), "C_max": float(mc.max()), "curve_C": mc.tolist()}
        res[name] = r
        print(name, json.dumps({k: {kk: v for kk, v in x.items() if kk != "curve_C"} for k, x in r.items()}), flush=True)
        json.dump(res, open(os.path.join(a.out, "cooperativity.json"), "w"), indent=1)
    if a.predict_only:
        return
    valid = [n for n in res if res[n]["AP1"]["valid"]]
    res["decision"] = {
        "valid_models": valid,
        "n_models_supported": sum(res[n]["AP1"]["wilcoxon_p_near_vs_far"] < 0.05 for n in valid),
        "negative_controls_passing": {n: [k for k in ("neg1", "neg2", "neg3") if res[n][k]["wilcoxon_p_near_vs_far"] < 0.05]
                                  for n in res if n != "decision"}}
    json.dump(res, open(os.path.join(a.out, "cooperativity.json"), "w"), indent=1)
    print(json.dumps(res["decision"], indent=1))


if __name__ == "__main__":
    main()
