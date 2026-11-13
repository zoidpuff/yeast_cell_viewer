"""Silent chromatin: the Sir-bound nucleosomes of the model against held-out data (IMPROVEMENTS Q2).

The model's Sir-bound nucleosomes come from Sir3 ChIP-seq (pipeline/silent_chromatin.py). Not used for that:
  - H4K16 acetylation per nucleosome in unstressed cells (Weiner et al. 2015 Mol Cell 58:371; log2 ChIP over input;
    Sir2 deacetylates H4K16, so Sir-bound nucleosomes should be hypoacetylated);
  - Sir4-Myc ChIP-seq with inputs and an untagged control (Goodnight & Rine 2020 eLife 9:e58910; W303, 37 C).
Compared: the Sir-bound set, the earlier rule (within 3 kb of an end, or HML/HMR), nucleosomes only one of the two
selects, and all nucleosomes. Output: work/validation/silent_chromatin.json
"""
import os, sys, json
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
from config import WORK, W3D, DATA  # noqa: E402

G = f"{DATA}/validation/genome/sir4_chipseq_goodnight2020"
SIR4 = [("GSM4675918_JRY12171", "GSM4675922_JRY12171-in"), ("GSM4675919_JRY12172", "GSM4675923_JRY12172-in")]
NOTAG = ("GSM4675940_JRY12269", "GSM4675941_JRY12269-in")


def ratio_at(ip, inp, chrom, bp, half=80):
    out = np.zeros(len(bp))
    for c in np.unique(chrom):
        a, b = ip[int(c)], inp[int(c)]; n = min(len(a), len(b))
        ca, cb = np.r_[0, np.cumsum(a[:n], dtype=np.float64)], np.r_[0, np.cumsum(b[:n], dtype=np.float64)]
        k = np.where(chrom == c)[0]
        lo = np.clip(bp[k] - half, 0, n); hi = np.clip(bp[k] + half + 1, 0, n)
        out[k] = (ca[hi] - ca[lo]) / np.maximum(cb[hi] - cb[lo], 0.2 * (hi - lo + 1))
    return out


def main():
    import silent_chromatin as sc, genome_map as gm, trna_telomeres as tt
    from wc3d import CHROM_KB
    z = np.load(f"{W3D}/prepass.npz"); segs = gm.segments()
    chrom = np.array([segs[s][0] + 1 for s in z["nuc_chain"]]); bp = z["nuc_bp"]
    new = sc.sir_bound(chrom, bp)
    L = np.array([CHROM_KB[c - 1] * 1000 for c in chrom]); dend = np.minimum(bp, L - bp)
    old = dend < 3000
    for c, a, b in tt.silent_loci():
        old |= (chrom == c) & (bp >= a) & (bp <= b)
    # H4K16ac of the nearest Weiner nucleosome within 40 bp
    h = pd.read_csv(gm.WEINER, header=None, nrows=2); ab, tm = h.iloc[0].values, h.iloc[1].values
    i16 = [i for i, (a, t) in enumerate(zip(ab, tm)) if a == "H4K16ac" and str(t) == "0"][0]
    w = pd.read_csv(gm.WEINER, header=None, skiprows=2, usecols=[1, 2, i16]); w.columns = ["chrom", "bp", "k16"]
    k16 = np.full(len(bp), np.nan)
    for c in np.unique(chrom):
        wc = w[w.chrom == c].sort_values("bp"); x = wc.bp.values; m = np.where(chrom == c)[0]
        j = np.clip(np.searchsorted(x, bp[m]), 1, len(x) - 1)
        jj = np.where(np.abs(x[j] - bp[m]) < np.abs(x[j - 1] - bp[m]), j, j - 1)
        ok = np.abs(x[jj] - bp[m]) <= 40
        k16[m[ok]] = wc.k16.values[jj[ok]]
    low = np.nanpercentile(k16, 10)
    rd = lambda n: sc.read_bedgraph(f"{G}/{n}_0_500_median_normalized.bedgraph.gz")
    sir4 = np.mean([ratio_at(rd(a), rd(b), chrom, bp) for a, b in SIR4], axis=0)
    notag = ratio_at(rd(NOTAG[0]), rd(NOTAG[1]), chrom, bp)
    sets = {"sir_bound (model)": new, "earlier rule (3 kb, HML/HMR)": old, "model only": new & ~old,
            "earlier rule only": old & ~new, "all nucleosomes": np.ones(len(bp), bool),
            "subtelomeric 30 kb, not Sir-bound": sc.silent_region(chrom, bp) & ~new}
    res = dict(threshold=sc.THRESH, subtel_kb=sc.SUBTEL_KB, h4k16ac_genome_10th_percentile=float(low), sets={})
    for k, m in sets.items():
        res["sets"][k] = dict(n=int(m.sum()), h4k16ac_median_log2=float(np.nanmedian(k16[m])),
                              h4k16ac_hypoacetylated_fraction=float(np.nanmean(k16[m] < low)),
                              sir4_chip_over_input_median=float(np.median(sir4[m])),
                              sir4_enriched_3x_fraction=float(np.mean(sir4[m] >= 3)), notag_median=float(np.median(notag[m])))
    tel = set((int(chrom[i]), bool(bp[i] < L[i] / 2)) for i in np.where(new & (dend < sc.SUBTEL_KB * 1000))[0])
    res["chromosome_ends_with_sir_bound_nucleosomes"] = len(tel)
    res["sir_bound_by_distance_kb"] = {"0-3": int((new & (dend < 3000)).sum()), "3-10": int((new & (dend >= 3000) & (dend < 10000)).sum()),
                                       "10-30": int((new & (dend >= 10000) & (dend < 30000)).sum())}
    json.dump(res, open(f"{WORK}/validation/silent_chromatin.json", "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
