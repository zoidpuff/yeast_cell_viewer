"""Measured nucleosome positions for the chromatin fibre (IMPROVEMENTS K2).

Dyads are called from the H3 chemical-cleavage dyad maps (H3Q85C; W303 background, YPD) of Chereji et al. 2018
Genome Biol 19:19 (GSE97290; replicates 1-3 summed; sacCer3 = R64, the coordinates of the model's chromosomes): the
summed dyad counts are smoothed (Gaussian, SIGMA_BP), and peaks are taken greedily from the highest down, each
excluding EXCL_BP on both sides, while the smoothed signal is above MIN_FRAC x its genome-wide median (below that the
region counts as nucleosome-free). Dyads inside promoter nucleosome-depleted regions are then removed (NDR_FILTER):
the regions between the upstream stable nucleosome and the +1 nucleosome of Rossi et al. 2021 Nature 592:309
(Supplementary Data 1; dyads more than 73 bp inside both boundary dyads; imputed or > 2 kb regions skipped).
Brogaard et al. 2012 (chemical map) is on the older sacCer2 coordinates (offsets of a few bp that grow along the
arms) and Jiang & Pugh 2009 on an intermediate pre-sacCer3 release, so neither is used as the input.
Output: data/genome/nucleosome_dyads.tsv.gz (chrom 1-16, dyad bp, smoothed signal), cached.
"""
import os
import numpy as np
import pandas as pd

from config import DATA

SRC = f"{DATA}/validation/genome/nucleosomes_chereji2018"
REPS = ["GSM2561057_Dyads_H3_CC_rep_1.bw", "GSM2561058_Dyads_H3_CC_rep_2.bw", "GSM2561059_Dyads_H3_CC_rep_3.bw"]
OUT = f"{DATA}/genome/nucleosome_dyads.tsv.gz"
SIGMA_BP, EXCL_BP, MIN_FRAC = 15.0, 135, 0.5
NDR_TABLE = f"{DATA}/validation/genome/chipexo_rossi2021/41586_2021_3314_MOESM3_ESM.xlsx"
NDR_FILTER = True
ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"]


def _call(s, excl, thr):
    o = np.argsort(-s); taken = np.zeros(len(s), bool); out = []
    for i in o:
        if s[i] < thr:
            break
        if taken[i]:
            continue
        out.append(i); taken[max(0, i - excl + 1):i + excl] = True
    return np.sort(np.array(out, np.int64))


def call(log=print):
    import pyBigWig
    from scipy.ndimage import gaussian_filter1d
    bw = [pyBigWig.open(f"{SRC}/{f}") for f in REPS]
    sig = {}
    for k, r in enumerate(ROMAN):
        c = "chr" + r; L = bw[0].chroms()[c]
        sig[k + 1] = gaussian_filter1d(sum(np.nan_to_num(np.array(b.values(c, 0, L))) for b in bw), SIGMA_BP)
    thr = MIN_FRAC * np.median(np.concatenate(list(sig.values())))
    rows = []
    for c, s in sig.items():
        p = _call(s, EXCL_BP, thr)
        rows.append(pd.DataFrame(dict(chrom=c, bp=p + 1, signal=s[p].round(3))))      # bigwig is 0-based
    d = pd.concat(rows, ignore_index=True)
    n0 = len(d)
    if NDR_FILTER:
        d = d[~in_ndr(d)].reset_index(drop=True)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    d.to_csv(OUT, sep="\t", index=False)
    sp = np.diff(d.bp.values)[np.diff(d.chrom.values) == 0]
    log(f"nucleosome map: {len(d)} dyads (Chereji 2018 H3 chemical map; {n0 - len(d)} removed in promoter NDRs); spacing median {np.median(sp):.0f} bp, "
        f"{(sp > 300).mean():.1%} gaps > 300 bp")
    return d


def in_ndr(d):
    x = pd.read_excel(NDR_TABLE, sheet_name="Supplementary Data 1")
    L = x[[c for c in x.columns if str(c).startswith("NFR/NDR_Dyad-to-Dyad_Length")][0]]
    x = x[(L != 9999) & (L <= 2000) & x["NFR/NDR_Left"].notna()]
    rm = np.zeros(len(d), bool)
    for c, g in x.groupby("Chrom"):
        m = np.where(d.chrom.values == int(c[3:]))[0]; p = d.bp.values[m]
        for a, b in zip(g["NFR/NDR_Left"].values + 73, g["NFR/NDR_Right"].values - 73):
            if b > a:
                rm[m[(p > a) & (p < b)]] = True
    return rm


def dyads():
    return pd.read_csv(OUT, sep="\t") if os.path.exists(OUT) else call()


def arm_count():
    """dyads on the chromosome arms (the reference rDNA repeats on XII, built separately, excluded)"""
    import genome_map
    return int(sum(len(p) for p in segment_dyads(genome_map.segments())))


def segment_dyads(segs):
    """per arm segment (genome_map.segments order: chromosome index, kind, start bp, end bp) the dyad positions in walk
    order (from the start of the walk towards its end)"""
    d = dyads()
    out = []
    for c, kind, s0, s1 in segs:
        p = d.bp.values[d.chrom.values == c + 1]
        lo, hi = min(s0, s1), max(s0, s1)
        if kind == "L":
            lo = -1                     # to the telomere (chromosome lengths are rounded to kb)
        if kind in ("R", "XIIR_out"):
            hi = 10 ** 9
        p = p[(p >= lo) & (p < hi)]
        out.append(p[::-1] if s1 < s0 else p)
    return out


if __name__ == "__main__":
    call()
