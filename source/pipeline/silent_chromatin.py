"""Silent chromatin from Sir3 ChIP-seq (IMPROVEMENTS Q2): which nucleosomes carry the Sir complex.

Before, Sir2/3/4 went on every nucleosome within 3 kb of a chromosome end and over HML and HMR (trna_telomeres.py), a
fixed extent. Here a nucleosome is Sir-bound where Sir3 ChIP-seq is enriched over input at its position:
anti-Sir3 ChIP of cycling W303 cells in YPD at 30 C, two replicates with inputs (Brothers & Rine 2022 eLife
11:e75653, GEO GSE189038; per-base coverage over the genome median, sacCer3 with the MAT locus removed). The strains
express Sir3 fused to the M.EcoGII methyltransferase, not untagged Sir3 (a caveat; the paper shows it silences).
Enrichment of a nucleosome: ChIP / input summed over +-HALF bp around its dyad; it is Sir-bound where the enrichment
is >= THRESH in both replicates. THRESH = 5: untagged and sir3-delta controls of the same laboratory's Sir4 ChIP
reach about 1 and 4-5 at silent loci (Goodnight & Rine 2020), and highly expressed genes are "hyper-ChIPable"
(Teytelman et al. 2013), so a low threshold would add active genes. For the same reason only nucleosomes within
SUBTEL_KB of a chromosome end or in HML/HMR can be Sir-bound: 33 nucleosomes elsewhere pass the threshold, at TDH3
and other highly expressed genes. Both were set before the test below was run.
Held out (validation/silent_chromatin_test.py): H4K16 acetylation of the same nucleosomes (Weiner et al. 2015; Sir2
deacetylates H4K16) and Sir4 ChIP-seq (Goodnight & Rine 2020).
"""
import os, gzip
import numpy as np
import pandas as pd

from config import DATA, WORK

DIR = f"{DATA}/validation/genome/sir3_chipseq_brothers2022"
REPS = [("GSM5823538_220115-JRY13670-IP", "GSM5823537_220115-JRY13670-IN"),
        ("GSM5823540_220115-JRY13671-IP", "GSM5823539_220115-JRY13671-IN")]
HALF, THRESH, SUBTEL_KB = 80, 5.0, 30.0
ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"]
CACHE = f"{WORK}/chromatin/sir3_tracks.npz"


def read_bedgraph(path):
    """per-base bedGraph of this laboratory (1-based start, value = coverage / genome median) -> {chrom 1-16: array}"""
    d = pd.read_csv(path, sep="\t", comment="#", skiprows=1, header=None, names=["c", "s", "e", "v"],
                    dtype={"c": str, "s": np.int64, "e": np.int64, "v": np.float32}, skip_blank_lines=True)
    out = {}
    for i, r in enumerate(ROMAN):
        x = d[d.c == r]
        a = np.zeros(int(x.s.max()) + 1, np.float32); a[x.s.values] = x.v.values
        out[i + 1] = a
    return out


def tracks():
    """{(replicate, 'ip'|'in', chrom): per-base array}, cached"""
    if os.path.exists(CACHE):
        z = np.load(CACHE)
        return {tuple(k.split("|")[:2]) + (int(k.split("|")[2]),): z[k] for k in z.files}
    out = {}
    for ri, (ip, inp) in enumerate(REPS):
        for kind, name in (("ip", ip), ("in", inp)):
            for c, a in read_bedgraph(f"{DIR}/{name}_median_normalized.bedGraph.gz").items():
                out[(str(ri), kind, c)] = a
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    np.savez_compressed(CACHE, **{f"{a}|{b}|{c}": v for (a, b, c), v in out.items()})
    return out


def enrichment(chrom, bp, T=None):
    """per position (chromosome 1-16, bp): Sir3 ChIP / input over +-HALF bp, one column per replicate"""
    T = tracks() if T is None else T
    chrom, bp = np.asarray(chrom), np.asarray(bp)
    out = np.zeros((len(bp), len(REPS)))
    for ri in range(len(REPS)):
        for c in np.unique(chrom):
            ip, inp = T[(str(ri), "ip", int(c))], T[(str(ri), "in", int(c))]
            n = min(len(ip), len(inp))
            ci, cn = np.r_[0, np.cumsum(ip[:n], dtype=np.float64)], np.r_[0, np.cumsum(inp[:n], dtype=np.float64)]
            k = np.where(chrom == c)[0]
            lo = np.clip(bp[k] - HALF, 0, n); hi = np.clip(bp[k] + HALF + 1, 0, n)
            out[k, ri] = (ci[hi] - ci[lo]) / np.maximum(cn[hi] - cn[lo], 0.2 * (hi - lo + 1))
    return out


def silent_region(chrom, bp):
    """boolean: within SUBTEL_KB of a chromosome end, or inside HML / HMR (SGD)"""
    import trna_telomeres as tt
    from wc3d import CHROM_KB
    chrom, bp = np.asarray(chrom), np.asarray(bp)
    L = np.array([CHROM_KB[c - 1] * 1000 for c in chrom])
    ok = np.minimum(bp, L - bp) < SUBTEL_KB * 1000
    for c, a, b in tt.silent_loci():
        ok |= (chrom == c) & (bp >= a) & (bp <= b)
    return ok


def sir_bound(chrom, bp, T=None):
    """boolean per nucleosome: Sir3 enriched >= THRESH in every replicate, in a subtelomere or silent cassette"""
    return (enrichment(chrom, bp, T) >= THRESH).all(1) & silent_region(chrom, bp)
