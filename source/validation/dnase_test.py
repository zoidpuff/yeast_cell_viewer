"""DNase I accessibility against the model's nucleosomes at base resolution (IMPROVEMENTS Q3; held-out test).

Data (not used by the model):
  - per-base DNase I cut counts, log-phase W303 in YPD, two replicates (Zhong et al. 2016 Genome Res 26:351, GEO
    GSE69651; sacCer2, lifted to sacCer3 here with the UCSC chain; GEO lists 37 C, possibly the digestion temperature);
  - DNase-hypersensitive sites of wild-type BY4741 in YPD (Zelin et al. 2012 Mol Cell 48:459, GSM989336; 3,260
    F-seq peaks, sacCer3).
Model: the placed nucleosomes (prepass nuc_bp; a core covers +-73 bp of its dyad). The nucleosome positions come from
a chemical-cleavage map (an input), so this tests that map and what the model kept of it (evictions at Pol II and
origins, promoter NDRs) against an independent assay, more finely than the ATAC regions of nucleosome_test.py.
Measures:
  cut profile     mean cut density (per bp, relative to the genome mean) by distance from a model dyad; a protected
                  core shows as a dip over +-60 bp against the linkers at +-(75-95) bp
  protected       share of nucleosomes with a lower cut density in their core (+-50 bp) than in their linkers
  bins            Spearman between model occupancy (share of bp inside a core) and cut density in 50 bp bins
  hypersensitive  model occupancy inside the DHS peaks against the genome average
Output: work/validation/dnase_test.json
"""
import os, sys, json, gzip
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
from config import WORK, W3D, DATA  # noqa: E402

G = f"{DATA}/validation/genome"
ZHONG = [f"{G}/dnase_zhong2016/GSM1705337_DNase-seq_W303_S_cerevisiae_Sample1.csv.gz",
         f"{G}/dnase_zhong2016/GSM1705338_DNase-seq_W303_S_cerevisiae_Sample2.csv.gz"]
ZELIN = f"{G}/dnase_zelin2012/GSM989336_BY4741_wild_type_DHS_peak.bed.gz"
CHAIN = f"{G}/liftover_ucsc/sacCer2ToSacCer3.over.chain.gz"
ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"]
CORE, W = 73, 150


def chain_map():
    """{chrom 1-16: array old 0-based position -> new 0-based position, -1 where unaligned} from a UCSC chain file"""
    out = {}
    with gzip.open(CHAIN, "rt") as f:
        cur = None
        for line in f:
            p = line.split()
            if not p:
                continue
            if p[0] == "chain":
                c = p[2].replace("chr", "")
                cur = ROMAN.index(c) + 1 if c in ROMAN and p[7] == p[2] and p[4] == "+" and p[9] == "+" else None
                if cur is not None:
                    m = out.setdefault(cur, np.full(int(p[3]), -1, np.int64)); t, q = int(p[5]), int(p[10])
            elif cur is not None:
                n = int(p[0]); m[t:t + n] = np.arange(q, q + n); t += n; q += n
                if len(p) == 3:
                    t += int(p[1]); q += int(p[2])
    return out


def cuts(lengths):
    """{chrom: per-base cut counts in sacCer3, replicates summed}"""
    cm = chain_map()
    out = {c: np.zeros(L, np.float64) for c, L in lengths.items()}
    for f in ZHONG:
        d = pd.read_csv(f, dtype={"chr": str})
        for c in range(1, 17):
            x = d[d.chr == str(c)]
            new = cm[c][np.clip(x.pos.values - 1, 0, len(cm[c]) - 1)]
            ok = (new >= 0) & (new < lengths[c])
            np.add.at(out[c], new[ok], x.total_count.values[ok])
    return out


def main():
    import genome_map as gm
    from wc3d import CHROM_KB
    from scipy.stats import spearmanr
    z = np.load(f"{W3D}/prepass.npz"); segs = gm.segments()
    chrom = np.array([segs[s][0] + 1 for s in z["nuc_chain"]]); bp = z["nuc_bp"]
    lengths = {c: CHROM_KB[c - 1] * 1000 + 2000 for c in range(1, 17)}
    C = cuts(lengths)
    tot = sum(v.sum() for v in C.values()); mean = tot / sum(len(v) for v in C.values())
    prof = np.zeros(2 * W + 1); n = 0; prot = []
    occ_bins, cut_bins = [], []
    occ = {}
    for c in range(1, 17):
        d = np.sort(bp[chrom == c]); d = d[(d > W) & (d < lengths[c] - W - 1)]
        a = C[c]
        idx = d[:, None] + np.arange(-W, W + 1)[None]
        m = a[idx]; prof += m.sum(0); n += len(d)
        core = m[:, W - 50:W + 51].mean(1); link = np.r_["1", m[:, W - 95:W - 74], m[:, W + 75:W + 96]].mean(1)
        prot.append(core < link)
        o = np.zeros(lengths[c], np.float64)
        for s in d:
            o[s - CORE:s + CORE + 1] = 1
        occ[c] = o
        k = (lengths[c] // 50) * 50
        occ_bins.append(o[:k].reshape(-1, 50).mean(1)); cut_bins.append(a[:k].reshape(-1, 50).mean(1))
    prof = prof / n / mean
    prot = np.concatenate(prot); ob, cb = np.concatenate(occ_bins), np.concatenate(cut_bins)
    pk = pd.read_csv(ZELIN, sep="\t", header=None, usecols=[0, 1, 2], names=["c", "s", "e"])
    inside = tot_bp = 0.0
    for r in pk.itertuples():
        cc = r.c.replace("chr", "")
        if cc in ROMAN:
            o = occ[ROMAN.index(cc) + 1][r.s:r.e]; inside += o.sum(); tot_bp += len(o)
    genome_occ = sum(o.sum() for o in occ.values()) / sum(len(o) for o in occ.values())
    res = dict(cuts_total=float(tot), nucleosomes=int(n),
               cut_profile_relative=dict(dyad=float(prof[W - 10:W + 11].mean()), core_50=float(prof[W - 50:W + 51].mean()),
                                         linker_75_95=float(np.r_[prof[W - 95:W - 74], prof[W + 75:W + 96]].mean())),
               core_over_linker=float(prof[W - 50:W + 51].mean() / np.r_[prof[W - 95:W - 74], prof[W + 75:W + 96]].mean()),
               nucleosomes_protected_fraction=float(prot.mean()),
               bins_50bp_spearman_occupancy_vs_cuts=float(spearmanr(ob, cb).correlation),
               cut_density_relative_by_occupancy=dict(free=float(cb[ob < 0.1].mean() / mean), mixed=float(cb[(ob >= 0.1) & (ob <= 0.9)].mean() / mean),
                                                      covered=float(cb[ob > 0.9].mean() / mean)),
               dhs_peaks=int(len(pk)), model_occupancy_in_dhs=float(inside / tot_bp), model_occupancy_genome=float(genome_occ),
               profile=[round(float(x), 3) for x in prof[::5]])
    json.dump(res, open(f"{WORK}/validation/dnase_test.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "profile"}, indent=1))


if __name__ == "__main__":
    main()
