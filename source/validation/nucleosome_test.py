"""Nucleosome positions of the model against independent maps (IMPROVEMENTS K2; all held out).

Model: the placed nucleosomes with their dyads (prepass nuc_bp); without a placed build, the input map
(nucleosome_map.py: Chereji 2018 chemical map, Rossi 2021 promoter NDRs removed).
  count            nucleosomes on the arms against 57,000-60,000 per cell from absolute occupancy (Oberbeckmann et
                   al. 2019 Genome Res 29:1996)
  occupancy        Oberbeckmann ORE-seq absolute occupancy at 13,807 restriction sites (BY4741, YPDA): mean occupancy
                   at sites inside a model nucleosome (dyad within 73 bp) and outside, and the AUROC of occupancy for
                   'covered'
  weiner_match     share of MNase nucleosomes of Weiner et al. 2015 (t = 0, unstressed; 66,359) with a model dyad within
                   40 bp, and the reverse
  atac             model dyads per kb inside the ATAC-accessible regions of Schep et al. 2015 (BY4741 YPD; merged
                   regions, median 1.3 kb), against the same regions shifted 2 kb and the genome average
Output: work/validation/nucleosome_test.json
"""
import os, sys, json, tarfile, io
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
from config import WORK, W3D, DATA

G = f"{DATA}/validation/genome"
ROMAN = {r: i + 1 for i, r in enumerate(["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"])}


def model_dyads():
    z = np.load(f"{W3D}/prepass.npz")
    if "nuc_bp" in z.files:
        import genome_map
        segs = genome_map.segments()
        chrom = np.array([segs[c][0] + 1 for c in z["nuc_chain"]])
        return pd.DataFrame(dict(chrom=chrom, bp=z["nuc_bp"])), "placed"
    import nucleosome_map
    return nucleosome_map.dyads()[["chrom", "bp"]], "input map (no placed build with nuc_bp)"


def nearest(a, b):
    """per position in a: distance to the nearest position in b (same chromosome); a, b: {chrom: sorted array}"""
    out = {}
    for c, x in a.items():
        y = b.get(c, np.array([], np.int64))
        if not len(y):
            out[c] = np.full(len(x), 10 ** 9); continue
        j = np.clip(np.searchsorted(y, x), 1, len(y) - 1) if len(y) > 1 else np.zeros(len(x), int)
        out[c] = np.minimum(np.abs(x - y[j]), np.abs(x - y[np.maximum(j - 1, 0)]))
    return out


def by_chrom(d):
    return {int(c): np.sort(g.bp.values.astype(np.int64)) for c, g in d.groupby("chrom")}


def main():
    m, src = model_dyads()
    M = by_chrom(m)
    res = dict(model=src, count=int(len(m)), count_oberbeckmann_range=[57000, 60000])
    o = pd.read_csv(f"{G}/nucleosome_occupancy_oberbeckmann2019/GSE141056_ORE-seq_map.tsv.gz", sep="\t")
    o = o[o.chr.str.match(r"^chr\d+$")]
    o["chrom"] = o.chr.str[3:].astype(int)
    dist = nearest(by_chrom(o.rename(columns={"pos": "bp"})), M)
    cov = np.concatenate([dist[c] <= 73 for c in sorted(dist)])
    occ = np.concatenate([g.sort_values("pos").occ.values for c, g in sorted(o.groupby("chrom"), key=lambda t: t[0])])
    au = mannwhitneyu(occ[cov], occ[~cov]).statistic / (cov.sum() * (~cov).sum())
    res["occupancy"] = dict(sites=int(len(occ)), covered_fraction=float(cov.mean()), mean_occ_covered=float(occ[cov].mean()),
                            mean_occ_uncovered=float(occ[~cov].mean()), auroc=float(au),
                            expected_coverage_from_occ=float(occ.mean()))
    w = pd.read_csv(f"{G}/histone_marks_weiner2015/GSE61888_nucs_normed.csv.gz", header=1, usecols=[1, 2])
    w.columns = ["chrom", "bp"]
    W = by_chrom(w.dropna().astype(int))
    dw = nearest(W, M); dm = nearest(M, W)
    res["weiner_match"] = dict(weiner_with_model_within_40bp=float(np.mean(np.concatenate(list(dw.values())) <= 40)),
                               model_with_weiner_within_40bp=float(np.mean(np.concatenate(list(dm.values())) <= 40)))
    with tarfile.open(f"{G}/atac_schep2015/GSE66386_all_nucleosome_calling.tar.gz") as t:
        f = [n for n in t.getnames() if n.endswith("lin_all.peaks.final.bed")][0]
        pk = pd.read_csv(io.BytesIO(t.extractfile(f).read()), sep="\t", header=None, usecols=[0, 1, 2], names=["chr", "a", "b"])
    pk["chrom"] = pk.chr.str.replace("chr", "").map(ROMAN)
    pk = pk.dropna()

    def density(shift):
        n, L = 0, 0
        for c, g in pk.groupby("chrom"):
            y = M.get(int(c), np.array([]))
            lo, hi = g.a.values + shift, g.b.values + shift
            n += int((np.searchsorted(y, hi) - np.searchsorted(y, lo)).sum()); L += int((hi - lo).sum())
        return n / L * 1000
    gen = len(m) / 12.07e3
    res["atac"] = dict(regions=int(len(pk)), median_width_bp=float((pk.b - pk.a).median()),
                       dyads_per_kb_inside=density(0), dyads_per_kb_shifted_2kb=density(2000), dyads_per_kb_genome=gen)
    json.dump(res, open(f"{WORK}/validation/nucleosome_test.json", "w"), indent=1)
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    main()
