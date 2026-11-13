"""Genome coordinates of the placed chromatin (IMPROVEMENTS K1).

Each nucleosome walk of the Rabl layout (wc3d.rabl_segments) is one arm segment, walked from the centromere to the
telomere (left arms towards position 0), except XII-R, which is split at the rDNA (XIIR_in: centromere to the rDNA;
XIIR_out: from the rDNA to the telomere). The walk order of nucleosomes is genomic order. Nucleosomes that did not fit
(about 2%) and the three evicted around each Pol II leave gaps, so positions are interpolated along each placed walk:
bp = segment start + (j + 0.5) / n_placed x segment length (accurate to about 2% of the arm). With the measured
nucleosome map (wc3d.NUC_MAP, nucleosome_map.py) each nucleosome carries its dyad position (prepass nuc_bp) instead.
Genes: SGD verified and uncharacterised ORFs (pipeline/data/sgd/SGD_features.tab; coordinates of the reference genome
R64, in which the rDNA array is collapsed to two repeats).
Histone marks (IMPROVEMENTS K7, annotation only): log2 ChIP/input at t = 0 (unstressed, YPD) of the nearest MNase
nucleosome within 40 bp in Weiner et al. 2015 Mol Cell 58:371 (GSE61888), for MARKS; NaN where none is that close.
Output: work/chromatin/nucleosome_genome.parquet (index into prepass nuc_pos, chromosome 1-16, bp, arm, ORF, gene,
marks)
"""
import os, json
import numpy as np
import pandas as pd

from config import W3D, WORK, DATA
from wc3d import CHROM_KB, CEN_KB, RDNA_KB

OUT = f"{WORK}/chromatin"
WEINER = f"{DATA}/validation/genome/histone_marks_weiner2015/GSE61888_nucs_normed.csv.gz"
MARKS = ["Htz1", "H3K4me3", "H3K36me3", "H3K79me3", "H3K56ac", "H3K9ac", "H4K16ac"]


def segments():
    """(chromosome index, kind, start bp, end bp) in the order of rabl_segments; start = where the walk begins"""
    segs = []
    for c, (kb, cen) in enumerate(zip(CHROM_KB, CEN_KB)):
        if c == 11:
            segs += [(c, "L", cen * 1000, 0), (c, "XIIR_in", cen * 1000, RDNA_KB[0] * 1000),
                     (c, "XIIR_out", RDNA_KB[1] * 1000, kb * 1000)]
        else:
            segs += [(c, "L", cen * 1000, 0), (c, "R", cen * 1000, kb * 1000)]
    return segs


def genes():
    f = pd.read_csv(f"{DATA}/sgd/SGD_features.tab", sep="\t", header=None, dtype=str)
    f = f[(f[1] == "ORF") & f[2].isin(["Verified", "Uncharacterized"])]
    roman = {"I": 1, "II": 2, "III": 3, "IV": 4, "V": 5, "VI": 6, "VII": 7, "VIII": 8, "IX": 9, "X": 10, "XI": 11,
             "XII": 12, "XIII": 13, "XIV": 14, "XV": 15, "XVI": 16}
    chrom = pd.to_numeric(f[8], errors="coerce")
    a, b = pd.to_numeric(f[9], errors="coerce"), pd.to_numeric(f[10], errors="coerce")
    g = pd.DataFrame(dict(chrom=chrom, start=np.minimum(a, b), end=np.maximum(a, b), orf=f[3], gene=f[4].fillna(f[3])))
    return g.dropna(subset=["chrom", "start", "end"]).astype({"chrom": int, "start": int, "end": int})


def build(log=print):
    os.makedirs(OUT, exist_ok=True)
    z = np.load(f"{W3D}/prepass.npz")
    ch = z["nuc_chain"]
    segs = segments()
    assert ch.max() + 1 == len(segs), f"{ch.max() + 1} walks vs {len(segs)} arm segments"
    rows = []
    for k, (c, kind, s0, s1) in enumerate(segs):
        idx = np.where(ch == k)[0]
        if "nuc_bp" in z.files:             # measured nucleosome map: exact dyad positions
            bp = z["nuc_bp"][idx]
        else:
            frac = (np.arange(len(idx)) + 0.5) / max(1, len(idx))
            bp = (s0 + frac * (s1 - s0)).astype(np.int64)
        rows.append(pd.DataFrame(dict(nuc=idx, chrom=c + 1, bp=bp, arm=kind)))
    d = pd.concat(rows, ignore_index=True)
    g = genes()
    d["orf"] = ""; d["gene"] = ""
    for c, gg in g.groupby("chrom"):
        m = d.chrom == c
        st, en = gg.start.values, gg.end.values
        o = np.argsort(st); st, en, orf, gn = st[o], en[o], gg.orf.values[o], gg.gene.values[o]
        pos = d.loc[m, "bp"].values
        j = np.searchsorted(st, pos, side="right") - 1
        hit = (j >= 0) & (pos <= en[np.clip(j, 0, None)])
        d.loc[m, "orf"] = np.where(hit, orf[np.clip(j, 0, None)], "")
        d.loc[m, "gene"] = np.where(hit, gn[np.clip(j, 0, None)], "")
    if os.path.exists(WEINER):
        d = d.join(marks(d))
    d.to_parquet(f"{OUT}/nucleosome_genome.parquet", index=False)
    log(f"genome map: {len(d)} nucleosomes on {d.chrom.nunique()} chromosomes; {int((d.orf != '').sum())} in ORFs "
        f"({(d.orf != '').mean():.0%}); {d[d.orf != ''].orf.nunique()} ORFs covered")
    return d


def marks(d, tol=40):
    """per nucleosome of d: the MARKS at t = 0 of the nearest Weiner 2015 nucleosome within tol bp"""
    h = pd.read_csv(WEINER, header=None, nrows=2)
    ab, tm = h.iloc[0].values, h.iloc[1].values
    cols = [1, 2] + [i for i, (a, t) in enumerate(zip(ab, tm)) if a in MARKS and str(t) == "0"]
    w = pd.read_csv(WEINER, header=None, skiprows=2, usecols=cols)
    w.columns = ["chrom", "bp"] + [ab[i] for i in cols[2:]]
    out = pd.DataFrame(np.nan, index=d.index, columns=MARKS)
    for c, g in d.groupby("chrom"):
        wc = w[w.chrom == c].sort_values("bp")
        if not len(wc):
            continue
        x = wc.bp.values; j = np.clip(np.searchsorted(x, g.bp.values), 1, len(x) - 1)
        jj = np.where(np.abs(x[j] - g.bp.values) < np.abs(x[j - 1] - g.bp.values), j, j - 1)
        ok = np.abs(x[jj] - g.bp.values) <= tol
        vals = wc[MARKS].values[jj]; vals[~ok] = np.nan
        out.loc[g.index] = vals
    return out


def write_loci(data_dir, d=None, log=print):
    """library/gene_loci.json for the viewer's search: per ORF its chromosome, coordinates and locus (centroid of the
    nucleosomes that overlap it, or the nucleosome nearest its midpoint for genes shorter than a nucleosome repeat)"""
    import json
    if d is None:
        d = pd.read_parquet(f"{OUT}/nucleosome_genome.parquet")
    P = np.load(f"{W3D}/prepass.npz")["nuc_pos"].astype(np.float64)
    g = genes()
    out = {}
    for c, gg in g.groupby("chrom"):
        dc = d[d.chrom == c].sort_values("bp")
        if len(dc) == 0:
            continue
        bp, idx = dc.bp.values, dc.nuc.values
        for r in gg.itertuples(index=False):
            m = (bp >= r.start) & (bp <= r.end)
            k = idx[m] if m.any() else idx[[int(np.argmin(np.abs(bp - (r.start + r.end) / 2)))]]
            x = P[k].mean(0)
            out[r.gene] = [round(float(x[0]), 1), round(float(x[1]), 1), round(float(x[2]), 1), int(c), int(r.start), int(r.end), r.orf]
    json.dump(dict(note="gene locus: centroid of the gene's nucleosomes (nm, cell frame); chrom, start, end (R64), ORF",
                   loci=out), open(f"{data_dir}/library/gene_loci.json", "w"), separators=(",", ":"))
    log(f"gene_loci.json: {len(out)} genes")
    return out


if __name__ == "__main__":
    from config import WEB_OUT
    write_loci(WEB_OUT, build())
