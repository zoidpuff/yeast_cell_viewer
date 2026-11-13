"""Chromatin folding of the packed cell against Hi-C (IMPROVEMENTS K4; held out: no Hi-C enters the model).

Model "Hi-C": nucleosome pairs closer than CONTACT_NM (centre distance) in the packed cell, binned to 5 kb with the
genome coordinates of genome_map.py (chromosome arms; the rDNA array is not compared).
Data: Lazar-Stefanita et al. 2017 EMBO J 36:2684 (G1-arrested cells, replicates I-III summed; 5 kb; the model cell is
G1-like) and Swygert et al. 2019 Mol Cell 73:533 (log phase; ICE-normalised 5 kb). sacCer3 bins.
Measures (each for model and data):
  ps_slope        log-log slope of the cis contact frequency over 10-200 kb
  cis_fraction    share of contacts within chromosomes
  cen_enrichment  trans contacts between centromere-proximal bins (+-CEN_KB) / mean trans contact (Rabl clustering)
  tel_enrichment  the same for telomere-proximal bins
  trna_enrichment trans contacts between bins holding tRNA genes (SGD; centromere and telomere bins +-50 kb left out)
                  / mean trans contact among the same non-centromeric, non-telomeric bins (tRNA-gene clustering;
                  Thompson et al. 2003 saw it by FISH)
  chrom_pair_r    Spearman correlation of the 16 x 16 trans contact matrix (per chromosome pair, normalised by
                  chromosome sizes) between model and data
Normalisation (2026-11-12): the measures depend on it (Swygert's ICE matrix re-balanced on the nuclear chromosomes
gives cis 0.51 instead of 0.41, telomere clustering 3.6x instead of 26x). Each map, and the model, is therefore also
scored after the same iterative correction (ICE: all bins to equal coverage, diagonal and the 2% least covered bins
removed); keys with the suffix "_ice". The two maps' chromosome-pair patterns agree with each other at Spearman
-0.04 as provided and 0.15 balanced, which bounds what that measure can show.
Output: work/validation/hic_test.json
"""
import os, sys, json, gzip
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
from config import WORK, W3D, DATA
from wc3d import CHROM_KB, CEN_KB

G = f"{DATA}/validation/genome"
BIN = 5000
CONTACT_NM = 30.0
CEN_W, TEL_W = 20000, 20000


def bins_genome():
    """bin index for (chromosome 1-16, bin) from chromosome sizes in 5 kb bins"""
    nb = [int(np.ceil(kb * 1000 / BIN)) for kb in CHROM_KB]
    off = np.r_[0, np.cumsum(nb)]
    return off


def model_matrix():
    z = np.load(f"{W3D}/prepass.npz"); P = z["nuc_pos"].astype(float)
    d = pd.read_parquet(f"{WORK}/chromatin/nucleosome_genome.parquet").sort_values("nuc")
    off = bins_genome()
    b = off[d.chrom.values - 1] + np.minimum(d.bp.values // BIN, np.diff(off)[d.chrom.values - 1] - 1)
    pairs = cKDTree(P[d.nuc.values]).query_pairs(CONTACT_NM, output_type="ndarray")
    n = off[-1]
    M = np.zeros((n, n))
    np.add.at(M, (b[pairs[:, 0]], b[pairs[:, 1]]), 1); np.add.at(M, (b[pairs[:, 1]], b[pairs[:, 0]]), 1)
    return M, off


def lazar_matrix():
    files = [f for f in os.listdir(f"{G}/hic_lazarstefanita2017") if "G1_" in f]
    M = None
    for f in files:
        d = pd.read_csv(f"{G}/hic_lazarstefanita2017/{f}", index_col=0)
        M = d.values.astype(float) if M is None else M + d.values.astype(float)
    labels = list(d.columns)
    keep = np.array([l.split("~")[0].replace("chr", "").isdigit() for l in labels])      # drop chrM
    chrom = np.array([int(l.split("~")[0].replace("chr", "")) for l in np.array(labels)[keep]])
    mid = np.array([int(l.split("~")[1]) for l in np.array(labels)[keep]])
    return M[np.ix_(keep, keep)], chrom, mid


def swygert_matrix():
    bed = pd.read_csv(f"{G}/hic_swygert2019_log/GSE120605_sacCer3.5000.bins.annot.bed.gz", sep="\t", header=None)
    roman = {r: i + 1 for i, r in enumerate(["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"])}
    keep = bed[0].str.replace("chr", "").isin(roman)
    chrom = bed[0].str.replace("chr", "").map(roman).values
    mid = ((bed[1] + bed[2]) // 2).values
    n = len(bed)
    M = np.zeros((n, n))
    t = pd.read_csv(f"{G}/hic_swygert2019_log/GSM3405625_Log.merged.5000.ice.matrix.txt.gz", sep=r"\s+", header=None,
                    names=["i", "j", "raw", "ice"], engine="python")
    t["v"] = t.ice.fillna(t.raw)            # a few rows carry a single value
    i, j, v = t.i.values.astype(int), t.j.values.astype(int), t.v.values.astype(float)
    np.add.at(M, (i, j), v)
    off = i != j
    np.add.at(M, (j[off], i[off]), v[off])
    k = np.where(keep.values)[0]
    return M[np.ix_(k, k)], chrom[k].astype(int), mid[k]


def ice(M, n=100):
    """iterative correction: diagonal and the 2% least covered bins removed, then all bins scaled to equal coverage"""
    M = np.nan_to_num(np.array(M, float)); np.fill_diagonal(M, 0)
    r = M.sum(1); keep = r > np.percentile(r[r > 0], 2)
    M[~keep, :] = 0; M[:, ~keep] = 0
    for _ in range(n):
        s = M.sum(1); keep = keep & (s > 0); f = np.where(keep, s / s[keep].mean(), 1.0)
        M /= np.outer(f, f)
    return M


def _trna_bins(chrom, mid):
    import trna_telomeres as tt
    g = tt.trna_genes(); t = np.zeros(len(chrom), bool)
    for c, b in zip(g.chrom, g.bp):
        t[(chrom == c) & (np.abs(mid - b) <= BIN / 2 + 1)] = True
    return t


def measures(M, chrom, mid):
    n = len(chrom)
    I, J = np.triu_indices(n, 1)
    v = M[I, J]
    cis = chrom[I] == chrom[J]
    s = np.abs(mid[I] - mid[J])
    out = dict(cis_fraction=float(v[cis].sum() / v.sum()))
    edges = np.logspace(np.log10(5000), np.log10(400000), 16); midb = np.sqrt(edges[1:] * edges[:-1])
    ps = np.array([v[cis & (s >= a) & (s < b)].mean() if (cis & (s >= a) & (s < b)).any() else np.nan for a, b in zip(edges[:-1], edges[1:])])
    ok = (midb > 10000) & (midb < 200000) & (ps > 0)
    out["ps_slope"] = float(np.polyfit(np.log(midb[ok]), np.log(ps[ok]), 1)[0])
    out["ps"] = dict(s_bp=midb.round(0).tolist(), p=(ps / np.nanmax(ps)).tolist())
    cen = np.array([abs(m - CEN_KB[c - 1] * 1000) <= CEN_W for c, m in zip(chrom, mid)])
    tel = np.array([min(m, CHROM_KB[c - 1] * 1000 - m) <= TEL_W for c, m in zip(chrom, mid)])
    tr = ~cis
    out["cen_enrichment"] = float(v[tr & cen[I] & cen[J]].mean() / v[tr].mean())
    out["tel_enrichment"] = float(v[tr & tel[I] & tel[J]].mean() / v[tr].mean())
    tg = _trna_bins(chrom, mid)
    far = np.array([abs(m - CEN_KB[c - 1] * 1000) > 50000 and min(m, CHROM_KB[c - 1] * 1000 - m) > 50000 for c, m in zip(chrom, mid)])
    base = tr & far[I] & far[J]
    out["trna_enrichment"] = float(v[base & tg[I] & tg[J]].mean() / v[base].mean())
    C = np.zeros((16, 16))
    np.add.at(C, (chrom[I][tr] - 1, chrom[J][tr] - 1), v[tr])
    C = C + C.T
    sz = np.bincount(chrom - 1, minlength=16).astype(float)
    out["trans_norm"] = (C / np.outer(sz, sz)).tolist()
    return out


def main():
    Mm, off = model_matrix()
    chrom_m = np.concatenate([np.full(off[c + 1] - off[c], c + 1) for c in range(16)])
    mid_m = np.concatenate([np.arange(off[c + 1] - off[c]) * BIN + BIN // 2 for c in range(16)])
    res = {"model": measures(Mm, chrom_m, mid_m), "model_ice": measures(ice(Mm), chrom_m, mid_m)}
    for name, (M, ch, mid) in (("lazarstefanita2017_G1", lazar_matrix()), ("swygert2019_log", swygert_matrix())):
        res[name] = measures(M, ch, mid)
        res[name + "_ice"] = measures(ice(M), ch, mid)
    iu = np.triu_indices(16, 1)
    for name in ("lazarstefanita2017_G1", "swygert2019_log"):
        for sfx in ("", "_ice"):
            a = np.array(res["model" + sfx]["trans_norm"])[iu]; b = np.array(res[name + sfx]["trans_norm"])[iu]
            res[name + sfx]["chrom_pair_r_vs_model"] = float(spearmanr(a, b).correlation)
    for sfx in ("", "_ice"):
        a = np.array(res["lazarstefanita2017_G1" + sfx]["trans_norm"])[iu]; b = np.array(res["swygert2019_log" + sfx]["trans_norm"])[iu]
        res["swygert2019_log" + sfx]["chrom_pair_r_vs_lazar"] = float(spearmanr(a, b).correlation)
    for k in res:
        res[k].pop("trans_norm", None)
    json.dump(res, open(f"{WORK}/validation/hic_test.json", "w"), indent=1)
    for k, v in res.items():
        print(k, {a: (round(b, 3) if isinstance(b, float) else None) for a, b in v.items() if a != "ps"})
    return res


if __name__ == "__main__":
    main()
