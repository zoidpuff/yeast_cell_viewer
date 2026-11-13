"""tRNA genes at the nucleolus (IMPROVEMENTS P6, test; nothing in the model was set from these data).

Data (FISH, log-phase BY4741 and W303 backgrounds):
  - Thompson et al. 2003 Science 302:1399 (doi 10.1126/science.1089814): the active SUP53 tRNA gene overlaps the
    nucleolar signal in 52% of 440 cells, the same gene with an inactivated promoter in 13% of 715 cells.
  - Haeusler et al. 2008 Genes Dev 22:2204 (doi 10.1101/gad.1675908): in wild type 83% of nuclei show a single
    tRNA-gene cluster, 76% of these at the nucleolus (so about 63% of nuclei have a nucleolar tRNA cluster).
Model: every nuclear tRNA gene of SGD (275) placed on the chromatin by its bp position between the two nearest
measured nucleosomes of its arm (the same interpolation as chromatin_sites.dna_point), its distance to the
nucleolus (the sphere of geometry G[7:11] inside the nucleus; negative inside), against all nucleosomes as the
random-locus baseline. "Overlap" is a centre within OVERLAP_NM of the nucleolus, about the radius of a
diffraction-limited FISH spot (set before comparing). One model cell stands for a population of loci, not of cells.
Output: work/validation/trna_nucleolus.json
"""
import os, sys, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
from config import WORK, W3D  # noqa: E402

OVERLAP_NM = 250.0


def locus_positions(z, chrom, bp):
    """positions of genome loci (chromosome 1-16, bp) between the two nearest placed nucleosomes of their arm"""
    import genome_map
    segs = genome_map.segments()
    nchrom = np.array([segs[s][0] + 1 for s in z["nuc_chain"]]); nbp = z["nuc_bp"]; P = z["nuc_pos"].astype(float)
    out = np.full((len(bp), 3), np.nan)
    for c in np.unique(chrom):
        k = np.where(nchrom == c)[0]
        if not len(k):
            continue
        o = np.argsort(nbp[k]); kb, kp = nbp[k][o], P[k][o]
        for i in np.where(chrom == c)[0]:
            j = np.searchsorted(kb, bp[i])
            if j == 0 or j == len(kb):
                out[i] = kp[min(j, len(kb) - 1)]
            else:
                f = (bp[i] - kb[j - 1]) / max(kb[j] - kb[j - 1], 1)
                out[i] = kp[j - 1] + f * (kp[j] - kp[j - 1])
    return out


def main():
    import cell3d as c3, trna_telomeres as ttl
    cell = c3.Cell.load(f"{W3D}/geometry.json")
    z = np.load(f"{W3D}/prepass.npz")
    G = cell.G; cen, r = np.asarray(G[7:10], float), float(G[10])
    t = ttl.trna_genes()
    X = locus_positions(z, t.chrom.values, t.bp.values)
    ok = ~np.isnan(X[:, 0])
    d_trna = np.linalg.norm(X[ok] - cen, axis=1) - r
    d_nuc = np.linalg.norm(z["nuc_pos"].astype(float) - cen, axis=1) - r
    import pandas as pd
    from config import DATA
    f = pd.read_csv(f"{DATA}/sgd/SGD_features.tab", sep="\t", header=None, dtype=str)
    sup53 = t.name.isin(f.loc[f[4] == "SUP53", 3]).values          # standard name -> systematic name
    res = dict(
        overlap_nm=OVERLAP_NM, trna_genes=int(ok.sum()),
        model_trna_within_overlap=float(np.mean(d_trna <= OVERLAP_NM)),
        model_nucleosomes_within_overlap=float(np.mean(d_nuc <= OVERLAP_NM)),
        model_trna_median_nm=float(np.median(d_trna)), model_nucleosome_median_nm=float(np.median(d_nuc)),
        model_sup53_nm=float(d_trna[np.where(sup53[ok])[0][0]]) if sup53[ok].any() else None,
        thompson2003_active_sup53=0.52, thompson2003_inactive_sup53=0.13,
        haeusler2008_nuclei_with_nucleolar_cluster=round(0.83 * 0.76, 2))
    print(json.dumps(res, indent=1))
    json.dump(res, open(f"{WORK}/validation/trna_nucleolus.json", "w"), indent=1)
    return res


if __name__ == "__main__":
    main()
