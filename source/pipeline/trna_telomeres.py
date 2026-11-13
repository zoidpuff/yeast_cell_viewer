"""Chromatin landmarks the ChIP-exo survey leaves out (IMPROVEMENTS N9): tRNA genes, telomeres and silent loci.

The peak files of Rossi et al. 2021 exclude tRNA genes and telomeres, so chromatin_sites.py placed no Pol III there and
no Rap1 or Sir proteins at telomeres. Here, with the same geometry (chromatin_sites.place, one copy per site):
  tRNA genes      the 275 nuclear tRNA genes of SGD (sacCer3), each with one engaged Pol III and one copy of each
                  TFIIIC subunit in the library (Tfc1, Tfc7; the other subunits have no library type). Pol III
                  occupies nearly all tRNA genes in growing cells; one polymerase per gene is an assumption (unchecked).
  telomeres       RAP1_PER_TEL Rap1 at the end of each of the 32 arms: the terminal TG1-3 tracts, about 300 bp, carry a
                  Rap1 site every 18 bp on average (Marcand et al. 1997 Science 275:986), so about 17 (was 15, unchecked)
  silent chromatin  per Sir-bound nucleosome: where Sir3 ChIP-seq is enriched 5-fold over input within 30 kb of a
                  chromosome end or in HML / HMR (silent_chromatin.py, 2026-11-12; 238 nucleosomes at 28 of the 32
                  ends; held-out H4K16 acetylation and Sir4 ChIP agree). Before, and still where the ChIP tracks are
                  missing: every nucleosome within SIR_KB of each telomere and over HML and HMR (529 nucleosomes, of
                  which only 28% are H4K16-hypoacetylated). SIR_PER_NUC copies per nucleosome, two Sir3 (Swygert et al. 2014 Nat Commun 5:4751) and one
                  Sir2-Sir4 (1-2 per nucleosome: Swygert et al. 2018 PNAS 115:12447). Both are saturating stoichiometries
                  on reconstituted arrays, so in cells they are upper bounds (before 2026-11-12: one of each, unchecked).
                  Sir2's nucleolar pool (RENT, rDNA) is not touched.
Copies come out of the free nuclear pools; the rest stay free. Nothing here is fitted; the expected consequences
(Pol III next to the nucleolus, Rap1/Sir foci at the envelope) are checks.
"""
import os
import numpy as np
import pandas as pd
from config import DATA

POL3, TFIIIC = "PolIII", ["TFC1", "TFC7"]
RAP1_PER_TEL = 17
SIR_PER_NUC = {"SIR3": 2, "SIR4": 1, "SIR2": 1}
SIR_KB = 3.0


def trna_genes():
    f = pd.read_csv(f"{DATA}/sgd/SGD_features.tab", sep="\t", header=None, dtype=str)
    t = f[(f[1] == "tRNA gene") & f[8].isin([str(i) for i in range(1, 17)])]
    a, b = t[9].astype(int), t[10].astype(int)
    return pd.DataFrame(dict(name=t[3].values, chrom=t[8].astype(int).values, bp=((a + b) // 2).values))


def silent_loci():
    f = pd.read_csv(f"{DATA}/sgd/SGD_features.tab", sep="\t", header=None, dtype=str)
    t = f[f[1] == "silent mating type cassette array"]
    return [(int(c), min(int(a), int(b)), max(int(a), int(b))) for c, a, b in zip(t[8], t[9], t[10])]


def sites(pre):
    """DataFrame gene, chrom, bp, score: one row per copy to place"""
    import genome_map
    from wc3d import CHROM_KB
    rows = []
    tg = trna_genes()
    for g in [POL3] + TFIIIC:
        rows.append(pd.DataFrame(dict(gene=g, chrom=tg.chrom, bp=tg.bp, score=1.0)))
    # telomere ends: bp 0 (left) and the chromosome end (right); several Rap1 per end at the same point
    ends = [(c, 0) for c in range(1, 17)] + [(c, CHROM_KB[c - 1] * 1000) for c in range(1, 17)]
    rows.append(pd.DataFrame([dict(gene="RAP1", chrom=c, bp=b, score=1.0) for c, b in ends for _ in range(RAP1_PER_TEL)]))
    # silent chromatin: placed nucleosomes within SIR_KB of a chromosome end, and on HML/HMR
    segs = genome_map.segments()
    chrom = np.array([segs[s][0] + 1 for s in pre["nuc_chain"]]); bp = pre["nuc_bp"]
    L = np.array([CHROM_KB[c - 1] * 1000 for c in chrom])
    import silent_chromatin as sc
    if os.path.exists(f"{sc.DIR}/{sc.REPS[0][0]}_median_normalized.bedGraph.gz"):
        sel = sc.sir_bound(chrom, bp)
    else:
        sel = (bp < SIR_KB * 1000) | (bp > L - SIR_KB * 1000)
        for c, a, b in silent_loci():
            sel |= (chrom == c) & (bp >= a) & (bp <= b)
    for g, n in SIR_PER_NUC.items():
        rows.append(pd.DataFrame(dict(gene=g, chrom=np.repeat(chrom[sel], n), bp=np.repeat(bp[sel], n), score=1.0)))
    return pd.concat(rows, ignore_index=True)
