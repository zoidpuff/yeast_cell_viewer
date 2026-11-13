"""Gene-specific mRNAs (IMPROVEMENTS M1-M2): every mRNA molecule is a transcript of a named gene, with the gene's
length (5' UTR + CDS + 3' UTR + poly(A)) and a ribosome load proportional to its CDS.

Inputs
  copies per gene   mRNA molecules per cell (RA) of Pelechano et al. 2010 PLoS One 5:e15442, Table S1 (compiled from
                    Nagalakshmi et al. 2008 and Miura et al. 2008); 5,085 genes, 25,900 mRNAs per cell
  CDS               summed CDS segments of each SGD ORF (R64), so introns are left out (until 2026-11-12 the ORF span,
                    which counted the introns of 208 genes, mostly ribosomal-protein genes, as coding sequence)
  UTRs, poly(A)     per gene: TIF-seq UTRs in YPD (Pelechano et al. 2013 Nature 497:127) and PAL-seq tails (Subtelny
                    et al. 2014 Nature 508:66); genes without data take UTR5_NT, UTR3_NT, POLYA_NT. TIF-seq UTRs
                    are genomic distances, so a 5' UTR intron (SGD five_prime_UTR_intron, 21 genes here) is
                    subtracted where the UTR is longer than the intron (the transcript starts upstream of it)
  ribosome load     ribosomes on an mRNA ~ Poisson(CDS / spacing_g), with one global spacing chosen so that the
                    translated mRNAs carry the polysome ribosome budget, modulated by the per-gene translational
                    efficiency from ribosome profiling (Subtelny et al. 2014). Per-gene polysome profiles (Arava et al.
                    2003) are not used: they remain a test. mRNAs drawing zero ribosomes are free (untranslated).
Output of transcripts(): DataFrame, one row per mRNA molecule: orf, gene, cds, utr5, utr3, polya, length, ribosomes.
"""
import os
import numpy as np
import pandas as pd

from config import DATA

PELECHANO = f"{DATA}/pelechano2010/files/pone.0015442.s002.xls"
SGD = f"{DATA}/sgd/SGD_features.tab"
UTR_TABLE = f"{DATA}/mrna/utr_per_gene.tsv"          # TIF-seq, YPD (Pelechano 2013); data/mrna/README.md
POLYA_TABLE = f"{DATA}/mrna/polya_per_gene.tsv"      # PAL-seq, total RNA (Subtelny 2014)
TE_TABLE = f"{DATA}/mrna/te_per_gene.tsv"            # ribosome profiling TE (Subtelny 2014)
UTR5_NT, UTR3_NT, POLYA_NT = 50, 120, 33     # where no per-gene value exists (near the data medians)
MAX_RIBO = 40             # cap on ribosomes per mRNA (the walk stops earlier when it runs out of space)
GENE_MRNA = True          # False: the earlier generic mRNAs (lognormal lengths, Poisson polysome sizes)


def gene_table():
    import warnings
    warnings.filterwarnings("ignore")
    p = pd.read_excel(PELECHANO, engine="xlrd")
    p = pd.DataFrame(dict(orf=p.iloc[:, 0].astype(str).str.upper(), gene=p.iloc[:, 1].astype(str),
                          copies=pd.to_numeric(p.iloc[:, 2], errors="coerce"))).dropna(subset=["copies"])
    sgd = pd.read_csv(SGD, sep="\t", header=None, dtype=str)

    def seg_nt(kind):          # summed length of a feature's segments per parent ORF
        x = sgd[sgd[1] == kind]
        n = (pd.to_numeric(x[10], errors="coerce") - pd.to_numeric(x[9], errors="coerce")).abs() + 1
        return n.groupby(x[6].str.upper().values).sum()
    f = sgd[sgd[1] == "ORF"]
    a, b = pd.to_numeric(f[9], errors="coerce"), pd.to_numeric(f[10], errors="coerce")
    span = pd.Series((abs(b - a) + 1).values, index=f[3].str.upper().values)
    cds = pd.DataFrame(dict(orf=span.index, cds=seg_nt("CDS").reindex(span.index).fillna(span).values))
    g = p.merge(cds, on="orf", how="inner")
    g["gene"] = np.where(g.gene.isin(["nan", ""]), g.orf, g.gene)
    g["utr5"], g["utr3"], g["polya"], g["te"] = float(UTR5_NT), float(UTR3_NT), float(POLYA_NT), 1.0
    for tab, cols in ((UTR_TABLE, ("utr5", "utr3")), (POLYA_TABLE, ("polya",)), (TE_TABLE, ("te",))):
        if tab is not None and os.path.exists(tab):
            t = pd.read_csv(tab, sep="\t")
            t["orf"] = t.orf.str.upper()
            g = g.drop(columns=list(cols)).merge(t[["orf", *cols]], on="orf", how="left")
            for c in cols:
                g[c] = g[c].fillna({"utr5": UTR5_NT, "utr3": UTR3_NT, "polya": POLYA_NT, "te": 1.0}[c])
    ui = g.orf.map(seg_nt("five_prime_UTR_intron")).fillna(0.0)
    g["utr5"] = np.where(g.utr5 > ui, g.utr5 - ui, g.utr5)
    return g


def transcripts(rng, n_ribo, n_total=None, log=print):
    """one row per mRNA molecule. n_ribo: ribosomes to put on mRNAs (the polysome budget); n_total: scale the
    per-gene copies to this many mRNAs (None: the data's total times inventory.CELL_SHARE, the size of this cell
    relative to the population-average cell the copies were measured for)"""
    import inventory
    g = gene_table()
    scale = inventory.CELL_SHARE if n_total is None else n_total / g.copies.sum()
    lam = g.copies.values * scale
    n = np.floor(lam).astype(np.int64) + (rng.random(len(lam)) < lam - np.floor(lam))      # stochastic rounding
    m = g.loc[np.repeat(np.arange(len(g)), n)].reset_index(drop=True)
    # one global spacing (nt of CDS per ribosome) so that the expected ribosomes on all mRNAs equal n_ribo
    w = (m.cds * m.te).values
    spacing = w.sum() / max(n_ribo, 1.0)
    m["ribosomes"] = np.minimum(MAX_RIBO, rng.poisson(w / spacing))
    m["length"] = (m.utr5 + m.cds + m.utr3 + m.polya).round().astype(np.int64)
    log(f"mRNA genes: {len(m)} mRNAs of {int((n > 0).sum())} genes (data total {g.copies.sum():.0f}, scale {scale:.2f}); "
        f"ribosome spacing {spacing:.0f} nt of CDS; {int((m.ribosomes > 0).sum())} translated, mean {m.ribosomes[m.ribosomes > 0].mean():.1f} "
        f"ribosomes; median length {int(m.length.median())} nt")
    return m
