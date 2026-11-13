"""mRNAs translated on the mitochondrial surface (IMPROVEMENTS M4).

Which mRNAs: the genes whose translating ribosomes are enriched next to the mitochondrial outer membrane in
proximity-specific ribosome profiling without cycloheximide (Williams, Jan & Weissman 2014 Science 346:748, Table S1
'minus.chx.enriched'; 223 genes; SD glucose, OD 0.4-0.6). Every polysome of these genes is tethered to the outer
membrane: its first ribosome starts in the cytosol MITO_SEED_NM from the OM and every ribosome centre stays within
MITO_TMAX_NM of it (polysome.walk). Free (untranslated) mRNAs of these genes are not tethered.
Held out: Saint-Georges et al. 2008 (fractionation, mRNA % in the mitochondrial fraction; validation/mrna_localization.py).

mRNAs translated on the ER (IMPROVEMENTS M5, 2026-11-12): the secretome genes (signal sequence or transmembrane
domain) whose translating ribosomes are enriched next to the Sec61 translocon in the same study design (Jan, Williams &
Weissman 2014 Science 346:1257521, Table S5: BirA-Sec63 with cycloheximide, log2 enrichment above the 95th percentile
of non-secretome genes, their definition of ER-enriched; 873 genes; SD glucose). Their polysomes are tethered to the
cortical ER the same way (first ribosome ER_SEED_NM from the ER membrane, all within MITO_TMAX_NM). The translocon-
docked single ribosomes of er_ribosomes.py stay as they are (not joined into these chains). Genes in both sets stay
mitochondrial. Held out: membrane-polysome enrichment per gene by fractionation in YPD (Chartron, Hunt & Frydman 2016
Nature 536:224, Supplementary Table 1).
"""
import numpy as np
import numba as nb
import pandas as pd
import cell3d as c3
from config import DATA

WILLIAMS = f"{DATA}/validation/mrna/williams2014_mito_prox/GSE61011_TableS1.txt.gz"
MITO_SEED_NM = (16.0, 30.0)     # first ribosome centre to OM midsurface (80S radius ~13 nm plus the TOM complex)
MITO_TMAX_NM = 45.0
MITO_TETHER = True
JAN = f"{DATA}/validation/mrna/jan2014_er_prox/GSE61012_TableS5.txt.gz"
ER_SEED_NM = (16.0, 30.0)
ER_TETHER = True


def mito_genes():
    w = pd.read_csv(WILLIAMS, sep="\t", skiprows=28)
    w.columns = [c.lstrip("#") for c in w.columns]
    return set(w.gene[w["minus.chx.enriched"] == "yes"].str.upper())


def er_genes():
    j = pd.read_csv(JAN, sep="\t", skiprows=33)
    j.columns = [c.lstrip("#") for c in j.columns]
    e = j["s63.7mchx.log2enrich"]
    thr = np.nanpercentile(e[j.secretome == "no"].dropna(), 95)
    return set(j.gene[(j.secretome == "yes") & (e > thr)].str.upper())


@nb.njit(cache=True, parallel=True)
def mem_dist(P, m, G, pores, pidx, PT, PQ, PC):
    d = np.empty(P.shape[0]); lab = np.empty(P.shape[0], np.int8)
    for i in nb.prange(P.shape[0]):
        out = np.zeros(c3.NMEM)
        mx = c3.mem_sdfs(P[i, 0], P[i, 1], P[i, 2], G, pores, pidx, PT, PQ, PC, out)
        d[i] = out[m]
        lab[i] = c3.label_from(P[i, 0], P[i, 1], P[i, 2], G, pores, out, mx)
    return d, lab


def surface_seeds(cell, n, rng, mem=c3.M_OM, prim_class=c3.C_OM, dist=MITO_SEED_NM, away=None):
    """n cytosolic points whose distance to membrane `mem` lies in `dist`, uniform over the membrane's surroundings
    (sampled in the bounding spheres of its primitives); away = (membrane, nm): also at least that far from it (the
    cortical ER's side facing the plasma membrane is a ~30 nm gap without room for ribosomes)"""
    k = np.where(cell.PC == prim_class)[0]
    if not len(k) or n == 0:
        return np.zeros((0, 3))
    G = cell.G; glo = np.full(3, -G[0] - 10.0)
    a = cell.args(glo, -glo, margin=50.0)
    w = (cell.br[k] + dist[1]) ** 3
    got = []
    while sum(len(g) for g in got) < n:
        i = k[rng.choice(len(k), 200000, p=w / w.sum())]
        u = rng.normal(size=(len(i), 3)); u /= np.linalg.norm(u, axis=1)[:, None]
        P = cell.bc[i] + u * ((cell.br[i] + dist[1]) * rng.random(len(i)) ** (1 / 3))[:, None]
        d, lab = mem_dist(P, mem, *a)
        ok = (lab == c3.CYTO) & (np.abs(d) >= dist[0]) & (np.abs(d) <= dist[1])
        if away is not None:
            d2, _ = mem_dist(P, away[0], *a)
            ok &= np.abs(d2) >= away[1]
        got.append(P[ok])
    P = np.concatenate(got)
    return P[rng.choice(len(P), n, replace=False)]


def tether(cell, chains, rng, log=print):
    """per polysome chain (mrna_genes.transcripts rows with ribosomes > 0): (membrane index or -1, seed point,
    max distance) for polysome.place"""
    tm = np.full(len(chains), -1, np.int64); seeds = np.zeros((len(chains), 3))
    on = np.zeros(len(chains), bool)
    if MITO_TETHER:
        on = chains.orf.str.upper().isin(mito_genes()).values
        tm[on] = c3.M_OM
        seeds[on] = surface_seeds(cell, int(on.sum()), rng)
        log(f"mRNA localization: {int(on.sum())} polysomes of {chains.orf[on].nunique()} genes tethered to the "
            f"mitochondrial outer membrane (Williams 2014, -CHX)")
    if ER_TETHER:
        er = chains.orf.str.upper().isin(er_genes()).values & ~on
        tm[er] = c3.M_ER
        seeds[er] = surface_seeds(cell, int(er.sum()), rng, mem=c3.M_ER, prim_class=c3.C_CER, dist=ER_SEED_NM, away=(c3.M_PM, 45.0))
        log(f"mRNA localization: {int(er.sum())} polysomes of {chains.orf[er].nunique()} genes tethered to the cortical "
            f"ER (Jan 2014, Sec63 proximity, +CHX)")
    return tm, seeds, MITO_TMAX_NM
