"""Chromatin-bound proteins at their measured genomic sites (IMPROVEMENTS K5).

Sites: ChExMix peaks of the ChIP-exo survey of about 400 chromatin proteins (Rossi et al. 2021 Nature 592:309; YPD,
log phase; sacCer3; filtered peak files of the paper's GitHub repository, 188,489 peaks). Each protein that is a
library type of its own in the nucleoplasm (293 proteins, 129,000 peaks; subunits of complexes the inventory does
not assemble, such as TFIIH or Mediator, are placed one by one, so they gather at the same promoters) puts
min(peaks, F_BOUND x nuclear copies) copies on its sites, drawn without replacement with probability proportional
to the peak's tag count. F_BOUND, the chromatin-bound share, is an assumption (single-molecule tracking finds bound
fractions of tens of percent for yeast TFs; unchecked). The rest stay free in the nucleoplasm.
Geometry: the DNA point of a site is interpolated between the two measured nucleosomes that flank it (by bp); the
protein (random orientation) comes within TOUCH_NM of it with one of its spheres, on the side away from the local
nucleosome mass, and overlaps nothing else. Requires the measured
nucleosome map (prepass nuc_bp).
"""
import os
import numpy as np
import pandas as pd
import pack3d as p3
import cell3d as c3
from config import DATA
from config import seed as _seed

PEAKS = f"{DATA}/validation/genome/chipexo_rossi2021/github_repo/04_ChExMix_Peaks"
F_BOUND = 0.3
TOUCH_NM, TRIES = 5.0, 60
SITES = True


def peaks():
    """DataFrame gene, chrom (1-16), bp, score over all peak files"""
    rows = []
    for f in sorted(os.listdir(PEAKS)):
        if not f.endswith(".bed") or os.path.getsize(f"{PEAKS}/{f}") == 0:
            continue
        gene = (f.split("_CX")[0] if "_CX" in f else f.split("_")[1]).upper()
        d = pd.read_csv(f"{PEAKS}/{f}", sep="\t", header=None, usecols=[0, 1, 4], names=["chr", "bp", "score"])
        d = d[d.chr.str.match(r"^chr\d+$")]
        rows.append(pd.DataFrame(dict(gene=gene, chrom=d.chr.str[3:].astype(int), bp=d.bp.astype(np.int64), score=d.score.astype(float))))
    return pd.concat(rows, ignore_index=True).drop_duplicates(["gene", "chrom", "bp"])


def place(cell, lp, pre, pools, obstacles, seed=53, log=print, sites=None, prefix="tf", what="ChIP-exo peaks"):
    """pools: {lib_id: free nuclear copies}. Returns pre-pass group `prefix` (tf_type, tf_pos, tf_q) and per-type
    counts. sites: DataFrame gene, chrom, bp, score of chosen sites (one copy per site, as many as the pool allows;
    trna_telomeres.py); None: the ChIP-exo peaks with the F_BOUND share."""
    from scipy.spatial import cKDTree
    from scipy.spatial.transform import Rotation
    import genome_map
    ids = list(lp["ids"])
    if not SITES or "nuc_bp" not in pre:
        return {}, {}
    seed = _seed(seed); rng = np.random.default_rng(seed)
    pk = peaks() if sites is None else sites
    pk = pk[pk.gene.isin([g for g in pools if g in ids])]
    # per segment: nucleosome indices sorted by bp
    pos, ch, bp = pre["nuc_pos"].astype(np.float64), pre["nuc_chain"], pre["nuc_bp"]
    segs = genome_map.segments()
    seg_of = {}
    for s, (c, kind, s0, s1) in enumerate(segs):
        I = np.where(ch == s)[0]; o = np.argsort(bp[I])
        seg_of.setdefault(c + 1, []).append((min(s0, s1), max(s0, s1), I[o], bp[I][o], kind))

    def dna_point(c, x):
        for lo, hi, I, g, kind in seg_of.get(c, ()):
            if not len(g):
                continue
            if lo <= x < hi or (kind == "L" and x < hi) or (kind in ("R", "XIIR_out") and x >= lo):
                j = np.searchsorted(g, x)
                if j == 0 or j == len(g):           # beyond the last nucleosome (chromosome ends): at that nucleosome
                    return pos[I[0 if j == 0 else -1]]
                f = (x - g[j - 1]) / max(1, g[j] - g[j - 1])
                return pos[I[j - 1]] + f * (pos[I[j]] - pos[I[j - 1]])
        return None

    cn = cell.G[3:6]; rn = float(cell.G[11:14].max()); glo, ghi = cn - rn - 80, cn + rn + 80
    a = cell.args(glo, ghi, margin=200.0)
    OC, OR = obstacles
    grid = p3.LocalGrid(glo, ghi, OC.astype(np.float64), OR.astype(np.float64), n_new_spheres=4_000_000)
    ntree = cKDTree(pos)
    out = np.zeros(c3.NMEM)

    types, P_, Q_, counts, wanted = [], [], [], {}, 0
    for gene, g in pk.groupby("gene"):
        if sites is None:
            n = int(min(len(g), np.floor(F_BOUND * pools[gene] + rng.random())))
        else:
            n = int(min(len(g), np.floor(pools[gene])))
        if n <= 0:
            continue
        wanted += n
        li = ids.index(gene)
        C = lp["STC"][lp["st_off"][li]:lp["st_off"][li + 1]].astype(np.float64)
        R = lp["STR"][lp["st_off"][li]:lp["st_off"][li + 1]].astype(np.float64)
        w = g.score.clip(lower=1.0).values          # a few peaks carry score 0
        pick = g.iloc[rng.choice(len(g), n, replace=False, p=w / w.sum())]
        got = 0
        for r in pick.itertuples():
            E = dna_point(r.chrom, r.bp)
            if E is None:
                continue
            # away from the local nucleosome mass (the site faces the open side of the fibre)
            nb = ntree.query_ball_point(E, 25.0)
            out_u = E - pos[nb].mean(0) if nb else np.zeros(3)
            out_u = out_u / (np.linalg.norm(out_u) + 1e-9)
            for _ in range(TRIES):
                Rm = Rotation.random(random_state=rng).as_matrix()
                S0 = C @ Rm.T; kk = rng.integers(len(S0))
                u = rng.normal(size=3) + 1.5 * out_u; u /= np.linalg.norm(u)
                t = E - S0[kk] + u * (R[kk] + 1.0 + rng.uniform(0, TOUCH_NM))      # DNA radius 1 nm
                S = S0 + t
                if grid.free(S, R) and p3.check_volume(S, R, t[0], t[1], t[2], float(lp["BR"][li]), 1 << c3.NUCP, *a, out):
                    grid.add(S, R)
                    types.append(li); P_.append(t); Q_.append(np.roll(Rotation.from_matrix(Rm).as_quat(), 1)); got += 1
                    break
        counts[gene] = got
    if sites is None:
        log(f"chromatin sites: {len(types)} of {wanted} copies of {len(counts)} proteins on their ChIP-exo peaks "
            f"(Rossi 2021; bound share {F_BOUND:.0%} of nuclear copies, at most one per peak)")
    else:
        log(f"{what}: {len(types)} of {wanted} copies (" + ", ".join(f"{k} {v}" for k, v in counts.items()) + ")")
    if not types:
        return {}, counts
    return {f"{prefix}_type": np.array(types, np.int32), f"{prefix}_pos": np.array(P_, np.float32),
            f"{prefix}_q": np.array(Q_, np.float32)}, counts
