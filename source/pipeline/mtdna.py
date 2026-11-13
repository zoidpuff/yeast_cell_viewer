"""Mitochondrial DNA nucleoids (literature item L16; sites and numbers revised 2026-11-13).

A haploid cell in glucose has 1.2-1.35 nucleoids per um of mitochondrial tubule, about 0.8 um apart, with one within
0.5 um of 69% of the tubule tips, and 1-2 copies of the 85.8 kb genome each (Osman et al. 2015 PNAS; Goke et al. 2020;
Seel et al. 2023). The sites are chosen with the geometry (cell3d: NUCLEOID_PER_UM, NUCLEOID_MIN_GAP), which keeps
cristae NUCLEOID_CLEAR away from them: a nucleoid is a compact body in crista-free matrix near the inner boundary
membrane (size: Goke 2020; position: by analogy with mammalian cryo-ET, Kukat et al. 2015; no yeast EM places nucleoids
relative to cristae). Until 2026-11-13 the 15 nucleoids were put anywhere in the matrix and their DNA threaded through
the gaps between cristae stacked every 42 nm. Each nucleoid is a cluster of LOOP_KB DNA loops: self-avoiding walks
(10 bp per 3.4 nm step) that start near the centre and are drawn back to it (pack3d.arm_walks). The DNA is drawn at
2 bp per bead, like the linker DNA. Abf2, which packages it, is packed free in the matrix.
"""
import numpy as np
from config import seed as _seed
import pack3d as p3
import cell3d as c3

N_NUCLEOIDS, GENOME_KB = 15, 85.8           # N_NUCLEOIDS: only for a geometry without nucleoid sites
GENOMES_PER_FL = 0.45                        # mtDNA copies per fL of cell: 21 +- 4 per haploid cell in YPD by qPCR (Goke
                                             # et al. 2020), about 0.45 per fL (Seel et al. 2023), both W303; was 35 per cell
LOOP_KB = 5.0
STEP, BP_PER_STEP = 3.4, 10
DNA_R = 1.0
DNA_STEP = 0.68                                    # 2 bp per drawn bead


def place(cell, obstacles=None, seed=81, log=print):
    seed = _seed(seed); rng = np.random.default_rng(seed)
    segs = np.array(cell.d["mito_segments"], np.float64)
    seg_len = np.linalg.norm(segs[:, 1] - segs[:, 0], axis=1)
    sites = cell.d.get("mito_nucleoids")
    if sites:
        # the geometry chose the sites and kept cristae away from them (cell3d, NUCLEOID_*)
        centres = [np.array(s_[0]) for s_ in sites if cell.label(np.array(s_[0]))[0] == c3.MATRIX]
    else:
        centres = []
        for _ in range(20000):
            if len(centres) >= N_NUCLEOIDS:
                break
            k = rng.choice(len(segs), p=seg_len / seg_len.sum())
            p = segs[k, 0] + (segs[k, 1] - segs[k, 0]) * rng.uniform(0.1, 0.9) + rng.normal(0, 25.0, 3)
            if cell.label(p)[0] != c3.MATRIX or any(np.linalg.norm(p - c) < 500.0 for c in centres):
                continue
            centres.append(p)
    centres = np.array(centres)
    import inventory as inv
    GENOMES = max(len(centres), int(round(GENOMES_PER_FL * inv.CELL_VOL)))
    genomes = np.bincount(rng.integers(0, len(centres), GENOMES - len(centres)), minlength=len(centres)) + 1
    loops = [(i, int(round(g * GENOME_KB / LOOP_KB))) for i, g in enumerate(genomes)]
    steps = int(LOOP_KB * 1000 / BP_PER_STEP)
    starts = np.array([centres[i] for i, n in loops for _ in range(n)])
    L = np.full(len(starts), steps, np.int64)
    G = cell.G; glo = centres.min(0) - 600; ghi = centres.max(0) + 600
    a = cell.args(glo, ghi, margin=100.0)
    OC, OR = obstacles if obstacles is not None else (np.zeros((0, 3)), np.zeros(0))
    k = np.all((OC > glo) & (OC < ghi), axis=1); OC, OR = OC[k].astype(np.float64), OR[k].astype(np.float64)
    nsph = int(L.sum()) + len(OR)
    pos, _, chain, brk = p3.arm_walks(seed, glo, ghi, *a, np.zeros((1, 3)), np.array([DNA_R]), DNA_R, 1 << c3.MATRIX,
                                      OC, OR, L, starts, starts, 40.0, STEP, int(nsph * 1.1) + 1000, int(nsph * 10) + 10000)
    # 2-bp beads along each unbroken piece of each loop
    pts, off = [], [0]
    cut = np.r_[np.where((np.diff(chain) != 0) | (brk[1:] == 1))[0] + 1, len(pos)]
    a0 = 0
    for b0 in cut:
        P = pos[a0:b0].astype(np.float64); a0 = b0
        if len(P) < 2:
            continue
        seg = np.linalg.norm(np.diff(P, axis=0), axis=1); s = np.r_[0, np.cumsum(seg)]
        t = np.arange(0, s[-1], DNA_STEP)
        Q = np.c_[[np.interp(t, s, P[:, j]) for j in range(3)]].T
        pts.append(Q); off.append(off[-1] + len(Q))
    pts = np.concatenate(pts).astype(np.float32)
    rad = [np.percentile(np.linalg.norm(pos[np.isin(chain, np.where(np.repeat(np.arange(len(loops)), [n for _, n in loops]) == i)[0])] - centres[i], axis=1), 90)
           for i in range(len(centres))]
    log(f"mtDNA: {len(centres)} nucleoids, {GENOMES} genomes ({genomes.min()}-{genomes.max()} per nucleoid), "
        f"{len(pos) * BP_PER_STEP / 1e6:.2f} of {L.sum() * BP_PER_STEP / 1e6:.2f} Mb placed in {len(off) - 1} pieces; "
        f"90% of each nucleoid's DNA within {np.median(rad):.0f} nm (median) of its centre")
    return dict(mtdna_pts=pts, mtdna_off=np.array(off, np.int64), mtdna_centres=centres.astype(np.float32),
                mtdna_genomes=genomes.astype(np.int32))
