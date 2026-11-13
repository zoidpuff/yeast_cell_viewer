"""Polysomes (inventory v3, item C1): cytosolic 80S ribosomes threaded on mRNAs.

Geometry: the mRNA channel of the 80S library structure (4V88, which has no mRNA) is taken from PDB 7ZW0 (yeast 80S
with 29 nt of mRNA, superposed on 4V88 by matching UniProt residues; median deviation 0.16 nm). The channel's 5' end is
the exit and its 3' end the entry. A polysome is grown as a self-avoiding chain: each next ribosome is placed so that
its mRNA entry lies a short linker (LINK_MIN..LINK_MAX nm) beyond the previous ribosome's exit, inside a cone around the
exit direction, with a random spin; ribosomes must lie in the cytosol and not clash. The mRNA polyline of each
ribosome is its channel segment plus the linker to the next ribosome (schematic beads, 0.5 nm apart, like linker DNA).
Polysome ribosomes are placed before packing (like NPCs and nucleosomes) and replace free cytosolic ribosomes.

Assumptions: POLY_FRACTION of cytosolic ribosomes are in polysomes; polysome sizes 1 + Poisson(MEAN_SIZE - 1)
(capped at MAX_SIZE); each polysome is one mRNA, so free mRNA threads are reduced by the number of polysomes.

Collided ribosomes (2026-11-12, IMPROVEMENTS P11): where the drawn free mRNA between two neighbours is at most
COLLIDE_NT nucleotides, the trailing (5') ribosome is placed in the pose it has in a collided disome, relative to the
leading one: the yeast disome 6I7O (Ikeuchi et al. 2019 EMBO J 38:e100276), both ribosomes superposed on the library
80S by sequence-aligned rRNA phosphates (disome_pose). The coarse sphere trees overlap in that pose, so the trailing
ribosome is moved away along the line between the centres until they clear (about 4 nm; centre distance 25.6 instead
of 21.6 nm); the relative orientation is the structure's. COLLIDE_NT = 5: disome footprints are 57-63 nt against 27-29
nt for monosomes (Meydan & Guydosh 2020 Mol Cell 79:588), so 2-5 nt lie between collided footprints. How often it
happens follows from the spacing model alone (independent loading); disome profiling (Zhao et al. 2021: 5.8% of
ribosomes in disomes) is the test (validation/ribosome_cryoet.py). If the pose clashes, the ribosome is placed as
before with a 5-nt linker.
"""
import os
import numpy as np
import numba as nb
import pack3d as p3
import cell3d as c3
import structlib as sl
from build_structs import rot_to_z
from config import LIB
from config import seed as _seed

RIBO = "ribosome_80S"
POLY_FRACTION = 0.7       # share of cytosolic ribosomes in polysomes (log phase, rich medium; assumption)
MEAN_SIZE = 5.0           # mean ribosomes per polysome (assumption), raised when needed so that the polysomes use at most
MRNA_SHARE = 0.85         # this share of the cell's mRNAs (each polysome is one mRNA)
MAX_SIZE = 20
LINK_MIN, LINK_MAX = 3.0, 12.0
CONE = 0.7                # rad, max deviation of the linker from the exit direction
# gene-specific polysomes (mrna_genes): the mRNA between ribosomes is a flexible single strand
SS_LP = 1.5               # nm, persistence length of single-stranded RNA (1-2 nm; MODEL.md 6.5)
SS_NM_PER_NT = 0.6        # nm of contour per nucleotide (as the mRNA threads)
FOOTPRINT_NT = 28         # nt covered by a translating ribosome (ribosome profiling footprint)
RANDOM_SPACING = True     # free nt between neighbouring ribosomes ~ exponential around the chain mean (2026-11-12);
                          # False: every linker of a chain holds the mean (the uniform spacing of builds K-N)
MAX_FAILS_FLEX = 1200      # failed tries before a chain with flexible linkers ends (400 before K_ORIENT)
K_ORIENT = 40             # orientations tried for one drawn linker before a new one is drawn (2026-11-12; was 1, which
                          # replaced every short linker that collided by a new, mostly longer one)
CONE_FLEX = 1.2           # rad: a flexible linker leaves the exit in a wide cone
H_WALK = 20.0             # nm, collision grid cell for the whole-cell walk
NT_STEP = 0.5             # nm between schematic mRNA beads (< 0.635 nm encoder delta limit)
CHANNEL_F = f"{LIB}/{RIBO}_mrna_channel.npy"
COLLIDE_NT = 5.0          # free nt at or below which neighbours are placed as a collided disome (0: off)
DISOME_PDB = "6I7O"
DISOME_F = f"{LIB}/{RIBO}_disome_{DISOME_PDB}.npy"


def channel_path():
    """mRNA channel of the 80S library structure, 5'->3', in the library frame (nm); computed once from 7ZW0."""
    if os.path.exists(CHANNEL_F):
        return np.load(CHANNEL_F)
    import gemmi, extend_complexes as xc, assemblies as A, membrane_orient as mo

    def ca(pid):
        st = gemmi.read_structure(xc.fetch_assembly(pid)); st.setup_entities()
        xyz, kind, acc, pos = mo.beads_with_positions(st, align=A.entity_align(pid))
        return st, xyz, acc, pos

    def kabsch(P, Q):
        pc, qc = P.mean(0), Q.mean(0); U, S, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
        d = np.sign(np.linalg.det(Vt.T @ U.T)); R = Vt.T @ np.diag([1, 1, d]) @ U.T
        return R, qc - R @ pc

    xs, _, _ = sl.beads_from_structure(gemmi.read_structure(xc.fetch_assembly("4V88")))
    _, x0, a0, p0 = ca("4V88")
    key0 = {}
    for i, k in enumerate(zip(a0, p0)):
        if k[0] and k[1] > 0:
            key0.setdefault(k, i)
    st1, x1, a1, p1 = ca("7ZW0")
    seen, pairs = set(), []
    for j, k in enumerate(zip(a1, p1)):
        if k in key0 and k not in seen:
            seen.add(k); pairs.append((key0[k], j))
    I, J = np.array(pairs).T
    P, Q = x1[J], x0[I]
    R, t = kabsch(P, Q)
    for _ in range(4):        # trim flexible parts (40S head, stalks)
        res = np.linalg.norm(P @ R.T + t - Q, axis=1)
        keep = res < max(0.3, np.percentile(res, 75)); R, t = kabsch(P[keep], Q[keep])
    desc = A.entity_desc("7ZW0"); sub2ent = {s: e.name for e in st1.entities for s in e.subchains}
    nts = sorted((r.seqid.num, (r.find_atom("P", "*") or r.find_atom("C4'", "*")).pos) for ch in st1[0] for r in ch
                 if desc.get(sub2ent.get(r.subchain, ""), "").lower() == "mrna" and (r.find_atom("P", "*") or r.find_atom("C4'", "*")))
    m = np.array([[p.x, p.y, p.z] for _, p in nts]) / 10.0 @ R.T + t - xs.mean(0)
    np.save(CHANNEL_F, m.astype(np.float32))
    return m


def disome_pose(stc_walk=None, str_walk=None, log=print):
    """pose of the trailing ribosome of a collided disome in the leading one's library frame: (quaternion w,x,y,z,
    translation nm, push nm), computed once from DISOME_PDB and cached"""
    if os.path.exists(DISOME_F):
        a = np.load(DISOME_F)
        return a[:4], a[4:7], float(a[7])
    import gemmi, extend_complexes as xc, rac, wc3d
    st4 = gemmi.read_structure(xc.fetch_assembly("4V88")); st4.remove_ligands_and_waters()
    xs, _, _ = sl.beads_from_structure(st4); off = xs.mean(0)
    r4 = rac._chains(st4, "P", True)
    st = gemmi.read_structure(xc.fetch_assembly(DISOME_PDB)); st.remove_ligands_and_waters()
    rd = rac._chains(st, "P", True)
    big = [c for c in sorted(rd, key=lambda c: -len(rd[c][0])) if len(rd[c][0]) > 2500][:2]     # the two 25S
    cen = {c: rd[c][1].mean(0) for c in rd}
    poses = []
    for c25 in big:            # each ribosome: the RNA chains closer to its 25S than to the other's
        mem = [c for c in rd if np.argmin([np.linalg.norm(cen[c] - cen[b]) for b in big]) == big.index(c25)]
        P, Q = [], []
        for c in mem:
            s, X = rd[c]
            best = max(((rac._pairs(s4, X4, s, X), k) for k, (s4, X4) in r4.items() if abs(len(s4) - len(s)) < 0.3 * len(s)),
                       key=lambda z: len(z[0][0]), default=None)
            if best:
                P.append(best[0][0]); Q.append(best[0][1])
        Rm, t, med, f = rac._trimmed(np.vstack(P), np.vstack(Q))
        poses.append((Rm, t)); log(f"  {DISOME_PDB} ribosome on {c25}: {sum(map(len, P))} phosphates, median deviation {med:.2f} nm")
    ch = channel_path()
    best = None
    for a_, b_ in ((0, 1), (1, 0)):        # leading, trailing: the trailing entry lies next to the leading exit
        (Ra, ta), (Rb, tb) = poses[a_], poses[b_]
        M = Ra.T @ Rb; t = Ra.T @ (Rb @ off + tb - ta) - off
        d = np.linalg.norm(M @ ch[-1] + t - ch[0])
        if best is None or d < best[0]:
            best = (d, M, t)
    d, M, t = best
    lp = wc3d.build_libpack([RIBO]); trees = [(lp["STC"].astype(float), lp["STR"].astype(float))]
    if stc_walk is not None:
        trees.append((stc_walk, str_walk))
    u = t / np.linalg.norm(t); push = 0.0
    while any((np.linalg.norm(C[:, None] - (C @ M.T + t + push * u)[None], axis=2) - R[:, None] - R[None]).min() < -p3.TOL
              for C, R in trees):
        push += 0.25
    t = t + push * u
    log(f"  collided disome {DISOME_PDB}: trailing entry {d:.2f} nm from the leading exit; pushed {push:.2f} nm to clear the "
        f"sphere trees (centre distance {np.linalg.norm(t):.1f} nm, entry-exit {np.linalg.norm(M @ ch[-1] + t - ch[0]):.2f} nm)")
    q = _mat2quat(M)
    np.save(DISOME_F, np.r_[q, t, push])
    return q, t, push


def _mat2quat(R):
    from scipy.spatial.transform import Rotation
    x, y, z, w = Rotation.from_matrix(R).as_quat()
    return np.array([w, x, y, z])          # pack3d convention (w, x, y, z)


@nb.njit(cache=True)
def _qmul(a, b, out):
    out[0] = a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3]
    out[1] = a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2]
    out[2] = a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1]
    out[3] = a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0]


@nb.njit(cache=True)
def walk(seed, glo, ghi, G, pores, pidx, PT, PQ, PC, stc, str_, br, mask, sizes, ex, tex, en, q_en, OC, OR, cap_s, cap_e,
         lrms, cone, max_fails, seeds, tmem, tmax, q_rel, t_rel, collide_nt):
    """Grow one chain per entry of sizes. ex/tex: exit point and outward exit direction (library frame);
    en: entry point; q_en: quaternion taking the outward entry direction to +z. lrms[c] > 0: exit-to-entry distance
    of chain c drawn from the end-to-end distribution of a Gaussian chain with that rms (the free mRNA between two
    ribosomes); lrms[c] == 0: uniform LINK_MIN..LINK_MAX. tmem[c] >= 0: chain c is tethered to membrane tmem[c]
    (cell3d M_* index): its first ribosome starts near seeds[c] and every ribosome centre stays within tmax of that
    membrane (after 300 failed starts the chain is placed untethered). q_rel, t_rel: pose of a collided trailing
    ribosome in the leading one's frame, used when the drawn free nucleotides are <= collide_nt (0: never).
    Returns pos, quat, chain index, collided flag (placed in the collided pose behind the previous ribosome)."""
    np.random.seed(seed)
    dims, head, nxt, sid, cen, rad, cnt = p3.grid_new(glo, ghi, H_WALK, cap_s, cap_e)
    if OC.shape[0]:
        p3.g_insert(glo, H_WALK, dims, head, nxt, sid, cen, rad, cnt, OC, OR)     # e.g. actin cables
    ntot = sizes.sum()
    pos = np.empty((ntot, 3), np.float32); qu = np.empty((ntot, 4), np.float32); chain = np.empty(ntot, np.int32)
    coll = np.zeros(ntot, np.bool_); qprev = np.empty(4); tprev = np.empty(3); Mp = np.empty((3, 3))
    out = np.zeros(c3.NMEM); M = np.empty((3, 3)); q = np.empty(4); q2 = np.empty(4); t = np.empty(3)
    P = np.empty((stc.shape[0], 3)); D = np.empty(3); negD = np.empty(3)
    exw = np.empty(3); texw = np.empty(3)
    n = 0
    R0 = G[2] - 30.0
    for c in range(sizes.shape[0]):
        placed = 0; fails = 0; tm = tmem[c]; cur_L = -1.0; tries_L = 0; x_pend = -1.0; coll_tried = False
        while placed < sizes[c]:
            coll_now = False
            if placed == 0 and tm >= 0:
                if fails > 300:
                    tm = -1
                for k in range(3):
                    t[k] = seeds[c, k] + (np.random.random() * 2 - 1) * 15.0
                p3.rand_quat(q)
            elif placed == 0:
                while True:
                    for k in range(3):
                        t[k] = (np.random.random() * 2 - 1) * R0
                    if t[0] * t[0] + t[1] * t[1] + t[2] * t[2] < R0 * R0:
                        break
                p3.rand_quat(q)
            elif collide_nt > 0.0 and lrms[c] > 0.0 and RANDOM_SPACING and cur_L < 0.0 and not coll_tried:
                # draw this link's free nucleotides; a collision places the trailing ribosome in the disome pose
                nt = lrms[c] ** 2 / (2.0 * SS_LP * SS_NM_PER_NT)
                x_pend = -np.log(1.0 - np.random.random()) * nt
                coll_tried = True
                if x_pend <= collide_nt:
                    coll_now = True
                    _qmul(qprev, q_rel, q); p3.quat_mat(qprev, Mp)
                    for k in range(3):
                        t[k] = tprev[k] + Mp[k, 0] * t_rel[0] + Mp[k, 1] * t_rel[1] + Mp[k, 2] * t_rel[2]
            if placed > 0 and not coll_now:
                while True:   # direction within CONE of the previous exit direction
                    for k in range(3):
                        D[k] = np.random.normal()
                    dn = np.sqrt(D[0] ** 2 + D[1] ** 2 + D[2] ** 2)
                    s = 0.0
                    for k in range(3):
                        D[k] = D[k] / dn
                        s += D[k] * texw[k]
                    if s > np.cos(cone):
                        break
                if cur_L > 0.0 and tries_L < K_ORIENT:
                    L = cur_L                    # the linker is set by the mRNA: try other orientations first
                elif lrms[c] > 0.0:
                    rms = lrms[c]
                    if RANDOM_SPACING:
                        # ribosomes load independently, so the free nucleotides between neighbours are exponential
                        # around the chain's mean (lrms^2 is proportional to the free nucleotides)
                        nt = lrms[c] ** 2 / (2.0 * SS_LP * SS_NM_PER_NT)
                        xd = x_pend if x_pend >= 0.0 else -np.log(1.0 - np.random.random()) * nt
                        x_pend = -1.0
                        rms = np.sqrt(2.0 * SS_LP * SS_NM_PER_NT * max(xd, 5.0))
                    sg = rms / np.sqrt(3.0)
                    a1 = np.random.normal() * sg; a2 = np.random.normal() * sg; a3 = np.random.normal() * sg
                    L = min(max(np.sqrt(a1 * a1 + a2 * a2 + a3 * a3), 2.0), 4.0 * rms)
                else:
                    L = LINK_MIN + np.random.random() * (LINK_MAX - LINK_MIN)
                if L != cur_L:
                    cur_L = L; tries_L = 0
                tries_L += 1
                for k in range(3):
                    negD[k] = -D[k]
                p3.quat_z_to(negD, np.random.random() * 2 * np.pi, q2)   # outward entry direction -> -D
                _qmul(q2, q_en, q)
                p3.quat_mat(q, M)
                for k in range(3):
                    t[k] = exw[k] + L * D[k] - (M[k, 0] * en[0] + M[k, 1] * en[1] + M[k, 2] * en[2])
            p3.quat_mat(q, M); p3.transform(stc, M, t, P)
            good = not p3.g_overlap(glo, H_WALK, dims, head, nxt, sid, cen, rad, P, str_, p3.TOL)
            if good:
                good = p3.check_volume(P, str_, t[0], t[1], t[2], br, mask, G, pores, pidx, PT, PQ, PC, out)
            if good and tm >= 0:
                c3.mem_sdfs(t[0], t[1], t[2], G, pores, pidx, PT, PQ, PC, out)
                good = abs(out[tm]) <= tmax
            if not good:
                fails += 1
                if placed > 0 and fails > max_fails:
                    break            # polysome ends shorter
                continue
            if not p3.g_insert(glo, H_WALK, dims, head, nxt, sid, cen, rad, cnt, P, str_):
                return pos[:n], qu[:n], chain[:n], coll[:n]
            for k in range(3):
                exw[k] = t[k] + M[k, 0] * ex[0] + M[k, 1] * ex[1] + M[k, 2] * ex[2]
                texw[k] = M[k, 0] * tex[0] + M[k, 1] * tex[1] + M[k, 2] * tex[2]
                pos[n, k] = t[k]
            for k in range(4):
                qu[n, k] = q[k]; qprev[k] = q[k]
            for k in range(3):
                tprev[k] = t[k]
            chain[n] = c; coll[n] = coll_now
            cur_L = -1.0; x_pend = -1.0; coll_tried = False
            n += 1; placed += 1; fails = 0
    return pos[:n], qu[:n], chain[:n], coll[:n]


def _resample(P, step=NT_STEP):
    seg = np.linalg.norm(np.diff(P, axis=0), axis=1)
    s = np.r_[0, np.cumsum(seg)]
    u = np.linspace(0, s[-1], max(2, int(np.ceil(s[-1] / step)) + 1))
    return np.c_[[np.interp(u, s, P[:, k]) for k in range(3)]].T


def _hermite(a, ta, b, tb):
    L = np.linalg.norm(b - a)
    s = np.linspace(0, 1, max(3, int(L * 2)))[:, None]
    h00 = 2 * s ** 3 - 3 * s ** 2 + 1; h10 = s ** 3 - 2 * s ** 2 + s; h01 = -2 * s ** 3 + 3 * s ** 2; h11 = s ** 3 - s ** 2
    return h00 * a + h10 * ta * L + h01 * b + h11 * tb * L


def place(cell, lp, n_ribo, seed=17, obstacles=None, log=print, n_mrna=None, sizes=None, linker_nt=None, tether=None):
    """Place polysomes holding ~n_ribo cytosolic ribosomes. Returns dict with ribosome instances, per-ribosome mRNA
    polylines (channel + linker to the next ribosome) and the linker beads as packing obstacles. sizes: ribosomes
    per chain, one chain per translated mRNA (mrna_genes.transcripts); None: 1 + Poisson(mean - 1) up to n_ribo.
    linker_nt: free mRNA nucleotides between neighbouring ribosomes per chain; the exit-to-entry distance is then the
    end-to-end distance of a single-stranded RNA chain of that contour (persistence SS_LP; Gaussian chain).
    tether: (membrane index per chain, -1 = free; seed point per chain; max ribosome-membrane distance), see walk."""
    from scipy.spatial.transform import Rotation
    li = int(np.where(lp["ids"] == RIBO)[0][0])
    z = np.load(f"{LIB}/{RIBO}.npz")
    stc, str_ = sl.sphere_tree(z["beads"].astype(np.float64), per=550)       # ~30 coarse spheres for the walk
    ch = channel_path().astype(np.float64)
    ex, en = ch[0], ch[-1]
    tex = ch[0] - ch[3]; tex /= np.linalg.norm(tex)
    ten = ch[-1] - ch[-4]; ten /= np.linalg.norm(ten)
    q_en = _mat2quat(rot_to_z(ten))
    seed = _seed(seed); rng = np.random.default_rng(seed)
    if sizes is not None:
        sizes = np.asarray(sizes, np.int64); mean = float(sizes.mean())
    else:
        mean = MEAN_SIZE if not n_mrna else max(MEAN_SIZE, n_ribo / (MRNA_SHARE * n_mrna))
        if mean > MEAN_SIZE:
            log(f"polysomes: mean size raised to {mean:.1f} ribosomes so that {n_ribo:.0f} ribosomes fit on {MRNA_SHARE:.0%} of {n_mrna} mRNAs")
        sizes = []
        while sum(sizes) < n_ribo:
            sizes.append(int(min(MAX_SIZE, 1 + rng.poisson(mean - 1))))
        sizes = np.array(sizes, np.int64)
    G = cell.G
    glo = np.full(3, -G[0] - 10.0); ghi = -glo
    a = cell.args(glo, ghi, margin=50.0)
    OC, OR = obstacles if obstacles is not None else (np.zeros((0, 3)), np.zeros(0))
    lrms = np.zeros(len(sizes)) if linker_nt is None else np.sqrt(2.0 * SS_LP * np.maximum(np.asarray(linker_nt, float), 5.0) * SS_NM_PER_NT)
    ns = len(str_) * len(sizes) * int(mean + 2) + len(OR)
    if tether is None:
        tether = (np.full(len(sizes), -1, np.int64), np.zeros((len(sizes), 3)), 0.0)
    if COLLIDE_NT > 0 and linker_nt is not None:
        q_rel, t_rel, _ = disome_pose(stc, str_, log=log)
    else:
        q_rel, t_rel = np.array([1.0, 0, 0, 0]), np.zeros(3)
    pos, qu, chn, coll = walk(seed, glo, ghi, *a, stc, str_, float(np.max(np.linalg.norm(stc, axis=1) + str_)), 1 << c3.CYTO,
                              sizes, ex, tex, en, q_en, OC.astype(np.float64), OR.astype(np.float64), ns + 1000, ns * 14 + 10000,
                              lrms, CONE_FLEX if linker_nt is not None else CONE, MAX_FAILS_FLEX if linker_nt is not None else 60,
                              np.asarray(tether[1], np.float64), np.asarray(tether[0], np.int64), float(tether[2]),
                              np.asarray(q_rel, np.float64), np.asarray(t_rel, np.float64),
                              COLLIDE_NT if linker_nt is not None else 0.0)
    # mRNA polylines
    Rm = Rotation.from_quat(np.c_[qu[:, 1:], qu[:, :1]].astype(np.float64)).as_matrix()
    pieces, links = [], []
    for i in range(len(pos)):
        P = ch @ Rm[i].T + pos[i]
        if i + 1 < len(pos) and chn[i + 1] == chn[i]:
            Pn = ch @ Rm[i + 1].T + pos[i + 1]
            lk = _hermite(P[0], Rm[i] @ tex, Pn[-1], -(Rm[i + 1] @ ten))
            links.append(lk[1:-1]); P = np.vstack([P[::-1], lk[1:]])
        else:
            P = P[::-1]
        pieces.append(_resample(P))
    off = np.r_[0, np.cumsum([len(p) for p in pieces])].astype(np.int64)
    n_chain = len(np.unique(chn))
    log(f"polysomes: {n_chain} chains, {len(pos)}/{sizes.sum()} ribosomes (mean {len(pos) / max(n_chain, 1):.1f} per mRNA), "
        f"{off[-1]} mRNA beads; {int(coll.sum())} ribosomes placed as collided behind their neighbour")
    return dict(poly_pos=pos, poly_q=qu, poly_chain=chn, poly_coll=coll, li_ribo=li, poly_mrna=np.concatenate(pieces).astype(np.float32),
                poly_mrna_off=off, poly_link=np.concatenate(links).astype(np.float32) if links else np.zeros((0, 3), np.float32),
                poly_chains=n_chain)
