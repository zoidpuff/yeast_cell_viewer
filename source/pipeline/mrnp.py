"""Proteins on the ends of the mRNA threads (IMPROVEMENTS M6): Pab1 on the poly(A) tail, eIF4E (Cdc33) on the cap.

Pab1 binds poly(A) with a footprint of about 27 nt (unchecked), so each 3' end carries max(1, poly(A) // 27) Pab1, at
most PAB_MAX (per-gene PAL-seq tails, mrna_genes.py; median 33 nt, so mostly one). Each mRNA's 5' end carries one
eIF4E. Pab1 and Cdc33 are the AlphaFold monomers of the library, placed at random orientations so that one of their
spheres touches the end bead (within TOUCH_NM) and nothing else overlaps. Ends: the 3' tail (mrna_kind 1) and the
5' tail (kind 2) of each polysome mRNA end at its 3' and 5' end; free mRNAs (kind 0) run 5' to 3'.
Copies come out of the free Pab1 and Cdc33 pools (142,600 and 41,000).
"""
import numpy as np
from config import seed as _seed
import pack3d as p3
import cell3d as c3

PAB_FOOTPRINT_NT, PAB_MAX = 27, 3
TOUCH_NM, TRIES = 1.5, 60
OWN_NT = 40                         # beads of the thread's own end that the protein may overlap (it binds them)


def ends(pre):
    """(3' end points, 3' genes, 5' end points) of all threads"""
    P, off, kind = pre["mrna_pts"], pre["mrna_off"], pre["mrna_kind"]
    g = pre.get("mrna_gene")
    e3, g3, e5, own3, own5 = [], [], [], [], []
    for i, k in enumerate(kind):
        a, b = off[i], off[i + 1]
        if b - a < 2:
            continue
        if k in (0, 1):
            e3.append(b - 1); g3.append(g[i] if g is not None else ""); own3.append((max(a, b - OWN_NT), b))
        if k == 0:
            e5.append(a); own5.append((a, min(b, a + OWN_NT)))
        if k == 2:
            e5.append(b - 1); own5.append((max(a, b - OWN_NT), b))
    return np.array(e3, np.int64), np.array(g3), np.array(e5, np.int64), own3, own5


def place(cell, lp, pre, n_pab, n_cap, obstacles, polya=None, seed=37, log=print):
    """obstacles: (centres, radii) of everything placed so far except the mRNA beads (added here).
    polya: {gene: poly(A) nt}. Returns pre-pass group 'rnp' (rnp_type, rnp_pos, rnp_q)."""
    from scipy.spatial import cKDTree
    from scipy.spatial.transform import Rotation
    import mrna
    ids = list(lp["ids"])
    if "mrna_pts" not in pre or "PAB1" not in ids:
        return {}
    seed = _seed(seed); rng = np.random.default_rng(seed)
    P = pre["mrna_pts"].astype(np.float64)
    e3, g3, e5, own3, own5 = ends(pre)
    OC, OR = obstacles
    AC = np.concatenate([OC, P]); AR = np.concatenate([OR, np.full(len(P), mrna.BEAD_R)])
    is_bead = np.r_[np.zeros(len(OC), bool), np.ones(len(P), bool)]
    bead_idx = np.r_[np.full(len(OC), -1), np.arange(len(P))]
    tree = cKDTree(AC); rmax = AR.max()
    G = cell.G; glo = np.full(3, -G[0] - 10.0); a = cell.args(glo, -glo, margin=50.0)
    out = np.zeros(c3.NMEM)
    placed_c, placed_r = [], []          # spheres of the proteins placed here
    st = dict(tree=None, n=0, R=None)      # tree over placed_c[:n], rebuilt every 200 placements

    def free(S, R, own):
        for nb, s, r in zip(tree.query_ball_point(S, R.max() + rmax), S, R):
            if not nb:
                continue
            nb = np.array(nb)
            nb = nb[~(is_bead[nb] & (bead_idx[nb] >= own[0]) & (bead_idx[nb] < own[1]))]
            if len(nb) and (np.linalg.norm(AC[nb] - s, axis=1) < AR[nb] + r - p3.TOL).any():
                return False
        if len(placed_c) - st["n"] > 200:      # was 2,000: the linear scan below over the placements not yet in
                                               # the tree took most of this step's 6 minutes; the test is the same
            st.update(tree=cKDTree(np.concatenate(placed_c)), n=len(placed_c), R=np.concatenate(placed_r))
        if st["tree"] is not None:
            for nb, s, r in zip(st["tree"].query_ball_point(S, R.max() + st["R"].max()), S, R):
                if nb and (np.linalg.norm(st["tree"].data[nb] - s, axis=1) < st["R"][nb] + r - p3.TOL).any():
                    return False
        for C_, R_ in zip(placed_c[st["n"]:], placed_r[st["n"]:]):
            if np.linalg.norm(C_.mean(0) - S.mean(0)) < 60 and (np.linalg.norm(S[:, None] - C_[None], axis=2) < R[:, None] + R_[None] - p3.TOL).any():
                return False
        return True

    types, pos, quat = [], [], []

    def put(li, E, own):
        C = lp["STC"][lp["st_off"][li]:lp["st_off"][li + 1]].astype(np.float64)
        R = lp["STR"][lp["st_off"][li]:lp["st_off"][li + 1]].astype(np.float64)
        for _ in range(TRIES):
            Rm = Rotation.random(random_state=rng).as_matrix()
            S0 = C @ Rm.T
            k = rng.integers(len(S0))                   # this sphere touches the end bead
            u = rng.normal(size=3); u /= np.linalg.norm(u)
            t = E - S0[k] + u * (R[k] + mrna.BEAD_R + rng.uniform(0, TOUCH_NM))
            S = S0 + t
            if not free(S, R, own):
                continue
            if not p3.check_volume(S, R, t[0], t[1], t[2], float(lp["BR"][li]), 1 << c3.CYTO, *a, out):
                continue
            placed_c.append(S); placed_r.append(R)
            types.append(li); pos.append(t); quat.append(np.roll(Rotation.from_matrix(Rm).as_quat(), 1))
            return True
        return False

    li_pab = ids.index("PAB1"); li_cap = ids.index("CDC33") if "CDC33" in ids else -1
    n3 = np.ones(len(e3), np.int64)
    if polya is not None:
        n3 = np.clip(np.array([polya.get(g, PAB_FOOTPRINT_NT) for g in g3]) // PAB_FOOTPRINT_NT, 1, PAB_MAX).astype(np.int64)
    got_p = 0
    for i in rng.permutation(len(e3)):
        for _ in range(n3[i]):
            if got_p >= n_pab:
                break
            if put(li_pab, P[e3[i]], own3[i]):
                got_p += 1
    got_c = 0
    if li_cap >= 0:
        for i in rng.permutation(len(e5)):
            if got_c >= n_cap:
                break
            if put(li_cap, P[e5[i]], own5[i]):
                got_c += 1
    log(f"mRNP ends: {got_p} Pab1 on {len(e3)} 3' ends ({int(n3.sum())} wanted), {got_c} eIF4E (Cdc33) on {len(e5)} 5' ends")
    if not types:
        return {}
    return dict(rnp_type=np.array(types, np.int32), rnp_pos=np.array(pos, np.float32), rnp_q=np.array(quat, np.float32),
                n_pab=got_p, n_cap=got_c)
