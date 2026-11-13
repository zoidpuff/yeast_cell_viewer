"""Kinetochores on the 16 kinetochore-microtubule plus ends (IMPROVEMENTS K6).

Budding yeast has point centromeres: one microtubule per chromosome. Each kinetochore here is one CCAN bound to its
centromeric DNA with Cbf1 (CPX-2533, 8OVW), centred CCAN_GAP nm beyond the plus end along the microtubule axis, and
N_MIND MIND (Mis12) complexes (CPX-1186, 9S5N) in a ring around the plus end, between the microtubule and the CCAN.
N_MIND = 6: quantitative imaging counts one or two copies of inner proteins (Mif2) up to 16 DAM/DASH per kinetochore
(Joglekar et al. 2006 Nat Cell Biol 8:581); the MIND number of about 6 is from the paper's table (unchecked). Not built:
the Cse4 nucleosome. The chromosome's two arms then start beyond the CCAN. The rest of the CCAN and MIND pools
stay free in the nucleoplasm.
Outer layer (IMPROVEMENTS N4; `outer`): per kinetochore N_DASH copies of each DASH/Dam1 subunit (16, the count of
Joglekar et al. 2006) in a ring around the kinetochore microtubule over its last DASH_SPAN nm, and N_NDC80 copies of
each Ndc80-complex subunit (about 8 in the same study; unchecked) in a wider band reaching from the lattice to just
beyond the tip. Both complexes are placed as their subunits (AlphaFold monomers; no assembled library types), so
the ring geometry is approximate.
"""
import numpy as np
from config import seed as _seed
import pack3d as p3
import cell3d as c3

CCAN_ID, MIND_ID = "CPX_2533_Kinetochore_CCAN_complex", "CPX_1186_Kinetochore_MIS12_complex"
N_MIND, CCAN_GAP, TRIES = 6, 4.0, 200
DASH = ["DAM1", "DUO1", "ASK1", "DAD1", "DAD2", "DAD3", "DAD4", "SPC19", "SPC34", "HSK3"]
NDC80C = ["NDC80", "NUF2", "SPC24", "SPC25"]
N_DASH, N_NDC80 = 16, 8
DASH_SPAN, DASH_R = (5.0, 60.0), (16.0, 32.0)          # nm behind the tip; nm from the microtubule axis
NDC_SPAN, NDC_R = (-15.0, 60.0), (14.0, 50.0)


def place(cell, lp, tips, dirs, obstacles, n_ccan, n_mind, seed=61, log=print):
    from scipy.spatial import cKDTree
    from scipy.spatial.transform import Rotation
    ids = list(lp["ids"])
    if CCAN_ID not in ids:
        return {}
    seed = _seed(seed); rng = np.random.default_rng(seed)
    OC, OR = obstacles
    tree = cKDTree(OC) if len(OC) else None; rmax = OR.max() if len(OR) else 0.0
    cn = cell.G[3:6]; rn = float(cell.G[11:14].max()); glo, ghi = cn - rn - 80, cn + rn + 80
    a = cell.args(glo, ghi, margin=200.0)
    out = np.zeros(c3.NMEM)
    placed = []

    def lib(li):
        return (lp["STC"][lp["st_off"][li]:lp["st_off"][li + 1]].astype(np.float64),
                lp["STR"][lp["st_off"][li]:lp["st_off"][li + 1]].astype(np.float64))

    def free(S, R, li, t):
        if tree is not None:
            for nb, s, r in zip(tree.query_ball_point(S, R.max() + rmax), S, R):
                if nb and (np.linalg.norm(OC[nb] - s, axis=1) < OR[nb] + r - p3.TOL).any():
                    return False
        for Q, QR in placed:
            if np.linalg.norm(Q.mean(0) - S.mean(0)) < 60 and (np.linalg.norm(S[:, None] - Q[None], axis=2) < R[:, None] + QR[None] - p3.TOL).any():
                return False
        return p3.check_volume(S, R, t[0], t[1], t[2], float(lp["BR"][li]), 1 << c3.NUCP, *a, out)

    def put(li, centre_fn):
        C, R = lib(li)
        for k in range(TRIES):
            Rm = Rotation.random(random_state=rng).as_matrix()
            t = centre_fn(k) - Rm @ C.mean(0)
            S = C @ Rm.T + t
            if free(S, R, li, t):
                placed.append((S, R))
                return t, np.roll(Rotation.from_matrix(Rm).as_quat(), 1), S
        return None

    li_c = ids.index(CCAN_ID); li_m = ids.index(MIND_ID) if MIND_ID in ids else -1
    rc = float(np.max(np.linalg.norm(lib(li_c)[0] - lib(li_c)[0].mean(0), axis=1)))
    types, pos, quat, starts = [], [], [], []
    nc = nm = 0
    for T, u in zip(tips, dirs):
        u = u / np.linalg.norm(u)
        e1 = np.cross(u, [0.3, 0.5, 0.8]); e1 /= np.linalg.norm(e1); e2 = np.cross(u, e1)
        start = T + u * 20.0
        if nc < n_ccan:
            r = put(li_c, lambda k: T + u * (CCAN_GAP + 0.6 * rc + 0.05 * k) + rng.normal(0, 1 + 0.05 * k, 3))
            if r is not None:
                types.append(li_c); pos.append(r[0]); quat.append(r[1]); nc += 1
                start = r[2][np.argmax(r[2] @ u)] + u * 12.0            # chromatin begins beyond the CCAN
        for j in range(N_MIND if li_m >= 0 else 0):
            if nm >= n_mind:
                break
            ph = 2 * np.pi * j / N_MIND
            r = put(li_m, lambda k: T + (12.0 + 0.05 * k) * (np.cos(ph) * e1 + np.sin(ph) * e2) + rng.normal(0, 1.5, 3))
            if r is not None:
                types.append(li_m); pos.append(r[0]); quat.append(r[1]); nm += 1
        starts.append(start)
    log(f"kinetochores: {nc} CCAN and {nm} MIND complexes on {len(tips)} kinetochore-microtubule plus ends")
    C_all = np.concatenate([S for S, _ in placed]) if placed else np.zeros((0, 3))
    R_all = np.concatenate([R for _, R in placed]) if placed else np.zeros(0)
    return dict(kt_type=np.array(types, np.int32), kt_pos=np.array(pos, np.float32).reshape(-1, 3),
                kt_q=np.array(quat, np.float32).reshape(-1, 4), starts=np.array(starts), spheres=(C_all, R_all), n_ccan=nc, n_mind=nm)


def outer(cell, lp, tips, dirs, obstacles, pools, seed=67, log=print):
    """DASH and Ndc80-complex subunits around the kinetochore-microtubule plus ends. pools: {gene: copies available}.
    Returns pre-pass group 'kto' (kto_type, kto_pos, kto_q), per-gene counts and the placed spheres."""
    from scipy.spatial import cKDTree
    from scipy.spatial.transform import Rotation
    ids = list(lp["ids"])
    seed = _seed(seed); rng = np.random.default_rng(seed)
    OC, OR = obstacles
    cn = cell.G[3:6]; rn = float(cell.G[11:14].max()); glo, ghi = cn - rn - 80, cn + rn + 80
    a = cell.args(glo, ghi, margin=200.0)
    out = np.zeros(c3.NMEM)
    placed_c, placed_r, types, pos, quat = [], [], [], [], []
    counts = {}

    grid = p3.LocalGrid(glo, ghi, OC.astype(np.float64), OR.astype(np.float64))

    plan = [(g, N_DASH, DASH_SPAN, DASH_R) for g in DASH] + [(g, N_NDC80, NDC_SPAN, NDC_R) for g in NDC80C]
    for T, u in zip(tips, dirs):
        u = u / np.linalg.norm(u)
        e1 = np.cross(u, [0.3, 0.5, 0.8]); e1 /= np.linalg.norm(e1); e2 = np.cross(u, e1)
        for g, n, span, rad in plan:
            if g not in ids:
                continue
            li = ids.index(g)
            C = lp["STC"][lp["st_off"][li]:lp["st_off"][li + 1]].astype(np.float64)
            R = lp["STR"][lp["st_off"][li]:lp["st_off"][li + 1]].astype(np.float64)
            fails = 0
            for _ in range(n):
                if counts.get(g, 0) >= pools.get(g, 0) or fails >= 3:
                    break
                fails += 1
                for k in range(TRIES):
                    ph = rng.uniform(0, 2 * np.pi); r = rng.uniform(*rad); z = rng.uniform(*span)
                    ctr = T - z * u + r * (np.cos(ph) * e1 + np.sin(ph) * e2)
                    Rm = Rotation.random(random_state=rng).as_matrix()
                    t = ctr - Rm @ C.mean(0); S = C @ Rm.T + t
                    if grid.free(S, R) and p3.check_volume(S, R, t[0], t[1], t[2], float(lp["BR"][li]), 1 << c3.NUCP, *a, out):
                        grid.add(S, R); placed_c.append(S); placed_r.append(R)
                        types.append(li); pos.append(t); quat.append(np.roll(Rotation.from_matrix(Rm).as_quat(), 1))
                        counts[g] = counts.get(g, 0) + 1; fails = 0
                        break
    want = sum(n for _, n, _, _ in plan) * len(tips)
    log(f"outer kinetochores: {len(types)} of {want} DASH and Ndc80-complex subunits on {len(tips)} plus ends")
    if not types:
        return {}, counts, (np.zeros((0, 3)), np.zeros(0))
    return (dict(kto_type=np.array(types, np.int32), kto_pos=np.array(pos, np.float32), kto_q=np.array(quat, np.float32)),
            counts, (np.concatenate(placed_c), np.concatenate(placed_r)))
