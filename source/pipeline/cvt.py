"""Cvt complex (literature item L18): the cytoplasm-to-vacuole targeting cargo of growing cells.

Precursor aminopeptidase I (prApe1) dodecamers assemble through their propeptides into one large complex per cell,
together with alpha-mannosidase Ams1, coated by the receptor Atg19 (which links it to Atg11 and Atg8). It sits at the
phagophore assembly site next to the vacuole and is delivered in a Cvt vesicle of ~150 nm (Baba et al. 1997 J Cell
Biol 139:1687; review: Lynch-Day & Klionsky 2010 FEBS Lett 584:1359).

Built here: one cluster in the cytosol, its edge CVT_GAP from the vacuole membrane, holding CVT_FRACTION of the Ape1
dodecamers (the library's 12-mer, 5JGF; the propeptides that glue dodecamers together are not resolved) and of Ams1,
packed at random to CVT_PACKING volume fraction, with Atg19 on the surface. The rest of Ape1 and Ams1 stays in the
vacuole as mature enzyme.
"""
import numpy as np
from config import seed as _seed
import pack3d as p3
import cell3d as c3

APE1, AMS1, ATG19 = "APE1_oligo", "AMS1_oligo", "ATG19"     # library ids (Ams1: AF dimer)
CVT_FRACTION = 0.5       # share of Ape1 (and Ams1) still precursor in the Cvt complex (assumption)
CVT_PACKING = 0.35       # protein volume fraction inside the complex (assumption; the complex is dense in EM)
CVT_GAP = 40.0           # nm between the complex and the vacuole membrane
H = 20.0


def place(cell, lp, n_ape1, n_ams1, n_atg19, obstacles=None, seed=61, log=print):
    """one Cvt complex: at most n_ape1 dodecamers, n_ams1 Ams1, n_atg19 Atg19 (instances)"""
    from scipy.spatial.transform import Rotation
    seed = _seed(seed); rng = np.random.default_rng(seed)
    ids = list(lp["ids"])
    li = {k: ids.index(k) for k in (APE1, AMS1, ATG19) if k in ids}
    sub = {k: (lp["STC"][lp["st_off"][v]:lp["st_off"][v + 1]], lp["STR"][lp["st_off"][v]:lp["st_off"][v + 1]], float(lp["BR"][v]))
           for k, v in li.items()}
    vol = {k: float((4 / 3 * np.pi * s[1] ** 3).sum()) for k, s in sub.items()}
    n_a = int(n_ape1 * CVT_FRACTION); n_m = int(n_ams1 * CVT_FRACTION) if AMS1 in li else 0
    Rb = ((n_a * vol[APE1] + n_m * vol.get(AMS1, 0)) / CVT_PACKING * 3 / (4 * np.pi)) ** (1 / 3)
    # site: next to a vacuole lobe, in the cytosol
    G = cell.G; vac = [i for i in range(len(cell.PT)) if cell.PC[i] == c3.C_VAC]
    out = np.zeros(c3.NMEM)
    centre = None
    vd = lambda x: (c3.mem_sdfs(x[0], x[1], x[2], *cell.args(x - 50, x + 50), out), out[c3.M_VAC])[1]
    dirs = rng.normal(size=(64, 3)); dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
    for _ in range(20000):
        p = rng.normal(size=3); p *= rng.uniform(0, G[0]) / np.linalg.norm(p)
        d = vd(p)
        if not (0.0 < d < 20.0):                       # just outside a vacuole lobe (within the field's band)
            continue
        n = np.array([vd(p + e) - vd(p - e) for e in np.eye(3)]); n /= np.linalg.norm(n)
        c = p + n * (Rb + CVT_GAP - d)
        shell = c + dirs * (Rb + 10.0)
        if cell.label(c)[0] == c3.CYTO and (cell.label(shell) == c3.CYTO).all():
            centre = c; break
    if centre is None:
        log("Cvt complex: no site next to the vacuole"); return {}, {APE1: 0, AMS1: 0, ATG19: 0}
    glo, ghi = centre - Rb - 60, centre + Rb + 60
    a = cell.args(glo, ghi, margin=50.0)
    cap = 200000 + (len(obstacles[1]) if obstacles is not None else 0)
    G1 = p3.grid_new(glo, ghi, H, cap, cap * 14)
    if obstacles is not None and len(obstacles[1]):
        k = np.all((obstacles[0] > glo) & (obstacles[0] < ghi), axis=1)
        p3.g_insert(glo, H, *G1, obstacles[0][k], obstacles[1][k])
    M = np.empty((3, 3)); q = np.empty(4)
    pos, rot, typ = [], [], []
    used = {APE1: 0, AMS1: 0, ATG19: 0}

    def put(lid, t):
        stc, str_, br = sub[lid]
        p3.rand_quat(q); p3.quat_mat(q, M); P = np.empty((len(stc), 3)); p3.transform(stc, M, t, P)
        if p3.g_overlap(glo, H, *G1[:6], P, str_, p3.TOL):
            return False
        if not p3.check_volume(P, str_, t[0], t[1], t[2], br, 1 << c3.CYTO, *a, out):
            return False
        p3.g_insert(glo, H, *G1, P, str_)
        pos.append(t.copy()); rot.append(q.copy()); typ.append(li[lid]); used[lid] += 1
        return True

    order = [APE1] * n_a + [AMS1] * n_m
    rng.shuffle(order)
    for lid in order:                                  # cargo: random sequential addition inside the ball
        r_in = Rb - sub[lid][2] * 0.5
        for _ in range(300):
            v = rng.normal(size=3); v *= r_in * rng.random() ** (1 / 3) / np.linalg.norm(v)
            if put(lid, centre + v):
                break
    if ATG19 in li:                                    # receptor coat on the surface
        for _ in range(n_atg19):
            for _ in range(200):
                v = rng.normal(size=3); v /= np.linalg.norm(v)
                if put(ATG19, centre + v * (Rb + sub[ATG19][2] * 0.6)):
                    break
    log(f"Cvt complex: radius {Rb:.0f} nm, {CVT_GAP:.0f} nm from the vacuole; {used[APE1]} of {n_a} Ape1 dodecamers, "
        f"{used[AMS1]} of {n_m} Ams1, {used[ATG19]} Atg19 on the surface")
    return dict(cvt_pos=np.array(pos, np.float32).reshape(-1, 3), cvt_q=np.array(rot, np.float32).reshape(-1, 4),
                cvt_type=np.array(typ, np.int32), cvt_centre=np.asarray(centre, np.float32), cvt_radius=np.float32(Rb)), used
