"""Nucleus-vacuole junction (literature item L19, the part this geometry supports).

Where a vacuole lobe lies against the nucleus, the outer nuclear membrane protein Nvj1 binds Vac8 on the vacuolar
membrane (Pan et al. 2000 Mol Biol Cell 11:2445); the structure of Vac8 with the Vac8-binding region of Nvj1 is 5XJG
(Jeong et al. 2017 PNAS 114:E4539). Built here: Vac8-Nvj1 complexes on the cytoplasmic face of the vacuolar membrane
wherever it is within CONTACT nm of the outer nuclear membrane, Vac8 against the membrane (its palmitoylated N-terminus
is not resolved) and the Nvj1 fragment pointing towards the nucleus. Nvj1's transmembrane anchor and linker are not
resolved and not drawn. Copies: all Nvj1 (each with one Vac8).
"""
import os, json
import numpy as np
import gemmi
import extend_complexes as xc
import structlib as sl
import build_structs as bs
import pack3d as p3
import cell3d as c3
from config import LIB
from config import seed as _seed

LIB_ID = "Vac8_Nvj1"
PDB = "5XJG"
CONTACT = 45.0        # nm, vacuolar membrane to outer nuclear membrane (midsurfaces)
GAP = 1.0             # nm between the complex and the vacuolar leaflet
H = 20.0
INFO_F = f"{LIB}/{LIB_ID}.json"


def model(log=print):
    """library entry in a frame with +z from the Vac8 N-terminus to the Nvj1 N-terminus (towards the nucleus)"""
    if os.path.exists(INFO_F) and os.path.exists(f"{LIB}/{LIB_ID}.npz"):
        return json.load(open(INFO_F))
    st = gemmi.read_structure(xc.fetch_assembly(PDB)); st.remove_ligands_and_waters()
    sub = gemmi.Structure(); sub.add_model(gemmi.Model("1"))
    for ch in st[0]:
        if ch.name in ("A", "B"):                         # one Vac8 (A) with its Nvj1 (B)
            sub[0].add_chain(ch)
    X, K, _ = sl.beads_from_structure(sub)
    ca = lambda c: np.array([r.find_atom("CA", "*").pos.tolist() for r in sub[0][c] if r.find_atom("CA", "*")]) / 10
    v8, nv = ca("A"), ca("B")
    z = nv[:5].mean(0) - v8[:5].mean(0); z /= np.linalg.norm(z)
    x = np.cross(z, [1.0, 0, 0]); x /= np.linalg.norm(x); y = np.cross(z, x)
    c = X.mean(0); Xl = (X - c) @ np.c_[x, y, z]
    bs.save(LIB_ID, Xl, K, extra=dict(pdb=PDB), per=30)
    info = dict(n_beads=int(len(Xl)), zmin=float(Xl[:, 2].min()), height=float(np.ptp(Xl[:, 2])))
    json.dump(info, open(INFO_F, "w"))
    log(f"  {LIB_ID} ({PDB}): {len(Xl)} beads, {info['height']:.1f} nm along the vacuole-to-nucleus axis")
    return info


def place(cell, lp, n_max, obstacles=None, seed=63, log=print):
    from scipy.spatial.transform import Rotation
    info = model(log=log)
    seed = _seed(seed); rng = np.random.default_rng(seed); n_max = int(round(n_max))
    li = int(list(lp["ids"]).index(LIB_ID))
    stc, str_ = lp["STC"][lp["st_off"][li]:lp["st_off"][li + 1]], lp["STR"][lp["st_off"][li]:lp["st_off"][li + 1]]
    z_lo = float((stc[:, 2] - str_).min())                   # lowest point of the sphere tree along +z
    G = cell.G; cn = G[3:6]; ext = float(G[11:14].max()) + 400
    glo, ghi = cn - ext, cn + ext
    a = cell.args(glo, ghi, margin=50.0); out = np.zeros(c3.NMEM)

    def sd(p):
        c3.mem_sdfs(p[0], p[1], p[2], *a, out); return out[c3.M_VAC], c3.nuc_sdf(p[0], p[1], p[2], G)
    # contact patch: vacuolar-membrane points close to the outer nuclear membrane
    pts, n_band = [], 0
    for p in cn + rng.uniform(-ext, ext, (400000, 3)):
        v, n = sd(p)
        if abs(v) < 20.0 and n < CONTACT + 20.0:             # near both: project onto the vacuolar midsurface
            n_band += abs(v) < 1.5
            nrm = np.array([sd(p + d)[0] - sd(p - d)[0] for d in np.eye(3) * 0.5]); nrm /= np.linalg.norm(nrm)
            q = p - nrm * v
            if abs(sd(q)[0]) < 0.5 and sd(q)[1] < CONTACT:
                pts.append(q)
    if not pts:
        log("nucleus-vacuole junction: no contact"); return {}, 0
    pts = np.array(pts)
    cap = 100000 + (len(obstacles[1]) if obstacles is not None else 0)
    G1 = p3.grid_new(glo, ghi, H, cap, cap * 14)
    if obstacles is not None and len(obstacles[1]):
        k = np.all((obstacles[0] > glo) & (obstacles[0] < ghi), axis=1)
        p3.g_insert(glo, H, *G1, obstacles[0][k], obstacles[1][k])
    M = np.empty((3, 3)); pos, rot = [], []
    for _ in range(int(n_max) * 30):
        if len(pos) >= n_max:
            break
        p = pts[rng.integers(len(pts))] + rng.normal(0, 8.0, 3)
        if sd(p)[1] > CONTACT + 5.0:
            continue
        e = 0.5
        nrm = np.array([sd(p + d)[0] - sd(p - d)[0] for d in np.eye(3) * e]); nrm /= np.linalg.norm(nrm)
        v0 = sd(p)[0]; base = p - nrm * v0                          # on the vacuolar midsurface; nrm points to the cytosol
        R = Rotation.from_rotvec(nrm * rng.uniform(0, 2 * np.pi)).as_matrix() @ _rot_z_to(nrm)
        t = base + nrm * (c3.HALF + GAP - z_lo)
        q = Rotation.from_matrix(R).as_quat(); q = np.array([q[3], q[0], q[1], q[2]])
        P = np.empty((len(stc), 3)); p3.quat_mat(q, M); p3.transform(stc, M, t, P)
        if p3.g_overlap(glo, H, *G1[:6], P, str_, p3.TOL):
            continue
        if not p3.check_volume(P, str_, t[0], t[1], t[2], float(lp["BR"][li]), 1 << c3.CYTO, *a, out):
            continue
        p3.g_insert(glo, H, *G1, P, str_); pos.append(t); rot.append(q)
    area = n_band / 400000 * (2 * ext) ** 3 / 3.0 / 1e6
    log(f"nucleus-vacuole junction: {len(pos)} of {int(n_max)} Vac8-Nvj1 on {area:.2f} um2 of vacuolar membrane within "
        f"{CONTACT:.0f} nm of the outer nuclear membrane")
    return dict(nvj_pos=np.array(pos, np.float32).reshape(-1, 3), nvj_q=np.array(rot, np.float32).reshape(-1, 4), li_nvj=li,
                nvj_centre=np.asarray(_on_patch(pos, pts), np.float32)), len(pos)


def _on_patch(pos, pts):
    """a point on the junction (the complex nearest the patch's mean; the mean of a curved patch lies off it)"""
    P = np.asarray(pos if len(pos) else pts, np.float64)
    return P[np.linalg.norm(P - P.mean(0), axis=1).argmin()]


def _rot_z_to(n):
    """rotation taking +z to unit n"""
    z = np.array([0, 0, 1.0]); v = np.cross(z, n); s = np.linalg.norm(v); c = float(z @ n)
    if s < 1e-9:
        return np.eye(3) if c > 0 else np.diag([1.0, -1, -1])
    K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]]) / s
    return np.eye(3) + s * K + (1 - c) * K @ K
