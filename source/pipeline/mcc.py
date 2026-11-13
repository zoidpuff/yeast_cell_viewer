"""MCC membrane proteins in the eisosome furrows (open item of L7).

The eisosome furrows are the membrane compartment of Can1 (MCC). Sur7 and Nce102, tetraspanning MCC markers, are
concentrated there (Malinska et al. 2003 Mol Biol Cell 14:4427; Grossmann et al. 2008 J Cell Biol 183:1075). Here
MCC_FRACTION of their copies is placed on the flat walls of the furrows (the half-cylinder bottom is under the Pil1/Lsp1
coat), oriented like any plasma-membrane protein (library z along the local normal, transmembrane residues in the
bilayer: pack3d.check_membrane), as many as fit; the rest stays in the rest of the plasma membrane. Can1 and other
nutrient transporters, which move in and out of the MCC with conditions, are not concentrated.
"""
import numpy as np
from config import seed as _seed
import pack3d as p3
import cell3d as c3

GENES = ("SUR7", "NCE102")
MCC_FRACTION = 0.8          # share of their copies in the furrows (assumption: both are MCC markers in normal growth)
H = 20.0


def place(cell, lp, want, obstacles=None, seed=73, log=print):
    """want: library id -> copies to place in the furrows"""
    seed = _seed(seed); rng = np.random.default_rng(seed)
    ids = list(lp["ids"])
    sub = {k: (ids.index(k), lp["STC"][lp["st_off"][ids.index(k)]:lp["st_off"][ids.index(k) + 1]],
               lp["STR"][lp["st_off"][ids.index(k)]:lp["st_off"][ids.index(k) + 1]]) for k in want if k in ids}
    G = cell.G; glo = np.full(3, -G[0] - 10.0); ghi = -glo
    a = cell.args(glo, ghi, margin=50.0); out = np.zeros(c3.NMEM)
    cap = 200000 + (len(obstacles[1]) if obstacles is not None else 0)
    G1 = p3.grid_new(glo, ghi, H, cap, cap * 14)
    if obstacles is not None and len(obstacles[1]):
        k = np.all((obstacles[0] > glo) & (obstacles[0] < ghi), axis=1)
        p3.g_insert(glo, H, *G1, obstacles[0][k], obstacles[1][k])
    furrows = [cell.PQ[i] for i in cell.furrows()]
    order = [k for k in sub for _ in range(int(round(want[k])))]
    rng.shuffle(order)
    M = np.empty((3, 3)); q = np.empty(4); t = np.empty(3); nrm = np.empty(3)
    pos, rot, typ = [], [], []
    for lid in order:
        li, stc, str_ = sub[lid]; ZL = np.ascontiguousarray(stc[:, 2])
        for _ in range(60):
            f = furrows[rng.integers(len(furrows))]
            A, B, r, h = f[:3], f[3:6], f[6], f[7]
            L = np.linalg.norm(B - A); e = (B - A) / L; U = (A + B) / np.linalg.norm(A + B); w = np.cross(e, U)
            t[:] = A + e * rng.uniform(0, L) + U * rng.uniform(4.0, h - 15.0) + w * r * rng.choice([-1.0, 1.0])
            if not p3.project_to_membrane(t, c3.M_PM, nrm, *a, out):
                continue
            p3.quat_z_to(nrm.copy(), rng.uniform(0, 2 * np.pi), q); p3.quat_mat(q, M)
            P = np.empty((len(stc), 3)); p3.transform(stc, M, t, P)
            if p3.g_overlap(glo, H, *G1[:6], P, str_, p3.TOL):
                continue
            if not p3.check_membrane(P, str_, ZL, c3.M_PM, *a, out):
                continue
            p3.g_insert(glo, H, *G1, P, str_); pos.append(t.copy()); rot.append(q.copy()); typ.append(li)
            break
    used = {k: int(sum(1 for x in typ if x == sub[k][0])) for k in sub}
    log("MCC proteins in the eisosome furrows: " + ", ".join(f"{k} {used[k]} of {int(round(want[k]))}" for k in sub))
    return dict(mcc_pos=np.array(pos, np.float32).reshape(-1, 3), mcc_q=np.array(rot, np.float32).reshape(-1, 4),
                mcc_type=np.array(typ, np.int32)), used
