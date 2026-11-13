"""Transport receptors in the nuclear pore channel (IMPROVEMENTS N5; the Jurkat model's 20 per pore).

Karyopherins and mRNA exporters spend a large share of their time inside the FG phase of the pores, which the packer
never fills (it places the FG mesh, then packs only the cytoplasm and nucleoplasm around it). Here N_PER_PORE receptors
per pore (an assumption, as in the Jurkat cell model; unchecked for yeast) are drawn from the free pools of the
receptor types below in proportion to their copies, and placed in the channel: centres within CH_R nm of the pore axis
and CH_Z nm of the envelope midplane, clear of the NPC scaffold. The FG domains are threaded afterwards around them
(fg_mesh.py). The rest of each pool stays in the cytoplasm and nucleoplasm.
"""
import numpy as np
from config import seed as _seed
import pack3d as p3
import cell3d as c3

RECEPTORS = ["KAP95", "SRP1", "KAP123", "PSE1", "CPX_26769_RAN_RANBP1_CRM1_nuclear_export_c",
             "CPX_1142_MEX67_MTR2_mRNA_nuclear_export_f", "NTF2_oligo"]
N_PER_PORE = 20
CH_R, CH_Z = 22.0, 25.0
TRIES = 300


def place(cell, lp, pools, obstacles, seed=73, log=print):
    """pools: {lib_id: free copies}. Returns pre-pass group 'kap' (kap_type, kap_pos, kap_q), per-type counts and
    the placed spheres."""
    from scipy.spatial.transform import Rotation
    ids = list(lp["ids"])
    seed = _seed(seed); rng = np.random.default_rng(seed)
    types_ok = [g for g in RECEPTORS if g in ids and pools.get(g, 0) >= 1]
    if not types_ok or not len(cell.pores):
        return {}, {}, (np.zeros((0, 3)), np.zeros(0))
    w = np.array([pools[g] for g in types_ok], float); w /= w.sum()
    G = cell.G; cn = G[3:6]; rn = float(G[11:14].max()); glo, ghi = cn - rn - 100, cn + rn + 100
    a = cell.args(glo, ghi, margin=200.0)
    OC, OR = obstacles
    grid = p3.LocalGrid(glo, ghi, OC.astype(np.float64), OR.astype(np.float64))
    out = np.zeros(c3.NMEM)
    mask = (1 << c3.PORE) | (1 << c3.NUCP) | (1 << c3.CYTO)
    types, pos, quat, counts, sph_c, sph_r = [], [], [], {}, [], []
    for pr in cell.pores:
        c0, u = pr[:3].astype(float), pr[3:6].astype(float); u /= np.linalg.norm(u)
        e1 = np.cross(u, [0.3, 0.5, 0.8]); e1 /= np.linalg.norm(e1); e2 = np.cross(u, e1)
        for g in rng.choice(types_ok, N_PER_PORE, p=w):
            if counts.get(g, 0) >= pools[g]:
                continue
            li = ids.index(g)
            C = lp["STC"][lp["st_off"][li]:lp["st_off"][li + 1]].astype(np.float64)
            R = lp["STR"][lp["st_off"][li]:lp["st_off"][li + 1]].astype(np.float64)
            for _ in range(TRIES):
                r = CH_R * np.sqrt(rng.random()); ph = rng.uniform(0, 2 * np.pi); z = rng.uniform(-CH_Z, CH_Z)
                ctr = c0 + z * u + r * (np.cos(ph) * e1 + np.sin(ph) * e2)
                Rm = Rotation.random(random_state=rng).as_matrix()
                t = ctr - Rm @ C.mean(0); S = C @ Rm.T + t
                if grid.free(S, R) and p3.check_volume(S, R, t[0], t[1], t[2], float(lp["BR"][li]), mask, *a, out):
                    grid.add(S, R); sph_c.append(S); sph_r.append(R)
                    types.append(li); pos.append(t); quat.append(np.roll(Rotation.from_matrix(Rm).as_quat(), 1))
                    counts[g] = counts.get(g, 0) + 1
                    break
    log(f"pore receptors: {len(types)} of {N_PER_PORE * len(cell.pores)} in the channels of {len(cell.pores)} pores ("
        + ", ".join(f"{g.split('_')[0] if not g.startswith('CPX') else g.split('_')[2]} {n}" for g, n in counts.items()) + ")")
    if not types:
        return {}, counts, (np.zeros((0, 3)), np.zeros(0))
    return (dict(kap_type=np.array(types, np.int32), kap_pos=np.array(pos, np.float32), kap_q=np.array(quat, np.float32)),
            counts, (np.concatenate(sph_c), np.concatenate(sph_r)))
