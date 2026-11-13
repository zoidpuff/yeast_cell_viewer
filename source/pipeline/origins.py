"""Licensed replication origins of the G1 cell (IMPROVEMENTS K6).

In G1, ORC sits on each origin and has loaded the MCM2-7 double hexamer around the origin DNA. Each confirmed origin
of OriDB (Siow et al. 2012 Nucleic Acids Res 40:D682; 410 'Confirmed', coordinates sacCer3) gets one MCM double
hexamer (CPX-2944, 7W8G) in the gap between the two measured nucleosomes that flank the origin's midpoint, its long
axis (the DNA channel) along the line joining them, and one ORC (CPX-768, 7TJH) next to it on the chromatin.
One double hexamer per origin is a lower bound: efficient origins load several (unchecked); the rest of the MCM and
ORC pools stay free in the nucleoplasm. Requires the measured nucleosome map (prepass nuc_bp).
The fibre walk packs nucleosomes too densely for a 21 nm double hexamer between most flanking pairs: it is moved up to
25 nm out of the fibre (the origin's free DNA looping out), and where that fails too the two flanking nucleosomes (last:
two on each side) are evicted (returned in `keep`, as for Pol II).
"""
import numpy as np
import pandas as pd
import pack3d as p3
import cell3d as c3
import rdna
from config import DATA
from config import seed as _seed

ORIDB = f"{DATA}/validation/genome/origins_oridb/oridb_all_origins.tsv"
MCM_ID, ORC_ID = "CPX_2944_MCM_complex", "CPX_768_Nuclear_origin_recognition_compl"


def confirmed():
    o = pd.read_csv(ORIDB, sep="\t")
    o = o[o.status == "Confirmed"]
    return pd.DataFrame(dict(chrom=o.chr.astype(int), mid=((o.start + o.end) // 2).astype(int), name=o.name))


def place(cell, lp, pre, n_mcm, n_orc, obstacles, seed=77, log=print):
    from scipy.spatial import cKDTree
    from scipy.spatial.transform import Rotation
    import genome_map
    ids = list(lp["ids"])
    if "nuc_bp" not in pre or MCM_ID not in ids:
        return dict(keep=np.ones(len(pre["nuc_pos"]), bool))
    seed = _seed(seed); rng = np.random.default_rng(seed)
    li_m, li_o, li_nuc = ids.index(MCM_ID), (ids.index(ORC_ID) if ORC_ID in ids else -1), int(pre["li_nuc"])
    lib = lambda li: (lp["STC"][lp["st_off"][li]:lp["st_off"][li + 1]].astype(np.float64), lp["STR"][lp["st_off"][li]:lp["st_off"][li + 1]].astype(np.float64))
    MC, MR = lib(li_m)
    u_m = np.linalg.eigh(np.cov(MC.T))[1][:, 2]          # double hexamer long axis = DNA channel
    pos, ch, bp = pre["nuc_pos"].astype(np.float64), pre["nuc_chain"], pre["nuc_bp"]
    segs = genome_map.segments()
    cn = cell.G[3:6]; rn = float(cell.G[11:14].max()); glo, ghi = cn - rn - 80, cn + rn + 80
    a = cell.args(glo, ghi, margin=200.0)
    NC, NR = p3.emit_spheres(np.full(len(pos), li_nuc, np.int32), pre["nuc_pos"], pre["nuc_q"], lp["st_off"], lp["STC"], lp["STR"],
                             lp["BR"], glo, ghi)
    OC, OR = obstacles
    k = np.all((OC > glo) & (OC < ghi), axis=1)
    AC = np.concatenate([NC, OC[k]]).astype(np.float64); AR = np.concatenate([NR, OR[k]]).astype(np.float64)
    tree = cKDTree(AC); rmax = AR.max()
    per = lp["st_off"][li_nuc + 1] - lp["st_off"][li_nuc]
    assert len(NC) == per * len(pos)
    owner = np.r_[np.repeat(np.arange(len(pos)), per), np.full(len(AC) - len(NC), -1)]
    keep = np.ones(len(pos), bool)
    placed = []                         # sphere sets placed here
    out = np.zeros(c3.NMEM)

    def free(P, R, t, li, skip=()):
        for nb, p_, r_ in zip(tree.query_ball_point(P, R.max() + rmax), P, R):
            if not nb:
                continue
            nb = np.array(nb); w = owner[nb]
            nb = nb[(w < 0) | (keep[np.maximum(w, 0)] & ~np.isin(w, skip))]
            if len(nb) and (np.linalg.norm(AC[nb] - p_, axis=1) < AR[nb] + r_ - p3.TOL).any():
                return False
        for Q, QR in placed[-50:]:
            if np.linalg.norm(Q.mean(0) - P.mean(0)) < 40 and (np.linalg.norm(P[:, None] - Q[None], axis=2) < R[:, None] + QR[None] - p3.TOL).any():
                return False
        return p3.check_volume(P, R, t[0], t[1], t[2], float(lp["BR"][li]), 1 << c3.NUCP, *a, out)

    o = confirmed()
    if len(o) > n_mcm:
        o = o.iloc[np.sort(rng.choice(len(o), int(n_mcm), replace=False))]
    m_pos, m_q, o_pos, o_q, names = [], [], [], [], []
    for r in o.itertuples():
        c = r.chrom - 1
        s = [i for i, (cc, kind, s0, s1) in enumerate(segs) if cc == c and min(s0, s1) <= r.mid < max(s0, s1)]
        if not s:
            continue
        I = np.where(ch == s[0])[0]
        g = bp[I]; o_ = np.argsort(g); I, g = I[o_], g[o_]
        j = np.searchsorted(g, r.mid)
        if j == 0 or j == len(g):
            continue
        up, dn = pos[I[j - 1]], pos[I[j]]
        u = dn - up; L = np.linalg.norm(u); u /= L
        # between the flanking nucleosomes, then further out of the fibre (the origin's free DNA looping out, up to
        # ~25 nm), then with the two flanking nucleosomes evicted, last with two on each side
        ok = False; flank = (I[j - 1], I[j]); flank2 = tuple(I[max(0, j - 2):j + 2])
        for k_ in range(96):
            skip = flank2 if k_ >= 72 else flank if k_ >= 48 else ()
            e = rng.normal(size=3); e -= (e @ u) * u; e /= np.linalg.norm(e)
            S = 0.5 * (up + dn) + (e * min(25.0, k_ * 0.6) if k_ < 48 else rng.normal(0, 2.0, 3))
            R = rdna._rot_u_to(u_m, u, rng.uniform(0, 2 * np.pi))
            t = S - R @ MC.mean(0)
            P = MC @ R.T + t
            if free(P, MR, t, li_m, skip):
                ok = True; break
        if not ok:
            continue
        if len(skip):
            keep[list(skip)] = False
        placed.append((P, MR)); m_pos.append(t); m_q.append(rdna._quat(R)); names.append(r.name)
        if li_o < 0 or len(o_pos) >= n_orc:
            continue
        OCc, ORr = lib(li_o)
        for _ in range(30):                 # ORC beside the double hexamer, against the origin DNA
            d = rng.normal(size=3); d -= (d @ u) * u; d /= np.linalg.norm(d)
            t = S + u * rng.uniform(-16, 16) + d * rng.uniform(12, 20)
            R2 = Rotation.random(random_state=rng).as_matrix()
            P2 = OCc @ R2.T + t
            if free(P2, ORr, t, li_o):
                placed.append((P2, ORr)); o_pos.append(t); o_q.append(rdna._quat(R2)); break
    log(f"origins: {len(m_pos)} MCM double hexamers and {len(o_pos)} ORC on {len(o)} confirmed origins (OriDB); "
        f"{(~keep).sum()} flanking nucleosomes evicted")
    f32 = lambda x, n: np.asarray(x, np.float32).reshape(-1, n)
    res = dict(keep=keep, mcm_pos=f32(m_pos, 3), mcm_q=f32(m_q, 4), li_mcm=li_m)
    if li_o >= 0:
        res.update(orc_pos=f32(o_pos, 3), orc_q=f32(o_q, 4), li_orc=li_o)
    return res
