"""Chromatin as a relaxed coarse polymer against the fibre walk (IMPROVEMENTS N2): reproduces the comparison behind
the decision to keep the fibre walk. Not part of validate_all (needs OpenMM on a GPU; about 1-2 min per replica).

  python validation/chromatin_polymer_scan.py [--reps 8] [--steps 400000] [--r-guide 40 60] [--gpu 0]

For each replica: chromatin_polymer.simulate (seed = replica) in the current cell with the current pre-pass obstacles
(pores, microtubules, kinetochores), the measured nucleosomes back-mapped along it (pack3d.guided_walks), then the
Hi-C measures of validation/hic_test.py on the ICE-balanced model matrix, per cell and for the summed ensemble, against
the ICE-balanced maps. R(s) at 20/100/200 kb is printed for the imaging envelopes (chromatin_stats.py).
Output: work/validation/chromatin_polymer_scan.json
"""
import os, sys, json, time, argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline")); sys.path.insert(0, HERE)
from config import WORK, W3D  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=8); ap.add_argument("--steps", type=int, default=400_000)
    ap.add_argument("--r-guide", type=float, nargs="+", default=[40.0, 60.0]); ap.add_argument("--gpu", default="0")
    a = ap.parse_args()
    os.environ.setdefault("CUDA_VISIBLE_DEVICES", a.gpu)
    import cell3d as c3, wc3d, pack3d as p3, genome_map, nucleosome_map
    import chromatin_polymer as cp, hic_test as ht
    from scipy.spatial import cKDTree
    from scipy.stats import spearmanr
    cell = c3.Cell.load(f"{W3D}/geometry.json")
    z = dict(np.load(f"{W3D}/prepass.npz"))
    ids = [str(x) for x in np.load(f"{W3D}/lib_ids.npy")]
    lp = wc3d.build_libpack(ids); li_nuc = ids.index("nucleosome")
    segs = genome_map.segments()
    ss, st, sk, sc = z["seg_start"], z["seg_target"], z["seg_kind"], z["seg_chrom"]
    starts = np.array([ss[[i for i in range(len(sk)) if sc[i] == c and sk[i] == "L"][0]] for c in range(16)])
    t_in = st[list(sk).index("XIIR_in")]; t_out = ss[list(sk).index("XIIR_out")]
    OCs, ORs = [], []
    for k, li in [("npc", "li_npc"), ("mt", None), ("kt", None), ("kto", None)]:
        if f"{k}_pos" not in z:
            continue
        ty = z[f"{k}_type"].astype(np.int32) if f"{k}_type" in z else np.full(len(z[f"{k}_pos"]), int(z[li]), np.int32)
        C_, R_ = p3.emit_spheres(ty, z[f"{k}_pos"].astype(np.float32), z[f"{k}_q"].astype(np.float32), lp["st_off"], lp["STC"],
                                 lp["STR"], lp["BR"], np.full(3, -1e4), np.full(3, 1e4))
        OCs.append(C_); ORs.append(R_)
    OC = np.concatenate(OCs); OR = np.concatenate(ORs)
    grid = cp.violation_grid(cell, (np.concatenate([OC, z["fg_pts"][::4].astype(float)]),
                                    np.concatenate([OR, np.full(len(z["fg_pts"][::4]), 1.0)])))
    cn = cell.G[3:6]; rn = float(cell.G[11:14].max()); glo = cn - rn - 80; ghi = cn + rn + 80
    k = np.all((OC > glo) & (OC < ghi), axis=1); OCw, ORw = OC[k].astype(float), OR[k].astype(float)
    args = cell.args(glo, ghi, margin=200.0)
    dy = nucleosome_map.segment_dyads(segs)
    sl = np.array([len(p) for p in dy], np.int64)
    _, lk = wc3d.fiber_steps(sl, np.concatenate([np.r_[wc3d.NRL_BP, np.abs(np.diff(p))] for p in dy]))
    s0, s1 = lp["st_off"][li_nuc], lp["st_off"][li_nuc + 1]; nsph = (s1 - s0) * sl.sum() + len(ORw)
    allbp = np.concatenate(dy); allchrom = np.concatenate([np.full(len(p), segs[i][0] + 1) for i, p in enumerate(dy)])
    off = ht.bins_genome(); n = off[-1]
    chrom_m = np.concatenate([np.full(off[c + 1] - off[c], c + 1) for c in range(16)])
    mid_m = np.concatenate([np.arange(off[c + 1] - off[c]) * ht.BIN + ht.BIN // 2 for c in range(16)])
    iu = np.triu_indices(16, 1)
    data = {}
    for nm, f in (("G1", ht.lazar_matrix), ("log", ht.swygert_matrix)):
        M_, c_, m_ = f(); data[nm] = ht.measures(ht.ice(M_), c_, m_)

    def score(M):
        m = ht.measures(ht.ice(M), chrom_m, mid_m)
        r = {k_: m[k_] for k_ in ("cis_fraction", "ps_slope", "cen_enrichment", "tel_enrichment")}
        for nm in data:
            r[f"pair_r_{nm}"] = float(spearmanr(np.array(m["trans_norm"])[iu], np.array(data[nm]["trans_norm"])[iu]).correlation)
        return r

    polys = [cp.simulate(cell, segs, starts, t_in, t_out, None, seed=s, steps=a.steps, grid=grid) for s in range(a.reps)]
    res = dict(steps=a.steps, reps=a.reps, params=dict(bead_bp=cp.BEAD_BP, bp_per_nm=cp.BP_PER_NM, lp_nm=cp.LP_NM,
               sigma_nm=cp.SIGMA_NM, e_rep=cp.E_REP), data_ice={nm: {k_: data[nm][k_] for k_ in
               ("cis_fraction", "ps_slope", "cen_enrichment", "tel_enrichment")} for nm in data}, runs={})
    for rg in a.r_guide:
        Msum = np.zeros((n, n)); cells, rs_all = [], []
        for X, bp, seg in polys:
            guide = cp.guide_points(X, bp, seg, dy)
            pos, qu, ch, brk, src = p3.guided_walks(7, glo, ghi, *args, lp["STC"][s0:s1], lp["STR"][s0:s1], lp["BR"][li_nuc],
                                                    1 << c3.NUCP, OCw, ORw, sl, guide, rg, lk, int(nsph * 1.1) + 1000, int(nsph * 10) + 10000)
            P = pos.astype(float)
            b = off[allchrom[src] - 1] + np.minimum(allbp[src] // ht.BIN, np.diff(off)[allchrom[src] - 1] - 1)
            pr = cKDTree(P).query_pairs(ht.CONTACT_NM, output_type="ndarray")
            M = np.zeros((n, n)); np.add.at(M, (b[pr[:, 0]], b[pr[:, 1]]), 1); np.add.at(M, (b[pr[:, 1]], b[pr[:, 0]]), 1)
            Msum += M; cells.append(score(M))
            rs = {}
            for s_kb in (20, 100, 200):
                d2 = []
                for c in range(len(sl)):
                    I = np.where(ch == c)[0]
                    q = allbp[src[I]].astype(float); o = np.argsort(q); q = q[o]; Pi = P[I][o]
                    j = np.searchsorted(q, q + s_kb * 1000); ok = j < len(q)
                    d2.append(np.sum((Pi[j[ok]] - Pi[ok]) ** 2, 1))
                rs[s_kb] = float(np.sqrt(np.mean(np.concatenate(d2))))
            rs_all.append(rs)
        res["runs"][str(rg)] = dict(cells=cells, ensemble=score(Msum),
                                    rms_nm={k_: float(np.mean([r[k_] for r in rs_all])) for k_ in (20, 100, 200)})
        print(f"r_guide {rg}: ensemble {res['runs'][str(rg)]['ensemble']}; R(s) {res['runs'][str(rg)]['rms_nm']}", flush=True)
    json.dump(res, open(f"{WORK}/validation/chromatin_polymer_scan.json", "w"), indent=1)
    return res


if __name__ == "__main__":
    main()
