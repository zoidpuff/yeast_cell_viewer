"""Chromatin folding of the packed cell against in-vivo measurements (IMPROVEMENTS section K).

Model: nucleosome walks of the 32 chromosome arms (prepass.npz: nuc_pos, nuc_chain; walk order = genomic order,
165 bp per nucleosome, or the measured dyad positions nuc_bp with the nucleosome map; MODEL.md 6.3).
  R(s)  root-mean-square 3D distance between nucleosomes s kb apart on the same arm
  P(s)  probability that two nucleosomes s kb apart are within 30 nm (centre distance), and its log-log slope
Data: Bystricky et al. 2004 PNAS 101:16495 (FISH at 14-100 kb and live imaging of whole arms) fitted a worm-like chain
with persistence length 170-220 nm and mass density 110-150 bp/nm; R(s) of that chain over the parameter ranges is the
comparison (a derived comparison: the paper's raw distances are not in its open abstract). Arbona et al. 2017 Genome
Biol 18:81 inferred 53-65 bp/nm and persistence 52-85 nm from imaging and Hi-C together; that chain's R(s) is the
second range (it used Hi-C data, so it is not independent of Hi-C tests).
Output: work/validation/chromatin_stats.json
"""
import os, sys, json, math
import numpy as np
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
from config import WORK, W3D

NRL = 165
S_KB = (2, 5, 10, 14, 20, 50, 100, 200)


def wlc_rms(s_kb, lp, bp_per_nm):
    L = s_kb * 1000 / bp_per_nm
    return math.sqrt(2 * lp * L * (1 - (lp / L) * (1 - math.exp(-L / lp))))


def main():
    z = np.load(f"{W3D}/prepass.npz")
    P, ch = z["nuc_pos"].astype(float), z["nuc_chain"]
    # genomic position along each walk: measured dyads (nuc_bp) or the uniform repeat
    gp = z["nuc_bp"].astype(float) if "nuc_bp" in z.files else np.concatenate([np.arange((ch == c).sum()) * NRL for c in np.unique(ch)])
    arms = np.unique(ch)
    R = {}
    for s in S_KB:
        d = []
        for c in arms:
            m = ch == c; x, g = P[m], gp[m]
            o = np.argsort(g); x, g = x[o], g[o]
            j = np.searchsorted(g, g + s * 1000)          # first nucleosome at least s kb further along
            ok = j < len(g)
            if ok.any():
                d.append(np.linalg.norm(x[j[ok]] - x[ok], axis=1))
        R[s] = float(np.sqrt((np.concatenate(d) ** 2).mean()))
    wlc = {s: [min(wlc_rms(s, lp, c) for lp in (170, 220) for c in (110, 150)),
               max(wlc_rms(s, lp, c) for lp in (170, 220) for c in (110, 150))] for s in S_KB if s >= 14}
    arb = {s: [min(wlc_rms(s, lp, c) for lp in (52, 85) for c in (53, 65)),
               max(wlc_rms(s, lp, c) for lp in (52, 85) for c in (53, 65))] for s in S_KB if s >= 5}
    pairs = cKDTree(P).query_pairs(30.0, output_type="ndarray")
    i, j = pairs[:, 0], pairs[:, 1]
    cis = ch[i] == ch[j]
    sep = np.abs(gp[i] - gp[j])[cis] / 1000
    bins = np.logspace(0, np.log10(500), 15); mid = np.sqrt(bins[1:] * bins[:-1])
    h, _ = np.histogram(sep, bins)
    npairs = np.array([sum(max(0, (ch == c).sum() - int(b * 1000 / NRL)) for c in arms) for b in mid])
    ps = h / np.diff(bins * 1000 / NRL) / np.maximum(npairs, 1)
    ok = (mid > 10) & (mid < 200) & (ps > 0)
    res = dict(arms=int(len(arms)), nucleosomes=int(len(P)), rms_distance_nm=R, bystricky2004_wlc_rms_range_nm=wlc, arbona2017_wlc_rms_range_nm=arb,
               contacts_30nm=int(len(pairs)), cis_arm_fraction=float(cis.mean()),
               ps_slope_10_200kb=float(np.polyfit(np.log(mid[ok]), np.log(ps[ok]), 1)[0]),
               ps=dict(s_kb=mid.round(2).tolist(), p=ps.tolist()))
    json.dump(res, open(f"{WORK}/validation/chromatin_stats.json", "w"), indent=1)
    print({k: v for k, v in res.items() if k != "ps"})
    return res


if __name__ == "__main__":
    main()
