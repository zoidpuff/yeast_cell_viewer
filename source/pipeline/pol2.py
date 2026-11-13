"""Transcribing RNA polymerase II on its genes (literature item L14).

Per-gene numbers of engaged Pol II from genomic run-on (Pelechano et al. 2010 PLoS One 5:e15442, Table S1, "RNA pol II
presence per gene"; 691 per cell over 4,670 genes, about 2% of all Pol II). Gene coordinates from SGD
(SGD_features.tab). Each gene gets Poisson(its mean) polymerases at uniform positions along it.

On the chromatin (wc3d.rabl_segments arms; the nucleosome nearest in bp with the measured map, else kb mapped
linearly along each arm walk): the nucleosome at the
polymerase's position and its two neighbours are removed (nucleosomes are evicted around elongating Pol II), the
polymerase takes the middle one's place with its DNA axis (elongation complex 1Y1W on the library's 1WCM) along the
chain, and DNA (2 bp per bead) runs from the upstream nucleosome through the cleft to the downstream one. The nascent
transcript (as long as the polymerase has transcribed) is threaded from the RNA exit, 2 nt per bead. Engaged polymerases
come out of the free Pol II pool.
"""
import numpy as np
import pandas as pd
import pack3d as p3
import cell3d as c3
import rdna
import mrna
from config import DATA as D
from config import seed as _seed

PELECHANO_F = f"{D}/pelechano2010/files/pone.0015442.s002.xls"
SGD_F = f"{D}/sgd/SGD_features.tab"
NT_PER_BEAD, RNA_STEP, RNA_R = 2, 0.6, 0.7
EVICT = 1                       # nucleosomes removed on each side of the polymerase's own


def engaged(rng):
    """(chromosome 0-15, bp position, nascent length nt) per engaged polymerase"""
    p = pd.read_excel(PELECHANO_F, header=0)
    p = p.iloc[:, [0, 9]]; p.columns = ["orf", "pol"]; p["pol"] = pd.to_numeric(p.pol, errors="coerce")
    p = p.dropna()
    s = pd.read_csv(SGD_F, sep="\t", header=None, usecols=[1, 3, 8, 9, 10, 11], names=["type", "orf", "chrom", "a", "b", "strand"])
    s = s[(s.type == "ORF") & s.chrom.astype(str).str.isdigit()].drop_duplicates("orf").set_index("orf")
    out = []
    for orf, lam in zip(p.orf, p.pol):
        if orf not in s.index:
            continue
        r = s.loc[orf]; lo, hi = sorted((int(r.a), int(r.b)))
        for _ in range(rng.poisson(lam)):
            x = rng.uniform(lo, hi)
            out.append((int(r.chrom) - 1, x, (x - lo) if r.strand == "W" else (hi - x)))
    return out


def _seg_index(pre, chrom_kb, cen_kb, rdna_kb):
    """per segment: (chromosome, kind, nucleosome indices in walk order) and a mapper (chrom, bp) -> (segment, j)"""
    ch, K, C = pre["nuc_chain"], pre["seg_kind"], pre["seg_chrom"]
    idx = {s: np.where(ch == s)[0] for s in range(len(K))}
    def where(c, bp):
        kb = bp / 1000
        for s in range(len(K)):
            if C[s] != c or not len(idx[s]):
                continue
            k = str(K[s])
            if k == "L" and kb <= cen_kb[c]:
                f = (cen_kb[c] - kb) / cen_kb[c]
            elif k == "R" and kb > cen_kb[c]:
                f = (kb - cen_kb[c]) / (chrom_kb[c] - cen_kb[c])
            elif k == "XIIR_in" and cen_kb[c] < kb < rdna_kb[0]:
                f = (kb - cen_kb[c]) / (rdna_kb[0] - cen_kb[c])
            elif k == "XIIR_out" and kb > rdna_kb[1]:
                f = (kb - rdna_kb[1]) / (chrom_kb[c] - rdna_kb[1])
            else:
                continue
            if "nuc_bp" in pre:         # measured nucleosome map: the nucleosome nearest the position
                return s, int(np.argmin(np.abs(pre["nuc_bp"][idx[s]] - bp)))
            return s, int(np.clip(round(f * (pre["chain_target"][s] - 1)), 0, len(idx[s]) - 1))
        return None
    return idx, where


def place(cell, lp, pre, n_pool, obstacles, seed=91, log=print):
    """returns (updated nucleosome arrays, pol2 instances, DNA pieces, nascent RNA threads)"""
    from scipy.spatial.transform import Rotation
    import wc3d
    seed = _seed(seed); rng = np.random.default_rng(seed)
    c_dna, u_dna, e_exit, d_exit = rdna.pol_frame(log=log, lib_id="PolII", lib_pdb="1WCM", ec_pdb="1Y1W")
    ids = list(lp["ids"]); li_pol, li_nuc = ids.index("PolII"), int(pre["li_nuc"])
    s0, s1 = lp["st_off"][li_pol], lp["st_off"][li_pol + 1]; STC, STR = lp["STC"][s0:s1], lp["STR"][s0:s1]
    want = engaged(rng)
    if len(want) > n_pool:
        want = [want[i] for i in rng.choice(len(want), int(n_pool), replace=False)]
    idx, where = _seg_index(pre, wc3d.CHROM_KB, wc3d.CEN_KB, wc3d.RDNA_KB)
    pos, q = pre["nuc_pos"].astype(np.float64), pre["nuc_q"]
    keep = np.ones(len(pos), bool)
    cn = cell.G[3:6]; rn = float(cell.G[11:14].max()); glo, ghi = cn - rn - 80, cn + rn + 80
    a = cell.args(glo, ghi, margin=200.0)
    NC, NR = p3.emit_spheres(np.full(len(pos), li_nuc, np.int32), pre["nuc_pos"], pre["nuc_q"], lp["st_off"], lp["STC"], lp["STR"],
                             lp["BR"], glo, ghi)
    per = len(lp["STC"][lp["st_off"][li_nuc]:lp["st_off"][li_nuc + 1]])
    owner = np.repeat(np.arange(len(pos)), per)          # emit_spheres keeps instance order (all inside the box)
    from scipy.spatial import cKDTree
    tree = cKDTree(NC)
    OC, OR = obstacles
    k = np.all((OC > glo) & (OC < ghi), axis=1); OC, OR = OC[k].astype(np.float64), OR[k].astype(np.float64)
    otree = cKDTree(OC) if len(OC) else None
    M = np.empty((3, 3)); out = np.zeros(c3.NMEM)
    pol_pos, pol_q, dna, rs, rd, rl = [], [], [], [], [], []
    placed_P = []
    for c, bp, nt in want:
        w = where(c, bp)
        if w is None:
            continue
        s, j = w; I = idx[s]
        if j - EVICT - 1 < 0 or j + EVICT + 1 >= len(I):
            continue
        ii = I[j - EVICT:j + EVICT + 1]
        if not keep[ii].all():
            continue
        up, dn = I[j - EVICT - 1], I[j + EVICT + 1]
        S = pos[I[j]]; u = pos[dn] - pos[up]; u /= np.linalg.norm(u)
        ok = False
        for _ in range(12):
            R = rdna._rot_u_to(u_dna, u, rng.uniform(0, 2 * np.pi))
            t = S - R @ c_dna
            P = STC @ R.T + t
            near = tree.query_ball_point(P, (STR + NR.max()).max())
            hit = False
            for pi, nb in enumerate(near):
                nb = [n for n in nb if keep[owner[n]] and owner[n] not in ii and np.linalg.norm(P[pi] - NC[n]) < STR[pi] + NR[n] - p3.TOL]
                if nb:
                    hit = True; break
            if not hit and otree is not None:
                hit = any(np.linalg.norm(P[pi] - OC[n]) < STR[pi] + OR[n] - p3.TOL for pi, nb in enumerate(otree.query_ball_point(P, (STR + OR.max()).max())) for n in nb)
            if not hit:
                hit = any((np.linalg.norm(P[:, None] - Q[None], axis=2) < STR[:, None] + STR[None] - p3.TOL).any() for Q in placed_P[-200:]
                          if np.linalg.norm(Q.mean(0) - P.mean(0)) < 40)
            if not hit and p3.check_volume(P, STR, t[0], t[1], t[2], float(lp["BR"][li_pol]), 1 << c3.NUCP, *a, out):
                ok = True; break
        if not ok:
            continue
        keep[ii] = False; placed_P.append(P)
        pol_pos.append(t); pol_q.append(rdna._quat(R))
        a_in, a_out = S - u * 6.0, S + u * 6.0               # DNA through the cleft to the flanking nucleosomes
        for p0, p1 in ((pos[up], a_in), (a_in, a_out), (a_out, pos[dn])):
            L = np.linalg.norm(p1 - p0); n_ = max(2, int(L / rdna.DNA_STEP))
            dna.append(p0 + (p1 - p0) * np.linspace(0, 1, n_)[:, None])
        if nt >= 2 * NT_PER_BEAD:
            rs.append(S - R @ c_dna + R @ e_exit); rd.append(R @ d_exit); rl.append(int(nt // NT_PER_BEAD))
    pol_pos, pol_q = np.array(pol_pos).reshape(-1, 3), np.array(pol_q).reshape(-1, 4)
    # nascent transcripts, threaded in the nucleoplasm clear of everything placed
    PC, PR = p3.emit_spheres(np.full(len(pol_pos), li_pol, np.int32), pol_pos.astype(np.float32), pol_q.astype(np.float32),
                             lp["st_off"], lp["STC"], lp["STR"], lp["BR"], glo, ghi)
    k2 = keep[owner]
    D_ = np.concatenate(dna) if dna else np.zeros((0, 3))
    OC2 = np.concatenate([NC[k2], PC, OC, D_]).astype(np.float64)
    OR2 = np.concatenate([NR[k2], PR, OR, np.full(len(D_), 1.0)]).astype(np.float64)
    lens = np.array(rl, np.int64)
    cap_s = len(OR2) + int(lens.sum()) + 100000
    if len(lens):
        rp, ro = mrna.thread_all(seed, glo, ghi, *a, OC2, OR2, np.array(rs), np.array(rd), lens, 0.0, cap_s,
                                 int(cap_s * 40), RNA_STEP, RNA_R, 0.8, 1 << c3.NUCP, 200.0, 30)
    else:
        rp, ro = np.zeros((0, 3), np.float32), np.zeros(1, np.int64)
    log(f"engaged Pol II: {len(pol_pos)} of {len(want)} on their genes (Pelechano 2010 per-gene counts); "
        f"{(~keep).sum()} nucleosomes evicted; nascent RNA {len(rp) * NT_PER_BEAD / 1e3:.0f} of {lens.sum() * NT_PER_BEAD / 1e3:.0f} knt")
    off = np.r_[0, np.cumsum([len(x) for x in dna])].astype(np.int64)
    return dict(keep=keep, pol2_pos=pol_pos.astype(np.float32), pol2_q=pol_q.astype(np.float32), li_pol2=li_pol,
                pol2dna_pts=D_.astype(np.float32), pol2dna_off=off, pol2rna_pts=np.asarray(rp, np.float32), pol2rna_off=np.asarray(ro, np.int64))
