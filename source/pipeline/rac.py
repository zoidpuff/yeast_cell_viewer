"""Ribosome-associated complex on the ribosome (literature item L15).

RAC (the Hsp40 Zuo1 with the Hsp70 Ssz1) sits on the 60S subunit at the polypeptide exit tunnel, Zuo1 reaching over
to the 40S, and recruits the Hsp70 Ssb to nascent chains. It is ribosome-bound in vivo (co-sedimentation; Gautschi et
al. 2001), so here it rides on ribosomes instead of floating free.

Structures: the yeast RAC of 7X3K (cryo-EM of the S. cerevisiae RNC-RAC complex; only RAC deposited) is placed on the
library 80S (4V88) through the C. thermophilum 80S-RAC complex 7Z3N (Kisonaite et al. 2023): 7Z3N's rRNA is superposed
on 4V88's by sequence-aligned phosphates, and 7X3K's Zuo1/Ssz1 on 7Z3N's by sequence-aligned C-alphas. RAC beads that
would sit inside the ribosome (flexible parts at the interface) are dropped.

Copies: every RAC complex goes on a ribosome (assumption), in polysomes and free 80S alike in proportion.
"""
import os, json
import numpy as np
import gemmi
import extend_complexes as xc
import structlib as sl
import build_structs as bs
from config import LIB
from config import seed as _seed

RIBO = "ribosome_80S"
LIB_ID = "ribosome_80S_RAC"
RAC_ING = "CPX_1743_Ribosome_associated_complex"
CT = ("7Z3N", {"C": "ZUO1", "D": "SSZ1"})          # C. thermophilum 80S-RAC: chain -> yeast gene
SC = ("7X3K", {"A": "ZUO1", "B": "SSZ1"})
INFO_F = f"{LIB}/{LIB_ID}.json"


def _kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0); U, S, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T)); R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


def _trimmed(P, Q, n=5):
    R, t = _kabsch(P, Q)
    for _ in range(n):
        r = np.linalg.norm(P @ R.T + t - Q, axis=1); k = r < max(0.3, np.percentile(r, 70)); R, t = _kabsch(P[k], Q[k])
    r = np.linalg.norm(P @ R.T + t - Q, axis=1)
    return R, t, float(np.median(r)), float((r < 0.5).mean())


def _chains(st, atom, rna):
    """chain -> (sequence, coordinates of `atom`) for RNA or protein chains"""
    out = {}
    for ch in st[0]:
        seq, X = [], []
        for r in ch:
            info = gemmi.find_tabulated_residue(r.name)
            if info is None or (info.is_nucleic_acid() if rna else info.is_amino_acid()) is False:
                continue
            a = r.find_atom(atom, "*")
            if a is None:
                continue
            seq.append(info.one_letter_code.upper()); X.append(a.pos.tolist())
        if len(seq) > 20:
            out[ch.name] = ("".join(seq), np.array(X) / 10)
    return out


def _pairs(s1, X1, s2, X2):
    """coordinate pairs of aligned identical residues"""
    al = gemmi.align_string_sequences(list(s1), list(s2), [])
    i = j = 0; A, B = [], []
    for n, op in al.cigar_tuples() if hasattr(al, "cigar_tuples") else _cigar(al.cigar_str()):
        for _ in range(n):
            if op == "M":
                if s1[i] == s2[j]:
                    A.append(X1[i]); B.append(X2[j])
                i += 1; j += 1
            elif op == "I":
                i += 1
            else:
                j += 1
    return np.array(A), np.array(B), al.calculate_identity()


def _cigar(s):
    import re
    return [(int(n), op) for n, op in re.findall(r"(\d+)([MID])", s)]


def build(log=print):
    if os.path.exists(INFO_F) and os.path.exists(f"{LIB}/{LIB_ID}.npz"):
        return json.load(open(INFO_F))
    st4 = gemmi.read_structure(xc.fetch_assembly("4V88")); st4.remove_ligands_and_waters()
    xs, ks, _ = sl.beads_from_structure(st4)
    lib = np.load(f"{LIB}/{RIBO}.npz"); off = xs.mean(0)
    assert np.allclose(lib["beads"], xs - off, atol=1e-3)
    stc = gemmi.read_structure(xc.fetch_assembly(CT[0])); stc.remove_ligands_and_waters()
    # 7Z3N rRNA -> 4V88 rRNA (each 7Z3N RNA chain to its best-matching 4V88 chain)
    r4, rc = _chains(st4, "P", True), _chains(stc, "P", True)
    P, Q = [], []
    for c, (s, X) in rc.items():
        best = max(((_pairs(s, X, s4, X4), k) for k, (s4, X4) in r4.items() if abs(len(s4) - len(s)) < 0.3 * len(s)),
                   key=lambda z: len(z[0][0]), default=None)
        if best is None:
            continue
        (A, B, ident), k = best
        P.append(A); Q.append(B)
        log(f"  7Z3N rRNA {c} ({len(s)} nt) ~ 4V88 {k}: {ident:.0f}% identity, {len(A)} identical aligned")
    R1, t1, d1, f1 = _trimmed(np.vstack(P), np.vstack(Q))
    log(f"  7Z3N -> 4V88: {sum(map(len, P))} phosphates, median deviation {d1:.2f} nm ({f1:.0%} within 0.5 nm)")
    # RAC beads: the ribosome-bound C. thermophilum RAC of 7Z3N itself (Zuo1 and Ssz1 47-58% identical to the yeast
    # proteins). The yeast RAC of 7X3K, superposed on it, deviates by ~2.4 nm (another state, 6 A, deposited without
    # its ribosome) and would float off the ribosome, so it only confirms the homology.
    sub = gemmi.Structure(); sub.add_model(gemmi.Model("1"))
    for ch in stc[0]:
        if ch.name in CT[1]:
            sub[0].add_chain(ch)
    xr, kr, _ = sl.beads_from_structure(sub)
    xr = (xr @ R1.T + t1) - off
    from scipy.spatial import cKDTree
    xb = lib["beads"].astype(np.float64)
    dmin = cKDTree(xb).query(xr)[0]
    keep = dmin > 0.5
    X = np.r_[xb, xr[keep]]; K = np.r_[lib["kind"], kr[keep]]
    bs.save(LIB_ID, X, K, extra=dict(pdb="7Z3N"), per=30)
    info = dict(rrna_dev=d1, rac_beads=int(keep.sum()), dropped=int((~keep).sum()),
                rac_contact_frac=float((dmin[keep] < 1.0).mean()))
    json.dump(info, open(INFO_F, "w"))
    log(f"  {LIB_ID}: RAC {keep.sum()} beads ({(~keep).sum()} inside the ribosome dropped), "
        f"{info['rac_contact_frac']:.0%} within 1 nm of ribosome beads")
    return info


def rac_spheres(pos, qu, typ, lp):
    """sphere tree of the RAC part of every decorated ribosome (types `typ`, ribosome frame), as obstacles for the
    mRNA threads that are placed after the decoration"""
    from scipy.spatial.transform import Rotation
    ids = list(lp["ids"])
    if LIB_ID not in ids:
        return np.zeros((0, 3)), np.zeros(0)
    k = np.asarray(typ) == ids.index(LIB_ID)
    z = np.load(f"{LIB}/{LIB_ID}.npz"); n0 = len(np.load(f"{LIB}/{RIBO}.npz")["beads"])
    rc, rr = sl.sphere_tree(z["beads"][n0:].astype(np.float64), per=30)
    q = np.asarray(qu, np.float64)[k]
    Rm = Rotation.from_quat(np.c_[q[:, 1:], q[:, :1]]).as_matrix()
    C = (np.einsum("nij,sj->nsi", Rm, rc) + np.asarray(pos, np.float64)[k][:, None, :]).reshape(-1, 3)
    return C, np.tile(rr, int(k.sum()))


def decorate(cell, lp, poly, n_target, obstacles=None, seed=71, log=print):
    """turn up to n_target polysome ribosomes into ribosome_80S_RAC where the RAC fits (its own sphere tree against
    the neighbouring ribosomes, the obstacles and the compartment); returns per-ribosome library types"""
    import pack3d as p3, cell3d as c3
    ids = list(lp["ids"]); li0, li1 = ids.index(RIBO), ids.index(LIB_ID)
    z = np.load(f"{LIB}/{LIB_ID}.npz"); n0 = len(np.load(f"{LIB}/{RIBO}.npz")["beads"])
    rc, rr = sl.sphere_tree(z["beads"][n0:].astype(np.float64), per=30)
    pos, qu = poly["poly_pos"].astype(np.float64), poly["poly_q"].astype(np.float64)
    typ = np.full(len(pos), li0, np.int32)
    G = cell.G; glo = np.full(3, -G[0] - 10.0); ghi = -glo
    # the RAC touches its own ribosome: test it against the obstacles (grid) and the other ribosomes (explicitly)
    s0, s1 = lp["st_off"][li0], lp["st_off"][li0 + 1]; SC, SR = lp["STC"][s0:s1], lp["STR"][s0:s1]
    OC, OR = obstacles if obstacles is not None else (np.zeros((0, 3)), np.zeros(0))
    cap = len(OR) + int(n_target) * len(rr) + 10000
    Gr = p3.grid_new(glo, ghi, 20.0, cap, cap * 14)
    p3.g_insert(glo, 20.0, *Gr, np.asarray(OC, np.float64), np.asarray(OR, np.float64))
    from scipy.spatial import cKDTree
    tree = cKDTree(pos); reach = float(lp["BR"][li0]) + float(np.linalg.norm(rc, axis=1).max() + rr.max())
    a = cell.args(glo, ghi, margin=50.0); out = np.zeros(c3.NMEM); M = np.empty((3, 3)); M2 = np.empty((3, 3))
    seed = _seed(seed); rng = np.random.default_rng(seed); done = 0
    for i in rng.permutation(len(pos)):
        if done >= n_target:
            break
        p3.quat_mat(qu[i], M); P = np.empty((len(rc), 3)); p3.transform(rc, M, pos[i], P)
        if p3.g_overlap(glo, 20.0, *Gr[:6], P, rr, p3.TOL):
            continue
        hit = False
        for j in tree.query_ball_point(pos[i], reach):
            if j == i:
                continue
            p3.quat_mat(qu[j], M2); Q = np.empty((len(SC), 3)); p3.transform(SC, M2, pos[j], Q)
            if (np.linalg.norm(P[:, None] - Q[None], axis=2) < rr[:, None] + SR[None] - p3.TOL).any():
                hit = True; break
        if hit:
            continue
        if not p3.check_volume(P, rr, P[:, 0].mean(), P[:, 1].mean(), P[:, 2].mean(), float(rr.max()) + 10.0, 1 << c3.CYTO, *a, out):
            continue
        p3.g_insert(glo, 20.0, *Gr, P, rr); typ[i] = li1; done += 1
    log(f"RAC on polysome ribosomes: {done} of {int(n_target)} wanted ({len(pos)} ribosomes)")
    return typ
