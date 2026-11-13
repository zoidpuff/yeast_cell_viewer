"""Eisosomes (literature item L7): Pil1/Lsp1 BAR-domain coats lining plasma-membrane furrows (MCC).

Furrows: ~30 per cell, about 300 nm long and 50 nm deep (Stradalova et al. 2009 J Cell Sci 122:2887; Karotki et al.
2011 J Cell Biol 195:889), with a half-pipe coat in vivo whose diameter is 30-60 nm in cryo-ET of vitreous sections
(Bharat et al. 2018 Structure 26:879). The furrows themselves are cell3d T_SLOT primitives (placed with the cell
geometry); this module builds the coat.

Coat: the near-native Pil1 (8QB7) and Lsp1 (8QB8) lattices of Kefauver et al. 2024 Nature 632:664, solved on eisosome
tubules isolated from yeast plasma membranes (coat diameter 31 nm, inside the measured in-situ range). The seven dimers
of 8QB7 are related by screw operations about the tube axis (0.03 nm rms), which give the lattice: two basis steps
(rotation about the axis, shift along it). The amphipathic helix (residues 39-48) sits in the cytoplasmic leaflet,
AH_DEPTH above the bilayer midplane, which sets the membrane radius of the furrow bottom (cell3d.EIS_R). The coat
covers the half-cylinder bottom of each furrow with whole dimers (the lattice cannot continue onto the flat walls
without bending the dimers). Each lattice site holds a Pil1 or an Lsp1 dimer in proportion to their pools (they
co-assemble; the two lattices are nearly identical).

Library frame of a dimer: x along the tube axis, z the outward radial direction through the dimer centroid (cytoplasmic
side, +z as for all membrane proteins), origin on the membrane midplane under the centroid.
"""
import os, json
import numpy as np
import gemmi
import extend_complexes as xc
import pack3d as p3
import cell3d as c3
import build_structs as bs
from config import LIB
from config import seed as _seed

LIB_IDS = ("Pil1_dimer", "Lsp1_dimer")
PDB = {"Pil1_dimer": "8QB7", "Lsp1_dimer": "8QB8"}
GENES = {"Pil1_dimer": "PIL1", "Lsp1_dimer": "LSP1"}
AH = (39, 48)            # amphipathic helix (Kefauver et al. 2024)
AH_DEPTH = 1.3           # nm from the bilayer midplane to the AH axis (interfacial helix; assumption)
H = 20.0
GEOM_F = f"{LIB}/eisosome_lattice.json"


def _ca(pid):
    st = gemmi.read_structure(xc.fetch_assembly(pid)); st.setup_entities()
    return {ch.name: (np.array([[r.find_atom("CA", "*").pos.x, r.find_atom("CA", "*").pos.y, r.find_atom("CA", "*").pos.z]
                                for r in ch if r.find_atom("CA", "*")]) / 10,
                      np.array([r.seqid.num for r in ch if r.find_atom("CA", "*")])) for ch in st[0]}


def _kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0); U, S, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T)); R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


def _rot(axis, ang):
    from scipy.spatial.transform import Rotation
    return Rotation.from_rotvec(np.asarray(axis) * ang).as_matrix()


def _lattice(ca):
    """dimers of the patch, screw operations from the first one, tube axis (direction, point)"""
    from scipy.spatial.transform import Rotation
    names = list(ca)
    X = {k: v[0] for k, v in ca.items()}
    # dimers: chain pairs with the most CA-CA contacts
    from scipy.spatial import cKDTree
    pairs, used = [], set()
    for a in names:
        if a in used:
            continue
        ta = cKDTree(X[a])
        b = max((n for n in names if n != a and n not in used), key=lambda n: sum(map(len, ta.query_ball_point(X[n], 0.8))))
        pairs.append((a, b)); used |= {a, b}
    X0 = np.vstack([X[pairs[0][0]], X[pairs[0][1]]])
    ops = []
    for a, b in pairs[1:]:
        best = None
        for Y in (np.vstack([X[a], X[b]]), np.vstack([X[b], X[a]])):     # the screw, not screw x dimer 2-fold
            R, t = _kabsch(X0, Y)
            rv = Rotation.from_matrix(R).as_rotvec(); ang = np.linalg.norm(rv)
            rms = np.sqrt(((X0 @ R.T + t - Y) ** 2).sum(1).mean())
            if best is None or ang < best[3]:
                best = (R, t, rms, ang, rv / max(ang, 1e-9))
        ops.append(best)
    axes = np.array([o[4] * np.sign(o[4] @ ops[0][4]) for o in ops])
    d = axes.mean(0); d /= np.linalg.norm(d)
    M = np.vstack([np.eye(3) - o[0] for o in ops]); y = np.concatenate([o[1] - (o[1] @ d) * d for o in ops])
    c = np.linalg.lstsq(M, y, rcond=None)[0]; c -= (c @ d) * d
    steps = [(float(np.degrees(o[3]) * np.sign(o[4] @ d)), float(o[1] @ d), float(o[2])) for o in ops]
    return pairs, steps, d, c


def model(log=print):
    """library entries Pil1_dimer, Lsp1_dimer and the lattice (eisosome_lattice.json)"""
    if os.path.exists(GEOM_F) and all(os.path.exists(f"{LIB}/{k}.npz") for k in LIB_IDS):
        return json.load(open(GEOM_F))
    geo = None
    for lid in LIB_IDS:
        ca = _ca(PDB[lid])
        pairs, steps, d, c = _lattice(ca)
        a, b = pairs[0]
        P = np.vstack([ca[a][0], ca[b][0]]); num = np.r_[ca[a][1], ca[b][1]]
        g = P.mean(0) - c; s0 = g @ d; rho = g - s0 * d; rho /= np.linalg.norm(rho)
        rad = lambda Q: np.linalg.norm((Q - c) - np.outer((Q - c) @ d, d), axis=1)
        ah = (num >= AH[0]) & (num <= AH[1])
        r_mem = float(rad(P[ah]).mean()) - AH_DEPTH
        o = c + s0 * d + r_mem * rho
        Rf = np.c_[d, np.cross(rho, d), rho]                   # columns: x (axis), y, z (radial)
        B = (P - o) @ Rf
        tm = np.hypot(B[:, 1], B[:, 2] + r_mem) - r_mem < c3.HALF  # residues inside the cytoplasmic leaflet (radius from the axis)
        bs.save(lid, B, np.zeros(len(B), np.int8), extra=dict(src="PDB", pdb=PDB[lid], membrane="TM", tm=tm), per=30)
        ang = np.degrees(np.arctan2(B[:, 1], B[:, 2] + r_mem))   # angular extent about the tube axis
        info = dict(r_mem=r_mem, half_angle=float(np.abs(ang).max()), half_len=float(np.abs(B[:, 0]).max()),
                    r_out=float(rad(P).max()), screw_rms=max(s[2] for s in steps), n_tm=int(tm.sum()))
        log(f"  {lid} ({PDB[lid]}): membrane midplane radius {r_mem:.2f} nm (AH at {r_mem + AH_DEPTH:.2f}), coat out to "
            f"{info['r_out']:.1f} nm, dimer spans +-{info['half_angle']:.0f} deg; screw operations fit to {info['screw_rms']:.3f} nm")
        if geo is None:
            # lattice basis: the two shortest independent steps
            arc = lambda s: np.array([np.radians(s[0]) * (r_mem + 3.0), s[1]])
            st = sorted(steps, key=lambda s: np.linalg.norm(arc(s)))
            a_ = st[0]
            b_ = next(s for s in st[1:] if abs(np.cross(arc(a_), arc(s))) > 0.3 * np.linalg.norm(arc(a_)) * np.linalg.norm(arc(s)))
            area = abs(np.radians(a_[0]) * b_[1] - a_[1] * np.radians(b_[0])) * (r_mem + 3.0)
            geo = dict(info, a=[a_[0], a_[1]], b=[b_[0], b_[1]], area_nm2=float(area), pdb=PDB[lid])
            log(f"  lattice: steps ({a_[0]:.2f} deg, {a_[1]:+.2f} nm) and ({b_[0]:.2f} deg, {b_[1]:+.2f} nm); "
                f"{area:.1f} nm2 per dimer at the coat mid-radius")
    assert abs(geo["r_mem"] - c3.EIS_R) < 0.3, f"cell3d.EIS_R {c3.EIS_R} vs coat membrane radius {geo['r_mem']:.2f}"
    json.dump(geo, open(GEOM_F, "w"))
    return geo


def place(cell, lp, n_pil1, n_lsp1, obstacles=None, seed=53, log=print):
    """coat dimers on the half-cylinder bottoms of the cell's eisosome furrows, at most n_pil1 / n_lsp1 dimers"""
    from scipy.spatial.transform import Rotation
    geo = model(log=log)
    seed = _seed(seed); rng = np.random.default_rng(seed)
    li = {k: int(np.where(lp["ids"] == k)[0][0]) for k in LIB_IDS}
    sub = {k: (lp["STC"][lp["st_off"][v]:lp["st_off"][v + 1]], lp["STR"][lp["st_off"][v]:lp["st_off"][v + 1]]) for k, v in li.items()}
    G = cell.G
    glo = np.full(3, -G[0] - 10.0); ghi = -glo
    slots = cell.furrows()
    cap = 200000 + (len(obstacles[1]) if obstacles is not None else 0)
    G1 = p3.grid_new(glo, ghi, H, cap, cap * 14)
    if obstacles is not None and len(obstacles[1]):
        k = np.all((obstacles[0] > glo) & (obstacles[0] < ghi), axis=1)
        p3.g_insert(glo, H, *G1, obstacles[0][k], obstacles[1][k])
    pa, pb = np.array(geo["a"]), np.array(geo["b"])
    phi_max = 90.0 - geo["half_angle"]
    left = {"Pil1_dimer": int(n_pil1), "Lsp1_dimer": int(n_lsp1)}
    M = np.empty((3, 3))
    pos, rot, typ, fid = [], [], [], []
    blocked = 0
    for j, i in enumerate(slots):
        q = cell.PQ[i]; a, b = q[:3], q[3:6]; r = q[6]
        L = np.linalg.norm(b - a); e1 = (b - a) / L; down = -(a + b) / np.linalg.norm(a + b)
        f0 = rng.uniform(0, abs(pa[0])); s0 = rng.uniform(0, abs(pb[1]) + abs(pa[1]))
        mm, nn = np.meshgrid(np.arange(-150, 151), np.arange(-60, 61), indexing="ij")
        phi = f0 + mm * pa[0] + nn * pb[0]; s = s0 + mm * pa[1] + nn * pb[1]
        ok = (np.abs(phi) <= phi_max) & (s >= geo["half_len"]) & (s <= L - geo["half_len"])
        for ph, sv in zip(phi[ok], s[ok]):
            tot = left["Pil1_dimer"] + left["Lsp1_dimer"]
            if tot <= 0:
                break
            lid = "Pil1_dimer" if rng.random() < left["Pil1_dimer"] / tot else "Lsp1_dimer"
            z = _rot(e1, np.radians(ph)) @ down
            Rm = np.c_[e1, np.cross(z, e1), z]
            o = a + sv * e1 + r * z
            qq = Rotation.from_matrix(Rm).as_quat(); qq = np.array([qq[3], qq[0], qq[1], qq[2]])
            stc, str_ = sub[lid]
            P = np.empty((len(stc), 3)); p3.quat_mat(qq, M); p3.transform(stc, M, o, P)
            if p3.g_overlap(glo, H, *G1[:6], P, str_, p3.TOL):
                blocked += 1
                continue
            pos.append(o); rot.append(qq); typ.append(li[lid]); fid.append(j); left[lid] -= 1
    log(f"eisosomes: {len(slots)} furrows, {len(pos)} coat dimers ({(np.array(typ) == li['Pil1_dimer']).sum()} Pil1, "
        f"{(np.array(typ) == li['Lsp1_dimer']).sum()} Lsp1; {len(pos) / max(len(slots), 1):.0f} per furrow, {blocked} sites blocked)")
    used = {k: int(v) - left[k] for k, v in (("Pil1_dimer", n_pil1), ("Lsp1_dimer", n_lsp1))}
    return dict(eis_pos=np.array(pos, np.float32).reshape(-1, 3), eis_q=np.array(rot, np.float32).reshape(-1, 4),
                eis_type=np.array(typ, np.int32), eis_id=np.array(fid, np.int32)), used
