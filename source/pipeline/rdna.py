"""rDNA in the nucleolus: actively transcribed rRNA genes ("Christmas trees") with RNA polymerase I and nascent pre-rRNA.

The rDNA is a tandem array of about 150 repeats of 9.1 kb on chromosome XII. In growing cells about half of the
repeats are active and largely free of nucleosomes (psoralen accessibility: Dammann et al. 1995; Hmo1-bound chromatin:
Merz et al. 2008; later ChIP work argues for unphased nucleosomes instead). Miller spreads give about 50 +- 20 Pol I
per active gene on the 6.7 kb 35S transcription unit and about 60 nt/s elongation (French et al. 2003); 75 genes x 50
polymerases then finish about 2,000 transcripts per minute, the ribosome production rate. The small-subunit processome (90S) assembles on the nascent
transcript; about 70% of transcripts are cleaved at site A2 while still being transcribed (Kos & Tollervey 2010),
which releases the 90S and leaves the large-subunit part on the polymerase.

Model:
  - Active genes: a persistent self-avoiding chain of Pol I sites inside the nucleolus, one site per 134 bp of DNA
    (6.7 kb / 50); the DNA (2 bp per bead) runs straight between sites and through each polymerase's cleft (Pol I
    oriented with the DNA axis of the elongation complex 5M3F mapped onto the library's 4C2M).
  - Nascent pre-rRNA from each polymerase's RNA exit channel, length set by its position along the gene; threaded
    self-avoiding at 2 nt per bead (0.6 nm per bead, about the rise of a base-paired transcript).
      * uncleaved transcripts beyond A2: the 5' end carries a 90S particle (5'ETS + 18S, 2.5 kb, inside the 90S);
      * cleaved transcripts: only the 3' part beyond A2 remains;
      * shorter transcripts: the whole transcript.
  - Inactive repeats and the non-transcribed spacers of active repeats are nucleosomal chromatin in the nucleolus
    (placed by s3 with the chromatin walk; their nucleosome count is nucleosomes_in_rdna()).
Polymerases and 90S particles come out of the PolI and preribo_90S pools.
"""
import numpy as np
import gemmi
import extend_complexes as xc
import assemblies as A
import membrane_orient as mo
import pack3d as p3
import cell3d as c3
from config import LIB
from config import seed as _seed

N_REPEATS = 150
ACTIVE_FRACTION = 0.5
REPEAT_BP, GENE_BP = 9100, 6700
POL_PER_GENE = 50
BP_NM = 0.34
SITE_STEP = GENE_BP / POL_PER_GENE * BP_NM      # 45.6 nm of DNA between neighbouring polymerases
DNA_STEP = 0.68                                  # 2 bp per bead (as linker DNA)
POL_R = 9.0                                      # exclusion radius of a polymerase site for the gene walk
A2_NT, SSU_NT, CO_CLEAVED = 2700, 2500, 0.7
NT_PER_BEAD, RNA_STEP, RNA_R = 2, 0.6, 0.7       # 0.3 nm per nt: about the rise of base-paired RNA
RNA_PERSIST = 0.8                                # lower values jam the walks in the crowded nucleolus
RNA_REACH = 200.0                                # nm: a nascent transcript stays near its polymerase
NRL_BP = 165


def n_active():
    return int(round(N_REPEATS * ACTIVE_FRACTION))


def nucleosomes_in_rdna():
    """nucleosomes on the inactive repeats and on the spacers of the active ones"""
    na = n_active()
    return int(round(((N_REPEATS - na) * REPEAT_BP + na * (REPEAT_BP - GENE_BP)) / NRL_BP))


def _kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0); U, S, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T)); R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


def pol_frame(log=print, lib_id="PolI", lib_pdb="4C2M", ec_pdb="5M3F"):
    """in the library frame of polymerase lib_id (built from lib_pdb): point on the DNA axis in the cleft, DNA axis
    (downstream), RNA exit point and direction, from the elongation complex ec_pdb"""
    def beads(pid):
        st = gemmi.read_structure(xc.fetch_assembly(pid)); st.remove_ligands_and_waters(); st.setup_entities()
        xyz, kind, acc, pos = mo.beads_with_positions(st, align=A.entity_align(pid))
        ptype = []
        for ch in st[0]:
            ent = st.get_entity_of(ch.get_polymer())
            for res in ch:
                info = gemmi.find_tabulated_residue(res.name)
                if info is None:
                    continue
                if (info.is_amino_acid() and res.find_atom("CA", "*")) or (info.is_nucleic_acid() and (res.find_atom("P", "*") or res.find_atom("C4'", "*"))):
                    ptype.append("RNA" if ent is not None and ent.polymer_type == gemmi.PolymerType.Rna else
                                 "DNA" if info.is_nucleic_acid() else "P")
        return xyz, acc, pos, np.array(ptype)
    lib = np.load(f"{LIB}/{lib_id}.npz")["beads"].astype(np.float64)
    x0, a0, p0, _ = beads(lib_pdb)
    assert len(x0) == len(lib) and np.allclose(x0 - x0.mean(0), lib, atol=1e-3)
    x0 = x0 - x0.mean(0)
    x1, a1, p1, t1 = beads(ec_pdb)
    key = {}
    for i, k in enumerate(zip(a0, p0)):
        if k[0] and k[1] > 0:
            key.setdefault(k, i)
    I, J = [], []
    for j, k in enumerate(zip(a1, p1)):
        if t1[j] == "P" and k in key:
            I.append(key[k]); J.append(j)
    R, t = _kabsch(x1[J], x0[I])
    # 4C2M's two largest subunits have no UniProt mapping and some mapped residues are numbered offset, so the
    # residue matches only seed the fit; it is refined on nearest protein residues (ICP)
    from scipy.spatial import cKDTree
    P1, tree = x1[t1 == "P"], cKDTree(x0)
    for _ in range(8):
        d, j = tree.query(P1 @ R.T + t)
        k = d < max(0.6, np.percentile(d, 70))
        R, t = _kabsch(P1[k], x0[j[k]])
    dev = tree.query(P1 @ R.T + t)[0]
    X = x1 @ R.T + t
    dna, rna = X[t1 == "DNA"], X[t1 == "RNA"]
    c = dna.mean(0); w, V = np.linalg.eigh(np.cov((dna - c).T)); u = V[:, 2]
    # RNA in file order runs 5' -> 3'; the 3' end sits at the active site, the 5' end at the exit channel
    e, r3 = rna[0], rna[-1]
    edir = e - r3; edir /= np.linalg.norm(edir)
    # downstream DNA lies ahead of the active site, away from the exit channel
    if np.dot(u, r3 - e) < 0:
        u = -u
    c = c + u * np.dot(r3 - c, u)                 # the axis point level with the active site
    log(f"  {lib_id} frame from {ec_pdb} on {lib_pdb}: {len(I)} seed residues, median nearest-residue distance {np.median(dev):.2f} nm; "
        f"{len(dna)} DNA and {len(rna)} RNA nucleotides; exit {np.linalg.norm(e - c):.1f} nm from the DNA axis")
    return c, u, e, edir


def _rot_u_to(u, t, roll):
    """rotation matrix taking unit u to unit t, then rolled by `roll` about t"""
    v = np.cross(u, t); s = np.linalg.norm(v); cth = float(np.dot(u, t))
    if s < 1e-9:              # parallel: identity, or a half turn about any axis perpendicular to u
        p = np.cross(u, [1.0, 0, 0] if abs(u[0]) < 0.9 else [0, 1.0, 0]); p /= np.linalg.norm(p)
        R = np.eye(3) if cth > 0 else 2 * np.outer(p, p) - np.eye(3)
    else:
        k = v / s; K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
        R = np.eye(3) + s * K + (1 - cth) * K @ K
    K = np.array([[0, -t[2], t[1]], [t[2], 0, -t[0]], [-t[1], t[0], 0]])
    Rr = np.eye(3) + np.sin(roll) * K + (1 - np.cos(roll)) * K @ K
    return Rr @ R


def _quat(R):
    from scipy.spatial.transform import Rotation
    x, y, z, w = Rotation.from_matrix(R).as_quat()
    return np.array([w, x, y, z])


def place(cell, lp, n_pol, n_90s, obstacles, seed=41, log=print):
    """Christmas trees in the nucleolus. n_pol, n_90s: available Pol I and 90S particles (nucleolar pools).
    obstacles: (centres, radii) of everything placed before. Returns pre-pass arrays:
    rdna_pts/off (DNA polylines), prerrna_pts/off (nascent RNA), rpol_pos/q + li_rpol, r90_pos/q + li_r90."""
    import mrna
    seed = _seed(seed); rng = np.random.default_rng(seed)
    li_pol = int(np.where(lp["ids"] == "PolI")[0][0]); li_90 = int(np.where(lp["ids"] == "preribo_90S")[0][0])
    c_dna, u_dna, e_exit, d_exit = pol_frame(log=log)
    G = cell.G
    cn = G[3:6]; rn = float(G[11:14].max())
    glo, ghi = cn - rn - 100, cn + rn + 100
    a = cell.args(glo, ghi, margin=100.0)
    OC, OR = obstacles
    k = np.all((OC > glo) & (OC < ghi), axis=1); OC, OR = OC[k].astype(np.float64), OR[k].astype(np.float64)
    nucl = 1 << c3.NUCL
    # random starts inside the nucleolus
    box = rng.uniform(glo, ghi, size=(20000, 3))
    lab = cell.label(box)
    pool = box[lab == c3.NUCL]
    ng = min(n_active(), int(n_pol // POL_PER_GENE))
    starts = pool[rng.choice(len(pool), ng, replace=False)]
    dirs = rng.normal(size=(ng, 3)); dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
    lens = np.full(ng, POL_PER_GENE + 2, np.int64)        # + promoter and terminator flanks
    cap_s = int(lens.sum()) + len(OR) + 1000
    sites, soff = mrna.thread_all(seed, glo, ghi, *a, OC, OR, starts, dirs, lens, 0.0, cap_s, int(cap_s * 40 + 8 * len(OR)),
                                  SITE_STEP, POL_R, 1.4, nucl, 0.0, 3)
    # polymerases on the interior sites, DNA through all sites
    pol_pos, pol_q, pol_len, dna, rna_start, rna_dir = [], [], [], [], [], []
    for g in range(ng):
        S = sites[soff[g]:soff[g + 1]].astype(np.float64)
        if len(S) < 3:
            continue
        # DNA: a Catmull-Rom spline through the sites (the polymerases sit on it, oriented along its tangent)
        m = max(2, int(np.ceil(SITE_STEP / DNA_STEP)))
        Sx = np.vstack([2 * S[0] - S[1], S, 2 * S[-1] - S[-2]])
        s = np.linspace(0, 1, m + 1)[1:, None]
        pts = [S[0]]
        for i in range(1, len(Sx) - 2):
            p0, p1, p2, p3_ = Sx[i - 1], Sx[i], Sx[i + 1], Sx[i + 2]
            pts.append(0.5 * (2 * p1 + (p2 - p0) * s + (2 * p0 - 5 * p1 + 4 * p2 - p3_) * s ** 2 + (3 * p1 - p0 - 3 * p2 + p3_) * s ** 3))
        dna.append(np.vstack(pts))
        for i in range(1, len(S) - 1):
            t = S[i + 1] - S[i - 1]; t /= np.linalg.norm(t)
            R = _rot_u_to(u_dna, t, rng.uniform(0, 2 * np.pi))
            pol_pos.append(S[i] - R @ c_dna); pol_q.append(_quat(R))
            pol_len.append(int(round((i - 0.5) * GENE_BP / POL_PER_GENE)))      # nt transcribed so far
            rna_start.append(S[i] - R @ c_dna + R @ e_exit); rna_dir.append(R @ d_exit)
    pol_pos = np.array(pol_pos); pol_q = np.array(pol_q); pol_len = np.array(pol_len)
    # what each transcript looks like
    cleaved = (pol_len > A2_NT) & (rng.random(len(pol_len)) < CO_CLEAVED)
    knob = (pol_len > A2_NT) & ~cleaved
    rna_nt = np.where(cleaved, pol_len - A2_NT, np.where(knob, pol_len - SSU_NT, pol_len))
    rna_beads = np.maximum(2, rna_nt // NT_PER_BEAD).astype(np.int64)
    # obstacles for the transcripts: polymerases, DNA and everything before
    PC, PR = p3.emit_spheres(np.full(len(pol_pos), li_pol, np.int32), pol_pos.astype(np.float32), pol_q.astype(np.float32),
                             lp["st_off"], lp["STC"], lp["STR"], lp["BR"], glo, ghi)
    D = np.vstack(dna) if dna else np.zeros((0, 3))
    OC2 = np.concatenate([OC, PC, D]); OR2 = np.concatenate([OR, PR, np.full(len(D), 1.0)])
    rs, rdir = np.array(rna_start), np.array(rna_dir)

    def thread(idx, OCx, ORx, sd):
        if not len(idx):
            return np.zeros((0, 3), np.float32), np.zeros(1, np.int64)
        cap_s = int(rna_beads[idx].sum()) + len(ORx) + 1000
        return mrna.thread_all(sd, glo, ghi, *a, OCx, ORx, rs[idx], rdir[idx], rna_beads[idx], 0.0, cap_s,
                               int(cap_s * 6 + 8 * len(ORx)), RNA_STEP, RNA_R, RNA_PERSIST, nucl, RNA_REACH, 25)
    # pass 1: transcripts that carry a 90S, each followed by its particle at the 5' end. The last END beads are left
    # out of the collision check: the 5' end of the pre-rRNA lies inside its 90S, and beads that end up inside the
    # particle are trimmed (they are part of it)
    s0, s1 = lp["st_off"][li_90], lp["st_off"][li_90 + 1]
    stc, strr = lp["STC"][s0:s1].astype(np.float64), lp["STR"][s0:s1].astype(np.float64); r90 = float(lp["BR"][li_90])
    END = 80
    kidx = rng.permutation(np.where(knob)[0]); oidx = rng.permutation(np.where(~knob)[0])
    rp, ro = thread(kidx, OC2, OR2, seed + 1)
    held = np.zeros(len(rp), bool)
    for j in range(len(kidx)):
        held[max(ro[j], ro[j + 1] - END):ro[j + 1]] = True
    dims, head, nxt, sid, cen, rad, cnt = p3.grid_new(glo, ghi, 10.0, len(OR2) + len(rp) + 200000, (len(OR2) + len(rp)) * 30 + 2000000)
    p3.g_insert(glo, 10.0, dims, head, nxt, sid, cen, rad, cnt, OC2, OR2)
    p3.g_insert(glo, 10.0, dims, head, nxt, sid, cen, rad, cnt, rp[~held].astype(np.float64), np.full(int((~held).sum()), RNA_R))
    k_pos, k_q, KC, KR = [], [], [], []
    keep1 = np.diff(ro).copy()
    out = np.zeros(c3.NMEM); br1 = np.full(1, RNA_R)
    for j in range(len(kidx)):
        seg = rp[ro[j]:ro[j + 1]].astype(np.float64)
        h0 = max(0, len(seg) - END)
        placed = None
        if len(seg) >= 2 and len(k_pos) < n_90s:
            end = seg[-1]; dirn = seg[-1] - seg[-2]; dirn /= np.linalg.norm(dirn)
            for attempt in range(40):
                v = dirn + 0.9 * rng.normal(size=3) * (attempt > 0); v /= np.linalg.norm(v)
                c = end + v * r90 * rng.uniform(0.35, 0.8)
                if (nucl >> c3.label_from(c[0], c[1], c[2], G, a[1], out, c3.mem_sdfs(c[0], c[1], c[2], *a, out))) & 1 == 0:
                    continue
                Rq = _rot_u_to(np.array([0, 0, 1.0]), v, rng.uniform(0, 2 * np.pi))
                Cs = stc @ Rq.T + c
                if p3.g_overlap(glo, 10.0, dims, head, nxt, sid, cen, rad, Cs, strr, 0.0):
                    continue
                p3.g_insert(glo, 10.0, dims, head, nxt, sid, cen, rad, cnt, Cs, strr)
                k_pos.append(c); k_q.append(_quat(Rq)); KC.append(Cs); KR.append(strr); placed = (Cs, strr); break
        # the held 5' end: beads inside the new particle, or clashing with anything placed meanwhile, end the thread
        n = len(seg)
        for b in range(h0, len(seg)):
            q = seg[b:b + 1]
            inside = placed is not None and bool(np.any(np.linalg.norm(placed[0] - q, axis=1) < placed[1] + RNA_R))
            if inside or p3.g_overlap(glo, 10.0, dims, head, nxt, sid, cen, rad, q, br1, 0.0):
                n = b; break
        if n > h0:
            p3.g_insert(glo, 10.0, dims, head, nxt, sid, cen, rad, cnt, seg[h0:n], np.full(n - h0, RNA_R))
        keep1[j] = n
    R1 = np.concatenate([rp[ro[j]:ro[j] + keep1[j]] for j in range(len(kidx))]) if len(kidx) else np.zeros((0, 3), np.float32)
    # pass 2: the other transcripts, around the first ones and their particles
    KC = np.concatenate(KC) if KC else np.zeros((0, 3)); KR = np.concatenate(KR) if KR else np.zeros(0)
    rp2, ro2 = thread(oidx, np.concatenate([OC2, R1.astype(np.float64), KC]), np.concatenate([OR2, np.full(len(R1), RNA_R), KR]), seed + 2)
    seg_of = {}
    for j, i in enumerate(kidx):
        seg_of[i] = rp[ro[j]:ro[j] + keep1[j]]
    for j, i in enumerate(oidx):
        seg_of[i] = rp2[ro2[j]:ro2[j + 1]]
    segs = [seg_of[i] for i in range(len(rna_beads))]
    got = np.array([len(x) for x in segs])
    log(f"rDNA: {ng} active genes ({len(dna)} placed), {len(pol_pos)} Pol I, {sum(len(x) for x in dna)} DNA beads; "
        f"{len(rna_beads)} nascent transcripts ({cleaved.sum()} cleaved at A2, {knob.sum()} uncleaved beyond A2 -> {len(k_pos)} with a 90S), "
        f"{got.sum() / 1e3:.0f}k of {rna_beads.sum() / 1e3:.0f}k RNA beads ({got.sum() / max(1, rna_beads.sum()):.0%})")
    doff = np.r_[0, np.cumsum([len(x) for x in dna])].astype(np.int64)
    rpts = np.concatenate(segs) if segs else np.zeros((0, 3), np.float32)
    roff = np.r_[0, np.cumsum(got)].astype(np.int64)
    return dict(rdna_pts=(D.astype(np.float32)), rdna_off=doff, prerrna_pts=rpts.astype(np.float32), prerrna_off=roff,
                rpol_pos=pol_pos.astype(np.float32), rpol_q=pol_q.astype(np.float32), li_rpol=li_pol,
                r90_pos=np.array(k_pos, np.float32).reshape(-1, 3), r90_q=np.array(k_q, np.float32).reshape(-1, 4), li_r90=li_90)
