"""Driver for the 3D whole-cell packing: library pack, expected counts per tile, NPC + chromatin pre-pass,
checkerboard-phased parallel tile packing."""
import os, re, json, glob, time
import numpy as np, pandas as pd
import numba as nb
import cell3d as c3, tiles3d as t3, pack3d as p3

from config import W3D as W, LIB, RAW, NPROC, DATA
from config import seed as _seed
CHROM_KB = [230, 813, 316, 1532, 577, 270, 1091, 563, 440, 746, 667, 1078, 924, 784, 1091, 948]  # chr I-XVI
# Rabl layout (literature item L13): centromeres at the spindle pole body (kinetochore MTs), each arm a walk to its
# telomere at the nuclear envelope. Centromere positions (S288C, kb), the rDNA array on XII-R
CEN_KB = [151, 238, 114, 450, 152, 148, 497, 106, 355, 436, 440, 151, 268, 628, 327, 556]
RDNA_KB = (451, 468)
# chromatin fibre (IMPROVEMENTS K3): worm-like backbone with the mass density and persistence length measured in living
# yeast (Bystricky et al. 2004 PNAS 101:16495: 110-150 bp/nm, 170-220 nm; mid values), nucleosomes in irregular
# clutches within FIBER_R of the backbone (assumption, sized so that the clutches can be packed: ~34% occupancy), each
# within LINK_MAX of the previous one (core radius 5.5 nm x 2 + a ~20 bp linker). FIBER = False: the earlier
# nucleosome-by-nucleosome walk (12.5 nm steps), too extended below ~50 kb.
FIBER = True
NUC_MAP = True            # fibre from the measured nucleosome map (K2); False: uniform NRL_BP repeat
NRL_BP = 165.0
FIBER_BP_PER_NM = 130.0
FIBER_LP = 195.0
FIBER_R = 60.0            # nucleosomes within this distance of the backbone: dispersed clutches (cryo-ET), R(s) inside the
                          # Arbona 2017 chain and at 100 kb inside Bystricky 2004; chosen among 9 variants with Hi-C in view
LINK_MAX = 18.0
TEL_CLUSTERS = 7               # telomeres gather in 6-8 foci at the nuclear envelope (Gotta et al. 1996)
TEL_REACH = 0.8               # telomere target within this share of the arm's fibre contour from its start
TEL_THETA = (35.0, 115.0)      # angle from the SPB axis of the telomere of the shortest / a >= 700 kb arm (assumption:
                               # the telomere-SPB distance grows with arm length; Therizols et al. 2010, Duan et al. 2010)
TEL_DEPTH = 25.0               # telomere targets this far inside the inner nuclear membrane
import inventory as _inv
N_TRNA = 3.3e6 * _inv.CELL_SHARE          # 3.3 million per population-average cell (assumption), about 11 per ribosome
N_MRNA = 3.0e4
MRNA_MEAN_NT = 1500


# ------------------------------------------------------------------ library
# Collision shapes of the packing. The library's sphere trees (about 30 residues per sphere; 40 for the curated
# complexes) are what the viewer's simulation and the contact statistics use. As shapes for packing they are fat:
# their union is 2.4 times the union of the residue beads for the most abundant matrix proteins, 1.9 times in the
# cytosol and 3.8 times in the nucleolus (2026-11-13), and that, not the algorithm, limited how much a crowded
# compartment held. Trees of PACK_PER residues per sphere are made from the beads for packing only (pack_trees).
# Six matrix-rich test tiles of build R: 87.9% of the matrix request placed with the library trees, 93.1% with 15
# residues per sphere (2.6 times the CPU), 94.9% with 8 (7 times); species under 80% placed: 13, 3, 9.
PACK_PER = 15           # residues (or nucleotides) per sphere; 0: the library's trees


def pack_trees(lib_ids, per, nproc=NPROC):
    """{lib id: (centres, radii)} at `per` residues per sphere, cached in W3D/pack_trees_<per>.npz and recomputed for
    entries whose library file changed (size or modification time)"""
    import structlib as sl
    from joblib import Parallel, delayed
    fn = f"{W}/pack_trees_{per}.npz"
    sig = {i: f"{os.path.getsize(f'{LIB}/{i}.npz')}:{os.path.getmtime(f'{LIB}/{i}.npz'):.0f}" for i in lib_ids}
    cache = {}                                       # id -> (signature, centres, radii), also for ids not asked for now
    if os.path.exists(fn):
        z = np.load(fn, allow_pickle=True); off = z["off"]
        cache = {str(i): (str(z["sig"][k]), z["C"][off[k]:off[k + 1]], z["R"][off[k]:off[k + 1]]) for k, i in enumerate(z["ids"])}
    todo = [i for i in lib_ids if i not in cache or cache[i][0] != sig[i]]

    def tree(i):
        c, r = sl.sphere_tree(np.load(f"{LIB}/{i}.npz")["beads"].astype(np.float64), per=per)
        return i, (np.asarray(c, np.float64), np.asarray(r, np.float64))
    if todo:
        for i, (c, r) in Parallel(n_jobs=min(nproc, 32))(delayed(tree)(i) for i in todo):
            cache[i] = (sig[i], c, r)
        ids = sorted(cache)
        off = np.r_[0, np.cumsum([len(cache[i][2]) for i in ids])]
        np.savez(fn, ids=np.array(ids), sig=np.array([cache[i][0] for i in ids]), off=off,
                 C=np.concatenate([cache[i][1] for i in ids]), R=np.concatenate([cache[i][2] for i in ids]))
    return {i: (cache[i][1], cache[i][2]) for i in lib_ids}


def build_libpack(lib_ids, per=None):
    per = PACK_PER if per is None else per
    fine = pack_trees([str(i) for i in lib_ids], per) if per else {}
    st_off = [0]; C, R, BR, MEM, NB = [], [], [], [], []
    for i in lib_ids:
        z = np.load(f"{LIB}/{i}.npz")
        c = z["st_c"].astype(np.float64); r = z["st_r"].astype(np.float64)
        if str(i) in fine:
            c, r = fine[str(i)]
        C.append(c); R.append(r); st_off.append(st_off[-1] + len(r))
        BR.append(float((np.linalg.norm(c, axis=1) + r).max()))
        MEM.append("membrane" in z.files); NB.append(len(z["beads"]))
    return dict(ids=np.array(lib_ids), st_off=np.array(st_off, np.int64), STC=np.concatenate(C), STR=np.concatenate(R),
                BR=np.array(BR), memb=np.array(MEM), nbeads=np.array(NB))


# ------------------------------------------------------------------ expected totals per (lib_id, region)
R = t3


SOLUBLE_ONLY = ("cytosol", "nucleoplasm", "nucleolus")
# Cvt cargo reaches the vacuole lumen by a selective-autophagy route, without a signal peptide (Lynch-Day & Klionsky
# 2010; Ape4: Yuga et al. 2011)
CVT_CARGO = {"APE1", "AMS1", "APE4"}
_SP = None


def _signal_and_location():
    """UniProt accession -> (has signal peptide, subcellular-location text)"""
    global _SP
    if _SP is None:
        f = pd.read_csv(f"{DATA}/validation/uniprot_features.tsv", sep="\t", dtype=str, usecols=["Entry", "Signal peptide"])
        u = pd.read_csv(f"{DATA}/uniprot_yeast.tsv", sep="\t", dtype=str,
                        usecols=["Entry", "Subcellular location [CC]", "Gene Names (primary)", "Gene Names (ordered locus)"])
        sig = dict(zip(f.Entry, f["Signal peptide"].fillna("").str.contains("SIGNAL")))
        loc = dict(zip(u.Entry, u["Subcellular location [CC]"].fillna("")))
        _SP = {a: (sig.get(a, False), loc.get(a, "")) for a in set(sig) | set(loc)}
        _SP["_gene"] = {}
        for a, g1, g2 in zip(u.Entry, u["Gene Names (primary)"].fillna(""), u["Gene Names (ordered locus)"].fillna("")):
            for g in (g1.split(";")[0].strip().upper(), g2.split(";")[0].strip().upper()):
                if g:
                    _SP["_gene"].setdefault(g, a)
    return _SP


def _accs(r):
    """UniProt accessions of an ingredient: its uniprot field, else its genes"""
    sp = _signal_and_location()
    a = re.findall(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}", str(r.uniprot))
    if not a:
        a = [sp["_gene"][g] for g in (x.strip().upper() for x in str(r.genes).split(";")) if g in sp["_gene"]]
    return a


def side_of(r, comp):
    """lumen or cytosolic face for a soluble protein (or complex) annotated to the ER, vacuole, Golgi/endosomes or
    peroxisomes. ER, vacuole and Golgi: the lumen is reached through the secretory pathway, so only proteins with a signal
    peptide (or annotated to the lumen, or Cvt cargo) are lumenal; the rest (V1 sector, EGO, TORC1, COPI, clathrin,
    Arf/Rab GTPases) sit on the cytosolic face. Peroxisomes import matrix proteins without signal peptides, so the
    matrix stays the default and only UniProt 'cytoplasmic side' proteins (Pex1/6 ...) go to the face.
    Complexes: by the majority of their subunits."""
    sp = _signal_and_location()
    info = [sp[a] for a in _accs(r) if a in sp]
    if comp == "peroxisome":
        face = [("ytoplasmic side" in l and "matrix" not in l.lower()) for _, l in info]
        return "face" if face and np.mean(face) > 0.5 else "lumen"
    genes = {g.strip() for g in str(r.genes).split(";")}
    if genes & CVT_CARGO and comp == "vacuole":
        return "lumen"
    lum = [(sg or "lumen" in l.lower()) for sg, l in info]
    return "lumen" if lum and np.mean(lum) > 0.5 else "face"


_GPI = None


def gpi_anchored(r):
    """a single protein with a UniProt GPI-anchor lipidation site"""
    global _GPI
    if _GPI is None:
        from config import DATA
        u = pd.read_csv(f"{DATA}/validation/uniprot_features.tsv", sep="\t", usecols=["Entry", "Lipidation"])
        _GPI = set(u.Entry[u.Lipidation.fillna("").str.contains("GPI-anchor")])
    a = _accs(r)
    return len(a) == 1 and a[0] in _GPI


def no_signal_peptide(r):
    sp = _signal_and_location()
    info = [sp[a][0] for a in _accs(r) if a in sp]
    return bool(info) and not any(info)
# outer nuclear membrane (half of the 18.3 um2 NE) / (cortical ER 51.8 um2 + ONM 9.1 um2), from tile_measure areas
ER_ON_NE = 0.15


def expected_totals(ing2, oligo_map):
    rows = []
    for r in ing2.itertuples():
        if r.ing_id in ("NPC", "nucleosome"):
            continue
        lib_id, n = r.ing_id, r.copies_per_cell
        if lib_id in oligo_map:
            lib_id, k = oligo_map[lib_id]; n = n / k
        if not os.path.exists(f"{LIB}/{lib_id}.npz"):
            rows.append((lib_id, r.ing_id, "no_structure", n)); continue
        ntm = float(r.n_tm) if r.n_tm == r.n_tm else 0.0
        tm = ntm > 0 or (isinstance(r.membrane, str) and r.membrane in ("PM", "NE", "VAC"))
        cw = json.loads(r.comp_weights)
        if not tm and cw.get("cell_wall", 0) > 0 and no_signal_peptide(r):
            # a soluble protein without a signal peptide cannot be secreted into the wall (e.g. Ssa1/2 from surface
            # proteomics): its wall share goes to its other compartments, or to the cytosol
            cw = {c: f for c, f in cw.items() if c != "cell_wall"} or {"cytosol": 1.0}
            z = sum(cw.values()); cw = {c: f / z for c, f in cw.items()}
        if not tm and gpi_anchored(r):
            # GPI-anchored proteins sit on the outer face of the plasma membrane or are cross-linked into the wall; the
            # peripheral rule had put their plasma-membrane share on the cytosolic face (64,000 copies of Gas1, Ecm33
            # and others; review comment 4). Their plasma-membrane and cytosol shares go to the wall region.
            w_ = sum(cw.pop(c, 0.0) for c in ("plasma_membrane", "cytosol"))
            if w_ > 0:
                cw["cell_wall"] = cw.get("cell_wall", 0.0) + w_
        if tm:   # an integral membrane protein listed in soluble and membrane compartments lives in the membranes
            mem = {c: f for c, f in cw.items() if c not in SOLUBLE_ONLY}
            if mem:
                z = sum(mem.values()); cw = {c: f / z for c, f in mem.items()}
        for c, f in cw.items():
            cp = n * f
            if cp <= 0:
                continue
            def add(reg, x=cp):
                rows.append((lib_id, r.ing_id, reg, x))
            if c == "cytosol":
                add("cytosol" if not tm else "dropped_TM_in_soluble_compartment")
            elif c == "nucleoplasm":
                add("nucleoplasm" if not tm else "dropped_TM_in_soluble_compartment")
            elif c == "nucleolus":
                add("nucleolus" if not tm else "dropped_TM_in_soluble_compartment")
            elif c == "mito_IM":
                add("mem_MITO_IM" if tm else "mito_matrix")     # peripheral IM proteins: mostly matrix-side (F1, bc1 core)
            elif c == "mito_OM":
                add("mem_MITO_OM" if tm else "cytosol")         # peripheral OM proteins face the cytosol
            elif c == "mito_IMS":
                add("mem_MITO_IM" if tm else "mito_IMS")
            elif c == "mito_matrix":
                add("mem_MITO_IM" if tm else "mito_matrix")
            elif c == "mitochondrion":
                if tm:
                    add("mem_MITO_IM", 0.75 * cp); add("mem_MITO_OM", 0.25 * cp)
                else:
                    add("mito_matrix", 0.9 * cp); add("mito_IMS", 0.1 * cp)
            elif c == "ER":
                if tm:   # the outer nuclear membrane is ER: share ER membrane proteins by membrane area
                    add("mem_ER", (1 - ER_ON_NE) * cp); add("mem_NE", ER_ON_NE * cp)
                elif r.ing_id == "ribosome_80S" or side_of(r, c) == "face":
                    # the lumen is reached co-translationally: soluble proteins without a signal peptide (or a lumen
                    # annotation) sit on the cytosolic face (COPII, Sar1, Cdc48, Get3, SRP, Scp160; 2026-11-12)
                    add("ER_peripheral")
                else:
                    add("ER_lumen")
            elif c == "plasma_membrane":
                add("mem_PM" if tm else "PM_peripheral")
            elif c == "nuclear_envelope":
                add("mem_NE" if tm else "NE_peripheral")
            elif c == "vacuole":
                add("mem_VAC" if tm else ("vacuole" if side_of(r, c) == "lumen" else "VAC_peripheral"))
            elif c == "golgi_endosome":
                add("mem_GOLGI" if tm else ("golgi_endosome" if side_of(r, c) == "lumen" else "GOLGI_peripheral"))
            elif c == "peroxisome":
                add("mem_PEROX" if tm else ("peroxisome" if side_of(r, c) == "lumen" else "PEROX_peripheral"))
            elif c == "lipid_droplet":
                add("mem_ER" if tm else "LD_peripheral")
            elif c == "cell_wall":
                add("cell_wall")
    import trna_species            # one type per isoacceptor where the predicted models exist (else one tRNA-Phe)
    rows += trna_species.expected_rows(N_TRNA)
    e = pd.DataFrame(rows, columns=["lib_id", "ing_id", "region", "expected"])
    return e.groupby(["lib_id", "ing_id", "region"], as_index=False).expected.sum()


def crowding_factors(e, M, reg_idx, logmod, vol, iters=30):
    """Per region, tile factors c_t for the modulated types such that the expected macromolecular volume per tile
    stays proportional to the region's measure in the tile (composition varies between tiles, crowding does not)"""
    out = {}
    for reg, g in e[e.lib_id.isin(logmod)].groupby("region"):
        ri = reg_idx[reg]; W = M[:, ri]; nz = np.nonzero(W)[0]
        if len(nz) == 0:
            continue
        A = np.stack([np.exp(logmod[l][nz] - logmod[l][nz].max()) for l in g.lib_id])         # rows x tiles
        V = (g.expected.values * np.array([vol[l] for l in g.lib_id]))[:, None]
        c = np.ones(len(nz))
        for _ in range(iters):
            Wc = A * (W[nz] * c)[None]
            occ = (V * Wc / Wc.sum(1, keepdims=True)).sum(0)
            ratio = (W[nz] / W[nz].sum()) / (occ / occ.sum())
            c *= ratio
        cf = np.zeros(len(W)); cf[nz] = c / c.mean(); out[ri] = cf
    return out


# ------------------------------------------------------------------ capacity of shells and crowded compartments
# Sphere volume fraction = summed sphere volume of the collision shapes per region volume (the spheres of one
# molecule overlap, so this is more than the fraction the shapes really fill: 1.6-1.8 times with the library trees,
# with which the cytosol of build R asked for 0.38). Thresholds are relative to the cytosol's own request.
SHELL_CAP = True        # a membrane-face shell asks for at most the cytosol's sphere volume fraction
THIN = True             # compartments that cannot hold their request are thinned uniformly (pilot_capacity)
THIN_MIN_REL = 1.1      # volume regions asking for more than this multiple of the cytosol's request are tested in the pilot
THIN_TILES = 16         # pilot tiles per region
THIN_SAFETY = 0.97      # the request is cut to this share of what the pilot tiles held


def _sphere_volumes(lp):
    return np.array([(4 / 3 * np.pi * lp["STR"][a:b] ** 3).sum() for a, b in zip(lp["st_off"][:-1], lp["st_off"][1:])])


def requested_phi(exp, measure, lp):
    """requested sphere volume fraction per volume region {name: phi}"""
    v = pd.Series(_sphere_volumes(lp), index=[str(x) for x in lp["ids"]])
    e = exp[exp.lib_id.isin(v.index)]
    req = (e.expected * v.reindex(e.lib_id).values).groupby(e.region).sum()
    tot = measure.sum()
    return {r: float(req[r] / tot[r]) for r in req.index if r in tot.index and not r.startswith("mem_") and tot[r] > 0}


def shell_capacity(exp, measure, lp, log=print):
    """Membrane-face shells (PM, envelope, ER, lipid-droplet, vacuole, Golgi/endosome, peroxisome faces) hold peripheral
    proteins up to the cytosol's own requested density; what a shell is asked to hold beyond that goes to the cytosol,
    every species of the shell by the same factor. A peripheral protein exchanges between its membrane and a soluble
    pool; the annotation gives no bound fraction, and until build R the whole pool was sent to the shell: the
    Golgi/endosome face was asked to hold 0.70 (cytosol 0.38), 84% of it fitted, 17% of the placed copies ended outside
    the shell and cytosolic proteins were ten times depleted next to these membranes. Returns the edited table."""
    phi = requested_phi(exp, measure, lp)
    cap = phi.get("cytosol")
    if cap is None:
        return exp
    exp = exp.copy(); moved = {}
    for r in t3.BANDS:
        name = t3.REGION_NAMES[r]
        if phi.get(name, 0.0) <= cap:
            continue
        f = cap / phi[name]
        k = exp.region == name
        extra = exp[k].assign(expected=exp.loc[k, "expected"] * (1 - f), region="cytosol")
        exp.loc[k, "expected"] *= f
        exp = pd.concat([exp, extra], ignore_index=True)
        moved[name] = (f, float(extra.expected.sum()))
        log(f"  shell {name}: asked for sphere volume fraction {phi[name]:.2f}, capped at the cytosol's {cap:.2f}: "
            f"{f:.2f} of every species stays, {extra.expected.sum():.0f} copies go to the cytosol")
    if moved:
        exp = exp.groupby(["lib_id", "ing_id", "region"], as_index=False).expected.sum()
    return exp


def pilot_capacity(counts, measure, lp, log=print, nproc=NPROC):
    """What crowded compartments can hold: pack THIN_TILES tiles per region that asks for more than THIN_MIN_REL x
    the cytosol's requested sphere volume fraction
    (drawn in proportion to the region's volume in them, not adjacent, no neighbours, one round) and return {region name: factor} with factor =
    THIN_SAFETY x placed / requested sphere volume in those tiles (regions that held >= 98% are left alone). The
    request of every species of such a region is then cut by that factor before the real packing, so that the shortfall
    is spread evenly instead of falling on mid-sized species (build R: Hsp60 49%, Sod2 37%, Ssc1 96% placed in the
    matrix). setup_shared must have been called with `counts`."""
    global RAW
    import multiprocessing as mp, shutil
    v = _sphere_volumes(lp)
    cv = (v[counts.lib.values] * counts.n.values); tot = measure.sum(); M = measure
    reqv = pd.Series(cv).groupby(counts.region.values).sum()
    phi_cyt = reqv.get(t3.R_CYTO, 0.0) / tot["cytosol"]
    regs = [int(r) for r in reqv.index if t3.IS_VOL[r] and tot[t3.REGION_NAMES[r]] > 0
            and reqv[r] / tot[t3.REGION_NAMES[r]] > THIN_MIN_REL * phi_cyt]
    if not regs:
        return {}
    tiles = _S["tiles"]; chosen = []; of = {}
    rng = np.random.default_rng(_seed(4242))
    for r in regs:
        # tiles drawn in proportion to the region's volume in them: a sample of where the region's copies go. (The
        # tiles with the most of the region, used at first, are its roomiest parts: they held 94% of a matrix request
        # of which the cell then held 81%, 2026-11-13.)
        w = M[t3.REGION_NAMES[r]].values.astype(np.float64); mine = []
        for ti in rng.choice(len(w), size=min(len(w), 40 * THIN_TILES), replace=False, p=w / w.sum()) if (w > 0).sum() > 40 * THIN_TILES \
                else rng.permutation(np.nonzero(w)[0]):
            if len(mine) == THIN_TILES:
                break
            if all(max(abs(a - b) for a, b in zip(tiles[ti], tiles[tj])) > 1 for tj in chosen):
                chosen.append(int(ti)); mine.append(int(ti))
        of[r] = mine
    raw0, edges0 = RAW, _S.get("round_edges")
    RAW = raw0 + "_pilot"; shutil.rmtree(RAW, ignore_errors=True); os.makedirs(RAW); _S["round_edges"] = None
    try:
        with mp.get_context("fork").Pool(min(nproc, len(chosen))) as pool:
            pool.map(_pw, [(ti, 0, _seed(500000 + ti), 0, True) for ti in chosen], chunksize=1)
        out = {}
        for r in regs:
            rq = pl = 0.0
            for ti in of[r]:
                z = np.load(f"{RAW}/{t3.tile_name(tiles[ti])}.npz"); k = z["job_region"] == r
                rq += (z["job_n"][k] * v[z["job_lib"][k]]).sum(); pl += (z["job_placed"][k] * v[z["job_lib"][k]]).sum()
            held = pl / max(rq, 1e-9)
            f = 1.0 if held >= 0.98 else THIN_SAFETY * held
            out[t3.REGION_NAMES[r]] = float(f)
            log(f"  pilot {t3.REGION_NAMES[r]}: asked for sphere volume fraction {reqv[r] / tot[t3.REGION_NAMES[r]]:.2f}; {len(of[r])} tiles "
                f"held {held:.3f} of the requested volume -> every species x {f:.3f}")
    finally:
        shutil.rmtree(RAW, ignore_errors=True); RAW = raw0; _S["round_edges"] = edges0
    return {k: f for k, f in out.items() if f < 1.0}


def thin(exp, factors):
    """the packing request: expected copies x the region's factor (the inventory table itself is not changed)"""
    exp = exp.copy()
    for name, f in factors.items():
        exp.loc[exp.region == name, "expected"] *= f
    return exp


def distribute(exp, measure, lib_index, seed=3, logmod=None, vol=None):
    """Multinomial split of each (lib, region) total over tiles proportional to the tile's region measure, times
    exp(logmod[lib]) per tile if given (organisation.tile_logmod: large-scale co-enrichment of STRING modules) and
    the crowding correction for those types (crowding_factors, needs vol: lib id -> volume).
    Returns long table tile_i, lib_i, region_i, n (cell totals conserved after stochastic rounding)."""
    seed = _seed(seed); rng = np.random.default_rng(seed)
    M = measure.values
    reg_idx = {n: i for i, n in enumerate(t3.REGION_NAMES)}
    e = exp[exp.region.isin(reg_idx)].groupby(["lib_id", "region"], as_index=False).expected.sum()
    cf = crowding_factors(e, M, reg_idx, logmod, vol) if logmod is not None else {}
    T, Li, Ri, Nn = [], [], [], []
    for r in e.itertuples():
        ri = reg_idx[r.region]
        w = M[:, ri]
        if w.sum() <= 0:
            continue
        N = int(np.floor(r.expected) + (rng.random() < r.expected - np.floor(r.expected)))
        if N == 0:
            continue
        nz = np.nonzero(w)[0]
        if logmod is not None and r.lib_id in logmod:
            lw = logmod[r.lib_id][nz]
            w = np.zeros_like(w); w[nz] = M[nz, ri] * np.exp(lw - lw.max()) * (cf[ri][nz] if ri in cf else 1.0)
        cnt = rng.multinomial(N, w[nz] / w[nz].sum())
        k = cnt > 0
        T.append(nz[k]); Li.append(np.full(k.sum(), lib_index[r.lib_id])); Ri.append(np.full(k.sum(), ri)); Nn.append(cnt[k])
    d = pd.DataFrame(dict(tile=np.concatenate(T).astype(np.int32), lib=np.concatenate(Li).astype(np.int32),
                          region=np.concatenate(Ri).astype(np.int16), n=np.concatenate(Nn).astype(np.int32)))
    return d.sort_values(["tile"]).reset_index(drop=True)


# ------------------------------------------------------------------ pre-pass: NPCs + chromatin
def prepass(cell, lp, n_nuc, seed=5, log=print, before_chromatin=None):
    """NPCs, then chromatin. before_chromatin(npc_spheres) may place more objects first and return a dict with
    'obstacles' (centres, radii) for the walks and optional 'starts' (one start point per chromosome)."""
    li_npc = int(np.where(lp["ids"] == "NPC")[0][0]); li_nuc = int(np.where(lp["ids"] == "nucleosome")[0][0])
    # NPC instances: library +z (cytoplasmic) along the pore normal
    npc_pos = cell.npc_pos
    q = np.empty(4); npc_q = []
    seed = _seed(seed); rng = np.random.default_rng(seed)
    for u in cell.pores[:, 3:6]:          # NE normals
        p3.quat_z_to(u.astype(np.float64), rng.uniform(0, 2 * np.pi), q); npc_q.append(q.copy())
    npc_q = np.array(npc_q)
    types = np.full(len(npc_pos), li_npc, np.int32)
    cn = cell.G[3:6]; rn = float(cell.G[11:14].max())
    glo = cn - rn - 80; ghi = cn + rn + 80
    OC, OR = p3.emit_spheres(types, npc_pos.astype(np.float32), npc_q.astype(np.float32), lp["st_off"], lp["STC"], lp["STR"], lp["BR"], glo, ghi)
    lens = np.array(CHROM_KB, float); lens = np.round(lens / lens.sum() * n_nuc).astype(np.int64)
    dirn = cell.G[7:10] - cn; dirn /= np.linalg.norm(dirn)
    spb = c3.nuc_surface(cell.G, -dirn, -(c3.PERI_W + 40))[0]
    starts, r_start = np.tile(spb, (len(lens), 1)), 150.0
    extra = {}
    if before_chromatin is not None:
        extra = before_chromatin(p3.emit_spheres(types, npc_pos.astype(np.float32), npc_q.astype(np.float32), lp["st_off"],
                                                 lp["STC"], lp["STR"], lp["BR"], np.full(3, -1e4), np.full(3, 1e4)))
        if "obstacles" in extra:
            C, R = extra["obstacles"]; k = np.all((C > glo) & (C < ghi), axis=1)
            OC = np.concatenate([OC, C[k]]); OR = np.concatenate([OR, R[k]])
        if len(extra.get("starts", ())) == len(lens):     # chromosome walks start at the kinetochore-MT plus ends
            starts, r_start = np.asarray(extra["starts"], np.float64), 20.0
    a = cell.args(glo, ghi, margin=200.0)
    s0, s1 = lp["st_off"][li_nuc], lp["st_off"][li_nuc + 1]
    mask = 1 << c3.NUCP
    nsph = (s1 - s0) * lens.sum() + len(OR)
    t0 = time.time()
    sl, ss, st, sc, sk = rabl_segments(cell, lens, starts, np.random.default_rng(seed + 100), log=log)
    dy = None
    if FIBER and NUC_MAP:
        # measured nucleosome positions (nucleosome_map.py): one nucleosome per dyad, in walk order; the spacing
        # sets the backbone advance and the reach to the previous nucleosome (wide across nucleosome-free regions)
        import nucleosome_map, genome_map
        dy = nucleosome_map.segment_dyads(genome_map.segments())
        sl = np.array([len(p) for p in dy], np.int64)
        log(f"nucleosome map: {sl.sum()} measured dyads on the arms (inventory {int(lens.sum())})")
    if FIBER:
        sp = None if dy is None else np.concatenate([np.r_[NRL_BP, np.abs(np.diff(p))] for p in dy])
        delta, lk = fiber_steps(sl, sp)
        pos, qu, chain, brk = p3.fiber_walks(seed, glo, ghi, *a, lp["STC"][s0:s1], lp["STR"][s0:s1], lp["BR"][li_nuc], mask,
                                             OC, OR, sl, ss, st, r_start, delta, FIBER_LP, FIBER_R,
                                             lk, int(nsph * 1.1) + 1000, int(nsph * 10) + 10000)
    else:
        pos, qu, chain, brk = p3.arm_walks(seed, glo, ghi, *a, lp["STC"][s0:s1], lp["STR"][s0:s1], lp["BR"][li_nuc], mask,
                                           OC, OR, sl, ss, st, r_start, 12.5, int(nsph * 1.1) + 1000, int(nsph * 10) + 10000)
    log(f"chromatin: {len(pos)}/{sl.sum()} nucleosomes, {int(brk.sum())} fragments, {time.time()-t0:.0f}s")
    if dy is not None:          # walks place their nucleosomes in order: the first n of each segment's dyads
        extra["nuc_bp"] = np.concatenate([dy[c][:int((chain == c).sum())] for c in range(len(dy))]).astype(np.int64)
        assert len(extra["nuc_bp"]) == len(pos) and np.all(np.diff(chain) >= 0)
    out = dict(npc_pos=npc_pos.astype(np.float32), npc_q=npc_q.astype(np.float32), li_npc=li_npc, li_nuc=li_nuc,
               nuc_pos=pos, nuc_q=qu, nuc_chain=chain, nuc_brk=brk, chain_target=sl,
               seg_start=ss, seg_target=st, seg_chrom=sc, seg_kind=sk.astype("U8"))
    out.update({k: v for k, v in extra.items() if k not in ("obstacles", "starts")})
    return out


def fiber_steps(seg_len, spacing_bp=None):
    """per nucleosome (segments concatenated): backbone advance (nm) and maximum centre distance from the previous
    nucleosome (nm). spacing_bp: dyad-to-dyad distance before each nucleosome (measured map); None: the uniform
    NRL_BP repeat. The extra linker DNA beyond a ~20 bp linker can stretch between neighbours (0.34 nm per bp,
    capped at 60 nm, e.g. across a nucleosome-free region)."""
    n = int(np.sum(seg_len))
    sp = np.full(n, NRL_BP) if spacing_bp is None else np.asarray(spacing_bp, np.float64)
    delta = sp / FIBER_BP_PER_NM
    lk = np.minimum(60.0, LINK_MAX + 0.34 * np.maximum(0.0, sp - NRL_BP))
    return delta.astype(np.float64), lk.astype(np.float64)


def rabl_segments(cell, lens, starts, rng, log=print):
    """chromosome-arm segments: (nucleosomes, start, target, chromosome, kind) with kind 'L'/'R' (centromere to
    telomere) or 'XIIR_in'/'XIIR_out' (chromosome XII-R up to the nucleolus, and from the nucleolus to its telomere;
    the rDNA itself is built in the nucleolus by rdna_chromatin)"""
    G = cell.G; cn = G[3:6]; ncl, rncl = G[7:10], G[10]
    spb_u = starts.mean(0) - cn; spb_u /= np.linalg.norm(spb_u)
    e1 = np.cross(spb_u, [0.3, 0.5, 0.8]); e1 /= np.linalg.norm(e1); e2 = np.cross(spb_u, e1)
    segs, tels = [], []
    for c, (kb, cen) in enumerate(zip(CHROM_KB, CEN_KB)):
        nn = lens[c] / kb
        if c == 11:
            segs += [[nn * cen, c, "L"], [nn * (RDNA_KB[0] - cen), c, "XIIR_in"], [nn * (kb - RDNA_KB[1]), c, "XIIR_out"]]
            tels += [(c, "L", cen), (c, "XIIR_out", kb - cen)]
        else:
            segs += [[nn * cen, c, "L"], [nn * (kb - cen), c, "R"]]
            tels += [(c, "L", cen), (c, "R", kb - cen)]
    # fibre contour of each telomere arm (nm), for the reach limit on its target
    step = NRL_BP / FIBER_BP_PER_NM if FIBER else 12.5
    contour = {(c, k): lens[c] / CHROM_KB[c] * arm * step for c, k, arm in tels}
    # telomere sites: angle from the SPB axis by arm length, gathered into TEL_CLUSTERS foci of similar arm length
    th = {(c, k): np.radians(TEL_THETA[0] + (TEL_THETA[1] - TEL_THETA[0]) * min(1.0, arm / 700.0)) for c, k, arm in tels}
    order = sorted(th, key=th.get)
    groups = np.array_split(np.arange(len(order)), TEL_CLUSTERS)
    dirs = lambda t, p: np.cos(t) * spb_u + np.sin(t) * (np.cos(p) * e1 + np.sin(p) * e2)
    target, centres = {}, []
    for g in groups:
        tg = float(np.mean([th[order[i]] for i in g]))
        for _ in range(500):
            ph = rng.uniform(0, 2 * np.pi); u = dirs(tg, ph)
            if all(np.degrees(np.arccos(np.clip(u @ v, -1, 1))) > 35.0 for v in centres):
                break
        while cell.label(c3.nuc_surface(G, u, -(c3.PERI_W + TEL_DEPTH))[0])[0] != c3.NUCP and tg > np.radians(20):
            tg -= np.radians(5); u = dirs(tg, ph)                    # keep the focus out of the nucleolus
        centres.append(u)
        for i in g:
            v = u + rng.normal(0, 0.04, 3); v /= np.linalg.norm(v)    # ~ +-60 nm around the focus
            t_ = c3.nuc_surface(G, v, -(c3.PERI_W + TEL_DEPTH))[0]
            # a telomere cannot sit further from its start than the arm's contour allows: move short arms' targets
            # along the great circle towards the pole until within TEL_REACH of the contour
            c_, k_ = order[i]
            st0 = starts[c_]; Lc = contour[(c_, k_)]
            for _ in range(60 if k_ != "XIIR_out" else 0):          # XII-R leaves from the nucleolus
                if np.linalg.norm(t_ - st0) <= TEL_REACH * Lc:
                    break
                v = v + 0.06 * spb_u; v /= np.linalg.norm(v)
                t_ = c3.nuc_surface(G, v, -(c3.PERI_W + TEL_DEPTH))[0]
            target[order[i]] = t_
    # chromosome XII-R: into the nucleolus (rDNA) and out again, at two points on its border facing the SPB
    t_in = ncl + (rncl + 40.0) * (_unit(spb_u + 0.35 * e1)); t_out = ncl + (rncl + 40.0) * _unit(spb_u - 0.35 * e1)
    L, S, T, C, K = [], [], [], [], []
    for n_, c, k in segs:
        L.append(int(round(n_))); C.append(c); K.append(k)
        S.append(t_out if k == "XIIR_out" else starts[c])
        T.append(t_in if k == "XIIR_in" else target[(c, k)])
    log(f"Rabl layout: {len(L)} arm segments, telomeres in {TEL_CLUSTERS} foci at the nuclear envelope "
        f"({TEL_THETA[0]:.0f}-{TEL_THETA[1]:.0f} deg from the SPB axis by arm length); XII-R through the nucleolus")
    return np.array(L, np.int64), np.array(S, np.float64), np.array(T, np.float64), np.array(C, np.int32), np.array(K)


def _unit(v):
    return v / np.linalg.norm(v)


def rdna_chromatin(cell, lp, n, obstacles, n_chains=8, seed=9, log=print):
    """nucleosomal rDNA (inactive repeats, spacers of active ones): persistent walks inside the nucleolus"""
    li_nuc = int(np.where(lp["ids"] == "nucleosome")[0][0])
    cn = cell.G[3:6]; rn = float(cell.G[11:14].max())
    glo = cn - rn - 80; ghi = cn + rn + 80
    a = cell.args(glo, ghi, margin=200.0)
    seed = _seed(seed); rng = np.random.default_rng(seed)
    box = rng.uniform(glo, ghi, size=(20000, 3)); pts = box[cell.label(box) == c3.NUCL]
    starts = pts[rng.choice(len(pts), n_chains, replace=False)].astype(np.float64)
    lens = np.full(n_chains, n // n_chains, np.int64); lens[: n - lens.sum()] += 1
    OC, OR = obstacles; k = np.all((OC > glo) & (OC < ghi), axis=1)
    OC, OR = OC[k].astype(np.float64), OR[k].astype(np.float64)
    s0, s1 = lp["st_off"][li_nuc], lp["st_off"][li_nuc + 1]
    nsph = (s1 - s0) * lens.sum() + len(OR)
    pos, qu, chain, brk = p3.chromatin_walk(seed, glo, ghi, *a, lp["STC"][s0:s1], lp["STR"][s0:s1], lp["BR"][li_nuc], 1 << c3.NUCL,
                                            OC, OR, lens, starts, 40.0, 12.5, int(nsph * 1.1) + 1000, int(nsph * 10) + 10000)
    log(f"rDNA chromatin: {len(pos)}/{n} nucleosomes in the nucleolus, {int(brk.sum())} fragments")
    return dict(rnuc_pos=pos, rnuc_q=qu, rnuc_chain=chain, rnuc_brk=brk)


# ------------------------------------------------------------------ tile worker
_S = {}   # shared state (set before fork)


PRE_GROUPS = [("npc", "li_npc"), ("nuc", "li_nuc"), ("poly", "li_ribo"), ("act", "li_actin"), ("tpm", "li_tpm"), ("mt", None), ("pat", None), ("atp", "li_atp"), ("eis", None), ("cvt", None), ("pol2", "li_pol2"), ("nvj", "li_nvj"), ("rvs", "li_rvs"), ("coat", None), ("mcc", None), ("kt", None), ("mcm", "li_mcm"), ("orc", "li_orc"), ("tf", None), ("rnp", None), ("kto", None), ("spb", None), ("lm", None), ("kap", None),
              ("rpol", "li_rpol"), ("r90", "li_r90"), ("rnuc", "li_nuc")]       # rDNA: Pol I, 90S on transcripts, nucleosomes


def pre_instances(pre):
    """types, positions, quaternions of all pre-placed instances: NPCs, nucleosomes, polysome ribosomes, actin,
    microtubules (+ gamma-TuRC caps), actin patches (actin, Arp2/3, capping protein); per-instance types in *_type"""
    parts = [(k, li) for k, li in PRE_GROUPS if f"{k}_pos" in pre]
    gt = np.concatenate([pre[f"{k}_type"].astype(np.int32) if f"{k}_type" in pre else np.full(len(pre[f"{k}_pos"]), int(pre[li]), np.int32)
                         for k, li in parts])
    gp = np.concatenate([pre[f"{k}_pos"] for k, _ in parts]); gq = np.concatenate([pre[f"{k}_q"] for k, _ in parts])
    return gt, gp, gq


def setup_shared(cell, lp, counts, tiles, pre, halo, partners=None, prm=None):
    """partners, prm: organisation.partners() CSR arrays and organisation.params() (None: uniform random placement)"""
    _S.update(cell=cell, lp=lp, tiles=tiles, halo=halo, tindex={t: i for i, t in enumerate(tiles)})
    _S["partners"] = partners; _S["prm"] = prm
    _S["counts"] = {k: g for k, g in counts.groupby("tile")}
    gt, gp, gq = pre_instances(pre)
    _S["glob"] = (gt, gp, gq)
    # pre-placed RNA beads (polysome linkers, full-length mRNA threads) as packing obstacles, bucketed by tile
    pb = [pre[k].astype(np.float32) for k in ("poly_link", "mrna_pts", "fg_pts", "rdna_pts", "prerrna_pts", "mtdna_pts", "pol2dna_pts", "pol2rna_pts") if k in pre and len(pre[k])]
    pb = np.concatenate(pb) if pb else np.zeros((0, 3), np.float32)
    key = np.floor(pb / t3.TILE).astype(np.int64)
    order = np.lexsort((key[:, 2], key[:, 1], key[:, 0])); pb = pb[order]; key = key[order]
    brk = np.r_[0, np.where(np.any(np.diff(key, axis=0) != 0, axis=1))[0] + 1, len(key)]
    _S["rna_beads"] = pb
    _S["rna_bucket"] = {tuple(key[a]): (a, b) for a, b in zip(brk[:-1], brk[1:]) if b > a}
    _S["mrna"] = None


def neighbours(t):
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                if dx or dy or dz:
                    yield (t[0] + dx, t[1] + dy, t[2] + dz)


def job_rounds(br, is_mem):
    """packing round of each job: membrane proteins and the largest volume species in round 0, smaller size classes
    in later rounds (_S["round_edges"]: descending bounding-radius thresholds; absent: one round)"""
    edges = _S.get("round_edges")
    if edges is None or len(edges) == 0:
        return np.zeros(len(br), np.int64)
    return np.where(is_mem, 0, np.searchsorted(-np.asarray(edges), -np.asarray(br), side="left")).astype(np.int64)


def round_edges(counts, lp, n_rounds):
    """bounding-radius thresholds that split the requested volume-region copies into n_rounds size classes of about
    equal summed sphere volume (largest first)"""
    if n_rounds <= 1:
        return np.zeros(0)
    c = counts[~t3.IS_MEM[counts.region.values]]
    br = lp["BR"][c.lib.values]
    v = np.array([(lp["STR"][a:b] ** 3).sum() for a, b in zip(lp["st_off"][:-1], lp["st_off"][1:])])[c.lib.values] * c.n.values
    o = np.argsort(-br); cum = np.cumsum(v[o]) / v.sum()
    return np.array([br[o][np.searchsorted(cum, k / n_rounds)] for k in range(1, n_rounds)])


def pack_one(ti, n_mrna, seed, rnd=0, last=True, max_tries=None, give_up=None):
    """pack the jobs of round `rnd` into tile ti, around what earlier rounds placed there (read back from its file);
    the squeeze runs in the last round, on all soluble molecules of the tile"""
    max_tries = PACK_TRIES if max_tries is None else max_tries; give_up = PACK_GIVE_UP if give_up is None else give_up
    cell, lp, halo = _S["cell"], _S["lp"], _S["halo"]
    t = _S["tiles"][ti]
    lo = np.array(t, float) * t3.TILE; hi = lo + t3.TILE
    glo = lo - halo; ghi = hi + halo
    a = cell.args(lo, hi, margin=80.0)
    g = np.load(f"{t3.GDIR}/{t3.tile_name(t)}.npz")
    # membrane triangles grouped by membrane id 0..7
    V, F, FM = g["verts"].astype(np.float64), g["faces"], g["face_mem"]
    order = np.argsort(FM, kind="stable")
    tri = V[F[order]].reshape(-1, 9) if len(F) else np.zeros((0, 9))
    fm = FM[order]
    tri_off = np.searchsorted(fm, np.arange(9)).astype(np.int64)
    ar = 0.5 * np.linalg.norm(np.cross(tri[:, 3:6] - tri[:, 0:3], tri[:, 6:9] - tri[:, 0:3]), axis=1) if len(tri) else np.zeros(0)
    cuma = np.cumsum(ar) if len(ar) else np.zeros(1)
    # voxel seeds per volume region
    n = int(t3.TILE / t3.VOX)
    band_m = np.array([t3.BANDS[r][0] for r in t3.BAND_ORDER], np.int64)
    band_w = np.array([t3.BANDS[r][1] for r in t3.BAND_ORDER], np.float64)
    lab, bits = t3._voxel_regions(lo, n, t3.VOX, *a, band_m, band_w, t3.WALL_DEPTH)
    masks = t3.region_masks(lab, bits)
    vox_off = [0]; VP = []
    for r in range(t3.NREG):
        if t3.IS_VOL[r]:
            ijk = np.argwhere(masks[r]); VP.append(lo + (ijk + 0.5) * t3.VOX)
            vox_off.append(vox_off[-1] + len(ijk))
        else:
            vox_off.append(vox_off[-1])
    VOXP = np.concatenate(VP) if VP else np.zeros((0, 3))
    # jobs: membranes first, then volume; each block largest bounding radius first
    c = _S["counts"].get(ti)
    if c is None:
        jl = np.zeros(0, np.int64); jr = np.zeros(0, np.int64); jn = np.zeros(0, np.int64)
    else:
        c = c.assign(br=lp["BR"][c.lib.values], mem=t3.IS_MEM[c.region.values])
        c = c[job_rounds(c.br.values, c.mem.values) == rnd]
        c = c.sort_values(["mem", "br"], ascending=[False, False])
        jl = c.lib.values.astype(np.int64); jr = c.region.values.astype(np.int64); jn = c.n.values.astype(np.int64)
    fn = f"{RAW}/{t3.tile_name(t)}.npz"
    own = np.load(fn) if rnd > 0 and os.path.exists(fn) else None
    if own is not None and len(jl) == 0 and not last:
        return ti, 0, 0, bool(own["overflow"]), 0.0, 0, 0, 0          # nothing to add in this round
    # obstacles: global (NPC, chromatin) + finished neighbour tiles
    gt, gp, gq = _S["glob"]
    OCs, ORs = [], []
    ATs, APs = [], []         # anchors for partner attraction: placed instances around the tile
    k = np.all((gp >= glo) & (gp <= ghi), axis=1); ATs.append(gt[k]); APs.append(gp[k])
    C0, R0 = p3.emit_spheres(gt, gp, gq, lp["st_off"], lp["STC"], lp["STR"], lp["BR"], glo, ghi)
    OCs.append(C0); ORs.append(R0)
    sl_ = [(_S["rna_bucket"][nt]) for nt in [t] + list(neighbours(t)) if nt in _S["rna_bucket"]]   # pre-placed RNA beads
    if sl_:
        pl = np.concatenate([_S["rna_beads"][a:b] for a, b in sl_]).astype(np.float64)
        k = np.all((pl >= glo - 1) & (pl <= ghi + 1), axis=1)
        OCs.append(pl[k]); ORs.append(np.full(k.sum(), 0.5))
    for nt in neighbours(t):
        f = f"{RAW}/{t3.tile_name(nt)}.npz"
        if nt in _S["tindex"] and os.path.exists(f):
            z = np.load(f)
            C1, R1 = p3.emit_spheres(z["types"], z["pos"], z["quat"], lp["st_off"], lp["STC"], lp["STR"], lp["BR"], glo, ghi)
            OCs.append(C1); ORs.append(R1)
            k = np.all((z["pos"] >= glo) & (z["pos"] <= ghi), axis=1); ATs.append(z["types"][k]); APs.append(z["pos"][k])
            mp_ = z["mrna_pts"]
            if len(mp_):
                k = np.all((mp_ >= glo - 1) & (mp_ <= ghi + 1), axis=1)
                OCs.append(mp_[k].astype(np.float64)); ORs.append(np.full(k.sum(), 0.5))
    OCx = np.concatenate(OCs); ORx = np.concatenate(ORs)          # everything outside this tile's own molecules
    if own is not None and len(own["types"]):
        C1, R1 = p3.emit_spheres(own["types"], own["pos"], own["quat"], lp["st_off"], lp["STC"], lp["STR"], lp["BR"], glo, ghi)
        OC = np.concatenate([OCx, C1]); OR = np.concatenate([ORx, R1])
        ATs.append(own["types"]); APs.append(own["pos"])
    else:
        OC, OR = OCx, ORx
    nsph = int((np.diff(lp["st_off"])[jl] * jn).sum()) if len(jl) else 0
    cap_s = len(OR) + nsph + n_mrna * 8000 + 1000
    cap_e = cap_s * 14
    t0 = time.time()
    pr = _S.get("partners")
    if pr is None:
        pr = dict(off=np.zeros(0, np.int64), idx=np.zeros(0, np.int64), strength=np.zeros(0), kind=np.zeros(0, np.int64))
    prm = _S.get("prm") if _S.get("prm") is not None else np.zeros(7)
    AT = np.concatenate(ATs).astype(np.int64); AP = np.concatenate(APs).astype(np.float64)
    res = p3.pack_tile(seed + 7919 * rnd, lo, t3.TILE, halo, *a, lp["st_off"], lp["STC"], lp["STR"], lp["BR"], t3.REG_MASK,
                       jl, jr, jn, t3.MEM_REGION0, np.array(vox_off, np.int64), VOXP, t3.VOX,
                       tri_off, tri, cuma, OC, OR, n_mrna, MRNA_MEAN_NT, t3.R_CYTO, max_tries, give_up, cap_s, cap_e,
                       AT, AP, pr["off"], pr["idx"], pr["strength"], pr["kind"], prm,
                       *thin_budget(max_tries, give_up), _flat_q(lp), *band_arrays())
    types, pos, quat, placed, mpts, moff, overflow, ns, ne, att = res
    n_req, n_new = int(jn.sum()), int(placed.sum())
    if own is not None:               # instances stay grouped by job: earlier rounds' jobs, then this round's
        types = np.concatenate([own["types"], types]); pos = np.concatenate([own["pos"], pos]); quat = np.concatenate([own["quat"], quat])
        jl = np.concatenate([own["job_lib"].astype(np.int64), jl]); jr = np.concatenate([own["job_region"].astype(np.int64), jr])
        jn = np.concatenate([own["job_n"].astype(np.int64), jn]); placed = np.concatenate([own["job_placed"].astype(np.int64), placed])
        overflow = bool(overflow) or bool(own["overflow"])
    if SQUEEZE and last and len(jl) and not overflow:
        cs = len(ORx) + int((np.diff(lp["st_off"])[jl] * jn).sum()) + 1000
        types, pos, quat, placed = squeeze_tile(seed, lo, halo, a, lp, jl, jr, jn, placed, types, pos, quat,
                                                np.array(vox_off, np.int64), VOXP, OCx, ORx, cs, cs * 14)
        n_new = int(placed.sum()) - (int(own["job_placed"].sum()) if own is not None else 0)
    np.savez(fn, types=types, pos=pos, quat=quat, job_lib=jl.astype(np.int32),
             job_region=jr.astype(np.int16), job_n=jn.astype(np.int32), job_placed=placed.astype(np.int32),
             mrna_pts=mpts, mrna_off=moff, overflow=overflow)
    return ti, n_req, n_new, bool(overflow), time.time() - t0, len(moff) - 1, int(att[0]), int(att[1])


SQUEEZE = True          # overlap relaxation after RSA (squeeze.py; literature item L1)

# Tiles are packed in 8 passes by tile parity (tiles of one pass never touch; later passes read finished neighbours).
# With one visit per tile, first-pass tiles over-filled the 5 nm at their faces (bodies protrude into still empty
# neighbours) and last-pass tiles under-filled it: sphere volume fraction 0.49 against 0.20 there, 0.40 inside
# (build R; validation/placement_checks.py). Now every tile is visited N_ROUNDS times, one size class per round
# (largest first; classes of equal requested volume), and the order of the 8 passes is reversed in every second
# round, so each tile is as often ahead of a neighbour as behind it.
N_ROUNDS = 6

# Attempt budget (IMPROVEMENTS N1, after the Jurkat model v7): with 200 tries per copy and 25 failures in a row, tiles
# stopped early; test tiles placed 2-10 points more of the Golgi face, envelope periphery, mitochondrial matrix and
# cytosol with 2,000 and 1,200. In thin lumens and membrane shells half of the tries align the molecule's flattest
# axis with the normal of the bounding membrane (small further gains in shells).
PACK_TRIES, PACK_GIVE_UP = 2000, 1200
FLAT = True
FLAT_MEM = {t3.R_ERL: c3.M_ER, t3.R_MATRIX: c3.M_IM, t3.R_IMS: c3.M_OM, t3.R_GOLGI: c3.M_GOLGI,
            **{r: m for r, (m, _) in t3.BANDS.items()}}


def thin_budget(max_tries, give_up):
    """per region: tries per copy, failures in a row (0: the global values) and membrane for flat alignment (-1: none)"""
    tr = np.zeros(t3.NREG, np.int64); gu = np.zeros(t3.NREG, np.int64); fm = np.full(t3.NREG, -1, np.int64)
    if FLAT:
        for r, m in FLAT_MEM.items():
            fm[r] = m
    return tr, gu, fm


SHELL_TOL = 4.0         # nm: a membrane-face protein may sit this far beyond its shell width


def band_arrays():
    """per region: membrane of its shell (-1: not a shell region) and the largest allowed distance of a molecule
    centre from that membrane (shell width + SHELL_TOL)"""
    bm = np.full(t3.NREG, -1, np.int64); bw = np.zeros(t3.NREG)
    for r, (m, w) in t3.BANDS.items():
        bm[r] = m; bw[r] = w + SHELL_TOL
    return bm, bw


def _flat_q(lp):
    """per library type: quaternion (w, x, y, z) taking its flattest principal axis (sphere tree, weighted by
    sphere volume) to +z; cached on lp"""
    if "flat_q" not in lp:
        q = np.zeros((len(lp["BR"]), 4)); q[:, 0] = 1.0
        for i in range(len(lp["BR"])):
            C = lp["STC"][lp["st_off"][i]:lp["st_off"][i + 1]].astype(np.float64)
            w = lp["STR"][lp["st_off"][i]:lp["st_off"][i + 1]].astype(np.float64) ** 3
            if len(C) < 3:
                continue
            X = C - (w[:, None] * C).sum(0) / w.sum()
            a = np.linalg.eigh((w[:, None] * X).T @ X)[1][:, 0]          # smallest-variance axis
            if a[2] < 0:
                a = -a
            # rotation a -> z: axis a x z, angle acos(a_z)
            ax = np.array([a[1], -a[0], 0.0]); s = np.linalg.norm(ax)
            if s > 1e-9:
                ang = np.arccos(np.clip(a[2], -1, 1))
                q[i] = [np.cos(ang / 2), *(np.sin(ang / 2) * ax / s)]
        lp["flat_q"] = q
    return lp["flat_q"]


def squeeze_tile(seed, lo, halo, a, lp, jl, jr, jn, placed, types, pos, quat, vox_off, VOXP, OC, OR, cap_s, cap_e):
    """insert the copies RSA could not place in rounds (at most squeeze.ROUND_FRAC more per round), relax the tile's
    soluble molecules, keep the overlap-free ones; a round that does not increase the count is undone. Returns types,
    pos, quat, placed (instances grouped by job, as pack_tile returns them)"""
    import squeeze as sq
    rng = np.random.default_rng(seed + 7)
    typ = types.astype(np.int64); P = pos.astype(np.float64); Q = quat.astype(np.float64)
    J = np.repeat(np.arange(len(jl)), placed)
    vol_job = np.array([t3.IS_VOL[jr[j]] and vox_off[jr[j] + 1] > vox_off[jr[j]] for j in range(len(jl))])
    for rnd in range(sq.N_ROUNDS):
        have = np.bincount(J, minlength=len(jl)); deficit = np.where(vol_job, jn - have, 0).clip(min=0)
        if deficit.sum() == 0:
            break
        n_mov = int(vol_job[J].sum())
        take = deficit * min(1.0, max(50, sq.ROUND_FRAC * n_mov) / deficit.sum())
        take = np.floor(take + rng.random(len(take))).astype(np.int64)
        nt, nj, npos = [], [], []
        for j in np.where(take > 0)[0]:
            r = jr[j]; v = rng.integers(vox_off[r], vox_off[r + 1], take[j])
            npos.append(VOXP[v] + (rng.random((take[j], 3)) - 0.5) * t3.VOX); nt.append(np.full(take[j], jl[j])); nj.append(np.full(take[j], j))
        if not nt:
            break
        nq = rng.normal(size=(sum(map(len, nt)), 4)); nq /= np.linalg.norm(nq, axis=1, keepdims=True)
        typ2 = np.concatenate([typ, np.concatenate(nt)]); P2 = np.concatenate([P, np.concatenate(npos)])
        Q2 = np.concatenate([Q, nq]); J2 = np.concatenate([J, np.concatenate(nj)])
        reg = jr[J2].astype(np.int64); movable = t3.IS_VOL[reg]
        sq.relax(seed + rnd, lo, t3.TILE, halo, *a, lp["st_off"], lp["STC"], lp["STR"], t3.REG_MASK, typ2, P2, Q2, reg, movable,
                 OC, OR, sq.N_ITER)
        keep = sq.accept(lo, t3.TILE, halo, *a, lp["st_off"], lp["STC"], lp["STR"], lp["BR"], t3.REG_MASK, typ2, P2, Q2, reg,
                         movable, OC, OR, cap_s + len(typ2) * 40, (cap_s + len(typ2) * 40) * 14, *band_arrays())
        if keep.sum() <= len(typ):
            break                                         # no gain: keep the previous configuration
        typ, P, Q, J = typ2[keep], P2[keep], Q2[keep], J2[keep]
    o = np.argsort(J, kind="stable")
    return (typ[o].astype(types.dtype), P[o].astype(pos.dtype), Q[o].astype(quat.dtype),
            np.bincount(J, minlength=len(jl)).astype(placed.dtype))


def _pw(args):
    nb.set_num_threads(1)
    return pack_one(*args)


def run_phase(tile_ids, mrna_counts, nproc=NPROC, seed0=1000, log=print, rnd=0, last=True):
    import multiprocessing as mp
    os.makedirs(RAW, exist_ok=True)
    jobs = [(ti, int(mrna_counts[ti]), _seed(seed0 + ti), rnd, last) for ti in tile_ids]
    # heavier tiles first for load balance
    jobs.sort(key=lambda j: -int(_S["counts"][j[0]].n.sum()) if j[0] in _S["counts"] else 0)
    out = []
    with mp.get_context("fork").Pool(nproc) as pool:
        for r in pool.imap_unordered(_pw, jobs, chunksize=1):
            out.append(r)
    return out


def phase_of(t):
    return (t[0] & 1) + 2 * (t[1] & 1) + 4 * (t[2] & 1)
