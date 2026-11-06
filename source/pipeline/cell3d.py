"""3D whole-cell geometry for a haploid S. cerevisiae cell as analytic signed-distance primitives.
Units nm, origin at cell centre. Organelle volumes/areas are data-driven (inventory geo dict); their
arrangement is illustrative (no tomogram).

Membrane SDF sign convention: negative on the *lumenal / non-cytoplasmic* side; the gradient therefore
points to the cytoplasmic (N-) side, which is +z in the membrane-protein library frame.
For the mitochondrial inner membrane the matrix is the N-side, so SDF_IM = -SDF_matrix."""
import json
import numpy as np
import numba as nb

# ----------------------------------------------------------------- labels / membranes
EXT, MANNAN, GLUCAN, CYTO, CER, PERI, NUCP, NUCL, VAC, IMS, MATRIX, LD, GOLGI, PEROX, PORE = range(15)
LABEL_NAMES = ["exterior", "wall_mannan", "wall_glucan", "cytosol", "cER_lumen", "perinuclear", "nucleoplasm",
               "nucleolus", "vacuole", "mito_IMS", "mito_matrix", "lipid_droplet", "golgi_endosome",
               "peroxisome", "pore_channel"]
M_PM, M_ER, M_NE, M_VAC, M_OM, M_IM, M_GOLGI, M_PEROX, M_LD, M_OUTER = range(10)
MEM_NAMES = ["PM", "ER", "NE", "VAC", "MITO_OM", "MITO_IM", "GOLGI", "PEROX", "LD", "CELL_OUTER"]
NMEM = 10
# primitive types and classes
T_SPHERE, T_CAPSULE, T_PILL, T_CERPATCH = range(4)
C_VAC, C_LD, C_GOLGI, C_PEROX, C_OM, C_IM, C_CRISTA, C_CER = range(8)
HALF = 2.0           # bilayer half thickness used for exclusion (as in the column model)
PORE_R = 45.0        # pore-membrane radius at the NE midplane
PERI_W = 34.0        # perinuclear space (ONM-INM midsurface distance)
IMS_W = 8.0          # OM-IM midsurface distance
MITO_R = 150.0
CRISTA_HALF = 9.0    # crista membrane half-separation = radius of the rim tube (fits ATP synthase dimers, see atp_dimers.py)
CER_GAP = 27.0       # PM midsurface -> outer cER membrane midsurface
CER_LUMEN = 20.0     # cER membrane midsurface separation


# ----------------------------------------------------------------- numba evaluators
@nb.njit(cache=True, fastmath=True)
def sdf_prim(px, py, pz, t, q):
    if t == T_SPHERE:
        dx = px - q[0]; dy = py - q[1]; dz = pz - q[2]
        return np.sqrt(dx * dx + dy * dy + dz * dz) - q[3]
    if t == T_CAPSULE:
        ax, ay, az = px - q[0], py - q[1], pz - q[2]
        bx, by, bz = q[3] - q[0], q[4] - q[1], q[5] - q[2]
        h = (ax * bx + ay * by + az * bz) / (bx * bx + by * by + bz * bz)
        h = min(max(h, 0.0), 1.0)
        dx = ax - bx * h; dy = ay - by * h; dz = az - bz * h
        return np.sqrt(dx * dx + dy * dy + dz * dz) - q[6]
    if t == T_PILL:   # points within q[7] of a flat disk (centre q0-2, normal q3-5, radius q6)
        vx, vy, vz = px - q[0], py - q[1], pz - q[2]
        z = vx * q[3] + vy * q[4] + vz * q[5]
        rx, ry, rz = vx - z * q[3], vy - z * q[4], vz - z * q[5]
        rho = np.sqrt(rx * rx + ry * ry + rz * rz)
        s = max(rho - q[6], 0.0)
        return np.sqrt(s * s + z * z) - q[7]
    # T_CERPATCH: spherical shell patch (axis q0-2, cos half-angle q3, half-angle q4, mid radius q5, half thick q6)
    r = np.sqrt(px * px + py * py + pz * pz)
    c = (px * q[0] + py * q[1] + pz * q[2]) / max(r, 1e-9)
    c = min(max(c, -1.0), 1.0)
    th = np.arccos(c)
    s = max(th - q[4], 0.0) * r
    dr = r - q[5]
    return np.sqrt(s * s + dr * dr) - q[6]


@nb.njit(cache=True, fastmath=True)
def mem_sdfs(px, py, pz, G, pores, pidx, PT, PQ, PC, out):
    """Fill out[NMEM] with signed distances to each membrane (negative = lumenal side).
    Returns the mito-matrix SDF (negative inside matrix)."""
    r = np.sqrt(px * px + py * py + pz * pz)
    out[M_PM] = G[2] - r
    out[M_OUTER] = G[0] - r
    # nuclear envelope lumen with toroidal pores
    vx, vy, vz = px - G[3], py - G[4], pz - G[5]
    dn = np.sqrt(vx * vx + vy * vy + vz * vz)
    hp = dn - (G[6] - 0.5 * PERI_W)
    d = abs(hp) - 0.5 * PERI_W
    rc = PORE_R + 0.5 * PERI_W
    for k in range(pores.shape[0]):
        z = vx * pores[k, 0] + vy * pores[k, 1] + vz * pores[k, 2]
        if z <= 0:
            continue
        rx, ry, rz = vx - z * pores[k, 0], vy - z * pores[k, 1], vz - z * pores[k, 2]
        rho = np.sqrt(rx * rx + ry * ry + rz * rz)
        if rho < rc:
            d = np.sqrt((rho - rc) ** 2 + hp * hp) - 0.5 * PERI_W
            break
    out[M_NE] = d
    for m in (M_ER, M_VAC, M_OM, M_GOLGI, M_PEROX, M_LD):
        out[m] = 1e9
    im = 1e9; cr = 1e9
    for ii in range(pidx.shape[0]):
        i = pidx[ii]
        d = sdf_prim(px, py, pz, PT[i], PQ[i])
        c = PC[i]
        if c == C_CER:
            out[M_ER] = min(out[M_ER], d)
        elif c == C_VAC:
            out[M_VAC] = min(out[M_VAC], d)
        elif c == C_OM:
            out[M_OM] = min(out[M_OM], d)
        elif c == C_IM:
            im = min(im, d)
        elif c == C_CRISTA:
            cr = min(cr, d)
        elif c == C_GOLGI:
            out[M_GOLGI] = min(out[M_GOLGI], d)
        elif c == C_PEROX:
            out[M_PEROX] = min(out[M_PEROX], d)
        elif c == C_LD:
            out[M_LD] = min(out[M_LD], d)
    matrix = max(im, -cr)
    out[M_IM] = -matrix
    return matrix


@nb.njit(cache=True, fastmath=True)
def label_from(px, py, pz, G, pores, out, matrix):
    """Compartment label given precomputed membrane SDFs (out) and matrix SDF."""
    r = np.sqrt(px * px + py * py + pz * pz)
    if r > G[0]:
        return EXT
    if r > G[1]:
        return MANNAN
    if r > G[2]:
        return GLUCAN
    if out[M_ER] < 0:
        return CER
    if out[M_NE] < 0:
        return PERI
    vx, vy, vz = px - G[3], py - G[4], pz - G[5]
    dn = np.sqrt(vx * vx + vy * vy + vz * vz)
    if dn < G[6] - 0.5 * PERI_W:
        hp = dn - (G[6] - 0.5 * PERI_W)
        rc = PORE_R + 0.5 * PERI_W
        if hp > -0.5 * PERI_W - 8.0:
            for k in range(pores.shape[0]):
                z = vx * pores[k, 0] + vy * pores[k, 1] + vz * pores[k, 2]
                if z <= 0:
                    continue
                rx, ry, rz = vx - z * pores[k, 0], vy - z * pores[k, 1], vz - z * pores[k, 2]
                if rx * rx + ry * ry + rz * rz < rc * rc:
                    return PORE
        wx, wy, wz = px - G[7], py - G[8], pz - G[9]
        if dn < G[6] - PERI_W and wx * wx + wy * wy + wz * wz < G[10] * G[10]:
            return NUCL
        return NUCP
    if dn < G[6]:      # pore channel between midplane and ONM
        return PORE
    if out[M_OM] < 0:
        return MATRIX if matrix < 0 else IMS
    if out[M_VAC] < 0:
        return VAC
    if out[M_LD] < 0:
        return LD
    if out[M_GOLGI] < 0:
        return GOLGI
    if out[M_PEROX] < 0:
        return PEROX
    return CYTO


@nb.njit(cache=True, fastmath=True)
def label_of(px, py, pz, G, pores, pidx, PT, PQ, PC, out):
    matrix = mem_sdfs(px, py, pz, G, pores, pidx, PT, PQ, PC, out)
    return label_from(px, py, pz, G, pores, out, matrix)


@nb.njit(cache=True, fastmath=True)
def min_mem_dist(out, skip):
    d = 1e9
    for m in range(NMEM):
        if m != skip:
            d = min(d, abs(out[m]))
    return d


@nb.njit(cache=True, parallel=True)
def label_grid(xs, ys, zs, G, pores, pidx, PT, PQ, PC):
    lab = np.zeros((xs.shape[0], ys.shape[0], zs.shape[0]), np.int8)
    for i in nb.prange(xs.shape[0]):
        out = np.zeros(NMEM)
        for j in range(ys.shape[0]):
            for k in range(zs.shape[0]):
                lab[i, j, k] = label_of(xs[i], ys[j], zs[k], G, pores, pidx, PT, PQ, PC, out)
    return lab


@nb.njit(cache=True, parallel=True)
def mem_grid(xs, ys, zs, G, pores, pidx, PT, PQ, PC):
    g = np.zeros((NMEM, xs.shape[0], ys.shape[0], zs.shape[0]), np.float32)
    for i in nb.prange(xs.shape[0]):
        out = np.zeros(NMEM)
        for j in range(ys.shape[0]):
            for k in range(zs.shape[0]):
                mem_sdfs(xs[i], ys[j], zs[k], G, pores, pidx, PT, PQ, PC, out)
                for m in range(NMEM):
                    g[m, i, j, k] = out[m]
    return g


@nb.njit(cache=True, parallel=True)
def label_points(P, G, pores, pidx, PT, PQ, PC):
    lab = np.zeros(P.shape[0], np.int8)
    for i in nb.prange(P.shape[0]):
        out = np.zeros(NMEM)
        lab[i] = label_of(P[i, 0], P[i, 1], P[i, 2], G, pores, pidx, PT, PQ, PC, out)
    return lab


# ----------------------------------------------------------------- python geometry container
class Cell:
    def __init__(self, d):
        self.d = d
        self.G = np.array(d["G"], np.float64)
        self.pores = np.array(d["pores"], np.float64).reshape(-1, 3)
        self.PT = np.array(d["PT"], np.int64)
        self.PQ = np.array(d["PQ"], np.float64).reshape(-1, 8)
        self.PC = np.array(d["PC"], np.int64)
        self.bc = np.array(d["bc"], np.float64).reshape(-1, 3)   # bounding spheres of prims
        self.br = np.array(d["br"], np.float64)
        self.npc_pos = np.array(d["npc_pos"], np.float64).reshape(-1, 3)

    def to_json(self, fn):
        json.dump(self.d, open(fn, "w"))

    @staticmethod
    def load(fn):
        return Cell(json.load(open(fn)))

    def select(self, lo, hi, margin=80.0):
        """Indices of prims / pores whose bounding sphere intersects box [lo,hi] grown by margin."""
        lo = np.asarray(lo) - margin; hi = np.asarray(hi) + margin
        q = np.clip(self.bc, lo, hi)
        pidx = np.where(((q - self.bc) ** 2).sum(1) <= self.br ** 2)[0].astype(np.int64)
        # pores: pore centre at NE midsurface
        pc = self.G[3:6] + self.pores * (self.G[6] - 0.5 * PERI_W)
        q = np.clip(pc, lo, hi)
        pk = np.where(((q - pc) ** 2).sum(1) <= (PORE_R + PERI_W + 40) ** 2)[0]
        return pidx, np.ascontiguousarray(self.pores[pk])

    def args(self, lo, hi, margin=80.0):
        pidx, pores = self.select(lo, hi, margin)
        return self.G, pores, pidx, self.PT, self.PQ, self.PC

    def label(self, P, margin=80.0):
        P = np.atleast_2d(P)
        lo, hi = P.min(0), P.max(0)
        a = self.args(lo, hi, margin)
        out = np.zeros(NMEM)
        return np.array([label_of(p[0], p[1], p[2], *a, out) for p in P], np.int8)


# ----------------------------------------------------------------- builder
def _rand_unit(rng, n=None):
    v = rng.normal(size=(n or 1, 3)); v /= np.linalg.norm(v, axis=1, keepdims=True)
    return v if n else v[0]


def build(geo, seed=11, log=print):
    rng = np.random.default_rng(seed)
    V = geo["volumes_um3"]
    Rcell = geo["cell_radius_um"] * 1000.0
    Rp = Rcell - geo["wall_thickness_um"] * 1000.0
    rn = geo["nucleus_radius_um"] * 1000.0
    rv = (V["vacuole"] * 3 / (4 * np.pi)) ** (1 / 3) * 1000.0
    cn = np.array([-990.0, 160.0, 90.0])
    cv = np.array([985.0, -170.0, -110.0])
    # nucleolus crescent: sphere centred beyond the NE on the side away from the vacuole, radius solved for
    # 20% of nuclear volume (MC)
    dirn = (cn - cv) / np.linalg.norm(cn - cv)
    cno = cn + dirn * rn * 0.80
    target = 0.2 * 4 / 3 * np.pi * rn ** 3
    pts = cn + _rand_unit(rng, 400000) * (rn - PERI_W) * rng.random((400000, 1)) ** (1 / 3)
    vin = 4 / 3 * np.pi * (rn - PERI_W) ** 3
    lo_, hi_ = 50.0, 2 * rn
    for _ in range(40):
        mid = 0.5 * (lo_ + hi_)
        f = (np.linalg.norm(pts - cno, axis=1) < mid).mean() * vin
        lo_, hi_ = (mid, hi_) if f < target else (lo_, mid)
    rno = 0.5 * (lo_ + hi_)
    G = [Rcell, Rcell - 35.0, Rp, *cn, rn, *cno, rno]
    PT, PQ, PC, bc, br = [], [], [], [], []

    def add(t, q, c, cen, rad):
        q = list(q) + [0.0] * (8 - len(q)); PT.append(t); PQ.append(q); PC.append(c); bc.append(list(cen)); br.append(rad)

    add(T_SPHERE, [*cv, rv], C_VAC, cv, rv)
    occupied = [(cn, rn), (cv, rv)]

    def free(c, r, gap):
        if np.linalg.norm(c) + r > Rp - CER_GAP - CER_LUMEN - 20 - gap:
            return False
        return all(np.linalg.norm(c - o) > r + ro + gap for o, ro in occupied)

    # ---- NPCs: 117 on the NE, >=150 nm apart (surface distance between centres)
    pores = []
    while len(pores) < geo["n_npc"]:
        u = _rand_unit(rng)
        if all(np.dot(u, w) < np.cos(150.0 / (rn - 17)) for w in pores):
            pores.append(u)
    pores = np.array(pores)
    npc_pos = cn + pores * (rn - 0.5 * PERI_W)

    # ---- mitochondrial network: cortical branched tubules, total length ~ V/(pi r^2)
    Vtarget = V["mitochondrion"] * 1e9
    sph = 4 / 3 * np.pi * MITO_R ** 3
    segs = []
    total = 0.0; vol = 0.0
    tries = 0
    while vol < Vtarget and tries < 5000:
        tries += 1
        parent = -1
        if segs and rng.random() < 0.4:   # branch from an existing segment
            parent = int(rng.integers(len(segs)))
            a, b = segs[parent]
            start = a + (b - a) * rng.uniform(0.3, 0.7)
        else:
            start = _rand_unit(rng) * rng.uniform(1560, 1700)
        d = _rand_unit(rng); d -= start * np.dot(d, start) / np.dot(start, start); d /= np.linalg.norm(d)
        p = start; nseg = 0
        for k in range(rng.integers(4, 12)):
            seglen = rng.uniform(250, 420)
            d = d + rng.normal(0, 0.3, 3); d -= p * np.dot(d, p) / np.dot(p, p); d /= np.linalg.norm(d)
            q = p + d * seglen
            q = q / np.linalg.norm(q) * np.clip(np.linalg.norm(q), 1520, 1720)
            samp = p + (q - p) * np.linspace(0.0, 1, 10)[:, None]
            okk = all(np.linalg.norm(s - cn) > rn + MITO_R + 40 and np.linalg.norm(s - cv) > rv + MITO_R + 40 for s in samp)
            ignore = {parent} if nseg == 0 else {len(segs) - 1}
            if okk and segs:
                A = np.array([s[0] for s in segs]); B = np.array([s[1] for s in segs])
                ab = B - A
                for s in samp[3:]:
                    h = np.clip(((s - A) * ab).sum(1) / (ab * ab).sum(1), 0, 1)
                    dd = np.linalg.norm(s - (A + ab * h[:, None]), axis=1)
                    dd[list(i for i in ignore if i >= 0)] = 1e9
                    if (dd < 2 * MITO_R + 60).any():
                        okk = False; break
            if not okk:
                break
            segs.append((p.copy(), q.copy())); total += np.linalg.norm(q - p); p = q; nseg += 1
            vol += np.pi * MITO_R ** 2 * np.linalg.norm(segs[-1][1] - segs[-1][0])
            if nseg == 1:
                vol += sph * (0.5 if parent < 0 else 0.0) + sph * 0.5   # free end(s) of the new chain/branch
            if vol >= Vtarget:
                break
    log(f"mito: {len(segs)} segments, length {total:.0f} nm, est. volume {vol/1e9:.3f} um3 (target {Vtarget/1e9:.3f})")
    n_cristae = 0
    for a, b in segs:
        cen = (a + b) / 2; rad = np.linalg.norm(b - a) / 2 + MITO_R
        add(T_CAPSULE, [*a, *b, MITO_R], C_OM, cen, rad)
        add(T_CAPSULE, [*a, *b, MITO_R - IMS_W], C_IM, cen, rad)
        ax = (b - a) / np.linalg.norm(b - a); Ls = np.linalg.norm(b - a)
        for s in np.arange(30.0, Ls - 30.0, 42.0):
            e = _rand_unit(rng); e -= ax * np.dot(e, ax); e /= np.linalg.norm(e)
            off = rng.uniform(30, 50); R = (MITO_R - IMS_W + 2.0) - off
            c = a + ax * (s + rng.uniform(-6, 6)) + e * off
            add(T_PILL, [*c, *ax, R, CRISTA_HALF], C_CRISTA, c, R + CRISTA_HALF + 1)
            n_cristae += 1
    for a, b in segs:
        for t in np.linspace(0, 1, max(2, int(np.linalg.norm(b - a) / 100))):
            occupied.append((a + (b - a) * t, MITO_R))
    log(f"cristae: {n_cristae}")

    def place_sphere(r, cls, gap=40.0, rmin=0, rmax=1900):
        for _ in range(20000):
            c = _rand_unit(rng) * rng.uniform(rmin, rmax)
            if free(c, r, gap):
                add(T_SPHERE, [*c, r], cls, c, r); occupied.append((c, r)); return c
        raise RuntimeError("could not place sphere")

    # ---- lipid droplets (often near the ER/NE; illustrative)
    nld = 10
    rld = (V["lipid_droplet"] * 1e9 * 3 / (4 * np.pi * nld)) ** (1 / 3)
    for _ in range(nld):
        place_sphere(rld, C_LD)
    # ---- peroxisomes
    npx = 4
    rpx = (V["peroxisome"] * 1e9 * 3 / (4 * np.pi * npx)) ** (1 / 3)
    for _ in range(npx):
        place_sphere(rpx, C_PEROX)
    # ---- Golgi cisternae (dispersed, unstacked in S. cerevisiae), endosomes, vesicles
    vg = V["golgi_endosome"] * 1e9
    n_cis, n_end, n_ves = 14, 8, 8
    vol_cis = 0.0
    for _ in range(n_cis):
        R = rng.uniform(170, 260); h = 10.0
        for _ in range(20000):
            c = _rand_unit(rng) * rng.uniform(300, 1850)
            if free(c, R + h, 40.0):
                break
        nrm = _rand_unit(rng)
        add(T_PILL, [*c, *nrm, R, h], C_GOLGI, c, R + h); occupied.append((c, R + h))
        vol_cis += 2 * np.pi * R * R * h + np.pi ** 2 * h * h * R + 4 / 3 * np.pi * h ** 3
    rves = 35.0
    vol_ves = n_ves * 4 / 3 * np.pi * rves ** 3
    rend = ((vg - vol_cis - vol_ves) * 3 / (4 * np.pi * n_end)) ** (1 / 3)
    for _ in range(n_end):
        place_sphere(rend, C_GOLGI)
    for _ in range(n_ves):
        place_sphere(rves, C_GOLGI, gap=30.0)
    log(f"LD r={rld:.0f}  perox r={rpx:.0f}  endosome r={rend:.0f}  golgi cisternae {n_cis}")

    # ---- cortical ER: spherical shell patches covering geo['cer_cover'] of the PM area
    rc = Rp - CER_GAP - 0.5 * CER_LUMEN
    test = _rand_unit(rng, 200000)
    covered = np.zeros(len(test), bool)
    caps = []
    while covered.mean() < geo["cer_cover"]:
        u = _rand_unit(rng); al = rng.uniform(0.12, 0.32)
        cand = test @ u > np.cos(al)
        # keep patches mostly separate (sheets with fenestrae-like gaps)
        if (cand & covered).mean() > 0.15 * cand.mean():
            continue
        caps.append((u, al)); covered |= cand
    for u, al in caps:
        add(T_CERPATCH, [*u, np.cos(al), al, rc, 0.5 * CER_LUMEN], C_CER, u * rc, al * rc + 20)
    log(f"cER: {len(caps)} patches, PM coverage {covered.mean():.3f}")

    d = dict(G=G, pores=pores.ravel().tolist(), PT=PT, PQ=np.array(PQ).ravel().tolist(), PC=PC,
             bc=np.array(bc).ravel().tolist(), br=br, npc_pos=npc_pos.ravel().tolist(),
             mito_segments=[[list(a), list(b)] for a, b in segs], seed=seed,
             consts=dict(HALF=HALF, PORE_R=PORE_R, PERI_W=PERI_W, IMS_W=IMS_W, MITO_R=MITO_R,
                         CER_GAP=CER_GAP, CER_LUMEN=CER_LUMEN),
             label_names=LABEL_NAMES, mem_names=MEM_NAMES)
    return Cell(d)
