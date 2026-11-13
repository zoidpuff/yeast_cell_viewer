"""3D whole-cell geometry for a haploid S. cerevisiae cell as analytic signed-distance primitives.
Units nm, origin at cell centre. Organelle volumes/areas are data-driven (inventory geo dict); nucleus shape, vacuole
lobes and lipid droplets match statistics of segmented cryo-FIB/SEM cells (organelle_stats.json); the arrangement
is illustrative (no single tomogram is copied).

G (cell constants): [0] cell outer radius, [1] mannan/glucan boundary, [2] plasma membrane, [3:6] nucleus centre,
[6] nucleus equivalent-sphere radius (ONM midsurface), [7:10] nucleolus sphere centre, [10] its radius,
[11:14] nucleus semi-axes (ONM midsurface), [14:23] nucleus rotation rows (local = R (p - centre)).
Pores: per pore the centre on the perinuclear midsurface and the outward NE normal (6 values).

Real shapes (cellshape.py): the cell outline, the nucleus and the vacuoles can instead be signed-distance fields,
appended to G after G_N: G[G_N] = number of fields, then 12 header values per field (grid origin xyz, spacing, bricks
along x/y/z, voxels per brick, clamp inside, clamp outside, offset of the brick index, offset of the brick data) and
the index and data sections (field_sdf). With a cell field, G[0] is the radius of a sphere enclosing the cell and
G[0] - G[1], G[0] - G[2] are the mannan-layer and wall depths below the outer surface.

Membrane SDF sign convention: negative on the *lumenal / non-cytoplasmic* side; the gradient therefore
points to the cytoplasmic (N-) side, which is +z in the membrane-protein library frame.
For the mitochondrial inner membrane the matrix is the N-side, so SDF_IM = -SDF_matrix."""
import json
import numpy as np
from config import seed as _seed
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
T_SPHERE, T_CAPSULE, T_PILL, T_CERPATCH, T_SLOT, T_TORUS = range(6)
C_VAC, C_LD, C_GOLGI, C_PEROX, C_OM, C_IM, C_CRISTA, C_CER, C_PMINV, C_SCAR = range(10)   # PMINV: PM invaginations (T_SLOT)
HALF = 2.0           # bilayer half thickness used for exclusion (as in the column model)
PORE_R = 45.0        # pore-membrane radius at the NE midplane
PERI_W = 34.0        # perinuclear space (ONM-INM midsurface distance)
# Mitochondria of a glucose-grown cell (2026-11-13; numbers and files: data/validation/literature/mitochondria/
# README_mitochondria_2026-11-13.md). Until then: straight 300 nm tubules, 8 nm between the membranes, and a disc
# across the whole tubule every 42 nm (crista membrane 1.4 x the outer membrane, the look of a respiring cell).
IMS_W = 14.0         # outer to inner boundary membrane, midsurfaces: 13-16 nm in cryo-tomograms of BY4741 in YPD
                     # (Chang et al. 2025 J Cell Biol, EMPIAR-12534; was 8)
MITO_R = 168.0       # outer radius: diameter 339 +- 5 nm in glucose (live 4Pi, Egner et al. 2002 PNAS 99:3370); median 330
                     # in the Chang 2025 meshes; was 150
MITO_R_SD, MITO_R_WALK, MITO_R_MIN, MITO_R_MAX = 18.0, 7.0, 125.0, 215.0   # it varies along a tubule (230-460 nm across in the meshes)
MITO_STEP = 110.0    # nm per capsule segment: short, so that tubules curve
MITO_TURN = 0.17     # direction noise per segment (rad)
MITO_IGNORE = 8      # a growing tubule does not test against its own last segments
NVJ_CONTACT = 45.0   # nm: the vacuolar membrane within this distance of the outer nuclear membrane is junction (as nvj.CONTACT)
NVJ_PORE_MARGIN = 30.0   # nm: no pore centre this close to the junction (a pore's membrane torus is 45 nm in radius)
CRISTA_HALF = 9.0    # crista membrane half-separation = radius of the rim tube (fits ATP synthase dimers, see atp_dimers.py)
CRISTA_PER_UM = 34.0 # cristae per um in the stretches that bear cristae; with the bare stretches around nucleoids and
                     # branch points about 18-20 per um of tubule overall (19-24 crista patches per um2 of outer
                     # membrane, i.e. about 20-25 per um, in the YPD meshes)
CRISTA_GAP_MIN = 22.0
# how far a crista reaches into the tubule: 60-110 nm, at most about 200, never across it ("almost none of them spanned
# the mitochondria", Stoldt et al. 2019 PNAS 116:9853; glucose cells have "relatively few, poorly defined cristae",
# Klecker & Westermann 2021). With these numbers the crista membrane is 0.2-0.3 x the outer membrane, as measured on
# the YPD meshes (0.72 in glycerol)
CRISTA_REACH, CRISTA_REACH_SHAPE, CRISTA_REACH_MIN, CRISTA_REACH_MAX = 85.0, 9.0, 50.0, 200.0
NUCLEOID_PER_UM = 1.27     # nucleoids per um of tubule (1.2-1.35; Osman et al. 2015 PNAS; Goke et al. 2020)
NUCLEOID_MIN_GAP = 550.0   # nm between nucleoids (mean spacing about 0.8 um, Osman 2015)
NUCLEOID_CLEAR = 110.0     # no crista within this distance of a nucleoid site: nucleoids are compact bodies (up to
                           # 0.3-0.4 um long, Goke 2020; one or two genomes are drawn about 130 nm across here) in
                           # crista-free matrix (placement after Kukat 2015, mammalian)
CER_GAP = 27.0       # PM midsurface -> outer cER membrane midsurface
CER_LUMEN = 20.0     # cER membrane midsurface separation
# eisosome furrows (MCC): ~30 per cell, ~300 nm long, ~50 nm deep (Stradalova 2009, Karotki 2011); half-pipe bottom with
# the membrane midplane radius of the native Pil1/Lsp1 coat (8QB7, eisosomes.py); EIS_K: smooth-union radius at the rim
EIS_N, EIS_LEN, EIS_DEPTH, EIS_R, EIS_K = 30, (300.0, 50.0), 50.0, 10.2, 8.0
# endocytic invaginations (one per actin patch): tubes with a rounded tip, growing to ~100 nm before scission (Kukulski
# et al. 2012 Cell 150:508); radius from the vesicles they release (6,400 nm2 mean area = 45 nm sphere; Kukulski 2016)
ENDO_R, ENDO_DEPTH = 20.0, (0.0, 100.0)         # depth uniform over the invagination phase (assumption)
ENDO_LEN = 0.01                                  # T_SLOT with a (near) zero-length bottom: a tube along the radial direction
# bud scar: chitin ring left on the mother's wall by each division, ~1 um across (AFM: Pillet et al. 2014); a torus
# (major radius SCAR_R, tube radius SCAR_r) whose tube centre lies on the outer surface, so the rim stands SCAR_r proud
# of it (height: assumption); smooth union with the outer surface only (the PM and the wall layers keep their depths)
BUD_SCARS, SCAR_R, SCAR_r, SCAR_K = 1, 450.0, 50.0, 30.0
SCAR_CHITIN = 0.5      # share of the cell's chitin in its bud scar(s) (assumption: "mostly in bud scars")


G_N = 23                 # analytic constants in G; distance fields follow
F_CELL, F_VAC, F_NUC = 0, 1, 2


# ----------------------------------------------------------------- numba evaluators
@nb.njit(cache=True, fastmath=True)
def n_fields(G):
    return int(G[G_N]) if G.shape[0] > G_N else 0


@nb.njit(cache=True, fastmath=True)
def field_sdf(px, py, pz, G, k):
    """signed distance field k (negative inside): trilinear within a stored brick, +-clamp elsewhere"""
    o = G_N + 1 + 12 * k
    h = G[o + 3]; B = int(G[o + 7])
    fx = (px - G[o]) / h; fy = (py - G[o + 1]) / h; fz = (pz - G[o + 2]) / h
    nbx = int(G[o + 4]); nby = int(G[o + 5]); nbz = int(G[o + 6])
    if fx < 0.0 or fy < 0.0 or fz < 0.0 or fx >= nbx * B or fy >= nby * B or fz >= nbz * B:
        return G[o + 9]
    ix = int(fx); iy = int(fy); iz = int(fz)
    bx = ix // B; by = iy // B; bz = iz // B
    bi = int(G[int(G[o + 10]) + (bx * nby + by) * nbz + bz])
    if bi == -1:
        return -G[o + 8]
    if bi == -2:
        return G[o + 9]
    n1 = B + 1
    tx = fx - ix; ty = fy - iy; tz = fz - iz
    b = int(G[o + 11]) + bi * n1 * n1 * n1 + ((ix - bx * B) * n1 + (iy - by * B)) * n1 + (iz - bz * B)
    sx = n1 * n1; sy = n1
    c00 = G[b] + (G[b + sx] - G[b]) * tx
    c01 = G[b + 1] + (G[b + sx + 1] - G[b + 1]) * tx
    c10 = G[b + sy] + (G[b + sx + sy] - G[b + sy]) * tx
    c11 = G[b + sy + 1] + (G[b + sx + sy + 1] - G[b + sy + 1]) * tx
    c0 = c00 + (c10 - c00) * ty; c1 = c01 + (c11 - c01) * ty
    return c0 + (c1 - c0) * tz


@nb.njit(cache=True, fastmath=True)
def outer_depth(px, py, pz, G):
    """depth below the outer cell surface (negative outside the cell)"""
    if n_fields(G) > F_CELL:
        return -field_sdf(px, py, pz, G, F_CELL)
    return G[0] - np.sqrt(px * px + py * py + pz * pz)


@nb.njit(cache=True, fastmath=True)
def nuc_sdf(px, py, pz, G):
    if n_fields(G) > F_NUC:
        return field_sdf(px, py, pz, G, F_NUC)
    return nuc_ellipsoid(px, py, pz, G)


@nb.njit(cache=True, fastmath=True)
def nuc_ellipsoid(px, py, pz, G):
    """signed distance to the ONM midsurface (ellipsoid; negative inside), first-order approximation
    k0 (k0 - 1) / k1: exact on the surface, error ~ d^2 x curvature difference nearby (< 0.5 nm across the NE)"""
    dx = px - G[3]; dy = py - G[4]; dz = pz - G[5]
    lx = (G[14] * dx + G[15] * dy + G[16] * dz) / G[11]
    ly = (G[17] * dx + G[18] * dy + G[19] * dz) / G[12]
    lz = (G[20] * dx + G[21] * dy + G[22] * dz) / G[13]
    k0 = np.sqrt(lx * lx + ly * ly + lz * lz)
    k1 = np.sqrt((lx / G[11]) ** 2 + (ly / G[12]) ** 2 + (lz / G[13]) ** 2)
    if k1 < 1e-12:
        return -min(G[11], min(G[12], G[13]))
    return k0 * (k0 - 1.0) / k1


@nb.njit(cache=True, fastmath=True)
def in_pore(px, py, pz, pores, rc):
    """index of the pore whose axis (through its centre along the NE normal) is within rc of p, or -1"""
    for k in range(pores.shape[0]):
        vx = px - pores[k, 0]; vy = py - pores[k, 1]; vz = pz - pores[k, 2]
        z = vx * pores[k, 3] + vy * pores[k, 4] + vz * pores[k, 5]
        if z < -150.0 or z > 150.0:
            continue
        rx = vx - z * pores[k, 3]; ry = vy - z * pores[k, 4]; rz = vz - z * pores[k, 5]
        if rx * rx + ry * ry + rz * rz < rc * rc:
            return k
    return -1


@nb.njit(cache=True, fastmath=True)
def sdf_cer(px, py, pz, q, pm, Rp):
    """cortical-ER patch: a sheet at constant depth (Rp - q[5]) below the plasma membrane (pm: signed PM distance,
    positive in the cytoplasm), within the angular cap q[0:3], cos q[3], q[4]; half-thickness q[6]"""
    r = np.sqrt(px * px + py * py + pz * pz)
    c = (px * q[0] + py * q[1] + pz * q[2]) / max(r, 1e-9)
    c = min(1.0, max(-1.0, c))
    s = max(np.arccos(c) - q[4], 0.0) * r
    dr = (Rp - q[5]) - pm
    return np.sqrt(s * s + dr * dr) - q[6]


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
    if t == T_TORUS:  # bud-scar rim: centre q0-2, axis q3-5, major radius q6, tube radius q7
        vx, vy, vz = px - q[0], py - q[1], pz - q[2]
        h = vx * q[3] + vy * q[4] + vz * q[5]
        rho = np.sqrt(max(vx * vx + vy * vy + vz * vz - h * h, 0.0))
        return np.sqrt((rho - q[6]) ** 2 + h * h) - q[7]
    if t == T_SLOT:   # eisosome furrow lumen: within q[6] of the rectangle swept from the bottom axis q0-2 -> q3-5 by
        # q[7] along the outward radial direction of its midpoint (axis perpendicular to it)
        mx, my, mz = 0.5 * (q[0] + q[3]), 0.5 * (q[1] + q[4]), 0.5 * (q[2] + q[5])
        m = np.sqrt(mx * mx + my * my + mz * mz); ux, uy, uz = mx / m, my / m, mz / m
        ex, ey, ez = q[3] - q[0], q[4] - q[1], q[5] - q[2]
        L = np.sqrt(ex * ex + ey * ey + ez * ez); ex /= L; ey /= L; ez /= L
        vx, vy, vz = px - q[0], py - q[1], pz - q[2]
        s = min(max(vx * ex + vy * ey + vz * ez, 0.0), L); h = min(max(vx * ux + vy * uy + vz * uz, 0.0), q[7])
        dx, dy, dz = vx - s * ex - h * ux, vy - s * ey - h * uy, vz - s * ez - h * uz
        return np.sqrt(dx * dx + dy * dy + dz * dz) - q[6]
    # T_CERPATCH: spherical shell patch (axis q0-2, cos half-angle q3, half-angle q4, mid radius q5, half thick q6)
    r = np.sqrt(px * px + py * py + pz * pz)
    c = (px * q[0] + py * q[1] + pz * q[2]) / max(r, 1e-9)
    c = min(max(c, -1.0), 1.0)
    th = np.arccos(c)
    s = max(th - q[4], 0.0) * r
    dr = r - q[5]
    return np.sqrt(s * s + dr * dr) - q[6]


@nb.njit(cache=True, fastmath=True)
def smin(a, b, k):
    """polynomial smooth minimum (rounds the union of two distance fields over ~k)"""
    h = max(k - abs(a - b), 0.0) / k
    return min(a, b) - h * h * k * 0.25


@nb.njit(cache=True, fastmath=True)
def mem_sdfs(px, py, pz, G, pores, pidx, PT, PQ, PC, out):
    """Fill out[NMEM] with signed distances to each membrane (negative = lumenal side).
    Returns the mito-matrix SDF (negative inside matrix)."""
    dep = outer_depth(px, py, pz, G)
    out[M_PM] = dep - (G[0] - G[2])
    out[M_OUTER] = dep
    # nuclear envelope lumen with toroidal pores (hp: offset from the perinuclear midsurface)
    hp = nuc_sdf(px, py, pz, G) + 0.5 * PERI_W
    d = abs(hp) - 0.5 * PERI_W
    rc = PORE_R + 0.5 * PERI_W
    if abs(hp) < rc + PERI_W:
        k = in_pore(px, py, pz, pores, rc)
        if k >= 0:
            vx = px - pores[k, 0]; vy = py - pores[k, 1]; vz = pz - pores[k, 2]
            z = vx * pores[k, 3] + vy * pores[k, 4] + vz * pores[k, 5]
            rho = np.sqrt(max(vx * vx + vy * vy + vz * vz - z * z, 0.0))
            d = np.sqrt((rho - rc) ** 2 + hp * hp) - 0.5 * PERI_W
    out[M_NE] = d
    for m in (M_ER, M_VAC, M_OM, M_GOLGI, M_PEROX, M_LD):
        out[m] = 1e9
    if n_fields(G) > F_VAC:
        out[M_VAC] = field_sdf(px, py, pz, G, F_VAC)
    im = 1e9; cr = 1e9; eis = 1e9; scar = 1e9
    for ii in range(pidx.shape[0]):
        i = pidx[ii]
        if PT[i] == T_CERPATCH:
            d = sdf_cer(px, py, pz, PQ[i], out[M_PM], G[2])
        else:
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
        elif c == C_PMINV:
            eis = min(eis, d)
        elif c == C_SCAR:
            scar = min(scar, d)
    if scar < -out[M_OUTER] + SCAR_K:  # bud scars: raised rims on the outer surface only
        out[M_OUTER] = -smin(-out[M_OUTER], scar, SCAR_K)
    if eis < out[M_PM] + EIS_K:        # PM invaginations: smooth union of their lumen with the outside of the PM
        out[M_PM] = smin(out[M_PM], eis, EIS_K)
    matrix = max(im, -cr)
    out[M_IM] = -matrix
    return matrix


@nb.njit(cache=True, fastmath=True)
def label_from(px, py, pz, G, pores, out, matrix):
    """Compartment label given precomputed membrane SDFs (out) and matrix SDF."""
    if out[M_OUTER] < 0:
        return EXT
    if out[M_OUTER] < G[0] - G[1]:
        return MANNAN
    if out[M_PM] < 0:
        return GLUCAN
    if out[M_ER] < 0:
        return CER
    if out[M_NE] < 0:
        return PERI
    ns = nuc_sdf(px, py, pz, G)
    if ns < -0.5 * PERI_W:
        if ns > -PERI_W - 8.0:
            if in_pore(px, py, pz, pores, PORE_R + 0.5 * PERI_W) >= 0:
                return PORE
        wx, wy, wz = px - G[7], py - G[8], pz - G[9]
        if ns < -PERI_W and wx * wx + wy * wy + wz * wz < G[10] * G[10]:
            return NUCL
        return NUCP
    if ns < 0:         # pore channel between midplane and ONM
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
def pack_fields(G, fields):
    """G (analytic constants) with distance fields appended; fields: list of cellshape.bricks() outputs, in the order
    F_CELL, F_VAC, F_NUC"""
    G = list(np.asarray(G, np.float64)[:G_N])
    head = len(G) + 1 + 12 * len(fields)
    secs, off = [], head
    hdr = []
    for h, idx, data in fields:
        hdr += list(h) + [off, off + len(idx)]
        secs += [idx, data]; off += len(idx) + len(data)
    return np.concatenate([np.array(G + [len(fields)] + hdr, np.float64)] + secs)


class Cell:
    def __init__(self, d, fields=None):
        self.d = d
        self.G = np.array(d["G"], np.float64)
        self.fields = fields                      # list of (header, index, data) or None (analytic cell)
        if fields:
            self.G = pack_fields(self.G, fields)
        self.pores = np.array(d["pores"], np.float64).reshape(-1, 6)    # centre (perinuclear midsurface), normal
        self.PT = np.array(d["PT"], np.int64)
        self.PQ = np.array(d["PQ"], np.float64).reshape(-1, 8)
        self.PC = np.array(d["PC"], np.int64)
        self.bc = np.array(d["bc"], np.float64).reshape(-1, 3)   # bounding spheres of prims
        self.br = np.array(d["br"], np.float64)
        self.npc_pos = np.array(d["npc_pos"], np.float64).reshape(-1, 3)

    def furrows(self):
        """prim indices of the eisosome furrows (T_SLOT with a long bottom)"""
        return [i for i in np.where(self.PT == T_SLOT)[0] if np.linalg.norm(self.PQ[i, 3:6] - self.PQ[i, :3]) > 1.0]

    def endocytic_sites(self):
        """(point 17.5 nm inside the PM on the invagination axis, outward radial direction, tip depth) per endocytic tube"""
        out = []
        for i in np.where(self.PT == T_SLOT)[0]:
            q = self.PQ[i]
            if np.linalg.norm(q[3:6] - q[:3]) <= 1.0:
                u = q[:3] / np.linalg.norm(q[:3]); mouth = u * ray_depth(self.G, u, self.G[0] - self.G[2])
                out.append((mouth - 17.5 * u, u, float(np.linalg.norm(mouth) - np.linalg.norm(q[:3]) + q[6])))
        return out

    def to_json(self, fn):
        import os
        d = dict(self.d)
        if self.fields:
            ff = os.path.splitext(fn)[0] + "_fields.npz"
            np.savez_compressed(ff, **{f"{n}{k}": a for k, (h, i, x) in enumerate(self.fields) for n, a in (("h", h), ("i", i), ("x", x))})
            d["fields"] = os.path.basename(ff)
        json.dump(d, open(fn, "w"))

    @staticmethod
    def load(fn):
        import os
        d = json.load(open(fn))
        fields = None
        if d.get("fields"):
            z = np.load(os.path.join(os.path.dirname(fn), d["fields"]))
            fields = [(z[f"h{k}"], z[f"i{k}"], z[f"x{k}"]) for k in range(len([f for f in z.files if f[0] == "h"]))]
        return Cell(d, fields)

    def select(self, lo, hi, margin=80.0):
        """Indices of prims / pores whose bounding sphere intersects box [lo,hi] grown by margin."""
        lo = np.asarray(lo) - margin; hi = np.asarray(hi) + margin
        q = np.clip(self.bc, lo, hi)
        pidx = np.where(((q - self.bc) ** 2).sum(1) <= self.br ** 2)[0].astype(np.int64)
        # pores: pore centre at NE midsurface
        pc = self.pores[:, :3]
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


@nb.njit(cache=True)
def _many(P, G, which):
    out = np.empty(P.shape[0])
    for n in range(P.shape[0]):
        if which == 0:
            out[n] = nuc_sdf(P[n, 0], P[n, 1], P[n, 2], G)
        else:
            out[n] = outer_depth(P[n, 0], P[n, 1], P[n, 2], G)
    return out


def nuc_sdf_np(P, G):
    """nuc_sdf for an array of points (n x 3)"""
    return _many(np.atleast_2d(np.asarray(P, np.float64)), np.asarray(G, np.float64), 0)


def ray_depth(G, u, depth, depth_fn=None):
    """radius along the unit direction u from the cell centre at which the depth below the outer cell surface is
    depth (exact for a spherical cell). The cell field in G is clamped (cellshape.BANDS); for deeper points pass
    depth_fn (e.g. the builder's dense grid)."""
    G = np.asarray(G, np.float64); u = np.asarray(u, np.float64) / np.linalg.norm(u)
    if n_fields(G) <= F_CELL:
        return G[0] - depth
    if depth_fn is None:
        assert depth < G[G_N + 1 + 12 * F_CELL + 8] - 5.0, "depth beyond the stored band: pass depth_fn"
        depth_fn = lambda p: outer_depth(p[0], p[1], p[2], G)
    lo_, hi_ = 0.0, G[0] + 50.0
    for _ in range(40):
        m = 0.5 * (lo_ + hi_)
        if float(np.ravel(depth_fn(m * u))[0]) > depth:
            lo_ = m
        else:
            hi_ = m
    return 0.5 * (lo_ + hi_)


def _grad(f, p, e=1.0):
    g = np.array([f(*(p + e * u)) - f(*(p - e * u)) for u in np.eye(3)])
    return g / max(np.linalg.norm(g), 1e-12)


def nuc_surface(G, u, offset=0.0):
    """point on the nuclear envelope along direction u from the nucleus centre, moved by offset (nm) along the
    outward normal (e.g. -PERI_W / 2: perinuclear midsurface), and that normal"""
    G = np.asarray(G, np.float64)
    if n_fields(G) > F_NUC:
        u = np.asarray(u, np.float64) / np.linalg.norm(u); c = G[3:6]
        lo_, hi_ = 0.0, 3000.0
        for _ in range(40):
            m = 0.5 * (lo_ + hi_)
            if nuc_sdf(*(c + m * u), G) < 0:
                lo_ = m
            else:
                hi_ = m
        p = c + 0.5 * (lo_ + hi_) * u
        n = _grad(lambda x, y, z: nuc_sdf(x, y, z, G), p)
        return p + offset * n, n
    R = G[14:23].reshape(3, 3); ax = G[11:14]
    ul = R @ (np.asarray(u, np.float64) / np.linalg.norm(u))
    pl = ul / np.sqrt(((ul / ax) ** 2).sum())                 # on the ellipsoid
    nl = pl / ax ** 2; nl /= np.linalg.norm(nl)
    n = R.T @ nl
    return G[3:6] + R.T @ pl + offset * n, n


def _sample_volumes(rng, pooled, target, jitter=0.25):
    """volumes drawn from the pooled real ones (with lognormal smoothing) until they add up to target; the last
    draw is rescaled so the sum is exact. Largest first."""
    out = []
    while sum(out) < target:
        out.append(float(rng.choice(pooled)) * float(np.exp(rng.normal(0.0, jitter))))
    out[-1] -= sum(out) - target
    if out[-1] < 0.25 * min(pooled):
        out.pop(); out = [v * target / sum(out) for v in out]
    return sorted(out, reverse=True)


def organelle_stats(fn=None):
    import os
    fn = fn or os.path.join(os.path.dirname(os.path.abspath(__file__)), "organelle_stats.json")
    st = json.load(open(fn))["cells"]
    frac = lambda key: np.concatenate([np.array(c[key]) / c["cell_volume_um3"] for c in st.values()])
    return dict(vacuole=frac("vacuole_volumes_um3"), lipid_droplet=frac("lipid_droplet_volumes_um3"),
                vac_rms_rel=float(np.median([c["vacuole_cluster_rms_um"] / (3 * c["cell_volume_um3"] / (4 * np.pi)) ** (1 / 3)
                                             for c in st.values()])))


NUC_AXES = (1.25, 1.10, 1.0)    # nucleus axis ratios: cell 05's nucleus (1.24 : 1.10 : 1) and the in-plane ratio of cell
                                # 04's (1.27; its short axis lies along the milling direction of a truncated cell)
NUC_OFFSET = 0.30               # nucleus centre offset / protoplast radius (two real nuclei: 0.13, 0.29)
VAC_OFFSET = 0.30               # vacuole cluster centre offset, opposite the nucleus (five cells: 0.04-0.33)
VAC_GAP = (30.0, 60.0)          # nm between touching vacuole lobes (real lobes are pressed together)


def build(geo, seed=11, log=print, shapes=None, n_endo=0, n_scars=0):
    """shapes: cellshape.build() output (real cell outline, nucleus and vacuoles as distance fields); None: an
    analytic cell (sphere, ellipsoidal nucleus, vacuole lobes resampled from the segmented cells); n_endo: endocytic
    invaginations (actin patch sites); n_scars: bud scars"""
    seed = _seed(seed); rng = np.random.default_rng(seed)
    V = geo["volumes_um3"]
    Rcell = geo["cell_radius_um"] * 1000.0
    rn = geo["nucleus_radius_um"] * 1000.0
    WALL = geo["wall_thickness_um"] * 1000.0
    stats = organelle_stats()
    PT, PQ, PC, bc, br = [], [], [], [], []

    def add(t, q, c, cen, rad):
        q = list(q) + [0.0] * (8 - len(q)); PT.append(t); PQ.append(q); PC.append(c); bc.append(list(cen)); br.append(rad)

    occupied = []
    fields = None
    if shapes is None:
        Rp = Rcell - WALL
        side = np.array([-990.0, 160.0, 90.0]); side /= np.linalg.norm(side)
        cn = side * NUC_OFFSET * Rp
        cv0 = -side * VAC_OFFSET * Rp
        # nucleus: ellipsoid of the inventory volume; long axis tangential, middle axis radial
        ax = np.array(NUC_AXES) * rn / np.prod(NUC_AXES) ** (1 / 3)
        e1 = np.cross(side, _rand_unit(rng)); e1 /= np.linalg.norm(e1); e3 = np.cross(e1, side)
        Rn = np.array([e1, side, e3])                       # rows: local x (long), y (middle, radial), z (short)
        G = [Rcell, Rcell - 35.0, Rp, *cn, rn, 0.0, 0.0, 0.0, 0.0, *ax, *Rn.ravel()]
        Gq = np.array(G)
        nuc_d = lambda P: nuc_sdf_np(P, Gq)
        depth = lambda P: Rcell - np.linalg.norm(np.atleast_2d(P), axis=1)
        nuc_ext = ax.max()
        nol_dir = -side                                     # towards the vacuole cluster
    else:
        cellg, nucg, vacg = shapes["cell"], shapes["nucleus"], shapes["vacuole"]
        Rcell = float(np.linalg.norm(shapes["meshes"]["cell"][0], axis=1).max()) + 20.0     # encloses the cell
        Rp = Rcell - WALL
        cn = np.asarray(shapes["nuc_centre"], np.float64); nuc_ext = shapes["nuc_extent"] + 40.0
        G = [Rcell, Rcell - 35.0, Rp, *cn, rn, 0.0, 0.0, 0.0, 0.0, nuc_ext, nuc_ext, nuc_ext, *np.eye(3).ravel()]
        fields = [shapes["bricks"]["cell"], shapes["bricks"]["vacuole"], shapes["bricks"]["nucleus"]]
        Gq = pack_fields(G, fields)
        nuc_d = nucg
        depth = lambda P: -cellg(P)
        vin_ = np.argwhere(vacg.sdf < 0)
        import cellshape
        cv = shapes["lo"] + vin_.mean(0) * cellshape.H              # vacuole centroid
        nol_dir = (cv - cn) / np.linalg.norm(cv - cn)
    # nucleolus crescent: sphere centred beyond the NE on the side facing the vacuole, radius solved for 20% of nuclear
    # volume (MC over the nucleoplasm). The nucleolus usually lies at the envelope next to the vacuole (Barbosa 2015 Mol
    # Biol Cell, citing Kvam & Goldfarb 2007; Santana-Sosa 2023 Int J Mol Sci); the spindle pole body, opposite the
    # nucleolus (Yang 1989 Chromosoma), then faces open cytoplasm. Up to build R the crescent was on the far side, which
    # with one large vacuole put the spindle pole body against the vacuole and left no room for cytoplasmic microtubules.
    cno = cn + 0.8 * (nuc_surface(Gq, nol_dir)[0] - cn)
    box = rng.uniform(-1, 1, (600000, 3)) * nuc_ext + cn
    nb_ = nuc_d(box)
    inn = box[nb_ < -PERI_W]
    vin = (nb_ < -PERI_W).mean() * (2 * nuc_ext) ** 3
    target = 0.2 * 4 / 3 * np.pi * rn ** 3
    lo_, hi_ = 50.0, 2 * rn
    for _ in range(40):
        mid = 0.5 * (lo_ + hi_)
        f = (np.linalg.norm(inn - cno, axis=1) < mid).mean() * vin
        lo_, hi_ = (mid, hi_) if f < target else (lo_, mid)
    rno = 0.5 * (lo_ + hi_)
    G[7:11] = [*cno, rno]
    Gq[7:11] = G[7:11]

    def free(c, r, gap):
        if depth(c)[0] - WALL < r + gap + CER_GAP + CER_LUMEN + 20:
            return False
        if nuc_d(c)[0] < r + gap:
            return False
        if vac_d(c)[0] < r + gap:
            return False
        return all(np.linalg.norm(c - o) > r + ro + gap for o, ro in occupied)

    if shapes is None:
        vac_d = lambda P: np.full(len(np.atleast_2d(P)), 1e9)
        # ---- vacuole: lobes with the size distribution of the segmented cells, grown as touching clusters around a
        # centre opposite the nucleus (real lobes are pressed together and spread over about half the cell radius)
        vols = _sample_volumes(rng, stats["vacuole"] * geo["cell_volume_um3"] * 1e9, V["vacuole"] * 1e9)
        lobes = []
        for v in vols:
            r = (3 * v / (4 * np.pi)) ** (1 / 3)
            for att in range(40000):
                if not lobes:      # the largest lobe: as close to the cluster centre as it fits
                    c = cv0 + _rand_unit(rng) * rng.uniform(0, 0.6 * Rp) * min(1.0, att / 2000)
                elif rng.random() < 0.8:
                    cj, rj = lobes[int(rng.integers(len(lobes)))]
                    c = cj + _rand_unit(rng) * (r + rj + rng.uniform(*VAC_GAP))
                else:
                    c = cv0 + rng.normal(0, 0.35 * Rp, 3)
                if free(c, r, VAC_GAP[0]):
                    break
            else:
                raise RuntimeError("could not place a vacuole lobe")
            lobes.append((c, r)); occupied.append((c, r)); add(T_SPHERE, [*c, r], C_VAC, c, r)
        lr = np.array([r for _, r in lobes]); lv = 4 / 3 * np.pi * lr ** 3; lc = np.array([c for c, _ in lobes])
        vcen = (lc * lv[:, None]).sum(0) / lv.sum()
        rms = np.sqrt((lv * ((lc - vcen) ** 2).sum(1)).sum() / lv.sum())
        log(f"vacuole: {len(lobes)} lobes, largest r {lobes[0][1]:.0f} nm, cluster rms {rms / Rp:.2f} x protoplast radius "
            f"(segmented cells {stats['vac_rms_rel']:.2f})")
        cv = vcen
        vac_d = lambda P: (np.linalg.norm(np.atleast_2d(P)[:, None] - lc[None], axis=2) - lr[None]).min(1)
    else:
        vac_d = vacg

    # ---- NPCs on the NE, >=150 nm apart, uniform over the envelope area outside the nucleus-vacuole junction (pores are
    # excluded from the junction: Severs 1976 and Pan 2000, as cited by Kohler & Buettner 2021 Contact 4 and by
    # Otto & Thumm 2020 Int J Mol Sci 21:4506)
    pores = []
    nvj_clear = 0.5 * PERI_W + NVJ_CONTACT + NVJ_PORE_MARGIN
    if shapes is None:
        amax = None
        while len(pores) < geo["n_npc"]:
            s, n = nuc_surface(G, _rand_unit(rng), -0.5 * PERI_W)
            # area element along a random direction ~ |p|^3 / (p . n) relative; accept proportionally
            p = s - cn; w = np.linalg.norm(p) ** 3 / max(np.dot(p, n), 1e-9)
            amax = amax or 1.3 * (ax.max() ** 3 / ax.min())
            if rng.random() * amax > w:
                continue
            if vac_d(s)[0] > nvj_clear and all(np.linalg.norm(s - q[:3]) > 150.0 for q in pores):
                pores.append(np.r_[s, n])
    else:
        nv, nf = shapes["meshes"]["nucleus"]
        tri = nv[nf]; ar = 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
        cum = np.cumsum(ar) / ar.sum()
        while len(pores) < geo["n_npc"]:
            t = tri[int(np.searchsorted(cum, rng.random()))]
            a_, b_ = rng.random(2)
            if a_ + b_ > 1:
                a_, b_ = 1 - a_, 1 - b_
            p = t[0] + a_ * (t[1] - t[0]) + b_ * (t[2] - t[0])
            n = _grad(lambda x, y, z: nuc_sdf(x, y, z, Gq), p)
            s = p - 0.5 * PERI_W * n
            if vac_d(s)[0] > nvj_clear and all(np.linalg.norm(s - q[:3]) > 150.0 for q in pores):
                pores.append(np.r_[s, n])
    pores = np.array(pores)
    npc_pos = pores[:, :3]

    # ---- mitochondrial network: cortical tubules with short cristae and crista-free nucleoid sites (rebuilt
    # 2026-11-13 against measurements of glucose-grown cells; constants and sources at MITO_* above)
    Vtarget = V["mitochondrion"] * 1e9
    segs, seg_r, seg_chain = [], [], []       # axis end points, outer radius and chain of every segment
    chains = []                               # per chain: dict(parent segment or -1, segment indices in order)
    branch_pts = []
    total = 0.0; vol = 0.0

    def seg_clear(samp, r, ignore):
        """sample points of a new segment keep their distance from the nucleus, the vacuole and other tubules"""
        if not ((nuc_d(samp) > r + 40).all() and (vac_d(samp) > r + 40).all()):
            return False
        if segs:
            A = np.array([s_[0] for s_ in segs]); B = np.array([s_[1] for s_ in segs]); ab = B - A
            lim = r + np.array(seg_r) + 60.0
            for s_ in samp:
                h = np.clip(((s_ - A) * ab).sum(1) / (ab * ab).sum(1), 0, 1)
                dd = np.linalg.norm(s_ - (A + ab * h[:, None]), axis=1) - lim
                if ignore:
                    dd[list(ignore)] = 1e9
                if (dd < 0).any():
                    return False
        return True

    n_pieces = 1 + int(rng.random() < 0.5)    # 1-3 separate mitochondria per cell (Osman 2015; Seel 2023: median 2)
    pieces = 0; tries = 0
    while vol < Vtarget and tries < 600:
        tries += 1
        new_piece = not segs or (pieces < n_pieces and rng.random() < 0.2)
        if new_piece:
            parent = -1
            u_ = _rand_unit(rng)                                   # 457-597 nm below the outer surface
            start = u_ * ray_depth(Gq, u_, rng.uniform(457, 597), depth)
            want = rng.uniform(3000, 6000)
        else:                                                      # a branch from the side of an existing tubule
            parent = int(rng.integers(len(segs)))
            a, b = segs[parent]
            start = a + (b - a) * rng.uniform(0.2, 0.8)
            if any(np.linalg.norm(start - bp) < 600.0 for bp in branch_pts):
                continue
            want = rng.uniform(1500, 4500)
        d = _rand_unit(rng); d -= start * np.dot(d, start) / np.dot(start, start); d /= np.linalg.norm(d)
        if parent >= 0:                                            # leave the parent sideways
            pa = (segs[parent][1] - segs[parent][0]); pa /= np.linalg.norm(pa)
            d -= pa * np.dot(d, pa) * 0.8; d /= np.linalg.norm(d)
        p = start; mine = []; r = float(np.clip(MITO_R + rng.normal(0, MITO_R_SD), MITO_R_MIN, MITO_R_MAX)); done = 0.0
        while done < want and vol < Vtarget:
            d = d + rng.normal(0, MITO_TURN, 3); d -= p * np.dot(d, p) / np.dot(p, p); d /= np.linalg.norm(d)
            q = p + d * MITO_STEP
            uq = q / np.linalg.norm(q)                             # keep within 437-637 nm of the outer surface
            q = uq * np.clip(np.linalg.norm(q), ray_depth(Gq, uq, 637.0, depth), ray_depth(Gq, uq, 437.0, depth))
            r = float(np.clip(r + rng.normal(0, MITO_R_WALK) + 0.15 * (MITO_R - r), MITO_R_MIN, MITO_R_MAX))
            samp = p + (q - p) * np.linspace(0.34, 1, 3)[:, None]
            ignore = set(mine[-MITO_IGNORE:])                      # its own last segments, and where it leaves its parent
            if parent >= 0 and len(mine) < MITO_IGNORE:
                pc = chains[seg_chain[parent]]["segs"]; k = pc.index(parent)
                ignore |= set(pc[max(0, k - 3):k + 4])
            if not seg_clear(samp, r, ignore):
                break
            segs.append((p.copy(), q.copy())); seg_r.append(r); seg_chain.append(len(chains)); mine.append(len(segs) - 1)
            L_ = np.linalg.norm(q - p); done += L_; total += L_; vol += np.pi * r * r * L_; p = q
        if len(mine) < 4:                                          # a stub: undo
            for _ in mine:
                a, b = segs.pop(); r_ = seg_r.pop(); seg_chain.pop()
                L_ = np.linalg.norm(b - a); total -= L_; vol -= np.pi * r_ * r_ * L_
            continue
        vol += 2 / 3 * np.pi * seg_r[mine[-1]] ** 3 + (2 / 3 * np.pi * seg_r[mine[0]] ** 3 if parent < 0 else 0.0)
        chains.append(dict(parent=parent, segs=mine))
        if parent < 0:
            pieces += 1
        else:
            branch_pts.append(start.copy())
    tips = [segs[c["segs"][-1]][1] for c in chains] + [segs[c["segs"][0]][0] for c in chains if c["parent"] < 0]
    log(f"mito: {pieces} piece(s), {len(chains)} chains, {len(branch_pts)} branch points, {len(tips)} tips; {len(segs)} segments, "
        f"length {total:.0f} nm, outer diameter {2 * np.min(seg_r):.0f}-{2 * np.max(seg_r):.0f} nm (median {2 * np.median(seg_r):.0f}), "
        f"est. volume {vol / 1e9:.3f} um3 (target {Vtarget / 1e9:.3f})")
    for (a, b), r in zip(segs, seg_r):
        cen = (a + b) / 2; rad = np.linalg.norm(b - a) / 2 + r
        add(T_CAPSULE, [*a, *b, r], C_OM, cen, rad)
        add(T_CAPSULE, [*a, *b, r - IMS_W], C_IM, cen, rad)

    # nucleoid sites: about NUCLEOID_PER_UM per um of tubule, 69% of the tips with one within 0.5 um (Osman 2015), at
    # least NUCLEOID_MIN_GAP apart; each keeps a stretch of matrix free of cristae and sits off the axis, towards the
    # inner boundary membrane
    def chain_point(ci, s):
        """axis point, direction, segment index at arclength s (nm) of chain ci"""
        for k in chains[ci]["segs"]:
            a, b = segs[k]; L_ = np.linalg.norm(b - a)
            if s <= L_:
                return a + (b - a) * (s / L_), (b - a) / L_, k
            s -= L_
        a, b = segs[chains[ci]["segs"][-1]]
        return b.copy(), (b - a) / np.linalg.norm(b - a), chains[ci]["segs"][-1]
    clen = [sum(np.linalg.norm(segs[k][1] - segs[k][0]) for k in c["segs"]) for c in chains]
    n_nucd = max(1, int(round(NUCLEOID_PER_UM * total / 1000.0)))
    nucd = []                                 # (axis point, axis direction, centre, segment)

    def add_nucleoid(ci, s):
        if not (NUCLEOID_CLEAR * 0.5 < s < clen[ci] - 60.0):
            return False
        pt, ax, k = chain_point(ci, s)
        if any(np.linalg.norm(pt - n_[0]) < NUCLEOID_MIN_GAP for n_ in nucd) or any(np.linalg.norm(pt - bp) < 250.0 for bp in branch_pts):
            return False
        e = _rand_unit(rng); e -= ax * np.dot(e, ax); e /= np.linalg.norm(e)
        nucd.append((pt, ax, pt + e * rng.uniform(0.25, 0.45) * (seg_r[k] - IMS_W), k))
        return True
    for ci, c in enumerate(chains):           # near the tips first
        ends = [clen[ci] - rng.uniform(150, 400)] + ([rng.uniform(150, 400)] if c["parent"] < 0 else [])
        for s in ends:
            if len(nucd) < n_nucd and rng.random() < 0.69:
                add_nucleoid(ci, s)
    for _ in range(4000):
        if len(nucd) >= n_nucd:
            break
        ci = int(rng.choice(len(chains), p=np.array(clen) / sum(clen)))
        add_nucleoid(ci, rng.uniform(0, clen[ci]))

    # cristae: short lamellae perpendicular to the tubule, opening into the intermembrane space where their far rim
    # crosses the inner boundary membrane; on two opposite sides that turn slowly along the tubule (Stoldt 2019)
    n_cristae = 0; cr_area = 0.0
    for ci, c in enumerate(chains):
        a, b = segs[c["segs"][0]]; ax0 = (b - a) / np.linalg.norm(b - a)
        e = _rand_unit(rng); e -= ax0 * np.dot(e, ax0); e /= np.linalg.norm(e)
        s = 70.0 if c["parent"] < 0 else seg_r[c["parent"]] + 60.0; kk = 0; phi = 0.0
        while s < clen[ci] - 70.0:
            pt, ax, k = chain_point(ci, s)
            e = e - ax * np.dot(e, ax); e /= np.linalg.norm(e)                    # carry the side along the tubule
            s_here = s
            s += CRISTA_GAP_MIN + rng.exponential(1000.0 / CRISTA_PER_UM - CRISTA_GAP_MIN)
            if any(np.linalg.norm(pt - n_[0]) < NUCLEOID_CLEAR for n_ in nucd) or any(np.linalg.norm(pt - bp) < 170.0 for bp in branch_pts):
                continue
            r_im = seg_r[k] - IMS_W
            reach = float(np.clip(rng.gamma(CRISTA_REACH_SHAPE, CRISTA_REACH / CRISTA_REACH_SHAPE), CRISTA_REACH_MIN, min(CRISTA_REACH_MAX, 2 * r_im - 40.0)))
            phi += rng.normal(0, 0.35)
            side = (np.cos(phi + np.pi * kk) * e + np.sin(phi + np.pi * kk) * np.cross(ax, e)); kk += 1
            half = reach / 2.0; Rd = half - CRISTA_HALF
            cen = pt + side * ((r_im + 2.0) - half)
            add(T_PILL, [*cen, *ax, Rd, CRISTA_HALF], C_CRISTA, cen, half + 1)
            n_cristae += 1; cr_area += 2 * np.pi * Rd * Rd + 2 * np.pi ** 2 * Rd * CRISTA_HALF + 4 * np.pi * CRISTA_HALF ** 2
    om_area = sum(2 * np.pi * r * np.linalg.norm(b - a) for (a, b), r in zip(segs, seg_r))
    for (a, b), r in zip(segs, seg_r):
        for t in np.linspace(0, 1, max(2, int(np.linalg.norm(b - a) / 100) + 1)):
            occupied.append((a + (b - a) * t, r))
    log(f"cristae: {n_cristae} ({n_cristae / (total / 1000):.1f} per um of tubule), membrane about {cr_area / 1e6:.1f} um2 = "
        f"{cr_area / om_area:.2f} x the outer membrane ({om_area / 1e6:.1f} um2); {len(nucd)} nucleoid sites")

    def place_sphere(r, cls, gap=40.0, rmin=0, rmax=None):
        for _ in range(20000):
            c = _rand_unit(rng) * rng.uniform(rmin, rmax or Rcell)
            if free(c, r, gap):
                add(T_SPHERE, [*c, r], cls, c, r); occupied.append((c, r)); return c
        raise RuntimeError("could not place sphere")

    # ---- lipid droplets: sizes from the segmented cells (depositor labels). Real cell: the segmented droplets as
    # volume-equivalent spheres at their positions (nudged where the reshaped vacuoles or nucleus now are)
    if shapes is None:
        ldr = [(3 * v / (4 * np.pi)) ** (1 / 3) for v in
               _sample_volumes(rng, stats["lipid_droplet"] * geo["cell_volume_um3"] * 1e9, V["lipid_droplet"] * 1e9)]
        for r in ldr:
            place_sphere(r, C_LD)
    else:
        ldr, moved = [], 0
        for c, r in shapes["droplets"]:
            c = np.asarray(c, np.float64)
            for t in range(4000):
                cc = c if t == 0 else c + _rand_unit(rng) * rng.uniform(0, 20 + 0.3 * t)
                if free(cc, r, 40.0):
                    break
            else:
                cc = None
            if cc is None:
                place_sphere(r, C_LD)
            else:
                add(T_SPHERE, [*cc, r], C_LD, cc, r); occupied.append((cc, r))
            moved += t > 0
            ldr.append(r)
        log(f"lipid droplets: {len(ldr)} segmented droplets, {moved} moved to free space")
    rld = np.median(ldr); nld = len(ldr)
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
            c = _rand_unit(rng) * rng.uniform(300, Rcell)
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
    log(f"LD n={nld} median r={rld:.0f}  perox r={rpx:.0f}  endosome r={rend:.0f}  golgi cisternae {n_cis}")

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
        rr = ray_depth(Gq, u, WALL + CER_GAP + 0.5 * CER_LUMEN, depth)    # the sheet's radius along the patch axis
        add(T_CERPATCH, [*u, np.cos(al), al, rc, 0.5 * CER_LUMEN], C_CER, u * rr, al * rr + 20 + abs(rr - rc) + (150 if fields else 0))
    log(f"cER: {len(caps)} patches, PM coverage {covered.mean():.3f}")

    # ---- eisosome furrows: where there is no cortical ER, the PM normal is close to radial (the slot opens along the
    # radial direction) and the nucleus, vacuole and other furrows are clear; random direction in the membrane plane
    rng2 = np.random.default_rng(seed + 1)
    eis, tries = [], 0
    cap_u = np.array([u for u, _ in caps]); cap_a = np.array([al for _, al in caps])
    while len(eis) < EIS_N and tries < 20000:
        tries += 1
        u = _rand_unit(rng2)
        L = float(np.clip(rng2.normal(*EIS_LEN), 150.0, 450.0))
        rp = ray_depth(Gq, u, WALL, depth)                               # PM midsurface along u
        if np.degrees(np.arccos(np.clip(np.dot(-_grad(lambda x, y, z: outer_depth(x, y, z, Gq), u * rp), u), -1, 1))) > 8.0:
            continue
        if len(cap_u) and (np.arccos(np.clip(cap_u @ u, -1, 1)) < cap_a + (0.5 * L + 80.0) / rp).any():
            continue
        if any(np.linalg.norm(u * rp - c) < 0.5 * (L + l) + 300.0 for c, l in eis):
            continue
        t = np.cross(u, _rand_unit(rng2)); t /= np.linalg.norm(t)
        bot = u * (rp - (EIS_DEPTH - EIS_R))                             # bottom axis midpoint
        a_, b_ = bot - 0.5 * L * t, bot + 0.5 * L * t
        probe = bot - u * 60.0 + np.linspace(-0.5, 0.5, 7)[:, None] * L * t
        if (nuc_d(probe) < 60.0).any() or (vac_d(probe) < 60.0).any():
            continue
        h = EIS_DEPTH - EIS_R + 15.0                                     # up to 15 nm above the PM midsurface
        add(T_SLOT, [*a_, *b_, EIS_R, h], C_PMINV, bot + 0.5 * h * u, np.hypot(0.5 * L, 0.5 * h) + EIS_R + EIS_K)
        eis.append((u * rp, L))
    log(f"eisosome furrows: {len(eis)}, length {np.mean([l for _, l in eis]):.0f} nm (mean), {EIS_DEPTH:.0f} nm deep")
    # ---- endocytic invaginations: same site rules, with room for the actin patch around them
    endo, tries = [], 0
    while len(endo) < n_endo and tries < 20000:
        tries += 1
        u = _rand_unit(rng2)
        rp = ray_depth(Gq, u, WALL, depth)
        if np.degrees(np.arccos(np.clip(np.dot(-_grad(lambda x, y, z: outer_depth(x, y, z, Gq), u * rp), u), -1, 1))) > 8.0:
            continue
        if len(cap_u) and (np.arccos(np.clip(cap_u @ u, -1, 1)) < cap_a + 100.0 / rp).any():
            continue
        if any(np.linalg.norm(u * rp - c) < 0.5 * l + 150.0 for c, l in eis) or any(np.linalg.norm(u * rp - c) < 300.0 for c in endo):
            continue
        probe = u * (rp - np.array([100.0, 200.0, 300.0]))[:, None]
        if (nuc_d(probe) < 40.0).any() or (vac_d(probe) < 40.0).any():
            continue
        D = rng2.uniform(*ENDO_DEPTH)
        bot = u * (rp - (D - ENDO_R)); t = np.cross(u, _rand_unit(rng2)); t /= np.linalg.norm(t)
        h = max(D - ENDO_R + 15.0, 0.0)                                 # shallow sites: a dimple of the tip sphere
        add(T_SLOT, [*bot, *(bot + ENDO_LEN * t), ENDO_R, h], C_PMINV, bot + 0.5 * h * u, 0.5 * abs(h) + ENDO_R + EIS_K)
        endo.append(u * rp)
    if n_endo:
        log(f"endocytic invaginations: {len(endo)} of {n_endo}, tubes {2 * ENDO_R:.0f} nm wide, {ENDO_DEPTH[0]:.0f}-{ENDO_DEPTH[1]:.0f} nm deep")
    # ---- bud scars: ring on the outer surface where the normal is close to radial and the rim stays inside the
    # enclosing sphere; the torus centre is lowered by the surface's sag under the ring (curved cell)
    scars, tries = [], 0
    th = np.linspace(0, 2 * np.pi, 48, endpoint=False)
    while len(scars) < n_scars and tries < 5000:
        tries += 1
        u = _rand_unit(rng2)
        rs = ray_depth(Gq, u, 0.0, depth)
        nrm = -_grad(lambda x, y, z: outer_depth(x, y, z, Gq), u * rs)
        if np.degrees(np.arccos(np.clip(nrm @ u, -1, 1))) > 10.0:
            continue
        c = u * rs; e0 = _rand_unit(rng2)
        for _ in range(6):           # ring plane fitted to the surface under the ring: tube centres on the surface
            e1 = np.cross(nrm, e0); e1 /= np.linalg.norm(e1); e2 = np.cross(nrm, e1)
            ring = c + SCAR_R * (np.cos(th)[:, None] * e1 + np.sin(th)[:, None] * e2)
            onsurf = ring + nrm * depth(ring)[:, None]
            c = onsurf.mean(0); nn = np.linalg.svd(onsurf - c)[2][2]; nrm = nn if nn @ u > 0 else -nn
        ring = c + SCAR_R * (np.cos(th)[:, None] * e1 + np.sin(th)[:, None] * e2)
        if np.abs(depth(ring)).max() > 10.0:          # surface too far from planar under the ring
            continue
        if np.linalg.norm(ring, axis=1).max() + SCAR_r + 10.0 > Gq[0] or any(np.linalg.norm(c - s) < 2 * SCAR_R + 300.0 for s in scars):
            continue
        add(T_TORUS, [*c, *nrm, SCAR_R, SCAR_r], C_SCAR, c, SCAR_R + SCAR_r + SCAR_K)
        scars.append(c)
    if n_scars:
        log(f"bud scars: {len(scars)}, rim {2 * (SCAR_R + SCAR_r) / 1000:.1f} um across, {SCAR_r:.0f} nm high")

    d = dict(G=G, pores=pores.ravel().tolist(), PT=PT, PQ=np.array(PQ).ravel().tolist(), PC=PC,
             bc=np.array(bc).ravel().tolist(), br=br, npc_pos=npc_pos.ravel().tolist(),
             mito_segments=[[list(a), list(b)] for a, b in segs], mito_radii=[float(r) for r in seg_r],
             mito_nucleoids=[[list(map(float, n_[2])), list(map(float, n_[1]))] for n_ in nucd], seed=seed,
             consts=dict(HALF=HALF, PORE_R=PORE_R, PERI_W=PERI_W, IMS_W=IMS_W, MITO_R=MITO_R,
                         CER_GAP=CER_GAP, CER_LUMEN=CER_LUMEN, EIS_K=EIS_K, SCAR_K=SCAR_K),
             label_names=LABEL_NAMES, mem_names=MEM_NAMES)
    return Cell(d, fields)
