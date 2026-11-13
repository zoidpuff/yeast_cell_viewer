"""Overlap relaxation after random sequential addition (literature item L1: densities beyond the RSA jamming limit).

Random sequential addition stops at ~35-40% sphere-tree volume fraction, short of the crowding in cells. After a tile
has been packed by RSA, the copies that did not fit are inserted anyway at random positions in their compartment, and
all soluble molecules of the tile are relaxed as rigid bodies: overlapping spheres of different molecules push each
other apart (translation and a small rotation per molecule, a position-based overlap removal), spheres that cross a
membrane or leave their compartment are pushed back along the membrane's distance gradient, and molecule centres stay
in the tile. Membrane proteins, pre-placed structures and finished neighbour tiles are fixed obstacles. Afterwards the
configuration is made exact again: molecules are accepted one by one (previously placed ones first) only if they
overlap nothing accepted so far and lie in their compartment, as in RSA.
"""
import numpy as np
import numba as nb
import cell3d as c3
import pack3d as p3

H = 8.0                 # neighbour-grid cell (nm); larger than any two touching sphere-tree radii in practice
N_ITER = 60
ROUND_FRAC = 0.15      # extra copies inserted per round, relative to the movable molecules
N_ROUNDS = 4
STEP_MAX = 0.6          # nm per iteration
ROT_MAX = 0.06          # rad per iteration
WALL_EVERY = 4


@nb.njit(cache=True)
def _grid(lo, dims, h, P, R):
    head = -np.ones(dims[0] * dims[1] * dims[2], np.int64); nxt = np.empty(P.shape[0], np.int64)
    for s in range(P.shape[0]):
        i = min(max(int((P[s, 0] - lo[0]) / h), 0), dims[0] - 1)
        j = min(max(int((P[s, 1] - lo[1]) / h), 0), dims[1] - 1)
        k = min(max(int((P[s, 2] - lo[2]) / h), 0), dims[2] - 1)
        c = (i * dims[1] + j) * dims[2] + k
        nxt[s] = head[c]; head[c] = s
    return head, nxt


@nb.njit(cache=True, fastmath=True)
def _rotq(q, wx, wy, wz):
    """q <- exp(w/2) * q for a small rotation vector w"""
    a = np.sqrt(wx * wx + wy * wy + wz * wz)
    if a < 1e-12:
        return
    s = np.sin(0.5 * a) / a; c = np.cos(0.5 * a)
    b0, b1, b2, b3 = c, wx * s, wy * s, wz * s
    a0, a1, a2, a3 = q[0], q[1], q[2], q[3]
    q[0] = b0 * a0 - b1 * a1 - b2 * a2 - b3 * a3
    q[1] = b0 * a1 + b1 * a0 + b2 * a3 - b3 * a2
    q[2] = b0 * a2 - b1 * a3 + b2 * a0 + b3 * a1
    q[3] = b0 * a3 + b1 * a2 - b2 * a1 + b3 * a0
    n = np.sqrt(q[0] ** 2 + q[1] ** 2 + q[2] ** 2 + q[3] ** 2)
    for t in range(4):
        q[t] /= n


@nb.njit(cache=True, fastmath=True)
def _wall(x, y, z, r, mask, G, pores, pidx, PT, PQ, PC, out, g):
    """push (g) for a sphere that crosses a membrane or sits outside its allowed compartments; returns violation"""
    m = c3.mem_sdfs(x, y, z, G, pores, pidx, PT, PQ, PC, out)
    l = c3.label_from(x, y, z, G, pores, out, m)
    dmin = 1e9; mi = -1
    for k in range(c3.NMEM):
        if abs(out[k]) < dmin:
            dmin = abs(out[k]); mi = k
    ok_label = (mask >> l) & 1
    viol = r + c3.HALF - dmin
    g[0] = 0.0; g[1] = 0.0; g[2] = 0.0
    if ok_label and viol <= 0:
        return 0.0
    e = 0.5
    for a in range(3):
        xp = x + (e if a == 0 else 0.0); yp = y + (e if a == 1 else 0.0); zp = z + (e if a == 2 else 0.0)
        c3.mem_sdfs(xp, yp, zp, G, pores, pidx, PT, PQ, PC, out); dp = out[mi]
        xm = x - (e if a == 0 else 0.0); ym = y - (e if a == 1 else 0.0); zm = z - (e if a == 2 else 0.0)
        c3.mem_sdfs(xm, ym, zm, G, pores, pidx, PT, PQ, PC, out); dm = out[mi]
        g[a] = (dp - dm) / (2 * e)
    gn = np.sqrt(g[0] ** 2 + g[1] ** 2 + g[2] ** 2) + 1e-12
    c3.mem_sdfs(x, y, z, G, pores, pidx, PT, PQ, PC, out); d = out[mi]
    # on the allowed side: move away from the membrane; on the wrong side: move across it
    sgn = 1.0 if d > 0 else -1.0
    if not ok_label:
        sgn = -sgn
        viol = abs(d) + r + c3.HALF
    for a in range(3):
        g[a] = g[a] / gn * sgn * viol
    return viol


@nb.njit(cache=True, fastmath=True)
def relax(seed, lo, size, halo, G, pores, pidx, PT, PQ, PC, st_off, STC, STR, regmask,
          typ, pos, quat, reg, movable, OC, OR, n_iter):
    """typ/pos/quat/reg: all tile instances (movable ones are relaxed in place); OC/OR: fixed spheres"""
    np.random.seed(seed)
    n = typ.shape[0]
    # sphere bookkeeping: movable instances' spheres, then fixed tile instances, then obstacles
    nsm = 0; nsf = 0
    for i in range(n):
        k = st_off[typ[i] + 1] - st_off[typ[i]]
        if movable[i]:
            nsm += k
        else:
            nsf += k
    NS = nsm + nsf + OR.shape[0]
    P = np.empty((NS, 3)); R = np.empty(NS); own = np.empty(NS, np.int64)
    s0 = np.empty(n, np.int64); s1 = np.empty(n, np.int64)
    s = 0
    for pas in range(2):
        for i in range(n):
            if (pas == 0) != movable[i]:
                continue
            a0 = st_off[typ[i]]; a1 = st_off[typ[i] + 1]
            s0[i] = s; s1[i] = s + a1 - a0
            for k in range(a1 - a0):
                R[s + k] = STR[a0 + k]; own[s + k] = i
            s += a1 - a0
    for k in range(OR.shape[0]):
        P[s + k, 0] = OC[k, 0]; P[s + k, 1] = OC[k, 1]; P[s + k, 2] = OC[k, 2]; R[s + k] = OR[k]; own[s + k] = -1
    M = np.empty((3, 3)); qq = np.empty(4)

    def place_spheres(i):
        a0 = st_off[typ[i]]
        for t in range(4):
            qq[t] = quat[i, t]
        p3.quat_mat(qq, M)
        for k in range(s1[i] - s0[i]):
            c = STC[a0 + k]
            for a in range(3):
                P[s0[i] + k, a] = M[a, 0] * c[0] + M[a, 1] * c[1] + M[a, 2] * c[2] + pos[i, a]

    for i in range(n):
        place_spheres(i)
    glo = lo - halo
    dims = np.empty(3, np.int64)
    for a in range(3):
        dims[a] = int(np.ceil((size + 2 * halo) / H)) + 1
    F = np.zeros((n, 3)); T = np.zeros((n, 3)); I = np.ones(n)
    for i in range(n):
        if movable[i]:
            acc = 1.0
            for k in range(s0[i], s1[i]):
                acc += (P[k, 0] - pos[i, 0]) ** 2 + (P[k, 1] - pos[i, 1]) ** 2 + (P[k, 2] - pos[i, 2]) ** 2
            I[i] = acc / max(s1[i] - s0[i], 1)
    out = np.zeros(c3.NMEM); g = np.zeros(3)
    near_mem = np.zeros(n, np.bool_)
    total = 0.0
    for it in range(n_iter):
        if it % WALL_EVERY == 0:          # wall forces only for molecules whose bounding sphere may reach a membrane
            for i in range(n):
                if movable[i]:
                    c3.mem_sdfs(pos[i, 0], pos[i, 1], pos[i, 2], G, pores, pidx, PT, PQ, PC, out)
                    br = 0.0
                    for k in range(s0[i], s1[i]):
                        br = max(br, np.sqrt((P[k, 0] - pos[i, 0]) ** 2 + (P[k, 1] - pos[i, 1]) ** 2 + (P[k, 2] - pos[i, 2]) ** 2) + R[k])
                    near_mem[i] = c3.min_mem_dist(out, -1) < br + c3.HALF + 2.0
        head, nxt = _grid(glo, dims, H, P, R)
        F[:, :] = 0.0; T[:, :] = 0.0; total = 0.0
        for i in range(n):
            if not movable[i]:
                continue
            for sp in range(s0[i], s1[i]):
                x, y, z, r = P[sp, 0], P[sp, 1], P[sp, 2], R[sp]
                ci = min(max(int((x - glo[0]) / H), 0), dims[0] - 1)
                cj = min(max(int((y - glo[1]) / H), 0), dims[1] - 1)
                ck = min(max(int((z - glo[2]) / H), 0), dims[2] - 1)
                for di in range(max(ci - 1, 0), min(ci + 2, dims[0])):
                    for dj in range(max(cj - 1, 0), min(cj + 2, dims[1])):
                        for dk in range(max(ck - 1, 0), min(ck + 2, dims[2])):
                            e = head[(di * dims[1] + dj) * dims[2] + dk]
                            while e >= 0:
                                if own[e] != i:
                                    dx = x - P[e, 0]; dy = y - P[e, 1]; dz = z - P[e, 2]
                                    d2 = dx * dx + dy * dy + dz * dz; rr = r + R[e]
                                    if d2 < rr * rr:
                                        d = np.sqrt(d2) + 1e-9; ov = rr - d
                                        w = 0.5 if (own[e] >= 0 and movable[own[e]]) else 1.0
                                        fx = w * ov * dx / d; fy = w * ov * dy / d; fz = w * ov * dz / d
                                        F[i, 0] += fx; F[i, 1] += fy; F[i, 2] += fz
                                        rx = x - pos[i, 0]; ry = y - pos[i, 1]; rz = z - pos[i, 2]
                                        T[i, 0] += ry * fz - rz * fy; T[i, 1] += rz * fx - rx * fz; T[i, 2] += rx * fy - ry * fx
                                        total += ov
                                e = nxt[e]
                if it % WALL_EVERY == 0 and near_mem[i]:
                    v = _wall(x, y, z, r, regmask[reg[i]], G, pores, pidx, PT, PQ, PC, out, g)
                    if v > 0:
                        F[i, 0] += g[0]; F[i, 1] += g[1]; F[i, 2] += g[2]; total += v
        if total < 1e-3 * max(nsm, 1):
            break
        for i in range(n):
            if not movable[i]:
                continue
            fn = np.sqrt(F[i, 0] ** 2 + F[i, 1] ** 2 + F[i, 2] ** 2)
            sc = 1.0 / max(1.0, s1[i] - s0[i]) ** 0.5
            stp = min(STEP_MAX, fn * sc) / (fn + 1e-12)
            for a in range(3):
                pos[i, a] += F[i, a] * stp
                pos[i, a] = min(max(pos[i, a], lo[a] + 1e-3), lo[a] + size - 1e-3)
            tn = np.sqrt(T[i, 0] ** 2 + T[i, 1] ** 2 + T[i, 2] ** 2)
            ang = min(ROT_MAX, 0.5 * tn / (I[i] * max(1.0, s1[i] - s0[i])))
            if tn > 1e-9:
                for t in range(4):
                    qq[t] = quat[i, t]
                _rotq(qq, T[i, 0] / tn * ang, T[i, 1] / tn * ang, T[i, 2] / tn * ang)
                for t in range(4):
                    quat[i, t] = qq[t]
            place_spheres(i)
    return total


@nb.njit(cache=True, fastmath=True)
def accept(lo, size, halo, G, pores, pidx, PT, PQ, PC, st_off, STC, STR, BR, regmask,
           typ, pos, quat, reg, movable, OC, OR, cap_s, cap_e, band_m, band_w):
    """keep[i]: fixed instances always; movable ones in order if they overlap nothing kept and lie in their compartment
    (and, for membrane-face regions, within band_w of the membrane band_m: the relaxation must not push them out)"""
    glo = lo - halo; ghi = lo + size + halo
    dims, head, nxt, sid, cen, rad, cnt = p3.grid_new(glo, ghi, p3.H_GRID, cap_s, cap_e)
    p3.g_insert(glo, p3.H_GRID, dims, head, nxt, sid, cen, rad, cnt, OC, OR)
    n = typ.shape[0]; keep = np.zeros(n, np.bool_)
    M = np.empty((3, 3)); q = np.empty(4); out = np.zeros(c3.NMEM)
    for pas in range(2):
        for i in range(n):
            if (pas == 0) == bool(movable[i]):
                continue
            a0 = st_off[typ[i]]; a1 = st_off[typ[i] + 1]
            for t in range(4):
                q[t] = quat[i, t]
            p3.quat_mat(q, M); Pi = np.empty((a1 - a0, 3)); p3.transform(STC[a0:a1], M, pos[i], Pi)
            if movable[i]:
                if p3.g_overlap(glo, p3.H_GRID, dims, head, nxt, sid, cen, rad, Pi, STR[a0:a1], p3.TOL):
                    continue
                if not p3.check_volume(Pi, STR[a0:a1], pos[i, 0], pos[i, 1], pos[i, 2], BR[typ[i]], regmask[reg[i]],
                                       G, pores, pidx, PT, PQ, PC, out):
                    continue
                if band_m[reg[i]] >= 0:
                    c3.mem_sdfs(pos[i, 0], pos[i, 1], pos[i, 2], G, pores, pidx, PT, PQ, PC, out)
                    if abs(out[band_m[reg[i]]]) > band_w[reg[i]]:
                        continue
            if not p3.g_insert(glo, p3.H_GRID, dims, head, nxt, sid, cen, rad, cnt, Pi, STR[a0:a1]):
                break
            keep[i] = True
    return keep
