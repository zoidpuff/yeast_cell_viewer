"""3D random-sequential-addition packing of the whole cell in 200 nm tiles (numba core).

Collision geometry: per-ingredient sphere trees (st_c, st_r) on a local uniform grid (tile + halo).
Region tests: exact analytic SDFs from cell3d (label + distance to every membrane).
Membrane proteins: sampled on the tile's marching-cubes mesh by area, projected onto the exact surface,
library +z aligned to the SDF gradient (cytoplasmic / N-side), random spin."""
import numpy as np
import numba as nb
import cell3d as c3

H_GRID = 6.0
TOL = 0.3


# ------------------------------------------------------------------ grid
@nb.njit(cache=True)
def grid_new(lo, hi, h, cap_s, cap_e):
    dims = np.empty(3, np.int64)
    for a in range(3):
        dims[a] = int(np.ceil((hi[a] - lo[a]) / h)) + 1
    head = -np.ones(dims[0] * dims[1] * dims[2], np.int64)
    nxt = np.empty(cap_e, np.int64); sid = np.empty(cap_e, np.int64)
    cen = np.empty((cap_s, 3), np.float64); rad = np.empty(cap_s, np.float64)
    cnt = np.zeros(2, np.int64)   # n_sph, n_ent
    return dims, head, nxt, sid, cen, rad, cnt


@nb.njit(cache=True, fastmath=True)
def _rng3(x, r, lo, h, d):
    i0 = int(np.floor((x - r - lo) / h)); i1 = int(np.floor((x + r - lo) / h))
    return max(i0, 0), min(i1, d - 1)


@nb.njit(cache=True, fastmath=True)
def g_overlap(lo, h, dims, head, nxt, sid, cen, rad, P, R, tol):
    for s in range(P.shape[0]):
        x, y, z, r = P[s, 0], P[s, 1], P[s, 2], R[s]
        i0, i1 = _rng3(x, r, lo[0], h, dims[0]); j0, j1 = _rng3(y, r, lo[1], h, dims[1]); k0, k1 = _rng3(z, r, lo[2], h, dims[2])
        for i in range(i0, i1 + 1):
            for j in range(j0, j1 + 1):
                for k in range(k0, k1 + 1):
                    e = head[(i * dims[1] + j) * dims[2] + k]
                    while e >= 0:
                        q = sid[e]
                        dx = x - cen[q, 0]; dy = y - cen[q, 1]; dz = z - cen[q, 2]
                        rr = r + rad[q] - tol
                        if rr > 0 and dx * dx + dy * dy + dz * dz < rr * rr:
                            return True
                        e = nxt[e]
    return False


@nb.njit(cache=True, fastmath=True)
def g_insert(lo, h, dims, head, nxt, sid, cen, rad, cnt, P, R):
    for s in range(P.shape[0]):
        x, y, z, r = P[s, 0], P[s, 1], P[s, 2], R[s]
        i0, i1 = _rng3(x, r, lo[0], h, dims[0]); j0, j1 = _rng3(y, r, lo[1], h, dims[1]); k0, k1 = _rng3(z, r, lo[2], h, dims[2])
        if i1 < i0 or j1 < j0 or k1 < k0:
            continue
        need = (i1 - i0 + 1) * (j1 - j0 + 1) * (k1 - k0 + 1)
        if cnt[0] + 1 > cen.shape[0] or cnt[1] + need > nxt.shape[0]:
            return False
        q = cnt[0]; cnt[0] += 1
        cen[q, 0] = x; cen[q, 1] = y; cen[q, 2] = z; rad[q] = r
        for i in range(i0, i1 + 1):
            for j in range(j0, j1 + 1):
                for k in range(k0, k1 + 1):
                    c = (i * dims[1] + j) * dims[2] + k
                    e = cnt[1]; cnt[1] += 1
                    sid[e] = q; nxt[e] = head[c]; head[c] = e
    return True


# ------------------------------------------------------------------ rotations
@nb.njit(cache=True, fastmath=True)
def quat_mat(q, M):
    a, b, c, d = q[0], q[1], q[2], q[3]
    M[0, 0] = a*a+b*b-c*c-d*d; M[0, 1] = 2*(b*c-a*d); M[0, 2] = 2*(b*d+a*c)
    M[1, 0] = 2*(b*c+a*d); M[1, 1] = a*a-b*b+c*c-d*d; M[1, 2] = 2*(c*d-a*b)
    M[2, 0] = 2*(b*d-a*c); M[2, 1] = 2*(c*d+a*b); M[2, 2] = a*a-b*b-c*c+d*d


@nb.njit(cache=True)
def rand_quat(q):
    s = 0.0
    for i in range(4):
        q[i] = np.random.normal(); s += q[i] * q[i]
    s = np.sqrt(s)
    for i in range(4):
        q[i] /= s
    if q[0] < 0:
        for i in range(4):
            q[i] = -q[i]


@nb.njit(cache=True)
def quat_z_to(n, spin, q):
    """quaternion = (rotation taking +z to unit n) * (rotation by spin about z)."""
    # q1: z -> n
    cx, cy, cz = -n[1], n[0], 0.0          # z x n
    w = 1.0 + n[2]
    if w < 1e-8:
        a1, b1, c1, d1 = 0.0, 1.0, 0.0, 0.0
    else:
        s = np.sqrt(w * w + cx * cx + cy * cy)
        a1, b1, c1, d1 = w / s, cx / s, cy / s, cz / s
    a2, b2, c2, d2 = np.cos(spin / 2), 0.0, 0.0, np.sin(spin / 2)
    q[0] = a1*a2 - b1*b2 - c1*c2 - d1*d2
    q[1] = a1*b2 + b1*a2 + c1*d2 - d1*c2
    q[2] = a1*c2 - b1*d2 + c1*a2 + d1*b2
    q[3] = a1*d2 + b1*c2 - c1*b2 + d1*a2
    if q[0] < 0:
        for i in range(4):
            q[i] = -q[i]


@nb.njit(cache=True, fastmath=True)
def transform(stc, M, t, P):
    for s in range(stc.shape[0]):
        for a in range(3):
            P[s, a] = M[a, 0] * stc[s, 0] + M[a, 1] * stc[s, 1] + M[a, 2] * stc[s, 2] + t[a]


# ------------------------------------------------------------------ region tests
@nb.njit(cache=True, fastmath=True)
def check_volume(P, R, cx, cy, cz, br, mask, G, pores, pidx, PT, PQ, PC, out):
    m = c3.mem_sdfs(cx, cy, cz, G, pores, pidx, PT, PQ, PC, out)
    l = c3.label_from(cx, cy, cz, G, pores, out, m)
    if not (mask >> l) & 1:
        return False
    if c3.min_mem_dist(out, -1) >= br + c3.HALF:
        return True
    for s in range(P.shape[0]):
        m = c3.mem_sdfs(P[s, 0], P[s, 1], P[s, 2], G, pores, pidx, PT, PQ, PC, out)
        l = c3.label_from(P[s, 0], P[s, 1], P[s, 2], G, pores, out, m)
        if not (mask >> l) & 1:
            return False
        if c3.min_mem_dist(out, -1) < R[s] + c3.HALF:
            return False
    return True


@nb.njit(cache=True, fastmath=True)
def check_membrane(P, R, ZL, own, G, pores, pidx, PT, PQ, PC, out):
    for s in range(P.shape[0]):
        c3.mem_sdfs(P[s, 0], P[s, 1], P[s, 2], G, pores, pidx, PT, PQ, PC, out)
        if c3.min_mem_dist(out, own) < R[s] + c3.HALF:
            return False
        if out[c3.M_OUTER] < R[s]:          # outside the cell
            return False
        z = ZL[s]
        if abs(z) >= R[s] + c3.HALF:
            d = out[own]
            if d * z <= 0 or abs(d) < R[s] + c3.HALF - 1.0:
                return False
    return True


@nb.njit(cache=True)
def quat_mul(a, b, q):
    """q = a * b (apply b, then a)"""
    q[0] = a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3]
    q[1] = a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2]
    q[2] = a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1]
    q[3] = a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0]
    if q[0] < 0:
        for i in range(4):
            q[i] = -q[i]


@nb.njit(cache=True, fastmath=True)
def mem_value(px, py, pz, m, G, pores, pidx, PT, PQ, PC, out):
    c3.mem_sdfs(px, py, pz, G, pores, pidx, PT, PQ, PC, out)
    return out[m]


@nb.njit(cache=True, fastmath=True)
def project_to_membrane(p, m, n, G, pores, pidx, PT, PQ, PC, out):
    """Newton projection onto the zero set of membrane m; n receives the unit gradient (cytoplasmic side)."""
    e = 0.4
    for it in range(4):
        d = mem_value(p[0], p[1], p[2], m, G, pores, pidx, PT, PQ, PC, out)
        gx = (mem_value(p[0] + e, p[1], p[2], m, G, pores, pidx, PT, PQ, PC, out) - mem_value(p[0] - e, p[1], p[2], m, G, pores, pidx, PT, PQ, PC, out)) / (2 * e)
        gy = (mem_value(p[0], p[1] + e, p[2], m, G, pores, pidx, PT, PQ, PC, out) - mem_value(p[0], p[1] - e, p[2], m, G, pores, pidx, PT, PQ, PC, out)) / (2 * e)
        gz = (mem_value(p[0], p[1], p[2] + e, m, G, pores, pidx, PT, PQ, PC, out) - mem_value(p[0], p[1], p[2] - e, m, G, pores, pidx, PT, PQ, PC, out)) / (2 * e)
        gn = np.sqrt(gx * gx + gy * gy + gz * gz)
        if gn < 1e-6:
            return False
        n[0] = gx / gn; n[1] = gy / gn; n[2] = gz / gn
        p[0] -= d * n[0]; p[1] -= d * n[1]; p[2] -= d * n[2]
        if abs(d) < 0.02:
            break
    return abs(mem_value(p[0], p[1], p[2], m, G, pores, pidx, PT, PQ, PC, out)) < 0.1


@nb.njit(cache=True, fastmath=True)
def near_pore(p, G, pores, lim):
    return c3.in_pore(p[0], p[1], p[2], pores, lim) >= 0


# ------------------------------------------------------------------ tile packer
KRES = 32     # anchors kept per library type (uniform reservoir sample of everything placed in the tile + halo)


@nb.njit(cache=True)
def _res_add(res_pos, res_n, ty, p):
    n = res_n[ty]
    slot = n if n < KRES else int(np.random.random() * (n + 1))
    if slot < KRES:
        for k in range(3):
            res_pos[ty, slot, k] = p[k]
    res_n[ty] = n + 1


@nb.njit(cache=True)
def _near(res_pos, res_n, pt, d, t):
    """t = a random kept anchor of type pt + d along a random direction"""
    slot = int(np.random.random() * min(res_n[pt], KRES))
    dx = np.random.normal(); dy = np.random.normal(); dz = np.random.normal()
    dn = np.sqrt(dx * dx + dy * dy + dz * dz) + 1e-12
    t[0] = res_pos[pt, slot, 0] + d * dx / dn; t[1] = res_pos[pt, slot, 1] + d * dy / dn; t[2] = res_pos[pt, slot, 2] + d * dz / dn


@nb.njit(cache=True)
def pack_tile(seed, lo, size, halo,
              G, pores, pidx, PT, PQ, PC,
              st_off, STC, STR, BR, regmask,
              jl, jr, jn, n_mem_regions0,
              vox_off, VOXP, vox_h,
              tri_off, TRI, CUMA,
              OC, OR,
              n_mrna, mrna_mean, cyto_region,
              max_tries, give_up, cap_s, cap_e,
              AT, AP, ptn_off, ptn_idx, ptn_str, ptn_kind, prm,
              reg_tries, reg_give_up, flat_mem, flat_q, band_m, band_w):
    """... AT, AP: types and positions of instances already placed around the tile (anchors); ptn_off/ptn_idx/ptn_str/ptn_kind:
    partner CSR over library types and prm (A_contact, A_loose, K_half, P_max, attract tries, loose gap min, max),
    see organisation.py. Empty ptn_off disables attraction.
    reg_tries, reg_give_up: per region, tries per copy and failures in a row before a job stops (0: max_tries, give_up).
    flat_mem: per volume region, the membrane (index into the SDF vector) whose normal half of the tries align the
    molecule's flattest axis with (thin lumens and shells), or -1; flat_q: per library type, the rotation taking its
    flattest axis to +z. band_m, band_w: per region, the membrane of its shell (-1: none) and the largest distance of
    the molecule centre from it: a copy proposed next to a partner must still lie in its membrane-face shell."""
    np.random.seed(seed)
    n_types = BR.shape[0]
    res_pos = np.empty((n_types, KRES, 3)); res_n = np.zeros(n_types, np.int64)
    for ia in range(AT.shape[0]):
        _res_add(res_pos, res_n, AT[ia], AP[ia])
    use_att = ptn_off.shape[0] == n_types + 1
    wk = np.empty(256); att_stats = np.zeros(2, np.int64)
    att_tries = int(prm[4]) if use_att else 0
    glo = lo - halo; ghi = lo + size + halo
    dims, head, nxt, sid, cen, rad, cnt = grid_new(glo, ghi, H_GRID, cap_s, cap_e)
    ok = g_insert(glo, H_GRID, dims, head, nxt, sid, cen, rad, cnt, OC, OR)
    ntot = 0
    for j in range(jn.shape[0]):
        ntot += jn[j]
    otype = np.empty(ntot, np.int32); opos = np.empty((ntot, 3), np.float32); oquat = np.empty((ntot, 4), np.float32)
    placed = np.zeros(jn.shape[0], np.int64)
    out = np.zeros(c3.NMEM)
    M = np.empty((3, 3)); q = np.empty(4); t = np.empty(3); nrm = np.empty(3)
    q1 = np.empty(4); tp = np.empty(3)
    no = 0
    overflow = not ok
    for j in range(jn.shape[0]):
        if overflow:
            break
        li = jl[j]; reg = jr[j]
        a0 = st_off[li]; a1 = st_off[li + 1]
        stc = STC[a0:a1]; str_ = STR[a0:a1]
        P = np.empty((a1 - a0, 3)); ZL = np.empty(a1 - a0)
        is_mem = reg >= n_mem_regions0 and reg < n_mem_regions0 + 8      # 8 membranes; codes above are volume shells
        mi = reg - n_mem_regions0
        if is_mem:
            t0 = tri_off[mi]; t1 = tri_off[mi + 1]
            if t1 <= t0:
                continue
            for s in range(a1 - a0):
                ZL[s] = stc[s, 2]
        else:
            v0 = vox_off[reg]; v1 = vox_off[reg + 1]
            if v1 <= v0:
                continue
        n_try = reg_tries[reg] if reg_tries[reg] > 0 else max_tries
        n_give = reg_give_up[reg] if reg_give_up[reg] > 0 else give_up
        fm = -1 if is_mem else flat_mem[reg]
        fails_in_row = 0
        for inst in range(jn[j]):
            success = False
            att_k = -1
            if use_att and ptn_off[li + 1] > ptn_off[li]:
                tot = 0.0
                for k in range(ptn_off[li], ptn_off[li + 1]):
                    nn = float(res_n[ptn_idx[k]])
                    a_ = ptn_str[k] * prm[ptn_kind[k]] * nn / (prm[2] + nn)
                    wk[k - ptn_off[li]] = a_; tot += a_
                if tot > 0 and np.random.random() < min(prm[3], tot):
                    x = np.random.random() * tot
                    for k in range(ptn_off[li], ptn_off[li + 1]):
                        x -= wk[k - ptn_off[li]]
                        if x <= 0:
                            att_k = k
                            break
                    if att_k < 0:
                        att_k = ptn_off[li + 1] - 1
                    att_stats[0] += 1
            for tr in range(n_try):
                near = att_k >= 0 and tr < att_tries
                if near:
                    pt = ptn_idx[att_k]
                    if ptn_kind[att_k] == 0:
                        dd = (BR[li] + BR[pt]) * (0.45 + 0.5 * np.random.random())
                    else:
                        dd = BR[li] + BR[pt] + prm[5] + (prm[6] - prm[5]) * np.random.random()
                    _near(res_pos, res_n, pt, dd, t)
                if is_mem and near:
                    pass
                elif is_mem:
                    u = np.random.random() * CUMA[t1 - 1] if t0 == 0 else CUMA[t0 - 1] + np.random.random() * (CUMA[t1 - 1] - CUMA[t0 - 1])
                    # binary search in [t0, t1)
                    lo_i = t0; hi_i = t1 - 1
                    while lo_i < hi_i:
                        mid = (lo_i + hi_i) // 2
                        if CUMA[mid] < u:
                            lo_i = mid + 1
                        else:
                            hi_i = mid
                    a = np.random.random(); b = np.random.random()
                    if a + b > 1:
                        a = 1 - a; b = 1 - b
                    for k in range(3):
                        t[k] = TRI[lo_i, k] + a * (TRI[lo_i, 3 + k] - TRI[lo_i, k]) + b * (TRI[lo_i, 6 + k] - TRI[lo_i, k])
                if is_mem:
                    if not project_to_membrane(t, mi, nrm, G, pores, pidx, PT, PQ, PC, out):
                        continue
                    if t[0] < lo[0] or t[1] < lo[1] or t[2] < lo[2] or t[0] >= lo[0] + size or t[1] >= lo[1] + size or t[2] >= lo[2] + size:
                        continue
                    if mi == c3.M_NE and near_pore(t, G, pores, c3.PORE_R + 0.5 * c3.PERI_W + 25.0):
                        continue
                    quat_z_to(nrm, np.random.random() * 2 * np.pi, q)
                    quat_mat(q, M)
                    transform(stc, M, t, P)
                    if g_overlap(glo, H_GRID, dims, head, nxt, sid, cen, rad, P, str_, TOL):
                        continue
                    if not check_membrane(P, str_, ZL, mi, G, pores, pidx, PT, PQ, PC, out):
                        continue
                else:
                    if near:
                        if t[0] < lo[0] or t[1] < lo[1] or t[2] < lo[2] or t[0] >= lo[0] + size or t[1] >= lo[1] + size or t[2] >= lo[2] + size:
                            continue
                        if band_m[reg] >= 0:
                            c3.mem_sdfs(t[0], t[1], t[2], G, pores, pidx, PT, PQ, PC, out)
                            if abs(out[band_m[reg]]) > band_w[reg]:
                                continue
                    else:
                        vi = v0 + int(np.random.random() * (v1 - v0))
                        for k in range(3):
                            t[k] = VOXP[vi, k] + (np.random.random() - 0.5) * vox_h
                    flat = False
                    if fm >= 0 and np.random.random() < 0.5:
                        for k in range(3):
                            tp[k] = t[k]
                        if project_to_membrane(tp, fm, nrm, G, pores, pidx, PT, PQ, PC, out):
                            quat_z_to(nrm, np.random.random() * 2 * np.pi, q1)
                            quat_mul(q1, flat_q[li], q)
                            flat = True
                    if not flat:
                        rand_quat(q)
                    quat_mat(q, M)
                    transform(stc, M, t, P)
                    if g_overlap(glo, H_GRID, dims, head, nxt, sid, cen, rad, P, str_, TOL):
                        continue
                    if not check_volume(P, str_, t[0], t[1], t[2], BR[li], regmask[reg], G, pores, pidx, PT, PQ, PC, out):
                        continue
                if not g_insert(glo, H_GRID, dims, head, nxt, sid, cen, rad, cnt, P, str_):
                    overflow = True
                    break
                otype[no] = li
                for k in range(3):
                    opos[no, k] = t[k]
                for k in range(4):
                    oquat[no, k] = q[k]
                no += 1; placed[j] += 1
                if use_att:
                    _res_add(res_pos, res_n, li, t)
                    if near:
                        att_stats[1] += 1
                success = True
                break
            if overflow:
                break
            if success:
                fails_in_row = 0
            else:
                fails_in_row += 1
                if fails_in_row >= n_give:
                    break
    # ---------------- mRNA threads (1 bead / nt, 0.6 nm rise) through free cytosol
    mr_pts = np.empty((0, 3), np.float32)
    mr_off = np.zeros(1, np.int64)
    if n_mrna > 0 and not overflow:
        v0 = vox_off[cyto_region]; v1 = vox_off[cyto_region + 1]
        maxb = 0
        lens = np.empty(n_mrna, np.int64)
        for c in range(n_mrna):
            lens[c] = max(300, int(np.random.exponential(mrna_mean)))
            maxb += lens[c]
        mr_pts = np.empty((maxb, 3), np.float32)
        mr_off = np.zeros(n_mrna + 1, np.int64)
        nb_ = 0
        bead = np.empty((1, 3)); br_ = np.full(1, 0.5)
        prev = np.empty((1, 3))
        mlo = lo - 10.0; mhi = lo + size + 10.0
        for c in range(n_mrna):
            if v1 <= v0:
                mr_off[c + 1] = nb_
                continue
            start_ok = False
            for at in range(200):
                vi = v0 + int(np.random.random() * (v1 - v0))
                for k in range(3):
                    bead[0, k] = VOXP[vi, k] + (np.random.random() - 0.5) * vox_h
                if not g_overlap(glo, H_GRID, dims, head, nxt, sid, cen, rad, bead, br_, 0.0):
                    m = c3.mem_sdfs(bead[0, 0], bead[0, 1], bead[0, 2], G, pores, pidx, PT, PQ, PC, out)
                    if c3.label_from(bead[0, 0], bead[0, 1], bead[0, 2], G, pores, out, m) == c3.CYTO and c3.min_mem_dist(out, -1) > 0.5 + c3.HALF:
                        start_ok = True
                        break
            if not start_ok:
                mr_off[c + 1] = nb_
                continue
            for k in range(3):
                mr_pts[nb_, k] = bead[0, k]
            nb_ += 1; nthis = 1
            dx = np.random.normal(); dy = np.random.normal(); dz = np.random.normal()
            dn = np.sqrt(dx * dx + dy * dy + dz * dz); dx /= dn; dy /= dn; dz /= dn
            stuck = 0
            pending = False
            while nthis < lens[c] and stuck < 200:
                ex = dx + np.random.normal() * 0.9; ey = dy + np.random.normal() * 0.9; ez = dz + np.random.normal() * 0.9
                en = np.sqrt(ex * ex + ey * ey + ez * ez); ex /= en; ey /= en; ez /= en
                bead[0, 0] = mr_pts[nb_ - 1, 0] + 0.6 * ex
                bead[0, 1] = mr_pts[nb_ - 1, 1] + 0.6 * ey
                bead[0, 2] = mr_pts[nb_ - 1, 2] + 0.6 * ez
                if bead[0, 0] < mlo[0] or bead[0, 1] < mlo[1] or bead[0, 2] < mlo[2] or bead[0, 0] > mhi[0] or bead[0, 1] > mhi[1] or bead[0, 2] > mhi[2]:
                    stuck += 1; continue
                if g_overlap(glo, H_GRID, dims, head, nxt, sid, cen, rad, bead, br_, 0.0):
                    stuck += 1; continue
                m = c3.mem_sdfs(bead[0, 0], bead[0, 1], bead[0, 2], G, pores, pidx, PT, PQ, PC, out)
                if c3.label_from(bead[0, 0], bead[0, 1], bead[0, 2], G, pores, out, m) != c3.CYTO or c3.min_mem_dist(out, -1) < 0.5 + c3.HALF:
                    stuck += 1; continue
                # insert the previous bead (lagged so consecutive beads do not self-collide)
                if pending:
                    if not g_insert(glo, H_GRID, dims, head, nxt, sid, cen, rad, cnt, prev, br_):
                        overflow = True; break
                for k in range(3):
                    prev[0, k] = mr_pts[nb_ - 1, k]
                pending = True
                for k in range(3):
                    mr_pts[nb_, k] = bead[0, k]
                nb_ += 1; nthis += 1
                dx, dy, dz = ex, ey, ez; stuck = 0
            if overflow:
                break
            # insert the last two beads
            for k in range(3):
                prev[0, k] = mr_pts[nb_ - 1, k]
            g_insert(glo, H_GRID, dims, head, nxt, sid, cen, rad, cnt, prev, br_)
            if nthis > 1:
                for k in range(3):
                    prev[0, k] = mr_pts[nb_ - 2, k]
                g_insert(glo, H_GRID, dims, head, nxt, sid, cen, rad, cnt, prev, br_)
            mr_off[c + 1] = nb_
        mr_pts = mr_pts[:nb_]
        mr_off = mr_off[:n_mrna + 1]
    return otype[:no], opos[:no], oquat[:no], placed, mr_pts, mr_off, overflow, cnt[0], cnt[1], att_stats


# ------------------------------------------------------------------ helpers for obstacles
@nb.njit(cache=True)
def emit_spheres(types, pos, quat, st_off, STC, STR, BR, lo, hi):
    """Sphere-tree spheres of instances whose bounding sphere intersects box [lo, hi]."""
    n = 0
    keep = np.zeros(types.shape[0], np.bool_)
    for i in range(types.shape[0]):
        r = BR[types[i]]
        ok = True
        for a in range(3):
            if pos[i, a] + r < lo[a] or pos[i, a] - r > hi[a]:
                ok = False
        if ok:
            keep[i] = True
            n += st_off[types[i] + 1] - st_off[types[i]]
    C = np.empty((n, 3)); R = np.empty(n)
    M = np.empty((3, 3)); q = np.empty(4); t = np.empty(3)
    k = 0
    for i in range(types.shape[0]):
        if not keep[i]:
            continue
        for a in range(4):
            q[a] = quat[i, a]
        quat_mat(q, M)
        for a in range(3):
            t[a] = pos[i, a]
        a0 = st_off[types[i]]; a1 = st_off[types[i] + 1]
        transform(STC[a0:a1], M, t, C[k:k + a1 - a0])
        R[k:k + a1 - a0] = STR[a0:a1]
        k += a1 - a0
    return C, R


# ------------------------------------------------------------------ chromatin pre-pass (global, nucleus only)
@nb.njit(cache=True)
def chromatin_walk(seed, glo, ghi, G, pores, pidx, PT, PQ, PC, stc, str_, br, mask,
                   OC, OR, chain_len, starts, r_start, step, cap_s, cap_e):
    """Persistent self-avoiding nucleosome walks, one per chromosome, first fragment starting within r_start of
    starts[c] (the SPB or the chromosome's kinetochore-MT tip: Rabl-like). Returns positions, quaternions, chain id
    and a 'break' flag (1 = fragment start)."""
    np.random.seed(seed)
    dims, head, nxt, sid, cen, rad, cnt = grid_new(glo, ghi, H_GRID, cap_s, cap_e)
    g_insert(glo, H_GRID, dims, head, nxt, sid, cen, rad, cnt, OC, OR)
    ntot = 0
    for c in range(chain_len.shape[0]):
        ntot += chain_len[c]
    pos = np.empty((ntot, 3), np.float32); qu = np.empty((ntot, 4), np.float32)
    chain = np.empty(ntot, np.int32); brk = np.zeros(ntot, np.int8)
    out = np.zeros(c3.NMEM); M = np.empty((3, 3)); q = np.empty(4); t = np.empty(3)
    P = np.empty((stc.shape[0], 3))
    n = 0
    for c in range(chain_len.shape[0]):
        placed = 0; have_prev = False
        dx, dy, dz = 0.0, 0.0, 0.0
        px, py, pz = 0.0, 0.0, 0.0
        fails = 0; restarts = 0
        while placed < chain_len[c] and restarts < 20000:
            if not have_prev:
                rr = r_start if placed == 0 else 900.0
                cxx = starts[c, 0] if placed == 0 else G[3]
                cyy = starts[c, 1] if placed == 0 else G[4]
                czz = starts[c, 2] if placed == 0 else G[5]
                t[0] = cxx + (np.random.random() * 2 - 1) * rr
                t[1] = cyy + (np.random.random() * 2 - 1) * rr
                t[2] = czz + (np.random.random() * 2 - 1) * rr
            else:
                ex = np.random.normal(); ey = np.random.normal(); ez = np.random.normal()
                en = np.sqrt(ex * ex + ey * ey + ez * ez)
                ex = 0.5 * dx + ex / en; ey = 0.5 * dy + ey / en; ez = 0.5 * dz + ez / en
                en = np.sqrt(ex * ex + ey * ey + ez * ez)
                ex /= en; ey /= en; ez /= en
                t[0] = px + step * ex; t[1] = py + step * ey; t[2] = pz + step * ez
            rand_quat(q); quat_mat(q, M); transform(stc, M, t, P)
            good = not g_overlap(glo, H_GRID, dims, head, nxt, sid, cen, rad, P, str_, TOL)
            if good:
                good = check_volume(P, str_, t[0], t[1], t[2], br, mask, G, pores, pidx, PT, PQ, PC, out)
            if not good:
                fails += 1
                if fails % 40 == 0 or not have_prev:
                    have_prev = False
                    restarts += 1
                continue
            if not g_insert(glo, H_GRID, dims, head, nxt, sid, cen, rad, cnt, P, str_):
                break
            if have_prev:
                dx = (t[0] - px) / step; dy = (t[1] - py) / step; dz = (t[2] - pz) / step
            else:
                brk[n] = 1
            px, py, pz = t[0], t[1], t[2]
            have_prev = True
            for k in range(3):
                pos[n, k] = t[k]
            for k in range(4):
                qu[n, k] = q[k]
            chain[n] = c
            n += 1; placed += 1; fails = 0
    return pos[:n], qu[:n], chain[:n], brk[:n]


@nb.njit(cache=True)
def arm_walks(seed, glo, ghi, G, pores, pidx, PT, PQ, PC, stc, str_, br, mask,
              OC, OR, seg_len, seg_start, seg_target, r_start, step, cap_s, cap_e):
    """Persistent self-avoiding nucleosome walks, one per chromosome-arm segment, from within r_start of seg_start
    (the kinetochore-MT tip) towards seg_target (a telomere site on the nuclear envelope, or where chromosome XII
    enters or leaves the nucleolus). The drift towards the target grows as the arm's remaining length runs short
    relative to its distance, so arms wander freely at first and end near their targets. A walk that gets stuck
    restarts near its last nucleosome (a break in the chain). Returns positions, quaternions, segment id, break flag."""
    np.random.seed(seed)
    dims, head, nxt, sid, cen, rad, cnt = grid_new(glo, ghi, H_GRID, cap_s, cap_e)
    g_insert(glo, H_GRID, dims, head, nxt, sid, cen, rad, cnt, OC, OR)
    ntot = 0
    for c in range(seg_len.shape[0]):
        ntot += seg_len[c]
    pos = np.empty((ntot, 3), np.float32); qu = np.empty((ntot, 4), np.float32)
    chain = np.empty(ntot, np.int32); brk = np.zeros(ntot, np.int8)
    out = np.zeros(c3.NMEM); M = np.empty((3, 3)); q = np.empty(4); t = np.empty(3)
    P = np.empty((stc.shape[0], 3))
    n = 0
    for c in range(seg_len.shape[0]):
        placed = 0; have_prev = False
        dx, dy, dz = 0.0, 0.0, 0.0
        px, py, pz = seg_start[c, 0], seg_start[c, 1], seg_start[c, 2]
        fails = 0; restarts = 0
        while placed < seg_len[c] and restarts < 20000:
            if not have_prev:
                rr = r_start if placed == 0 else 30.0 + 2.0 * min(restarts, 200)
                t[0] = px + (np.random.random() * 2 - 1) * rr
                t[1] = py + (np.random.random() * 2 - 1) * rr
                t[2] = pz + (np.random.random() * 2 - 1) * rr
            else:
                gx = seg_target[c, 0] - px; gy = seg_target[c, 1] - py; gz = seg_target[c, 2] - pz
                d = np.sqrt(gx * gx + gy * gy + gz * gz) + 1e-9
                rem = seg_len[c] - placed
                w = d / (2.0 * step * np.sqrt(rem + 1.0)) - 0.25      # 0 while the target is easily within reach
                w = min(max(w, 0.0), 1.5)
                ex = np.random.normal(); ey = np.random.normal(); ez = np.random.normal()
                en = np.sqrt(ex * ex + ey * ey + ez * ez)
                ex = 0.5 * dx + ex / en + w * gx / d; ey = 0.5 * dy + ey / en + w * gy / d; ez = 0.5 * dz + ez / en + w * gz / d
                en = np.sqrt(ex * ex + ey * ey + ez * ez)
                ex /= en; ey /= en; ez /= en
                t[0] = px + step * ex; t[1] = py + step * ey; t[2] = pz + step * ez
            rand_quat(q); quat_mat(q, M); transform(stc, M, t, P)
            good = not g_overlap(glo, H_GRID, dims, head, nxt, sid, cen, rad, P, str_, TOL)
            if good:
                good = check_volume(P, str_, t[0], t[1], t[2], br, mask, G, pores, pidx, PT, PQ, PC, out)
            if not good:
                fails += 1
                if fails % 40 == 0 or not have_prev:
                    have_prev = False
                    restarts += 1
                continue
            if not g_insert(glo, H_GRID, dims, head, nxt, sid, cen, rad, cnt, P, str_):
                break
            if have_prev:
                dx = (t[0] - px) / step; dy = (t[1] - py) / step; dz = (t[2] - pz) / step
            else:
                brk[n] = 1
            px, py, pz = t[0], t[1], t[2]
            have_prev = True
            for k in range(3):
                pos[n, k] = t[k]
            for k in range(4):
                qu[n, k] = q[k]
            chain[n] = c
            n += 1; placed += 1; fails = 0
    return pos[:n], qu[:n], chain[:n], brk[:n]


W_MAX = 5.0      # cap on the fibre backbone's drift weight towards its target (1.5 left the final approach to noise)


@nb.njit(cache=True)
def guided_walks(seed, glo, ghi, G, pores, pidx, PT, PQ, PC, stc, str_, br, mask,
                 OC, OR, seg_len, guide, r_guide, link_max, cap_s, cap_e):
    """Nucleosomes along a given path (IMPROVEMENTS N2): nucleosome i (segments concatenated, walk order) goes within
    r_guide of its guide point (the relaxed coarse polymer at its bp, chromatin_polymer.py) and within link_max[i]
    of the previous nucleosome. Where it does not fit, the radius grows (up to 3 r_guide) and the reach is relaxed
    (up to 2 link_max); failing that the chain restarts near the guide (a break). Returns positions, quaternions,
    segment id, break flag, and the index (into the concatenated nucleosome list) of each placed nucleosome."""
    np.random.seed(seed)
    dims, head, nxt, sid, cen, rad, cnt = grid_new(glo, ghi, H_GRID, cap_s, cap_e)
    g_insert(glo, H_GRID, dims, head, nxt, sid, cen, rad, cnt, OC, OR)
    ntot = 0
    for c in range(seg_len.shape[0]):
        ntot += seg_len[c]
    pos = np.empty((ntot, 3), np.float32); qu = np.empty((ntot, 4), np.float32)
    chain = np.empty(ntot, np.int32); brk = np.zeros(ntot, np.int8); src = np.empty(ntot, np.int64)
    out = np.zeros(c3.NMEM); M = np.empty((3, 3)); q = np.empty(4); t = np.empty(3)
    P = np.empty((stc.shape[0], 3))
    n = 0
    o = 0
    for c in range(seg_len.shape[0]):
        if c > 0:
            o += seg_len[c - 1]
        have_prev = False
        px, py, pz = 0.0, 0.0, 0.0
        i = 0
        while i < seg_len[c]:
            gx, gy, gz = guide[o + i, 0], guide[o + i, 1], guide[o + i, 2]
            ok = False
            for rnd in range(6):
                r = r_guide * (1.0 + 0.4 * rnd)
                lk = link_max[o + i] * (1.0 + 0.2 * rnd)
                for k in range(40):
                    ex = np.random.normal(); ey = np.random.normal(); ez = np.random.normal()
                    en = np.sqrt(ex * ex + ey * ey + ez * ez) + 1e-9
                    if have_prev:
                        # within reach of the previous nucleosome, and within r of the guide point
                        rr = lk * np.random.random() ** (1.0 / 3.0)
                        t[0] = px + rr * ex / en; t[1] = py + rr * ey / en; t[2] = pz + rr * ez / en
                        lx = t[0] - gx; ly = t[1] - gy; lz = t[2] - gz
                        if lx * lx + ly * ly + lz * lz > r * r:
                            continue
                    else:
                        rr = r * np.random.random() ** (1.0 / 3.0)
                        t[0] = gx + rr * ex / en; t[1] = gy + rr * ey / en; t[2] = gz + rr * ez / en
                    rand_quat(q); quat_mat(q, M); transform(stc, M, t, P)
                    if g_overlap(glo, H_GRID, dims, head, nxt, sid, cen, rad, P, str_, TOL):
                        continue
                    if not check_volume(P, str_, t[0], t[1], t[2], br, mask, G, pores, pidx, PT, PQ, PC, out):
                        continue
                    ok = True
                    break
                if ok:
                    break
            if not ok:
                if have_prev:
                    have_prev = False          # restart this nucleosome without the reach constraint
                    continue
                i += 1                         # no room near the guide at all: skip this nucleosome
                continue
            if not g_insert(glo, H_GRID, dims, head, nxt, sid, cen, rad, cnt, P, str_):
                break
            if not have_prev:
                brk[n] = 1
            px, py, pz = t[0], t[1], t[2]
            have_prev = True
            for k in range(3):
                pos[n, k] = t[k]
            for k in range(4):
                qu[n, k] = q[k]
            chain[n] = c; src[n] = o + i
            n += 1; i += 1
    return pos[:n], qu[:n], chain[:n], brk[:n], src[:n]


@nb.njit(cache=True)
def fiber_walks(seed, glo, ghi, G, pores, pidx, PT, PQ, PC, stc, str_, br, mask,
                OC, OR, seg_len, seg_start, seg_target, r_start, delta, lp, r_fiber, link_max, cap_s, cap_e):
    # delta, link_max: per nucleosome (all segments concatenated in order): backbone advance before it and its
    # maximum centre distance from the previous nucleosome (both grow with the linker DNA in between)
    """Chromatin as a compact, irregular fibre (IMPROVEMENTS K3): a worm-like backbone per arm segment that advances
    `delta` nm per nucleosome (165 bp / delta = mass density in bp/nm) with persistence length `lp`, and nucleosomes
    placed at random within `r_fiber` of the backbone point, each within `link_max` (centre to centre) of the
    previous one. The backbone drifts towards seg_target as the remaining contour runs short (as in arm_walks); when
    nucleosomes repeatedly fail to fit, the backbone turns; when that fails too, the chain restarts nearby (a break).
    Returns positions, quaternions, segment id, break flag."""
    np.random.seed(seed)
    dims, head, nxt, sid, cen, rad, cnt = grid_new(glo, ghi, H_GRID, cap_s, cap_e)
    g_insert(glo, H_GRID, dims, head, nxt, sid, cen, rad, cnt, OC, OR)
    ntot = 0
    for c in range(seg_len.shape[0]):
        ntot += seg_len[c]
    pos = np.empty((ntot, 3), np.float32); qu = np.empty((ntot, 4), np.float32)
    chain = np.empty(ntot, np.int32); brk = np.zeros(ntot, np.int8)
    out = np.zeros(c3.NMEM); M = np.empty((3, 3)); q = np.empty(4); t = np.empty(3)
    P = np.empty((stc.shape[0], 3))
    n = 0
    o = 0
    for c in range(seg_len.shape[0]):
        if c > 0:
            o += seg_len[c - 1]
        placed = 0; have_prev = False
        bx, by, bz = seg_start[c, 0], seg_start[c, 1], seg_start[c, 2]
        px, py, pz = bx, by, bz
        ux, uy, uz = seg_target[c, 0] - bx, seg_target[c, 1] - by, seg_target[c, 2] - bz
        un = np.sqrt(ux * ux + uy * uy + uz * uz) + 1e-9; ux /= un; uy /= un; uz /= un
        restarts = 0; rounds = 0
        while placed < seg_len[c] and restarts < 20000:
            dl = delta[o + placed]; lk = link_max[o + placed]
            sig = np.sqrt(dl / lp)          # tangent noise per step: <theta^2> = 2 dl / lp
            if have_prev:
                # backbone step: worm-like tangent noise plus drift towards the target
                gx = seg_target[c, 0] - bx; gy = seg_target[c, 1] - by; gz = seg_target[c, 2] - bz
                d = np.sqrt(gx * gx + gy * gy + gz * gz) + 1e-9
                rem = seg_len[c] - placed
                w = d / (2.0 * np.sqrt(2.0 * lp * dl * (rem + 1.0)) + 1e-9) - 0.25
                w = min(max(w, 0.0), W_MAX)
                ux += sig * np.random.normal() + w * (dl / 20.0) * gx / d
                uy += sig * np.random.normal() + w * (dl / 20.0) * gy / d
                uz += sig * np.random.normal() + w * (dl / 20.0) * gz / d
                un = np.sqrt(ux * ux + uy * uy + uz * uz); ux /= un; uy /= un; uz /= un
                bx += dl * ux; by += dl * uy; bz += dl * uz
            ok = False
            for k in range(40):
                if have_prev:
                    # random point in the disc of radius r_fiber around the backbone, normal to the tangent
                    ex = np.random.normal(); ey = np.random.normal(); ez = np.random.normal()
                    dot = ex * ux + ey * uy + ez * uz
                    ex -= dot * ux; ey -= dot * uy; ez -= dot * uz
                    en = np.sqrt(ex * ex + ey * ey + ez * ez) + 1e-9
                    r = r_fiber * np.sqrt(np.random.random())
                    a = (np.random.random() - 0.5) * dl * 4.0
                    t[0] = bx + r * ex / en + a * ux; t[1] = by + r * ey / en + a * uy; t[2] = bz + r * ez / en + a * uz
                    lx = t[0] - px; ly = t[1] - py; lz = t[2] - pz
                    if lx * lx + ly * ly + lz * lz > lk * lk:
                        continue
                else:
                    rr = r_start if placed == 0 else 30.0 + 2.0 * min(restarts, 200)
                    t[0] = px + (np.random.random() * 2 - 1) * rr
                    t[1] = py + (np.random.random() * 2 - 1) * rr
                    t[2] = pz + (np.random.random() * 2 - 1) * rr
                rand_quat(q); quat_mat(q, M); transform(stc, M, t, P)
                if g_overlap(glo, H_GRID, dims, head, nxt, sid, cen, rad, P, str_, TOL):
                    continue
                if not check_volume(P, str_, t[0], t[1], t[2], br, mask, G, pores, pidx, PT, PQ, PC, out):
                    continue
                ok = True
                break
            if not ok:
                rounds += 1
                if have_prev:
                    # turn the backbone back to the last nucleosome and point it somewhere new, keeping a pull
                    # towards the target (a fully random turn erased the drift: arms ended ~800 nm from their targets)
                    bx, by, bz = px, py, pz
                    gx = seg_target[c, 0] - bx; gy = seg_target[c, 1] - by; gz = seg_target[c, 2] - bz
                    d = np.sqrt(gx * gx + gy * gy + gz * gz) + 1e-9
                    wt = 0.5 + 2.0 * min(1.0, d / (dl * (seg_len[c] - placed) + 1e-9))     # stronger as the contour runs short
                    ux += np.random.normal() + wt * gx / d
                    uy += np.random.normal() + wt * gy / d
                    uz += np.random.normal() + wt * gz / d
                    un = np.sqrt(ux * ux + uy * uy + uz * uz); ux /= un; uy /= un; uz /= un
                if rounds % 30 == 0 or not have_prev:
                    have_prev = False
                    restarts += 1
                continue
            if not g_insert(glo, H_GRID, dims, head, nxt, sid, cen, rad, cnt, P, str_):
                break
            if not have_prev:
                brk[n] = 1
                bx, by, bz = t[0], t[1], t[2]
            px, py, pz = t[0], t[1], t[2]
            have_prev = True; rounds = 0
            for k in range(3):
                pos[n, k] = t[k]
            for k in range(4):
                qu[n, k] = q[k]
            chain[n] = c
            n += 1; placed += 1
    return pos[:n], qu[:n], chain[:n], brk[:n]


class LocalGrid:
    """overlap test for pre-pass placements in a box (spatial hash of pack3d): obstacles inside the box, then every
    placed molecule; much faster than pairwise checks in Python"""

    def __init__(self, glo, ghi, OC, OR, n_new_spheres=200000):
        self.lo = np.asarray(glo, np.float64); hi = np.asarray(ghi, np.float64)
        k = np.all((OC > self.lo) & (OC < hi), axis=1)
        OC, OR = np.ascontiguousarray(OC[k], np.float64), np.ascontiguousarray(OR[k], np.float64)
        cap_s = len(OR) + n_new_spheres
        self.g = grid_new(self.lo, hi, H_GRID, cap_s, cap_s * 14)
        if len(OR):
            g_insert(self.lo, H_GRID, *self.g, OC, OR)

    def free(self, P, R):
        dims, head, nxt, sid, cen, rad, cnt = self.g
        return not g_overlap(self.lo, H_GRID, dims, head, nxt, sid, cen, rad, np.ascontiguousarray(P, np.float64),
                             np.ascontiguousarray(R, np.float64), TOL)

    def add(self, P, R):
        return g_insert(self.lo, H_GRID, *self.g, np.ascontiguousarray(P, np.float64), np.ascontiguousarray(R, np.float64))
