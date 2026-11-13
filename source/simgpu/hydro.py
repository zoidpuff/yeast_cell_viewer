"""Many-body hydrodynamics for the GPU Brownian dynamics: the Rotne-Prager-Yamakawa (RPY) mobility of N spheres of
different radii in a periodic cube, Ewald-summed, as a dense 3N x 3N matrix on the GPU.

RPY for non-overlapping spheres i, j (radii a_i, a_j, sigma^2 = a_i^2 + a_j^2) in units c = kT / (8 pi eta):
    M_ij = c [ G(r) + sigma^2/6 lap G(r) ],   G = (I + rr)/r  (Stokeslet),  lap G = 2 (I - 3 rr) / r^3
Overlapping pairs (r < a_i + a_j) use the generalised RPY of Zuk et al. 2014 (J Fluid Mech 741:R5). Self term
kT / (6 pi eta a_i) I. Periodic sum by Hasimoto's splitting, H(k) = (1 + k^2/4xi^2) exp(-k^2/4xi^2):
  real space   G_R = (erfc/r - 2xi/sqrt(pi) e) I + (erfc/r + 2xi/sqrt(pi) e) rr,       e = exp(-xi^2 r^2)
               lapG_R = 2/r^3 [erfc + 2xi r/sqrt(pi) (1 + 4 xi^2 r^2 - 2 xi^4 r^4) e] I
                      + 2/r^3 [-3 erfc + 2xi r/sqrt(pi) (2 xi^4 r^4 - 2 xi^2 r^2 - 3) e] rr
  reciprocal   (1/V) sum_k!=0 8 pi (I - kk) / k^2 H(k) (1 - sigma^2 k^2 / 6) cos(k.r)
  self         M_ii = kT/(6 pi eta a_i) I + reciprocal(r = 0) - c [4 xi/sqrt(pi) - (2 a_i^2/6) 40 xi^3 / (3 sqrt(pi))] I
(derived with sympy from the biharmonic splitting; the total is independent of xi, which test_hydro.py checks).
The minimum image is used in real space, so the real-space cutoff must stay below L/2.
"""
import math
import numpy as np
import cupy as cp
from cupyx.scipy.special import erfc as cerfc

SQPI = math.sqrt(math.pi)


def choose_xi(L, tol=1e-5):
    """splitting parameter and reciprocal cutoff for a real-space cutoff of L/2 at relative accuracy tol"""
    s = math.sqrt(-math.log(tol))
    xi = s / (0.5 * L)
    kmax = 2 * xi * s
    return xi, kmax


def kvectors(L, kmax):
    """half-space k vectors (k != 0, one of each +-k pair) inside the cutoff, as float64 numpy (Nk, 3)"""
    n = int(math.ceil(kmax * L / (2 * math.pi)))
    r = np.arange(-n, n + 1)
    g = np.stack(np.meshgrid(r, r, r, indexing="ij"), -1).reshape(-1, 3)
    half = (g[:, 0] > 0) | ((g[:, 0] == 0) & (g[:, 1] > 0)) | ((g[:, 0] == 0) & (g[:, 1] == 0) & (g[:, 2] > 0))
    k = g[half] * (2 * math.pi / L)
    k2 = (k ** 2).sum(1)
    return k[k2 <= kmax ** 2]


def rpy_free(r, ai, aj):
    """free-space RPY coefficients (A, B): M = c (A I + B rr), for any r, including overlaps (Zuk et al. 2014)"""
    s2 = ai ** 2 + aj ** 2
    A = 1 / r + s2 / (3 * r ** 3)
    B = 1 / r - s2 / r ** 3
    # partial overlap |ai - aj| < r < ai + aj: M = 1/(6 pi eta ai aj) [...]; c = 1/(8 pi eta) -> factor 8/6 = 4/3
    ov = r < ai + aj
    if bool(cp.any(ov)):
        d2 = (ai - aj) ** 2
        Ao = (4 / 3) / (ai * aj) * (16 * r ** 3 * (ai + aj) - (d2 + 3 * r ** 2) ** 2) / (32 * r ** 3)
        Bo = (4 / 3) / (ai * aj) * 3 * (d2 - r ** 2) ** 2 / (32 * r ** 3)
        inside = r <= cp.abs(ai - aj)
        Ai = (4 / 3) / cp.maximum(ai, aj)
        A = cp.where(ov, cp.where(inside, Ai, Ao), A)
        B = cp.where(ov, cp.where(inside, 0.0, Bo), B)
    return A, B


def mobility(x, a, L, c, tol=1e-5, dtype=cp.float64):
    """dense periodic RPY mobility (3N, 3N), units of c (nm^3/ns per kT gives nm^2/ns per kT): x (N,3), a (N,) radii"""
    x = cp.asarray(x, dtype); a = cp.asarray(a, dtype)
    N = x.shape[0]
    xi, kmax = choose_xi(L, tol)
    rc = 0.5 * L
    # real space, minimum image
    d = x[:, None, :] - x[None, :, :]
    d -= L * cp.rint(d / L)
    r = cp.sqrt((d ** 2).sum(-1))
    eye = cp.eye(N, dtype=bool)
    r_safe = cp.where(eye, 1.0, r)
    rh = d / r_safe[..., None]
    ai, aj = a[:, None], a[None, :]
    s2 = ai ** 2 + aj ** 2
    xr = xi * r_safe
    e = cp.exp(-xr ** 2)
    ec = cerfc(xr)
    r3 = r_safe ** 3
    GA = ec / r_safe - 2 * xi / SQPI * e
    GB = ec / r_safe + 2 * xi / SQPI * e
    LA = 2 / r3 * (ec + 2 * xr / SQPI * (1 + 4 * xr ** 2 - 2 * xr ** 4) * e)
    LB = 2 / r3 * (-3 * ec + 2 * xr / SQPI * (2 * xr ** 4 - 2 * xr ** 2 - 3) * e)
    A = GA + s2 / 6 * LA
    B = GB + s2 / 6 * LB
    # overlapping pairs: replace the free-space part of the real-space kernel by the overlap form
    ov = (r_safe < ai + aj) & ~eye
    if bool(cp.any(ov)):
        Af, Bf = rpy_free(r_safe, ai, aj)
        A0 = 1 / r_safe + s2 / (3 * r3)
        B0 = 1 / r_safe - s2 / r3
        A = cp.where(ov, A + (Af - A0), A)
        B = cp.where(ov, B + (Bf - B0), B)
    keep = (r < rc) & ~eye
    A = cp.where(keep, A, 0.0); B = cp.where(keep, B, 0.0)
    M = cp.zeros((N, 3, N, 3), dtype)
    for p in range(3):
        for q in range(3):
            M[:, p, :, q] = (A * (p == q) + B * rh[..., p] * rh[..., q])
    # reciprocal space
    k = cp.asarray(kvectors(L, kmax), dtype)
    k2 = (k ** 2).sum(1)
    Hk = (1 + k2 / (4 * xi ** 2)) * cp.exp(-k2 / (4 * xi ** 2))
    W = 2 * 8 * math.pi * Hk / k2 / L ** 3                    # factor 2: half space
    kh = k / cp.sqrt(k2)[:, None]
    ph = x @ k.T
    C, S = cp.cos(ph), cp.sin(ph)
    a2 = a ** 2
    for p in range(3):
        for q in range(p, 3):
            P = W * ((p == q) - kh[:, p] * kh[:, q])
            X = (C * P) @ C.T + (S * P) @ S.T
            Y = (C * (P * k2)) @ C.T + (S * (P * k2)) @ S.T
            blk = X - (a2[:, None] * Y + Y * a2[None, :]) / 6
            M[:, p, :, q] += blk
            if q != p:
                M[:, q, :, p] += blk
    # self terms
    self_c = (4 / 3) / a - (4 * xi / SQPI - (2 * a2 / 6) * 40 * xi ** 3 / (3 * SQPI))
    idx = cp.arange(N)
    for p in range(3):
        M[idx, p, idx, p] += self_c
    return c * M.reshape(3 * N, 3 * N)
