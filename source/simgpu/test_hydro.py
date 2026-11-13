"""Checks of the Ewald-summed RPY mobility (hydro.py):
1. the matrix does not depend on the splitting parameter xi (real/reciprocal split is exact);
2. one sphere in a periodic box: D(L)/D0 = 1 - 2.837297 a/L + 4/3 pi (a/L)^3 (Hasimoto 1959; Yeh & Hummer 2004 use
   the leading term for finite-size corrections);
3. the matrix is symmetric positive definite for a random non-overlapping configuration.
usage: python simgpu/test_hydro.py"""
import math, sys, os
import numpy as np
import cupy as cp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hydro as H

c = 0.20895        # kT/(8 pi eta) in nm^3/ns at 30 C (water)
L = 60.0
rng = np.random.default_rng(1)
# random non-overlapping polydisperse spheres
a, x = [], []
while len(x) < 120:
    ai = rng.uniform(1.5, 6.0); xi_ = rng.uniform(0, L, 3)
    ok = True
    for xj, aj in zip(x, a):
        d = xi_ - xj; d -= L * np.rint(d / L)
        if np.linalg.norm(d) < ai + aj + 0.2:
            ok = False; break
    if ok:
        a.append(ai); x.append(xi_)
a, x = np.array(a), np.array(x)

def mob_with(scale):
    """mobility with xi scaled by `scale` (monkeypatched choose_xi)"""
    orig = H.choose_xi
    def ch(L_, tol=1e-5):
        xi, kmax = orig(L_, tol)
        xi2 = xi * scale
        s = math.sqrt(-math.log(tol))
        return xi2, 2 * xi2 * s * 1.15
    H.choose_xi = ch
    try:
        return H.mobility(x, a, L, c, tol=1e-7)
    finally:
        H.choose_xi = orig

M1, M2 = mob_with(1.0), mob_with(1.4)
err = float(cp.abs(M1 - M2).max() / cp.abs(M1).max())
print(f"xi independence: max rel diff {err:.2e}")
assert err < 1e-4, err
sym = float(cp.abs(M1 - M1.T).max() / cp.abs(M1).max())
ev = float(cp.linalg.eigvalsh(M1).min())
print(f"symmetry {sym:.1e}; min eigenvalue {ev:.3e} (self ~{c*4/3/a.max():.3e})")
assert sym < 1e-10 and ev > 0
# Hasimoto finite-size correction for one sphere
for ar in (0.02, 0.05, 0.1):
    a1 = ar * L
    M = H.mobility(np.zeros((1, 3)), np.array([a1]), L, c, tol=1e-8)
    D0 = c * 4 / 3 / a1
    got = float(M[0, 0]) / D0
    exp = 1 - 2.837297 * ar + 4 / 3 * math.pi * ar ** 3
    print(f"a/L={ar}: D/D0 {got:.6f} vs Hasimoto {exp:.6f}")
    assert abs(got - exp) < 2e-3 * ar + 1e-5
print("ok")
