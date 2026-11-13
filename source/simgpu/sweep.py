"""Diffusion sweeps for the diffusion model (IMPROVEMENTS H4): per-body D vs hydrodynamic radius in packed boxes, and
40 nm GEM tracers (Delarue et al. 2018). One job per GPU.
usage: python simgpu/sweep.py JOB OUT.json      JOB in: cyto120_{none,mf,rpy}, nuc100_mf, gem200_mf
"""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import boxes as bx
from engine import BD, msd_d

CENTERS = {"cyto": [-100, 500, -1900], "nuc": [-300, -700, 500]}


def gem_shell(R=18.5, r=3.2, n_sub=2):
    """GEM (40 nm PfV encapsulin with fluorescent cargo): spheres on an icosphere of radius R (162 for n_sub=2)"""
    t = (1 + 5 ** 0.5) / 2
    v = np.array([[-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0], [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t],
                  [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]], float)
    f = [[0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11], [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6],
         [7, 1, 8], [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9], [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1]]
    verts = [x / np.linalg.norm(x) for x in v]
    for _ in range(n_sub):
        cache, nf = {}, []
        def mid(a, b):
            k = (min(a, b), max(a, b))
            if k not in cache:
                m = verts[a] + verts[b]; verts.append(m / np.linalg.norm(m)); cache[k] = len(verts) - 1
            return cache[k]
        for a, b, c in f:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            nf += [[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]
        f = nf
    p = np.array(verts) * R
    return np.hstack([p, np.full((len(p), 1), r)]).astype(np.float32)


def run(job):
    lib = bx.Library()
    kind, rest = job.split("_")
    hyd = rest
    if kind.startswith("cyto"):
        size = float(kind[4:]); box = bx.cut_box(lib, CENTERS["cyto"], size); extra = None; total = 2000.0
    elif kind.startswith("nuc"):
        size = float(kind[3:]); box = bx.cut_box(lib, CENTERS["nuc"], size); extra = None; total = 2000.0
    elif kind.startswith("gem"):
        size = float(kind[3:]); box = bx.cut_box(lib, CENTERS["cyto"], size)
        sh = gem_shell()
        extra = [(sh, 20.0, np.array(p) * size) for p in ([0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.25, 0.75, 0.75], [0.75, 0.25, 0.75])]
        total = 3000.0
    sim = BD(box, lib, hydro=hyd, seed=7, extra_bodies=extra)
    t0 = time.time()
    times, tr = sim.run(total, 2.0)
    idx = np.where(sim.mobile)[0]
    D = msd_d(times, tr, idx, 0.1 * total, 0.5 * total) * 1e3
    D0 = bx.d_water(sim.rh[idx]) * 1e3
    return dict(job=job, n=int(sim.n), phi=box["phi"], wall_s=time.time() - t0, rh=sim.rh[idx].tolist(), D=D.tolist(),
                D0=D0.tolist(), tracer=sim.tracer[idx].tolist(),
                ids=[lib.types[t]["id"] if t >= 0 else "GEM40" for t in sim.types[idx]])


if __name__ == "__main__":
    out = run(sys.argv[1])
    json.dump(out, open(sys.argv[2], "w"))
    print(sys.argv[1], "done", out["wall_s"])
