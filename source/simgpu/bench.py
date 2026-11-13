"""Speed and diffusion check of the GPU engine on a cytoplasm box.
usage: python simgpu/bench.py HYDRO [SIZE_NM] [TOTAL_NS] [CENTER x,y,z]"""
import sys, os, time, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cupy as cp
import boxes as bx
from engine import BD, msd_d

hyd = sys.argv[1]; size = float(sys.argv[2]) if len(sys.argv) > 2 else 120.0
total = float(sys.argv[3]) if len(sys.argv) > 3 else 1000.0
center = [float(v) for v in sys.argv[4].split(",")] if len(sys.argv) > 4 else [-100, 500, -1900]
lib = bx.Library()
box = bx.cut_box(lib, center, size)
sim = BD(box, lib, hydro=hyd, seed=1)
print(f"{box['n']} bodies, {sim.ns} spheres, phi {box['phi']:.3f}, cells {sim.nc}^3", flush=True)
sim.step(20); cp.cuda.Device().synchronize()
t0 = time.time(); sim.step(200); cp.cuda.Device().synchronize(); dtw = time.time() - t0
print(f"{200 / dtw:.0f} steps/s ({200 * sim.dt / 1e3 / dtw:.3f} us/s)", flush=True)
times, tr = sim.run(total, 2.0)
rh = sim.rh
out = {}
for name, m in (("GFP-sized (R_h<3)", rh < 3), ("R_h 3-6", (rh >= 3) & (rh < 6)), ("R_h>6", rh >= 6)):
    idx = np.where(m & sim.mobile)[0]
    D = msd_d(times, tr, idx, 0.1 * total, 0.7 * total) * 1e3          # um^2/s
    D0 = bx.d_water(rh[idx]) * 1e3
    out[name] = dict(n=len(idx), D=float(np.mean(D)), D_over_D0=float(np.mean(D / D0)))
    print(f"{name}: n={len(idx)} D={np.mean(D):.1f} um^2/s, D/D0={np.mean(D / D0):.3f}", flush=True)
print(json.dumps(out))
