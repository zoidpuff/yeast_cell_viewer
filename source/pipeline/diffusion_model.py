"""Predicted in-vivo diffusion coefficient of every soluble species (MODEL.md 14.5, IMPROVEMENTS H4).

D = D0(R_h) * exp(a_c + k_c R_h) * s_c
  R_h      Kirkwood radius of the library beads (as the simulation's hydro()), times the factor for trimmed IDRs
           (library/rh_scale.json, simdata.rh_scale)
  a_c, k_c crowding (excluded volume + mean-field hydrodynamics) from GPU Brownian dynamics of packed boxes
           (simgpu/sweep.py: 120 nm cytoplasm box, 100 nm nucleoplasm box with fixed chromatin; 2 us each)
  s_c      remaining slowdown fitted on the free tracers of the training papers of the diffusion benchmark
           (validation/diffusion.py; held-out papers: within-paper ratio RMSE 0.11 log10 units)
Writes library/diffusion.json: {type index: [D um^2/s, compartment]} for species placed mostly in the cytosol or the
nucleus (membrane proteins and other compartments are left out: not calibrated).
"""
import os, sys, json, math
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "simgpu"))

CROWD = {"cytosol": (-1.0423, -0.09880), "nucleus": (-1.5500, -0.14172)}
S = {"cytosol": 0.3645, "nucleus": 0.2579}
KT, ETA = 1.380649e-23 * 303.15, 0.797e-3
REGION = {0: "cytosol", 1: "nucleus", 2: "nucleus", 10: "cytosol", 11: "cytosol", 12: "cytosol", 13: "cytosol",
          24: "cytosol", 25: "cytosol", 26: "cytosol"}      # cytosolic shells next to membranes count as cytosol


def predict(rh_nm, comp):
    a, k = CROWD[comp]
    return KT / (6 * math.pi * ETA * rh_nm * 1e-9) * 1e12 * math.exp(a + k * rh_nm) * S[comp]


def write(data_dir, log=print):
    import boxes
    lib = boxes.Library(data_dir)
    sc = {}
    p = os.path.join(data_dir, "library", "rh_scale.json")
    if os.path.exists(p):
        sc = {int(k): v for k, v in json.load(open(p))["scale"].items()}
    out = {}
    for ti, t in enumerate(lib.types):
        reg = {int(r): c for r, c in t.get("regions", {}).items()}
        if not reg:
            continue
        r0 = max(reg, key=reg.get)
        comp = REGION.get(r0)
        if comp is None or reg[r0] < 0.5 * sum(reg.values()):
            continue
        rh = lib.rh(ti) * sc.get(ti, 1.0)
        out[str(ti)] = [round(predict(rh, comp), 3), comp]
    json.dump(dict(note="predicted in-vivo D (um^2/s) at 30 C; MODEL.md 14.5", types=out),
              open(os.path.join(data_dir, "library", "diffusion.json"), "w"), separators=(",", ":"))
    v = [x[0] for x in out.values()]
    log(f"diffusion.json: {len(out)} soluble species, median D {np.median(v):.1f} um^2/s")
