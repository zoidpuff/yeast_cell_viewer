"""Encounter rates in Brownian dynamics of a packed box vs Smoluchowski (IMPROVEMENTS I5).

An encounter of bodies i and j is counted when their sphere-tree surface gap falls to <= G1 after having been
>= G2 since their last encounter (hysteresis, so rattling inside one encounter counts once). For two spheres
diffusing with relative diffusion D = D_i + D_j the steady rate of such events per pair is
    k = 4 pi D R1 R2 / (R2 - R1) / V,   R1 = a_i + a_j + G1,  R2 = a_i + a_j + G2
(flux between two absorbing shells; V box volume). The test uses each body's long-time D measured in the same run
(MSD) and its hydrodynamic radius a, so a ratio observed / predicted of ~1 means a kinetic model can use
Smoluchowski rates with in-vivo D (encounters.py); systematic deviations measure what shape, crowding (caging,
correlated re-collisions) and the anisotropic contact geometry add.

  python simgpu/encounter_bd.py [cyto|nuc] [total_ns] OUT.json
"""
import sys, os, json, time
import numpy as np
import cupy as cp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import boxes as bx
from engine import BD, msd_d
from sweep import CENTERS

G1 = float(os.environ.get("ENC_G1", 1.0)); G2 = float(os.environ.get("ENC_G2", 5.0))

_SRC = r"""
extern "C" __global__ void near_pbc(const float* __restrict__ w, const float* __restrict__ rad, const int* __restrict__ sb,
    const int* __restrict__ order, const int* __restrict__ cstart, const int* __restrict__ cend, const int* __restrict__ cid,
    int ns, int nc, float L, float gmax, int* __restrict__ cnt, int cap, int* __restrict__ ra, int* __restrict__ rb,
    float* __restrict__ rg) {
  int s = blockIdx.x * blockDim.x + threadIdx.x; if (s >= ns) return;
  int b = sb[s];
  float x = w[3*s], y = w[3*s+1], z = w[3*s+2], r = rad[s];
  int c = cid[s]; int cz = c % nc, cy = (c / nc) % nc, cx = c / (nc * nc);
  for (int dx = -1; dx <= 1; dx++) for (int dy = -1; dy <= 1; dy++) for (int dz = -1; dz <= 1; dz++) {
    int cc = (((cx + dx + nc) % nc) * nc + (cy + dy + nc) % nc) * nc + (cz + dz + nc) % nc;
    for (int k = cstart[cc]; k < cend[cc]; k++) {
      int t = order[k]; int bt = sb[t]; if (bt <= b) continue;
      float ex = x - w[3*t], ey = y - w[3*t+1], ez = z - w[3*t+2];
      ex -= L * rintf(ex / L); ey -= L * rintf(ey / L); ez -= L * rintf(ez / L);
      float g = sqrtf(ex*ex + ey*ey + ez*ez) - r - rad[t];
      if (g > gmax) continue;
      int i = atomicAdd(cnt, 1);
      if (i < cap) { ra[i] = b; rb[i] = bt; rg[i] = g; }
    }
  }
}
"""


def pair_gaps(sim, kern, buf, gmax=G2):
    """(keys b*n+bt, min gap) of body pairs with a sphere gap <= gmax, from the engine's current world spheres"""
    ns, L = sim.ns, sim.L
    rmax = float(sim.rad.max())
    nc = max(3, int(L // (2 * rmax + gmax))); h = L / nc
    w = sim.w
    cid = cp.minimum((w / h).astype(cp.int32), nc - 1)
    cid = (cid[:, 0] * nc + cid[:, 1]) * nc + cid[:, 2]
    order = cp.argsort(cid).astype(cp.int32); sc = cid[order]
    cells = cp.arange(nc ** 3, dtype=cp.int32)
    cs = cp.searchsorted(sc, cells, side="left").astype(cp.int32); ce = cp.searchsorted(sc, cells, side="right").astype(cp.int32)
    while True:
        cnt, ra, rb, rg = buf
        cnt.fill(0)
        kern(((ns + 127) // 128,), (128,), (w, sim.rad, sim.sb, order, cs, ce, cid.astype(cp.int32), np.int32(ns), np.int32(nc),
                                            np.float32(L), np.float32(gmax), cnt, np.int32(len(ra)), ra, rb, rg))
        n = int(cnt.get()[0])
        if n <= len(ra):
            break
        m = int(n * 1.3)                                  # grow the pair buffer and redo this step's search
        buf[1:] = [cp.empty(m, cp.int32), cp.empty(m, cp.int32), cp.empty(m, cp.float32)]
    key = ra[:n].astype(cp.int64) * sim.n + rb[:n]
    o = cp.lexsort(cp.stack([rg[:n], key])); key, g = key[o], rg[:n][o]
    first = cp.concatenate([cp.ones(1, bool), key[1:] != key[:-1]]) if n else cp.zeros(0, bool)
    return cp.asnumpy(key[first]), cp.asnumpy(g[first])


def run(kind="cyto", total=1000.0, size=120.0, seed=11):
    lib = bx.Library()
    box = bx.cut_box(lib, CENTERS[kind], size)
    sim = BD(box, lib, hydro="mf", seed=seed)
    kern = cp.RawModule(code=_SRC).get_function("near_pbc")
    cap = 4_000_000
    buf = [cp.zeros(1, cp.int32), cp.empty(cap, cp.int32), cp.empty(cap, cp.int32), cp.empty(cap, cp.float32)]
    nsteps = int(round(total / sim.dt)); rec = max(1, int(round(2.0 / sim.dt)))
    traj = [cp.asnumpy(sim.pos).copy()]
    contact = np.zeros(0, np.int64); events = []
    t0 = time.time()
    for s in range(nsteps):
        sim.step(1)                                      # sim.w: world spheres at the start of this step
        key, g = pair_gaps(sim, kern, buf)
        now = key[g <= G1]
        if s == 0:                                       # pairs already within G2 have not been separated: not eligible
            contact = key.copy()
        new = np.setdiff1d(now, contact, assume_unique=True)
        if len(new):
            events.append(new)
        contact = np.union1d(np.intersect1d(contact, key, assume_unique=True), now)
        if (s + 1) % rec == 0:
            traj.append(cp.asnumpy(sim.pos).copy())
    wall = time.time() - t0
    tr = np.array(traj); times = np.arange(len(tr)) * rec * sim.dt
    mob = np.where(sim.mobile)[0]
    D = np.zeros(sim.n)
    D[mob] = msd_d(times, tr, mob, 0.1 * total, 0.5 * total)            # nm^2/ns
    ev = np.concatenate(events) if events else np.zeros(0, np.int64)
    bi, bj = ev // sim.n, ev % sim.n
    T = (nsteps - 1) * sim.dt; V = size ** 3
    a = sim.rh
    m = sim.mobile
    # predicted events over all mobile pairs, binned by a_i + a_j
    ii, jj = np.triu_indices(sim.n, 1)
    keep = m[ii] & m[jj]; ii, jj = ii[keep], jj[keep]
    Dij = np.clip(D[ii], 0, None) + np.clip(D[jj], 0, None)
    R1 = a[ii] + a[jj] + G1; R2 = a[ii] + a[jj] + G2
    pred = 4 * np.pi * Dij * R1 * R2 / (R2 - R1) / V * T
    Ds = sim.dt0                                          # short-time (input, mean-field) D, nm^2/ns
    pred_s = 4 * np.pi * (Ds[ii] + Ds[jj]) * R1 * R2 / (R2 - R1) / V * T
    obs_m = m[bi] & m[bj]
    edges = np.array([0, 4, 6, 8, 10, 14, 20, 40])
    sb_obs = np.histogram((a[bi] + a[bj])[obs_m], edges)[0]
    sb_pred = np.histogram(a[ii] + a[jj], edges, weights=pred)[0]
    out = dict(kind=kind, size=size, total_ns=total, n=int(sim.n), mobile=int(m.sum()), phi=box["phi"], wall_s=wall,
               G1=G1, G2=G2, events_observed=int(obs_m.sum()), events_predicted=float(pred.sum()),
               ratio=float(obs_m.sum() / pred.sum()), events_predicted_shorttime_D=float(pred_s.sum()),
               ratio_shorttime_D=float(obs_m.sum() / pred_s.sum()),
               D_long_over_short=float(np.median(D[mob] / sim.dt0[mob])), bins_a_sum_nm=edges.tolist(), obs_by_bin=sb_obs.tolist(),
               pred_by_bin=sb_pred.tolist(), median_D_um2_s=float(np.median(D[mob]) * 1e3))
    return out


if __name__ == "__main__":
    kind = sys.argv[1] if len(sys.argv) > 1 else "cyto"
    total = float(sys.argv[2]) if len(sys.argv) > 2 else 1000.0
    o = run(kind, total)
    print(json.dumps({k: v for k, v in o.items() if not isinstance(v, list)}, indent=1))
    if len(sys.argv) > 3:
        json.dump(o, open(sys.argv[3], "w"), indent=1)
