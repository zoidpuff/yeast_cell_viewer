"""Rigid-body Brownian dynamics of a periodic box of the packed cell on a CUDA GPU (CuPy), for parameter sweeps and
calibration runs on the cluster. The physics follows the browser engine (web/app/src/sim, MODEL.md section 14) at
the schematic resolution: sphere-tree bodies, harmonic overlap penalty, overdamped Euler-Maruyama with the
deterministic displacement capped at 0.15 nm per step. It adds what the browser engine lacks:

- hydro='rpy': many-body hydrodynamics. Translational mobility of the mobile bodies from the Ewald-summed
  Rotne-Prager-Yamakawa tensor of spheres with the bodies' hydrodynamic radii (hydro.py), as in Ando & Skolnick 2010.
  Displacements dX = M F dt + sqrt(2 dt) L xi with L L^T = M (Cholesky); M and L are rebuilt every `hi_every` steps.
  Rotations keep the free-particle rotational diffusion (no rotational coupling).
- hydro='mf': the browser's mean-field factors (Cichocki et al. 1999) at the box's volume fraction.
- hydro='none': free-draining (no hydrodynamic interactions).
- per-body friction factors (`dscale`), used for per-protein transient-binding models.

Units nm, ns, kT. Positions are kept unwrapped for MSDs; the force kernel uses the minimum image.
"""
import math
import numpy as np
import cupy as cp

from boxes import d_water, dr_water
import hydro as HY

C_RPY = 0.20895        # kT / (8 pi eta_water) at 303.15 K in nm^3/ns

_SRC = r"""
extern "C" {
__device__ inline float wrapd(float d, float L) { return d - L * rintf(d / L); }

__global__ void world(const float* __restrict__ pos, const float* __restrict__ quat, const float* __restrict__ loc,
                      const int* __restrict__ sb, float* __restrict__ w, float* __restrict__ rel, int ns, float L) {
  int s = blockIdx.x * blockDim.x + threadIdx.x; if (s >= ns) return;
  int b = sb[s];
  float qw = quat[4*b], qx = quat[4*b+1], qy = quat[4*b+2], qz = quat[4*b+3];
  float vx = loc[4*s], vy = loc[4*s+1], vz = loc[4*s+2];
  // v' = v + 2w (u x v) + 2 u x (u x v)
  float tx = 2*(qy*vz - qz*vy), ty = 2*(qz*vx - qx*vz), tz = 2*(qx*vy - qy*vx);
  float rx = vx + qw*tx + (qy*tz - qz*ty), ry = vy + qw*ty + (qz*tx - qx*tz), rz = vz + qw*tz + (qx*ty - qy*tx);
  rel[3*s] = rx; rel[3*s+1] = ry; rel[3*s+2] = rz;
  float x = pos[3*b] + rx, y = pos[3*b+1] + ry, z = pos[3*b+2] + rz;
  w[3*s] = x - L * floorf(x / L); w[3*s+1] = y - L * floorf(y / L); w[3*s+2] = z - L * floorf(z / L);
}

__global__ void cellid(const float* __restrict__ w, int* __restrict__ cid, int ns, int nc, float h) {
  int s = blockIdx.x * blockDim.x + threadIdx.x; if (s >= ns) return;
  int cx = min(nc - 1, (int)(w[3*s] / h)), cy = min(nc - 1, (int)(w[3*s+1] / h)), cz = min(nc - 1, (int)(w[3*s+2] / h));
  cid[s] = (cx * nc + cy) * nc + cz;
}

// one thread per sphere: overlap forces against spheres of other bodies in the 27 neighbouring cells
__global__ void forces(const float* __restrict__ w, const float* __restrict__ rel, const float* __restrict__ rad,
                       const int* __restrict__ sb, const unsigned char* __restrict__ mob,
                       const int* __restrict__ order, const int* __restrict__ cstart, const int* __restrict__ cend,
                       const int* __restrict__ cid, float* __restrict__ F, float* __restrict__ T,
                       int ns, int nc, float h, float L, float krep, float eps, float range) {
  int s = blockIdx.x * blockDim.x + threadIdx.x; if (s >= ns) return;
  int b = sb[s];
  if (!mob[b]) return;
  float x = w[3*s], y = w[3*s+1], z = w[3*s+2], r = rad[s];
  int c = cid[s]; int cz = c % nc, cy = (c / nc) % nc, cx = c / (nc * nc);
  float fx = 0, fy = 0, fz = 0;
  for (int dx = -1; dx <= 1; dx++) for (int dy = -1; dy <= 1; dy++) for (int dz = -1; dz <= 1; dz++) {
    int ox = (cx + dx + nc) % nc, oy = (cy + dy + nc) % nc, oz = (cz + dz + nc) % nc;
    int cc = (ox * nc + oy) * nc + oz;
    for (int k = cstart[cc]; k < cend[cc]; k++) {
      int t = order[k]; if (sb[t] == b) continue;
      float ex = wrapd(x - w[3*t], L), ey = wrapd(y - w[3*t+1], L), ez = wrapd(z - w[3*t+2], L);
      float d2 = ex*ex + ey*ey + ez*ez, rr = r + rad[t];
      float lim = rr + range;
      if (d2 >= lim * lim || d2 < 1e-12f) continue;
      float d = sqrtf(d2), f = 0.f;
      if (d < rr) f = krep * (rr - d);
      else if (eps > 0.f) { float u = 1.5707963f * (d - rr) / range; f = -eps * 1.5707963f / range * sinf(2.f * u); }
      f /= d; fx += f * ex; fy += f * ey; fz += f * ez;
    }
  }
  if (fx == 0.f && fy == 0.f && fz == 0.f) return;
  atomicAdd(&F[3*b], fx); atomicAdd(&F[3*b+1], fy); atomicAdd(&F[3*b+2], fz);
  float px = rel[3*s], py = rel[3*s+1], pz = rel[3*s+2];
  atomicAdd(&T[3*b], py*fz - pz*fy); atomicAdd(&T[3*b+1], pz*fx - px*fz); atomicAdd(&T[3*b+2], px*fy - py*fx);
}
}
"""
_mod = None


def _qrot(q, v):
    """rotate vectors v by unit quaternions q (w, x, y, z)"""
    w, u = q[:, :1], q[:, 1:]
    t = 2 * cp.cross(u, v)
    return v + w * t + cp.cross(u, t)


def _kern(name):
    global _mod
    if _mod is None:
        _mod = cp.RawModule(code=_SRC, options=("--use_fast_math",))
    return _mod.get_function(name)


class BD:
    def __init__(self, box, lib, hydro="rpy", dt=0.2, krep=20.0, eps=0.0, viscosity=1.0, hi_every=10, seed=0,
                 dscale=None, cap=0.15, extra_bodies=None, links=None):
        """box: boxes.cut_box output (periodic cube). extra_bodies: list of (sphere array (k,4) local, rh, pos) tracers
        inserted after removing overlapping molecules (e.g. GEM nanoparticles). links: the mRNA channel path of the
        80S library frame (polysome.channel_path) to tie neighbouring ribosomes of one mRNA together with entropic
        springs (boxes.polysome_bonds; the mRNA itself is not a body); None: polysome ribosomes diffuse freely."""
        self.L = float(box["size"]); self.dt = dt; self.krep = krep; self.eps = eps; self.range = 0.6
        self.hydro = hydro; self.hi_every = hi_every; self.cap = cap
        self.rng = cp.random.default_rng(seed)
        types = list(box["type"]); pos = np.array(box["pos"]) - box["lo"]; quat = np.array(box["quat"])
        ptx = np.array(box.get("poly_tx", np.full(len(types), -1))); prk = np.array(box.get("poly_rank", np.full(len(types), -1)))
        rh = list(box["rh"]); mob = list(box["mobile"]); sph = [lib.spheres(t) for t in types]
        self.tracer = [False] * len(types)
        if extra_bodies:
            keep = np.ones(len(types), bool)
            for (s_loc, r_h, p) in extra_bodies:
                p = np.asarray(p, float); R = float(np.max(np.linalg.norm(s_loc[:, :3], axis=1) + s_loc[:, 3]))
                for i in range(len(types)):
                    d = pos[i] - p; d -= self.L * np.rint(d / self.L)
                    if np.linalg.norm(d) < R + box["rh"][i]:
                        keep[i] = False
            idx = np.where(keep)[0]
            types = [types[i] for i in idx]; pos = pos[idx]; quat = quat[idx]; rh = [rh[i] for i in idx]
            mob = [mob[i] for i in idx]; sph = [sph[i] for i in idx]; self.tracer = [False] * len(idx)
            ptx = np.r_[ptx[idx], np.full(len(extra_bodies), -1)]; prk = np.r_[prk[idx], np.full(len(extra_bodies), -1)]
            for (s_loc, r_h, p) in extra_bodies:
                types.append(-1); pos = np.vstack([pos, p]); quat = np.vstack([quat, [1, 0, 0, 0]]); rh.append(r_h)
                mob.append(True); sph.append(np.asarray(s_loc, np.float32)); self.tracer.append(True)
        self.types = np.array(types); self.n = n = len(types)
        self.rh = np.array(rh, np.float64); self.mobile = np.array(mob, bool); self.tracer = np.array(self.tracer)
        self.dscale = np.ones(n) if dscale is None else np.asarray(dscale, float)
        self.dt0 = d_water(self.rh) / viscosity * self.dscale
        self.dr0 = dr_water(self.rh) / viscosity * self.dscale
        if hydro == "mf":
            phi = box["phi"]
            self.dt0 = self.dt0 * max(0.05, 1 - 1.832 * phi - 0.219 * phi ** 2)
            self.dr0 = self.dr0 * max(0.05, 1 - 0.631 * phi - 0.726 * phi ** 2)
        self.viscosity = viscosity
        cnt = np.array([len(s) for s in sph]); off = np.concatenate([[0], np.cumsum(cnt)])
        self.ns = int(off[-1])
        loc = np.concatenate(sph).astype(np.float32)
        self.sb = cp.asarray(np.repeat(np.arange(n), cnt).astype(np.int32))
        self.loc = cp.asarray(loc); self.rad = cp.asarray(loc[:, 3].copy())
        self.pos = cp.asarray(pos, cp.float64)               # unwrapped body centres
        self.quat = cp.asarray(quat, cp.float64)
        self.mob_u8 = cp.asarray(self.mobile.astype(np.uint8))
        rmax = float(loc[:, 3].max())
        self.h = max(2 * rmax + self.range, self.L / 64)
        self.nc = max(3, int(self.L // self.h)); self.h = self.L / self.nc
        self.w = cp.zeros((self.ns, 3), cp.float32); self.rel = cp.zeros((self.ns, 3), cp.float32)
        self.cid = cp.zeros(self.ns, cp.int32)
        self.F = cp.zeros((n, 3), cp.float32); self.T = cp.zeros((n, 3), cp.float32)
        self.mi = cp.asarray(np.where(self.mobile)[0])
        self.Dt = cp.asarray(self.dt0); self.Dr = cp.asarray(self.dr0)
        self.step_n = 0; self.Lchol = None; self.M = None
        self.bonds = None
        if links is not None:
            from boxes import polysome_bonds
            bi, bj, pi, pj, k = polysome_bonds(ptx, prk, pos, quat, np.asarray(links, float))
            if len(bi):
                self.bonds = (cp.asarray(bi), cp.asarray(bj), cp.asarray(pi), cp.asarray(pj), cp.asarray(k)[:, None])
            self.n_bonds = int(len(bi))

    # ------------------------------------------------------------------ forces
    def _forces(self):
        ns, nc, bs = self.ns, self.nc, 128
        grid = ((ns + bs - 1) // bs,)
        _kern("world")(grid, (bs,), (self.pos.astype(cp.float32), self.quat.astype(cp.float32), self.loc, self.sb,
                                     self.w, self.rel, np.int32(ns), np.float32(self.L)))
        _kern("cellid")(grid, (bs,), (self.w, self.cid, np.int32(ns), np.int32(nc), np.float32(self.h)))
        order = cp.argsort(self.cid).astype(cp.int32)
        sc = self.cid[order]
        cells = cp.arange(nc ** 3, dtype=cp.int32)
        cstart = cp.searchsorted(sc, cells, side="left").astype(cp.int32)
        cend = cp.searchsorted(sc, cells, side="right").astype(cp.int32)
        self.F.fill(0); self.T.fill(0)
        _kern("forces")(grid, (bs,), (self.w, self.rel, self.rad, self.sb, self.mob_u8, order, cstart, cend, self.cid,
                                      self.F, self.T, np.int32(ns), np.int32(nc), np.float32(self.h),
                                      np.float32(self.L), np.float32(self.krep), np.float32(self.eps),
                                      np.float32(self.range)))
        if self.bonds is not None:          # entropic springs between polysome neighbours (forces and torques, kT units)
            bi, bj, pi, pj, k = self.bonds
            ri, rj = _qrot(self.quat[bi], pi), _qrot(self.quat[bj], pj)
            f = k * (self.pos[bj] + rj - self.pos[bi] - ri)               # on i, towards j's entry
            f32 = f.astype(cp.float32)
            cp.add.at(self.F, bi, f32); cp.add.at(self.F, bj, -f32)
            cp.add.at(self.T, bi, cp.cross(ri, f).astype(cp.float32)); cp.add.at(self.T, bj, cp.cross(rj, -f).astype(cp.float32))

    def _hi_update(self):
        mi = self.mi
        x = self.pos[mi] % self.L
        a = cp.asarray(self.rh)[mi]
        M = HY.mobility(x, a, self.L, C_RPY / self.viscosity, tol=1e-4, dtype=cp.float32)
        ds = cp.repeat(cp.sqrt(cp.asarray(self.dscale)[mi]), 3).astype(cp.float32)
        M = M * ds[:, None] * ds[None, :]                    # per-body friction factors keep M symmetric
        try:
            self.Lchol = cp.linalg.cholesky(M)
        except Exception:
            jitter = 1e-4 * float(cp.diag(M).mean())
            self.Lchol = cp.linalg.cholesky(M + jitter * cp.eye(M.shape[0], dtype=M.dtype))
        if not bool(cp.all(cp.isfinite(self.Lchol))):
            w_, V = cp.linalg.eigh(M)
            self.Lchol = V * cp.sqrt(cp.clip(w_, 1e-6 * float(w_.max()), None))[None, :]
        self.M = M

    # ------------------------------------------------------------------ integration
    def step(self, nsteps=1):
        dt, mi = self.dt, self.mi
        nm = len(mi)
        for _ in range(nsteps):
            self._forces()
            F = self.F[mi].astype(cp.float64)
            if self.hydro == "rpy":
                if self.Lchol is None or self.step_n % self.hi_every == 0:
                    self._hi_update()
                drift = (self.M @ F.reshape(-1).astype(cp.float32)).reshape(nm, 3).astype(cp.float64) * dt
                xi = self.rng.standard_normal(3 * nm, dtype=cp.float32)
                noise = (self.Lchol @ xi).reshape(nm, 3).astype(cp.float64) * math.sqrt(2 * dt)
            else:
                Dt = self.Dt[mi][:, None]
                drift = Dt * F * dt
                noise = cp.sqrt(2 * Dt * dt) * self.rng.standard_normal((nm, 3))
            dn = cp.linalg.norm(drift, axis=1, keepdims=True)
            drift *= cp.minimum(1.0, self.cap / cp.maximum(dn, 1e-12))
            self.pos[mi] += drift + noise
            # rotation: free rotational diffusion + torque
            Dr = self.Dr[mi][:, None]
            dth = Dr * self.T[mi].astype(cp.float64) * dt
            dthn = cp.linalg.norm(dth, axis=1, keepdims=True)
            rmax = cp.asarray(self.rh)[mi][:, None]
            dth *= cp.minimum(1.0, self.cap / cp.maximum(dthn * rmax, 1e-12))
            dth += cp.sqrt(2 * Dr * dt) * self.rng.standard_normal((nm, 3))
            ang = cp.linalg.norm(dth, axis=1, keepdims=True)
            axis = dth / cp.maximum(ang, 1e-12)
            dq = cp.concatenate([cp.cos(ang / 2), axis * cp.sin(ang / 2)], 1)
            q = self.quat[mi]
            w1, v1 = dq[:, :1], dq[:, 1:]; w2, v2 = q[:, :1], q[:, 1:]
            qn = cp.concatenate([w1 * w2 - (v1 * v2).sum(1, keepdims=True), w1 * v2 + w2 * v1 + cp.cross(v1, v2)], 1)
            self.quat[mi] = qn / cp.linalg.norm(qn, axis=1, keepdims=True)
            self.step_n += 1

    def run(self, total_ns, record_every_ns):
        """run and record unwrapped centres; returns (times ns, positions (nframes, n, 3) numpy)"""
        ev = max(1, int(round(record_every_ns / self.dt)))
        nfr = int(round(total_ns / (ev * self.dt)))
        out = np.zeros((nfr + 1, self.n, 3), np.float64)
        out[0] = cp.asnumpy(self.pos)
        for f in range(nfr):
            self.step(ev)
            out[f + 1] = cp.asnumpy(self.pos)
        return np.arange(nfr + 1) * ev * self.dt, out


def msd_d(times, traj, idx, lag_lo_ns, lag_hi_ns):
    """D (nm^2/ns) per body from the time-averaged MSD, least-squares slope of MSD vs lag over [lag_lo, lag_hi]"""
    dtf = times[1] - times[0]
    lags = np.arange(max(1, int(round(lag_lo_ns / dtf))), int(round(lag_hi_ns / dtf)) + 1)
    m = np.array([((traj[l:, idx] - traj[:-l, idx]) ** 2).sum(-1).mean(0) for l in lags])   # (nlag, nidx)
    t = lags * dtf
    # MSD = 6 D t + b
    A = np.vstack([t, np.ones_like(t)]).T
    coef = np.linalg.lstsq(A, m, rcond=None)[0]
    return coef[0] / 6.0
