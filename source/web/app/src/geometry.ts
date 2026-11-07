// JavaScript port of the compartment geometry (same as the WGSL fill shader) for the minimap.
import { F_VAC, F_NUC, nFields, fieldSdf, outerDepth } from './fields';

export class Geometry {
  G: number[]; PT: number[]; PQ: number[]; PC: number[]; pores: number[]; c: any; grid: any; n3: number;
  F: Float32Array | null;
  constructor(geo: any) {
    this.F = geo.F ?? null;
    this.G = geo.G; this.PT = geo.PT; this.PQ = geo.PQ; this.PC = geo.PC; this.pores = geo.pores; this.c = geo.consts; this.grid = geo.grid;
    this.n3 = geo.grid.n ** 3;
  }
  cell(x: number, y: number, z: number) {
    const n = this.grid.n, ext = this.grid.ext, h = 2 * ext / n;
    const f = (v: number) => Math.min(n - 1, Math.max(0, Math.floor((v + ext) / h)));
    return (f(x) * n + f(y)) * n + f(z);
  }
  /** signed distance to the ONM midsurface (ellipsoid), first-order approximation as in cell3d.nuc_sdf */
  nucSdf(x: number, y: number, z: number) {
    if (nFields(this.F) > F_NUC) return fieldSdf(this.F!, x, y, z, F_NUC);
    const G = this.G, dx = x - G[3], dy = y - G[4], dz = z - G[5];
    const lx = (G[14] * dx + G[15] * dy + G[16] * dz) / G[11], ly = (G[17] * dx + G[18] * dy + G[19] * dz) / G[12],
      lz = (G[20] * dx + G[21] * dy + G[22] * dz) / G[13];
    const k0 = Math.hypot(lx, ly, lz), k1 = Math.hypot(lx / G[11], ly / G[12], lz / G[13]);
    return k1 < 1e-12 ? -Math.min(G[11], G[12], G[13]) : k0 * (k0 - 1) / k1;
  }
  /** cortical-ER patch: a sheet at constant depth below the plasma membrane (cell3d.sdf_cer) */
  sdfCer(i: number, x: number, y: number, z: number, pm: number) {
    const q = this.PQ, o = 8 * i, r = Math.hypot(x, y, z) || 1e-9;
    const c = Math.min(1, Math.max(-1, (x * q[o] + y * q[o + 1] + z * q[o + 2]) / r));
    const s = Math.max(Math.acos(c) - q[o + 4], 0) * r, dr = (this.G[2] - q[o + 5]) - pm;
    return Math.sqrt(s * s + dr * dr) - q[o + 6];
  }
  sdfPrim(i: number, x: number, y: number, z: number) {
    const t = this.PT[i], q = this.PQ, o = 8 * i;
    if (t === 0) return Math.hypot(x - q[o], y - q[o + 1], z - q[o + 2]) - q[o + 3];
    if (t === 1) {
      const ax = x - q[o], ay = y - q[o + 1], az = z - q[o + 2];
      const bx = q[o + 3] - q[o], by = q[o + 4] - q[o + 1], bz = q[o + 5] - q[o + 2];
      const hh = Math.min(1, Math.max(0, (ax * bx + ay * by + az * bz) / (bx * bx + by * by + bz * bz)));
      return Math.hypot(ax - bx * hh, ay - by * hh, az - bz * hh) - q[o + 6];
    }
    if (t === 2) {
      const vx = x - q[o], vy = y - q[o + 1], vz = z - q[o + 2];
      const nx = q[o + 3], ny = q[o + 4], nz = q[o + 5];
      const zz = vx * nx + vy * ny + vz * nz;
      const rho = Math.hypot(vx - zz * nx, vy - zz * ny, vz - zz * nz);
      const s = Math.max(rho - q[o + 6], 0);
      return Math.sqrt(s * s + zz * zz) - q[o + 7];
    }
    const r = Math.hypot(x, y, z) || 1e-9;
    const c = Math.min(1, Math.max(-1, (x * q[o] + y * q[o + 1] + z * q[o + 2]) / r));
    const th = Math.acos(c);
    const s = Math.max(th - q[o + 4], 0) * r;
    const dr = r - q[o + 5];
    return Math.sqrt(s * s + dr * dr) - q[o + 6];
  }
  label(x: number, y: number, z: number): number {
    const G = this.G, PERI = this.c.PERI_W, PORE_R = this.c.PORE_R;
    const dep = outerDepth(this.F, G, x, y, z), pm = dep - (G[0] - G[2]);
    if (dep < 0) return 0;
    if (dep < G[0] - G[1]) return 1;
    if (pm < 0) return 2;
    const cell = this.cell(x, y, z);
    const g = this.grid;
    const d = [1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1e9];
    if (nFields(this.F) > F_VAC) d[3] = fieldSdf(this.F!, x, y, z, F_VAC);
    let im = 1e9, cr = 1e9;
    for (let k = g.prim_off[cell]; k < g.prim_off[cell + 1]; k++) {
      const i = g.prim_idx[k], dd = this.PT[i] === 3 ? this.sdfCer(i, x, y, z, pm) : this.sdfPrim(i, x, y, z), c = this.PC[i];
      if (c === 7) d[1] = Math.min(d[1], dd); else if (c === 0) d[3] = Math.min(d[3], dd); else if (c === 4) d[4] = Math.min(d[4], dd);
      else if (c === 5) im = Math.min(im, dd); else if (c === 6) cr = Math.min(cr, dd); else if (c === 2) d[6] = Math.min(d[6], dd);
      else if (c === 3) d[7] = Math.min(d[7], dd); else if (c === 1) d[8] = Math.min(d[8], dd);
    }
    const matrix = Math.max(im, -cr);
    if (d[1] < 0) return 4;
    const ns = this.nucSdf(x, y, z), hp = ns + 0.5 * PERI;
    const rc = PORE_R + 0.5 * PERI;
    let poreHit = false;
    for (let k = g.pore_off[cell]; k < g.pore_off[cell + 1]; k++) {
      const o = 6 * g.pore_idx[k], P = this.pores;
      const vx = x - P[o], vy = y - P[o + 1], vz = z - P[o + 2];
      const zz = vx * P[o + 3] + vy * P[o + 4] + vz * P[o + 5]; if (Math.abs(zz) > 150) continue;
      const rho = Math.hypot(vx - zz * P[o + 3], vy - zz * P[o + 4], vz - zz * P[o + 5]);
      if (rho < rc) { poreHit = true; break; }
    }
    const inPeri = Math.abs(hp) - 0.5 * PERI < 0 && !poreHit;
    if (inPeri) return 5;
    if (ns < -0.5 * PERI) {
      if (poreHit && ns > -PERI - 8) return 14;
      const wx = x - G[7], wy = y - G[8], wz = z - G[9];
      if (ns < -PERI && wx * wx + wy * wy + wz * wz < G[10] * G[10]) return 7;
      return 6;
    }
    if (ns < 0) return 14;
    if (d[4] < 0) return matrix < 0 ? 10 : 9;
    if (d[3] < 0) return 8;
    if (d[8] < 0) return 11;
    if (d[6] < 0) return 12;
    if (d[7] < 0) return 13;
    return 3;
  }
}
