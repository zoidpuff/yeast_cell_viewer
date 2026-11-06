// JavaScript port of the analytic compartment geometry (same as the WGSL fill shader) for the minimap.
export class Geometry {
  G: number[]; PT: number[]; PQ: number[]; PC: number[]; pores: number[]; c: any; grid: any; n3: number;
  constructor(geo: any) {
    this.G = geo.G; this.PT = geo.PT; this.PQ = geo.PQ; this.PC = geo.PC; this.pores = geo.pores; this.c = geo.consts; this.grid = geo.grid;
    this.n3 = geo.grid.n ** 3;
  }
  cell(x: number, y: number, z: number) {
    const n = this.grid.n, ext = this.grid.ext, h = 2 * ext / n;
    const f = (v: number) => Math.min(n - 1, Math.max(0, Math.floor((v + ext) / h)));
    return (f(x) * n + f(y)) * n + f(z);
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
    const r = Math.hypot(x, y, z);
    if (r > G[0]) return 0;
    if (r > G[1]) return 1;
    if (r > G[2]) return 2;
    const cell = this.cell(x, y, z);
    const g = this.grid;
    const d = [1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1e9];
    let im = 1e9, cr = 1e9;
    for (let k = g.prim_off[cell]; k < g.prim_off[cell + 1]; k++) {
      const i = g.prim_idx[k], dd = this.sdfPrim(i, x, y, z), c = this.PC[i];
      if (c === 7) d[1] = Math.min(d[1], dd); else if (c === 0) d[3] = Math.min(d[3], dd); else if (c === 4) d[4] = Math.min(d[4], dd);
      else if (c === 5) im = Math.min(im, dd); else if (c === 6) cr = Math.min(cr, dd); else if (c === 2) d[6] = Math.min(d[6], dd);
      else if (c === 3) d[7] = Math.min(d[7], dd); else if (c === 1) d[8] = Math.min(d[8], dd);
    }
    const matrix = Math.max(im, -cr);
    if (d[1] < 0) return 4;
    const vx = x - G[3], vy = y - G[4], vz = z - G[5], dn = Math.hypot(vx, vy, vz);
    const hp = dn - (G[6] - 0.5 * PERI);
    const rc = PORE_R + 0.5 * PERI;
    let poreHit = false;
    for (let k = g.pore_off[cell]; k < g.pore_off[cell + 1]; k++) {
      const p = g.pore_idx[k], px = this.pores[3 * p], py = this.pores[3 * p + 1], pz = this.pores[3 * p + 2];
      const zz = vx * px + vy * py + vz * pz; if (zz <= 0) continue;
      const rho = Math.hypot(vx - zz * px, vy - zz * py, vz - zz * pz);
      if (rho < rc) { poreHit = true; break; }
    }
    const inPeri = Math.abs(hp) - 0.5 * PERI < 0 && !poreHit;
    if (inPeri) return 5;
    if (dn < G[6] - 0.5 * PERI) {
      if (poreHit && hp > -0.5 * PERI - 8) return 14;
      const wx = x - G[7], wy = y - G[8], wz = z - G[9];
      if (dn < G[6] - PERI && wx * wx + wy * wy + wz * wz < G[10] * G[10]) return 7;
      return 6;
    }
    if (dn < G[6]) return 14;
    if (d[4] < 0) return matrix < 0 ? 10 : 9;
    if (d[3] < 0) return 8;
    if (d[8] < 0) return 11;
    if (d[6] < 0) return 12;
    if (d[7] < 0) return 13;
    return 3;
  }
}
