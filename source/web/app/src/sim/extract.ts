// Cut a box out of the packed cell and turn its molecules into rigid bodies for the simulation.
import type { TileManager, LoadedTile } from '../tiles';
import { tkey, type Library, PAL, fetchBin, fetchJson } from '../data';
import { packColor, type V3 } from '../math';
import type { SimSystem } from './engine';
import { hydro, HYD_VOL_PER_DA, saffmanDelbruck, NT_CHARGE, MRNA_LINK_NT, MRNA_NM_PER_NT, MRNA_LP, RIBO_CHANNEL, type SimSettings, type SimMode, type TypeHydro } from './params';

// structures built explicitly (filaments, lattices, assemblies on DNA or membranes): held fixed. Newer builds flag them
// (instance flag bit 1); this list covers older data.
const FIXED_IDS = new Set(['F_actin', 'Arp23', 'CP_cap', 'TPM1_coil', 'MT_tubulin', 'gTuRC_nuc', 'gTuRC_cyto', 'ATP_dimer',
  'Pil1_dimer', 'Lsp1_dimer', 'Vac8_Nvj1', 'RVS161_RVS167_het', 'NPC', 'nucleosome']);
const R_MEM0 = 14, R_MEM1 = 21, R_NUCLEOSOME = 22, R_NPC = 23;
export const RESIDUE_R = { protein: 0.36, nucleic: 0.6 };

/** residue beads of a species that another molecule's residue can touch: a probe bead (protein radius) placed in
 *  contact along any of 26 directions overlaps none of the species' other beads. Buried beads are left out of the
 *  residue and mid collision spheres; core spheres (below) keep molecules from passing through each other. */
const exposedCache = new Map<number, Uint32Array>();
function exposedBeads(lib: Library, t: number): Uint32Array {
  const c = exposedCache.get(t);
  if (c) return c;
  const T = lib.types[t], B = lib.beads, BU = new Uint32Array(B.buffer), n = T.nb, o = T.bo;
  const rad = (i: number) => ((BU[4 * (o + i) + 3] & 255) === 1 ? RESIDUE_R.nucleic : RESIDUE_R.protein);
  const h = 1.6, key = (x: number, y: number, z: number) => `${Math.floor(x / h)},${Math.floor(y / h)},${Math.floor(z / h)}`;
  const grid = new Map<string, number[]>();
  for (let i = 0; i < n; i++) { const q = 4 * (o + i), k = key(B[q], B[q + 1], B[q + 2]); (grid.get(k) ?? grid.set(k, []).get(k)!).push(i); }
  const dirs: number[][] = [];
  for (let a = -1; a <= 1; a++) for (let b = -1; b <= 1; b++) for (let c2 = -1; c2 <= 1; c2++) {
    if (a || b || c2) { const l = Math.hypot(a, b, c2); dirs.push([a / l, b / l, c2 / l]); }
  }
  const out: number[] = [], rp = RESIDUE_R.protein;
  for (let i = 0; i < n; i++) {
    const q = 4 * (o + i), x = B[q], y = B[q + 1], z = B[q + 2], dc = rad(i) + rp;
    let exposed = false;
    for (const d of dirs) {
      const px = x + d[0] * dc, py = y + d[1] * dc, pz = z + d[2] * dc;
      let free = true;
      const cx = Math.floor(px / h), cy = Math.floor(py / h), cz = Math.floor(pz / h);
      for (let ax = -1; ax <= 1 && free; ax++) for (let ay = -1; ay <= 1 && free; ay++) for (let az = -1; az <= 1 && free; az++) {
        const lst = grid.get(`${cx + ax},${cy + ay},${cz + az}`);
        if (!lst) continue;
        for (const j of lst) {
          if (j === i) continue;
          const r = 4 * (o + j), m = rad(j) + rp - 0.1;          // 0.1 nm tolerance: probes may graze neighbours
          const ex = B[r] - px, ey = B[r + 1] - py, ez = B[r + 2] - pz;
          if (ex * ex + ey * ey + ez * ez < m * m) { free = false; break; }
        }
      }
      if (free) { exposed = true; break; }
    }
    if (exposed) out.push(i);
  }
  const res = new Uint32Array(out);
  exposedCache.set(t, res);
  return res;
}

/** core spheres: the species' sphere tree shrunk by CORE_INSET (inside the residue surface by about 1 nm). With
 *  surface spheres only, molecules are hollow shells, and a strong attraction pulled them through each other (the mid
 *  model at 1.5 kT: cells packed far beyond capacity, missed neighbours, molecules faster than in water). Uncharged,
 *  no native contacts. */
const CORE_INSET = 1.0;
function coreSpheres(lib: Library, t: number): Float32Array {
  const T = lib.types[t], out: number[] = [];
  for (let j = 0; j < T.sn; j++) {
    const q = 4 * (T.so + j), r = lib.st[q + 3] - CORE_INSET;
    if (r >= 0.4) out.push(lib.st[q], lib.st[q + 1], lib.st[q + 2], r);
  }
  return new Float32Array(out);
}

// ---------------------------------------------------------------- optional per-bead data for the simulation
/** charges.bin: int8 per library bead (+1 Lys/Arg, -1 Asp/Glu, -1 per nucleotide phosphate before condensation, 0 else);
 *  contacts.json + contacts.bin: native residue contacts of predicted complexes between library types (AlphaFold DB
 *  dimers), uint16 bead-index pairs. Built by the pipeline (pipeline/simdata.py); absent in older data. */
export interface SimData {
  charges: Int8Array | null;
  rhScale: Map<number, number>;          // per type: R_h factor for trimmed disordered regions (library/rh_scale.json)
  memScale: number[];                     // per membrane (PM, ER, NE, VAC, OM, IM, GOLGI, PEROX): factor on the
                                          // Saffman-Delbrueck D, calibrated on yeast measurements (library/membrane_diffusion.json)
  contacts: Map<number, { other: number; mine: Uint16Array; theirs: Uint16Array }[]>;
  nPairs: number;
}
let simData: Promise<SimData> | null = null;
export function loadSimData(): Promise<SimData> {
  simData ??= (async () => {
    const out: SimData = { charges: null, rhScale: new Map(), memScale: new Array(8).fill(1), contacts: new Map(), nPairs: 0 };
    try {
      const j = await fetchJson('library/rh_scale.json');
      for (const [t, v] of Object.entries(j.scale as Record<string, number>)) out.rhScale.set(Number(t), v);
    } catch { /* not in this build */ }
    try {
      const j = await fetchJson('library/membrane_diffusion.json');
      (j.scale as number[]).forEach((v, m) => { if (m < 8 && v > 0) out.memScale[m] = v; });
    } catch { /* not in this build */ }
    try { out.charges = new Int8Array(await fetchBin('library/charges.bin')); } catch { /* not in this build */ }
    try {
      const j = await fetchJson('library/contacts.json'), bin = new Uint16Array(await fetchBin('library/contacts.bin'));
      for (const [a, b, off, cnt] of j.pairs as number[][]) {
        const A = bin.subarray(2 * off, 2 * (off + cnt));
        const ia = new Uint16Array(cnt), ib = new Uint16Array(cnt);
        for (let k = 0; k < cnt; k++) { ia[k] = A[2 * k]; ib[k] = A[2 * k + 1]; }
        const add = (t: number, o: number, m: Uint16Array, th: Uint16Array) => {
          (out.contacts.get(t) ?? out.contacts.set(t, []).get(t)!).push({ other: o, mine: m, theirs: th });
        };
        if (a === b) {      // homodimer: each copy's interface beads meet the other copy's
          const m = new Uint16Array(2 * cnt), th = new Uint16Array(2 * cnt);
          m.set(ia); m.set(ib, cnt); th.set(ib); th.set(ia, cnt);
          add(a, a, m, th);
        } else { add(a, b, ia, ib); add(b, a, ib, ia); }
        out.nPairs++;
      }
    } catch { /* not in this build */ }
    return out;
  })();
  return simData;
}

// ---------------------------------------------------------------- collision spheres of a type, per model
/** spheres of one library type in one model: local xyzr; for every residue bead the sphere it belongs to (its charge
 *  and native contacts go there); per sphere the summed charge */
interface TypeSph { local: Float32Array; beadTo: Int32Array; n: number; nReg: number }   // spheres from nReg on are cores
const sphCache = new Map<string, TypeSph>();
const MID_PER = [4, 2];                 // residues / nucleotides per mid-model cluster
const MID_RMS_F = 0.8;                  // mid cluster radius = bead radius + 0.8 x rms distance of its members from the centre:
                                        // the cluster surface lies 0.13-0.22 nm (mean) outside the residue surface for
                                        // typical proteins, 0.3 nm for the ribosome (rays toward the centre, 13 types)

function nearestOf(P: Float32Array, nP: number, x: number, y: number, z: number, h: number, grid: Map<string, number[]>, scale?: Float32Array): number {
  const cx = Math.floor(x / h), cy = Math.floor(y / h), cz = Math.floor(z / h);
  let best = -1, bd = Infinity;
  for (let ring = 1; ring <= 4 && best < 0; ring++) {
    for (let a = -ring; a <= ring; a++) for (let b = -ring; b <= ring; b++) for (let c = -ring; c <= ring; c++) {
      const l = grid.get(`${cx + a},${cy + b},${cz + c}`);
      if (!l) continue;
      for (const j of l) {
        const dx = P[4 * j] - x, dy = P[4 * j + 1] - y, dz = P[4 * j + 2] - z;
        let d = dx * dx + dy * dy + dz * dz; if (scale) d /= scale[j] * scale[j];
        if (d < bd) { bd = d; best = j; }
      }
    }
  }
  if (best < 0) for (let j = 0; j < nP; j++) {
    const dx = P[4 * j] - x, dy = P[4 * j + 1] - y, dz = P[4 * j + 2] - z;
    let d = dx * dx + dy * dy + dz * dz; if (scale) d /= scale[j] * scale[j];
    if (d < bd) { bd = d; best = j; }
  }
  return best;
}
function hashOf(P: Float32Array, n: number, h: number): Map<string, number[]> {
  const g = new Map<string, number[]>();
  for (let j = 0; j < n; j++) {
    const k = `${Math.floor(P[4 * j] / h)},${Math.floor(P[4 * j + 1] / h)},${Math.floor(P[4 * j + 2] / h)}`;
    (g.get(k) ?? g.set(k, []).get(k)!).push(j);
  }
  return g;
}

function typeSpheres(lib: Library, t: number, mode: SimMode): TypeSph {
  const key = mode + ':' + t, c = sphCache.get(key);
  if (c) return c;
  const T = lib.types[t], B = lib.beads, BU = new Uint32Array(B.buffer), n = T.nb, o = T.bo;
  const kindOf = (i: number) => BU[4 * (o + i) + 3] & 255;
  const beadTo = new Int32Array(n);
  let local: Float32Array;
  if (mode === 'schematic') {
    local = lib.st.slice(4 * T.so, 4 * (T.so + T.sn));
    const radius = new Float32Array(T.sn); for (let j = 0; j < T.sn; j++) radius[j] = local[4 * j + 3];
    const g = hashOf(local, T.sn, 3);
    for (let i = 0; i < n; i++) beadTo[i] = nearestOf(local, T.sn, B[4 * (o + i)], B[4 * (o + i) + 1], B[4 * (o + i) + 2], 3, g, radius);
  } else if (mode === 'residue') {
    const ex = exposedBeads(lib, t), core = coreSpheres(lib, t);
    local = new Float32Array(4 * ex.length + core.length);
    ex.forEach((b, j) => { const q = 4 * (o + b); local.set([B[q], B[q + 1], B[q + 2], kindOf(b) === 1 ? RESIDUE_R.nucleic : RESIDUE_R.protein], 4 * j); });
    local.set(core, 4 * ex.length);
    beadTo.fill(-1); ex.forEach((b, j) => { beadTo[b] = j; });
    const g = hashOf(local, ex.length, 1.5);       // buried beads: their nearest exposed bead carries their charge
    for (let i = 0; i < n; i++) if (beadTo[i] < 0) beadTo[i] = nearestOf(local, ex.length, B[4 * (o + i)], B[4 * (o + i) + 1], B[4 * (o + i) + 2], 1.5, g);
  } else {
    // mid: k-means clusters of ~4 residues (2 nucleotides) per kind, seeded by runs of the Morton-ordered beads
    const cl: number[] = [], assign = new Int32Array(n).fill(-1);
    for (const kind of [0, 1]) {
      const idx: number[] = []; for (let i = 0; i < n; i++) if (kindOf(i) === kind) idx.push(i);
      if (!idx.length) continue;
      let k = Math.ceil(idx.length / MID_PER[kind]);
      let lab = idx.map((_, j) => Math.floor(j * k / idx.length));
      let C = new Float32Array(4 * k);
      for (let it = 0; it < 5; it++) {
        const acc = new Float64Array(4 * k);
        idx.forEach((b, j) => { const q = 4 * (o + b), c4 = 4 * lab[j]; acc[c4] += B[q]; acc[c4 + 1] += B[q + 1]; acc[c4 + 2] += B[q + 2]; acc[c4 + 3]++; });
        const keep: number[] = []; for (let j = 0; j < k; j++) if (acc[4 * j + 3] > 0) keep.push(j);
        C = new Float32Array(4 * keep.length);
        keep.forEach((j, m) => { C[4 * m] = acc[4 * j] / acc[4 * j + 3]; C[4 * m + 1] = acc[4 * j + 1] / acc[4 * j + 3]; C[4 * m + 2] = acc[4 * j + 2] / acc[4 * j + 3]; });
        k = keep.length;
        const g = hashOf(C, k, 1.5);
        lab = idx.map(b => nearestOf(C, k, B[4 * (o + b)], B[4 * (o + b) + 1], B[4 * (o + b) + 2], 1.5, g));
      }
      // final centres and radii
      const acc = new Float64Array(4 * k), ss = new Float64Array(k);
      idx.forEach((b, j) => { const q = 4 * (o + b), c4 = 4 * lab[j]; acc[c4] += B[q]; acc[c4 + 1] += B[q + 1]; acc[c4 + 2] += B[q + 2]; acc[c4 + 3]++; });
      const map = new Int32Array(k).fill(-1);
      for (let j = 0; j < k; j++) if (acc[4 * j + 3] > 0) {
        map[j] = (cl.length / 4); const m = acc[4 * j + 3];
        cl.push(acc[4 * j] / m, acc[4 * j + 1] / m, acc[4 * j + 2] / m, 0);
      }
      idx.forEach((b, j) => {
        const c = map[lab[j]], q = 4 * (o + b);
        const dx = B[q] - cl[4 * c], dy = B[q + 1] - cl[4 * c + 1], dz = B[q + 2] - cl[4 * c + 2];
        ss[lab[j]] += dx * dx + dy * dy + dz * dz; assign[b] = c;
      });
      for (let j = 0; j < k; j++) if (map[j] >= 0) {
        cl[4 * map[j] + 3] = (kind === 1 ? RESIDUE_R.nucleic : RESIDUE_R.protein) + MID_RMS_F * Math.sqrt(ss[j] / acc[4 * j + 3]);
      }
    }
    // keep the clusters that hold an exposed residue (the others' beads go to the nearest kept cluster), plus cores
    const ex = exposedBeads(lib, t), exposed = new Uint8Array(cl.length / 4), core = coreSpheres(lib, t);
    for (const b of ex) exposed[assign[b]] = 1;
    const newIdx = new Int32Array(cl.length / 4).fill(-1); let m = 0;
    for (let j = 0; j < cl.length / 4; j++) if (exposed[j]) newIdx[j] = m++;
    local = new Float32Array(4 * m + core.length);
    for (let j = 0; j < cl.length / 4; j++) if (newIdx[j] >= 0) local.set(cl.slice(4 * j, 4 * j + 4), 4 * newIdx[j]);
    local.set(core, 4 * m);
    const g = hashOf(local, m, 1.5);
    for (let i = 0; i < n; i++) {
      const c = newIdx[assign[i]];
      beadTo[i] = c >= 0 ? c : nearestOf(local, m, B[4 * (o + i)], B[4 * (o + i) + 1], B[4 * (o + i) + 2], 1.5, g);
    }
  }
  const nCore = mode === 'schematic' ? 0 : coreSpheres(lib, t).length / 4;
  const out = { local, beadTo, n: local.length / 4, nReg: local.length / 4 - nCore };
  sphCache.set(key, out);
  return out;
}
/** summed charge per sphere (e; nucleotides scaled by counterion condensation) */
const chargeCache = new Map<string, Float32Array>();
function sphereCharges(lib: Library, t: number, mode: SimMode, ts: TypeSph, charges: Int8Array): Float32Array {
  const key = mode + ':' + t, c = chargeCache.get(key);
  if (c) return c;
  const T = lib.types[t], BU = new Uint32Array(lib.beads.buffer), q = new Float32Array(ts.n);
  for (let i = 0; i < T.nb; i++) {
    const z = charges[T.bo + i]; if (!z) continue;
    q[ts.beadTo[i]] += (BU[4 * (T.bo + i) + 3] & 255) === 1 ? -NT_CHARGE * z : z;
  }
  chargeCache.set(key, q);
  return q;
}

export interface SimBox { center: V3; size: number }
export interface SimMeta {
  types: number[]; region: Uint8Array; hydro: TypeHydro[]; box: SimBox; nMobile: number; nFixed: number; nShell: number;
  boundary: 'periodic' | 'walls' | 'embedded'; hasMembrane: boolean; nMem2d: number;
  charges: boolean; nCharged: number; netCharge: number;          // electrostatics data, charged spheres, net charge of the box (e)
  contacts: boolean; nNativePairs: number; nNativeBodies: number;  // native contacts: type pairs present, bodies with a partner type in the box
  nLinks: number; nPolyRibo: number;                              // polysome mRNA linkers inside the box, ribosomes in them
}

const i16f = (x: number) => (((x << 16) >> 16) / 32767) * 0.70710678;
function quatOf(w: Uint32Array, i: number): [number, number, number, number] {
  const w0 = w[4 * i], w2 = w[4 * i + 2], w3 = w[4 * i + 3];
  const a = i16f(w2 >>> 16), b = i16f(w3 & 0xFFFF), c = i16f(w3 >>> 16);
  const big = Math.sqrt(Math.max(0, 1 - a * a - b * b - c * c));
  const qi = (w0 >>> 14) & 3;
  const q = qi === 0 ? [big, a, b, c] : qi === 1 ? [a, big, b, c] : qi === 2 ? [a, b, big, c] : [a, b, c, big];
  return [q[1], q[2], q[3], q[0]];        // (w,x,y,z) -> (x,y,z,w)
}

/** tile indices covering the box */
export function boxTiles(T: TileManager, box: SimBox): number[] {
  const h = box.size / 2, out: number[] = [];
  const lo = box.center.map(c => Math.floor((c - h) / 200)), hi = box.center.map(c => Math.floor((c + h) / 200));
  for (let x = lo[0]; x <= hi[0]; x++) for (let y = lo[1]; y <= hi[1]; y++) for (let z = lo[2]; z <= hi[2]; z++) {
    const i = T.idx.key.get(tkey(x, y, z));
    if (i !== undefined) out.push(i);
  }
  return out;
}

export function extract(T: TileManager, lib: Library, box: SimBox, set: SimSettings, cache: Map<number, TypeHydro>, data: SimData): { sys: SimSystem; meta: SimMeta } {
  const h = box.size / 2, c = box.center;
  const lo: V3 = [c[0] - h, c[1] - h, c[2] - h], hi: V3 = [c[0] + h, c[1] + h, c[2] + h];
  const tu = new Uint32Array(lib.typeInfo);
  const picks: { t: LoadedTile; i: number; type: number; region: number; flags: number; p: V3 }[] = [];
  for (const ti of boxTiles(T, box)) {
    const t = T.loaded.get(ti);
    if (!t) continue;
    for (let i = 0; i < t.nInst; i++) {
      const d = T.decode(t, i);
      if (d.p[0] < lo[0] || d.p[1] < lo[1] || d.p[2] < lo[2] || d.p[0] >= hi[0] || d.p[1] >= hi[1] || d.p[2] >= hi[2]) continue;
      picks.push({ t, i, ...d });
    }
  }
  const n = picks.length, mode = set.mode;
  const tsOf = new Map<number, TypeSph>();
  for (const k of picks) if (!tsOf.has(k.type)) tsOf.set(k.type, typeSpheres(lib, k.type, mode));
  let nSph = 0, nBead = 0;
  for (const k of picks) { nSph += tsOf.get(k.type)!.n; nBead += lib.types[k.type].nb; }
  const sys: SimSystem = {
    n, type: new Uint32Array(n), pos: new Float32Array(3 * n), quat: new Float32Array(4 * n), mobile: new Uint8Array(n),
    dt: new Float32Array(n), dr: new Float32Array(n), wall: new Uint8Array(n),
    sphOff: new Uint32Array(n), sphCnt: new Uint32Array(n), sphLocal: new Float32Array(4 * nSph), sphProp: new Uint32Array(4 * nSph),
    native: new Uint32Array(0), bondOff: new Uint32Array(n), bondCnt: new Uint32Array(n), bonds: new Float32Array(0),
    beadOff: new Uint32Array(n), beadCnt: new Uint32Array(n), beadLocal: new Float32Array(4 * nBead), beadKind: new Uint32Array(nBead),
    color: new Uint32Array(2 * n), lo, hi, phi: 0, tether: new Uint8Array(n), mem2d: new Uint8Array(n), boundary: 'walls',
  };
  // boundary: periodic suits bulk compartments; a membrane crossing the box does not continue across its faces, so
  // boxes with membranes default to the embedded boundary
  const hasMembrane = picks.some(k => k.region >= R_MEM0 && k.region <= R_MEM1);
  const boundary = set.boundary === 'auto' ? (hasMembrane ? 'embedded' : 'periodic') : set.boundary;
  sys.boundary = boundary;
  // volume fraction: hydrated volume of every molecule in the box (moving or not) over the box volume (membranes and
  // other compartments inside the box are counted as free volume, so mixed boxes get a low estimate)
  let vol = 0;
  for (const k of picks) vol += lib.types[k.type].mass * HYD_VOL_PER_DA;
  sys.phi = Math.min(0.5, vol / box.size ** 3);
  // native contacts between the types present: per type, one shared list per sphere (partner key, energy share)
  const present = new Set(picks.map(k => k.type));
  const natOff = new Map<number, { off: Uint32Array; cnt: Uint32Array }>();
  const natWords: number[] = [];
  let nNativePairs = 0;
  const f32u = new Float32Array(1), u32f = new Uint32Array(f32u.buffer);
  if (data.contacts.size) {
    for (const t of present) {
      const groups = (data.contacts.get(t) ?? []).filter(g => present.has(g.other));
      if (!groups.length) continue;
      const mine = tsOf.get(t)!, lists = new Map<number, Map<number, number>>();     // sphere -> partner key -> share
      for (const g of groups) {
        const theirs = tsOf.get(g.other)!;
        const pairs = new Set<number>();
        for (let k = 0; k < g.mine.length; k++) {
          const a = mine.beadTo[g.mine[k]], b = theirs.beadTo[g.theirs[k]];
          if (a >= 0 && b >= 0) pairs.add(a * 65536 + b);
        }
        if (g.other >= t) nNativePairs++;
        // a homodimer's list holds both directions of each pair; its interface energy is shared over the distinct pairs
        const nUnique = g.other === t ? new Set([...pairs].map(v => { const a = Math.floor(v / 65536), b = v % 65536; return a < b ? a * 65536 + b : b * 65536 + a; })).size : pairs.size;
        for (const v of pairs) {
          const a = Math.floor(v / 65536), b = v % 65536;
          const l = lists.get(a) ?? lists.set(a, new Map()).get(a)!;
          l.set((g.other << 16) | b, 1 / Math.max(1, nUnique));
        }
      }
      const off = new Uint32Array(mine.n), cnt = new Uint32Array(mine.n);
      for (const [a, l] of lists) {
        off[a] = natWords.length / 2; cnt[a] = l.size;
        for (const [key, share] of l) { f32u[0] = share; natWords.push(key >>> 0, u32f[0]); }
      }
      natOff.set(t, { off, cnt });
    }
  }
  sys.native = new Uint32Array(natWords);
  const region = new Uint8Array(n), hyd: TypeHydro[] = [];
  const B = lib.beads, BU = new Uint32Array(B.buffer);
  let so = 0, bo = 0, nMob = 0, nShell = 0, nMem2d = 0, nCharged = 0, netCharge = 0, nNativeBodies = 0;
  picks.forEach((k, j) => {
    const Ty = lib.types[k.type];
    sys.type[j] = k.type; region[j] = k.region;
    sys.pos.set(k.p, 3 * j); sys.quat.set(quatOf(k.t.inst, k.i), 4 * j);
    const membrane = k.region >= R_MEM0 && k.region <= R_MEM1;
    const scaffold = (k.flags & 2) !== 0 || FIXED_IDS.has(Ty.id) || k.region === R_NUCLEOSOME || k.region === R_NPC;
    const inShell = boundary === 'embedded' && k.p.some((v, a) => v < lo[a] + set.shell || v > hi[a] - set.shell);
    const mem2d = membrane && !scaffold && set.memFluid;     // membrane proteins: lateral diffusion in their membrane
    const mob = (!membrane || mem2d) && !scaffold;
    sys.mobile[j] = mob ? 1 : 0; sys.wall[j] = mob && !membrane ? 1 : 0; sys.tether[j] = mob && inShell ? 1 : 0; sys.mem2d[j] = mem2d ? 1 : 0;
    if (mem2d) nMem2d++;
    if (mob && !inShell) nMob++; else if (mob) nShell++;
    let hy = hydro(lib, k.type, cache);
    const sc = data.rhScale.get(k.type);
    if (sc) hy = { rh: hy.rh * sc, dt: hy.dt / sc, dr: hy.dr / (sc * sc * sc) };   // trimmed IDRs restored as coils
    hyd.push(hy);
    if (mem2d) {           // fluid-bilayer value times the membrane's calibration factor (crowding, domains, wall)
      const sd = saffmanDelbruck(hy.rh, set.memVisc), f = set.memCalibrated ? data.memScale[k.region - R_MEM0] : 1;
      sys.dt[j] = sd.dt * f; sys.dr[j] = sd.dr * f;
    }
    else { sys.dt[j] = hy.dt; sys.dr[j] = hy.dr; }
    // collision spheres and their properties: charge, native-contact list, type << 16 | local id
    const ts = tsOf.get(k.type)!, q = data.charges ? sphereCharges(lib, k.type, mode, ts, data.charges) : null, nat = natOff.get(k.type);
    if (nat) nNativeBodies++;
    sys.sphOff[j] = so;
    sys.sphLocal.set(ts.local, 4 * so);
    for (let m = 0; m < ts.n; m++, so++) {
      const z = q ? q[m] : 0;
      if (z) { nCharged++; netCharge += z; }
      f32u[0] = z; sys.sphProp[4 * so] = u32f[0];
      if (nat) { sys.sphProp[4 * so + 1] = nat.off[m]; sys.sphProp[4 * so + 2] = nat.cnt[m]; }
      // key of native contacts (no type has them past 65534 spheres); cores: 0xFFFFFFFE
      sys.sphProp[4 * so + 3] = m >= ts.nReg ? 0xFFFFFFFE : m < 0xFFFF ? ((k.type << 16) | m) >>> 0 : 0xFFFFFFFF;
    }
    sys.sphCnt[j] = so - sys.sphOff[j];
    // drawing beads
    sys.beadOff[j] = bo;
    for (let b = 0; b < Ty.nb; b++, bo++) {
      const q4 = 4 * (Ty.bo + b);
      sys.beadLocal[4 * bo] = B[q4]; sys.beadLocal[4 * bo + 1] = B[q4 + 1]; sys.beadLocal[4 * bo + 2] = B[q4 + 2];
      sys.beadKind[bo] = BU[q4 + 3];
    }
    sys.beadCnt[j] = Ty.nb;
    const c0 = tu[16 * k.type + 8], c1 = tu[16 * k.type + 9];
    sys.color[2 * j] = c0 || packColor(PAL[k.region] ?? '#7EA7D3'); sys.color[2 * j + 1] = c1;
  });
  const links = set.mrna ? polysomeLinks(lib, sys) : { nLinks: 0, nPolyRibo: 0 };
  return {
    sys,
    meta: { types: Array.from(sys.type), region, hydro: hyd, box, nMobile: nMob, nFixed: n - nMob - nShell, nShell, boundary, hasMembrane, nMem2d,
      charges: !!data.charges, nCharged, netCharge, contacts: data.contacts.size > 0, nNativePairs, nNativeBodies, ...links },
  };
}

/** Polysome mRNA linkers: the packer placed each next ribosome of a polysome so that its mRNA entry lies 3-12 nm beyond
 *  the previous ribosome's exit, within 0.7 rad of the exit direction, with its outward entry direction pointing back
 *  along the linker. Pairs of ribosomes in the box meeting all three conditions (to the tile quantisation) are linked
 *  by a worm-like chain of the mean linker length (params: MRNA_*). Bonds are stored on both bodies. */
function polysomeLinks(lib: Library, sys: SimSystem): { nLinks: number; nPolyRibo: number } {
  const rt = lib.byId.get('ribosome_80S');
  const ribos: number[] = [];
  if (rt !== undefined) for (let i = 0; i < sys.n; i++) if (sys.type[i] === rt && sys.mobile[i]) ribos.push(i);
  const C = RIBO_CHANNEL, w = (i: number, v: number[]) => qrotJs(sys.quat.subarray(4 * i, 4 * i + 4), v[0], v[1], v[2]);
  const ex = ribos.map(i => { const v = w(i, C.exit); return [sys.pos[3 * i] + v[0], sys.pos[3 * i + 1] + v[1], sys.pos[3 * i + 2] + v[2]]; });
  const en = ribos.map(i => { const v = w(i, C.entry); return [sys.pos[3 * i] + v[0], sys.pos[3 * i + 1] + v[1], sys.pos[3 * i + 2] + v[2]]; });
  const exd = ribos.map(i => w(i, C.exitDir)), end = ribos.map(i => w(i, C.entryDir));
  const pairs: [number, number][] = [];
  for (let a = 0; a < ribos.length; a++) for (let b = 0; b < ribos.length; b++) {
    if (a === b) continue;
    const d = [en[b][0] - ex[a][0], en[b][1] - ex[a][1], en[b][2] - ex[a][2]], L = Math.hypot(d[0], d[1], d[2]);
    if (L < C.linkMin - 0.05 || L > C.linkMax + 0.05) continue;
    const u = d.map(v => v / L);
    if (u[0] * exd[a][0] + u[1] * exd[a][1] + u[2] * exd[a][2] < Math.cos(C.cone + 0.02)) continue;
    if (u[0] * end[b][0] + u[1] * end[b][1] + u[2] * end[b][2] > -0.995) continue;
    pairs.push([ribos[a], ribos[b]]);
  }
  const per = sys.boundary === 'periodic' ? sys.hi[0] - sys.lo[0] : 0;
  const Lc = MRNA_LINK_NT * MRNA_NM_PER_NT;
  const lists: number[][] = Array.from({ length: sys.n }, () => []);
  for (const [i, j] of pairs) {
    lists[i].push(...C.exit, Lc, ...C.entry, MRNA_LP, j, per, i, 1);
    lists[j].push(...C.entry, Lc, ...C.exit, MRNA_LP, i, per, j, 0);
  }
  const out: number[] = [];
  for (let i = 0; i < sys.n; i++) { sys.bondOff[i] = out.length / 12; sys.bondCnt[i] = lists[i].length / 12; out.push(...lists[i]); }
  sys.bonds = new Float32Array(out);
  const inPoly = new Set<number>(); for (const [i, j] of pairs) { inPoly.add(i); inPoly.add(j); }
  return { nLinks: pairs.length, nPolyRibo: inPoly.size };
}

function qrotJs(q: ArrayLike<number>, x: number, y: number, z: number): [number, number, number] {
  const qx = q[0], qy = q[1], qz = q[2], qw = q[3];
  const tx = 2 * (qy * z - qz * y), ty = 2 * (qz * x - qx * z), tz = 2 * (qx * y - qy * x);
  return [x + qw * tx + (qy * tz - qz * ty), y + qw * ty + (qz * tx - qx * tz), z + qw * tz + (qx * ty - qy * tx)];
}
