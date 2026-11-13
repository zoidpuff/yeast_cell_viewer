// Rigid-body Brownian dynamics of a box of the cell on the GPU (WebGPU compute).
//
// Each molecule is a rigid body: a set of spheres in its own frame (schematic mode: the library sphere tree, ~30
// residues per sphere; mid mode: clusters of ~4 residues / 2 nucleotides; residue mode: one sphere per exposed residue
// bead). Per step:
//   pose:      sphere world positions from the body poses
//   bin:       spheres into two uniform grids (atomic counters, capacity from the starting occupancy): grid A for the
//              ordinary spheres (radius up to the 90th percentile), grid B, with larger cells, for the few big ones, so
//              that the common small spheres need not scan far for the rare big ones; each cell also keeps the lowest
//              and highest body id binned into it, so cells holding only the scanning sphere's own body are skipped
//   force:     per sphere, against spheres of other bodies: overlap repulsion (kRep x overlap), a short-range attraction
//              (nonspecific depth eps over ATT_W nm; native contacts of known complexes: depth epsGo), screened
//              electrostatics between charged spheres (Debye-Hueckel); membrane walls from a precomputed distance grid.
//              Cells farther than the interaction range from the sphere are skipped
//   bonds:     per body, entropic springs to linked bodies (polysome mRNA linkers: worm-like chain)
//   integrate: one workgroup per body sums force and torque over its spheres (parallel reduction), then the
//              overdamped Brownian update
//              dx = D F dt + sqrt(2 D dt) xi,  dtheta = Dr tau dt + sqrt(2 Dr dt) xi   (kT = 1)
//              with D, Dr the Stokes-Einstein(-Debye) values in water / viscosity, times the mean-field
//              hydrodynamic factors of the box's volume fraction when enabled (params.hiTrans / hiRot)
// Fixed bodies (filaments, chromatin, nuclear pores) only act as obstacles. Box faces: periodic (minimum image,
// wrapped centres), closed (reflecting centres) or embedded (closed, with the outer shell tethered to its start).
import type { Renderer } from '../renderer';
import { FILL } from '../shaders';
import { hiTrans, hiRot, debyeKappa, BJERRUM_NM, type SimSettings } from './params';

export interface SimSystem {
  n: number;                       // bodies
  type: Uint32Array; pos: Float32Array; quat: Float32Array;   // quat (x,y,z,w)
  mobile: Uint8Array; dt: Float32Array; dr: Float32Array;     // per body: mobile flag, D_t (nm^2/ns), D_r (1/ns) in water
  wall: Uint8Array;                // 1: soluble body that membranes repel
  sphOff: Uint32Array; sphCnt: Uint32Array; sphLocal: Float32Array;   // per body sphere ranges, local xyzr
  sphProp: Uint32Array;            // per sphere: charge (f32 bits, e), native-contact list offset, count, type << 16 | local id
  native: Uint32Array;             // native contacts, 2 words each: partner key (type << 16 | local id), share of the interface energy (f32 bits)
  bondOff: Uint32Array; bondCnt: Uint32Array;   // per body: its bonds (each bond is stored on both bodies)
  bonds: Float32Array;             // 12 per bond: own anchor xyz, contour length (nm); other anchor xyz, persistence length (nm); other body, periodic box length (0: none), self body, drawn (1 on one of the two copies)
  beadOff: Uint32Array; beadCnt: Uint32Array; beadLocal: Float32Array; beadKind: Uint32Array;  // drawing
  color: Uint32Array;              // per body: packed colour0, colour1 (2 words)
  lo: [number, number, number]; hi: [number, number, number];
  tether: Uint8Array;              // 1: mobile but tethered to its start (embedded boundary shell)
  mem2d: Uint8Array;               // 1: membrane protein diffusing laterally in its membrane (dt, dr: Saffman-Delbrueck)
  boundary: 'periodic' | 'walls' | 'embedded';
  phi: number;                     // macromolecular volume fraction of the box (mean-field hydrodynamics)
}

const CAP_MIN = 16;       // spheres per grid cell: at least this, else 1.5x the fullest cell of the starting configuration
                          // (overflowing spheres would be skipped as neighbours: asymmetric, self-propelling forces)
const WG = 64;
const ATT_W = 0.6;        // nm, range of the attraction beyond contact
const MAX_DRIFT = 0.15;   // nm, cap on the deterministic displacement per step (body centre, and body surface by rotation):
                          // contact forces summed over many sphere pairs are stiff; uncapped they overshoot the 0.6 nm
                          // attraction well and heat the system
export const ELEC_RANGE = 3;   // electrostatics cut off this many Debye lengths beyond contact (force shifted to 0 there)
const LINK_SEG = 40;      // drawn beads per mRNA linker

// geometry of the cell (membrane distance functions) from the fill shader, rebound for compute
const GEOM = FILL.slice(FILL.indexOf('struct Prim'), FILL.indexOf('fn tileLoaded'));

function qrotJs(q: ArrayLike<number>, x: number, y: number, z: number): [number, number, number] {
  const qx = q[0], qy = q[1], qz = q[2], qw = q[3];
  const tx = 2 * (qy * z - qz * y), ty = 2 * (qz * x - qx * z), tz = 2 * (qx * y - qy * x);
  return [x + qw * tx + (qy * tz - qz * ty), y + qw * ty + (qz * tx - qx * tz), z + qw * tz + (qx * ty - qy * tx)];
}

const SIM_COMMON = /* wgsl */ `
struct SimU {
  lo: vec4f,          // box min, w = dt (ns)
  hi: vec4f,          // box max, w = eps (kT)
  grid: vec4u,        // grid A dims xyz, cap
  cell: vec4f,        // grid A cell size, kRep, wall stiffness, max drift (nm)
  counts: vec4u,      // bodies, spheres, step, seed
  misc: vec4f,        // diffusion scale (1/viscosity), attraction width, membrane half thickness, max angle (rad)
  wgrid: vec4f,       // wall grid: origin xyz, voxel size
  wdims: vec4u,       // wall grid dims, enabled
  hydro: vec4f,       // mean-field hydrodynamic factors: translational, rotational (1: off); volume fraction; w: largest sphere radius
  bc: vec4f,          // boundary: x 0 closed walls, 1 periodic, 2 embedded (tethered shell); y tether stiffness (kT/nm^2)
  gridB: vec4u,       // grid B (spheres larger than the split radius) dims xyz, cap
  cellB: vec4f,       // grid B cell size, split radius, cores in B, largest surface sphere radius
  offs: vec4u,        // grid B's first cell counter, its first item, total cells (owner words follow the counters), surface spheres in B
  elec: vec4f,        // electrostatics: Debye kappa (1/nm), Bjerrum length (nm), cutoff beyond contact (nm), on
  spec: vec4f,        // native contacts: energy of a whole interface (kT; each contact pair gets its share), on; z: exp(-kappa x electrostatic cutoff)
};
@group(0) @binding(0) var<uniform> su: SimU;
fn hash32(x0: u32) -> u32 {
  var x = x0; x ^= x >> 16u; x *= 0x7feb352du; x ^= x >> 15u; x *= 0x846ca68bu; x ^= x >> 16u; return x;
}
fn qrot(q: vec4f, v: vec3f) -> vec3f { let t = 2.0 * cross(q.xyz, v); return v + q.w * t + cross(q.xyz, t); }
fn periodic() -> bool { return su.bc.x == 1.0; }
fn boxL() -> vec3f { return su.hi.xyz - su.lo.xyz; }
fn wrapP(p: vec3f) -> vec3f { let L = boxL(); let q = p - su.lo.xyz; return su.lo.xyz + q - L * floor(q / L); }
fn minImage(v: vec3f) -> vec3f { let L = boxL(); return v - L * round(v / L); }
fn cellOfG(p: vec3f, h: f32, G: vec3u) -> vec3i {
  // periodic: cells tile the box exactly; otherwise one padding cell on each side (positions outside are clamped)
  if (periodic()) { return min(vec3i(floor((wrapP(p) - su.lo.xyz) / h)), vec3i(G) - vec3i(1)); }
  return clamp(vec3i(floor((p - su.lo.xyz) / h)) + vec3i(1), vec3i(0), vec3i(G) - vec3i(1));
}
fn cellIdxG(c: vec3i, G: vec3u) -> u32 { return (u32(c.x) * G.y + u32(c.y)) * G.z + u32(c.z); }
`;

const SIM_CS = SIM_COMMON + /* wgsl */ `
@group(0) @binding(1) var<storage, read_write> bodyPos: array<vec4f>;    // xyz, w = D_t
@group(0) @binding(2) var<storage, read_write> bodyQuat: array<vec4f>;
@group(0) @binding(3) var<storage, read> bodyInfo: array<vec4u>;         // sphere offset, count, flags (1 mobile, 2 wall, 4 tether, 8 membrane; radius << 8), D_r bits
@group(0) @binding(4) var<storage, read> sphLocal: array<vec4f>;
@group(0) @binding(5) var<storage, read> sphBody: array<u32>;
@group(0) @binding(6) var<storage, read_write> sphWorld: array<vec4f>;
@group(0) @binding(7) var<storage, read_write> sphForce: array<vec4f>;
@group(0) @binding(8) var<storage, read_write> cellCount: array<atomic<u32>>;   // per cell (A, then B): count; then per cell: lowest, highest body
@group(0) @binding(9) var<storage, read_write> cellItems: array<vec2u>;        // sphere, body | charged << 31 | core << 30
@group(0) @binding(10) var<storage, read_write> wallGrid: array<vec4f>;   // signed distance to the nearest membrane, gradient
@group(0) @binding(11) var<storage, read_write> ctr: array<u32>;         // [0]: step counter (RNG stream)
@group(0) @binding(12) var<storage, read_write> bodyAnchor: array<vec4f>; // start: xyz position (embedded tether), w depth to the membrane (lateral diffusion)
@group(0) @binding(13) var<storage, read> sphProp: array<vec4u>;
@group(0) @binding(14) var<storage, read> native: array<vec2u>;        // partner key (type << 16 | local id), share of the interface energy (f32 bits)
@group(0) @binding(15) var<storage, read> bonds: array<vec4f>;
@group(0) @binding(16) var<storage, read_write> bodyExt: array<vec4f>;   // per body: bond force, bond torque
@group(0) @binding(17) var<storage, read> bodyBond: array<vec2u>;       // per body: first bond, count

@compute @workgroup_size(${WG})
fn clearGrid(@builtin(global_invocation_id) g: vec3u) {
  let i = g.x; if (i >= su.offs.z) { return; }
  atomicStore(&cellCount[i], 0u);
  atomicStore(&cellCount[su.offs.z + 2u * i], 0xFFFFFFFFu); atomicStore(&cellCount[su.offs.z + 2u * i + 1u], 0u);
}
@compute @workgroup_size(1)
fn tick() { ctr[0] = ctr[0] + 1u; }

@compute @workgroup_size(${WG})
fn pose(@builtin(global_invocation_id) g: vec3u) {
  let s = g.x; if (s >= su.counts.y) { return; }
  let b = sphBody[s]; let l = sphLocal[s];
  sphWorld[s] = vec4f(bodyPos[b].xyz + qrot(bodyQuat[b], l.xyz), l.w);
}

@compute @workgroup_size(${WG})
fn bin(@builtin(global_invocation_id) g: vec3u) {
  let s = g.x; if (s >= su.counts.y) { return; }
  let w = sphWorld[s]; let b = sphBody[s];
  let pr = sphProp[s];
  let tag = b | select(0u, 0x80000000u, bitcast<f32>(pr.x) != 0.0) | select(0u, 0x40000000u, pr.w == 0xFFFFFFFEu);   // bit 31 charged, bit 30 core
  var ci = 0u;
  if (w.w <= su.cellB.y) {
    let c = cellIdxG(cellOfG(w.xyz, su.cell.x, su.grid.xyz), su.grid.xyz);
    let k = atomicAdd(&cellCount[c], 1u);
    if (k < su.grid.w) { cellItems[c * su.grid.w + k] = vec2u(s, tag); }
    ci = c;
  } else {
    let c = cellIdxG(cellOfG(w.xyz, su.cellB.x, su.gridB.xyz), su.gridB.xyz);
    let k = atomicAdd(&cellCount[su.offs.x + c], 1u);
    if (k < su.gridB.w) { cellItems[su.offs.y + c * su.gridB.w + k] = vec2u(s, tag); }
    ci = su.offs.x + c;
  }
  atomicMin(&cellCount[su.offs.z + 2u * ci], b); atomicMax(&cellCount[su.offs.z + 2u * ci + 1u], b);
}

fn wallAt(p: vec3f) -> vec4f {
  let f = (p - su.wgrid.xyz) / su.wgrid.w - vec3f(0.5);
  let d = vec3i(su.wdims.xyz);
  let i0 = clamp(vec3i(floor(f)), vec3i(0), d - vec3i(2));
  let t = clamp(f - vec3f(i0), vec3f(0.0), vec3f(1.0));
  var acc = vec4f(0.0);
  for (var k = 0u; k < 8u; k++) {
    let o = vec3i(i32(k & 1u), i32((k >> 1u) & 1u), i32((k >> 2u) & 1u));
    let w = select(1.0 - t.x, t.x, o.x == 1) * select(1.0 - t.y, t.y, o.y == 1) * select(1.0 - t.z, t.z, o.z == 1);
    let c = i0 + o;
    acc += w * wallGrid[(u32(c.x) * su.wdims.y + u32(c.y)) * su.wdims.z + u32(c.z)];
  }
  return acc;
}

// force on sphere (p, r, properties me) from sphere q of another body
fn pairF(p: vec3f, r: f32, me: vec4u, q: u32, core: bool) -> vec3f {
  let o = sphWorld[q];
  var v = p - o.xyz; if (periodic()) { v = minImage(v); }
  let d2 = dot(v, v); let rr = r + o.w; let aw = su.misc.y;
  if (core) {      // cores: overlap repulsion only (they keep molecules from passing through each other)
    if (d2 >= rr * rr || d2 < 1e-12) { return vec3f(0.0); }
    let d = sqrt(d2); return su.cell.y * (rr - d) * v / d;
  }
  let myQ = bitcast<f32>(me.x);
  let charged = su.elec.w > 0.0 && myQ != 0.0;
  var lim = rr + aw;
  if (charged) { lim = max(lim, rr + su.elec.z); }
  if (d2 >= lim * lim || d2 < 1e-12) { return vec3f(0.0); }
  let d = sqrt(d2);
  var f = 0.0;
  if (d < rr) { f += su.cell.y * (rr - d); }
  var qp = vec4u(0u);
  if (charged || me.z > 0u) { qp = sphProp[q]; }
  var depth = su.hi.w;
  if (me.z > 0u && su.spec.y > 0.0) {
    for (var k = me.y; k < me.y + me.z; k++) { let e = native[k]; if (e.x == qp.w) { depth = max(depth, su.spec.x * bitcast<f32>(e.y)); break; } }
  }
  if (depth > 0.0 && d >= rr && d < rr + aw) { f -= depth * 1.5707963 / aw * sin(3.1415927 * (d - rr) / aw); }
  if (charged) {
    let oq = bitcast<f32>(qp.x);
    let rc = rr + su.elec.z;
    if (oq != 0.0 && d < rc) {
      // Debye-Hueckel between charged spheres of radii r, o.w: U = lB q1 q2 exp(-k(d - r - o.w)) / ((1 + k r)(1 + k o.w) d),
      // force shifted to zero at the cutoff; inside contact the contact value
      let k = su.elec.x; let re = max(d, rr);
      let A = su.elec.y * myQ * oq / ((1.0 + k * r) * (1.0 + k * o.w));
      f += A * (exp(k * (rr - re)) * (1.0 + k * re) / (re * re) - su.spec.z * (1.0 + k * rc) / (rc * rc));
    }
  }
  return f * v / d;
}

// all spheres of other bodies in one grid within cut of p (cells whose box is farther than cut are skipped, as are
// cells holding only body b)
fn scan(p: vec3f, r: f32, b: u32, me: vec4u, core: bool, cut: f32, cutContact: f32, h: f32, G: vec3u, cntOff: u32, itemOff: u32, cap: u32) -> vec3f {
  var f = vec3f(0.0);
  let per = periodic(); let Gi = vec3i(G);
  let c0 = cellOfG(p, h, G);
  let R = i32(ceil(cut / h));
  var c1 = c0 - vec3i(R); var nc = vec3i(2 * R + 1);
  var full = vec3<bool>(false, false, false);
  if (per) {
    if (nc.x >= Gi.x) { c1.x = 0; nc.x = Gi.x; full.x = true; }
    if (nc.y >= Gi.y) { c1.y = 0; nc.y = Gi.y; full.y = true; }
    if (nc.z >= Gi.z) { c1.z = 0; nc.z = Gi.z; full.z = true; }
  }
  // p relative to the lower corner of its own cell
  let corner = select(su.lo.xyz + (vec3f(c0) - vec3f(1.0)) * h, su.lo.xyz + vec3f(c0) * h, per);
  let pl = select(p, wrapP(p), per) - corner;
  let cut2 = cut * cut;
  for (var ix = 0; ix < nc.x; ix++) {
    for (var iy = 0; iy < nc.y; iy++) {
      for (var iz = 0; iz < nc.z; iz++) {
        var c = c1 + vec3i(ix, iy, iz);
        let o = vec3f(c - c0);
        var dd = max(vec3f(0.0), max(o * h - pl, pl - (o + 1.0) * h));
        if (per) { dd = select(dd, vec3f(0.0), full); c = ((c % Gi) + Gi) % Gi; }
        else {
          if (any(c < vec3i(0)) || any(c >= Gi)) { continue; }
          dd = select(dd, vec3f(0.0), (c == vec3i(0)) | (c == Gi - vec3i(1)));     // padding cells hold clamped outsiders
        }
        let cd2 = dot(dd, dd);
        if (cd2 >= cut2) { continue; }
        let elecOnly = cd2 >= cutContact * cutContact;     // beyond contact range: only charged spheres matter
        let ci = cntOff + cellIdxG(c, G);
        let n = min(atomicLoad(&cellCount[ci]), cap);
        if (n == 0u) { continue; }
        if (atomicLoad(&cellCount[su.offs.z + 2u * ci]) == b && atomicLoad(&cellCount[su.offs.z + 2u * ci + 1u]) == b) { continue; }
        let base = itemOff + (ci - cntOff) * cap;
        for (var k = 0u; k < n; k++) {
          let it = cellItems[base + k];
          // cores meet only cores, surface spheres only surface spheres
          if ((it.y & 0x3FFFFFFFu) == b || (elecOnly && it.y < 0x80000000u) || (((it.y & 0x40000000u) != 0u) != core)) { continue; }
          f += pairF(p, r, me, it.x, core);
        }
      }
    }
  }
  return f;
}

@compute @workgroup_size(${WG})
fn force(@builtin(global_invocation_id) g: vec3u) {
  let s = g.x; if (s >= su.counts.y) { return; }
  let b = sphBody[s];
  let info = bodyInfo[b];
  if ((info.z & 1u) == 0u) { sphForce[s] = vec4f(0.0); return; }
  let p = sphWorld[s].xyz; let r = sphWorld[s].w;
  let me = sphProp[s];
  var ext = su.misc.y;
  if (su.elec.w > 0.0 && bitcast<f32>(me.x) != 0.0) { ext = max(ext, su.elec.z); }
  // grid A holds spheres up to the split radius, grid B the larger ones
  let aw = su.misc.y;
  let core = me.w == 0xFFFFFFFEu;
  if (core) { ext = 0.0; }
  var f = scan(p, r, b, me, core, r + su.cellB.y + ext, r + su.cellB.y + aw, su.cell.x, su.grid.xyz, 0u, 0u, su.grid.w);
  // grid B: surface spheres scan it only if it holds larger surface spheres (schematic, mid), cores if it holds cores
  if ((!core && su.offs.w > 0u) || (core && su.cellB.z > 0.0)) {
    let rb = select(su.cellB.w, su.hydro.w, core);      // largest surface sphere, or largest sphere
    f += scan(p, r, b, me, core, r + rb + ext, r + rb + aw, su.cellB.x, su.gridB.xyz, su.offs.x, su.offs.y, su.gridB.w);
  }
  // membranes: soluble molecules stay HALF + r away from every bilayer, on their own side
  if ((info.z & 2u) != 0u && su.wdims.w == 1u) {
    let w = wallAt(p);
    let gap = abs(w.x) - (r + su.misc.z);
    if (gap < 0.0) {
      let n = normalize(w.yzw + vec3f(1e-6)) * sign(w.x);
      f += su.cell.z * (-gap) * n;
    }
  }
  // the box itself acts on molecule centres (integrate): wrapped (periodic) or reflected (closed, embedded)
  sphForce[s] = vec4f(f, 0.0);
}

// entropic springs between linked bodies (worm-like chain, Marko-Siggia interpolation, extension capped at 0.95 L_c)
@compute @workgroup_size(${WG})
fn bondForce(@builtin(global_invocation_id) g: vec3u) {
  let b = g.x; if (b >= su.counts.x) { return; }
  let bb = bodyBond[b];
  var F = vec3f(0.0); var T = vec3f(0.0);
  let x = bodyPos[b].xyz; let q = bodyQuat[b];
  for (var k = bb.x; k < bb.x + bb.y; k++) {
    let A = bonds[3u * k]; let B = bonds[3u * k + 1u]; let C = bonds[3u * k + 2u];
    let o = u32(C.x);
    let la = qrot(q, A.xyz);
    var dv = bodyPos[o].xyz + qrot(bodyQuat[o], B.xyz) - (x + la);
    if (periodic()) { dv = minImage(dv); }
    let L = length(dv);
    if (L < 1e-6) { continue; }
    let u = min(L / A.w, 0.95);
    let fm = (0.25 / ((1.0 - u) * (1.0 - u)) - 0.25 + u) / B.w;
    let f = fm * dv / L;
    F += f; T += cross(la, f);
  }
  bodyExt[2u * b] = vec4f(F, 0.0); bodyExt[2u * b + 1u] = vec4f(T, 0.0);
}

fn gauss(seed: u32) -> vec3f {
  let a = (f32(hash32(seed) & 0xFFFFFFu) + 0.5) / 16777216.0;
  let b = f32(hash32(seed ^ 0x68E31DA4u) & 0xFFFFFFu) / 16777216.0;
  let c = (f32(hash32(seed ^ 0xB5297A4Du) & 0xFFFFFFu) + 0.5) / 16777216.0;
  let d = f32(hash32(seed ^ 0x1B56C4E9u) & 0xFFFFFFu) / 16777216.0;
  let r1 = sqrt(-2.0 * log(a)); let r2 = sqrt(-2.0 * log(c));
  return vec3f(r1 * cos(6.2831853 * b), r1 * sin(6.2831853 * b), r2 * cos(6.2831853 * d));
}
fn qmul(a: vec4f, b: vec4f) -> vec4f { return vec4f(a.w * b.xyz + b.w * a.xyz + cross(a.xyz, b.xyz), a.w * b.w - dot(a.xyz, b.xyz)); }

// membrane proteins: remember the signed distance of their centre to the nearest membrane midplane (kept constant)
@compute @workgroup_size(${WG})
fn anchorInit(@builtin(global_invocation_id) g: vec3u) {
  let b = g.x; if (b >= su.counts.x) { return; }
  if ((bodyInfo[b].z & 8u) != 0u) { bodyAnchor[b].w = wallAt(bodyPos[b].xyz).x; }
}
fn boxFaces(x: vec3f) -> vec3f {
  if (periodic()) { return wrapP(x); }
  var xn = select(x, 2.0 * su.lo.xyz - x, x < su.lo.xyz);       // reflecting faces for the molecule centre
  xn = select(xn, 2.0 * su.hi.xyz - xn, xn > su.hi.xyz);
  return clamp(xn, su.lo.xyz, su.hi.xyz);
}

var<workgroup> redF: array<vec3f, ${WG}>;
var<workgroup> redT: array<vec3f, ${WG}>;
@compute @workgroup_size(${WG})
fn integrate(@builtin(workgroup_id) wg: vec3u, @builtin(local_invocation_index) li: u32) {
  let b = wg.x + wg.y * 65535u;
  if (b >= su.counts.x) { return; }
  let info = bodyInfo[b];
  if ((info.z & 1u) == 0u) { return; }
  // force and torque over the body's spheres: strided partial sums, then a tree reduction in workgroup memory
  let x = bodyPos[b].xyz;
  var Fp = vec3f(0.0); var Tp = vec3f(0.0);
  for (var s = info.x + li; s < info.x + info.y; s += ${WG}u) {
    let f = sphForce[s].xyz; Fp += f; Tp += cross(sphWorld[s].xyz - x, f);
  }
  redF[li] = Fp; redT[li] = Tp;
  workgroupBarrier();
  for (var st = ${WG / 2}u; st > 0u; st >>= 1u) {
    if (li < st) { redF[li] += redF[li + st]; redT[li] += redT[li + st]; }
    workgroupBarrier();
  }
  if (li != 0u) { return; }
  var F = redF[0] + bodyExt[2u * b].xyz; var T = redT[0] + bodyExt[2u * b + 1u].xyz;
  // embedded boundary: shell molecules are held near their place in the cell by a harmonic tether
  if ((info.z & 4u) != 0u) { F -= su.bc.y * (x - bodyAnchor[b].xyz); }
  let dt = su.lo.w; let ds = su.misc.x;
  if ((info.z & 8u) != 0u) {
    // membrane protein: Brownian motion in the tangent plane of its membrane (Saffman-Delbrueck D, D_r about the normal),
    // then back to its starting depth along the membrane's distance gradient
    let seed2 = hash32((b * 0x9E3779B9u) ^ hash32((ctr[0] * 0x85EBCA6Bu) ^ su.counts.w));
    let nrm = normalize(wallAt(x).yzw + vec3f(1e-6));
    let D2 = bodyPos[b].w; let Dr2 = bitcast<f32>(info.w);
    var dr2 = D2 * (F - nrm * dot(nrm, F)) * dt;
    let l2 = length(dr2); if (l2 > su.cell.w) { dr2 *= su.cell.w / l2; }
    var gn = gauss(seed2); gn -= nrm * dot(nrm, gn);
    var x2 = x + dr2 + sqrt(2.0 * D2 * dt) * gn;
    let w2 = wallAt(x2); x2 += normalize(w2.yzw + vec3f(1e-6)) * (bodyAnchor[b].w - w2.x);
    let rb2 = max(0.5, f32(info.z >> 8u) * 0.01);
    var ang = Dr2 * dot(T, nrm) * dt; ang = clamp(ang, -su.cell.w / rb2, su.cell.w / rb2);
    ang += sqrt(2.0 * Dr2 * dt) * gauss(seed2 ^ 0x27D4EB2Fu).x;
    bodyPos[b] = vec4f(boxFaces(x2), bodyPos[b].w);
    bodyQuat[b] = normalize(qmul(vec4f(nrm * sin(0.5 * ang), cos(0.5 * ang)), bodyQuat[b]));
    return;
  }
  let D = bodyPos[b].w * ds * su.hydro.x; let Dr = bitcast<f32>(info.w) * ds * su.hydro.y;
  let seed = hash32((b * 0x9E3779B9u) ^ hash32((ctr[0] * 0x85EBCA6Bu) ^ su.counts.w));
  var drift = D * F * dt;
  let dl = length(drift); if (dl > su.cell.w) { drift *= su.cell.w / dl; }
  let xn = boxFaces(x + drift + sqrt(2.0 * D * dt) * gauss(seed));
  // deterministic rotation capped so that the body surface moves at most the drift cap per step (stability of stiff
  // contact forces on large bodies, whose torques have long lever arms)
  var th = Dr * T * dt;
  let rb = max(0.5, f32(info.z >> 8u) * 0.01);
  let amax = min(su.misc.w, su.cell.w / rb);
  let tl = length(th); if (tl > amax) { th *= amax / tl; }
  th += sqrt(2.0 * Dr * dt) * gauss(seed ^ 0x27D4EB2Fu);
  let a = length(th);
  var q = bodyQuat[b];
  if (a > 1e-9) { q = normalize(qmul(vec4f(th / a * sin(0.5 * a), cos(0.5 * a)), q)); }
  bodyPos[b] = vec4f(xn, bodyPos[b].w);
  bodyQuat[b] = q;
}
`;

// wall grid: signed distance to the nearest membrane (sign: side of that membrane) and its gradient
const WALL_CS = SIM_COMMON + GEOM + /* wgsl */ `
@group(2) @binding(0) var<storage, read_write> wallOut: array<vec4f>;
fn nearest(p: vec3f) -> vec2f {
  let m = memSdfs(p);
  var best = 1e9; var k = 0u;
  for (var i = 0u; i < 9u; i++) { if (abs(m.d[i]) < abs(best)) { best = m.d[i]; k = i; } }
  return vec2f(best, f32(k));
}
fn memDist(p: vec3f, k: u32) -> f32 { let m = memSdfs(p); return m.d[k]; }
@compute @workgroup_size(${WG})
fn wall(@builtin(global_invocation_id) g: vec3u) {
  let n = su.wdims.x * su.wdims.y * su.wdims.z;
  let i = g.x; if (i >= n) { return; }
  let iz = i % su.wdims.z; let iy = (i / su.wdims.z) % su.wdims.y; let ix = i / (su.wdims.y * su.wdims.z);
  let p = su.wgrid.xyz + (vec3f(f32(ix), f32(iy), f32(iz)) + vec3f(0.5)) * su.wgrid.w;
  let nk = nearest(p); let k = u32(nk.y);
  let e = 0.5;
  let gr = vec3f(memDist(p + vec3f(e, 0.0, 0.0), k) - memDist(p - vec3f(e, 0.0, 0.0), k),
                 memDist(p + vec3f(0.0, e, 0.0), k) - memDist(p - vec3f(0.0, e, 0.0), k),
                 memDist(p + vec3f(0.0, 0.0, e), k) - memDist(p - vec3f(0.0, 0.0, e), k)) / (2.0 * e);
  wallOut[i] = vec4f(nk.x, gr);
}
`;

interface GridDef { h: number; dims: number[]; cap: number; cells: number }
// bindings of each entry point ('auto' layouts only contain the bindings an entry point uses)
const USED: Record<string, number[]> = {
  clearGrid: [0, 8], pose: [0, 1, 2, 4, 5, 6], bin: [0, 5, 6, 8, 9, 13], force: [0, 3, 5, 6, 7, 8, 9, 10, 13, 14],
  bondForce: [0, 1, 2, 15, 16, 17], integrate: [0, 1, 2, 3, 6, 7, 10, 11, 12, 16], tick: [11], anchorInit: [0, 1, 3, 10, 12] };

export class SimEngine {
  r: Renderer; d: GPUDevice; sys: SimSystem; set: SimSettings;
  uni: GPUBuffer; bufs: Record<string, GPUBuffer> = {};
  pipes: Record<string, GPUComputePipeline> = {};
  nSph = 0; rmax = 0.5; rmaxReg = 0.5; rSplit = 0.5; nLarge = 0; nLargeCore = 0; nBonds = 0;
  A!: GridDef; B!: GridDef;
  gridInfo = { h: 0, rmax: 0, cap: 0, cells: 0, reach: 0, hB: 0, capB: 0, cellsB: 0, split: 0, large: 0 };
  stepNo = 0; timeNs = 0; seed = (Math.random() * 1e9) >>> 0;
  drawPipe!: GPURenderPipeline; drawBG!: GPUBindGroup; drawG0!: GPUBindGroup; nBeads = 0;
  linkPipe: GPURenderPipeline | null = null; linkBG: GPUBindGroup | null = null;
  constructor(r: Renderer, sys: SimSystem, set: SimSettings) {
    this.r = r; this.d = r.device; this.sys = sys; this.set = set;
    this.uni = this.d.createBuffer({ size: 256, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
  }

  private buf(name: string, data: ArrayBufferView | number, usage = GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST | GPUBufferUsage.COPY_SRC) {
    const size = typeof data === 'number' ? data : data.byteLength;
    const b = this.d.createBuffer({ size: Math.max(16, Math.ceil(size / 16) * 16), usage, label: 'sim.' + name });
    if (typeof data !== 'number' && data.byteLength) this.d.queue.writeBuffer(b, 0, data as any);
    this.bufs[name]?.destroy(); this.bufs[name] = b;
    return b;
  }

  /** a uniform grid over the box for spheres of the given radii: cell = 2 r + attraction range (periodic: whole cells
   *  tiling the box, else one padding cell per side); capacity from the occupancy of the starting configuration */
  private grid(rCell: number, pick: (k: number) => boolean, world: Float32Array): GridDef {
    const s = this.sys, per = s.boundary === 'periodic', L = s.hi[0] - s.lo[0];
    const h0 = 2 * rCell + ATT_W + 0.05;
    const h = per ? L / Math.max(1, Math.floor(L / h0)) : h0;
    const dims = [0, 1, 2].map(a => per ? Math.max(1, Math.round((s.hi[a] - s.lo[a]) / h)) : Math.ceil((s.hi[a] - s.lo[a]) / h) + 2);
    const cells = dims[0] * dims[1] * dims[2];
    const cellOf = (x: number, a: number) => {
      if (per) { const w = x - s.lo[a] - L * Math.floor((x - s.lo[a]) / L); return Math.min(dims[a] - 1, Math.floor(w / h)); }
      return Math.min(dims[a] - 1, Math.max(0, Math.floor((x - s.lo[a]) / h) + 1));
    };
    const occ = new Uint32Array(cells);
    for (let k = 0; k < this.nSph; k++) {
      if (!pick(k)) continue;
      occ[(cellOf(world[3 * k], 0) * dims[1] + cellOf(world[3 * k + 1], 1)) * dims[2] + cellOf(world[3 * k + 2], 2)]++;
    }
    let omax = 0; for (const o of occ) omax = Math.max(omax, o);
    return { h, dims, cap: Math.max(CAP_MIN, Math.ceil(omax * 1.5) + 4), cells };
  }

  async init() {
    const d = this.d, s = this.sys, n = s.n;
    const pos4 = new Float32Array(4 * n), info = new Uint32Array(4 * n), infF = new Float32Array(info.buffer);
    for (let i = 0; i < n; i++) {
      pos4[4 * i] = s.pos[3 * i]; pos4[4 * i + 1] = s.pos[3 * i + 1]; pos4[4 * i + 2] = s.pos[3 * i + 2]; pos4[4 * i + 3] = s.dt[i];
      info[4 * i] = s.sphOff[i]; info[4 * i + 1] = s.sphCnt[i];
      // body radius (outermost collision sphere surface) in bits 8+ (0.01 nm): caps the rotation per step
      let rb = 0.5;
      for (let k = s.sphOff[i]; k < s.sphOff[i] + s.sphCnt[i]; k++)
        rb = Math.max(rb, Math.hypot(s.sphLocal[4 * k], s.sphLocal[4 * k + 1], s.sphLocal[4 * k + 2]) + s.sphLocal[4 * k + 3]);
      info[4 * i + 2] = (s.mobile[i] ? 1 : 0) | (s.wall[i] ? 2 : 0) | (s.tether[i] ? 4 : 0) | (s.mem2d[i] ? 8 : 0) | (Math.min(0xFFFFFF, Math.round(rb * 100)) << 8); infF[4 * i + 3] = s.dr[i];
    }
    this.nSph = s.sphLocal.length / 4;
    const sphBody = new Uint32Array(this.nSph);
    for (let i = 0; i < n; i++) sphBody.fill(i, s.sphOff[i], s.sphOff[i] + s.sphCnt[i]);
    this.buf('bodyPos', pos4); this.buf('bodyQuat', s.quat); this.buf('bodyInfo', info); this.buf('bodyAnchor', pos4.slice());
    this.buf('sphLocal', s.sphLocal); this.buf('sphBody', sphBody);
    this.buf('sphWorld', 16 * this.nSph); this.buf('sphForce', 16 * this.nSph);
    this.buf('sphProp', s.sphProp); this.buf('native', s.native.length ? s.native : new Uint32Array(4));
    // bonds (each stored on both bodies): per body first bond and count
    this.nBonds = s.bonds.length / 12;
    const bb = new Uint32Array(2 * Math.max(1, n));
    for (let i = 0; i < n; i++) { bb[2 * i] = s.bondOff[i]; bb[2 * i + 1] = s.bondCnt[i]; }
    this.buf('bonds', this.nBonds ? s.bonds : new Float32Array(12)); this.buf('bodyBond', bb); this.buf('bodyExt', 32 * Math.max(1, n));
    // two neighbour grids: A for spheres up to the 90th-percentile radius (most of them), B for the larger ones
    const world = new Float32Array(3 * this.nSph), radii = new Float32Array(this.nSph);
    for (let i = 0; i < n; i++) {
      const q = s.quat.subarray(4 * i, 4 * i + 4);
      for (let k = s.sphOff[i]; k < s.sphOff[i] + s.sphCnt[i]; k++) {
        const v = qrotJs(q, s.sphLocal[4 * k], s.sphLocal[4 * k + 1], s.sphLocal[4 * k + 2]);
        world[3 * k] = s.pos[3 * i] + v[0]; world[3 * k + 1] = s.pos[3 * i + 1] + v[1]; world[3 * k + 2] = s.pos[3 * i + 2] + v[2];
        radii[k] = s.sphLocal[4 * k + 3];
      }
    }
    const sorted = radii.slice().sort();
    this.rmax = sorted[this.nSph - 1] ?? 0.5;
    this.rSplit = sorted[Math.floor(0.9 * (this.nSph - 1))] ?? 0.5;
    if (this.rmax <= 1.05 * this.rSplit) this.rSplit = this.rmax;        // nearly uniform: one grid
    const split = this.rSplit;
    this.A = this.grid(split, k => radii[k] <= split, world);
    this.nLarge = 0; this.nLargeCore = 0; this.rmaxReg = 0;
    for (let k = 0; k < this.nSph; k++) if (s.sphProp[4 * k + 3] !== 0xFFFFFFFE) this.rmaxReg = Math.max(this.rmaxReg, radii[k]);
    for (let k = 0; k < this.nSph; k++) if (radii[k] > split) { if (s.sphProp[4 * k + 3] === 0xFFFFFFFE) this.nLargeCore++; else this.nLarge++; }
    this.B = this.nLarge + this.nLargeCore ? this.grid(this.rmax, k => radii[k] > split, world) : { h: 1, dims: [1, 1, 1], cap: 1, cells: 1 };
    const cellsAll = this.A.cells + this.B.cells;
    this.gridInfo = { h: this.A.h, rmax: this.rmax, cap: this.A.cap, cells: this.A.cells, reach: Math.ceil((2 * split + ATT_W) / this.A.h),
      hB: this.B.h, capB: this.B.cap, cellsB: this.nLarge + this.nLargeCore ? this.B.cells : 0, split, large: this.nLarge + this.nLargeCore };
    this.buf('cellCount', 4 * 3 * cellsAll); this.buf('cellItems', 8 * (this.A.cells * this.A.cap + this.B.cells * this.B.cap));
    // wall grid over the box
    const box = Math.max(...[0, 1, 2].map(a => s.hi[a] - s.lo[a]));
    const vox = Math.max(1, box / 96);
    const wd = [0, 1, 2].map(a => Math.max(2, Math.ceil((s.hi[a] - s.lo[a]) / vox)));
    this.buf('wallGrid', 16 * wd[0] * wd[1] * wd[2]);
    this.wall = { vox, wd };
    this.writeUniforms();
    // pipelines
    const mod = async (code: string, label: string) => {
      const m = d.createShaderModule({ code, label });
      const info = await m.getCompilationInfo();
      const errs = info.messages.filter(x => x.type === 'error');
      for (const e of errs) console.error(`[${label}] ${e.lineNum}:${e.linePos} ${e.message}`);
      if (errs.length) throw new Error(`WGSL error in ${label}: ${errs[0].message}`);
      return m;
    };
    const mSim = await mod(SIM_CS, 'sim');
    this.buf('ctr', new Uint32Array([this.stepNo, 0, 0, 0]));
    for (const ep of Object.keys(USED)) this.pipes[ep] = await d.createComputePipelineAsync({ layout: 'auto', compute: { module: mSim, entryPoint: ep } });
    this.makeBindGroups();
    await this.buildWalls(wd);
    {   // depth of each membrane protein below its membrane's midplane
      const enc = d.createCommandEncoder(), cp = enc.beginComputePass();
      cp.setPipeline(this.pipes.anchorInit); cp.setBindGroup(0, this.bgFor.anchorInit); cp.dispatchWorkgroups(Math.ceil(n / WG)); cp.end();
      d.queue.submit([enc.finish()]);
    }
    await this.initDraw();
  }
  bgFor!: Record<string, GPUBindGroup>;
  wall!: { vox: number; wd: number[] };
  private makeBindGroups() {
    const B = this.bufs;
    const all = [this.uni, B.bodyPos, B.bodyQuat, B.bodyInfo, B.sphLocal, B.sphBody, B.sphWorld, B.sphForce, B.cellCount, B.cellItems,
      B.wallGrid, B.ctr, B.bodyAnchor, B.sphProp, B.native, B.bonds, B.bodyExt, B.bodyBond];
    this.bgFor = {} as any;
    for (const ep of Object.keys(USED))
      this.bgFor[ep] = this.d.createBindGroup({ layout: this.pipes[ep].getBindGroupLayout(0), entries: USED[ep].map(i => ({ binding: i, resource: { buffer: all[i] } })) });
  }
  /** more room per grid cell (call between submissions): molecules that cluster tightly can fill cells beyond the
   *  capacity set from the starting configuration, and overflowing spheres would be missed as neighbours */
  growGrid(f = 2) {
    this.A.cap = Math.min(4096, Math.ceil(this.A.cap * f)); this.B.cap = Math.min(4096, Math.ceil(this.B.cap * f));
    this.buf('cellItems', 8 * (this.A.cells * this.A.cap + this.B.cells * this.B.cap));
    this.gridInfo.cap = this.A.cap; this.gridInfo.capB = this.B.cap; this.grows++;
    this.writeUniforms(); this.makeBindGroups();
  }
  grows = 0;
  /** wall-grid layout (export of the membrane surface) */
  get uniState() { return { vox: this.wall.vox, wd: this.wall.wd }; }

  private writeUniforms() {
    const s = this.sys, set = this.set, A = this.A, B = this.B, wd = this.wall.wd;
    const ab = new ArrayBuffer(256), f = new Float32Array(ab), u = new Uint32Array(ab);
    f[0] = s.lo[0]; f[1] = s.lo[1]; f[2] = s.lo[2]; f[3] = set.dt;
    f[4] = s.hi[0]; f[5] = s.hi[1]; f[6] = s.hi[2]; f[7] = set.eps;
    u[8] = A.dims[0]; u[9] = A.dims[1]; u[10] = A.dims[2]; u[11] = A.cap;
    f[12] = A.h; f[13] = set.kRep; f[14] = 4 * set.kRep; f[15] = MAX_DRIFT;
    u[16] = s.n; u[17] = this.nSph; u[18] = this.stepNo; u[19] = this.seed;
    f[20] = 1 / set.viscosity; f[21] = ATT_W; f[22] = 2.0; f[23] = 0.15;
    f[24] = s.lo[0]; f[25] = s.lo[1]; f[26] = s.lo[2]; f[27] = this.wall.vox;
    u[28] = wd[0]; u[29] = wd[1]; u[30] = wd[2]; u[31] = 1;
    f[32] = set.hydro ? hiTrans(s.phi) : 1; f[33] = set.hydro ? hiRot(s.phi) : 1; f[34] = s.phi; f[35] = this.rmax;
    f[36] = s.boundary === 'periodic' ? 1 : s.boundary === 'embedded' ? 2 : 0; f[37] = set.tether;
    u[40] = B.dims[0]; u[41] = B.dims[1]; u[42] = B.dims[2]; u[43] = B.cap;
    f[44] = B.h; f[45] = this.rSplit; f[46] = this.nLargeCore; f[47] = this.rmaxReg;
    u[48] = A.cells; u[49] = A.cells * A.cap; u[50] = A.cells + B.cells; u[51] = this.nLarge;
    const kappa = debyeKappa(set.ionic);
    f[52] = kappa; f[53] = BJERRUM_NM; f[54] = ELEC_RANGE / kappa; f[55] = set.elec ? 1 : 0;
    f[56] = set.epsGo; f[57] = set.specific ? 1 : 0; f[58] = Math.exp(-ELEC_RANGE);
    this.d.queue.writeBuffer(this.uni, 0, ab);
  }

  setSettings(set: SimSettings) { this.set = set; this.writeUniforms(); }

  private async buildWalls(wd: number[]) {
    const d = this.d, r = this.r;
    const m = d.createShaderModule({ code: WALL_CS, label: 'simWall' });
    const info = await m.getCompilationInfo();
    for (const e of info.messages.filter(x => x.type === 'error')) console.error(`[simWall] ${e.lineNum}:${e.linePos} ${e.message}`);
    const p = await d.createComputePipelineAsync({ layout: 'auto', compute: { module: m, entryPoint: 'wall' } });
    const fb = r.fillBufs;
    const geo = [fb[0], fb[1], fb[2], fb[3], fb[4], (r as any).tileState as GPUBuffer, fb[5]];
    const g0 = d.createBindGroup({ layout: p.getBindGroupLayout(0), entries: [{ binding: 0, resource: { buffer: this.uni } }] });
    const lay1 = p.getBindGroupLayout(1);
    const used1 = [0, 1, 2, 3, 4, 6];      // memSdfs needs prims, pores, grid offsets/indices, geo uniform, fields
    const g1 = d.createBindGroup({ layout: lay1, entries: used1.map(i => ({ binding: i, resource: { buffer: geo[i] } })) });
    const g2 = d.createBindGroup({ layout: p.getBindGroupLayout(2), entries: [{ binding: 0, resource: { buffer: this.bufs.wallGrid } }] });
    const enc = d.createCommandEncoder(); const cp = enc.beginComputePass();
    cp.setPipeline(p); cp.setBindGroup(0, g0); cp.setBindGroup(1, g1); cp.setBindGroup(2, g2);
    cp.dispatchWorkgroups(Math.ceil(wd[0] * wd[1] * wd[2] / WG)); cp.end();
    d.queue.submit([enc.finish()]);
    await d.queue.onSubmittedWorkDone();
  }

  /** the dispatches of one step: [entry point, workgroups x, y] in order */
  plan(): [string, number, number][] {
    const gS = Math.ceil(this.nSph / WG), gB = Math.ceil(this.sys.n / WG), gC = Math.ceil((this.A.cells + this.B.cells) / WG);
    const n = this.sys.n, gx = Math.min(n, 65535), gy = Math.ceil(n / 65535);
    const p: [string, number, number][] = [['clearGrid', gC, 1], ['pose', gS, 1], ['bin', gS, 1], ['force', gS, 1]];
    if (this.nBonds) p.push(['bondForce', gB, 1]);
    p.push(['integrate', gx, gy], ['tick', 1, 1]);
    return p;
  }

  /** n Brownian steps into the encoder, all in one compute pass (dispatches in a pass are ordered and see each other's
   *  writes): clear grids, pose, bin, force, bonds, integrate, then advance the step counter (RNG stream) on the GPU */
  step(enc: GPUCommandEncoder, n: number) {
    const P = this.pipes, G = this.bgFor, plan = this.plan();
    const cp = enc.beginComputePass();
    for (let k = 0; k < n; k++) {
      for (const [name, gx, gy] of plan) { cp.setPipeline(P[name]); cp.setBindGroup(0, G[name]); cp.dispatchWorkgroups(gx, gy); }
    }
    cp.end();
    this.stepNo += n;
    this.timeNs += n * this.set.dt;
  }

  // ---------------------------------------------------------------- drawing: residue beads of every body, mRNA linkers
  private async initDraw() {
    const d = this.d, r = this.r, s = this.sys;
    this.nBeads = s.beadLocal.length / 4;
    const beadBody = new Uint32Array(this.nBeads);
    for (let i = 0; i < s.n; i++) beadBody.fill(i, s.beadOff[i], s.beadOff[i] + s.beadCnt[i]);
    const loc = new Float32Array(s.beadLocal), lu = new Uint32Array(loc.buffer);
    for (let i = 0; i < this.nBeads; i++) lu[4 * i + 3] = s.beadKind[i];
    this.buf('beadLocal', loc); this.buf('beadBody', beadBody);
    const bodyCol = new Uint32Array(4 * s.n);
    for (let i = 0; i < s.n; i++) { bodyCol[4 * i] = s.color[2 * i]; bodyCol[4 * i + 1] = s.color[2 * i + 1]; bodyCol[4 * i + 2] = s.mobile[i]; bodyCol[4 * i + 3] = s.type[i]; }
    this.buf('bodyCol', bodyCol);
    const m = d.createShaderModule({ code: SIM_DRAW, label: 'simDraw' });
    const info = await m.getCompilationInfo();
    for (const e of info.messages.filter(x => x.type === 'error')) console.error(`[simDraw] ${e.lineNum}:${e.linePos} ${e.message}`);
    const ro = (b: number): GPUBindGroupLayoutEntry => ({ binding: b, visibility: GPUShaderStage.VERTEX | (b === 4 ? GPUShaderStage.FRAGMENT : 0), buffer: { type: 'read-only-storage' } });
    const l1 = d.createBindGroupLayout({ entries: [ro(0), ro(1), ro(2), ro(3), ro(4)] });
    const g0 = slimG0(r);
    this.drawG0 = g0.bg;
    const desc = (entry: string): GPURenderPipelineDescriptor => ({
      layout: d.createPipelineLayout({ bindGroupLayouts: [g0.layout, l1] }),
      vertex: { module: m, entryPoint: entry },
      fragment: { module: m, entryPoint: entry === 'vsBead' ? 'fsBead' : 'fsMain', targets: [{ format: 'rgba8unorm' }, { format: 'r32uint' }] },
      primitive: { topology: 'triangle-list', cullMode: 'none' },
      depthStencil: { format: 'depth32float', depthWriteEnabled: true, depthCompare: 'greater' },
    });
    this.drawPipe = await d.createRenderPipelineAsync(desc('vsBead'));
    const B = this.bufs;
    this.drawBG = d.createBindGroup({ layout: l1, entries: [B.bodyPos, B.bodyQuat, B.beadLocal, B.beadBody, B.bodyCol].map((b, i) => ({ binding: i, resource: { buffer: b } })) });
    if (this.nBonds) {
      this.linkPipe = await d.createRenderPipelineAsync(desc('vsLink'));
      this.linkBG = d.createBindGroup({ layout: l1, entries: [B.bodyPos, B.bodyQuat, B.bonds, B.beadBody, B.bodyCol].map((b, i) => ({ binding: i, resource: { buffer: b } })) });
    }
  }

  /** tracked molecules: drawn as usual, every other one faded (none tracked: all normal) */
  setTracked(ids: number[]) {
    const s = this.sys, set = new Set(ids);
    const col = new Uint32Array(4 * s.n);
    for (let i = 0; i < s.n; i++) {
      col[4 * i] = s.color[2 * i]; col[4 * i + 1] = s.color[2 * i + 1]; col[4 * i + 3] = s.type[i];
      col[4 * i + 2] = (s.mobile[i] ? 1 : 0) | (set.size ? (set.has(i) ? 4 : 2) : 0);
    }
    this.d.queue.writeBuffer(this.bufs.bodyCol, 0, col);
  }

  draw(rp: GPURenderPassEncoder) {
    rp.setPipeline(this.drawPipe); rp.setBindGroup(0, this.drawG0); rp.setBindGroup(1, this.drawBG);
    rp.draw(6, this.nBeads);
    if (this.linkPipe && this.linkBG) { rp.setPipeline(this.linkPipe); rp.setBindGroup(1, this.linkBG); rp.draw(6, this.nBonds * LINK_SEG); }
  }

  /** fullest grid cells of the last step and the share of spheres that did not fit (should be 0) */
  async gridLoad(): Promise<{ max: number; cap: number; dropped: number }> {
    const nc = this.A.cells + this.B.cells;
    const d = this.d, rb = d.createBuffer({ size: 4 * nc, usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST });
    const enc = d.createCommandEncoder(); enc.copyBufferToBuffer(this.bufs.cellCount, 0, rb, 0, 4 * nc); d.queue.submit([enc.finish()]);
    await rb.mapAsync(GPUMapMode.READ); const c = new Uint32Array(rb.getMappedRange().slice(0)); rb.unmap(); rb.destroy();
    let max = 0, drop = 0;
    for (let i = 0; i < nc; i++) {
      const cap = i < this.A.cells ? this.A.cap : this.B.cap, v = c[i];
      max = Math.max(max, v / cap * this.A.cap); if (v > cap) drop += v - cap;
    }
    return { max: Math.round(max), cap: this.A.cap, dropped: drop / Math.max(1, this.nSph) };
  }

  /** queue a copy of the body poses (position+D, quaternion: 32 bytes per body) into dst */
  copyPoses(enc: GPUCommandEncoder, dst: GPUBuffer) {
    const n = this.sys.n;
    enc.copyBufferToBuffer(this.bufs.bodyPos, 0, dst, 0, 16 * n); enc.copyBufferToBuffer(this.bufs.bodyQuat, 0, dst, 16 * n, 16 * n);
  }
  unpackPoses(a: Float32Array): { pos: Float32Array; quat: Float32Array } {
    const n = this.sys.n, pos = new Float32Array(3 * n);
    for (let i = 0; i < n; i++) { pos[3 * i] = a[4 * i]; pos[3 * i + 1] = a[4 * i + 1]; pos[3 * i + 2] = a[4 * i + 2]; }
    return { pos, quat: a.slice(4 * n, 8 * n) };
  }
  /** current body poses (for recording and export) */
  async readPoses(): Promise<{ pos: Float32Array; quat: Float32Array }> {
    const d = this.d, rb = d.createBuffer({ size: 32 * this.sys.n, usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST });
    const enc = d.createCommandEncoder(); this.copyPoses(enc, rb); d.queue.submit([enc.finish()]);
    await rb.mapAsync(GPUMapMode.READ);
    const a = new Float32Array(rb.getMappedRange().slice(0)); rb.unmap(); rb.destroy();
    return this.unpackPoses(a);
  }

  /** upload poses (playback of a recorded frame) */
  writePoses(pos: Float32Array, quat: Float32Array) {
    const n = this.sys.n, p4 = new Float32Array(4 * n);
    for (let i = 0; i < n; i++) { p4[4 * i] = pos[3 * i]; p4[4 * i + 1] = pos[3 * i + 1]; p4[4 * i + 2] = pos[3 * i + 2]; p4[4 * i + 3] = this.sys.dt[i]; }
    this.d.queue.writeBuffer(this.bufs.bodyPos, 0, p4); this.d.queue.writeBuffer(this.bufs.bodyQuat, 0, quat as any);
  }

  destroy() { for (const b of Object.values(this.bufs)) b.destroy(); this.uni.destroy(); }
}

/** group 0 for the simulation's draws: the frame uniforms and the palette only (the renderer's group 0 also binds
 *  its instance pools, which would exceed the vertex-stage storage-buffer limit together with the simulation's) */
export function slimG0(r: Renderer): { layout: GPUBindGroupLayout; bg: GPUBindGroup } {
  const d = r.device, VF = GPUShaderStage.VERTEX | GPUShaderStage.FRAGMENT;
  const layout = d.createBindGroupLayout({ entries: [
    { binding: 0, visibility: VF, buffer: { type: 'uniform' } }, { binding: 5, visibility: VF, buffer: { type: 'read-only-storage' } }] });
  const bg = d.createBindGroup({ layout, entries: [{ binding: 0, resource: { buffer: (r as any).uni } }, { binding: 5, resource: { buffer: (r as any).palette } }] });
  return { layout, bg };
}

const SIM_DRAW = FILL.slice(0, FILL.indexOf('@group(0) @binding(5)')) + /* wgsl */ `
@group(0) @binding(5) var<storage, read> palette: array<vec4f>;
` + (() => {
  // impostor helpers from the sphere shader (VOut/FOut, emitQuad, shadeFrag)
  const i0 = DRAW_SPHERES_SRC.indexOf('struct VOut'), i1 = DRAW_SPHERES_SRC.indexOf('@group(1) @binding(0)');
  return DRAW_SPHERES_SRC.slice(i0, i1);
})() + /* wgsl */ `
@group(1) @binding(0) var<storage, read> bodyPos: array<vec4f>;
@group(1) @binding(1) var<storage, read> bodyQuat: array<vec4f>;
@group(1) @binding(2) var<storage, read> beadLocal: array<vec4f>;    // xyz, w = bits(kind | plddt << 8)
@group(1) @binding(3) var<storage, read> beadBody: array<u32>;
@group(1) @binding(4) var<storage, read> bodyCol: array<vec4u>;      // colour0, colour1, flags (1 mobile, 2 faded while tracking others, 4 tracked), type
@vertex
fn vsBead(@builtin(vertex_index) vi: u32, @builtin(instance_index) ii: u32) -> VOut {
  let b = beadBody[ii]; let l = beadLocal[ii];
  let kind = bitcast<u32>(l.w) & 255u;
  let p = bodyPos[b].xyz + qrot(bodyQuat[b], l.xyz);
  let bc = bodyCol[b];
  var c = unpackColor(bc.x);
  if (kind == 1u && bc.y != 0u) { c = unpackColor(bc.y); }
  if ((bc.z & 1u) == 0u) { c = mix(c, vec3f(0.80, 0.80, 0.78), 0.55); }    // fixed bodies greyed
  var r = u.lod.w; if (kind == 1u) { r = u.misc.y; }
  return emitQuad(p - u.eye.xyz, r, vi, c, ID_INST | b);
}
// mRNA linkers: beads along the line between the two anchors (minimum image in periodic boxes); each bond is stored on
// both bodies and drawn from the copy flagged in its w
@vertex
fn vsLink(@builtin(vertex_index) vi: u32, @builtin(instance_index) ii: u32) -> VOut {
  let k = ii / ${LINK_SEG}u; let t = (f32(ii % ${LINK_SEG}u) + 0.5) / ${LINK_SEG}.0;
  let A = beadLocal[3u * k]; let B = beadLocal[3u * k + 1u]; let C = beadLocal[3u * k + 2u];
  let me = u32(C.z); let o = u32(C.x);
  let pa = bodyPos[me].xyz + qrot(bodyQuat[me], A.xyz);
  var dv = bodyPos[o].xyz + qrot(bodyQuat[o], B.xyz) - pa;
  let L = C.y;                // periodic box length (0: not periodic)
  if (L > 0.0) { dv -= L * round(dv / L); }
  var r = 0.32; if (C.w == 0.0) { r = 0.0; }
  return emitQuad(pa + t * dv - u.eye.xyz, r, vi, vec3f(0.36, 0.22, 0.14), 0u);     // dark brown: apart from the orange box outline and the path colours
}
@fragment
fn fsMain(in: VOut) -> FOut { return shadeFrag(in); }
// molecules faded while others are tracked: screen-door transparency (one pixel in four drawn), no sorting needed
@fragment
fn fsBead(in: VOut) -> FOut {
  if ((bodyCol[in.id & 0x1FFFFFFFu].z & 2u) != 0u) {
    let q = vec2u(in.pos.xy);
    if (((q.x & 1u) | (q.y & 1u)) != 0u) { discard; }
    var o = shadeFrag(in);
    o.color = vec4f(mix(o.color.xyz, vec3f(0.965, 0.945, 0.906), 0.3), 1.0);
    o.id = in.id | 0x20000000u;          // ghost: still pickable, no outlines (composite)
    return o;
  }
  return shadeFrag(in);
}
`;
import { DRAW_SPHERES as DRAW_SPHERES_SRC } from '../shaders';
