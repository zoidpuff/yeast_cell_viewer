(function(){const e=document.createElement("link").relList;if(e&&e.supports&&e.supports("modulepreload"))return;for(const n of document.querySelectorAll('link[rel="modulepreload"]'))i(n);new MutationObserver(n=>{for(const r of n)if(r.type==="childList")for(const o of r.addedNodes)o.tagName==="LINK"&&o.rel==="modulepreload"&&i(o)}).observe(document,{childList:!0,subtree:!0});function t(n){const r={};return n.integrity&&(r.integrity=n.integrity),n.referrerPolicy&&(r.referrerPolicy=n.referrerPolicy),n.crossOrigin==="use-credentials"?r.credentials="include":n.crossOrigin==="anonymous"?r.credentials="omit":r.credentials="same-origin",r}function i(n){if(n.ep)return;n.ep=!0;const r=t(n);fetch(n.href,r)}})();const H=1024,Y=4096,ye=`
const PAGE: u32 = ${H}u;
const SPAGE: u32 = ${Y}u;
const INV_SQRT2: f32 = 0.70710678;
struct U {
  proj: mat4x4f,
  view: mat4x4f,          // rotation only (camera-relative)
  invViewProj: mat4x4f,   // camera-relative clip -> world-relative
  eye: vec4f,             // world eye, w = px per nm at unit view depth
  plane: vec4f,           // world n.xyz, d
  planeV: vec4f,          // view-space n.xyz, sEye = dot(n,eye)-d
  cut: vec4f,             // x enabled, y slab thickness, z fog range, w time
  lod: vec4f,             // x L0 bead px, y L1 instance px, z min px for free spheres/lipids, w protein bead radius
  vp: vec4f,              // width, height, 1/w, 1/h
  frustum: array<vec4f, 6>,  // camera-relative planes, inside: dot(n,p)+w >= 0
  modes: vec4u,           // x colorMode, y regionMask, z unused, w highlightActive
  sel: vec4u,             // x selected id, y selected type+1, z shaded(1), w organelle mode
  misc: vec4f,            // x QSTEP, y nucleic bead radius, z outline strength, w dpr
  cell: vec4f,            // x R outer, y R mannan inner, z R PM, w near
  caps: vec4u,            // x capL0 (offset of L1 items in list01)
};
@group(0) @binding(0) var<uniform> u: U;

fn sPlaneW(p: vec3f) -> f32 { return dot(u.plane.xyz, p + u.eye.xyz) - u.plane.w; }   // p camera-relative

fn hash32(x0: u32) -> u32 {
  var x = x0; x ^= x >> 16u; x *= 0x7feb352du; x ^= x >> 15u; x *= 0x846ca68bu; x ^= x >> 16u; return x;
}
fn unpackColor(c: u32) -> vec3f {
  return vec3f(f32(c & 255u), f32((c >> 8u) & 255u), f32((c >> 16u) & 255u)) / 255.0;
}
fn qrot(q: vec4f, v: vec3f) -> vec3f {   // q = (x,y,z,w)
  let t = 2.0 * cross(q.xyz, v);
  return v + q.w * t + cross(q.xyz, t);
}
fn sphereInFrustum(c: vec3f, r: f32) -> bool {
  for (var k = 0; k < 6; k++) { if (dot(u.frustum[k].xyz, c) + u.frustum[k].w < -r) { return false; } }
  return true;
}
fn i16f(x: u32) -> f32 { return f32(bitcast<i32>(x << 16u) >> 16u) / 32767.0 * INV_SQRT2; }
`,vt=`
struct TypeInfo {
  beadOff: u32, beadCount: u32, clusterOff: u32, clusterCount: u32,
  stOff: u32, stCount: u32, br: f32, req: f32,
  color0: u32, color1: u32, flags: u32, plddt: f32,
};
struct Inst { t: u32, region: u32, flags: u32, p: vec3f, q: vec4f };
@group(0) @binding(1) var<storage, read> inst: array<vec4u>;
@group(0) @binding(2) var<storage, read> pageOrg: array<vec4f>;       // xyz origin (camera-relative is computed), w = bitcast count
@group(0) @binding(3) var<storage, read> types: array<TypeInfo>;
@group(0) @binding(4) var<storage, read> typeFlags: array<u32>;       // bit0 highlight, bit1 hidden

fn decodeInst(i: u32) -> Inst {
  let w = inst[i];
  let org = pageOrg[i / PAGE].xyz;
  var o: Inst;
  o.t = w.x & 0x3FFFu;
  let qi = (w.x >> 14u) & 3u;
  o.region = (w.x >> 16u) & 31u;
  o.flags = w.x >> 21u;
  o.p = org + vec3f(f32(w.y & 0xFFFFu), f32(w.y >> 16u), f32(w.z & 0xFFFFu)) * u.misc.x - u.eye.xyz;
  let a = i16f(w.z >> 16u); let b = i16f(w.w & 0xFFFFu); let c = i16f(w.w >> 16u);
  let big = sqrt(max(0.0, 1.0 - a * a - b * b - c * c));
  // stored order (w,x,y,z) with component qi replaced by 'big'
  var q4: vec4f;
  if (qi == 0u) { q4 = vec4f(big, a, b, c); }
  else if (qi == 1u) { q4 = vec4f(a, big, b, c); }
  else if (qi == 2u) { q4 = vec4f(a, b, big, c); }
  else { q4 = vec4f(a, b, c, big); }
  o.q = vec4f(q4.y, q4.z, q4.w, q4.x);    // -> (x,y,z,w)
  return o;
}
`,Gt=ye+vt+`
struct Counters { l0: atomic<u32>, l1: atomic<u32>, l2: atomic<u32>, lip: atomic<u32>, fre: atomic<u32>, ninst: atomic<u32>, p0: u32, p1: u32 };
struct CullParams { numPages: u32, strideX: u32, capL0: u32, capL1: u32, capL2: u32, capLip: u32, capFree: u32, numSPages: u32 };
struct Cluster { c: vec4f, info: vec4u };   // centre xyz, radius; bead offset, count
@group(1) @binding(0) var<storage, read_write> counters: Counters;
@group(1) @binding(1) var<storage, read_write> list01: array<vec2u>;   // [0,capL0) L0 items, [capL0, capL0+capL1) L1 items
@group(1) @binding(2) var<storage, read_write> listL2: array<u32>;
@group(1) @binding(3) var<uniform> cp: CullParams;
@group(1) @binding(4) var<storage, read> clusters: array<Cluster>;

fn visibleClass(t: TypeInfo, region: u32, tf: u32) -> bool {
  if ((tf & 2u) != 0u) { return false; }
  if (((u.modes.y >> region) & 1u) == 0u) { return false; }
  return true;
}

@compute @workgroup_size(64)
fn cullInstances(@builtin(global_invocation_id) gid: vec3u) {
  let i = gid.x + gid.y * cp.strideX;
  let page = i / PAGE;
  if (page >= cp.numPages) { return; }
  let cnt = bitcast<u32>(pageOrg[page].w);
  if ((i % PAGE) >= cnt) { return; }
  let I = decodeInst(i);
  let T = types[I.t];
  if (!visibleClass(T, I.region, typeFlags[I.t])) { return; }
  let br = T.br;
  if (u.cut.x > 0.5) {
    let s = sPlaneW(I.p);
    if (s < -br || s - br > u.cut.y) { return; }
  }
  if (!sphereInFrustum(I.p, br)) { return; }
  atomicAdd(&counters.ninst, 1u);
  let dist = max(length(I.p), 0.01);
  let pxPerNm = u.eye.w / dist;
  let beadPx = 2.0 * u.lod.w * pxPerNm;
  if (beadPx >= u.lod.x && T.beadCount > 0u) {
    for (var c = 0u; c < T.clusterCount; c++) {
      let CL = clusters[T.clusterOff + c];
      let cs = CL.c;
      let cc = I.p + qrot(I.q, cs.xyz);
      let r = cs.w;
      if (u.cut.x > 0.5) {
        let s = sPlaneW(cc);
        if (s < -r || s - r > u.cut.y) { continue; }
      }
      if (!sphereInFrustum(cc, r)) { continue; }
      let ci = CL.info;
      let k = atomicAdd(&counters.l0, 1u);
      if (k < cp.capL0) { list01[k] = vec2u(i, ci.x | ((ci.y - 1u) << 22u)); }
    }
  } else if (br * pxPerNm >= u.lod.y && T.stCount > 0u) {
    let n = T.stCount;
    var o = 0u;
    loop {
      if (o >= n) { break; }
      let m = min(16u, n - o);
      let k = atomicAdd(&counters.l1, 1u);
      if (k < cp.capL1) { list01[cp.capL0 + k] = vec2u(i, (T.stOff + o) | ((m - 1u) << 22u)); }
      o += 16u;
    }
  } else {
    let k = atomicAdd(&counters.l2, 1u);
    if (k < cp.capL2) { listL2[k] = i; }
  }
}
`,Nt=ye+`
struct Counters { l0: atomic<u32>, l1: atomic<u32>, l2: atomic<u32>, lip: atomic<u32>, fre: atomic<u32>, ninst: atomic<u32>, p0: u32, p1: u32 };
struct SParams { numPages: u32, strideX: u32, cap: u32, isLipid: u32 };
@group(1) @binding(0) var<storage, read_write> counters: Counters;
@group(1) @binding(1) var<storage, read> spool: array<vec4f>;
@group(1) @binding(2) var<storage, read> spages: array<vec4u>;      // x count, y tile slot, z visible flag
@group(1) @binding(3) var<storage, read_write> slist: array<u32>;
@group(1) @binding(4) var<uniform> sp: SParams;

fn freeRadius(kind: u32) -> f32 {
  if (kind == 0u) { return 0.6; }   // mRNA nucleotide
  if (kind == 1u) { return 1.0; }   // linker DNA (2 bp per bead)
  if (kind == 2u) { return 0.75; }  // beta-glucan triple helix
  return 0.6;                        // mannan
}

@compute @workgroup_size(64)
fn cullSpheres(@builtin(global_invocation_id) gid: vec3u) {
  let i = gid.x + gid.y * sp.strideX;
  let page = i / SPAGE;
  if (page >= sp.numPages) { return; }
  let pg = spages[page];
  if (pg.z == 0u || (i % SPAGE) >= pg.x) { return; }
  let rec = spool[i];
  let p = rec.xyz - u.eye.xyz;
  var r = 2.6;
  if (sp.isLipid == 0u) {
    let kind = bitcast<u32>(rec.w) & 15u;
    if (((u.modes.z >> kind) & 1u) != 0u) { return; }
    r = freeRadius(kind);
  }
  if (u.cut.x > 0.5) {
    let s = sPlaneW(p);
    if (s < -r || s - r > u.cut.y) { return; }
  }
  if (!sphereInFrustum(p, r)) { return; }
  let px = r * u.eye.w / max(length(p), 0.01);
  if (px < u.lod.z) { return; }
  if (sp.isLipid == 1u) {
    let k = atomicAdd(&counters.lip, 1u);
    if (k < sp.cap) { slist[k] = i; }
  } else {
    let k = atomicAdd(&counters.fre, 1u);
    if (k < sp.cap) { slist[k] = i; }
  }
}
`,$t=`
struct Counters { l0: u32, l1: u32, l2: u32, lip: u32, fre: u32, ninst: u32, p0: u32, p1: u32 };
struct Caps { capL0: u32, capL1: u32, capL2: u32, capLip: u32, capFree: u32, a: u32, b: u32, c: u32 };
@group(0) @binding(0) var<storage, read> counters: Counters;
@group(0) @binding(1) var<storage, read_write> args: array<u32>;
@group(0) @binding(2) var<uniform> caps: Caps;
@compute @workgroup_size(1)
fn main() {
  // drawIndirect args: vertexCount, instanceCount, firstVertex, firstInstance
  args[0] = 384u; args[1] = min(counters.l0, caps.capL0); args[2] = 0u; args[3] = 0u;
  args[4] = 96u;  args[5] = min(counters.l1, caps.capL1); args[6] = 0u; args[7] = 0u;
  args[8] = 6u;   args[9] = min(counters.l2, caps.capL2); args[10] = 0u; args[11] = 0u;
  args[12] = 18u; args[13] = min(counters.lip, caps.capLip); args[14] = 0u; args[15] = 0u;
  args[16] = 6u;  args[17] = min(counters.fre, caps.capFree); args[18] = 0u; args[19] = 0u;
}
`,wt=`
struct VOut {
  @builtin(position) pos: vec4f,
  @location(0) vpos: vec3f,
  @location(1) @interpolate(flat) center: vec3f,
  @location(2) @interpolate(flat) radius: f32,
  @location(3) @interpolate(flat) color: vec3f,
  @location(4) @interpolate(flat) id: u32,
};
struct FOut {
  @location(0) color: vec4f,
  @location(1) id: u32,
  @builtin(frag_depth) depth: f32,
};
const ID_INST: u32 = 0x40000000u;
const ID_LIPID: u32 = 0x80000000u;
const ID_FREE: u32 = 0xC0000000u;

fn emitQuad(cw: vec3f, r: f32, corner: u32, color: vec3f, id: u32) -> VOut {
  var o: VOut;
  let cv = (u.view * vec4f(cw, 1.0)).xyz;
  let d = length(cv);
  o.center = cv; o.radius = r; o.color = color; o.id = id;
  if (r <= 0.0 || d <= r * 1.02 || cv.z > -u.cell.w) {
    o.pos = vec4f(2.0, 2.0, 2.0, 1.0); o.vpos = vec3f(0.0); return o;
  }
  let k = d / sqrt(d * d - r * r);
  var off = vec2f(-1.0, -1.0);
  if (corner == 1u || corner == 4u) { off = vec2f(1.0, -1.0); }
  else if (corner == 2u || corner == 3u) { off = vec2f(-1.0, 1.0); }
  else if (corner == 5u) { off = vec2f(1.0, 1.0); }
  // billboard perpendicular to the view ray through the centre
  let fwd = cv / d;
  var up0 = vec3f(0.0, 1.0, 0.0);
  if (abs(fwd.y) > 0.99) { up0 = vec3f(1.0, 0.0, 0.0); }
  let right = normalize(cross(fwd, up0));
  let up = cross(right, fwd);
  let vp = cv + (right * off.x + up * off.y) * r * k * 1.05;
  o.vpos = vp;
  o.pos = u.proj * vec4f(vp, 1.0);
  return o;
}

fn shadeFrag(in: VOut) -> FOut {
  var o: FOut;
  let rd = normalize(in.vpos);
  let b = dot(rd, in.center);
  let c = dot(in.center, in.center) - in.radius * in.radius;
  var h = b * b - c;
  if (h < 0.0) { discard; }
  h = sqrt(h);
  var t = b - h;
  let t1 = b + h;
  var nrm = (rd * t - in.center) / in.radius;
  var cap = false;
  if (u.cut.x > 0.5) {
    let dn = dot(u.planeV.xyz, rd);
    let s0 = u.planeV.w + t * dn;
    if (s0 < 0.0) {
      if (dn <= 0.0) { discard; }
      let tp = -u.planeV.w / dn;
      if (tp > t1) { discard; }
      t = tp; nrm = -u.planeV.xyz; cap = true;
    }
    if (u.planeV.w + t * dn > u.cut.y) { discard; }
  }
  if (t < u.cell.w) { discard; }
  let p = rd * t;
  let clip = u.proj * vec4f(p, 1.0);
  o.depth = clip.z / clip.w;
  var shade = 1.0;
  if (u.sel.z == 1u) {
    let l = normalize(vec3f(0.35, 0.45, 1.0));
    shade = 0.45 + 0.55 * max(dot(nrm, l), 0.0) + 0.12 * pow(max(nrm.z, 0.0), 8.0);
  } else {
    shade = 0.93 + 0.07 * max(nrm.z, 0.0);
  }
  if (cap) { shade = 0.9; }
  o.color = vec4f(in.color * shade, 1.0);
  o.id = in.id;
  return o;
}
`,qt=`
@group(0) @binding(5) var<storage, read> palette: array<vec4f>;   // 0..31 region colours, 32.. special

fn hsv2rgb(h: f32, s: f32, v: f32) -> vec3f {
  let k = vec3f(1.0, 2.0 / 3.0, 1.0 / 3.0);
  let p = abs(fract(vec3f(h) + k) * 6.0 - 3.0);
  return v * mix(vec3f(1.0), clamp(p - 1.0, vec3f(0.0), vec3f(1.0)), s);
}
fn viridis(t0: f32) -> vec3f {
  let t = clamp(t0, 0.0, 1.0);
  let c0 = vec3f(0.2777, 0.0054, 0.3340); let c1 = vec3f(0.1050, 1.4046, 1.3846);
  let c2 = vec3f(-0.3308, 0.2148, 0.0950); let c3 = vec3f(-4.6342, -5.7991, -19.3324);
  let c4 = vec3f(6.2282, 14.1799, 56.6906); let c5 = vec3f(4.7763, -13.7451, -65.3530);
  let c6 = vec3f(-5.4355, 4.6459, 26.3124);
  return c0 + t * (c1 + t * (c2 + t * (c3 + t * (c4 + t * (c5 + t * c6)))));
}
fn plddtColor(p: f32) -> vec3f {
  if (p > 90.0) { return vec3f(0.0, 0.325, 0.839); }
  if (p > 70.0) { return vec3f(0.396, 0.796, 0.953); }
  if (p > 50.0) { return vec3f(1.0, 0.859, 0.075); }
  return vec3f(1.0, 0.490, 0.271);
}
// kind: 0 protein bead, 1 nucleic-acid bead, 2 coarse (sphere tree / single sphere)
fn instColor(I: Inst, T: TypeInfo, kind: u32, plddt: f32, id: u32) -> vec3f {
  let mode = u.modes.x;
  let hj = f32(hash32(I.t * 2654435761u + 17u) & 1023u) / 1023.0;
  var c: vec3f;
  if (mode == 0u) {          // Goodsell compartment palette with per-type lightness jitter
    if (T.color0 != 0u) {
      c = unpackColor(T.color0);
      if (kind == 1u && T.color1 != 0u) { c = unpackColor(T.color1); }
    } else { c = palette[I.region].xyz; }
    c = c * (0.93 + 0.14 * hj);
  } else if (mode == 1u) {   // structure source
    if ((T.flags & 1u) != 0u) { c = vec3f(0.30, 0.55, 0.85); } else { c = vec3f(0.96, 0.66, 0.30); }
    if (kind == 1u) { c = c * 0.8; }
  } else if (mode == 2u) {   // pLDDT
    var p = plddt;
    if (kind == 2u) { p = T.plddt; }
    if (p > 100.5 || p < 0.0) { c = vec3f(0.72, 0.72, 0.72); } else { c = plddtColor(p); }
  } else if (mode == 3u) {   // abundance (log10 copies, 1 .. 1e6)
    c = viridis(log(max(f32(T.flags >> 8u), 1.0)) / log(10.0) / 6.0);
  } else if (mode == 4u) {   // random per type
    c = hsv2rgb(hj, 0.55, 0.92);
  } else {                    // compartment only
    c = palette[I.region].xyz;
  }
  if (u.modes.w == 1u && (typeFlags[I.t] & 1u) == 0u) {
    let g = dot(c, vec3f(0.3, 0.59, 0.11));
    c = mix(vec3f(g), vec3f(0.93), 0.72);
  }
  if (u.sel.x == id && id != 0u) { c = mix(c, vec3f(1.0, 0.92, 0.2), 0.55); }
  return c;
}
`,ot=ye+vt+wt+qt,Vt=ot+`
@group(1) @binding(0) var<storage, read> beads: array<vec4f>;     // xyz, w = bitcast(kind | plddt<<8)
@group(1) @binding(1) var<storage, read> list01: array<vec2u>;
@vertex
fn vsL0(@builtin(vertex_index) vi: u32, @builtin(instance_index) ii: u32) -> VOut {
  let item = list01[ii];
  let j = vi / 6u; let corner = vi % 6u;
  let cnt = (item.y >> 22u) + 1u;
  if (j >= cnt) { var o: VOut; o.pos = vec4f(2.0, 2.0, 2.0, 1.0); return o; }
  let I = decodeInst(item.x);
  let b = beads[(item.y & 0x3FFFFFu) + j];
  let bw = bitcast<u32>(b.w);
  let kind = bw & 255u;
  let T = types[I.t];
  var r = u.lod.w;
  if (kind == 1u) { r = u.misc.y; }
  let id = ID_INST | item.x;
  let col = instColor(I, T, kind, f32((bw >> 8u) & 255u), id);
  return emitQuad(I.p + qrot(I.q, b.xyz), r, corner, col, id);
}
@fragment
fn fsMain(in: VOut) -> FOut { return shadeFrag(in); }
`,Wt=ot+`
@group(1) @binding(0) var<storage, read> stSpheres: array<vec4f>; // xyz r
@group(1) @binding(1) var<storage, read> list01: array<vec2u>;
@vertex
fn vsL1(@builtin(vertex_index) vi: u32, @builtin(instance_index) ii: u32) -> VOut {
  let item = list01[u.caps.x + ii];
  let j = vi / 6u; let corner = vi % 6u;
  let cnt = (item.y >> 22u) + 1u;
  if (j >= cnt) { var o: VOut; o.pos = vec4f(2.0, 2.0, 2.0, 1.0); return o; }
  let I = decodeInst(item.x);
  let s = stSpheres[(item.y & 0x3FFFFFu) + j];
  let T = types[I.t];
  let id = ID_INST | item.x;
  // sphere-tree spheres are padded for collision; shrink slightly for display
  return emitQuad(I.p + qrot(I.q, s.xyz), s.w * 0.9, corner, instColor(I, T, 2u, -1.0, id), id);
}
@fragment
fn fsMain(in: VOut) -> FOut { return shadeFrag(in); }
`,jt=ot+`
@group(1) @binding(0) var<storage, read> listL2: array<u32>;
@vertex
fn vsL2(@builtin(vertex_index) vi: u32, @builtin(instance_index) ii: u32) -> VOut {
  let i = listL2[ii];
  let I = decodeInst(i);
  let T = types[I.t];
  let id = ID_INST | i;
  return emitQuad(I.p, T.req, vi, instColor(I, T, 2u, -1.0, id), id);
}
@fragment
fn fsMain(in: VOut) -> FOut { return shadeFrag(in); }
`,Ht=ye+wt+`
@group(1) @binding(0) var<storage, read> spool: array<vec4f>;
@group(1) @binding(1) var<storage, read> spages: array<vec4u>;
@group(1) @binding(2) var<storage, read> slist: array<u32>;
@group(0) @binding(5) var<storage, read> palette: array<vec4f>;

fn octDecode(e: vec2f) -> vec3f {
  var v = vec3f(e.x, e.y, 1.0 - abs(e.x) - abs(e.y));
  if (v.z < 0.0) { let xy = (1.0 - abs(v.yx)) * select(vec2f(-1.0), vec2f(1.0), v.xy >= vec2f(0.0)); v = vec3f(xy, v.z); }
  return normalize(v);
}

// lipid record: xyz = midsurface point, w = bits: oct normal (8+8), class (3) <<16, leaflet sign <<19, membrane (4) <<20
@vertex
fn vsLipid(@builtin(vertex_index) vi: u32, @builtin(instance_index) ii: u32) -> VOut {
  let i = slist[ii];
  let rec = spool[i];
  let w = bitcast<u32>(rec.w);
  let n0 = octDecode(vec2f(f32(w & 255u), f32((w >> 8u) & 255u)) / 127.5 - 1.0);
  let cls = (w >> 16u) & 7u;
  let sgn = select(-1.0, 1.0, ((w >> 19u) & 1u) == 1u);
  let mem = (w >> 20u) & 15u;
  let j = vi / 6u; let corner = vi % 6u;
  var off = 1.95; var r = 0.55;
  var col = palette[40u].xyz;                  // head
  if (cls == 5u) {                              // ergosterol: two beads, rigid ring system
    col = palette[42u].xyz;
    if (j == 0u) { off = 1.55; r = 0.45; } else if (j == 1u) { off = 0.75; r = 0.45; } else { r = 0.0; }
  } else {
    col = palette[32u + cls].xyz;
    if (j == 1u) { off = 1.25; r = 0.45; col = palette[41u].xyz; }
    else if (j == 2u) { off = 0.55; r = 0.45; col = palette[41u].xyz; }
  }
  let slot = spages[i / SPAGE].y;
  let id = ID_LIPID | (slot << 5u) | (mem << 1u) | select(0u, 1u, sgn > 0.0);
  let p = rec.xyz - u.eye.xyz + n0 * sgn * off;
  return emitQuad(p, r, corner, col, id);
}

// free sphere record: xyz, w = bits: kind (4) | chain (16) << 4
@vertex
fn vsFree(@builtin(vertex_index) vi: u32, @builtin(instance_index) ii: u32) -> VOut {
  let i = slist[ii];
  let rec = spool[i];
  let w = bitcast<u32>(rec.w);
  let kind = w & 15u;
  let chain = (w >> 4u) & 0xFFFFu;
  var r = 0.6;
  if (kind == 1u) { r = 1.0; } else if (kind == 2u) { r = 0.75; }
  let hj = f32(hash32(chain * 747796405u + kind) & 255u) / 255.0;
  let col = palette[48u + kind].xyz * (0.94 + 0.12 * hj);
  let slot = spages[i / SPAGE].y;
  let id = ID_FREE | ((slot & 0xFFFu) << 18u) | ((kind & 3u) << 16u) | chain;
  return emitQuad(rec.xyz - u.eye.xyz, r, vi, col, id);
}

@fragment
fn fsMain(in: VOut) -> FOut { return shadeFrag(in); }
`,Yt=ye+`
@group(0) @binding(5) var<storage, read> palette: array<vec4f>;
struct MP { slot: u32, organelle: u32, offset0: f32, offset1: f32 };
@group(1) @binding(0) var<uniform> mp: MP;
struct VOut {
  @builtin(position) pos: vec4f,
  @location(0) prel: vec3f,
  @location(1) nrm: vec3f,
  @location(2) @interpolate(flat) mem: u32,
  @location(3) @interpolate(flat) leaf: u32,
};
@vertex
fn vsMesh(@location(0) p: vec4f, @location(1) n: vec4f, @builtin(instance_index) ii: u32) -> VOut {
  var o: VOut;
  let mem = bitcast<u32>(p.w);
  var off = select(mp.offset0, mp.offset1, ii == 1u);
  if (mem == 8u) { off = select(-0.8, 0.8, ii == 1u); }   // lipid droplet monolayer
  let pr = p.xyz + n.xyz * off - u.eye.xyz;
  o.prel = pr; o.nrm = n.xyz; o.mem = mem; o.leaf = ii;
  o.pos = u.proj * (u.view * vec4f(pr, 1.0));
  return o;
}
struct FOut { @location(0) color: vec4f, @location(1) id: u32 };
@fragment
fn fsMesh(in: VOut, @builtin(front_facing) ff: bool) -> FOut {
  var o: FOut;
  if (u.cut.x > 0.5) {
    let s = sPlaneW(in.prel);
    if (s < 0.0 || (mp.organelle != 1u && s > u.cut.y)) { discard; }
  }
  let vn = normalize((u.view * vec4f(in.nrm, 0.0)).xyz);
  var c = palette[40u].xyz;
  if (mp.organelle == 1u) { c = palette[56u + in.mem].xyz; } else if (mp.organelle == 2u) { c = palette[41u].xyz; }
  let l = abs(vn.z);
  var shade = 0.82 + 0.18 * l;
  if (mp.organelle == 1u) { shade = 0.55 + 0.45 * l; if (!ff) { shade *= 0.8; } }
  o.color = vec4f(c * shade, 1.0);
  o.id = 0x80000000u | (mp.slot << 5u) | (in.mem << 1u) | in.leaf;
  return o;
}
`,Xt=ye+`
@group(0) @binding(5) var<storage, read> palette: array<vec4f>;
struct Prim { q0: vec4f, q1: vec4f, tc: vec4u };      // type, class
struct Geo { G0: vec4f, G1: vec4f, G2: vec4f, consts: vec4f, grid: vec4f };   // consts: PORE_R, PERI_W, IMS_W, nPores
@group(1) @binding(0) var<storage, read> prims: array<Prim>;
@group(1) @binding(1) var<storage, read> pores: array<vec4f>;
@group(1) @binding(2) var<storage, read> gridOff: array<u32>;   // prim list offsets (n^3+1) then pore list offsets
@group(1) @binding(3) var<storage, read> gridIdx: array<u32>;
@group(1) @binding(4) var<uniform> geo: Geo;
@group(1) @binding(5) var<storage, read> tileState: array<u32>;   // 22^3, bit0 = molecules loaded

const NMEM: u32 = 10u;
fn sdfPrim(p: vec3f, pr: Prim) -> f32 {
  let t = pr.tc.x; let q0 = pr.q0; let q1 = pr.q1;
  if (t == 0u) { return length(p - q0.xyz) - q0.w; }
  if (t == 1u) {
    let a = p - q0.xyz; let b = vec3f(q0.w, q1.x, q1.y) - q0.xyz;
    let h = clamp(dot(a, b) / dot(b, b), 0.0, 1.0);
    return length(a - b * h) - q1.z;
  }
  if (t == 2u) {
    let v = p - q0.xyz; let nn = vec3f(q0.w, q1.x, q1.y);
    let z = dot(v, nn); let rho = length(v - z * nn);
    let s = max(rho - q1.z, 0.0);
    return sqrt(s * s + z * z) - q1.w;
  }
  let r = length(p);
  let c = clamp(dot(p, q0.xyz) / max(r, 1e-9), -1.0, 1.0);
  let th = acos(c);
  let s = max(th - q1.x, 0.0) * r;
  let dr = r - q1.y;
  return sqrt(s * s + dr * dr) - q1.z;
}
struct Mem { d: array<f32, 10>, matrix: f32, cell: u32 };
fn gridCell(p: vec3f) -> u32 {
  let n = u32(geo.grid.x); let ext = geo.grid.y; let h = 2.0 * ext / f32(n);
  let g = clamp(vec3i(floor((p + vec3f(ext)) / h)), vec3i(0), vec3i(i32(n) - 1));
  return (u32(g.x) * n + u32(g.y)) * n + u32(g.z);
}
fn memSdfs(p: vec3f) -> Mem {
  var m: Mem;
  let G = geo;
  let PERI = G.consts.y; let PORE_R = G.consts.x;
  let r = length(p);
  m.d[0] = G.G0.z - r;     // PM
  m.d[9] = G.G0.x - r;     // outer
  let v = p - vec3f(G.G0.w, G.G1.x, G.G1.y);
  let dn = length(v);
  let hp = dn - (G.G1.z - 0.5 * PERI);
  var d = abs(hp) - 0.5 * PERI;
  let rc = PORE_R + 0.5 * PERI;
  let cell = gridCell(p); m.cell = cell;
  let n3 = u32(geo.grid.x * geo.grid.x * geo.grid.x);
  let po0 = gridOff[n3 + 1u + cell]; let po1 = gridOff[n3 + 1u + cell + 1u];
  for (var k = po0; k < po1; k++) {
    let pd = pores[gridIdx[k]].xyz;
    let z = dot(v, pd);
    if (z <= 0.0) { continue; }
    let rho = length(v - z * pd);
    if (rho < rc) { d = sqrt((rho - rc) * (rho - rc) + hp * hp) - 0.5 * PERI; break; }
  }
  m.d[2] = d;
  m.d[1] = 1e9; m.d[3] = 1e9; m.d[4] = 1e9; m.d[6] = 1e9; m.d[7] = 1e9; m.d[8] = 1e9;
  var im = 1e9; var cr = 1e9;
  let o0 = gridOff[cell]; let o1 = gridOff[cell + 1u];
  for (var k = o0; k < o1; k++) {
    let pr = prims[gridIdx[k]];
    let dd = sdfPrim(p, pr);
    let c = pr.tc.y;
    if (c == 7u) { m.d[1] = min(m.d[1], dd); }
    else if (c == 0u) { m.d[3] = min(m.d[3], dd); }
    else if (c == 4u) { m.d[4] = min(m.d[4], dd); }
    else if (c == 5u) { im = min(im, dd); }
    else if (c == 6u) { cr = min(cr, dd); }
    else if (c == 2u) { m.d[6] = min(m.d[6], dd); }
    else if (c == 3u) { m.d[7] = min(m.d[7], dd); }
    else if (c == 1u) { m.d[8] = min(m.d[8], dd); }
  }
  m.matrix = max(im, -cr);
  m.d[5] = -m.matrix;
  return m;
}
// labels: 0 EXT,1 MANNAN,2 GLUCAN,3 CYTO,4 CER,5 PERI,6 NUCP,7 NUCL,8 VAC,9 IMS,10 MATRIX,11 LD,12 GOLGI,13 PEROX,14 PORE
fn labelOf(p: vec3f, m: Mem) -> u32 {
  let G = geo; let PERI = G.consts.y; let PORE_R = G.consts.x;
  let r = length(p);
  if (r > G.G0.x) { return 0u; }
  if (r > G.G0.y) { return 1u; }
  if (r > G.G0.z) { return 2u; }
  if (m.d[1] < 0.0) { return 4u; }
  if (m.d[2] < 0.0) { return 5u; }
  let v = p - vec3f(G.G0.w, G.G1.x, G.G1.y);
  let dn = length(v);
  if (dn < G.G1.z - 0.5 * PERI) {
    let hp = dn - (G.G1.z - 0.5 * PERI);
    let rc = PORE_R + 0.5 * PERI;
    if (hp > -0.5 * PERI - 8.0) {
      let n3 = u32(geo.grid.x * geo.grid.x * geo.grid.x);
      let po0 = gridOff[n3 + 1u + m.cell]; let po1 = gridOff[n3 + 1u + m.cell + 1u];
      for (var k = po0; k < po1; k++) {
        let pd = pores[gridIdx[k]].xyz;
        let z = dot(v, pd);
        if (z <= 0.0) { continue; }
        let rr = v - z * pd;
        if (dot(rr, rr) < rc * rc) { return 14u; }
      }
    }
    let w = p - vec3f(G.G1.w, G.G2.x, G.G2.y);
    if (dn < G.G1.z - PERI && dot(w, w) < G.G2.z * G.G2.z) { return 7u; }
    return 6u;
  }
  if (dn < G.G1.z) { return 14u; }
  if (m.d[4] < 0.0) { if (m.matrix < 0.0) { return 10u; } return 9u; }
  if (m.d[3] < 0.0) { return 8u; }
  if (m.d[8] < 0.0) { return 11u; }
  if (m.d[6] < 0.0) { return 12u; }
  if (m.d[7] < 0.0) { return 13u; }
  return 3u;
}
fn tileLoaded(p: vec3f) -> bool {
  let g = vec3i(floor(p / 200.0)) + vec3i(11);
  if (any(g < vec3i(0)) || any(g > vec3i(21))) { return false; }
  return (tileState[(u32(g.x) * 22u + u32(g.y)) * 22u + u32(g.z)] & 1u) != 0u;
}

struct VO { @builtin(position) pos: vec4f, @location(0) ndc: vec2f };
@vertex
fn vsFill(@builtin(vertex_index) vi: u32) -> VO {
  var o: VO;
  let xy = vec2f(f32((vi << 1u) & 2u), f32(vi & 2u)) * 2.0 - 1.0;
  o.pos = vec4f(xy, 0.0, 1.0); o.ndc = xy;
  return o;
}
struct FOut { @location(0) color: vec4f, @location(1) id: u32, @builtin(frag_depth) depth: f32 };

fn shadeLabel(p: vec3f, o: ptr<function, FOut>, bands: bool) {
  let m = memSdfs(p);
  let lab = labelOf(p, m);
  var c = palette[80u + lab].xyz;
  var id = 0x100u + lab;
  if (bands) {
    for (var k = 0u; k < 9u; k++) {
      let d = m.d[k];
      if (k == 8u) {           // lipid droplet: phospholipid monolayer on the outside of the neutral-lipid core
        if (d > -2.0 && d < 0.0) { c = select(palette[41u].xyz, palette[40u].xyz, d > -0.9); id = 0x200u + k; }
        continue;
      }
      let a = abs(d);
      if (a < 2.5) { c = select(palette[41u].xyz, palette[40u].xyz, a > 1.35); id = 0x200u + k; }
    }
  }
  (*o).color = vec4f(c, 1.0);
  (*o).id = id;
}

@fragment
fn fsFill(in: VO) -> FOut {
  var o: FOut;
  let pw = u.invViewProj * vec4f(in.ndc, 1.0, 1.0);    // near plane point (camera-relative)
  let rd = normalize(pw.xyz / pw.w);
  let eye = u.eye.xyz;
  let R = u.cell.x;
  // ray vs outer cell sphere (centre at origin in world)
  let b = dot(rd, eye);
  let c = dot(eye, eye) - R * R;
  let h = b * b - c;
  o.color = vec4f(palette[80u].xyz, 1.0); o.id = 0u; o.depth = 0.0;
  if (h < 0.0) { return o; }
  let sh = sqrt(h);
  var t0 = max(-b - sh, 0.0);
  let t1 = -b + sh;
  if (t1 <= 0.0) { return o; }
  var onCut = false;
  if (u.cut.x > 0.5) {
    let dn = dot(u.plane.xyz, rd);
    let sE = u.planeV.w;
    if (abs(dn) < 1e-6) { if (sE < 0.0) { return o; } }
    else {
      let tp = -sE / dn;
      if (dn > 0.0) { if (tp > t1) { return o; } if (tp > t0) { t0 = tp; onCut = true; } }
      else { if (tp < t0) { return o; } }
    }
    if (u.sel.w == 1u) {
      // organelle mode: cut face only (meshes provide the rest)
      if (!onCut) { o.color = vec4f(palette[80u].xyz, 1.0); return o; }
    }
  } else if (u.sel.w == 0u) {
    if (c < 0.0) { return o; }      // camera inside the cell without a cut: nothing analytic to show
  }
  var t = t0 + select(0.5, 0.0, onCut);
  var p = eye + rd * t;
  var bands = true;
  if (u.sel.w == 0u && tileLoaded(p)) {
    // molecules are drawn here: show the compartment backdrop behind the visible slab
    if (onCut) { t = t + u.cut.y / max(dot(u.plane.xyz, rd), 0.05); }
    else { t = t + 140.0; }
    if (t >= t1) { o.color = vec4f(palette[80u].xyz, 1.0); return o; }
    p = eye + rd * t;
  }
  shadeLabel(p, &o, bands);
  let pv = (u.view * vec4f(rd * t, 1.0)).xyz;
  let clip = u.proj * vec4f(pv, 1.0);
  o.depth = clamp(clip.z / clip.w, 0.0, 1.0);
  if (!onCut && u.cut.x < 0.5) {
    // outer surface: gentle sphere shading
    let nn = normalize(p);
    let vn = (u.view * vec4f(nn, 0.0)).xyz;
    o.color = vec4f(o.color.xyz * (0.75 + 0.25 * max(vn.z, 0.0)), 1.0);
  }
  return o;
}
`,Qt=`
@group(0) @binding(0) var idTex: texture_2d<u32>;
@group(0) @binding(1) var depthTex: texture_depth_2d;
@group(0) @binding(2) var<uniform> px: vec4u;
@group(0) @binding(3) var<storage, read_write> outp: array<u32, 4>;
@compute @workgroup_size(1)
fn main() {
  let p = vec2i(px.xy);
  outp[0] = textureLoad(idTex, p, 0).x;
  outp[1] = bitcast<u32>(textureLoad(depthTex, p, 0));
}
`,Kt=ye+`
@group(1) @binding(0) var colorTex: texture_2d<f32>;
@group(1) @binding(1) var idTex: texture_2d<u32>;
@group(1) @binding(2) var depthTex: texture_depth_2d;
struct VO { @builtin(position) pos: vec4f };
@vertex
fn vsFull(@builtin(vertex_index) vi: u32) -> VO {
  var o: VO;
  let xy = vec2f(f32((vi << 1u) & 2u), f32(vi & 2u)) * 2.0 - 1.0;
  o.pos = vec4f(xy, 0.0, 1.0);
  return o;
}
fn viewPos(px: vec2i, d: f32) -> vec3f {
  let ndc = vec2f((f32(px.x) + 0.5) * u.vp.z * 2.0 - 1.0, 1.0 - (f32(px.y) + 0.5) * u.vp.w * 2.0);
  let w = u.invViewProj * vec4f(ndc, max(d, 1e-7), 1.0);
  return w.xyz / w.w;
}
@fragment
fn fsComposite(in: VO) -> @location(0) vec4f {
  let px = vec2i(in.pos.xy);
  let dims = vec2i(textureDimensions(colorTex));
  var c = textureLoad(colorTex, px, 0).xyz;
  let id = textureLoad(idTex, px, 0).x;
  let d = textureLoad(depthTex, px, 0);
  let cat = id >> 30u;
  if (d <= 0.0) { return vec4f(c, 1.0); }
  let prel = viewPos(px, d);
  let dist = length(prel);
  let pxPerNm = u.eye.w / max(dist, 1e-3);
  // depth cue behind the cut face (Goodsell-style darkening with depth)
  if (u.cut.x > 0.5 && cat != 0u && u.sel.w == 0u) {
    let s = sPlaneW(prel);
    let f = smoothstep(0.0, u.cut.z, s);
    c = mix(c, c * vec3f(0.55, 0.56, 0.62), f * 0.8);
  }
  // outlines at molecule (id) boundaries; strength fades when molecules are only a few pixels wide
  let rad = i32(max(1.0, round(u.misc.w)));
  var edge = 0.0;
  var selEdge = 0.0;
  let offs = array<vec2i, 4>(vec2i(rad, 0), vec2i(-rad, 0), vec2i(0, rad), vec2i(0, -rad));
  for (var k = 0; k < 4; k++) {
    let q = clamp(px + offs[k], vec2i(0), dims - vec2i(1));
    let id2 = textureLoad(idTex, q, 0).x;
    if (id2 != id) {
      let d2 = textureLoad(depthTex, q, 0);
      if (cat == 0u) { if ((id2 >> 30u) == 0u) { edge = 1.0; } }
      else if (d2 <= d * 1.0005 || (id2 >> 30u) == 0u) { edge = 1.0; }
      if (u.sel.x != 0u && (id2 == u.sel.x || id == u.sel.x)) { selEdge = 1.0; }
    }
  }
  var strength = u.misc.z;
  if (cat == 1u) { strength *= smoothstep(0.6, 2.5, pxPerNm); }
  else if (cat == 2u) { strength *= 0.55 * smoothstep(1.5, 6.0, pxPerNm); }
  else if (cat == 3u) { strength *= 0.7 * smoothstep(1.5, 5.0, pxPerNm); }
  else { strength *= 0.35; }
  c = mix(c, c * 0.18, edge * strength);
  if (selEdge > 0.0) { c = vec3f(1.0, 0.85, 0.1); }
  return vec4f(c, 1.0);
}
`,_=(s,e)=>[s[0]+e[0],s[1]+e[1],s[2]+e[2]],ae=(s,e)=>[s[0]-e[0],s[1]-e[1],s[2]-e[2]],L=(s,e)=>[s[0]*e,s[1]*e,s[2]*e],$=(s,e)=>s[0]*e[0]+s[1]*e[1]+s[2]*e[2],de=(s,e)=>[s[1]*e[2]-s[2]*e[1],s[2]*e[0]-s[0]*e[2],s[0]*e[1]-s[1]*e[0]],ge=s=>Math.hypot(s[0],s[1],s[2]),X=s=>{const e=ge(s)||1;return[s[0]/e,s[1]/e,s[2]/e]},it=(s,e,t)=>[s[0]+(e[0]-s[0])*t,s[1]+(e[1]-s[1])*t,s[2]+(e[2]-s[2])*t];function Be(s,e,t){const i=Math.cos(t),n=Math.sin(t),r=de(e,s),o=$(e,s)*(1-i);return[s[0]*i+r[0]*n+e[0]*o,s[1]*i+r[1]*n+e[1]*o,s[2]*i+r[2]*n+e[2]*o]}function Jt(s,e){const t=new Float32Array(16);for(let i=0;i<4;i++)for(let n=0;n<4;n++){let r=0;for(let o=0;o<4;o++)r+=s[o*4+n]*e[i*4+o];t[i*4+n]=r}return t}function Zt(s){const e=Array.from(s),t=new Float32Array(16),[i,n,r,o,c,d,a,u,m,l,h,f,M,v,p,b]=e,P=i*d-n*c,x=i*a-r*c,A=i*u-o*c,k=n*a-r*d,F=n*u-o*d,C=r*u-o*a,B=m*v-l*M,R=m*p-h*M,q=m*b-f*M,Q=l*p-h*v,K=l*b-f*v,V=h*b-f*p;let D=P*V-x*K+A*Q+k*q-F*R+C*B;return D&&(D=1/D,t[0]=(d*V-a*K+u*Q)*D,t[1]=(r*K-n*V-o*Q)*D,t[2]=(v*C-p*F+b*k)*D,t[3]=(h*F-l*C-f*k)*D,t[4]=(a*q-c*V-u*R)*D,t[5]=(i*V-r*q+o*R)*D,t[6]=(p*A-M*C-b*x)*D,t[7]=(m*C-h*A+f*x)*D,t[8]=(c*K-d*q+u*B)*D,t[9]=(n*q-i*K-o*B)*D,t[10]=(M*F-v*A+b*P)*D,t[11]=(l*A-m*F-f*P)*D,t[12]=(d*R-c*Q-a*B)*D,t[13]=(i*Q-n*R+r*B)*D,t[14]=(v*x-M*k-p*P)*D,t[15]=(m*k-l*x+h*P)*D),t}function ut(s){const e=parseInt(s.replace("#",""),16);return(e>>16&255|(e>>8&255)<<8|(e&255)<<16|4278190080)>>>0}function ei(s){return s<.5?4*s*s*s:1-Math.pow(-2*s+2,3)/2}const xt=new URL("data/",document.baseURI).href,ti=fetch(xt+"version.json",{cache:"no-store"}).then(s=>s.ok?s.json():{v:""}).then(s=>String(s.v??"")).catch(()=>""),ii=s=>s?`?v=${encodeURIComponent(s)}`:"";function Mt(s){const e=atob(s),t=new Uint8Array(e.length);for(let i=0;i<e.length;i++)t[i]=e.charCodeAt(i);return t}async function Xe(s){const e=new Uint8Array(s,0,Math.min(2,s.byteLength));if(e[0]===31&&e[1]===139){const t=new DecompressionStream("gzip"),i=new Blob([s]).stream().pipeThrough(t);return await new Response(i).arrayBuffer()}return s}async function Ie(s,e){const t=window.__YC_EMBED__;if(t&&t[s]!==void 0)return Xe(Mt(t[s]).buffer);const i=await fetch(xt+s+ii(await ti));if(!i.ok)throw new Error(`HTTP ${i.status} for ${s}`);if(!e||!i.body)return Xe(await i.arrayBuffer());const n=Number(i.headers.get("content-length")||0),r=i.body.getReader(),o=[];let c=0;for(;;){const{done:u,value:m}=await r.read();if(u)break;o.push(m),c+=m.length,e(c,n)}const d=new Uint8Array(c);let a=0;for(const u of o)d.set(u,a),a+=u.length;return Xe(d.buffer)}async function qe(s){const e=await Ie(s);return JSON.parse(new TextDecoder().decode(e))}const Ye=["cytosol","nucleoplasm","nucleolus","ER lumen","vacuole lumen","mitochondrial matrix","mitochondrial IMS","Golgi / endosome","peroxisome","cell wall","PM peripheral","NE peripheral","ER peripheral","lipid droplet surface","plasma membrane","ER membrane","nuclear envelope","vacuolar membrane","mito outer membrane","mito inner membrane","Golgi membrane","peroxisomal membrane","chromatin (nucleosome)","nuclear pore"],dt=["plasma membrane","ER membrane","nuclear envelope","vacuolar membrane","mito outer membrane","mito inner membrane","Golgi/endosome membrane","peroxisomal membrane","lipid droplet monolayer","cell surface"],si=["exterior","cell wall (mannoprotein layer)","cell wall (glucan/chitin layer)","cytosol","cortical ER lumen","perinuclear space","nucleoplasm","nucleolus","vacuole","mito intermembrane space","mito matrix","lipid droplet","Golgi / endosome","peroxisome","nuclear pore channel"],N={0:"#7EA7D3",1:"#B993C9",2:"#D39BC0",3:"#A3C48C",4:"#9FC7C0",5:"#D9946B",6:"#E8B48F",7:"#B7C77E",8:"#C5B0E0",9:"#9A8352",10:"#7EA7D3",11:"#8FA3D0",12:"#7EA7D3",13:"#7EA7D3",14:"#5FAE75",15:"#7DBE86",16:"#6BB596",17:"#5BA8A0",18:"#C27C55",19:"#B5653F",20:"#93B060",21:"#A68BC9",22:"#F4B26A",23:"#3E9C9A",32:"#F3C643",33:"#F0BC3C",34:"#F6D15E",35:"#EAB542",36:"#EFA83A",37:"#E7A93C",38:"#F2CA55",40:"#F3C643",41:"#F8E3A0",42:"#E7A93C",48:"#B0437A",49:"#E07636",50:"#CDA866",51:"#C3C27A",56:"#5FAE75",57:"#7DBE86",58:"#6BB596",59:"#5BA8A0",60:"#C27C55",61:"#B5653F",62:"#93B060",63:"#A68BC9",64:"#E9C46A",65:"#9A8352",80:"#F6F1E7",81:"#F1E9D3",82:"#EEE3CC",83:"#E4ECF5",84:"#E8F0DF",85:"#E8F0DF",86:"#F1E7F1",87:"#EBD9EA",88:"#E3F0EE",89:"#F6E6DA",90:"#F2DCCB",91:"#F7EDC4",92:"#EEF1DC",93:"#EEE8F6",94:"#ECE8F2"};function ni(){const s=new Float32Array(384);for(const[e,t]of Object.entries(N)){const i=parseInt(t.slice(1),16),n=Number(e)*4;s[n]=(i>>16&255)/255,s[n+1]=(i>>8&255)/255,s[n+2]=(i&255)/255,s[n+3]=1}return s}const Pt={ribosome_80S:{c0:"#EE9BC0",c1:"#C9558F",cls:"ribosome"},tRNA:{c0:"#E27FA8",c1:"#E27FA8",cls:"tRNA"},nucleosome:{c0:"#F4B26A",c1:"#E07636",cls:"nucleosome"},NPC:{c0:"#3E9C9A",cls:"NPC"},PolI:{c0:"#9A72B8",c1:"#9A72B8",cls:"polymerase"},PolII:{c0:"#9A72B8",c1:"#9A72B8",cls:"polymerase"},PolIII:{c0:"#9A72B8",c1:"#9A72B8",cls:"polymerase"},FAS:{c0:"#8B88C9",cls:"complex"},proteasome_26S:{c0:"#8B88C9",cls:"complex"},CCT:{c0:"#8B88C9",cls:"complex"},F_actin:{c0:"#E0893A",c1:"#E0893A",cls:"actin"},MT_tubulin:{c0:"#74C476",c1:"#74C476",cls:"microtubule"},gTuRC_nuc:{c0:"#2E8B57",c1:"#2E8B57",cls:"microtubule"},gTuRC_cyto:{c0:"#2E8B57",c1:"#2E8B57",cls:"microtubule"},ribosome_80S_Ssh1:{c0:"#EE9BC0",c1:"#C9558F",cls:"ribosome"},ribosome_80S_Sec61:{c0:"#EE9BC0",c1:"#C9558F",cls:"ribosome"}};async function ri(s){s("Loading structure index",.02);const e=await qe("library/types.json"),t=e.types,i=e.bead_q;s("Loading residue coordinates (2.5 M beads)",.06);const n=await Ie("library/beads.bin",(x,A)=>s("Loading residue coordinates (2.5 M beads)",.06+.5*(A?x/A:0))),r=new Float32Array(await Ie("library/spheretree.bin"));s("Building level-of-detail clusters",.6);const o=n.byteLength/8,c=new Int16Array(n),d=new Uint8Array(n),a=new Float32Array(o*4),u=new Uint32Array(a.buffer);for(let x=0;x<o;x++)a[4*x]=c[4*x]*i,a[4*x+1]=c[4*x+1]*i,a[4*x+2]=c[4*x+2]*i,u[4*x+3]=d[8*x+6]|d[8*x+7]<<8;let m=0;for(const x of t)m+=Math.ceil(x.nb/64);const l=new ArrayBuffer(m*32),h=new Float32Array(l),f=new Uint32Array(l),M=new ArrayBuffer(t.length*48),v=new Uint32Array(M),p=new Float32Array(M);let b=0;t.forEach((x,A)=>{const k=b;for(let B=0;B<x.nb;B+=64){const R=Math.min(64,x.nb-B),q=x.bo+B;let Q=1e9,K=1e9,V=1e9,D=-1e9,Le=-1e9,Ee=-1e9;for(let ce=0;ce<R;ce++){const W=4*(q+ce),se=a[W],G=a[W+1],j=a[W+2];se<Q&&(Q=se),se>D&&(D=se),G<K&&(K=G),G>Le&&(Le=G),j<V&&(V=j),j>Ee&&(Ee=j)}const _e=(Q+D)/2,ze=(K+Le)/2,Re=(V+Ee)/2;let Te=0;for(let ce=0;ce<R;ce++){const W=4*(q+ce);Te=Math.max(Te,Math.hypot(a[W]-_e,a[W+1]-ze,a[W+2]-Re))}h[8*b]=_e,h[8*b+1]=ze,h[8*b+2]=Re,h[8*b+3]=Te+.85,f[8*b+4]=q,f[8*b+5]=R,b++}const F=Pt[x.id],C=x.src==="PDB";v[12*A+0]=x.bo,v[12*A+1]=x.nb,v[12*A+2]=k,v[12*A+3]=b-k,v[12*A+4]=x.so,v[12*A+5]=x.sn,p[12*A+6]=x.br,p[12*A+7]=Math.sqrt(5/3)*x.rg,v[12*A+8]=F?ut(F.c0):0,v[12*A+9]=F&&F.c1?ut(F.c1):0,v[12*A+10]=((C?1:0)|Math.min(x.copies,16777215)<<8)>>>0,p[12*A+11]=C?255:x.plddt_mean??255});const P=new Map;return t.forEach((x,A)=>P.set(x.id,A)),{types:t,beads:a,st:r,clusters:l,typeInfo:M,nClusters:m,byId:P}}const st=(s,e,t)=>((s+11)*22+(e+11))*22+(t+11);async function oi(){const s=await qe("tiles/index.json"),e=s.tile_columns,t=s.tiles,i=t.length,n=c=>e.indexOf(c),r=(c,d)=>{const a=new c(i),u=n(d);for(let m=0;m<i;m++)a[m]=t[m][u];return a},o={n:i,ix:r(Int32Array,"ix"),iy:r(Int32Array,"iy"),iz:r(Int32Array,"iz"),nInst:r(Uint32Array,"n_inst"),nGhost:r(Uint32Array,"n_ghost"),nMr:r(Uint32Array,"n_mrna_beads"),nLk:r(Uint32Array,"n_link_beads"),nVert:r(Uint32Array,"n_vert"),nTri:r(Uint32Array,"n_tri"),gz:r(Uint32Array,"gz_bytes"),labelMask:r(Uint32Array,"label_mask"),volGlucan:r(Float32Array,"vol_glucan_nm3"),volMannan:r(Float32Array,"vol_mannan_nm3"),key:new Map,wall:s.wall,regionNames:s.region_names};for(let c=0;c<i;c++)o.key.set(st(o.ix[c],o.iy[c],o.iz[c]),c);return o}const ai=1073741824,O=GPUShaderStage,E=GPUBufferUsage;class li{device;ctx;format;canvas;limits={};adapterInfo="";instPages=0;lipPages=0;freePages=0;capL0=0;capL1=0;capL2=0;capLip=0;capFree=0;maxSlots=4096;uni;inst;pageOrg;typesBuf;typeFlags;palette;beads;st;clusters;counters;list01;listL2;args;cullParams;caps;lipPool;lipPagesBuf;lipList;lipParams;freePool;freePagesBuf;freeList;freeParams;meshUni;tileState;fillBufs=[];pickBuf;pickParams;pickOut;statsBuf;pipes={};bg={};bgl={};colorTex;idTex;depthTex;width=0;height=0;usedInstPages=0;usedLipPages=0;usedFreePages=0;pickPending=null;pickBusy=!1;statsBusy=!1;stats={l0:0,l1:0,l2:0,lip:0,fre:0,ninst:0};lastStatsT=0;organelleMeshes=[];onLost=()=>{};constructor(e){this.canvas=e}async init(e,t){if(!("gpu"in navigator))throw new Error("WebGPU is not available in this browser.");const i=await navigator.gpu.requestAdapter({powerPreference:"high-performance"});if(!i)throw new Error("No WebGPU adapter found (GPU blocked or unsupported).");try{const m=i.info||{};this.adapterInfo=[m.vendor,m.architecture,m.device,m.description].filter(Boolean).join(" ")}catch{}const n=i.limits,r={maxStorageBufferBindingSize:Math.min(n.maxStorageBufferBindingSize,2**31-4),maxBufferSize:Math.min(n.maxBufferSize,2**31-4),maxComputeWorkgroupsPerDimension:n.maxComputeWorkgroupsPerDimension};this.device=await i.requestDevice({requiredLimits:r}),this.device.lost.then(m=>this.onLost(`GPU device lost: ${m.message}`)),this.device.addEventListener("uncapturederror",m=>console.error("WebGPU error:",m.error?.message)),this.limits={...r};const o=this.device.limits.maxStorageBufferBindingSize,c=m=>Math.min(m,Math.floor(o/16));this.instPages=Math.floor(c(8*1024*1024)/H),this.lipPages=Math.floor(c(8*1024*1024)/Y),this.freePages=Math.floor(c(10*1024*1024)/Y),this.capL0=Math.min(4*1024*1024,Math.floor(o/8)-2*1024*1024),this.capL1=2*1024*1024,this.capL2=Math.min(8*1024*1024,Math.floor(o/4)),this.capLip=this.lipPages*Y,this.capFree=this.freePages*Y,this.ctx=this.canvas.getContext("webgpu"),this.format=navigator.gpu.getPreferredCanvasFormat(),this.ctx.configure({device:this.device,format:this.format,alphaMode:"opaque",usage:GPUTextureUsage.RENDER_ATTACHMENT|GPUTextureUsage.COPY_SRC});const d=this.device,a=(m,l,h)=>d.createBuffer({size:Math.max(16,Math.ceil(m/16)*16),usage:l,label:h}),u=E.STORAGE|E.COPY_DST;this.uni=a(512,E.UNIFORM|E.COPY_DST,"uniforms"),this.inst=a(this.instPages*H*16,u,"instances"),this.pageOrg=a(this.instPages*16,u,"pageOrg"),this.typesBuf=a(e.typeInfo.byteLength,u,"types"),this.typeFlags=a(e.types.length*4,u,"typeFlags"),this.palette=a(1536,u,"palette"),this.beads=a(e.beads.byteLength,u,"beads"),this.st=a(e.st.byteLength,u,"spheretree"),this.clusters=a(e.clusters.byteLength,u,"clusters"),d.queue.writeBuffer(this.typesBuf,0,e.typeInfo),d.queue.writeBuffer(this.palette,0,ni()),d.queue.writeBuffer(this.beads,0,e.beads),d.queue.writeBuffer(this.st,0,e.st),d.queue.writeBuffer(this.clusters,0,e.clusters),this.counters=a(32,E.STORAGE|E.COPY_DST|E.COPY_SRC,"counters"),this.list01=a((this.capL0+this.capL1)*8,E.STORAGE,"list01"),this.listL2=a(this.capL2*4,E.STORAGE,"listL2"),this.args=a(80,E.STORAGE|E.INDIRECT,"args"),this.cullParams=a(32,E.UNIFORM|E.COPY_DST,"cullParams"),this.caps=a(32,E.UNIFORM|E.COPY_DST,"caps"),d.queue.writeBuffer(this.caps,0,new Uint32Array([this.capL0,this.capL1,this.capL2,this.capLip,this.capFree,0,0,0])),this.lipPool=a(this.lipPages*Y*16,u,"lipidPool"),this.lipPagesBuf=a(this.lipPages*16,u,"lipidPages"),this.lipList=a(this.capLip*4,E.STORAGE,"lipidList"),this.lipParams=a(16,E.UNIFORM|E.COPY_DST,"lipParams"),this.freePool=a(this.freePages*Y*16,u,"freePool"),this.freePagesBuf=a(this.freePages*16,u,"freePages"),this.freeList=a(this.capFree*4,E.STORAGE,"freeList"),this.freeParams=a(16,E.UNIFORM|E.COPY_DST,"freeParams"),this.meshUni=a(256*(this.maxSlots+16),E.UNIFORM|E.COPY_DST,"meshUniforms"),this.tileState=a(484*22*4,u,"tileState"),this.pickBuf=d.createBuffer({size:16,usage:E.MAP_READ|E.COPY_DST}),this.pickParams=a(16,E.UNIFORM|E.COPY_DST,"pickParams"),this.pickOut=a(16,E.STORAGE|E.COPY_SRC,"pickOut"),this.statsBuf=d.createBuffer({size:32,usage:E.MAP_READ|E.COPY_DST}),this.initFill(t),await this.initPipelines(),this.resize()}initFill(e){const t=this.device,i=e.PT.length,n=new ArrayBuffer(i*48),r=new Float32Array(n),o=new Uint32Array(n);for(let p=0;p<i;p++){for(let b=0;b<8;b++)r[12*p+b]=e.PQ[8*p+b];o[12*p+8]=e.PT[p],o[12*p+9]=e.PC[p]}const c=e.pores.length/3,d=new Float32Array(Math.max(1,c)*4);for(let p=0;p<c;p++)d[4*p]=e.pores[3*p],d[4*p+1]=e.pores[3*p+1],d[4*p+2]=e.pores[3*p+2];const a=e.grid,u=a.prim_idx.length,m=new Uint32Array(a.prim_off.length+a.pore_off.length);m.set(a.prim_off,0);for(let p=0;p<a.pore_off.length;p++)m[a.prim_off.length+p]=a.pore_off[p]+u;const l=new Uint32Array(u+a.pore_idx.length+1);l.set(a.prim_idx,0),l.set(a.pore_idx,u);const h=e.G,f=e.consts,M=new Float32Array([h[0],h[1],h[2],h[3],h[4],h[5],h[6],h[7],h[8],h[9],h[10],0,f.PORE_R,f.PERI_W,f.IMS_W,c,a.n,a.ext,0,0]),v=(p,b)=>{const P=t.createBuffer({size:Math.ceil(p.byteLength/16)*16+16,usage:b|E.COPY_DST});return t.queue.writeBuffer(P,0,p),P};this.fillBufs=[v(n,E.STORAGE),v(d,E.STORAGE),v(m,E.STORAGE),v(l,E.STORAGE),v(M,E.UNIFORM)]}async initPipelines(){const e=this.device,t=(G,j)=>({binding:G,visibility:j,buffer:{type:"read-only-storage"}}),i=G=>({binding:G,visibility:O.COMPUTE,buffer:{type:"storage"}}),n=(G,j,ne=!1)=>({binding:G,visibility:j,buffer:{type:"uniform",hasDynamicOffset:ne}}),r=O.COMPUTE|O.VERTEX|O.FRAGMENT,o=O.COMPUTE|O.VERTEX,c=O.VERTEX|O.FRAGMENT,d=G=>e.createBindGroupLayout({entries:G});this.bgl.g0=d([n(0,r),t(1,o),t(2,o),t(3,o),t(4,o),t(5,c)]),this.bgl.cull=d([i(0),i(1),i(2),n(3,O.COMPUTE),t(4,O.COMPUTE)]),this.bgl.cullS=d([i(0),t(1,O.COMPUTE),t(2,O.COMPUTE),i(3),n(4,O.COMPUTE)]),this.bgl.fin=d([t(0,O.COMPUTE),i(1),n(2,O.COMPUTE)]),this.bgl.l0=d([t(0,O.VERTEX),t(1,O.VERTEX)]),this.bgl.l2=d([t(0,O.VERTEX)]),this.bgl.sph=d([t(0,O.VERTEX),t(1,O.VERTEX),t(2,O.VERTEX)]),this.bgl.mesh=d([n(0,c,!0)]),this.bgl.fill=d([t(0,O.FRAGMENT),t(1,O.FRAGMENT),t(2,O.FRAGMENT),t(3,O.FRAGMENT),n(4,O.FRAGMENT),t(5,O.FRAGMENT)]),this.bgl.comp=d([{binding:0,visibility:O.FRAGMENT,texture:{sampleType:"unfilterable-float"}},{binding:1,visibility:O.FRAGMENT,texture:{sampleType:"uint"}},{binding:2,visibility:O.FRAGMENT,texture:{sampleType:"depth"}}]),this.bgl.pick=d([{binding:0,visibility:O.COMPUTE,texture:{sampleType:"uint"}},{binding:1,visibility:O.COMPUTE,texture:{sampleType:"depth"}},n(2,O.COMPUTE),i(3)]);const a=(...G)=>e.createPipelineLayout({bindGroupLayouts:G}),u=async(G,j)=>{const ne=e.createShaderModule({code:G,label:j}),ve=await ne.getCompilationInfo(),te=ve.messages.filter(we=>we.type==="error");for(const we of ve.messages)(we.type==="error"?console.error:console.warn)(`[${j}] ${we.lineNum}:${we.linePos} ${we.message}`);if(te.length)throw new Error(`WGSL compile error in ${j}: ${te[0].lineNum}:${te[0].linePos} ${te[0].message}`);return ne},m=await u(Gt,"cull"),l=await u(Nt,"cullSpheres"),h=await u($t,"finalize"),f=await u(Vt,"drawL0"),M=await u(Wt,"drawL1"),v=await u(jt,"drawL2"),p=await u(Ht,"drawSpheres"),b=await u(Yt,"mesh"),P=await u(Xt,"fill"),x=await u(Kt,"composite"),A=await u(Qt,"pick"),k=(G,j,ne)=>e.createComputePipelineAsync({layout:ne,compute:{module:G,entryPoint:j}}),F=[{format:"rgba8unorm"},{format:"r32uint"}],C=(G="greater")=>({format:"depth32float",depthWriteEnabled:!0,depthCompare:G}),B=(G,j,ne,ve,te={})=>e.createRenderPipelineAsync({layout:ve,vertex:{module:G,entryPoint:j,...te.vertex?{buffers:te.vertex.buffers}:{}},fragment:{module:G,entryPoint:ne,targets:F},primitive:{topology:"triangle-list",cullMode:"none"},depthStencil:te.depthStencil??C()}),R=this.bgl.g0,[q,Q,K,V,D,Le,Ee,_e,ze,Re,Te,ce]=await Promise.all([k(m,"cullInstances",a(R,this.bgl.cull)),k(l,"cullSpheres",a(R,this.bgl.cullS)),k(h,"main",a(this.bgl.fin)),B(f,"vsL0","fsMain",a(R,this.bgl.l0)),B(M,"vsL1","fsMain",a(R,this.bgl.l0)),B(v,"vsL2","fsMain",a(R,this.bgl.l2)),B(p,"vsLipid","fsMain",a(R,this.bgl.sph)),B(p,"vsFree","fsMain",a(R,this.bgl.sph)),B(b,"vsMesh","fsMesh",a(R,this.bgl.mesh),{vertex:{buffers:[{arrayStride:32,attributes:[{shaderLocation:0,offset:0,format:"float32x4"},{shaderLocation:1,offset:16,format:"float32x4"}]}]}}),B(P,"vsFill","fsFill",a(R,this.bgl.fill),{depthStencil:C("always")}),e.createRenderPipelineAsync({layout:a(R,this.bgl.comp),vertex:{module:x,entryPoint:"vsFull"},fragment:{module:x,entryPoint:"fsComposite",targets:[{format:this.format}]},primitive:{topology:"triangle-list"}}),k(A,"main",a(this.bgl.pick))]);Object.assign(this.pipes,{cull:q,cullS:Q,fin:K,l0:V,l1:D,l2:Le,lip:Ee,fre:_e,mesh:ze,fill:Re,comp:Te,pick:ce});const W=(G,j,ne=[])=>e.createBindGroup({layout:G,entries:j.map((ve,te)=>({binding:te,resource:{buffer:ve,size:ne[te]}}))});this.bg.g0=W(R,[this.uni,this.inst,this.pageOrg,this.typesBuf,this.typeFlags,this.palette]),this.bg.cull=W(this.bgl.cull,[this.counters,this.list01,this.listL2,this.cullParams,this.clusters]),this.bg.cullLip=W(this.bgl.cullS,[this.counters,this.lipPool,this.lipPagesBuf,this.lipList,this.lipParams]),this.bg.cullFree=W(this.bgl.cullS,[this.counters,this.freePool,this.freePagesBuf,this.freeList,this.freeParams]),this.bg.fin=W(this.bgl.fin,[this.counters,this.args,this.caps]),this.bg.l0=W(this.bgl.l0,[this.beads,this.list01]),this.bg.l1=W(this.bgl.l0,[this.st,this.list01]),this.bg.l2=W(this.bgl.l2,[this.listL2]),this.bg.lip=W(this.bgl.sph,[this.lipPool,this.lipPagesBuf,this.lipList]),this.bg.fre=W(this.bgl.sph,[this.freePool,this.freePagesBuf,this.freeList]),this.bg.mesh=e.createBindGroup({layout:this.bgl.mesh,entries:[{binding:0,resource:{buffer:this.meshUni,size:16}}]});const se=this.fillBufs;this.bg.fill=W(this.bgl.fill,[se[0],se[1],se[2],se[3],se[4],this.tileState])}resize(e=1){const t=Math.min(window.devicePixelRatio||1,2)*e,i=Math.max(1,Math.floor(this.canvas.clientWidth*t)),n=Math.max(1,Math.floor(this.canvas.clientHeight*t));if(i===this.width&&n===this.height&&this.colorTex)return;this.canvas.width=i,this.canvas.height=n,this.width=i,this.height=n,this.colorTex?.destroy(),this.idTex?.destroy(),this.depthTex?.destroy();const r=GPUTextureUsage;this.colorTex=this.device.createTexture({size:[i,n],format:"rgba8unorm",usage:r.RENDER_ATTACHMENT|r.TEXTURE_BINDING}),this.idTex=this.device.createTexture({size:[i,n],format:"r32uint",usage:r.RENDER_ATTACHMENT|r.TEXTURE_BINDING|r.COPY_SRC}),this.depthTex=this.device.createTexture({size:[i,n],format:"depth32float",usage:r.RENDER_ATTACHMENT|r.TEXTURE_BINDING|r.COPY_SRC}),this.bg.comp=this.device.createBindGroup({layout:this.bgl.comp,entries:[{binding:0,resource:this.colorTex.createView()},{binding:1,resource:this.idTex.createView()},{binding:2,resource:this.depthTex.createView()}]}),this.bg.pick=this.device.createBindGroup({layout:this.bgl.pick,entries:[{binding:0,resource:this.idTex.createView()},{binding:1,resource:this.depthTex.createView()},{binding:2,resource:{buffer:this.pickParams}},{binding:3,resource:{buffer:this.pickOut}}]})}writeUniforms(e,t,i){const n=this.width/this.height,r=e.viewRot(),o=e.proj(n),c=Jt(o,r),d=Zt(c),a=e.eye,u=new ArrayBuffer(512),m=new Float32Array(u),l=new Uint32Array(u);m.set(o,0),m.set(r,16),m.set(d,32);const h=e.pxPerNmAtUnit(this.height);m.set([a[0],a[1],a[2],h],48);const f=[t.plane[0],t.plane[1],t.plane[2]];m.set([f[0],f[1],f[2],t.plane[3]],52);const M=[r[0]*f[0]+r[4]*f[1]+r[8]*f[2],r[1]*f[0]+r[5]*f[1]+r[9]*f[2],r[2]*f[0]+r[6]*f[1]+r[10]*f[2]];m.set([M[0],M[1],M[2],$(f,a)-t.plane[3]],56),m.set([t.cutOn?1:0,t.slab,Math.min(400,Math.max(20,t.slab)),performance.now()/1e3],60),m.set([t.l0px,t.l1px,t.minpx,.6],64),m.set([this.width,this.height,1/this.width,1/this.height],68),m.set(e.frustum(n),72),l.set([t.colorMode,t.regionMask>>>0,t.freeHide,t.highlight?1:0],96),l.set([t.selId>>>0,0,t.shaded?1:0,t.organelle?1:0],100),m.set([360/65535,.8,t.outline,Math.min(window.devicePixelRatio||1,2)*t.renderScale],104),m.set([i[0],i[1],i[2],e.near()],108),l.set([this.capL0,0,0,0],112),this.device.queue.writeBuffer(this.uni,0,u,0,464)}frame(e,t,i,n,r){const o=this.device;this.writeUniforms(e,t,i);const c=o.createCommandEncoder();c.clearBuffer(this.counters);const d=!t.organelle;if(d){const f=c.beginComputePass();f.setBindGroup(0,this.bg.g0);const M=(v,p,b)=>{const P=Math.ceil(v/64),x=Math.min(P,65535),A=Math.ceil(P/x);return o.queue.writeBuffer(p,0,b(x*64)),[x,A]};if(this.usedInstPages>0){const[v,p]=M(this.usedInstPages*H,this.cullParams,b=>new Uint32Array([this.usedInstPages,b,this.capL0,this.capL1,this.capL2,this.capLip,this.capFree,0]).buffer);f.setPipeline(this.pipes.cull),f.setBindGroup(1,this.bg.cull),f.dispatchWorkgroups(v,p)}if(this.usedLipPages>0&&t.showLipids){const[v,p]=M(this.usedLipPages*Y,this.lipParams,b=>new Uint32Array([this.usedLipPages,b,this.capLip,1]).buffer);f.setPipeline(this.pipes.cullS),f.setBindGroup(1,this.bg.cullLip),f.dispatchWorkgroups(v,p)}if(this.usedFreePages>0){const[v,p]=M(this.usedFreePages*Y,this.freeParams,b=>new Uint32Array([this.usedFreePages,b,this.capFree,0]).buffer);f.setPipeline(this.pipes.cullS),f.setBindGroup(1,this.bg.cullFree),f.dispatchWorkgroups(v,p)}f.setPipeline(this.pipes.fin),f.setBindGroup(0,this.bg.fin),f.dispatchWorkgroups(1),f.end()}const a=c.beginRenderPass({colorAttachments:[{view:this.colorTex.createView(),loadOp:"clear",storeOp:"store",clearValue:{r:.965,g:.945,b:.906,a:1}},{view:this.idTex.createView(),loadOp:"clear",storeOp:"store",clearValue:{r:0,g:0,b:0,a:0}}],depthStencilAttachment:{view:this.depthTex.createView(),depthLoadOp:"clear",depthStoreOp:"store",depthClearValue:0}});if(a.setBindGroup(0,this.bg.g0),a.setPipeline(this.pipes.fill),a.setBindGroup(1,this.bg.fill),a.draw(3),t.organelle){a.setPipeline(this.pipes.mesh);const f=this.maxSlots;o.queue.writeBuffer(this.meshUni,f*256,new Uint32Array([0,1,0,0]));for(const M of this.organelleMeshes)a.setBindGroup(1,this.bg.mesh,[f*256]),a.setVertexBuffer(0,M.vb),a.setIndexBuffer(M.ib,"uint32"),a.drawIndexed(M.nIdx,1)}else if(t.showMeshes&&n.length){a.setPipeline(this.pipes.mesh);for(const f of n)a.setBindGroup(1,this.bg.mesh,[f.slot*256]),a.setVertexBuffer(0,f.vb),a.setIndexBuffer(f.ib,"uint32"),a.drawIndexed(f.nIdx,2)}d&&(a.setPipeline(this.pipes.l2),a.setBindGroup(1,this.bg.l2),a.drawIndirect(this.args,32),a.setPipeline(this.pipes.l1),a.setBindGroup(1,this.bg.l1),a.drawIndirect(this.args,16),a.setPipeline(this.pipes.l0),a.setBindGroup(1,this.bg.l0),a.drawIndirect(this.args,0),t.showLipids&&(a.setPipeline(this.pipes.lip),a.setBindGroup(1,this.bg.lip),a.drawIndirect(this.args,48)),a.setPipeline(this.pipes.fre),a.setBindGroup(1,this.bg.fre),a.drawIndirect(this.args,64)),a.end();const u=this.ctx.getCurrentTexture(),m=c.beginRenderPass({colorAttachments:[{view:u.createView(),loadOp:"clear",storeOp:"store",clearValue:{r:1,g:1,b:1,a:1}}]});m.setPipeline(this.pipes.comp),m.setBindGroup(0,this.bg.g0),m.setBindGroup(1,this.bg.comp),m.draw(3),m.end();let l=null;if(this.pickPending&&!this.pickBusy){l=this.pickPending,this.pickPending=null;const f=Math.max(0,Math.min(this.width-1,Math.floor(l.x))),M=Math.max(0,Math.min(this.height-1,Math.floor(l.y)));o.queue.writeBuffer(this.pickParams,0,new Uint32Array([f,M,0,0]));const v=c.beginComputePass();v.setPipeline(this.pipes.pick),v.setBindGroup(0,this.bg.pick),v.dispatchWorkgroups(1),v.end(),c.copyBufferToBuffer(this.pickOut,0,this.pickBuf,0,16)}let h=!1;if(r&&!this.statsBusy&&d&&(c.copyBufferToBuffer(this.counters,0,this.statsBuf,0,32),h=!0),o.queue.submit([c.finish()]),l){this.pickBusy=!0;const f=l;this.pickBuf.mapAsync(GPUMapMode.READ).then(()=>{const M=this.pickBuf.getMappedRange(0,8),v=new Uint32Array(M)[0],p=new Float32Array(M)[1];this.pickBuf.unmap(),this.pickBusy=!1,f.resolve({id:v>>>0,depth:p})}).catch(()=>{this.pickBusy=!1})}h&&(this.statsBusy=!0,this.statsBuf.mapAsync(GPUMapMode.READ).then(()=>{const f=new Uint32Array(this.statsBuf.getMappedRange());this.stats={l0:f[0],l1:f[1],l2:f[2],lip:f[3],fre:f[4],ninst:f[5]},this.statsBuf.unmap(),this.statsBusy=!1}).catch(()=>{this.statsBusy=!1}))}pick(e,t){return new Promise(i=>{this.pickPending={x:e,y:t,resolve:i}})}setMeshSlot(e,t){const i=new ArrayBuffer(16),n=new Uint32Array(i),r=new Float32Array(i);n[0]=e,n[1]=t?2:0,r[2]=t?-1:-1.9,r[3]=t?1:1.9,this.device.queue.writeBuffer(this.meshUni,e*256,i)}uploadOrganelleMeshes(e){const t=this.device,i=new Uint32Array(e,0,1)[0],n=new Uint32Array(e,4,i*3);let r=4+i*12;for(const o of this.organelleMeshes)o.vb.destroy(),o.ib.destroy();this.organelleMeshes=[];for(let o=0;o<i;o++){const c=n[3*o],d=n[3*o+1],a=n[3*o+2],u=new Float32Array(e.slice(r,r+d*12));r+=d*12;const m=new Uint32Array(e.slice(r,r+a*12));r+=a*12;const l=new Float32Array(d*3);for(let p=0;p<a;p++){const b=m[3*p],P=m[3*p+1],x=m[3*p+2],A=u[3*P]-u[3*b],k=u[3*P+1]-u[3*b+1],F=u[3*P+2]-u[3*b+2],C=u[3*x]-u[3*b],B=u[3*x+1]-u[3*b+1],R=u[3*x+2]-u[3*b+2],q=k*R-F*B,Q=F*C-A*R,K=A*B-k*C;for(const V of[b,P,x])l[3*V]+=q,l[3*V+1]+=Q,l[3*V+2]+=K}const h=new Float32Array(d*8),f=new Uint32Array(h.buffer);for(let p=0;p<d;p++){const b=Math.hypot(l[3*p],l[3*p+1],l[3*p+2])||1;h[8*p]=u[3*p],h[8*p+1]=u[3*p+1],h[8*p+2]=u[3*p+2],f[8*p+3]=c,h[8*p+4]=l[3*p]/b,h[8*p+5]=l[3*p+1]/b,h[8*p+6]=l[3*p+2]/b}const M=t.createBuffer({size:h.byteLength,usage:E.VERTEX|E.COPY_DST});t.queue.writeBuffer(M,0,h);const v=t.createBuffer({size:Math.ceil(m.byteLength/4)*4,usage:E.INDEX|E.COPY_DST});t.queue.writeBuffer(v,0,m),this.organelleMeshes.push({vb:M,ib:v,nIdx:a*3,mem:c})}}gpuMemoryMB(){return[this.inst,this.pageOrg,this.beads,this.st,this.clusters,this.list01,this.listL2,this.lipPool,this.lipList,this.freePool,this.freeList,this.meshUni].reduce((t,i)=>t+i.size,0)/1048576+this.width*this.height*12/1048576}}class ci{target=[0,0,0];dist=9e3;fwd=[0,0,-1];up=[0,1,0];fovy=35*Math.PI/180;minDist=1.5;maxDist=14e3;anim=null;get right(){return X(de(this.fwd,this.up))}get eye(){return ae(this.target,L(this.fwd,this.dist))}orthonormalize(){this.fwd=X(this.fwd);const e=X(de(this.fwd,this.up));this.up=X(de(e,this.fwd))}orbit(e,t){const i=this.right;this.fwd=Be(this.fwd,this.up,-e),this.fwd=Be(this.fwd,i,-t),this.up=Be(this.up,i,-t),this.orthonormalize()}roll(e){this.up=Be(this.up,this.fwd,e),this.orthonormalize()}pan(e,t,i){const n=1/i;this.target=_(this.target,_(L(this.right,-e*n),L(this.up,t*n)))}zoom(e){this.dist=Math.min(this.maxDist,Math.max(this.minDist,this.dist*e))}fly(e){this.target=_(this.target,e)}pxPerNmAtUnit(e){return e/(2*Math.tan(this.fovy/2))}distForWidth(e,t){return e/(2*Math.tan(this.fovy/2)*Math.min(t,1))}near(){return Math.min(5,Math.max(.02,this.dist*.004))}viewRot(){const e=this.right,t=this.up,i=this.fwd,n=new Float32Array(16);return n[0]=e[0],n[4]=e[1],n[8]=e[2],n[1]=t[0],n[5]=t[1],n[9]=t[2],n[2]=-i[0],n[6]=-i[1],n[10]=-i[2],n[15]=1,n}proj(e){const t=1/Math.tan(this.fovy/2),i=new Float32Array(16);return i[0]=t/e,i[5]=t,i[11]=-1,i[14]=this.near(),i}frustum(e){const t=Math.tan(this.fovy/2),i=t*e,n=this.right,r=this.up,o=this.fwd,c=(a,u=0)=>{const m=X(a);return[m[0],m[1],m[2],u]},d=[...c(_(n,L(o,i))),...c(_(L(n,-1),L(o,i))),...c(_(r,L(o,t))),...c(_(L(r,-1),L(o,t))),o[0],o[1],o[2],-this.near(),0,0,0,1];return new Float32Array(d)}ray(e,t,i){const n=Math.tan(this.fovy/2),r=n*i,o=(e*2-1)*r,c=(1-t*2)*n;return X(_(this.fwd,_(L(this.right,o),L(this.up,c))))}state(){return{target:[...this.target],dist:this.dist,fwd:[...this.fwd],up:[...this.up]}}set(e){this.target=[...e.target],this.dist=e.dist,this.fwd=[...e.fwd],this.up=[...e.up],this.orthonormalize()}flyTo(e,t=1800,i,n){this.anim={t0:performance.now(),dur:t,a:this.state(),b:e,onDone:i,onStep:n}}step(e){if(!this.anim)return!1;const t=this.anim;let i=Math.min(1,(e-t.t0)/t.dur);const n=ei(i),r=Math.log(t.a.dist),o=Math.log(t.b.dist),c=Math.hypot(...ae(t.b.target,t.a.target)),d=Math.log(Math.max(t.a.dist,t.b.dist,c*1.2)),a=Math.max(0,d-Math.max(r,o)),u=r+(o-r)*n+a*Math.sin(Math.PI*n);this.target=it(t.a.target,t.b.target,n),this.dist=Math.exp(u),this.fwd=ui(t.a.fwd,t.b.fwd,n);let m=t.a.up;return this.up=X(it(m,t.b.up,n)),Math.abs($(this.up,this.fwd))>.98&&(this.up=X(de(de(this.fwd,t.b.up),this.fwd))),this.orthonormalize(),t.onStep?.(n),i>=1&&(this.anim=null,t.onDone?.()),!0}}function ui(s,e,t){const i=Math.max(-1,Math.min(1,$(s,e))),n=Math.acos(i);if(n<1e-4)return X(it(s,e,t));if(n>Math.PI-.001){const o=X(de(s,Math.abs(s[0])<.9?[1,0,0]:[0,1,0]));return Be(s,o,Math.PI*t)}const r=Math.sin(n);return X(_(L(s,Math.sin((1-t)*n)/r),L(e,Math.sin(t*n)/r)))}const At='const vt=.005493247882810712,cn=.005,ln=.01,fn=.70710678;let Ct="",Ot={},kt,bt,Ut,W,jt=null,lt=[],St={glucanPerNm3:.018147,mannanPerProtein:1908.6},zt="";const nt=new Map,hn=48;async function Lt(s){if(Ot[s])return Xt(Ot[s].slice().buffer);const t=await fetch(Ct+s+(zt?`?v=${encodeURIComponent(zt)}`:""));if(!t.ok)throw new Error(`HTTP ${t.status} ${s}`);return Xt(await t.arrayBuffer())}async function Xt(s){const t=new Uint8Array(s,0,Math.min(2,s.byteLength));return t[0]===31&&t[1]===139?await new Response(new Blob([s]).stream().pipeThrough(new DecompressionStream("gzip"))).arrayBuffer():s}async function yn(s){Ct=s.base,lt=s.G;try{const c=await fetch(Ct+"version.json",{cache:"no-store"});c.ok&&(zt=String((await c.json()).v??""))}catch{zt=""}s.wall&&(St=s.wall),s.embedded&&(Ot=s.embedded);const t=JSON.parse(new TextDecoder().decode(await Lt("library/types.json"))),n=t.types,e=n.length,i=t.bead_q;kt=new Uint32Array(e),bt=new Uint32Array(e),Ut=new Float32Array(e);for(let c=0;c<e;c++)kt[c]=n[c].bo,bt[c]=n[c].nb,Ut[c]=n[c].br;const r=await Lt("library/beads.bin"),y=new Int16Array(r),f=r.byteLength/8;W=new Float32Array(f*3);for(let c=0;c<f;c++)W[3*c]=y[4*c]*i,W[3*c+1]=y[4*c+1]*i,W[3*c+2]=y[4*c+2]*i}const Nt=s=>s+3&-4;function Jt(s,t,n,e,i,r,y,f,c,a){const M=new Uint32Array(s,t,n);t+=4*n;const p=new Int32Array(s,t,3*n);t+=12*n;const U=new Int8Array(s,t,3*(e-n));t+=3*(e-n),t=Nt(t);let o=0,u=a;for(let d=0;d<n;d++){let h=p[3*d],g=p[3*d+1],A=p[3*d+2];const w=M[d];for(let I=0;I<w;I++)I>0&&(h+=U[3*o],g+=U[3*o+1],A+=U[3*o+2],o++),f[4*u]=i[0]+h*r,f[4*u+1]=i[1]+g*r,f[4*u+2]=i[2]+A*r,c[4*u+3]=(y|(d&65535)<<4)>>>0,u++}return{off:t,o:u}}function Yt(s){const t=new Uint32Array(s,0,16);if(t[0]!==827605849)throw new Error("bad tile magic");const n=new Int32Array(s,4,3),e=[n[0]*200,n[1]*200,n[2]*200],[i,r,y,f,c,a,M,p,U]=[t[4],t[5],t[6],t[7],t[8],t[9],t[10],t[11],t[12]];let o=64;const u=new Uint32Array(s.slice(o,o+16*i));o+=16*i;const d=new Uint32Array(s.slice(o,o+16*r));o+=16*r;const h=new Float32Array((f+a)*4),g=new Uint32Array(h.buffer);let A=0;if(y){const m=Jt(s,o,y,f,e,cn,0,h,g,A);o=m.off,A=m.o}if(c){const m=Jt(s,o,c,a,e,ln,1,h,g,A);o=m.off,A=m.o}const w=new Float32Array(M*3);if(M){const m=new Uint16Array(s,o,3*M);for(let l=0;l<M;l++)for(let x=0;x<3;x++)w[3*l+x]=e[x]-80+m[3*l+x]*vt;o+=Nt(6*M)}let I=new Uint32Array(p*3);p&&(U&1?I=new Uint32Array(s.slice(o,o+12*p)):I=Uint32Array.from(new Uint16Array(s,o,3*p)),o+=Nt((U&1?4:2)*3*p));const ot=p?new Uint8Array(s.slice(o,o+p)):new Uint8Array(0),b=new Float32Array(M*3),P=new Uint8Array(M);for(let m=0;m<p;m++){const l=I[3*m],x=I[3*m+1],z=I[3*m+2],F=w[3*x]-w[3*l],H=w[3*x+1]-w[3*l+1],v=w[3*x+2]-w[3*l+2],k=w[3*z]-w[3*l],T=w[3*z+1]-w[3*l+1],Y=w[3*z+2]-w[3*l+2],B=H*Y-v*T,K=v*k-F*Y,X=F*T-H*k;for(const O of[l,x,z])b[3*O]+=B,b[3*O+1]+=K,b[3*O+2]+=X,P[O]=ot[m]}for(let m=0;m<M;m++){const l=Math.hypot(b[3*m],b[3*m+1],b[3*m+2])||1;b[3*m]/=l,b[3*m+1]/=l,b[3*m+2]/=l}return{lo:e,inst:u,ghost:d,free:h,verts:w,normals:b,tris:I,fmem:ot,vmem:P,nI:i,nG:r}}const Gt=s=>(s<<16>>16)/32767*fn;function Rt(s,t,n,e){const i=s[4*t],r=s[4*t+1],y=s[4*t+2],f=s[4*t+3];e[0]=n[0]-80+(r&65535)*vt,e[1]=n[1]-80+(r>>>16)*vt,e[2]=n[2]-80+(y&65535)*vt;const c=Gt(y>>>16),a=Gt(f&65535),M=Gt(f>>>16),p=Math.sqrt(Math.max(0,1-c*c-a*a-M*M)),U=i>>>14&3,o=U===0?[p,c,a,M]:U===1?[c,p,a,M]:U===2?[c,a,p,M]:[c,a,M,p],[u,d,h,g]=o;return e[3]=u*u+d*d-h*h-g*g,e[4]=2*(d*h-u*g),e[5]=2*(d*g+u*h),e[6]=2*(d*h+u*g),e[7]=u*u-d*d+h*h-g*g,e[8]=2*(h*g-u*d),e[9]=2*(d*g-u*h),e[10]=2*(h*g+u*d),e[11]=u*u-d*d-h*h+g*g,{type:i&16383,region:i>>>16&31,flags:i>>>21}}class Zt{h;tbl;mask;start;pts;inst;constructor(t,n,e,i){this.h=i;let r=1;for(;r<e*2;)r<<=1;r=Math.max(r,1024),this.tbl=r,this.mask=r-1;const y=new Uint32Array(e),f=new Uint32Array(r+1);for(let a=0;a<e;a++)y[a]=this.k(Math.floor(t[3*a]/i),Math.floor(t[3*a+1]/i),Math.floor(t[3*a+2]/i)),f[y[a]+1]++;for(let a=0;a<r;a++)f[a+1]+=f[a];this.start=f.slice(),this.pts=new Float32Array(e*3),this.inst=new Int32Array(e);const c=f.slice(0,r);for(let a=0;a<e;a++){const M=c[y[a]]++;this.pts[3*M]=t[3*a],this.pts[3*M+1]=t[3*a+1],this.pts[3*M+2]=t[3*a+2],this.inst[M]=n[a]}}k(t,n,e){return(Math.imul(t,73856093)^Math.imul(n,19349663)^Math.imul(e,83492791))&this.mask}near(t,n,e,i,r=-1){const y=this.h,f=Math.floor(t/y),c=Math.floor(n/y),a=Math.floor(e/y),M=i*i;for(let p=-1;p<=1;p++)for(let U=-1;U<=1;U++)for(let o=-1;o<=1;o++){const u=this.k(f+p,c+U,a+o);for(let d=this.start[u],h=this.start[u+1];d<h;d++){const g=this.pts[3*d]-t,A=this.pts[3*d+1]-n,w=this.pts[3*d+2]-e;if(g*g+A*A+w*w<M&&this.inst[d]!==r)return!0}}return!1}}function Wt(s){let t=0;const n=new Float64Array(12);for(const y of s)t+=bt[y.w[4*y.i]&16383];const e=new Float32Array(t*3),i=new Int32Array(t);let r=0;for(const y of s){const{type:f}=Rt(y.w,y.i,y.lo,n),c=kt[f],a=bt[f];for(let M=0;M<a;M++){const p=W[3*(c+M)],U=W[3*(c+M)+1],o=W[3*(c+M)+2];e[3*r]=n[0]+n[3]*p+n[4]*U+n[5]*o,e[3*r+1]=n[1]+n[6]*p+n[7]*U+n[8]*o,e[3*r+2]=n[2]+n[9]*p+n[10]*U+n[11]*o,i[r]=y.id,r++}}return{xyz:e,owner:i,n:t}}function Qt(s){let t=s>>>0;return()=>{t=t+1831565813>>>0;let n=t;return n=Math.imul(n^n>>>15,n|1),n^=n+Math.imul(n^n>>>7,n|61),((n^n>>>14)>>>0)/4294967296}}const $t=(s,t,n)=>{let e=Math.imul(s^2654435769,2246822507)^Math.imul(t+1663821227,3266489909)^Math.imul(n+668265263,374761393);return e^=e>>>16,e=Math.imul(e,2146121005),e^=e>>>15,e=Math.imul(e,2221713035),e^=e>>>16,e>>>0};function R(s){return Math.sqrt(-2*Math.log(s()+1e-12))*Math.cos(6.2831853*s())}function dn(s,t,n){const e=Math.abs(s)+Math.abs(t)+Math.abs(n)||1;let i=s/e,r=t/e;if(n<0){const c=i;i=(1-Math.abs(r))*(c>=0?1:-1),r=(1-Math.abs(c))*(r>=0?1:-1)}const y=Math.max(0,Math.min(255,Math.round((i+1)*127.5))),f=Math.max(0,Math.min(255,Math.round((r+1)*127.5)));return y|f<<8}const un=(()=>{const s=[20.3,14.9,14.3,8.3,1.8,12,28.4],t=s.reduce((e,i)=>e+i);let n=0;return s.map(e=>n+=e/t)})(),wn=.65,mn=.95;function ct(s,t,n,e){return s.n[0]*t+s.n[1]*n+s.n[2]*e-s.d}function Mn(s,t,n,e){const i=[],r=[],y={lipCand:0,lipMasked:0,glucan:0,mannan:0},f=new Float64Array(12),c=(o,u,d,h)=>{if(n.mode==="shell")return Math.hypot(o,u,d)>lt[0]-70-h;const g=ct(n,o,u,d);return g>=n.lo-h&&g<=n.hi+h};if(e.lipids&&t.tris.length&&n.mode==="slab"){const o=[],u=t.inst.length/4,d=t.ghost.length/4;for(let b=0;b<u;b++)t.inst[4*b]>>>21&1&&o.push({w:t.inst,i:b,lo:t.lo,id:b});for(let b=0;b<d;b++)o.push({w:t.ghost,i:b,lo:t.lo,id:u+b});const h=Wt(o),g=new Zt(h.xyz,h.owner,h.n,1),A=t.verts,w=t.normals,I=t.tris,ot=I.length/3;for(let b=0;b<ot;b++){const P=I[3*b],m=I[3*b+1],l=I[3*b+2],x=A[3*P],z=A[3*P+1],F=A[3*P+2],H=A[3*m],v=A[3*m+1],k=A[3*m+2],T=A[3*l],Y=A[3*l+1],B=A[3*l+2],K=ct(n,x,z,F),X=ct(n,H,v,k),O=ct(n,T,Y,B);if(Math.max(K,X,O)<n.lo-3||Math.min(K,X,O)>n.hi+3)continue;const L=H-x,_=v-z,E=k-F,j=T-x,S=Y-z,q=B-F,G=.5*Math.hypot(_*q-E*S,E*j-L*q,L*S-_*j),Q=t.fmem[b],J=Q===8?[1]:[-1,1];for(const $ of J){const C=Qt($t(s,b,$+7)),N=G/wn,rt=Math.floor(N)+(C()<N-Math.floor(N)?1:0);if(!rt)continue;const D=Math.ceil(Math.sqrt(rt)),tt=D*D,at=new Uint16Array(tt);for(let V=0;V<tt;V++)at[V]=V;for(let V=0;V<rt;V++){const ft=V+Math.floor(C()*(tt-V)),Z=at[V];at[V]=at[ft],at[ft]=Z}for(let V=0;V<rt;V++){const ft=at[V],Z=Math.floor(Math.sqrt(ft)),Vt=ft-Z*Z,it=Vt>>1,Dt=(Vt&1)===0;let ht,yt,Tt,It,_t,Et;Dt?(ht=Z,yt=it,Tt=Z+1,It=it,_t=Z+1,Et=it+1):(ht=Z,yt=it,Tt=Z+1,It=it+1,_t=Z,Et=it+1);let dt=C(),ut=C();dt+ut>1&&(dt=1-dt,ut=1-ut);const tn=ht+dt*(Tt-ht)+ut*(_t-ht),Ht=yt+dt*(It-yt)+ut*(Et-yt),et=(tn-Ht)/D,st=Ht/D,qt=1-et-st,nn=x+et*L+st*j,en=z+et*_+st*S,sn=F+et*E+st*q;let wt=qt*w[3*P]+et*w[3*m]+st*w[3*l],mt=qt*w[3*P+1]+et*w[3*m+1]+st*w[3*l+1],Mt=qt*w[3*P+2]+et*w[3*m+2]+st*w[3*l+2];const Ft=Math.hypot(wt,mt,Mt)||1;wt/=Ft,mt/=Ft,Mt/=Ft;const on=C();let pt=0;for(;pt<6&&on>un[pt];)pt++;y.lipCand++;let gt=nn,At=en,xt=sn;Q===8&&(gt-=wt*2.4,At-=mt*2.4,xt-=Mt*2.4);const Pt=ct(n,gt,At,xt);if(Pt<n.lo-2.5||Pt>n.hi+2.5)continue;const rn=pt===5?[1.55,.75]:[1.95,1.25,.55];let Kt=!1;for(const Bt of rn)if(g.near(gt+wt*$*Bt,At+mt*$*Bt,xt+Mt*$*Bt,mn)){Kt=!0;break}if(Kt){y.lipMasked++;continue}const an=(dn(wt,mt,Mt)|pt<<16|($>0?1:0)<<19|(Q&15)<<20)>>>0;i.push(gt,At,xt,an)}}}}if(e.wall&&lt.length>0){const o=lt[0],u=lt[1],d=lt[2],h=t.lo,g=t.inst.length/4,A=[];for(let l=0;l<g;l++)(t.inst[4*l]>>>16&31)===9&&A.push({w:t.inst,i:l,lo:h,id:l});const w=Wt(A),I=new Zt(w.xyz,w.owner,w.n,1.4),ot=(l,x,z)=>l>=h[0]&&l<h[0]+200&&x>=h[1]&&x<h[1]+200&&z>=h[2]&&z<h[2]+200,b=self.__vol?.[s]??0,P=Math.round(b*St.glucanPerNm3);if(P>0){const l=Qt($t(s,1234567,1)),z=Math.ceil(P/48);for(let F=0;F<z;F++){const H=24+Math.floor(l()*49);let v=0,k=0,T=0,Y=!1;for(let q=0;q<40&&!Y;q++){v=h[0]+l()*200,k=h[1]+l()*200,T=h[2]+l()*200;const G=Math.hypot(v,k,T);Y=G>d+2&&G<u-1}if(!Y)continue;let B=Math.hypot(v,k,T),K=v/B,X=k/B,O=T/B,L=R(l),_=R(l),E=R(l),j=L*K+_*X+E*O;L-=j*K,_-=j*X,E-=j*O;let S=Math.hypot(L,_,E)||1;L/=S,_/=S,E/=S;for(let q=0;q<H&&ot(v,k,T);q++){c(v,k,T,1)&&!I.near(v,k,T,1.35)&&(r.push(v,k,T,(2|(F&65535)<<4)>>>0),y.glucan++),L+=.12*R(l),_+=.12*R(l),E+=.12*R(l),B=Math.hypot(v,k,T),K=v/B,X=k/B,O=T/B,j=L*K+_*X+E*O;const G=B,Q=G<d+3&&j<0||G>u-2&&j>0?-1:.6;L-=(1-Q)*j*K,_-=(1-Q)*j*X,E-=(1-Q)*j*O,S=Math.hypot(L,_,E)||1,L/=S,_/=S,E/=S,v+=1.2*L,k+=1.2*_,T+=1.2*E}}}let m=0;for(const l of A){const{type:x}=Rt(l.w,l.i,h,f),z=f[0],F=f[1],H=f[2],v=Math.hypot(z,F,H);if(v<u-5)continue;if(n.mode==="slab"){const L=ct(n,z,F,H);if(L<n.lo-45-Ut[x]||L>n.hi+45+Ut[x]){m+=80;continue}}else if(v<o-120)continue;const k=Qt($t(s,l.i,99)),T=kt[x],Y=bt[x];let B=St.mannanPerProtein;const K=z/v,X=F/v,O=H/v;for(;B>0;){const L=Math.min(Math.ceil(B),12+Math.floor(k()*27));B-=L;const _=T+Math.floor(k()*Y),E=W[3*_],j=W[3*_+1],S=W[3*_+2];let q=z+f[3]*E+f[4]*j+f[5]*S,G=F+f[6]*E+f[7]*j+f[8]*S,Q=H+f[9]*E+f[10]*j+f[11]*S,J=q-z,$=G-F,C=Q-H,N=Math.hypot(J,$,C)||1;J=J/N+K,$=$/N+X,C=C/N+O,N=Math.hypot(J,$,C)||1,J/=N,$/=N,C/=N,q+=1.3*J,G+=1.3*$,Q+=1.3*C;const rt=m++;for(let D=0;D<L;D++){const tt=Math.hypot(q,G,Q);if(tt>o-.6)break;c(q,G,Q,1)&&!I.near(q,G,Q,1.15,D<2?l.i:-1)&&(r.push(q,G,Q,(3|(rt&65535)<<4)>>>0),y.mannan++),J+=.35*R(k)+.18*q/tt,$+=.35*R(k)+.18*G/tt,C+=.35*R(k)+.18*Q/tt,N=Math.hypot(J,$,C)||1,J/=N,$/=N,C/=N,q+=1.05*J,G+=1.05*$,Q+=1.05*C}}}}const a=new Float32Array(i.length),M=new Uint32Array(a.buffer);for(let o=0;o<i.length;o+=4)a[o]=i[o],a[o+1]=i[o+1],a[o+2]=i[o+2],M[o+3]=i[o+3];const p=new Float32Array(r.length),U=new Uint32Array(p.buffer);for(let o=0;o<r.length;o+=4)p[o]=r[o],p[o+1]=r[o+1],p[o+2]=r[o+2],U[o+3]=r[o+3];return{lipids:a,free:p,stats:y}}self.onmessage=async s=>{const t=s.data;try{if(t.cmd==="init"){jt=yn(t),self.__vol=t.volGlucan,await jt,self.postMessage({cmd:"ready"});return}if(await jt,t.cmd==="load"){const n=await Lt(`tiles/${t.ix}_${t.iy}_${t.iz}.bin.gz`),e=Yt(n),i=e.verts.length/3,r=new Float32Array(i*8),y=new Uint32Array(r.buffer);for(let a=0;a<i;a++)r[8*a]=e.verts[3*a],r[8*a+1]=e.verts[3*a+1],r[8*a+2]=e.verts[3*a+2],y[8*a+3]=e.vmem[a],r[8*a+4]=e.normals[3*a],r[8*a+5]=e.normals[3*a+1],r[8*a+6]=e.normals[3*a+2];for(nt.set(t.idx,{lo:e.lo,inst:e.inst,ghost:e.ghost,verts:e.verts,normals:e.normals,tris:e.tris,fmem:e.fmem,vmem:e.vmem});nt.size>hn;)nt.delete(nt.keys().next().value);const f=e.inst.slice(),c=e.tris.slice();self.postMessage({cmd:"loaded",idx:t.idx,inst:f,free:e.free,vb:r,ib:c,lo:e.lo},[f.buffer,e.free.buffer,r.buffer,c.buffer])}else if(t.cmd==="detail"){let n=nt.get(t.idx);if(n)nt.delete(t.idx),nt.set(t.idx,n);else{const r=Yt(await Lt(`tiles/${t.ix}_${t.iy}_${t.iz}.bin.gz`));n={lo:r.lo,inst:r.inst,ghost:r.ghost,verts:r.verts,normals:r.normals,tris:r.tris,fmem:r.fmem,vmem:r.vmem},nt.set(t.idx,n)}const e=performance.now(),i=Mn(t.idx,n,t.range,t.want);self.postMessage({cmd:"detail",idx:t.idx,key:t.key,lipids:i.lipids,free:i.free,stats:i.stats,ms:performance.now()-e},[i.lipids.buffer,i.free.buffer])}}catch(n){self.postMessage({cmd:"error",idx:t.idx,key:t.key,req:t.cmd,error:String(n&&n.stack||n)})}};\n',ht=typeof self<"u"&&self.Blob&&new Blob(["URL.revokeObjectURL(import.meta.url);",At],{type:"text/javascript;charset=utf-8"});function di(s){let e;try{if(e=ht&&(self.URL||self.webkitURL).createObjectURL(ht),!e)throw"";const t=new Worker(e,{type:"module",name:s?.name});return t.addEventListener("error",()=>{(self.URL||self.webkitURL).revokeObjectURL(e)}),t}catch{return new Worker("data:text/javascript;charset=utf-8,"+encodeURIComponent(At),{type:"module",name:s?.name})}}const Qe=200,Ke=173.3,ft=70;class Je{used;n;count=0;high=0;constructor(e){this.n=e,this.used=new Uint8Array(e)}alloc(e){if(e<=0)return[];if(this.n-this.count<e)return null;const t=[];for(let i=0;i<this.n&&t.length<e;i++)this.used[i]||(this.used[i]=1,t.push(i));return this.count+=e,this.high=Math.max(this.high,t[t.length-1]+1),t}free(e){for(const t of e)this.used[t]&&(this.used[t]=0,this.count--);for(;this.high>0&&!this.used[this.high-1];)this.high--}get freeCount(){return this.n-this.count}}class hi{idx;r;lib;G;workers=[];state;loaded=new Map;instA;lipA;freeA;slots=[];pageTile;pageLocal;tileState=new Uint32Array(484*22);desired=[];desiredSet=new Set;inflight=0;maxInflight=12;detInflight=0;maxDetInflight=3;frameNo=0;onChange=()=>{};onTileLoaded=[];bytesLoaded=0;tilesLoadedTotal=0;starved=!1;detailStats={lipids:0,glucan:0,mannan:0};detailEnabled=!0;lastLoadErr="";constructor(e,t,i,n,r,o){this.idx=e,this.r=t,this.lib=i,this.G=n.G,this.state=new Uint8Array(e.n),this.instA=new Je(t.instPages),this.lipA=new Je(t.lipPages),this.freeA=new Je(t.freePages);for(let a=t.maxSlots-1;a>=0;a--)this.slots.push(a);this.pageTile=new Int32Array(t.instPages).fill(-1),this.pageLocal=new Uint32Array(t.instPages);const c=Array.from(e.volGlucan),d=new URL("data/",document.baseURI).href;for(let a=0;a<r;a++){const u=new di;u.onmessage=m=>this.onMessage(m.data),u.postMessage({cmd:"init",base:d,G:n.G,volGlucan:c,embedded:o,wall:{glucanPerNm3:e.wall.glucan_beads_per_nm3,mannanPerProtein:e.wall.mannan_beads_per_protein}}),this.workers.push(u)}}tileCenter(e){return[this.idx.ix[e]*Qe+100,this.idx.iy[e]*Qe+100,this.idx.iz[e]*Qe+100]}update(e,t,i,n,r,o,c){if(this.frameNo++,o)return;const d=e.eye,a=e.frustum(t),u=this.idx.n,m=this.G[0],l=[n[0],n[1],n[2]],h=100*(Math.abs(l[0])+Math.abs(l[1])+Math.abs(l[2])),f=[];for(let b=0;b<u;b++){const P=this.tileCenter(b),x=ae(P,d);let A=!0;for(let k=0;k<5;k++)if(a[4*k]*x[0]+a[4*k+1]*x[1]+a[4*k+2]*x[2]+a[4*k+3]<-243.3){A=!1;break}if(A){if(i){const k=$(l,P)-n[3];if(k+h+ft<0||k-h-ft>r)continue}else{const k=ge(P);if(k<m-150-Ke||k>m+Ke)continue;const F=ae(d,P);if($(P,F)<-k*260)continue}f.push({i:b,d:ge(x)})}}f.sort((b,P)=>b.d-P.d);let M=this.r.instPages*H*.94,v=this.r.maxSlots-8;const p=[];for(const b of f){const P=this.idx.nInst[b.i];if(P>M||v<=0){this.starved=!0;break}M-=Math.ceil(P/H)*H,v--,p.push(b.i)}p.length===f.length&&(this.starved=!1),this.desired=p,this.desiredSet=new Set(p);for(const b of p){const P=this.loaded.get(b);P&&(P.lastUsed=this.frameNo)}for(const b of p){if(this.inflight>=this.maxInflight)break;this.state[b]===0&&(this.state[b]=1,this.inflight++,this.workers[b%this.workers.length].postMessage({cmd:"load",idx:b,ix:this.idx.ix[b],iy:this.idx.iy[b],iz:this.idx.iz[b]}))}this.updateDetail(e,i,n,r,c)}updateDetail(e,t,i,n,r){if(!this.detailEnabled)return;const o=e.eye,c=e.pxPerNmAtUnit(r),d=[i[0],i[1],i[2]],a=Math.min(n,250),u=[];for(const l of this.desired){const h=this.loaded.get(l);if(!h)continue;const f=this.idx.nTri[l]>0,M=this.idx.volGlucan[l]>0||(this.idx.labelMask[l]&2)!==0;if(!f&&!M)continue;const v=this.tileCenter(l),p=Math.max(1,ge(ae(v,o))-Ke);c/p*.55<.8||!t&&!M||u.push({t:h,d:p})}u.sort((l,h)=>l.d-h.d);const m=new Set(u.slice(0,40).map(l=>l.t.idx));for(const l of this.loaded.values())l.det&&!m.has(l.idx)&&!l.detPending&&this.dropDetail(l);for(const{t:l}of u.slice(0,40)){if(this.detInflight>=this.maxDetInflight)break;if(l.detPending)continue;const h=l.det;let f=!1;if(h&&(t?f=!h.shell&&$(h.n,d)>.99985&&h.a<=i[3]+.5&&h.b>=i[3]+a-.5:f=h.shell),f)continue;const M=Math.max(10,.35*a),v=t?{n:d,d:i[3],lo:-M,hi:a+M,mode:"slab"}:{n:d,d:i[3],lo:0,hi:0,mode:"shell"};l.detPending=!0,this.detInflight++;const p=t?`${i[3]-M}|${i[3]+a+M}`:"shell";l._req={n:d,a:i[3]-M,b:i[3]+a+M,shell:!t};const b=this.idx.nTri[l.idx]>0,P=this.idx.volGlucan[l.idx]>0||(this.idx.labelMask[l.idx]&2)!==0;this.workers[l.idx%this.workers.length].postMessage({cmd:"detail",idx:l.idx,key:p,ix:this.idx.ix[l.idx],iy:this.idx.iy[l.idx],iz:this.idx.iz[l.idx],range:v,want:{lipids:b&&t,wall:P}})}}dropDetail(e){this.lipA.free(e.lipPages),this.freeA.free(e.detPages),this.clearSPages(this.r.lipPagesBuf,e.lipPages),this.clearSPages(this.r.freePagesBuf,e.detPages),e.lipPages=[],e.detPages=[],e.nLip=0,e.nDet=0,e.det=null,e.underMesh&&(e.underMesh=!1,this.r.setMeshSlot(e.slot,!1)),this.syncHigh()}clearSPages(e,t){const i=new Uint32Array(4);for(const n of t)this.r.device.queue.writeBuffer(e,n*16,i)}syncHigh(){this.r.usedInstPages=this.instA.high,this.r.usedLipPages=this.lipA.high,this.r.usedFreePages=this.freeA.high}evictOne(e){let t=null;for(const i of this.loaded.values())i.idx===e||this.desiredSet.has(i.idx)||(!t||i.lastUsed<t.lastUsed)&&(t=i);return t?(this.evict(t),!0):!1}evict(e){const t=this.r.device.queue,i=new Float32Array(4);for(const n of e.pages)t.writeBuffer(this.r.pageOrg,n*16,i),this.pageTile[n]=-1;this.instA.free(e.pages),this.freeA.free(e.basePages),this.clearSPages(this.r.freePagesBuf,e.basePages),this.lipA.free(e.lipPages),this.clearSPages(this.r.lipPagesBuf,e.lipPages),this.freeA.free(e.detPages),this.clearSPages(this.r.freePagesBuf,e.detPages),e.vb?.destroy(),e.ib?.destroy(),this.slots.push(e.slot),this.loaded.delete(e.idx),this.state[e.idx]=0,this.tileState[st(this.idx.ix[e.idx],this.idx.iy[e.idx],this.idx.iz[e.idx])]=0,t.writeBuffer(this.r.tileState,0,this.tileState),this.syncHigh()}allocWithEvict(e,t,i){let n=e.alloc(t);for(;!n&&this.evictOne(i);)n=e.alloc(t);return n}onMessage(e){if(e.cmd!=="ready"){if(e.cmd==="error"){if(console.warn("worker error",e.req,e.idx,e.error),this.lastLoadErr=String(e.error).split(`
`)[0],e.req==="load"&&(this.state[e.idx]=0,this.inflight--),e.req==="detail"){this.detInflight--;const t=this.loaded.get(e.idx);t&&(t.detPending=!1)}return}e.cmd==="loaded"&&(this.inflight--,this.onLoaded(e)),e.cmd==="detail"&&(this.detInflight--,this.onDetail(e))}}onLoaded(e){const t=e.idx;if(!this.desiredSet.has(t)){this.state[t]=0,this.onChange();return}const i=e.inst,n=i.length/4,r=this.r.device.queue;this.slots.length||this.evictOne(t);const o=this.allocWithEvict(this.instA,Math.ceil(n/H),t);if(!o||!this.slots.length){this.state[t]=0,this.starved=!0,o&&this.instA.free(o);return}const c=this.slots.pop(),d=e.lo,a=new ArrayBuffer(16),u=new Float32Array(a),m=new Uint32Array(a);u[0]=d[0]-80,u[1]=d[1]-80,u[2]=d[2]-80,o.forEach((x,A)=>{const k=Math.min(H,n-A*H);r.writeBuffer(this.r.inst,x*H*16,i.buffer,A*H*16,k*16),m[3]=k,r.writeBuffer(this.r.pageOrg,x*16,a),this.pageTile[x]=t,this.pageLocal[x]=A*H});const l=e.free,h=l.length/4;let f=[],M=!1;if(h){const x=this.allocWithEvict(this.freeA,Math.ceil(h/Y),t);x?(f=x,this.writeSpheres(this.r.freePool,this.r.freePagesBuf,x,l,c)):M=!0}let v=null,p=null;const b=e.ib.length;b&&(v=this.r.device.createBuffer({size:e.vb.byteLength,usage:GPUBufferUsage.VERTEX|GPUBufferUsage.COPY_DST}),r.writeBuffer(v,0,e.vb),p=this.r.device.createBuffer({size:Math.ceil(e.ib.byteLength/4)*4,usage:GPUBufferUsage.INDEX|GPUBufferUsage.COPY_DST}),r.writeBuffer(p,0,e.ib),this.r.setMeshSlot(c,!1));const P={idx:t,slot:c,lo:d,inst:i,nInst:n,pages:o,basePages:f,nBase:h,lipPages:[],detPages:[],nLip:0,nDet:0,vb:v,ib:p,nIdx:b,lastUsed:this.frameNo,det:null,detPending:!1,underMesh:!1,baseFreeMissing:M};this.loaded.set(t,P),this.state[t]=2,this.tileState[st(this.idx.ix[t],this.idx.iy[t],this.idx.iz[t])]=1,r.writeBuffer(this.r.tileState,0,this.tileState),this.bytesLoaded+=this.idx.gz[t],this.tilesLoadedTotal++,this.syncHigh();for(const x of this.onTileLoaded.slice())x(P);this.onChange()}writeSpheres(e,t,i,n,r){const o=this.r.device.queue,c=n.length/4;i.forEach((d,a)=>{const u=Math.min(Y,c-a*Y);o.writeBuffer(e,d*Y*16,n.buffer,n.byteOffset+a*Y*16,u*16),o.writeBuffer(t,d*16,new Uint32Array([u,r,1,0]))})}onDetail(e){const t=this.loaded.get(e.idx);if(!t)return;t.detPending=!1;const i=t._req;this.lipA.free(t.lipPages),this.clearSPages(this.r.lipPagesBuf,t.lipPages),this.freeA.free(t.detPages),this.clearSPages(this.r.freePagesBuf,t.detPages),t.lipPages=[],t.detPages=[];const n=e.lipids,r=e.free,o=n.length/4,c=r.length/4;if(o){const a=this.lipA.alloc(Math.ceil(o/Y))??this.freeDetailSpace(this.lipA,Math.ceil(o/Y),t.idx);a&&(t.lipPages=a,this.writeSpheres(this.r.lipPool,this.r.lipPagesBuf,a,n,t.slot))}if(c){const a=this.freeA.alloc(Math.ceil(c/Y))??this.freeDetailSpace(this.freeA,Math.ceil(c/Y),t.idx);a&&(t.detPages=a,this.writeSpheres(this.r.freePool,this.r.freePagesBuf,a,r,t.slot))}t.nLip=t.lipPages.length?o:0,t.nDet=t.detPages.length?c:0,t.det={n:i.n,a:i.a,b:i.b,shell:i.shell};const d=t.nLip>0;d!==t.underMesh&&t.nIdx&&(t.underMesh=d,this.r.setMeshSlot(t.slot,d)),this.detailStats.lipids+=e.stats.lipCand,this.syncHigh(),this.onChange()}freeDetailSpace(e,t,i){const n=[...this.loaded.values()].filter(r=>r.idx!==i&&(r.lipPages.length||r.detPages.length)&&!this.desiredSet.has(r.idx));for(const r of n){this.dropDetail(r);const o=e.alloc(t);if(o)return o}return null}meshDraws(){const e=[];for(const t of this.loaded.values())t.nIdx&&t.vb&&t.ib&&e.push({slot:t.slot,vb:t.vb,ib:t.ib,nIdx:t.nIdx,under:t.underMesh});return e}instanceAt(e){const t=Math.floor(e/H),i=this.pageTile[t];if(i<0)return null;const n=this.loaded.get(i);if(!n)return null;const r=this.pageLocal[t]+e%H;return{t:n,local:r,...this.decode(n,r)}}decode(e,t){const i=e.inst,n=i[4*t],r=i[4*t+1],o=i[4*t+2],c=360/65535,d=[e.lo[0]-80+(r&65535)*c,e.lo[1]-80+(r>>>16)*c,e.lo[2]-80+(o&65535)*c];return{type:n&16383,region:n>>>16&31,flags:n>>>21,p:d}}poolIndex(e,t){return e.pages[Math.floor(t/H)]*H+t%H}tileBySlot(e){for(const t of this.loaded.values())if(t.slot===e)return t;return null}countLoaded(e){let t=0;for(const i of this.loaded.values()){const n=i.inst;for(let r=0;r<i.nInst;r++)(n[4*r]&16383)===e&&t++}return t}nearestInTile(e,t,i){let n=-1,r=1e18;for(let o=0;o<e.nInst;o++){if((e.inst[4*o]&16383)!==t)continue;const c=this.decode(e,o),d=(c.p[0]-i[0])**2+(c.p[1]-i[1])**2+(c.p[2]-i[2])**2;d<r&&(r=d,n=o)}return n}totals(){let e=0,t=0,i=0,n=0;for(const r of this.loaded.values())e+=r.nInst,t+=r.nLip,i+=r.nBase+r.nDet,n+=r.nIdx/3;return{tiles:this.loaded.size,inst:e,lip:t,free:i,mesh:n,inflight:this.inflight+this.detInflight,desired:this.desired.length}}pendingForView(){let e=0;for(const i of this.desired)this.state[i]!==2&&e++;let t=0;for(const i of this.loaded.values())i.detPending&&t++;return e+t}}class fi{G;PT;PQ;PC;pores;c;grid;n3;constructor(e){this.G=e.G,this.PT=e.PT,this.PQ=e.PQ,this.PC=e.PC,this.pores=e.pores,this.c=e.consts,this.grid=e.grid,this.n3=e.grid.n**3}cell(e,t,i){const n=this.grid.n,r=this.grid.ext,o=2*r/n,c=d=>Math.min(n-1,Math.max(0,Math.floor((d+r)/o)));return(c(e)*n+c(t))*n+c(i)}sdfPrim(e,t,i,n){const r=this.PT[e],o=this.PQ,c=8*e;if(r===0)return Math.hypot(t-o[c],i-o[c+1],n-o[c+2])-o[c+3];if(r===1){const h=t-o[c],f=i-o[c+1],M=n-o[c+2],v=o[c+3]-o[c],p=o[c+4]-o[c+1],b=o[c+5]-o[c+2],P=Math.min(1,Math.max(0,(h*v+f*p+M*b)/(v*v+p*p+b*b)));return Math.hypot(h-v*P,f-p*P,M-b*P)-o[c+6]}if(r===2){const h=t-o[c],f=i-o[c+1],M=n-o[c+2],v=o[c+3],p=o[c+4],b=o[c+5],P=h*v+f*p+M*b,x=Math.hypot(h-P*v,f-P*p,M-P*b),A=Math.max(x-o[c+6],0);return Math.sqrt(A*A+P*P)-o[c+7]}const d=Math.hypot(t,i,n)||1e-9,a=Math.min(1,Math.max(-1,(t*o[c]+i*o[c+1]+n*o[c+2])/d)),u=Math.acos(a),m=Math.max(u-o[c+4],0)*d,l=d-o[c+5];return Math.sqrt(m*m+l*l)-o[c+6]}label(e,t,i){const n=this.G,r=this.c.PERI_W,o=this.c.PORE_R,c=Math.hypot(e,t,i);if(c>n[0])return 0;if(c>n[1])return 1;if(c>n[2])return 2;const d=this.cell(e,t,i),a=this.grid,u=[1e9,1e9,1e9,1e9,1e9,1e9,1e9,1e9,1e9];let m=1e9,l=1e9;for(let k=a.prim_off[d];k<a.prim_off[d+1];k++){const F=a.prim_idx[k],C=this.sdfPrim(F,e,t,i),B=this.PC[F];B===7?u[1]=Math.min(u[1],C):B===0?u[3]=Math.min(u[3],C):B===4?u[4]=Math.min(u[4],C):B===5?m=Math.min(m,C):B===6?l=Math.min(l,C):B===2?u[6]=Math.min(u[6],C):B===3?u[7]=Math.min(u[7],C):B===1&&(u[8]=Math.min(u[8],C))}const h=Math.max(m,-l);if(u[1]<0)return 4;const f=e-n[3],M=t-n[4],v=i-n[5],p=Math.hypot(f,M,v),b=p-(n[6]-.5*r),P=o+.5*r;let x=!1;for(let k=a.pore_off[d];k<a.pore_off[d+1];k++){const F=a.pore_idx[k],C=this.pores[3*F],B=this.pores[3*F+1],R=this.pores[3*F+2],q=f*C+M*B+v*R;if(q<=0)continue;if(Math.hypot(f-q*C,M-q*B,v-q*R)<P){x=!0;break}}if(Math.abs(b)-.5*r<0&&!x)return 5;if(p<n[6]-.5*r){if(x&&b>-.5*r-8)return 14;const k=e-n[7],F=t-n[8],C=i-n[9];return p<n[6]-r&&k*k+F*F+C*C<n[10]*n[10]?7:6}return p<n[6]?14:u[4]<0?h<0?10:9:u[3]<0?8:u[8]<0?11:u[6]<0?12:u[7]<0?13:3}}const pi=`
<h2>About this model</h2>
<p>This is one complete haploid <i>Saccharomyces cerevisiae</i> cell (strain BY4741, log phase in YPD), built as a
structural whole-cell model in the spirit of Maritan <i>et al.</i> (J. Mol. Biol. 2022; whole <i>Mycoplasma genitalium</i> cell built with cellPACK)
and drawn in the style of David Goodsell. Every macromolecule is drawn from its experimental (PDB) or predicted
(AlphaFold DB) coordinates, at one bead per amino-acid residue or nucleotide when it is large enough on screen.
Further away, molecules switch to a sphere-tree stand-in and then to one sphere each. This only happens once
the detail is smaller than a pixel. Click any molecule to inspect it and to load its full-atom model in Mol*.</p>

<h3>What is in the cell</h3>
<ul>
<li><b>21.6 million macromolecule instances (most multi-chain assemblies), 5,642 molecular species.</b> Protein abundances are from the PaxDb integrated dataset rescaled to
4.2 × 10⁷ protein molecules per cell (Ho <i>et al.</i>), in a 42 µm³ cell. There are 169,000 80S ribosomes, 3.3 million tRNAs, 73,150 nucleosomes and 117 nuclear pore complexes.</li>
<li><b>Structures.</b> 55% of protein mass is drawn from PDB entries. These include the ribosome 4V88, nuclear pore 7N9F,
nucleosome 1ID3, 26S proteasome 4CR2, FAS 2UV8, Pol I/II/III, V-ATPase 6O7V, Pma1 7VH5, CCT, about 100 further Complex Portal
complexes with a yeast PDB assembly (e.g. phosphofructokinase, eIF2/eIF3, the mitochondrial ribosome, mediator, exosome, snoRNPs with
their RNA) and 133 homo-oligomers. Membrane complexes (ATP synthase, respiratory complexes II-IV, Sec translocon, OST, EMC, glucan
synthase, flippases, ...) are oriented in their membrane automatically from UniProt transmembrane and topology annotations.
15% are AlphaFold DB predicted complexes: homodimers and 511 heterodimers (ipTM &ge; 0.5, partners sharing a
compartment; e.g. actin-cofilin, actin-profilin, RNR2-RNR4, SOD1-CCS1, eIF4G-eIF4A); 30% are AlphaFold monomers.
AlphaFold residues with pLDDT &lt; 50 are trimmed. Outer-membrane beta-barrels (TOM40, SAM50, MDM10) are oriented by their strands.</li>
<li><b>Polysomes:</b> 23,000 mRNAs carry 93,000 of the 80S ribosomes (about 4 per mRNA), threaded through the mRNA
channel (path from PDB 7ZW0) and linked by short mRNA segments; the remaining 7,000 mRNAs are free nucleotide-resolution
threads. Linker DNA joins consecutive nucleosomes.</li>
<li><b>ER-bound ribosomes:</b> 4,300 ribosomes sit on translocons in the ER and outer nuclear membrane: on the Ssh1 complex
(yeast ribosome-Ssh1 cryo-EM structure 2WW9) or on Sec61 (7KAH, placed through the shared Sss1 subunit). Their number is
limited by the free translocon pools (80% assumed ribosome-bound); the PM-facing side of the cortical ER is too narrow for ribosomes.</li>
<li><b>Microtubules:</b> 21 microtubules (11,500 &alpha;/&beta;-tubulin dimers) grow from &gamma;-tubulin rings on the spindle pole body,
opposite the nucleolus: 16 short kinetochore microtubules whose plus ends anchor the 16 chromosomes (Rabl configuration),
2 more nuclear ones and 3 cytoplasmic ones about 1 &micro;m long that curve around organelles. The 13-protofilament lattice
and the ring come from the native yeast &gamma;-TuRC-capped microtubule 8QV3; the lattice dimer is 5W3F (8.3 nm repeat).</li>
<li><b>Actin cables:</b> 10 bundles of 8 filaments (33,000 actin subunits, 92 &micro;m) run 90 nm beneath the plasma
membrane, built from the yeast F-actin cryo-EM structure 9GO5 (167&deg; and 2.76 nm per subunit).</li>
<li><b>Membranes</b> (plasma membrane, cortical and perinuclear ER, nuclear envelope with pores, vacuole, a reticular
mitochondrion, Golgi/endosomes, peroxisomes, lipid droplets) are represented by analytic signed-distance fields.
Near the camera every lipid is drawn individually. Class composition (Ejsing <i>et al.</i> 2009): PI 20%, PE 15%, PC 14%, PA 8%, PS 2%,
ergosterol 12%, 0.65 nm² per lipid. Lipids are removed wherever they would overlap a membrane protein.</li>
<li><b>Cell wall:</b> an inner β-1,3-glucan layer and an outer mannoprotein layer with mannan brushes. Bead counts follow a
wall mass of about 20% of dry mass.</li>
</ul>

<h3>How it was built</h3>
<p>Proteins were assigned to compartments from UniProt subcellular-location annotations, with GO as fallback. Membrane proteins were oriented by their
transmembrane segments. Everything was then packed into 200 nm tiles with collision detection on residue-level
sphere trees. Validation on the final model:</p>
<ul>
<li>Bead clash rate 6.8×10⁻⁷.</li>
<li>Containment violations 6.4×10⁻⁵.</li>
<li>Protein densities: cytosol 89 mg/mL, nucleoplasm 62 mg/mL, mitochondrial matrix 117 mg/mL.</li>
<li>99.8% of requested copies were placed outside mitochondria.</li>
</ul>
<p>Positions are stored to 0.0055 nm and orientations to 0.004°.</p>

<h3>Honest caveats</h3>
<ul>
<li><b>Layout.</b> The organelle shapes and positions are a plausible illustration, not a tomogram of one specific cell.</li>
<li><b>Abundances.</b> Absolute abundances are uncertain to about 2-fold.</li>
<li><b>Mitochondria.</b> The mitochondrion saturated: only 69% of matrix and 68% of IMS proteins fit, so their copies are under-represented.</li>
<li><b>Chromatin</b> is a coarse Rabl-like arrangement of nucleosome chains in 59 fragments, each chromosome starting at its kinetochore microtubule. rDNA repeats are not modelled.
Linker DNA and mRNA paths are schematic.</li>
<li><b>Procedural content.</b> Lipids, β-glucan and mannan are generated procedurally in your browser from deterministic seeds. They are not stored
coordinates, and polysaccharides are only partly clash-checked.</li>
<li><b>Spindle pole body.</b> Only the &gamma;-tubulin rings are modelled, not the SPB plaques; PaxDb lists fewer SPC98 copies than 21 rings need.</li>
<li><b>Not included:</b> actin patches, septins, metabolites and water. The FG-repeat mesh of the pore is not modelled.</li>
<li><b>Dropped copies.</b> 9,946 transmembrane-protein copies assigned to soluble compartments were dropped, and 6,888 copies had no structure.</li>
</ul>
<p style="color:#a9a69e">Renderer: WebGPU compute culling and level-of-detail selection, sphere impostors with exact cut caps, analytic
compartment backgrounds and Goodsell-style outlines. Full-atom models are fetched live from RCSB PDB and AlphaFold DB
when you request them.</p>
<p style="color:#a9a69e">Data sources: RCSB Protein Data Bank (CC0) and the individual entries named above; AlphaFold Protein Structure
Database monomers and predicted complexes (CC-BY 4.0, Jumper <i>et al.</i> 2021; Varadi <i>et al.</i> 2024); PaxDb protein abundances;
Complex Portal; UniProt and Gene Ontology annotations; IntAct interactions; Alliance of Genome Resources paralogs;
lipid composition from Ejsing <i>et al.</i> 2009; total protein from Ho <i>et al.</i> 2018.</p>
`,g=s=>document.getElementById(s),re=s=>s>=1e6?(s/1e6).toFixed(s>=1e7?1:2)+" M":s>=1e4?(s/1e3).toFixed(0)+" k":s.toLocaleString("en-US"),kt=s=>s>=1e3?(s/1e3).toFixed(s>=1e4?0:2).replace(/\.?0+$/,"")+" µm":(s>=10?s.toFixed(0):s.toFixed(1))+" nm";let J,Se,xe,Lt,fe,$e,z,I;const y=new ci,S=g("gpu"),w={cutOn:!0,plane:[0,0,1,0],slab:60,colorMode:0,regionMask:16777215,freeHide:0,highlight:!1,selId:0,shaded:!0,organelle:!1,l0px:1.6,l1px:2.5,minpx:.3,outline:1,showMeshes:!0,showLipids:!0,renderScale:1};let pe=!1,T=!0,at=0,Ze=0,Et="auto",le,me=-1,Tt=null,U=null,oe=null;const Z=[];function Ue(s,e){g("loading").querySelector(".m").textContent=s,g("loading").querySelector(".pb div").style.width=`${Math.round(e*100)}%`}function pt(s){const e=g("loading").querySelector(".err");e.style.display="block",e.innerHTML=s,g("loading").querySelector(".t").textContent="Could not start the viewer"}function Me(s,e=2200){const t=g("toast");t.textContent=s,t.style.display="block",clearTimeout(Me._t),Me._t=setTimeout(()=>t.style.display="none",e)}async function mi(){if(!("gpu"in navigator)){pt(`This browser does not expose <b>WebGPU</b>. Use a current Chrome or Edge (version 113 or later), or Safari 26 or later.
      On Linux Chrome, enable <code>chrome://flags/#enable-unsafe-webgpu</code> and <code>#enable-vulkan</code>.
      Firefox 141+ supports WebGPU on Windows.<br><br>Also check that the page is served over <code>http://localhost</code>
      (for example <code>python -m http.server</code>). Opening the file directly with <code>file://</code> does not work.`);return}try{Ue("Loading cell geometry",.01),[xe,fe]=await Promise.all([qe("geometry.json"),qe("bookmarks.json")]),Lt=new fi(xe),Se=await oi(),J=await ri(Ue);const s=await Ie("types_tiles.bin"),e=J.types.length;$e={off:new Uint32Array(s,0,e+1),tiles:new Uint16Array(s,(e+1)*4)},Ue("Starting WebGPU",.7),z=new li(S),z.onLost=o=>{Me(o,1e4)},await z.init(J,xe),le=new Uint32Array(e),Ue("Compiling shaders and starting tile streaming",.9);let t;if(window.__YC_EMBED__){t={};for(const[o,c]of Object.entries(window.__YC_EMBED__))(o.startsWith("tiles/")||o.startsWith("library/"))&&(t[o]=Mt(c))}const i=Math.max(2,Math.min(4,(navigator.hardwareConcurrency||4)-1));I=new hi(Se,z,J,xe,i,t),I.onChange=()=>{T=!0},xi();const n=new URLSearchParams(location.search).get("view"),r=n===null?0:isNaN(Number(n))?Math.max(0,fe.findIndex(o=>o.name.toLowerCase().startsWith(n.toLowerCase()))):Math.min(fe.length-1,Math.max(0,Number(n)));Ae(r,!1),g("loading").style.display="none",requestAnimationFrame(Dt),window.yc={R:z,T:I,cam:y,fs:w,lib:J,idx:Se,applyBookmark:Ae,setType:ke,pickAt:We,bookmarks:fe,setQuality:Ot,describe:lt,ready:!0,state:()=>({pending:I.pendingForView(),detailPending:I.detInflight,loaded:I.loaded.size,totals:I.totals(),stats:z.stats,cam:{t:y.target,d:y.dist},plane:w.plane,gpuMs:zt,gpuBusy:He,dirty:T}),render:()=>{T=!0}}}catch(s){console.error(s),pt(`${String(s?.message||s)}<br><br>See the browser console for details.`)}}function ie(){return[w.plane[0],w.plane[1],w.plane[2]]}function ee(s,e){s=X(s),w.plane=[s[0],s[1],s[2],e],he(),T=!0}function Pe(){const s=y.fwd;ee(s,$(s,y.target))}function gi(s){const e=S.clientWidth/S.clientHeight,t=X(s.normal);let i=X(ae(s.up,L(t,$(s.up,t))));return{target:[...s.target],dist:y.distForWidth(s.width,e),fwd:t,up:i}}function Ae(s,e=!0,t){const i=fe[s];if(!i)return;w.organelle&&s!==0&&Oe(!1);const n=gi(i),r=X(i.normal),o=$(r,i.target),c=ie(),d=w.plane[3],a=w.slab;if(w.cutOn=!0,g("cutOn").checked=!0,pe=!1,g("followView").classList.remove("on"),!e){y.set(n),ee(r,o),w.slab=i.slab,he(),T=!0,t?.();return}y.flyTo(n,2600,()=>{ee(r,o),w.slab=i.slab,he(),t?.()},u=>{let m=X(_(L(c,1-u),L(r,u)));ge(_(L(c,1-u),L(r,u)))<.001&&(m=r),w.plane=[m[0],m[1],m[2],d+(o-d)*u],w.slab=Math.exp(Math.log(a)+(Math.log(i.slab)-Math.log(a))*u),he()}),Ve(i),T=!0}function Ve(s){const e=g("caption");if(!s){e.style.display="none";return}e.innerHTML=`<b>${s.name}</b>${s.caption}`,e.style.display="block",clearTimeout(Ve._t),U||(Ve._t=setTimeout(()=>e.style.display="none",9e3))}function Ct(){U={i:0,phase:"fly",t:performance.now()},g("tourBtn").textContent="■ Stop",g("tourBtn").classList.add("on"),Ae(0,!0,()=>{U&&(U.phase="wait",U.t=performance.now())})}function be(){U=null,g("tourBtn").textContent="▶ Tour",g("tourBtn").classList.remove("on"),Ve(null)}function bi(s){if(U){if(U.phase==="wait"&&(I.pendingForView()===0||s-U.t>14e3))U.phase="hold",U.t=s;else if(U.phase==="hold"&&s-U.t>6500){U.i=(U.i+1)%fe.length,U.phase="fly";const e=U.i;Ae(e,!0,()=>{U&&U.i===e&&(U.phase="wait",U.t=performance.now())})}}}function ke(s){me>=0&&(le[me]&=-2),me=s,s>=0&&(le[s]|=1),w.highlight=s>=0,z.device.queue.writeBuffer(z.typeFlags,0,le),Ft(),T=!0}function Ft(){const s=g("selType");if(me<0){s.style.display="none";return}const e=J.types[me],t=I.countLoaded(me);s.style.display="block",s.innerHTML=`<div><b>${e.gene}</b> ${e.name}</div>
    <div class="sub" style="margin:3px 0 6px">${re(e.copies)} copies in the cell · ${re(t)} in loaded tiles · ${St(e)}</div>
    <div class="grid2"><button class="small" id="flyNearest">Fly to nearest</button><button class="small" id="clearType">Clear highlight</button></div>`,g("flyNearest").onclick=()=>Bt(me),g("clearType").onclick=()=>ke(-1)}function yi(s){if(s=s.trim().toLowerCase(),!s)return[];const e=[];return J.types.forEach((t,i)=>{const n=t.gene.toLowerCase(),r=t.id.toLowerCase(),o=(t.uniprot||"").toLowerCase(),c=(t.name||"").toLowerCase(),d=(t.genes||"").toLowerCase();let a=0;n===s||r===s||o===s?a=100:n.startsWith(s)||r.startsWith(s)?a=60:d.split(/[;,\s]+/).includes(s)?a=55:c.includes(s)?a=30+(c.startsWith(s)?10:0):(d.includes(s)||o.startsWith(s)||(t.pdb||"").toLowerCase()===s)&&(a=20),a&&e.push({i,s:a+Math.log10(t.copies+1)})}),e.sort((t,i)=>i.s-t.s),e.slice(0,14).map(t=>t.i)}function Bt(s){const e=$e.off[s],t=$e.off[s+1];if(t<=e){Me("No copies of this molecule were placed.");return}let i=-1,n=1e18;for(let m=e;m<t;m++){const l=$e.tiles[m],h=I.tileCenter(l),f=ge(ae(h,y.target));f<n&&(n=f,i=l)}const r=I.tileCenter(i),o=J.types[s],c=Math.max(90,o.br*9),d=m=>{const l=I.nearestInTile(m,s,r);if(l<0)return;const h=I.decode(m,l),f=ie(),M=y.state();M.target=h.p,M.dist=y.distForWidth(c,S.clientWidth/S.clientHeight),y.flyTo(M,900,()=>{ee(f,$(f,h.p)-o.br*.4)}),w.plane[3]=$(f,h.p)-o.br*.4;const v=I.poolIndex(m,l);je({kind:"inst",id:(ai|v)>>>0,info:I.instanceAt(v)})};I.loaded.get(i);const a=ie(),u=y.state();u.target=r,u.dist=y.distForWidth(c*2.5,S.clientWidth/S.clientHeight),w.cutOn=!0,g("cutOn").checked=!0,w.organelle&&Oe(!1),y.flyTo(u,2200,()=>{ee(a,$(a,r)-20);const m=I.loaded.get(i);if(m)d(m);else{const l=h=>{h.idx===i&&(I.onTileLoaded.splice(I.onTileLoaded.indexOf(l),1),d(h))};I.onTileLoaded.push(l)}}),ee(a,$(a,r)-20)}async function We(s,e){const t=z.width/S.clientWidth;return z.pick(s*t,e*t)}function St(s){return s.src==="PDB"?"PDB "+s.pdb:s.src==="AF-multimer"?"AlphaFold complex":"AlphaFold"}function lt(s){const e=s>>>30;if(s===0)return{kind:"none",title:"Outside the cell"};if(e===1){const t=I.instanceAt(s&1073741823);return t?{kind:"inst",id:s,info:t,title:J.types[t.type].gene,sub:Ye[t.region]}:{kind:"none",title:""}}if(e===2){const t=s>>>1&15,i=s&1,n=dt[t]||"membrane";return{kind:"lipid",id:s,mem:t,leaf:i,title:"Lipid bilayer",sub:`${n} · ${t===8?"phospholipid monolayer":i?"cytoplasmic leaflet":t===0?"exoplasmic (outer) leaflet":"lumenal leaflet"}`}}if(e===3){const t=s>>>16&3,i=s&65535,n=["mRNA","Linker DNA","β-1,3-glucan fibril","Mannan chain"][t],r=["nucleotide-resolution thread (schematic path)","between consecutive nucleosomes, 2 bp per bead","procedural, 36 glucose per bead","procedural α-mannan on a cell-wall mannoprotein"][t];return{kind:"free",id:s,fk:t,chain:i,title:n,sub:r}}return s>=512&&s<528?{kind:"band",id:s,title:"Membrane (cut face)",sub:dt[s-512]}:s>=256&&s<272?{kind:"label",id:s,title:si[s-256],sub:"compartment (not loaded at this distance, or beyond the slab)"}:{kind:"none",title:""}}function je(s){Tt=s,w.selId=s&&s.kind==="inst"?s.id>>>0:0,It(),T=!0}function It(){const s=g("inspect"),e=Tt;if(!e||e.kind==="none"){s.style.display="none";return}s.style.display="block";let t='<button class="close" id="inspClose">×</button>';if(e.kind==="inst"){const i=J.types[e.info.type],n=i.src==="PDB",r=i.nb,o=i.nnt;t+=`<h2>${i.gene}</h2><div>${i.name}</div><div style="margin-top:6px"><span class="badge ${n?"pdb":"af"}">${St(i)}</span>
      ${i.oligomer>1?`<span class="badge" style="background:#666">${i.oligomer}-mer</span>`:""}</div>
      <div class="kv">
      ${i.genes&&i.genes!==i.gene?`<div>Subunits</div><div style="word-break:break-word">${i.genes.split(/;\s*/).slice(0,40).join(", ")}${i.genes.split(";").length>40?" …":""}</div>`:""}
      ${i.uniprot?`<div>UniProt</div><div>${i.uniprot}</div>`:""}
      <div>Location</div><div>${Ye[e.info.region]}</div>
      <div>Copies / cell</div><div>${i.copies.toLocaleString("en-US")}</div>
      <div>Mass</div><div>${(i.mass/1e3).toFixed(1)} kDa</div>
      <div>Residues drawn</div><div>${(r-o).toLocaleString("en-US")} aa${o?` + ${o.toLocaleString("en-US")} nt`:""}</div>
      <div>Radius</div><div>${i.br.toFixed(1)} nm (Rg ${i.rg.toFixed(1)} nm)</div>
      ${!n&&i.plddt_mean!=null?`<div>Mean pLDDT</div><div>${i.plddt_mean.toFixed(1)} (${Math.round((i.plddt_frac_kept??1)*100)}% of residues kept, pLDDT ≥ 50)</div>`:""}
      <div>Position</div><div>${e.info.p.map(c=>c.toFixed(1)).join(", ")} nm</div>
      </div>
      <div class="grid2"><button class="small" id="loadAtoms">Load atomic model</button><button class="small" id="hlType">Highlight all</button>
      <button class="small" id="centerOn">Center view</button><button class="small" id="hideType">${le[e.info.type]&2?"Show":"Hide"} this type</button></div>
      <div style="margin-top:9px">${i.uniprot?`<a href="https://www.uniprot.org/uniprotkb/${i.uniprot}" target="_blank">UniProt</a>`:""}
      ${n?`<a href="https://www.rcsb.org/structure/${i.pdb}" target="_blank">RCSB ${i.pdb}</a>`:""}
      ${(i.model||i.uniprot)&&!n?`<a href="https://alphafold.ebi.ac.uk/entry/${i.model||i.uniprot}" target="_blank">AlphaFold DB</a>`:""}
      ${/^[A-Z0-9-]+$/.test(i.gene)&&!Pt[i.id]&&!i.id.startsWith("CPX")?`<a href="https://www.yeastgenome.org/locus/${i.gene}" target="_blank">SGD</a>`:""}</div>`}else t+=`<h2>${e.title}</h2><div class="sub" style="font-size:12.5px">${e.sub||""}</div>`,e.kind==="lipid"&&(t+='<p class="sub">Lipids are generated in your browser on the mesh of this membrane. Classes are PI, PE, PC, PA and PS, with ergosterol drawn as two beads. They sit at 0.65 nm² per lipid and are removed where they would overlap protein.</p>');s.innerHTML=t,g("inspClose").onclick=()=>je(null),e.kind==="inst"&&(g("loadAtoms").onclick=()=>wi(e.info.type),g("hlType").onclick=()=>ke(e.info.type),g("centerOn").onclick=()=>{const i=y.state();i.target=e.info.p,y.flyTo(i,700)},g("hideType").onclick=()=>{le[e.info.type]^=2,z.device.queue.writeBuffer(z.typeFlags,0,le),It(),T=!0})}let De=null,Ce=null;function vi(){if(Ce)return Ce;const s="https://cdn.jsdelivr.net/npm/molstar@4.18.0/build/viewer/",e=t=>new Promise((i,n)=>{const r=document.createElement("link");r.rel="stylesheet",r.href=t+"molstar.css",document.head.appendChild(r);const o=document.createElement("script");o.src=t+"molstar.js",o.onload=()=>i(),o.onerror=()=>{r.remove(),o.remove(),n(new Error("Mol* failed to load"))},document.head.appendChild(o)});return Ce=(window.__YC_EMBED__?e(s):e("molstar/").catch(()=>e(s))).catch(t=>{throw Ce=null,t}),Ce}async function wi(s){const e=J.types[s];g("molModal").style.display="flex",g("molTitle").textContent=`${e.gene} · ${e.name} · loading…`;try{await vi();const t=window.molstar;De?await De.plugin.clear():De=await t.Viewer.create("molHost",{layoutIsExpanded:!1,layoutShowControls:!1,layoutShowRemoteState:!1,layoutShowSequence:!0,layoutShowLog:!1,layoutShowLeftPanel:!1,viewportShowExpand:!1,viewportShowSelectionMode:!1,viewportShowAnimation:!1,collapseLeftPanel:!0,pdbProvider:"rcsb",emdbProvider:"rcsb"});let i="",n="";if(e.src==="PDB"&&e.pdb){i=`https://files.rcsb.org/download/${e.pdb}-assembly1.cif`,n=`https://www.rcsb.org/structure/${e.pdb}`;try{(await fetch(i,{method:"HEAD"})).ok||(i=`https://files.rcsb.org/download/${e.pdb}.cif`)}catch{i=`https://files.rcsb.org/download/${e.pdb}.cif`}}else{const o=e.model||e.uniprot;n=`https://alphafold.ebi.ac.uk/entry/${o}`;try{const d=await(await fetch(`https://alphafold.ebi.ac.uk/api/prediction/${o}`)).json();i=(Array.isArray(d)?d[0]:d)?.cifUrl}catch{}i||(i=e.model?`https://alphafold.ebi.ac.uk/files/${e.model}-model_v1.cif`:`https://alphafold.ebi.ac.uk/files/AF-${e.uniprot}-F1-model_v6.cif`)}await De.loadStructureFromUrl(i,"mmcif",!1),g("molTitle").textContent=`${e.gene} · ${e.name} · ${e.src==="PDB"?"PDB "+e.pdb+" (biological assembly)":e.model?`AlphaFold DB predicted ${e.oligomer}-mer (full length, incl. low-pLDDT regions)`:"AlphaFold DB model (full length, incl. low-pLDDT regions)"}`;const r=g("molLink");r.href=n,r.textContent="Open source page"}catch(t){g("molTitle").textContent=`Could not load the atomic model: ${t?.message||t}. Internet access to RCSB/AlphaFold DB is needed.`}}function Oe(s){w.organelle=s,g("modeSection").classList.toggle("on",!s),g("modeOrg").classList.toggle("on",s),s&&!z.organelleMeshes.length&&(Me("Loading organelle surfaces (6.7 MB)…"),Ie("organelle_meshes.bin.gz").then(e=>{z.uploadOrganelleMeshes(e),T=!0})),T=!0}function Ot(s){Et=s;const e={low:[2.6,.75],med:[1.6,1],high:[1,1],ultra:[.6,1],auto:[1.6,1]},[t,i]=e[s];w.l0px=t,w.renderScale=i,z.resize(i),T=!0}function mt(){const s=g("legend"),e=(i,n)=>`<div><span class="sw" style="background:${i}"></span>${n}</div>`,t=w.colorMode;t===0?s.innerHTML=[e(N[0],"cytosolic protein"),e(N[1],"nuclear protein"),e("#EE9BC0","ribosome protein"),e("#C9558F","rRNA"),e("#E27FA8","tRNA"),e("#B0437A","mRNA"),e("#F4B26A","histone"),e("#E07636","DNA"),e(N[23],"nuclear pore"),e("#E0893A","actin filament"),e("#74C476","microtubule"),e("#2E8B57","γ-tubulin ring"),e(N[14],"membrane protein"),e(N[5],"mitochondrial"),e(N[3],"ER lumen"),e("#F3C643","lipid head"),e("#F8E3A0","lipid tail"),e("#CDA866","β-glucan"),e("#C3C27A","mannan"),e(N[9],"wall protein"),e("#8B88C9","large complex")].join(""):t===1?s.innerHTML=e("#4C8CD9","experimental (PDB)")+e("#F5A84C","AlphaFold model"):t===2?s.innerHTML=e("#0053D6","pLDDT > 90")+e("#65CBF3","70–90")+e("#FFDB13","50–70")+e("#B8B8B8","experimental (PDB)"):t===3?s.innerHTML='<div style="grid-column: 1 / 3"><div style="height:8px;border-radius:3px;background:linear-gradient(90deg,#440154,#3b528b,#21918c,#5ec962,#fde725)"></div><div style="display:flex;justify-content:space-between"><span>1</span><span>10³</span><span>10⁶ copies</span></div></div>':t===5?s.innerHTML=[0,1,2,3,4,5,6,7,8,9,14,15,18,19].map(i=>e(N[i],Ye[i])).join(""):s.innerHTML=""}function he(){const s=w.plane[3];g("planeD").value=String(Math.round(s)),g("planeDv").textContent=Math.abs(s-Math.round(s))<.05?`${Math.round(s)} nm`:`${s.toFixed(1)} nm`;const e=w.slab;g("slab").value=String(Math.round(Math.log(e/5)/Math.log(4400/5)*1e3)),g("slabv").textContent=kt(e),["axX","axY","axZ"].forEach((t,i)=>g(t).classList.toggle("on",Math.abs(w.plane[i])>.999)),nt=!0}function xi(){const s=J.types.length;let e=0;for(const l of J.types)e+=l.copies;if(g("subCount").textContent=`${(e/1e6).toFixed(1)} M macromolecules of ${s.toLocaleString("en-US")} species, each a real PDB or AlphaFold structure`,window.__YC_PREVIEW__){const l=Se.nInst.reduce((h,f)=>h+f,0);g("subCount").innerHTML=`<b>Preview:</b> ${Se.n} of 6,416 tiles (${(l/1e3).toFixed(0)} k molecules) from the full ${(e/1e6).toFixed(1)} M-molecule cell. The full cell is in yeast_cell_viewer.zip.`,g("modeOrg").disabled=!0,g("modeOrg").title="Organelle surfaces are only in the full app"}const t=g("views");fe.forEach((l,h)=>{const f=document.createElement("button");f.textContent=`${(h+1)%10===0&&h<10?0:h<9?h+1:""}${h<10?" · ":""}${l.name}`,f.title=l.caption,f.onclick=()=>{be(),Ae(h)},t.appendChild(f)}),g("tourBtn").onclick=()=>U?be():Ct();const i=g("search"),n=g("results");let r=0,o=[];const c=()=>{n.innerHTML=o.map((l,h)=>{const f=J.types[l];return`<div data-i="${l}" class="${h===r?"act":""}"><span class="nm"><b>${f.gene}</b> ${f.name}</span><span class="n">${re(f.copies)}</span></div>`}).join(""),n.querySelectorAll("div").forEach(l=>l.addEventListener("click",()=>{ke(Number(l.dataset.i)),o=[],c(),i.value=""}))};i.addEventListener("input",()=>{o=yi(i.value),r=0,c()}),i.addEventListener("keydown",l=>{l.key==="ArrowDown"?(r=Math.min(o.length-1,r+1),c(),l.preventDefault()):l.key==="ArrowUp"?(r=Math.max(0,r-1),c(),l.preventDefault()):l.key==="Enter"&&o.length?(ke(o[r]),Bt(o[r]),o=[],c(),i.value="",i.blur()):l.key==="Escape"&&(o=[],c(),i.blur()),l.stopPropagation()}),g("cutOn").onchange=l=>{w.cutOn=l.target.checked,T=!0};const d=l=>()=>{const h=[0,0,0];h[l]=$(y.fwd,[l===0?1:0,l===1?1:0,l===2?1:0])>=0?1:-1,w.cutOn=!0,g("cutOn").checked=!0,ee(h,$(h,y.target))};g("axX").onclick=d(0),g("axY").onclick=d(1),g("axZ").onclick=d(2),g("alignView").onclick=()=>{w.cutOn=!0,g("cutOn").checked=!0,Pe()},g("followView").onclick=()=>{pe=!pe,g("followView").classList.toggle("on",pe),pe&&Pe()},g("planeD").oninput=l=>{const h=Number(l.target.value),f=ie(),M=h-w.plane[3];y.target=_(y.target,L(f,M)),ee(f,h)},g("slab").oninput=l=>{const h=Number(l.target.value)/1e3;w.slab=5*Math.pow(4400/5,h),he(),T=!0},g("sweepBtn").onclick=()=>{if(oe){oe=null,g("sweepBtn").classList.remove("on");return}const l=ie();w.cutOn=!0,g("cutOn").checked=!0;const h=xe.G[0];oe={t0:performance.now(),d0:-h,d1:h,dur:3e4},y.target=_(y.target,L(l,-h-w.plane[3])),g("sweepBtn").classList.add("on")},g("flipBtn").onclick=()=>{const l=ie(),h=y.state(),f=y.target,M=L(l,-1);ee(M,-w.plane[3]),h.fwd=L(h.fwd,-1),h.target=f,y.flyTo(h,1200)},g("modeSection").onclick=()=>Oe(!1),g("modeOrg").onclick=()=>Oe(!0),g("colorMode").onchange=l=>{w.colorMode=Number(l.target.value),mt(),T=!0},g("shaded").onchange=l=>{w.shaded=l.target.checked,T=!0},g("outlines").onchange=l=>{w.outline=l.target.checked?1:0,T=!0},mt();const a=g("showList"),u=[["Ribosomes","#EE9BC0",l=>Ge(["ribosome_80S"],!l)],["tRNA","#E27FA8",l=>Ge(["tRNA"],!l)],["Nucleosomes","#F4B26A",l=>Ge(["nucleosome"],!l)],["Nuclear pores",N[23],l=>Ge(["NPC"],!l)],["mRNA","#B0437A",l=>Ne(0,l)],["Linker DNA","#E07636",l=>Ne(1,l)],["β-Glucan","#CDA866",l=>Ne(2,l)],["Mannan","#C3C27A",l=>Ne(3,l)],["Lipids","#F3C643",l=>{w.showLipids=l,T=!0}],["Membranes","#F8E3A0",l=>{w.showMeshes=l,T=!0}],["Soluble proteins",N[0],l=>et([0,1,2,3,4,5,6,7,8,10,11,12,13],l)],["Membrane proteins",N[14],l=>et([14,15,16,17,18,19,20,21],l)]];for(const[l,h,f]of u){const M=document.createElement("label");M.className="chk",M.innerHTML=`<input type="checkbox" checked><span class="sw" style="background:${h}"></span>${l}`,M.querySelector("input").onchange=v=>f(v.target.checked),a.appendChild(M)}const m=g("regionList");Ye.forEach((l,h)=>{const f=document.createElement("label");f.className="chk",f.innerHTML=`<input type="checkbox" checked data-r="${h}"><span class="sw" style="background:${N[h]}"></span>${l}`,f.querySelector("input").onchange=M=>et([h],M.target.checked),m.appendChild(f)}),g("quality").onchange=l=>Ot(l.target.value),g("detailOn").onchange=l=>{if(I.detailEnabled=l.target.checked,!I.detailEnabled)for(const h of I.loaded.values())h.det&&I.dropDetail(h);T=!0},g("shotBtn").onclick=_t,g("helpBtn").onclick=()=>g("helpModal").style.display="flex",g("aboutBtn").onclick=()=>g("aboutModal").style.display="flex",g("aboutBody").innerHTML=pi,document.querySelectorAll(".modal").forEach(l=>l.addEventListener("click",h=>{(h.target===l||h.target.classList.contains("close"))&&(l.style.display="none")})),g("molClose").onclick=()=>g("molModal").style.display="none",g("panelToggle").onclick=()=>{g("panel").classList.toggle("collapsed"),g("panelToggle").textContent=g("panel").classList.contains("collapsed")?"›":"‹"},he(),ct(),Pi()}function Ge(s,e){for(const t of s){const i=J.byId.get(t);i!==void 0&&(e?le[i]|=2:le[i]&=-3)}z.device.queue.writeBuffer(z.typeFlags,0,le),T=!0}function Ne(s,e){e?w.freeHide&=~(1<<s):w.freeHide|=1<<s,T=!0}function et(s,e){for(const t of s)e?w.regionMask|=1<<t:w.regionMask&=~(1<<t);w.regionMask>>>=0,document.querySelectorAll("#regionList input").forEach(t=>{t.checked=(w.regionMask>>>Number(t.dataset.r)&1)===1}),T=!0}function ct(){let s=null;const e=new Set;S.addEventListener("contextmenu",i=>i.preventDefault()),S.addEventListener("pointerdown",i=>{S.setPointerCapture(i.pointerId),S.focus(),s={x:i.clientX,y:i.clientY,btn:i.button,moved:0,shift:i.shiftKey},S.classList.add("dragging"),U&&be()}),S.addEventListener("pointermove",i=>{if(!s){Mi(i.clientX,i.clientY);return}const n=i.clientX-s.x,r=i.clientY-s.y;if(s.x=i.clientX,s.y=i.clientY,s.moved+=Math.abs(n)+Math.abs(r),s.btn===2||s.btn===1||s.shift){const o=y.pxPerNmAtUnit(S.clientHeight)/y.dist;y.pan(n,r,o)}else y.orbit(n*.005,r*.005),pe&&Pe();y.anim=null,T=!0});const t=i=>{if(!s)return;const n=s.moved<4&&s.btn===0;s=null,S.classList.remove("dragging"),n&&We(i.offsetX,i.offsetY).then(r=>{const o=lt(r.id);je(o.kind==="inst"||o.kind==="lipid"||o.kind==="free"||o.kind==="band"?o:null)})};S.addEventListener("pointerup",t),S.addEventListener("pointercancel",t),S.addEventListener("dblclick",i=>{We(i.offsetX,i.offsetY).then(n=>{if(n.depth<=0)return;const r=y.near(),o=r/n.depth,c=S.clientWidth/S.clientHeight,d=y.ray(i.offsetX/S.clientWidth,i.offsetY/S.clientHeight,c),a=o/$(d,y.fwd),u=_(y.eye,L(d,a)),m=y.state();m.target=u,m.dist=Math.max(y.minDist*4,y.dist*.55),y.flyTo(m,900)})}),S.addEventListener("wheel",i=>{i.preventDefault();const n=i.deltaY||i.deltaX,r=i.deltaMode===1?n*30:n;if(i.altKey){const o=ie(),c=Math.sign(r)*Math.min(2,Math.abs(r)*.005);y.target=_(y.target,L(o,c)),ee(o,w.plane[3]+c)}else if(i.shiftKey){const o=ie(),c=Math.max(1,Math.min(w.slab,y.dist*.05))*Math.sign(r)*Math.min(4,Math.abs(r)/60+.5);y.target=_(y.target,L(o,c)),ee(o,w.plane[3]+c)}else{const o=Math.exp(r*.0015*(i.metaKey?.125:1)),c=S.clientWidth/S.clientHeight,d=y.ray(i.offsetX/S.clientWidth,i.offsetY/S.clientHeight,c),a=y.dist,u=_(y.eye,L(d,a/$(d,y.fwd))),m=y.dist;y.zoom(o);const l=1-y.dist/m;y.target=_(y.target,L(ae(u,y.target),l))}y.anim=null,U&&be(),T=!0},{passive:!1}),window.addEventListener("keydown",i=>{if(i.target.tagName==="INPUT"||i.target.tagName==="SELECT")return;const n=i.key.toLowerCase();if(e.add(n),n==="c")w.cutOn=!w.cutOn,g("cutOn").checked=w.cutOn,T=!0;else if(n==="v")w.cutOn=!0,g("cutOn").checked=!0,Pe();else if(n==="[")w.slab=Math.max(5,w.slab/1.25),he(),T=!0;else if(n==="]")w.slab=Math.min(4400,w.slab*1.25),he(),T=!0;else if(n===","||n==="."||i.code==="Comma"||i.code==="Period"){const r=n==="."||i.code==="Period",o=ie(),c=(r?1:-1)*(i.altKey?.5:Math.max(2,w.slab*.25));y.target=_(y.target,L(o,c)),ee(o,w.plane[3]+c)}else if(n==="o")Oe(!w.organelle);else if(n==="t")U?be():Ct();else if(n==="p")_t();else if(n==="escape")je(null),ke(-1);else if(/^[0-9]$/.test(n)){const r=n==="0"?9:Number(n)-1;r<fe.length&&(be(),Ae(r))}}),window.addEventListener("keyup",i=>e.delete(i.key.toLowerCase())),window.addEventListener("blur",()=>e.clear()),ct.keys=e,window.addEventListener("resize",()=>{z.resize(w.renderScale),T=!0})}let gt=0,tt=!1;function Mi(s,e){const t=performance.now();tt||t-gt<90||(gt=t,tt=!0,T=!0,We(s,e).then(i=>{tt=!1;const n=lt(i.id),r=g("hover");if(!n.title||n.kind==="none"){r.style.display="none";return}let o=`<b>${n.title}</b>`;if(n.kind==="inst"){const c=J.types[n.info.type];o+=` · ${c.name} <span style="color:#aaa">· ${n.sub}</span>`}else o+=` <span style="color:#aaa">· ${n.sub}</span>`;r.innerHTML=o,r.style.display="block",r.style.left=`${s+14}px`,r.style.top=`${e+12}px`}))}S.addEventListener("pointerleave",()=>g("hover").style.display="none");let nt=!0,bt=0,ue=null;const Fe={data:null};function Pi(){const s=g("minimap");s.addEventListener("click",e=>{if(!ue)return;const t=s.getBoundingClientRect(),i=((e.clientX-t.left)/t.width*2-1)*ue.ext,n=(1-(e.clientY-t.top)/t.height*2)*ue.ext,r=_(ue.c,_(L(ue.u,i),L(ue.v,n))),o=y.state();o.target=r,y.flyTo(o,1200),U&&be()})}function rt(s){const e=g("minimap"),t=e.getContext("2d"),i=e.width,n=e.height,r=w.cutOn?ie():y.fwd,o=w.cutOn?w.plane[3]:$(y.fwd,y.target);if(nt&&s-bt>250){nt=!1,bt=s;let v=ae(y.right,L(r,$(y.right,r)));ge(v)<.001&&(v=X(de(r,[0,0,1]))),v=X(v);let p=de(r,v);$(p,y.up)<0&&(p=L(p,-1));const b=2300,P=L(r,o);ue={u:v,v:p,c:P,ext:b};const x=120,A=t.createImageData(x,x),k=[N[80],N[9],"#D2BE8E",N[83],N[3],N[3],N[1],N[2],N[4],N[6],N[5],"#E9C46A",N[7],N[8],"#C9C2DC"].map(F=>{const C=parseInt(F.slice(1),16);return[C>>16&255,C>>8&255,C&255]});k[3]=[196,214,235];for(let F=0;F<x;F++)for(let C=0;C<x;C++){const B=((C+.5)/x*2-1)*b,R=(1-(F+.5)/x*2)*b,q=_(P,_(L(v,B),L(p,R))),Q=Lt.label(q[0],q[1],q[2]),K=k[Q],V=4*(F*x+C);A.data[V]=K[0],A.data[V+1]=K[1],A.data[V+2]=K[2],A.data[V+3]=Q===0?0:255}Fe.data=A}if(!Fe.data||!ue)return;const c=rt.tmp||(rt.tmp=document.createElement("canvas"));c.width=Fe.data.width,c.height=Fe.data.height,c.getContext("2d").putImageData(Fe.data,0,0),t.clearRect(0,0,i,n),t.imageSmoothingEnabled=!0,t.drawImage(c,0,0,i,n);const d=ue,a=v=>{const p=ae(v,d.c);return[($(p,d.u)/d.ext+1)/2*i,(1-$(p,d.v)/d.ext)/2*n]},u=S.clientWidth/S.clientHeight,m=y.dist*Math.tan(y.fovy/2),l=m*u,h=[[-1,-1],[1,-1],[1,1],[-1,1]].map(([v,p])=>a(_(y.target,_(L(y.right,v*l),L(y.up,p*m)))));t.strokeStyle="#FF9400",t.lineWidth=3,t.beginPath(),h.forEach(([v,p],b)=>b?t.lineTo(v,p):t.moveTo(v,p)),t.closePath(),t.stroke();const[f,M]=a(y.target);t.fillStyle="#FF9400",t.beginPath(),t.arc(f,M,5,0,7),t.fill(),g("mmLabel").textContent=w.cutOn?`plane ${Math.round(o)} nm`:"view plane"}function _t(){const s=w.renderScale;z.resize(2*s),w.renderScale=2*s,Rt(!0),S.toBlob(e=>{if(e){const t=document.createElement("a");t.href=URL.createObjectURL(e),t.download=`yeast_cell_${new Date().toISOString().replace(/[:.]/g,"-").slice(0,19)}.png`,t.click(),setTimeout(()=>URL.revokeObjectURL(t.href),5e3)}w.renderScale=s,z.resize(s),T=!0},"image/png"),Me("Screenshot saved (2× resolution)")}function Ai(){const s=I.totals(),e=z.stats,t=Z.length>3?1e3/(Z.reduce((a,u)=>a+u)/Z.length):0,i=I.pendingForView();g("hud").innerHTML=`<b>${re(s.inst)}</b> molecules in ${s.tiles} tiles${i?` · <span style="color:#FF9400">streaming ${i}</span>`:""}${I.starved?' · <span style="color:#FF9400">memory budget reached (thinner slab loads more)</span>':""}<br>
    drawn: <b>${re(e.ninst)}</b> in view · ${re(e.l0*64)} residue beads (≤) · ${re(e.l1*16)} LOD1 · ${re(e.l2)} LOD2 · ${re(e.lip)} lipids · ${re(e.fre)} polymer beads<br>
    ${t?`${t.toFixed(0)} fps · `:""}LOD bead ≥ ${w.l0px.toFixed(1)} px · GPU pools ${z.gpuMemoryMB().toFixed(0)} MB · ${(I.bytesLoaded/1048576).toFixed(0)} MB streamed${z.adapterInfo?" · "+z.adapterInfo:""}`;const n=z.width/S.clientWidth,r=y.pxPerNmAtUnit(z.height)/y.dist/n,o=120/r,c=Math.pow(10,Math.floor(Math.log10(o))),d=[1,2,5,10].map(a=>a*c).reduce((a,u)=>Math.abs(u-o)<Math.abs(a-o)?u:a);g("scalebar").querySelector(".bar").style.width=`${d*r}px`,g("scalebar").querySelector("span").textContent=kt(d)}let zt=0,He=!1;function Rt(s=!1){const e=performance.now();z.frame(y,w,xe.G,I.meshDraws(),s||e-z.lastStatsT>400),z.lastStatsT=performance.now(),He=!0,z.device.queue.onSubmittedWorkDone().then(()=>{zt=performance.now()-e,He=!1})}let yt=0;function Ut(){I.update(y,S.clientWidth/S.clientHeight,w.cutOn,w.plane,w.slab,w.organelle,z.height),at=performance.now(),(I.inflight||I.detInflight)&&(T=!0)}setInterval(()=>{I&&performance.now()-at>150&&Ut()},150);function Dt(s){requestAnimationFrame(Dt);const e=Ze?s-Ze:16;Ze=s;const t=ct.keys;if(t&&t.size){const i=y.dist*.9*e/1e3;let n=[0,0,0];t.has("w")&&(n=_(n,L(y.fwd,i))),t.has("s")&&(n=_(n,L(y.fwd,-i))),t.has("a")&&(n=_(n,L(y.right,-i))),t.has("d")&&(n=_(n,L(y.right,i))),t.has("q")&&(n=_(n,L(y.up,-i))),t.has("e")&&(n=_(n,L(y.up,i))),t.has("z")&&y.roll(-e*.0015),t.has("x")&&y.roll(e*.0015),(ge(n)>0||t.has("z")||t.has("x"))&&(y.fly(n),T=!0,pe&&Pe())}if(y.step(s)&&(T=!0,pe&&!U&&Pe()),oe){const i=Math.min(1,(s-oe.t0)/oe.dur),n=ie(),r=oe.d0+(oe.d1-oe.d0)*i;y.target=_(y.target,L(n,r-w.plane[3])),w.plane[3]=r,he(),i>=1&&(oe=null,g("sweepBtn").classList.remove("on")),T=!0}if(bi(s),s-at>90&&Ut(),T&&!He){if(T=!1,Rt(),Z.push(e),Z.length>30&&Z.shift(),Et==="auto"&&Z.length>=20){const i=Z.reduce((n,r)=>n+r)/Z.length;i>45&&w.l0px<3.5?(w.l0px=Math.min(3.5,w.l0px*1.15),Z.length=0):i<22&&w.l0px>.9&&(w.l0px=Math.max(.9,w.l0px/1.1),Z.length=0)}setTimeout(()=>{T=T||z.statsBusy},0)}else T||(Z.length=0);s-yt>300&&(yt=s,Ai(),rt(s),me>=0&&I.inflight&&Ft())}mi();
