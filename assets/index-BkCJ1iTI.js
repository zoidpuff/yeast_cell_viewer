(function(){const e=document.createElement("link").relList;if(e&&e.supports&&e.supports("modulepreload"))return;for(const s of document.querySelectorAll('link[rel="modulepreload"]'))n(s);new MutationObserver(s=>{for(const o of s)if(o.type==="childList")for(const r of o.addedNodes)r.tagName==="LINK"&&r.rel==="modulepreload"&&n(r)}).observe(document,{childList:!0,subtree:!0});function t(s){const o={};return s.integrity&&(o.integrity=s.integrity),s.referrerPolicy&&(o.referrerPolicy=s.referrerPolicy),s.crossOrigin==="use-credentials"?o.credentials="include":s.crossOrigin==="anonymous"?o.credentials="omit":o.credentials="same-origin",o}function n(s){if(s.ep)return;s.ep=!0;const o=t(s);fetch(s.href,o)}})();const X=1024,Q=4096,ve=`
const PAGE: u32 = ${X}u;
const SPAGE: u32 = ${Q}u;
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
  caps: vec4u,            // x capL0 (offset of L1 items in list01), y motion amount x 1000, z motion cap nm x 1000
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
`,Tt=`
struct TypeInfo {
  beadOff: u32, beadCount: u32, clusterOff: u32, clusterCount: u32,
  stOff: u32, stCount: u32, br: f32, req: f32,
  color0: u32, color1: u32, flags: u32, plddt: f32,
  modeOff: u32, coefOff: u32, nConf: u32, pad: u32,   // pool-ensemble motion modes (nConf 0: none)
};
struct Inst { t: u32, region: u32, flags: u32, p: vec3f, q: vec4f, seed: u32 };
@group(0) @binding(1) var<storage, read> inst: array<vec4u>;
@group(0) @binding(2) var<storage, read> pageOrg: array<vec4f>;       // xyz origin (camera-relative is computed), w = bitcast count
@group(0) @binding(3) var<storage, read> types: array<TypeInfo>;
@group(0) @binding(4) var<storage, read> typeFlags: array<u32>;       // bit0 highlight, bit1 hidden, bit2 own class, bit3 no motion

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
  o.seed = hash32(w.y ^ hash32(w.z ^ hash32(w.w)));   // stable per instance (tile-relative pose bits)
  return o;
}

// ---- illustrative motion (u.caps.y = amount x 1000; time u.cut.w in s). Not physics: smooth pseudo-random
// rigid-body wobble per instance plus flexing of low-confidence regions, slowed to about 1 Hz.
fn rnd(s: u32, k: u32) -> f32 { return f32(hash32(s ^ (k * 0x9E3779B9u)) & 0xFFFFu) / 65535.0; }
fn wob(s: u32, k: u32, t: f32) -> f32 {          // smooth signal in [-1, 1]
  return 0.6 * sin(t * (0.5 + rnd(s, k)) + 6.2832 * rnd(s, k + 101u))
       + 0.4 * sin(t * (1.1 + 1.3 * rnd(s, k + 211u)) + 6.2832 * rnd(s, k + 307u));
}
fn qmulq(a: vec4f, b: vec4f) -> vec4f {
  return vec4f(a.w * b.xyz + b.w * a.xyz + cross(a.xyz, b.xyz), a.w * b.w - dot(a.xyz, b.xyz));
}
fn animInst(I: Inst, br: f32) -> Inst {
  if (u.caps.y == 0u || (typeFlags[I.t] & 8u) != 0u) { return I; }
  var o = I;
  let a = f32(u.caps.y) / 1000.0; let t = u.cut.w; let s = I.seed;
  let sz = max(br, 2.0);
  var dp = a * 0.35 * sqrt(2.0 / sz) * vec3f(wob(s, 1u, t), wob(s, 2u, t), wob(s, 3u, t));   // nm; larger = slower
  var ax = normalize(vec3f(wob(s, 4u, 0.3 * t), wob(s, 5u, 0.3 * t), wob(s, 6u, 0.3 * t)) + vec3f(1e-3));
  if (I.region >= 14u && I.region <= 21u) {       // membrane protein: stay in the bilayer
    let nrm = qrot(I.q, vec3f(0.0, 0.0, 1.0));
    dp -= nrm * dot(dp, nrm); ax = nrm;
  }
  let h = 0.5 * a * 0.35 / sz * wob(s, 7u, 0.8 * t);  // half-angle: the rim moves about 0.35 nm
  o.p = I.p + dp;
  o.q = qmulq(vec4f(ax * sin(h), cos(h)), I.q);
  return o;
}
// per-bead flexing in the molecule frame, spatially coherent (~5 nm wavelength); amplitude from pLDDT
fn flexBead(b: vec3f, plddt: f32, s: u32) -> vec3f {
  if (u.caps.y == 0u) { return b; }
  var amp = 0.04;
  if (plddt < 254.5) { amp = mix(0.45, 0.04, clamp((plddt - 50.0) / 40.0, 0.0, 1.0)); }
  let t = u.cut.w * 0.9; let k = 1.25;
  let d = vec3f(sin(dot(b, vec3f(k, 0.4, -0.3)) + t * (0.8 + rnd(s, 11u)) + 6.2832 * rnd(s, 12u)),
                sin(dot(b, vec3f(-0.3, k, 0.5)) + t * (0.7 + rnd(s, 13u)) + 6.2832 * rnd(s, 14u)),
                sin(dot(b, vec3f(0.5, -0.2, k)) + t * (0.9 + rnd(s, 15u)) + 6.2832 * rnd(s, 16u)));
  return b + (f32(u.caps.y) / 1000.0) * amp * d;
}
`,Yt=ve+Tt+`
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
  // bit2: molecule classes with their own Show switch (ribosomes, cytoskeleton, ...) ignore the compartment filter
  if ((tf & 4u) == 0u && ((u.modes.y >> region) & 1u) == 0u) { return false; }
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
`,Xt=ve+`
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
  if (kind == 5u) { return 0.45; }  // disordered protein residue (FG repeats)
  return 0.6;                        // mannan, pre-rRNA
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
`,Qt=`
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
`,Ct=`
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
`,Kt=`
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
`,ut=ve+Tt+Ct+Kt,Jt=ut+`
@group(1) @binding(0) var<storage, read> beads: array<vec4f>;     // xyz, w = bitcast(kind | plddt<<8)
@group(1) @binding(1) var<storage, read> list01: array<vec2u>;
@group(1) @binding(2) var<storage, read> modes: array<u32>;      // per bead 4 x (int8 xyz), then f32 conformer paths
const N_MODES: u32 = 4u;
fn i8x(w: u32, sh: u32) -> f32 { return f32(bitcast<i32>(w << (24u - sh)) >> 24u); }
fn coefAt(T: TypeInfo, k: u32, j: u32) -> f32 { return bitcast<f32>(modes[T.coefOff + (k % T.nConf) * N_MODES + j]); }
// displacement of bead 'local' along the protein's pool-ensemble path: Catmull-Rom through the sampled conformers'
// mode coefficients, about 2 s per conformer, each instance at its own phase
fn modeDisp(T: TypeInfo, local: u32, s: u32) -> vec3f {
  let n = T.nConf;
  let ph = u.cut.w * 0.5 + rnd(s, 21u) * f32(n);
  let k = u32(floor(ph)) % n; let f = fract(ph);
  let f2 = f * f; let f3 = f2 * f;
  var d = vec3f(0.0);
  for (var j = 0u; j < N_MODES; j++) {
    let p0 = coefAt(T, k + n - 1u, j); let p1 = coefAt(T, k, j); let p2 = coefAt(T, k + 1u, j); let p3 = coefAt(T, k + 2u, j);
    let c = 0.5 * (2.0 * p1 + (p2 - p0) * f + (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * f2 + (3.0 * p1 - p0 - 3.0 * p2 + p3) * f3);
    let w = modes[T.modeOff + local * N_MODES + j];
    d += c * vec3f(i8x(w, 0u), i8x(w, 8u), i8x(w, 16u));
  }
  let cap = f32(u.caps.z) / 1000.0;     // nm, largest excursion of one bead (modes.json "cap"; 0 = none)
  let L = length(d);
  if (cap > 0.0 && L > cap) { d *= cap / L; }
  return d * f32(u.caps.y) / 1000.0;
}
@vertex
fn vsL0(@builtin(vertex_index) vi: u32, @builtin(instance_index) ii: u32) -> VOut {
  let item = list01[ii];
  let j = vi / 6u; let corner = vi % 6u;
  let cnt = (item.y >> 22u) + 1u;
  if (j >= cnt) { var o: VOut; o.pos = vec4f(2.0, 2.0, 2.0, 1.0); return o; }
  let I0 = decodeInst(item.x); let T = types[I0.t];
  let I = animInst(I0, T.br);
  let b = beads[(item.y & 0x3FFFFFu) + j];
  let bw = bitcast<u32>(b.w);
  let kind = bw & 255u;
  var r = u.lod.w;
  if (kind == 1u) { r = u.misc.y; }
  let id = ID_INST | item.x;
  let pl = f32((bw >> 8u) & 255u);
  let col = instColor(I, T, kind, pl, id);
  var bp = b.xyz;
  if (u.caps.y > 0u && T.nConf > 0u) { bp = b.xyz + modeDisp(T, (item.y & 0x3FFFFFu) + j - T.beadOff, I.seed); }
  else if ((typeFlags[I.t] & 8u) == 0u) { bp = flexBead(b.xyz, select(255.0, pl, kind == 0u), I.seed); }
  return emitQuad(I.p + qrot(I.q, bp), r, corner, col, id);
}
@fragment
fn fsMain(in: VOut) -> FOut { return shadeFrag(in); }
`,Zt=ut+`
@group(1) @binding(0) var<storage, read> stSpheres: array<vec4f>; // xyz r
@group(1) @binding(1) var<storage, read> list01: array<vec2u>;
@vertex
fn vsL1(@builtin(vertex_index) vi: u32, @builtin(instance_index) ii: u32) -> VOut {
  let item = list01[u.caps.x + ii];
  let j = vi / 6u; let corner = vi % 6u;
  let cnt = (item.y >> 22u) + 1u;
  if (j >= cnt) { var o: VOut; o.pos = vec4f(2.0, 2.0, 2.0, 1.0); return o; }
  let I0 = decodeInst(item.x); let T = types[I0.t];
  let I = animInst(I0, T.br);
  let s = stSpheres[(item.y & 0x3FFFFFu) + j];
  let id = ID_INST | item.x;
  // sphere-tree spheres are padded for collision; shrink slightly for display
  return emitQuad(I.p + qrot(I.q, s.xyz), s.w * 0.9, corner, instColor(I, T, 2u, -1.0, id), id);
}
@fragment
fn fsMain(in: VOut) -> FOut { return shadeFrag(in); }
`,en=ut+`
@group(1) @binding(0) var<storage, read> listL2: array<u32>;
@vertex
fn vsL2(@builtin(vertex_index) vi: u32, @builtin(instance_index) ii: u32) -> VOut {
  let i = listL2[ii];
  let I0 = decodeInst(i); let T = types[I0.t];
  let I = animInst(I0, T.br);
  let id = ID_INST | i;
  return emitQuad(I.p, T.req, vi, instColor(I, T, 2u, -1.0, id), id);
}
@fragment
fn fsMain(in: VOut) -> FOut { return shadeFrag(in); }
`,tn=ve+Ct+`
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

// free sphere record: xyz, w = bits: kind (4) | chain (15) << 4
@vertex
fn vsFree(@builtin(vertex_index) vi: u32, @builtin(instance_index) ii: u32) -> VOut {
  let i = slist[ii];
  let rec = spool[i];
  let w = bitcast<u32>(rec.w);
  let kind = w & 15u;
  let chain = (w >> 4u) & 0x7FFFu;
  var r = 0.6;
  if (kind == 1u) { r = 1.0; } else if (kind == 2u) { r = 0.75; } else if (kind == 5u) { r = 0.45; }
  let hj = f32(hash32(chain * 747796405u + kind) & 255u) / 255.0;
  let col = palette[48u + kind].xyz * (0.94 + 0.12 * hj);
  let slot = spages[i / SPAGE].y;
  let id = ID_FREE | ((slot & 0xFFFu) << 18u) | ((kind & 7u) << 15u) | chain;
  var p = rec.xyz;
  if (u.caps.y > 0u) {     // illustrative motion: FG repeats and mRNA most, wall polysaccharides and DNA little
    var amps = array<f32, 6>(0.25, 0.08, 0.06, 0.25, 0.25, 0.5);
    let amp = amps[min(kind, 5u)] * f32(u.caps.y) / 1000.0;
    let t = u.cut.w; let ph = 6.2832 * f32(hash32(chain * 2654435761u + kind) & 1023u) / 1023.0;
    p += amp * vec3f(sin(dot(p, vec3f(0.9, 0.3, -0.2)) + 1.1 * t + ph), sin(dot(p, vec3f(-0.2, 0.9, 0.4)) + 0.9 * t + 2.0 * ph),
                     sin(dot(p, vec3f(0.3, -0.3, 0.9)) + 1.3 * t + 3.0 * ph));
  }
  return emitQuad(p - u.eye.xyz, r, vi, col, id);
}

@fragment
fn fsMain(in: VOut) -> FOut { return shadeFrag(in); }
`,nn=ve+`
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
`,sn=ve+`
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
`,on=`
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
`,rn=ve+`
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
`,D=(i,e)=>[i[0]+e[0],i[1]+e[1],i[2]+e[2]],ce=(i,e)=>[i[0]-e[0],i[1]-e[1],i[2]-e[2]],E=(i,e)=>[i[0]*e,i[1]*e,i[2]*e],$=(i,e)=>i[0]*e[0]+i[1]*e[1]+i[2]*e[2],pe=(i,e)=>[i[1]*e[2]-i[2]*e[1],i[2]*e[0]-i[0]*e[2],i[0]*e[1]-i[1]*e[0]],ye=i=>Math.hypot(i[0],i[1],i[2]),K=i=>{const e=ye(i)||1;return[i[0]/e,i[1]/e,i[2]/e]},at=(i,e,t)=>[i[0]+(e[0]-i[0])*t,i[1]+(e[1]-i[1])*t,i[2]+(e[2]-i[2])*t];function De(i,e,t){const n=Math.cos(t),s=Math.sin(t),o=pe(e,i),r=$(e,i)*(1-n);return[i[0]*n+o[0]*s+e[0]*r,i[1]*n+o[1]*s+e[1]*r,i[2]*n+o[2]*s+e[2]*r]}function an(i,e){const t=new Float32Array(16);for(let n=0;n<4;n++)for(let s=0;s<4;s++){let o=0;for(let r=0;r<4;r++)o+=i[r*4+s]*e[n*4+r];t[n*4+s]=o}return t}function ln(i){const e=Array.from(i),t=new Float32Array(16),[n,s,o,r,c,u,l,d,p,m,y,g,P,x,a,h]=e,b=n*u-s*c,A=n*l-o*c,C=n*d-r*c,T=s*l-o*u,M=s*d-r*u,k=o*d-r*l,N=p*x-m*P,O=p*a-y*P,q=p*h-g*P,j=m*a-y*x,H=m*h-g*x,W=y*h-g*a;let U=b*W-A*H+C*j+T*q-M*O+k*N;return U&&(U=1/U,t[0]=(u*W-l*H+d*j)*U,t[1]=(o*H-s*W-r*j)*U,t[2]=(x*k-a*M+h*T)*U,t[3]=(y*M-m*k-g*T)*U,t[4]=(l*q-c*W-d*O)*U,t[5]=(n*W-o*q+r*O)*U,t[6]=(a*C-P*k-h*A)*U,t[7]=(p*k-y*C+g*A)*U,t[8]=(c*H-u*q+d*N)*U,t[9]=(s*q-n*H-r*N)*U,t[10]=(P*M-x*C+h*b)*U,t[11]=(m*C-p*M-g*b)*U,t[12]=(u*O-c*j-l*N)*U,t[13]=(n*j-s*O+o*N)*U,t[14]=(x*A-P*T-a*b)*U,t[15]=(p*T-m*A+y*b)*U),t}function ft(i){const e=parseInt(i.replace("#",""),16);return(e>>16&255|(e>>8&255)<<8|(e&255)<<16|4278190080)>>>0}function cn(i){return i<.5?4*i*i*i:1-Math.pow(-2*i+2,3)/2}const Et=new URL("data/",document.baseURI).href,dn=fetch(Et+"version.json",{cache:"no-store"}).then(i=>i.ok?i.json():{v:""}).then(i=>String(i.v??"")).catch(()=>""),un=i=>i?`?v=${encodeURIComponent(i)}`:"";function Lt(i){const e=atob(i),t=new Uint8Array(e.length);for(let n=0;n<e.length;n++)t[n]=e.charCodeAt(n);return t}async function Ze(i){const e=new Uint8Array(i,0,Math.min(2,i.byteLength));if(e[0]===31&&e[1]===139){const t=new DecompressionStream("gzip"),n=new Blob([i]).stream().pipeThrough(t);return await new Response(n).arrayBuffer()}return i}async function Pe(i,e){const t=window.__YC_EMBED__;if(t&&t[i]!==void 0)return Ze(Lt(t[i]).buffer);const n=await fetch(Et+i+un(await dn));if(!n.ok)throw new Error(`HTTP ${n.status} for ${i}`);if(!e||!n.body)return Ze(await n.arrayBuffer());const s=Number(n.headers.get("content-length")||0),o=n.body.getReader(),r=[];let c=0;for(;;){const{done:d,value:p}=await o.read();if(d)break;r.push(p),c+=p.length,e(c,s)}const u=new Uint8Array(c);let l=0;for(const d of r)u.set(d,l),l+=d.length;return Ze(u.buffer)}async function ze(i){const e=await Pe(i);return JSON.parse(new TextDecoder().decode(e))}const Ue=["cytosol","nucleoplasm","nucleolus","ER lumen","vacuole lumen","mitochondrial matrix","mitochondrial IMS","Golgi / endosome","peroxisome","cell wall","PM peripheral","NE peripheral","ER peripheral","lipid droplet surface","plasma membrane","ER membrane","nuclear envelope","vacuolar membrane","mito outer membrane","mito inner membrane","Golgi membrane","peroxisomal membrane","chromatin (nucleosome)","nuclear pore"],gt=["plasma membrane","ER membrane","nuclear envelope","vacuolar membrane","mito outer membrane","mito inner membrane","Golgi/endosome membrane","peroxisomal membrane","lipid droplet monolayer","cell surface"],hn=["exterior","cell wall (mannoprotein layer)","cell wall (glucan/chitin layer)","cytosol","cortical ER lumen","perinuclear space","nucleoplasm","nucleolus","vacuole","mito intermembrane space","mito matrix","lipid droplet","Golgi / endosome","peroxisome","nuclear pore channel"],_={0:"#7EA7D3",1:"#B993C9",2:"#D39BC0",3:"#A3C48C",4:"#9FC7C0",5:"#D9946B",6:"#E8B48F",7:"#B7C77E",8:"#C5B0E0",9:"#9A8352",10:"#7EA7D3",11:"#8FA3D0",12:"#7EA7D3",13:"#7EA7D3",14:"#5FAE75",15:"#7DBE86",16:"#6BB596",17:"#5BA8A0",18:"#C27C55",19:"#B5653F",20:"#93B060",21:"#A68BC9",22:"#F4B26A",23:"#3E9C9A",32:"#F3C643",33:"#F0BC3C",34:"#F6D15E",35:"#EAB542",36:"#EFA83A",37:"#E7A93C",38:"#F2CA55",40:"#F3C643",41:"#F8E3A0",42:"#E7A93C",48:"#B0437A",49:"#E07636",50:"#CDA866",51:"#C3C27A",52:"#C9558F",53:"#7CC6C2",56:"#5FAE75",57:"#7DBE86",58:"#6BB596",59:"#5BA8A0",60:"#C27C55",61:"#B5653F",62:"#93B060",63:"#A68BC9",64:"#E9C46A",65:"#9A8352",80:"#F6F1E7",81:"#F1E9D3",82:"#EEE3CC",83:"#E4ECF5",84:"#E8F0DF",85:"#E8F0DF",86:"#F1E7F1",87:"#EBD9EA",88:"#E3F0EE",89:"#F6E6DA",90:"#F2DCCB",91:"#F7EDC4",92:"#EEF1DC",93:"#EEE8F6",94:"#ECE8F2"};function pn(){const i=new Float32Array(384);for(const[e,t]of Object.entries(_)){const n=parseInt(t.slice(1),16),s=Number(e)*4;i[s]=(n>>16&255)/255,i[s+1]=(n>>8&255)/255,i[s+2]=(n&255)/255,i[s+3]=1}return i}const Ft={ribosome_80S:{c0:"#EE9BC0",c1:"#C9558F",cls:"ribosome"},tRNA:{c0:"#E27FA8",c1:"#E27FA8",cls:"tRNA"},nucleosome:{c0:"#F4B26A",c1:"#E07636",cls:"nucleosome"},NPC:{c0:"#3E9C9A",cls:"NPC"},PolI:{c0:"#9A72B8",c1:"#9A72B8",cls:"polymerase"},PolII:{c0:"#9A72B8",c1:"#9A72B8",cls:"polymerase"},PolIII:{c0:"#9A72B8",c1:"#9A72B8",cls:"polymerase"},FAS:{c0:"#8B88C9",cls:"complex"},proteasome_26S:{c0:"#8B88C9",cls:"complex"},CCT:{c0:"#8B88C9",cls:"complex"},F_actin:{c0:"#E0893A",c1:"#E0893A",cls:"actin"},Arp23:{c0:"#D9534F",c1:"#D9534F",cls:"actin"},TPM1_coil:{c0:"#2FA4A0",c1:"#2FA4A0",cls:"actin"},CP_cap:{c0:"#8E6CB0",c1:"#8E6CB0",cls:"actin"},MT_tubulin:{c0:"#74C476",c1:"#74C476",cls:"microtubule"},gTuRC_nuc:{c0:"#2E8B57",c1:"#2E8B57",cls:"microtubule"},gTuRC_cyto:{c0:"#2E8B57",c1:"#2E8B57",cls:"microtubule"},ribosome_80S_Ssh1:{c0:"#EE9BC0",c1:"#C9558F",cls:"ribosome"},ribosome_80S_Sec61:{c0:"#EE9BC0",c1:"#C9558F",cls:"ribosome"}};async function mn(i){i("Loading structure index",.02);const e=await ze("library/types.json"),t=e.types,n=e.bead_q;i("Loading residue coordinates (2.5 M beads)",.06);const s=await Pe("library/beads.bin",(M,k)=>i("Loading residue coordinates (2.5 M beads)",.06+.5*(k?M/k:0))),o=new Float32Array(await Pe("library/spheretree.bin"));let r=new Uint32Array(4),c={},u=0;try{const M=await ze("library/modes.json");r=new Uint32Array(await Pe("library/modes.bin")),c=M.types,u=M.cap??1.5}catch{}i("Building level-of-detail clusters",.6);const l=s.byteLength/8,d=new Int16Array(s),p=new Uint8Array(s),m=new Float32Array(l*4),y=new Uint32Array(m.buffer);for(let M=0;M<l;M++)m[4*M]=d[4*M]*n,m[4*M+1]=d[4*M+1]*n,m[4*M+2]=d[4*M+2]*n,y[4*M+3]=p[8*M+6]|p[8*M+7]<<8;let g=0;for(const M of t)g+=Math.ceil(M.nb/64);const P=new ArrayBuffer(g*32),x=new Float32Array(P),a=new Uint32Array(P),h=new ArrayBuffer(t.length*64),b=new Uint32Array(h),A=new Float32Array(h);let C=0;t.forEach((M,k)=>{const N=C;for(let H=0;H<M.nb;H+=64){const W=Math.min(64,M.nb-H),U=M.bo+H;let Ce=1e9,Ee=1e9,Le=1e9,Fe=-1e9,Se=-1e9,Ie=-1e9;for(let V=0;V<W;V++){const Y=4*(U+V),ie=m[Y],Z=m[Y+1],se=m[Y+2];ie<Ce&&(Ce=ie),ie>Fe&&(Fe=ie),Z<Ee&&(Ee=Z),Z>Se&&(Se=Z),se<Le&&(Le=se),se>Ie&&(Ie=se)}const Ge=(Ce+Fe)/2,ee=(Ee+Se)/2,me=(Le+Ie)/2;let G=0;for(let V=0;V<W;V++){const Y=4*(U+V);G=Math.max(G,Math.hypot(m[Y]-Ge,m[Y+1]-ee,m[Y+2]-me))}x[8*C]=Ge,x[8*C+1]=ee,x[8*C+2]=me,x[8*C+3]=G+.85+(c[k]?.[3]??0),a[8*C+4]=U,a[8*C+5]=W,C++}const O=Ft[M.id],q=M.src==="PDB";b[16*k+0]=M.bo,b[16*k+1]=M.nb,b[16*k+2]=N,b[16*k+3]=C-N,b[16*k+4]=M.so,b[16*k+5]=M.sn,A[16*k+6]=M.br,A[16*k+7]=Math.sqrt(5/3)*M.rg,b[16*k+8]=O?ft(O.c0):0,b[16*k+9]=O&&O.c1?ft(O.c1):0,b[16*k+10]=((q?1:0)|Math.min(M.copies,16777215)<<8)>>>0,A[16*k+11]=q?255:M.plddt_mean??255;const j=c[k];j&&(b[16*k+12]=j[0],b[16*k+13]=j[1],b[16*k+14]=j[2])});const T=new Map;return t.forEach((M,k)=>T.set(M.id,k)),{types:t,beads:m,st:o,clusters:P,typeInfo:h,modes:r,modeCap:u,nClusters:g,byId:T}}const lt=(i,e,t)=>((i+11)*22+(e+11))*22+(t+11);async function fn(){const i=await ze("tiles/index.json"),e=i.tile_columns,t=i.tiles,n=t.length,s=c=>e.indexOf(c),o=(c,u)=>{const l=new c(n),d=s(u);for(let p=0;p<n;p++)l[p]=t[p][d];return l},r={n,ix:o(Int32Array,"ix"),iy:o(Int32Array,"iy"),iz:o(Int32Array,"iz"),nInst:o(Uint32Array,"n_inst"),nGhost:o(Uint32Array,"n_ghost"),nMr:o(Uint32Array,"n_mrna_beads"),nLk:o(Uint32Array,"n_link_beads"),nVert:o(Uint32Array,"n_vert"),nTri:o(Uint32Array,"n_tri"),gz:o(Uint32Array,"gz_bytes"),labelMask:o(Uint32Array,"label_mask"),volGlucan:o(Float32Array,"vol_glucan_nm3"),volMannan:o(Float32Array,"vol_mannan_nm3"),key:new Map,wall:i.wall,regionNames:i.region_names};for(let c=0;c<n;c++)r.key.set(lt(r.ix[c],r.iy[c],r.iz[c]),c);return r}const gn=1073741824,B=GPUShaderStage,F=GPUBufferUsage;class bn{device;ctx;format;canvas;limits={};adapterInfo="";instPages=0;lipPages=0;freePages=0;capL0=0;capL1=0;capL2=0;capLip=0;capFree=0;maxSlots=4096;uni;inst;pageOrg;typesBuf;typeFlags;palette;beads;st;clusters;modes;modeCap=0;counters;list01;listL2;args;cullParams;caps;lipPool;lipPagesBuf;lipList;lipParams;freePool;freePagesBuf;freeList;freeParams;meshUni;tileState;fillBufs=[];pickBuf;pickParams;pickOut;statsBuf;pipes={};bg={};bgl={};colorTex;idTex;depthTex;width=0;height=0;usedInstPages=0;usedLipPages=0;usedFreePages=0;pickPending=null;pickBusy=!1;statsBusy=!1;stats={l0:0,l1:0,l2:0,lip:0,fre:0,ninst:0};lastStatsT=0;organelleMeshes=[];onLost=()=>{};constructor(e){this.canvas=e}async init(e,t){if(!("gpu"in navigator))throw new Error("WebGPU is not available in this browser.");const n=await navigator.gpu.requestAdapter({powerPreference:"high-performance"});if(!n)throw new Error("No WebGPU adapter found (GPU blocked or unsupported).");try{const p=n.info||{};this.adapterInfo=[p.vendor,p.architecture,p.device,p.description].filter(Boolean).join(" ")}catch{}const s=n.limits,o={maxStorageBufferBindingSize:Math.min(s.maxStorageBufferBindingSize,2**31-4),maxBufferSize:Math.min(s.maxBufferSize,2**31-4),maxComputeWorkgroupsPerDimension:s.maxComputeWorkgroupsPerDimension};this.device=await n.requestDevice({requiredLimits:o}),this.device.lost.then(p=>this.onLost(`GPU device lost: ${p.message}`)),this.device.addEventListener("uncapturederror",p=>console.error("WebGPU error:",p.error?.message)),this.limits={...o};const r=this.device.limits.maxStorageBufferBindingSize,c=p=>Math.min(p,Math.floor(r/16));this.instPages=Math.floor(c(8*1024*1024)/X),this.lipPages=Math.floor(c(8*1024*1024)/Q),this.freePages=Math.floor(c(10*1024*1024)/Q),this.capL0=Math.min(4*1024*1024,Math.floor(r/8)-2*1024*1024),this.capL1=2*1024*1024,this.capL2=Math.min(8*1024*1024,Math.floor(r/4)),this.capLip=this.lipPages*Q,this.capFree=this.freePages*Q,this.ctx=this.canvas.getContext("webgpu"),this.format=navigator.gpu.getPreferredCanvasFormat(),this.ctx.configure({device:this.device,format:this.format,alphaMode:"opaque",usage:GPUTextureUsage.RENDER_ATTACHMENT|GPUTextureUsage.COPY_SRC});const u=this.device,l=(p,m,y)=>u.createBuffer({size:Math.max(16,Math.ceil(p/16)*16),usage:m,label:y}),d=F.STORAGE|F.COPY_DST;this.uni=l(512,F.UNIFORM|F.COPY_DST,"uniforms"),this.inst=l(this.instPages*X*16,d,"instances"),this.pageOrg=l(this.instPages*16,d,"pageOrg"),this.typesBuf=l(e.typeInfo.byteLength,d,"types"),this.typeFlags=l(e.types.length*4,d,"typeFlags"),this.palette=l(1536,d,"palette"),this.beads=l(e.beads.byteLength,d,"beads"),this.modes=l(Math.max(16,e.modes.byteLength),d,"modes"),this.modeCap=e.modeCap,this.st=l(e.st.byteLength,d,"spheretree"),this.clusters=l(e.clusters.byteLength,d,"clusters"),u.queue.writeBuffer(this.typesBuf,0,e.typeInfo),u.queue.writeBuffer(this.palette,0,pn()),u.queue.writeBuffer(this.beads,0,e.beads),u.queue.writeBuffer(this.modes,0,e.modes),u.queue.writeBuffer(this.st,0,e.st),u.queue.writeBuffer(this.clusters,0,e.clusters),this.counters=l(32,F.STORAGE|F.COPY_DST|F.COPY_SRC,"counters"),this.list01=l((this.capL0+this.capL1)*8,F.STORAGE,"list01"),this.listL2=l(this.capL2*4,F.STORAGE,"listL2"),this.args=l(80,F.STORAGE|F.INDIRECT,"args"),this.cullParams=l(32,F.UNIFORM|F.COPY_DST,"cullParams"),this.caps=l(32,F.UNIFORM|F.COPY_DST,"caps"),u.queue.writeBuffer(this.caps,0,new Uint32Array([this.capL0,this.capL1,this.capL2,this.capLip,this.capFree,0,0,0])),this.lipPool=l(this.lipPages*Q*16,d,"lipidPool"),this.lipPagesBuf=l(this.lipPages*16,d,"lipidPages"),this.lipList=l(this.capLip*4,F.STORAGE,"lipidList"),this.lipParams=l(16,F.UNIFORM|F.COPY_DST,"lipParams"),this.freePool=l(this.freePages*Q*16,d,"freePool"),this.freePagesBuf=l(this.freePages*16,d,"freePages"),this.freeList=l(this.capFree*4,F.STORAGE,"freeList"),this.freeParams=l(16,F.UNIFORM|F.COPY_DST,"freeParams"),this.meshUni=l(256*(this.maxSlots+16),F.UNIFORM|F.COPY_DST,"meshUniforms"),this.tileState=l(484*22*4,d,"tileState"),this.pickBuf=u.createBuffer({size:16,usage:F.MAP_READ|F.COPY_DST}),this.pickParams=l(16,F.UNIFORM|F.COPY_DST,"pickParams"),this.pickOut=l(16,F.STORAGE|F.COPY_SRC,"pickOut"),this.statsBuf=u.createBuffer({size:32,usage:F.MAP_READ|F.COPY_DST}),this.initFill(t),await this.initPipelines(),this.resize()}initFill(e){const t=this.device,n=e.PT.length,s=new ArrayBuffer(n*48),o=new Float32Array(s),r=new Uint32Array(s);for(let a=0;a<n;a++){for(let h=0;h<8;h++)o[12*a+h]=e.PQ[8*a+h];r[12*a+8]=e.PT[a],r[12*a+9]=e.PC[a]}const c=e.pores.length/3,u=new Float32Array(Math.max(1,c)*4);for(let a=0;a<c;a++)u[4*a]=e.pores[3*a],u[4*a+1]=e.pores[3*a+1],u[4*a+2]=e.pores[3*a+2];const l=e.grid,d=l.prim_idx.length,p=new Uint32Array(l.prim_off.length+l.pore_off.length);p.set(l.prim_off,0);for(let a=0;a<l.pore_off.length;a++)p[l.prim_off.length+a]=l.pore_off[a]+d;const m=new Uint32Array(d+l.pore_idx.length+1);m.set(l.prim_idx,0),m.set(l.pore_idx,d);const y=e.G,g=e.consts,P=new Float32Array([y[0],y[1],y[2],y[3],y[4],y[5],y[6],y[7],y[8],y[9],y[10],0,g.PORE_R,g.PERI_W,g.IMS_W,c,l.n,l.ext,0,0]),x=(a,h)=>{const b=t.createBuffer({size:Math.ceil(a.byteLength/16)*16+16,usage:h|F.COPY_DST});return t.queue.writeBuffer(b,0,a),b};this.fillBufs=[x(s,F.STORAGE),x(u,F.STORAGE),x(p,F.STORAGE),x(m,F.STORAGE),x(P,F.UNIFORM)]}async initPipelines(){const e=this.device,t=(G,V)=>({binding:G,visibility:V,buffer:{type:"read-only-storage"}}),n=G=>({binding:G,visibility:B.COMPUTE,buffer:{type:"storage"}}),s=(G,V,Y=!1)=>({binding:G,visibility:V,buffer:{type:"uniform",hasDynamicOffset:Y}}),o=B.COMPUTE|B.VERTEX|B.FRAGMENT,r=B.COMPUTE|B.VERTEX,c=B.VERTEX|B.FRAGMENT,u=G=>e.createBindGroupLayout({entries:G});this.bgl.g0=u([s(0,o),t(1,r),t(2,r),t(3,r),t(4,r),t(5,c)]),this.bgl.cull=u([n(0),n(1),n(2),s(3,B.COMPUTE),t(4,B.COMPUTE)]),this.bgl.cullS=u([n(0),t(1,B.COMPUTE),t(2,B.COMPUTE),n(3),s(4,B.COMPUTE)]),this.bgl.fin=u([t(0,B.COMPUTE),n(1),s(2,B.COMPUTE)]),this.bgl.l0=u([t(0,B.VERTEX),t(1,B.VERTEX),t(2,B.VERTEX)]),this.bgl.l2=u([t(0,B.VERTEX)]),this.bgl.sph=u([t(0,B.VERTEX),t(1,B.VERTEX),t(2,B.VERTEX)]),this.bgl.mesh=u([s(0,c,!0)]),this.bgl.fill=u([t(0,B.FRAGMENT),t(1,B.FRAGMENT),t(2,B.FRAGMENT),t(3,B.FRAGMENT),s(4,B.FRAGMENT),t(5,B.FRAGMENT)]),this.bgl.comp=u([{binding:0,visibility:B.FRAGMENT,texture:{sampleType:"unfilterable-float"}},{binding:1,visibility:B.FRAGMENT,texture:{sampleType:"uint"}},{binding:2,visibility:B.FRAGMENT,texture:{sampleType:"depth"}}]),this.bgl.pick=u([{binding:0,visibility:B.COMPUTE,texture:{sampleType:"uint"}},{binding:1,visibility:B.COMPUTE,texture:{sampleType:"depth"}},s(2,B.COMPUTE),n(3)]);const l=(...G)=>e.createPipelineLayout({bindGroupLayouts:G}),d=async(G,V)=>{const Y=e.createShaderModule({code:G,label:V}),ie=await Y.getCompilationInfo(),Z=ie.messages.filter(se=>se.type==="error");for(const se of ie.messages)(se.type==="error"?console.error:console.warn)(`[${V}] ${se.lineNum}:${se.linePos} ${se.message}`);if(Z.length)throw new Error(`WGSL compile error in ${V}: ${Z[0].lineNum}:${Z[0].linePos} ${Z[0].message}`);return Y},p=await d(Yt,"cull"),m=await d(Xt,"cullSpheres"),y=await d(Qt,"finalize"),g=await d(Jt,"drawL0"),P=await d(Zt,"drawL1"),x=await d(en,"drawL2"),a=await d(tn,"drawSpheres"),h=await d(nn,"mesh"),b=await d(sn,"fill"),A=await d(rn,"composite"),C=await d(on,"pick"),T=(G,V,Y)=>e.createComputePipelineAsync({layout:Y,compute:{module:G,entryPoint:V}}),M=[{format:"rgba8unorm"},{format:"r32uint"}],k=(G="greater")=>({format:"depth32float",depthWriteEnabled:!0,depthCompare:G}),N=(G,V,Y,ie,Z={})=>e.createRenderPipelineAsync({layout:ie,vertex:{module:G,entryPoint:V,...Z.vertex?{buffers:Z.vertex.buffers}:{}},fragment:{module:G,entryPoint:Y,targets:M},primitive:{topology:"triangle-list",cullMode:"none"},depthStencil:Z.depthStencil??k()}),O=this.bgl.g0,[q,j,H,W,U,Ce,Ee,Le,Fe,Se,Ie,Ge]=await Promise.all([T(p,"cullInstances",l(O,this.bgl.cull)),T(m,"cullSpheres",l(O,this.bgl.cullS)),T(y,"main",l(this.bgl.fin)),N(g,"vsL0","fsMain",l(O,this.bgl.l0)),N(P,"vsL1","fsMain",l(O,this.bgl.l0)),N(x,"vsL2","fsMain",l(O,this.bgl.l2)),N(a,"vsLipid","fsMain",l(O,this.bgl.sph)),N(a,"vsFree","fsMain",l(O,this.bgl.sph)),N(h,"vsMesh","fsMesh",l(O,this.bgl.mesh),{vertex:{buffers:[{arrayStride:32,attributes:[{shaderLocation:0,offset:0,format:"float32x4"},{shaderLocation:1,offset:16,format:"float32x4"}]}]}}),N(b,"vsFill","fsFill",l(O,this.bgl.fill),{depthStencil:k("always")}),e.createRenderPipelineAsync({layout:l(O,this.bgl.comp),vertex:{module:A,entryPoint:"vsFull"},fragment:{module:A,entryPoint:"fsComposite",targets:[{format:this.format}]},primitive:{topology:"triangle-list"}}),T(C,"main",l(this.bgl.pick))]);Object.assign(this.pipes,{cull:q,cullS:j,fin:H,l0:W,l1:U,l2:Ce,lip:Ee,fre:Le,mesh:Fe,fill:Se,comp:Ie,pick:Ge});const ee=(G,V,Y=[])=>e.createBindGroup({layout:G,entries:V.map((ie,Z)=>({binding:Z,resource:{buffer:ie,size:Y[Z]}}))});this.bg.g0=ee(O,[this.uni,this.inst,this.pageOrg,this.typesBuf,this.typeFlags,this.palette]),this.bg.cull=ee(this.bgl.cull,[this.counters,this.list01,this.listL2,this.cullParams,this.clusters]),this.bg.cullLip=ee(this.bgl.cullS,[this.counters,this.lipPool,this.lipPagesBuf,this.lipList,this.lipParams]),this.bg.cullFree=ee(this.bgl.cullS,[this.counters,this.freePool,this.freePagesBuf,this.freeList,this.freeParams]),this.bg.fin=ee(this.bgl.fin,[this.counters,this.args,this.caps]),this.bg.l0=ee(this.bgl.l0,[this.beads,this.list01,this.modes]),this.bg.l1=ee(this.bgl.l0,[this.st,this.list01,this.modes]),this.bg.l2=ee(this.bgl.l2,[this.listL2]),this.bg.lip=ee(this.bgl.sph,[this.lipPool,this.lipPagesBuf,this.lipList]),this.bg.fre=ee(this.bgl.sph,[this.freePool,this.freePagesBuf,this.freeList]),this.bg.mesh=e.createBindGroup({layout:this.bgl.mesh,entries:[{binding:0,resource:{buffer:this.meshUni,size:16}}]});const me=this.fillBufs;this.bg.fill=ee(this.bgl.fill,[me[0],me[1],me[2],me[3],me[4],this.tileState])}resize(e=1){const t=Math.min(window.devicePixelRatio||1,2)*e,n=Math.max(1,Math.floor(this.canvas.clientWidth*t)),s=Math.max(1,Math.floor(this.canvas.clientHeight*t));if(n===this.width&&s===this.height&&this.colorTex)return;this.canvas.width=n,this.canvas.height=s,this.width=n,this.height=s,this.colorTex?.destroy(),this.idTex?.destroy(),this.depthTex?.destroy();const o=GPUTextureUsage;this.colorTex=this.device.createTexture({size:[n,s],format:"rgba8unorm",usage:o.RENDER_ATTACHMENT|o.TEXTURE_BINDING}),this.idTex=this.device.createTexture({size:[n,s],format:"r32uint",usage:o.RENDER_ATTACHMENT|o.TEXTURE_BINDING|o.COPY_SRC}),this.depthTex=this.device.createTexture({size:[n,s],format:"depth32float",usage:o.RENDER_ATTACHMENT|o.TEXTURE_BINDING|o.COPY_SRC}),this.bg.comp=this.device.createBindGroup({layout:this.bgl.comp,entries:[{binding:0,resource:this.colorTex.createView()},{binding:1,resource:this.idTex.createView()},{binding:2,resource:this.depthTex.createView()}]}),this.bg.pick=this.device.createBindGroup({layout:this.bgl.pick,entries:[{binding:0,resource:this.idTex.createView()},{binding:1,resource:this.depthTex.createView()},{binding:2,resource:{buffer:this.pickParams}},{binding:3,resource:{buffer:this.pickOut}}]})}writeUniforms(e,t,n){const s=this.width/this.height,o=e.viewRot(),r=e.proj(s),c=an(r,o),u=ln(c),l=e.eye,d=new ArrayBuffer(512),p=new Float32Array(d),m=new Uint32Array(d);p.set(r,0),p.set(o,16),p.set(u,32);const y=e.pxPerNmAtUnit(this.height);p.set([l[0],l[1],l[2],y],48);const g=[t.plane[0],t.plane[1],t.plane[2]];p.set([g[0],g[1],g[2],t.plane[3]],52);const P=[o[0]*g[0]+o[4]*g[1]+o[8]*g[2],o[1]*g[0]+o[5]*g[1]+o[9]*g[2],o[2]*g[0]+o[6]*g[1]+o[10]*g[2]];p.set([P[0],P[1],P[2],$(g,l)-t.plane[3]],56),p.set([t.cutOn?1:0,t.slab,Math.min(400,Math.max(20,t.slab)),performance.now()/1e3],60),p.set([t.l0px,t.l1px,t.minpx,.6],64),p.set([this.width,this.height,1/this.width,1/this.height],68),p.set(e.frustum(s),72),m.set([t.colorMode,t.regionMask>>>0,t.freeHide,t.highlight?1:0],96),m.set([t.selId>>>0,0,t.shaded?1:0,t.organelle?1:0],100),p.set([360/65535,.8,t.outline,Math.min(window.devicePixelRatio||1,2)*t.renderScale],104),p.set([n[0],n[1],n[2],e.near()],108),m.set([this.capL0,Math.round(t.motion*1e3),Math.round(this.modeCap*1e3),0],112),this.device.queue.writeBuffer(this.uni,0,d,0,464)}frame(e,t,n,s,o){const r=this.device;this.writeUniforms(e,t,n);const c=r.createCommandEncoder();c.clearBuffer(this.counters);const u=!t.organelle;if(u){const g=c.beginComputePass();g.setBindGroup(0,this.bg.g0);const P=(x,a,h)=>{const b=Math.ceil(x/64),A=Math.min(b,65535),C=Math.ceil(b/A);return r.queue.writeBuffer(a,0,h(A*64)),[A,C]};if(this.usedInstPages>0){const[x,a]=P(this.usedInstPages*X,this.cullParams,h=>new Uint32Array([this.usedInstPages,h,this.capL0,this.capL1,this.capL2,this.capLip,this.capFree,0]).buffer);g.setPipeline(this.pipes.cull),g.setBindGroup(1,this.bg.cull),g.dispatchWorkgroups(x,a)}if(this.usedLipPages>0&&t.showLipids){const[x,a]=P(this.usedLipPages*Q,this.lipParams,h=>new Uint32Array([this.usedLipPages,h,this.capLip,1]).buffer);g.setPipeline(this.pipes.cullS),g.setBindGroup(1,this.bg.cullLip),g.dispatchWorkgroups(x,a)}if(this.usedFreePages>0){const[x,a]=P(this.usedFreePages*Q,this.freeParams,h=>new Uint32Array([this.usedFreePages,h,this.capFree,0]).buffer);g.setPipeline(this.pipes.cullS),g.setBindGroup(1,this.bg.cullFree),g.dispatchWorkgroups(x,a)}g.setPipeline(this.pipes.fin),g.setBindGroup(0,this.bg.fin),g.dispatchWorkgroups(1),g.end()}const l=c.beginRenderPass({colorAttachments:[{view:this.colorTex.createView(),loadOp:"clear",storeOp:"store",clearValue:{r:.965,g:.945,b:.906,a:1}},{view:this.idTex.createView(),loadOp:"clear",storeOp:"store",clearValue:{r:0,g:0,b:0,a:0}}],depthStencilAttachment:{view:this.depthTex.createView(),depthLoadOp:"clear",depthStoreOp:"store",depthClearValue:0}});if(l.setBindGroup(0,this.bg.g0),l.setPipeline(this.pipes.fill),l.setBindGroup(1,this.bg.fill),l.draw(3),t.organelle){l.setPipeline(this.pipes.mesh);const g=this.maxSlots;r.queue.writeBuffer(this.meshUni,g*256,new Uint32Array([0,1,0,0]));for(const P of this.organelleMeshes)l.setBindGroup(1,this.bg.mesh,[g*256]),l.setVertexBuffer(0,P.vb),l.setIndexBuffer(P.ib,"uint32"),l.drawIndexed(P.nIdx,1)}else if(t.showMeshes&&s.length){l.setPipeline(this.pipes.mesh);for(const g of s)l.setBindGroup(1,this.bg.mesh,[g.slot*256]),l.setVertexBuffer(0,g.vb),l.setIndexBuffer(g.ib,"uint32"),l.drawIndexed(g.nIdx,2)}u&&(l.setPipeline(this.pipes.l2),l.setBindGroup(1,this.bg.l2),l.drawIndirect(this.args,32),l.setPipeline(this.pipes.l1),l.setBindGroup(1,this.bg.l1),l.drawIndirect(this.args,16),l.setPipeline(this.pipes.l0),l.setBindGroup(1,this.bg.l0),l.drawIndirect(this.args,0),t.showLipids&&(l.setPipeline(this.pipes.lip),l.setBindGroup(1,this.bg.lip),l.drawIndirect(this.args,48)),l.setPipeline(this.pipes.fre),l.setBindGroup(1,this.bg.fre),l.drawIndirect(this.args,64)),l.end();const d=this.ctx.getCurrentTexture(),p=c.beginRenderPass({colorAttachments:[{view:d.createView(),loadOp:"clear",storeOp:"store",clearValue:{r:1,g:1,b:1,a:1}}]});p.setPipeline(this.pipes.comp),p.setBindGroup(0,this.bg.g0),p.setBindGroup(1,this.bg.comp),p.draw(3),p.end();let m=null;if(this.pickPending&&!this.pickBusy){m=this.pickPending,this.pickPending=null;const g=Math.max(0,Math.min(this.width-1,Math.floor(m.x))),P=Math.max(0,Math.min(this.height-1,Math.floor(m.y)));r.queue.writeBuffer(this.pickParams,0,new Uint32Array([g,P,0,0]));const x=c.beginComputePass();x.setPipeline(this.pipes.pick),x.setBindGroup(0,this.bg.pick),x.dispatchWorkgroups(1),x.end(),c.copyBufferToBuffer(this.pickOut,0,this.pickBuf,0,16)}let y=!1;if(o&&!this.statsBusy&&u&&(c.copyBufferToBuffer(this.counters,0,this.statsBuf,0,32),y=!0),r.queue.submit([c.finish()]),m){this.pickBusy=!0;const g=m;this.pickBuf.mapAsync(GPUMapMode.READ).then(()=>{const P=this.pickBuf.getMappedRange(0,8),x=new Uint32Array(P)[0],a=new Float32Array(P)[1];this.pickBuf.unmap(),this.pickBusy=!1,g.resolve({id:x>>>0,depth:a})}).catch(()=>{this.pickBusy=!1})}y&&(this.statsBusy=!0,this.statsBuf.mapAsync(GPUMapMode.READ).then(()=>{const g=new Uint32Array(this.statsBuf.getMappedRange());this.stats={l0:g[0],l1:g[1],l2:g[2],lip:g[3],fre:g[4],ninst:g[5]},this.statsBuf.unmap(),this.statsBusy=!1}).catch(()=>{this.statsBusy=!1}))}pick(e,t){return new Promise(n=>{this.pickPending={x:e,y:t,resolve:n}})}setMeshSlot(e,t){const n=new ArrayBuffer(16),s=new Uint32Array(n),o=new Float32Array(n);s[0]=e,s[1]=t?2:0,o[2]=t?-1:-1.9,o[3]=t?1:1.9,this.device.queue.writeBuffer(this.meshUni,e*256,n)}uploadOrganelleMeshes(e){const t=this.device,n=new Uint32Array(e,0,1)[0],s=new Uint32Array(e,4,n*3);let o=4+n*12;for(const r of this.organelleMeshes)r.vb.destroy(),r.ib.destroy();this.organelleMeshes=[];for(let r=0;r<n;r++){const c=s[3*r],u=s[3*r+1],l=s[3*r+2],d=new Float32Array(e.slice(o,o+u*12));o+=u*12;const p=new Uint32Array(e.slice(o,o+l*12));o+=l*12;const m=new Float32Array(u*3);for(let a=0;a<l;a++){const h=p[3*a],b=p[3*a+1],A=p[3*a+2],C=d[3*b]-d[3*h],T=d[3*b+1]-d[3*h+1],M=d[3*b+2]-d[3*h+2],k=d[3*A]-d[3*h],N=d[3*A+1]-d[3*h+1],O=d[3*A+2]-d[3*h+2],q=T*O-M*N,j=M*k-C*O,H=C*N-T*k;for(const W of[h,b,A])m[3*W]+=q,m[3*W+1]+=j,m[3*W+2]+=H}const y=new Float32Array(u*8),g=new Uint32Array(y.buffer);for(let a=0;a<u;a++){const h=Math.hypot(m[3*a],m[3*a+1],m[3*a+2])||1;y[8*a]=d[3*a],y[8*a+1]=d[3*a+1],y[8*a+2]=d[3*a+2],g[8*a+3]=c,y[8*a+4]=m[3*a]/h,y[8*a+5]=m[3*a+1]/h,y[8*a+6]=m[3*a+2]/h}const P=t.createBuffer({size:y.byteLength,usage:F.VERTEX|F.COPY_DST});t.queue.writeBuffer(P,0,y);const x=t.createBuffer({size:Math.ceil(p.byteLength/4)*4,usage:F.INDEX|F.COPY_DST});t.queue.writeBuffer(x,0,p),this.organelleMeshes.push({vb:P,ib:x,nIdx:l*3,mem:c})}}gpuMemoryMB(){return[this.inst,this.pageOrg,this.beads,this.st,this.clusters,this.list01,this.listL2,this.lipPool,this.lipList,this.freePool,this.freeList,this.meshUni].reduce((t,n)=>t+n.size,0)/1048576+this.width*this.height*12/1048576}}class yn{target=[0,0,0];dist=9e3;fwd=[0,0,-1];up=[0,1,0];fovy=35*Math.PI/180;minDist=1.5;maxDist=14e3;anim=null;get right(){return K(pe(this.fwd,this.up))}get eye(){return ce(this.target,E(this.fwd,this.dist))}orthonormalize(){this.fwd=K(this.fwd);const e=K(pe(this.fwd,this.up));this.up=K(pe(e,this.fwd))}orbit(e,t){const n=this.right;this.fwd=De(this.fwd,this.up,-e),this.fwd=De(this.fwd,n,-t),this.up=De(this.up,n,-t),this.orthonormalize()}roll(e){this.up=De(this.up,this.fwd,e),this.orthonormalize()}pan(e,t,n){const s=1/n;this.target=D(this.target,D(E(this.right,-e*s),E(this.up,t*s)))}zoom(e){this.dist=Math.min(this.maxDist,Math.max(this.minDist,this.dist*e))}fly(e){this.target=D(this.target,e)}pxPerNmAtUnit(e){return e/(2*Math.tan(this.fovy/2))}distForWidth(e,t){return e/(2*Math.tan(this.fovy/2)*Math.min(t,1))}near(){return Math.min(5,Math.max(.02,this.dist*.004))}viewRot(){const e=this.right,t=this.up,n=this.fwd,s=new Float32Array(16);return s[0]=e[0],s[4]=e[1],s[8]=e[2],s[1]=t[0],s[5]=t[1],s[9]=t[2],s[2]=-n[0],s[6]=-n[1],s[10]=-n[2],s[15]=1,s}proj(e){const t=1/Math.tan(this.fovy/2),n=new Float32Array(16);return n[0]=t/e,n[5]=t,n[11]=-1,n[14]=this.near(),n}frustum(e){const t=Math.tan(this.fovy/2),n=t*e,s=this.right,o=this.up,r=this.fwd,c=(l,d=0)=>{const p=K(l);return[p[0],p[1],p[2],d]},u=[...c(D(s,E(r,n))),...c(D(E(s,-1),E(r,n))),...c(D(o,E(r,t))),...c(D(E(o,-1),E(r,t))),r[0],r[1],r[2],-this.near(),0,0,0,1];return new Float32Array(u)}ray(e,t,n){const s=Math.tan(this.fovy/2),o=s*n,r=(e*2-1)*o,c=(1-t*2)*s;return K(D(this.fwd,D(E(this.right,r),E(this.up,c))))}state(){return{target:[...this.target],dist:this.dist,fwd:[...this.fwd],up:[...this.up]}}set(e){this.target=[...e.target],this.dist=e.dist,this.fwd=[...e.fwd],this.up=[...e.up],this.orthonormalize()}flyTo(e,t=1800,n,s){this.anim={t0:performance.now(),dur:t,a:this.state(),b:e,onDone:n,onStep:s}}step(e){if(!this.anim)return!1;const t=this.anim;let n=Math.min(1,(e-t.t0)/t.dur);const s=cn(n),o=Math.log(t.a.dist),r=Math.log(t.b.dist),c=Math.hypot(...ce(t.b.target,t.a.target)),u=Math.log(Math.max(t.a.dist,t.b.dist,c*1.2)),l=Math.max(0,u-Math.max(o,r)),d=o+(r-o)*s+l*Math.sin(Math.PI*s);this.target=at(t.a.target,t.b.target,s),this.dist=Math.exp(d),this.fwd=wn(t.a.fwd,t.b.fwd,s);let p=t.a.up;return this.up=K(at(p,t.b.up,s)),Math.abs($(this.up,this.fwd))>.98&&(this.up=K(pe(pe(this.fwd,t.b.up),this.fwd))),this.orthonormalize(),t.onStep?.(s),n>=1&&(this.anim=null,t.onDone?.()),!0}}function wn(i,e,t){const n=Math.max(-1,Math.min(1,$(i,e))),s=Math.acos(n);if(s<1e-4)return K(at(i,e,t));if(s>Math.PI-.001){const r=K(pe(i,Math.abs(i[0])<.9?[1,0,0]:[0,1,0]));return De(i,r,Math.PI*t)}const o=Math.sin(s);return K(D(E(i,Math.sin((1-t)*s)/o),E(e,Math.sin(t*s)/o)))}const St='const kt=.005493247882810712,Jt=.005,ln=.01,fn=.70710678;let St="",Nt={},Ut,gt,zt,R,Gt=null,lt=[],Vt={glucanPerNm3:.018147,mannanPerProtein:1908.6},Lt="";const et=new Map,hn=48;async function Tt(s){if(Nt[s])return Yt(Nt[s].slice().buffer);const t=await fetch(St+s+(Lt?`?v=${encodeURIComponent(Lt)}`:""));if(!t.ok)throw new Error(`HTTP ${t.status} ${s}`);return Yt(await t.arrayBuffer())}async function Yt(s){const t=new Uint8Array(s,0,Math.min(2,s.byteLength));return t[0]===31&&t[1]===139?await new Response(new Blob([s]).stream().pipeThrough(new DecompressionStream("gzip"))).arrayBuffer():s}async function yn(s){St=s.base,lt=s.G;try{const l=await fetch(St+"version.json",{cache:"no-store"});l.ok&&(Lt=String((await l.json()).v??""))}catch{Lt=""}s.wall&&(Vt=s.wall),s.embedded&&(Nt=s.embedded);const t=JSON.parse(new TextDecoder().decode(await Tt("library/types.json"))),n=t.types,e=n.length,a=t.bead_q;Ut=new Uint32Array(e),gt=new Uint32Array(e),zt=new Float32Array(e);for(let l=0;l<e;l++)Ut[l]=n[l].bo,gt[l]=n[l].nb,zt[l]=n[l].br;const o=await Tt("library/beads.bin"),d=new Int16Array(o),h=o.byteLength/8;R=new Float32Array(h*3);for(let l=0;l<h;l++)R[3*l]=d[4*l]*a,R[3*l+1]=d[4*l+1]*a,R[3*l+2]=d[4*l+2]*a}const bt=s=>s+3&-4;function Qt(s,t,n,e,a,o,d,h,l,r,m){const p=new Uint32Array(s,t,n);t+=4*n;const I=new Int32Array(s,t,3*n);t+=12*n;const c=new Int8Array(s,t,3*(e-n));t+=3*(e-n),t=bt(t);let u=0,i=r;for(let f=0;f<n;f++){let b=I[3*f],A=I[3*f+1],v=I[3*f+2];const O=p[f];for(let L=0;L<O;L++)L>0&&(b+=c[3*u],A+=c[3*u+1],v+=c[3*u+2],u++),h[4*i]=a[0]+b*o,h[4*i+1]=a[1]+A*o,h[4*i+2]=a[2]+v*o,l[4*i+3]=((m?m[f]:d)|(f&32767)<<4)>>>0,i++}return{off:t,o:i}}function Zt(s){const t=new Uint32Array(s,0,16);if(t[0]!==827605849)throw new Error("bad tile magic");const n=new Int32Array(s,4,3),e=[n[0]*200,n[1]*200,n[2]*200],[a,o,d,h,l,r,m,p,I,c,u]=[t[4],t[5],t[6],t[7],t[8],t[9],t[10],t[11],t[12],t[13],t[14]];let i=64;const f=new Uint32Array(s.slice(i,i+16*a));i+=16*a;const b=new Uint32Array(s.slice(i,i+16*o));i+=16*o;const A=new Float32Array((h+r+u)*4),v=new Uint32Array(A.buffer);let O=0;if(d){const y=Qt(s,i,d,h,e,Jt,0,A,v,O);i=y.off,O=y.o}if(l){const y=Qt(s,i,l,r,e,ln,1,A,v,O);i=y.off,O=y.o}const L=new Float32Array(m*3);if(m){const y=new Uint16Array(s,i,3*m);for(let M=0;M<m;M++)for(let x=0;x<3;x++)L[3*M+x]=e[x]-80+y[3*M+x]*kt;i+=bt(6*m)}let g=new Uint32Array(p*3);p&&(I&1?g=new Uint32Array(s.slice(i,i+12*p)):g=Uint32Array.from(new Uint16Array(s,i,3*p)),i+=bt((I&1?4:2)*3*p));const P=p?new Uint8Array(s.slice(i,i+p)):new Uint8Array(0);if(p&&(i+=bt(p)),c){const y=new Uint8Array(s.slice(i,i+c));i+=bt(c),Qt(s,i,c,u,e,Jt,0,A,v,O,y)}const _=new Float32Array(m*3),w=new Uint8Array(m);for(let y=0;y<p;y++){const M=g[3*y],x=g[3*y+1],S=g[3*y+2],k=L[3*x]-L[3*M],U=L[3*x+1]-L[3*M+1],E=L[3*x+2]-L[3*M+2],Z=L[3*S]-L[3*M],q=L[3*S+1]-L[3*M+1],H=L[3*S+2]-L[3*M+2],K=U*H-E*q,J=E*Z-k*H,T=k*q-U*Z;for(const z of[M,x,S])_[3*z]+=K,_[3*z+1]+=J,_[3*z+2]+=T,w[z]=P[y]}for(let y=0;y<m;y++){const M=Math.hypot(_[3*y],_[3*y+1],_[3*y+2])||1;_[3*y]/=M,_[3*y+1]/=M,_[3*y+2]/=M}return{lo:e,inst:f,ghost:b,free:A,verts:L,normals:_,tris:g,fmem:P,vmem:w,nI:a,nG:o}}const $t=s=>(s<<16>>16)/32767*fn;function Dt(s,t,n,e){const a=s[4*t],o=s[4*t+1],d=s[4*t+2],h=s[4*t+3];e[0]=n[0]-80+(o&65535)*kt,e[1]=n[1]-80+(o>>>16)*kt,e[2]=n[2]-80+(d&65535)*kt;const l=$t(d>>>16),r=$t(h&65535),m=$t(h>>>16),p=Math.sqrt(Math.max(0,1-l*l-r*r-m*m)),I=a>>>14&3,c=I===0?[p,l,r,m]:I===1?[l,p,r,m]:I===2?[l,r,p,m]:[l,r,m,p],[u,i,f,b]=c;return e[3]=u*u+i*i-f*f-b*b,e[4]=2*(i*f-u*b),e[5]=2*(i*b+u*f),e[6]=2*(i*f+u*b),e[7]=u*u-i*i+f*f-b*b,e[8]=2*(f*b-u*i),e[9]=2*(i*b-u*f),e[10]=2*(f*b+u*i),e[11]=u*u-i*i-f*f+b*b,{type:a&16383,region:a>>>16&31,flags:a>>>21}}class Wt{h;tbl;mask;start;pts;inst;constructor(t,n,e,a){this.h=a;let o=1;for(;o<e*2;)o<<=1;o=Math.max(o,1024),this.tbl=o,this.mask=o-1;const d=new Uint32Array(e),h=new Uint32Array(o+1);for(let r=0;r<e;r++)d[r]=this.k(Math.floor(t[3*r]/a),Math.floor(t[3*r+1]/a),Math.floor(t[3*r+2]/a)),h[d[r]+1]++;for(let r=0;r<o;r++)h[r+1]+=h[r];this.start=h.slice(),this.pts=new Float32Array(e*3),this.inst=new Int32Array(e);const l=h.slice(0,o);for(let r=0;r<e;r++){const m=l[d[r]]++;this.pts[3*m]=t[3*r],this.pts[3*m+1]=t[3*r+1],this.pts[3*m+2]=t[3*r+2],this.inst[m]=n[r]}}k(t,n,e){return(Math.imul(t,73856093)^Math.imul(n,19349663)^Math.imul(e,83492791))&this.mask}near(t,n,e,a,o=-1){const d=this.h,h=Math.floor(t/d),l=Math.floor(n/d),r=Math.floor(e/d),m=a*a;for(let p=-1;p<=1;p++)for(let I=-1;I<=1;I++)for(let c=-1;c<=1;c++){const u=this.k(h+p,l+I,r+c);for(let i=this.start[u],f=this.start[u+1];i<f;i++){const b=this.pts[3*i]-t,A=this.pts[3*i+1]-n,v=this.pts[3*i+2]-e;if(b*b+A*A+v*v<m&&this.inst[i]!==o)return!0}}return!1}}function Rt(s){let t=0;const n=new Float64Array(12);for(const d of s)t+=gt[d.w[4*d.i]&16383];const e=new Float32Array(t*3),a=new Int32Array(t);let o=0;for(const d of s){const{type:h}=Dt(d.w,d.i,d.lo,n),l=Ut[h],r=gt[h];for(let m=0;m<r;m++){const p=R[3*(l+m)],I=R[3*(l+m)+1],c=R[3*(l+m)+2];e[3*o]=n[0]+n[3]*p+n[4]*I+n[5]*c,e[3*o+1]=n[1]+n[6]*p+n[7]*I+n[8]*c,e[3*o+2]=n[2]+n[9]*p+n[10]*I+n[11]*c,a[o]=d.id,o++}}return{xyz:e,owner:a,n:t}}function Ct(s){let t=s>>>0;return()=>{t=t+1831565813>>>0;let n=t;return n=Math.imul(n^n>>>15,n|1),n^=n+Math.imul(n^n>>>7,n|61),((n^n>>>14)>>>0)/4294967296}}const Ot=(s,t,n)=>{let e=Math.imul(s^2654435769,2246822507)^Math.imul(t+1663821227,3266489909)^Math.imul(n+668265263,374761393);return e^=e>>>16,e=Math.imul(e,2146121005),e^=e>>>15,e=Math.imul(e,2221713035),e^=e>>>16,e>>>0};function D(s){return Math.sqrt(-2*Math.log(s()+1e-12))*Math.cos(6.2831853*s())}function dn(s,t,n){const e=Math.abs(s)+Math.abs(t)+Math.abs(n)||1;let a=s/e,o=t/e;if(n<0){const l=a;a=(1-Math.abs(o))*(l>=0?1:-1),o=(1-Math.abs(l))*(o>=0?1:-1)}const d=Math.max(0,Math.min(255,Math.round((a+1)*127.5))),h=Math.max(0,Math.min(255,Math.round((o+1)*127.5)));return d|h<<8}const wn=(()=>{const s=[20.3,14.9,14.3,8.3,1.8,12,28.4],t=s.reduce((e,a)=>e+a);let n=0;return s.map(e=>n+=e/t)})(),mn=.65,un=.95;function ct(s,t,n,e){return s.n[0]*t+s.n[1]*n+s.n[2]*e-s.d}function Mn(s,t,n,e){const a=[],o=[],d={lipCand:0,lipMasked:0,glucan:0,mannan:0},h=new Float64Array(12),l=(c,u,i,f)=>{if(n.mode==="shell")return Math.hypot(c,u,i)>lt[0]-70-f;const b=ct(n,c,u,i);return b>=n.lo-f&&b<=n.hi+f};if(e.lipids&&t.tris.length&&n.mode==="slab"){const c=[],u=t.inst.length/4,i=t.ghost.length/4;for(let g=0;g<u;g++)t.inst[4*g]>>>21&1&&c.push({w:t.inst,i:g,lo:t.lo,id:g});for(let g=0;g<i;g++)c.push({w:t.ghost,i:g,lo:t.lo,id:u+g});const f=Rt(c),b=new Wt(f.xyz,f.owner,f.n,1),A=t.verts,v=t.normals,O=t.tris,L=O.length/3;for(let g=0;g<L;g++){const P=O[3*g],_=O[3*g+1],w=O[3*g+2],y=A[3*P],M=A[3*P+1],x=A[3*P+2],S=A[3*_],k=A[3*_+1],U=A[3*_+2],E=A[3*w],Z=A[3*w+1],q=A[3*w+2],H=ct(n,y,M,x),K=ct(n,S,k,U),J=ct(n,E,Z,q);if(Math.max(H,K,J)<n.lo-3||Math.min(H,K,J)>n.hi+3)continue;const T=S-y,z=k-M,B=U-x,j=E-y,N=Z-M,F=q-x,G=.5*Math.hypot(z*F-B*N,B*j-T*F,T*N-z*j),Q=t.fmem[g],Y=Q===8?[1]:[-1,1];for(const $ of Y){const C=Ct(Ot(s,g,$+7)),V=G/mn,rt=Math.floor(V)+(C()<V-Math.floor(V)?1:0);if(!rt)continue;const tt=Math.ceil(Math.sqrt(rt)),nt=tt*tt,it=new Uint16Array(nt);for(let X=0;X<nt;X++)it[X]=X;for(let X=0;X<rt;X++){const ft=X+Math.floor(C()*(nt-X)),W=it[X];it[X]=it[ft],it[ft]=W}for(let X=0;X<rt;X++){const ft=it[X],W=Math.floor(Math.sqrt(ft)),Xt=ft-W*W,at=Xt>>1,tn=(Xt&1)===0;let ht,yt,It,_t,Et,qt;tn?(ht=W,yt=at,It=W+1,_t=at,Et=W+1,qt=at+1):(ht=W,yt=at,It=W+1,_t=at+1,Et=W,qt=at+1);let dt=C(),wt=C();dt+wt>1&&(dt=1-dt,wt=1-wt);const nn=ht+dt*(It-ht)+wt*(Et-ht),Ht=yt+dt*(_t-yt)+wt*(qt-yt),st=(nn-Ht)/tt,ot=Ht/tt,Bt=1-st-ot,en=y+st*T+ot*j,sn=M+st*z+ot*N,on=x+st*B+ot*F;let mt=Bt*v[3*P]+st*v[3*_]+ot*v[3*w],ut=Bt*v[3*P+1]+st*v[3*_+1]+ot*v[3*w+1],Mt=Bt*v[3*P+2]+st*v[3*_+2]+ot*v[3*w+2];const Ft=Math.hypot(mt,ut,Mt)||1;mt/=Ft,ut/=Ft,Mt/=Ft;const rn=C();let pt=0;for(;pt<6&&rn>wn[pt];)pt++;d.lipCand++;let At=en,xt=sn,vt=on;Q===8&&(At-=mt*2.4,xt-=ut*2.4,vt-=Mt*2.4);const Pt=ct(n,At,xt,vt);if(Pt<n.lo-2.5||Pt>n.hi+2.5)continue;const an=pt===5?[1.55,.75]:[1.95,1.25,.55];let Kt=!1;for(const jt of an)if(b.near(At+mt*$*jt,xt+ut*$*jt,vt+Mt*$*jt,un)){Kt=!0;break}if(Kt){d.lipMasked++;continue}const cn=(dn(mt,ut,Mt)|pt<<16|($>0?1:0)<<19|(Q&15)<<20)>>>0;a.push(At,xt,vt,cn)}}}}if(e.wall&&lt.length>0){const c=lt[0],u=lt[1],i=lt[2],f=t.lo,b=t.inst.length/4,A=[];for(let w=0;w<b;w++)(t.inst[4*w]>>>16&31)===9&&A.push({w:t.inst,i:w,lo:f,id:w});const v=Rt(A),O=new Wt(v.xyz,v.owner,v.n,1.4),L=(w,y,M)=>w>=f[0]&&w<f[0]+200&&y>=f[1]&&y<f[1]+200&&M>=f[2]&&M<f[2]+200,g=self.__vol?.[s]??0,P=Math.round(g*Vt.glucanPerNm3);if(P>0){const w=Ct(Ot(s,1234567,1)),M=Math.ceil(P/48);for(let x=0;x<M;x++){const S=24+Math.floor(w()*49);let k=0,U=0,E=0,Z=!1;for(let F=0;F<40&&!Z;F++){k=f[0]+w()*200,U=f[1]+w()*200,E=f[2]+w()*200;const G=Math.hypot(k,U,E);Z=G>i+2&&G<u-1}if(!Z)continue;let q=Math.hypot(k,U,E),H=k/q,K=U/q,J=E/q,T=D(w),z=D(w),B=D(w),j=T*H+z*K+B*J;T-=j*H,z-=j*K,B-=j*J;let N=Math.hypot(T,z,B)||1;T/=N,z/=N,B/=N;for(let F=0;F<S&&L(k,U,E);F++){l(k,U,E,1)&&!O.near(k,U,E,1.35)&&(o.push(k,U,E,(2|(x&65535)<<4)>>>0),d.glucan++),T+=.12*D(w),z+=.12*D(w),B+=.12*D(w),q=Math.hypot(k,U,E),H=k/q,K=U/q,J=E/q,j=T*H+z*K+B*J;const G=q,Q=G<i+3&&j<0||G>u-2&&j>0?-1:.6;T-=(1-Q)*j*H,z-=(1-Q)*j*K,B-=(1-Q)*j*J,N=Math.hypot(T,z,B)||1,T/=N,z/=N,B/=N,k+=1.2*T,U+=1.2*z,E+=1.2*B}}}let _=0;for(const w of A){const{type:y}=Dt(w.w,w.i,f,h),M=h[0],x=h[1],S=h[2],k=Math.hypot(M,x,S);if(k<u-5)continue;if(n.mode==="slab"){const T=ct(n,M,x,S);if(T<n.lo-45-zt[y]||T>n.hi+45+zt[y]){_+=80;continue}}else if(k<c-120)continue;const U=Ct(Ot(s,w.i,99)),E=Ut[y],Z=gt[y];let q=Vt.mannanPerProtein;const H=M/k,K=x/k,J=S/k;for(;q>0;){const T=Math.min(Math.ceil(q),12+Math.floor(U()*27));q-=T;const z=E+Math.floor(U()*Z),B=R[3*z],j=R[3*z+1],N=R[3*z+2];let F=M+h[3]*B+h[4]*j+h[5]*N,G=x+h[6]*B+h[7]*j+h[8]*N,Q=S+h[9]*B+h[10]*j+h[11]*N,Y=F-M,$=G-x,C=Q-S,V=Math.hypot(Y,$,C)||1;Y=Y/V+H,$=$/V+K,C=C/V+J,V=Math.hypot(Y,$,C)||1,Y/=V,$/=V,C/=V,F+=1.3*Y,G+=1.3*$,Q+=1.3*C;const rt=_++;for(let tt=0;tt<T;tt++){const nt=Math.hypot(F,G,Q);if(nt>c-.6)break;l(F,G,Q,1)&&!O.near(F,G,Q,1.15,tt<2?w.i:-1)&&(o.push(F,G,Q,(3|(rt&65535)<<4)>>>0),d.mannan++),Y+=.35*D(U)+.18*F/nt,$+=.35*D(U)+.18*G/nt,C+=.35*D(U)+.18*Q/nt,V=Math.hypot(Y,$,C)||1,Y/=V,$/=V,C/=V,F+=1.05*Y,G+=1.05*$,Q+=1.05*C}}}}const r=new Float32Array(a.length),m=new Uint32Array(r.buffer);for(let c=0;c<a.length;c+=4)r[c]=a[c],r[c+1]=a[c+1],r[c+2]=a[c+2],m[c+3]=a[c+3];const p=new Float32Array(o.length),I=new Uint32Array(p.buffer);for(let c=0;c<o.length;c+=4)p[c]=o[c],p[c+1]=o[c+1],p[c+2]=o[c+2],I[c+3]=o[c+3];return{lipids:r,free:p,stats:d}}self.onmessage=async s=>{const t=s.data;try{if(t.cmd==="init"){Gt=yn(t),self.__vol=t.volGlucan,await Gt,self.postMessage({cmd:"ready"});return}if(await Gt,t.cmd==="load"){const n=await Tt(`tiles/${t.ix}_${t.iy}_${t.iz}.bin.gz`),e=Zt(n),a=e.verts.length/3,o=new Float32Array(a*8),d=new Uint32Array(o.buffer);for(let r=0;r<a;r++)o[8*r]=e.verts[3*r],o[8*r+1]=e.verts[3*r+1],o[8*r+2]=e.verts[3*r+2],d[8*r+3]=e.vmem[r],o[8*r+4]=e.normals[3*r],o[8*r+5]=e.normals[3*r+1],o[8*r+6]=e.normals[3*r+2];for(et.set(t.idx,{lo:e.lo,inst:e.inst,ghost:e.ghost,verts:e.verts,normals:e.normals,tris:e.tris,fmem:e.fmem,vmem:e.vmem});et.size>hn;)et.delete(et.keys().next().value);const h=e.inst.slice(),l=e.tris.slice();self.postMessage({cmd:"loaded",idx:t.idx,inst:h,free:e.free,vb:o,ib:l,lo:e.lo},[h.buffer,e.free.buffer,o.buffer,l.buffer])}else if(t.cmd==="detail"){let n=et.get(t.idx);if(n)et.delete(t.idx),et.set(t.idx,n);else{const o=Zt(await Tt(`tiles/${t.ix}_${t.iy}_${t.iz}.bin.gz`));n={lo:o.lo,inst:o.inst,ghost:o.ghost,verts:o.verts,normals:o.normals,tris:o.tris,fmem:o.fmem,vmem:o.vmem},et.set(t.idx,n)}const e=performance.now(),a=Mn(t.idx,n,t.range,t.want);self.postMessage({cmd:"detail",idx:t.idx,key:t.key,lipids:a.lipids,free:a.free,stats:a.stats,ms:performance.now()-e},[a.lipids.buffer,a.free.buffer])}}catch(n){self.postMessage({cmd:"error",idx:t.idx,key:t.key,req:t.cmd,error:String(n&&n.stack||n)})}};\n',bt=typeof self<"u"&&self.Blob&&new Blob(["URL.revokeObjectURL(import.meta.url);",St],{type:"text/javascript;charset=utf-8"});function vn(i){let e;try{if(e=bt&&(self.URL||self.webkitURL).createObjectURL(bt),!e)throw"";const t=new Worker(e,{type:"module",name:i?.name});return t.addEventListener("error",()=>{(self.URL||self.webkitURL).revokeObjectURL(e)}),t}catch{return new Worker("data:text/javascript;charset=utf-8,"+encodeURIComponent(St),{type:"module",name:i?.name})}}const et=200,tt=173.3,yt=70;class nt{used;n;count=0;high=0;constructor(e){this.n=e,this.used=new Uint8Array(e)}alloc(e){if(e<=0)return[];if(this.n-this.count<e)return null;const t=[];for(let n=0;n<this.n&&t.length<e;n++)this.used[n]||(this.used[n]=1,t.push(n));return this.count+=e,this.high=Math.max(this.high,t[t.length-1]+1),t}free(e){for(const t of e)this.used[t]&&(this.used[t]=0,this.count--);for(;this.high>0&&!this.used[this.high-1];)this.high--}get freeCount(){return this.n-this.count}}class xn{idx;r;lib;G;workers=[];state;loaded=new Map;instA;lipA;freeA;slots=[];pageTile;pageLocal;tileState=new Uint32Array(484*22);desired=[];desiredSet=new Set;inflight=0;maxInflight=12;detInflight=0;maxDetInflight=3;frameNo=0;onChange=()=>{};onTileLoaded=[];bytesLoaded=0;tilesLoadedTotal=0;starved=!1;detailStats={lipids:0,glucan:0,mannan:0};detailEnabled=!0;lastLoadErr="";constructor(e,t,n,s,o,r){this.idx=e,this.r=t,this.lib=n,this.G=s.G,this.state=new Uint8Array(e.n),this.instA=new nt(t.instPages),this.lipA=new nt(t.lipPages),this.freeA=new nt(t.freePages);for(let l=t.maxSlots-1;l>=0;l--)this.slots.push(l);this.pageTile=new Int32Array(t.instPages).fill(-1),this.pageLocal=new Uint32Array(t.instPages);const c=Array.from(e.volGlucan),u=new URL("data/",document.baseURI).href;for(let l=0;l<o;l++){const d=new vn;d.onmessage=p=>this.onMessage(p.data),d.postMessage({cmd:"init",base:u,G:s.G,volGlucan:c,embedded:r,wall:{glucanPerNm3:e.wall.glucan_beads_per_nm3,mannanPerProtein:e.wall.mannan_beads_per_protein}}),this.workers.push(d)}}tileCenter(e){return[this.idx.ix[e]*et+100,this.idx.iy[e]*et+100,this.idx.iz[e]*et+100]}update(e,t,n,s,o,r,c){if(this.frameNo++,r)return;const u=e.eye,l=e.frustum(t),d=this.idx.n,p=this.G[0],m=[s[0],s[1],s[2]],y=100*(Math.abs(m[0])+Math.abs(m[1])+Math.abs(m[2])),g=[];for(let h=0;h<d;h++){const b=this.tileCenter(h),A=ce(b,u);let C=!0;for(let T=0;T<5;T++)if(l[4*T]*A[0]+l[4*T+1]*A[1]+l[4*T+2]*A[2]+l[4*T+3]<-243.3){C=!1;break}if(C){if(n){const T=$(m,b)-s[3];if(T+y+yt<0||T-y-yt>o)continue}else{const T=ye(b);if(T<p-150-tt||T>p+tt)continue;const M=ce(u,b);if($(b,M)<-T*260)continue}g.push({i:h,d:ye(A)})}}g.sort((h,b)=>h.d-b.d);let P=this.r.instPages*X*.94,x=this.r.maxSlots-8;const a=[];for(const h of g){const b=this.idx.nInst[h.i];if(b>P||x<=0){this.starved=!0;break}P-=Math.ceil(b/X)*X,x--,a.push(h.i)}a.length===g.length&&(this.starved=!1),this.desired=a,this.desiredSet=new Set(a);for(const h of a){const b=this.loaded.get(h);b&&(b.lastUsed=this.frameNo)}for(const h of a){if(this.inflight>=this.maxInflight)break;this.state[h]===0&&(this.state[h]=1,this.inflight++,this.workers[h%this.workers.length].postMessage({cmd:"load",idx:h,ix:this.idx.ix[h],iy:this.idx.iy[h],iz:this.idx.iz[h]}))}this.updateDetail(e,n,s,o,c)}updateDetail(e,t,n,s,o){if(!this.detailEnabled)return;const r=e.eye,c=e.pxPerNmAtUnit(o),u=[n[0],n[1],n[2]],l=Math.min(s,250),d=[];for(const m of this.desired){const y=this.loaded.get(m);if(!y)continue;const g=this.idx.nTri[m]>0,P=this.idx.volGlucan[m]>0||(this.idx.labelMask[m]&2)!==0;if(!g&&!P)continue;const x=this.tileCenter(m),a=Math.max(1,ye(ce(x,r))-tt);c/a*.55<.8||!t&&!P||d.push({t:y,d:a})}d.sort((m,y)=>m.d-y.d);const p=new Set(d.slice(0,40).map(m=>m.t.idx));for(const m of this.loaded.values())m.det&&!p.has(m.idx)&&!m.detPending&&this.dropDetail(m);for(const{t:m}of d.slice(0,40)){if(this.detInflight>=this.maxDetInflight)break;if(m.detPending)continue;const y=m.det;let g=!1;if(y&&(t?g=!y.shell&&$(y.n,u)>.99985&&y.a<=n[3]+.5&&y.b>=n[3]+l-.5:g=y.shell),g)continue;const P=Math.max(10,.35*l),x=t?{n:u,d:n[3],lo:-P,hi:l+P,mode:"slab"}:{n:u,d:n[3],lo:0,hi:0,mode:"shell"};m.detPending=!0,this.detInflight++;const a=t?`${n[3]-P}|${n[3]+l+P}`:"shell";m._req={n:u,a:n[3]-P,b:n[3]+l+P,shell:!t};const h=this.idx.nTri[m.idx]>0,b=this.idx.volGlucan[m.idx]>0||(this.idx.labelMask[m.idx]&2)!==0;this.workers[m.idx%this.workers.length].postMessage({cmd:"detail",idx:m.idx,key:a,ix:this.idx.ix[m.idx],iy:this.idx.iy[m.idx],iz:this.idx.iz[m.idx],range:x,want:{lipids:h&&t,wall:b}})}}dropDetail(e){this.lipA.free(e.lipPages),this.freeA.free(e.detPages),this.clearSPages(this.r.lipPagesBuf,e.lipPages),this.clearSPages(this.r.freePagesBuf,e.detPages),e.lipPages=[],e.detPages=[],e.nLip=0,e.nDet=0,e.det=null,e.underMesh&&(e.underMesh=!1,this.r.setMeshSlot(e.slot,!1)),this.syncHigh()}clearSPages(e,t){const n=new Uint32Array(4);for(const s of t)this.r.device.queue.writeBuffer(e,s*16,n)}syncHigh(){this.r.usedInstPages=this.instA.high,this.r.usedLipPages=this.lipA.high,this.r.usedFreePages=this.freeA.high}evictOne(e){let t=null;for(const n of this.loaded.values())n.idx===e||this.desiredSet.has(n.idx)||(!t||n.lastUsed<t.lastUsed)&&(t=n);return t?(this.evict(t),!0):!1}evict(e){const t=this.r.device.queue,n=new Float32Array(4);for(const s of e.pages)t.writeBuffer(this.r.pageOrg,s*16,n),this.pageTile[s]=-1;this.instA.free(e.pages),this.freeA.free(e.basePages),this.clearSPages(this.r.freePagesBuf,e.basePages),this.lipA.free(e.lipPages),this.clearSPages(this.r.lipPagesBuf,e.lipPages),this.freeA.free(e.detPages),this.clearSPages(this.r.freePagesBuf,e.detPages),e.vb?.destroy(),e.ib?.destroy(),this.slots.push(e.slot),this.loaded.delete(e.idx),this.state[e.idx]=0,this.tileState[lt(this.idx.ix[e.idx],this.idx.iy[e.idx],this.idx.iz[e.idx])]=0,t.writeBuffer(this.r.tileState,0,this.tileState),this.syncHigh()}allocWithEvict(e,t,n){let s=e.alloc(t);for(;!s&&this.evictOne(n);)s=e.alloc(t);return s}onMessage(e){if(e.cmd!=="ready"){if(e.cmd==="error"){if(console.warn("worker error",e.req,e.idx,e.error),this.lastLoadErr=String(e.error).split(`
`)[0],e.req==="load"&&(this.state[e.idx]=0,this.inflight--),e.req==="detail"){this.detInflight--;const t=this.loaded.get(e.idx);t&&(t.detPending=!1)}return}e.cmd==="loaded"&&(this.inflight--,this.onLoaded(e)),e.cmd==="detail"&&(this.detInflight--,this.onDetail(e))}}onLoaded(e){const t=e.idx;if(!this.desiredSet.has(t)){this.state[t]=0,this.onChange();return}const n=e.inst,s=n.length/4,o=this.r.device.queue;this.slots.length||this.evictOne(t);const r=this.allocWithEvict(this.instA,Math.ceil(s/X),t);if(!r||!this.slots.length){this.state[t]=0,this.starved=!0,r&&this.instA.free(r);return}const c=this.slots.pop(),u=e.lo,l=new ArrayBuffer(16),d=new Float32Array(l),p=new Uint32Array(l);d[0]=u[0]-80,d[1]=u[1]-80,d[2]=u[2]-80,r.forEach((A,C)=>{const T=Math.min(X,s-C*X);o.writeBuffer(this.r.inst,A*X*16,n.buffer,C*X*16,T*16),p[3]=T,o.writeBuffer(this.r.pageOrg,A*16,l),this.pageTile[A]=t,this.pageLocal[A]=C*X});const m=e.free,y=m.length/4;let g=[],P=!1;if(y){const A=this.allocWithEvict(this.freeA,Math.ceil(y/Q),t);A?(g=A,this.writeSpheres(this.r.freePool,this.r.freePagesBuf,A,m,c)):P=!0}let x=null,a=null;const h=e.ib.length;h&&(x=this.r.device.createBuffer({size:e.vb.byteLength,usage:GPUBufferUsage.VERTEX|GPUBufferUsage.COPY_DST}),o.writeBuffer(x,0,e.vb),a=this.r.device.createBuffer({size:Math.ceil(e.ib.byteLength/4)*4,usage:GPUBufferUsage.INDEX|GPUBufferUsage.COPY_DST}),o.writeBuffer(a,0,e.ib),this.r.setMeshSlot(c,!1));const b={idx:t,slot:c,lo:u,inst:n,nInst:s,pages:r,basePages:g,nBase:y,lipPages:[],detPages:[],nLip:0,nDet:0,vb:x,ib:a,nIdx:h,lastUsed:this.frameNo,det:null,detPending:!1,underMesh:!1,baseFreeMissing:P};this.loaded.set(t,b),this.state[t]=2,this.tileState[lt(this.idx.ix[t],this.idx.iy[t],this.idx.iz[t])]=1,o.writeBuffer(this.r.tileState,0,this.tileState),this.bytesLoaded+=this.idx.gz[t],this.tilesLoadedTotal++,this.syncHigh();for(const A of this.onTileLoaded.slice())A(b);this.onChange()}writeSpheres(e,t,n,s,o){const r=this.r.device.queue,c=s.length/4;n.forEach((u,l)=>{const d=Math.min(Q,c-l*Q);r.writeBuffer(e,u*Q*16,s.buffer,s.byteOffset+l*Q*16,d*16),r.writeBuffer(t,u*16,new Uint32Array([d,o,1,0]))})}onDetail(e){const t=this.loaded.get(e.idx);if(!t)return;t.detPending=!1;const n=t._req;this.lipA.free(t.lipPages),this.clearSPages(this.r.lipPagesBuf,t.lipPages),this.freeA.free(t.detPages),this.clearSPages(this.r.freePagesBuf,t.detPages),t.lipPages=[],t.detPages=[];const s=e.lipids,o=e.free,r=s.length/4,c=o.length/4;if(r){const l=this.lipA.alloc(Math.ceil(r/Q))??this.freeDetailSpace(this.lipA,Math.ceil(r/Q),t.idx);l&&(t.lipPages=l,this.writeSpheres(this.r.lipPool,this.r.lipPagesBuf,l,s,t.slot))}if(c){const l=this.freeA.alloc(Math.ceil(c/Q))??this.freeDetailSpace(this.freeA,Math.ceil(c/Q),t.idx);l&&(t.detPages=l,this.writeSpheres(this.r.freePool,this.r.freePagesBuf,l,o,t.slot))}t.nLip=t.lipPages.length?r:0,t.nDet=t.detPages.length?c:0,t.det={n:n.n,a:n.a,b:n.b,shell:n.shell};const u=t.nLip>0;u!==t.underMesh&&t.nIdx&&(t.underMesh=u,this.r.setMeshSlot(t.slot,u)),this.detailStats.lipids+=e.stats.lipCand,this.syncHigh(),this.onChange()}freeDetailSpace(e,t,n){const s=[...this.loaded.values()].filter(o=>o.idx!==n&&(o.lipPages.length||o.detPages.length)&&!this.desiredSet.has(o.idx));for(const o of s){this.dropDetail(o);const r=e.alloc(t);if(r)return r}return null}meshDraws(){const e=[];for(const t of this.loaded.values())t.nIdx&&t.vb&&t.ib&&e.push({slot:t.slot,vb:t.vb,ib:t.ib,nIdx:t.nIdx,under:t.underMesh});return e}instanceAt(e){const t=Math.floor(e/X),n=this.pageTile[t];if(n<0)return null;const s=this.loaded.get(n);if(!s)return null;const o=this.pageLocal[t]+e%X;return{t:s,local:o,...this.decode(s,o)}}decode(e,t){const n=e.inst,s=n[4*t],o=n[4*t+1],r=n[4*t+2],c=360/65535,u=[e.lo[0]-80+(o&65535)*c,e.lo[1]-80+(o>>>16)*c,e.lo[2]-80+(r&65535)*c];return{type:s&16383,region:s>>>16&31,flags:s>>>21,p:u}}poolIndex(e,t){return e.pages[Math.floor(t/X)]*X+t%X}tileBySlot(e){for(const t of this.loaded.values())if(t.slot===e)return t;return null}countLoaded(e){let t=0;for(const n of this.loaded.values()){const s=n.inst;for(let o=0;o<n.nInst;o++)(s[4*o]&16383)===e&&t++}return t}nearestInTile(e,t,n){let s=-1,o=1e18;for(let r=0;r<e.nInst;r++){if((e.inst[4*r]&16383)!==t)continue;const c=this.decode(e,r),u=(c.p[0]-n[0])**2+(c.p[1]-n[1])**2+(c.p[2]-n[2])**2;u<o&&(o=u,s=r)}return s}totals(){let e=0,t=0,n=0,s=0;for(const o of this.loaded.values())e+=o.nInst,t+=o.nLip,n+=o.nBase+o.nDet,s+=o.nIdx/3;return{tiles:this.loaded.size,inst:e,lip:t,free:n,mesh:s,inflight:this.inflight+this.detInflight,desired:this.desired.length}}pendingForView(){let e=0;for(const n of this.desired)this.state[n]!==2&&e++;let t=0;for(const n of this.loaded.values())n.detPending&&t++;return e+t}}class Pn{G;PT;PQ;PC;pores;c;grid;n3;constructor(e){this.G=e.G,this.PT=e.PT,this.PQ=e.PQ,this.PC=e.PC,this.pores=e.pores,this.c=e.consts,this.grid=e.grid,this.n3=e.grid.n**3}cell(e,t,n){const s=this.grid.n,o=this.grid.ext,r=2*o/s,c=u=>Math.min(s-1,Math.max(0,Math.floor((u+o)/r)));return(c(e)*s+c(t))*s+c(n)}sdfPrim(e,t,n,s){const o=this.PT[e],r=this.PQ,c=8*e;if(o===0)return Math.hypot(t-r[c],n-r[c+1],s-r[c+2])-r[c+3];if(o===1){const y=t-r[c],g=n-r[c+1],P=s-r[c+2],x=r[c+3]-r[c],a=r[c+4]-r[c+1],h=r[c+5]-r[c+2],b=Math.min(1,Math.max(0,(y*x+g*a+P*h)/(x*x+a*a+h*h)));return Math.hypot(y-x*b,g-a*b,P-h*b)-r[c+6]}if(o===2){const y=t-r[c],g=n-r[c+1],P=s-r[c+2],x=r[c+3],a=r[c+4],h=r[c+5],b=y*x+g*a+P*h,A=Math.hypot(y-b*x,g-b*a,P-b*h),C=Math.max(A-r[c+6],0);return Math.sqrt(C*C+b*b)-r[c+7]}const u=Math.hypot(t,n,s)||1e-9,l=Math.min(1,Math.max(-1,(t*r[c]+n*r[c+1]+s*r[c+2])/u)),d=Math.acos(l),p=Math.max(d-r[c+4],0)*u,m=u-r[c+5];return Math.sqrt(p*p+m*m)-r[c+6]}label(e,t,n){const s=this.G,o=this.c.PERI_W,r=this.c.PORE_R,c=Math.hypot(e,t,n);if(c>s[0])return 0;if(c>s[1])return 1;if(c>s[2])return 2;const u=this.cell(e,t,n),l=this.grid,d=[1e9,1e9,1e9,1e9,1e9,1e9,1e9,1e9,1e9];let p=1e9,m=1e9;for(let T=l.prim_off[u];T<l.prim_off[u+1];T++){const M=l.prim_idx[T],k=this.sdfPrim(M,e,t,n),N=this.PC[M];N===7?d[1]=Math.min(d[1],k):N===0?d[3]=Math.min(d[3],k):N===4?d[4]=Math.min(d[4],k):N===5?p=Math.min(p,k):N===6?m=Math.min(m,k):N===2?d[6]=Math.min(d[6],k):N===3?d[7]=Math.min(d[7],k):N===1&&(d[8]=Math.min(d[8],k))}const y=Math.max(p,-m);if(d[1]<0)return 4;const g=e-s[3],P=t-s[4],x=n-s[5],a=Math.hypot(g,P,x),h=a-(s[6]-.5*o),b=r+.5*o;let A=!1;for(let T=l.pore_off[u];T<l.pore_off[u+1];T++){const M=l.pore_idx[T],k=this.pores[3*M],N=this.pores[3*M+1],O=this.pores[3*M+2],q=g*k+P*N+x*O;if(q<=0)continue;if(Math.hypot(g-q*k,P-q*N,x-q*O)<b){A=!0;break}}if(Math.abs(h)-.5*o<0&&!A)return 5;if(a<s[6]-.5*o){if(A&&h>-.5*o-8)return 14;const T=e-s[7],M=t-s[8],k=n-s[9];return a<s[6]-o&&T*T+M*M+k*k<s[10]*s[10]?7:6}return a<s[6]?14:d[4]<0?y<0?10:9:d[3]<0?8:d[8]<0?11:d[6]<0?12:d[7]<0?13:3}}const Mn=`
<h2>About this model</h2>
<p>This is one complete haploid <i>Saccharomyces cerevisiae</i> cell (strain BY4741, log phase in YPD), built as a
structural whole-cell model in the spirit of Maritan <i>et al.</i> (J. Mol. Biol. 2022; whole <i>Mycoplasma genitalium</i> cell built with cellPACK)
and drawn in the style of David Goodsell. Every macromolecule is drawn from its experimental (PDB) or predicted
(AlphaFold DB) coordinates, at one bead per amino-acid residue or nucleotide when it is large enough on screen.
Further away, molecules switch to a sphere-tree stand-in and then to one sphere each. This only happens once
the detail is smaller than a pixel. Click any molecule to inspect it and to load its full-atom model in Mol*.</p>
<p>This page is a summary. The full description of the model (data sources, assumptions, design decisions, validation,
densities and limitations) is under <a href="#" id="aboutToModel" style="color:#8fc3ff">Model details</a>.</p>

<h3>What is in the cell</h3>
<ul>
<li><b>21.6 million macromolecule instances (most multi-chain assemblies), 5,638 molecular species.</b> Protein abundances are from the PaxDb integrated dataset rescaled to
4.2 × 10⁷ protein molecules per cell (Ho <i>et al.</i>), in a 42 µm³ cell. There are 169,000 80S ribosomes, 3.3 million tRNAs, 73,150 nucleosomes and 117 nuclear pore complexes.</li>
<li><b>Structures.</b> 55% of protein mass is drawn from PDB entries. These include the ribosome 4V88, nuclear pore 7N9F,
nucleosome 1ID3, 26S proteasome 4CR2, FAS 2UV8, Pol I/II/III, V-ATPase 6O7V, Pma1 7VH5, CCT, about 100 further Complex Portal
complexes with a yeast PDB assembly (e.g. phosphofructokinase, eIF2/eIF3, the mitochondrial ribosome, mediator, exosome, snoRNPs with
their RNA) and 133 homo-oligomers. Membrane complexes (ATP synthase, respiratory complexes II-IV, Sec translocon, OST, EMC, glucan
synthase, flippases, ...) are oriented in their membrane automatically from UniProt transmembrane and topology annotations.
14% are AlphaFold DB predicted complexes: 452 homodimers and 505 heterodimers (ipTM &ge; 0.5, partners sharing a
compartment; e.g. actin-cofilin, actin-profilin, RNR2-RNR4, SOD1-CCS1, eIF4G-eIF4A); 30% are AlphaFold monomers.
AlphaFold residues with pLDDT &lt; 50 are trimmed. Outer-membrane beta-barrels (TOM40, SAM50, MDM10) are oriented by their strands.</li>
<li><b>Nuclear pore FG mesh:</b> the disordered FG-repeat domains of NSP1, NUP49, NUP57 and NUP159 (128 per pore, 6.1 million
residues in all) are threaded as self-avoiding residue chains from where the 7N9F model ends, filling the pore channel and
reaching into the cytoplasm (checkbox "FG repeats").</li>
<li><b>mRNA and polysomes:</b> 30,000 full-length mRNAs (lognormal lengths, mean 1.5 kb including UTRs and poly(A);
40 million nucleotides, one bead per nucleotide) are threaded through the cytosol before the proteins are packed.
23,000 of them carry 93,000 of the 80S ribosomes (about 4 per mRNA): the mRNA runs through each ribosome's channel
(path from PDB 7ZW0), short linkers join neighbouring ribosomes, and the rest of the transcript leaves the first and
last ribosome as 3&prime; and 5&prime; tails. The other 7,000 mRNAs carry no ribosomes. Threads that run out of room end
early (86% of the target length is reached overall). Linker DNA joins consecutive nucleosomes.</li>
<li><b>ER-bound ribosomes:</b> 4,300 ribosomes sit on translocons in the ER and outer nuclear membrane: on the Ssh1 complex
(yeast ribosome-Ssh1 cryo-EM structure 2WW9) or on Sec61 (7KAH, placed through the shared Sss1 subunit). Their number is
limited by the free translocon pools (80% assumed ribosome-bound); the PM-facing side of the cortical ER is too narrow for ribosomes.</li>
<li><b>Microtubules:</b> 21 microtubules (11,500 &alpha;/&beta;-tubulin dimers) grow from &gamma;-tubulin rings on the spindle pole body,
opposite the nucleolus: 16 short kinetochore microtubules whose plus ends anchor the 16 chromosomes (Rabl configuration),
2 more nuclear ones and 3 cytoplasmic ones about 1 &micro;m long that curve around organelles. The 13-protofilament lattice
and the ring come from the native yeast &gamma;-TuRC-capped microtubule 8QV3; the lattice dimer is 5W3F (8.3 nm repeat).</li>
<li><b>Mitochondrial cristae:</b> all 11,800 ATP synthase copies form 5,900 dimers in rows along the rims of the lamellar
cristae, every 13 nm. Each dimer is two copies of the yeast monomer 8F39 placed in the dimer arrangement of 4B2Q
(cryo-ET), whose 84&deg; angle matches the curvature of the rim (crista lumen 18 nm wide). Complex IV is entirely in
1,560 respiratory supercomplexes III<sub>2</sub>IV<sub>2</sub> (6HU9) in the inner membrane, next to 1,400 free complex III dimers.</li>
<li><b>Actin cables:</b> 8 bundles of 8 filaments, 12 nm apart (33,000 actin subunits, 92 &micro;m), run 90 nm beneath the
plasma membrane, built from the yeast F-actin cryo-EM structure 9GO5 (167&deg; and 2.76 nm per subunit). Both strands of
every filament are coated head-to-tail by 6,600 tropomyosin TPM1 dimers (200 residues, one per 5 actins), following the
tropomyosin path of the mammalian actin-tropomyosin structure 5JLF wound as an ideal coiled coil.</li>
<li><b>Actin patches:</b> 22 endocytic patches under the plasma membrane, each a branched network of about 1,500 actin
subunits: short filaments branch at 73&deg; from Arp2/3 complexes (2,200 branches; geometry from the branch-junction
structure 8UXW) and are capped by capping protein (8F8Q), barbed ends toward the membrane. The yeast Arp2/3 complex
(all 7 subunits) and capping protein are AlphaFold models superposed on those templates. A further 4,800 Arp2/3
complexes are free in the cytosol. Assumption: 30% of actin in cables, 30% in patches, the rest monomeric
(much of it bound to profilin, cofilin or Srv2).</li>
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
<li>Bead clash rate 1.4×10⁻⁶.</li>
<li>Containment violations 4.3×10⁻⁵.</li>
<li>Macromolecule densities (protein plus the RNA of ribosomes and tRNAs and the DNA of nucleosomes; lipids and small
molecules not counted): cytosol 89 mg/mL (67 protein + 22 RNA), nucleoplasm 62 mg/mL, mitochondrial matrix 112 mg/mL.</li>
<li>99.8% of requested copies were placed outside mitochondria.</li>
</ul>
<p>Positions are stored to 0.0055 nm and orientations to 0.004°.</p>

<h3>Honest caveats</h3>
<ul>
<li><b>Layout.</b> The organelle shapes and positions are a plausible illustration, not a tomogram of one specific cell.</li>
<li><b>Abundances.</b> Absolute abundances are uncertain to about 2-fold.</li>
<li><b>Mitochondria.</b> The mitochondrion saturated: only 65% of matrix and 86% of IMS proteins fit, so their copies are under-represented.</li>
<li><b>Chromatin</b> is a coarse Rabl-like arrangement of nucleosome chains in 59 fragments, each chromosome starting at its kinetochore microtubule. rDNA repeats are not modelled.
Linker DNA and mRNA paths are schematic.</li>
<li><b>Procedural content.</b> Lipids, β-glucan and mannan are generated procedurally in your browser from deterministic seeds. They are not stored
coordinates, and polysaccharides are only partly clash-checked.</li>
<li><b>Spindle pole body.</b> Only the &gamma;-tubulin rings are modelled, not the SPB plaques; PaxDb lists fewer SPC98 copies than 21 rings need.</li>
<li><b>Not included:</b> septins, fimbrin cross-links, ribosome-bound chaperones (NAC, RAC/Ssb), metabolites and water.
The FG nucleoporins absent from 7N9F (NUP100, NUP116, NUP145N, NUP42 and the basket NUP1, NUP60, NUP2) have no FG domains.</li>
<li><b>Dropped copies.</b> 9,946 transmembrane-protein copies assigned to soluble compartments were dropped, and 6,888 copies had no structure.</li>
</ul>
<p style="color:#a9a69e">Renderer: WebGPU compute culling and level-of-detail selection, sphere impostors with exact cut caps, analytic
compartment backgrounds and Goodsell-style outlines. Full-atom models are fetched live from RCSB PDB and AlphaFold DB
when you request them.</p>
<p style="color:#a9a69e">Data sources: RCSB Protein Data Bank (CC0) and the individual entries named above; AlphaFold Protein Structure
Database monomers and predicted complexes (CC-BY 4.0, Jumper <i>et al.</i> 2021; Varadi <i>et al.</i> 2024); PaxDb protein abundances;
Complex Portal; UniProt and Gene Ontology annotations; IntAct interactions; Alliance of Genome Resources paralogs;
lipid composition from Ejsing <i>et al.</i> 2009; total protein from Ho <i>et al.</i> 2018.</p>
`,wt=i=>i.replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;");function $e(i){const e=[];return i=i.replace(/`([^`]+)`/g,(t,n)=>`\0${e.push(wt(n))-1}\0`),i=wt(i).replace(/\[([^\]]+)\]\(([^)\s]+)\)/g,(t,n,s)=>`<a href="${s}" target="_blank" rel="noopener">${n}</a>`).replace(/\*\*([^*]+)\*\*/g,"<b>$1</b>").replace(/\*([^*\s][^*]*)\*/g,"<i>$1</i>"),i.replace(/\u0000(\d+)\u0000/g,(t,n)=>`<code>${e[+n]}</code>`)}const it=/^([-*]|\d+\.)\s+/,st=i=>i.length-i.trimStart().length,vt=i=>i.trim().replace(/^\||\|$/g,"").split("|").map(e=>e.trim());function It(i){let e="",t=0;for(;t<i.length;){const n=i[t],s=n.trim();if(!s){t++;continue}let o;if(o=s.match(/^(#{1,6})\s+(.*)$/)){e+=`<h${o[1].length}>${$e(o[2])}</h${o[1].length}>`,t++;continue}if(/^-{3,}$/.test(s)){e+="<hr>",t++;continue}if(s.startsWith("|")&&t+1<i.length&&/^\|[\s|:-]+\|$/.test(i[t+1].trim())){const c=vt(s);t+=2;let u="";for(;t<i.length&&i[t].trim().startsWith("|");)u+=`<tr>${vt(i[t]).map(l=>`<td>${$e(l)}</td>`).join("")}</tr>`,t++;e+=`<table><thead><tr>${c.map(l=>`<th>${$e(l)}</th>`).join("")}</tr></thead><tbody>${u}</tbody></table>`;continue}if(o=s.match(it)){const c=/\d/.test(o[1]),u=st(n);let l="";for(;t<i.length;){const d=i[t],p=d.trim().match(it);if(!p||st(d)!==u||/\d/.test(p[1])!==c)break;const m=u+p[0].length,y=[d.slice(m)];for(t++;t<i.length&&(!i[t].trim()||st(i[t])>=m);)y.push(i[t].slice(m)),t++;for(;y.length&&!y[y.length-1].trim();)y.pop();const g=It(y),P=g.match(/^<p>((?:(?!<p>)[\s\S])*)<\/p>$/);l+=`<li>${P?P[1]:g}</li>`}e+=c?`<ol>${l}</ol>`:`<ul>${l}</ul>`;continue}const r=[];for(;t<i.length&&i[t].trim()&&!/^(#{1,6}\s|\||-{3,}$)/.test(i[t].trim())&&!it.test(i[t].trim());)r.push(i[t++].trim());e+=`<p>${$e(r.join(" "))}</p>`}return e}function An(i){return It(i.replace(/\r\n/g,`
`).split(`
`))}const kn=`# Whole yeast cell: model description

This document describes how the structural model shown in the viewer was built. It covers what the model contains, which data and structures it comes from, the assumptions and design decisions behind it, how it was validated, and where it is known to be wrong or incomplete. It is written for readers who want to judge the model scientifically or reuse its numbers.

All numbers refer to the current build (5,638 molecular species, 21.6 million placed instances) unless stated otherwise.

## 1. What the model is

- **The cell.** One haploid *Saccharomyces cerevisiae* cell (strain BY4741 / S288C background) growing exponentially in rich glucose medium (YPD). The cell is unbudded, and the microtubules and chromatin are arranged as in G1.
- **What a molecule is.** Every macromolecule is an individual copy of an experimental (PDB) or predicted (AlphaFold DB) structure. It is drawn at one bead per amino-acid residue (Cα) or per nucleotide (P or C4′). There are 14.2 billion beads in all.
- **How the cell was assembled.** Copy numbers come from proteomics. Compartment volumes and membrane areas come from imaging literature. The organelles are built as analytic geometry, and molecules are packed into them with residue-level collision checks.
- **What it is not:**
  - **Not a dynamic or physical simulation.** There are no forces, no diffusion and no binding equilibria. Except where stated below, molecules are placed at random positions and orientations by random sequential addition. A snapshot from a real cell would show interaction-driven clustering that the model lacks.
  - **Not a reconstruction of a particular cell.** The organelle volumes and areas are taken from measurements, but their shapes and positions are a plausible arrangement, not a tomogram.
- **Lineage.** The approach follows the cellPACK-style whole-cell models of *Mycoplasma genitalium* (Maritan *et al.* 2022). The visual style follows David Goodsell.

## 2. Composition

### 2.1 Proteins

- **Abundances:** PaxDb integrated whole-organism dataset for *S. cerevisiae* (taxon 4932), in ppm, matched to UniProt by ORF.
- **Absolute scale:** ppm values are rescaled to 4.2 × 10⁷ protein molecules per cell, the unified total of Ho, Baryshnikova and Brown (2018), who used the same 42 µm³ cell volume. This gives 2.7 pg of protein per cell (abundance-weighted mean mass 39 kDa).
- **Uncertainty:**
  - PaxDb integrated values are not specific to BY4741 in YPD. Per-protein uncertainty is roughly 2-fold.
  - The absolute total is the largest single uncertainty in the model. 4.2 × 10⁷ molecules in 42 µm³ is 1.0 × 10⁶ per µm³, at the low end of published normalised estimates (1–3 × 10⁶ per µm³; Dolgalev *et al.* 2023). If the true total is higher, every protein density in this document scales up with it.

### 2.2 Nucleic acids and fixed counts

| Component | Copies | Basis |
|---|---|---|
| 80S ribosomes | 169,000 | Ribosomal protein abundances (median over the r-protein families) |
| tRNA | 3.3 × 10⁶ | Assumption, literature order of magnitude |
| mRNA | 30,000 | Assumption; mean length 1,500 nt including UTRs and poly(A) |
| Nucleosomes | 73,150 | 12.07 Mb nuclear genome (rDNA repeats not expanded) ÷ 165 bp repeat length |
| Nuclear pore complexes | 117 | Winey *et al.* 1997 |

The ribosome count is probably low. A growth-law estimate for YPD is about 3 × 10⁵ (Kostinski & Reuveni 2021). Because the count follows the ribosomal-protein abundances, it inherits the uncertainty in the proteome total.

### 2.3 Not included

- **Small molecules and ions:** metabolites, free amino acids, nucleotides, ions and water.
- **Storage carbohydrates:** glycogen and trehalose.
- **Polyphosphate.**
- **rDNA and pre-rRNA:** the rDNA array and nascent pre-rRNA in the nucleolus, and most pre-ribosomal particles.
- **Proteins without a structure:** 6,888 copies.
- **Transmembrane proteins annotated only to soluble compartments:** 9,946 copies (see 4.3).

Lipids and the cell-wall polysaccharides are drawn in the viewer but generated procedurally (section 7).

## 3. Geometry

### 3.1 Volumes

The cell is a sphere of 42 µm³ (mean log-phase volume; Jorgensen *et al.* 2002) with a 110 nm cell wall, the lower bound of 110–200 nm (Orlean 2012).

Organelle volume fractions:
- **Nucleus, vacuole and lipid droplets:** soft X-ray tomography (Uchida *et al.* 2011): nucleus 7%, vacuole 7%, lipid droplets 1%. The nucleolus is 20% of the nucleus.
- **Mitochondria:** raised from the measured 1% to 2%. The assumption is that X-ray tomography under-segments thin tubules. At 1%, the mitochondrial proteome would imply about 480 mg/mL of protein.
- **Golgi, endosomes and peroxisomes:** 1% together (assumption).
- **Cortical ER:** covers 48% of the plasma membrane (Wolf *et al.* 2012).

Measured on the final voxelised geometry (4 nm grid):

| Compartment | Volume (µm³) | Share of cell |
|---|---|---|
| Cytosol | 27.81 | 66.2% |
| Cell wall (glucan 4.09, mannan 2.01) | 6.10 | 14.5% |
| Vacuole | 2.94 | 7.0% |
| Nucleoplasm | 2.02 | 4.8% |
| Nucleolus | 0.59 | 1.4% |
| Mitochondrial matrix | 0.61 | 1.5% |
| Mitochondrial intermembrane space (incl. crista lumen) | 0.24 | 0.6% |
| ER lumen: cortical ER | 0.50 | 1.2% |
| ER lumen: perinuclear space | 0.29 | 0.7% |
| Lipid droplets | 0.42 | 1.0% |
| Golgi / endosomes | 0.34 | 0.8% |
| Peroxisomes | 0.08 | 0.2% |
| **Total** | **41.98** | |

### 3.2 Membranes

Membrane mid-surface areas, measured on the final geometry:

| Membrane | Area (µm²) |
|---|---|
| Plasma membrane | 52.6 |
| ER (cortical sheets and perinuclear) | 51.8 |
| Nuclear envelope | 18.3 |
| Mitochondrial inner membrane | 28.4 |
| Mitochondrial outer membrane | 12.3 |
| Vacuole | 9.9 |
| Golgi / endosomes | 8.9 |
| Peroxisomes | 1.5 |

### 3.3 Organelle shapes

All membranes are analytic signed-distance fields (spheres, capsules, pills and cortical-ER patches), so containment and membrane distance are exact everywhere. Shapes and fixed dimensions:

| Element | Shape and dimensions |
|---|---|
| Nucleus | A sphere, with the nucleolus as a crescent on the side away from the vacuole. Nuclear envelope: perinuclear space 34 nm between the membrane mid-surfaces; 117 pores of 45 nm radius at the mid-plane. |
| Cortical ER | Sheets 27 nm beneath the plasma membrane, with a 20 nm lumen. |
| Mitochondria | One reticulum of tubules with 150 nm radius. Intermembrane space 8 nm. Lamellar cristae with an 18 nm lumen (rim-tube radius 9 nm, chosen to fit ATP synthase dimers; section 5.6). |
| Other organelles | Vacuole, lipid droplets, Golgi/endosome cisternae and peroxisomes are simple primitives. |

**Arrangement.** Organelle positions follow a plausible layout:
- The nucleus is off-centre.
- The spindle pole body faces away from the nucleolus.
- The mitochondria form branched tubules along the cortex.

This layout is illustrative; it is not a reconstruction of a particular cell.

**Known geometric deviations from EM.** Serial-section EM (Yamaguchi *et al.* 2011) gives a nucleus of 10.5% of cell volume; this model has about 7%. Mitochondria are 1.7% there and about 2% here.

## 4. Localization

### 4.1 Source annotations

**UniProt subcellular-location text** is mapped to compartments by ordered regular-expression rules. Examples:
- "nucleolus"
- "mitochondrion inner membrane"
- "endoplasmic reticulum"
- "cell cortex" → plasma membrane

Annotations are filtered by their evidence:
- **GFP data preferred:** where a location token is backed by the genome-wide GFP study of Huh *et al.* (2003), only GFP-backed tokens are kept. This gives the bulk steady-state localization.
- **Proteomics-only calls dropped:** tokens supported only by two high-throughput mitochondrial proteomics papers (Sickmann 2003, Reinders 2006) are discarded. Those studies assign many abundant cytosolic enzymes to mitochondria.

**Fallbacks:**
- **GO cellular component:** used when UniProt has no location. Complex terms are skipped, except pores, ribosomes and proteasomes. If GO lists the cytosol, only cytosolic and nuclear terms are kept.
- **Cytosol:** the default for proteins with no annotation at all (625 proteins).

### 4.2 Weights across compartments

- **Rank order:** a protein's compartments get weights 1, ½, ¼, … in annotation order, with the first-listed location taken as primary.
- **Cytosol and nucleoplasm:** their shares are pooled and re-split by volume, i.e. equal concentration across the nuclear pore.
- **Secondary mitochondrial locations:** a mitochondrial location listed after a cytosolic or nuclear primary one (e.g. SOD1, ADK1 in the intermembrane space) gets 0.1× weight.
- **Mitochondrial sub-compartments:** taken from UniProt or GO where named (inner membrane, outer membrane, intermembrane space, matrix). Otherwise transmembrane proteins are split 75/25 between the inner and outer membranes, and soluble proteins 90/10 between matrix and intermembrane space.

### 4.3 Membrane proteins

- **Definition:** proteins with UniProt TRANSMEM segments. β-barrels (TOM40, SAM50, MDM10, MDM34) have no TRANSMEM annotation and are added explicitly.
- **Assigned to soluble compartments only:** these take their membrane from GO, reading "cell periphery" as plasma membrane. Copies still without a membrane are dropped (9,946).
- **ER membrane proteins:** 15% go to the nuclear envelope, in proportion to area, because the outer nuclear membrane is continuous with the ER.
- **Lipid-droplet transmembrane proteins:** placed in the ER membrane.

### 4.4 Peripheral proteins

Soluble proteins annotated to a membrane are placed in thin "peripheral" shells next to it:
- **Plasma membrane, nuclear envelope, ER and lipid-droplet surface:** on the cytosolic side.
- **Inner-membrane-associated mitochondrial proteins:** in the matrix.

## 5. Structures

### 5.1 Bead representation and mass

- **Beads:** the Cα for amino acids, and P (or C4′) for nucleotides.
- **Mass:** for densities, a residue counts as 110 Da and a nucleotide as 330 Da.
- **What is missing from the beads:** residues missing from a structure (unresolved in the PDB entry, or trimmed from an AlphaFold model) are absent from the model and from its mass.

### 5.2 Choice of structure (in order of preference)

1. **Curated PDB complexes:**

   | Complex | PDB |
   |---|---|
   | 80S ribosome | 4V88 |
   | Nuclear pore complex | 7N9F |
   | Nucleosome | 1ID3 |
   | 26S proteasome | 4CR2 |
   | Fatty acid synthase | 2UV8 |
   | CCT | 4V94 |
   | RNA polymerase I | 4C2M |
   | RNA polymerase II | 1WCM |
   | RNA polymerase III | 5FJ8 |
   | Pma1 hexamer | 7VH5 |
   | V-ATPase | 6O7V |

2. **Complex Portal complexes** (126 added, 22 of them membrane complexes) with a *S. cerevisiae* RCSB assembly. Rules:
   - **Acceptance:** the assembly must cover at least 80% of the subunit mass, and at least 60% of the residues must be modelled.
   - **RNA chains** are kept only where they belong to the complex: snoRNAs, snRNAs, telomerase RNA, and mitoribosomal rRNA.
   - **Copy number:** the median subunit pool, capped at 2× the limiting subunit. Complexes with four or more subunit types are capped at 2× the 25th-percentile subunit.
   - **Mitochondrially encoded subunits** (COB, COX1–3, ATP6, ATP8, OLI1, VAR1) are ignored when setting copy numbers, because mass spectrometry badly under-detects them.
   - **Subunit bookkeeping:** subunits are consumed from the monomer pool, so total protein copies are conserved.
3. **Homo-oligomers:**
   - **PDB:** about 130 homo-oligomeric assemblies.
   - **AlphaFold:** 452 AlphaFold-predicted homodimers with ipTM ≥ 0.5.
4. **Heterodimers:** 505 AlphaFold-predicted heterodimers, from 1,250 models. Requirements:
   - ipTM ≥ 0.5.
   - The two partners share a compartment. A membrane protein may pair with a soluble partner on its cytosolic or matrix side.
   - Not both partners are subunits of an already-modelled complex.
   - Neither partner forms a larger homo-oligomer.

   **Copy allocation:**
   - Each protein's free copies are split evenly among its accepted partners.
   - Pairs supported by Complex Portal or a direct IntAct interaction bind all available copies; the others bind 50%. This is an assumption: without affinity data, stable and transient interactions cannot be told apart.
5. **AlphaFold monomers** for everything else. Residues with pLDDT < 50 are trimmed, which removes 3.4% of placed protein mass overall. The loss is concentrated in nuclear and nucleolar proteins: 15–16% of mass there, measured on an earlier build.

**Result by mass of placed protein:** PDB 55.3%, AlphaFold complexes 14.3%, AlphaFold monomers 30.4%.

### 5.3 Membrane orientation

- **Method:** UniProt TRANSMEM residues are mapped onto the structure, through RCSB's UniProt alignment for PDB chains. A search finds the membrane normal and offset that put the most TM residues, and the fewest other residues, inside a ±1.6 nm slab.
- **Side:** from UniProt topology domains. Without them, the side carrying more soluble mass is taken as cytoplasmic.
- **β-barrels:** oriented by their β-strands, found from Cα geometry.
- **Check:** the method reproduces hand-oriented Pma1 and V-ATPase (z-coordinate correlation 1.000 and 0.9996).
- **Result:** in the packed cell, the median distance of TM residues from the membrane mid-surface is 0.73 nm, and 97.2% lie within 2 nm.

## 6. Supramolecular organisation built explicitly

These structures are built before random packing. They then act as fixed obstacles for everything packed afterwards.

### 6.1 Nuclear pores and FG mesh

- **Pores:** 117 copies of 7N9F, with the pore axis along the local nuclear-envelope normal.
- **FG mesh:** the disordered FG-repeat domains of the FG-Nups modelled in 7N9F are threaded as self-avoiding residue chains (0.38 nm per residue, low persistence) from the first modelled residue of each chain:

  | Nup | Residues | Copies per pore |
  |---|---|---|
  | NSP1 | 1–636 | 48 |
  | NUP57 | 1–286 | 32 |
  | NUP49 | 1–269 | 32 |
  | NUP159 | 388–1102 | 16 |

  - **Confinement:** the chains stay within 60 nm of their anchor, in the pore channel, cytoplasm or nucleoplasm.
  - **Result:** 6.1 of 6.9 million residues placed (89%).
  - **Missing:** FG-Nups absent from 7N9F (NUP100, NUP116, NUP145N, NUP42, NUP1, NUP60, NUP2), so the nuclear basket side is sparse.
  - **Caveat:** the chain conformations are excluded-volume random walks, not a physical model of FG cohesion.

### 6.2 Microtubules and the spindle pole body

- **Number:** 21 microtubules, made of 11,505 αβ-tubulin dimers:
  - 16 kinetochore microtubules, 229 nm mean length;
  - 2 other nuclear microtubules, 311 nm;
  - 3 cytoplasmic microtubules, 1.0 µm.
- **Lattice:** an ideal 13-protofilament, 3-start lattice with a seam. The lateral geometry comes from 8QV3 (native yeast γ-TuRC-capped microtubule, cryo-ET). The dimer is 5W3F, with an 8.3 nm axial repeat.
- **Caps:** each minus end has a γ-tubulin ring of 7 γ-TuSCs from 8QV3.
- **Growth:** from 90 nm plaques at the spindle pole body. Microtubules bend away from membranes, with a minimum radius of curvature of 300 nm.
- **Assumptions:** the lengths (G1-like), and the 70% cap on the use of the tubulin pool.
- **Not built:** the spindle pole body's own plaques (Spc42, Spc72 layers).

### 6.3 Chromatin

- **Walks:** 16 chromosomes are self-avoiding, persistent nucleosome walks (12.5 nm step), with lengths proportional to chromosome size.
- **Rabl-like layout:** each chromosome starts within 20 nm of its kinetochore microtubule tip.
- **Result:** 73,150 nucleosomes in 59 fragments; walks that get trapped restart elsewhere.
- **Linker DNA:** joins consecutive nucleosomes, with a schematic path.
- **Not modelled:** chromosome territories, telomere anchoring, the rDNA array, nucleosome-depleted regions and transcribing Pol II.

### 6.4 Actin

The ACT1 pool (111,000) is split by assumption: 30% in cables, 30% in patches, and the rest monomeric or in profilin, cofilin and Srv2 heterodimers.

**Cables:**
- **Number:** 8 cables of 8 filaments, 33,300 subunits and 92 µm of filament in total.
- **Position:** 90 nm beneath the plasma membrane, along great circles.
- **Bundle:** filaments sit on a hexagonal lattice 12 nm apart (assumed fimbrin bundle spacing).
- **Filament structure:** yeast F-actin 9GO5, 167° and 2.76 nm per subunit.
- **Tropomyosin:** 6,610 TPM1 coiled-coil dimers, one per 5 actins per strand. The path is transferred from the mammalian actin–tropomyosin structure 5JLF, because no yeast structure exists. The chains are wound as an ideal coiled coil.

**Endocytic patches:**
- **Number:** 22 patches of about 1,500 actin each.
- **Branching:** dendritic networks branch at 73° from Arp2/3 complexes (2,169 branches). The branch-junction geometry is 8UXW. The yeast Arp2/3 subunits are AlphaFold models superposed on the fission-yeast template, 32–74% identity.
- **Capping:** capping protein (from 8F8Q) caps filaments after 20 ± 6 subunits. Filaments that end against an obstacle stay uncapped.
- **Orientation:** barbed ends face the membrane.
- **Collision checks:** residue by residue inside each patch.

Not built: fimbrin cross-links, septins, formins.

### 6.5 Ribosomes and mRNA

- **Polysomes:** 70% of cytosolic ribosomes are in polysomes (assumption), with sizes 1 + Poisson(4).
  - **Chain geometry:** each ribosome's mRNA entry lies 3–12 nm beyond the previous one's exit, within a 40° cone. Neighbouring ribosomes end up 25–31 nm apart, as in cryo-ET polysomes.
  - **mRNA channel:** the path through each ribosome comes from 7ZW0, superposed on 4V88.
  - **Result:** 93,000 ribosomes on 23,000 mRNAs.
- **Full-length mRNA:** all 30,000 mRNAs are threaded at full length through the cytosol before protein packing.
  - **Lengths:** lognormal, mean 1,500 nt, σ = 0.5 in log space, minimum 300 nt.
  - **Representation:** one bead per nucleotide, 0.6 nm apart.
  - **Polysome tails:** the untranslated remainder leaves the first and last ribosome as 3′ and 5′ tails. The 5′ side gets 10–40% of it.
  - **Result:** 86% of the target nucleotides were placed. Threads that run out of room end early.
- **ER-bound ribosomes:** 4,264 ribosomes are docked on translocons in the ER and the outer nuclear membrane.
  - **Templates:** 3,349 on Ssh1 (2WW9, yeast 80S–Ssh1) and 915 on Sec61 (7KAH, placed through the shared Sss1 subunit).
  - **Assumption:** 80% of free translocons are ribosome-bound.
  - **Not modelled:** ER-bound polysomes and their mRNA.
- **Not placed on ribosomes:** the ribosome-associated factors NAC, RAC/Ssb and eEF1A are free in the cytosol.

### 6.6 Mitochondrial inner membrane

- **ATP synthase dimers:** all 11,762 ATP synthase copies form 5,881 dimers. Each dimer is two copies of the yeast monomer 8F39 arranged as in the cryo-ET dimer 4B2Q (dimer angle 83.5°).
  - **Placement:** in rows along the crista rims, every 13 nm (assumed spacing).
  - **Missing subunits:** 8F39 lacks the dimer-specific subunits e and k.
- **Respiratory supercomplexes:** complex IV is entirely in 1,561 III₂IV₂ supercomplexes (6HU9). The remaining 1,422 complex III dimers are free.
- **Not modelled:** III₂IV₁ supercomplexes, and confinement of supercomplexes to the flat crista membranes.

## 7. Procedural content (generated in the browser)

These are generated deterministically from fixed seeds when you zoom in. They are not stored coordinates and do not count in the densities.

- **Lipids:**
  - **Packing:** 0.65 nm² per lipid per leaflet.
  - **Composition:** whole-cell, from Ejsing *et al.* (2009): PI 20.3%, PE 14.9%, PC 14.3%, PA 8.3%, PS 1.8%, ergosterol 12%, other 28.4%. The same composition is used for every membrane, with no leaflet asymmetry and no organelle-specific lipids (e.g. cardiolipin).
  - **Masking:** lipids within 0.95 nm of a protein bead are removed.
- **Cell wall:**
  - **Layers:** an inner β-1,3-glucan layer (semi-flexible chains, 0.018 beads/nm³) and an outer mannan layer (branched brushes on the mannoproteins).
  - **Amount:** bead counts are set so the wall is about 20% of cell dry mass.
  - **Collisions:** only partly checked.
- **Motion (Display → Motion; off by default):** illustrative animation, not molecular dynamics. Nothing is simulated and neighbouring molecules are not checked for collisions while they move.
  - **Conformational motion from AlphaFold 3 ensembles:**
    - **Coverage:** 3,692 of the 4,403 single-chain AlphaFold species (54.5% of all placed molecules, 29.7% of the mass).
    - **Ensembles:** each of these proteins was predicted by AlphaFold 3 in many random multi-protein "pools", about 600 models per protein (in-house predictions, unpublished).
    - **Fit:** for each protein, 64 pool models are sampled, keeping those with mean pLDDT ≥ 60 over the residues the viewer draws (at least 16). Their Cα atoms are superposed on the library structure, and a principal component analysis gives 4 modes, which capture a median 77% of the variance.
    - **Playback:** the sampled conformers are ordered into a closed nearest-neighbour path in mode space. Each molecule follows that path at its own phase, about 2 s per conformer, interpolated smoothly between conformers.
    - **Size:** no amplitude cap. Typical bead excursion is 0.74 nm per protein (median), and floppy termini can swing several nm.
    - **Skipped proteins:** 707 have too few confident models, and 4 have too few pool models. These get only the generic motion below.
  - **Limits of the conformational motion:** the modes are linear, so large swings of termini or domains stretch the chain. Along the paths, a median 3.7% of neighbouring Cα–Cα distances leave 0.30–0.46 nm at "On", and 10.6% at "Exaggerated" (2×). Proteins also move independently of their neighbours and can briefly interpenetrate. Motions within complexes, dimers and PDB structures are not represented.
  - **Generic motion for everything else:** a small rigid-body wobble whose amplitude falls with molecular size (about 0.35 nm translation for small proteins). Membrane proteins move only within the bilayer and rotate only about the membrane normal. Residues with low pLDDT flex by up to 0.45 nm. Free polymers (mRNA, linker DNA, wall polysaccharides, FG repeats) wiggle along their length.
  - **Static:** actin filaments, microtubules, nuclear pores and nucleosomes stay still, so that their assemblies stay intact.

## 8. Packing

### 8.1 Spatial decomposition

- **Tiles:** the cell is cut into 200 nm cubic tiles; 6,416 of them intersect the cell.
- **Copies per tile:** each species' expected copy number per compartment is split over tiles in proportion to the compartment's volume (or membrane area) in each tile. A multinomial draw with stochastic rounding conserves the cell totals.
- **Concentrations:** uniform within each compartment.

### 8.2 Algorithm

Random sequential addition, run independently in each tile.

**Placement:**
- **Order:** membrane proteins first, then the rest by decreasing size.
- **Soluble molecules:** a uniformly random position in the target compartment and a uniformly random orientation.
- **Membrane proteins:** projected onto the exact membrane surface, with their membrane normal along the local surface normal and a random spin about it.
- **Retries:** a bounded number of trials per copy. Copies that don't fit after that are left out.

**Collision checks:**
- **Sphere trees:** each structure is approximated by k-means clusters of about 30 residues, each sphere padded by a 0.42 nm bead radius. Spheres may interpenetrate by up to 0.3 nm.
- **Compartments:** a copy must lie entirely within its compartment and clear of membranes, which are treated as 2 nm half-thickness slabs.

**Parallelism:**
- **Phases:** tiles are packed in 8 checkerboard phases, so neighbouring tiles are never packed simultaneously.
- **Halos:** each tile sees the already placed molecules of its neighbours through a halo.
- **Obstacles:** the pre-placed structures of section 6 are fixed obstacles.

### 8.3 Consequences

- **The density is not an equilibrium value.** Random sequential addition jams at lower density than real cytoplasm, and the outcome depends on insertion order. Where the requested copies exceed what random sequential addition can fit, the shortfall is reported, not hidden.
- **The shortfall is concentrated in mitochondria.** Placed copies are 97.9% of requested overall and 99.8% outside mitochondria. In the mitochondrial matrix only 65% fit, so matrix protein density is under-represented.
- **No interaction-driven organisation.** There are no metabolons and no condensates; enzymes are not co-localised beyond what the explicit structures in section 6 provide.

## 9. Validation

| Check | Result |
|---|---|
| Residue-bead clashes (beads of different molecules closer than 0.3 nm) | 1.4 × 10⁻⁶ of 2.9 M sampled beads in 64 test cubes; 0 between tiles |
| Containment (beads outside their compartment) | 4.3 × 10⁻⁵ of 17 M beads |
| Membrane-protein orientation vs local surface normal | median 0.0°, 100% within 1° |
| TM residues within 2 nm of the membrane mid-surface | 97.2% (median 0.73 nm) |
| ATP synthase dimer TM residues within 2 nm of the mid-surface | 98.8% |
| Placed / requested copies | 97.9% overall; 99.8% outside mitochondria; matrix 65%; nucleolus 95% |
| Inventory protein mass vs proteome | within 0.4% before packing |

## 10. Densities, mass and volume fractions

### 10.1 Macromolecular densities

Densities are the mass of placed macromolecules per compartment volume. They count protein plus the RNA of ribosomes, tRNAs and complexes, the DNA of nucleosomes, and mRNA. Lipids, polysaccharides and small molecules are not counted.

| Compartment | Protein (mg/mL) | RNA/DNA (mg/mL) | Total (mg/mL) | Dry volume fraction |
|---|---|---|---|---|
| Cytosol | 67 | 22 | 89 | 6.1% |
| Nucleoplasm | 56 | 6 | 62 | 4.4% |
| Nucleolus | 87 | 6 | 94 | 6.7% |
| Mitochondrial matrix | 108 | 4 | 113 | 8.1% |
| Mitochondrial intermembrane space | 27 | 0 | 27 | 1.9% |
| ER lumen | 32 | 0 | 32 | 2.3% |
| Golgi / endosomes | 64 | 0 | 64 | 4.7% |
| Peroxisomes | 83 | 0 | 83 | 6.1% |
| Vacuole | 5 | 0 | 5 | 0.4% |

- **Dry volume fraction:** computed with partial specific volumes of 0.73 mL/g for protein and 0.55 mL/g for nucleic acid.
- **Hydrated volume fraction:** with a hydration shell of about 0.3–0.4 g water per g, the hydrated volume fraction is roughly 1.5× the dry one: about 9–10% in the cytosol.
- **Excluded volume:** the volume excluded to a diffusing particle is larger still, and grows with the particle's size. It can be computed from the coordinates as accessible volume versus probe radius; this has not been done yet.

### 10.2 Whole-cell totals

| Quantity | Value |
|---|---|
| Protein in the inventory | 2.7 pg |
| Protein placed | 2.2 pg |
| RNA and DNA placed | 0.64 pg |
| mRNA | 0.02 pg |

Placed protein is lower than the inventory for four reasons:
- residues not resolved in PDB entries;
- AlphaFold residues trimmed for pLDDT < 50 (3.4% of mass);
- mitochondrial copies that did not fit;
- dropped and structure-less copies.

With full-length sequences, cytosolic protein would be about 80–85 mg/mL.

### 10.3 Comparison with dry-mass measurements

Measured cytoplasmic dry-mass densities are around 200 mg/mL, more than twice the macromolecular density here. These are rough estimates of the missing components, averaged over the cytoplasm (cytosol plus organelles, excluding the vacuole):

| Component | Estimate (mg/mL) | Basis |
|---|---|---|
| Macromolecules as placed | about 90 | This model |
| Full-length protein correction | about +15 | Section 10.2 |
| Lipids of internal membranes and lipid droplets | about 20–25 | About 130 µm² of internal membrane plus 0.4 µm³ of lipid droplets in this model, about 0.8 pg of lipid |
| Small molecules and ions | probably 30–40 | Free amino acids, nucleotides, glutathione, sugar phosphates, phosphate and K⁺, a few hundred mM in total |
| Glycogen and trehalose | low in log phase on glucose | Can exceed 10% of dry weight after the diauxic shift |

Together these give about 150–160 mg/mL. The remainder is within the uncertainty of the absolute proteome size (section 2.1). A higher protein total than 4.2 × 10⁷ molecules would close most of the gap.

## 11. Known limitations

| Area | Limitation |
|---|---|
| Physics | Static random sequential addition: no interactions, no dynamics, random orientations, and densities below the jamming limit of real cytoplasm. The optional motion is illustrative (section 7). |
| Abundances | About 2-fold per-protein uncertainty. The proteome total is probably low. Ribosomes are probably low (169,000 vs about 3 × 10⁵ estimated). Fixed counts for tRNA and mRNA are assumptions. |
| Geometry | Illustrative layout. Nucleus smaller than in EM. A single mitochondrial reticulum. No bud, so no cell-cycle-dependent structures. |
| Mitochondria | Matrix only 65% filled; mitochondrially encoded subunits are inferred from their complexes. |
| Disorder | AlphaFold residues with pLDDT < 50 are deleted rather than represented, except the FG domains of 4 Nups. |
| Localization | Bulk steady state from annotations. 625 proteins default to the cytosol. Dual localizations are weighted by rank, not measured. |
| Nucleolus | No rDNA, Pol I transcription units or pre-rRNA; few pre-ribosomes. |
| Cytoskeleton | No septins, fimbrin, formins or spindle-pole-body plaques. Actin partition and filament lengths are assumptions. |
| Lipids and wall | One lipid composition for all membranes. Procedural cell wall. |
| Missing classes | Metabolites, ions, water, glycogen, trehalose, polyphosphate. |

## 12. Data sources and references

### 12.1 Data

| Source | Use | Terms |
|---|---|---|
| RCSB Protein Data Bank | Experimental structures | CC0 |
| AlphaFold Protein Structure Database | Monomers and predicted homo- and heterodimers (Jumper *et al.* 2021, Nature 596:583; Varadi *et al.* 2024, Nucleic Acids Res. 52:D368) | CC-BY 4.0 |
| PaxDb | Integrated whole-organism protein abundances for *S. cerevisiae* | |
| UniProt | Proteome UP000002311: sequences, subcellular location, transmembrane segments, topology, subunit structure | |
| Gene Ontology | Cellular-component annotations | |
| Complex Portal (EMBL-EBI) | Curated complexes for taxon 559292 | |
| IntAct | Physical interactions, to decide heterodimer binding | |
| Alliance of Genome Resources | Paralog pairs | |
| AlphaFold 3 pool predictions (in-house, unpublished) | Conformational ensembles behind the optional motion | not public |

### 12.2 Literature behind the parameters

| Reference | Used for |
|---|---|
| Ho B, Baryshnikova A, Brown GW (2018). *Cell Systems* 6:192–205 | Total protein molecules per cell |
| Jorgensen P *et al.* (2002). *Science* 297:395–400 | Cell volume |
| Uchida M *et al.* (2011). *Yeast* 28:227–236 | Organelle volume fractions (soft X-ray tomography) |
| Yamaguchi M *et al.* (2011). *Microscopy* 60:321–335, doi:10.1093/jmicro/dfr052 | Serial-section EM comparison |
| Orlean P (2012). *Genetics* 192:775–818 | Cell-wall thickness and composition |
| Wolf W *et al.* (2012). *PLoS ONE* 7:e39703 | Cortical ER coverage |
| Winey M *et al.* (1997). *Mol. Biol. Cell* 8:2119–2132 | NPC number |
| Ejsing CS *et al.* (2009). *PNAS* 106:2136–2141 | Lipid composition |
| Huh WK *et al.* (2003). *Nature* 425:686–691 | GFP localization |
| Sickmann A *et al.* (2003), *PNAS*; Reinders J *et al.* (2006), *J. Proteome Res.* | Excluded high-throughput organelle calls |
| Dolgalev G *et al.* (2023). *Int. J. Mol. Sci.*, doi:10.3390/ijms24032081 | Protein density comparison |
| Kostinski S, Reuveni S (2021). *Phys. Rev. Research*, doi:10.1103/physrevresearch.3.013020 | Ribosome count comparison |
| Maritan M *et al.* (2022). *J. Mol. Biol.* 434:167351 | Whole-cell modelling approach |

### 12.3 PDB entries used for assemblies and templates

| Assembly or template | PDB entries |
|---|---|
| Ribosome and polysomes | 4V88, 7ZW0 |
| ER translocons | 2WW9, 7KAH |
| Nuclear pore and nucleosome | 7N9F, 1ID3 |
| Microtubules | 8QV3, 5W3F |
| Actin and actin-binding proteins | 9GO5, 5JLF, 8UXW, 8F8Q |
| Mitochondrial inner membrane | 8F39, 4B2Q, 6HU9 |
| Other curated complexes | 4CR2, 2UV8, 4V94, 4C2M, 1WCM, 5FJ8, 7VH5, 6O7V |

The other complexes and homo-oligomers are listed per species in the viewer's inspect panel: click any molecule.

## 13. Reproducibility and data access

- **Pipeline.** The model is built by a Python pipeline (numba-accelerated) in five stages: inventory, geometry, packing, validation and encoding. A full rebuild takes roughly 10–15 minutes on a recent multi-core desktop or laptop, with fixed seeds throughout. An optional sixth stage fits the motion modes (section 7); it needs the pool predictions, which are not public. The shipped modes are in \`data/library/modes.bin\`.
- **Source code.** The pipeline and the viewer source are in \`source/\`; \`source/README.md\` explains how to rebuild.
- **This repository** holds the built viewer and the encoded data.
  - The tiles under \`data/tiles/\` store every instance as a species index, a position quantised to about 0.0055 nm and an orientation quantised to about 0.004°.
  - The residue-bead library is in \`data/library/\`.
- **Full-atom coordinates.** The viewer fetches full-atom models live from RCSB and AlphaFold DB for display in Mol*.
`,f=i=>document.getElementById(i),ae=i=>i>=1e6?(i/1e6).toFixed(i>=1e7?1:2)+" M":i>=1e4?(i/1e3).toFixed(0)+" k":i.toLocaleString("en-US"),Bt=i=>i>=1e3?(i/1e3).toFixed(i>=1e4?0:2).replace(/\.?0+$/,"")+" µm":(i>=10?i.toFixed(0):i.toFixed(1))+" nm";let J,Oe,xe,Rt,fe,je,R,I;const w=new yn,S=f("gpu"),v={cutOn:!0,plane:[0,0,1,0],slab:60,colorMode:0,regionMask:16777215,freeHide:0,highlight:!1,selId:0,shaded:!0,organelle:!1,motion:0,l0px:1.6,l1px:2.5,minpx:.3,outline:1,showMeshes:!0,showLipids:!0,renderScale:1};let ge=!1,L=!0,ht=0,ot=0,Nt="auto",te;const He={actin:["F_actin","Arp23","CP_cap","TPM1_coil"],microtubules:["MT_tubulin","gTuRC_nuc","gTuRC_cyto"],ribosomes:["ribosome_80S","ribosome_80S_Ssh1","ribosome_80S_Sec61"],tRNA:["tRNA"],nucleosomes:["nucleosome"],pores:["NPC"]},Tn=[...He.actin,...He.microtubules,"NPC","nucleosome"];let be=-1,Dt=null,z=null,oe=null;const ne=[];function qe(i,e){f("loading").querySelector(".m").textContent=i,f("loading").querySelector(".pb div").style.width=`${Math.round(e*100)}%`}function xt(i){const e=f("loading").querySelector(".err");e.style.display="block",e.innerHTML=i,f("loading").querySelector(".t").textContent="Could not start the viewer"}function Me(i,e=2200){const t=f("toast");t.textContent=i,t.style.display="block",clearTimeout(Me._t),Me._t=setTimeout(()=>t.style.display="none",e)}async function Cn(){if(!("gpu"in navigator)){xt(`This browser does not expose <b>WebGPU</b>. Use a current Chrome or Edge (version 113 or later), or Safari 26 or later.
      On Linux Chrome, enable <code>chrome://flags/#enable-unsafe-webgpu</code> and <code>#enable-vulkan</code>.
      Firefox 141+ supports WebGPU on Windows.<br><br>Also check that the page is served over <code>http://localhost</code>
      (for example <code>python -m http.server</code>). Opening the file directly with <code>file://</code> does not work.`);return}try{qe("Loading cell geometry",.01),[xe,fe]=await Promise.all([ze("geometry.json"),ze("bookmarks.json")]),Rt=new Pn(xe),Oe=await fn(),J=await mn(qe);const i=await Pe("types_tiles.bin"),e=J.types.length;je={off:new Uint32Array(i,0,e+1),tiles:new Uint16Array(i,(e+1)*4)},qe("Starting WebGPU",.7),R=new bn(S),R.onLost=r=>{Me(r,1e4)},await R.init(J,xe),te=new Uint32Array(e);for(const r of Object.values(He))for(const c of r){const u=J.byId.get(c);u!==void 0&&(te[u]|=4)}for(const r of Tn){const c=J.byId.get(r);c!==void 0&&(te[c]|=8)}R.device.queue.writeBuffer(R.typeFlags,0,te),qe("Compiling shaders and starting tile streaming",.9);let t;if(window.__YC_EMBED__){t={};for(const[r,c]of Object.entries(window.__YC_EMBED__))(r.startsWith("tiles/")||r.startsWith("library/"))&&(t[r]=Lt(c))}const n=Math.max(2,Math.min(4,(navigator.hardwareConcurrency||4)-1));I=new xn(Oe,R,J,xe,n,t),I.onChange=()=>{L=!0},Bn();const s=new URLSearchParams(location.search).get("view"),o=s===null?0:isNaN(Number(s))?Math.max(0,fe.findIndex(r=>r.name.toLowerCase().startsWith(s.toLowerCase()))):Math.min(fe.length-1,Math.max(0,Number(s)));ke(o,!1),f("loading").style.display="none",requestAnimationFrame(Ht),window.yc={R,T:I,cam:w,fs:v,lib:J,idx:Oe,applyBookmark:ke,setType:Te,pickAt:Qe,bookmarks:fe,setQuality:$t,describe:pt,ready:!0,state:()=>({pending:I.pendingForView(),detailPending:I.detInflight,loaded:I.loaded.size,totals:I.totals(),stats:R.stats,cam:{t:w.target,d:w.dist},plane:v.plane,gpuMs:Vt,gpuBusy:Je,dirty:L}),render:()=>{L=!0}}}catch(i){console.error(i),xt(`${String(i?.message||i)}<br><br>See the browser console for details.`)}}function de(){return[v.plane[0],v.plane[1],v.plane[2]]}function re(i,e){i=K(i),v.plane=[i[0],i[1],i[2],e],ue(),L=!0}function Ae(){const i=w.fwd;re(i,$(i,w.target))}let le=.5;try{const i=Number(localStorage.getItem("yc.nudgeStep"));i>0&&(le=i)}catch{}function Ye(i){const e=de();w.target=D(w.target,E(e,i)),re(e,v.plane[3]+i)}function En(i){const e=S.clientWidth/S.clientHeight,t=K(i.normal);let n=K(ce(i.up,E(t,$(i.up,t))));return{target:[...i.target],dist:w.distForWidth(i.width,e),fwd:t,up:n}}function ke(i,e=!0,t){const n=fe[i];if(!n)return;v.organelle&&i!==0&&_e(!1);const s=En(n),o=K(n.normal),r=$(o,n.target),c=de(),u=v.plane[3],l=v.slab;if(v.cutOn=!0,f("cutOn").checked=!0,ge=!1,f("followView").classList.remove("on"),!e){w.set(s),re(o,r),v.slab=n.slab,ue(),L=!0,t?.();return}w.flyTo(s,2600,()=>{re(o,r),v.slab=n.slab,ue(),t?.()},d=>{let p=K(D(E(c,1-d),E(o,d)));ye(D(E(c,1-d),E(o,d)))<.001&&(p=o),v.plane=[p[0],p[1],p[2],u+(r-u)*d],v.slab=Math.exp(Math.log(l)+(Math.log(n.slab)-Math.log(l))*d),ue()}),Xe(n),L=!0}function Xe(i){const e=f("caption");if(!i){e.style.display="none";return}e.innerHTML=`<b>${i.name}</b>${i.caption}`,e.style.display="block",clearTimeout(Xe._t),z||(Xe._t=setTimeout(()=>e.style.display="none",9e3))}function Ot(){z={i:0,phase:"fly",t:performance.now()},f("tourBtn").textContent="■ Stop",f("tourBtn").classList.add("on"),ke(0,!0,()=>{z&&(z.phase="wait",z.t=performance.now())})}function we(){z=null,f("tourBtn").textContent="▶ Tour",f("tourBtn").classList.remove("on"),Xe(null)}function Ln(i){if(z){if(z.phase==="wait"&&(I.pendingForView()===0||i-z.t>14e3))z.phase="hold",z.t=i;else if(z.phase==="hold"&&i-z.t>6500){z.i=(z.i+1)%fe.length,z.phase="fly";const e=z.i;ke(e,!0,()=>{z&&z.i===e&&(z.phase="wait",z.t=performance.now())})}}}function Te(i){be>=0&&(te[be]&=-2),be=i,i>=0&&(te[i]|=1),v.highlight=i>=0,R.device.queue.writeBuffer(R.typeFlags,0,te),zt(),L=!0}function zt(){const i=f("selType");if(be<0){i.style.display="none";return}const e=J.types[be],t=I.countLoaded(be);i.style.display="block",i.innerHTML=`<div><b>${e.gene}</b> ${e.name}</div>
    <div class="sub" style="margin:3px 0 6px">${ae(e.copies)} copies in the cell · ${ae(t)} in loaded tiles · ${_t(e)}</div>
    <div class="grid2"><button class="small" id="flyNearest">Fly to nearest</button><button class="small" id="clearType">Clear highlight</button></div>`,f("flyNearest").onclick=()=>Ut(be),f("clearType").onclick=()=>Te(-1)}function Fn(i){if(i=i.trim().toLowerCase(),!i)return[];const e=[];return J.types.forEach((t,n)=>{const s=t.gene.toLowerCase(),o=t.id.toLowerCase(),r=(t.uniprot||"").toLowerCase(),c=(t.name||"").toLowerCase(),u=(t.genes||"").toLowerCase();let l=0;s===i||o===i||r===i?l=100:s.startsWith(i)||o.startsWith(i)?l=60:u.split(/[;,\s]+/).includes(i)?l=55:c.includes(i)?l=30+(c.startsWith(i)?10:0):(u.includes(i)||r.startsWith(i)||(t.pdb||"").toLowerCase()===i)&&(l=20),l&&e.push({i:n,s:l+Math.log10(t.copies+1)})}),e.sort((t,n)=>n.s-t.s),e.slice(0,14).map(t=>t.i)}function Ut(i){const e=je.off[i],t=je.off[i+1];if(t<=e){Me("No copies of this molecule were placed.");return}let n=-1,s=1e18;for(let p=e;p<t;p++){const m=je.tiles[p],y=I.tileCenter(m),g=ye(ce(y,w.target));g<s&&(s=g,n=m)}const o=I.tileCenter(n),r=J.types[i],c=Math.max(90,r.br*9),u=p=>{const m=I.nearestInTile(p,i,o);if(m<0)return;const y=I.decode(p,m),g=de(),P=w.state();P.target=y.p,P.dist=w.distForWidth(c,S.clientWidth/S.clientHeight),w.flyTo(P,900,()=>{re(g,$(g,y.p)-r.br*.4)}),v.plane[3]=$(g,y.p)-r.br*.4;const x=I.poolIndex(p,m);Ke({kind:"inst",id:(gn|x)>>>0,info:I.instanceAt(x)})};I.loaded.get(n);const l=de(),d=w.state();d.target=o,d.dist=w.distForWidth(c*2.5,S.clientWidth/S.clientHeight),v.cutOn=!0,f("cutOn").checked=!0,v.organelle&&_e(!1),w.flyTo(d,2200,()=>{re(l,$(l,o)-20);const p=I.loaded.get(n);if(p)u(p);else{const m=y=>{y.idx===n&&(I.onTileLoaded.splice(I.onTileLoaded.indexOf(m),1),u(y))};I.onTileLoaded.push(m)}}),re(l,$(l,o)-20)}async function Qe(i,e){const t=R.width/S.clientWidth;return R.pick(i*t,e*t)}function _t(i){return i.src==="PDB"?"PDB "+i.pdb:i.src==="AF-multimer"?"AlphaFold complex":"AlphaFold"}function pt(i){const e=i>>>30;if(i===0)return{kind:"none",title:"Outside the cell"};if(e===1){const t=I.instanceAt(i&1073741823);return t?{kind:"inst",id:i,info:t,title:J.types[t.type].gene,sub:Ue[t.region]}:{kind:"none",title:""}}if(e===2){const t=i>>>1&15,n=i&1,s=gt[t]||"membrane";return{kind:"lipid",id:i,mem:t,leaf:n,title:"Lipid bilayer",sub:`${s} · ${t===8?"phospholipid monolayer":n?"cytoplasmic leaflet":t===0?"exoplasmic (outer) leaflet":"lumenal leaflet"}`}}if(e===3){const t=i>>>15&7,n=i&32767,s=["mRNA","Linker DNA","β-1,3-glucan fibril","Mannan chain","pre-rRNA","FG repeats"][t]??"polymer",o=["nucleotide-resolution thread (schematic path)","between consecutive nucleosomes, 2 bp per bead","procedural, 36 glucose per bead","procedural α-mannan on a cell-wall mannoprotein","nascent ribosomal RNA","disordered domain of an FG nucleoporin (NSP1, NUP49, NUP57, NUP159), one bead per residue"][t]??"";return{kind:"free",id:i,fk:t,chain:n,title:s,sub:o}}return i>=512&&i<528?{kind:"band",id:i,title:"Membrane (cut face)",sub:gt[i-512]}:i>=256&&i<272?{kind:"label",id:i,title:hn[i-256],sub:"compartment (not loaded at this distance, or beyond the slab)"}:{kind:"none",title:""}}function Ke(i){Dt=i,v.selId=i&&i.kind==="inst"?i.id>>>0:0,Gt(),L=!0}function Gt(){const i=f("inspect"),e=Dt;if(!e||e.kind==="none"){i.style.display="none";return}i.style.display="block";let t='<button class="close" id="inspClose">×</button>';if(e.kind==="inst"){const n=J.types[e.info.type],s=n.src==="PDB",o=n.nb,r=n.nnt;t+=`<h2>${n.gene}</h2><div>${n.name}</div><div style="margin-top:6px"><span class="badge ${s?"pdb":"af"}">${_t(n)}</span>
      ${n.oligomer>1?`<span class="badge" style="background:#666">${n.oligomer}-mer</span>`:""}</div>
      <div class="kv">
      ${n.genes&&n.genes!==n.gene?`<div>Subunits</div><div style="word-break:break-word">${n.genes.split(/;\s*/).slice(0,40).join(", ")}${n.genes.split(";").length>40?" …":""}</div>`:""}
      ${n.uniprot?`<div>UniProt</div><div>${n.uniprot}</div>`:""}
      <div>Location</div><div>${Ue[e.info.region]}</div>
      <div>Copies / cell</div><div>${n.copies.toLocaleString("en-US")}</div>
      <div>Mass</div><div>${(n.mass/1e3).toFixed(1)} kDa</div>
      <div>Residues drawn</div><div>${(o-r).toLocaleString("en-US")} aa${r?` + ${r.toLocaleString("en-US")} nt`:""}</div>
      <div>Radius</div><div>${n.br.toFixed(1)} nm (Rg ${n.rg.toFixed(1)} nm)</div>
      ${!s&&n.plddt_mean!=null?`<div>Mean pLDDT</div><div>${n.plddt_mean.toFixed(1)} (${Math.round((n.plddt_frac_kept??1)*100)}% of residues kept, pLDDT ≥ 50)</div>`:""}
      <div>Position</div><div>${e.info.p.map(c=>c.toFixed(1)).join(", ")} nm</div>
      </div>
      <div class="grid2"><button class="small" id="loadAtoms">Load atomic model</button><button class="small" id="hlType">Highlight all</button>
      <button class="small" id="centerOn">Center view</button><button class="small" id="hideType">${te[e.info.type]&2?"Show":"Hide"} this type</button></div>
      <div style="margin-top:9px">${n.uniprot?`<a href="https://www.uniprot.org/uniprotkb/${n.uniprot}" target="_blank">UniProt</a>`:""}
      ${s?`<a href="https://www.rcsb.org/structure/${n.pdb}" target="_blank">RCSB ${n.pdb}</a>`:""}
      ${(n.model||n.uniprot)&&!s?`<a href="https://alphafold.ebi.ac.uk/entry/${n.model||n.uniprot}" target="_blank">AlphaFold DB</a>`:""}
      ${/^[A-Z0-9-]+$/.test(n.gene)&&!Ft[n.id]&&!n.id.startsWith("CPX")?`<a href="https://www.yeastgenome.org/locus/${n.gene}" target="_blank">SGD</a>`:""}</div>`}else t+=`<h2>${e.title}</h2><div class="sub" style="font-size:12.5px">${e.sub||""}</div>`,e.kind==="lipid"&&(t+='<p class="sub">Lipids are generated in your browser on the mesh of this membrane. Classes are PI, PE, PC, PA and PS, with ergosterol drawn as two beads. They sit at 0.65 nm² per lipid and are removed where they would overlap protein.</p>');i.innerHTML=t,f("inspClose").onclick=()=>Ke(null),e.kind==="inst"&&(f("loadAtoms").onclick=()=>In(e.info.type),f("hlType").onclick=()=>Te(e.info.type),f("centerOn").onclick=()=>{const n=w.state();n.target=e.info.p,w.flyTo(n,700)},f("hideType").onclick=()=>{te[e.info.type]^=2,R.device.queue.writeBuffer(R.typeFlags,0,te),Gt(),L=!0})}let Ve=null,Be=null;function Sn(){if(Be)return Be;const i="https://cdn.jsdelivr.net/npm/molstar@4.18.0/build/viewer/",e=t=>new Promise((n,s)=>{const o=document.createElement("link");o.rel="stylesheet",o.href=t+"molstar.css",document.head.appendChild(o);const r=document.createElement("script");r.src=t+"molstar.js",r.onload=()=>n(),r.onerror=()=>{o.remove(),r.remove(),s(new Error("Mol* failed to load"))},document.head.appendChild(r)});return Be=(window.__YC_EMBED__?e(i):e("molstar/").catch(()=>e(i))).catch(t=>{throw Be=null,t}),Be}async function In(i){const e=J.types[i];f("molModal").style.display="flex",f("molTitle").textContent=`${e.gene} · ${e.name} · loading…`;try{await Sn();const t=window.molstar;Ve?await Ve.plugin.clear():Ve=await t.Viewer.create("molHost",{layoutIsExpanded:!1,layoutShowControls:!1,layoutShowRemoteState:!1,layoutShowSequence:!0,layoutShowLog:!1,layoutShowLeftPanel:!1,viewportShowExpand:!1,viewportShowSelectionMode:!1,viewportShowAnimation:!1,collapseLeftPanel:!0,pdbProvider:"rcsb",emdbProvider:"rcsb"});let n="",s="";if(e.src==="PDB"&&e.pdb){n=`https://files.rcsb.org/download/${e.pdb}-assembly1.cif`,s=`https://www.rcsb.org/structure/${e.pdb}`;try{(await fetch(n,{method:"HEAD"})).ok||(n=`https://files.rcsb.org/download/${e.pdb}.cif`)}catch{n=`https://files.rcsb.org/download/${e.pdb}.cif`}}else{const r=e.model||e.uniprot;s=`https://alphafold.ebi.ac.uk/entry/${r}`;try{const u=await(await fetch(`https://alphafold.ebi.ac.uk/api/prediction/${r}`)).json();n=(Array.isArray(u)?u[0]:u)?.cifUrl}catch{}n||(n=e.model?`https://alphafold.ebi.ac.uk/files/${e.model}-model_v1.cif`:`https://alphafold.ebi.ac.uk/files/AF-${e.uniprot}-F1-model_v6.cif`)}await Ve.loadStructureFromUrl(n,"mmcif",!1),f("molTitle").textContent=`${e.gene} · ${e.name} · ${e.src==="PDB"?"PDB "+e.pdb+" (biological assembly)":e.model?`AlphaFold DB predicted ${e.oligomer}-mer (full length, incl. low-pLDDT regions)`:"AlphaFold DB model (full length, incl. low-pLDDT regions)"}`;const o=f("molLink");o.href=s,o.textContent="Open source page"}catch(t){f("molTitle").textContent=`Could not load the atomic model: ${t?.message||t}. Internet access to RCSB/AlphaFold DB is needed.`}}function _e(i){v.organelle=i,f("modeSection").classList.toggle("on",!i),f("modeOrg").classList.toggle("on",i),i&&!R.organelleMeshes.length&&(Me("Loading organelle surfaces (6.7 MB)…"),Pe("organelle_meshes.bin.gz").then(e=>{R.uploadOrganelleMeshes(e),L=!0})),L=!0}function $t(i){Nt=i;const e={low:[2.6,.75],med:[1.6,1],high:[1,1],ultra:[.6,1],auto:[1.6,1]},[t,n]=e[i];v.l0px=t,v.renderScale=n,R.resize(n),L=!0}function Pt(){const i=f("legend"),e=(n,s)=>`<div><span class="sw" style="background:${n}"></span>${s}</div>`,t=v.colorMode;t===0?i.innerHTML=[e(_[0],"cytosolic protein"),e(_[1],"nuclear protein"),e("#EE9BC0","ribosome protein"),e("#C9558F","rRNA"),e("#E27FA8","tRNA"),e("#B0437A","mRNA"),e("#F4B26A","histone"),e("#E07636","DNA"),e(_[23],"nuclear pore"),e("#7CC6C2","FG repeats"),e("#E0893A","actin filament"),e("#D9534F","Arp2/3 complex"),e("#8E6CB0","capping protein"),e("#2FA4A0","tropomyosin"),e("#74C476","microtubule"),e("#2E8B57","γ-tubulin ring"),e(_[14],"membrane protein"),e(_[5],"mitochondrial"),e(_[3],"ER lumen"),e("#F3C643","lipid head"),e("#F8E3A0","lipid tail"),e("#CDA866","β-glucan"),e("#C3C27A","mannan"),e(_[9],"wall protein"),e("#8B88C9","large complex")].join(""):t===1?i.innerHTML=e("#4C8CD9","experimental (PDB)")+e("#F5A84C","AlphaFold model"):t===2?i.innerHTML=e("#0053D6","pLDDT > 90")+e("#65CBF3","70–90")+e("#FFDB13","50–70")+e("#B8B8B8","experimental (PDB)"):t===3?i.innerHTML='<div style="grid-column: 1 / 3"><div style="height:8px;border-radius:3px;background:linear-gradient(90deg,#440154,#3b528b,#21918c,#5ec962,#fde725)"></div><div style="display:flex;justify-content:space-between"><span>1</span><span>10³</span><span>10⁶ copies</span></div></div>':t===5?i.innerHTML=[0,1,2,3,4,5,6,7,8,9,14,15,18,19].map(n=>e(_[n],Ue[n])).join(""):i.innerHTML=""}function ue(){const i=v.plane[3];f("planeD").value=String(Math.round(i)),f("planeDv").textContent=`${Number(i.toFixed(le<.5?2:1))} nm`;const e=v.slab;f("slab").value=String(Math.round(Math.log(e/5)/Math.log(4400/5)*1e3)),f("slabv").textContent=Bt(e),["axX","axY","axZ"].forEach((t,n)=>f(t).classList.toggle("on",Math.abs(v.plane[n])>.999)),ct=!0}function Bn(){const i=J.types.length;let e=0;for(const a of J.types)e+=a.copies;if(f("subCount").textContent=`${(e/1e6).toFixed(1)} M macromolecules of ${i.toLocaleString("en-US")} species, each a real PDB or AlphaFold structure`,window.__YC_PREVIEW__){const a=Oe.nInst.reduce((h,b)=>h+b,0);f("subCount").innerHTML=`<b>Preview:</b> ${Oe.n} of 6,416 tiles (${(a/1e3).toFixed(0)} k molecules) from the full ${(e/1e6).toFixed(1)} M-molecule cell. The full cell is in yeast_cell_viewer.zip.`,f("modeOrg").disabled=!0,f("modeOrg").title="Organelle surfaces are only in the full app"}const t=f("views");fe.forEach((a,h)=>{const b=document.createElement("button");b.textContent=`${(h+1)%10===0&&h<10?0:h<9?h+1:""}${h<10?" · ":""}${a.name}`,b.title=a.caption,b.onclick=()=>{we(),ke(h)},t.appendChild(b)}),f("tourBtn").onclick=()=>z?we():Ot();const n=f("search"),s=f("results");let o=0,r=[];const c=()=>{s.innerHTML=r.map((a,h)=>{const b=J.types[a];return`<div data-i="${a}" class="${h===o?"act":""}"><span class="nm"><b>${b.gene}</b> ${b.name}</span><span class="n">${ae(b.copies)}</span></div>`}).join(""),s.querySelectorAll("div").forEach(a=>a.addEventListener("click",()=>{Te(Number(a.dataset.i)),r=[],c(),n.value=""}))};n.addEventListener("input",()=>{r=Fn(n.value),o=0,c()}),n.addEventListener("keydown",a=>{a.key==="ArrowDown"?(o=Math.min(r.length-1,o+1),c(),a.preventDefault()):a.key==="ArrowUp"?(o=Math.max(0,o-1),c(),a.preventDefault()):a.key==="Enter"&&r.length?(Te(r[o]),Ut(r[o]),r=[],c(),n.value="",n.blur()):a.key==="Escape"&&(r=[],c(),n.blur()),a.stopPropagation()}),f("cutOn").onchange=a=>{v.cutOn=a.target.checked,L=!0};const u=a=>()=>{const h=[0,0,0];h[a]=$(w.fwd,[a===0?1:0,a===1?1:0,a===2?1:0])>=0?1:-1,v.cutOn=!0,f("cutOn").checked=!0,re(h,$(h,w.target))};f("axX").onclick=u(0),f("axY").onclick=u(1),f("axZ").onclick=u(2),f("alignView").onclick=()=>{v.cutOn=!0,f("cutOn").checked=!0,Ae()},f("followView").onclick=()=>{ge=!ge,f("followView").classList.toggle("on",ge),ge&&Ae()},f("planeD").oninput=a=>{const h=Number(a.target.value),b=de(),A=h-v.plane[3];w.target=D(w.target,E(b,A)),re(b,h)};const l=f("nudgeStep");[...l.options].some(a=>Number(a.value)===le)&&(l.value=String(le)),l.onchange=()=>{le=Number(l.value),ue(),l.blur();try{localStorage.setItem("yc.nudgeStep",String(le))}catch{}};for(const[a,h]of[["nudgeBack",-1],["nudgeFwd",1]]){const b=f(a);let A=0,C=0;const T=()=>{clearTimeout(A),clearInterval(C)};b.addEventListener("pointerdown",M=>{M.preventDefault(),T(),oe=null,f("sweepBtn").classList.remove("on"),Ye(h*le),A=window.setTimeout(()=>{C=window.setInterval(()=>Ye(h*le),60)},400)});for(const M of["pointerup","pointerleave","pointercancel"])b.addEventListener(M,T)}f("slab").oninput=a=>{const h=Number(a.target.value)/1e3;v.slab=5*Math.pow(4400/5,h),ue(),L=!0},f("sweepBtn").onclick=()=>{if(oe){oe=null,f("sweepBtn").classList.remove("on");return}const a=de();v.cutOn=!0,f("cutOn").checked=!0;const h=xe.G[0];oe={t0:performance.now(),d0:-h,d1:h,dur:3e4},w.target=D(w.target,E(a,-h-v.plane[3])),f("sweepBtn").classList.add("on")},f("flipBtn").onclick=()=>{const a=de(),h=w.state(),b=w.target,A=E(a,-1);re(A,-v.plane[3]),h.fwd=E(h.fwd,-1),h.target=b,w.flyTo(h,1200)},f("modeSection").onclick=()=>_e(!1),f("modeOrg").onclick=()=>_e(!0),f("colorMode").onchange=a=>{v.colorMode=Number(a.target.value),Pt(),L=!0},f("shaded").onchange=a=>{v.shaded=a.target.checked,L=!0},f("motion").onchange=a=>{v.motion=Number(a.target.value),a.target.blur(),L=!0},f("outlines").onchange=a=>{v.outline=a.target.checked?1:0,L=!0},Pt();const d=f("showList"),p=a=>h=>Rn(He[a],!h),m=[["Actin","#E0893A",p("actin")],["Microtubules","#74C476",p("microtubules")],["Ribosomes","#EE9BC0",p("ribosomes")],["tRNA","#E27FA8",p("tRNA")],["Nucleosomes","#F4B26A",p("nucleosomes")],["Nuclear pores",_[23],p("pores")],["mRNA","#B0437A",a=>Re(0,a)],["Linker DNA","#E07636",a=>Re(1,a)],["FG repeats","#7CC6C2",a=>Re(5,a)],["β-Glucan","#CDA866",a=>Re(2,a)],["Mannan","#C3C27A",a=>Re(3,a)],["Lipids","#F3C643",a=>{v.showLipids=a,L=!0}],["Membranes","#F8E3A0",a=>{v.showMeshes=a,L=!0}],["Soluble proteins",_[0],a=>We([0,1,2,3,4,5,6,7,8,10,11,12,13],a)],["Membrane proteins",_[14],a=>We([14,15,16,17,18,19,20,21],a)]],y=[],g=Ue.map((a,h)=>h),P=a=>{m.forEach(([,,b],A)=>{y[A].checked=a(A),b(a(A))});const h=m.every((b,A)=>a(A));We(g,h)};m.forEach(([a,h,b],A)=>{const C=document.createElement("label");C.className="chk",C.innerHTML=`<input type="checkbox" checked><span class="sw" style="background:${h}"></span><span class="nm">${a}</span><a class="only" title="Show only ${a.toLowerCase()}">only</a>`;const T=C.querySelector("input");y.push(T),T.onchange=()=>b(T.checked),C.querySelector(".only").onclick=M=>{M.preventDefault(),M.stopPropagation(),P(k=>k===A)},d.appendChild(C)}),f("showAll").onclick=()=>P(()=>!0),f("hideAll").onclick=()=>P(()=>!1);const x=f("regionList");Ue.forEach((a,h)=>{const b=document.createElement("label");b.className="chk",b.innerHTML=`<input type="checkbox" checked data-r="${h}"><span class="sw" style="background:${_[h]}"></span>${a}`,b.querySelector("input").onchange=A=>We([h],A.target.checked),x.appendChild(b)}),f("quality").onchange=a=>$t(a.target.value),f("detailOn").onchange=a=>{if(I.detailEnabled=a.target.checked,!I.detailEnabled)for(const h of I.loaded.values())h.det&&I.dropDetail(h);L=!0},f("shotBtn").onclick=qt,f("helpBtn").onclick=()=>f("helpModal").style.display="flex",f("aboutBtn").onclick=()=>f("aboutModal").style.display="flex",f("aboutBody").innerHTML=Mn,f("modelBtn").onclick=()=>{f("modelBody").innerHTML||(f("modelBody").innerHTML=An(kn)),f("modelModal").style.display="flex"},f("aboutModal").addEventListener("click",a=>{a.target.id==="aboutToModel"&&(a.preventDefault(),f("aboutModal").style.display="none",f("modelBtn").click())}),document.querySelectorAll(".modal").forEach(a=>a.addEventListener("click",h=>{(h.target===a||h.target.classList.contains("close"))&&(a.style.display="none")})),f("molClose").onclick=()=>f("molModal").style.display="none",f("panelToggle").onclick=()=>{f("panel").classList.toggle("collapsed"),f("panelToggle").textContent=f("panel").classList.contains("collapsed")?"›":"‹"},ue(),mt(),Dn()}function Rn(i,e){for(const t of i){const n=J.byId.get(t);n!==void 0&&(e?te[n]|=2:te[n]&=-3)}R.device.queue.writeBuffer(R.typeFlags,0,te),L=!0}function Re(i,e){e?v.freeHide&=~(1<<i):v.freeHide|=1<<i,L=!0}function We(i,e){for(const t of i)e?v.regionMask|=1<<t:v.regionMask&=~(1<<t);v.regionMask>>>=0,document.querySelectorAll("#regionList input").forEach(t=>{t.checked=(v.regionMask>>>Number(t.dataset.r)&1)===1}),L=!0}function mt(){let i=null;const e=new Set;S.addEventListener("contextmenu",n=>n.preventDefault()),S.addEventListener("pointerdown",n=>{S.setPointerCapture(n.pointerId),S.focus(),i={x:n.clientX,y:n.clientY,btn:n.button,moved:0,shift:n.shiftKey},S.classList.add("dragging"),z&&we()}),S.addEventListener("pointermove",n=>{if(!i){Nn(n.clientX,n.clientY);return}const s=n.clientX-i.x,o=n.clientY-i.y;if(i.x=n.clientX,i.y=n.clientY,i.moved+=Math.abs(s)+Math.abs(o),i.btn===2||i.btn===1||i.shift){const r=w.pxPerNmAtUnit(S.clientHeight)/w.dist;w.pan(s,o,r)}else w.orbit(s*.005,o*.005),ge&&Ae();w.anim=null,L=!0});const t=n=>{if(!i)return;const s=i.moved<4&&i.btn===0;i=null,S.classList.remove("dragging"),s&&Qe(n.offsetX,n.offsetY).then(o=>{const r=pt(o.id);Ke(r.kind==="inst"||r.kind==="lipid"||r.kind==="free"||r.kind==="band"?r:null)})};S.addEventListener("pointerup",t),S.addEventListener("pointercancel",t),S.addEventListener("dblclick",n=>{Qe(n.offsetX,n.offsetY).then(s=>{if(s.depth<=0)return;const o=w.near(),r=o/s.depth,c=S.clientWidth/S.clientHeight,u=w.ray(n.offsetX/S.clientWidth,n.offsetY/S.clientHeight,c),l=r/$(u,w.fwd),d=D(w.eye,E(u,l)),p=w.state();p.target=d,p.dist=Math.max(w.minDist*4,w.dist*.55),w.flyTo(p,900)})}),S.addEventListener("wheel",n=>{n.preventDefault();const s=n.deltaY||n.deltaX,o=n.deltaMode===1?s*30:s;if(n.altKey)Ye(Math.sign(o)*le*Math.min(4,Math.abs(o)/100));else if(n.shiftKey){const r=de(),c=Math.max(1,Math.min(v.slab,w.dist*.05))*Math.sign(o)*Math.min(4,Math.abs(o)/60+.5);w.target=D(w.target,E(r,c)),re(r,v.plane[3]+c)}else{const r=Math.exp(o*.0015*(n.metaKey?.125:1)),c=S.clientWidth/S.clientHeight,u=w.ray(n.offsetX/S.clientWidth,n.offsetY/S.clientHeight,c),l=w.dist,d=D(w.eye,E(u,l/$(u,w.fwd))),p=w.dist;w.zoom(r);const m=1-w.dist/p;w.target=D(w.target,E(ce(d,w.target),m))}w.anim=null,z&&we(),L=!0},{passive:!1}),window.addEventListener("keydown",n=>{if(n.target.tagName==="INPUT"||n.target.tagName==="SELECT")return;const s=n.key.toLowerCase();if(e.add(s),s==="c")v.cutOn=!v.cutOn,f("cutOn").checked=v.cutOn,L=!0;else if(s==="v")v.cutOn=!0,f("cutOn").checked=!0,Ae();else if(s==="[")v.slab=Math.max(5,v.slab/1.25),ue(),L=!0;else if(s==="]")v.slab=Math.min(4400,v.slab*1.25),ue(),L=!0;else if(s===","||s==="."||n.code==="Comma"||n.code==="Period"){const o=s==="."||n.code==="Period";Ye((o?1:-1)*(n.altKey?le:Math.max(2,v.slab*.25)))}else if(s==="o")_e(!v.organelle);else if(s==="t")z?we():Ot();else if(s==="p")qt();else if(s==="escape")Ke(null),Te(-1);else if(/^[0-9]$/.test(s)){const o=s==="0"?9:Number(s)-1;o<fe.length&&(we(),ke(o))}}),window.addEventListener("keyup",n=>e.delete(n.key.toLowerCase())),window.addEventListener("blur",()=>e.clear()),mt.keys=e,window.addEventListener("resize",()=>{R.resize(v.renderScale),L=!0})}let Mt=0,rt=!1;function Nn(i,e){const t=performance.now();rt||t-Mt<90||(Mt=t,rt=!0,L=!0,Qe(i,e).then(n=>{rt=!1;const s=pt(n.id),o=f("hover");if(!s.title||s.kind==="none"){o.style.display="none";return}let r=`<b>${s.title}</b>`;if(s.kind==="inst"){const c=J.types[s.info.type];r+=` · ${c.name} <span style="color:#aaa">· ${s.sub}</span>`}else r+=` <span style="color:#aaa">· ${s.sub}</span>`;o.innerHTML=r,o.style.display="block",o.style.left=`${i+14}px`,o.style.top=`${e+12}px`}))}S.addEventListener("pointerleave",()=>f("hover").style.display="none");let ct=!0,At=0,he=null;const Ne={data:null};function Dn(){const i=f("minimap");i.addEventListener("click",e=>{if(!he)return;const t=i.getBoundingClientRect(),n=((e.clientX-t.left)/t.width*2-1)*he.ext,s=(1-(e.clientY-t.top)/t.height*2)*he.ext,o=D(he.c,D(E(he.u,n),E(he.v,s))),r=w.state();r.target=o,w.flyTo(r,1200),z&&we()})}function dt(i){const e=f("minimap"),t=e.getContext("2d"),n=e.width,s=e.height,o=v.cutOn?de():w.fwd,r=v.cutOn?v.plane[3]:$(w.fwd,w.target);if(ct&&i-At>250){ct=!1,At=i;let x=ce(w.right,E(o,$(w.right,o)));ye(x)<.001&&(x=K(pe(o,[0,0,1]))),x=K(x);let a=pe(o,x);$(a,w.up)<0&&(a=E(a,-1));const h=2300,b=E(o,r);he={u:x,v:a,c:b,ext:h};const A=120,C=t.createImageData(A,A),T=[_[80],_[9],"#D2BE8E",_[83],_[3],_[3],_[1],_[2],_[4],_[6],_[5],"#E9C46A",_[7],_[8],"#C9C2DC"].map(M=>{const k=parseInt(M.slice(1),16);return[k>>16&255,k>>8&255,k&255]});T[3]=[196,214,235];for(let M=0;M<A;M++)for(let k=0;k<A;k++){const N=((k+.5)/A*2-1)*h,O=(1-(M+.5)/A*2)*h,q=D(b,D(E(x,N),E(a,O))),j=Rt.label(q[0],q[1],q[2]),H=T[j],W=4*(M*A+k);C.data[W]=H[0],C.data[W+1]=H[1],C.data[W+2]=H[2],C.data[W+3]=j===0?0:255}Ne.data=C}if(!Ne.data||!he)return;const c=dt.tmp||(dt.tmp=document.createElement("canvas"));c.width=Ne.data.width,c.height=Ne.data.height,c.getContext("2d").putImageData(Ne.data,0,0),t.clearRect(0,0,n,s),t.imageSmoothingEnabled=!0,t.drawImage(c,0,0,n,s);const u=he,l=x=>{const a=ce(x,u.c);return[($(a,u.u)/u.ext+1)/2*n,(1-$(a,u.v)/u.ext)/2*s]},d=S.clientWidth/S.clientHeight,p=w.dist*Math.tan(w.fovy/2),m=p*d,y=[[-1,-1],[1,-1],[1,1],[-1,1]].map(([x,a])=>l(D(w.target,D(E(w.right,x*m),E(w.up,a*p)))));t.strokeStyle="#FF9400",t.lineWidth=3,t.beginPath(),y.forEach(([x,a],h)=>h?t.lineTo(x,a):t.moveTo(x,a)),t.closePath(),t.stroke();const[g,P]=l(w.target);t.fillStyle="#FF9400",t.beginPath(),t.arc(g,P,5,0,7),t.fill(),f("mmLabel").textContent=v.cutOn?`plane ${Math.round(r)} nm`:"view plane"}function qt(){const i=v.renderScale;R.resize(2*i),v.renderScale=2*i,Wt(!0),S.toBlob(e=>{if(e){const t=document.createElement("a");t.href=URL.createObjectURL(e),t.download=`yeast_cell_${new Date().toISOString().replace(/[:.]/g,"-").slice(0,19)}.png`,t.click(),setTimeout(()=>URL.revokeObjectURL(t.href),5e3)}v.renderScale=i,R.resize(i),L=!0},"image/png"),Me("Screenshot saved (2× resolution)")}function On(){const i=I.totals(),e=R.stats,t=ne.length>3?1e3/(ne.reduce((l,d)=>l+d)/ne.length):0,n=I.pendingForView();f("hud").innerHTML=`<b>${ae(i.inst)}</b> molecules in ${i.tiles} tiles${n?` · <span style="color:#FF9400">streaming ${n}</span>`:""}${I.starved?' · <span style="color:#FF9400">memory budget reached (thinner slab loads more)</span>':""}<br>
    drawn: <b>${ae(e.ninst)}</b> in view · ${ae(e.l0*64)} residue beads (≤) · ${ae(e.l1*16)} LOD1 · ${ae(e.l2)} LOD2 · ${ae(e.lip)} lipids · ${ae(e.fre)} polymer beads<br>
    ${t?`${t.toFixed(0)} fps · `:""}LOD bead ≥ ${v.l0px.toFixed(1)} px · GPU pools ${R.gpuMemoryMB().toFixed(0)} MB · ${(I.bytesLoaded/1048576).toFixed(0)} MB streamed${R.adapterInfo?" · "+R.adapterInfo:""}`;const s=R.width/S.clientWidth,o=w.pxPerNmAtUnit(R.height)/w.dist/s,r=120/o,c=Math.pow(10,Math.floor(Math.log10(r))),u=[1,2,5,10].map(l=>l*c).reduce((l,d)=>Math.abs(d-r)<Math.abs(l-r)?d:l);f("scalebar").querySelector(".bar").style.width=`${u*o}px`,f("scalebar").querySelector("span").textContent=Bt(u)}let Vt=0,Je=!1;function Wt(i=!1){const e=performance.now();R.frame(w,v,xe.G,I.meshDraws(),i||e-R.lastStatsT>400),R.lastStatsT=performance.now(),Je=!0,R.device.queue.onSubmittedWorkDone().then(()=>{Vt=performance.now()-e,Je=!1})}let kt=0;function jt(){I.update(w,S.clientWidth/S.clientHeight,v.cutOn,v.plane,v.slab,v.organelle,R.height),ht=performance.now(),(I.inflight||I.detInflight)&&(L=!0)}setInterval(()=>{I&&performance.now()-ht>150&&jt()},150);function Ht(i){requestAnimationFrame(Ht);const e=ot?i-ot:16;ot=i;const t=mt.keys;if(t&&t.size){const n=w.dist*.9*e/1e3;let s=[0,0,0];t.has("w")&&(s=D(s,E(w.fwd,n))),t.has("s")&&(s=D(s,E(w.fwd,-n))),t.has("a")&&(s=D(s,E(w.right,-n))),t.has("d")&&(s=D(s,E(w.right,n))),t.has("q")&&(s=D(s,E(w.up,-n))),t.has("e")&&(s=D(s,E(w.up,n))),t.has("z")&&w.roll(-e*.0015),t.has("x")&&w.roll(e*.0015),(ye(s)>0||t.has("z")||t.has("x"))&&(w.fly(s),L=!0,ge&&Ae())}if(w.step(i)&&(L=!0,ge&&!z&&Ae()),oe){const n=Math.min(1,(i-oe.t0)/oe.dur),s=de(),o=oe.d0+(oe.d1-oe.d0)*n;w.target=D(w.target,E(s,o-v.plane[3])),v.plane[3]=o,ue(),n>=1&&(oe=null,f("sweepBtn").classList.remove("on")),L=!0}if(Ln(i),v.motion>0&&(L=!0),i-ht>90&&jt(),L&&!Je){if(L=!1,Wt(),ne.push(e),ne.length>30&&ne.shift(),Nt==="auto"&&ne.length>=20){const n=ne.reduce((s,o)=>s+o)/ne.length;n>45&&v.l0px<3.5?(v.l0px=Math.min(3.5,v.l0px*1.15),ne.length=0):n<22&&v.l0px>.9&&(v.l0px=Math.max(.9,v.l0px/1.1),ne.length=0)}setTimeout(()=>{L=L||R.statsBusy},0)}else L||(ne.length=0);i-kt>300&&(kt=i,On(),dt(i),be>=0&&I.inflight&&zt())}Cn();
