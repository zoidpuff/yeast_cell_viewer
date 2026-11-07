(function(){const e=document.createElement("link").relList;if(e&&e.supports&&e.supports("modulepreload"))return;for(const s of document.querySelectorAll('link[rel="modulepreload"]'))n(s);new MutationObserver(s=>{for(const o of s)if(o.type==="childList")for(const r of o.addedNodes)r.tagName==="LINK"&&r.rel==="modulepreload"&&n(r)}).observe(document,{childList:!0,subtree:!0});function t(s){const o={};return s.integrity&&(o.integrity=s.integrity),s.referrerPolicy&&(o.referrerPolicy=s.referrerPolicy),s.crossOrigin==="use-credentials"?o.credentials="include":s.crossOrigin==="anonymous"?o.credentials="omit":o.credentials="same-origin",o}function n(s){if(s.ep)return;s.ep=!0;const o=t(s);fetch(s.href,o)}})();const X=1024,Q=4096,ke=`
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
  modes: vec4u,           // x colorMode, y regionMask, z hidden free-sphere kinds, w highlightActive
  sel: vec4u,             // x selected id, y selected type+1, z shaded(1), w organelle mode
  misc: vec4f,            // x QSTEP, y nucleic bead radius, z outline strength, w dpr
  cell: vec4f,            // x R outer, y R mannan inner, z R PM, w near
  caps: vec4u,            // x capL0 (offset of L1 items in list01), y motion amount x 1000, z motion cap nm x 1000,
                          // w membranes-only view (whole-cell surfaces, palette 66 + membrane, no cut)
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
`,Ut=`
struct TypeInfo {
  beadOff: u32, beadCount: u32, clusterOff: u32, clusterCount: u32,
  stOff: u32, stCount: u32, br: f32, req: f32,
  color0: u32, color1: u32, flags: u32, plddt: f32,
  modeOff: u32, coefOff: u32, nConf: u32,             // pool-ensemble motion modes (nConf 0: none)
  modColor: u32,                                       // packed colour of the STRING interaction module (0: none)
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
`,un=ke+Ut+`
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
`,hn=ke+`
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
`,pn=`
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
`,Gt=`
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
  // spheres a few pixels wide show only their bright, front-facing centre: fade to the sphere's mean shade there,
  // so residue beads (L0) are not lighter than the larger spheres of the coarser levels at the LOD boundary
  let rpx = in.radius * u.eye.w / max(length(in.center), 1e-3);
  let sharp = smoothstep(0.75, 3.0, rpx);
  if (u.sel.z == 1u) {
    let l = normalize(vec3f(0.35, 0.45, 1.0));
    shade = 0.45 + 0.55 * max(dot(nrm, l), 0.0) + 0.12 * pow(max(nrm.z, 0.0), 8.0);
    shade = mix(0.80, shade, sharp);       // mean over the projected disc: n.l ~ 2/3 l.z, pow(n.z, 8) ~ 1/5
  } else {
    shade = mix(0.977, 0.93 + 0.07 * max(nrm.z, 0.0), sharp);
  }
  if (cap) { shade = 0.9; }
  o.color = vec4f(in.color * shade, 1.0);
  o.id = in.id;
  return o;
}
`,mn=`
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
  } else if (mode == 6u) {   // STRING interaction module
    if (T.modColor != 0u) { c = unpackColor(T.modColor); } else { c = vec3f(0.82, 0.82, 0.80); }
    if (kind == 1u) { c = c * 0.8; }
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
`,Pt=ke+Ut+Gt+mn,fn=Pt+`
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
`,gn=Pt+`
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
`,bn=Pt+`
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
`,yn=ke+Gt+`
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
`,wn=ke+`
@group(0) @binding(5) var<storage, read> palette: array<vec4f>;
struct MP { slot: u32, organelle: u32, offset0: f32, offset1: f32 };   // organelle: 0 tile, 1 organelle mode, 2 under lipids, 3 membranes only
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
  let outside = !ff;      // membrane triangles are wound inward (tiles3d marching cubes): the outside is the back face
  // membranes only: the plasma membrane is seen from inside (its outside dropped) so the organelles show through it
  if (mp.organelle == 3u && ((in.mem == 0u && outside) || in.mem == 9u)) { discard; }
  if (u.cut.x > 0.5) {
    let s = sPlaneW(in.prel);
    if (s < 0.0 || (mp.organelle != 1u && s > u.cut.y)) { discard; }
  }
  let vn = normalize((u.view * vec4f(in.nrm, 0.0)).xyz);
  var c = palette[41u].xyz;
  if (mp.organelle == 1u) { c = palette[56u + in.mem].xyz; } else if (mp.organelle == 3u) { c = palette[66u + in.mem].xyz; }
  let l = abs(vn.z);
  var shade = 0.82 + 0.18 * l;
  if (mp.organelle == 1u || mp.organelle == 3u) { shade = 0.55 + 0.45 * l; if (!outside) { shade *= 0.8; } }
  o.color = vec4f(c * shade, 1.0);
  o.id = 0x80000000u | (mp.slot << 5u) | (in.mem << 1u) | in.leaf;
  return o;
}
`,vn=ke+`
@group(0) @binding(5) var<storage, read> palette: array<vec4f>;
struct Prim { q0: vec4f, q1: vec4f, tc: vec4u };      // type, class
struct Geo { G0: vec4f, G1: vec4f, G2: vec4f, consts: vec4f, grid: vec4f,       // consts: PORE_R, PERI_W, IMS_W, nPores
             nax: vec4f, nr0: vec4f, nr1: vec4f, nr2: vec4f };                    // nucleus semi-axes, rotation rows
@group(1) @binding(0) var<storage, read> prims: array<Prim>;
@group(1) @binding(1) var<storage, read> pores: array<vec4f>;      // per pore: centre (perinuclear midsurface), NE normal
@group(1) @binding(2) var<storage, read> gridOff: array<u32>;   // prim list offsets (n^3+1) then pore list offsets
@group(1) @binding(3) var<storage, read> gridIdx: array<u32>;
@group(1) @binding(4) var<uniform> geo: Geo;
@group(1) @binding(5) var<storage, read> tileState: array<u32>;   // 26^3 (data.TGRID), bit0 = molecules loaded
@group(1) @binding(6) var<storage, read> fields: array<f32>;       // real-shape distance fields (fields.ts layout)

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
// distance field k (cell3d.field_sdf): trilinear inside a stored brick, +-clamp elsewhere
fn fieldSdf(p: vec3f, k: u32) -> f32 {
  let o = 1u + 12u * k; let h = fields[o + 3u]; let B = u32(fields[o + 7u]);
  let f = (p - vec3f(fields[o], fields[o + 1u], fields[o + 2u])) / h;
  let nb = vec3u(u32(fields[o + 4u]), u32(fields[o + 5u]), u32(fields[o + 6u]));
  if (any(f < vec3f(0.0)) || any(f >= vec3f(nb * B))) { return fields[o + 9u]; }
  let ii = vec3u(floor(f)); let bb = ii / B;
  let bi = i32(fields[u32(fields[o + 10u]) + (bb.x * nb.y + bb.y) * nb.z + bb.z]);
  if (bi == -1) { return -fields[o + 8u]; }
  if (bi == -2) { return fields[o + 9u]; }
  let n1 = B + 1u; let sx = n1 * n1; let sy = n1;
  let t = f - vec3f(ii); let l = ii - bb * B;
  let b = u32(fields[o + 11u]) + u32(bi) * n1 * n1 * n1 + (l.x * n1 + l.y) * n1 + l.z;
  let c00 = mix(fields[b], fields[b + sx], t.x); let c01 = mix(fields[b + 1u], fields[b + sx + 1u], t.x);
  let c10 = mix(fields[b + sy], fields[b + sx + sy], t.x); let c11 = mix(fields[b + sy + 1u], fields[b + sx + sy + 1u], t.x);
  return mix(mix(c00, c10, t.y), mix(c01, c11, t.y), t.z);
}
// depth below the outer cell surface (negative outside)
fn outerDepth(p: vec3f) -> f32 {
  if (fields[0] > 0.5) { return -fieldSdf(p, 0u); }
  return geo.G0.x - length(p);
}
// cortical-ER patch: a sheet at constant depth below the plasma membrane (pm: signed PM distance, cell3d.sdf_cer)
fn sdfCer(p: vec3f, pr: Prim, pm: f32) -> f32 {
  let r = length(p);
  let c = clamp(dot(p, pr.q0.xyz) / max(r, 1e-9), -1.0, 1.0);
  let s = max(acos(c) - pr.q1.x, 0.0) * r;
  let dr = (geo.G0.z - pr.q1.y) - pm;
  return sqrt(s * s + dr * dr) - pr.q1.z;
}
fn nucSdf(p: vec3f) -> f32 {
  if (fields[0] > 2.5) { return fieldSdf(p, 2u); }
  return nucEllipsoid(p);
}
// signed distance to the ONM midsurface (ellipsoid), first-order approximation as in cell3d.nuc_sdf
fn nucEllipsoid(p: vec3f) -> f32 {
  let v = p - vec3f(geo.G0.w, geo.G1.x, geo.G1.y);
  let l = vec3f(dot(geo.nr0.xyz, v), dot(geo.nr1.xyz, v), dot(geo.nr2.xyz, v)) / geo.nax.xyz;
  let k0 = length(l); let k1 = length(l / geo.nax.xyz);
  if (k1 < 1e-12) { return -min(geo.nax.x, min(geo.nax.y, geo.nax.z)); }
  return k0 * (k0 - 1.0) / k1;
}
// pore whose axis (through its centre along the NE normal) passes within rc of p: returns rho, or -1
fn poreRho(p: vec3f, cell: u32, rc: f32) -> f32 {
  let n3 = u32(geo.grid.x * geo.grid.x * geo.grid.x);
  let po0 = gridOff[n3 + 1u + cell]; let po1 = gridOff[n3 + 1u + cell + 1u];
  for (var k = po0; k < po1; k++) {
    let i = gridIdx[k];
    let v = p - pores[2u * i].xyz; let nn = pores[2u * i + 1u].xyz;
    let z = dot(v, nn);
    if (abs(z) > 150.0) { continue; }
    let rho = length(v - z * nn);
    if (rho < rc) { return rho; }
  }
  return -1.0;
}
fn gridCell(p: vec3f) -> u32 {
  let n = u32(geo.grid.x); let ext = geo.grid.y; let h = 2.0 * ext / f32(n);
  let g = clamp(vec3i(floor((p + vec3f(ext)) / h)), vec3i(0), vec3i(i32(n) - 1));
  return (u32(g.x) * n + u32(g.y)) * n + u32(g.z);
}
fn memSdfs(p: vec3f) -> Mem {
  var m: Mem;
  let G = geo;
  let PERI = G.consts.y; let PORE_R = G.consts.x;
  let dep = outerDepth(p);
  m.d[0] = dep - (G.G0.x - G.G0.z);     // PM
  m.d[9] = dep;                         // outer
  let hp = nucSdf(p) + 0.5 * PERI;
  var d = abs(hp) - 0.5 * PERI;
  let rc = PORE_R + 0.5 * PERI;
  let cell = gridCell(p); m.cell = cell;
  if (abs(hp) < rc + PERI) {
    let rho = poreRho(p, cell, rc);
    if (rho >= 0.0) { d = sqrt((rho - rc) * (rho - rc) + hp * hp) - 0.5 * PERI; }
  }
  m.d[2] = d;
  m.d[1] = 1e9; m.d[3] = 1e9; m.d[4] = 1e9; m.d[6] = 1e9; m.d[7] = 1e9; m.d[8] = 1e9;
  if (fields[0] > 1.5) { m.d[3] = fieldSdf(p, 1u); }
  var im = 1e9; var cr = 1e9;
  let o0 = gridOff[cell]; let o1 = gridOff[cell + 1u];
  for (var k = o0; k < o1; k++) {
    let pr = prims[gridIdx[k]];
    var dd: f32;
    if (pr.tc.x == 3u) { dd = sdfCer(p, pr, m.d[0]); } else { dd = sdfPrim(p, pr); }
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
  if (m.d[9] < 0.0) { return 0u; }
  if (m.d[9] < G.G0.x - G.G0.y) { return 1u; }
  if (m.d[0] < 0.0) { return 2u; }
  if (m.d[1] < 0.0) { return 4u; }
  if (m.d[2] < 0.0) { return 5u; }
  let ns = nucSdf(p);
  if (ns < -0.5 * PERI) {
    if (ns > -PERI - 8.0 && poreRho(p, m.cell, PORE_R + 0.5 * PERI) >= 0.0) { return 14u; }
    let w = p - vec3f(G.G1.w, G.G2.x, G.G2.y);
    if (ns < -PERI && dot(w, w) < G.G2.z * G.G2.z) { return 7u; }
    return 6u;
  }
  if (ns < 0.0) { return 14u; }
  if (m.d[4] < 0.0) { if (m.matrix < 0.0) { return 10u; } return 9u; }
  if (m.d[3] < 0.0) { return 8u; }
  if (m.d[8] < 0.0) { return 11u; }
  if (m.d[6] < 0.0) { return 12u; }
  if (m.d[7] < 0.0) { return 13u; }
  return 3u;
}
fn tileLoaded(p: vec3f) -> bool {
  let g = vec3i(floor(p / 200.0)) + vec3i(13);
  if (any(g < vec3i(0)) || any(g > vec3i(25))) { return false; }
  return (tileState[(u32(g.x) * 26u + u32(g.y)) * 26u + u32(g.z)] & 1u) != 0u;
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
        if (d > -2.0 && d < 0.0) { c = palette[41u].xyz; id = 0x200u + k; }
        continue;
      }
      if (abs(d) < 2.5) { c = palette[41u].xyz; id = 0x200u + k; }
    }
  }
  (*o).color = vec4f(c, 1.0);
  (*o).id = id;
}

fn fillDepth(pv0: vec3f) -> f32 {
  let clip = u.proj * vec4f((u.view * vec4f(pv0, 1.0)).xyz, 1.0);
  return clamp(clip.z / clip.w, 0.0, 1.0);
}

@fragment
fn fsFill(in: VO) -> FOut {
  var o: FOut;
  let pw = u.invViewProj * vec4f(in.ndc, 1.0, 1.0);    // near plane point (camera-relative)
  let rd = normalize(pw.xyz / pw.w);
  let eye = u.eye.xyz;
  let R = u.cell.x;
  // ray vs the enclosing sphere of the cell (centre at origin in world)
  let b = dot(rd, eye);
  let c = dot(eye, eye) - R * R;
  let h = b * b - c;
  o.color = vec4f(palette[80u].xyz, 1.0); o.id = 0u; o.depth = 0.0;
  if (h < 0.0 || u.caps.w == 1u) { return o; }      // no ray hit, or the membranes-only view (meshes only)
  let sh = sqrt(h);
  var t0 = max(-b - sh, 0.0);
  let t1 = -b + sh;
  if (t1 <= 0.0) { return o; }
  var onCut = false;
  var tEnd = t1;          // the exterior trace ends here (where the ray leaves the cell or enters the removed half)
  if (u.cut.x > 0.5) {
    let dn = dot(u.plane.xyz, rd);
    let sE = u.planeV.w;
    if (abs(dn) < 1e-6) { if (sE < 0.0) { return o; } }
    else {
      let tp = -sE / dn;
      if (dn > 0.0) { if (tp > t1) { return o; } if (tp > t0) { t0 = tp; onCut = true; } }
      else { if (tp < t0) { return o; } tEnd = min(t1, tp); }
    }
    if (u.sel.w == 1u) {
      // organelle mode: cut face only (meshes provide the rest)
      if (!onCut) { o.color = vec4f(palette[80u].xyz, 1.0); return o; }
    }
  }
  var t = t0;
  var p = eye + rd * t;
  if (onCut && outerDepth(p) >= 0.0) {
    if (u.sel.w == 0u && tileLoaded(p)) {
      // molecules are drawn here: show the compartment backdrop behind the visible slab
      t = t + u.cut.y / max(dot(u.plane.xyz, rd), 0.05);
      if (t >= t1) { return o; }
      p = eye + rd * t;
    }
    shadeLabel(p, &o, true);
    o.depth = fillDepth(rd * t);
    return o;
  }
  // exterior: the outer cell surface, found by sphere tracing the (clamped) distance field; with a cut it is also
  // what lies beyond a cut point outside the real outline (but inside the enclosing sphere)
  if (onCut && u.sel.w == 1u) { return o; }      // organelle mode: the meshes show it
  if (!onCut && outerDepth(eye) > 0.0) { return o; }   // camera inside the cell without a cut: nothing analytic to show
  for (var k = 0; k < 128; k++) {
    let d = -outerDepth(p);
    if (d < 0.25) { break; }
    t = t + max(d, 0.5);
    if (t >= tEnd) { return o; }
    p = eye + rd * t;
  }
  let ps = p;
  if (u.sel.w == 0u && tileLoaded(p) && (u.cut.x < 0.5 || dot(u.plane.xyz, p) - u.plane.w < u.cut.y)) {
    // wall molecules (mannoproteins, mannan) are drawn: the backdrop is the top of the glucan layer, at constant depth
    // below the surface, so the view does not depend on the viewing angle; grazing rays that miss it see background
    let D = u.cell.x - u.cell.y;
    for (var k = 0; k < 64; k++) {
      let d = D - outerDepth(p);
      if (d < 0.25) { break; }
      t = t + max(d, 0.5);
      if (t >= tEnd) { return o; }
      p = eye + rd * t;
    }
  }
  // one wall colour, lit like the molecules by the smoothed surface normal (no compartment labels or membrane
  // bands: those are cut-face features, and on a traced surface they would trace contours of the field's noise)
  let e = 10.0;
  let g = vec3f(outerDepth(ps + vec3f(e, 0.0, 0.0)) - outerDepth(ps - vec3f(e, 0.0, 0.0)),
                outerDepth(ps + vec3f(0.0, e, 0.0)) - outerDepth(ps - vec3f(0.0, e, 0.0)),
                outerDepth(ps + vec3f(0.0, 0.0, e)) - outerDepth(ps - vec3f(0.0, 0.0, e)));
  var nn = normalize(ps);
  if (dot(g, g) > 1e-6) { nn = -normalize(g); }
  let vn = (u.view * vec4f(nn, 0.0)).xyz;
  var shade = 0.93 + 0.07 * max(vn.z, 0.0);
  if (u.sel.z == 1u) {
    let l = normalize(vec3f(0.35, 0.45, 1.0));
    shade = 0.5 + 0.5 * max(dot(vn, l), 0.0) + 0.06 * pow(max(vn.z, 0.0), 8.0);
  }
  o.color = vec4f(palette[82u].xyz * shade, 1.0);
  o.id = 0x101u;
  o.depth = fillDepth(rd * t);
  return o;
}
`,xn=`
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
`,Mn=ke+`
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
`,_=(i,e)=>[i[0]+e[0],i[1]+e[1],i[2]+e[2]],ue=(i,e)=>[i[0]-e[0],i[1]-e[1],i[2]-e[2]],L=(i,e)=>[i[0]*e,i[1]*e,i[2]*e],V=(i,e)=>i[0]*e[0]+i[1]*e[1]+i[2]*e[2],ge=(i,e)=>[i[1]*e[2]-i[2]*e[1],i[2]*e[0]-i[0]*e[2],i[0]*e[1]-i[1]*e[0]],Ae=i=>Math.hypot(i[0],i[1],i[2]),J=i=>{const e=Ae(i)||1;return[i[0]/e,i[1]/e,i[2]/e]},ft=(i,e,t)=>[i[0]+(e[0]-i[0])*t,i[1]+(e[1]-i[1])*t,i[2]+(e[2]-i[2])*t];function Ue(i,e,t){const n=Math.cos(t),s=Math.sin(t),o=ge(e,i),r=V(e,i)*(1-n);return[i[0]*n+o[0]*s+e[0]*r,i[1]*n+o[1]*s+e[1]*r,i[2]*n+o[2]*s+e[2]*r]}function Pn(i,e){const t=new Float32Array(16);for(let n=0;n<4;n++)for(let s=0;s<4;s++){let o=0;for(let r=0;r<4;r++)o+=i[r*4+s]*e[n*4+r];t[n*4+s]=o}return t}function An(i){const e=Array.from(i),t=new Float32Array(16),[n,s,o,r,a,d,l,c,u,m,f,T,w,M,p,x]=e,h=n*d-s*a,g=n*l-o*a,b=n*c-r*a,A=s*l-o*d,C=s*c-r*d,F=o*c-r*l,O=u*M-m*w,k=u*p-f*w,E=u*x-T*w,Y=m*p-f*M,$=m*x-T*M,W=f*x-T*p;let z=h*W-g*$+b*Y+A*E-C*k+F*O;return z&&(z=1/z,t[0]=(d*W-l*$+c*Y)*z,t[1]=(o*$-s*W-r*Y)*z,t[2]=(M*F-p*C+x*A)*z,t[3]=(f*C-m*F-T*A)*z,t[4]=(l*E-a*W-c*k)*z,t[5]=(n*W-o*E+r*k)*z,t[6]=(p*b-w*F-x*g)*z,t[7]=(u*F-f*b+T*g)*z,t[8]=(a*$-d*E+c*O)*z,t[9]=(s*E-n*$-r*O)*z,t[10]=(w*C-M*b+x*h)*z,t[11]=(m*b-u*C-T*h)*z,t[12]=(d*k-a*Y-l*O)*z,t[13]=(n*Y-s*k+o*O)*z,t[14]=(M*g-w*A-p*h)*z,t[15]=(u*A-m*g+f*h)*z),t}function at(i){const e=parseInt(i.replace("#",""),16);return(e>>16&255|(e>>8&255)<<8|(e&255)<<16|4278190080)>>>0}function kn(i){return i<.5?4*i*i*i:1-Math.pow(-2*i+2,3)/2}const qt=new URL("data/",document.baseURI).href,Tn=fetch(qt+"version.json",{cache:"no-store"}).then(i=>i.ok?i.json():{v:""}).then(i=>String(i.v??"")).catch(()=>""),Cn=i=>i?`?v=${encodeURIComponent(i)}`:"";function $t(i){const e=atob(i),t=new Uint8Array(e.length);for(let n=0;n<e.length;n++)t[n]=e.charCodeAt(n);return t}async function lt(i){const e=new Uint8Array(i,0,Math.min(2,i.byteLength));if(e[0]===31&&e[1]===139){const t=new DecompressionStream("gzip"),n=new Blob([i]).stream().pipeThrough(t);return await new Response(n).arrayBuffer()}return i}async function Me(i,e){const t=window.__YC_EMBED__;if(t&&t[i]!==void 0)return lt($t(t[i]).buffer);const n=await fetch(qt+i+Cn(await Tn));if(!n.ok)throw new Error(`HTTP ${n.status} for ${i}`);if(!e||!n.body)return lt(await n.arrayBuffer());const s=Number(n.headers.get("content-length")||0),o=n.body.getReader(),r=[];let a=0;for(;;){const{done:c,value:u}=await o.read();if(c)break;r.push(u),a+=u.length,e(a,s)}const d=new Uint8Array(a);let l=0;for(const c of r)d.set(c,l),l+=c.length;return lt(d.buffer)}async function Ce(i){const e=await Me(i);return JSON.parse(new TextDecoder().decode(e))}const $e=["cytosol","nucleoplasm","nucleolus","ER lumen","vacuole lumen","mitochondrial matrix","mitochondrial IMS","Golgi / endosome","peroxisome","cell wall","PM peripheral","NE peripheral","ER peripheral","lipid droplet surface","plasma membrane","ER membrane","nuclear envelope","vacuolar membrane","mito outer membrane","mito inner membrane","Golgi membrane","peroxisomal membrane","chromatin (nucleosome)","nuclear pore"],St=["plasma membrane","ER membrane","nuclear envelope","vacuolar membrane","mito outer membrane","mito inner membrane","Golgi/endosome membrane","peroxisomal membrane","lipid droplet monolayer","cell surface"],En=["exterior","cell wall (mannoprotein layer)","cell wall (glucan/chitin layer)","cytosol","cortical ER lumen","perinuclear space","nucleoplasm","nucleolus","vacuole","mito intermembrane space","mito matrix","lipid droplet","Golgi / endosome","peroxisome","nuclear pore channel"],G={0:"#7EA7D3",1:"#B993C9",2:"#D39BC0",3:"#A3C48C",4:"#9FC7C0",5:"#D9946B",6:"#E8B48F",7:"#B7C77E",8:"#C5B0E0",9:"#9A8352",10:"#7EA7D3",11:"#8FA3D0",12:"#7EA7D3",13:"#7EA7D3",14:"#5FAE75",15:"#7DBE86",16:"#6BB596",17:"#5BA8A0",18:"#C27C55",19:"#B5653F",20:"#93B060",21:"#A68BC9",22:"#F4B26A",23:"#3E9C9A",32:"#F3C643",33:"#F0BC3C",34:"#F6D15E",35:"#EAB542",36:"#EFA83A",37:"#E7A93C",38:"#F2CA55",40:"#F3C643",41:"#F8E3A0",42:"#E7A93C",48:"#B0437A",49:"#E07636",50:"#CDA866",51:"#C3C27A",52:"#C9558F",53:"#7CC6C2",56:"#5FAE75",57:"#7DBE86",58:"#6BB596",59:"#5BA8A0",60:"#C27C55",61:"#B5653F",62:"#93B060",63:"#A68BC9",64:"#E9C46A",65:"#9A8352",66:"#A9C3E3",67:"#A3C48C",68:"#B993C9",69:"#9FC7C0",70:"#E8B48F",71:"#D9946B",72:"#B7C77E",73:"#C5B0E0",74:"#E9C46A",75:"#9A8352",80:"#F6F1E7",81:"#F1E9D3",82:"#EEE3CC",83:"#E4ECF5",84:"#E8F0DF",85:"#E8F0DF",86:"#F1E7F1",87:"#EBD9EA",88:"#E3F0EE",89:"#F6E6DA",90:"#F2DCCB",91:"#F7EDC4",92:"#EEF1DC",93:"#EEE8F6",94:"#ECE8F2"};function Sn(){const i=new Float32Array(384);for(const[e,t]of Object.entries(G)){const n=parseInt(t.slice(1),16),s=Number(e)*4;i[s]=(n>>16&255)/255,i[s+1]=(n>>8&255)/255,i[s+2]=(n&255)/255,i[s+3]=1}return i}const Vt={ribosome_80S:{c0:"#EE9BC0",c1:"#C9558F",cls:"ribosome"},tRNA:{c0:"#E27FA8",c1:"#E27FA8",cls:"tRNA"},nucleosome:{c0:"#F4B26A",c1:"#E07636",cls:"nucleosome"},NPC:{c0:"#3E9C9A",cls:"NPC"},PolI:{c0:"#9A72B8",c1:"#9A72B8",cls:"polymerase"},PolII:{c0:"#9A72B8",c1:"#9A72B8",cls:"polymerase"},PolIII:{c0:"#9A72B8",c1:"#9A72B8",cls:"polymerase"},FAS:{c0:"#8B88C9",cls:"complex"},proteasome_26S:{c0:"#8B88C9",cls:"complex"},CCT:{c0:"#8B88C9",cls:"complex"},F_actin:{c0:"#E0893A",c1:"#E0893A",cls:"actin"},Arp23:{c0:"#D9534F",c1:"#D9534F",cls:"actin"},TPM1_coil:{c0:"#2FA4A0",c1:"#2FA4A0",cls:"actin"},CP_cap:{c0:"#8E6CB0",c1:"#8E6CB0",cls:"actin"},MT_tubulin:{c0:"#74C476",c1:"#74C476",cls:"microtubule"},gTuRC_nuc:{c0:"#2E8B57",c1:"#2E8B57",cls:"microtubule"},gTuRC_cyto:{c0:"#2E8B57",c1:"#2E8B57",cls:"microtubule"},ribosome_80S_Ssh1:{c0:"#EE9BC0",c1:"#C9558F",cls:"ribosome"},ribosome_80S_Sec61:{c0:"#EE9BC0",c1:"#C9558F",cls:"ribosome"}};function Wt(i){const e=i*137.508%360,t=.62,n=[.58,.44,.7][i%3],s=o=>{const r=(o+e/30)%12;return n-t*Math.min(n,1-n)*Math.max(-1,Math.min(r-3,9-r,1))};return"#"+[s(0),s(8),s(4)].map(o=>Math.round(o*255).toString(16).padStart(2,"0")).join("")}async function Ln(i){i("Loading structure index",.02);const e=await Ce("library/types.json"),t=e.types,n=e.bead_q;i("Loading residue coordinates (2.5 M beads)",.06);const s=await Me("library/beads.bin",(k,E)=>i("Loading residue coordinates (2.5 M beads)",.06+.5*(E?k/E:0))),o=new Float32Array(await Me("library/spheretree.bin"));let r=new Uint32Array(4),a={},d=0;try{const k=await Ce("library/modes.json");r=new Uint32Array(await Me("library/modes.bin")),a=k.types,d=k.cap??1.5}catch{}let l=[],c={};try{const k=await Ce("library/organisation.json");l=k.modules,c=k.types}catch{}i("Building level-of-detail clusters",.6);const u=s.byteLength/8,m=new Int16Array(s),f=new Uint8Array(s),T=new Float32Array(u*4),w=new Uint32Array(T.buffer);for(let k=0;k<u;k++)T[4*k]=m[4*k]*n,T[4*k+1]=m[4*k+1]*n,T[4*k+2]=m[4*k+2]*n,w[4*k+3]=f[8*k+6]|f[8*k+7]<<8;let M=0;for(const k of t)M+=Math.ceil(k.nb/64);const p=new ArrayBuffer(M*32),x=new Float32Array(p),h=new Uint32Array(p),g=new ArrayBuffer(t.length*64),b=new Uint32Array(g),A=new Float32Array(g);let C=0;t.forEach((k,E)=>{const Y=C;for(let me=0;me<k.nb;me+=64){const ve=Math.min(64,k.nb-me),Fe=k.bo+me;let Be=1e9,Re=1e9,Ne=1e9,De=-1e9,ee=-1e9,ae=-1e9;for(let H=0;H<ve;H++){const te=4*(Fe+H),We=T[te],je=T[te+1],He=T[te+2];We<Be&&(Be=We),We>De&&(De=We),je<Re&&(Re=je),je>ee&&(ee=je),He<Ne&&(Ne=He),He>ae&&(ae=He)}const q=(Be+De)/2,Z=(Re+ee)/2,ne=(Ne+ae)/2;let le=0;for(let H=0;H<ve;H++){const te=4*(Fe+H);le=Math.max(le,Math.hypot(T[te]-q,T[te+1]-Z,T[te+2]-ne))}x[8*C]=q,x[8*C+1]=Z,x[8*C+2]=ne,x[8*C+3]=le+.85+(a[E]?.[3]??0),h[8*C+4]=Fe,h[8*C+5]=ve,C++}const $=Vt[k.id],W=k.src==="PDB";b[16*E+0]=k.bo,b[16*E+1]=k.nb,b[16*E+2]=Y,b[16*E+3]=C-Y,b[16*E+4]=k.so,b[16*E+5]=k.sn,A[16*E+6]=k.br,A[16*E+7]=Math.sqrt(5/3)*k.rg,b[16*E+8]=$?at($.c0):0,b[16*E+9]=$&&$.c1?at($.c1):0,b[16*E+10]=((W?1:0)|Math.min(k.copies,16777215)<<8)>>>0,A[16*E+11]=W?255:k.plddt_mean??255;const z=a[E];z&&(b[16*E+12]=z[0],b[16*E+13]=z[1],b[16*E+14]=z[2]),b[16*E+15]=c[E]!==void 0?at(Wt(c[E])):0});const F=new Map;t.forEach((k,E)=>F.set(k.id,E));const O=new Int32Array(t.length).fill(-1);for(const[k,E]of Object.entries(c))O[Number(k)]=E;return{types:t,beads:T,st:o,clusters:p,typeInfo:g,modes:r,modeCap:d,modules:l,moduleOf:O,nClusters:M,byId:F}}const Ge=26,gt=(i,e,t)=>((i+13)*Ge+(e+13))*Ge+(t+13);async function In(){const i=await Ce("tiles/index.json"),e=i.tile_columns,t=i.tiles,n=t.length,s=a=>e.indexOf(a),o=(a,d)=>{const l=new a(n),c=s(d);for(let u=0;u<n;u++)l[u]=t[u][c];return l},r={n,ix:o(Int32Array,"ix"),iy:o(Int32Array,"iy"),iz:o(Int32Array,"iz"),nInst:o(Uint32Array,"n_inst"),nGhost:o(Uint32Array,"n_ghost"),nMr:o(Uint32Array,"n_mrna_beads"),nLk:o(Uint32Array,"n_link_beads"),nVert:o(Uint32Array,"n_vert"),nTri:o(Uint32Array,"n_tri"),gz:o(Uint32Array,"gz_bytes"),labelMask:o(Uint32Array,"label_mask"),volGlucan:o(Float32Array,"vol_glucan_nm3"),volMannan:o(Float32Array,"vol_mannan_nm3"),key:new Map,wall:i.wall,regionNames:i.region_names};for(let a=0;a<n;a++)r.key.set(gt(r.ix[a],r.iy[a],r.iz[a]),a);return r}const Fn=1073741824,B=GPUShaderStage,I=GPUBufferUsage;class Bn{device;ctx;format;canvas;limits={};adapterInfo="";instPages=0;lipPages=0;freePages=0;capL0=0;capL1=0;capL2=0;capLip=0;capFree=0;maxSlots=4096;uni;inst;pageOrg;typesBuf;typeFlags;palette;beads;st;clusters;modes;modeCap=0;counters;list01;listL2;args;cullParams;caps;lipPool;lipPagesBuf;lipList;lipParams;freePool;freePagesBuf;freeList;freeParams;meshUni;tileState;fillBufs=[];pickBuf;pickParams;pickOut;statsBuf;pipes={};bg={};bgl={};colorTex;idTex;depthTex;width=0;height=0;usedInstPages=0;usedLipPages=0;usedFreePages=0;pickPending=null;pickBusy=!1;statsBusy=!1;stats={l0:0,l1:0,l2:0,lip:0,fre:0,ninst:0};lastStatsT=0;organelleMeshes=[];onLost=()=>{};constructor(e){this.canvas=e}async init(e,t){if(!("gpu"in navigator))throw new Error("WebGPU is not available in this browser.");const n=await navigator.gpu.requestAdapter({powerPreference:"high-performance"});if(!n)throw new Error("No WebGPU adapter found (GPU blocked or unsupported).");try{const u=n.info||{};this.adapterInfo=[u.vendor,u.architecture,u.device,u.description].filter(Boolean).join(" ")}catch{}const s=n.limits,o={maxStorageBufferBindingSize:Math.min(s.maxStorageBufferBindingSize,2**31-4),maxBufferSize:Math.min(s.maxBufferSize,2**31-4),maxComputeWorkgroupsPerDimension:s.maxComputeWorkgroupsPerDimension};this.device=await n.requestDevice({requiredLimits:o}),this.device.lost.then(u=>this.onLost(`GPU device lost: ${u.message}`)),this.device.addEventListener("uncapturederror",u=>console.error("WebGPU error:",u.error?.message)),this.limits={...o};const r=this.device.limits.maxStorageBufferBindingSize,a=u=>Math.min(u,Math.floor(r/16));this.instPages=Math.floor(a(8*1024*1024)/X),this.lipPages=Math.floor(a(8*1024*1024)/Q),this.freePages=Math.floor(a(10*1024*1024)/Q),this.capL0=Math.min(4*1024*1024,Math.floor(r/8)-2*1024*1024),this.capL1=2*1024*1024,this.capL2=Math.min(8*1024*1024,Math.floor(r/4)),this.capLip=this.lipPages*Q,this.capFree=this.freePages*Q,this.ctx=this.canvas.getContext("webgpu"),this.format=navigator.gpu.getPreferredCanvasFormat(),this.ctx.configure({device:this.device,format:this.format,alphaMode:"opaque",usage:GPUTextureUsage.RENDER_ATTACHMENT|GPUTextureUsage.COPY_SRC});const d=this.device,l=(u,m,f)=>d.createBuffer({size:Math.max(16,Math.ceil(u/16)*16),usage:m,label:f}),c=I.STORAGE|I.COPY_DST;this.uni=l(512,I.UNIFORM|I.COPY_DST,"uniforms"),this.inst=l(this.instPages*X*16,c,"instances"),this.pageOrg=l(this.instPages*16,c,"pageOrg"),this.typesBuf=l(e.typeInfo.byteLength,c,"types"),this.typeFlags=l(e.types.length*4,c,"typeFlags"),this.palette=l(1536,c,"palette"),this.beads=l(e.beads.byteLength,c,"beads"),this.modes=l(Math.max(16,e.modes.byteLength),c,"modes"),this.modeCap=e.modeCap,this.st=l(e.st.byteLength,c,"spheretree"),this.clusters=l(e.clusters.byteLength,c,"clusters"),d.queue.writeBuffer(this.typesBuf,0,e.typeInfo),d.queue.writeBuffer(this.palette,0,Sn()),d.queue.writeBuffer(this.beads,0,e.beads),d.queue.writeBuffer(this.modes,0,e.modes),d.queue.writeBuffer(this.st,0,e.st),d.queue.writeBuffer(this.clusters,0,e.clusters),this.counters=l(32,I.STORAGE|I.COPY_DST|I.COPY_SRC,"counters"),this.list01=l((this.capL0+this.capL1)*8,I.STORAGE,"list01"),this.listL2=l(this.capL2*4,I.STORAGE,"listL2"),this.args=l(80,I.STORAGE|I.INDIRECT,"args"),this.cullParams=l(32,I.UNIFORM|I.COPY_DST,"cullParams"),this.caps=l(32,I.UNIFORM|I.COPY_DST,"caps"),d.queue.writeBuffer(this.caps,0,new Uint32Array([this.capL0,this.capL1,this.capL2,this.capLip,this.capFree,0,0,0])),this.lipPool=l(this.lipPages*Q*16,c,"lipidPool"),this.lipPagesBuf=l(this.lipPages*16,c,"lipidPages"),this.lipList=l(this.capLip*4,I.STORAGE,"lipidList"),this.lipParams=l(16,I.UNIFORM|I.COPY_DST,"lipParams"),this.freePool=l(this.freePages*Q*16,c,"freePool"),this.freePagesBuf=l(this.freePages*16,c,"freePages"),this.freeList=l(this.capFree*4,I.STORAGE,"freeList"),this.freeParams=l(16,I.UNIFORM|I.COPY_DST,"freeParams"),this.meshUni=l(256*(this.maxSlots+16),I.UNIFORM|I.COPY_DST,"meshUniforms"),this.tileState=l(676*26*4,c,"tileState"),this.pickBuf=d.createBuffer({size:16,usage:I.MAP_READ|I.COPY_DST}),this.pickParams=l(16,I.UNIFORM|I.COPY_DST,"pickParams"),this.pickOut=l(16,I.STORAGE|I.COPY_SRC,"pickOut"),this.statsBuf=d.createBuffer({size:32,usage:I.MAP_READ|I.COPY_DST}),this.initFill(t),await this.initPipelines(),this.resize()}initFill(e){const t=this.device,n=e.PT.length,s=new ArrayBuffer(n*48),o=new Float32Array(s),r=new Uint32Array(s);for(let p=0;p<n;p++){for(let x=0;x<8;x++)o[12*p+x]=e.PQ[8*p+x];r[12*p+8]=e.PT[p],r[12*p+9]=e.PC[p]}const a=e.pores.length/6,d=new Float32Array(Math.max(1,a)*8);for(let p=0;p<a;p++)for(let x=0;x<3;x++)d[8*p+x]=e.pores[6*p+x],d[8*p+4+x]=e.pores[6*p+3+x];const l=e.grid,c=l.prim_idx.length,u=new Uint32Array(l.prim_off.length+l.pore_off.length);u.set(l.prim_off,0);for(let p=0;p<l.pore_off.length;p++)u[l.prim_off.length+p]=l.pore_off[p]+c;const m=new Uint32Array(c+l.pore_idx.length+1);m.set(l.prim_idx,0),m.set(l.pore_idx,c);const f=e.G,T=e.consts,w=new Float32Array([f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],0,T.PORE_R,T.PERI_W,T.IMS_W,a,l.n,l.ext,0,0,f[11],f[12],f[13],0,f[14],f[15],f[16],0,f[17],f[18],f[19],0,f[20],f[21],f[22],0]),M=(p,x)=>{const h=t.createBuffer({size:Math.ceil(p.byteLength/16)*16+16,usage:x|I.COPY_DST});return t.queue.writeBuffer(h,0,p),h};this.fillBufs=[M(s,I.STORAGE),M(d,I.STORAGE),M(u,I.STORAGE),M(m,I.STORAGE),M(w,I.UNIFORM),M(e.F??new Float32Array(4),I.STORAGE)]}async initPipelines(){const e=this.device,t=(q,Z)=>({binding:q,visibility:Z,buffer:{type:"read-only-storage"}}),n=q=>({binding:q,visibility:B.COMPUTE,buffer:{type:"storage"}}),s=(q,Z,ne=!1)=>({binding:q,visibility:Z,buffer:{type:"uniform",hasDynamicOffset:ne}}),o=B.COMPUTE|B.VERTEX|B.FRAGMENT,r=B.COMPUTE|B.VERTEX,a=B.VERTEX|B.FRAGMENT,d=q=>e.createBindGroupLayout({entries:q});this.bgl.g0=d([s(0,o),t(1,r),t(2,r),t(3,r),t(4,r),t(5,a)]),this.bgl.cull=d([n(0),n(1),n(2),s(3,B.COMPUTE),t(4,B.COMPUTE)]),this.bgl.cullS=d([n(0),t(1,B.COMPUTE),t(2,B.COMPUTE),n(3),s(4,B.COMPUTE)]),this.bgl.fin=d([t(0,B.COMPUTE),n(1),s(2,B.COMPUTE)]),this.bgl.l0=d([t(0,B.VERTEX),t(1,B.VERTEX),t(2,B.VERTEX)]),this.bgl.l2=d([t(0,B.VERTEX)]),this.bgl.sph=d([t(0,B.VERTEX),t(1,B.VERTEX),t(2,B.VERTEX)]),this.bgl.mesh=d([s(0,a,!0)]),this.bgl.fill=d([t(0,B.FRAGMENT),t(1,B.FRAGMENT),t(2,B.FRAGMENT),t(3,B.FRAGMENT),s(4,B.FRAGMENT),t(5,B.FRAGMENT),t(6,B.FRAGMENT)]),this.bgl.comp=d([{binding:0,visibility:B.FRAGMENT,texture:{sampleType:"unfilterable-float"}},{binding:1,visibility:B.FRAGMENT,texture:{sampleType:"uint"}},{binding:2,visibility:B.FRAGMENT,texture:{sampleType:"depth"}}]),this.bgl.pick=d([{binding:0,visibility:B.COMPUTE,texture:{sampleType:"uint"}},{binding:1,visibility:B.COMPUTE,texture:{sampleType:"depth"}},s(2,B.COMPUTE),n(3)]);const l=(...q)=>e.createPipelineLayout({bindGroupLayouts:q}),c=async(q,Z)=>{const ne=e.createShaderModule({code:q,label:Z}),le=await ne.getCompilationInfo(),H=le.messages.filter(te=>te.type==="error");for(const te of le.messages)(te.type==="error"?console.error:console.warn)(`[${Z}] ${te.lineNum}:${te.linePos} ${te.message}`);if(H.length)throw new Error(`WGSL compile error in ${Z}: ${H[0].lineNum}:${H[0].linePos} ${H[0].message}`);return ne},u=await c(un,"cull"),m=await c(hn,"cullSpheres"),f=await c(pn,"finalize"),T=await c(fn,"drawL0"),w=await c(gn,"drawL1"),M=await c(bn,"drawL2"),p=await c(yn,"drawSpheres"),x=await c(wn,"mesh"),h=await c(vn,"fill"),g=await c(Mn,"composite"),b=await c(xn,"pick"),A=(q,Z,ne)=>e.createComputePipelineAsync({layout:ne,compute:{module:q,entryPoint:Z}}),C=[{format:"rgba8unorm"},{format:"r32uint"}],F=(q="greater")=>({format:"depth32float",depthWriteEnabled:!0,depthCompare:q}),O=(q,Z,ne,le,H={})=>e.createRenderPipelineAsync({layout:le,vertex:{module:q,entryPoint:Z,...H.vertex?{buffers:H.vertex.buffers}:{}},fragment:{module:q,entryPoint:ne,targets:C},primitive:{topology:"triangle-list",cullMode:"none"},depthStencil:H.depthStencil??F()}),k=this.bgl.g0,[E,Y,$,W,z,me,ve,Fe,Be,Re,Ne,De]=await Promise.all([A(u,"cullInstances",l(k,this.bgl.cull)),A(m,"cullSpheres",l(k,this.bgl.cullS)),A(f,"main",l(this.bgl.fin)),O(T,"vsL0","fsMain",l(k,this.bgl.l0)),O(w,"vsL1","fsMain",l(k,this.bgl.l0)),O(M,"vsL2","fsMain",l(k,this.bgl.l2)),O(p,"vsLipid","fsMain",l(k,this.bgl.sph)),O(p,"vsFree","fsMain",l(k,this.bgl.sph)),O(x,"vsMesh","fsMesh",l(k,this.bgl.mesh),{vertex:{buffers:[{arrayStride:32,attributes:[{shaderLocation:0,offset:0,format:"float32x4"},{shaderLocation:1,offset:16,format:"float32x4"}]}]}}),O(h,"vsFill","fsFill",l(k,this.bgl.fill),{depthStencil:F("always")}),e.createRenderPipelineAsync({layout:l(k,this.bgl.comp),vertex:{module:g,entryPoint:"vsFull"},fragment:{module:g,entryPoint:"fsComposite",targets:[{format:this.format}]},primitive:{topology:"triangle-list"}}),A(b,"main",l(this.bgl.pick))]);Object.assign(this.pipes,{cull:E,cullS:Y,fin:$,l0:W,l1:z,l2:me,lip:ve,fre:Fe,mesh:Be,fill:Re,comp:Ne,pick:De});const ee=(q,Z,ne=[])=>e.createBindGroup({layout:q,entries:Z.map((le,H)=>({binding:H,resource:{buffer:le,size:ne[H]}}))});this.bg.g0=ee(k,[this.uni,this.inst,this.pageOrg,this.typesBuf,this.typeFlags,this.palette]),this.bg.cull=ee(this.bgl.cull,[this.counters,this.list01,this.listL2,this.cullParams,this.clusters]),this.bg.cullLip=ee(this.bgl.cullS,[this.counters,this.lipPool,this.lipPagesBuf,this.lipList,this.lipParams]),this.bg.cullFree=ee(this.bgl.cullS,[this.counters,this.freePool,this.freePagesBuf,this.freeList,this.freeParams]),this.bg.fin=ee(this.bgl.fin,[this.counters,this.args,this.caps]),this.bg.l0=ee(this.bgl.l0,[this.beads,this.list01,this.modes]),this.bg.l1=ee(this.bgl.l0,[this.st,this.list01,this.modes]),this.bg.l2=ee(this.bgl.l2,[this.listL2]),this.bg.lip=ee(this.bgl.sph,[this.lipPool,this.lipPagesBuf,this.lipList]),this.bg.fre=ee(this.bgl.sph,[this.freePool,this.freePagesBuf,this.freeList]),this.bg.mesh=e.createBindGroup({layout:this.bgl.mesh,entries:[{binding:0,resource:{buffer:this.meshUni,size:16}}]});const ae=this.fillBufs;this.bg.fill=ee(this.bgl.fill,[ae[0],ae[1],ae[2],ae[3],ae[4],this.tileState,ae[5]])}resize(e=1){const t=Math.min(window.devicePixelRatio||1,2)*e,n=Math.max(1,Math.floor(this.canvas.clientWidth*t)),s=Math.max(1,Math.floor(this.canvas.clientHeight*t));if(n===this.width&&s===this.height&&this.colorTex)return;this.canvas.width=n,this.canvas.height=s,this.width=n,this.height=s,this.colorTex?.destroy(),this.idTex?.destroy(),this.depthTex?.destroy();const o=GPUTextureUsage;this.colorTex=this.device.createTexture({size:[n,s],format:"rgba8unorm",usage:o.RENDER_ATTACHMENT|o.TEXTURE_BINDING}),this.idTex=this.device.createTexture({size:[n,s],format:"r32uint",usage:o.RENDER_ATTACHMENT|o.TEXTURE_BINDING|o.COPY_SRC}),this.depthTex=this.device.createTexture({size:[n,s],format:"depth32float",usage:o.RENDER_ATTACHMENT|o.TEXTURE_BINDING|o.COPY_SRC}),this.bg.comp=this.device.createBindGroup({layout:this.bgl.comp,entries:[{binding:0,resource:this.colorTex.createView()},{binding:1,resource:this.idTex.createView()},{binding:2,resource:this.depthTex.createView()}]}),this.bg.pick=this.device.createBindGroup({layout:this.bgl.pick,entries:[{binding:0,resource:this.idTex.createView()},{binding:1,resource:this.depthTex.createView()},{binding:2,resource:{buffer:this.pickParams}},{binding:3,resource:{buffer:this.pickOut}}]})}writeUniforms(e,t,n){const s=this.width/this.height,o=e.viewRot(),r=e.proj(s),a=Pn(r,o),d=An(a),l=e.eye,c=new ArrayBuffer(512),u=new Float32Array(c),m=new Uint32Array(c),f=t.memColors&&!t.organelle;u.set(r,0),u.set(o,16),u.set(d,32);const T=e.pxPerNmAtUnit(this.height);u.set([l[0],l[1],l[2],T],48);const w=[t.plane[0],t.plane[1],t.plane[2]];u.set([w[0],w[1],w[2],t.plane[3]],52);const M=[o[0]*w[0]+o[4]*w[1]+o[8]*w[2],o[1]*w[0]+o[5]*w[1]+o[9]*w[2],o[2]*w[0]+o[6]*w[1]+o[10]*w[2]];u.set([M[0],M[1],M[2],V(w,l)-t.plane[3]],56),u.set([t.cutOn&&!f?1:0,t.slab,Math.min(400,Math.max(20,t.slab)),performance.now()/1e3],60),u.set([t.l0px,t.l1px,t.minpx,.6],64),u.set([this.width,this.height,1/this.width,1/this.height],68),u.set(e.frustum(s),72),m.set([t.colorMode,t.regionMask>>>0,t.freeHide,t.highlight?1:0],96),m.set([t.selId>>>0,0,t.shaded?1:0,t.organelle?1:0],100),u.set([360/65535,.8,t.outline,Math.min(window.devicePixelRatio||1,2)*t.renderScale],104),u.set([n[0],n[1],n[2],e.near()],108),m.set([this.capL0,Math.round(t.motion*1e3),Math.round(this.modeCap*1e3),f?1:0],112),this.device.queue.writeBuffer(this.uni,0,c,0,464)}frame(e,t,n,s,o){const r=this.device;this.writeUniforms(e,t,n);const a=r.createCommandEncoder();a.clearBuffer(this.counters);const d=t.memColors&&!t.organelle,l=!t.organelle&&!d;if(l){const w=a.beginComputePass();w.setBindGroup(0,this.bg.g0);const M=(p,x,h)=>{const g=Math.ceil(p/64),b=Math.min(g,65535),A=Math.ceil(g/b);return r.queue.writeBuffer(x,0,h(b*64)),[b,A]};if(this.usedInstPages>0){const[p,x]=M(this.usedInstPages*X,this.cullParams,h=>new Uint32Array([this.usedInstPages,h,this.capL0,this.capL1,this.capL2,this.capLip,this.capFree,0]).buffer);w.setPipeline(this.pipes.cull),w.setBindGroup(1,this.bg.cull),w.dispatchWorkgroups(p,x)}if(this.usedLipPages>0&&t.showLipids){const[p,x]=M(this.usedLipPages*Q,this.lipParams,h=>new Uint32Array([this.usedLipPages,h,this.capLip,1]).buffer);w.setPipeline(this.pipes.cullS),w.setBindGroup(1,this.bg.cullLip),w.dispatchWorkgroups(p,x)}if(this.usedFreePages>0){const[p,x]=M(this.usedFreePages*Q,this.freeParams,h=>new Uint32Array([this.usedFreePages,h,this.capFree,0]).buffer);w.setPipeline(this.pipes.cullS),w.setBindGroup(1,this.bg.cullFree),w.dispatchWorkgroups(p,x)}w.setPipeline(this.pipes.fin),w.setBindGroup(0,this.bg.fin),w.dispatchWorkgroups(1),w.end()}const c=a.beginRenderPass({colorAttachments:[{view:this.colorTex.createView(),loadOp:"clear",storeOp:"store",clearValue:{r:.965,g:.945,b:.906,a:1}},{view:this.idTex.createView(),loadOp:"clear",storeOp:"store",clearValue:{r:0,g:0,b:0,a:0}}],depthStencilAttachment:{view:this.depthTex.createView(),depthLoadOp:"clear",depthStoreOp:"store",depthClearValue:0}});if(c.setBindGroup(0,this.bg.g0),c.setPipeline(this.pipes.fill),c.setBindGroup(1,this.bg.fill),c.draw(3),t.organelle||d){c.setPipeline(this.pipes.mesh);const w=this.maxSlots;r.queue.writeBuffer(this.meshUni,w*256,new Uint32Array([0,t.organelle?1:3,0,0]));for(const M of this.organelleMeshes)c.setBindGroup(1,this.bg.mesh,[w*256]),c.setVertexBuffer(0,M.vb),c.setIndexBuffer(M.ib,"uint32"),c.drawIndexed(M.nIdx,1)}else if(t.showMeshes&&s.length){c.setPipeline(this.pipes.mesh);for(const w of s)c.setBindGroup(1,this.bg.mesh,[w.slot*256]),c.setVertexBuffer(0,w.vb),c.setIndexBuffer(w.ib,"uint32"),c.drawIndexed(w.nIdx,2)}l&&(c.setPipeline(this.pipes.l2),c.setBindGroup(1,this.bg.l2),c.drawIndirect(this.args,32),c.setPipeline(this.pipes.l1),c.setBindGroup(1,this.bg.l1),c.drawIndirect(this.args,16),c.setPipeline(this.pipes.l0),c.setBindGroup(1,this.bg.l0),c.drawIndirect(this.args,0),t.showLipids&&(c.setPipeline(this.pipes.lip),c.setBindGroup(1,this.bg.lip),c.drawIndirect(this.args,48)),c.setPipeline(this.pipes.fre),c.setBindGroup(1,this.bg.fre),c.drawIndirect(this.args,64)),c.end();const u=this.ctx.getCurrentTexture(),m=a.beginRenderPass({colorAttachments:[{view:u.createView(),loadOp:"clear",storeOp:"store",clearValue:{r:1,g:1,b:1,a:1}}]});m.setPipeline(this.pipes.comp),m.setBindGroup(0,this.bg.g0),m.setBindGroup(1,this.bg.comp),m.draw(3),m.end();let f=null;if(this.pickPending&&!this.pickBusy){f=this.pickPending,this.pickPending=null;const w=Math.max(0,Math.min(this.width-1,Math.floor(f.x))),M=Math.max(0,Math.min(this.height-1,Math.floor(f.y)));r.queue.writeBuffer(this.pickParams,0,new Uint32Array([w,M,0,0]));const p=a.beginComputePass();p.setPipeline(this.pipes.pick),p.setBindGroup(0,this.bg.pick),p.dispatchWorkgroups(1),p.end(),a.copyBufferToBuffer(this.pickOut,0,this.pickBuf,0,16)}let T=!1;if(o&&!this.statsBusy&&l&&(a.copyBufferToBuffer(this.counters,0,this.statsBuf,0,32),T=!0),r.queue.submit([a.finish()]),f){this.pickBusy=!0;const w=f;this.pickBuf.mapAsync(GPUMapMode.READ).then(()=>{const M=this.pickBuf.getMappedRange(0,8),p=new Uint32Array(M)[0],x=new Float32Array(M)[1];this.pickBuf.unmap(),this.pickBusy=!1,w.resolve({id:p>>>0,depth:x})}).catch(()=>{this.pickBusy=!1})}T&&(this.statsBusy=!0,this.statsBuf.mapAsync(GPUMapMode.READ).then(()=>{const w=new Uint32Array(this.statsBuf.getMappedRange());this.stats={l0:w[0],l1:w[1],l2:w[2],lip:w[3],fre:w[4],ninst:w[5]},this.statsBuf.unmap(),this.statsBusy=!1}).catch(()=>{this.statsBusy=!1}))}pick(e,t){return new Promise(n=>{this.pickPending={x:e,y:t,resolve:n}})}setMeshSlot(e,t){const n=new ArrayBuffer(16),s=new Uint32Array(n),o=new Float32Array(n);s[0]=e,s[1]=t?2:0,o[2]=t?-1:-1.9,o[3]=t?1:1.9,this.device.queue.writeBuffer(this.meshUni,e*256,n)}uploadOrganelleMeshes(e){const t=this.device,n=new Uint32Array(e,0,1)[0],s=new Uint32Array(e,4,n*3);let o=4+n*12;for(const r of this.organelleMeshes)r.vb.destroy(),r.ib.destroy();this.organelleMeshes=[];for(let r=0;r<n;r++){const a=s[3*r],d=s[3*r+1],l=s[3*r+2],c=new Float32Array(e.slice(o,o+d*12));o+=d*12;const u=new Uint32Array(e.slice(o,o+l*12));o+=l*12;const m=new Float32Array(d*3);for(let p=0;p<l;p++){const x=u[3*p],h=u[3*p+1],g=u[3*p+2],b=c[3*h]-c[3*x],A=c[3*h+1]-c[3*x+1],C=c[3*h+2]-c[3*x+2],F=c[3*g]-c[3*x],O=c[3*g+1]-c[3*x+1],k=c[3*g+2]-c[3*x+2],E=A*k-C*O,Y=C*F-b*k,$=b*O-A*F;for(const W of[x,h,g])m[3*W]+=E,m[3*W+1]+=Y,m[3*W+2]+=$}const f=new Float32Array(d*8),T=new Uint32Array(f.buffer);for(let p=0;p<d;p++){const x=Math.hypot(m[3*p],m[3*p+1],m[3*p+2])||1;f[8*p]=c[3*p],f[8*p+1]=c[3*p+1],f[8*p+2]=c[3*p+2],T[8*p+3]=a,f[8*p+4]=m[3*p]/x,f[8*p+5]=m[3*p+1]/x,f[8*p+6]=m[3*p+2]/x}const w=t.createBuffer({size:f.byteLength,usage:I.VERTEX|I.COPY_DST});t.queue.writeBuffer(w,0,f);const M=t.createBuffer({size:Math.ceil(u.byteLength/4)*4,usage:I.INDEX|I.COPY_DST});t.queue.writeBuffer(M,0,u),this.organelleMeshes.push({vb:w,ib:M,nIdx:l*3,mem:a})}}gpuMemoryMB(){return[this.inst,this.pageOrg,this.beads,this.st,this.clusters,this.list01,this.listL2,this.lipPool,this.lipList,this.freePool,this.freeList,this.meshUni].reduce((t,n)=>t+n.size,0)/1048576+this.width*this.height*12/1048576}}class Rn{target=[0,0,0];dist=9e3;fwd=[0,0,-1];up=[0,1,0];fovy=35*Math.PI/180;minDist=1.5;maxDist=14e3;anim=null;get right(){return J(ge(this.fwd,this.up))}get eye(){return ue(this.target,L(this.fwd,this.dist))}orthonormalize(){this.fwd=J(this.fwd);const e=J(ge(this.fwd,this.up));this.up=J(ge(e,this.fwd))}orbit(e,t){const n=this.right;this.fwd=Ue(this.fwd,this.up,-e),this.fwd=Ue(this.fwd,n,-t),this.up=Ue(this.up,n,-t),this.orthonormalize()}roll(e){this.up=Ue(this.up,this.fwd,e),this.orthonormalize()}pan(e,t,n){const s=1/n;this.target=_(this.target,_(L(this.right,-e*s),L(this.up,t*s)))}zoom(e){this.dist=Math.min(this.maxDist,Math.max(this.minDist,this.dist*e))}fly(e){this.target=_(this.target,e)}pxPerNmAtUnit(e){return e/(2*Math.tan(this.fovy/2))}distForWidth(e,t){return e/(2*Math.tan(this.fovy/2)*Math.min(t,1))}near(){return Math.min(5,Math.max(.02,this.dist*.004))}viewRot(){const e=this.right,t=this.up,n=this.fwd,s=new Float32Array(16);return s[0]=e[0],s[4]=e[1],s[8]=e[2],s[1]=t[0],s[5]=t[1],s[9]=t[2],s[2]=-n[0],s[6]=-n[1],s[10]=-n[2],s[15]=1,s}proj(e){const t=1/Math.tan(this.fovy/2),n=new Float32Array(16);return n[0]=t/e,n[5]=t,n[11]=-1,n[14]=this.near(),n}frustum(e){const t=Math.tan(this.fovy/2),n=t*e,s=this.right,o=this.up,r=this.fwd,a=(l,c=0)=>{const u=J(l);return[u[0],u[1],u[2],c]},d=[...a(_(s,L(r,n))),...a(_(L(s,-1),L(r,n))),...a(_(o,L(r,t))),...a(_(L(o,-1),L(r,t))),r[0],r[1],r[2],-this.near(),0,0,0,1];return new Float32Array(d)}ray(e,t,n){const s=Math.tan(this.fovy/2),o=s*n,r=(e*2-1)*o,a=(1-t*2)*s;return J(_(this.fwd,_(L(this.right,r),L(this.up,a))))}state(){return{target:[...this.target],dist:this.dist,fwd:[...this.fwd],up:[...this.up]}}set(e){this.target=[...e.target],this.dist=e.dist,this.fwd=[...e.fwd],this.up=[...e.up],this.orthonormalize()}flyTo(e,t=1800,n,s){this.anim={t0:performance.now(),dur:t,a:this.state(),b:e,onDone:n,onStep:s}}step(e){if(!this.anim)return!1;const t=this.anim;let n=Math.min(1,(e-t.t0)/t.dur);const s=kn(n),o=Math.log(t.a.dist),r=Math.log(t.b.dist),a=Math.hypot(...ue(t.b.target,t.a.target)),d=Math.log(Math.max(t.a.dist,t.b.dist,a*1.2)),l=Math.max(0,d-Math.max(o,r)),c=o+(r-o)*s+l*Math.sin(Math.PI*s);this.target=ft(t.a.target,t.b.target,s),this.dist=Math.exp(c),this.fwd=Nn(t.a.fwd,t.b.fwd,s);let u=t.a.up;return this.up=J(ft(u,t.b.up,s)),Math.abs(V(this.up,this.fwd))>.98&&(this.up=J(ge(ge(this.fwd,t.b.up),this.fwd))),this.orthonormalize(),t.onStep?.(s),n>=1&&(this.anim=null,t.onDone?.()),!0}}function Nn(i,e,t){const n=Math.max(-1,Math.min(1,V(i,e))),s=Math.acos(n);if(s<1e-4)return J(ft(i,e,t));if(s>Math.PI-.001){const r=J(ge(i,Math.abs(i[0])<.9?[1,0,0]:[0,1,0]));return Ue(i,r,Math.PI*t)}const o=Math.sin(s);return J(_(L(i,Math.sin((1-t)*s)/o),L(e,Math.sin(t*s)/o)))}const bt=0,Lt=1,It=2;function Ze(i){return i&&i.length?i[0]:0}function yt(i,e,t,n,s){const o=1+12*s,r=i[o+3],a=i[o+7],d=(e-i[o])/r,l=(t-i[o+1])/r,c=(n-i[o+2])/r,u=i[o+4],m=i[o+5],f=i[o+6];if(d<0||l<0||c<0||d>=u*a||l>=m*a||c>=f*a)return i[o+9];const T=Math.floor(d),w=Math.floor(l),M=Math.floor(c),p=Math.floor(T/a),x=Math.floor(w/a),h=Math.floor(M/a),g=i[i[o+10]+(p*m+x)*f+h];if(g===-1)return-i[o+8];if(g===-2)return i[o+9];const b=a+1,A=b*b,C=b,F=d-T,O=l-w,k=c-M,E=i[o+11]+g*b*b*b+((T-p*a)*b+(w-x*a))*b+(M-h*a),Y=i[E]+(i[E+A]-i[E])*F,$=i[E+1]+(i[E+A+1]-i[E+1])*F,W=i[E+C]+(i[E+A+C]-i[E+C])*F,z=i[E+C+1]+(i[E+A+C+1]-i[E+C+1])*F,me=Y+(W-Y)*O,ve=$+(z-$)*O;return me+(ve-me)*k}function At(i,e,t,n,s){return Ze(i)>bt?-yt(i,t,n,s,bt):e[0]-Math.hypot(t,n,s)}function Dn(i,e,t,n,s,o=2){if(Ze(i)<=bt){const u=Math.hypot(t,n,s)||1;return[t/u,n/u,s/u]}const r=(u,m,f)=>At(i,e,u,m,f),a=r(t-o,n,s)-r(t+o,n,s),d=r(t,n-o,s)-r(t,n+o,s),l=r(t,n,s-o)-r(t,n,s+o),c=Math.hypot(a,d,l);return c>1e-9?[a/c,d/c,l/c]:[0,0,0]}async function On(i,e){const t=i.fields;if(!t)return null;const n=await e(t.file),s=t.head.length;let o=1+12*s;for(let d=0;d<s;d++)o+=t.n_index[d]+t.n_data[d];const r=new Float32Array(o);r[0]=s;let a=0;for(let d=0;d<s;d++){r.set(t.head[d],1+12*d);const l=new Int32Array(n,a,t.n_index[d]);a+=4*t.n_index[d],r.set(l,t.head[d][10]);const c=new Int16Array(n,a,t.n_data[d]);a+=2*t.n_data[d],a=Math.ceil(a/4)*4;const u=t.head[d][11];for(let m=0;m<c.length;m++)r[u+m]=c[m]*t.q}return r}const jt='function sn(n){return n&&n.length?n[0]:0}function wn(n,t,e,o,a){const s=1+12*a,h=n[s+3],l=n[s+7],c=(t-n[s])/h,i=(e-n[s+1])/h,d=(o-n[s+2])/h,m=n[s+4],k=n[s+5],f=n[s+6];if(c<0||i<0||d<0||c>=m*l||i>=k*l||d>=f*l)return n[s+9];const p=Math.floor(c),r=Math.floor(i),M=Math.floor(d),g=Math.floor(p/l),x=Math.floor(r/l),U=Math.floor(M/l),C=n[n[s+10]+(g*k+x)*f+U];if(C===-1)return-n[s+8];if(C===-2)return n[s+9];const A=l+1,b=A*A,I=A,u=c-p,B=i-r,w=d-M,y=n[s+11]+C*A*A*A+((p-g*l)*A+(r-x*l))*A+(M-U*l),v=n[y]+(n[y+b]-n[y])*u,z=n[y+1]+(n[y+b+1]-n[y+1])*u,L=n[y+I]+(n[y+b+I]-n[y+I])*u,T=n[y+I+1]+(n[y+b+I+1]-n[y+I+1])*u,H=v+(L-v)*B,S=z+(T-z)*B;return H+(S-H)*w}function lt(n,t,e,o,a){return sn(n)>0?-wn(n,e,o,a,0):t[0]-Math.hypot(e,o,a)}function Qt(n,t,e,o,a,s=2){if(sn(n)<=0){const m=Math.hypot(e,o,a)||1;return[e/m,o/m,a/m]}const h=(m,k,f)=>lt(n,t,m,k,f),l=h(e-s,o,a)-h(e+s,o,a),c=h(e,o-s,a)-h(e,o+s,a),i=h(e,o,a-s)-h(e,o,a+s),d=Math.hypot(l,c,i);return d>1e-9?[l/d,c/d,i/d]:[0,0,0]}const nt=200,zt=80,Lt=360/65535,Rt=.005,mn=.01,Mn=.70710678;let Kt="",Jt={},_t,At,Bt,D,Vt=null,J=[],tt=null,Wt={glucanPerNm3:.018147,mannanPerProtein:1908.6},qt="";const ot=new Map,pn=48;async function Tt(n){if(Jt[n])return tn(Jt[n].slice().buffer);const t=await fetch(Kt+n+(qt?`?v=${encodeURIComponent(qt)}`:""));if(!t.ok)throw new Error(`HTTP ${t.status} ${n}`);return tn(await t.arrayBuffer())}async function tn(n){const t=new Uint8Array(n,0,Math.min(2,n.byteLength));return t[0]===31&&t[1]===139?await new Response(new Blob([n]).stream().pipeThrough(new DecompressionStream("gzip"))).arrayBuffer():n}async function bn(n){Kt=n.base,J=n.G,tt=n.F??null;try{const c=await fetch(Kt+"version.json",{cache:"no-store"});c.ok&&(qt=String((await c.json()).v??""))}catch{qt=""}n.wall&&(Wt=n.wall),n.embedded&&(Jt=n.embedded);const t=JSON.parse(new TextDecoder().decode(await Tt("library/types.json"))),e=t.types,o=e.length,a=t.bead_q;_t=new Uint32Array(o),At=new Uint32Array(o),Bt=new Float32Array(o);for(let c=0;c<o;c++)_t[c]=e[c].bo,At[c]=e[c].nb,Bt[c]=e[c].br;const s=await Tt("library/beads.bin"),h=new Int16Array(s),l=s.byteLength/8;D=new Float32Array(l*3);for(let c=0;c<l;c++)D[3*c]=h[4*c]*a,D[3*c+1]=h[4*c+1]*a,D[3*c+2]=h[4*c+2]*a}const gt=n=>n+3&-4;function Xt(n,t,e,o,a,s,h,l,c,i,d){const m=new Uint32Array(n,t,e);t+=4*e;const k=new Int32Array(n,t,3*e);t+=12*e;const f=new Int8Array(n,t,3*(o-e));t+=3*(o-e),t=gt(t);let p=0,r=i;for(let M=0;M<e;M++){let g=k[3*M],x=k[3*M+1],U=k[3*M+2];const C=m[M];for(let A=0;A<C;A++)A>0&&(g+=f[3*p],x+=f[3*p+1],U+=f[3*p+2],p++),l[4*r]=a[0]+g*s,l[4*r+1]=a[1]+x*s,l[4*r+2]=a[2]+U*s,c[4*r+3]=((d?d[M]:h)|(M&32767)<<4)>>>0,r++}return{off:t,o:r}}function nn(n){const t=new Uint32Array(n,0,16);if(t[0]!==827605849)throw new Error("bad tile magic");const e=new Int32Array(n,4,3),o=[e[0]*nt,e[1]*nt,e[2]*nt],[a,s,h,l,c,i,d,m,k,f,p]=[t[4],t[5],t[6],t[7],t[8],t[9],t[10],t[11],t[12],t[13],t[14]];let r=64;const M=new Uint32Array(n.slice(r,r+16*a));r+=16*a;const g=new Uint32Array(n.slice(r,r+16*s));r+=16*s;const x=new Float32Array((l+i+p)*4),U=new Uint32Array(x.buffer);let C=0;if(h){const w=Xt(n,r,h,l,o,Rt,0,x,U,C);r=w.off,C=w.o}if(c){const w=Xt(n,r,c,i,o,mn,1,x,U,C);r=w.off,C=w.o}const A=new Float32Array(d*3);if(d){const w=new Uint16Array(n,r,3*d);for(let y=0;y<d;y++)for(let v=0;v<3;v++)A[3*y+v]=o[v]-zt+w[3*y+v]*Lt;r+=gt(6*d)}let b=new Uint32Array(m*3);m&&(k&1?b=new Uint32Array(n.slice(r,r+12*m)):b=Uint32Array.from(new Uint16Array(n,r,3*m)),r+=gt((k&1?4:2)*3*m));const I=m?new Uint8Array(n.slice(r,r+m)):new Uint8Array(0);if(m&&(r+=gt(m)),f){const w=new Uint8Array(n.slice(r,r+f));r+=gt(f),Xt(n,r,f,p,o,Rt,0,x,U,C,w)}const u=new Float32Array(d*3),B=new Uint8Array(d);for(let w=0;w<m;w++){const y=b[3*w],v=b[3*w+1],z=b[3*w+2],L=A[3*v]-A[3*y],T=A[3*v+1]-A[3*y+1],H=A[3*v+2]-A[3*y+2],S=A[3*z]-A[3*y],O=A[3*z+1]-A[3*y+1],P=A[3*z+2]-A[3*y+2],E=T*P-H*O,q=H*S-L*P,j=L*O-T*S;for(const _ of[y,v,z])u[3*_]+=E,u[3*_+1]+=q,u[3*_+2]+=j,B[_]=I[w]}for(let w=0;w<d;w++){const y=Math.hypot(u[3*w],u[3*w+1],u[3*w+2])||1;u[3*w]/=y,u[3*w+1]/=y,u[3*w+2]/=y}return{lo:o,inst:M,ghost:g,free:x,verts:A,normals:u,tris:b,fmem:I,vmem:B,nI:a,nG:s}}const Ht=n=>(n<<16>>16)/32767*Mn;function rn(n,t,e,o){const a=n[4*t],s=n[4*t+1],h=n[4*t+2],l=n[4*t+3];o[0]=e[0]-zt+(s&65535)*Lt,o[1]=e[1]-zt+(s>>>16)*Lt,o[2]=e[2]-zt+(h&65535)*Lt;const c=Ht(h>>>16),i=Ht(l&65535),d=Ht(l>>>16),m=Math.sqrt(Math.max(0,1-c*c-i*i-d*d)),k=a>>>14&3,f=k===0?[m,c,i,d]:k===1?[c,m,i,d]:k===2?[c,i,m,d]:[c,i,d,m],[p,r,M,g]=f;return o[3]=p*p+r*r-M*M-g*g,o[4]=2*(r*M-p*g),o[5]=2*(r*g+p*M),o[6]=2*(r*M+p*g),o[7]=p*p-r*r+M*M-g*g,o[8]=2*(M*g-p*r),o[9]=2*(r*g-p*M),o[10]=2*(M*g+p*r),o[11]=p*p-r*r-M*M+g*g,{type:a&16383,region:a>>>16&31,flags:a>>>21}}class en{h;tbl;mask;start;pts;inst;constructor(t,e,o,a){this.h=a;let s=1;for(;s<o*2;)s<<=1;s=Math.max(s,1024),this.tbl=s,this.mask=s-1;const h=new Uint32Array(o),l=new Uint32Array(s+1);for(let i=0;i<o;i++)h[i]=this.k(Math.floor(t[3*i]/a),Math.floor(t[3*i+1]/a),Math.floor(t[3*i+2]/a)),l[h[i]+1]++;for(let i=0;i<s;i++)l[i+1]+=l[i];this.start=l.slice(),this.pts=new Float32Array(o*3),this.inst=new Int32Array(o);const c=l.slice(0,s);for(let i=0;i<o;i++){const d=c[h[i]]++;this.pts[3*d]=t[3*i],this.pts[3*d+1]=t[3*i+1],this.pts[3*d+2]=t[3*i+2],this.inst[d]=e[i]}}k(t,e,o){return(Math.imul(t,73856093)^Math.imul(e,19349663)^Math.imul(o,83492791))&this.mask}near(t,e,o,a,s=-1){const h=this.h,l=Math.floor(t/h),c=Math.floor(e/h),i=Math.floor(o/h),d=a*a;for(let m=-1;m<=1;m++)for(let k=-1;k<=1;k++)for(let f=-1;f<=1;f++){const p=this.k(l+m,c+k,i+f);for(let r=this.start[p],M=this.start[p+1];r<M;r++){const g=this.pts[3*r]-t,x=this.pts[3*r+1]-e,U=this.pts[3*r+2]-o;if(g*g+x*x+U*U<d&&this.inst[r]!==s)return!0}}return!1}}function on(n){let t=0;const e=new Float64Array(12);for(const h of n)t+=At[h.w[4*h.i]&16383];const o=new Float32Array(t*3),a=new Int32Array(t);let s=0;for(const h of n){const{type:l}=rn(h.w,h.i,h.lo,e),c=_t[l],i=At[l];for(let d=0;d<i;d++){const m=D[3*(c+d)],k=D[3*(c+d)+1],f=D[3*(c+d)+2];o[3*s]=e[0]+e[3]*m+e[4]*k+e[5]*f,o[3*s+1]=e[1]+e[6]*m+e[7]*k+e[8]*f,o[3*s+2]=e[2]+e[9]*m+e[10]*k+e[11]*f,a[s]=h.id,s++}}return{xyz:o,owner:a,n:t}}function Ot(n){let t=n>>>0;return()=>{t=t+1831565813>>>0;let e=t;return e=Math.imul(e^e>>>15,e|1),e^=e+Math.imul(e^e>>>7,e|61),((e^e>>>14)>>>0)/4294967296}}const Pt=(n,t,e)=>{let o=Math.imul(n^2654435769,2246822507)^Math.imul(t+1663821227,3266489909)^Math.imul(e+668265263,374761393);return o^=o>>>16,o=Math.imul(o,2146121005),o^=o>>>15,o=Math.imul(o,2221713035),o^=o>>>16,o>>>0};function R(n){return Math.sqrt(-2*Math.log(n()+1e-12))*Math.cos(6.2831853*n())}function gn(n,t,e){const o=Math.abs(n)+Math.abs(t)+Math.abs(e)||1;let a=n/o,s=t/o;if(e<0){const c=a;a=(1-Math.abs(s))*(c>=0?1:-1),s=(1-Math.abs(c))*(s>=0?1:-1)}const h=Math.max(0,Math.min(255,Math.round((a+1)*127.5))),l=Math.max(0,Math.min(255,Math.round((s+1)*127.5)));return h|l<<8}const An=(()=>{const n=[20.3,14.9,14.3,8.3,1.8,12,28.4],t=n.reduce((o,a)=>o+a);let e=0;return n.map(o=>e+=o/t)})(),xn=.65,vn=.95;function ct(n,t,e,o){return n.n[0]*t+n.n[1]*e+n.n[2]*o-n.d}function kn(n,t,e,o){const a=[],s=[],h={lipCand:0,lipMasked:0,glucan:0,mannan:0},l=new Float64Array(12),c=(f,p,r,M)=>{if(e.mode==="shell")return lt(tt,J,f,p,r)<70+M;const g=ct(e,f,p,r);return g>=e.lo-M&&g<=e.hi+M};if(o.lipids&&t.tris.length&&e.mode==="slab"){const f=[],p=t.inst.length/4,r=t.ghost.length/4;for(let b=0;b<p;b++)t.inst[4*b]>>>21&1&&f.push({w:t.inst,i:b,lo:t.lo,id:b});for(let b=0;b<r;b++)f.push({w:t.ghost,i:b,lo:t.lo,id:p+b});const M=on(f),g=new en(M.xyz,M.owner,M.n,1),x=t.verts,U=t.normals,C=t.tris,A=C.length/3;for(let b=0;b<A;b++){const I=C[3*b],u=C[3*b+1],B=C[3*b+2],w=x[3*I],y=x[3*I+1],v=x[3*I+2],z=x[3*u],L=x[3*u+1],T=x[3*u+2],H=x[3*B],S=x[3*B+1],O=x[3*B+2],P=ct(e,w,y,v),E=ct(e,z,L,T),q=ct(e,H,S,O);if(Math.max(P,E,q)<e.lo-3||Math.min(P,E,q)>e.hi+3)continue;const j=z-w,_=L-y,Q=T-v,W=H-w,$=S-y,V=O-v,F=.5*Math.hypot(_*V-Q*$,Q*W-j*V,j*$-_*W),K=t.fmem[b],Y=K===8?[1]:[-1,1];for(const G of Y){const N=Ot(Pt(n,b,G+7)),ft=F/xn,et=Math.floor(ft)+(N()<ft-Math.floor(ft)?1:0);if(!et)continue;const xt=Math.ceil(Math.sqrt(et)),jt=xt*xt,it=new Uint16Array(jt);for(let X=0;X<jt;X++)it[X]=X;for(let X=0;X<et;X++){const ht=X+Math.floor(N()*(jt-X)),Z=it[X];it[X]=it[ht],it[ht]=Z}for(let X=0;X<et;X++){const ht=it[X],Z=Math.floor(Math.sqrt(ht)),Yt=ht-Z*Z,at=Yt>>1,an=(Yt&1)===0;let ut,yt,It,Ct,Et,St;an?(ut=Z,yt=at,It=Z+1,Ct=at,Et=Z+1,St=at+1):(ut=Z,yt=at,It=Z+1,Ct=at+1,Et=Z,St=at+1);let dt=N(),wt=N();dt+wt>1&&(dt=1-dt,wt=1-wt);const cn=ut+dt*(It-ut)+wt*(Et-ut),Zt=yt+dt*(Ct-yt)+wt*(St-yt),st=(cn-Zt)/xt,rt=Zt/xt,$t=1-st-rt,ln=w+st*j+rt*W,fn=y+st*_+rt*$,hn=v+st*Q+rt*V;let mt=$t*U[3*I]+st*U[3*u]+rt*U[3*B],Mt=$t*U[3*I+1]+st*U[3*u+1]+rt*U[3*B+1],pt=$t*U[3*I+2]+st*U[3*u+2]+rt*U[3*B+2];const Nt=Math.hypot(mt,Mt,pt)||1;mt/=Nt,Mt/=Nt,pt/=Nt;const un=N();let bt=0;for(;bt<6&&un>An[bt];)bt++;h.lipCand++;let vt=ln,kt=fn,Ut=hn;K===8&&(vt-=mt*2.4,kt-=Mt*2.4,Ut-=pt*2.4);const Dt=ct(e,vt,kt,Ut);if(Dt<e.lo-2.5||Dt>e.hi+2.5)continue;const yn=bt===5?[1.55,.75]:[1.95,1.25,.55];let Ft=!1;for(const Gt of yn)if(g.near(vt+mt*G*Gt,kt+Mt*G*Gt,Ut+pt*G*Gt,vn)){Ft=!0;break}if(Ft){h.lipMasked++;continue}const dn=(gn(mt,Mt,pt)|bt<<16|(G>0?1:0)<<19|(K&15)<<20)>>>0;a.push(vt,kt,Ut,dn)}}}}if(o.wall&&J.length>0){const f=J[0]-J[2],p=J[0]-J[1],r=t.lo,M=t.inst.length/4,g=[];for(let u=0;u<M;u++)(t.inst[4*u]>>>16&31)===9&&g.push({w:t.inst,i:u,lo:r,id:u});const x=on(g),U=new en(x.xyz,x.owner,x.n,1.4),C=(u,B,w)=>u>=r[0]&&u<r[0]+nt&&B>=r[1]&&B<r[1]+nt&&w>=r[2]&&w<r[2]+nt,A=self.__vol?.[n]??0,b=Math.round(A*Wt.glucanPerNm3);if(b>0){const u=Ot(Pt(n,1234567,1)),w=Math.ceil(b/48);for(let y=0;y<w;y++){const v=24+Math.floor(u()*49);let z=0,L=0,T=0,H=!1;for(let W=0;W<40&&!H;W++){z=r[0]+u()*nt,L=r[1]+u()*nt,T=r[2]+u()*nt;const $=lt(tt,J,z,L,T);H=$<f-2&&$>p+1}if(!H)continue;let[S,O,P]=Qt(tt,J,z,L,T),E=R(u),q=R(u),j=R(u),_=E*S+q*O+j*P;E-=_*S,q-=_*O,j-=_*P;let Q=Math.hypot(E,q,j)||1;E/=Q,q/=Q,j/=Q;for(let W=0;W<v&&C(z,L,T);W++){c(z,L,T,1)&&!U.near(z,L,T,1.35)&&(s.push(z,L,T,(2|(y&65535)<<4)>>>0),h.glucan++),E+=.12*R(u),q+=.12*R(u),j+=.12*R(u),[S,O,P]=Qt(tt,J,z,L,T),_=E*S+q*O+j*P;const $=lt(tt,J,z,L,T),V=$>f-3&&_<0||$<p+2&&_>0?-1:.6;E-=(1-V)*_*S,q-=(1-V)*_*O,j-=(1-V)*_*P,Q=Math.hypot(E,q,j)||1,E/=Q,q/=Q,j/=Q,z+=1.2*E,L+=1.2*q,T+=1.2*j}}}let I=0;for(const u of g){const{type:B}=rn(u.w,u.i,r,l),w=l[0],y=l[1],v=l[2],z=lt(tt,J,w,y,v);if(z>p+5)continue;if(e.mode==="slab"){const q=ct(e,w,y,v);if(q<e.lo-45-Bt[B]||q>e.hi+45+Bt[B]){I+=80;continue}}else if(z>120)continue;const L=Ot(Pt(n,u.i,99)),T=_t[B],H=At[B];let S=Wt.mannanPerProtein;const[O,P,E]=Qt(tt,J,w,y,v);for(;S>0;){const q=Math.min(Math.ceil(S),12+Math.floor(L()*27));S-=q;const j=T+Math.floor(L()*H),_=D[3*j],Q=D[3*j+1],W=D[3*j+2];let $=w+l[3]*_+l[4]*Q+l[5]*W,V=y+l[6]*_+l[7]*Q+l[8]*W,F=v+l[9]*_+l[10]*Q+l[11]*W,K=$-w,Y=V-y,G=F-v,N=Math.hypot(K,Y,G)||1;K=K/N+O,Y=Y/N+P,G=G/N+E,N=Math.hypot(K,Y,G)||1,K/=N,Y/=N,G/=N,$+=1.3*K,V+=1.3*Y,F+=1.3*G;const ft=I++;for(let et=0;et<q&&!(lt(tt,J,$,V,F)<.6);et++)c($,V,F,1)&&!U.near($,V,F,1.15,et<2?u.i:-1)&&(s.push($,V,F,(3|(ft&65535)<<4)>>>0),h.mannan++),K+=.35*R(L)+.18*O,Y+=.35*R(L)+.18*P,G+=.35*R(L)+.18*E,N=Math.hypot(K,Y,G)||1,K/=N,Y/=N,G/=N,$+=1.05*K,V+=1.05*Y,F+=1.05*G}}}const i=new Float32Array(a.length),d=new Uint32Array(i.buffer);for(let f=0;f<a.length;f+=4)i[f]=a[f],i[f+1]=a[f+1],i[f+2]=a[f+2],d[f+3]=a[f+3];const m=new Float32Array(s.length),k=new Uint32Array(m.buffer);for(let f=0;f<s.length;f+=4)m[f]=s[f],m[f+1]=s[f+1],m[f+2]=s[f+2],k[f+3]=s[f+3];return{lipids:i,free:m,stats:h}}self.onmessage=async n=>{const t=n.data;try{if(t.cmd==="init"){Vt=bn(t),self.__vol=t.volGlucan,await Vt,self.postMessage({cmd:"ready"});return}if(await Vt,t.cmd==="load"){const e=await Tt(`tiles/${t.ix}_${t.iy}_${t.iz}.bin.gz`),o=nn(e),a=o.verts.length/3,s=new Float32Array(a*8),h=new Uint32Array(s.buffer);for(let i=0;i<a;i++)s[8*i]=o.verts[3*i],s[8*i+1]=o.verts[3*i+1],s[8*i+2]=o.verts[3*i+2],h[8*i+3]=o.vmem[i],s[8*i+4]=o.normals[3*i],s[8*i+5]=o.normals[3*i+1],s[8*i+6]=o.normals[3*i+2];for(ot.set(t.idx,{lo:o.lo,inst:o.inst,ghost:o.ghost,verts:o.verts,normals:o.normals,tris:o.tris,fmem:o.fmem,vmem:o.vmem});ot.size>pn;)ot.delete(ot.keys().next().value);const l=o.inst.slice(),c=o.tris.slice();self.postMessage({cmd:"loaded",idx:t.idx,inst:l,free:o.free,vb:s,ib:c,lo:o.lo},[l.buffer,o.free.buffer,s.buffer,c.buffer])}else if(t.cmd==="detail"){let e=ot.get(t.idx);if(e)ot.delete(t.idx),ot.set(t.idx,e);else{const s=nn(await Tt(`tiles/${t.ix}_${t.iy}_${t.iz}.bin.gz`));e={lo:s.lo,inst:s.inst,ghost:s.ghost,verts:s.verts,normals:s.normals,tris:s.tris,fmem:s.fmem,vmem:s.vmem},ot.set(t.idx,e)}const o=performance.now(),a=kn(t.idx,e,t.range,t.want);self.postMessage({cmd:"detail",idx:t.idx,key:t.key,lipids:a.lipids,free:a.free,stats:a.stats,ms:performance.now()-o},[a.lipids.buffer,a.free.buffer])}}catch(e){self.postMessage({cmd:"error",idx:t.idx,key:t.key,req:t.cmd,error:String(e&&e.stack||e)})}};\n',Ft=typeof self<"u"&&self.Blob&&new Blob(["URL.revokeObjectURL(import.meta.url);",jt],{type:"text/javascript;charset=utf-8"});function _n(i){let e;try{if(e=Ft&&(self.URL||self.webkitURL).createObjectURL(Ft),!e)throw"";const t=new Worker(e,{type:"module",name:i?.name});return t.addEventListener("error",()=>{(self.URL||self.webkitURL).revokeObjectURL(e)}),t}catch{return new Worker("data:text/javascript;charset=utf-8,"+encodeURIComponent(jt),{type:"module",name:i?.name})}}const Te=200,zn=173.3,ct=70,Un=60;class dt{used;n;count=0;high=0;constructor(e){this.n=e,this.used=new Uint8Array(e)}alloc(e){if(e<=0)return[];if(this.n-this.count<e)return null;const t=[];for(let n=0;n<this.n&&t.length<e;n++)this.used[n]||(this.used[n]=1,t.push(n));return this.count+=e,this.high=Math.max(this.high,t[t.length-1]+1),t}free(e){for(const t of e)this.used[t]&&(this.used[t]=0,this.count--);for(;this.high>0&&!this.used[this.high-1];)this.high--}get freeCount(){return this.n-this.count}}class Gn{idx;r;lib;G;F=null;workers=[];state;loaded=new Map;instA;lipA;freeA;slots=[];pageTile;pageLocal;tileState=new Uint32Array(Ge*Ge*Ge);desired=[];desiredSet=new Set;inflight=0;maxInflight=12;detInflight=0;maxDetInflight=3;frameNo=0;onChange=()=>{};onTileLoaded=[];bytesLoaded=0;tilesLoadedTotal=0;starved=!1;detailStats={lipids:0,glucan:0,mannan:0};detailEnabled=!0;lastLoadErr="";shellRange=null;constructor(e,t,n,s,o,r){this.idx=e,this.r=t,this.lib=n,this.G=s.G,this.F=s.F??null,this.state=new Uint8Array(e.n),this.instA=new dt(t.instPages),this.lipA=new dt(t.lipPages),this.freeA=new dt(t.freePages);for(let l=t.maxSlots-1;l>=0;l--)this.slots.push(l);this.pageTile=new Int32Array(t.instPages).fill(-1),this.pageLocal=new Uint32Array(t.instPages);const a=Array.from(e.volGlucan),d=new URL("data/",document.baseURI).href;for(let l=0;l<o;l++){const c=new _n;c.onmessage=u=>this.onMessage(u.data),c.postMessage({cmd:"init",base:d,G:s.G,F:s.F,volGlucan:a,embedded:r,wall:{glucanPerNm3:e.wall.glucan_beads_per_nm3,mannanPerProtein:e.wall.mannan_beads_per_protein}}),this.workers.push(c)}}tileCenter(e){return[this.idx.ix[e]*Te+100,this.idx.iy[e]*Te+100,this.idx.iz[e]*Te+100]}update(e,t,n,s,o,r,a){if(this.frameNo++,r){this.desired=[],this.desiredSet.clear();return}const d=e.eye,l=e.frustum(t),c=this.idx.n,u=[s[0],s[1],s[2]],m=100*(Math.abs(u[0])+Math.abs(u[1])+Math.abs(u[2])),f=[];for(let p=0;p<c;p++){const x=this.tileCenter(p),h=ue(x,d);let g=!0;for(let b=0;b<5;b++)if(l[4*b]*h[0]+l[4*b+1]*h[1]+l[4*b+2]*h[2]+l[4*b+3]<-243.3){g=!1;break}if(g){if(n){const b=V(u,x)-s[3];if(b+m+ct<0||b-m-ct>o)continue}else{const b=this.tileShellRange();if(b[5*p]>Un+87||b[5*p+1]<-87-ct)continue;const A=ue(d,x),C=5*p+2;if(b[C]*A[0]+b[C+1]*A[1]+b[C+2]*A[2]<-300)continue}f.push({i:p,d:Ae(h)})}}f.sort((p,x)=>p.d-x.d);let T=this.r.instPages*X*.94,w=this.r.maxSlots-8;const M=[];for(const p of f){const x=this.idx.nInst[p.i];if(x>T||w<=0){this.starved=!0;break}T-=Math.ceil(x/X)*X,w--,M.push(p.i)}M.length===f.length&&(this.starved=!1),this.desired=M,this.desiredSet=new Set(M);for(const p of M){const x=this.loaded.get(p);x&&(x.lastUsed=this.frameNo)}for(const p of M){if(this.inflight>=this.maxInflight)break;this.state[p]===0&&(this.state[p]=1,this.inflight++,this.workers[p%this.workers.length].postMessage({cmd:"load",idx:p,ix:this.idx.ix[p],iy:this.idx.iy[p],iz:this.idx.iz[p]}))}this.updateDetail(e,n,s,o,a)}tileShellRange(){if(this.shellRange)return this.shellRange;const e=this.idx.n,t=new Float32Array(5*e);for(let n=0;n<e;n++){const s=this.idx.ix[n]*Te,o=this.idx.iy[n]*Te,r=this.idx.iz[n]*Te;let a=1/0,d=-1/0;for(let c=0;c<=2;c++)for(let u=0;u<=2;u++)for(let m=0;m<=2;m++){const f=At(this.F,this.G,s+100*c,o+100*u,r+100*m);a=Math.min(a,f),d=Math.max(d,f)}const l=Dn(this.F,this.G,s+100,o+100,r+100,60);t[5*n]=a,t[5*n+1]=d,t.set(l,5*n+2)}return this.shellRange=t}updateDetail(e,t,n,s,o){if(!this.detailEnabled)return;const r=e.eye,a=e.pxPerNmAtUnit(o),d=[n[0],n[1],n[2]],l=Math.min(s,250),c=[];for(const m of this.desired){const f=this.loaded.get(m);if(!f)continue;const T=this.idx.nTri[m]>0,w=this.idx.volGlucan[m]>0||(this.idx.labelMask[m]&2)!==0;if(!T&&!w)continue;const M=this.tileCenter(m),p=Math.max(1,Ae(ue(M,r))-zn);a/p*.55<.8||!t&&!w||c.push({t:f,d:p})}c.sort((m,f)=>m.d-f.d);const u=new Set(c.slice(0,40).map(m=>m.t.idx));for(const m of this.loaded.values())m.det&&!u.has(m.idx)&&!m.detPending&&this.dropDetail(m);for(const{t:m}of c.slice(0,40)){if(this.detInflight>=this.maxDetInflight)break;if(m.detPending)continue;const f=m.det;let T=!1;if(f&&(t?T=!f.shell&&V(f.n,d)>.99985&&f.a<=n[3]+.5&&f.b>=n[3]+l-.5:T=f.shell),T)continue;const w=Math.max(10,.35*l),M=t?{n:d,d:n[3],lo:-w,hi:l+w,mode:"slab"}:{n:d,d:n[3],lo:0,hi:0,mode:"shell"};m.detPending=!0,this.detInflight++;const p=t?`${n[3]-w}|${n[3]+l+w}`:"shell";m._req={n:d,a:n[3]-w,b:n[3]+l+w,shell:!t};const x=this.idx.nTri[m.idx]>0,h=this.idx.volGlucan[m.idx]>0||(this.idx.labelMask[m.idx]&2)!==0;this.workers[m.idx%this.workers.length].postMessage({cmd:"detail",idx:m.idx,key:p,ix:this.idx.ix[m.idx],iy:this.idx.iy[m.idx],iz:this.idx.iz[m.idx],range:M,want:{lipids:x&&t,wall:h}})}}dropDetail(e){this.lipA.free(e.lipPages),this.freeA.free(e.detPages),this.clearSPages(this.r.lipPagesBuf,e.lipPages),this.clearSPages(this.r.freePagesBuf,e.detPages),e.lipPages=[],e.detPages=[],e.nLip=0,e.nDet=0,e.det=null,e.underMesh&&(e.underMesh=!1,this.r.setMeshSlot(e.slot,!1)),this.syncHigh()}clearSPages(e,t){const n=new Uint32Array(4);for(const s of t)this.r.device.queue.writeBuffer(e,s*16,n)}syncHigh(){this.r.usedInstPages=this.instA.high,this.r.usedLipPages=this.lipA.high,this.r.usedFreePages=this.freeA.high}evictOne(e){let t=null;for(const n of this.loaded.values())n.idx===e||this.desiredSet.has(n.idx)||(!t||n.lastUsed<t.lastUsed)&&(t=n);return t?(this.evict(t),!0):!1}evict(e){const t=this.r.device.queue,n=new Float32Array(4);for(const s of e.pages)t.writeBuffer(this.r.pageOrg,s*16,n),this.pageTile[s]=-1;this.instA.free(e.pages),this.freeA.free(e.basePages),this.clearSPages(this.r.freePagesBuf,e.basePages),this.lipA.free(e.lipPages),this.clearSPages(this.r.lipPagesBuf,e.lipPages),this.freeA.free(e.detPages),this.clearSPages(this.r.freePagesBuf,e.detPages),e.vb?.destroy(),e.ib?.destroy(),this.slots.push(e.slot),this.loaded.delete(e.idx),this.state[e.idx]=0,this.tileState[gt(this.idx.ix[e.idx],this.idx.iy[e.idx],this.idx.iz[e.idx])]=0,t.writeBuffer(this.r.tileState,0,this.tileState),this.syncHigh()}allocWithEvict(e,t,n){let s=e.alloc(t);for(;!s&&this.evictOne(n);)s=e.alloc(t);return s}onMessage(e){if(e.cmd!=="ready"){if(e.cmd==="error"){if(console.warn("worker error",e.req,e.idx,e.error),this.lastLoadErr=String(e.error).split(`
`)[0],e.req==="load"&&(this.state[e.idx]=0,this.inflight--),e.req==="detail"){this.detInflight--;const t=this.loaded.get(e.idx);t&&(t.detPending=!1)}return}e.cmd==="loaded"&&(this.inflight--,this.onLoaded(e)),e.cmd==="detail"&&(this.detInflight--,this.onDetail(e))}}onLoaded(e){const t=e.idx;if(!this.desiredSet.has(t)){this.state[t]=0,this.onChange();return}const n=e.inst,s=n.length/4,o=this.r.device.queue;this.slots.length||this.evictOne(t);const r=this.allocWithEvict(this.instA,Math.ceil(s/X),t);if(!r||!this.slots.length){this.state[t]=0,this.starved=!0,r&&this.instA.free(r);return}const a=this.slots.pop(),d=e.lo,l=new ArrayBuffer(16),c=new Float32Array(l),u=new Uint32Array(l);c[0]=d[0]-80,c[1]=d[1]-80,c[2]=d[2]-80,r.forEach((g,b)=>{const A=Math.min(X,s-b*X);o.writeBuffer(this.r.inst,g*X*16,n.buffer,b*X*16,A*16),u[3]=A,o.writeBuffer(this.r.pageOrg,g*16,l),this.pageTile[g]=t,this.pageLocal[g]=b*X});const m=e.free,f=m.length/4;let T=[],w=!1;if(f){const g=this.allocWithEvict(this.freeA,Math.ceil(f/Q),t);g?(T=g,this.writeSpheres(this.r.freePool,this.r.freePagesBuf,g,m,a)):w=!0}let M=null,p=null;const x=e.ib.length;x&&(M=this.r.device.createBuffer({size:e.vb.byteLength,usage:GPUBufferUsage.VERTEX|GPUBufferUsage.COPY_DST}),o.writeBuffer(M,0,e.vb),p=this.r.device.createBuffer({size:Math.ceil(e.ib.byteLength/4)*4,usage:GPUBufferUsage.INDEX|GPUBufferUsage.COPY_DST}),o.writeBuffer(p,0,e.ib),this.r.setMeshSlot(a,!1));const h={idx:t,slot:a,lo:d,inst:n,nInst:s,pages:r,basePages:T,nBase:f,lipPages:[],detPages:[],nLip:0,nDet:0,vb:M,ib:p,nIdx:x,lastUsed:this.frameNo,det:null,detPending:!1,underMesh:!1,baseFreeMissing:w};this.loaded.set(t,h),this.state[t]=2,this.tileState[gt(this.idx.ix[t],this.idx.iy[t],this.idx.iz[t])]=1,o.writeBuffer(this.r.tileState,0,this.tileState),this.bytesLoaded+=this.idx.gz[t],this.tilesLoadedTotal++,this.syncHigh();for(const g of this.onTileLoaded.slice())g(h);this.onChange()}writeSpheres(e,t,n,s,o){const r=this.r.device.queue,a=s.length/4;n.forEach((d,l)=>{const c=Math.min(Q,a-l*Q);r.writeBuffer(e,d*Q*16,s.buffer,s.byteOffset+l*Q*16,c*16),r.writeBuffer(t,d*16,new Uint32Array([c,o,1,0]))})}onDetail(e){const t=this.loaded.get(e.idx);if(!t)return;t.detPending=!1;const n=t._req;this.lipA.free(t.lipPages),this.clearSPages(this.r.lipPagesBuf,t.lipPages),this.freeA.free(t.detPages),this.clearSPages(this.r.freePagesBuf,t.detPages),t.lipPages=[],t.detPages=[];const s=e.lipids,o=e.free,r=s.length/4,a=o.length/4;if(r){const l=this.lipA.alloc(Math.ceil(r/Q))??this.freeDetailSpace(this.lipA,Math.ceil(r/Q),t.idx);l&&(t.lipPages=l,this.writeSpheres(this.r.lipPool,this.r.lipPagesBuf,l,s,t.slot))}if(a){const l=this.freeA.alloc(Math.ceil(a/Q))??this.freeDetailSpace(this.freeA,Math.ceil(a/Q),t.idx);l&&(t.detPages=l,this.writeSpheres(this.r.freePool,this.r.freePagesBuf,l,o,t.slot))}t.nLip=t.lipPages.length?r:0,t.nDet=t.detPages.length?a:0,t.det={n:n.n,a:n.a,b:n.b,shell:n.shell};const d=t.nLip>0;d!==t.underMesh&&t.nIdx&&(t.underMesh=d,this.r.setMeshSlot(t.slot,d)),this.detailStats.lipids+=e.stats.lipCand,this.syncHigh(),this.onChange()}freeDetailSpace(e,t,n){const s=[...this.loaded.values()].filter(o=>o.idx!==n&&(o.lipPages.length||o.detPages.length)&&!this.desiredSet.has(o.idx));for(const o of s){this.dropDetail(o);const r=e.alloc(t);if(r)return r}return null}meshDraws(){const e=[];for(const t of this.loaded.values())t.nIdx&&t.vb&&t.ib&&e.push({slot:t.slot,vb:t.vb,ib:t.ib,nIdx:t.nIdx,under:t.underMesh});return e}instanceAt(e){const t=Math.floor(e/X),n=this.pageTile[t];if(n<0)return null;const s=this.loaded.get(n);if(!s)return null;const o=this.pageLocal[t]+e%X;return{t:s,local:o,...this.decode(s,o)}}decode(e,t){const n=e.inst,s=n[4*t],o=n[4*t+1],r=n[4*t+2],a=360/65535,d=[e.lo[0]-80+(o&65535)*a,e.lo[1]-80+(o>>>16)*a,e.lo[2]-80+(r&65535)*a];return{type:s&16383,region:s>>>16&31,flags:s>>>21,p:d}}poolIndex(e,t){return e.pages[Math.floor(t/X)]*X+t%X}tileBySlot(e){for(const t of this.loaded.values())if(t.slot===e)return t;return null}countLoaded(e){let t=0;for(const n of this.loaded.values()){const s=n.inst;for(let o=0;o<n.nInst;o++)(s[4*o]&16383)===e&&t++}return t}nearestInTile(e,t,n){let s=-1,o=1e18;for(let r=0;r<e.nInst;r++){if((e.inst[4*r]&16383)!==t)continue;const a=this.decode(e,r),d=(a.p[0]-n[0])**2+(a.p[1]-n[1])**2+(a.p[2]-n[2])**2;d<o&&(o=d,s=r)}return s}totals(){let e=0,t=0,n=0,s=0;for(const o of this.loaded.values())e+=o.nInst,t+=o.nLip,n+=o.nBase+o.nDet,s+=o.nIdx/3;return{tiles:this.loaded.size,inst:e,lip:t,free:n,mesh:s,inflight:this.inflight+this.detInflight,desired:this.desired.length}}pendingForView(){let e=0;for(const n of this.desired)this.state[n]!==2&&e++;let t=0;for(const n of this.loaded.values())n.detPending&&t++;return e+t}}class qn{G;PT;PQ;PC;pores;c;grid;n3;F;constructor(e){this.F=e.F??null,this.G=e.G,this.PT=e.PT,this.PQ=e.PQ,this.PC=e.PC,this.pores=e.pores,this.c=e.consts,this.grid=e.grid,this.n3=e.grid.n**3}cell(e,t,n){const s=this.grid.n,o=this.grid.ext,r=2*o/s,a=d=>Math.min(s-1,Math.max(0,Math.floor((d+o)/r)));return(a(e)*s+a(t))*s+a(n)}nucSdf(e,t,n){if(Ze(this.F)>It)return yt(this.F,e,t,n,It);const s=this.G,o=e-s[3],r=t-s[4],a=n-s[5],d=(s[14]*o+s[15]*r+s[16]*a)/s[11],l=(s[17]*o+s[18]*r+s[19]*a)/s[12],c=(s[20]*o+s[21]*r+s[22]*a)/s[13],u=Math.hypot(d,l,c),m=Math.hypot(d/s[11],l/s[12],c/s[13]);return m<1e-12?-Math.min(s[11],s[12],s[13]):u*(u-1)/m}sdfCer(e,t,n,s,o){const r=this.PQ,a=8*e,d=Math.hypot(t,n,s)||1e-9,l=Math.min(1,Math.max(-1,(t*r[a]+n*r[a+1]+s*r[a+2])/d)),c=Math.max(Math.acos(l)-r[a+4],0)*d,u=this.G[2]-r[a+5]-o;return Math.sqrt(c*c+u*u)-r[a+6]}sdfPrim(e,t,n,s){const o=this.PT[e],r=this.PQ,a=8*e;if(o===0)return Math.hypot(t-r[a],n-r[a+1],s-r[a+2])-r[a+3];if(o===1){const f=t-r[a],T=n-r[a+1],w=s-r[a+2],M=r[a+3]-r[a],p=r[a+4]-r[a+1],x=r[a+5]-r[a+2],h=Math.min(1,Math.max(0,(f*M+T*p+w*x)/(M*M+p*p+x*x)));return Math.hypot(f-M*h,T-p*h,w-x*h)-r[a+6]}if(o===2){const f=t-r[a],T=n-r[a+1],w=s-r[a+2],M=r[a+3],p=r[a+4],x=r[a+5],h=f*M+T*p+w*x,g=Math.hypot(f-h*M,T-h*p,w-h*x),b=Math.max(g-r[a+6],0);return Math.sqrt(b*b+h*h)-r[a+7]}const d=Math.hypot(t,n,s)||1e-9,l=Math.min(1,Math.max(-1,(t*r[a]+n*r[a+1]+s*r[a+2])/d)),c=Math.acos(l),u=Math.max(c-r[a+4],0)*d,m=d-r[a+5];return Math.sqrt(u*u+m*m)-r[a+6]}label(e,t,n){const s=this.G,o=this.c.PERI_W,r=this.c.PORE_R,a=At(this.F,s,e,t,n),d=a-(s[0]-s[2]);if(a<0)return 0;if(a<s[0]-s[1])return 1;if(d<0)return 2;const l=this.cell(e,t,n),c=this.grid,u=[1e9,1e9,1e9,1e9,1e9,1e9,1e9,1e9,1e9];Ze(this.F)>Lt&&(u[3]=yt(this.F,e,t,n,Lt));let m=1e9,f=1e9;for(let g=c.prim_off[l];g<c.prim_off[l+1];g++){const b=c.prim_idx[g],A=this.PT[b]===3?this.sdfCer(b,e,t,n,d):this.sdfPrim(b,e,t,n),C=this.PC[b];C===7?u[1]=Math.min(u[1],A):C===0?u[3]=Math.min(u[3],A):C===4?u[4]=Math.min(u[4],A):C===5?m=Math.min(m,A):C===6?f=Math.min(f,A):C===2?u[6]=Math.min(u[6],A):C===3?u[7]=Math.min(u[7],A):C===1&&(u[8]=Math.min(u[8],A))}const T=Math.max(m,-f);if(u[1]<0)return 4;const w=this.nucSdf(e,t,n),M=w+.5*o,p=r+.5*o;let x=!1;for(let g=c.pore_off[l];g<c.pore_off[l+1];g++){const b=6*c.pore_idx[g],A=this.pores,C=e-A[b],F=t-A[b+1],O=n-A[b+2],k=C*A[b+3]+F*A[b+4]+O*A[b+5];if(Math.abs(k)>150)continue;if(Math.hypot(C-k*A[b+3],F-k*A[b+4],O-k*A[b+5])<p){x=!0;break}}if(Math.abs(M)-.5*o<0&&!x)return 5;if(w<-.5*o){if(x&&w>-o-8)return 14;const g=e-s[7],b=t-s[8],A=n-s[9];return w<-o&&g*g+b*b+A*A<s[10]*s[10]?7:6}return w<0?14:u[4]<0?T<0?10:9:u[3]<0?8:u[8]<0?11:u[6]<0?12:u[7]<0?13:3}}const $n=`
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
23,000 of them carry 91,000 of the 80S ribosomes (about 4 per mRNA): the mRNA runs through each ribosome's channel
(path from PDB 7ZW0), short linkers join neighbouring ribosomes, and the rest of the transcript leaves the first and
last ribosome as 3&prime; and 5&prime; tails. The other 7,000 mRNAs carry no ribosomes. Threads that run out of room end
early (86% of the target length is reached overall). Linker DNA joins consecutive nucleosomes.</li>
<li><b>ER-bound ribosomes:</b> 4,300 ribosomes sit on translocons in the ER and outer nuclear membrane: on the Ssh1 complex
(yeast ribosome-Ssh1 cryo-EM structure 2WW9) or on Sec61 (7KAH, placed through the shared Sss1 subunit). Their number is
limited by the free translocon pools (80% assumed ribosome-bound); the PM-facing side of the cortical ER is too narrow for ribosomes.</li>
<li><b>Microtubules:</b> 21 microtubules (11,700 &alpha;/&beta;-tubulin dimers) grow from &gamma;-tubulin rings on the spindle pole body,
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
mitochondrion, Golgi/endosomes, peroxisomes, lipid droplets) are represented by signed-distance fields. The cell
outline, vacuoles and lipid droplets are those of a cell segmented from a cryo-FIB/SEM volume (EMPIAR-11416, Dumoux
<i>et al.</i> 2023), the nucleus that of a second cell; the shapes are smoothed, and the vacuoles' thickness along the
milling direction, which the segmentation loses, is reconstructed.
Near the camera every lipid is drawn individually. Class composition (Ejsing <i>et al.</i> 2009): PI 20%, PE 15%, PC 14%, PA 8%, PS 2%,
ergosterol 12%, 0.65 nm² per lipid. Lipids are removed wherever they would overlap a membrane protein.</li>
<li><b>Cell wall:</b> an inner β-1,3-glucan layer and an outer mannoprotein layer with mannan brushes. Bead counts follow a
wall mass of about 20% of dry mass.</li>
</ul>

<h3>How it was built</h3>
<p>Proteins were assigned to compartments from UniProt subcellular-location annotations, with GO as fallback. Membrane proteins were oriented by their
transmembrane segments. Everything was then packed into 200 nm tiles with collision detection on residue-level
sphere trees. Placement is not purely random: proteins of the same STRING interaction module are co-enriched in
domains a few hundred nm across, and about a fifth of all molecules were placed next to an interaction partner predicted
by AlphaFold 3 or linked in STRING (an assumed, tunable organisation; see Model details). Validation on the final model:</p>
<ul>
<li>Bead clash rate 1.9×10⁻⁶.</li>
<li>Containment violations 5.6×10⁻⁵.</li>
<li>Macromolecule densities (protein plus the RNA of ribosomes and tRNAs and the DNA of nucleosomes; lipids and small
molecules not counted): cytosol 93 mg/mL (70 protein + 23 RNA), nucleoplasm 56 mg/mL, mitochondrial matrix 111 mg/mL.</li>
<li>99.9% of requested copies were placed outside mitochondria.</li>
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
Complex Portal; UniProt and Gene Ontology annotations; IntAct interactions; STRING v12 associations (CC-BY 4.0); EMPIAR-11416 cryo-FIB/SEM volume (CC0); Alliance of Genome Resources paralogs;
lipid composition from Ejsing <i>et al.</i> 2009; total protein from Ho <i>et al.</i> 2018.</p>
`,Bt=i=>i.replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;");function Ye(i){const e=[];return i=i.replace(/`([^`]+)`/g,(t,n)=>`\0${e.push(Bt(n))-1}\0`),i=Bt(i).replace(/\[([^\]]+)\]\(([^)\s]+)\)/g,(t,n,s)=>`<a href="${s}" target="_blank" rel="noopener">${n}</a>`).replace(/\*\*([^*]+)\*\*/g,"<b>$1</b>").replace(/\*([^*\s][^*]*)\*/g,"<i>$1</i>"),i.replace(/\u0000(\d+)\u0000/g,(t,n)=>`<code>${e[+n]}</code>`)}const ut=/^([-*]|\d+\.)\s+/,ht=i=>i.length-i.trimStart().length,Rt=i=>i.trim().replace(/^\||\|$/g,"").split("|").map(e=>e.trim());function Ht(i){let e="",t=0;for(;t<i.length;){const n=i[t],s=n.trim();if(!s){t++;continue}let o;if(o=s.match(/^(#{1,6})\s+(.*)$/)){e+=`<h${o[1].length}>${Ye(o[2])}</h${o[1].length}>`,t++;continue}if(/^-{3,}$/.test(s)){e+="<hr>",t++;continue}if(s.startsWith("|")&&t+1<i.length&&/^\|[\s|:-]+\|$/.test(i[t+1].trim())){const a=Rt(s);t+=2;let d="";for(;t<i.length&&i[t].trim().startsWith("|");)d+=`<tr>${Rt(i[t]).map(l=>`<td>${Ye(l)}</td>`).join("")}</tr>`,t++;e+=`<table><thead><tr>${a.map(l=>`<th>${Ye(l)}</th>`).join("")}</tr></thead><tbody>${d}</tbody></table>`;continue}if(o=s.match(ut)){const a=/\d/.test(o[1]),d=ht(n);let l="";for(;t<i.length;){const c=i[t],u=c.trim().match(ut);if(!u||ht(c)!==d||/\d/.test(u[1])!==a)break;const m=d+u[0].length,f=[c.slice(m)];for(t++;t<i.length&&(!i[t].trim()||ht(i[t])>=m);)f.push(i[t].slice(m)),t++;for(;f.length&&!f[f.length-1].trim();)f.pop();const T=Ht(f),w=T.match(/^<p>((?:(?!<p>)[\s\S])*)<\/p>$/);l+=`<li>${w?w[1]:T}</li>`}e+=a?`<ol>${l}</ol>`:`<ul>${l}</ul>`;continue}const r=[];for(;t<i.length&&i[t].trim()&&!/^(#{1,6}\s|\||-{3,}$)/.test(i[t].trim())&&!ut.test(i[t].trim());)r.push(i[t++].trim());e+=`<p>${Ye(r.join(" "))}</p>`}return e}function Vn(i){return Ht(i.replace(/\r\n/g,`
`).split(`
`))}const Wn=`# Whole yeast cell: model description

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

The cell is 42 µm³ (mean log-phase volume; Jorgensen *et al.* 2002) with a 110 nm cell wall, the lower bound of 110–200 nm (Orlean 2012). Its outline, nucleus and vacuoles are real segmented shapes (section 3.3).

Organelle volume fractions:
- **Nucleus, vacuole and lipid droplets:** from segmented cryo-FIB/SEM cells (EMPIAR-11416; Dumoux *et al.* 2023; section 3.3), set against the literature:
  - **Nucleus 8%.** The two segmented nuclei are 8.0% and 9.7% of their cells (both from cells cut by the edge of the imaged volume, which biases the fraction up). Fluorescence gives about 7% (Jorgensen *et al.* 2007) and serial-section EM 10.5% (Yamaguchi *et al.* 2011). The nucleolus is 20% of the nucleus (Uchida *et al.* 2011).
  - **Vacuole 9%.** The median of five segmented cells (6–13%). The segmentation misses faint vacuoles, so this is a lower bound; soft X-ray tomography gave 7% (Uchida *et al.* 2011).
  - **Lipid droplets 0.34%.** The 14 segmented droplets of the cell whose outline is used (the five cells range 0.3–0.8%, from the depositors' hand-corrected labels); soft X-ray tomography gave 1%.
- **Mitochondria:** raised from the measured 1% to 2%. The assumption is that X-ray tomography under-segments thin tubules. At 1%, the mitochondrial proteome would imply about 480 mg/mL of protein.
- **Golgi, endosomes and peroxisomes:** 1% together (assumption).
- **Cortical ER:** covers 48% of the plasma membrane (Wolf *et al.* 2012).

Measured on the final voxelised geometry (4 nm grid):

| Compartment | Volume (µm³) | Share of cell |
|---|---|---|
| Cytosol | 26.62 | 63.5% |
| Cell wall (glucan 4.23, mannan 2.09) | 6.32 | 15.1% |
| Vacuole (13 compartments) | 3.75 | 8.9% |
| Nucleoplasm | 2.26 | 5.4% |
| Nucleolus | 0.67 | 1.6% |
| Mitochondrial matrix | 0.59 | 1.4% |
| Mitochondrial intermembrane space (incl. crista lumen) | 0.24 | 0.6% |
| ER lumen: cortical ER | 0.53 | 1.3% |
| ER lumen: perinuclear space | 0.36 | 0.9% |
| Lipid droplets (14) | 0.14 | 0.3% |
| Golgi / endosomes | 0.34 | 0.8% |
| Peroxisomes | 0.08 | 0.2% |
| **Total** (incl. 0.04 µm³ of pore channels) | **41.93** | |

### 3.2 Membranes

Membrane mid-surface areas, measured on the final geometry:

| Membrane | Area (µm²) |
|---|---|
| Plasma membrane | 54.3 |
| ER (cortical sheets and perinuclear) | 54.2 |
| Nuclear envelope | 22.5 |
| Mitochondrial inner membrane | 28.0 |
| Mitochondrial outer membrane | 11.7 |
| Vacuole (all lobes) | 29.2 |
| Golgi / endosomes | 9.5 |
| Peroxisomes | 1.5 |

### 3.3 Organelle shapes

**A composite of real cells.** The cell outline, the nucleus and the vacuoles are real shapes from a cryo-plasma FIB/SEM volume of plunge-frozen log-phase yeast (EMPIAR-11416; Dumoux *et al.* 2023), segmented into meshes by a separate pipeline (Cellpose plus manual curation):
- **Outline, vacuoles and lipid droplets:** cell 08, the only segmented cell not cut by the edge of the imaged volume (44.2 µm³, scaled by 0.983 to 42 µm³).
- **Nucleus:** cell 05's. Cell 08's nucleus was not segmented; cell 05's is one of two real ones, and the other (cell 04's) has its short axis along the milling direction of a cut cell. It is scaled to 8% of the cell volume and placed where the segmentation pipeline put cell 08's stand-in nucleus (the cell's most open region, where its unsegmented nucleus probably is), oriented along that stand-in's principal axes, and moved 155 nm inward to stay 260 nm below the outer surface.
- **Not used:** the bundle's procedural mitochondria, cortical ER and peroxisomes.

**Cleaning the meshes** (\`pipeline/cellshape.py\`; the meshes are in \`pipeline/cell_meshes.npz\`, extracted by \`scripts/cell_meshes.py\`):
- **Smoothing:** each part is voxelised on a 20 nm grid, closed morphologically (radius 100 nm for the outline, 60 nm for the rest), which fills grooves a curtaining stripe of the FIB data cut into the outline and some vacuoles. It is then Gaussian-smoothed (σ 50 nm outline, 40 nm nucleus, 30 nm vacuoles), which removes the mesh facets and the 50 nm slice steps.
- **Vacuole z-correction:** the segmentation loses the tops and bottoms of vacuoles. Across the 107 segmented vacuoles, z extent divided by in-plane extent is 0.55–0.61, against 1.02 for the cells, and cell 08's are mostly 0.2–0.4. The in-plane outlines are kept, and each lobe is corrected in four steps:
  - it is stretched along z about its centroid until its z extent matches its smaller in-plane extent (stretch 1.0–3.0, median 2.0);
  - it is smoothed more along z (σ 90 nm);
  - its distance field is blended 1:1 with that of its moment-matched ellipsoid, which rounds the missing caps;
  - one thin streak along the curtaining stripe (in-plane aspect 4.9) is dropped as an artefact.
- **Vacuole volume and spacing:** the remaining 15 lobes are scaled together to 9% of the cell volume. They are carved apart where they come within 30 nm of each other or of the nucleus, or closer than 220 nm to the outer surface, and smoothed again (σ 20 nm). Several lobes press against each other and against the nucleus, as in the segmented cells; the result is 13 separate vacuole compartments with a surface-to-volume ratio of 7.8 µm⁻¹ (segmented cells: 6.8–9.9 µm⁻¹).
- **Lipid droplets:** the 14 segmented droplets are kept as volume-equivalent spheres at their positions (radii 64–171 nm); 5 that now overlap a vacuole or the nucleus are moved to the nearest free space. Spheres also remove the droplets' z-blur (z extent / in-plane extent 1.25).

**Distance fields.** Each cleaned surface is re-meshed (marching cubes) and the exact distance to it is computed at the grid nodes near the surface (closest triangle). Between nodes the distance is interpolated linearly, accurate to about 0.2 nm, so membranes, offsets (the 35 nm mannan layer, the 110 nm wall, the cortical-ER sheets) and protein orientations are as precise as with the analytic shapes. The fields are stored in sparse 8³-voxel bricks within a band around each surface, and the viewer uses the same fields (re-bricked, 7.5 MB compressed). \`scripts/s2_geometry.py --sphere\` builds the previous analytic cell instead.

Remaining shapes and fixed dimensions:

| Element | Shape and dimensions |
|---|---|
| Nuclear envelope | Perinuclear space 34 nm between the membrane mid-surfaces, following the real nucleus. 117 pores of 45 nm radius at the mid-plane, spread uniformly over the envelope area and set along the local normal. The nucleolus is a crescent on the side away from the vacuoles (20% of the nucleus). |
| Cortical ER | Sheets 27 nm beneath the plasma membrane, following its real shape, with a 20 nm lumen. |
| Mitochondria | One reticulum of tubules with 150 nm radius, 440–640 nm below the outer surface. Intermembrane space 8 nm. Lamellar cristae with an 18 nm lumen (rim-tube radius 9 nm, chosen to fit ATP synthase dimers; section 5.6). |
| Other organelles | Golgi/endosome cisternae and peroxisomes are simple primitives. |

**Arrangement.** The outline, vacuoles and lipid droplets keep their arrangement in cell 08. The nucleus sits in that cell's open region. The spindle pole body faces away from the nucleolus, and the mitochondria form branched tubules along the cortex; both are procedural.

**Shape statistics** of all five segmented cells (\`scripts/organelle_stats.py\` → \`pipeline/organelle_stats.json\`) were used to choose the composite's parts and to check it. The cell 08 outline is mildly non-spherical: its volume-moment semi-axes are 2.41, 2.22 and 2.00 µm.

**Known geometric deviations from EM.** Serial-section EM (Yamaguchi *et al.* 2011) gives a nucleus of 10.5% of cell volume; this model has 8%. Mitochondria are 1.7% there and about 2% here.

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

- **Number:** 21 microtubules, made of 11,661 αβ-tubulin dimers:
  - 16 kinetochore microtubules, 229 nm mean length;
  - 2 other nuclear microtubules, 361 nm;
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
  - **Result:** 91,000 ribosomes on 23,000 mRNAs.
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
- **Copies per tile:** each species' expected copy number per compartment is split over tiles in proportion to the compartment's volume (or membrane area) in each tile, modulated by the large-scale organisation of section 8.4. A multinomial draw with stochastic rounding conserves the cell totals.
- **Concentrations:** total macromolecular crowding is uniform within each compartment; the composition varies (section 8.4).

### 8.2 Algorithm

Random sequential addition, run independently in each tile.

**Placement:**
- **Order:** membrane proteins first, then the rest by decreasing size.
- **Soluble molecules:** a uniformly random position in the target compartment and a uniformly random orientation, or, for a share of copies, a position next to an interaction partner (section 8.4).
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
- **The shortfall is concentrated in mitochondria.** Placed copies are 97.8% of requested overall and 99.9% outside mitochondria. In the mitochondrial matrix only 62% fit, so matrix protein density is under-represented.
- **Interaction-driven organisation is a hypothesis layer.** Section 8.4 places partners near each other and co-enriches network modules, but with assumed strengths. There are no explicit condensates, and the layer is not fitted to any measured protein positions.

### 8.4 Organisation from association data

Uniform random placement is itself an assumption: that proteins have no spatial organisation beyond their compartment. The model replaces it with a weak, data-driven organisation at two scales. Neither scale is calibrated against measured protein positions, which hardly exist at this resolution. The strengths are tuned so that the packing stays random-looking but is measurably organised, and total crowding stays uniform. Code: \`pipeline/organisation.py\`; \`scripts/s3_pack.py --random\` builds the same cell without this layer.

**Data:**
- **STRING v12** (Szklarczyk *et al.* 2023, Nucleic Acids Res. 51:D638): 102,672 links with combined score ≥ 700 between yeast proteins. STRING includes indirect associations (shared pathways, co-expression, text mining) that a structure prediction cannot see.
- **AlphaFold 3 pool predictions** (in-house, unpublished): every protein was predicted in many random multi-protein pools, giving a chain-pair score for about 10.3 million protein pairs. The best (lowest) chain-pair PAE min-min per pair is used.
  - **Calibration against known interactions:** pairs at ≤ 2 Å are 73% Complex Portal pairs and 72% have IntAct evidence; pairs at 2–4 Å, 43% and 48%; pairs at 6–8 Å are still 50 and 10 times above the background of 0.06% and 0.6%.
  - **Shared localization:** confident pairs also share compartments much more often than random pairs (overlap of their localization weights 0.56–0.62 vs 0.22).
  - **Used:** 15,262 pairs at ≤ 8 Å after filtering. A pair is dropped when a paralog of one partner scores better with the other, which removes homologous-interface artefacts such as Pol II subunits scored against their Pol III paralogs.

**Large scale (hundreds of nm): network-informed sampling of where each species goes.**
- **Modules:** STRING (≥ 700) is split by Louvain clustering (resolution 2) into 31 modules of ≥ 10 proteins, covering 5,459 proteins. Examples: glycolysis and translation factors; actin and cortical proteins; trafficking GTPases and coat proteins; ubiquitin–proteasome; mRNA-binding proteins; ergosterol synthesis.
- **Membership:** a protein's membership in each module is the share of its STRING weight that links into that module. A complex takes the mean over its subunits. 4,607 of the 5,638 species are modulated.
- **Fields:** each module has a smooth Gaussian random field over the cell (correlation length 500 nm, unit variance). A species' copies are distributed over tiles in proportion to region volume × exp(0.8 × Σ membership × field).
- **Crowding:** an iterative correction keeps each tile's total macromolecular volume proportional to its compartment volume, so composition varies between regions but crowding does not.
- **Effect:** in fully cytosolic tiles, module densities vary across tiles with a coefficient of variation of 0.44, against 0.10 for random placement, and form domains a few hundred nm across with roughly 0.3–3 times the mean density. Total crowding varies as little as with random placement (CV 0.042 vs 0.041).

**Meso and small scale (nm to tens of nm): partners placed next to each other.**
- **Partners:** each species gets up to 16 partners. Pool pairs are "contact" partners; STRING links that are not pool pairs are "loose" partners. Strength is scaled from the PAE (1 at ≤ 2 Å, falling to 0 at 8 Å) or the STRING score (0.2 at 700 to 1 at 1000), and divided by √(degree/20) for proteins with more than 20 links, so hubs do not dominate.
- **Attraction:** when a copy is placed, each partner p contributes a_p = strength × A × n_p / (8 + n_p). Here n_p is the number of copies of p already placed in the tile and its halo, and A is 0.6 for contact and 0.25 for loose partners. The pull saturates with the partner's local copy number, like mass action.
  - With probability min(0.6, Σ a_p), the copy is proposed next to a randomly chosen placed copy of a partner (chosen in proportion to a_p): in contact for contact partners, with a 2–15 nm gap for loose partners.
  - The copy gets up to 40 such proposals; if none passes the collision and compartment checks, it is placed at random as before.
- **Effect:** 5.8 M copies drew a partner and 4.2 M of them (about 20% of all placed molecules) were placed next to one. A molecule's contact partner lies within contact distance 6.4% of the time (random placement: 0.8%), and a loose partner within 15 nm 19.6% of the time (random: 14.6%). Placement completeness is unchanged (97.9%).

**What it is not:** strengths and bound fractions are assumptions, as no affinities are used; the partners are not bound in a specific orientation (pre-docked pairs are a separate future step); and smaller species are pulled towards larger ones that are already placed, since species are placed in order of decreasing size.

**In the viewer:** Display → colour "Interaction module (STRING)" colours each species by its dominant module. Clicking a module in the legend highlights it.

## 9. Validation

| Check | Result |
|---|---|
| Residue-bead clashes (beads of different molecules closer than 0.3 nm) | 1.9 × 10⁻⁶ of 3.1 M sampled beads in 64 test cubes; 4 between tiles |
| Containment (beads outside their compartment) | 5.6 × 10⁻⁵ of 17 M beads |
| Membrane-protein orientation vs local surface normal | median 0.0°, 100% within 1° |
| TM residues within 2 nm of the membrane mid-surface | 97.3% (median 0.73 nm) |
| ATP synthase dimer TM residues within 2 nm of the mid-surface | 98.8% |
| Placed / requested copies | 97.8% overall; 99.9% outside mitochondria; matrix 62%; nucleolus 99% |
| Inventory protein mass vs proteome | within 0.4% before packing |
| Organisation layer (section 8.4), against a random repack | contact partners in contact 6.4% vs 0.8%; loose partners within 15 nm 19.6% vs 14.6%; module density CV 0.44 vs 0.10; total crowding CV 0.042 vs 0.041 |

## 10. Densities, mass and volume fractions

### 10.1 Macromolecular densities

Densities are the mass of placed macromolecules per compartment volume. They count protein plus the RNA of ribosomes, tRNAs and complexes, the DNA of nucleosomes, and mRNA. Lipids, polysaccharides and small molecules are not counted.

| Compartment | Protein (mg/mL) | RNA/DNA (mg/mL) | Total (mg/mL) | Dry volume fraction |
|---|---|---|---|---|
| Cytosol | 70 | 23 | 93 | 6.4% |
| Nucleoplasm | 50 | 5 | 56 | 4.0% |
| Nucleolus | 79 | 6 | 85 | 6.1% |
| Mitochondrial matrix | 107 | 4 | 111 | 8.0% |
| Mitochondrial intermembrane space | 27 | 0 | 27 | 2.0% |
| ER lumen | 28 | 0 | 28 | 2.1% |
| Golgi / endosomes | 64 | 0 | 64 | 4.7% |
| Peroxisomes | 83 | 0 | 83 | 6.1% |
| Vacuole | 4 | 0 | 4 | 0.3% |

- **Dry volume fraction:** computed with partial specific volumes of 0.73 mL/g for protein and 0.55 mL/g for nucleic acid.
- **Hydrated volume fraction:** with a hydration shell of about 0.3–0.4 g water per g, the hydrated volume fraction is roughly 1.5× the dry one: about 10% in the cytosol.
- **Excluded volume:** the volume excluded to a diffusing particle is larger still, and grows with the particle's size. It can be computed from the coordinates as accessible volume versus probe radius; this has not been done yet.

### 10.2 Whole-cell totals

| Quantity | Value |
|---|---|
| Protein in the inventory | 2.7 pg |
| Protein placed | 2.2 pg |
| RNA and DNA placed | 0.62 pg |
| mRNA | 0.02 pg |

Placed protein is lower than the inventory for four reasons:
- residues not resolved in PDB entries;
- AlphaFold residues trimmed for pLDDT < 50 (3.4% of mass);
- mitochondrial copies that did not fit;
- dropped and structure-less copies.

With full-length sequences, cytosolic protein would be about 85–90 mg/mL.

### 10.3 Comparison with dry-mass measurements

Measured cytoplasmic dry-mass densities are around 200 mg/mL, more than twice the macromolecular density here. These are rough estimates of the missing components, averaged over the cytoplasm (cytosol plus organelles, excluding the vacuole):

| Component | Estimate (mg/mL) | Basis |
|---|---|---|
| Macromolecules as placed | about 90 | This model (cytosol 93, nucleus about 60) |
| Full-length protein correction | about +15 | Section 10.2 |
| Lipids of internal membranes and lipid droplets | about 20–25 | About 155 µm² of internal membrane plus 0.14 µm³ of lipid droplets in this model, about 0.8 pg of lipid |
| Small molecules and ions | probably 30–40 | Free amino acids, nucleotides, glutathione, sugar phosphates, phosphate and K⁺, a few hundred mM in total |
| Glycogen and trehalose | low in log phase on glucose | Can exceed 10% of dry weight after the diauxic shift |

Together these give about 150–160 mg/mL. The remainder is within the uncertainty of the absolute proteome size (section 2.1). A higher protein total than 4.2 × 10⁷ molecules would close most of the gap.

## 11. Known limitations

| Area | Limitation |
|---|---|
| Physics | Static random sequential addition: no physical interactions or dynamics, random orientations, and densities below the jamming limit of real cytoplasm. The optional motion is illustrative (section 7). |
| Organisation | Partner proximity and module domains (section 8.4) follow assumed strengths, not measurements. Partners are near each other, not bound in specific orientations. |
| Abundances | About 2-fold per-protein uncertainty. The proteome total is probably low. Ribosomes are probably low (169,000 vs about 3 × 10⁵ estimated). Fixed counts for tRNA and mRNA are assumptions. |
| Geometry | A composite of two segmented cells; vacuole thickness along z is reconstructed (the segmentation flattens vacuoles). Nucleus smaller than in EM. A single procedural mitochondrial reticulum. No bud, so no cell-cycle-dependent structures. |
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
| AlphaFold 3 pool predictions (in-house, unpublished) | Conformational ensembles behind the optional motion; pair scores for partner placement | not public |
| STRING v12 | Functional and physical associations, for interaction modules and loose partners | CC-BY 4.0 |
| EMPIAR-11416 (cryo-plasma FIB/SEM of yeast), segmented into cell, vacuole, lipid-droplet and nucleus meshes | The cell outline, vacuoles and lipid droplets (cell 08) and the nucleus (cell 05) of the model cell (section 3.3) | CC0 |

### 12.2 Literature behind the parameters

| Reference | Used for |
|---|---|
| Ho B, Baryshnikova A, Brown GW (2018). *Cell Systems* 6:192–205 | Total protein molecules per cell |
| Jorgensen P *et al.* (2002). *Science* 297:395–400 | Cell volume |
| Uchida M *et al.* (2011). *Yeast* 28:227–236 | Organelle volume fractions (soft X-ray tomography) |
| Dumoux M *et al.* (2023). *eLife* 12:e83623; EMPIAR-11416 | Cryo-plasma FIB/SEM volume of yeast; segmented cells used for nucleus, vacuole and lipid-droplet characteristics |
| Jorgensen P *et al.* (2007). *Mol. Biol. Cell* 18:3523–3532 | Nuclear-to-cell volume ratio |
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

- **Pipeline.** The model is built by a Python pipeline (numba-accelerated) in five stages: inventory, geometry, packing, validation and encoding. A full rebuild takes roughly 10–15 minutes on a recent multi-core desktop or laptop, with fixed seeds throughout. The real shapes come from \`pipeline/cell_meshes.npz\` (in the repository). The organisation layer of the packing (section 8.4) needs STRING (public) and the pool pair scores (not public; the derived partner lists are reproducible only with them). An optional sixth stage fits the motion modes (section 7); it also needs the pool predictions. The shipped modes are in \`data/library/modes.bin\`.
- **Source code.** The pipeline and the viewer source are in \`source/\`; \`source/README.md\` explains how to rebuild.
- **This repository** holds the built viewer and the encoded data.
  - The tiles under \`data/tiles/\` store every instance as a species index, a position quantised to about 0.0055 nm and an orientation quantised to about 0.004°.
  - The residue-bead library is in \`data/library/\`.
- **Full-atom coordinates.** The viewer fetches full-atom models live from RCSB and AlphaFold DB for display in Mol*.
`,y=i=>document.getElementById(i),ce=i=>i>=1e6?(i/1e6).toFixed(i>=1e7?1:2)+" M":i>=1e4?(i/1e3).toFixed(0)+" k":i.toLocaleString("en-US"),Yt=i=>i>=1e3?(i/1e3).toFixed(i>=1e4?0:2).replace(/\.?0+$/,"")+" µm":(i>=10?i.toFixed(0):i.toFixed(1))+" nm";let j,qe,be,Xt,ye,Ke,D,N;const P=new Rn,R=y("gpu"),v={cutOn:!0,plane:[0,0,1,0],slab:60,colorMode:0,regionMask:16777215,freeHide:0,highlight:!1,selId:0,shaded:!0,organelle:!1,motion:0,l0px:1.6,l1px:2.5,minpx:.3,outline:1,showMeshes:!0,showLipids:!0,renderScale:1,memColors:!1};let we=!1,S=!0,kt=0,pt=0,Qt="auto",K;const et={actin:["F_actin","Arp23","CP_cap","TPM1_coil"],microtubules:["MT_tubulin","gTuRC_nuc","gTuRC_cyto"],ribosomes:["ribosome_80S","ribosome_80S_Ssh1","ribosome_80S_Sec61"],tRNA:["tRNA"],nucleosomes:["nucleosome"],pores:["NPC"]},jn=[...et.actin,...et.microtubules,"NPC","nucleosome"];let oe=-1,Jt=null,U=null,se=null;const ie=[];function Xe(i,e){y("loading").querySelector(".m").textContent=i,y("loading").querySelector(".pb div").style.width=`${Math.round(e*100)}%`}function Nt(i){const e=y("loading").querySelector(".err");e.style.display="block",e.innerHTML=i,y("loading").querySelector(".t").textContent="Could not start the viewer"}function Ee(i,e=2200){const t=y("toast");t.textContent=i,t.style.display="block",clearTimeout(Ee._t),Ee._t=setTimeout(()=>t.style.display="none",e)}async function Hn(){if(!("gpu"in navigator)){Nt(`This browser does not expose <b>WebGPU</b>. Use a current Chrome or Edge (version 113 or later), or Safari 26 or later.
      On Linux Chrome, enable <code>chrome://flags/#enable-unsafe-webgpu</code> and <code>#enable-vulkan</code>.
      Firefox 141+ supports WebGPU on Windows.<br><br>Also check that the page is served over <code>http://localhost</code>
      (for example <code>python -m http.server</code>). Opening the file directly with <code>file://</code> does not work.`);return}try{Xe("Loading cell geometry",.01),[be,ye]=await Promise.all([Ce("geometry.json"),Ce("bookmarks.json")]),be.F=await On(be,Me),Xt=new qn(be),qe=await In(),j=await Ln(Xe);const i=await Me("types_tiles.bin"),e=j.types.length;Ke={off:new Uint32Array(i,0,e+1),tiles:new Uint16Array(i,(e+1)*4)},Xe("Starting WebGPU",.7),D=new Bn(R),D.onLost=r=>{Ee(r,1e4)},await D.init(j,be),K=new Uint32Array(e);for(const r of Object.values(et))for(const a of r){const d=j.byId.get(a);d!==void 0&&(K[d]|=4)}for(const r of jn){const a=j.byId.get(r);a!==void 0&&(K[a]|=8)}D.device.queue.writeBuffer(D.typeFlags,0,K),Xe("Compiling shaders and starting tile streaming",.9);let t;if(window.__YC_EMBED__){t={};for(const[r,a]of Object.entries(window.__YC_EMBED__))(r.startsWith("tiles/")||r.startsWith("library/"))&&(t[r]=$t(a))}const n=Math.max(2,Math.min(4,(navigator.hardwareConcurrency||4)-1));N=new Gn(qe,D,j,be,n,t),N.onChange=()=>{S=!0},ei();const s=new URLSearchParams(location.search).get("view"),o=s===null?0:isNaN(Number(s))?Math.max(0,ye.findIndex(r=>r.name.toLowerCase().startsWith(s.toLowerCase()))):Math.min(ye.length-1,Math.max(0,Number(s)));Le(o,!1),y("loading").style.display="none",requestAnimationFrame(dn),window.yc={R:D,T:N,cam:P,fs:v,lib:j,idx:qe,applyBookmark:Le,setType:Ie,setModule:rt,pickAt:it,bookmarks:ye,setQuality:on,describe:Ct,ready:!0,state:()=>({pending:N.pendingForView(),detailPending:N.detInflight,loaded:N.loaded.size,totals:N.totals(),stats:D.stats,cam:{t:P.target,d:P.dist},plane:v.plane,gpuMs:an,gpuBusy:ot,dirty:S}),render:()=>{S=!0}}}catch(i){console.error(i),Nt(`${String(i?.message||i)}<br><br>See the browser console for details.`)}}function he(){return[v.plane[0],v.plane[1],v.plane[2]]}function re(i,e){i=J(i),v.plane=[i[0],i[1],i[2],e],pe(),S=!0}function Se(){const i=P.fwd;re(i,V(i,P.target))}let de=.5;try{const i=Number(localStorage.getItem("yc.nudgeStep"));i>0&&(de=i)}catch{}function tt(i){const e=he();P.target=_(P.target,L(e,i)),re(e,v.plane[3]+i)}function Yn(i){const e=R.clientWidth/R.clientHeight,t=J(i.normal);let n=J(ue(i.up,L(t,V(i.up,t))));return{target:[...i.target],dist:P.distForWidth(i.width,e),fwd:t,up:n}}function Le(i,e=!0,t){const n=ye[i];if(!n)return;(v.organelle||v.memColors)&&i!==0&&Pe("molecules");const s=Yn(n),o=J(n.normal),r=V(o,n.target),a=he(),d=v.plane[3],l=v.slab;if(v.cutOn=!0,y("cutOn").checked=!0,we=!1,y("followView").classList.remove("on"),!e){P.set(s),re(o,r),v.slab=n.slab,pe(),S=!0,t?.();return}P.flyTo(s,2600,()=>{re(o,r),v.slab=n.slab,pe(),t?.()},c=>{let u=J(_(L(a,1-c),L(o,c)));Ae(_(L(a,1-c),L(o,c)))<.001&&(u=o),v.plane=[u[0],u[1],u[2],d+(r-d)*c],v.slab=Math.exp(Math.log(l)+(Math.log(n.slab)-Math.log(l))*c),pe()}),nt(n),S=!0}function nt(i){const e=y("caption");if(!i){e.style.display="none";return}e.innerHTML=`<b>${i.name}</b>${i.caption}`,e.style.display="block",clearTimeout(nt._t),U||(nt._t=setTimeout(()=>e.style.display="none",9e3))}function Kt(){U={i:0,phase:"fly",t:performance.now()},y("tourBtn").textContent="■ Stop",y("tourBtn").classList.add("on"),Le(0,!0,()=>{U&&(U.phase="wait",U.t=performance.now())})}function xe(){U=null,y("tourBtn").textContent="▶ Tour",y("tourBtn").classList.remove("on"),nt(null)}function Xn(i){if(U){if(U.phase==="wait"&&(N.pendingForView()===0||i-U.t>14e3))U.phase="hold",U.t=i;else if(U.phase==="hold"&&i-U.t>6500){U.i=(U.i+1)%ye.length,U.phase="fly";const e=U.i;Le(e,!0,()=>{U&&U.i===e&&(U.phase="wait",U.t=performance.now())})}}}let Ve=-1;function rt(i){oe>=0&&(K[oe]&=-2,oe=-1,Tt()),Ve=i;for(let e=0;e<K.length;e++)i>=0&&j.moduleOf[e]===i?K[e]|=1:K[e]&=-2;v.highlight=i>=0,D.device.queue.writeBuffer(D.typeFlags,0,K),vt(),S=!0}function Ie(i){Ve>=0&&rt(-1),oe>=0&&(K[oe]&=-2),oe=i,i>=0&&(K[i]|=1),v.highlight=i>=0,D.device.queue.writeBuffer(D.typeFlags,0,K),Tt(),S=!0}function Tt(){const i=y("selType");if(oe<0){i.style.display="none";return}const e=j.types[oe],t=N.countLoaded(oe);i.style.display="block",i.innerHTML=`<div><b>${e.gene}</b> ${e.name}</div>
    <div class="sub" style="margin:3px 0 6px">${ce(e.copies)} copies in the cell · ${ce(t)} in loaded tiles · ${en(e)}</div>
    <div class="grid2"><button class="small" id="flyNearest">Fly to nearest</button><button class="small" id="clearType">Clear highlight</button></div>`,y("flyNearest").onclick=()=>Zt(oe),y("clearType").onclick=()=>Ie(-1)}function Qn(i){if(i=i.trim().toLowerCase(),!i)return[];const e=[];return j.types.forEach((t,n)=>{const s=t.gene.toLowerCase(),o=t.id.toLowerCase(),r=(t.uniprot||"").toLowerCase(),a=(t.name||"").toLowerCase(),d=(t.genes||"").toLowerCase();let l=0;s===i||o===i||r===i?l=100:s.startsWith(i)||o.startsWith(i)?l=60:d.split(/[;,\s]+/).includes(i)?l=55:a.includes(i)?l=30+(a.startsWith(i)?10:0):(d.includes(i)||r.startsWith(i)||(t.pdb||"").toLowerCase()===i)&&(l=20),l&&e.push({i:n,s:l+Math.log10(t.copies+1)})}),e.sort((t,n)=>n.s-t.s),e.slice(0,14).map(t=>t.i)}function Zt(i){const e=Ke.off[i],t=Ke.off[i+1];if(t<=e){Ee("No copies of this molecule were placed.");return}let n=-1,s=1e18;for(let u=e;u<t;u++){const m=Ke.tiles[u],f=N.tileCenter(m),T=Ae(ue(f,P.target));T<s&&(s=T,n=m)}const o=N.tileCenter(n),r=j.types[i],a=Math.max(90,r.br*9),d=u=>{const m=N.nearestInTile(u,i,o);if(m<0)return;const f=N.decode(u,m),T=he(),w=P.state();w.target=f.p,w.dist=P.distForWidth(a,R.clientWidth/R.clientHeight),P.flyTo(w,900,()=>{re(T,V(T,f.p)-r.br*.4)}),v.plane[3]=V(T,f.p)-r.br*.4;const M=N.poolIndex(u,m);st({kind:"inst",id:(Fn|M)>>>0,info:N.instanceAt(M)})};N.loaded.get(n);const l=he(),c=P.state();c.target=o,c.dist=P.distForWidth(a*2.5,R.clientWidth/R.clientHeight),v.cutOn=!0,y("cutOn").checked=!0,(v.organelle||v.memColors)&&Pe("molecules"),P.flyTo(c,2200,()=>{re(l,V(l,o)-20);const u=N.loaded.get(n);if(u)d(u);else{const m=f=>{f.idx===n&&(N.onTileLoaded.splice(N.onTileLoaded.indexOf(m),1),d(f))};N.onTileLoaded.push(m)}}),re(l,V(l,o)-20)}async function it(i,e){const t=D.width/R.clientWidth;return D.pick(i*t,e*t)}function en(i){return i.src==="PDB"?"PDB "+i.pdb:i.src==="AF-multimer"?"AlphaFold complex":"AlphaFold"}function Ct(i){const e=i>>>30;if(i===0)return{kind:"none",title:"Outside the cell"};if(e===1){const t=N.instanceAt(i&1073741823);return t?{kind:"inst",id:i,info:t,title:j.types[t.type].gene,sub:$e[t.region]}:{kind:"none",title:""}}if(e===2){const t=i>>>1&15,n=i&1,s=St[t]||"membrane";return{kind:"lipid",id:i,mem:t,leaf:n,title:"Lipid bilayer",sub:`${s} · ${t===8?"phospholipid monolayer":n?"cytoplasmic leaflet":t===0?"exoplasmic (outer) leaflet":"lumenal leaflet"}`}}if(e===3){const t=i>>>15&7,n=i&32767,s=["mRNA","Linker DNA","β-1,3-glucan fibril","Mannan chain","pre-rRNA","FG repeats"][t]??"polymer",o=["nucleotide-resolution thread (schematic path)","between consecutive nucleosomes, 2 bp per bead","procedural, 36 glucose per bead","procedural α-mannan on a cell-wall mannoprotein","nascent ribosomal RNA","disordered domain of an FG nucleoporin (NSP1, NUP49, NUP57, NUP159), one bead per residue"][t]??"";return{kind:"free",id:i,fk:t,chain:n,title:s,sub:o}}return i>=512&&i<528?{kind:"band",id:i,title:"Membrane (cut face)",sub:St[i-512]}:i>=256&&i<272?{kind:"label",id:i,title:En[i-256],sub:"compartment (not loaded at this distance, or beyond the slab)"}:{kind:"none",title:""}}function st(i){Jt=i,v.selId=i&&i.kind==="inst"?i.id>>>0:0,tn(),S=!0}function tn(){const i=y("inspect"),e=Jt;if(!e||e.kind==="none"){i.style.display="none";return}i.style.display="block";let t='<button class="close" id="inspClose">×</button>';if(e.kind==="inst"){const n=j.types[e.info.type],s=n.src==="PDB",o=n.nb,r=n.nnt;t+=`<h2>${n.gene}</h2><div>${n.name}</div><div style="margin-top:6px"><span class="badge ${s?"pdb":"af"}">${en(n)}</span>
      ${n.oligomer>1?`<span class="badge" style="background:#666">${n.oligomer}-mer</span>`:""}</div>
      <div class="kv">
      ${n.genes&&n.genes!==n.gene?`<div>Subunits</div><div style="word-break:break-word">${n.genes.split(/;\s*/).slice(0,40).join(", ")}${n.genes.split(";").length>40?" …":""}</div>`:""}
      ${n.uniprot?`<div>UniProt</div><div>${n.uniprot}</div>`:""}
      <div>Location</div><div>${$e[e.info.region]}</div>
      <div>Copies / cell</div><div>${n.copies.toLocaleString("en-US")}</div>
      <div>Mass</div><div>${(n.mass/1e3).toFixed(1)} kDa</div>
      <div>Residues drawn</div><div>${(o-r).toLocaleString("en-US")} aa${r?` + ${r.toLocaleString("en-US")} nt`:""}</div>
      <div>Radius</div><div>${n.br.toFixed(1)} nm (Rg ${n.rg.toFixed(1)} nm)</div>
      ${!s&&n.plddt_mean!=null?`<div>Mean pLDDT</div><div>${n.plddt_mean.toFixed(1)} (${Math.round((n.plddt_frac_kept??1)*100)}% of residues kept, pLDDT ≥ 50)</div>`:""}
      <div>Position</div><div>${e.info.p.map(a=>a.toFixed(1)).join(", ")} nm</div>
      </div>
      <div class="grid2"><button class="small" id="loadAtoms">Load atomic model</button><button class="small" id="hlType">Highlight all</button>
      <button class="small" id="centerOn">Center view</button><button class="small" id="hideType">${K[e.info.type]&2?"Show":"Hide"} this type</button></div>
      <div style="margin-top:9px">${n.uniprot?`<a href="https://www.uniprot.org/uniprotkb/${n.uniprot}" target="_blank">UniProt</a>`:""}
      ${s?`<a href="https://www.rcsb.org/structure/${n.pdb}" target="_blank">RCSB ${n.pdb}</a>`:""}
      ${(n.model||n.uniprot)&&!s?`<a href="https://alphafold.ebi.ac.uk/entry/${n.model||n.uniprot}" target="_blank">AlphaFold DB</a>`:""}
      ${/^[A-Z0-9-]+$/.test(n.gene)&&!Vt[n.id]&&!n.id.startsWith("CPX")?`<a href="https://www.yeastgenome.org/locus/${n.gene}" target="_blank">SGD</a>`:""}</div>`}else t+=`<h2>${e.title}</h2><div class="sub" style="font-size:12.5px">${e.sub||""}</div>`,e.kind==="lipid"&&(t+='<p class="sub">Lipids are generated in your browser on the mesh of this membrane. Classes are PI, PE, PC, PA and PS, with ergosterol drawn as two beads. They sit at 0.65 nm² per lipid and are removed where they would overlap protein.</p>');i.innerHTML=t,y("inspClose").onclick=()=>st(null),e.kind==="inst"&&(y("loadAtoms").onclick=()=>Kn(e.info.type),y("hlType").onclick=()=>Ie(e.info.type),y("centerOn").onclick=()=>{const n=P.state();n.target=e.info.p,P.flyTo(n,700)},y("hideType").onclick=()=>{K[e.info.type]^=2,D.device.queue.writeBuffer(D.typeFlags,0,K),tn(),S=!0})}let Qe=null,Oe=null;function Jn(){if(Oe)return Oe;const i="https://cdn.jsdelivr.net/npm/molstar@4.18.0/build/viewer/",e=t=>new Promise((n,s)=>{const o=document.createElement("link");o.rel="stylesheet",o.href=t+"molstar.css",document.head.appendChild(o);const r=document.createElement("script");r.src=t+"molstar.js",r.onload=()=>n(),r.onerror=()=>{o.remove(),r.remove(),s(new Error("Mol* failed to load"))},document.head.appendChild(r)});return Oe=(window.__YC_EMBED__?e(i):e("molstar/").catch(()=>e(i))).catch(t=>{throw Oe=null,t}),Oe}async function Kn(i){const e=j.types[i];y("molModal").style.display="flex",y("molTitle").textContent=`${e.gene} · ${e.name} · loading…`;try{await Jn();const t=window.molstar;Qe?await Qe.plugin.clear():Qe=await t.Viewer.create("molHost",{layoutIsExpanded:!1,layoutShowControls:!1,layoutShowRemoteState:!1,layoutShowSequence:!0,layoutShowLog:!1,layoutShowLeftPanel:!1,viewportShowExpand:!1,viewportShowSelectionMode:!1,viewportShowAnimation:!1,collapseLeftPanel:!0,pdbProvider:"rcsb",emdbProvider:"rcsb"});let n="",s="";if(e.src==="PDB"&&e.pdb){n=`https://files.rcsb.org/download/${e.pdb}-assembly1.cif`,s=`https://www.rcsb.org/structure/${e.pdb}`;try{(await fetch(n,{method:"HEAD"})).ok||(n=`https://files.rcsb.org/download/${e.pdb}.cif`)}catch{n=`https://files.rcsb.org/download/${e.pdb}.cif`}}else{const r=e.model||e.uniprot;s=`https://alphafold.ebi.ac.uk/entry/${r}`;try{const d=await(await fetch(`https://alphafold.ebi.ac.uk/api/prediction/${r}`)).json();n=(Array.isArray(d)?d[0]:d)?.cifUrl}catch{}n||(n=e.model?`https://alphafold.ebi.ac.uk/files/${e.model}-model_v1.cif`:`https://alphafold.ebi.ac.uk/files/AF-${e.uniprot}-F1-model_v6.cif`)}await Qe.loadStructureFromUrl(n,"mmcif",!1),y("molTitle").textContent=`${e.gene} · ${e.name} · ${e.src==="PDB"?"PDB "+e.pdb+" (biological assembly)":e.model?`AlphaFold DB predicted ${e.oligomer}-mer (full length, incl. low-pLDDT regions)`:"AlphaFold DB model (full length, incl. low-pLDDT regions)"}`;const o=y("molLink");o.href=s,o.textContent="Open source page"}catch(t){y("molTitle").textContent=`Could not load the atomic model: ${t?.message||t}. Internet access to RCSB/AlphaFold DB is needed.`}}let Dt=!1;function nn(){Dt||(Dt=!0,Ee("Loading organelle surfaces (6.7 MB)…"),Me("organelle_meshes.bin.gz").then(i=>{D.uploadOrganelleMeshes(i),S=!0}))}let wt=()=>{};function sn(){const i=v.memColors&&!v.organelle;y("modeSection").classList.toggle("on",!v.organelle&&!i),y("modeOrg").classList.toggle("on",v.organelle),y("modeMem").classList.toggle("on",i)}function Zn(i){v.organelle=i,i&&nn(),sn(),S=!0}function Pe(i){i!=="membranes"&&v.memColors&&wt(!1),Zn(i==="organelles"),i==="membranes"&&wt(!0)}function on(i){Qt=i;const e={low:[2.6,.75],med:[1.6,1],high:[1,1],ultra:[.6,1],auto:[1.6,1]},[t,n]=e[i];v.l0px=t,v.renderScale=n,D.resize(n),S=!0}function vt(){const i=y("legend"),e=(n,s)=>`<div><span class="sw" style="background:${n}"></span>${s}</div>`,t=v.colorMode;t===0?i.innerHTML=[e(G[0],"cytosolic protein"),e(G[1],"nuclear protein"),e("#EE9BC0","ribosome protein"),e("#C9558F","rRNA"),e("#E27FA8","tRNA"),e("#B0437A","mRNA"),e("#F4B26A","histone"),e("#E07636","DNA"),e(G[23],"nuclear pore"),e("#7CC6C2","FG repeats"),e("#E0893A","actin filament"),e("#D9534F","Arp2/3 complex"),e("#8E6CB0","capping protein"),e("#2FA4A0","tropomyosin"),e("#74C476","microtubule"),e("#2E8B57","γ-tubulin ring"),e(G[14],"membrane protein"),e(G[5],"mitochondrial"),e(G[3],"ER lumen"),e("#F3C643","lipid head"),e("#F8E3A0","lipid tail"),e("#CDA866","β-glucan"),e("#C3C27A","mannan"),e(G[9],"wall protein"),e("#8B88C9","large complex")].join(""):t===1?i.innerHTML=e("#4C8CD9","experimental (PDB)")+e("#F5A84C","AlphaFold model"):t===2?i.innerHTML=e("#0053D6","pLDDT > 90")+e("#65CBF3","70–90")+e("#FFDB13","50–70")+e("#B8B8B8","experimental (PDB)"):t===3?i.innerHTML='<div style="grid-column: 1 / 3"><div style="height:8px;border-radius:3px;background:linear-gradient(90deg,#440154,#3b528b,#21918c,#5ec962,#fde725)"></div><div style="display:flex;justify-content:space-between"><span>1</span><span>10³</span><span>10⁶ copies</span></div></div>':t===6?(i.innerHTML=j.modules.length?'<div class="sub" style="grid-column: 1 / 3">Click a module to highlight it</div>'+j.modules.map((n,s)=>`<div class="modsw${Ve===s?" on":""}" data-m="${s}" title="${n.n} proteins"><span class="sw" style="background:${Wt(s)}"></span>${n.label}</div>`).join("")+e("#D1D1CC","no module"):'<div style="grid-column: 1 / 3">No interaction modules in this build</div>',i.querySelectorAll(".modsw").forEach(n=>n.onclick=()=>{const s=Number(n.dataset.m);rt(Ve===s?-1:s)})):t===5?i.innerHTML=[0,1,2,3,4,5,6,7,8,9,14,15,18,19].map(n=>e(G[n],$e[n])).join(""):i.innerHTML=""}function pe(){const i=v.plane[3];y("planeD").value=String(Math.round(i)),y("planeDv").textContent=`${Number(i.toFixed(de<.5?2:1))} nm`;const e=v.slab;y("slab").value=String(Math.round(Math.log(e/5)/Math.log(4400/5)*1e3)),y("slabv").textContent=Yt(e),["axX","axY","axZ"].forEach((t,n)=>y(t).classList.toggle("on",Math.abs(v.plane[n])>.999)),xt=!0}function ei(){const i=j.types.length;let e=0;for(const h of j.types)e+=h.copies;if(y("subCount").textContent=`${(e/1e6).toFixed(1)} M macromolecules of ${i.toLocaleString("en-US")} species, each a real PDB or AlphaFold structure`,window.__YC_PREVIEW__){const h=qe.nInst.reduce((g,b)=>g+b,0);y("subCount").innerHTML=`<b>Preview:</b> ${qe.n} of 6,416 tiles (${(h/1e3).toFixed(0)} k molecules) from the full ${(e/1e6).toFixed(1)} M-molecule cell. The full cell is in yeast_cell_viewer.zip.`,y("modeOrg").disabled=!0,y("modeOrg").title="Organelle surfaces are only in the full app",y("modeMem").disabled=!0,y("modeMem").title="Membrane surfaces are only in the full app"}const t=y("views");ye.forEach((h,g)=>{const b=document.createElement("button");b.textContent=`${(g+1)%10===0&&g<10?0:g<9?g+1:""}${g<10?" · ":""}${h.name}`,b.title=h.caption,b.onclick=()=>{xe(),Le(g)},t.appendChild(b)}),y("tourBtn").onclick=()=>U?xe():Kt();const n=y("search"),s=y("results");let o=0,r=[];const a=()=>{s.innerHTML=r.map((h,g)=>{const b=j.types[h];return`<div data-i="${h}" class="${g===o?"act":""}"><span class="nm"><b>${b.gene}</b> ${b.name}</span><span class="n">${ce(b.copies)}</span></div>`}).join(""),s.querySelectorAll("div").forEach(h=>h.addEventListener("click",()=>{Ie(Number(h.dataset.i)),r=[],a(),n.value=""}))};n.addEventListener("input",()=>{r=Qn(n.value),o=0,a()}),n.addEventListener("keydown",h=>{h.key==="ArrowDown"?(o=Math.min(r.length-1,o+1),a(),h.preventDefault()):h.key==="ArrowUp"?(o=Math.max(0,o-1),a(),h.preventDefault()):h.key==="Enter"&&r.length?(Ie(r[o]),Zt(r[o]),r=[],a(),n.value="",n.blur()):h.key==="Escape"&&(r=[],a(),n.blur()),h.stopPropagation()}),y("cutOn").onchange=h=>{v.cutOn=h.target.checked,S=!0};const d=h=>()=>{const g=[0,0,0];g[h]=V(P.fwd,[h===0?1:0,h===1?1:0,h===2?1:0])>=0?1:-1,v.cutOn=!0,y("cutOn").checked=!0,re(g,V(g,P.target))};y("axX").onclick=d(0),y("axY").onclick=d(1),y("axZ").onclick=d(2),y("alignView").onclick=()=>{v.cutOn=!0,y("cutOn").checked=!0,Se()},y("followView").onclick=()=>{we=!we,y("followView").classList.toggle("on",we),we&&Se()},y("planeD").oninput=h=>{const g=Number(h.target.value),b=he(),A=g-v.plane[3];P.target=_(P.target,L(b,A)),re(b,g)};const l=y("nudgeStep");[...l.options].some(h=>Number(h.value)===de)&&(l.value=String(de)),l.onchange=()=>{de=Number(l.value),pe(),l.blur();try{localStorage.setItem("yc.nudgeStep",String(de))}catch{}};for(const[h,g]of[["nudgeBack",-1],["nudgeFwd",1]]){const b=y(h);let A=0,C=0;const F=()=>{clearTimeout(A),clearInterval(C)};b.addEventListener("pointerdown",O=>{O.preventDefault(),F(),se=null,y("sweepBtn").classList.remove("on"),tt(g*de),A=window.setTimeout(()=>{C=window.setInterval(()=>tt(g*de),60)},400)});for(const O of["pointerup","pointerleave","pointercancel"])b.addEventListener(O,F)}y("slab").oninput=h=>{const g=Number(h.target.value)/1e3;v.slab=5*Math.pow(4400/5,g),pe(),S=!0},y("sweepBtn").onclick=()=>{if(se){se=null,y("sweepBtn").classList.remove("on");return}const h=he();v.cutOn=!0,y("cutOn").checked=!0;const g=be.G[0];se={t0:performance.now(),d0:-g,d1:g,dur:3e4},P.target=_(P.target,L(h,-g-v.plane[3])),y("sweepBtn").classList.add("on")},y("flipBtn").onclick=()=>{const h=he(),g=P.state(),b=P.target,A=L(h,-1);re(A,-v.plane[3]),g.fwd=L(g.fwd,-1),g.target=b,P.flyTo(g,1200)},y("modeSection").onclick=()=>Pe("molecules"),y("modeOrg").onclick=()=>Pe("organelles"),y("modeMem").onclick=()=>Pe("membranes"),y("colorMode").onchange=h=>{v.colorMode=Number(h.target.value),v.colorMode!==6&&Ve>=0&&rt(-1),vt(),S=!0},y("shaded").onchange=h=>{v.shaded=h.target.checked,S=!0},y("motion").onchange=h=>{v.motion=Number(h.target.value),h.target.blur(),S=!0},y("outlines").onchange=h=>{v.outline=h.target.checked?1:0,S=!0},vt();const c=y("showList"),u=h=>g=>ti(et[h],!g),m=[["Actin","#E0893A",u("actin")],["Microtubules","#74C476",u("microtubules")],["Ribosomes","#EE9BC0",u("ribosomes")],["tRNA","#E27FA8",u("tRNA")],["Nucleosomes","#F4B26A",u("nucleosomes")],["Nuclear pores",G[23],u("pores")],["mRNA","#B0437A",h=>_e(0,h)],["Linker DNA","#E07636",h=>_e(1,h)],["FG repeats","#7CC6C2",h=>_e(5,h)],["β-Glucan","#CDA866",h=>_e(2,h)],["Mannan","#C3C27A",h=>_e(3,h)],["Lipids","#F3C643",h=>{v.showLipids=h,S=!0}],["Membranes","#F8E3A0",h=>{v.showMeshes=h,S=!0}],["Soluble proteins",G[0],h=>Je([0,1,2,3,4,5,6,7,8,10,11,12,13],h)],["Membrane proteins",G[14],h=>Je([14,15,16,17,18,19,20,21],h)]],f=[],T=$e.map((h,g)=>g),w=m.findIndex(([h])=>h==="Membranes"),M=()=>{v.memColors=f.every((h,g)=>h.checked===(g===w))&&!window.__YC_PREVIEW__,v.memColors&&nn(),sn(),S=!0},p=h=>{m.forEach(([,,b],A)=>{f[A].checked=h(A),b(h(A))});const g=m.every((b,A)=>h(A));Je(T,g),M()};m.forEach(([h,g,b],A)=>{const C=document.createElement("label");C.className="chk",C.innerHTML=`<input type="checkbox" checked><span class="sw" style="background:${g}"></span><span class="nm">${h}</span><a class="only" title="Show only ${h.toLowerCase()}">only</a>`;const F=C.querySelector("input");f.push(F),F.onchange=()=>{b(F.checked),M()},C.querySelector(".only").onclick=O=>{O.preventDefault(),O.stopPropagation(),p(k=>k===A)},c.appendChild(C)}),wt=h=>p(h?(g=>g===w):(()=>!0)),y("showAll").onclick=()=>p(()=>!0),y("hideAll").onclick=()=>p(()=>!1);const x=y("regionList");$e.forEach((h,g)=>{const b=document.createElement("label");b.className="chk",b.innerHTML=`<input type="checkbox" checked data-r="${g}"><span class="sw" style="background:${G[g]}"></span>${h}`,b.querySelector("input").onchange=A=>Je([g],A.target.checked),x.appendChild(b)}),y("quality").onchange=h=>on(h.target.value),y("detailOn").onchange=h=>{if(N.detailEnabled=h.target.checked,!N.detailEnabled)for(const g of N.loaded.values())g.det&&N.dropDetail(g);S=!0},y("shotBtn").onclick=rn,y("helpBtn").onclick=()=>y("helpModal").style.display="flex",y("aboutBtn").onclick=()=>y("aboutModal").style.display="flex",y("aboutBody").innerHTML=$n,y("modelBtn").onclick=()=>{y("modelBody").innerHTML||(y("modelBody").innerHTML=Vn(Wn)),y("modelModal").style.display="flex"},y("aboutModal").addEventListener("click",h=>{h.target.id==="aboutToModel"&&(h.preventDefault(),y("aboutModal").style.display="none",y("modelBtn").click())}),document.querySelectorAll(".modal").forEach(h=>h.addEventListener("click",g=>{(g.target===h||g.target.classList.contains("close"))&&(h.style.display="none")})),y("molClose").onclick=()=>y("molModal").style.display="none",y("panelToggle").onclick=()=>{y("panel").classList.toggle("collapsed"),y("panelToggle").textContent=y("panel").classList.contains("collapsed")?"›":"‹"},pe(),Et(),ii()}function ti(i,e){for(const t of i){const n=j.byId.get(t);n!==void 0&&(e?K[n]|=2:K[n]&=-3)}D.device.queue.writeBuffer(D.typeFlags,0,K),S=!0}function _e(i,e){e?v.freeHide&=~(1<<i):v.freeHide|=1<<i,S=!0}function Je(i,e){for(const t of i)e?v.regionMask|=1<<t:v.regionMask&=~(1<<t);v.regionMask>>>=0,document.querySelectorAll("#regionList input").forEach(t=>{t.checked=(v.regionMask>>>Number(t.dataset.r)&1)===1}),S=!0}function Et(){let i=null;const e=new Set;R.addEventListener("contextmenu",n=>n.preventDefault()),R.addEventListener("pointerdown",n=>{R.setPointerCapture(n.pointerId),R.focus(),i={x:n.clientX,y:n.clientY,btn:n.button,moved:0,shift:n.shiftKey},R.classList.add("dragging"),U&&xe()}),R.addEventListener("pointermove",n=>{if(!i){ni(n.clientX,n.clientY);return}const s=n.clientX-i.x,o=n.clientY-i.y;if(i.x=n.clientX,i.y=n.clientY,i.moved+=Math.abs(s)+Math.abs(o),i.btn===2||i.btn===1||i.shift){const r=P.pxPerNmAtUnit(R.clientHeight)/P.dist;P.pan(s,o,r)}else P.orbit(s*.005,o*.005),we&&Se();P.anim=null,S=!0});const t=n=>{if(!i)return;const s=i.moved<4&&i.btn===0;i=null,R.classList.remove("dragging"),s&&it(n.offsetX,n.offsetY).then(o=>{const r=Ct(o.id);st(r.kind==="inst"||r.kind==="lipid"||r.kind==="free"||r.kind==="band"?r:null)})};R.addEventListener("pointerup",t),R.addEventListener("pointercancel",t),R.addEventListener("dblclick",n=>{it(n.offsetX,n.offsetY).then(s=>{if(s.depth<=0)return;const o=P.near(),r=o/s.depth,a=R.clientWidth/R.clientHeight,d=P.ray(n.offsetX/R.clientWidth,n.offsetY/R.clientHeight,a),l=r/V(d,P.fwd),c=_(P.eye,L(d,l)),u=P.state();u.target=c,u.dist=Math.max(P.minDist*4,P.dist*.55),P.flyTo(u,900)})}),R.addEventListener("wheel",n=>{n.preventDefault();const s=n.deltaY||n.deltaX,o=n.deltaMode===1?s*30:s;if(n.altKey)tt(Math.sign(o)*de*Math.min(4,Math.abs(o)/100));else if(n.shiftKey){const r=he(),a=Math.max(1,Math.min(v.slab,P.dist*.05))*Math.sign(o)*Math.min(4,Math.abs(o)/60+.5);P.target=_(P.target,L(r,a)),re(r,v.plane[3]+a)}else{const r=Math.exp(o*.0015*(n.metaKey?.125:1)),a=R.clientWidth/R.clientHeight,d=P.ray(n.offsetX/R.clientWidth,n.offsetY/R.clientHeight,a),l=P.dist,c=_(P.eye,L(d,l/V(d,P.fwd))),u=P.dist;P.zoom(r);const m=1-P.dist/u;P.target=_(P.target,L(ue(c,P.target),m))}P.anim=null,U&&xe(),S=!0},{passive:!1}),window.addEventListener("keydown",n=>{if(n.target.tagName==="INPUT"||n.target.tagName==="SELECT")return;const s=n.key.toLowerCase();if(e.add(s),s==="c")v.cutOn=!v.cutOn,y("cutOn").checked=v.cutOn,S=!0;else if(s==="v")v.cutOn=!0,y("cutOn").checked=!0,Se();else if(s==="[")v.slab=Math.max(5,v.slab/1.25),pe(),S=!0;else if(s==="]")v.slab=Math.min(4400,v.slab*1.25),pe(),S=!0;else if(s===","||s==="."||n.code==="Comma"||n.code==="Period"){const o=s==="."||n.code==="Period";tt((o?1:-1)*(n.altKey?de:Math.max(2,v.slab*.25)))}else if(s==="o")Pe(v.organelle?"molecules":"organelles");else if(s==="m")Pe(v.memColors&&!v.organelle?"molecules":"membranes");else if(s==="t")U?xe():Kt();else if(s==="p")rn();else if(s==="escape")st(null),Ie(-1);else if(/^[0-9]$/.test(s)){const o=s==="0"?9:Number(s)-1;o<ye.length&&(xe(),Le(o))}}),window.addEventListener("keyup",n=>e.delete(n.key.toLowerCase())),window.addEventListener("blur",()=>e.clear()),Et.keys=e,window.addEventListener("resize",()=>{D.resize(v.renderScale),S=!0})}let Ot=0,mt=!1;function ni(i,e){const t=performance.now();mt||t-Ot<90||(Ot=t,mt=!0,S=!0,it(i,e).then(n=>{mt=!1;const s=Ct(n.id),o=y("hover");if(!s.title||s.kind==="none"){o.style.display="none";return}let r=`<b>${s.title}</b>`;if(s.kind==="inst"){const a=j.types[s.info.type];r+=` · ${a.name} <span style="color:#aaa">· ${s.sub}</span>`}else r+=` <span style="color:#aaa">· ${s.sub}</span>`;o.innerHTML=r,o.style.display="block",o.style.left=`${i+14}px`,o.style.top=`${e+12}px`}))}R.addEventListener("pointerleave",()=>y("hover").style.display="none");let xt=!0,_t=0,fe=null;const ze={data:null};function ii(){const i=y("minimap");i.addEventListener("click",e=>{if(!fe)return;const t=i.getBoundingClientRect(),n=((e.clientX-t.left)/t.width*2-1)*fe.ext,s=(1-(e.clientY-t.top)/t.height*2)*fe.ext,o=_(fe.c,_(L(fe.u,n),L(fe.v,s))),r=P.state();r.target=o,P.flyTo(r,1200),U&&xe()})}function Mt(i){const e=y("minimap"),t=e.getContext("2d"),n=e.width,s=e.height,o=v.cutOn?he():P.fwd,r=v.cutOn?v.plane[3]:V(P.fwd,P.target);if(xt&&i-_t>250){xt=!1,_t=i;let M=ue(P.right,L(o,V(P.right,o)));Ae(M)<.001&&(M=J(ge(o,[0,0,1]))),M=J(M);let p=ge(o,M);V(p,P.up)<0&&(p=L(p,-1));const x=2300,h=L(o,r);fe={u:M,v:p,c:h,ext:x};const g=120,b=t.createImageData(g,g),A=[G[80],G[9],"#D2BE8E",G[83],G[3],G[3],G[1],G[2],G[4],G[6],G[5],"#E9C46A",G[7],G[8],"#C9C2DC"].map(C=>{const F=parseInt(C.slice(1),16);return[F>>16&255,F>>8&255,F&255]});A[3]=[196,214,235];for(let C=0;C<g;C++)for(let F=0;F<g;F++){const O=((F+.5)/g*2-1)*x,k=(1-(C+.5)/g*2)*x,E=_(h,_(L(M,O),L(p,k))),Y=Xt.label(E[0],E[1],E[2]),$=A[Y],W=4*(C*g+F);b.data[W]=$[0],b.data[W+1]=$[1],b.data[W+2]=$[2],b.data[W+3]=Y===0?0:255}ze.data=b}if(!ze.data||!fe)return;const a=Mt.tmp||(Mt.tmp=document.createElement("canvas"));a.width=ze.data.width,a.height=ze.data.height,a.getContext("2d").putImageData(ze.data,0,0),t.clearRect(0,0,n,s),t.imageSmoothingEnabled=!0,t.drawImage(a,0,0,n,s);const d=fe,l=M=>{const p=ue(M,d.c);return[(V(p,d.u)/d.ext+1)/2*n,(1-V(p,d.v)/d.ext)/2*s]},c=R.clientWidth/R.clientHeight,u=P.dist*Math.tan(P.fovy/2),m=u*c,f=[[-1,-1],[1,-1],[1,1],[-1,1]].map(([M,p])=>l(_(P.target,_(L(P.right,M*m),L(P.up,p*u)))));t.strokeStyle="#FF9400",t.lineWidth=3,t.beginPath(),f.forEach(([M,p],x)=>x?t.lineTo(M,p):t.moveTo(M,p)),t.closePath(),t.stroke();const[T,w]=l(P.target);t.fillStyle="#FF9400",t.beginPath(),t.arc(T,w,5,0,7),t.fill(),y("mmLabel").textContent=v.cutOn?`plane ${Math.round(r)} nm`:"view plane"}function rn(){const i=v.renderScale;D.resize(2*i),v.renderScale=2*i,ln(!0),R.toBlob(e=>{if(e){const t=document.createElement("a");t.href=URL.createObjectURL(e),t.download=`yeast_cell_${new Date().toISOString().replace(/[:.]/g,"-").slice(0,19)}.png`,t.click(),setTimeout(()=>URL.revokeObjectURL(t.href),5e3)}v.renderScale=i,D.resize(i),S=!0},"image/png"),Ee("Screenshot saved (2× resolution)")}function si(){const i=N.totals(),e=D.stats,t=ie.length>3?1e3/(ie.reduce((l,c)=>l+c)/ie.length):0,n=N.pendingForView();y("hud").innerHTML=`<b>${ce(i.inst)}</b> molecules in ${i.tiles} tiles${n?` · <span style="color:#FF9400">streaming ${n}</span>`:""}${N.starved?' · <span style="color:#FF9400">memory budget reached (thinner slab loads more)</span>':""}<br>
    drawn: <b>${ce(e.ninst)}</b> in view · ${ce(e.l0*64)} residue beads (≤) · ${ce(e.l1*16)} LOD1 · ${ce(e.l2)} LOD2 · ${ce(e.lip)} lipids · ${ce(e.fre)} polymer beads<br>
    ${t?`${t.toFixed(0)} fps · `:""}LOD bead ≥ ${v.l0px.toFixed(1)} px · GPU pools ${D.gpuMemoryMB().toFixed(0)} MB · ${(N.bytesLoaded/1048576).toFixed(0)} MB streamed${D.adapterInfo?" · "+D.adapterInfo:""}`;const s=D.width/R.clientWidth,o=P.pxPerNmAtUnit(D.height)/P.dist/s,r=120/o,a=Math.pow(10,Math.floor(Math.log10(r))),d=[1,2,5,10].map(l=>l*a).reduce((l,c)=>Math.abs(c-r)<Math.abs(l-r)?c:l);y("scalebar").querySelector(".bar").style.width=`${d*o}px`,y("scalebar").querySelector("span").textContent=Yt(d)}let an=0,ot=!1;function ln(i=!1){const e=performance.now();D.frame(P,v,be.G,N.meshDraws(),i||e-D.lastStatsT>400),D.lastStatsT=performance.now(),ot=!0,D.device.queue.onSubmittedWorkDone().then(()=>{an=performance.now()-e,ot=!1})}let zt=0;function cn(){N.update(P,R.clientWidth/R.clientHeight,v.cutOn,v.plane,v.slab,v.organelle||v.memColors,D.height),kt=performance.now(),(N.inflight||N.detInflight)&&(S=!0)}setInterval(()=>{N&&performance.now()-kt>150&&cn()},150);function dn(i){requestAnimationFrame(dn);const e=pt?i-pt:16;pt=i;const t=Et.keys;if(t&&t.size){const n=P.dist*.9*e/1e3;let s=[0,0,0];t.has("w")&&(s=_(s,L(P.fwd,n))),t.has("s")&&(s=_(s,L(P.fwd,-n))),t.has("a")&&(s=_(s,L(P.right,-n))),t.has("d")&&(s=_(s,L(P.right,n))),t.has("q")&&(s=_(s,L(P.up,-n))),t.has("e")&&(s=_(s,L(P.up,n))),t.has("z")&&P.roll(-e*.0015),t.has("x")&&P.roll(e*.0015),(Ae(s)>0||t.has("z")||t.has("x"))&&(P.fly(s),S=!0,we&&Se())}if(P.step(i)&&(S=!0,we&&!U&&Se()),se){const n=Math.min(1,(i-se.t0)/se.dur),s=he(),o=se.d0+(se.d1-se.d0)*n;P.target=_(P.target,L(s,o-v.plane[3])),v.plane[3]=o,pe(),n>=1&&(se=null,y("sweepBtn").classList.remove("on")),S=!0}if(Xn(i),v.motion>0&&(S=!0),i-kt>90&&cn(),S&&!ot){if(S=!1,ln(),ie.push(e),ie.length>30&&ie.shift(),Qt==="auto"&&ie.length>=20){const n=ie.reduce((s,o)=>s+o)/ie.length;n>45&&v.l0px<3.5?(v.l0px=Math.min(3.5,v.l0px*1.15),ie.length=0):n<22&&v.l0px>.9&&(v.l0px=Math.max(.9,v.l0px/1.1),ie.length=0)}setTimeout(()=>{S=S||D.statsBusy},0)}else S||(ie.length=0);i-zt>300&&(zt=i,si(),Mt(i),oe>=0&&N.inflight&&Tt())}Hn();
