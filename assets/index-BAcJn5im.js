(function(){const e=document.createElement("link").relList;if(e&&e.supports&&e.supports("modulepreload"))return;for(const s of document.querySelectorAll('link[rel="modulepreload"]'))n(s);new MutationObserver(s=>{for(const r of s)if(r.type==="childList")for(const a of r.addedNodes)a.tagName==="LINK"&&a.rel==="modulepreload"&&n(a)}).observe(document,{childList:!0,subtree:!0});function t(s){const r={};return s.integrity&&(r.integrity=s.integrity),s.referrerPolicy&&(r.referrerPolicy=s.referrerPolicy),s.crossOrigin==="use-credentials"?r.credentials="include":s.crossOrigin==="anonymous"?r.credentials="omit":r.credentials="same-origin",r}function n(s){if(s.ep)return;s.ep=!0;const r=t(s);fetch(s.href,r)}})();const Vo="modulepreload",Yo=function(o,e){return new URL(o,e).href},In={},Ko=function(e,t,n){let s=Promise.resolve();if(t&&t.length>0){let a=function(c){return Promise.all(c.map(d=>Promise.resolve(d).then(f=>({status:"fulfilled",value:f}),f=>({status:"rejected",reason:f}))))};const i=document.getElementsByTagName("link"),h=document.querySelector("meta[property=csp-nonce]"),l=h?.nonce||h?.getAttribute("nonce");s=a(t.map(c=>{if(c=Yo(c,n),c in In)return;In[c]=!0;const d=c.endsWith(".css"),f=d?'[rel="stylesheet"]':"";if(!!n)for(let m=i.length-1;m>=0;m--){const u=i[m];if(u.href===c&&(!d||u.rel==="stylesheet"))return}else if(document.querySelector(`link[href="${c}"]${f}`))return;const b=document.createElement("link");if(b.rel=d?"stylesheet":Vo,d||(b.as="script"),b.crossOrigin="",b.href=c,l&&b.setAttribute("nonce",l),document.head.appendChild(b),d)return new Promise((m,u)=>{b.addEventListener("load",m),b.addEventListener("error",()=>u(new Error(`Unable to preload CSS for ${c}`)))})}))}function r(a){const i=new Event("vite:preloadError",{cancelable:!0});if(i.payload=a,window.dispatchEvent(i),!i.defaultPrevented)throw a}return s.then(a=>{for(const i of a||[])i.status==="rejected"&&r(i.reason);return e().catch(r)})},ce=1024,de=4096,Qe=`
const PAGE: u32 = ${ce}u;
const SPAGE: u32 = ${de}u;
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
  caps: vec4u,            // x capL0 (offset of L1 items in list01), y motion amount x 1000, z unused,
                          // w membranes-only view (whole-cell surfaces, palette 66 + membrane, no cut)
  simBox: vec4f,          // simulation mode: box centre (world), half size; w = 0 off. Membrane meshes are clipped to it
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
`,ao=`
struct TypeInfo {
  beadOff: u32, beadCount: u32, clusterOff: u32, clusterCount: u32,
  stOff: u32, stCount: u32, br: f32, req: f32,
  color0: u32, color1: u32, flags: u32, plddt: f32,
  pad0: u32, pad1: u32, pad2: u32,                    // unused (were the pool-ensemble motion modes)
  modColor: u32,                                       // packed network-neighbourhood colour (0: no network profile)
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
`,Jo=Qe+ao+`
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
`,Xo=Qe+`
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
  if (kind == 4u) { return 0.7; }   // nascent pre-rRNA, 2 nt per bead
  if (kind == 5u) { return 0.45; }  // disordered protein residue (FG repeats)
  if (kind == 6u) { return 1.2; }   // chitin microfibril (~18 chains)
  if (kind == 7u) { return 0.45; }  // polyphosphate, 4 phosphate units per bead
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
`,Qo=`
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
`,ro=`
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
`,Zo=`
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
  } else if (mode == 6u) {   // network neighbourhood (STRING personalized PageRank)
    if (T.modColor != 0u) { c = unpackColor(T.modColor); } else { c = vec3f(0.82, 0.82, 0.80); }
    if (kind == 1u) { c = c * 0.8; }
  } else if (mode >= 10u && mode < 17u) {   // histone mark: nucleosomes (region 22) whose instance flag bit (mode - 8) is set
    if (I.region == 22u) {
      if (((I.flags >> (mode - 8u)) & 1u) != 0u) { c = vec3f(0.85, 0.20, 0.35); } else { c = vec3f(0.80, 0.78, 0.72); }
      if (kind == 1u) { c = c * 0.85; }
    } else { c = vec3f(0.90, 0.90, 0.89); }
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
`,wn=Qe+ao+ro+Zo,es=wn+`
@group(1) @binding(0) var<storage, read> beads: array<vec4f>;     // xyz, w = bitcast(kind | plddt<<8)
@group(1) @binding(1) var<storage, read> list01: array<vec2u>;
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
  if ((typeFlags[I.t] & 8u) == 0u) { bp = flexBead(b.xyz, select(255.0, pl, kind == 0u), I.seed); }
  return emitQuad(I.p + qrot(I.q, bp), r, corner, col, id);
}
@fragment
fn fsMain(in: VOut) -> FOut { return shadeFrag(in); }
`,ts=wn+`
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
`,ns=wn+`
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
`,Ee=Qe+ro+`
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
  if (kind == 1u) { r = 1.0; } else if (kind == 2u) { r = 0.75; } else if (kind == 4u) { r = 0.7; } else if (kind == 5u) { r = 0.45; } else if (kind == 6u) { r = 1.2; } else if (kind == 7u) { r = 0.45; }
  let hj = f32(hash32(chain * 747796405u + kind) & 255u) / 255.0;
  let col = palette[48u + kind].xyz * (0.94 + 0.12 * hj);
  let slot = spages[i / SPAGE].y;
  let id = ID_FREE | ((slot & 0xFFFu) << 18u) | ((kind & 7u) << 15u) | chain;
  var p = rec.xyz;
  if (u.caps.y > 0u) {     // illustrative motion: FG repeats and mRNA most, wall polysaccharides and DNA little
    var amps = array<f32, 8>(0.25, 0.08, 0.06, 0.25, 0.25, 0.5, 0.03, 0.3);
    let amp = amps[min(kind, 7u)] * f32(u.caps.y) / 1000.0;
    let t = u.cut.w; let ph = 6.2832 * f32(hash32(chain * 2654435761u + kind) & 1023u) / 1023.0;
    p += amp * vec3f(sin(dot(p, vec3f(0.9, 0.3, -0.2)) + 1.1 * t + ph), sin(dot(p, vec3f(-0.2, 0.9, 0.4)) + 0.9 * t + 2.0 * ph),
                     sin(dot(p, vec3f(0.3, -0.3, 0.9)) + 1.3 * t + 3.0 * ph));
  }
  return emitQuad(p - u.eye.xyz, r, vi, col, id);
}

@fragment
fn fsMain(in: VOut) -> FOut { return shadeFrag(in); }
`,os=Qe+`
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
  if (u.simBox.w > 0.0) {
    let q = abs(in.prel + u.eye.xyz - u.simBox.xyz);
    if (max(q.x, max(q.y, q.z)) > u.simBox.w) { discard; }
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
`,_e=Qe+`
@group(0) @binding(5) var<storage, read> palette: array<vec4f>;
struct Prim { q0: vec4f, q1: vec4f, tc: vec4u };      // type, class
struct Geo { G0: vec4f, G1: vec4f, G2: vec4f, consts: vec4f, grid: vec4f,       // consts: PORE_R, PERI_W, IMS_W, nPores; grid.z: bud scars
             nax: vec4f, nr0: vec4f, nr1: vec4f, nr2: vec4f,                      // nucleus semi-axes, rotation rows
             scars: array<vec4f, 8> };                                            // per scar: centre, major R; axis, tube r
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
  if (t == 4u) {   // eisosome furrow lumen (cell3d T_SLOT)
    let a = q0.xyz; let b = vec3f(q0.w, q1.x, q1.y); let u = normalize(0.5 * (a + b));
    let L = length(b - a); let e = (b - a) / L; let v = p - a;
    let s = clamp(dot(v, e), 0.0, L); let h = clamp(dot(v, u), 0.0, q1.w);
    return length(v - s * e - h * u) - q1.z;
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
// the same with the bud-scar rims (cell3d T_TORUS, smooth union with the outside surface only)
fn outerDepthS(p: vec3f) -> f32 {
  let raw = outerDepth(p);
  var s = 1e9;
  for (var k = 0u; k < u32(geo.grid.z); k++) {
    let a = geo.scars[2u * k]; let b = geo.scars[2u * k + 1u];
    let v = p - a.xyz; let h = dot(v, b.xyz); let rho = sqrt(max(dot(v, v) - h * h, 0.0));
    s = min(s, length(vec2f(rho - a.w, h)) - b.w);
  }
  let K = 30.0; let hh = max(K - abs(-raw - s), 0.0) / K;
  return -(min(-raw, s) - hh * hh * K * 0.25);
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
  m.d[9] = outerDepthS(p);              // outer (with bud scars)
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
  var im = 1e9; var cr = 1e9; var eis = 1e9;
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
    else if (c == 8u) { eis = min(eis, dd); }
  }
  let K = 8.0; let hh = max(K - abs(m.d[0] - eis), 0.0) / K;      // eisosome furrows: smooth union (cell3d.EIS_K)
  m.d[0] = min(m.d[0], eis) - hh * hh * K * 0.25;
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
  if (onCut && outerDepthS(p) >= 0.0) {
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
  if (!onCut && outerDepthS(eye) > 0.0) { return o; }   // camera inside the cell without a cut: nothing analytic to show
  for (var k = 0; k < 128; k++) {
    let d = -outerDepthS(p);
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
  let g = vec3f(outerDepthS(ps + vec3f(e, 0.0, 0.0)) - outerDepthS(ps - vec3f(e, 0.0, 0.0)),
                outerDepthS(ps + vec3f(0.0, e, 0.0)) - outerDepthS(ps - vec3f(0.0, e, 0.0)),
                outerDepthS(ps + vec3f(0.0, 0.0, e)) - outerDepthS(ps - vec3f(0.0, 0.0, e)));
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
`,ss=`
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
`,is=Qe+`
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
  // simulation: molecules faded while others are tracked (instance ids with bit 29) take no part in outlines
  let ghost = cat == 1u && (id & 0x20000000u) != 0u;
  for (var k = 0; k < 4; k++) {
    let q = clamp(px + offs[k], vec2i(0), dims - vec2i(1));
    var id2 = textureLoad(idTex, q, 0).x;
    if ((id2 >> 30u) == 1u && (id2 & 0x20000000u) != 0u) { id2 = 0u; }
    if (id2 != id && !ghost) {
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
`,J=(o,e)=>[o[0]+e[0],o[1]+e[1],o[2]+e[2]],Fe=(o,e)=>[o[0]-e[0],o[1]-e[1],o[2]-e[2]],O=(o,e)=>[o[0]*e,o[1]*e,o[2]*e],oe=(o,e)=>o[0]*e[0]+o[1]*e[1]+o[2]*e[2],Le=(o,e)=>[o[1]*e[2]-o[2]*e[1],o[2]*e[0]-o[0]*e[2],o[0]*e[1]-o[1]*e[0]],Je=o=>Math.hypot(o[0],o[1],o[2]),he=o=>{const e=Je(o)||1;return[o[0]/e,o[1]/e,o[2]/e]},on=(o,e,t)=>[o[0]+(e[0]-o[0])*t,o[1]+(e[1]-o[1])*t,o[2]+(e[2]-o[2])*t];function ut(o,e,t){const n=Math.cos(t),s=Math.sin(t),r=Le(e,o),a=oe(e,o)*(1-n);return[o[0]*n+r[0]*s+e[0]*a,o[1]*n+r[1]*s+e[1]*a,o[2]*n+r[2]*s+e[2]*a]}function as(o,e){const t=new Float32Array(16);for(let n=0;n<4;n++)for(let s=0;s<4;s++){let r=0;for(let a=0;a<4;a++)r+=o[a*4+s]*e[n*4+a];t[n*4+s]=r}return t}function rs(o){const e=Array.from(o),t=new Float32Array(16),[n,s,r,a,i,h,l,c,d,f,w,b,m,u,p,x]=e,y=n*h-s*i,A=n*l-r*i,M=n*c-a*i,_=s*l-r*h,S=s*c-a*h,g=r*c-a*l,k=d*u-f*m,P=d*p-w*m,T=d*x-b*m,E=f*p-w*u,I=f*x-b*u,F=w*x-b*p;let N=y*F-A*I+M*E+_*T-S*P+g*k;return N&&(N=1/N,t[0]=(h*F-l*I+c*E)*N,t[1]=(r*I-s*F-a*E)*N,t[2]=(u*g-p*S+x*_)*N,t[3]=(w*S-f*g-b*_)*N,t[4]=(l*T-i*F-c*P)*N,t[5]=(n*F-r*T+a*P)*N,t[6]=(p*M-m*g-x*A)*N,t[7]=(d*g-w*M+b*A)*N,t[8]=(i*I-h*T+c*k)*N,t[9]=(s*T-n*I-a*k)*N,t[10]=(m*S-u*M+x*y)*N,t[11]=(f*M-d*S-b*y)*N,t[12]=(h*P-i*E-l*k)*N,t[13]=(n*E-s*P+r*k)*N,t[14]=(u*A-m*_-p*y)*N,t[15]=(d*_-f*A+w*y)*N),t}function Bt(o){const e=parseInt(o.replace("#",""),16);return(e>>16&255|(e>>8&255)<<8|(e&255)<<16|4278190080)>>>0}function ls(o){return o<.5?4*o*o*o:1-Math.pow(-2*o+2,3)/2}const vn=new URL("data/",document.baseURI).href,cs=fetch(vn+"version.json",{cache:"no-store"}).then(o=>o.ok?o.json():{v:""}).then(o=>String(o.v??"")).catch(()=>""),ds=o=>o?`?v=${encodeURIComponent(o)}`:"";function lo(o){const e=atob(o),t=new Uint8Array(e.length);for(let n=0;n<e.length;n++)t[n]=e.charCodeAt(n);return t}async function Ht(o){const e=new Uint8Array(o,0,Math.min(2,o.byteLength));if(e[0]===31&&e[1]===139){const t=new DecompressionStream("gzip"),n=new Blob([o]).stream().pipeThrough(t);return await new Response(n).arrayBuffer()}return o}async function He(o,e){const t=window.__YC_EMBED__;if(t&&t[o]!==void 0)return Ht(lo(t[o]).buffer);const n=await fetch(vn+o+ds(await cs));if(!n.ok)throw new Error(`HTTP ${n.status} for ${o}`);if(!e||!n.body)return Ht(await n.arrayBuffer());const s=Number(n.headers.get("content-length")||0),r=n.body.getReader(),a=[];let i=0;for(;;){const{done:c,value:d}=await r.read();if(c)break;a.push(d),i+=d.length,e(i,s)}const h=new Uint8Array(i);let l=0;for(const c of a)h.set(c,l),l+=c.length;return Ht(h.buffer)}async function Ae(o){const e=await He(o);return JSON.parse(new TextDecoder().decode(e))}const bt=["cytosol","nucleoplasm","nucleolus","ER lumen","vacuole lumen","mitochondrial matrix","mitochondrial IMS","Golgi / endosome","peroxisome","cell wall","PM peripheral","NE peripheral","ER peripheral","lipid droplet surface","plasma membrane","ER membrane","nuclear envelope","vacuolar membrane","mito outer membrane","mito inner membrane","Golgi membrane","peroxisomal membrane","chromatin (nucleosome)","nuclear pore","vacuole surface (cytosolic face)","Golgi surface (cytosolic face)","peroxisome surface (cytosolic face)"],Dn=["plasma membrane","ER membrane","nuclear envelope","vacuolar membrane","mito outer membrane","mito inner membrane","Golgi/endosome membrane","peroxisomal membrane","lipid droplet monolayer","cell surface"],hs=["exterior","cell wall (mannoprotein layer)","cell wall (glucan/chitin layer)","cytosol","cortical ER lumen","perinuclear space","nucleoplasm","nucleolus","vacuole","mito intermembrane space","mito matrix","lipid droplet","Golgi / endosome","peroxisome","nuclear pore channel"],ee={0:"#7EA7D3",1:"#B993C9",2:"#D39BC0",3:"#A3C48C",4:"#9FC7C0",5:"#D9946B",6:"#E8B48F",7:"#B7C77E",8:"#C5B0E0",9:"#9A8352",10:"#7EA7D3",11:"#8FA3D0",12:"#7EA7D3",13:"#7EA7D3",14:"#5FAE75",15:"#7DBE86",16:"#6BB596",17:"#5BA8A0",18:"#C27C55",19:"#B5653F",20:"#93B060",21:"#A68BC9",22:"#F4B26A",23:"#3E9C9A",24:"#7EA7D3",25:"#7EA7D3",26:"#7EA7D3",32:"#F3C643",33:"#F0BC3C",34:"#F6D15E",35:"#EAB542",36:"#EFA83A",37:"#E7A93C",38:"#F2CA55",39:"#E8873A",40:"#F3C643",41:"#F8E3A0",42:"#E7A93C",48:"#B0437A",49:"#E07636",50:"#CDA866",51:"#C3C27A",52:"#C9558F",53:"#7CC6C2",54:"#E2D6A8",55:"#9B8BC4",56:"#5FAE75",57:"#7DBE86",58:"#6BB596",59:"#5BA8A0",60:"#C27C55",61:"#B5653F",62:"#93B060",63:"#A68BC9",64:"#E9C46A",65:"#9A8352",66:"#A9C3E3",67:"#A3C48C",68:"#B993C9",69:"#9FC7C0",70:"#E8B48F",71:"#D9946B",72:"#B7C77E",73:"#C5B0E0",74:"#E9C46A",75:"#9A8352",80:"#F6F1E7",81:"#F1E9D3",82:"#EEE3CC",83:"#E4ECF5",84:"#E8F0DF",85:"#E8F0DF",86:"#F1E7F1",87:"#EBD9EA",88:"#E3F0EE",89:"#F6E6DA",90:"#F2DCCB",91:"#F7EDC4",92:"#EEF1DC",93:"#EEE8F6",94:"#ECE8F2"};function co(o){const e=[];for(let t=0;t<o.PT.length&&e.length<32;t++){if(o.PC[t]!==9)continue;const n=o.PQ.slice(8*t,8*t+8);e.push(n[0],n[1],n[2],n[6],n[3],n[4],n[5],n[7])}return e}function us(){const o=new Float32Array(384);for(const[e,t]of Object.entries(ee)){const n=parseInt(t.slice(1),16),s=Number(e)*4;o[s]=(n>>16&255)/255,o[s+1]=(n>>8&255)/255,o[s+2]=(n&255)/255,o[s+3]=1}return o}const sn={ribosome_80S:{c0:"#EE9BC0",c1:"#C9558F",cls:"ribosome"},tRNA:{c0:"#E27FA8",c1:"#E27FA8",cls:"tRNA"},nucleosome:{c0:"#F4B26A",c1:"#E07636",cls:"nucleosome"},NPC:{c0:"#3E9C9A",cls:"NPC"},PolI:{c0:"#9A72B8",c1:"#9A72B8",cls:"polymerase"},PolII:{c0:"#9A72B8",c1:"#9A72B8",cls:"polymerase"},PolIII:{c0:"#9A72B8",c1:"#9A72B8",cls:"polymerase"},FAS:{c0:"#8B88C9",cls:"complex"},proteasome_26S:{c0:"#8B88C9",cls:"complex"},CCT:{c0:"#8B88C9",cls:"complex"},F_actin:{c0:"#E0893A",c1:"#E0893A",cls:"actin"},Arp23:{c0:"#D9534F",c1:"#D9534F",cls:"actin"},TPM1_coil:{c0:"#2FA4A0",c1:"#2FA4A0",cls:"actin"},CP_cap:{c0:"#8E6CB0",c1:"#8E6CB0",cls:"actin"},MT_tubulin:{c0:"#74C476",c1:"#74C476",cls:"microtubule"},Pil1_dimer:{c0:"#3F8FCB",c1:"#3F8FCB",cls:"eisosome"},APE1_oligo:{c0:"#E0A43A",cls:"complex"},AMS1_oligo:{c0:"#B9792E",cls:"complex"},ATG19:{c0:"#C0473F",cls:"complex"},Vac8_Nvj1:{c0:"#2E9C8F",cls:"complex"},RVS161_RVS167_het:{c0:"#9B4DCA",cls:"complex"},SLA2:{c0:"#E3B23C",cls:"complex"},SLA1:{c0:"#CFA040",cls:"complex"},END3:{c0:"#BB8A35",cls:"complex"},Lsp1_dimer:{c0:"#7FB6E0",c1:"#7FB6E0",cls:"eisosome"},gTuRC_nuc:{c0:"#2E8B57",c1:"#2E8B57",cls:"microtubule"},gTuRC_cyto:{c0:"#2E8B57",c1:"#2E8B57",cls:"microtubule"},ribosome_80S_Ssh1:{c0:"#EE9BC0",c1:"#C9558F",cls:"ribosome"},ribosome_80S_RAC:{c0:"#EE9BC0",c1:"#C9558F",cls:"ribosome"},ribosome_80S_Sec61:{c0:"#EE9BC0",c1:"#C9558F",cls:"ribosome"},preribo_90S:{c0:"#EF9A7E",c1:"#C9558F",cls:"pre-ribosome"},preribo_60S_nucleolar:{c0:"#EF9A7E",c1:"#C9558F",cls:"pre-ribosome"},preribo_60S_nucleoplasm:{c0:"#EF9A7E",c1:"#C9558F",cls:"pre-ribosome"},preribo_60S_cytoplasm:{c0:"#EF9A7E",c1:"#C9558F",cls:"pre-ribosome"},preribo_40S:{c0:"#EF9A7E",c1:"#C9558F",cls:"pre-ribosome"}};async function ms(o){o("Loading structure index",.02);const e=await Ae("library/types.json"),t=e.types,n=e.bead_q;o("Loading residue coordinates (2.5 M beads)",.06);const s=await He("library/beads.bin",(S,g)=>o("Loading residue coordinates (2.5 M beads)",.06+.5*(g?S/g:0))),r=new Float32Array(await He("library/spheretree.bin"));let a={},i={},h=!1;try{const S=await Ae("library/organisation.json");a=S.colors||{},i=S.neighbours||{},h=!!S.colors}catch{}o("Building level-of-detail clusters",.6);const l=s.byteLength/8,c=new Int16Array(s),d=new Uint8Array(s),f=new Float32Array(l*4),w=new Uint32Array(f.buffer);for(let S=0;S<l;S++)f[4*S]=c[4*S]*n,f[4*S+1]=c[4*S+1]*n,f[4*S+2]=c[4*S+2]*n,w[4*S+3]=d[8*S+6]|d[8*S+7]<<8;let b=0;for(const S of t)b+=Math.ceil(S.nb/64);const m=new ArrayBuffer(b*32),u=new Float32Array(m),p=new Uint32Array(m),x=new ArrayBuffer(t.length*64),y=new Uint32Array(x),A=new Float32Array(x);let M=0;t.forEach((S,g)=>{const k=M;for(let E=0;E<S.nb;E+=64){const I=Math.min(64,S.nb-E),F=S.bo+E;let N=1e9,G=1e9,$=1e9,Z=-1e9,V=-1e9,ie=-1e9;for(let z=0;z<I;z++){const j=4*(F+z),X=f[j],pe=f[j+1],se=f[j+2];X<N&&(N=X),X>Z&&(Z=X),pe<G&&(G=pe),pe>V&&(V=pe),se<$&&($=se),se>ie&&(ie=se)}const me=(N+Z)/2,be=(G+V)/2,ae=($+ie)/2;let L=0;for(let z=0;z<I;z++){const j=4*(F+z);L=Math.max(L,Math.hypot(f[j]-me,f[j+1]-be,f[j+2]-ae))}u[8*M]=me,u[8*M+1]=be,u[8*M+2]=ae,u[8*M+3]=L+.85,p[8*M+4]=F,p[8*M+5]=I,M++}const P=sn[S.id]??(S.id.startsWith("tRNA_")?sn.tRNA:void 0),T=S.src==="PDB";y[16*g+0]=S.bo,y[16*g+1]=S.nb,y[16*g+2]=k,y[16*g+3]=M-k,y[16*g+4]=S.so,y[16*g+5]=S.sn,A[16*g+6]=S.br,A[16*g+7]=Math.sqrt(5/3)*S.rg,y[16*g+8]=P?Bt(P.c0):0,y[16*g+9]=P&&P.c1?Bt(P.c1):0,y[16*g+10]=((T?1:0)|Math.min(S.copies,16777215)<<8)>>>0,A[16*g+11]=T?255:S.plddt_mean??255,y[16*g+15]=a[g]?Bt(a[g]):0});const _=new Map;return t.forEach((S,g)=>_.set(S.id,g)),{types:t,beads:f,st:r,clusters:m,typeInfo:x,netOn:h,netNb:i,nClusters:b,byId:_}}const pt=26,_t=(o,e,t)=>((o+13)*pt+(e+13))*pt+(t+13);async function ps(){const o=await Ae("tiles/index.json"),e=o.tile_columns,t=o.tiles,n=t.length,s=i=>e.indexOf(i),r=(i,h)=>{const l=new i(n),c=s(h);for(let d=0;d<n;d++)l[d]=t[d][c];return l},a={n,ix:r(Int32Array,"ix"),iy:r(Int32Array,"iy"),iz:r(Int32Array,"iz"),nInst:r(Uint32Array,"n_inst"),nGhost:r(Uint32Array,"n_ghost"),nMr:r(Uint32Array,"n_mrna_beads"),nLk:r(Uint32Array,"n_link_beads"),nVert:r(Uint32Array,"n_vert"),nTri:r(Uint32Array,"n_tri"),gz:r(Uint32Array,"gz_bytes"),labelMask:r(Uint32Array,"label_mask"),volGlucan:r(Float32Array,"vol_glucan_nm3"),volMannan:r(Float32Array,"vol_mannan_nm3"),key:new Map,wall:o.wall,regionNames:o.region_names};for(let i=0;i<n;i++)a.key.set(_t(a.ix[i],a.iy[i],a.iz[i]),i);return a}const fs=1073741824,Y=GPUShaderStage,U=GPUBufferUsage;class bs{device;ctx;format;canvas;limits={};adapterInfo="";instPages=0;lipPages=0;freePages=0;capL0=0;capL1=0;capL2=0;capLip=0;capFree=0;maxSlots=4096;uni;inst;pageOrg;typesBuf;typeFlags;palette;beads;st;clusters;counters;list01;listL2;args;cullParams;caps;lipPool;lipPagesBuf;lipList;lipParams;freePool;freePagesBuf;freeList;freeParams;meshUni;tileState;fillBufs=[];pickBuf;pickParams;pickOut;statsBuf;pipes={};bg={};bgl={};sim=null;colorTex;idTex;depthTex;width=0;height=0;usedInstPages=0;usedLipPages=0;usedFreePages=0;pickPending=null;pickBusy=!1;statsBusy=!1;stats={l0:0,l1:0,l2:0,lip:0,fre:0,ninst:0};lastStatsT=0;organelleMeshes=[];onLost=()=>{};constructor(e){this.canvas=e}async init(e,t){if(!("gpu"in navigator))throw new Error("WebGPU is not available in this browser.");const n=await navigator.gpu.requestAdapter({powerPreference:"high-performance"});if(!n)throw new Error("No WebGPU adapter found (GPU blocked or unsupported).");try{const d=n.info||{};this.adapterInfo=[d.vendor,d.architecture,d.device,d.description].filter(Boolean).join(" ")}catch{}const s=n.limits,r={maxStorageBufferBindingSize:Math.min(s.maxStorageBufferBindingSize,2**31-4),maxBufferSize:Math.min(s.maxBufferSize,2**31-4),maxComputeWorkgroupsPerDimension:s.maxComputeWorkgroupsPerDimension,maxStorageBuffersPerShaderStage:Math.min(s.maxStorageBuffersPerShaderStage,16)};this.device=await n.requestDevice({requiredLimits:r}),this.device.lost.then(d=>this.onLost(`GPU device lost: ${d.message}`)),this.device.addEventListener("uncapturederror",d=>console.error("WebGPU error:",d.error?.message)),this.limits={...r};const a=this.device.limits.maxStorageBufferBindingSize,i=d=>Math.min(d,Math.floor(a/16));this.instPages=Math.floor(i(8*1024*1024)/ce),this.lipPages=Math.floor(i(8*1024*1024)/de),this.freePages=Math.floor(i(10*1024*1024)/de),this.capL0=Math.min(4*1024*1024,Math.floor(a/8)-2*1024*1024),this.capL1=2*1024*1024,this.capL2=Math.min(8*1024*1024,Math.floor(a/4)),this.capLip=this.lipPages*de,this.capFree=this.freePages*de,this.ctx=this.canvas.getContext("webgpu"),this.format=navigator.gpu.getPreferredCanvasFormat(),this.ctx.configure({device:this.device,format:this.format,alphaMode:"opaque",usage:GPUTextureUsage.RENDER_ATTACHMENT|GPUTextureUsage.COPY_SRC});const h=this.device,l=(d,f,w)=>h.createBuffer({size:Math.max(16,Math.ceil(d/16)*16),usage:f,label:w}),c=U.STORAGE|U.COPY_DST;this.uni=l(512,U.UNIFORM|U.COPY_DST,"uniforms"),this.inst=l(this.instPages*ce*16,c,"instances"),this.pageOrg=l(this.instPages*16,c,"pageOrg"),this.typesBuf=l(e.typeInfo.byteLength,c,"types"),this.typeFlags=l(e.types.length*4,c,"typeFlags"),this.palette=l(1536,c,"palette"),this.beads=l(e.beads.byteLength,c,"beads"),this.st=l(e.st.byteLength,c,"spheretree"),this.clusters=l(e.clusters.byteLength,c,"clusters"),h.queue.writeBuffer(this.typesBuf,0,e.typeInfo),h.queue.writeBuffer(this.palette,0,us()),h.queue.writeBuffer(this.beads,0,e.beads),h.queue.writeBuffer(this.st,0,e.st),h.queue.writeBuffer(this.clusters,0,e.clusters),this.counters=l(32,U.STORAGE|U.COPY_DST|U.COPY_SRC,"counters"),this.list01=l((this.capL0+this.capL1)*8,U.STORAGE,"list01"),this.listL2=l(this.capL2*4,U.STORAGE,"listL2"),this.args=l(80,U.STORAGE|U.INDIRECT,"args"),this.cullParams=l(32,U.UNIFORM|U.COPY_DST,"cullParams"),this.caps=l(32,U.UNIFORM|U.COPY_DST,"caps"),h.queue.writeBuffer(this.caps,0,new Uint32Array([this.capL0,this.capL1,this.capL2,this.capLip,this.capFree,0,0,0])),this.lipPool=l(this.lipPages*de*16,c,"lipidPool"),this.lipPagesBuf=l(this.lipPages*16,c,"lipidPages"),this.lipList=l(this.capLip*4,U.STORAGE,"lipidList"),this.lipParams=l(16,U.UNIFORM|U.COPY_DST,"lipParams"),this.freePool=l(this.freePages*de*16,c,"freePool"),this.freePagesBuf=l(this.freePages*16,c,"freePages"),this.freeList=l(this.capFree*4,U.STORAGE,"freeList"),this.freeParams=l(16,U.UNIFORM|U.COPY_DST,"freeParams"),this.meshUni=l(256*(this.maxSlots+16),U.UNIFORM|U.COPY_DST,"meshUniforms"),this.tileState=l(676*26*4,c,"tileState"),this.pickBuf=h.createBuffer({size:16,usage:U.MAP_READ|U.COPY_DST}),this.pickParams=l(16,U.UNIFORM|U.COPY_DST,"pickParams"),this.pickOut=l(16,U.STORAGE|U.COPY_SRC,"pickOut"),this.statsBuf=h.createBuffer({size:32,usage:U.MAP_READ|U.COPY_DST}),this.initFill(t),await this.initPipelines(),this.resize()}initFill(e){const t=this.device,n=e.PT.length,s=new ArrayBuffer(n*48),r=new Float32Array(s),a=new Uint32Array(s);for(let x=0;x<n;x++){for(let y=0;y<8;y++)r[12*x+y]=e.PQ[8*x+y];a[12*x+8]=e.PT[x],a[12*x+9]=e.PC[x]}const i=e.pores.length/6,h=new Float32Array(Math.max(1,i)*8);for(let x=0;x<i;x++)for(let y=0;y<3;y++)h[8*x+y]=e.pores[6*x+y],h[8*x+4+y]=e.pores[6*x+3+y];const l=e.grid,c=l.prim_idx.length,d=new Uint32Array(l.prim_off.length+l.pore_off.length);d.set(l.prim_off,0);for(let x=0;x<l.pore_off.length;x++)d[l.prim_off.length+x]=l.pore_off[x]+c;const f=new Uint32Array(c+l.pore_idx.length+1);f.set(l.prim_idx,0),f.set(l.pore_idx,c);const w=e.G,b=e.consts,m=co(e),u=new Float32Array([w[0],w[1],w[2],w[3],w[4],w[5],w[6],w[7],w[8],w[9],w[10],0,b.PORE_R,b.PERI_W,b.IMS_W,i,l.n,l.ext,m.length/8,0,w[11],w[12],w[13],0,w[14],w[15],w[16],0,w[17],w[18],w[19],0,w[20],w[21],w[22],0,...m,...new Array(32-m.length).fill(0)]),p=(x,y)=>{const A=t.createBuffer({size:Math.ceil(x.byteLength/16)*16+16,usage:y|U.COPY_DST});return t.queue.writeBuffer(A,0,x),A};this.fillBufs=[p(s,U.STORAGE),p(h,U.STORAGE),p(d,U.STORAGE),p(f,U.STORAGE),p(u,U.UNIFORM),p(e.F??new Float32Array(4),U.STORAGE)]}async initPipelines(){const e=this.device,t=(z,j)=>({binding:z,visibility:j,buffer:{type:"read-only-storage"}}),n=z=>({binding:z,visibility:Y.COMPUTE,buffer:{type:"storage"}}),s=(z,j,X=!1)=>({binding:z,visibility:j,buffer:{type:"uniform",hasDynamicOffset:X}}),r=Y.COMPUTE|Y.VERTEX|Y.FRAGMENT,a=Y.COMPUTE|Y.VERTEX,i=Y.VERTEX|Y.FRAGMENT,h=z=>e.createBindGroupLayout({entries:z});this.bgl.g0=h([s(0,r),t(1,a),t(2,a),t(3,a),t(4,a),t(5,i)]),this.bgl.cull=h([n(0),n(1),n(2),s(3,Y.COMPUTE),t(4,Y.COMPUTE)]),this.bgl.cullS=h([n(0),t(1,Y.COMPUTE),t(2,Y.COMPUTE),n(3),s(4,Y.COMPUTE)]),this.bgl.fin=h([t(0,Y.COMPUTE),n(1),s(2,Y.COMPUTE)]),this.bgl.l0=h([t(0,Y.VERTEX),t(1,Y.VERTEX)]),this.bgl.l2=h([t(0,Y.VERTEX)]),this.bgl.sph=h([t(0,Y.VERTEX),t(1,Y.VERTEX),t(2,Y.VERTEX)]),this.bgl.mesh=h([s(0,i,!0)]),this.bgl.fill=h([t(0,Y.FRAGMENT),t(1,Y.FRAGMENT),t(2,Y.FRAGMENT),t(3,Y.FRAGMENT),s(4,Y.FRAGMENT),t(5,Y.FRAGMENT),t(6,Y.FRAGMENT)]),this.bgl.comp=h([{binding:0,visibility:Y.FRAGMENT,texture:{sampleType:"unfilterable-float"}},{binding:1,visibility:Y.FRAGMENT,texture:{sampleType:"uint"}},{binding:2,visibility:Y.FRAGMENT,texture:{sampleType:"depth"}}]),this.bgl.pick=h([{binding:0,visibility:Y.COMPUTE,texture:{sampleType:"uint"}},{binding:1,visibility:Y.COMPUTE,texture:{sampleType:"depth"}},s(2,Y.COMPUTE),n(3)]);const l=(...z)=>e.createPipelineLayout({bindGroupLayouts:z}),c=async(z,j)=>{const X=e.createShaderModule({code:z,label:j}),pe=await X.getCompilationInfo(),se=pe.messages.filter(ne=>ne.type==="error");for(const ne of pe.messages)(ne.type==="error"?console.error:console.warn)(`[${j}] ${ne.lineNum}:${ne.linePos} ${ne.message}`);if(se.length)throw new Error(`WGSL compile error in ${j}: ${se[0].lineNum}:${se[0].linePos} ${se[0].message}`);return X},d=await c(Jo,"cull"),f=await c(Xo,"cullSpheres"),w=await c(Qo,"finalize"),b=await c(es,"drawL0"),m=await c(ts,"drawL1"),u=await c(ns,"drawL2"),p=await c(Ee,"drawSpheres"),x=await c(os,"mesh"),y=await c(_e,"fill"),A=await c(is,"composite"),M=await c(ss,"pick"),_=(z,j,X)=>e.createComputePipelineAsync({layout:X,compute:{module:z,entryPoint:j}}),S=[{format:"rgba8unorm"},{format:"r32uint"}],g=(z="greater")=>({format:"depth32float",depthWriteEnabled:!0,depthCompare:z}),k=(z,j,X,pe,se={})=>e.createRenderPipelineAsync({layout:pe,vertex:{module:z,entryPoint:j,...se.vertex?{buffers:se.vertex.buffers}:{}},fragment:{module:z,entryPoint:X,targets:S},primitive:{topology:"triangle-list",cullMode:"none"},depthStencil:se.depthStencil??g()}),P=this.bgl.g0,[T,E,I,F,N,G,$,Z,V,ie,me,be]=await Promise.all([_(d,"cullInstances",l(P,this.bgl.cull)),_(f,"cullSpheres",l(P,this.bgl.cullS)),_(w,"main",l(this.bgl.fin)),k(b,"vsL0","fsMain",l(P,this.bgl.l0)),k(m,"vsL1","fsMain",l(P,this.bgl.l0)),k(u,"vsL2","fsMain",l(P,this.bgl.l2)),k(p,"vsLipid","fsMain",l(P,this.bgl.sph)),k(p,"vsFree","fsMain",l(P,this.bgl.sph)),k(x,"vsMesh","fsMesh",l(P,this.bgl.mesh),{vertex:{buffers:[{arrayStride:32,attributes:[{shaderLocation:0,offset:0,format:"float32x4"},{shaderLocation:1,offset:16,format:"float32x4"}]}]}}),k(y,"vsFill","fsFill",l(P,this.bgl.fill),{depthStencil:g("always")}),e.createRenderPipelineAsync({layout:l(P,this.bgl.comp),vertex:{module:A,entryPoint:"vsFull"},fragment:{module:A,entryPoint:"fsComposite",targets:[{format:this.format}]},primitive:{topology:"triangle-list"}}),_(M,"main",l(this.bgl.pick))]);Object.assign(this.pipes,{cull:T,cullS:E,fin:I,l0:F,l1:N,l2:G,lip:$,fre:Z,mesh:V,fill:ie,comp:me,pick:be});const ae=(z,j,X=[])=>e.createBindGroup({layout:z,entries:j.map((pe,se)=>({binding:se,resource:{buffer:pe,size:X[se]}}))});this.bg.g0=ae(P,[this.uni,this.inst,this.pageOrg,this.typesBuf,this.typeFlags,this.palette]),this.bg.cull=ae(this.bgl.cull,[this.counters,this.list01,this.listL2,this.cullParams,this.clusters]),this.bg.cullLip=ae(this.bgl.cullS,[this.counters,this.lipPool,this.lipPagesBuf,this.lipList,this.lipParams]),this.bg.cullFree=ae(this.bgl.cullS,[this.counters,this.freePool,this.freePagesBuf,this.freeList,this.freeParams]),this.bg.fin=ae(this.bgl.fin,[this.counters,this.args,this.caps]),this.bg.l0=ae(this.bgl.l0,[this.beads,this.list01]),this.bg.l1=ae(this.bgl.l0,[this.st,this.list01]),this.bg.l2=ae(this.bgl.l2,[this.listL2]),this.bg.lip=ae(this.bgl.sph,[this.lipPool,this.lipPagesBuf,this.lipList]),this.bg.fre=ae(this.bgl.sph,[this.freePool,this.freePagesBuf,this.freeList]),this.bg.mesh=e.createBindGroup({layout:this.bgl.mesh,entries:[{binding:0,resource:{buffer:this.meshUni,size:16}}]});const L=this.fillBufs;this.bg.fill=ae(this.bgl.fill,[L[0],L[1],L[2],L[3],L[4],this.tileState,L[5]])}resize(e=1){const t=Math.min(window.devicePixelRatio||1,2)*e,n=Math.max(1,Math.floor(this.canvas.clientWidth*t)),s=Math.max(1,Math.floor(this.canvas.clientHeight*t));if(n===this.width&&s===this.height&&this.colorTex)return;this.canvas.width=n,this.canvas.height=s,this.width=n,this.height=s,this.colorTex?.destroy(),this.idTex?.destroy(),this.depthTex?.destroy();const r=GPUTextureUsage;this.colorTex=this.device.createTexture({size:[n,s],format:"rgba8unorm",usage:r.RENDER_ATTACHMENT|r.TEXTURE_BINDING}),this.idTex=this.device.createTexture({size:[n,s],format:"r32uint",usage:r.RENDER_ATTACHMENT|r.TEXTURE_BINDING|r.COPY_SRC}),this.depthTex=this.device.createTexture({size:[n,s],format:"depth32float",usage:r.RENDER_ATTACHMENT|r.TEXTURE_BINDING|r.COPY_SRC}),this.bg.comp=this.device.createBindGroup({layout:this.bgl.comp,entries:[{binding:0,resource:this.colorTex.createView()},{binding:1,resource:this.idTex.createView()},{binding:2,resource:this.depthTex.createView()}]}),this.bg.pick=this.device.createBindGroup({layout:this.bgl.pick,entries:[{binding:0,resource:this.idTex.createView()},{binding:1,resource:this.depthTex.createView()},{binding:2,resource:{buffer:this.pickParams}},{binding:3,resource:{buffer:this.pickOut}}]})}writeUniforms(e,t,n){const s=this.width/this.height,r=e.viewRot(),a=e.proj(s),i=as(a,r),h=rs(i),l=e.eye,c=new ArrayBuffer(512),d=new Float32Array(c),f=new Uint32Array(c),w=t.memColors&&!t.organelle;d.set(a,0),d.set(r,16),d.set(h,32);const b=e.pxPerNmAtUnit(this.height);d.set([l[0],l[1],l[2],b],48);const m=[t.plane[0],t.plane[1],t.plane[2]];d.set([m[0],m[1],m[2],t.plane[3]],52);const u=[r[0]*m[0]+r[4]*m[1]+r[8]*m[2],r[1]*m[0]+r[5]*m[1]+r[9]*m[2],r[2]*m[0]+r[6]*m[1]+r[10]*m[2]];d.set([u[0],u[1],u[2],oe(m,l)-t.plane[3]],56),d.set([t.cutOn&&!w?1:0,t.slab,Math.min(400,Math.max(20,t.slab)),performance.now()/1e3],60),d.set([t.l0px,t.l1px,t.minpx,.6],64),d.set([this.width,this.height,1/this.width,1/this.height],68),d.set(e.frustum(s),72),f.set([t.colorMode,t.regionMask>>>0,t.freeHide,t.highlight?1:0],96),f.set([t.selId>>>0,0,t.shaded?1:0,t.organelle?1:0],100),d.set([360/65535,.8,t.outline,Math.min(window.devicePixelRatio||1,2)*t.renderScale],104),d.set([n[0],n[1],n[2],e.near()],108),f.set([this.capL0,Math.round(t.motion*1e3),0,w?1:0],112),d.set(t.simMode&&this.sim?.clipBox()||[0,0,0,0],116),this.device.queue.writeBuffer(this.uni,0,c,0,480)}frame(e,t,n,s,r){const a=this.device;this.writeUniforms(e,t,n);const i=a.createCommandEncoder();i.clearBuffer(this.counters);const h=t.memColors&&!t.organelle,l=!!t.simMode&&!!this.sim,c=!t.organelle&&!h&&!l;if(l&&this.sim.compute(i),c){const u=i.beginComputePass();u.setBindGroup(0,this.bg.g0);const p=(x,y,A)=>{const M=Math.ceil(x/64),_=Math.min(M,65535),S=Math.ceil(M/_);return a.queue.writeBuffer(y,0,A(_*64)),[_,S]};if(this.usedInstPages>0){const[x,y]=p(this.usedInstPages*ce,this.cullParams,A=>new Uint32Array([this.usedInstPages,A,this.capL0,this.capL1,this.capL2,this.capLip,this.capFree,0]).buffer);u.setPipeline(this.pipes.cull),u.setBindGroup(1,this.bg.cull),u.dispatchWorkgroups(x,y)}if(this.usedLipPages>0&&t.showLipids){const[x,y]=p(this.usedLipPages*de,this.lipParams,A=>new Uint32Array([this.usedLipPages,A,this.capLip,1]).buffer);u.setPipeline(this.pipes.cullS),u.setBindGroup(1,this.bg.cullLip),u.dispatchWorkgroups(x,y)}if(this.usedFreePages>0){const[x,y]=p(this.usedFreePages*de,this.freeParams,A=>new Uint32Array([this.usedFreePages,A,this.capFree,0]).buffer);u.setPipeline(this.pipes.cullS),u.setBindGroup(1,this.bg.cullFree),u.dispatchWorkgroups(x,y)}u.setPipeline(this.pipes.fin),u.setBindGroup(0,this.bg.fin),u.dispatchWorkgroups(1),u.end()}const d=i.beginRenderPass({colorAttachments:[{view:this.colorTex.createView(),loadOp:"clear",storeOp:"store",clearValue:{r:.965,g:.945,b:.906,a:1}},{view:this.idTex.createView(),loadOp:"clear",storeOp:"store",clearValue:{r:0,g:0,b:0,a:0}}],depthStencilAttachment:{view:this.depthTex.createView(),depthLoadOp:"clear",depthStoreOp:"store",depthClearValue:0}});if(d.setBindGroup(0,this.bg.g0),d.setPipeline(this.pipes.fill),d.setBindGroup(1,this.bg.fill),d.draw(3),t.organelle||h){d.setPipeline(this.pipes.mesh);const u=this.maxSlots;a.queue.writeBuffer(this.meshUni,u*256,new Uint32Array([0,t.organelle?1:3,0,0]));for(const p of this.organelleMeshes)d.setBindGroup(1,this.bg.mesh,[u*256]),d.setVertexBuffer(0,p.vb),d.setIndexBuffer(p.ib,"uint32"),d.drawIndexed(p.nIdx,1)}else if(t.showMeshes&&s.length){d.setPipeline(this.pipes.mesh);for(const u of s)d.setBindGroup(1,this.bg.mesh,[u.slot*256]),d.setVertexBuffer(0,u.vb),d.setIndexBuffer(u.ib,"uint32"),d.drawIndexed(u.nIdx,2)}l&&this.sim.draw(d),c&&(d.setPipeline(this.pipes.l2),d.setBindGroup(1,this.bg.l2),d.drawIndirect(this.args,32),d.setPipeline(this.pipes.l1),d.setBindGroup(1,this.bg.l1),d.drawIndirect(this.args,16),d.setPipeline(this.pipes.l0),d.setBindGroup(1,this.bg.l0),d.drawIndirect(this.args,0),t.showLipids&&(d.setPipeline(this.pipes.lip),d.setBindGroup(1,this.bg.lip),d.drawIndirect(this.args,48)),d.setPipeline(this.pipes.fre),d.setBindGroup(1,this.bg.fre),d.drawIndirect(this.args,64)),this.sim&&this.sim.drawOverlay(d),d.end();const f=this.ctx.getCurrentTexture(),w=i.beginRenderPass({colorAttachments:[{view:f.createView(),loadOp:"clear",storeOp:"store",clearValue:{r:1,g:1,b:1,a:1}}]});w.setPipeline(this.pipes.comp),w.setBindGroup(0,this.bg.g0),w.setBindGroup(1,this.bg.comp),w.draw(3),w.end();let b=null;if(this.pickPending&&!this.pickBusy){b=this.pickPending,this.pickPending=null;const u=Math.max(0,Math.min(this.width-1,Math.floor(b.x))),p=Math.max(0,Math.min(this.height-1,Math.floor(b.y)));a.queue.writeBuffer(this.pickParams,0,new Uint32Array([u,p,0,0]));const x=i.beginComputePass();x.setPipeline(this.pipes.pick),x.setBindGroup(0,this.bg.pick),x.dispatchWorkgroups(1),x.end(),i.copyBufferToBuffer(this.pickOut,0,this.pickBuf,0,16)}let m=!1;if(r&&!this.statsBusy&&c&&(i.copyBufferToBuffer(this.counters,0,this.statsBuf,0,32),m=!0),a.queue.submit([i.finish()]),b){this.pickBusy=!0;const u=b;this.pickBuf.mapAsync(GPUMapMode.READ).then(()=>{const p=this.pickBuf.getMappedRange(0,8),x=new Uint32Array(p)[0],y=new Float32Array(p)[1];this.pickBuf.unmap(),this.pickBusy=!1,u.resolve({id:x>>>0,depth:y})}).catch(()=>{this.pickBusy=!1})}m&&(this.statsBusy=!0,this.statsBuf.mapAsync(GPUMapMode.READ).then(()=>{const u=new Uint32Array(this.statsBuf.getMappedRange());this.stats={l0:u[0],l1:u[1],l2:u[2],lip:u[3],fre:u[4],ninst:u[5]},this.statsBuf.unmap(),this.statsBusy=!1}).catch(()=>{this.statsBusy=!1}))}pick(e,t){return new Promise(n=>{this.pickPending={x:e,y:t,resolve:n}})}setMeshSlot(e,t){const n=new ArrayBuffer(16),s=new Uint32Array(n),r=new Float32Array(n);s[0]=e,s[1]=t?2:0,r[2]=t?-1:-1.9,r[3]=t?1:1.9,this.device.queue.writeBuffer(this.meshUni,e*256,n)}uploadOrganelleMeshes(e){const t=this.device,n=new Uint32Array(e,0,1)[0],s=new Uint32Array(e,4,n*3);let r=4+n*12;for(const a of this.organelleMeshes)a.vb.destroy(),a.ib.destroy();this.organelleMeshes=[];for(let a=0;a<n;a++){const i=s[3*a],h=s[3*a+1],l=s[3*a+2],c=new Float32Array(e.slice(r,r+h*12));r+=h*12;const d=new Uint32Array(e.slice(r,r+l*12));r+=l*12;const f=new Float32Array(h*3);for(let p=0;p<l;p++){const x=d[3*p],y=d[3*p+1],A=d[3*p+2],M=c[3*y]-c[3*x],_=c[3*y+1]-c[3*x+1],S=c[3*y+2]-c[3*x+2],g=c[3*A]-c[3*x],k=c[3*A+1]-c[3*x+1],P=c[3*A+2]-c[3*x+2],T=_*P-S*k,E=S*g-M*P,I=M*k-_*g;for(const F of[x,y,A])f[3*F]+=T,f[3*F+1]+=E,f[3*F+2]+=I}const w=new Float32Array(h*8),b=new Uint32Array(w.buffer);for(let p=0;p<h;p++){const x=Math.hypot(f[3*p],f[3*p+1],f[3*p+2])||1;w[8*p]=c[3*p],w[8*p+1]=c[3*p+1],w[8*p+2]=c[3*p+2],b[8*p+3]=i,w[8*p+4]=f[3*p]/x,w[8*p+5]=f[3*p+1]/x,w[8*p+6]=f[3*p+2]/x}const m=t.createBuffer({size:w.byteLength,usage:U.VERTEX|U.COPY_DST});t.queue.writeBuffer(m,0,w);const u=t.createBuffer({size:Math.ceil(d.byteLength/4)*4,usage:U.INDEX|U.COPY_DST});t.queue.writeBuffer(u,0,d),this.organelleMeshes.push({vb:m,ib:u,nIdx:l*3,mem:i})}}gpuMemoryMB(){return[this.inst,this.pageOrg,this.beads,this.st,this.clusters,this.list01,this.listL2,this.lipPool,this.lipList,this.freePool,this.freeList,this.meshUni].reduce((t,n)=>t+n.size,0)/1048576+this.width*this.height*12/1048576}}const xn=303.15,gt=1380649e-29*xn,yt=797e-6,kn=.3;function ho(o,e,t){const n=t.get(e);if(n)return n;const s=o.types[e],r=o.beads,a=s.nb,i=s.bo;let h=.5;if(a>=2){const w=Math.min(a,700),b=a/w;let m=0,u=0;for(let p=0;p<w;p++){const x=4*(i+Math.floor(p*b));for(let y=p+1;y<w;y++){const A=4*(i+Math.floor(y*b)),M=Math.hypot(r[x]-r[A],r[x+1]-r[A+1],r[x+2]-r[A+2]);m+=1/Math.max(M,.38),u++}}h=u?u/m:.5}const l=Math.max(.5,1.2*h+kn),c=gt/(6*Math.PI*yt*l*1e-9)*1e9,d=gt/(8*Math.PI*yt*(l*1e-9)**3)*1e-9,f={rh:l,dt:c,dr:d};return t.set(e,f),f}const Gt=o=>Math.max(.05,1-1.832*o-.219*o*o),Ot=o=>Math.max(.05,1-.631*o-.726*o*o),uo=(.73+.35)/602.214,An=.1,mo=4e-9;function po(o,e=An){const t=Math.min(5,Math.max(1,.5*o))*1e-9,n=e*mo,s=gt/(4*Math.PI*n)*(Math.log(n/(yt*t))-.5772)*1e9,r=gt/(4*Math.PI*n*t*t)*1e-9;return{dt:s,dr:r}}const Mn=.72,Ft=-.25,Pn=o=>3.3*Math.sqrt(Math.max(1e-4,o)),fo=127,bo=.49,an=1.5,go={exit:[-6.3208,-3.1478,5.8007],entry:[-11.2279,-1.9972,-2.6486],exitDir:[-.1301,.3672,.921],entryDir:[-.6272,-.7241,-.2869],linkMin:3,linkMax:12,cone:.7},Sn={schematic:{label:"Schematic (rigid sphere trees)",dt:.2,steps:200,kRep:20},mid:{label:"Mid (clusters of ~4 residues)",dt:.2,steps:60,kRep:30},residue:{label:"Residue detail (rigid, residue-bead contacts)",dt:.1,steps:20,kRep:50}};function Cn(o,e){const t=Sn[o].dt;return e>0?Math.min(t,yo[o]/e):t}const yo={schematic:.05,mid:.01,residue:.0125},wo={schematic:2,mid:1,residue:.5},vo={schematic:"Stickiness reaches the measured value at 1.5–2 kT, but by joining most molecules into one network within a microsecond.",mid:"Stickiness of 1 kT makes most molecules aggregate (about 12 µm²/s); stronger values are unstable and are not offered.",residue:"Stickiness of 0.5 kT makes most molecules aggregate (25–32 µm²/s)."},xo=2.75;function ko(o="schematic"){const e=Sn[o];return{mode:o,dt:e.dt,stepsPerFrame:e.steps,eps:0,viscosity:xo,hydro:!0,boundary:"auto",tether:1,memFluid:!0,memVisc:An,memCalibrated:!0,kRep:e.kRep,shell:8,recordEvery:200,elec:!0,ionic:.3,specific:!0,epsGo:20,mrna:!0}}const gs=Object.freeze(Object.defineProperty({__proto__:null,BJERRUM_NM:Mn,ETA_WATER:yt,HYDRATION:kn,HYD_VOL_PER_DA:uo,KT:gt,MEM_H:mo,MEM_VISC:An,MODES:Sn,MRNA_LINK_NT:fo,MRNA_LP:an,MRNA_NM_PER_NT:bo,NT_CHARGE:Ft,RIBO_CHANNEL:go,STICKY_DT:yo,STICKY_MAX:wo,STICKY_NOTE:vo,T_K:xn,VISC_CALIBRATED:xo,debyeKappa:Pn,defaultSettings:ko,dtMax:Cn,hiRot:Ot,hiTrans:Gt,hydro:ho,saffmanDelbruck:po},Symbol.toStringTag,{value:"Module"})),ys=16,fe=64,Wt=.6,ws=.15,Ln=3,Nt=40,vs=_e.slice(_e.indexOf("struct Prim"),_e.indexOf("fn tileLoaded"));function xs(o,e,t,n){const s=o[0],r=o[1],a=o[2],i=o[3],h=2*(r*n-a*t),l=2*(a*e-s*n),c=2*(s*t-r*e);return[e+i*h+(r*c-a*l),t+i*l+(a*h-s*c),n+i*c+(s*l-r*h)]}const Ao=`
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
`,ks=Ao+`
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

@compute @workgroup_size(${fe})
fn clearGrid(@builtin(global_invocation_id) g: vec3u) {
  let i = g.x; if (i >= su.offs.z) { return; }
  atomicStore(&cellCount[i], 0u);
  atomicStore(&cellCount[su.offs.z + 2u * i], 0xFFFFFFFFu); atomicStore(&cellCount[su.offs.z + 2u * i + 1u], 0u);
}
@compute @workgroup_size(1)
fn tick() { ctr[0] = ctr[0] + 1u; }

@compute @workgroup_size(${fe})
fn pose(@builtin(global_invocation_id) g: vec3u) {
  let s = g.x; if (s >= su.counts.y) { return; }
  let b = sphBody[s]; let l = sphLocal[s];
  sphWorld[s] = vec4f(bodyPos[b].xyz + qrot(bodyQuat[b], l.xyz), l.w);
}

@compute @workgroup_size(${fe})
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

@compute @workgroup_size(${fe})
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
@compute @workgroup_size(${fe})
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
@compute @workgroup_size(${fe})
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

var<workgroup> redF: array<vec3f, ${fe}>;
var<workgroup> redT: array<vec3f, ${fe}>;
@compute @workgroup_size(${fe})
fn integrate(@builtin(workgroup_id) wg: vec3u, @builtin(local_invocation_index) li: u32) {
  let b = wg.x + wg.y * 65535u;
  if (b >= su.counts.x) { return; }
  let info = bodyInfo[b];
  if ((info.z & 1u) == 0u) { return; }
  // force and torque over the body's spheres: strided partial sums, then a tree reduction in workgroup memory
  let x = bodyPos[b].xyz;
  var Fp = vec3f(0.0); var Tp = vec3f(0.0);
  for (var s = info.x + li; s < info.x + info.y; s += ${fe}u) {
    let f = sphForce[s].xyz; Fp += f; Tp += cross(sphWorld[s].xyz - x, f);
  }
  redF[li] = Fp; redT[li] = Tp;
  workgroupBarrier();
  for (var st = ${fe/2}u; st > 0u; st >>= 1u) {
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
`,As=Ao+vs+`
@group(2) @binding(0) var<storage, read_write> wallOut: array<vec4f>;
fn nearest(p: vec3f) -> vec2f {
  let m = memSdfs(p);
  var best = 1e9; var k = 0u;
  for (var i = 0u; i < 9u; i++) { if (abs(m.d[i]) < abs(best)) { best = m.d[i]; k = i; } }
  return vec2f(best, f32(k));
}
fn memDist(p: vec3f, k: u32) -> f32 { let m = memSdfs(p); return m.d[k]; }
@compute @workgroup_size(${fe})
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
`,jt={clearGrid:[0,8],pose:[0,1,2,4,5,6],bin:[0,5,6,8,9,13],force:[0,3,5,6,7,8,9,10,13,14],bondForce:[0,1,2,15,16,17],integrate:[0,1,2,3,6,7,10,11,12,16],tick:[11],anchorInit:[0,1,3,10,12]};class Ms{r;d;sys;set;uni;bufs={};pipes={};nSph=0;rmax=.5;rmaxReg=.5;rSplit=.5;nLarge=0;nLargeCore=0;nBonds=0;A;B;gridInfo={h:0,rmax:0,cap:0,cells:0,reach:0,hB:0,capB:0,cellsB:0,split:0,large:0};stepNo=0;timeNs=0;seed=Math.random()*1e9>>>0;drawPipe;drawBG;drawG0;nBeads=0;linkPipe=null;linkBG=null;constructor(e,t,n){this.r=e,this.d=e.device,this.sys=t,this.set=n,this.uni=this.d.createBuffer({size:256,usage:GPUBufferUsage.UNIFORM|GPUBufferUsage.COPY_DST})}buf(e,t,n=GPUBufferUsage.STORAGE|GPUBufferUsage.COPY_DST|GPUBufferUsage.COPY_SRC){const s=typeof t=="number"?t:t.byteLength,r=this.d.createBuffer({size:Math.max(16,Math.ceil(s/16)*16),usage:n,label:"sim."+e});return typeof t!="number"&&t.byteLength&&this.d.queue.writeBuffer(r,0,t),this.bufs[e]?.destroy(),this.bufs[e]=r,r}grid(e,t,n){const s=this.sys,r=s.boundary==="periodic",a=s.hi[0]-s.lo[0],i=2*e+Wt+.05,h=r?a/Math.max(1,Math.floor(a/i)):i,l=[0,1,2].map(b=>r?Math.max(1,Math.round((s.hi[b]-s.lo[b])/h)):Math.ceil((s.hi[b]-s.lo[b])/h)+2),c=l[0]*l[1]*l[2],d=(b,m)=>{if(r){const u=b-s.lo[m]-a*Math.floor((b-s.lo[m])/a);return Math.min(l[m]-1,Math.floor(u/h))}return Math.min(l[m]-1,Math.max(0,Math.floor((b-s.lo[m])/h)+1))},f=new Uint32Array(c);for(let b=0;b<this.nSph;b++)t(b)&&f[(d(n[3*b],0)*l[1]+d(n[3*b+1],1))*l[2]+d(n[3*b+2],2)]++;let w=0;for(const b of f)w=Math.max(w,b);return{h,dims:l,cap:Math.max(ys,Math.ceil(w*1.5)+4),cells:c}}async init(){const e=this.d,t=this.sys,n=t.n,s=new Float32Array(4*n),r=new Uint32Array(4*n),a=new Float32Array(r.buffer);for(let y=0;y<n;y++){s[4*y]=t.pos[3*y],s[4*y+1]=t.pos[3*y+1],s[4*y+2]=t.pos[3*y+2],s[4*y+3]=t.dt[y],r[4*y]=t.sphOff[y],r[4*y+1]=t.sphCnt[y];let A=.5;for(let M=t.sphOff[y];M<t.sphOff[y]+t.sphCnt[y];M++)A=Math.max(A,Math.hypot(t.sphLocal[4*M],t.sphLocal[4*M+1],t.sphLocal[4*M+2])+t.sphLocal[4*M+3]);r[4*y+2]=(t.mobile[y]?1:0)|(t.wall[y]?2:0)|(t.tether[y]?4:0)|(t.mem2d[y]?8:0)|Math.min(16777215,Math.round(A*100))<<8,a[4*y+3]=t.dr[y]}this.nSph=t.sphLocal.length/4;const i=new Uint32Array(this.nSph);for(let y=0;y<n;y++)i.fill(y,t.sphOff[y],t.sphOff[y]+t.sphCnt[y]);this.buf("bodyPos",s),this.buf("bodyQuat",t.quat),this.buf("bodyInfo",r),this.buf("bodyAnchor",s.slice()),this.buf("sphLocal",t.sphLocal),this.buf("sphBody",i),this.buf("sphWorld",16*this.nSph),this.buf("sphForce",16*this.nSph),this.buf("sphProp",t.sphProp),this.buf("native",t.native.length?t.native:new Uint32Array(4)),this.nBonds=t.bonds.length/12;const h=new Uint32Array(2*Math.max(1,n));for(let y=0;y<n;y++)h[2*y]=t.bondOff[y],h[2*y+1]=t.bondCnt[y];this.buf("bonds",this.nBonds?t.bonds:new Float32Array(12)),this.buf("bodyBond",h),this.buf("bodyExt",32*Math.max(1,n));const l=new Float32Array(3*this.nSph),c=new Float32Array(this.nSph);for(let y=0;y<n;y++){const A=t.quat.subarray(4*y,4*y+4);for(let M=t.sphOff[y];M<t.sphOff[y]+t.sphCnt[y];M++){const _=xs(A,t.sphLocal[4*M],t.sphLocal[4*M+1],t.sphLocal[4*M+2]);l[3*M]=t.pos[3*y]+_[0],l[3*M+1]=t.pos[3*y+1]+_[1],l[3*M+2]=t.pos[3*y+2]+_[2],c[M]=t.sphLocal[4*M+3]}}const d=c.slice().sort();this.rmax=d[this.nSph-1]??.5,this.rSplit=d[Math.floor(.9*(this.nSph-1))]??.5,this.rmax<=1.05*this.rSplit&&(this.rSplit=this.rmax);const f=this.rSplit;this.A=this.grid(f,y=>c[y]<=f,l),this.nLarge=0,this.nLargeCore=0,this.rmaxReg=0;for(let y=0;y<this.nSph;y++)t.sphProp[4*y+3]!==4294967294&&(this.rmaxReg=Math.max(this.rmaxReg,c[y]));for(let y=0;y<this.nSph;y++)c[y]>f&&(t.sphProp[4*y+3]===4294967294?this.nLargeCore++:this.nLarge++);this.B=this.nLarge+this.nLargeCore?this.grid(this.rmax,y=>c[y]>f,l):{h:1,dims:[1,1,1],cap:1,cells:1};const w=this.A.cells+this.B.cells;this.gridInfo={h:this.A.h,rmax:this.rmax,cap:this.A.cap,cells:this.A.cells,reach:Math.ceil((2*f+Wt)/this.A.h),hB:this.B.h,capB:this.B.cap,cellsB:this.nLarge+this.nLargeCore?this.B.cells:0,split:f,large:this.nLarge+this.nLargeCore},this.buf("cellCount",12*w),this.buf("cellItems",8*(this.A.cells*this.A.cap+this.B.cells*this.B.cap));const b=Math.max(...[0,1,2].map(y=>t.hi[y]-t.lo[y])),m=Math.max(1,b/96),u=[0,1,2].map(y=>Math.max(2,Math.ceil((t.hi[y]-t.lo[y])/m)));this.buf("wallGrid",16*u[0]*u[1]*u[2]),this.wall={vox:m,wd:u},this.writeUniforms();const x=await(async(y,A)=>{const M=e.createShaderModule({code:y,label:A}),S=(await M.getCompilationInfo()).messages.filter(g=>g.type==="error");for(const g of S)console.error(`[${A}] ${g.lineNum}:${g.linePos} ${g.message}`);if(S.length)throw new Error(`WGSL error in ${A}: ${S[0].message}`);return M})(ks,"sim");this.buf("ctr",new Uint32Array([this.stepNo,0,0,0]));for(const y of Object.keys(jt))this.pipes[y]=await e.createComputePipelineAsync({layout:"auto",compute:{module:x,entryPoint:y}});this.makeBindGroups(),await this.buildWalls(u);{const y=e.createCommandEncoder(),A=y.beginComputePass();A.setPipeline(this.pipes.anchorInit),A.setBindGroup(0,this.bgFor.anchorInit),A.dispatchWorkgroups(Math.ceil(n/fe)),A.end(),e.queue.submit([y.finish()])}await this.initDraw()}bgFor;wall;makeBindGroups(){const e=this.bufs,t=[this.uni,e.bodyPos,e.bodyQuat,e.bodyInfo,e.sphLocal,e.sphBody,e.sphWorld,e.sphForce,e.cellCount,e.cellItems,e.wallGrid,e.ctr,e.bodyAnchor,e.sphProp,e.native,e.bonds,e.bodyExt,e.bodyBond];this.bgFor={};for(const n of Object.keys(jt))this.bgFor[n]=this.d.createBindGroup({layout:this.pipes[n].getBindGroupLayout(0),entries:jt[n].map(s=>({binding:s,resource:{buffer:t[s]}}))})}growGrid(e=2){this.A.cap=Math.min(4096,Math.ceil(this.A.cap*e)),this.B.cap=Math.min(4096,Math.ceil(this.B.cap*e)),this.buf("cellItems",8*(this.A.cells*this.A.cap+this.B.cells*this.B.cap)),this.gridInfo.cap=this.A.cap,this.gridInfo.capB=this.B.cap,this.grows++,this.writeUniforms(),this.makeBindGroups()}grows=0;get uniState(){return{vox:this.wall.vox,wd:this.wall.wd}}writeUniforms(){const e=this.sys,t=this.set,n=this.A,s=this.B,r=this.wall.wd,a=new ArrayBuffer(256),i=new Float32Array(a),h=new Uint32Array(a);i[0]=e.lo[0],i[1]=e.lo[1],i[2]=e.lo[2],i[3]=t.dt,i[4]=e.hi[0],i[5]=e.hi[1],i[6]=e.hi[2],i[7]=t.eps,h[8]=n.dims[0],h[9]=n.dims[1],h[10]=n.dims[2],h[11]=n.cap,i[12]=n.h,i[13]=t.kRep,i[14]=4*t.kRep,i[15]=ws,h[16]=e.n,h[17]=this.nSph,h[18]=this.stepNo,h[19]=this.seed,i[20]=1/t.viscosity,i[21]=Wt,i[22]=2,i[23]=.15,i[24]=e.lo[0],i[25]=e.lo[1],i[26]=e.lo[2],i[27]=this.wall.vox,h[28]=r[0],h[29]=r[1],h[30]=r[2],h[31]=1,i[32]=t.hydro?Gt(e.phi):1,i[33]=t.hydro?Ot(e.phi):1,i[34]=e.phi,i[35]=this.rmax,i[36]=e.boundary==="periodic"?1:e.boundary==="embedded"?2:0,i[37]=t.tether,h[40]=s.dims[0],h[41]=s.dims[1],h[42]=s.dims[2],h[43]=s.cap,i[44]=s.h,i[45]=this.rSplit,i[46]=this.nLargeCore,i[47]=this.rmaxReg,h[48]=n.cells,h[49]=n.cells*n.cap,h[50]=n.cells+s.cells,h[51]=this.nLarge;const l=Pn(t.ionic);i[52]=l,i[53]=Mn,i[54]=Ln/l,i[55]=t.elec?1:0,i[56]=t.epsGo,i[57]=t.specific?1:0,i[58]=Math.exp(-Ln),this.d.queue.writeBuffer(this.uni,0,a)}setSettings(e){this.set=e,this.writeUniforms()}async buildWalls(e){const t=this.d,n=this.r,s=t.createShaderModule({code:As,label:"simWall"}),r=await s.getCompilationInfo();for(const u of r.messages.filter(p=>p.type==="error"))console.error(`[simWall] ${u.lineNum}:${u.linePos} ${u.message}`);const a=await t.createComputePipelineAsync({layout:"auto",compute:{module:s,entryPoint:"wall"}}),i=n.fillBufs,h=[i[0],i[1],i[2],i[3],i[4],n.tileState,i[5]],l=t.createBindGroup({layout:a.getBindGroupLayout(0),entries:[{binding:0,resource:{buffer:this.uni}}]}),c=a.getBindGroupLayout(1),d=[0,1,2,3,4,6],f=t.createBindGroup({layout:c,entries:d.map(u=>({binding:u,resource:{buffer:h[u]}}))}),w=t.createBindGroup({layout:a.getBindGroupLayout(2),entries:[{binding:0,resource:{buffer:this.bufs.wallGrid}}]}),b=t.createCommandEncoder(),m=b.beginComputePass();m.setPipeline(a),m.setBindGroup(0,l),m.setBindGroup(1,f),m.setBindGroup(2,w),m.dispatchWorkgroups(Math.ceil(e[0]*e[1]*e[2]/fe)),m.end(),t.queue.submit([b.finish()]),await t.queue.onSubmittedWorkDone()}plan(){const e=Math.ceil(this.nSph/fe),t=Math.ceil(this.sys.n/fe),n=Math.ceil((this.A.cells+this.B.cells)/fe),s=this.sys.n,r=Math.min(s,65535),a=Math.ceil(s/65535),i=[["clearGrid",n,1],["pose",e,1],["bin",e,1],["force",e,1]];return this.nBonds&&i.push(["bondForce",t,1]),i.push(["integrate",r,a],["tick",1,1]),i}step(e,t){const n=this.pipes,s=this.bgFor,r=this.plan(),a=e.beginComputePass();for(let i=0;i<t;i++)for(const[h,l,c]of r)a.setPipeline(n[h]),a.setBindGroup(0,s[h]),a.dispatchWorkgroups(l,c);a.end(),this.stepNo+=t,this.timeNs+=t*this.set.dt}async initDraw(){const e=this.d,t=this.r,n=this.sys;this.nBeads=n.beadLocal.length/4;const s=new Uint32Array(this.nBeads);for(let m=0;m<n.n;m++)s.fill(m,n.beadOff[m],n.beadOff[m]+n.beadCnt[m]);const r=new Float32Array(n.beadLocal),a=new Uint32Array(r.buffer);for(let m=0;m<this.nBeads;m++)a[4*m+3]=n.beadKind[m];this.buf("beadLocal",r),this.buf("beadBody",s);const i=new Uint32Array(4*n.n);for(let m=0;m<n.n;m++)i[4*m]=n.color[2*m],i[4*m+1]=n.color[2*m+1],i[4*m+2]=n.mobile[m],i[4*m+3]=n.type[m];this.buf("bodyCol",i);const h=e.createShaderModule({code:Ps,label:"simDraw"}),l=await h.getCompilationInfo();for(const m of l.messages.filter(u=>u.type==="error"))console.error(`[simDraw] ${m.lineNum}:${m.linePos} ${m.message}`);const c=m=>({binding:m,visibility:GPUShaderStage.VERTEX|(m===4?GPUShaderStage.FRAGMENT:0),buffer:{type:"read-only-storage"}}),d=e.createBindGroupLayout({entries:[c(0),c(1),c(2),c(3),c(4)]}),f=rn(t);this.drawG0=f.bg;const w=m=>({layout:e.createPipelineLayout({bindGroupLayouts:[f.layout,d]}),vertex:{module:h,entryPoint:m},fragment:{module:h,entryPoint:m==="vsBead"?"fsBead":"fsMain",targets:[{format:"rgba8unorm"},{format:"r32uint"}]},primitive:{topology:"triangle-list",cullMode:"none"},depthStencil:{format:"depth32float",depthWriteEnabled:!0,depthCompare:"greater"}});this.drawPipe=await e.createRenderPipelineAsync(w("vsBead"));const b=this.bufs;this.drawBG=e.createBindGroup({layout:d,entries:[b.bodyPos,b.bodyQuat,b.beadLocal,b.beadBody,b.bodyCol].map((m,u)=>({binding:u,resource:{buffer:m}}))}),this.nBonds&&(this.linkPipe=await e.createRenderPipelineAsync(w("vsLink")),this.linkBG=e.createBindGroup({layout:d,entries:[b.bodyPos,b.bodyQuat,b.bonds,b.beadBody,b.bodyCol].map((m,u)=>({binding:u,resource:{buffer:m}}))}))}setTracked(e){const t=this.sys,n=new Set(e),s=new Uint32Array(4*t.n);for(let r=0;r<t.n;r++)s[4*r]=t.color[2*r],s[4*r+1]=t.color[2*r+1],s[4*r+3]=t.type[r],s[4*r+2]=(t.mobile[r]?1:0)|(n.size?n.has(r)?4:2:0);this.d.queue.writeBuffer(this.bufs.bodyCol,0,s)}draw(e){e.setPipeline(this.drawPipe),e.setBindGroup(0,this.drawG0),e.setBindGroup(1,this.drawBG),e.draw(6,this.nBeads),this.linkPipe&&this.linkBG&&(e.setPipeline(this.linkPipe),e.setBindGroup(1,this.linkBG),e.draw(6,this.nBonds*Nt))}async gridLoad(){const e=this.A.cells+this.B.cells,t=this.d,n=t.createBuffer({size:4*e,usage:GPUBufferUsage.MAP_READ|GPUBufferUsage.COPY_DST}),s=t.createCommandEncoder();s.copyBufferToBuffer(this.bufs.cellCount,0,n,0,4*e),t.queue.submit([s.finish()]),await n.mapAsync(GPUMapMode.READ);const r=new Uint32Array(n.getMappedRange().slice(0));n.unmap(),n.destroy();let a=0,i=0;for(let h=0;h<e;h++){const l=h<this.A.cells?this.A.cap:this.B.cap,c=r[h];a=Math.max(a,c/l*this.A.cap),c>l&&(i+=c-l)}return{max:Math.round(a),cap:this.A.cap,dropped:i/Math.max(1,this.nSph)}}copyPoses(e,t){const n=this.sys.n;e.copyBufferToBuffer(this.bufs.bodyPos,0,t,0,16*n),e.copyBufferToBuffer(this.bufs.bodyQuat,0,t,16*n,16*n)}unpackPoses(e){const t=this.sys.n,n=new Float32Array(3*t);for(let s=0;s<t;s++)n[3*s]=e[4*s],n[3*s+1]=e[4*s+1],n[3*s+2]=e[4*s+2];return{pos:n,quat:e.slice(4*t,8*t)}}async readPoses(){const e=this.d,t=e.createBuffer({size:32*this.sys.n,usage:GPUBufferUsage.MAP_READ|GPUBufferUsage.COPY_DST}),n=e.createCommandEncoder();this.copyPoses(n,t),e.queue.submit([n.finish()]),await t.mapAsync(GPUMapMode.READ);const s=new Float32Array(t.getMappedRange().slice(0));return t.unmap(),t.destroy(),this.unpackPoses(s)}writePoses(e,t){const n=this.sys.n,s=new Float32Array(4*n);for(let r=0;r<n;r++)s[4*r]=e[3*r],s[4*r+1]=e[3*r+1],s[4*r+2]=e[3*r+2],s[4*r+3]=this.sys.dt[r];this.d.queue.writeBuffer(this.bufs.bodyPos,0,s),this.d.queue.writeBuffer(this.bufs.bodyQuat,0,t)}destroy(){for(const e of Object.values(this.bufs))e.destroy();this.uni.destroy()}}function rn(o){const e=o.device,t=GPUShaderStage.VERTEX|GPUShaderStage.FRAGMENT,n=e.createBindGroupLayout({entries:[{binding:0,visibility:t,buffer:{type:"uniform"}},{binding:5,visibility:t,buffer:{type:"read-only-storage"}}]}),s=e.createBindGroup({layout:n,entries:[{binding:0,resource:{buffer:o.uni}},{binding:5,resource:{buffer:o.palette}}]});return{layout:n,bg:s}}const Ps=_e.slice(0,_e.indexOf("@group(0) @binding(5)"))+`
@group(0) @binding(5) var<storage, read> palette: array<vec4f>;
`+(()=>{const o=Ee.indexOf("struct VOut"),e=Ee.indexOf("@group(1) @binding(0)");return Ee.slice(o,e)})()+`
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
  let k = ii / ${Nt}u; let t = (f32(ii % ${Nt}u) + 0.5) / ${Nt}.0;
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
`,Ss=new Set(["F_actin","Arp23","CP_cap","TPM1_coil","MT_tubulin","gTuRC_nuc","gTuRC_cyto","ATP_dimer","Pil1_dimer","Lsp1_dimer","Vac8_Nvj1","RVS161_RVS167_het","NPC","nucleosome"]),Vt=14,zn=21,Cs=22,Rs=23,We={protein:.36,nucleic:.6},Gn=new Map;function On(o,e){const t=Gn.get(e);if(t)return t;const n=o.types[e],s=o.beads,r=new Uint32Array(s.buffer),a=n.nb,i=n.bo,h=u=>(r[4*(i+u)+3]&255)===1?We.nucleic:We.protein,l=1.6,c=(u,p,x)=>`${Math.floor(u/l)},${Math.floor(p/l)},${Math.floor(x/l)}`,d=new Map;for(let u=0;u<a;u++){const p=4*(i+u),x=c(s[p],s[p+1],s[p+2]);(d.get(x)??d.set(x,[]).get(x)).push(u)}const f=[];for(let u=-1;u<=1;u++)for(let p=-1;p<=1;p++)for(let x=-1;x<=1;x++)if(u||p||x){const y=Math.hypot(u,p,x);f.push([u/y,p/y,x/y])}const w=[],b=We.protein;for(let u=0;u<a;u++){const p=4*(i+u),x=s[p],y=s[p+1],A=s[p+2],M=h(u)+b;let _=!1;for(const S of f){const g=x+S[0]*M,k=y+S[1]*M,P=A+S[2]*M;let T=!0;const E=Math.floor(g/l),I=Math.floor(k/l),F=Math.floor(P/l);for(let N=-1;N<=1&&T;N++)for(let G=-1;G<=1&&T;G++)for(let $=-1;$<=1&&T;$++){const Z=d.get(`${E+N},${I+G},${F+$}`);if(Z)for(const V of Z){if(V===u)continue;const ie=4*(i+V),me=h(V)+b-.1,be=s[ie]-g,ae=s[ie+1]-k,L=s[ie+2]-P;if(be*be+ae*ae+L*L<me*me){T=!1;break}}}if(T){_=!0;break}}_&&w.push(u)}const m=new Uint32Array(w);return Gn.set(e,m),m}const Ts=1;function Yt(o,e){const t=o.types[e],n=[];for(let s=0;s<t.sn;s++){const r=4*(t.so+s),a=o.st[r+3]-Ts;a>=.4&&n.push(o.st[r],o.st[r+1],o.st[r+2],a)}return new Float32Array(n)}let Un=null;function Mo(){return Un??=(async()=>{const o={charges:null,rhScale:new Map,memScale:new Array(8).fill(1),contacts:new Map,nPairs:0};try{const e=await Ae("library/rh_scale.json");for(const[t,n]of Object.entries(e.scale))o.rhScale.set(Number(t),n)}catch{}try{(await Ae("library/membrane_diffusion.json")).scale.forEach((t,n)=>{n<8&&t>0&&(o.memScale[n]=t)})}catch{}try{o.charges=new Int8Array(await He("library/charges.bin"))}catch{}try{const e=await Ae("library/contacts.json"),t=new Uint16Array(await He("library/contacts.bin"));for(const[n,s,r,a]of e.pairs){const i=t.subarray(2*r,2*(r+a)),h=new Uint16Array(a),l=new Uint16Array(a);for(let d=0;d<a;d++)h[d]=i[2*d],l[d]=i[2*d+1];const c=(d,f,w,b)=>{(o.contacts.get(d)??o.contacts.set(d,[]).get(d)).push({other:f,mine:w,theirs:b})};if(n===s){const d=new Uint16Array(2*a),f=new Uint16Array(2*a);d.set(h),d.set(l,a),f.set(l),f.set(h,a),c(n,n,d,f)}else c(n,s,h,l),c(s,n,l,h);o.nPairs++}}catch{}return o})(),Un}const $n=new Map,Bs=[4,2],Ns=.8;function xt(o,e,t,n,s,r,a,i){const h=Math.floor(t/r),l=Math.floor(n/r),c=Math.floor(s/r);let d=-1,f=1/0;for(let w=1;w<=4&&d<0;w++)for(let b=-w;b<=w;b++)for(let m=-w;m<=w;m++)for(let u=-w;u<=w;u++){const p=a.get(`${h+b},${l+m},${c+u}`);if(p)for(const x of p){const y=o[4*x]-t,A=o[4*x+1]-n,M=o[4*x+2]-s;let _=y*y+A*A+M*M;i&&(_/=i[x]*i[x]),_<f&&(f=_,d=x)}}if(d<0)for(let w=0;w<e;w++){const b=o[4*w]-t,m=o[4*w+1]-n,u=o[4*w+2]-s;let p=b*b+m*m+u*u;i&&(p/=i[w]*i[w]),p<f&&(f=p,d=w)}return d}function kt(o,e,t){const n=new Map;for(let s=0;s<e;s++){const r=`${Math.floor(o[4*s]/t)},${Math.floor(o[4*s+1]/t)},${Math.floor(o[4*s+2]/t)}`;(n.get(r)??n.set(r,[]).get(r)).push(s)}return n}function Es(o,e,t){const n=t+":"+e,s=$n.get(n);if(s)return s;const r=o.types[e],a=o.beads,i=new Uint32Array(a.buffer),h=r.nb,l=r.bo,c=m=>i[4*(l+m)+3]&255,d=new Int32Array(h);let f;if(t==="schematic"){f=o.st.slice(4*r.so,4*(r.so+r.sn));const m=new Float32Array(r.sn);for(let p=0;p<r.sn;p++)m[p]=f[4*p+3];const u=kt(f,r.sn,3);for(let p=0;p<h;p++)d[p]=xt(f,r.sn,a[4*(l+p)],a[4*(l+p)+1],a[4*(l+p)+2],3,u,m)}else if(t==="residue"){const m=On(o,e),u=Yt(o,e);f=new Float32Array(4*m.length+u.length),m.forEach((x,y)=>{const A=4*(l+x);f.set([a[A],a[A+1],a[A+2],c(x)===1?We.nucleic:We.protein],4*y)}),f.set(u,4*m.length),d.fill(-1),m.forEach((x,y)=>{d[x]=y});const p=kt(f,m.length,1.5);for(let x=0;x<h;x++)d[x]<0&&(d[x]=xt(f,m.length,a[4*(l+x)],a[4*(l+x)+1],a[4*(l+x)+2],1.5,p))}else{const m=[],u=new Int32Array(h).fill(-1);for(const S of[0,1]){const g=[];for(let N=0;N<h;N++)c(N)===S&&g.push(N);if(!g.length)continue;let k=Math.ceil(g.length/Bs[S]),P=g.map((N,G)=>Math.floor(G*k/g.length)),T=new Float32Array(4*k);for(let N=0;N<5;N++){const G=new Float64Array(4*k);g.forEach((V,ie)=>{const me=4*(l+V),be=4*P[ie];G[be]+=a[me],G[be+1]+=a[me+1],G[be+2]+=a[me+2],G[be+3]++});const $=[];for(let V=0;V<k;V++)G[4*V+3]>0&&$.push(V);T=new Float32Array(4*$.length),$.forEach((V,ie)=>{T[4*ie]=G[4*V]/G[4*V+3],T[4*ie+1]=G[4*V+1]/G[4*V+3],T[4*ie+2]=G[4*V+2]/G[4*V+3]}),k=$.length;const Z=kt(T,k,1.5);P=g.map(V=>xt(T,k,a[4*(l+V)],a[4*(l+V)+1],a[4*(l+V)+2],1.5,Z))}const E=new Float64Array(4*k),I=new Float64Array(k);g.forEach((N,G)=>{const $=4*(l+N),Z=4*P[G];E[Z]+=a[$],E[Z+1]+=a[$+1],E[Z+2]+=a[$+2],E[Z+3]++});const F=new Int32Array(k).fill(-1);for(let N=0;N<k;N++)if(E[4*N+3]>0){F[N]=m.length/4;const G=E[4*N+3];m.push(E[4*N]/G,E[4*N+1]/G,E[4*N+2]/G,0)}g.forEach((N,G)=>{const $=F[P[G]],Z=4*(l+N),V=a[Z]-m[4*$],ie=a[Z+1]-m[4*$+1],me=a[Z+2]-m[4*$+2];I[P[G]]+=V*V+ie*ie+me*me,u[N]=$});for(let N=0;N<k;N++)F[N]>=0&&(m[4*F[N]+3]=(S===1?We.nucleic:We.protein)+Ns*Math.sqrt(I[N]/E[4*N+3]))}const p=On(o,e),x=new Uint8Array(m.length/4),y=Yt(o,e);for(const S of p)x[u[S]]=1;const A=new Int32Array(m.length/4).fill(-1);let M=0;for(let S=0;S<m.length/4;S++)x[S]&&(A[S]=M++);f=new Float32Array(4*M+y.length);for(let S=0;S<m.length/4;S++)A[S]>=0&&f.set(m.slice(4*S,4*S+4),4*A[S]);f.set(y,4*M);const _=kt(f,M,1.5);for(let S=0;S<h;S++){const g=A[u[S]];d[S]=g>=0?g:xt(f,M,a[4*(l+S)],a[4*(l+S)+1],a[4*(l+S)+2],1.5,_)}}const w=t==="schematic"?0:Yt(o,e).length/4,b={local:f,beadTo:d,n:f.length/4,nReg:f.length/4-w};return $n.set(n,b),b}const qn=new Map;function _s(o,e,t,n,s){const r=t+":"+e,a=qn.get(r);if(a)return a;const i=o.types[e],h=new Uint32Array(o.beads.buffer),l=new Float32Array(n.n);for(let c=0;c<i.nb;c++){const d=s[i.bo+c];d&&(l[n.beadTo[c]]+=(h[4*(i.bo+c)+3]&255)===1?-Ft*d:d)}return qn.set(r,l),l}const Kt=o=>(o<<16>>16)/32767*.70710678;function Fs(o,e){const t=o[4*e],n=o[4*e+2],s=o[4*e+3],r=Kt(n>>>16),a=Kt(s&65535),i=Kt(s>>>16),h=Math.sqrt(Math.max(0,1-r*r-a*a-i*i)),l=t>>>14&3,c=l===0?[h,r,a,i]:l===1?[r,h,a,i]:l===2?[r,a,h,i]:[r,a,i,h];return[c[1],c[2],c[3],c[0]]}function Po(o,e){const t=e.size/2,n=[],s=e.center.map(a=>Math.floor((a-t)/200)),r=e.center.map(a=>Math.floor((a+t)/200));for(let a=s[0];a<=r[0];a++)for(let i=s[1];i<=r[1];i++)for(let h=s[2];h<=r[2];h++){const l=o.idx.key.get(_t(a,i,h));l!==void 0&&n.push(l)}return n}function Is(o,e,t,n,s,r){const a=t.size/2,i=t.center,h=[i[0]-a,i[1]-a,i[2]-a],l=[i[0]+a,i[1]+a,i[2]+a],c=new Uint32Array(e.typeInfo),d=[];for(const L of Po(o,t)){const z=o.loaded.get(L);if(z)for(let j=0;j<z.nInst;j++){const X=o.decode(z,j);X.p[0]<h[0]||X.p[1]<h[1]||X.p[2]<h[2]||X.p[0]>=l[0]||X.p[1]>=l[1]||X.p[2]>=l[2]||d.push({t:z,i:j,...X})}}const f=d.length,w=n.mode,b=new Map;for(const L of d)b.has(L.type)||b.set(L.type,Es(e,L.type,w));let m=0,u=0;for(const L of d)m+=b.get(L.type).n,u+=e.types[L.type].nb;const p={n:f,type:new Uint32Array(f),pos:new Float32Array(3*f),quat:new Float32Array(4*f),mobile:new Uint8Array(f),dt:new Float32Array(f),dr:new Float32Array(f),wall:new Uint8Array(f),sphOff:new Uint32Array(f),sphCnt:new Uint32Array(f),sphLocal:new Float32Array(4*m),sphProp:new Uint32Array(4*m),native:new Uint32Array(0),bondOff:new Uint32Array(f),bondCnt:new Uint32Array(f),bonds:new Float32Array(0),beadOff:new Uint32Array(f),beadCnt:new Uint32Array(f),beadLocal:new Float32Array(4*u),beadKind:new Uint32Array(u),color:new Uint32Array(2*f),lo:h,hi:l,phi:0,tether:new Uint8Array(f),mem2d:new Uint8Array(f),boundary:"walls"},x=d.some(L=>L.region>=Vt&&L.region<=zn),y=n.boundary==="auto"?x?"embedded":"periodic":n.boundary;p.boundary=y;let A=0;for(const L of d)A+=e.types[L.type].mass*uo;p.phi=Math.min(.5,A/t.size**3);const M=new Set(d.map(L=>L.type)),_=new Map,S=[];let g=0;const k=new Float32Array(1),P=new Uint32Array(k.buffer);if(r.contacts.size)for(const L of M){const z=(r.contacts.get(L)??[]).filter(ne=>M.has(ne.other));if(!z.length)continue;const j=b.get(L),X=new Map;for(const ne of z){const Re=b.get(ne.other),ge=new Set;for(let ye=0;ye<ne.mine.length;ye++){const ve=j.beadTo[ne.mine[ye]],xe=Re.beadTo[ne.theirs[ye]];ve>=0&&xe>=0&&ge.add(ve*65536+xe)}ne.other>=L&&g++;const Te=ne.other===L?new Set([...ge].map(ye=>{const ve=Math.floor(ye/65536),xe=ye%65536;return ve<xe?ve*65536+xe:xe*65536+ve})).size:ge.size;for(const ye of ge){const ve=Math.floor(ye/65536),xe=ye%65536;(X.get(ve)??X.set(ve,new Map).get(ve)).set(ne.other<<16|xe,1/Math.max(1,Te))}}const pe=new Uint32Array(j.n),se=new Uint32Array(j.n);for(const[ne,Re]of X){pe[ne]=S.length/2,se[ne]=Re.size;for(const[ge,Te]of Re)k[0]=Te,S.push(ge>>>0,P[0])}_.set(L,{off:pe,cnt:se})}p.native=new Uint32Array(S);const T=new Uint8Array(f),E=[],I=e.beads,F=new Uint32Array(I.buffer);let N=0,G=0,$=0,Z=0,V=0,ie=0,me=0,be=0;d.forEach((L,z)=>{const j=e.types[L.type];p.type[z]=L.type,T[z]=L.region,p.pos.set(L.p,3*z),p.quat.set(Fs(L.t.inst,L.i),4*z);const X=L.region>=Vt&&L.region<=zn,pe=(L.flags&2)!==0||Ss.has(j.id)||L.region===Cs||L.region===Rs,se=y==="embedded"&&L.p.some((le,we)=>le<h[we]+n.shell||le>l[we]-n.shell),ne=X&&!pe&&n.memFluid,Re=(!X||ne)&&!pe;p.mobile[z]=Re?1:0,p.wall[z]=Re&&!X?1:0,p.tether[z]=Re&&se?1:0,p.mem2d[z]=ne?1:0,ne&&V++,Re&&!se?$++:Re&&Z++;let ge=ho(e,L.type,s);const Te=r.rhScale.get(L.type);if(Te&&(ge={rh:ge.rh*Te,dt:ge.dt/Te,dr:ge.dr/(Te*Te*Te)}),E.push(ge),ne){const le=po(ge.rh,n.memVisc),we=n.memCalibrated?r.memScale[L.region-Vt]:1;p.dt[z]=le.dt*we,p.dr[z]=le.dr*we}else p.dt[z]=ge.dt,p.dr[z]=ge.dr;const ye=b.get(L.type),ve=r.charges?_s(e,L.type,w,ye,r.charges):null,xe=_.get(L.type);xe&&be++,p.sphOff[z]=N,p.sphLocal.set(ye.local,4*N);for(let le=0;le<ye.n;le++,N++){const we=ve?ve[le]:0;we&&(ie++,me+=we),k[0]=we,p.sphProp[4*N]=P[0],xe&&(p.sphProp[4*N+1]=xe.off[le],p.sphProp[4*N+2]=xe.cnt[le]),p.sphProp[4*N+3]=le>=ye.nReg?4294967294:le<65535?(L.type<<16|le)>>>0:4294967295}p.sphCnt[z]=N-p.sphOff[z],p.beadOff[z]=G;for(let le=0;le<j.nb;le++,G++){const we=4*(j.bo+le);p.beadLocal[4*G]=I[we],p.beadLocal[4*G+1]=I[we+1],p.beadLocal[4*G+2]=I[we+2],p.beadKind[G]=F[we+3]}p.beadCnt[z]=j.nb;const Fn=c[16*L.type+8],jo=c[16*L.type+9];p.color[2*z]=Fn||Bt(ee[L.region]??"#7EA7D3"),p.color[2*z+1]=jo});const ae=n.mrna?Ds(e,p):{nLinks:0,nPolyRibo:0};return{sys:p,meta:{types:Array.from(p.type),region:T,hydro:E,box:t,nMobile:$,nFixed:f-$-Z,nShell:Z,boundary:y,hasMembrane:x,nMem2d:V,charges:!!r.charges,nCharged:ie,netCharge:me,contacts:r.contacts.size>0,nNativePairs:g,nNativeBodies:be,...ae}}}function Ds(o,e){const t=o.byId.get("ribosome_80S"),n=[];if(t!==void 0)for(let u=0;u<e.n;u++)e.type[u]===t&&e.mobile[u]&&n.push(u);const s=go,r=(u,p)=>Ls(e.quat.subarray(4*u,4*u+4),p[0],p[1],p[2]),a=n.map(u=>{const p=r(u,s.exit);return[e.pos[3*u]+p[0],e.pos[3*u+1]+p[1],e.pos[3*u+2]+p[2]]}),i=n.map(u=>{const p=r(u,s.entry);return[e.pos[3*u]+p[0],e.pos[3*u+1]+p[1],e.pos[3*u+2]+p[2]]}),h=n.map(u=>r(u,s.exitDir)),l=n.map(u=>r(u,s.entryDir)),c=[];for(let u=0;u<n.length;u++)for(let p=0;p<n.length;p++){if(u===p)continue;const x=[i[p][0]-a[u][0],i[p][1]-a[u][1],i[p][2]-a[u][2]],y=Math.hypot(x[0],x[1],x[2]);if(y<s.linkMin-.05||y>s.linkMax+.05)continue;const A=x.map(M=>M/y);A[0]*h[u][0]+A[1]*h[u][1]+A[2]*h[u][2]<Math.cos(s.cone+.02)||A[0]*l[p][0]+A[1]*l[p][1]+A[2]*l[p][2]>-.995||c.push([n[u],n[p]])}const d=e.boundary==="periodic"?e.hi[0]-e.lo[0]:0,f=fo*bo,w=Array.from({length:e.n},()=>[]);for(const[u,p]of c)w[u].push(...s.exit,f,...s.entry,an,p,d,u,1),w[p].push(...s.entry,f,...s.exit,an,u,d,p,0);const b=[];for(let u=0;u<e.n;u++)e.bondOff[u]=b.length/12,e.bondCnt[u]=w[u].length/12,b.push(...w[u]);e.bonds=new Float32Array(b);const m=new Set;for(const[u,p]of c)m.add(u),m.add(p);return{nLinks:c.length,nPolyRibo:m.size}}function Ls(o,e,t,n){const s=o[0],r=o[1],a=o[2],i=o[3],h=2*(r*n-a*t),l=2*(a*e-s*n),c=2*(s*t-r*e);return[e+i*h+(r*c-a*l),t+i*l+(a*h-s*c),n+i*c+(s*l-r*h)]}const zs=`#!/usr/bin/env python
"""Run a box exported from the yeast whole-cell viewer in OpenMM, as a coarse-grained model.

Files (from the viewer's simulation mode, "Export MD files"):
  yeast_box.cif    one pseudo-atom per residue (CA) / nucleotide (P) bead, one chain per molecule; chain ZZZZ: membrane
                   midplane beads; occupancy 1 = mobile, 0 = fixed (membrane proteins, filaments, chromatin)
  yeast_box.json   manifest: molecule sources (PDB / AlphaFold DB), poses, hydrodynamic radii, settings of the browser run

Model (deliberately simple; edit to taste):
  - each molecule: elastic network over its beads (pairs closer than ENM_CUT, force constant ENM_K), so it stays close
    to its structure; inter-molecular: Weeks-Chandler-Andersen excluded volume (EPS_REP, sigma = bead diameters) plus
    the optional short-range attraction of the browser run (depth "eps" kT over 0.6 nm, cos^2 well)
  - fixed beads (occupancy 0, membrane beads) have zero mass, which OpenMM keeps in place
  - box faces as in the browser run: periodic boundary conditions, or soft walls keeping mobile beads inside the box;
    "embedded" boxes also tether the shell molecules (harmonic restraint of each bead to its start, the body's tether
    shared over its beads)
  - interactions of the browser run, when exported: screened electrostatics between charged beads (Debye-Hueckel,
    manifest "interactions.charges"), polysome mRNA links (worm-like chain between the beads at the mRNA exit of one
    ribosome and the entry of the next), native contacts of predicted complexes (cos^2 wells between residue pairs)
  - Langevin dynamics at 303 K after an energy minimisation

  python run_openmm.py [--steps 50000] [--dt 0.02] [--out traj.dcd] [--platform OpenCL]
Requires: openmm (pip install openmm), numpy.
"""
import argparse, json, sys
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as u

ENM_CUT, ENM_K = 1.2, 500.0          # nm, kJ/mol/nm^2
ENM_MIN = 0.1                         # nm; shorter pairs get no bond (bonds of near-zero length blew up on CUDA)
WCA_CORE = 0.8                        # below this fraction of sigma the excluded volume continues linearly (capped force)
EPS_REP = 2.5                         # kJ/mol (~1 kT)
SIGMA = {"CA": 0.50, "P": 0.80, "M": 3.0}   # nm, bead diameters (membrane midplane bead: most of the bilayer)
MASS = {"CA": 110.0, "P": 330.0, "M": 0.0}  # amu


def pairs_within(x, cut):
    """index pairs (i < j) closer than cut, with their distances (scipy k-d tree if available, else a cell list)"""
    try:
        from scipy.spatial import cKDTree
        p = cKDTree(x).query_pairs(cut, output_type="ndarray")
        return np.c_[p, np.linalg.norm(x[p[:, 0]] - x[p[:, 1]], axis=1)] if len(p) else np.zeros((0, 3))
    except ImportError:
        pass
    keys = np.floor(x / cut).astype(np.int64)
    cells = {}
    for i, k in enumerate(map(tuple, keys)):
        cells.setdefault(k, []).append(i)
    out = []
    for k, idx in cells.items():
        idx = np.array(idx)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    nb = cells.get((k[0] + dx, k[1] + dy, k[2] + dz))
                    if nb is None:
                        continue
                    nb = np.array(nb)
                    d = np.linalg.norm(x[idx][:, None] - x[nb][None], axis=2)
                    a, b = np.nonzero((d < cut) & (idx[:, None] < nb[None]))
                    out.append(np.c_[idx[a], nb[b], d[a, b]])
    return np.vstack(out) if out else np.zeros((0, 3))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif", default="yeast_box.cif"); ap.add_argument("--json", default="yeast_box.json")
    ap.add_argument("--steps", type=int, default=50000); ap.add_argument("--dt", type=float, default=0.02, help="ps")
    ap.add_argument("--out", default="yeast_box.dcd"); ap.add_argument("--every", type=int, default=500)
    ap.add_argument("--platform", default=None)
    a = ap.parse_args()
    meta = json.load(open(a.json))
    eps_att = float(meta["settings"].get("eps", 0.0)) * 2.519            # kT -> kJ/mol at 303 K
    cif = app.PDBxFile(a.cif)
    top = cif.topology
    x = np.array(cif.getPositions(asNumpy=True).value_in_unit(u.nanometer))
    atoms = list(top.atoms())
    kind = np.array(["M" if at.residue.chain.id == "ZZZZ" else at.name for at in atoms])
    chain = np.array([at.residue.chain.index for at in atoms])
    mobile_chain = {b["chain"]: b["mobile"] for b in meta["bodies"]}
    mobile = np.array([kind[i] != "M" and mobile_chain.get(at.residue.chain.id, False) for i, at in enumerate(atoms)])
    print(f"{len(atoms)} beads in {top.getNumChains()} chains; {mobile.sum()} mobile", flush=True)

    system = mm.System()
    for i in range(len(atoms)):
        system.addParticle(MASS[kind[i]] if mobile[i] else 0.0)
    # elastic network within each mobile molecule
    bonds = mm.HarmonicBondForce(); nb_ = 0
    for c in np.unique(chain[mobile]):
        idx = np.where((chain == c) & mobile)[0]
        if len(idx) < 2:
            continue
        p = pairs_within(x[idx], ENM_CUT)
        p = p[p[:, 2] >= ENM_MIN] if len(p) else p      # near-duplicate beads (overlapping RNA/protein beads of
        for i, j, d in p:                               # ribosome models) would get zero-length bonds
            bonds.addBond(int(idx[int(i)]), int(idx[int(j)]), float(d), ENM_K); nb_ += 1
    system.addForce(bonds)
    # excluded volume (+ attraction) between molecules
    # pairs inside one molecule are left to the elastic network; pairs of two fixed beads are switched off (m = mobile)
    inter = meta.get("interactions") or {}
    el = inter.get("electrostatics")
    q = np.zeros(len(atoms))
    for i, z in (inter.get("charges") or []) if el else []:
        q[int(i)] = z
    kappa = float(el["kappa_per_nm"]) if el else 1.0
    kt = 2.519                                                           # kJ/mol at 303 K
    expr = ("(1 - delta(c1 - c2)) * step(m1 + m2 - 0.5) * (wca + att + elec);"
            "wca = step(rc - r) * (4 * er * ((s / rw)^12 - (s / rw)^6) + er + fw * (rw - r));"
            "fw = 4 * er * (12 * (s / rw)^12 - 6 * (s / rw)^6) / rw; rw = max(r, core * s);"
            "att = -ea * step(r - s) * step(s + 0.6 - r) * cos(1.5707963 * (r - s) / 0.6)^2;"
            "elec = lb * q1 * q2 * (exp(-kappa * rr) / rr - exp(-kappa * rcut) / rcut) * step(rcut - r); rr = max(r, s);"
            "rc = s * 1.122462; s = 0.5 * (s1 + s2)")
    nbf = mm.CustomNonbondedForce(expr)
    nbf.addGlobalParameter("er", EPS_REP); nbf.addGlobalParameter("ea", eps_att); nbf.addGlobalParameter("core", WCA_CORE)
    nbf.addGlobalParameter("lb", kt * float(el["bjerrum_nm"]) if el else 0.0); nbf.addGlobalParameter("kappa", kappa)
    nbf.addPerParticleParameter("c"); nbf.addPerParticleParameter("s"); nbf.addPerParticleParameter("m"); nbf.addPerParticleParameter("q")
    for i in range(len(atoms)):
        nbf.addParticle([float(chain[i]), SIGMA[kind[i]], 1.0 if mobile[i] else 0.0, float(q[i])])
    lo, hi = np.array(meta["box"]["lo"]), np.array(meta["box"]["hi"])
    boundary = meta["box"].get("boundary", "walls")
    if boundary == "periodic":
        L = hi - lo
        system.setDefaultPeriodicBoxVectors(mm.Vec3(L[0], 0, 0), mm.Vec3(0, L[1], 0), mm.Vec3(0, 0, L[2]))
        nbf.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    else:
        nbf.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
    s_max = 0.5 * (max(SIGMA["CA"], SIGMA["P"]) + max(SIGMA.values()))    # mobile bead against anything
    rcut = max(s_max * 1.123, s_max + 0.6) + 0.05
    if el:
        rcut = max(rcut, 2 * 0.4 + 3.0 / kappa)                              # three Debye lengths beyond contact
    nbf.addGlobalParameter("rcut", rcut)
    nbf.setCutoffDistance(rcut)
    nbf.setUseLongRangeCorrection(False)
    mob_ix = [int(i) for i in np.where(mobile)[0]]
    system.addForce(nbf)
    # box walls (closed and embedded boxes)
    if boundary != "periodic":
        wall = mm.CustomExternalForce("kw * (max(0, xl - x)^2 + max(0, x - xh)^2 + max(0, yl - y)^2 + max(0, y - yh)^2 + "
                                      "max(0, zl - z)^2 + max(0, z - zh)^2)")
        wall.addGlobalParameter("kw", 500.0)
        for nm_, v in zip(("xl", "yl", "zl", "xh", "yh", "zh"), list(lo) + list(hi)):
            wall.addGlobalParameter(nm_, float(v))
        for i in mob_ix:
            wall.addParticle(i, [])
        system.addForce(wall)
    # embedded boundary: shell molecules tethered to where they started
    tk = float(meta["box"].get("tether_kT_per_nm2", 0.0)) * 2.519        # kT/nm^2 -> kJ/mol/nm^2
    teth = {b["chain"] for b in meta["bodies"] if b.get("tethered")}
    if tk > 0 and teth:
        tf = mm.CustomExternalForce("0.5 * kt * ((x - x0)^2 + (y - y0)^2 + (z - z0)^2)")
        tf.addPerParticleParameter("kt"); tf.addPerParticleParameter("x0"); tf.addPerParticleParameter("y0"); tf.addPerParticleParameter("z0")
        nb_body = {}
        for at in atoms:
            nb_body[at.residue.chain.id] = nb_body.get(at.residue.chain.id, 0) + 1
        for i, at in enumerate(atoms):
            c = at.residue.chain.id
            if c in teth and mobile[i]:
                tf.addParticle(i, [tk / nb_body[c], *map(float, x[i])])
        system.addForce(tf)
    # polysome mRNA links: worm-like chain (Marko-Siggia), extension capped at 0.95 of the contour
    links = inter.get("mrna_links") or []
    if links:
        wlc = mm.CustomBondForce("kt * lc / lp * (0.5 * u^2 - 0.25 * u + 0.25 / (1 - u) - 0.25); u = min(r / lc, 0.95)")
        wlc.addGlobalParameter("kt", kt); wlc.addPerBondParameter("lc"); wlc.addPerBondParameter("lp")
        wlc.setUsesPeriodicBoundaryConditions(boundary == "periodic")
        for l in links:
            wlc.addBond(int(l["atoms"][0]), int(l["atoms"][1]), [float(l["contour_nm"]), float(l["persistence_nm"])])
        system.addForce(wlc)
    # native contacts of predicted complexes: each residue pair a cos^2 well of its share of the interface energy
    nat = inter.get("native_contacts") or []
    if nat:
        e_if = float(inter["native_contact_model"]["interface_energy_kT"]) * kt
        go = mm.CustomBondForce("-e * select(step(0.5 - r), 1, step(1.1 - r) * cos(1.5707963 * (r - 0.5) / 0.6)^2)")
        go.addPerBondParameter("e"); go.setUsesPeriodicBoundaryConditions(boundary == "periodic")
        for i, j, share in nat:
            if mobile[int(i)] or mobile[int(j)]:
                go.addBond(int(i), int(j), [e_if * float(share)])
        system.addForce(go)
    print(f"elastic-network bonds {nb_}; attraction {eps_att:.2f} kJ/mol; boundary {boundary}; charged beads {int((q != 0).sum())}; "
          f"mRNA links {len(links)}; native contact pairs {len(nat)}", flush=True)

    integ = mm.LangevinMiddleIntegrator(meta.get("temperature_K", 303.15) * u.kelvin, 1.0 / u.picosecond, a.dt * u.picoseconds)
    plat = mm.Platform.getPlatformByName(a.platform) if a.platform else None
    sim = app.Simulation(top, system, integ, plat) if plat else app.Simulation(top, system, integ)
    sim.context.setPositions(x * u.nanometer)
    print("platform", sim.context.getPlatform().getName(), flush=True)
    sim.minimizeEnergy(maxIterations=500)
    sim.reporters.append(app.DCDReporter(a.out, a.every))
    sim.reporters.append(app.StateDataReporter(sys.stdout, a.every * 10, step=True, time=True, potentialEnergy=True,
                                               temperature=True, speed=True))
    sim.step(a.steps)
    print("done:", a.out)


if __name__ == "__main__":
    main()
`,Gs=`#!/usr/bin/env python
"""Expand the rigid-body trajectory of a browser simulation into a bead-level DCD (VMD, ChimeraX, MDAnalysis, mdtraj).

Files (from the viewer's simulation mode, "Export MD files"): yeast_box.cif, yeast_box.json, yeast_box_trajectory.json.
Each molecule's beads are moved rigidly with its recorded pose (position, quaternion x,y,z,w); membrane beads stay put.

  python trajectory_to_dcd.py [--out yeast_box_browser.dcd]
  then e.g.  vmd yeast_box.cif yeast_box_browser.dcd   or   ChimeraX: open yeast_box.cif; open yeast_box_browser.dcd
Requires: openmm (for the mmCIF reader and DCD writer), numpy.
"""
import argparse, json
import numpy as np
import openmm.app as app
import openmm.unit as u


def rotmat(q):
    """rotation matrices from quaternions (..., 4) as x, y, z, w"""
    x, y, z, w = q[..., 0], q[..., 1], q[..., 2], q[..., 3]
    return np.stack([np.stack([1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)], -1),
                     np.stack([2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)], -1),
                     np.stack([2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)], -1)], -2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif", default="yeast_box.cif"); ap.add_argument("--json", default="yeast_box.json")
    ap.add_argument("--traj", default="yeast_box_trajectory.json"); ap.add_argument("--out", default="yeast_box_browser.dcd")
    a = ap.parse_args()
    cif = app.PDBxFile(a.cif)
    x = np.array(cif.getPositions(asNumpy=True).value_in_unit(u.nanometer))
    meta = json.load(open(a.json)); tr = json.load(open(a.traj))["frames"]
    bodies = meta["bodies"]
    chain_of = np.array([at.residue.chain.id for at in cif.topology.atoms()])
    # bead -> body (manifest order = chain order); local coordinates from the exported pose
    idx = {b["chain"]: i for i, b in enumerate(bodies)}
    body = np.array([idx.get(c, -1) for c in chain_of])
    P = np.array([b["pos_nm"] for b in bodies]); Q = np.array([b["quat_xyzw"] for b in bodies])
    R0 = rotmat(Q)
    m = body >= 0
    local = np.zeros_like(x)
    local[m] = np.einsum("nji,nj->ni", R0[body[m]], x[m] - P[body[m]])      # R^T (x - p)
    lo, hi = np.array(meta["box"]["lo"]), np.array(meta["box"]["hi"])
    dt_ps = (tr[1]["t_ns"] - tr[0]["t_ns"]) * 1e3 if len(tr) > 1 else 1.0
    with open(a.out, "wb") as fh:
        dcd = app.DCDFile(fh, cif.topology, dt_ps * u.picoseconds)
        for f in tr:
            pf = np.array(f["pos"]).reshape(-1, 3); qf = np.array(f["quat"]).reshape(-1, 4)
            Rf = rotmat(qf)
            y = x.copy()
            y[m] = pf[body[m]] + np.einsum("nij,nj->ni", Rf[body[m]], local[m])
            dcd.writeModel(y * u.nanometer, periodicBoxVectors=None)
    print(f"{len(tr)} frames, {len(x)} beads ({m.sum()} in {len(bodies)} molecules), "
          f"{tr[0]['t_ns']:.1f}-{tr[-1]['t_ns']:.1f} ns -> {a.out}; box {lo} to {hi} nm")


if __name__ == "__main__":
    main()
`;function At(o){return o.toString(36).toUpperCase().padStart(2,"0")}function rt(o,e){const t=typeof e=="string"?new Blob([e],{type:"text/plain"}):e,n=document.createElement("a");n.href=URL.createObjectURL(t),n.download=o,document.body.appendChild(n),n.click(),n.remove(),setTimeout(()=>URL.revokeObjectURL(n.href),4e3)}function Os(o,e,t){const n=o[e],s=o[e+1],r=o[e+2],a=o[e+3],i=2*(s*t[2]-r*t[1]),h=2*(r*t[0]-n*t[2]),l=2*(n*t[1]-s*t[0]);return[t[0]+a*i+(s*l-r*h),t[1]+a*h+(r*i-n*l),t[2]+a*l+(n*h-s*i)]}async function Us(o){const e=o.d,t=o.uniState,n=t.wd[0]*t.wd[1]*t.wd[2],s=e.createBuffer({size:16*n,usage:GPUBufferUsage.MAP_READ|GPUBufferUsage.COPY_DST}),r=e.createCommandEncoder();r.copyBufferToBuffer(o.bufs.wallGrid,0,s,0,16*n),e.queue.submit([r.finish()]),await s.mapAsync(GPUMapMode.READ);const a=new Float32Array(s.getMappedRange().slice(0));s.unmap(),s.destroy();const i=[],h=t.vox,l=o.sys.lo;for(let c=0;c<t.wd[0];c++)for(let d=0;d<t.wd[1];d++)for(let f=0;f<t.wd[2];f++){const w=(c*t.wd[1]+d)*t.wd[2]+f,b=a[4*w];if(Math.abs(b)>.5*h)continue;const m=[a[4*w+1],a[4*w+2],a[4*w+3]],u=Math.hypot(m[0],m[1],m[2])||1,p=[l[0]+(c+.5)*h,l[1]+(d+.5)*h,l[2]+(f+.5)*h];i.push(p[0]-b*m[0]/u,p[1]-b*m[1]/u,p[2]-b*m[2]/u)}return new Float32Array(i)}async function $s(o,e,t,n,s,r,a){const i=o.sys,h=i.n,l=["data_yeast_cell_box","#",`_cell.length_a ${(10*(i.hi[0]-i.lo[0])).toFixed(3)}`,`_cell.length_b ${(10*(i.hi[1]-i.lo[1])).toFixed(3)}`,`_cell.length_c ${(10*(i.hi[2]-i.lo[2])).toFixed(3)}`,"_cell.angle_alpha 90","_cell.angle_beta 90","_cell.angle_gamma 90","#","loop_","_atom_site.group_PDB","_atom_site.id","_atom_site.type_symbol","_atom_site.label_atom_id","_atom_site.label_comp_id","_atom_site.label_asym_id","_atom_site.label_entity_id","_atom_site.label_seq_id","_atom_site.Cartn_x","_atom_site.Cartn_y","_atom_site.Cartn_z","_atom_site.occupancy","_atom_site.B_iso_or_equiv","_atom_site.auth_seq_id","_atom_site.auth_asym_id","_atom_site.pdbx_PDB_model_num"];let c=1;const d=i.beadLocal,f=i.beadKind;for(let g=0;g<h;g++){const k=At(g),P=i.mobile[g]?"1.00":"0.00";for(let T=0;T<i.beadCnt[g];T++){const E=i.beadOff[g]+T,I=(f[E]&255)===1,F=Os(r,4*g,[d[4*E],d[4*E+1],d[4*E+2]]),N=10*(s[3*g]+F[0]),G=10*(s[3*g+1]+F[1]),$=10*(s[3*g+2]+F[2]);l.push(`ATOM ${c++} ${I?"P":"C"} ${I?"P":"CA"} ${I?"N":"UNK"} ${k} 1 ${T+1} ${N.toFixed(2)} ${G.toFixed(2)} ${$.toFixed(2)} ${P} 0 ${T+1} ${k} 1`)}}const w=await Mo(),b=new Uint32Array(h);{let g=0;for(let k=0;k<h;k++)b[k]=g,g+=i.beadCnt[k]}const m=[];if(w.charges)for(let g=0;g<h;g++){const k=t.types[i.type[g]];for(let P=0;P<i.beadCnt[g];P++){const T=w.charges[k.bo+P];T&&m.push([b[g]+P,(f[i.beadOff[g]+P]&255)===1?-Ft*T:T])}}const u=(g,k,P,T)=>{let E=0,I=1/0;for(let F=0;F<i.beadCnt[g];F++){const N=4*(i.beadOff[g]+F),G=(d[N]-k)**2+(d[N+1]-P)**2+(d[N+2]-T)**2;G<I&&(I=G,E=F)}return b[g]+E},p=[];for(let g=0;g<h;g++)for(let k=i.bondOff[g];k<i.bondOff[g]+i.bondCnt[g];k++){const P=i.bonds.subarray(12*k,12*k+12);if(P[11]!==1)continue;const T=P[8];p.push({chains:[At(g),At(T)],atoms:[u(g,P[0],P[1],P[2]),u(T,P[4],P[5],P[6])],anchors_local_nm:[[P[0],P[1],P[2]].map(E=>+E.toFixed(3)),[P[4],P[5],P[6]].map(E=>+E.toFixed(3))],contour_nm:+P[3].toFixed(2),persistence_nm:P[7]})}const x=[];if(w.contacts.size){const g=new Map;for(let k=0;k<h;k++)i.mobile[k]&&(g.get(i.type[k])??g.set(i.type[k],[]).get(i.type[k])).push(k);for(const[k,P]of g)for(const T of w.contacts.get(k)??[]){if(T.other<k)continue;const E=g.get(T.other);if(!E)continue;const I=T.other===k?T.mine.length/2:T.mine.length,F=1/I;for(const N of P)for(const G of E)if(!(N===G||T.other===k&&G<N))for(let $=0;$<I&&x.length<3e5;$++)x.push([b[N]+T.mine[$],b[G]+T.theirs[$],+F.toFixed(5)])}}const y=await Us(o);for(let g=0;g<y.length/3;g++)l.push(`ATOM ${c++} C M MEM ZZZZ 2 ${g+1} ${(10*y[3*g]).toFixed(2)} ${(10*y[3*g+1]).toFixed(2)} ${(10*y[3*g+2]).toFixed(2)} 0.00 0 ${g+1} ZZZZ 1`);l.push("#");const A=t.types,M={format:"yeast-cell box export v1",units:{length:"nm",time:"ns",diffusion:"um^2/s (in water at the given temperature; divide by the viscosity factor)"},temperature_K:xn,water_viscosity_Pa_s:yt,hydration_shell_nm:kn,box:{lo:i.lo,hi:i.hi,boundary:i.boundary,shell_nm:i.boundary==="embedded"?n.shell:0,tether_kT_per_nm2:i.boundary==="embedded"?n.tether:0},settings:n,volume_fraction_hydrated:+i.phi.toFixed(4),mean_field_hydrodynamics:n.hydro?{D_t_factor:+Gt(i.phi).toFixed(4),D_r_factor:+Ot(i.phi).toFixed(4),model:"short-time self-diffusion of hard spheres at the box volume fraction, Cichocki et al. 1999 J Chem Phys 111:3265"}:null,note:`Coarse-grained export of the yeast whole-cell model. Each molecule is one chain of pseudo-atoms, one per residue (CA) or nucleotide (P) bead of its structure, posed as in the cell (or the selected simulation frame). Chain ZZZZ holds fixed membrane surface beads (bilayer midplane). Mobile molecules have occupancy 1, fixed ones (membrane proteins, filaments, chromatin) 0. box.boundary: periodic (molecule centres wrapped into the box), walls, or embedded (molecules flagged tethered are held near their position by a harmonic tether of box.tether_kT_per_nm2 on their centre). To build an atomistic or Martini system, take each molecule's structure from "source" and apply its rotation (quaternion x,y,z,w) and translation (nm) to the structure's bead frame (CA/P centred as in the library).`,bodies:Array.from({length:h},(g,k)=>{const P=A[i.type[k]],T=e.hydro[k];return{chain:At(k),type:P.id,name:P.name,genes:P.genes,uniprot:P.uniprot||void 0,source:P.src==="PDB"?`PDB ${P.pdb}`:P.model?`AlphaFold DB ${P.model}`:P.uniprot?`AlphaFold DB AF-${P.uniprot}-F1`:P.src,region:e.region[k],mobile:!!i.mobile[k],tethered:!!i.tether[k],beads:i.beadCnt[k],R_h_nm:+T.rh.toFixed(3),D_t_um2_s:+(T.dt*1e3).toFixed(3),D_r_per_us:+(T.dr*1e3).toFixed(4),pos_nm:Array.from(s.subarray(3*k,3*k+3),E=>+E.toFixed(3)),quat_xyzw:Array.from(r.subarray(4*k,4*k+4),E=>+E.toFixed(6))}}),membrane_beads:y.length/3,interactions:{electrostatics:n.elec&&w.charges?{model:"Debye-Hueckel between charged beads (kT units)",bjerrum_nm:Mn,ionic_strength_M:n.ionic,kappa_per_nm:+Pn(n.ionic).toFixed(4),nucleotide_charge:Ft,note:"charges: [atom index (0-based, CIF order), charge in e]; Lys/Arg +1, Asp/Glu -1 (pH 7.2), nucleotides scaled for counterion condensation"}:null,charges:n.elec?m:[],mrna_links:n.mrna?p:[],mrna_link_model:"worm-like chain (Marko-Siggia) between the beads nearest the mRNA exit of one ribosome and the entry of the next",native_contacts:n.specific?x:[],native_contact_model:{interface_energy_kT:n.epsGo,note:"[atom i, atom j, share]: residue contacts of a confident AlphaFold DB dimer between two molecules of the partner types; each pair is a well of depth share x interface_energy_kT (cos^2, 0.6 nm wide, beyond contact)"}}},_=Math.max(1,Math.ceil(a.length/100)),S={note:"Rigid-body trajectory of the browser simulation: per frame, time (ns) and for every molecule (manifest order) position (nm) and quaternion (x,y,z,w).",frames:a.filter((g,k)=>k%_===0).map(g=>({t_ns:+g.tNs.toFixed(3),pos:Array.from(g.pos,k=>+k.toFixed(2)),quat:Array.from(g.quat,k=>+k.toFixed(4))}))};rt("yeast_box.cif",l.join(`
`)+`
`),rt("yeast_box.json",JSON.stringify(M)),rt("yeast_box_trajectory.json",JSON.stringify(S)),rt("run_openmm.py",zs),rt("trajectory_to_dcd.py",Gs)}const qs=_e.slice(0,_e.indexOf("@group(0) @binding(5)"))+Ee.slice(Ee.indexOf("struct VOut"),Ee.indexOf("@group(1) @binding(0)"))+`
@group(1) @binding(0) var<uniform> box: vec4f;      // centre xyz, half size
const PER: u32 = 480u;
@vertex
fn vsEdge(@builtin(vertex_index) vi: u32, @builtin(instance_index) ii: u32) -> VOut {
  let e = ii / PER; let t = f32(ii % PER) / f32(PER - 1u) * 2.0 - 1.0;
  // 12 edges: axis a = e / 4, the other two coordinates at +-1 from the edge number
  let a = e / 4u; let s1 = select(-1.0, 1.0, (e & 1u) == 1u); let s2 = select(-1.0, 1.0, (e & 2u) == 2u);
  var v = vec3f(0.0);
  if (a == 0u) { v = vec3f(t, s1, s2); } else if (a == 1u) { v = vec3f(s1, t, s2); } else { v = vec3f(s1, s2, t); }
  let p = box.xyz + v * box.w;
  // thickness in pixels rather than nm, so the outline reads at any zoom
  let pe = p - u.eye.xyz;
  return emitQuad(pe, max(0.3, 1.6 * length(pe) / u.eye.w), vi, vec3f(1.0, 0.55, 0.1), 0u);
}
// outline on top of everything: flat disc test, no cut-plane test (the slab would hide most of the box)
@fragment
fn fsMain(in: VOut) -> FOut {
  var o: FOut;
  let rd = normalize(in.vpos); let b = dot(rd, in.center);
  if (b * b - (dot(in.center, in.center) - in.radius * in.radius) < 0.0) { discard; }
  o.color = vec4f(in.color, 1.0); o.id = 0u; o.depth = 1.0;
  return o;
}
`,Ze=640,Hn=12,Wn=8,So=["#e6194b","#3cb44b","#4363d8","#f58231","#911eb4","#42d4f4","#f032e6","#9a6324"],Hs=_e.slice(0,_e.indexOf("@group(0) @binding(5)"))+Ee.slice(Ee.indexOf("struct VOut"),Ee.indexOf("@group(1) @binding(0)"))+`
@group(1) @binding(0) var<uniform> pu: vec4u;                 // segments drawn per molecule, stride, beads per segment, molecules
@group(1) @binding(1) var<storage, read> pts: array<vec4f>;   // per molecule, per frame: centre, w = periodic box length (0: none)
const COLS = array<vec3f, 8>(${So.map(o=>`vec3f(${[1,3,5].map(e=>(parseInt(o.slice(e,e+2),16)/255).toFixed(3)).join(", ")})`).join(", ")});
@vertex
fn vsPath(@builtin(vertex_index) vi: u32, @builtin(instance_index) ii: u32) -> VOut {
  let S = pu.z; let nSeg = pu.x;
  let k = ii / (S * nSeg); let rem = ii % (S * nSeg); let j = rem / S;
  let t = (f32(rem % S) + 0.5) / f32(S);
  let a = pts[k * pu.y + j]; let b = pts[k * pu.y + j + 1u];
  let dv = b.xyz - a.xyz;
  let p = a.xyz + t * dv; let pe = p - u.eye.xyz;
  var r = max(0.06, 1.6 * length(pe) / u.eye.w);
  if (a.w > 0.0 && any(abs(dv) > vec3f(0.5 * a.w))) { r = 0.0; }
  return emitQuad(pe, r, vi, COLS[k % 8u], 0u);
}
@fragment
fn fsPath(in: VOut) -> FOut { var o = shadeFrag(in); o.id = 0u; return o; }
`;class Ws{R;T;lib;box=null;set=ko();eng=null;meta=null;running=!1;active=!1;replayIdx=-1;frames=[];recording=!0;maxFrames=600;stepsSinceRec=0;recEvery=200;computing=!1;abort=!1;aheadRate=0;playing=!1;playFps=24;playStart=0;playFrom=0;hydroCache=new Map;boxPipe=null;boxBG=null;boxG0=null;boxBuf;onUpdate=()=>{};stepsPerSec=0;lastT=0;lastSteps=0;totalSteps=0;tracked=[];trackPipe=null;trackBG=null;trackUni;trackPts;constructor(e,t,n){this.R=e,this.T=t,this.lib=n,this.boxBuf=e.device.createBuffer({size:16,usage:GPUBufferUsage.UNIFORM|GPUBufferUsage.COPY_DST}),this.trackUni=e.device.createBuffer({size:16,usage:GPUBufferUsage.UNIFORM|GPUBufferUsage.COPY_DST}),this.trackPts=e.device.createBuffer({size:16*Ze*Wn,usage:GPUBufferUsage.STORAGE|GPUBufferUsage.COPY_DST})}setTracked(e){this.tracked=[...new Set(e)].slice(0,Wn),this.eng?.setTracked(this.tracked),this.uploadTracks(),this.onUpdate()}toggleTrack(e){!this.eng||e<0||e>=this.eng.sys.n||this.setTracked(this.tracked.includes(e)?this.tracked.filter(t=>t!==e):[...this.tracked,e])}uploadTracks(){const e=this.eng;if(!e||!this.tracked.length)return;const t=Math.min(this.frames.length,Ze),n=e.sys,s=n.boundary==="periodic"?n.hi[0]-n.lo[0]:0,r=new Float32Array(4*Ze*this.tracked.length);this.tracked.forEach((a,i)=>{for(let h=0;h<t;h++){const l=this.frames[h].pos,c=4*(i*Ze+h);r[c]=l[3*a],r[c+1]=l[3*a+1],r[c+2]=l[3*a+2],r[c+3]=s}}),this.R.device.queue.writeBuffer(this.trackPts,0,r)}async ensureTrackPipe(){if(this.trackPipe)return;const e=this.R.device,t=e.createShaderModule({code:Hs,label:"simTrack"}),n=await t.getCompilationInfo();for(const a of n.messages.filter(i=>i.type==="error"))console.error(`[simTrack] ${a.lineNum}:${a.linePos} ${a.message}`);const s=e.createBindGroupLayout({entries:[{binding:0,visibility:GPUShaderStage.VERTEX,buffer:{type:"uniform"}},{binding:1,visibility:GPUShaderStage.VERTEX,buffer:{type:"read-only-storage"}}]}),r=rn(this.R);this.trackPipe=await e.createRenderPipelineAsync({layout:e.createPipelineLayout({bindGroupLayouts:[r.layout,s]}),vertex:{module:t,entryPoint:"vsPath"},fragment:{module:t,entryPoint:"fsPath",targets:[{format:"rgba8unorm"},{format:"r32uint"}]},primitive:{topology:"triangle-list",cullMode:"none"},depthStencil:{format:"depth32float",depthWriteEnabled:!1,depthCompare:"always"}}),this.trackBG=e.createBindGroup({layout:s,entries:[{binding:0,resource:{buffer:this.trackUni}},{binding:1,resource:{buffer:this.trackPts}}]}),this.trackG0=r.bg}trackG0=null;trackStats(){const e=this.frames,t=this.replayIdx>=0?this.replayIdx:e.length-1;return this.tracked.map(n=>{if(t<1)return{b:n,disp:0,d:NaN};const s=l=>l.upos??l.pos;let r=0,a=0;for(let l=1;l<=t;l++){const c=s(e[l-1]),d=s(e[l]);r+=(d[3*n]-c[3*n])**2+(d[3*n+1]-c[3*n+1])**2+(d[3*n+2]-c[3*n+2])**2,a+=e[l].tNs-e[l-1].tNs}const i=s(e[0]),h=s(e[t]);return{b:n,disp:Math.hypot(h[3*n]-i[3*n],h[3*n+1]-i[3*n+1],h[3*n+2]-i[3*n+2]),d:a>0?r/(6*a)*1e3:NaN}})}setPlayFps(e){this.playing&&this.replayIdx>=0&&(this.playFrom=this.replayIdx,this.playStart=performance.now()),this.playFps=Math.max(.5,e)}async showBox(e,t){if(this.box={center:[...e],size:t},this.R.device.queue.writeBuffer(this.boxBuf,0,new Float32Array([...e,t/2])),!this.boxPipe){const n=this.R.device,s=n.createShaderModule({code:qs,label:"simBox"}),r=n.createBindGroupLayout({entries:[{binding:0,visibility:GPUShaderStage.VERTEX,buffer:{type:"uniform"}}]}),a=rn(this.R);this.boxG0=a.bg,this.boxPipe=await n.createRenderPipelineAsync({layout:n.createPipelineLayout({bindGroupLayouts:[a.layout,r]}),vertex:{module:s,entryPoint:"vsEdge"},fragment:{module:s,entryPoint:"fsMain",targets:[{format:"rgba8unorm"},{format:"r32uint",writeMask:0}]},primitive:{topology:"triangle-list",cullMode:"none"},depthStencil:{format:"depth32float",depthWriteEnabled:!1,depthCompare:"always"}}),this.boxBG=n.createBindGroup({layout:r,entries:[{binding:0,resource:{buffer:this.boxBuf}}]})}this.R.sim=this}hideBox(){this.box=null,this.active||(this.R.sim=null)}async start(e){if(!this.box)throw new Error("no box");e("Loading the tiles of the box");try{await this.T.ensureLoaded(Po(this.T,this.box))}catch(r){throw this.T.unpin(),r}e("Loading charges and native contacts");const t=await Mo();e("Building rigid bodies and diffusion coefficients");const{sys:n,meta:s}=Is(this.T,this.lib,this.box,this.set,this.hydroCache,t);if(!n.n)throw this.T.unpin(),new Error("the box is empty");e(`Compiling the simulation (${n.n} molecules, ${n.sphLocal.length/4} collision spheres)`),this.eng?.destroy(),this.eng=new Ms(this.R,n,this.set),await this.eng.init(),this.meta=s,this.frames=[],this.replayIdx=-1,this.totalSteps=0,this.stepsSinceRec=0,this.grid=null,this.tracked=[],await this.ensureTrackPipe(),this.recEvery=this.set.recordEvery,this.playing=!1,this.maxFrames=Math.max(100,Math.min(600,Math.floor(6e7/(28*n.n)))),this.active=!0,this.running=!0,this.R.sim=this,this.frames.push({tNs:0,pos:n.pos.slice(),quat:n.quat.slice()}),this.lastT=performance.now(),this.lastSteps=0}stop(){this.running=!1,this.active=!1,this.playing=!1,this.abort=!0,this.T.unpin(),this.eng?.destroy(),this.eng=null,this.meta=null,this.frames=[],this.R.sim=this.box?this:null}applySettings(e){this.set={...this.set,...e},this.set.dt=Math.min(this.set.dt,Cn(this.set.mode,this.set.eps)),this.eng?.setSettings(this.set)}grid=null;polling=!1;pollGrid(){if(!this.eng||this.polling)return;this.polling=!0;const e=this.eng;e.gridLoad().then(t=>{e===this.eng&&(this.grid=t,t.dropped>0&&e.growGrid())}).catch(()=>{}).finally(()=>{this.polling=!1})}clipBox(){return this.active&&this.box?[...this.box.center,this.box.size/2]:null}compute(e){if(!this.eng)return;if(this.playing){this.advancePlayback();return}if(!this.running||this.replayIdx>=0||this.computing)return;const t=Math.max(1,Math.min(512,Math.round(this.set.stepsPerFrame)));this.eng.step(e,t),this.totalSteps+=t,this.stepsSinceRec+=t;const n=performance.now();if(n-this.lastT>1e3&&(this.stepsPerSec=(this.totalSteps-this.lastSteps)*1e3/(n-this.lastT),this.lastT=n,this.lastSteps=this.totalSteps),this.recording&&this.stepsSinceRec>=this.recEvery&&!this.reading){this.stepsSinceRec=0,this.reading=!0;const s=this.eng.timeNs,r=this.eng;queueMicrotask(async()=>{try{const a=await r.readPoses();r===this.eng&&this.record({tNs:s,...a})}catch{}this.reading=!1})}}reading=!1;record(e){const t=this.frames[this.frames.length-1],n=this.eng?.sys;if(t&&n&&n.boundary==="periodic"){const s=n.hi[0]-n.lo[0],r=t.upos??t.pos,a=new Float32Array(e.pos.length);for(let i=0;i<a.length;i++){let h=e.pos[i]-t.pos[i];h-=s*Math.round(h/s),a[i]=r[i]+h}e.upos=a}this.frames.push(e),this.frames.length>this.maxFrames&&(this.frames=this.frames.filter((s,r)=>r%2===0),this.recEvery*=2),this.tracked.length&&this.uploadTracks(),this.onUpdate()}async computeAhead(e,t,n){const s=this.eng;if(!s||this.computing)return;const r=this.R.device,a=32*s.sys.n;for(this.playing=!1,this.replayIdx>=0&&this.resumeLive(),this.running=!1,this.computing=!0,this.abort=!1;this.reading;)await new Promise(x=>setTimeout(x,5));const i=Math.max(1,Math.round(e/this.set.dt));this.recEvery=Math.max(1,Math.round(i/t)),this.frames=[{tNs:s.timeNs,...await s.readPoses()}],this.uploadTracks(),this.onUpdate();const h=[],l=[],c=()=>{const x=h.pop()??r.createBuffer({size:a,usage:GPUBufferUsage.MAP_READ|GPUBufferUsage.COPY_DST});return l.includes(x)||l.push(x),x},d=[],f=async x=>{for(;d.length>x;){const y=d.shift();if(await y.p,!y.buf)continue;const A=new Float32Array(y.buf.getMappedRange().slice(0));y.buf.unmap(),h.push(y.buf),s===this.eng&&this.record({tNs:y.t,...s.unpackPoses(A)})}};let w=0,b=0,m=0,u=16;const p=performance.now();try{for(;w<i&&!this.abort&&s===this.eng;){if(m>0&&m%8===0){if(await f(0),this.abort||s!==this.eng)break;const M=await s.gridLoad();this.grid=M,M.dropped>0&&s.growGrid();const _=w/Math.max(.001,(performance.now()-p)/1e3);if(u=Math.max(4,Math.min(256,Math.round(_*.1))),this.abort||s!==this.eng)break}m++;const x=r.createCommandEncoder(),y=[];let A=0;for(;A<u&&w<i;){const M=Math.min(this.recEvery-b,i-w,u-A);if(s.step(x,M),w+=M,b+=M,A+=M,this.totalSteps+=M,b>=this.recEvery||w>=i){b=0;const _=c();s.copyPoses(x,_),y.push({buf:_,t:s.timeNs})}}r.queue.submit([x.finish()]);for(const M of y)d.push({...M,p:M.buf.mapAsync(GPUMapMode.READ)});if(y.length||d.push({buf:null,t:0,p:r.queue.onSubmittedWorkDone()}),await f(8),this.abort||s!==this.eng)break;this.aheadRate=w*this.set.dt/1e3/Math.max(.001,(performance.now()-p)/1e3),n(w/i)}await f(0).catch(()=>{})}finally{this.computing=!1;for(const x of l)x.destroy()}}play(e){e&&this.frames.length>1?(this.running=!1,this.playing=!0,this.playFrom=this.replayIdx>=0&&this.replayIdx<this.frames.length-1?this.replayIdx:0,this.playStart=performance.now()):this.playing=!1}advancePlayback(){const e=this.frames.length,t=(this.playFrom+Math.floor((performance.now()-this.playStart)*this.playFps/1e3))%e;if(t!==this.replayIdx){this.replayIdx=t;const n=this.frames[t];this.eng.writePoses(n.pos,n.quat),this.onUpdate()}}draw(e){if(!this.eng||!this.active)return;this.eng.draw(e);const t=Math.min(this.frames.length,Ze),n=Math.min(t-1,this.replayIdx>=0?this.replayIdx:t-1);this.tracked.length&&n>=1&&this.trackPipe&&this.trackBG&&this.trackG0&&(this.R.device.queue.writeBuffer(this.trackUni,0,new Uint32Array([n,Ze,Hn,this.tracked.length])),e.setPipeline(this.trackPipe),e.setBindGroup(0,this.trackG0),e.setBindGroup(1,this.trackBG),e.draw(6,this.tracked.length*n*Hn))}drawOverlay(e){this.box&&this.boxPipe&&this.boxBG&&(e.setPipeline(this.boxPipe),e.setBindGroup(0,this.boxG0),e.setBindGroup(1,this.boxBG),e.draw(6,12*480))}replay(e){if(!this.eng||!this.frames.length)return;if(this.playing=!1,e<0||e>=this.frames.length-1&&this.running){this.replayIdx=-1;return}this.replayIdx=Math.min(e,this.frames.length-1);const t=this.frames[this.replayIdx];this.eng.writePoses(t.pos,t.quat)}resumeLive(){if(!this.eng)return;this.playing=!1;const e=this.frames[this.frames.length-1];this.replayIdx>=0&&e&&this.eng.writePoses(e.pos,e.quat),this.replayIdx=-1}diffusion(){if(!this.eng||!this.meta)return[];const e=this.replayIdx>=0&&!this.playing?this.replayIdx:this.frames.length-1;if(e<1)return[];const t=this.frames,n=this.eng.sys,s=this.meta.hydro,r=Math.max(5,this.meta.box.size-(n.boundary==="embedded"?2*this.set.shell:0)),a=(r/4)**2,i=c=>c.upos??c.pos,h=[{label:"R_h < 3 nm (GFP-sized)",lo:0,hi:3},{label:"R_h 3-6 nm",lo:3,hi:6},{label:"R_h > 6 nm (ribosomes, large complexes)",lo:6,hi:1e9}],l=[];for(let c=1;c<=e;c=Math.max(c+1,Math.round(c*1.5)))l.push(c);return h.map(c=>{const d=[];let f=0;for(let b=0;b<n.n;b++)n.mobile[b]&&!n.tether[b]&&s[b].rh>=c.lo&&s[b].rh<c.hi&&(d.push(b),f+=s[b].dt);if(!d.length)return null;let w={d:NaN,lag:0};for(const b of l){const m=e-b+1,u=Math.max(1,Math.floor(m/40));let p=0,x=0,y=0;for(let M=0;M+b<=e;M+=u){const _=t[M],S=t[M+b],g=i(_),k=i(S);y+=S.tNs-_.tNs;for(const P of d){const T=k[3*P]-g[3*P],E=k[3*P+1]-g[3*P+1],I=k[3*P+2]-g[3*P+2];p+=T*T+E*E+I*I}x++}const A=p/(x*d.length);if(y/=x,!(y<=0)){if(w.lag&&A>a)break;w={d:A/(6*y)*1e3,lag:y}}}return{label:c.label,d:w.d,n:d.length,d0:f/d.length*1e3/this.set.viscosity,lagNs:w.lag}}).filter(c=>!!c)}async export(){if(!this.eng||!this.meta)return;const e=this.replayIdx>=0?this.frames[this.replayIdx]:await this.eng.readPoses();await $s(this.eng,this.meta,this.lib,this.set,e.pos,e.quat,this.frames)}}class js{target=[0,0,0];dist=9e3;fwd=[0,0,-1];up=[0,1,0];fovy=35*Math.PI/180;minDist=1.5;maxDist=14e3;anim=null;get right(){return he(Le(this.fwd,this.up))}get eye(){return Fe(this.target,O(this.fwd,this.dist))}orthonormalize(){this.fwd=he(this.fwd);const e=he(Le(this.fwd,this.up));this.up=he(Le(e,this.fwd))}orbit(e,t){const n=this.right;this.fwd=ut(this.fwd,this.up,-e),this.fwd=ut(this.fwd,n,-t),this.up=ut(this.up,n,-t),this.orthonormalize()}roll(e){this.up=ut(this.up,this.fwd,e),this.orthonormalize()}pan(e,t,n){const s=1/n;this.target=J(this.target,J(O(this.right,-e*s),O(this.up,t*s)))}zoom(e){this.dist=Math.min(this.maxDist,Math.max(this.minDist,this.dist*e))}fly(e){this.target=J(this.target,e)}pxPerNmAtUnit(e){return e/(2*Math.tan(this.fovy/2))}distForWidth(e,t){return e/(2*Math.tan(this.fovy/2)*Math.min(t,1))}near(){return Math.min(5,Math.max(.02,this.dist*.004))}viewRot(){const e=this.right,t=this.up,n=this.fwd,s=new Float32Array(16);return s[0]=e[0],s[4]=e[1],s[8]=e[2],s[1]=t[0],s[5]=t[1],s[9]=t[2],s[2]=-n[0],s[6]=-n[1],s[10]=-n[2],s[15]=1,s}proj(e){const t=1/Math.tan(this.fovy/2),n=new Float32Array(16);return n[0]=t/e,n[5]=t,n[11]=-1,n[14]=this.near(),n}frustum(e){const t=Math.tan(this.fovy/2),n=t*e,s=this.right,r=this.up,a=this.fwd,i=(l,c=0)=>{const d=he(l);return[d[0],d[1],d[2],c]},h=[...i(J(s,O(a,n))),...i(J(O(s,-1),O(a,n))),...i(J(r,O(a,t))),...i(J(O(r,-1),O(a,t))),a[0],a[1],a[2],-this.near(),0,0,0,1];return new Float32Array(h)}ray(e,t,n){const s=Math.tan(this.fovy/2),r=s*n,a=(e*2-1)*r,i=(1-t*2)*s;return he(J(this.fwd,J(O(this.right,a),O(this.up,i))))}state(){return{target:[...this.target],dist:this.dist,fwd:[...this.fwd],up:[...this.up]}}set(e){this.target=[...e.target],this.dist=e.dist,this.fwd=[...e.fwd],this.up=[...e.up],this.orthonormalize()}flyTo(e,t=1800,n,s){this.anim={t0:performance.now(),dur:t,a:this.state(),b:e,onDone:n,onStep:s}}step(e){if(!this.anim)return!1;const t=this.anim;let n=Math.min(1,(e-t.t0)/t.dur);const s=ls(n),r=Math.log(t.a.dist),a=Math.log(t.b.dist),i=Math.hypot(...Fe(t.b.target,t.a.target)),h=Math.log(Math.max(t.a.dist,t.b.dist,i*1.2)),l=Math.max(0,h-Math.max(r,a)),c=r+(a-r)*s+l*Math.sin(Math.PI*s);this.target=on(t.a.target,t.b.target,s),this.dist=Math.exp(c),this.fwd=Vs(t.a.fwd,t.b.fwd,s);let d=t.a.up;return this.up=he(on(d,t.b.up,s)),Math.abs(oe(this.up,this.fwd))>.98&&(this.up=he(Le(Le(this.fwd,t.b.up),this.fwd))),this.orthonormalize(),t.onStep?.(s),n>=1&&(this.anim=null,t.onDone?.()),!0}}function Vs(o,e,t){const n=Math.max(-1,Math.min(1,oe(o,e))),s=Math.acos(n);if(s<1e-4)return he(on(o,e,t));if(s>Math.PI-.001){const a=he(Le(o,Math.abs(o[0])<.9?[1,0,0]:[0,1,0]));return ut(o,a,Math.PI*t)}const r=Math.sin(s);return he(J(O(o,Math.sin((1-t)*s)/r),O(e,Math.sin(t*s)/r)))}const ln=0,jn=1,Vn=2;function It(o){return o&&o.length?o[0]:0}function cn(o,e,t,n,s){const r=1+12*s,a=o[r+3],i=o[r+7],h=(e-o[r])/a,l=(t-o[r+1])/a,c=(n-o[r+2])/a,d=o[r+4],f=o[r+5],w=o[r+6];if(h<0||l<0||c<0||h>=d*i||l>=f*i||c>=w*i)return o[r+9];const b=Math.floor(h),m=Math.floor(l),u=Math.floor(c),p=Math.floor(b/i),x=Math.floor(m/i),y=Math.floor(u/i),A=o[o[r+10]+(p*f+x)*w+y];if(A===-1)return-o[r+8];if(A===-2)return o[r+9];const M=i+1,_=M*M,S=M,g=h-b,k=l-m,P=c-u,T=o[r+11]+A*M*M*M+((b-p*i)*M+(m-x*i))*M+(u-y*i),E=o[T]+(o[T+_]-o[T])*g,I=o[T+1]+(o[T+_+1]-o[T+1])*g,F=o[T+S]+(o[T+_+S]-o[T+S])*g,N=o[T+S+1]+(o[T+_+S+1]-o[T+S+1])*g,G=E+(F-E)*k,$=I+(N-I)*k;return G+($-G)*P}function Rn(o,e,t,n,s){return It(o)>ln?-cn(o,t,n,s,ln):e[0]-Math.hypot(t,n,s)}function Ys(o,e,t,n,s,r=2){if(It(o)<=ln){const d=Math.hypot(t,n,s)||1;return[t/d,n/d,s/d]}const a=(d,f,w)=>Rn(o,e,d,f,w),i=a(t-r,n,s)-a(t+r,n,s),h=a(t,n-r,s)-a(t,n+r,s),l=a(t,n,s-r)-a(t,n,s+r),c=Math.hypot(i,h,l);return c>1e-9?[i/c,h/c,l/c]:[0,0,0]}async function Ks(o,e){const t=o.fields;if(!t)return null;const n=await e(t.file),s=t.head.length;let r=1+12*s;for(let h=0;h<s;h++)r+=t.n_index[h]+t.n_data[h];const a=new Float32Array(r);a[0]=s;let i=0;for(let h=0;h<s;h++){a.set(t.head[h],1+12*h);const l=new Int32Array(n,i,t.n_index[h]);i+=4*t.n_index[h],a.set(l,t.head[h][10]);const c=new Int16Array(n,i,t.n_data[h]);i+=2*t.n_data[h],i=Math.ceil(i/4)*4;const d=t.head[h][11];for(let f=0;f<c.length;f++)a[d+f]=c[f]*t.q}return a}const Co='function nn(n){return n&&n.length?n[0]:0}function fn(n,t,e,o,l){const s=1+12*l,y=n[s+3],u=n[s+7],f=(t-n[s])/y,i=(e-n[s+1])/y,b=(o-n[s+2])/y,g=n[s+4],j=n[s+5],c=n[s+6];if(f<0||i<0||b<0||f>=g*u||i>=j*u||b>=c*u)return n[s+9];const m=Math.floor(f),r=Math.floor(i),d=Math.floor(b),p=Math.floor(m/u),z=Math.floor(r/u),q=Math.floor(d/u),v=n[n[s+10]+(p*j+z)*c+q];if(v===-1)return-n[s+8];if(v===-2)return n[s+9];const A=u+1,x=A*A,G=A,B=f-m,w=i-r,h=b-d,a=n[s+11]+v*A*A*A+((m-p*u)*A+(r-z*u))*A+(d-q*u),M=n[a]+(n[a+x]-n[a])*B,U=n[a+1]+(n[a+x+1]-n[a+1])*B,T=n[a+G]+(n[a+x+G]-n[a+G])*B,L=n[a+G+1]+(n[a+x+G+1]-n[a+G+1])*B,C=M+(T-M)*w,X=U+(L-U)*w;return C+(X-C)*h}function dt(n,t,e,o,l){return nn(n)>0?-fn(n,e,o,l,0):t[0]-Math.hypot(e,o,l)}function Ht(n,t,e,o,l,s=2){if(nn(n)<=0){const g=Math.hypot(e,o,l)||1;return[e/g,o/g,l/g]}const y=(g,j,c)=>dt(n,t,g,j,c),u=y(e-s,o,l)-y(e+s,o,l),f=y(e,o-s,l)-y(e,o+s,l),i=y(e,o,l-s)-y(e,o,l+s),b=Math.hypot(u,f,i);return b>1e-9?[u/b,f/b,i/b]:[0,0,0]}const Ft=[[[16.8,20.3,17.7,33.6,.2,3.9,6.9],.4/.6],[[51.3,33.4,7.5,6.6,.4,.2,.5],.26],[[44.6,26.9,15.1,5.9,.5,2.2,4.3],.64],[[46.5,19.4,18.3,4.4,1.6,2.1,7.7],.18],[[45.6,32.6,10.2,1.2,5.9,4.4,0],.01],[[38.4,24,16.2,3.8,16.1,1.5,0],.32],[[35,22.3,19.1,12.9,.7,1.2,8.8],.75],[[48.2,22.9,15.8,4.5,7,1.6,0],.39],[[51.3,33.4,7.5,6.6,.4,.2,.5],0]],rn={4:1,1:.8,6:0};function xn(n,t){const[e,o]=Ft[Math.min(n,Ft.length-1)],[l,s,y,u,f,i,b]=e,g=e.reduce((r,d)=>r+d),j=1/(1+o),c=o/(1+o);let m=[y,s,l,i,u,0,b,f].map(r=>r/g*j*100);if(m[5]=c*100,n===0){m=m.map((d,p)=>p===5?d:d*2*(t?rn[p]??.5:1-(rn[p]??.5)));const r=m.reduce((d,p,z)=>d+(z===5?0:p),0);m=m.map((d,p)=>p===5?d:d/r*(100-m[5]))}return m}const An=Ft.map((n,t)=>[0,1].map(e=>{const o=xn(t,e),l=o.reduce((y,u)=>y+u);let s=0;return o.map(y=>s+=y/l)})),vn=300,an=4,N=200,Tt=80,Ct=360/65535,cn=.005,zn=.01,kn=.70710678;let Rt="",Dt={},It,_t,Nt,rt,Yt=null,tt=[],et=null,Et={glucanPerNm3:.018147,mannanPerProtein:1908.6},tn=[],jt="";const ct=new Map,Un=48;async function St(n){if(Dt[n])return ln(Dt[n].slice().buffer);const t=await fetch(Rt+n+(jt?`?v=${encodeURIComponent(jt)}`:""));if(!t.ok)throw new Error(`HTTP ${t.status} ${n}`);return ln(await t.arrayBuffer())}async function ln(n){const t=new Uint8Array(n,0,Math.min(2,n.byteLength));return t[0]===31&&t[1]===139?await new Response(new Blob([n]).stream().pipeThrough(new DecompressionStream("gzip"))).arrayBuffer():n}async function _n(n){Rt=n.base,tt=n.G,et=n.F??null,tn=n.scars??[];try{const f=await fetch(Rt+"version.json",{cache:"no-store"});f.ok&&(jt=String((await f.json()).v??""))}catch{jt=""}n.wall&&(Et=n.wall),n.embedded&&(Dt=n.embedded);const t=JSON.parse(new TextDecoder().decode(await St("library/types.json"))),e=t.types,o=e.length,l=t.bead_q;It=new Uint32Array(o),_t=new Uint32Array(o),Nt=new Float32Array(o);for(let f=0;f<o;f++)It[f]=e[f].bo,_t[f]=e[f].nb,Nt[f]=e[f].br;const s=await St("library/beads.bin"),y=new Int16Array(s),u=s.byteLength/8;rt=new Float32Array(u*3);for(let f=0;f<u;f++)rt[3*f]=y[4*f]*l,rt[3*f+1]=y[4*f+1]*l,rt[3*f+2]=y[4*f+2]*l}const Ut=n=>n+3&-4;function Kt(n,t,e,o,l,s,y,u,f,i,b,g=0){const j=new Uint32Array(n,t,e);t+=4*e;const c=new Int32Array(n,t,3*e);t+=12*e;const m=new Int8Array(n,t,3*(o-e));t+=3*(o-e),t=Ut(t);let r=0,d=i;for(let p=0;p<e;p++){let z=c[3*p],q=c[3*p+1],v=c[3*p+2];const A=j[p];for(let x=0;x<A;x++)x>0&&(z+=m[3*r],q+=m[3*r+1],v+=m[3*r+2],r++),u[4*d]=l[0]+z*s,u[4*d+1]=l[1]+q*s,u[4*d+2]=l[2]+v*s,f[4*d+3]=((b?b[p]:y)|(g+p&32767)<<4)>>>0,d++}return{off:t,o:d}}function hn(n){const t=new Uint32Array(n,0,16);if(t[0]!==827605849)throw new Error("bad tile magic");const e=new Int32Array(n,4,3),o=[e[0]*N,e[1]*N,e[2]*N],[l,s,y,u,f,i,b,g,j,c,m]=[t[4],t[5],t[6],t[7],t[8],t[9],t[10],t[11],t[12],t[13],t[14]];let r=64;const d=new Uint32Array(n.slice(r,r+16*l));r+=16*l;const p=new Uint32Array(n.slice(r,r+16*s));r+=16*s;const z=new Float32Array((u+i+m)*4),q=new Uint32Array(z.buffer);let v=0;if(y){const h=Kt(n,r,y,u,o,cn,0,z,q,v);r=h.off,v=h.o}if(f){const h=Kt(n,r,f,i,o,zn,1,z,q,v);r=h.off,v=h.o}const A=new Float32Array(b*3);if(b){const h=new Uint16Array(n,r,3*b);for(let a=0;a<b;a++)for(let M=0;M<3;M++)A[3*a+M]=o[M]-Tt+h[3*a+M]*Ct;r+=Ut(6*b)}let x=new Uint32Array(g*3);g&&(j&1?x=new Uint32Array(n.slice(r,r+12*g)):x=Uint32Array.from(new Uint16Array(n,r,3*g)),r+=Ut((j&1?4:2)*3*g));const G=g?new Uint8Array(n.slice(r,r+g)):new Uint8Array(0);if(g&&(r+=Ut(g)),c){const h=new Uint8Array(n.slice(r,r+c));r+=Ut(c),Kt(n,r,c,m,o,cn,0,z,q,v,h,y+f)}const B=new Float32Array(b*3),w=new Uint8Array(b);for(let h=0;h<g;h++){const a=x[3*h],M=x[3*h+1],U=x[3*h+2],T=A[3*M]-A[3*a],L=A[3*M+1]-A[3*a+1],C=A[3*M+2]-A[3*a+2],X=A[3*U]-A[3*a],Q=A[3*U+1]-A[3*a+1],H=A[3*U+2]-A[3*a+2],F=L*H-C*Q,E=C*X-T*H,k=T*Q-L*X;for(const I of[a,M,U])B[3*I]+=F,B[3*I+1]+=E,B[3*I+2]+=k,w[I]=G[h]}for(let h=0;h<b;h++){const a=Math.hypot(B[3*h],B[3*h+1],B[3*h+2])||1;B[3*h]/=a,B[3*h+1]/=a,B[3*h+2]/=a}return{lo:o,inst:d,ghost:p,free:z,verts:A,normals:B,tris:x,fmem:G,vmem:w,nI:l,nG:s}}const Zt=n=>(n<<16>>16)/32767*kn;function un(n,t,e,o){const l=n[4*t],s=n[4*t+1],y=n[4*t+2],u=n[4*t+3];o[0]=e[0]-Tt+(s&65535)*Ct,o[1]=e[1]-Tt+(s>>>16)*Ct,o[2]=e[2]-Tt+(y&65535)*Ct;const f=Zt(y>>>16),i=Zt(u&65535),b=Zt(u>>>16),g=Math.sqrt(Math.max(0,1-f*f-i*i-b*b)),j=l>>>14&3,c=j===0?[g,f,i,b]:j===1?[f,g,i,b]:j===2?[f,i,g,b]:[f,i,b,g],[m,r,d,p]=c;return o[3]=m*m+r*r-d*d-p*p,o[4]=2*(r*d-m*p),o[5]=2*(r*p+m*d),o[6]=2*(r*d+m*p),o[7]=m*m-r*r+d*d-p*p,o[8]=2*(d*p-m*r),o[9]=2*(r*p-m*d),o[10]=2*(d*p+m*r),o[11]=m*m-r*r-d*d+p*p,{type:l&16383,region:l>>>16&31,flags:l>>>21}}class Jt{h;tbl;mask;start;pts;inst;constructor(t,e,o,l){this.h=l;let s=1;for(;s<o*2;)s<<=1;s=Math.max(s,1024),this.tbl=s,this.mask=s-1;const y=new Uint32Array(o),u=new Uint32Array(s+1);for(let i=0;i<o;i++)y[i]=this.k(Math.floor(t[3*i]/l),Math.floor(t[3*i+1]/l),Math.floor(t[3*i+2]/l)),u[y[i]+1]++;for(let i=0;i<s;i++)u[i+1]+=u[i];this.start=u.slice(),this.pts=new Float32Array(o*3),this.inst=new Int32Array(o);const f=u.slice(0,s);for(let i=0;i<o;i++){const b=f[y[i]]++;this.pts[3*b]=t[3*i],this.pts[3*b+1]=t[3*i+1],this.pts[3*b+2]=t[3*i+2],this.inst[b]=e[i]}}k(t,e,o){return(Math.imul(t,73856093)^Math.imul(e,19349663)^Math.imul(o,83492791))&this.mask}near(t,e,o,l,s=-1){const y=this.h,u=Math.floor(t/y),f=Math.floor(e/y),i=Math.floor(o/y),b=l*l;for(let g=-1;g<=1;g++)for(let j=-1;j<=1;j++)for(let c=-1;c<=1;c++){const m=this.k(u+g,f+j,i+c);for(let r=this.start[m],d=this.start[m+1];r<d;r++){const p=this.pts[3*r]-t,z=this.pts[3*r+1]-e,q=this.pts[3*r+2]-o;if(p*p+z*z+q*q<b&&this.inst[r]!==s)return!0}}return!1}}function Wt(n){let t=0;const e=new Float64Array(12);for(const y of n)t+=_t[y.w[4*y.i]&16383];const o=new Float32Array(t*3),l=new Int32Array(t);let s=0;for(const y of n){const{type:u}=un(y.w,y.i,y.lo,e),f=It[u],i=_t[u];for(let b=0;b<i;b++){const g=rt[3*(f+b)],j=rt[3*(f+b)+1],c=rt[3*(f+b)+2];o[3*s]=e[0]+e[3]*g+e[4]*j+e[5]*c,o[3*s+1]=e[1]+e[6]*g+e[7]*j+e[8]*c,o[3*s+2]=e[2]+e[9]*g+e[10]*j+e[11]*c,l[s]=y.id,s++}}return{xyz:o,owner:l,n:t}}function zt(n){let t=n>>>0;return()=>{t=t+1831565813>>>0;let e=t;return e=Math.imul(e^e>>>15,e|1),e^=e+Math.imul(e^e>>>7,e|61),((e^e>>>14)>>>0)/4294967296}}const kt=(n,t,e)=>{let o=Math.imul(n^2654435769,2246822507)^Math.imul(t+1663821227,3266489909)^Math.imul(e+668265263,374761393);return o^=o>>>16,o=Math.imul(o,2146121005),o^=o>>>15,o=Math.imul(o,2221713035),o^=o>>>16,o>>>0};function W(n){return Math.sqrt(-2*Math.log(n()+1e-12))*Math.cos(6.2831853*n())}function Ln(n,t,e){const o=Math.abs(n)+Math.abs(t)+Math.abs(e)||1;let l=n/o,s=t/o;if(e<0){const f=l;l=(1-Math.abs(s))*(f>=0?1:-1),s=(1-Math.abs(f))*(s>=0?1:-1)}const y=Math.max(0,Math.min(255,Math.round((l+1)*127.5))),u=Math.max(0,Math.min(255,Math.round((s+1)*127.5)));return y|u<<8}const Bn=.65,qn=.95;function yt(n,t,e,o){return n.n[0]*t+n.n[1]*e+n.n[2]*o-n.d}function Tn(n,t,e,o){const l=[],s=[],y={lipCand:0,lipMasked:0,glucan:0,mannan:0,chitin:0,polyp:0},u=new Float64Array(12),f=(c,m,r,d)=>{if(e.mode==="shell")return dt(et,tt,c,m,r)<70+d;const p=yt(e,c,m,r);return p>=e.lo-d&&p<=e.hi+d};if(o.lipids&&t.tris.length&&e.mode==="slab"){const c=[],m=t.inst.length/4,r=t.ghost.length/4;for(let x=0;x<m;x++)t.inst[4*x]>>>21&1&&c.push({w:t.inst,i:x,lo:t.lo,id:x});for(let x=0;x<r;x++)c.push({w:t.ghost,i:x,lo:t.lo,id:m+x});const d=Wt(c),p=new Jt(d.xyz,d.owner,d.n,1),z=t.verts,q=t.normals,v=t.tris,A=v.length/3;for(let x=0;x<A;x++){const G=v[3*x],B=v[3*x+1],w=v[3*x+2],h=z[3*G],a=z[3*G+1],M=z[3*G+2],U=z[3*B],T=z[3*B+1],L=z[3*B+2],C=z[3*w],X=z[3*w+1],Q=z[3*w+2],H=yt(e,h,a,M),F=yt(e,U,T,L),E=yt(e,C,X,Q);if(Math.max(H,F,E)<e.lo-3||Math.min(H,F,E)>e.hi+3)continue;const k=U-h,I=T-a,Y=L-M,Z=C-h,R=X-a,_=Q-M,P=.5*Math.hypot(I*_-Y*R,Y*Z-k*_,k*R-I*Z),$=t.fmem[x],K=$===8?[1]:[-1,1];for(const V of K){const O=zt(kt(n,x,V+7)),S=P/Bn,it=Math.floor(S)+(O()<S-Math.floor(S)?1:0);if(!it)continue;const nt=Math.ceil(Math.sqrt(it)),ft=nt*nt,ot=new Uint16Array(ft);for(let J=0;J<ft;J++)ot[J]=J;for(let J=0;J<it;J++){const at=J+Math.floor(O()*(ft-J)),D=ot[J];ot[J]=ot[at],ot[at]=D}for(let J=0;J<it;J++){const at=ot[J],D=Math.floor(Math.sqrt(at)),ut=at-D*D,st=ut>>1,Mt=(ut&1)===0;let mt,pt,Pt,$t,Vt,Gt;Mt?(mt=D,pt=st,Pt=D+1,$t=st,Vt=D+1,Gt=st+1):(mt=D,pt=st,Pt=D+1,$t=st+1,Vt=D,Gt=st+1);let wt=O(),bt=O();wt+bt>1&&(wt=1-wt,bt=1-bt);const yn=mt+wt*(Pt-mt)+bt*(Vt-mt),en=pt+wt*($t-pt)+bt*(Gt-pt),lt=(yn-en)/nt,ht=en/nt,Ot=1-lt-ht,dn=h+lt*k+ht*Z,Mn=a+lt*I+ht*R,mn=M+lt*Y+ht*_;let gt=Ot*q[3*G]+lt*q[3*B]+ht*q[3*w],xt=Ot*q[3*G+1]+lt*q[3*B+1]+ht*q[3*w+1],At=Ot*q[3*G+2]+lt*q[3*B+2]+ht*q[3*w+2];const Qt=Math.hypot(gt,xt,At)||1;gt/=Qt,xt/=Qt,At/=Qt;const pn=O(),wn=An[Math.min($,8)][V>0?1:0];let vt=0;for(;vt<7&&pn>wn[vt];)vt++;y.lipCand++;let Lt=dn,Bt=Mn,qt=mn;$===8&&(Lt-=gt*2.4,Bt-=xt*2.4,qt-=At*2.4);const on=yt(e,Lt,Bt,qt);if(on<e.lo-2.5||on>e.hi+2.5)continue;const bn=vt===5?[1.55,.75]:[1.95,1.25,.55];let sn=!1;for(const Xt of bn)if(p.near(Lt+gt*V*Xt,Bt+xt*V*Xt,qt+At*V*Xt,qn)){sn=!0;break}if(sn){y.lipMasked++;continue}const gn=(Ln(gt,xt,At)|vt<<16|(V>0?1:0)<<19|($&15)<<20)>>>0;l.push(Lt,Bt,qt,gn)}}}}if(o.wall&&tt.length>0){const c=tt[0]-tt[2],m=tt[0]-tt[1],r=t.lo,d=t.inst.length/4,p=[];for(let w=0;w<d;w++)(t.inst[4*w]>>>16&31)===9&&p.push({w:t.inst,i:w,lo:r,id:w});const z=Wt(p),q=new Jt(z.xyz,z.owner,z.n,1.4),v=(w,h,a)=>w>=r[0]&&w<r[0]+N&&h>=r[1]&&h<r[1]+N&&a>=r[2]&&a<r[2]+N,A=self.__vol?.[n]??0,x=Math.round(A*Et.glucanPerNm3);if(x>0){const w=zt(kt(n,1234567,1)),a=Math.ceil(x/48);for(let M=0;M<a;M++){const U=24+Math.floor(w()*49);let T=0,L=0,C=0,X=!1;for(let R=0;R<40&&!X;R++){T=r[0]+w()*N,L=r[1]+w()*N,C=r[2]+w()*N;const _=dt(et,tt,T,L,C);X=_<c-2&&_>m+1}if(!X)continue;let[Q,H,F]=Ht(et,tt,T,L,C),E=W(w),k=W(w),I=W(w),Y=E*Q+k*H+I*F;E-=Y*Q,k-=Y*H,I-=Y*F;let Z=Math.hypot(E,k,I)||1;E/=Z,k/=Z,I/=Z;for(let R=0;R<U&&v(T,L,C);R++){f(T,L,C,1)&&!q.near(T,L,C,1.35)&&(s.push(T,L,C,(2|(M&65535)<<4)>>>0),y.glucan++),E+=.12*W(w),k+=.12*W(w),I+=.12*W(w),[Q,H,F]=Ht(et,tt,T,L,C),Y=E*Q+k*H+I*F;const _=dt(et,tt,T,L,C),P=_>c-3&&Y<0||_<m+2&&Y>0?-1:.6;E-=(1-P)*Y*Q,k-=(1-P)*Y*H,I-=(1-P)*Y*F,Z=Math.hypot(E,k,I)||1,E/=Z,k/=Z,I/=Z,T+=1.2*E,L+=1.2*k,C+=1.2*I}}}const G=Et.chitinPerNm3??0;for(let w=0;G>0&&w+8<=tn.length;w+=8){const[h,a,M,U,T,L,C,X]=tn.slice(w,w+8),Q=Math.min(Math.max(h,r[0]),r[0]+N),H=Math.min(Math.max(a,r[1]),r[1]+N),F=Math.min(Math.max(M,r[2]),r[2]+N);if(Math.hypot(Q-h,H-a,F-M)>U+X)continue;const E=(_,P,$)=>{const K=_-h,V=P-a,O=$-M,S=K*T+V*L+O*C;return Math.hypot(Math.sqrt(Math.max(K*K+V*V+O*O-S*S,0))-U,S)<X-1.2},k=zt(kt(n,7654321,w));let I=0;for(let _=0;_<4e3;_++)E(r[0]+k()*N,r[1]+k()*N,r[2]+k()*N)&&I++;const Y=Math.round(I/4e3*N**3*G);let Z=0,R=0;for(;Z<Y&&R<1e5;){let _=0,P=0,$=0,K=!1;for(let S=0;S<60&&!K;S++)_=r[0]+k()*N,P=r[1]+k()*N,$=r[2]+k()*N,K=E(_,P,$);if(!K)break;const V=30+Math.floor(k()*50),O=R++;for(let S=0;S<V&&Z<Y;S++){const it=_-h,nt=P-a,ft=$-M;let ot=L*ft-C*nt,J=C*it-T*ft,at=T*nt-L*it;const D=Math.hypot(ot,J,at)||1;ot/=D,J/=D,at/=D;const ut=_+1.04*(ot+.08*W(k)),st=P+1.04*(J+.08*W(k)),Mt=$+1.04*(at+.08*W(k));if(!E(ut,st,Mt)||!v(ut,st,Mt))break;_=ut,P=st,$=Mt,Z++,f(_,P,$,1)&&!q.near(_,P,$,1.6)&&(s.push(_,P,$,(6|(O&32767)<<4)>>>0),y.chitin++)}}}let B=0;for(const w of p){const{type:h}=un(w.w,w.i,r,u),a=u[0],M=u[1],U=u[2],T=dt(et,tt,a,M,U);if(T>m+5)continue;if(e.mode==="slab"){const k=yt(e,a,M,U);if(k<e.lo-45-Nt[h]||k>e.hi+45+Nt[h]){B+=80;continue}}else if(T>120)continue;const L=zt(kt(n,w.i,99)),C=It[h],X=_t[h];let Q=Et.mannanPerProtein;const[H,F,E]=Ht(et,tt,a,M,U);for(;Q>0;){const k=Math.min(Math.ceil(Q),12+Math.floor(L()*27));Q-=k;const I=C+Math.floor(L()*X),Y=rt[3*I],Z=rt[3*I+1],R=rt[3*I+2];let _=a+u[3]*Y+u[4]*Z+u[5]*R,P=M+u[6]*Y+u[7]*Z+u[8]*R,$=U+u[9]*Y+u[10]*Z+u[11]*R,K=_-a,V=P-M,O=$-U,S=Math.hypot(K,V,O)||1;K=K/S+H,V=V/S+F,O=O/S+E,S=Math.hypot(K,V,O)||1,K/=S,V/=S,O/=S,_+=1.3*K,P+=1.3*V,$+=1.3*O;const it=B++;for(let nt=0;nt<k&&!(dt(et,tt,_,P,$)<.6);nt++)f(_,P,$,1)&&!q.near(_,P,$,1.15,nt<2?w.i:-1)&&(s.push(_,P,$,(3|(it&65535)<<4)>>>0),y.mannan++),K+=.35*W(L)+.18*H,V+=.35*W(L)+.18*F,O+=.35*W(L)+.18*E,S=Math.hypot(K,V,O)||1,K/=S,V/=S,O/=S,_+=1.05*K,P+=1.05*V,$+=1.05*O}}}if(o.vac&&e.mode==="slab"&&nn(et)>1){const c=t.lo,m=[],r=t.inst.length/4;for(let h=0;h<r;h++){const a=t.inst[4*h]>>>16&31;(a===4||a===17)&&m.push({w:t.inst,i:h,lo:c,id:h})}const d=Wt(m),p=new Jt(d.xyz,d.owner,d.n,1.2),z=(h,a,M)=>fn(et,h,a,M,1)<-3,q=(h,a,M)=>h>=c[0]&&h<c[0]+N&&a>=c[1]&&a<c[1]+N&&M>=c[2]&&M<c[2]+N,v=zt(kt(n,4242421,7)),A=[],x=8e3;for(let h=0;h<x;h++){const a=c[0]+v()*N,M=c[1]+v()*N,U=c[2]+v()*N;z(a,M,U)&&f(a,M,U,0)&&A.push(a,M,U)}const G=Math.round(A.length/3/x*N**3*vn*6022e-7/an);let B=0,w=0;for(;B<G&&A.length&&w<2e5;){const h=3*Math.floor(v()*A.length/3);let a=A[h]+(v()-.5)*6,M=A[h+1]+(v()-.5)*6,U=A[h+2]+(v()-.5)*6;const T=Math.max(4,Math.round(Math.exp(Math.log(24)-.18+.6*W(v)))),L=Math.max(1,Math.round(T/an)),C=w++;let X=W(v),Q=W(v),H=W(v);for(let F=0;F<L&&B<G&&!(!z(a,M,U)||!q(a,M,U));F++){B++,f(a,M,U,1)&&!p.near(a,M,U,1.1)&&(s.push(a,M,U,(7|(C&32767)<<4)>>>0),y.polyp++),X+=.6*W(v),Q+=.6*W(v),H+=.6*W(v);const E=Math.hypot(X,Q,H)||1;a+=1.2*X/E,M+=1.2*Q/E,U+=1.2*H/E}}}const i=new Float32Array(l.length),b=new Uint32Array(i.buffer);for(let c=0;c<l.length;c+=4)i[c]=l[c],i[c+1]=l[c+1],i[c+2]=l[c+2],b[c+3]=l[c+3];const g=new Float32Array(s.length),j=new Uint32Array(g.buffer);for(let c=0;c<s.length;c+=4)g[c]=s[c],g[c+1]=s[c+1],g[c+2]=s[c+2],j[c+3]=s[c+3];return{lipids:i,free:g,stats:y}}self.onmessage=async n=>{const t=n.data;try{if(t.cmd==="init"){Yt=_n(t),self.__vol=t.volGlucan,await Yt,self.postMessage({cmd:"ready"});return}if(await Yt,t.cmd==="load"){const e=await St(`tiles/${t.ix}_${t.iy}_${t.iz}.bin.gz`),o=hn(e),l=o.verts.length/3,s=new Float32Array(l*8),y=new Uint32Array(s.buffer);for(let i=0;i<l;i++)s[8*i]=o.verts[3*i],s[8*i+1]=o.verts[3*i+1],s[8*i+2]=o.verts[3*i+2],y[8*i+3]=o.vmem[i],s[8*i+4]=o.normals[3*i],s[8*i+5]=o.normals[3*i+1],s[8*i+6]=o.normals[3*i+2];for(ct.set(t.idx,{lo:o.lo,inst:o.inst,ghost:o.ghost,verts:o.verts,normals:o.normals,tris:o.tris,fmem:o.fmem,vmem:o.vmem});ct.size>Un;)ct.delete(ct.keys().next().value);const u=o.inst.slice(),f=o.tris.slice();self.postMessage({cmd:"loaded",idx:t.idx,inst:u,free:o.free,vb:s,ib:f,lo:o.lo},[u.buffer,o.free.buffer,s.buffer,f.buffer])}else if(t.cmd==="detail"){let e=ct.get(t.idx);if(e)ct.delete(t.idx),ct.set(t.idx,e);else{const s=hn(await St(`tiles/${t.ix}_${t.iy}_${t.iz}.bin.gz`));e={lo:s.lo,inst:s.inst,ghost:s.ghost,verts:s.verts,normals:s.normals,tris:s.tris,fmem:s.fmem,vmem:s.vmem},ct.set(t.idx,e)}const o=performance.now(),l=Tn(t.idx,e,t.range,t.want);self.postMessage({cmd:"detail",idx:t.idx,key:t.key,lipids:l.lipids,free:l.free,stats:l.stats,ms:performance.now()-o},[l.lipids.buffer,l.free.buffer])}}catch(e){self.postMessage({cmd:"error",idx:t.idx,key:t.key,req:t.cmd,error:String(e&&e.stack||e)})}};\n',Yn=typeof self<"u"&&self.Blob&&new Blob(["URL.revokeObjectURL(import.meta.url);",Co],{type:"text/javascript;charset=utf-8"});function Js(o){let e;try{if(e=Yn&&(self.URL||self.webkitURL).createObjectURL(Yn),!e)throw"";const t=new Worker(e,{type:"module",name:o?.name});return t.addEventListener("error",()=>{(self.URL||self.webkitURL).revokeObjectURL(e)}),t}catch{return new Worker("data:text/javascript;charset=utf-8,"+encodeURIComponent(Co),{type:"module",name:o?.name})}}const et=200,Kn=173.3,Jt=70,Xs=60;class Xt{used;n;count=0;high=0;constructor(e){this.n=e,this.used=new Uint8Array(e)}alloc(e){if(e<=0)return[];if(this.n-this.count<e)return null;const t=[];for(let n=0;n<this.n&&t.length<e;n++)this.used[n]||(this.used[n]=1,t.push(n));return this.count+=e,this.high=Math.max(this.high,t[t.length-1]+1),t}free(e){for(const t of e)this.used[t]&&(this.used[t]=0,this.count--);for(;this.high>0&&!this.used[this.high-1];)this.high--}get freeCount(){return this.n-this.count}}class Qs{idx;r;lib;G;F=null;workers=[];state;loaded=new Map;instA;lipA;freeA;slots=[];pageTile;pageLocal;tileState=new Uint32Array(pt*pt*pt);desired=[];desiredSet=new Set;inflight=0;maxInflight=12;detInflight=0;maxDetInflight=3;frameNo=0;onChange=()=>{};onTileLoaded=[];bytesLoaded=0;tilesLoadedTotal=0;starved=!1;detailStats={lipids:0,glucan:0,mannan:0};detailEnabled=!0;lastLoadErr="";shellRange=null;constructor(e,t,n,s,r,a){this.idx=e,this.r=t,this.lib=n,this.G=s.G,this.F=s.F??null,this.state=new Uint8Array(e.n),this.instA=new Xt(t.instPages),this.lipA=new Xt(t.lipPages),this.freeA=new Xt(t.freePages);for(let l=t.maxSlots-1;l>=0;l--)this.slots.push(l);this.pageTile=new Int32Array(t.instPages).fill(-1),this.pageLocal=new Uint32Array(t.instPages);const i=Array.from(e.volGlucan),h=new URL("data/",document.baseURI).href;for(let l=0;l<r;l++){const c=new Js;c.onmessage=d=>this.onMessage(d.data),c.postMessage({cmd:"init",base:h,G:s.G,F:s.F,volGlucan:i,embedded:a,scars:co(s),wall:{glucanPerNm3:e.wall.glucan_beads_per_nm3,mannanPerProtein:e.wall.mannan_beads_per_protein,chitinPerNm3:e.wall.chitin_beads_per_nm3_scar??0}}),this.workers.push(c)}}tileCenter(e){return[this.idx.ix[e]*et+100,this.idx.iy[e]*et+100,this.idx.iz[e]*et+100]}update(e,t,n,s,r,a,i){if(this.frameNo++,a){this.desired=[],this.desiredSet.clear();return}const h=e.eye,l=e.frustum(t),c=this.idx.n,d=[s[0],s[1],s[2]],f=100*(Math.abs(d[0])+Math.abs(d[1])+Math.abs(d[2])),w=[];for(let p=0;p<c;p++){const x=this.tileCenter(p),y=Fe(x,h);let A=!0;for(let M=0;M<5;M++)if(l[4*M]*y[0]+l[4*M+1]*y[1]+l[4*M+2]*y[2]+l[4*M+3]<-243.3){A=!1;break}if(A){if(n){const M=oe(d,x)-s[3];if(M+f+Jt<0||M-f-Jt>r)continue}else{const M=this.tileShellRange();if(M[5*p]>Xs+87||M[5*p+1]<-87-Jt)continue;const _=Fe(h,x),S=5*p+2;if(M[S]*_[0]+M[S+1]*_[1]+M[S+2]*_[2]<-300)continue}w.push({i:p,d:Je(y)})}}w.sort((p,x)=>p.d-x.d);let b=this.r.instPages*ce*.94,m=this.r.maxSlots-8;const u=[];for(const p of w){const x=this.idx.nInst[p.i];if(x>b||m<=0){this.starved=!0;break}b-=Math.ceil(x/ce)*ce,m--,u.push(p.i)}u.length===w.length&&(this.starved=!1),this.desired=u,this.desiredSet=new Set(u);for(const p of u){const x=this.loaded.get(p);x&&(x.lastUsed=this.frameNo)}for(const p of u){if(this.inflight>=this.maxInflight)break;this.state[p]===0&&(this.state[p]=1,this.inflight++,this.workers[p%this.workers.length].postMessage({cmd:"load",idx:p,ix:this.idx.ix[p],iy:this.idx.iy[p],iz:this.idx.iz[p]}))}this.updateDetail(e,n,s,r,i)}tileShellRange(){if(this.shellRange)return this.shellRange;const e=this.idx.n,t=new Float32Array(5*e);for(let n=0;n<e;n++){const s=this.idx.ix[n]*et,r=this.idx.iy[n]*et,a=this.idx.iz[n]*et;let i=1/0,h=-1/0;for(let c=0;c<=2;c++)for(let d=0;d<=2;d++)for(let f=0;f<=2;f++){const w=Rn(this.F,this.G,s+100*c,r+100*d,a+100*f);i=Math.min(i,w),h=Math.max(h,w)}const l=Ys(this.F,this.G,s+100,r+100,a+100,60);t[5*n]=i,t[5*n+1]=h,t.set(l,5*n+2)}return this.shellRange=t}updateDetail(e,t,n,s,r){if(!this.detailEnabled)return;const a=e.eye,i=e.pxPerNmAtUnit(r),h=[n[0],n[1],n[2]],l=Math.min(s,250),c=[];for(const f of this.desired){const w=this.loaded.get(f);if(!w)continue;const b=this.idx.nTri[f]>0,m=this.idx.volGlucan[f]>0||(this.idx.labelMask[f]&2)!==0,u=(this.idx.labelMask[f]&256)!==0;if(!b&&!m&&!(u&&t))continue;const p=this.tileCenter(f),x=Math.max(1,Je(Fe(p,a))-Kn);i/(t?Math.max(1,e.dist-Kn):x)*.55<.8||!t&&!m||c.push({t:w,d:x})}c.sort((f,w)=>f.d-w.d);const d=new Set(c.slice(0,40).map(f=>f.t.idx));for(const f of this.loaded.values())f.det&&!d.has(f.idx)&&!f.detPending&&this.dropDetail(f);for(const{t:f}of c.slice(0,40)){if(this.detInflight>=this.maxDetInflight)break;if(f.detPending)continue;const w=f.det;let b=!1;if(w&&(t?b=!w.shell&&oe(w.n,h)>.99985&&w.a<=n[3]+.5&&w.b>=n[3]+l-.5:b=w.shell),b)continue;const m=Math.max(10,.35*l),u=t?{n:h,d:n[3],lo:-m,hi:l+m,mode:"slab"}:{n:h,d:n[3],lo:0,hi:0,mode:"shell"};f.detPending=!0,this.detInflight++;const p=t?`${n[3]-m}|${n[3]+l+m}`:"shell";f._req={n:h,a:n[3]-m,b:n[3]+l+m,shell:!t};const x=this.idx.nTri[f.idx]>0,y=this.idx.volGlucan[f.idx]>0||(this.idx.labelMask[f.idx]&2)!==0,A=(this.idx.labelMask[f.idx]&256)!==0;this.workers[f.idx%this.workers.length].postMessage({cmd:"detail",idx:f.idx,key:p,ix:this.idx.ix[f.idx],iy:this.idx.iy[f.idx],iz:this.idx.iz[f.idx],range:u,want:{lipids:x&&t,wall:y,vac:A&&t}})}}dropDetail(e){this.lipA.free(e.lipPages),this.freeA.free(e.detPages),this.clearSPages(this.r.lipPagesBuf,e.lipPages),this.clearSPages(this.r.freePagesBuf,e.detPages),e.lipPages=[],e.detPages=[],e.nLip=0,e.nDet=0,e.det=null,e.underMesh&&(e.underMesh=!1,this.r.setMeshSlot(e.slot,!1)),this.syncHigh()}clearSPages(e,t){const n=new Uint32Array(4);for(const s of t)this.r.device.queue.writeBuffer(e,s*16,n)}syncHigh(){this.r.usedInstPages=this.instA.high,this.r.usedLipPages=this.lipA.high,this.r.usedFreePages=this.freeA.high}evictOne(e){let t=null;for(const n of this.loaded.values())n.idx===e||this.desiredSet.has(n.idx)||this.pinned.has(n.idx)||(!t||n.lastUsed<t.lastUsed)&&(t=n);return t?(this.evict(t),!0):!1}evict(e){const t=this.r.device.queue,n=new Float32Array(4);for(const s of e.pages)t.writeBuffer(this.r.pageOrg,s*16,n),this.pageTile[s]=-1;this.instA.free(e.pages),this.freeA.free(e.basePages),this.clearSPages(this.r.freePagesBuf,e.basePages),this.lipA.free(e.lipPages),this.clearSPages(this.r.lipPagesBuf,e.lipPages),this.freeA.free(e.detPages),this.clearSPages(this.r.freePagesBuf,e.detPages),e.vb?.destroy(),e.ib?.destroy(),this.slots.push(e.slot),this.loaded.delete(e.idx),this.state[e.idx]=0,this.tileState[_t(this.idx.ix[e.idx],this.idx.iy[e.idx],this.idx.iz[e.idx])]=0,t.writeBuffer(this.r.tileState,0,this.tileState),this.syncHigh()}allocWithEvict(e,t,n){let s=e.alloc(t);for(;!s&&this.evictOne(n);)s=e.alloc(t);return s}onMessage(e){if(e.cmd!=="ready"){if(e.cmd==="error"){if(console.warn("worker error",e.req,e.idx,e.error),this.lastLoadErr=String(e.error).split(`
`)[0],e.req==="load"&&(this.state[e.idx]=0,this.inflight--),e.req==="detail"){this.detInflight--;const t=this.loaded.get(e.idx);t&&(t.detPending=!1)}return}e.cmd==="loaded"&&(this.inflight--,this.onLoaded(e)),e.cmd==="detail"&&(this.detInflight--,this.onDetail(e))}}onLoaded(e){const t=e.idx;if(!this.desiredSet.has(t)&&!this.pinned.has(t)){this.state[t]=0,this.onChange();return}const n=e.inst,s=n.length/4,r=this.r.device.queue;this.slots.length||this.evictOne(t);const a=this.allocWithEvict(this.instA,Math.ceil(s/ce),t);if(!a||!this.slots.length){this.state[t]=0,this.starved=!0,a&&this.instA.free(a);return}const i=this.slots.pop(),h=e.lo,l=new ArrayBuffer(16),c=new Float32Array(l),d=new Uint32Array(l);c[0]=h[0]-80,c[1]=h[1]-80,c[2]=h[2]-80,a.forEach((A,M)=>{const _=Math.min(ce,s-M*ce);r.writeBuffer(this.r.inst,A*ce*16,n.buffer,M*ce*16,_*16),d[3]=_,r.writeBuffer(this.r.pageOrg,A*16,l),this.pageTile[A]=t,this.pageLocal[A]=M*ce});const f=e.free,w=f.length/4;let b=[],m=!1;if(w){const A=this.allocWithEvict(this.freeA,Math.ceil(w/de),t);A?(b=A,this.writeSpheres(this.r.freePool,this.r.freePagesBuf,A,f,i)):m=!0}let u=null,p=null;const x=e.ib.length;x&&(u=this.r.device.createBuffer({size:e.vb.byteLength,usage:GPUBufferUsage.VERTEX|GPUBufferUsage.COPY_DST}),r.writeBuffer(u,0,e.vb),p=this.r.device.createBuffer({size:Math.ceil(e.ib.byteLength/4)*4,usage:GPUBufferUsage.INDEX|GPUBufferUsage.COPY_DST}),r.writeBuffer(p,0,e.ib),this.r.setMeshSlot(i,!1));const y={idx:t,slot:i,lo:h,inst:n,nInst:s,pages:a,basePages:b,nBase:w,lipPages:[],detPages:[],nLip:0,nDet:0,vb:u,ib:p,nIdx:x,lastUsed:this.frameNo,det:null,detPending:!1,underMesh:!1,baseFreeMissing:m};this.loaded.set(t,y),this.state[t]=2,this.tileState[_t(this.idx.ix[t],this.idx.iy[t],this.idx.iz[t])]=1,r.writeBuffer(this.r.tileState,0,this.tileState),this.bytesLoaded+=this.idx.gz[t],this.tilesLoadedTotal++,this.syncHigh();for(const A of this.onTileLoaded.slice())A(y);this.onChange()}writeSpheres(e,t,n,s,r){const a=this.r.device.queue,i=s.length/4;n.forEach((h,l)=>{const c=Math.min(de,i-l*de);a.writeBuffer(e,h*de*16,s.buffer,s.byteOffset+l*de*16,c*16),a.writeBuffer(t,h*16,new Uint32Array([c,r,1,0]))})}onDetail(e){const t=this.loaded.get(e.idx);if(!t)return;t.detPending=!1;const n=t._req;this.lipA.free(t.lipPages),this.clearSPages(this.r.lipPagesBuf,t.lipPages),this.freeA.free(t.detPages),this.clearSPages(this.r.freePagesBuf,t.detPages),t.lipPages=[],t.detPages=[];const s=e.lipids,r=e.free,a=s.length/4,i=r.length/4;if(a){const l=this.lipA.alloc(Math.ceil(a/de))??this.freeDetailSpace(this.lipA,Math.ceil(a/de),t.idx);l&&(t.lipPages=l,this.writeSpheres(this.r.lipPool,this.r.lipPagesBuf,l,s,t.slot))}if(i){const l=this.freeA.alloc(Math.ceil(i/de))??this.freeDetailSpace(this.freeA,Math.ceil(i/de),t.idx);l&&(t.detPages=l,this.writeSpheres(this.r.freePool,this.r.freePagesBuf,l,r,t.slot))}t.nLip=t.lipPages.length?a:0,t.nDet=t.detPages.length?i:0,t.det={n:n.n,a:n.a,b:n.b,shell:n.shell};const h=t.nLip>0;h!==t.underMesh&&t.nIdx&&(t.underMesh=h,this.r.setMeshSlot(t.slot,h)),this.detailStats.lipids+=e.stats.lipCand,this.syncHigh(),this.onChange()}freeDetailSpace(e,t,n){const s=[...this.loaded.values()].filter(r=>r.idx!==n&&(r.lipPages.length||r.detPages.length)&&!this.desiredSet.has(r.idx));for(const r of s){this.dropDetail(r);const a=e.alloc(t);if(a)return a}return null}meshDraws(){const e=[];for(const t of this.loaded.values())t.nIdx&&t.vb&&t.ib&&e.push({slot:t.slot,vb:t.vb,ib:t.ib,nIdx:t.nIdx,under:t.underMesh});return e}instanceAt(e){const t=Math.floor(e/ce),n=this.pageTile[t];if(n<0)return null;const s=this.loaded.get(n);if(!s)return null;const r=this.pageLocal[t]+e%ce;return{t:s,local:r,...this.decode(s,r)}}decode(e,t){const n=e.inst,s=n[4*t],r=n[4*t+1],a=n[4*t+2],i=360/65535,h=[e.lo[0]-80+(r&65535)*i,e.lo[1]-80+(r>>>16)*i,e.lo[2]-80+(a&65535)*i];return{type:s&16383,region:s>>>16&31,flags:s>>>21,p:h}}pinned=new Set;unpin(){this.pinned.clear()}ensureLoaded(e,t=9e4){for(const a of e)this.pinned.add(a);const n=()=>e.filter(a=>!this.loaded.has(a)),s=()=>{for(const a of n())this.state[a]===0&&(this.state[a]=1,this.inflight++,this.workers[a%this.workers.length].postMessage({cmd:"load",idx:a,ix:this.idx.ix[a],iy:this.idx.iy[a],iz:this.idx.iz[a]}))},r=performance.now();return new Promise((a,i)=>{const h=()=>{if(!n().length)return a();if(performance.now()-r>t)return i(new Error(`${n().length} tiles of the box did not load${this.lastLoadErr?": "+this.lastLoadErr:""}`));s(),setTimeout(h,50)};h()})}poolIndex(e,t){return e.pages[Math.floor(t/ce)]*ce+t%ce}tileBySlot(e){for(const t of this.loaded.values())if(t.slot===e)return t;return null}countLoaded(e){let t=0;for(const n of this.loaded.values()){const s=n.inst;for(let r=0;r<n.nInst;r++)(s[4*r]&16383)===e&&t++}return t}nearestInTile(e,t,n){let s=-1,r=1e18;for(let a=0;a<e.nInst;a++){if((e.inst[4*a]&16383)!==t)continue;const i=this.decode(e,a),h=(i.p[0]-n[0])**2+(i.p[1]-n[1])**2+(i.p[2]-n[2])**2;h<r&&(r=h,s=a)}return s}totals(){let e=0,t=0,n=0,s=0;for(const r of this.loaded.values())e+=r.nInst,t+=r.nLip,n+=r.nBase+r.nDet,s+=r.nIdx/3;return{tiles:this.loaded.size,inst:e,lip:t,free:n,mesh:s,inflight:this.inflight+this.detInflight,desired:this.desired.length}}pendingForView(){let e=0;for(const n of this.desired)this.state[n]!==2&&e++;let t=0;for(const n of this.loaded.values())n.detPending&&t++;return e+t}}class Zs{G;PT;PQ;PC;pores;c;grid;n3;F;constructor(e){this.F=e.F??null,this.G=e.G,this.PT=e.PT,this.PQ=e.PQ,this.PC=e.PC,this.pores=e.pores,this.c=e.consts,this.grid=e.grid,this.n3=e.grid.n**3}cell(e,t,n){const s=this.grid.n,r=this.grid.ext,a=2*r/s,i=h=>Math.min(s-1,Math.max(0,Math.floor((h+r)/a)));return(i(e)*s+i(t))*s+i(n)}nucSdf(e,t,n){if(It(this.F)>Vn)return cn(this.F,e,t,n,Vn);const s=this.G,r=e-s[3],a=t-s[4],i=n-s[5],h=(s[14]*r+s[15]*a+s[16]*i)/s[11],l=(s[17]*r+s[18]*a+s[19]*i)/s[12],c=(s[20]*r+s[21]*a+s[22]*i)/s[13],d=Math.hypot(h,l,c),f=Math.hypot(h/s[11],l/s[12],c/s[13]);return f<1e-12?-Math.min(s[11],s[12],s[13]):d*(d-1)/f}sdfCer(e,t,n,s,r){const a=this.PQ,i=8*e,h=Math.hypot(t,n,s)||1e-9,l=Math.min(1,Math.max(-1,(t*a[i]+n*a[i+1]+s*a[i+2])/h)),c=Math.max(Math.acos(l)-a[i+4],0)*h,d=this.G[2]-a[i+5]-r;return Math.sqrt(c*c+d*d)-a[i+6]}sdfPrim(e,t,n,s){const r=this.PT[e],a=this.PQ,i=8*e;if(r===0)return Math.hypot(t-a[i],n-a[i+1],s-a[i+2])-a[i+3];if(r===1){const w=t-a[i],b=n-a[i+1],m=s-a[i+2],u=a[i+3]-a[i],p=a[i+4]-a[i+1],x=a[i+5]-a[i+2],y=Math.min(1,Math.max(0,(w*u+b*p+m*x)/(u*u+p*p+x*x)));return Math.hypot(w-u*y,b-p*y,m-x*y)-a[i+6]}if(r===5){const w=t-a[i],b=n-a[i+1],m=s-a[i+2],u=w*a[i+3]+b*a[i+4]+m*a[i+5],p=Math.sqrt(Math.max(w*w+b*b+m*m-u*u,0));return Math.hypot(p-a[i+6],u)-a[i+7]}if(r===4){const w=.5*(a[i]+a[i+3]),b=.5*(a[i+1]+a[i+4]),m=.5*(a[i+2]+a[i+5]),u=Math.hypot(w,b,m),p=w/u,x=b/u,y=m/u;let A=a[i+3]-a[i],M=a[i+4]-a[i+1],_=a[i+5]-a[i+2];const S=Math.hypot(A,M,_);A/=S,M/=S,_/=S;const g=t-a[i],k=n-a[i+1],P=s-a[i+2],T=Math.min(Math.max(g*A+k*M+P*_,0),S),E=Math.min(Math.max(g*p+k*x+P*y,0),a[i+7]);return Math.hypot(g-T*A-E*p,k-T*M-E*x,P-T*_-E*y)-a[i+6]}if(r===2){const w=t-a[i],b=n-a[i+1],m=s-a[i+2],u=a[i+3],p=a[i+4],x=a[i+5],y=w*u+b*p+m*x,A=Math.hypot(w-y*u,b-y*p,m-y*x),M=Math.max(A-a[i+6],0);return Math.sqrt(M*M+y*y)-a[i+7]}const h=Math.hypot(t,n,s)||1e-9,l=Math.min(1,Math.max(-1,(t*a[i]+n*a[i+1]+s*a[i+2])/h)),c=Math.acos(l),d=Math.max(c-a[i+4],0)*h,f=h-a[i+5];return Math.sqrt(d*d+f*f)-a[i+6]}label(e,t,n){const s=this.G,r=this.c.PERI_W,a=this.c.PORE_R,i=Rn(this.F,s,e,t,n),h=i-(s[0]-s[2]),l=this.cell(e,t,n),c=this.grid;let d=1e9;for(let E=c.prim_off[l];E<c.prim_off[l+1];E++){const I=c.prim_idx[E];this.PC[I]===9&&(d=Math.min(d,this.sdfPrim(I,e,t,n)))}const f=this.c.SCAR_K??30,w=Math.max(f-Math.abs(-i-d),0)/f,b=-(Math.min(-i,d)-w*w*f*.25);if(b<0)return 0;if(b<s[0]-s[1])return 1;if(h<0)return 2;let m=1e9;const u=this.grid,p=[1e9,1e9,1e9,1e9,1e9,1e9,1e9,1e9,1e9];It(this.F)>jn&&(p[3]=cn(this.F,e,t,n,jn));let x=1e9,y=1e9;for(let E=u.prim_off[l];E<u.prim_off[l+1];E++){const I=u.prim_idx[E],F=this.PT[I]===3?this.sdfCer(I,e,t,n,h):this.sdfPrim(I,e,t,n),N=this.PC[I];N===7?p[1]=Math.min(p[1],F):N===0?p[3]=Math.min(p[3],F):N===4?p[4]=Math.min(p[4],F):N===5?x=Math.min(x,F):N===6?y=Math.min(y,F):N===2?p[6]=Math.min(p[6],F):N===3?p[7]=Math.min(p[7],F):N===1?p[8]=Math.min(p[8],F):N===8&&(m=Math.min(m,F))}const A=this.c.EIS_K??8,M=Math.max(A-Math.abs(h-m),0)/A;if(Math.min(h,m)-M*M*A*.25<0)return 2;const _=Math.max(x,-y);if(p[1]<0)return 4;const S=this.nucSdf(e,t,n),g=S+.5*r,k=a+.5*r;let P=!1;for(let E=u.pore_off[l];E<u.pore_off[l+1];E++){const I=6*u.pore_idx[E],F=this.pores,N=e-F[I],G=t-F[I+1],$=n-F[I+2],Z=N*F[I+3]+G*F[I+4]+$*F[I+5];if(Math.abs(Z)>150)continue;if(Math.hypot(N-Z*F[I+3],G-Z*F[I+4],$-Z*F[I+5])<k){P=!0;break}}if(Math.abs(g)-.5*r<0&&!P)return 5;if(S<-.5*r){if(P&&S>-r-8)return 14;const E=e-s[7],I=t-s[8],F=n-s[9];return S<-r&&E*E+I*I+F*F<s[10]*s[10]?7:6}return S<0?14:p[4]<0?_<0?10:9:p[3]<0?8:p[8]<0?11:p[6]<0?12:p[7]<0?13:3}}const ei=["PI","PE","PC","PA","PS","ergosterol","other (sphingolipids, PG, DAG, ...)","cardiolipin"],dn=[[[16.8,20.3,17.7,33.6,.2,3.9,6.9],.4/.6],[[51.3,33.4,7.5,6.6,.4,.2,.5],.26],[[44.6,26.9,15.1,5.9,.5,2.2,4.3],.64],[[46.5,19.4,18.3,4.4,1.6,2.1,7.7],.18],[[45.6,32.6,10.2,1.2,5.9,4.4,0],.01],[[38.4,24,16.2,3.8,16.1,1.5,0],.32],[[35,22.3,19.1,12.9,.7,1.2,8.8],.75],[[48.2,22.9,15.8,4.5,7,1.6,0],.39],[[51.3,33.4,7.5,6.6,.4,.2,.5],0]],Jn={4:1,1:.8,6:0};function Ro(o,e){const[t,n]=dn[Math.min(o,dn.length-1)],[s,r,a,i,h,l,c]=t,d=t.reduce((m,u)=>m+u),f=1/(1+n),w=n/(1+n);let b=[a,r,s,l,i,0,c,h].map(m=>m/d*f*100);if(b[5]=w*100,o===0){b=b.map((u,p)=>p===5?u:u*2*(e?Jn[p]??.5:1-(Jn[p]??.5)));const m=b.reduce((u,p,x)=>u+(x===5?0:p),0);b=b.map((u,p)=>p===5?u:u/m*(100-b[5]))}return b}dn.map((o,e)=>[0,1].map(t=>{const n=Ro(e,t),s=n.reduce((a,i)=>a+i);let r=0;return n.map(a=>r+=a/s)}));const ti=`
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
<li><b>35.8 million macromolecule instances (most multi-chain assemblies), 5,615 molecular species.</b> Protein abundances are from the PaxDb integrated dataset rescaled to
5.0 pg of protein per cell (measured 4–6 pg), 7.8 × 10⁷ protein molecules in a 42 µm³ cell. There are 301,000 80S ribosomes, 28,000 pre-ribosomes, 3.3 million tRNAs, 78,000 nucleosomes and 117 nuclear pore complexes.</li>
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
25,600 of them carry 121,000 of the 80S ribosomes (about 5 per mRNA): the mRNA runs through each ribosome's channel
(path from PDB 7ZW0), short linkers join neighbouring ribosomes, and the rest of the transcript leaves the first and
last ribosome as 3&prime; and 5&prime; tails. The other 4,400 mRNAs carry no ribosomes. Threads that run out of room end
early (86% of the target length is reached overall). Linker DNA joins consecutive nucleosomes.</li>
<li><b>ER-bound ribosomes:</b> 8,100 ribosomes sit on translocons in the ER and outer nuclear membrane: on the Ssh1 complex
(yeast ribosome-Ssh1 cryo-EM structure 2WW9) or on Sec61 (7KAH, placed through the shared Sss1 subunit). Their number is
limited by the free translocon pools (80% assumed ribosome-bound); the PM-facing side of the cortical ER is too narrow for ribosomes.
A further 15,000 free ribosomes sit next to the ER and outer nuclear membrane, bringing the density to the 550&ndash;900 per &micro;m&sup2;
seen in tomograms (West <i>et al.</i> 2011).</li>
<li><b>RAC on ribosomes:</b> all 91,000 ribosome-associated complexes (Hsp40 Zuo1 with Hsp70 Ssz1) sit on ribosomes at the
tunnel exit, posed as in the 80S&ndash;RAC structure 7Z3N (<i>C. thermophilum</i> RAC standing in for the yeast one).</li>
<li><b>Chromosomes:</b> the Rabl configuration: 16 centromeres at the spindle pole body, each of the 32 arms a nucleosome chain
running to its telomere at the nuclear envelope; telomeres gather in 7 foci, longer arms reaching farther from the pole
(telomere distance follows arm length, Spearman 0.71); chromosome XII runs through the nucleolus. 550 transcribing Pol II sit
on their genes (per-gene counts from run-on, Pelechano <i>et al.</i> 2010; nucleosomes evicted around each), DNA through
the cleft and nascent mRNA from the exit.</li>
<li><b>mtDNA:</b> 15 nucleoids in the matrix holding 35 copies of the 85.8 kb mitochondrial genome (3 Mb of DNA).</li>
<li><b>Cvt complex:</b> one cluster of precursor aminopeptidase I dodecamers, Ams1 and the receptor Atg19, 180 nm across, next
to the vacuole (the cytoplasm-to-vacuole targeting cargo).</li>
<li><b>Ribosome assembly:</b> 15,000 pre-ribosomes (about 9% of ribosomes) follow the assembly pathway from the nucleolus
to the cytoplasm: 90S particles (SSU processome, 5WLC) and nucleolar pre-60S (7R7A) in the nucleolus, Rix1&ndash;Rea1
pre-60S (6YLH) in the nucleoplasm, and Nmd3&ndash;Lsg1 pre-60S (5T62) and pre-40S (6FAI) in the cytoplasm, each with its
pre-rRNA. Their numbers follow the abundances of the assembly factors specific to each stage. In the nucleolus, 75 active
rRNA genes carry 50 RNA polymerase I each on nucleosome-free DNA, with nascent pre-rRNA growing longer along each gene and
90S particles on uncleaved transcripts (the "Christmas trees" of Miller spreads); the inactive repeats are nucleosomal.
The <i>rRNA genes</i> guided view shows them on their own.</li>
<li><b>Microtubules:</b> 21 microtubules (11,700 &alpha;/&beta;-tubulin dimers) grow from &gamma;-tubulin rings on the spindle pole body,
opposite the nucleolus: 16 short kinetochore microtubules whose plus ends anchor the 16 chromosomes (Rabl configuration),
2 more nuclear ones and 3 cytoplasmic ones about 1 &micro;m long that curve around organelles. The 13-protofilament lattice
and the ring come from the native yeast &gamma;-TuRC-capped microtubule 8QV3; the lattice dimer is 5W3F (8.3 nm repeat).</li>
<li><b>Mitochondrial cristae:</b> 11,900 of the 21,800 ATP synthase copies form 5,900 dimers in rows along the rims of the lamellar
cristae, every 13 nm. Each dimer is two copies of the yeast monomer 8F39 placed in the dimer arrangement of 4B2Q
(cryo-ET), whose 84&deg; angle matches the curvature of the rim (crista lumen 18 nm wide). Complex IV is entirely in
2,900 respiratory supercomplexes III<sub>2</sub>IV<sub>2</sub> (6HU9) in the inner membrane, next to 2,600 free complex III dimers.</li>
<li><b>Actin cables:</b> 15 bundles of filaments 12 nm apart (62,000 actin subunits, 171 &micro;m), run 90 nm beneath the
plasma membrane, built from the yeast F-actin cryo-EM structure 9GO5 (167&deg; and 2.76 nm per subunit). Both strands of
every filament are coated head-to-tail by 11,650 tropomyosin TPM1 dimers (200 residues, one per 5 actins), following the
tropomyosin path of the mammalian actin-tropomyosin structure 5JLF wound as an ideal coiled coil.</li>
<li><b>Actin patches:</b> 41 endocytic sites, each a plasma-membrane invagination (a tube about 40 nm wide, 0&ndash;100 nm deep
at random stages; Kukulski <i>et al.</i> 2012) inside a branched network of about 1,500 actin subunits: short filaments branch at 73&deg; from Arp2/3 complexes (4,200 branches; geometry from the branch-junction
structure 8UXW) and are capped by capping protein (8F8Q), barbed ends toward the membrane. The yeast Arp2/3 complex
(all 7 subunits) and capping protein are AlphaFold models superposed on those templates. A further 4,800 Arp2/3
complexes are free in the cytosol. Assumption: 30% of actin in cables, 30% in patches, the rest monomeric
(much of it bound to profilin, cofilin or Srv2).</li>
<li><b>Endocytic sites</b> also carry their coat at the tip (Sla2, Sla1, End3; numbers per site from Picco <i>et al.</i> 2015)
and amphiphysin Rvs161/167 on the tubes deeper than 50 nm.</li>
<li><b>Nucleus&ndash;vacuole junction:</b> 758 Vac8&ndash;Nvj1 complexes (5XJG) where a vacuole lobe lies against the nucleus.</li>
<li><b>Eisosomes:</b> 27 plasma-membrane furrows, about 300 nm long and 50 nm deep, whose half-pipe bottoms are lined by
4,700 Pil1 and Lsp1 BAR-domain dimers in their native lattice (8QB7, 8QB8; Kefauver <i>et al.</i> 2024), amphipathic
helices in the inner leaflet. The coats hold 6.7% of the cell's Pil1 and Lsp1; the rest stay free near the membrane.
The MCC markers Sur7 and Nce102 crowd the furrow walls.</li>
<li><b>Membranes</b> (plasma membrane, cortical and perinuclear ER, nuclear envelope with pores, vacuole, a reticular
mitochondrion, Golgi/endosomes, peroxisomes, lipid droplets) are represented by signed-distance fields. The cell
outline, vacuoles and lipid droplets are those of a cell segmented from a cryo-FIB/SEM volume (EMPIAR-11416, Dumoux
<i>et al.</i> 2023), the nucleus that of a second cell; the shapes are smoothed, and the vacuoles' thickness along the
milling direction, which the segmentation loses, is reconstructed.
Near the camera every lipid is drawn individually, with the composition of each membrane from isolated organelle
membranes (Zinser <i>et al.</i> 1991), cardiolipin in the inner mitochondrial membrane and PS/PE in the cytoplasmic leaflet
of the plasma membrane; 0.65 nm² per lipid. Lipids are removed wherever they would overlap a membrane protein.</li>
<li><b>Vacuolar polyphosphate:</b> the vacuole lumen holds little protein but is packed with short polyphosphate chains
(about 300 mM phosphate units, an assumption within the reported "hundreds of millimolar"), the cell's phosphate store,
drawn procedurally when you zoom in.</li>
<li><b>Cell wall:</b> an inner β-1,3-glucan layer and an outer mannoprotein layer with mannan brushes. Bead counts follow a
wall mass of about 20% of dry mass. The cell is drawn as a mother that has budded once: one bud scar, a chitin ring 1 &micro;m
across whose rim (half the wall's chitin, as microfibrils) stands about 50 nm above the surface.</li>
</ul>

<h3>How it was built</h3>
<p>Proteins were assigned to compartments from UniProt subcellular-location annotations, with GO as fallback, and measured shares from the yeast GFP collection (CYCLoPs) and mitochondrial proteomics (Vögtle 2017, Morgenstern 2017). Membrane proteins were oriented by their
transmembrane segments. Everything was then packed into 200 nm tiles with collision detection on residue-level
sphere trees. Placement is not purely random: proteins with similar STRING network neighbourhoods (personalized
PageRank profiles, no clustering) are co-enriched in regions a few hundred nm across, and about a quarter of all
molecules were placed next to an interaction partner predicted by AlphaFold 3 or a close network neighbour (an assumed,
tunable organisation within each compartment; see Model details). Random addition jams below cellular crowding, so
tiles that could not fit everything were relaxed: overlapping molecules pushed apart as rigid bodies and the missing copies
added in rounds. Validation on the final model:</p>
<ul>
<li>Bead clash rate 9×10⁻⁶.</li>
<li>Containment violations 1.4×10⁻⁴.</li>
<li>Macromolecule densities (protein plus RNA and DNA; lipids and small molecules not counted): cytosol 163 mg/mL
(127 protein + 36 RNA), nucleoplasm 103 mg/mL, nucleolus 137 mg/mL, mitochondrial matrix 107 mg/mL. With lipids and
small molecules this comes to roughly the measured dry mass of the cytoplasm.</li>
<li>98.2% of requested copies were placed outside mitochondria (94.7% overall).</li>
</ul>
<p>Positions are stored to 0.0055 nm and orientations to 0.004°.</p>

<h3>Simulating a piece of the cell</h3>
<p><b>Simulate a box</b> in the side panel cuts a 30–200 nm cube out of the packed cell and runs Brownian dynamics of
its molecules on your GPU, in WebGPU compute shaders written for this viewer. The rest of the cell stops rendering meanwhile.
Each molecule is a rigid body made of collision spheres at one of three resolutions: its sphere tree (<i>schematic</i>, about
30 residues per sphere, microseconds), clusters of about 4 residues (<i>mid</i>), or every exposed residue (<i>residue</i>,
nanoseconds). Diffusion comes from its shape (Kirkwood hydrodynamic radius and Stokes–Einstein), with mean-field
hydrodynamics for the crowding. Molecules repel on overlap and interact by:</p>
<ul>
<li><b>Screened electrostatics</b> between charged residues and nucleotides, at the cytosol's salt and pH.</li>
<li><b>Specific binding</b> at the predicted interfaces of 371 AlphaFold dimers.</li>
<li><b>mRNA links</b> between the ribosomes of a polysome.</li>
<li><b>An optional nonspecific stickiness.</b> Without it, proteins diffuse about 5× faster than measured in yeast.
Stickiness can close that gap, as in bacterial cytoplasm models, but here it does so by clustering the molecules into
networks. A higher viscosity, standing in for the many-body hydrodynamics the model lacks, is the cleaner option.</li>
</ul>
<p>Membranes act as walls, and membrane proteins diffuse within them. Filaments and chromatin stay fixed. The box is
periodic, or embedded in the cell by a tethered outer shell when it contains membranes.</p>
<ul>
<li><b>Running:</b> watch it live, or compute a stretch first and play it back at a chosen frame rate.</li>
<li><b>Tracking:</b> click molecules to track them; the others fade and the tracked paths are drawn.</li>
<li><b>Readout:</b> compares the simulated diffusion with yeast measurements.</li>
<li><b>Export MD files:</b> saves
  <ul>
  <li>a coarse-grained mmCIF;</li>
  <li>a manifest with each molecule's structure source, pose, charges, links and contacts;</li>
  <li>the trajectory, with a script that turns it into a DCD for VMD or ChimeraX;</li>
  <li>a ready-to-run OpenMM script with the same interactions.</li>
  </ul></li>
</ul>
<p>See Model details, section 14.</p>

<h3>Honest caveats</h3>
<ul>
<li><b>Layout.</b> The organelle shapes and positions are a plausible illustration, not a tomogram of one specific cell.</li>
<li><b>Abundances.</b> Absolute abundances are uncertain to about 2-fold.</li>
<li><b>Full compartments.</b> Small, crowded compartments saturated: only 40% of the mitochondrial matrix proteins, 80% of the IMS, 77% of the nucleolus and 81% of the Golgi/endosomes fit, so their copies are under-represented.</li>
<li><b>Histones</b> beyond the fixed number of nucleosomes are not removed: 80,000 surplus H2B are left out and about 50,000 free H4 and histone dimers sit in the nucleoplasm, although cells keep free histones low.</li>
<li><b>Chromatin</b> follows the Rabl layout at the level of arms (kb mapped linearly along each walk); there are no
sequence-level nucleosome positions or inter-chromosomal contacts. Linker DNA and mRNA paths are schematic.</li>
<li><b>Procedural content.</b> Lipids, β-glucan and mannan are generated procedurally in your browser from deterministic seeds. They are not stored
coordinates, and polysaccharides are only partly clash-checked.</li>
<li><b>Spindle pole body.</b> Only the &gamma;-tubulin rings are modelled, not the SPB plaques; PaxDb lists fewer SPC98 copies than 21 rings need.</li>
<li><b>Not included:</b> septins, fimbrin cross-links, the ribosome-bound chaperones NAC and Ssb, metabolites and water.
The FG nucleoporins absent from 7N9F (NUP100, NUP116, NUP145N, NUP42 and the basket NUP1, NUP60, NUP2) have no FG domains.</li>
<li><b>Dropped copies.</b> 15,982 transmembrane-protein copies assigned to soluble compartments were dropped, and 90,403 copies had no structure (80,000 of them the surplus histone H2B).</li>
</ul>
<p style="color:#a9a69e">Renderer: WebGPU compute culling and level-of-detail selection, sphere impostors with exact cut caps, analytic
compartment backgrounds and Goodsell-style outlines. Full-atom models are fetched live from RCSB PDB and AlphaFold DB
when you request them.</p>
<p style="color:#a9a69e">Data sources: RCSB Protein Data Bank (CC0) and the individual entries named above; AlphaFold Protein Structure
Database monomers and predicted complexes (CC-BY 4.0, Jumper <i>et al.</i> 2021; Varadi <i>et al.</i> 2024); PaxDb protein abundances;
Complex Portal; UniProt and Gene Ontology annotations; IntAct interactions; STRING v12 associations (CC-BY 4.0); EMPIAR-11416 cryo-FIB/SEM volume (CC0); Alliance of Genome Resources paralogs;
lipid composition from Ejsing <i>et al.</i> 2009; total protein from Ho <i>et al.</i> 2018.</p>
`,hn=o=>o.replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;");function Mt(o){const e=[];return o=o.replace(/`([^`]+)`/g,(t,n)=>`\0${e.push(hn(n))-1}\0`),o=hn(o).replace(/\[([^\]]+)\]\(([^)\s]+)\)/g,(t,n,s)=>`<a href="${s}" target="_blank" rel="noopener">${n}</a>`).replace(/\*\*([^*]+)\*\*/g,"<b>$1</b>").replace(/\*([^*\s][^*]*)\*/g,"<i>$1</i>"),o.replace(/\u0000(\d+)\u0000/g,(t,n)=>`<code>${e[+n]}</code>`)}const Qt=/^([-*]|\d+\.)\s+/,Pt=o=>o.length-o.trimStart().length,Xn=o=>o.trim().replace(/^\||\|$/g,"").split("|").map(e=>e.trim());function To(o){let e="",t=0;for(;t<o.length;){const n=o[t],s=n.trim();if(!s){t++;continue}let r;if(r=s.match(/^(#{1,6})\s+(.*)$/)){e+=`<h${r[1].length}>${Mt(r[2])}</h${r[1].length}>`,t++;continue}if(/^-{3,}$/.test(s)){e+="<hr>",t++;continue}if(s.startsWith("```")){const i=[];for(t++;t<o.length&&!o[t].trim().startsWith("```");)i.push(o[t++]);t++;const h=Math.min(...i.filter(l=>l.trim()).map(Pt),99);e+=`<pre><code>${hn(i.map(l=>l.slice(h)).join(`
`))}</code></pre>`;continue}if(s.startsWith("|")&&t+1<o.length&&/^\|[\s|:-]+\|$/.test(o[t+1].trim())){const i=Xn(s);t+=2;let h="";for(;t<o.length&&o[t].trim().startsWith("|");)h+=`<tr>${Xn(o[t]).map(l=>`<td>${Mt(l)}</td>`).join("")}</tr>`,t++;e+=`<table><thead><tr>${i.map(l=>`<th>${Mt(l)}</th>`).join("")}</tr></thead><tbody>${h}</tbody></table>`;continue}if(r=s.match(Qt)){const i=/\d/.test(r[1]),h=Pt(n);let l="";for(;t<o.length;){const c=o[t],d=c.trim().match(Qt);if(!d||Pt(c)!==h||/\d/.test(d[1])!==i)break;const f=h+d[0].length,w=[c.slice(f)];for(t++;t<o.length&&(!o[t].trim()||Pt(o[t])>=f);)w.push(o[t].slice(f)),t++;for(;w.length&&!w[w.length-1].trim();)w.pop();const b=To(w),m=b.match(/^<p>((?:(?!<p>)[\s\S])*)<\/p>$/);l+=`<li>${m?m[1]:b}</li>`}e+=i?`<ol>${l}</ol>`:`<ul>${l}</ul>`;continue}const a=[];for(;t<o.length&&o[t].trim()&&!/^(#{1,6}\s|\||-{3,}$|```)/.test(o[t].trim())&&!Qt.test(o[t].trim());)a.push(o[t++].trim());e+=`<p>${Mt(a.join(" "))}</p>`}return e}function ni(o){return To(o.replace(/\r\n/g,`
`).split(`
`))}const oi=`# Whole yeast cell: model description

This document describes how the structural model shown in the viewer was built. It covers what the model contains, which data and structures it comes from, the assumptions and design decisions behind it, how it was validated, and where it is known to be wrong or incomplete. It is written for readers who want to judge the model scientifically or reuse its numbers.

All numbers refer to the current build (build S of 2026-11-13, candidate build 3 of round S: commit 9568da3, master seed 0, build id \`20261113-155407-9568da3-s0\`; 5,562 molecular species, 24.9 million packed instances, 25.3 million with the pre-placed structures) unless stated otherwise. Build S keeps the cell and the protein copy numbers of build R (3.5 pg of protein in 34 µm³ instead of 5.0 pg in 42 µm³; section 2.1); its inventory differs only in the nucleus/cytosol split, which now uses the compartment volumes as built (section 4.2), and in free H2B, which now has a structure (section 2.3). It changes the vacuole (one near-spherical body instead of 12 lobes), the mitochondria, the side of the nucleolus, the pores at the nucleus–vacuole junction and the packing algorithm (sections 3.3, 6.7, 6.10 and 8), and it is the first build made in an isolated run folder with a build stamp (section 13). Where a build R number is informative it is given as "(build R: x)"; numbers marked "build R" alone were not re-measured on build S. Numbers quoted for builds up to Q refer to the larger cell. Build S is one build at one seed: the spread between replicate builds is not yet measured (section 11).

## 1. What the model is

- **The cell.** One haploid *Saccharomyces cerevisiae* cell (strain BY4741 / S288C background) growing exponentially in rich glucose medium (YPD). The cell is unbudded, and the microtubules and chromatin are arranged as in G1. Since build R its size and protein content are also those of one G1 cell with one genome, not of the average cell of a growing culture (section 2.1).
- **What a molecule is.** Every macromolecule is an individual copy of an experimental (PDB) or predicted (AlphaFold DB) structure. It is drawn at one bead per amino-acid residue (Cα) or per nucleotide (P or C4′). There are 18.3 billion beads in all (build R: 17.9 billion).
- **How the cell was assembled.** Copy numbers come from proteomics. Compartment volumes and membrane areas come from imaging literature. The organelles are built as analytic geometry, and molecules are packed into them with residue-level collision checks.
- **What it is not:**
  - **Not a dynamic or physical simulation.** There are no forces, no diffusion and no binding equilibria. Except where stated below, molecules are placed at random positions and orientations by random sequential addition. A snapshot from a real cell would show interaction-driven clustering that the model lacks.
  - **Not a reconstruction of a particular cell.** The organelle volumes and areas are taken from measurements, but their shapes and positions are a plausible arrangement, not a tomogram.
- **Lineage.** The approach follows the cellPACK-style whole-cell models of *Mycoplasma genitalium* (Maritan *et al.* 2022). The visual style follows David Goodsell.

## 2. Composition

### 2.1 Proteins

- **Abundances:** PaxDb integrated whole-organism dataset for *S. cerevisiae* (taxon 4932), in ppm, matched to UniProt by ORF.
- **Absolute scale (2026-11-12, build R):** ppm values are scaled so that the proteins weigh 3.5 pg. With the abundance-weighted mean mass of 38.7 kDa this is 5.45 × 10⁷ protein molecules.
  - **Why not the bulk value:** bulk measurements give 4–5 pg (Carroll *et al.* 2011) and 6 pg (Sherman 2002) per haploid cell, and builds up to Q used 5.0 pg in the culture's mean cell volume of 42 µm³. Those are averages over a growing culture, which includes budded cells with up to two genomes. The model cell is one unbudded G1 cell with one genome, a fixed number of pores and one spindle pole body.
  - **Anchors (\`validation/anchors.py\`):** structures counted per cell say how much protein this cell holds. At 5.0 pg the core histones were 1.3 times what one genome's nucleosomes need (so 3.9 pg), and the scaffold and FG nucleoporins 1.5 times what the pores need (3.2 pg). 3.5 pg lies between the two.
  - **Check, not used:** the scale gives 211,500 ribosomes. Counts from rRNA are about 200,000 (Warner 1999) and 187,000 ± 56,000 (von der Haar 2008); a growth-law estimate for YPD is about 3 × 10⁵ (Kostinski & Reuveni 2021).
  - **At 3.5 pg:** core histones are 0.9 times the nucleosomes' need (median of four; H3 0.70, H2B 1.31), nucleoporins 1.36 times the need of 93 pores.
  - **The cell shrinks with it** (section 3.1), so that concentrations stay where they were: the whole-cell protein concentration is 103 mg/mL (119 before), near the 100 mg/mL commonly used (Ohsawa *et al.* 2026 and references there), above a recent Bradford-based 71 mg/mL (Kunzi *et al.* 2026, preprint).
- **Earlier total:** earlier versions used the unified total of Ho, Baryshnikova and Brown (2018), 4.2 × 10⁷ molecules. That comes to only 2.7 pg of protein and left the cytosol at about 90 mg/mL of macromolecules.
- **Uncertainty:**
  - PaxDb integrated values are not specific to BY4741 in YPD. Per-protein uncertainty is roughly 2-fold.
  - The absolute total is still the largest single uncertainty: the anchors span 3.2–3.9 pg, and bulk measurements of average cells 4–6 pg.
- **Proteins this cell does not have (2026-11-12; \`inventory.NOT_IN_THIS_CELL\`):** PaxDb mixes strains, media and stages. Curated lists set 73 proteins to zero before scaling: the genes deleted in BY4741 (His3, Leu2, Met17, Ura3), α-specific proteins in this MATa cell (8), meiosis and sporulation proteins (42), glucose-repressed proteins (11: Gal1, Pck1, Mls1, Icl1, Cyc7 and others) and proteins degraded before G1 (8: Mcd1, Cdc5, Sgo1, Pds1, mitotic cyclins). Together 0.27% of the abundance signal. The lists treat a symptom; a condition-matched abundance set would replace the last three.

### 2.2 Nucleic acids and fixed counts

| Component | Copies | Basis |
|---|---|---|
| 80S ribosomes | 211,500 | Ribosomal protein abundances (median over the r-protein families, with Asc1 and Stm1, section 6.5), less the r-protein content of the pre-ribosomes |
| Pre-ribosomes | 19,500 | Abundances of stage-specific assembly factors (section 5.2) |
| tRNA | 2.31 × 10⁶ | 3.3 × 10⁶ per average cell (assumption, literature order of magnitude) × 0.7, the size of this cell relative to the average cell: about 11 per ribosome; split over 42 isoacceptors by tRNA gene copies (Percudani *et al.* 1997; section 6.5) |
| mRNA | 18,100 | Per gene from Pelechano *et al.* 2010 (Table S1: molecules per cell compiled from Nagalakshmi *et al.* 2008 and Miura *et al.* 2008; 25,900 per average cell over 5,085 genes) × 0.7; each mRNA is a transcript of a named gene (section 6.5) |
| Nucleosomes | 69,714 | 64,487 measured dyads on the chromosome arms (section 6.3; before the measured map: genome ÷ 165 bp = 73,150), plus 5,227 on the nucleosomal rDNA (section 6.6) |
| Nuclear pore complexes | 93 | 10.8 per µm² of nuclear surface: the density of 117 pores (Winey *et al.* 1997, mean over the cell cycle; 65–182, rising from G1) on the earlier, larger nucleus. An assumption, checked 2026-11-13 against Table 1 of that study: 10 G1 nuclei have 86 ± 16 pores (65–113), 11.7 per µm², so 93 is inside the G1 range |

Because the ribosome count follows the ribosomal-protein abundances, it inherits the uncertainty in the proteome total (section 2.1 for the measurements).

### 2.3 Not included

- **Small molecules and ions:** metabolites, free amino acids, nucleotides, ions and water.
- **Storage carbohydrates:** glycogen and trehalose.
- **Polyphosphate.**
- **Proteins without a structure:** 7,200 copies (build R: 40,300; 33,100 of those were free histone H2B, whose library entries HTB1 and HTB2 were missing until stage 0 rebuilt the library for build S, section 13).
- **Transmembrane proteins annotated only to soluble compartments:** 11,200 copies (see 4.3).
- **Surplus nucleoporins:** 13,000 copies beyond a free pool of 25% of the copies in pores (section 6.1).
- **Proteins this cell does not have:** section 2.1.

Lipids and the cell-wall polysaccharides are drawn in the viewer but generated procedurally (section 7).

## 3. Geometry

### 3.1 Volumes

The cell is 34 µm³ with a 110 nm cell wall, the lower bound of 110–200 nm (Orlean 2012). Its outline and nucleus are real segmented shapes; since build S the vacuole is one near-spherical body with the shape statistics measured by soft X-ray tomography of G1 cells (section 3.3).

**Cell volume (2026-11-12, build R; was 42 µm³, the mean log-phase volume of Jorgensen *et al.* 2002).** With 3.5 pg of protein (section 2.1) and the fractions below, 34 µm³ is the volume at which the cytosol keeps the macromolecular density it had (about 170 mg/mL), so the crowding and the diffusion calibration of section 14.5 stay valid. It lies between the culture mean (42 µm³) and unbudded G1 cells of the same strain and medium by soft X-ray tomography (24.9 µm³; Chen *et al.* 2025). The measurements do not agree on it: the dry-mass density of the cytoplasm (section 10.3) alone would put the cell at 36–38 µm³, and the ribosome concentration by cryo-ET (23 µM; the model has 17 µM) at about 26 µm³. The volume is a choice inside that range, not a measurement.

Organelle volume fractions:
- **Nucleus, vacuole and lipid droplets:** set against segmented cryo-FIB/SEM cells (EMPIAR-11416; Dumoux *et al.* 2023; section 3.3) and the literature:
  - **Nucleus 7% (build R; was 8%).** The nucleus-to-cell ratio is constant across cell sizes at about 7% (fluorescence; Jorgensen *et al.* 2007). Soft X-ray tomography of unbudded G1 BY4741 in YPD gives 5.9 ± 1.3% (Chen *et al.* 2025) and 5.3–6.3% in control strains (Ohsawa *et al.* 2026, a preprint); the two segmented nuclei are 8.0% and 9.7% of their cells (both from cells cut by the edge of the imaged volume, which biases the fraction up); serial-section EM gives 10.5% (Yamaguchi *et al.* 2011). The nucleus is 2.0 µm³ inside the inner membrane (2.9 before). The nucleolus is 20% of the nucleus (Uchida *et al.* 2011).
  - **Vacuole 17.2% (build R; was 9%).** Chen *et al.* 2025 (Mol Biol Cell 36:ar132): 17.2 ± 6.0% in 287 unbudded G1 cells of BY4741 in YPD by soft X-ray tomography, the measurement matched in strain, medium and stage. The five segmented cells give 6–13% (a lower bound: the segmentation misses faint vacuoles), an older soft X-ray study 7% (Uchida *et al.* 2011), fluorescence 25% (Joyner *et al.* 2016).
  - **Lipid droplets 0.3%.** The 14 segmented droplets of the cell whose outline is used (the five cells range 0.3–0.8%, from the depositors' hand-corrected labels); soft X-ray tomography gave 1%.
- **Mitochondria 3% (2026-11-11; was 2%, and 1% before that).** Measurements span 1.7% (serial-section EM; Yamaguchi *et al.* 2011) to 4.77% (3D fluorescence in glucose; Di Bartolomeo *et al.* 2020 PNAS 117:7524; CEN.PK in minimal medium, and diffraction inflates 300 nm tubules); soft X-ray tomography gave 1% (Uchida *et al.* 2011). At 2% the matrix asked for about 270 mg/mL of protein and only 40% fitted. 3% is where the requested matrix protein reaches 1.1× the cytosol's density, the ratio of soft X-ray absorption of mitochondria to cytosol (Chen *et al.* 2025): a fit to that ratio, inside the measured range. In build S, with the rebuilt mitochondria (section 3.3), the placed matrix density is 1.05 times the cytosol's (build R: 0.89).
- **Golgi, endosomes and peroxisomes:** 1% together (assumption).
- **Cortical ER:** covers 48% of the plasma membrane (Wolf *et al.* 2012).

Measured on the final voxelised geometry (Monte Carlo volumes of stage 2, build S; \`work/s2_geometry.log\`):

| Compartment | Volume (µm³) | Share of cell | Build R (µm³) |
|---|---|---|---|
| Cytosol | 18.38 | 54.1% | 18.42 |
| Cell wall (glucan 3.66, mannan 1.82) | 5.48 | 16.1% | 5.48 |
| Vacuole (1 compartment) | 5.83 | 17.2% | 5.81 (12 compartments) |
| Nucleoplasm | 1.56 | 4.6% | 1.56 |
| Nucleolus | 0.47 | 1.4% | 0.47 |
| Mitochondrial matrix | 0.82 | 2.4% | 0.70 |
| Mitochondrial intermembrane space (incl. crista lumen) | 0.18 | 0.5% | 0.27 |
| ER lumen: cortical ER | 0.45 | 1.3% | 0.45 |
| ER lumen: perinuclear space | 0.29 | 0.9% | 0.29 |
| Lipid droplets (14) | 0.12 | 0.3% | 0.12 |
| Golgi / endosomes | 0.27 | 0.8% | 0.26 |
| Peroxisomes | 0.07 | 0.2% | 0.07 |
| **Total** (incl. 0.03 µm³ of pore channels) | **33.95** | | 33.95 |

The tile measures used for packing differ from these by up to 0.01 µm³ (cytosol 18.37, vacuole 5.84, matrix 0.82 µm³).

### 3.2 Membranes

Membrane mid-surface areas, measured on the final geometry:

| Membrane | Area, build S (µm²) | Build R (µm²) |
|---|---|---|
| Plasma membrane (with furrows and invaginations) | 47.4 | 47.4 |
| ER (cortical sheets) | 46.6 | 46.5 |
| Nuclear envelope (both membranes) | 17.7 | 17.7 |
| Mitochondrial inner membrane (boundary membrane and cristae) | 13.1 | 31.5 |
| Mitochondrial outer membrane | 11.7 | 14.3 |
| Vacuole | 16.0 | 40.1 (12 lobes) |
| Golgi / endosomes | 8.2 | 8.6 |
| Peroxisomes | 1.3 | 1.3 |

The internal membranes total 115 µm² in build S against 160 µm² in build R: one round vacuole has 40% of the membrane of the carved lobes, and the rebuilt mitochondria have short cristae (section 3.3).

**Plasma-membrane invaginations** (\`cell3d.py\`, primitive \`T_SLOT\`). The plasma membrane is otherwise a fixed depth under the outer surface; these pockets are added to it by a smooth union (rim radius about 8 nm), and their lumen counts as cell wall.
- **Eisosome furrows (MCC):** 23 in build S (22 in build R; of 30 tried; they keep 300 nm apart, and the smaller cell has less free plasma membrane), about 300 nm long (lengths drawn from 300 ± 50 nm; mean 268 nm in build S, 272 nm in build R) and 50 nm deep (Strádalová *et al.* 2009; Karotki *et al.* 2011: about 30 per cell). The cross-section is a U whose bottom is a half-cylinder with the membrane radius of the native Pil1/Lsp1 coat (10.2 nm midplane, 31 nm coat; in situ 30–60 nm, Bharat *et al.* 2018). The coat is in section 6.8.
- **Endocytic invaginations:** one per actin patch (29). Each is a tube with a rounded tip, 40 nm wide at the membrane midplane (the 6,400 nm² mean area of the vesicles they release is a 45 nm sphere; Kukulski *et al.* 2012, 2016), 0–100 nm deep, drawn uniformly because the patches are at random stages and invaginations reach about 100 nm at scission. Measured tips are narrower and invaginations deeper than this (section 11).
- **Sites:** where there is no cortical ER, the membrane normal is within 8° of radial, and the nucleus and vacuole are clear. Furrows run in random directions in the membrane plane.

**Bud scar** (\`cell3d.py\`, primitive \`T_TORUS\`). The model cell is taken to be a mother that has budded once (assumption), so it carries one bud scar: the chitin ring left on the mother's wall by cytokinesis.
- **Shape:** a ring 1.0 µm across (AFM: about 1 µm; Pillet *et al.* 2014), modelled as a torus whose tube (radius 50 nm) is centred on the outer surface, so the rim stands 50 nm proud of it (the height is an assumption). The torus is lowered by the surface's sag under the ring, and it is joined to the outer surface only by a smooth union (30 nm), so the plasma membrane and the wall layers under it keep their depths. The AFM profile's convex centre and the depression inside the rim are not modelled.
- **Chitin:** the rim is filled with procedural chitin microfibrils (section 7).
- **Birth scar:** not modelled. A daughter's birth scar has no chitin rim and is only a faint feature of the surface.

### 3.3 Organelle shapes

**A composite of real cells.** The cell outline, the nucleus and the lipid droplets are real shapes from a cryo-plasma FIB/SEM volume of plunge-frozen log-phase yeast (EMPIAR-11416; Dumoux *et al.* 2023), segmented into meshes by a separate pipeline (Cellpose plus manual curation). Up to build R the vacuoles were too; since build S the vacuole is one ellipsoid placed near them (below):
- **Outline and lipid droplets:** cell 08, the only segmented cell not cut by the edge of the imaged volume (44.2 µm³, scaled by 0.916 in length to 34 µm³; 0.983 to 42 µm³ before build R). Its segmented vacuoles decide between equally good positions of the vacuole.
- **Nucleus:** cell 05's. Cell 08's nucleus was not segmented; cell 05's is one of two real ones, and the other (cell 04's) has its short axis along the milling direction of a cut cell. It is scaled to 7% of the cell volume and placed where the segmentation pipeline put cell 08's stand-in nucleus (the cell's most open region, where its unsegmented nucleus probably is), oriented along that stand-in's principal axes, and moved 118 nm from there to stay 260 nm below the outer surface.
- **Not used:** the bundle's procedural mitochondria, cortical ER and peroxisomes.

**Cleaning the meshes** (\`pipeline/cellshape.py\`; the meshes are in \`pipeline/cell_meshes.npz\`, extracted by \`scripts/cell_meshes.py\`):
- **Smoothing:** each part is voxelised on a 20 nm grid, closed morphologically (radius 100 nm for the outline, 60 nm for the rest), which fills grooves a curtaining stripe of the FIB data cut into the outline and some vacuoles. It is then Gaussian-smoothed (σ 50 nm outline, 40 nm nucleus, 30 nm vacuoles), which removes the mesh facets and the 50 nm slice steps.
- **One vacuole (build S, 2026-11-13; \`cellshape.VAC_MODE = "single"\`).** Cryo soft X-ray tomography of unbudded G1 BY4741 in YPD (Chen *et al.* 2025, Supplementary file S8: 490 cells at 30 nm voxels) finds one vacuole object in 488 cells, with a corrected sphericity of 0.98 (median) and a best-fit ellipsoid of elongation 1.10 and flatness 1.07 (medians); the same study gives the 17.2% volume fraction used since build R. The vacuole is therefore an ellipsoid with these axis ratios and 17.2% of the cell volume. Of 81 orientations (80 random), it is put where it needs the least carving to stay 220 nm below the outer surface and 30 nm from the nucleus (ties: nearest to the volume-weighted centroid of cell 08's segmented vacuoles), then carved where it still violates either, grown back to the target volume and smoothed. Result: one compartment of 5.84 µm³ (semi-axes 1.22, 1.11 and 1.04 µm before carving; 1.8% of it carved by the nucleus and the cortex; centre 1.2 µm from the centroid of the segmented vacuoles), 16.0 µm² of membrane, sphericity 0.98. 30 nm voxels cannot tell one vacuole from lobes pressed flat against each other, so the tomography does not exclude lobes; the carved lobes of build R (sphericity 0.39, below) are far from what it measures.
- **The segmented lobes (builds up to R; \`VAC_MODE = "lobes"\`, kept in the code):** the segmentation loses the tops and bottoms of vacuoles (z extent over in-plane extent 0.55–0.61 across the 107 segmented vacuoles, against 1.02 for the cells), so each lobe of cell 08 was stretched along z to its smaller in-plane extent (stretch 1.0–3.0, median 2.0), smoothed more along z (σ 90 nm), blended 1:1 with its moment-matched ellipsoid, and one streak along the curtaining stripe was dropped. The 15 lobes were scaled to 17.2% of the cell (build R: their segmented size, in an outline shrunk by 0.916), carved 30 nm apart from each other and from the nucleus and 220 nm from the outer surface, and smoothed (σ 20 nm). In build R that gave 12 compartments (5.81 µm³) with 40.1 µm² of membrane and a sphericity of 0.39: at 17% of a 34 µm³ cell the lobes collided and were carved into boxy pieces with slits between them, and one lobe wrapped another (visual check of build R, 2026-11-13).
- **Lipid droplets:** the 14 segmented droplets are kept as volume-equivalent spheres at their positions (median radius 124 nm); 13 that overlap the vacuole or the nucleus are moved to the nearest free space (build R: 8). Spheres also remove the droplets' z-blur (z extent / in-plane extent 1.25).

**Distance fields.** Each cleaned surface is re-meshed (marching cubes) and the exact distance to it is computed at the grid nodes near the surface (closest triangle). Between nodes the distance is interpolated linearly, accurate to about 0.2 nm, so membranes, offsets (the 35 nm mannan layer, the 110 nm wall, the cortical-ER sheets) and protein orientations are as precise as with the analytic shapes. The fields are stored in sparse 8³-voxel bricks within a band around each surface, and the viewer uses the same fields (re-bricked, 7.5 MB compressed). \`scripts/s2_geometry.py --sphere\` builds the previous analytic cell instead.

Remaining shapes and fixed dimensions:

| Element | Shape and dimensions |
|---|---|
| Nuclear envelope | Perinuclear space 34 nm between the membrane mid-surfaces, following the real nucleus. 93 pores of 45 nm radius at the mid-plane, at least 150 nm apart, spread uniformly over the envelope area outside the nucleus–vacuole junction (section 6.10) and set along the local normal. The nucleolus is a crescent on the side facing the vacuole (20% of the nucleus; below). |
| Cortical ER | Sheets 27 nm beneath the plasma membrane, following its real shape, with a 20 nm lumen. |
| Mitochondria (rebuilt for build S, 2026-11-13) | Curved tubules of capsule segments 110 nm long, started 457–597 nm and kept 437–637 nm below the outer surface, 1–2 separate networks per cell with side branches (build S: 2 networks, 7 chains, 5 branch points, 9 tips; 10.4 µm of tubule). Outer diameter varies along a tubule, 250–414 nm in build S (median 338 nm; 339 ± 5 nm in glucose by live 4Pi microscopy, Egner *et al.* 2002; median 330 nm in cryo-tomograms of BY4741 in YPD, Chang *et al.* 2025). Intermembrane space 14 nm between the membrane mid-surfaces (13–16 nm in those tomograms). 169 short lamellar cristae (16 per µm of tubule) with an 18 nm lumen (rim-tube radius 9 nm, chosen to fit ATP synthase dimers; section 6.7), each reaching 50–200 nm into the tubule (gamma-distributed around 85 nm; never across it) from alternating sides that turn slowly along the tubule. Crista membrane about 2.6 µm², 0.24 times the outer membrane (0.28 measured on the YPD meshes of Chang *et al.* 2025; 0.72 in glycerol). 11 crista-free nucleoid sites (section 6.7). Until build R: a dozen straight tubules of 150 nm radius, an 8 nm intermembrane space and a disc across the tubule every 42 nm (239 cristae; crista membrane 1.4 times the outer membrane, the look of a respiring cell). |
| Other organelles | Golgi/endosome cisternae and peroxisomes are simple primitives. |

**Arrangement.** The outline and lipid droplets keep their arrangement in cell 08 (13 of the 14 droplets moved out of the vacuole's way), and the vacuole sits where it needs the least carving, nearest the centroid of cell 08's vacuoles among such positions (1.2 µm from it). The nucleus sits in that cell's open region. The nucleolus faces the vacuole, and the spindle pole body, opposite the nucleolus, faces open cytoplasm; the mitochondria are tubules along the cortex. The last three are procedural.
- **Nucleolus next to the vacuole (build S, 2026-11-13).** The nucleolus usually lies at the nuclear envelope next to the vacuole (Barbosa *et al.* 2015; Santana-Sosa *et al.* 2023), and the spindle pole body lies opposite the nucleolus (Yang *et al.* 1989). The crescent is a sphere centred beyond the envelope in the direction of the vacuole centroid, sized to 20% of the nucleus. Up to build R it was on the side away from the vacuoles; with one large vacuole that put the spindle pole body against the vacuole, and no cytoplasmic microtubule could be built (round S candidate builds 1 and 2; build R had 3).

**Mitochondrial sources.** Measurements of glucose-grown cells collected for the rebuild are in \`pipeline/data/validation/literature/mitochondria/\` (README of 2026-11-13). Crista shape: "most cristae were rather sheet-like", "almost none of them spanned the mitochondria", mostly on two opposing sides of the tubule (FIB-SEM, galactose; Stoldt *et al.* 2019); glucose cells have relatively few, poorly defined cristae (Klecker & Westermann 2021). The crista density and reach were measured on the YPD meshes of EMPIAR-12534 (Chang *et al.* 2025): 24 crista patches per µm² of outer membrane, median depth 85 nm, maximum about 200 nm. Number of separate networks: median 2 per haploid cell (Seel *et al.* 2023).

**Shape statistics** of all five segmented cells (\`scripts/organelle_stats.py\` → \`pipeline/organelle_stats.json\`) were used to choose the composite's parts and to check it. The cell 08 outline is mildly non-spherical: its volume-moment semi-axes are 2.41, 2.22 and 2.00 µm.

**Known geometric deviations from EM.** Serial-section EM (Yamaguchi *et al.* 2011) gives a nucleus of 10.5% of cell volume; this model has 7%. Mitochondria are 1.7% there and 2.9% here (matrix plus intermembrane space, 0.99 µm³ in build S; 0.98 in build R).

## 4. Localization

### 4.1 Source annotations

**UniProt subcellular-location text** is mapped to compartments by ordered regular-expression rules. Examples:
- "nucleolus"
- "mitochondrion inner membrane"
- "endoplasmic reticulum"
- "cell cortex" → plasma membrane

Annotations are filtered by their evidence:
- **GFP data preferred:** where a location token is backed by the genome-wide GFP study of Huh *et al.* (2003), the tokens kept are the GFP-backed ones plus those with experimental evidence from a dedicated study. A paper cited as location evidence in 40 or more UniProt entries counts as high-throughput. Tokens are kept in UniProt's order, giving the bulk steady-state localization without proteomics-only calls. Examples: GPD1 keeps its cytoplasm, so it is two-thirds cytosolic; PNC1 is mostly cytosolic and nuclear with 14% in peroxisomes.
- **Proteomics-only calls dropped:** tokens supported only by two high-throughput mitochondrial proteomics papers (Sickmann 2003, Reinders 2006) are discarded. Those studies assign many abundant cytosolic enzymes to mitochondria.

**Fallbacks:**
- **GO cellular component:** used when UniProt has no location. Complex terms are skipped, except pores, ribosomes and proteasomes. If GO lists the cytosol, only cytosolic and nuclear terms are kept.
- **Cytosol:** the default for proteins with no annotation at all (625 proteins).

**Curated fix:** Pil1 goes to the plasma membrane. Its UniProt location is only "lipid droplet", from an organelle proteomics screen, while GO and the eisosome structures place it in the eisosome coat (section 6.8).

### 4.2 Weights across compartments

**Measured shares (CYCLoPs).** For 2,930 proteins, the share of each compartment comes from imaging: the CYCLoPs LOC-scores of the three wild-type screens of the yeast GFP collection (Koh *et al.* 2015 G3 5:1223; fraction of cells showing the protein in each of 16 compartments, averaged over the screens).
- **Groups:** the scores set the share of each compartment group (nucleocytoplasm, ER with nuclear envelope, mitochondria, Golgi/endosomes, vacuole, peroxisomes, plasma membrane). Within a group the annotation-based rule below decides, because imaging is less reliable there: the classifier calls many abundant cytosolic enzymes "Nucleus" (ERG10, TKL1), and LOC-scores are fractions of cells, not concentrations.
- **Periphery:** "cell periphery" and bud calls count as plasma membrane for membrane proteins and proteins annotated there, and as cytosol otherwise (abundant cytosolic proteins look peripheral because the vacuole presses the cytoplasm against the cortex). Cortical patches count as plasma membrane.
- **Not used** when the scores and the annotation share no group (534 proteins; C-terminal GFP mislocalises some secreted, lumenal and imported proteins, e.g. Hsp60, cytochrome c), for membrane proteins without a membrane call, for secreted cell-wall proteins and for ribosomal proteins (GFP-tagged r-proteins accumulate in the nucleus unassembled).
- **New groups:** a group the annotation does not list is added only from a LOC-score of 0.5, never the nucleocytoplasm (where fusions that fail to import end up), and mitochondria only for proteins of the high-confidence mitochondrial proteome.
- **Nucleus and cytosol (2026-11-11).** Within the nucleocytoplasm group the annotation rule used to split the two by volume (equal concentration), so proteins that imaging shows enriched in the nucleus had 92% of their copies in the cytosol. The split now follows the CYCLoPs Nucleus − Cytoplasm LOC-score difference *d*: the nuclear/cytosolic concentration ratio is exp(5*d*) (12× at *d* = 0.5, 1/12 at *d* = −0.5). The factor 5 was fitted on the dev half of an independent comparison (section 9.1) and raised agreement on the held-out half with N-terminal SWAT calls from 64.0% to 66.8% (paired difference +2.9 points, 95% CI 2.1–3.7; 49 proteins fixed, none broken).
- **Effect:** 1.25 M protein copies change compartment (counted in the build before the proteome was rescaled to 5.0 pg, when it held 42 M copies; not recounted since). Examples: proteasome RPN1 and WTM1 are nuclear, COPI coat proteins Golgi, actin-patch proteins (ABP1, SAC6) cortical, profilin cytosolic.
- **Complexes inherit the split (2026-11-11):** the curated 26S proteasome took a fixed 50/50 nucleus/cytosol split; it now takes the copy-weighted split of its subunits' measured weights (81% of RPN1-containing copies placed in the nucleus). The proteasome is assembled in, and stays mostly in, the nucleus of growing cells (Laporte *et al.* 2008 J Cell Biol 181:737).

**Annotation-based rule** (proteins without usable CYCLoPs data, and the split within groups):
- **Rank order:** a protein's compartments get weights 1, ½, ¼, … in annotation order, with the first-listed location taken as primary.
- **Cytosol and nucleoplasm:** their shares are pooled and re-split by volume, i.e. equal concentration across the nuclear pore. Since build S (2026-11-13) the volumes are those of the built cell (nucleoplasm 1.65 and cytosol 18.44 µm³, against 1.56 and 18.42 as built; the rest is pore channels and the smoothing of the nuclear surface). Until build R the inventory used 1.90 µm³ of nucleoplasm, which includes the lumen of the nuclear envelope, so "equal concentration" came out 1.24 times enriched in the nucleus. The fitted nucleus rule below uses the same volumes.
- **Secondary mitochondrial locations:** a mitochondrial location listed after a cytosolic or nuclear primary one (e.g. SOD1, ADK1 in the intermembrane space) gets 0.1× weight.
- **High-confidence mitochondrial proteins:** for proteins in the high-confidence mitochondrial proteome whose first location is mitochondrial, any other location gets 0.1× weight. This is the reverse of the rule for cytosolic proteins with a minor mitochondrial pool. Example: PET9, the major ADP/ATP carrier, is 5% on peroxisomes rather than a third.

**Mitochondria.**
- **Sub-compartments, measured:** for 810 proteins the sub-compartment (outer membrane, inner membrane, intermembrane space, matrix) comes from SILAC-based submitochondrial fractionation (Vögtle *et al.* 2017 Nat Commun 8:290); it replaces annotated or generic locations for 252 inventory proteins. For proteins Vögtle does not cover, the submitochondrial classes of Morgenstern *et al.* 2017 (Cell Rep 19:2836; Matrix/IM, IMS/IM, OM) replace a generic "mitochondrion" (soluble proteins to the matrix or intermembrane space, membrane proteins to the inner membrane).
- **False mitochondrial assignments removed:** proteins outside Morgenstern's high-confidence mitochondrial proteome (901 proteins) that are depleted from pure mitochondria (their Class 3; 28 proteins, e.g. Glt1, Lys21, COPI coat) lose their mitochondrial share.
- **Remaining generic locations** (now 0.02 M copies, was 0.52 M): transmembrane proteins split 75/25 between the inner and outer membranes, soluble proteins 90/10 between matrix and intermembrane space.
- **Copy numbers:** Morgenstern's absolute copy numbers in glucose correlate with PaxDb (Spearman 0.86 for 759 mitochondrial proteins) and give the same mitochondrial share of protein mass (27% vs 25% over the 4,137 shared proteins); their total is 1.5× higher (63 M molecules), a different proteome scale. PaxDb is kept, so that the whole proteome stays on one scale.

### 4.3 Membrane proteins

- **Definition:** proteins with UniProt TRANSMEM segments. β-barrels (TOM40, SAM50, MDM10, MDM34) have no TRANSMEM annotation and are added explicitly.
- **Assigned to soluble compartments only:** these take their membrane from GO, reading "cell periphery" as plasma membrane. Copies still without a membrane are dropped (11,200).
- **ER membrane proteins:** 15% go to the nuclear envelope, in proportion to area, because the outer nuclear membrane is continuous with the ER.
- **Lipid-droplet transmembrane proteins:** placed in the ER membrane.

### 4.4 Peripheral proteins

Soluble proteins annotated to a membrane are placed in thin "peripheral" shells next to it:
- **Plasma membrane, nuclear envelope and lipid-droplet surface:** on the cytosolic side.
- **ER (2026-11-12):** a soluble ER-annotated protein goes into the lumen only with a signal peptide or a lumen annotation (Kar2, Pdi1, Ero1 ...); the rest sit in a 25 nm shell on the cytosolic face (COPII coat and Sar1, Cdc48, Ypt1, Get3, SRP subunits, Scp160, Acc1). Before this fix 181,000 copies of 72 such proteins were packed inside the ER lumen, the same error the rule below had removed from the other organelles.
- **Inner-membrane-associated mitochondrial proteins:** in the matrix.
- **Vacuole, Golgi/endosomes and peroxisomes (2026-11-11):** a soluble protein annotated to one of these organelles goes into its lumen only if it can get there:
  - vacuole and Golgi: proteins with a signal peptide (UniProt), annotated to the lumen, or Cvt cargo (Ape1, Ams1, Ape4, which reach the vacuole by a selective-autophagy route; Lynch-Day & Klionsky 2010; Yuga *et al.* 2011);
  - peroxisomes: the matrix stays the default (matrix proteins are imported without signal peptides); only proteins UniProt places on the cytoplasmic side (Pex1, Pex6, Pex4) move out;
  - all others go into a 20 nm shell on the cytosolic face of the organelle's membrane: free V-ATPase V1 subunits, the EGO complex, Tco89 and part of TORC1, Vac8, COPI, clathrin, COPII, Arf1/2 and Rab GTPases. Before this rule they were packed inside the organelles (680,000 copies; the Golgi lumen held mostly coat proteins).
  - Complexes follow the majority of their subunits.
- **No cell-wall share without a signal peptide (2026-11-11):** soluble proteins with a cell-wall annotation but no signal peptide (e.g. Ssa1/2 from surface proteomics, Fus3, Kss1) lose the wall share (132,000 Ssa1/2 copies).
- **GPI-anchored proteins (2026-11-12):** proteins with a GPI anchor (UniProt lipidation; Gas1, Ecm33, Yps1, Plb1 and others) are soluble by the transmembrane definition and annotated to the cell membrane, so the peripheral rule had put 64,000 copies on the cytosolic face of the plasma membrane. Their plasma-membrane and cytosol shares now go to the wall region; the model has no separate layer between the bilayer and the wall.
- **"Spindle pole" (2026-11-12):** the UniProt term used to send any protein to the nuclear envelope (43,000 calmodulin were spread over the whole envelope). It now does so only for the 16 components of the spindle pole body itself, which section 6.2 then draws; for all other proteins the term is ignored and their other annotations decide.
- **Capacity of the shells (build S, 2026-11-13; \`wc3d.shell_capacity\`):** a membrane-face shell is asked to hold at most the sphere volume fraction that the cytosol is asked to hold (0.41 in build S); beyond that, every species of the shell is reduced by the same factor and the rest of its copies go to the cytosol as the protein's soluble pool. In build S this affects two shells: the Golgi/endosome face asked for 0.67, so 0.61 of every species stays and 141,000 copies go to the cytosol; the nuclear-envelope face asked for 0.51, so 0.80 stays and 8,000 copies go to the cytosol. A peripheral protein exchanges between its membrane and a soluble pool, and the annotation gives no bound fraction; much of the Arf1 and coatomer pool is probably cytosolic at any moment, which the imaging-based Golgi weight cannot separate.
- **Result (build S):** the vacuole face holds 116,200 copies (99.2% of its request), the Golgi face 211,600 (97.7% of the capped request of 216,700). Build R, without the cap: the Golgi face was asked for 358,000, of which 84% fit in a 20 nm shell around 8.6 µm² of membrane, 17% of the placed copies ended outside the shell, and the vacuole face held 117,000 (99.9%).
- **In the viewer (fixed 2026-11-13):** the default region mask had 24 bits, so the vacuole, Golgi/endosome and peroxisome faces (regions 24–26, added later) were packed but not drawn, about 300,000 proteins in build R; this, not the packing, was the empty-looking layer around Golgi and endosome membranes in the published viewer of build R.

## 5. Structures

### 5.1 Bead representation and mass

- **Beads:** the Cα for amino acids, and P (or C4′) for nucleotides.
- **Mass:** for densities, a residue counts as 110 Da and a nucleotide as 330 Da.
- **What is missing from the beads:** residues missing from a structure (unresolved in the PDB entry, or trimmed from an AlphaFold model) are absent from the model and from its mass.

### 5.2 Choice of structure (in order of preference)

1. **Curated PDB complexes:**

   | Complex | PDB |
   |---|---|
   | 80S ribosome (with Asc1 and Stm1, section 6.5) | 4V88 |
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

2. **Pre-ribosomes** (ribosome assembly intermediates; section 6.6), one yeast cryo-EM structure per stage, with all their RNA:

   | Stage | PDB | Compartment | Particles |
   |---|---|---|---|
   | 90S pre-ribosome (SSU processome) | 5WLC | nucleolus | 2,699 |
   | Nucleolar pre-60S (state E1) | 7R7A | nucleolus | 6,176 |
   | Nucleoplasmic pre-60S (Rix1–Rea1) | 6YLH | nucleoplasm | 1,886 |
   | Cytoplasmic pre-60S (Nmd3, Lsg1, Tif6) | 5T62 | cytosol | 6,126 |
   | Pre-40S | 6FAI | cytosol (85%), nucleoplasm (15%) | 2,654 |

   - **Copy number:** set by the stage's own assembly factors, as for complexes (below): the median over factors found only in that stage's structure, capped by all of its factors. The particles consume their factors before the Complex Portal step, so subcomplexes that are mostly particle-bound (UTP-B, PeBoW, Noc2–Noc3, Nog1–Tif6) are no longer placed free.
   - **Ribosomal proteins** in the particles are not taken from the pools; the mature 80S count drops by the particles' r-protein content instead (8,292 ribosome equivalents), so r-protein mass is conserved.
   - **Total:** 19,500 particles, 8.9% of ribosomes. At a production rate of about 2,000 ribosomes per minute for an average cell (Warner 1999), scaled to this cell, this is 2–4 minutes per stage.
   - **Placed (build S):** the nucleolus cannot hold its request, and its request is thinned by 0.907 before packing (section 8.2), so 4,801 of the 6,176 nucleolar pre-60S and 2,160 of the 2,699 90S particles are placed (329 of those on nascent transcripts); the other stages fit 99–100%. Build R, without thinning: 5,079 and 2,565 (349 on transcripts); the other stages 98–99%.
3. **Complex Portal complexes** (118 added) with a *S. cerevisiae* RCSB assembly. Rules:
   - **Acceptance:** the assembly must cover at least 80% of the subunit mass, and at least 60% of the residues must be modelled.
   - **RNA chains** are kept only where they belong to the complex: snoRNAs, snRNAs, telomerase RNA, and mitoribosomal rRNA.
   - **Copy number:** the median subunit pool, capped at 2× the limiting subunit. Complexes with four or more subunit types are capped at 2× the 25th-percentile subunit.
   - **Mitochondrially encoded subunits** (COB, COX1–3, ATP6, ATP8, OLI1, VAR1) are ignored when setting copy numbers, because mass spectrometry badly under-detects them.
   - **Subunit bookkeeping:** subunits are consumed from the monomer pool, so total protein copies are conserved.
   - **Location (2026-11-11):** the mean of the subunits' compartment weights, each subunit weighted by 1 / (number of Complex Portal complexes it belongs to): a subunit shared between complexes (Tor2 and Lst8 in TORC1 and TORC2, the Lsm ring of U6 snRNP and the cytoplasmic Lsm1-7) reports a mixture of their locations. Five of 118 complexes moved by more than 10% (U5 and U6 snRNPs more nuclear, CORVET towards endosomes, the G protein and Dnf1–Lem3 towards the plasma membrane).
4. **Homo-oligomers:**
   - **PDB:** 134 homo-oligomeric assemblies (112 found from the UniProt subunit text, 22 curated earlier).
   - **Curated from a PDB audit (2026-11-12; \`assemblies.CURATED_HOMO\`):** abundant proteins whose UniProt text does not name the oligomer and that were therefore monomers: Cdc48 hexamer (8U9P, substrate peptide left out), Tsa1 decamer (3SBC; a C47S mutant; the decamer exchanges with dimers), Gln1 decamer (3FKY; first 18 residues missing), Mmf1 and Hmf1 trimers (3QUW, 1JD1), Pst2 tetramer (5MP4), Ynk1 hexamer (3B54, an AlphaFold dimer before). Sah1 and Gdh1, also oligomers, have no *S. cerevisiae* entry in the PDB and stay monomers.
   - **AlphaFold:** 443 AlphaFold-predicted homodimers with ipTM ≥ 0.5.
5. **Heterodimers:** 476 AlphaFold-predicted heterodimers, from 1,250 models (build R: 475). Requirements:
   - ipTM ≥ 0.5.
   - The two partners share a compartment. A membrane protein may pair with a soluble partner on its cytosolic or matrix side.
   - Not both partners are subunits of an already-modelled complex.
   - Neither partner forms a larger homo-oligomer.
   - A pair with a histone is a native histone pair (H2A–H2B, H3–H4, H2A.Z–H2B, Cse4–H4): histone folds pair promiscuously in AlphaFold (H4–H2B, H2A.Z–TAF9 and others were accepted before 2026-11-12).

   **Copy allocation:**
   - Each protein offers its free copies to its accepted partners in proportion to the partners' copies, and a pair forms the smaller of the two offers (2026-11-12). Before, copies were split evenly over the partners: the rare paralog Btt1 took half of Egd2, and the obligate NAC heterodimer Egd1–Egd2 formed 61,000 dimers beside 123,000 free subunits (in the 5 pg cell). Now it forms 85,500 dimers with 1,550 free subunits.
   - Pairs supported by Complex Portal or a direct IntAct interaction bind all available copies; the others bind 50%. This is an assumption: without affinity data, stable and transient interactions cannot be told apart.
6. **AlphaFold monomers** for everything else. Residues with pLDDT < 50 are trimmed, which removes 3.4% of placed protein mass overall. The loss is concentrated in nuclear and nucleolar proteins: 15–16% of mass there, measured on an earlier build.

**Result by placed mass (build S):** PDB 52.0%, AlphaFold complexes 15.4%, AlphaFold monomers 30.1%, tRNA models predicted with Boltz-2 2.6% (build R: 52.4, 15.0, 30.1 and 2.6%).

**A systematic audit of the PDB (2026-11-12)** screened all 7,067 *S. cerevisiae* entries against the library. Of 869 distinct multi-protein subunit sets, 496 have at least 80% of their genes inside one PDB type of the model. What it found and fixed is in this section and in 6.5; the open candidates (a translating ribosome with tRNAs, eEF3, Ssb and eIF5A on the ribosome, 43S/48S complexes, remodelers on nucleosomes, eEF1A–eEF1Bα and others) are listed in the project's \`pdb_audit_2026-11-12/README.md\`.

### 5.3 Membrane orientation

- **Method:** UniProt TRANSMEM residues are mapped onto the structure, through RCSB's UniProt alignment for PDB chains. A search finds the membrane normal and offset that put the most TM residues, and the fewest other residues, inside a ±1.6 nm slab.
- **Side:** from UniProt topology domains: the N side (cytoplasm; matrix for the inner mitochondrial membrane; nucleoplasm for the nuclear envelope) faces +z, the P side (lumen, intermembrane space, vacuole interior, outside) −z. Without topology domains, the side carrying more soluble mass is taken as cytoplasmic.
  - **Fix (2026-11):** single-chain AlphaFold models had their side set only from "Cytoplasmic" domains, so 66 proteins annotated with other sides (mostly inner-membrane proteins: "Mitochondrial matrix" / "Mitochondrial intermembrane") were oriented by chance; 39 of them, e.g. TIM50, TIM21 and Mgm1, had their intermembrane-space domains in the matrix. They now use the same side rule as everything else.
- **β-barrels:** oriented by their β-strands, found from Cα geometry.
- **Library entries older than the side rule (found 2026-11-13):** when stage 0 rebuilt the structure library from the pipeline (section 13), 82 transmembrane proteins came out rotated by 180° against the library on disk, whose entries predated the side rule above; 61 of them were packed in build R (42,800 copies; Kre2, Ssh1, Ato3, Erv29, Ctr1 and others) upside down. Build S uses the rebuilt entries.
- **Check:** the method reproduces hand-oriented Pma1 and V-ATPase (z-coordinate correlation 1.000 and 0.9996).
- **Result:** in the packed cell, the median distance of TM residues from the membrane mid-surface is 0.72 nm, and 97.4% lie within 2 nm (build S; build R 0.73 nm and 97.3%).

## 6. Supramolecular organisation built explicitly

These structures are built before random packing. They then act as fixed obstacles for everything packed afterwards.

### 6.1 Nuclear pores and FG mesh

- **Pores:** 93 copies of 7N9F (the in-situ yeast pore), with the pore axis along the local nuclear-envelope normal.
- **None in the nucleus–vacuole junction (build S, 2026-11-13):** pores are absent from the envelope where it is closely apposed to the vacuole (Severs *et al.* 1976 and Pan *et al.* 2000, as cited by Kohler & Büttner 2021 and Otto & Thumm 2020). No pore centre is placed where the vacuole is closer than 92 nm (half the perinuclear space, the 45 nm junction contact distance of section 6.10 and a 30 nm margin; a pore's membrane torus has a 45 nm radius). Up to build R pores were spread over the whole envelope, including the junction.
- **Nucleoporin bookkeeping (2026-11-12):** each pore takes the measured copies of every nucleoporin: 8 to 48 per pore, about 550 proteins and 52 MDa in all (Kim *et al.* 2018 Nature 555:475, Supplementary Table 2A). Until then the subunit list came from Complex Portal, which gives no stoichiometry, so one copy of each nucleoporin per pore was booked and about 100,000 nucleoporin copies were placed as free proteins, mostly at the nuclear envelope.
  - **Free pool:** PaxDb still holds more scaffold and FG nucleoporins than the pores need (median 1.36 times at 3.5 pg). These exist almost only in pores, so their free pool is capped at 25% of the copies in pores (an assumption) and the surplus, 13,000 copies, is not placed. Sec13, Seh1, Dyn2 and the basket proteins Mlp1/2, which have other homes, are not capped.
- **FG mesh:** the disordered FG-repeat domains of the FG-Nups modelled in 7N9F are threaded as self-avoiding residue chains (0.38 nm per residue, low persistence) from the first modelled residue of each chain:

  | Nup | Residues | Copies per pore |
  |---|---|---|
  | NSP1 | 1–636 | 48 |
  | NUP57 | 1–286 | 32 |
  | NUP49 | 1–269 | 32 |
  | NUP159 | 388–1102 | 16 |

  - **Confinement:** the chains stay within 60 nm of their anchor, in the pore channel, cytoplasm or nucleoplasm.
  - **Result:** 4.9 of 5.5 million residues placed (89%).
  - **Missing:** FG-Nups absent from 7N9F (NUP100, NUP116, NUP145N, NUP42, NUP1, NUP60, NUP2), so the nuclear basket side is sparse.
  - **Caveat:** the chain conformations are excluded-volume random walks, not a physical model of FG cohesion.
- **Transport receptors in the channel (2026-11-12; \`pipeline/pore_receptors.py\`):** 20 receptors per pore (the Jurkat cell model's number; unchecked for yeast), drawn from the free pools of Kap95, Srp1, Kap123, Pse1, the Crm1 export complex, Mex67–Mtr2 and Ntf2 in proportion to their copies, are placed in the channel (within 22 nm of the pore axis and 25 nm of the envelope midplane) before the FG domains are threaded around them. Before, the packer never put a receptor inside a pore. All 1,860 fit (builds R and S).

### 6.2 Microtubules and the spindle pole body

- **Number:** 20 microtubules, made of 12,324 αβ-tubulin dimers (53% of the tubulin pool; build R: 12,831, 55%):
  - 15 kinetochore microtubules, 262 nm mean length (16 are asked for; the 16th did not fit at the plaque in builds R and S; its kinetochore sits 40 nm from the plaque in build R, not re-measured);
  - 2 other nuclear microtubules (4 asked for, about 400 nm): mean length 208 nm in build S (357 nm in build R);
  - 3 cytoplasmic microtubules, 1.2 µm (mean 1,173 nm).
- **Nucleation sites (build S, 2026-11-13):** up to 400 minus-end sites on the 90 nm plaque are tried per microtubule (80 before), since the rings nearly fill the plaque. With the nucleolus moved to the vacuole side (section 3.3) the spindle pole body faces open cytoplasm again, and the 3 cytoplasmic microtubules are built (round S candidate builds 1 and 2, with the spindle pole body against the vacuole, had none).
- **Lattice:** an ideal 13-protofilament, 3-start lattice with a seam. The lateral geometry comes from 8QV3 (native yeast γ-TuRC-capped microtubule, cryo-ET). The dimer is 5W3F, with an 8.3 nm axial repeat.
- **Caps:** each minus end has a γ-tubulin ring of 7 γ-TuSCs from 8QV3.
- **Growth:** from 90 nm plaques at the spindle pole body. Microtubules bend away from membranes, with a minimum radius of curvature of 300 nm.
- **Assumptions:** the lengths (G1-like), and the 70% cap on the use of the tubulin pool.
- **Spindle pole body plaques (2026-11-12; \`pipeline/spb.py\`):** the SPB is a layered cylinder in the envelope (cryo-EM of isolated SPBs, Bullitt *et al.* 1997 Cell 89:1077; integrative model of the core, Viswanath *et al.* 2017 Mol Biol Cell 28:3298). Its proteins had been spread through the nucleoplasm, so a simulated Spc42-GFP image showed no spot. Now each SPB protein goes into its layer, a disc of the plaque radius (90 nm) on the SPB axis: outer plaque (Spc72, Nud1; 45–75 nm above the outer membrane), Cnm67 (20–50 nm), central plaque (Spc42, Spc29, calmodulin; 3–25 nm), inner plaque (Spc110; 3–60 nm below the inner membrane), half-bridge (Sfi1, Kar1; beside the plaques).
  - **Core stoichiometry (2026-11-12):** Spc42 : Spc29 : Spc110 : Cmd1 : Cnm67 = 2 : 1 : 1 : 1 : 1 (Muller *et al.* 2005 Mol Biol Cell 16:3341, as used by Viswanath *et al.* 2017). The inner plaque is filled first; the Spc110 at the SPB is then that layer plus the 14 per nuclear γ-TuRC, calmodulin takes one per Spc110, and Spc42, Spc29 and Cnm67 are added in units of 2 : 1 : 1 until units stop fitting. Before, each protein was filled from its whole pool, which gave Spc42 : Spc29 : Cnm67 about 1 : 1 : 1 and calmodulin only on the inner-plaque Spc110.
  - **Result (build S):** Spc110 270 (238 in rings, 32 in the inner plaque), calmodulin 270, Spc42 297, Spc29 152, Cnm67 152, Nud1 536, Spc72 184, Sfi1 93, Kar1 38 (build R: Spc110 274, calmodulin 274, Spc42 239, Spc29 124, Cnm67 124, Nud1 413, Spc72 111, Sfi1 74, Kar1 38). The 90 nm disc runs out of room at about 1.1 Spc42 per Spc110 instead of 2 (build R: 0.9), so the core is short of Spc42, Spc29 and Cnm67 relative to Spc110. The absolute number per SPB is not measured here. Layer heights are approximate (unchecked). Simplification: the envelope has no fenestra, so the central plaque sits on its cytoplasmic face. The inventory copies exceed one SPB (573 Spc42), so what does not fit stays free; the simulated Spc42-GFP image shows one spot.

### 6.3 Chromatin

- **Rabl configuration:** each chromosome is two arms, both starting at its centromere, within 20 nm of its kinetochore-microtubule tip at the spindle pole body, and running to a telomere at the nuclear envelope. Arm lengths come from the S288C centromere positions (32 arms; chromosome XII-R is split where the rDNA array sits, so it runs into the nucleolus and out again, and the rDNA itself is in section 6.6).
- **Fibre (2026-11-11):** each arm is a worm-like backbone with the mass density and persistence length measured for interphase chromatin in living yeast (Bystricky *et al.* 2004 PNAS 101:16495: 110–150 bp/nm and 170–220 nm; the middle values 130 bp/nm and 195 nm are used): the backbone advances 1.27 nm per 165 bp nucleosome repeat. Nucleosomes sit at random within 60 nm of the backbone, each within reach (18 nm centre to centre, more across long linkers) of the previous one: dispersed, irregular clutches, not a regular 30 nm fibre, as in-situ cryo-ET of yeast nuclei suggests. The backbone drifts towards the arm's target as its remaining contour runs short, so arms wander freely at first and end at their telomeres; where nucleosomes repeatedly fail to fit, the backbone turns, and as a last resort the chain restarts nearby (a break).
  - **Why:** the earlier walk placed one nucleosome per 12.5 nm step (about 13 bp/nm), so short stretches were far too extended (rms 214 nm at 20 kb against 118–160 nm in vivo) while long stretches were compressed by the nucleus.
  - **Fibre radius (2026-11-11, late):** a first version with a 20 nm radius gave 122–146 nm at 14–20 kb (Bystricky's range) but failed the Hi-C test: 62% of contacts within chromosomes (Hi-C 37–41%), P(s) slope −1.6 (−1.09 to −1.27) and centromere clustering 13× (4–5×). Nine variants were built as pre-passes only (compaction 60–130 bp/nm, persistence 70–195 nm, radius 20–120 nm, and the earlier walk), and the 60 nm radius was chosen: it is the only variant whose R(s) lies inside the chain inferred by Arbona *et al.* 2017 (Genome Biol 18:81; 53–65 bp/nm, persistence 52–85 nm, from imaging and Hi-C together) at both 20 and 100 kb (197 and 411 nm against 163–223 and 393–553 nm), and inside Bystricky's at 100 kb (above it at 20 kb). The Hi-C test is therefore partly tuned: it was looked at in this choice.
  - **Telomere steering (fixed 2026-11-11, late):** the fibre walk turned the backbone in a fully random direction after every failed placement, which erased the drift towards the telomere focus, and the foci of short arms lay beyond their own fibre contour. Arms ended a median 790 nm from their targets (the earlier walk: 37 nm). Now a failed placement turns the backbone with a pull towards the target that grows as the contour runs short, the drift weight is capped at 5 (was 1.5), and each telomere target is moved towards the pole until it lies within 80% of its arm's contour of the start (\`TEL_REACH\`). In a pre-pass the median miss is 160 nm; the last nucleosome lies up to 60 nm off the backbone. The variant scan above was repeated after the fix (cis 0.47, slope −0.96 for the chosen radius) and the choice held.
- **Measured nucleosome positions (2026-11-11; \`pipeline/nucleosome_map.py\`):** the fibre is built from a measured nucleosome map: one nucleosome per measured dyad, in genomic order, so every nucleosome carries its genome position (\`nuc_bp\` in the pre-pass).
  - **Map:** dyads called from the H3 chemical-cleavage dyad maps (H3Q85C; W303 background, YPD) of Chereji *et al.* 2018 (Genome Biol 19:19; three replicates summed; sacCer3, the coordinates of the model's chromosomes). The signal is smoothed (Gaussian, 15 bp), and peaks are taken from the highest down, each excluding 135 bp on both sides, while the signal is above half its genome-wide median; below that the DNA counts as nucleosome-free. That gives 66,245 dyads. Dyads inside promoter nucleosome-depleted regions are then removed: the regions between the upstream stable nucleosome and the +1 nucleosome in Rossi *et al.* 2021 (Nature 592:309, Supplementary Data 1; 10,313 regions; dyads more than 73 bp inside both boundaries). That removes 1,758, leaving 64,487 nucleosomes on the arms (spacing median 169 bp, 4.3% of gaps over 300 bp). The chemical map of Brogaard *et al.* 2012 (Nature 486:496; 67,548 unique nucleosomes) is on the older sacCer2 coordinates (offsets of a few bp that grow along the arms), and Jiang & Pugh 2009 on an intermediate pre-sacCer3 release, so neither is used.
  - **Checks:** a nucleosome-free region 200–100 bp upstream of start codons (0.4 times the flanking density) and a phased array downstream. Confirmed origins sit in wider gaps (median 236 bp between the flanking dyads, against 174 bp at random positions), as expected from origin nucleosome positioning (Eaton *et al.* 2010 Genes Dev 24:748).
  - **Geometry:** the backbone advances by each nucleosome's spacing over 130 bp/nm. The reach from the previous nucleosome grows with the extra linker, 18 nm plus 0.34 nm per bp beyond 165 bp, at most 60 nm (across nucleosome-free regions).
  - **Count:** the inventory now takes the map's count (64,487 plus the rDNA) instead of the genome length over a 165 bp repeat (73,047). That is still 7–13% more than the 57,000–60,000 per cell measured by absolute occupancy (Oberbeckmann *et al.* 2019 Genome Res 29:1996; held out, not tuned to): a population map holds positions that a single cell does not all occupy.
  - **Held-out tests (\`validation/nucleosome_test.py\`):** the absolute occupancy at 13,807 restriction sites (Oberbeckmann 2019) is 0.75 at sites inside a placed model nucleosome and 0.57 outside (AUROC 0.75; 79% of sites covered against a mean occupancy of 0.72; build M). Of the MNase nucleosomes of Weiner *et al.* 2015 (unstressed), 79% have a model dyad within 40 bp. ATAC-accessible regions (Schep *et al.* 2015) hold 10% fewer model dyads per kb than the genome average; the merged regions are 1.3 kb wide, so the test is coarse.
  - **DNase I at base resolution (2026-11-12; \`validation/dnase_test.py\`; held out):** 83 million DNase I cuts in log-phase cells (Zhong *et al.* 2016 Genome Res 26:351; W303 in YPD, two replicates, lifted from sacCer2 to sacCer3). Around the model's dyads the cut density dips to 0.62 of the genome mean at the dyad and rises to 1.04 in the linkers 75–95 bp away; the core (±50 bp) is cut 0.66 times as often as the linkers, and 76% of nucleosomes are cut less in their core than in their own linkers. In 50 bp bins, stretches without a model nucleosome are cut 2.3 times the mean and fully covered ones 0.71 times. DNase-hypersensitive sites of BY4741 (Zelin *et al.* 2012 Mol Cell 48:459; 3,260 peaks) have a model occupancy of 0.62 against 0.76 genome-wide.
  - **Fit test:** in test pre-passes all 64,487 nucleosomes are placed (20 and 60 nm radius).
- **Telomeres:** gathered in 7 foci at the nuclear envelope (6–8 seen by immunofluorescence; Gotta *et al.* 1996), 25 nm inside the inner nuclear membrane. A focus's angle from the spindle-pole axis grows with the arm length of its telomeres, 35° for the shortest to 115° for arms of 700 kb or more. That angle range is an assumption; the trend is from Hi-C and imaging (Duan *et al.* 2010; Therizols *et al.* 2010). The foci stay out of the nucleolus.
- **Result (build S, in the nucleus of 2.0 µm³, with the nucleolus and the spindle pole body turned to the other side; section 3.3):** 62,882 of 64,487 nucleosomes placed (the rest evicted by Pol II and origins). All 32 telomeres are within 150 nm of the nuclear envelope (median 94 nm; build R 112 nm), a median 164 nm from their targets. Telomere-to-pole distance follows arm length (Spearman 0.92; build R 0.85). Against the contact maps (\`validation/hic_test.py\`), as provided and ICE-balanced: 44% / 50% of contacts within chromosomes (Hi-C 37–41% / 39–51%), P(s) slope −0.94 / −1.01 over 10–200 kb (−1.09 to −1.27 / −1.07 to −1.35), centromere clustering 6.5× / 5.9× (4.0–5.3× / 4.0–7.2×), telomere clustering 3.3× balanced (1.7–3.6×). Balanced, the cis share and both clustering measures are inside the range of the two maps; the contact decay is flatter than both. Build R had 51% cis, −1.04, 6.0× and 1.9× (balanced), build Q, in the 2.9 µm³ nucleus, 55%, −1.03, 7.1× and 2.3×. One chromatin realisation moves these measures by about this much (builds O and P differed as much with identical rules), so differences between builds R and S are not evidence for or against the new layout. Rms distances along the chromosome: 191, 289, 398 and 516 nm at 20, 50, 100 and 200 kb (build R: 195, 284, 376, 467; build Q: 193, 294, 403, 536). At 20, 50 and 100 kb that is inside the chain inferred by Arbona *et al.* 2017 (393–553 nm at 100 kb, inside which the fibre radius had been chosen; build R was 4% below it), at 200 kb below it (561–792 nm). Which chromosome pairs touch is not reproduced (balanced Spearman 0.14 with both maps), but the two maps also disagree with each other on that (Spearman −0.04 / 0.15).
- **Linker DNA:** joins consecutive nucleosomes, with a schematic path.
- **In the viewer (2026-11-12; \`pipeline/polymer_annot.py\`):** hovering or clicking chromatin gives its genome position (R64 coordinates):
  - **Nucleosomes:** the chromosome and the measured dyad, the gene or element there (SGD features: ORFs, RNA genes, centromeres, ARS, telomeres, Ty elements; otherwise the neighbouring genes and the distance upstream of a gene start), and the enriched histone marks.
  - **Linker DNA:** the bp range between the two wrapped cores.
  - **Pol II and the DNA through it:** the midpoint of its two flanking nucleosomes, and the gene there. Near gene ends that can be the neighbour of the transcribed gene. Nascent transcripts carry the same position.
  - **Without positions:** rDNA (nucleosomal or active) and mtDNA are named by kind only.
- **Transcribing Pol II:** per-gene numbers of engaged Pol II from genomic run-on (Pelechano *et al.* 2010 PLoS One 5:e15442, Table S1; 691 per cell over 4,670 genes, about 2% of all Pol II), at uniform positions along each gene (SGD coordinates).
  - **On the arms:** with the measured map, a polymerase replaces the nucleosome nearest its position in bp (earlier builds mapped kb linearly along each arm walk). The nucleosome at a polymerase's position and its two neighbours are removed (nucleosomes are evicted around elongating Pol II), and the polymerase takes the middle one's place. Its DNA axis comes from the elongation complex 1Y1W, superposed on the library 1WCM to 0.01 nm, and lies along the chain. DNA runs from the upstream nucleosome through the cleft to the downstream one.
  - **Nascent RNA:** threaded from the RNA exit, as long as the polymerase has transcribed, 2 nt per bead.
  - **Source pool:** engaged polymerases come out of the free Pol II pool.
- **Licensed replication origins (2026-11-11; \`pipeline/origins.py\`):** the model cell is G1-like, so each confirmed origin of OriDB (Siow *et al.* 2012 Nucleic Acids Res 40:D682; 410 'Confirmed', sacCer3) gets one MCM2–7 double hexamer (7W8G) and one ORC (7TJH).
  - **Placement:** the double hexamer goes in the gap between the two measured nucleosomes that flank the origin, its long axis (the DNA channel) along the line joining them. ORC sits beside it.
  - **Where it does not fit:** the fibre packs nucleosomes too densely for a 21 nm double hexamer between most flanking pairs. It is then moved up to 25 nm out of the fibre (the origin's free DNA looping out). Where that fails too, the two flanking nucleosomes are evicted (last resort: two on each side).
  - **Result (build S):** 365 of 410 origins licensed, with 356 ORC and 90 nucleosomes evicted (build R: 347, 340, 86; a test pre-pass before the fibre's final settings: 317, 296, 236).
  - **Limits:** one double hexamer per origin is a lower bound, since efficient origins load several (unchecked). The rest of the MCM (1,383 double hexamers in build R; not re-counted) and ORC (506 in build R) pools stay free in the nucleoplasm.
- **Kinetochores (2026-11-11; \`pipeline/kinetochores.py\`):** one per chromosome on its kinetochore-microtubule plus end: a CCAN bound to centromeric DNA with Cbf1 (8OVW) just beyond the plus end along the microtubule axis, and 6 MIND (Mis12) complexes (9S5N) in a ring around the plus end. Quantitative imaging counts per kinetochore in metaphase two Cse4, one or two copies of the inner proteins, 6–7 MIND, 8 Ndc80 complexes and 16–20 DAM/DASH (Joglekar *et al.* 2006 Nat Cell Biol 8:581, Table 1; checked 2026-11-12 against the paper, so the model's 6 MIND, 8 Ndc80 and 16 DASH are within the measured counts). The chromosome's two arms start beyond the CCAN. All 16 CCAN and 96 MIND fit (builds R and S).
  - **Outer kinetochore (2026-11-12):** per kinetochore 16 copies of each DASH/Dam1 subunit (Joglekar's count) in a ring 16–32 nm from the microtubule axis over its last 5–60 nm, and 8 copies of each Ndc80-complex subunit (Ndc80, Nuf2, Spc24, Spc25; 8 in Table 1 of the same study) in a band from the lattice to just beyond the tip; 3,050 of 3,072 subunits fit (build S; build R 2,908). Both complexes are placed as their AlphaFold subunits (no assembled library types), so the ring is approximate. Not built: the Cse4 nucleosome.
- **Chromatin proteins at their sites (2026-11-11; \`pipeline/chromatin_sites.py\`):** ChExMix peaks of the ChIP-exo survey of Rossi *et al.* 2021 (about 400 chromatin proteins, YPD, sacCer3; filtered peak files of the paper's repository, 188,489 peaks).
  - **Which proteins:** each protein that is a library type of its own in the nucleoplasm (about 280; subunits of complexes the inventory does not assemble, such as TFIIH, Mediator or SAGA, are placed one by one, so they gather at the same promoters).
  - **How many:** min(peaks, 30% of its nuclear copies) copies go on sites, drawn without replacement with probability proportional to the peak's tag count. The 30% bound share is an assumption (unchecked); the rest stay free.
  - **Geometry:** the DNA point of a site is interpolated between the two measured nucleosomes that flank it. The protein (random orientation) comes within 5 nm of it with one of its spheres, on the side away from the local nucleosome mass, and overlaps nothing else.
  - **Result (build S):** 20,226 of 50,931 wanted copies of 268 proteins fit (40%; build R 19,238 of 51,073 copies of 270 proteins, 38%); the fibre leaves little room next to most linkers. The peak calling excluded the rDNA, tRNA genes and telomeres, so Pol III at tRNAs and Sir/Rap1 at telomeres are under-represented.
- **tRNA genes, telomeres and silent chromatin (2026-11-12; \`pipeline/trna_telomeres.py\`):** the ChIP-exo peak files exclude tRNA genes and telomeres, so these sites are placed with the same geometry (\`chromatin_sites.place\`, one copy per site): one Pol III and one copy of each TFIIIC subunit in the library (Tfc1, Tfc7) on each of the 275 nuclear tRNA genes (SGD); 17 Rap1 at each of the 32 chromosome ends (TG1–3 tracts of about 300 bp with a Rap1 site every 18 bp on average: Marcand *et al.* 1997 Science 275:986; 15 before 2026-11-12); per Sir-bound nucleosome, two Sir3 and one each of Sir2 and Sir4 (saturating stoichiometries on reconstituted arrays, Swygert *et al.* 2014 Nat Commun 5:4751 and 2018 PNAS 115:12447, so upper bounds in cells; one each before). The chromatin is crowded, so a share does not fit (build S: 811 of 2,305 copies: Pol III on 102 of 275 tRNA genes, Rap1 89, Sir2 123, Sir3 82, Sir4 53; build R: 842 of 2,317, with Rap1 90, Sir2 135, Sir3 95, Sir4 57).
  - **Which nucleosomes are Sir-bound (2026-11-12; \`pipeline/silent_chromatin.py\`):** those where Sir3 ChIP-seq is enriched at least 5-fold over input in both replicates (anti-Sir3 ChIP of cycling W303 cells in YPD; Brothers & Rine 2022 eLife 11:e75653), within 30 kb of a chromosome end or in HML or HMR: 234 placed nucleosomes at 28 of the 32 ends in build S (109 within 3 kb, 70 at 3–10 kb, 55 at 10–30 kb; build R 237). Before, every nucleosome within 3 kb of an end and in HML/HMR was taken (529).
    - **Threshold and region, set before the test:** the same laboratory's untagged and *sir3Δ* controls reach about 1 and 4–5 at silent loci, and highly expressed genes are "hyper-ChIPable" (Teytelman *et al.* 2013 PNAS 110:18602): 33 nucleosomes at TDH3 and other such genes pass the threshold and are excluded by the 30 kb limit.
    - **Caveat:** the strains express Sir3 fused to a methyltransferase, not untagged Sir3.
    - **Held-out test (\`validation/silent_chromatin_test.py\`):** 77% of the Sir-bound nucleosomes are H4K16-hypoacetylated (below the genome's 10th percentile; Weiner *et al.* 2015; Sir2 deacetylates H4K16), against 28% under the earlier rule and 14% of the nucleosomes only the earlier rule took. Sir4 ChIP-seq (Goodnight & Rine 2020 eLife 9:e58910; with an untagged control at 0.9) is enriched 36-fold on them, against 4-fold under the earlier rule (build S; build R 78% and 37-fold).
  - **Check (conflicting evidence):** by FISH, tRNA genes cluster at the nucleolus (Thompson *et al.* 2003 Science 302:1399), an arrangement the model has no rule for; the balanced contact maps show no excess of trans contacts between tRNA-gene bins (1.05× and 1.00×; the model 1.14× in build S, 0.94× in build R), so no rule was added.
  - **Test (2026-11-12, \`validation/trna_nucleolus.py\`):** FISH finds the active SUP53 tRNA gene overlapping the nucleolus in 52% of 440 cells and the same gene with an inactive promoter in 13% of 715 (Thompson *et al.* 2003); in wild type 76% of the nuclei with a single tRNA-gene cluster have it at the nucleolus (Haeusler *et al.* 2008 Genes Dev 22:2204). In the model (build S), 23% of the 275 tRNA genes lie within 250 nm of the nucleolus (a FISH spot's radius, set before comparing), against 33% of nucleosomes, and SUP53 is 520 nm away: tRNA genes are no closer to the nucleolus than chromatin in general, so the test fails. (Build R: 35% against 43%, SUP53 740 nm; build Q, in the larger nucleus: 33% against 36%.)
- **A relaxed polymer instead of the fibre walk (2026-11-12, tested, not adopted; \`pipeline/chromatin_polymer.py\`, \`validation/chromatin_polymer_scan.py\`):** after the Jurkat cell model, every arm was simulated as a coarse polymer (2 kb beads; Arbona's 60 bp/nm and 70 nm persistence; soft excluded volume; confined to the nucleoplasm and clear of the nucleolus, pores and microtubules; centromeres tethered at their kinetochores, telomeres within 60 nm of the envelope, no foci imposed) with Langevin dynamics in OpenMM, and the measured nucleosomes were placed along the relaxed chain (\`pack3d.guided_walks\`). Its R(s) lies in the Arbona envelope. Scored on ICE-balanced Hi-C, it improved the contact decay (P(s) slope −1.24 to −1.42 against −1.07 and −1.35; the fibre walk −1.03) but made the chromosomes more territorial (cis share 0.63–0.77 against 0.39–0.51; fibre 0.57) and the centromeres more clustered (8–20× against 4.0–7.2×; fibre 6.4×), in all variants tried (longer runs, weaker repulsion, Bystricky's 130 bp/nm and 195 nm, looser telomere tethers; 4–8 replicas each). The fibre walk stays.
- **Not modelled:** inter-chromosomal contacts beyond the Rabl layout, histone variants and marks as structure (H2A.Z, H3K4me3, H3K36me3, H3K79me3, H3K56ac, H3K9ac and H4K16ac of the nearest unstressed MNase nucleosome within 40 bp in Weiner *et al.* 2015 are attached to 83% of nucleosomes in \`work/chromatin/nucleosome_genome.parquet\`, as annotation only; readers are not placed by them), cohesin and condensin.

### 6.4 Actin

The ACT1 pool (144,600 copies at 3.5 pg; 206,500 in builds J to Q) is split by assumption: 30% in cables, 30% in patches, and the rest monomeric or in profilin, cofilin and Srv2 heterodimers.

**Cables:**
- **Number:** 11 cables, 81 filaments, 43,400 subunits and 120 µm of filament in total (build S; build R 82 filaments).
- **Position:** 90 nm beneath the plasma membrane, along great circles.
- **Bundle:** filaments sit on a hexagonal lattice 12 nm apart (assumed fimbrin bundle spacing).
- **Filament structure:** yeast F-actin 9GO5, 167° and 2.76 nm per subunit.
- **Tropomyosin:** 8,269 TPM1 coiled-coil dimers (build S; 8,206 in build R), one per 5 actins per strand. The path is transferred from the mammalian actin–tropomyosin structure 5JLF, because no yeast structure exists. The chains are wound as an ideal coiled coil.

**Endocytic patches:**
- **Number:** 29 patches of about 1,500 actin each, one at each endocytic invagination (section 3.2; the number of sites follows the plasma-membrane area). The network grows in the cytosol around the tube.
- **Amphiphysin on the tubes:** Rvs161–Rvs167 heterodimers (AlphaFold model) on the tubular part of every invagination deeper than 50 nm, between the neck and the coated tip region (the last 35 nm). Their number rises linearly to 100 at 100 nm depth: Rvs assembles from about 50 nm and peaks at about 100 Rvs167 at scission (Picco *et al.* 2015 eLife 4:e04535). Each dimer lies around the tube with its concave face on the membrane. 422 of the 695 called for fit (build S; build R 386 of 748); random addition jams before the dense lattice that would cover the tube. The library entry of the heterodimer is made by the builder itself since build S (\`amphiphysin.model\`); before, the only copy was a leftover of 2026-11-07 that no stage wrote, found when a build started from an empty library.
- **Endocytic coat:** Sla2, Sla1 and End3 on the tip of each invagination (the last 35 nm; at flat sites, the membrane patch under the actin), with numbers per site from Picco *et al.* 2015. Sla2 is about 40 throughout. Sla1 (about 90) and End3 (about 60) are present at the start and disassemble during invagination, here linearly to none at 100 nm depth (assumption). Orientations are random, since there is no structure of the coat. Each molecule sits as close to the membrane as it fits.
- **Not built:** clathrin, Pan1, the epsins Ent1/2. Picco *et al.* also count about 3,000 actin per site at peak, twice the 1,500 used here.
- **Branching:** dendritic networks branch at 73° from Arp2/3 complexes (2,998 branches in build S; 3,071 in build R). The branch-junction geometry is 8UXW. The yeast Arp2/3 subunits are AlphaFold models superposed on the fission-yeast template, 32–74% identity.
- **Capping:** capping protein (from 8F8Q) caps filaments after 20 ± 6 subunits. Filaments that end against an obstacle stay uncapped.
- **Orientation:** barbed ends face the membrane.
- **Collision checks:** residue by residue inside each patch. Some overlaps remain inside the built structures: 0.25–0.5% of their beads lie within 0.3 nm of a bead of another built actin component (tropomyosin on its filament, filaments within a bundle, Arp2/3 at its branch), in every build (81 of 17,154 sampled beads in build S, 101 of 36,786 in build R). The overlap check of section 9 reports them separately.

Not built: fimbrin cross-links, septins, formins.

### 6.5 Ribosomes and mRNA

- **Gene-specific mRNAs (2026-11-11; \`pipeline/mrna_genes.py\`):** every mRNA is a transcript of a named gene, drawn from per-gene molecules per cell (Pelechano *et al.* 2010, Table S1, a population mean of 25,900; times 0.7 for the G1 cell like the proteome, section 2.1: 18,100 mRNAs of about 3,600 genes in one draw). Its length is 5′ UTR + CDS + 3′ UTR + poly(A), with the CDS the summed coding segments of the SGD ORF (introns left out; see the fix below); UTRs and poly(A) are per gene where measured: TIF-seq UTRs in YPD for 4,298 genes (Pelechano *et al.* 2013 Nature 497:127) and PAL-seq total-RNA tails for 3,454 genes (Subtelny *et al.* 2014 Nature 508:66); other genes take 50, 120 and 33 nt (near the data medians). Ribosome load is weighted by per-gene translational efficiency from the same ribosome profiling (Subtelny *et al.* 2014; clipped to 0.1–10 times the median). Tables and their derivation: \`pipeline/data/mrna/README.md\`.
- **Polysomes:** 70% of cytosolic ribosomes are on mRNAs (assumption). Each mRNA carries Poisson(CDS × TE / s) ribosomes (TE: translational efficiency, above), with one global s chosen so that the total matches (167 nt of CDS per ribosome at TE 1 in builds R and S); mRNAs drawing none are free. This is close to the average density from polysome profiling quoted in section 14.3 (one ribosome per 156 nt; Arava *et al.* 2003), which was not used; per-gene densities from that study can serve as a test.
  - **Chain geometry:** the mRNA between two ribosomes is a flexible single strand of (spacing − 28 nt footprint) nucleotides; the distance from one ribosome's exit to the next one's entry is drawn from the end-to-end distribution of a Gaussian chain of that contour (persistence 1.5 nm, 0.6 nm per nt; rms about 14 nm), within a 70° cone around the exit direction. Earlier builds used a fixed 3–12 nm linker in a 40° cone and gave up after 60 failed tries, which ended many chains early; the polysome builder now keeps 95% of the requested ribosomes in an otherwise empty cell (60% before).
  - **Spacing (2026-11-12):** the free nucleotides between neighbouring ribosomes are drawn from an exponential distribution around the mRNA's mean (ribosomes load independently; before, every linker of an mRNA held the mean), and the walk keeps a drawn linker while it tries up to 40 orientations of the next ribosome (before, a short linker that collided was replaced by a new, mostly longer one), with up to 1,200 failed tries per chain. In an empty-cell test this raised ribosomes in tight pairs (exit to entry < 8 nm) from 2.8% to 13.8% of chain ribosomes and placement from 90% to 94%.
  - **Not modelled:** the stacked relative orientations of ribosomes in tightly packed polysomes (cryo-ET).
  - **mRNA channel:** the path through each ribosome comes from 7ZW0, superposed on 4V88.
  - **Collided ribosomes (2026-11-12, \`polysome.disome_pose\`):** where the drawn free mRNA between two neighbours is at most 5 nt, the trailing ribosome takes its pose in the yeast collided disome 6I7O (Ikeuchi *et al.* 2019 EMBO J 38:e100276; both ribosomes superposed on the library 80S by sequence-aligned rRNA phosphates, median deviation 0.10 nm). The cut-off is from disome footprints (57–63 nt against 27–29 nt for monosomes; Meydan & Guydosh 2020 Mol Cell 79:588). The coarse sphere trees overlap in the structure's pose, so the trailing ribosome is moved 4.75 nm apart along the line between the centres (26.4 instead of 21.6 nm; exit-to-entry 8.2 instead of 4.8 nm), keeping the structure's relative orientation; the poly(A)-stalled disome 6T83 (Tesina *et al.* 2020) lies within 15° and 2.9 nm of 6I7O. How often it happens follows from independent loading alone. If the pose clashes, the ribosome is placed as before with a 5-nt linker.
  - **Result (build S):** 116,723 ribosomes on 17,566 mRNAs (mean 6.6 per mRNA; 7.6 requested). Chains stop where they run out of space, so 88% of the 132,700 ribosomes meant for polysomes are in them (62% of cytosolic ribosomes, against the assumed 70%; build R 115,139 and 61%; the generic mRNAs before 2026-11-11 reached 44%); the rest are free. 3,704 ribosomes were placed collided behind their neighbour, so 3.8% of cytosolic ribosomes are in collided pairs, against 5.8% by disome profiling (Zhao *et al.* 2021 Genome Biol 22:16; test; build O had none). Exit-to-entry distances: median 13.7 nm (5–95% 7.8–31.9 nm); with the collided pairs, 10.2% of cytosolic ribosomes are in tight pairs (under 8 nm), against 14.5–18.6% by cryo-ET (build R 10.0%, build O 7.2%, build M 1.4%). These shares are the same as in builds Q and R: cytosolic crowding was kept when the cell was made smaller.
- **Full-length mRNA:** all mRNAs are threaded at full length through the cytosol before protein packing.
  - **Lengths:** gene-specific (above). The 5′ tail is the 5′ UTR plus the part of the CDS upstream of the ribosomes (assumed centred on the CDS); the 3′ tail the rest. Each thread records its gene (\`mrna_gene\` in the pre-pass).
  - **Representation:** one bead per nucleotide, 0.6 nm apart.
  - **Polysome tails:** the untranslated remainder leaves the first and last ribosome as 3′ and 5′ tails. The 5′ side gets 10–40% of it.
  - **Result (build S):** 85% of the target nucleotides were placed (as in build R). Threads that run out of room end early. In all, 19.2 million nucleotides of mRNA are drawn, including the ribosome channels (27.7 million in build P, with 1.4 times the mRNAs).
  - **Proteins on the ends (2026-11-11; \`pipeline/mrnp.py\`):** Pab1 on each 3′ end, max(1, poly(A) // 27 nt) copies up to 3 (a footprint of about 27 nt, unchecked; per-gene tails, so mostly one), and one eIF4E (Cdc33) on each 5′ cap. Each is an AlphaFold monomer at a random orientation with one of its spheres touching the end bead. Copies come out of the free pools (Pab1 100,100, Cdc33 28,800). In build S, 15,299 Pab1 (84% of the 18,158 wanted) and 17,529 eIF4E (97% of 5′ ends) fit (build R: 15,046 and 17,425). eIF4G, eIF4A and the closed loop are not modelled.
- **ER-bound ribosomes:** 5,675 ribosomes are docked on translocons in the ER and the outer nuclear membrane (5,586 fit in build S; 5,507 in build R).
  - **Templates:** 4,449 on Ssh1 (2WW9, yeast 80S–Ssh1) and 1,226 on Sec61 (7KAH, placed through the shared Sss1 subunit).
  - **Assumption:** 80% of free translocons are ribosome-bound.
  - **ER-associated ribosomes:** tomograms show 550–900 ribosomes per µm² on the cytoplasmic face of the cortical ER (West *et al.* 2011 J Cell Biol 193:333). The translocon pools fall far short of that, so more free ribosomes are packed in the band next to the ER and the outer nuclear membrane, up to 725 per µm² with the docked ones (the outer nuclear membrane is assumed to match). 11,000 and 5,400 are requested; the envelope face is capped since build S (section 4.4), which leaves 4,300 for it. In build S, 6,924 of 10,953 fit next to the ER and 3,378 of 4,341 next to the outer nuclear membrane (build R: 6,121 of 10,952 and 3,263 of 5,394); the density reached with the tethered polysomes is measured below.
  - **Not modelled:** ER-bound polysomes and their mRNA. The docked ribosomes are single ribosomes, not members of the chains of secretory mRNAs.
- **mRNAs at the mitochondrial surface (2026-11-11; \`pipeline/mrna_localization.py\`):** polysomes of genes whose translating ribosomes are enriched next to the mitochondrial outer membrane are tethered to it. The genes are the 223 −CHX calls of proximity-specific ribosome profiling in glucose (Williams *et al.* 2014 Science 346:748; 147–196 of them are expressed in one draw, about 1.5% of mRNAs).
  - **Tether:** the first ribosome starts in the cytosol 16–30 nm from the outer membrane, and every ribosome of the chain stays within 45 nm of it. After 300 failed starts a chain is placed untethered. Free (untranslated) mRNAs of these genes are not tethered.
  - **Test pre-pass:** the tethered ribosomes sit a median 30 nm from the outer membrane (99% within 50 nm, against 3% of the others), but only 69% of the requested ribosomes fit on the tethered chains.
  - **Held-out test:** the mitochondrial-fraction mRNA shares of Saint-Georges *et al.* 2008 (PLoS One 3:e2293; galactose) in \`validation/mrna_localization.py\`. Build S: 275 polysomes of 147 genes tethered, Spearman 0.45 and AUROC 0.66 over 561 genes. Build R, with the same mRNAs and tethers: 0.35 and 0.61. Builds N to Q, with 1.4 times the mRNAs: 0.41–0.44 and 0.66–0.68. With fewer mRNAs per gene, the share of a gene's mRNAs near mitochondria is a coarser number (often one or two molecules), and fewer class I genes are expressed in the draw; the rule is unchanged. Builds R and S share the mRNA draw and the tethers and differ in geometry and packing, so the difference between them (0.35 against 0.45) is a measure of how much one realisation moves this test. Build M: 386 polysomes of 170 genes tethered; over 610 genes Spearman 0.44 between each gene's share of mRNAs with a ribosome within 50 nm of the outer membrane and its mitochondrial-fraction share, AUROC 0.66 for class I against classes III–IV (untethered build: −0.05 and 0.47). The model is all-or-none per gene, while fractionation finds a median 43% of class I mRNAs on mitochondria.
- **mRNAs translated on the ER (2026-11-12, IMPROVEMENTS M5):** polysomes of the secretome genes (signal sequence or transmembrane domain) whose translating ribosomes are enriched next to the translocon are tethered to the cortical ER the same way. The genes come from ER proximity-specific ribosome profiling in glucose (Jan, Williams & Weissman 2014 Science 346:1257521, Table S5: BirA–Sec63 with cycloheximide, enrichment above the 95th percentile of non-secretome genes, their definition; 873 genes). Seeds lie 16–30 nm from the ER membrane and at least 45 nm from the plasma membrane: the cortical ER's face towards the plasma membrane is a ~30 nm gap without room for ribosomes. The translocon-docked single ribosomes stay as before; joining them into these chains is not done.
  - **Held-out test:** membrane-polysome enrichment per gene by fractionation in YPD (Chartron, Hunt & Frydman 2016 Nature 536:224, Supplementary Table 1) in \`validation/mrna_localization.py\`. Build M, without the ER tether: Spearman −0.08 between the share of a gene's mRNAs with a ribosome within 50 nm of the ER and its log2 enrichment (AUROC 0.49 for the top decile). Build S: 1,741 polysomes of 556 genes tethered; Spearman 0.35, AUROC 0.89 (build R the same; builds N to P 0.35 and 0.90–0.91). Both methods ultimately reflect signal sequences and transmembrane domains, so a sequence rule would also score well.
  - **ER ribosome density (\`validation/er_ribosome_density.py\`):** with the tethered polysomes, 19,300 ribosome centres lie within 30 nm of the cortical ER, 827 per µm² of its cytoplasmic face (23.3 µm²; West *et al.* 2011: 550–900; build R 796, build M about 470). Builds N to R agree within 2%; build S is 4% above build R; round N quoted about 790 from a one-off area estimate. Partly tuned: the ER tether was chosen with this number in view.
- **RAC on ribosomes:** the ribosome-associated complex (Zuo1–Ssz1; 63,600 copies) sits on the 60S subunit at the tunnel exit and is ribosome-bound in vivo, so every copy rides on a ribosome. In build S that is 39,200 on polysome ribosomes, where its own sphere tree must clear the neighbouring ribosomes, and 24,400 on free ribosomes (build R: 38,700 and 25,000).
  - **Pose:** from the *C. thermophilum* 80S–RAC structure 7Z3N (RAC-1). Its rRNA is superposed on the library 4V88 by 4,558 sequence-aligned phosphates (median deviation 0.27 nm).
  - **Stand-in:** its Zuo1 and Ssz1 (58% and 47% identical to the yeast proteins) stand in for the yeast RAC. The yeast RAC structure 7X3K was deposited without its ribosome and, superposed on 7Z3N, would float off it (2.4 nm deviation; another state).
- **Asc1 and Stm1 (2026-11-12, found by the PDB audit of section 5.2):** both are chains of the 80S structure 4V88 and were also placed as free proteins (316,000 and 188,000 copies in build P). They are now booked to the ribosome (\`ingredients.RIBO_EXTRA\`): Stm1 has fewer copies than there are ribosomes, so none is free; 21,300 Asc1 remain free. Stm1 sits in the mRNA channel of idle ribosomes (4V88 is from starved cells), so translating ribosomes are drawn with a protein they do not carry; a structure of a translating 80S is an open item of the audit.
- **Not placed on ribosomes:** NAC (Egd1–Egd2), Ssb, eEF3, eIF5A and eEF1A are free in the cytosol, and no tRNA is on a ribosome. Yeast structures exist for all but NAC (PDB audit, \`pdb_audit_2026-11-12/\`); not built yet.
- **tRNAs (2026-11-12; \`pipeline/trna_species.py\`):** one type per isoacceptor instead of one tRNA-Phe model (1EHZ) for all of them.
  - **Sequences:** the 275 nuclear tRNA genes of S288C in GtRNAdb (mature, introns removed; Chan & Lowe 2016 Nucleic Acids Res 44:D184) form 42 isoacceptors (anticodon families, initiator and elongator Met apart; two genes of undetermined identity are left out). Each is drawn as its most common sequence, with the 3′ CCA that yeast nuclear tRNA genes do not encode: 74–78 nt, and 85–87 nt for the class II tRNAs (Leu, Ser) with their long variable arm.
  - **Abundance:** the 2.31 million tRNAs (3.3 million for the population-mean cell, an assumption, times 0.7; section 2.2) are split by gene copies (tRNA levels follow gene copy number in growing yeast; Percudani *et al.* 1997 J Mol Biol 268:322). Build S places 2,309,616, all in the cytosol (build R 2,309,466).
  - **Structures:** predicted with Boltz-2 (Passaro *et al.* 2025, bioRxiv; an AlphaFold3-like model; best of 5 samples), confidence 0.73–0.89. The predicted tRNA-Phe is within 1.2 Å P-atom RMSD of the crystal structure 1EHZ (not an independent test: 1EHZ is surely in Boltz-2's training data). Modified nucleosides are not modelled. The viewer opens each model with "Load atomic model" (\`data/models/\`).
- **Introns (found and fixed 2026-11-12, build O):** through build N, the CDS length of each gene was its SGD ORF span, which includes introns.
  - **Affected:** the 208 intron-containing genes, mostly ribosomal-protein genes. RPL25, for example, got 843 nt instead of 429.
  - **Share of build N:** 18% of mRNA molecules and 17.5% of polysome ribosomes. Intron sequence was 5.2% of all mRNA nucleotides. These mRNAs were too long and, since ribosome load scales with CDS length, carried too many ribosomes.
  - **Fix:** the CDS is now the sum of the ORF's SGD CDS segments. TIF-seq 5′ UTRs are genomic distances, so a 5′ UTR intron (SGD; 21 of the genes here, mostly ribosomal-protein genes) is subtracted where the measured UTR is longer than the intron. RPS8B's 5′ UTR, for example, goes from 389 to 29 nt.
  - **Effect in build O:** 5.8% fewer mRNA nucleotides; ribosomes on mRNAs 168,017 instead of 171,872 (median 5 per mRNA instead of 6); the mRNA tests are unchanged within noise (above).
- **In the viewer (2026-11-12; \`pipeline/polymer_annot.py\`):** hovering or clicking an mRNA names its gene and the part it belongs to: the 5′ tail, the stretch through ribosome *k* of *n* (counted from the 5′ end), the 3′ tail, or a free mRNA.
  - **Positions:** nucleotides are counted from the end of the tail, the cap or the 3′ end, where eIF4E and Pab1 sit. A thread that ran out of room ends early, and its drawn end counts as its end. The regions (5′ UTR, CDS, 3′ UTR, poly(A)) follow from the gene's layout.
  - **Between ribosomes:** the mRNA there is drawn shorter than the nucleotides it stands for, since each linker is its end-to-end distance. Those pieces are named by their ribosome only.
  - **Ribosomes:** polysome ribosomes name the mRNA they translate and their rank on it.

### 6.6 Ribosome assembly: pre-ribosomes and the rDNA

- **Pre-ribosomes:** 19,500 assembly intermediates in five stages, from the nucleolus to the cytoplasm (structures and counts in section 5.2). They are packed like other complexes in their compartments, except the 90S particles still attached to nascent transcripts (below).
- **rDNA:** about 150 tandem repeats of 9.1 kb on chromosome XII. Half of them (75) are active (psoralen accessibility; Dammann *et al.* 1995).
  - **Active genes:** nucleosome-depleted (Merz *et al.* 2008; later ChIP studies argue for unphased nucleosomes instead), with 50 RNA polymerase I complexes on the 6.7 kb 35S transcription unit, one every 134 bp (Miller spreads: 50 ± 20; French *et al.* 2003). At 60 nt/s this finishes about 2,000 transcripts per minute, the ribosome production rate.
  - **Geometry:** each gene is a persistent self-avoiding walk of polymerase sites, 45.6 nm of DNA apart, inside the nucleolus. The DNA (2 bp per bead) is a smooth spline through the sites and passes through each polymerase's cleft: Pol I is oriented by the DNA axis of the yeast Pol I elongation complex 5M3F, superposed on the library's 4C2M. The walks are not connected into the real tandem array.
  - **Nascent pre-rRNA:** one transcript per polymerase from its RNA exit, as long as the polymerase has transcribed (0–6.7 kb), threaded self-avoiding at 2 nt per bead (0.6 nm per bead, about the rise of base-paired RNA). About 70% of transcripts are cleaved at site A2 while still being made (Kos & Tollervey 2010): those keep only the part beyond A2. Uncleaved transcripts beyond A2 carry a 90S particle on their 5′ end (5′ETS and 18S inside it).
  - **Inactive repeats and the spacers of active ones:** nucleosomal chromatin in the nucleolus, 5,227 nucleosomes in persistent walks.
- **Result (build S):** 73 of the 75 genes fit, with 3,650 Pol I (of 10,000; the rest are free in the nucleolus), 329 transcripts with a 90S (of 666 that should have one), and 57% of the nascent RNA placed (build R: 3,601 Pol I, 349 of 649, 57%). The number of active genes and their polymerases are per genome and do not scale with the cell, so the smaller nucleolus is more crowded: 170 µm of transcribed DNA in 0.47 µm³ (0.67 before build R) puts gene axes about 50 nm apart, and threads that run out of room end early. Real nascent transcripts are compact RNPs, so they need less room than these threads.

### 6.7 Mitochondrial inner membrane

- **ATP synthase dimers:** 3,880 of the 15,300 ATP synthase copies form 1,940 dimers (build S), which fill the rim sites that fit (72% of 2,703 tried); the rest stay monomers in the flat inner membrane. Build R, with a crista across the tubule every 42 nm, had 6,257 dimers (52% of 11,935 sites). The short cristae of glucose-grown cells (section 3.3) have far less rim, so most ATP synthase is now monomeric in the inner boundary membrane; in glucose cells ATP synthase is 2–3 times less abundant per mitochondrial mass than in respiring cells (Klecker & Westermann 2021, citing Di Bartolomeo *et al.* 2020); the model's copy number is PaxDb's (section 2.1). Each dimer is two copies of the yeast monomer 8F39 arranged as in the cryo-ET dimer 4B2Q (dimer angle 83.5°).
  - **Placement:** in rows along the crista rims, every 13 nm (assumed spacing). The residue test of each rim site uses only the geometry near its crista and is compiled (build S; same tests and the same result on build R's geometry).
  - **Missing subunits:** 8F39 lacks the dimer-specific subunits e and k.
- **Respiratory supercomplexes:** complex IV is entirely in 2,031 III₂IV₂ supercomplexes (6HU9). The remaining 1,850 complex III dimers are free (builds R and S).
- **A crowded inner membrane (build S):** the inner membrane has 13.1 µm² instead of 31.5 for 309,800 requested membrane-protein copies (build R: 301,200; more ATP synthase stays monomeric now), and holds 87.4% of them (build R: 98.2%). Carriers lose most: Yhm2 40%, Ggc1 42%, Odc1 and Odc2 37–38% placed. Membrane regions are not thinned (section 8.2).
- **Mitoribosomes on the membrane:** in yeast the vast majority of mitoribosomes are tethered to the inner membrane, exit tunnel towards it, through Mba1 and the 21S rRNA expansion segment 96-ES1 (cryo-ET; Pfeffer *et al.* 2015 Nat Commun 6:6019).
  - **Structure:** the 74S monosome of 9HM0 (yeast mitoribosome with mRNA; the translational activator Aep3 left out).
  - **Frame:** the axis from the large subunit's centre to the exit-site proteins (uL23m MRP20, uL24m MRPL40, uL29m MRPL4) is the membrane normal, pointing into the matrix; the monosome is lowered until its nearest bead is 1 nm above the matrix-side leaflet, which puts the exit site 5.1 nm above it (the gap Mba1 and Oxa1 bridge). No residue enters the bilayer.
  - **Copies:** 90% of the smaller subunit pool (assumption from "vast majority"): 1,367 monosomes, of which 1,364 fit in build S (1,311 in build R), along the inner boundary membrane and at crista tips. In build R, with cristae 24 nm apart along the whole tubule, the 37 nm tall monosome had no room between them. The remaining 152 large and 169 small subunits are free in the matrix.
- **mtDNA nucleoids (revised for build S, 2026-11-13):** a haploid cell in glucose has 1.2–1.35 nucleoids per µm of tubule, about 0.8 µm apart, with one within 0.5 µm of 69% of the tubule tips (Osman *et al.* 2015), and about 21 ± 4 copies of the 85.8 kb genome per haploid cell in YPD (Göke *et al.* 2020), about 0.45 per fL of cell (Seel *et al.* 2023; both W303). The geometry chooses the sites (1.27 per µm of tubule, near tips first, at least 550 nm apart and 250 nm from branch points, off the axis towards the inner boundary membrane), and no crista starts within 110 nm of one: a nucleoid is a compact body in crista-free matrix (size from Göke *et al.* 2020; position by analogy with mammalian cryo-ET, Kukat *et al.* 2015; no yeast EM places nucleoids relative to cristae). Build S: 11 nucleoids holding 15 genome copies (0.45 per fL × 34 fL; 1–3 per nucleoid), 1.27 Mb of DNA drawn at 2 bp per bead. Up to build R: 15 nucleoids anywhere in the matrix with 35 copies (haploid cells: 10–20 nucleoids, 25–50 genome copies, as summarised by Zelenaya-Troitskaya *et al.* 1998 Genetics 148:1763), 3.0 Mb.
  - **Walks:** each nucleoid is a cluster of 5 kb loops, self-avoiding walks (10 bp per 3.4 nm step) that start near the nucleoid centre and are drawn back to it. Ninety percent of each nucleoid's DNA lies within 53 nm (median) of its centre (build R, threaded between cristae every 42 nm: 63 nm).
  - **Abf2:** the packaging protein (35,600 copies) is packed free in the matrix.
- **Not modelled:** III₂IV₁ supercomplexes, confinement of supercomplexes to the flat crista membranes, and mitochondrial polysomes.

### 6.8 Eisosomes

- **Structure:** the near-native Pil1 and Lsp1 lattices of Kefauver *et al.* 2024 (8QB7, 8QB8). These were solved on eisosome tubules isolated from yeast plasma membranes, with the membrane lipids still bound.
- **Lattice:** the seven dimers of each entry are related by screw operations about the tube axis (0.03 nm rms). Two of them span the lattice: 13.6° with −4.0 nm along the axis, and −37.7° with −3.1 nm. That gives 44.8 nm² per dimer at the coat's mid-radius.
- **Membrane:** the amphipathic helix (residues 39–48) sits 11.5 nm from the tube axis. Taking it as 1.3 nm above the bilayer midplane (assumption: an interfacial helix in the cytoplasmic leaflet) puts the midplane at 10.2 nm, which is the furrow-bottom radius in section 3.2.
- **Coat:** whole dimers on the half-cylinder bottom of each furrow, within ±60° of its lowest line, so that no dimer runs onto the flat walls. The lattice phase is random per furrow. Each site holds a Pil1 or an Lsp1 dimer in proportion to the pools; the two lattices are nearly identical.
- **Result (build S):** 3,757 dimers on 23 furrows (2,013 Pil1, 1,744 Lsp1; 163 per furrow), 7.6% of the cell's Pil1 and Lsp1 (build R: 3,660 on 22 furrows, 7.4%). The rest stays near the membrane as free monomers and dimers. At about 2,500 molecules per eisosome (unchecked), the measured copy numbers imply coats much larger than one half-pipe 300 nm long: longer furrows, coated walls, or a cytosolic pool.
- **MCC membrane proteins:** the furrows are the membrane compartment of Can1 (MCC). Its tetraspanning markers Sur7 and Nce102 are concentrated there (Malínská *et al.* 2003; Grossmann *et al.* 2008). 80% of their copies (assumption) go onto the flat furrow walls, oriented like other plasma-membrane proteins, as many as fit. The rest stay in the rest of the plasma membrane.
- **Not built:** Can1 and other transporters that move in and out of the MCC with conditions; PI(4,5)P₂ or sterol enrichment in the furrows; the stretched state of the lattice.

### 6.9 Cvt complex

- **What it is:** the cytoplasm-to-vacuole targeting cargo of growing cells. Precursor aminopeptidase I (prApe1) dodecamers assemble through their propeptides into one large complex with α-mannosidase Ams1, coated by the receptor Atg19. It waits at the vacuole to be wrapped into a Cvt vesicle of about 150 nm (Baba *et al.* 1997; Lynch-Day & Klionsky 2010).
- **Built:** one cluster in the cytosol, 40 nm from the vacuole. It holds half of the Ape1 dodecamers (the library 12-mer 5JGF; the propeptides are not resolved) and half of the Ams1 dimers. Both halves are an assumption for the share still in precursor form. They are packed at random to 35% volume fraction, with all Atg19 on the surface.
- **Result (build S):** radius 84 nm (168 nm across, close to a 150 nm Cvt vesicle; build R 162 nm, 180 nm at 5.0 pg), 229 Ape1 dodecamers, 137 Ams1 dimers, 254 Atg19. The other halves of Ape1 and Ams1 stay in the vacuole as mature enzymes.

### 6.10 Nucleus–vacuole junction

- **What it is:** where the vacuole lies against the nucleus, the outer-nuclear-membrane protein Nvj1 binds Vac8 on the vacuolar membrane (Pan *et al.* 2000 Mol Biol Cell 11:2445).
- **Built:** all 532 Nvj1, each with one Vac8, as the Vac8–Nvj1 crystal structure 5XJG (Jeong *et al.* 2017 PNAS 114:E4539). They sit on the cytoplasmic face of the vacuolar membrane wherever it is within 45 nm of the outer nuclear membrane, which in build S is 1.01 µm² (build R, with 12 lobes: 0.86 µm²; 0.46 µm² in the larger cell of the builds up to Q).
- **No pores (build S):** the envelope of the junction carries no nuclear pore (section 6.1). Vac8 is against the membrane (its palmitoylated N-terminus is not resolved), with the Nvj1 fragment towards the nucleus.
- **Not drawn:** Nvj1's transmembrane anchor and linker, which are not resolved.

## 7. Procedural content (generated in the browser)

These are generated deterministically from fixed seeds when you zoom in. They are not stored coordinates and do not count in the densities.

- **Lipids:**
  - **Packing:** 0.65 nm² per lipid per leaflet.
  - **Composition per membrane** (Zinser *et al.* 1991 J Bacteriol 173:2026, isolated organelle membranes): phospholipid classes (PC, PE, PI, PS, PA, cardiolipin, other) from their Table 5 and ergosterol from the ergosterol : phospholipid ratios of their Table 4. ER: light microsomes; Golgi/endosomes: secretory vesicles; nuclear envelope: nuclei; lipid-droplet monolayer: ER phospholipids without sterol. Examples: inner mitochondrial membrane 12% cardiolipin, outer 6%; plasma membrane 33% of phospholipid PS; ER 41% PC and 21% ergosterol.
  - **Plasma membrane:** ergosterol 40 mol% (assumption; the measured ratio of 3.3 would mean 77 mol%, more than a bilayer holds). Leaflets (consensus): PS in the cytoplasmic leaflet, PE 80% there, sphingolipids (in "other") in the exoplasmic leaflet. Other membranes are drawn symmetric.
  - **Caveats:** mitochondria, microsomes and peroxisomes were isolated from lactate- or oleate-grown cells, the others from YPD; fractions carry some cross-contamination.
  - **Masking:** lipids within 0.95 nm of a protein bead are removed.
  - **In the viewer:** clicking a membrane lists the composition of that membrane (and leaflet). In cut views the zoom level at which lipids are drawn is tested at the cut slab for every tile (2026-11-13); tested per tile, it drew lipids only on the tiles nearest the eye, a wedge on membranes seen face-on.
- **Cell wall:**
  - **Layers:** an inner β-1,3-glucan layer (semi-flexible chains, 0.018 beads/nm³) and an outer mannan layer (branched brushes on the mannoproteins).
  - **Amount:** bead counts are set so the wall is about 20% of cell dry mass.
  - **Chitin:** 1.5% of the wall (Klis *et al.* 2002). Half of it (assumption: chitin is mostly in bud scars) fills the bud-scar rim as microfibrils running around the ring, one bead per 1.04 nm of an ~18-chain α-chitin microfibril (36 GlcNAc per bead, like the glucan beads).
  - **Collisions:** only partly checked.
- **Vacuolar polyphosphate:** the vacuole is the cell's phosphate store. The VTC complex in its membrane synthesises polyphosphate (polyP) straight into the lumen.
  - **Amount:** reported as "hundreds of millimolar" of phosphate units (Bru *et al.* 2025 eLife 14:RP108181). Here it is 300 mM in the lumen (assumption), which averages to about 27 mM over the cell. Whole-cell measurements range from about 23 to 250 mM; the location is disputed between older fractionations and the current vacuolar view.
  - **Chains:** short and linear, lognormal lengths around 24 units (the mean measured in polyP-rich yeast, Christ *et al.* 2020).
  - **Drawing:** one bead per 4 phosphate units (1.2 nm of chain), generated only inside the cut slab, kept 3 nm inside the vacuolar membrane and clear of vacuolar proteins.
  - **Not drawn:** counter-ions (Mg²⁺, K⁺, basic amino acids). It does not count in the protein densities.
- **Motion (Display → Motion; off by default):** illustrative animation, not molecular dynamics. Nothing is simulated and neighbouring molecules are not checked for collisions while they move.
  - **Wobble and flexing:** a small rigid-body wobble whose amplitude falls with molecular size (about 0.35 nm translation for small proteins). Membrane proteins move only within the bilayer and rotate only about the membrane normal. Residues with low pLDDT flex by up to 0.45 nm. Free polymers (mRNA, linker DNA, wall polysaccharides, FG repeats) wiggle along their length.
  - **Static:** actin filaments, microtubules, nuclear pores and nucleosomes stay still, so that their assemblies stay intact.

## 8. Packing

### 8.1 Spatial decomposition

- **Tiles:** the cell is cut into 200 nm cubic tiles; 5,912 of them intersect the cell (builds R and S; 6,416 in the larger cell).
- **Copies per tile:** each species' expected copy number per compartment is split over tiles in proportion to the compartment's volume (or membrane area) in each tile, modulated by the large-scale organisation of section 8.4. A multinomial draw with stochastic rounding conserves the cell totals.
- **Concentrations:** total macromolecular crowding is uniform within each compartment; the composition varies (section 8.4).

### 8.2 Algorithm

Random sequential addition in each tile, in rounds by size class (since build S).

**Before packing: what each region is asked to hold (build S, 2026-11-13; \`pipeline/wc3d.py\`).** The measure is the requested sphere volume fraction of a region: the summed sphere volume of the packing shapes of its copies per region volume (the spheres of one molecule overlap, so this exceeds the fraction they really fill). The cytosol is asked for 0.41.
- **Membrane-face shells** are capped at the cytosol's fraction; what a shell is asked for beyond it goes to the cytosol, every species by the same factor (section 4.4: Golgi/endosome face ×0.61, nuclear-envelope face ×0.80).
- **Over-full volume regions are thinned uniformly.** Regions asking for more than 1.1 times the cytosol's fraction are packed first in a pilot: 16 tiles per region, drawn in proportion to the region's volume in them, not adjacent, without neighbours, one round. A region whose pilot tiles held less than 98% of their requested sphere volume has the request of every one of its species multiplied by 0.97 times the share held (\`work/yc3d/thinning.json\`). In build S the nucleolus asked for 0.64 and its pilot tiles held 0.935, so every nucleolar species is packed at 0.907 of its inventory copies; the matrix asked for 0.56 and held 0.991, so it is not thinned. This replaces a shortfall that fell on mid-sized species (build R: Hsp60 49%, Sod2 37% and Ssc1 96% placed in the matrix) by one reported factor. The inventory itself is not changed, and the validation reports placed copies against the inventory as well as against the request. Membrane regions are not thinned. (Candidate build 1 drew the pilot tiles with the most of the region, its roomiest parts; they overstated what the matrix holds.)

**Placement:**
- **Rounds (build S):** every tile is visited in 6 rounds. Round 0 places the membrane proteins and the largest soluble species; later rounds place smaller size classes, split by bounding radius so that each class holds about the same requested sphere volume (thresholds 15.6, 8.1, 7.1, 5.6 and 5.0 nm in build S). Each round packs the tiles in the 8 checkerboard passes below, in the reverse pass order in every second round, so that each tile is as often ahead of a neighbour as behind it. Within a round, molecules go in by decreasing size. Up to build R each tile was packed once, membrane proteins first and then the rest by decreasing size: first-pass tiles over-filled the 5 nm at their faces (molecules protrude into still empty neighbours) and last-pass tiles under-filled it (section 9).
- **Soluble molecules:** a uniformly random position in the target compartment and a uniformly random orientation, or, for a share of copies, a position next to an interaction partner (section 8.4).
- **Membrane proteins:** projected onto the exact membrane surface, with their membrane normal along the local surface normal and a random spin about it.
- **Retries (2026-11-12, after the Jurkat cell model):** up to 2,000 trials per copy, and a job stops after 1,200 copies in a row fail (was 200 and 25). With the old budget tiles stopped early: on test tiles the Golgi face rose from 57% to 66% placed, the envelope periphery from 43% to 51%, the vacuole face from 57% to 63%, the cytosol of crowded tiles from 92% to 100%, the mitochondrial matrix only from 47% to 51% (space-limited).
- **Flat alignment:** in thin lumens and membrane shells (ER lumen, matrix, intermembrane space, Golgi lumen and the peripheral shells) half of the trials align the molecule's flattest principal axis with the normal of the bounding membrane; a small further gain in shells (vacuole face 63% to 67% on a test tile).
- **Membrane-face shells (build S):** a membrane-face protein proposed next to a partner (section 8.4), or moved by the squeeze below, must still have its centre within its shell width plus 4 nm of its membrane. In build R 4–17% of the placed copies of these shells had ended outside them (Golgi face 17%, envelope 11%, lipid droplets 7%, vacuole 7%, ER 5%, plasma membrane 4%).
- **Squeeze relaxation:** random sequential addition jams at roughly 35–40% sphere-tree volume fraction, short of cellular crowding. In the last round (since build S; before, after the single pass), in a tile where copies did not fit, those copies are inserted at random positions in their compartment, in rounds of at most 15% of the tile's soluble molecules (4 rounds).
  - **Relaxation:** after each insertion round, all soluble molecules of the tile are relaxed as rigid bodies for 60 iterations. Overlapping spheres of different molecules push each other apart (translation and a small rotation per molecule, at most 0.6 nm and 0.06 rad per iteration). Spheres crossing a membrane or leaving their compartment are pushed back along the membrane's distance gradient, and centres stay in the tile.
  - **Fixed:** membrane proteins, pre-placed structures and finished neighbour tiles do not move.
  - **Acceptance:** the configuration is then made exact again. Molecules are accepted one by one, previously placed ones first, only if they overlap nothing accepted so far and lie in their compartment, as in RSA. A round that does not increase the count is undone.
  - **Effect** on a test tile of the mitochondrial matrix: 72% → 89% placed at the earlier protein total, and 49% → 56% at the current one (with the library sphere trees; not re-measured with the packing shapes below). Cytosol and nucleoplasm tiles fit everything by RSA alone.
  - Copies that still don't fit are left out.

**Collision checks:**
- **Packing shapes (build S, 2026-11-13; \`wc3d.PACK_PER = 15\`):** for packing, each structure is approximated by k-means clusters of about 15 residues (or nucleotides) per sphere, made from the beads and cached in \`work/yc3d/pack_trees_15.npz\`; each sphere is padded by a 0.42 nm bead radius, and spheres may interpenetrate by up to 0.3 nm. The library's sphere trees of about 30 residues per sphere (40 for the curated complexes) are still what the viewer's simulation and the contact statistics use. As packing shapes they were fat: their union is 2.4 times the union of the residue beads for the most abundant matrix proteins, 1.9 times in the cytosol and 3.8 times in the nucleolus. On six matrix-rich test tiles of build R, 87.9% of the matrix request was placed with the library trees, 93.1% with 15 residues per sphere (2.6 times the CPU time) and 94.9% with 8 (7 times); species under 80% placed: 13, 3 and 9.
- **Compartments:** a copy must lie entirely within its compartment and clear of membranes, which are treated as 2 nm half-thickness slabs.

**Parallelism:**
- **Phases:** in each round, tiles are packed in 8 checkerboard passes, so neighbouring tiles are never packed simultaneously. A tile packed in a later round reads back what earlier rounds placed in it and its neighbours.
- **Halos:** each tile sees the already placed molecules of its neighbours through a halo.
- **Obstacles:** the pre-placed structures of section 6 are fixed obstacles.

### 8.3 Consequences

- **The density is not an equilibrium value.** Random sequential addition jams at lower density than real cytoplasm, and the outcome depends on insertion order; the squeeze relaxation (section 8.2) pushes past the jamming limit but is not an equilibrium either. Where the requested copies exceed what fits, the shortfall is reported, not hidden; where a compartment is thinned before packing (the nucleolus, by 0.907), the factor is reported.
- **Tile seams are gone (build S).** In 918 fully cytosolic tiles the sphere volume fraction within 5 nm of a tile face is 0.41 in first-pass and 0.40 in last-pass tiles, against 0.41 in the interior (15–35 nm from the faces). Build R, packed in one visit per tile: 0.49 and 0.20 against 0.40.
- **Membrane-face shells (build S):** none of the sampled membrane-face copies lies outside its shell (build R: 4–17%). Cytosolic proteins are still sparse in the 4–20 nm next to the Golgi/endosome membranes, where the shell is (0.16 of their density further out; vacuole 0.56, nuclear envelope 0.68), while all molecules together are at 1.25 times the density further out there (\`placement_checks.json\`, halos).
- **The shortfall is smaller and spread more evenly (build S).** Placed copies are 99.3% of the packing request and 99.1% of the inventory, 99.5% outside mitochondria, and 98.1% of the inventory by sphere volume (build R: 97.5% and 98.9% of the request, 96.0% by volume; builds N to Q 96.5–96.7% and 98.5–98.6%; build M 95.9 / 97.9% before the larger attempt budget). The cytosol (99.7%) and nucleoplasm (99.8%) fit nearly everything. The matrix holds 98.5% of its copies at 178 mg/mL (build R: 75% at 150 mg/mL; 132 in build P, 107 with mitochondria at 2% of the cell), with the finer packing shapes and a larger matrix (0.82 against 0.70 µm³; round S candidate build 1, with the old mitochondria and the library trees, held 73.6% of its inventory). The nucleolus holds 92.7% of its thinned request, 84.1% of its inventory (build R 80%). Shortfalls that remain: the mitochondrial inner membrane (87.4%; section 6.7), the cytosolic faces of the lipid droplets (95.0%), the ER (95.3%) and the nuclear envelope (96.8%); the Golgi face holds 97.7% of a capped request (build R 84% of the full one). Species placed at under half of their inventory copies: 10 of 4,469 (35,000 requested copies; build R 76); under 80%: 60 (build R 238; \`validation/placement_checks.py\`). Counted per region, 12 species fall below half of the region's share: six inner-membrane carriers (Yhm2, Ggc1, Odc1, Odc2, Mrs3, Yea6), three nucleolar proteins (Tif3, Esf1, Rrp12), two Golgi-face complexes (AP-3, Ric1–Rgp1) and the ER membrane protein Frt2.
- **Interaction-driven organisation is a hypothesis layer.** Section 8.4 places partners near each other and co-enriches network neighbours, but with assumed strengths. There are no explicit condensates, and the layer is not fitted to any measured protein positions.

### 8.4 Organisation from association data

Uniform random placement is itself an assumption: that proteins have no spatial organisation beyond their compartment. The model replaces it with a weak, data-driven organisation at two scales. It never moves a protein out of its compartment: it only changes where in its compartments a protein goes, and next to whom. Neither scale is calibrated against measured protein positions, which hardly exist at this resolution. The strengths are tuned so that the packing stays random-looking but is measurably organised, and total crowding stays uniform. Code: \`pipeline/organisation.py\`; \`scripts/s3_pack.py --random\` builds the same cell without this layer.

**Data:**
- **STRING v12** (Szklarczyk *et al.* 2023, Nucleic Acids Res. 51:D638): 271,240 links with combined score ≥ 400 (medium confidence) between 5,881 yeast proteins, weighted by score. STRING includes indirect associations (shared pathways, co-expression, text mining) that a structure prediction cannot see.
- **AlphaFold 3 pool predictions** (in-house, unpublished): every protein was predicted in many random multi-protein pools, giving a chain-pair score for about 10.3 million protein pairs. The best (lowest) chain-pair PAE min-min per pair is used.
  - **Calibration against known interactions:** pairs at ≤ 2 Å are 73% Complex Portal pairs and 72% have IntAct evidence; pairs at 2–4 Å, 43% and 48%; pairs at 6–8 Å are still 50 and 10 times above the background of 0.06% and 0.6%.
  - **Used:** 15,262 pairs at ≤ 8 Å after filtering. A pair is dropped when a paralog of one partner scores better with the other, which removes homologous-interface artefacts such as Pol II subunits scored against their Pol III paralogs.

**Network neighbourhoods (no clustering).** An earlier version split STRING into modules by Louvain clustering; the module resolution and minimum size were arbitrary and changed the result. The current version uses network propagation instead:
- **Profiles:** each protein gets its personalized PageRank vector on the weighted STRING graph: the stationary distribution of a random walk that continues along an edge with probability 0.5 (the damping) and otherwise restarts at the protein. A low damping keeps the profile local, concentrated on the protein's near neighbourhood. The protein's own entry is dropped and the vector renormalised. A complex takes the mean profile of its subunits; 5,452 of the 5,562 species of build S have a profile (the rest are prebuilt assemblies, tRNA and proteins without STRING links; 5,542 of 5,605 in an earlier build).
- **Similarity:** the Pearson correlation *r* of two profiles. It is high for proteins that share a network neighbourhood, whether or not they are linked directly. Typical neighbourhoods are small: 1.0% of species pairs have *r* ≥ 0.4. Examples: TDH3's closest neighbours are TDH1 (0.92), TDH2 (0.91), PGK1 (0.78), TPI1, ENO1, GPM1, PGI1 and FBA1; Sec61's are Sec63, Sec62, the SRP receptor, Ssh1 and OST.
- **Choice of parameters** (\`scripts/organisation_sweep.py\`): STRING thresholds 400 and 700, damping 0.3–0.85 and four profile variants were scored by how well *r* predicts co-complex membership (Complex Portal co-members against random pairs) and shared localization (same coarse compartment, for proteins with a dominant location). Every variant separates co-complex pairs almost perfectly (AUROC 0.993–0.998); shared localization is predicted only moderately (AUROC 0.56–0.61), much as for pairs that are not directly linked. Log-transformed profiles scored slightly higher on localization but are not local (a typical protein has thousands of partners at *r* > 0.3). Among the local variants, STRING ≥ 400 with damping 0.5 scored best (localization 0.594, co-complex 0.996).

**Large scale (hundreds of nm): network-informed sampling of where each species goes.**
- **Correlated fields:** the species correlation matrix *R* is used as the covariance of the species' log-density fields over the cell. It is approximated as *R* ≈ *E Eᵀ* from its 128 leading eigenvectors (54% of the trace). Each of the 128 components gets an independent smooth Gaussian random field *F*ₖ over the cell (correlation length 500 nm, unit variance). A species' copies are distributed over tiles in proportion to region volume × exp(0.8 × Σₖ *E*ᵢₖ *F*ₖ).
- **Interpretation:** the log-densities of two species co-vary across the cell by 0.8² × *r*: network neighbours are co-enriched in the same places, unrelated species vary independently, and species without a profile are spread uniformly. Nothing is clustered and no module sizes are chosen.
- **Crowding:** an iterative correction keeps each tile's total macromolecular volume proportional to its compartment volume, so composition varies between regions but crowding does not.
- **Effect** (abundant species in fully cytosolic tiles): the correlation of two species' tile densities follows their network similarity (Spearman 0.58 over 86 species, 0.71 over the 24 most abundant; the module version gave 0.44 on the same measure). Pairs with *r* ≥ 0.5 have tile densities correlated at +0.56; unrelated abundant pairs are slightly anti-correlated (−0.18), because the crowding correction makes abundant species compete for volume within a tile. Total crowding varies as little as with random placement (CV 0.040 vs 0.041).

**Meso and small scale (nm to tens of nm): partners placed next to each other.**
- **Contact partners:** pool pairs. Strength falls from 1 at PAE ≤ 2 Å to 0 at 8 Å, is multiplied by (1 + *r*⁺)/2 so that pairs that also share a network neighbourhood pull fully and network-unsupported pairs at half strength, and is divided by √(degree/20) for proteins with more than 20 pool partners, so hubs do not dominate.
- **Loose partners:** the species most correlated with it (*r* ≥ 0.4) that are not contact partners. Strength rises from 0.2 at *r* = 0.4 to 1 at *r* = 1. Each species keeps up to 16 partners, contact partners first.
- **Attraction:** when a copy is placed, each partner p contributes a_p = strength × A × n_p / (8 + n_p). Here n_p is the number of copies of p already placed in the tile and its halo, and A is 0.6 for contact and 0.25 for loose partners. The pull saturates with the partner's local copy number, like mass action.
  - With probability min(0.6, Σ a_p), the copy is proposed next to a randomly chosen placed copy of a partner (chosen in proportion to a_p): in contact for contact partners, with a 2–15 nm gap for loose partners.
  - The copy gets up to 40 such proposals; if none passes the collision and compartment checks, it is placed at random as before.
- **Effect:** 7.1 M copies drew a partner and 5.1 M of them (about 24% of all placed molecules) were placed next to one. A molecule's contact partner lies within contact distance 8.4% of the time, against 1.3% for a control species of the same local abundance, and a loose partner within 15 nm 19.9% of the time (control 14.0%). Placement completeness is unchanged (97.8%). Since build S a proposal next to a partner must, for a membrane-face protein, also lie in its shell (section 8.2). These shares were measured on an earlier build; in build S, 8.7 M copies tried for a partner and 4.9 M (20% of placed molecules) were placed next to one (build R: 9.5 M and 4.7 M, 19%).

**What it is not:** strengths and bound fractions are assumptions, as no affinities are used; the partners are not bound in a specific orientation (pre-docked pairs are a separate future step); and smaller species are pulled towards larger ones that are already placed, since species are placed in order of decreasing size. The same profiles could inform any later packing method (for example a denser, relaxed packing): the layer only supplies a covariance for where species go and a list of partners.

**In the viewer:** Display → colour "Network neighbourhood (STRING)" colours each species from a two-dimensional t-SNE map of the profiles (distance 1 − *r*), so network neighbours get similar colours. Clicking a molecule and then "Network neighbours" highlights it with its 30 most correlated species (*r* ≥ 0.2) and lists them.

## 9. Validation

| Check | Result |
|---|---|
| Residue-bead clashes (beads of different molecules closer than 0.3 nm) | Packed molecules: 3.4 × 10⁻⁶ of 5.3 M sampled beads in 64 test cubes (18 beads, none between tiles; build S; build R 6.4 × 10⁻⁶, builds N to Q 3.7–9.3 × 10⁻⁶), under the 10⁻⁵ acceptance limit; 8 of the 18 are Pab1 or eIF4E on the mRNA end it binds, which they are allowed to overlap. Including the built actin structures: 1.87 × 10⁻⁵ (99 beads; build R 2.65 × 10⁻⁵, 133 beads), above the limit. The other 81 beads are overlaps inside those structures (tropomyosin on its filament 33, filament on filament 18; section 6.4), 0.25–0.5% of their beads in every build; the smaller cell has the same cortex in fewer cubes' worth of volume, so the sample holds more of them. The two rates are reported separately since build R (\`clash_rate_without_actin_structures\`; four other cube samples of build R gave 3.8–6.6 × 10⁻⁶). The squeeze relaxation presses molecules to the 0.3 nm sphere-tree tolerance far more often than random addition alone (6.3 × 10⁻⁷ before it) |
| Containment (beads outside their compartment) | 6.7 × 10⁻⁵ of 19 M beads; 0.047% of 30,000 sampled instances have a bead outside (build S; build R 4.1 × 10⁻⁵ and 0.04%; round S candidate build 1 9.1 × 10⁻⁵). Relaxed molecules sit against compartment walls, and the packer tests packing shapes, not residue beads. Why the rate rose from build R is not yet analysed (section 11) |
| Membrane-protein orientation vs local surface normal | median 0.0°, all of 10,000 sampled within 1° |
| TM residues within 2 nm of the membrane mid-surface | 97.4% (median 0.72 nm; build R 97.3%, 0.73 nm) |
| Telomeres (Rabl layout) | 32 of 32 within 150 nm of the nuclear envelope, median 164 nm from their targets; telomere–pole distance vs arm length Spearman 0.92 (build S; build R 0.85, earlier builds 0.86–0.90; section 6.3) |
| Eisosome coat: membrane-inserted residues within 2 nm of the furrow mid-surface | all sampled (median 1.36 nm) and all other residues on the cytoplasmic side (builds R and S; an earlier build: 99.2% and 98.9%) |
| ATP synthase dimer TM residues within 2 nm of the mid-surface | 98.8% (builds R and S) |
| Placed / requested copies | 99.3% of the request, 99.1% of the inventory; 99.5% outside mitochondria, above the 99% acceptance target (build S). Lowest: mitochondrial inner membrane 87.4%, nucleolus 92.7% of its request (84.1% of its inventory after thinning by 0.907), LD cytosolic face 95.0%, ER cytosolic face 95.3%, nuclear-envelope periphery 96.8%, Golgi cytosolic face 97.7% (of a capped request); matrix 98.5%, intermembrane space 99.0%. Build R: 97.5% overall, 98.9% outside mitochondria (nucleolus 80%, Golgi face 84%, matrix 75%, intermembrane space 87%) |
| Placement checks (\`validation/placement_checks.py\`, since build S) | placed 99.1% of the inventory by copies and 98.1% by sphere volume (build R 97.5% and 96.0%); species under half placed 10, under 80% 60 (build R 76 and 238); tile seams: sphere volume fraction within 5 nm of a tile face 0.41 / 0.40 (first / last pass) against 0.41 inside (build R 0.49 / 0.20 against 0.40); membrane-face copies outside their shell 0% (build R 4–17%) |
| Inventory protein mass vs proteome | within 0.4% before packing |
| Protein scale against counted structures (\`validation/anchors.py\`; check) | At 3.5 pg: core histones 0.90 of what the nucleosomes need (H3 0.70, H2A 0.77, H4 1.04, H2B 1.31); scaffold and FG nucleoporins a median 1.36 times what 93 pores need (10th–90th percentile 0.99–2.15); ribosomes 1.06–1.13 times the 187,000–200,000 counted from rRNA. Each anchor alone would put the protein total at 3.9, 2.6 and 3.1–3.3 pg. The outer kinetochore (Ndc80 complex 6–10 times its 8 per kinetochore) is not an anchor: these proteins have soluble pools |
| Organisation layer (section 8.4) | contact partners in contact 8.4% vs 1.3% (abundance-matched control); loose partners within 15 nm 19.9% vs 14.0%; tile-density correlation vs network similarity, Spearman 0.58; total crowding CV 0.040 (random repack 0.041) (an earlier build; not re-measured) |

### 9.1 Validation against independent data (2026-11-11)

These checks use data that never entered the model, with held-out splits fixed before any tuning (\`validation/splits.py\`: proteins by a hash of the ORF name, half dev and half test; diffusion measurements by a hash of the paper's DOI). Changes were tuned on dev only; test numbers were computed once per finished change. Code and outputs: \`validation/\`, \`work/validation/\`.

**Localization** (\`validation/localization.py\`). The model's compartment shares are compared with genome-wide localization that does not come from the C-terminal GFP collection behind CYCLoPs and Huh *et al.* (2003): the N-terminally tagged SWAT library under its native promoter and under the NOP1 promoter (Weill *et al.* 2018) and organelle spatial proteomics (hyperLOPIT; Nightingale *et al.* 2019). Classes: cytosol, nucleus, nucleolus, ER and nuclear envelope, mitochondrion, Golgi/endosome, vacuole, cell periphery, peroxisome, lipid droplet; SWAT "punctate" counts as Golgi/endosome, peroxisome or lipid droplet. A protein agrees when the class holding most of its copies is among the reference calls (top-1).

| Reference (test half) | Proteins | Top-1 before | Top-1 after the nucleus rule (4.2) | Macro-F1 after |
|---|---|---|---|---|
| SWAT, native promoter | 1,715 | 64.0% | 66.8% | 0.48 |
| SWAT, NOP1 promoter | 2,614 | 72.6% | 73.8% | 0.50 |
| hyperLOPIT (assigned proteins only) | 468 | 84.6% | 84.6% | 0.75 |

- **Placed copies.** The table scores the inventory's compartment shares. Scored on the copies actually placed in the rebuilt cell, test top-1 is 66.6% (SWAT native), 74.2% (SWAT NOP1) and 84.8% (hyperLOPIT) in build S (build R 67.4, 74.6 and 84.8%; 67.0, 74.1 and 84.5% when the table was made): packing shortfalls, the thinning of the nucleolus and the shell caps (section 8.2) hardly change the comparison.
- **Ceiling.** The references disagree with each other: SWAT (native) and CYCLoPs share a class for only 74% of the proteins both cover. Agreement with CYCLoPs, an input, is a fit (78% top-1), not a test.
- **Where the errors are** (dev half): proteins with measured shares agree far more often (73%) than those placed from annotation alone (55%), and the 625 proteins without any annotation, defaulted to the cytosol, agree least (40%). Peroxisome and lipid-droplet calls are rarely matched (F1 0.10 and 0.03), Golgi/endosome poorly (0.28).
- **Concentration instead of copies.** Scoring by the compartment of highest concentration (shares divided by volumes) agreed worse with every reference (it over-calls the small compartments), so copies are the primary measure.
- **Not adopted:** membrane association of soluble proteins from lipid anchors (myristoyl, palmitoyl, prenyl; UniProt) or membrane-binding domains (PH, PX, FYVE, C2, BAR, ENTH; InterPro, and lipid-array hits of Gallego *et al.* 2010). A bound fraction for domain proteins made dev agreement worse: most of them image as cytosolic. Lipid anchors helped dev by 0.1 point and not test, so the rule is switched off (\`membrane_anchor.py\`, \`BOUND_F = 0\`).

**Nucleus against cytoplasm from raw images** (2026-11-12; \`validation/image_nc.py\`; the first test on pixels rather than on published calls). The C-terminal SWAT library re-imaged for YeastRGB (Dubreuil *et al.* 2019; library: Meurer *et al.* 2018): mNeonGreen fusions under the native promoter, every strain with the envelope marker Nup49-mCherry, 6,145 strains with up to 225 single-cell crops each (9.7 GB of 16-bit images). Per cell the nucleus is segmented from the marker and the nuclear/cytoplasmic intensity ratio measured; the segmentation and thresholds were fixed on four strains before the comparison. Blur compresses the measured ratios, so ranks are compared.
- **Result (test half, 531 proteins, build S):** Spearman 0.81 between the model's and the measured nuclear/cytoplasmic concentration ratio; proteins the model puts in the nucleus are separated from those it puts in the cytosol with AUROC 0.98. Build R: 0.80 and 0.98 (538 proteins); build Q: 0.81 and 0.98 (562 proteins).
- **How independent:** the images are not an input. The nucleus rule (section 4.2) was tuned with CYCLoPs scores of the other half of the proteins, also C-terminal fusions, so tag artefacts shared by C-terminal fusions would not show.
- **Abundance (check):** the mean cellular intensity of 4,886 strains against the model's copies per cell, Spearman 0.74 (0.81 over the brightest 70%).

**Copy numbers** (\`validation/copy_numbers.py\`; check; the same in builds R and S, which share the copy numbers). At 3.5 pg the model's copies are a median 1.23 times those of the mass-spectrometric proteome of Das *et al.* 2026 in glucose (Spearman 0.81; 46% within 2-fold; 1.75 times at 5.0 pg) and 0.37 times the unified data set of Ho *et al.* 2018 (Spearman 0.93; 0.53 at 5.0 pg). The two data sets differ from each other by a median factor of 3.3, and both are culture averages. The check reads the per-protein feature table, which had not been regenerated since 2026-11-11; the build script now refreshes it (\`scripts/run_all.sh\`).

**Diffusion** (\`validation/diffusion.py\`): a literature table of 93 measurements from 25 papers, every value quoted from the paper (\`pipeline/data/validation/diffusion_measurements.tsv\`); 35 in-vivo values of soluble molecules in unstressed cells form the benchmark (free GFP and tandem GFPs, tagged cytoplasmic and nuclear proteins, 40 nm GEM particles). Free GFP alone ranges from 8.7 to 21 µm²/s between laboratories and methods, so the benchmark has absolute tracer anchors and, separately, within-paper ratios of a molecule to the same paper's free tracer. Section 14.5 describes the model; on the held-out papers:

| Model | Ratio RMSE (log10) | Absolute RMSE (log10) |
|---|---|---|
| Stokes–Einstein × one constant | 0.20 | 0.48 |
| + size-dependent crowding from Brownian dynamics of packed boxes | **0.09** | **0.21** |
| + binding propensity (RNA-interactome hits, pool partners) | 0.17 (build R 0.18) | 0.25 |

Build S, 15 absolute values and 5 ratios in the test split, which since 2026-11-12 includes the 40 nm GEM values (section 14.5; they are why the size-blind null is worse than the 0.31 quoted before). Build S gives the same values to two decimals as build R (ratio 0.09 against 0.09, absolute 0.21 against 0.21). The crowding term comes from Brownian dynamics of boxes cut from an earlier build; the cytosol's density is unchanged (170 mg/mL in build S, 168 in build R, 170 in the build the boxes came from), so it was not re-run. Two free-histone rows (H2B) use the H2B–H2A.Z dimer of the library as their species since build R, because the non-native H4–H2B pair they used before is no longer built (section 5.2) and the library has no free H2A–H2B dimer.

**Condensates** (\`pipeline/condensates.py\`, IMPROVEMENTS H5). Can membrane-less compartments emerge from generic interaction data at the model's concentrations? Interaction energies came only from the AlphaFold 3 pool scores (PAE converted to the likelihood ratio of a known IntAct interaction, including weak pairs at PAE 8–25 Å, which are still 1.5–10× enriched; pairs within one stable complex excluded), from RNA-interactome capture (Beckmann *et al.* 2015; Matia-González *et al.* 2015; Mitchell *et al.* 2013) and from IDR sticker residues. A multicomponent second-virial free energy over all species gives the stability of the mixture and, by linear response, the enrichment of each protein in an RNA-rich phase. The weights were chosen on the nucleus only (nucleolar members, CD-CODE: AUROC 0.68) and then frozen. In the cytosol, CD-CODE P-body members were enriched (AUROC 0.73, against 0.57 for permuted RNA-binding labels), but:
- the pool network added nothing (the same 0.73 with RNA binding alone and with a permuted network);
- the negative controls failed: stress-granule-only members (0.72) and the hypoxia-induced glycolytic G-body enzymes (0.91) were enriched just as much in log phase;
- the cytosol was stable at the calibrated interaction scale, in log phase and after an *in-silico* stress that releases all mRNA.

Re-run on build S as part of the validation run (\`validate_all.sh\`): P-body members AUROC 0.74, stress-granule-only members 0.70, G-body enzymes 0.89, the same with RNA binding alone (0.74) and with a permuted pool network (0.74), 0.59–0.64 with permuted RNA-binding labels; the cytosol stays stable in log phase and after the in-silico stress, and the stress changes none of the enrichments. So this shows that RNA-binding proteins co-enrich with RNA, not that P-bodies emerge. The model places no condensates; specific P-body assembly needs information that generic pair scores do not carry (valence, the decapping/Lsm network, regulation).

## 10. Densities, mass and volume fractions

### 10.1 Macromolecular densities

Densities are the mass of placed macromolecules per compartment volume. They count protein plus the RNA of ribosomes, pre-ribosomes, tRNAs and complexes, the DNA of nucleosomes and active rRNA genes, nascent pre-rRNA, and mRNA. Lipids, polysaccharides and small molecules are not counted.

| Compartment | Total, build S (mg/mL) | Build R (mg/mL) | Build P, 5.0 pg in 42 µm³ (mg/mL) | Protein / RNA+DNA split (earlier build, mg/mL) |
|---|---|---|---|---|
| Cytosol | 170 | 168 | 170 | 129 / 36 |
| Nucleoplasm | 111 | 110 | 106 | 99 / 7 |
| Nucleolus | 167 | 164 | 166 (136 before the larger attempt budget) | 108 / 28 |
| Mitochondrial matrix | 178 | 150 | 132 | 101 / 6 |
| Mitochondrial intermembrane space | 25 | 14 | 16 | 21 / 0 |
| ER lumen | 21 | 21 | 25 (44 before the ER topology fix moved 181,000 cytosolic-face copies out) | 44 / 0 |
| Golgi / endosome lumen | 2 | 2 | 3 | 2 / 0 |
| Peroxisomes | 97 | 97 | 120 | 120 / 0 |
| Vacuole lumen | 1 | 1 | 3 | 3 / 0 |

Builds R and S kept the densities of the large compartments by design: the cell volume was chosen so that the cytosol stays near 170 mg/mL (section 3.1). The vacuole lumen is more dilute than before build R because the vacuole is 17% of the cell instead of 9%. In build S the matrix and the intermembrane space are denser: the matrix now holds 98.5% of its proteins (finer packing shapes, section 8.2), and the intermembrane space is smaller (0.18 against 0.27 µm³) and holds more of its proteins (99.0% against 87%).

The protein / nucleic-acid split and the dry volume fractions (cytosol 11.4%, nucleoplasm 7.6%) were computed on an earlier build and are not recomputed by s4.

- **Organelle lumens:** since the 2026-11-11 rule (section 4.4) the vacuole and Golgi/endosome lumens hold only proteins that can reach them (signal peptide or Cvt cargo); their coat and peripheral proteins are counted in the cytosol, whose shells on the organelle faces they occupy. The low lumenal densities do not include glycans, cargo in transit or vacuolar polyphosphate.
- **Dry volume fraction:** computed with partial specific volumes of 0.73 mL/g for protein and 0.55 mL/g for nucleic acid.
- **Hydrated volume fraction:** with a hydration shell of about 0.3–0.4 g water per g, the hydrated volume fraction is roughly 1.5× the dry one: about 17% in the cytosol.
- **Protein and RNA/DNA split:** protein is the placed species' mass minus their nucleotides (330 Da each); RNA/DNA is the rest of the validated total, which also holds the mRNA and other threads. Where volume rounding made the remainder slightly negative it is shown as 0.
- **Excluded volume:** the volume excluded to a diffusing particle is larger still, and grows with the particle's size. It can be computed from the coordinates as accessible volume versus probe radius; this has not been done yet.

### 10.2 Whole-cell totals

| Quantity | Value |
|---|---|
| Protein in the inventory | 3.5 pg (5.0 pg through build Q) |
| Protein placed | about 2.9 pg (build R; not re-measured; 4.2 of 5.0 pg before) |
| RNA and DNA placed | about 0.8 pg (build R; not re-measured; 1.0 pg before): rRNA of ribosomes 0.64, tRNA 0.10, pre-ribosomes about 0.05, DNA 0.013 |
| mRNA | 0.011 pg (19.2 million nucleotides; builds R and S) |

The build R values of placed protein and RNA are estimates from the compartment densities and volumes above and the placed counts of section 2.2, not a separate tally. Build S places more of the matrix (section 8.3) and fewer nucleolar copies (thinning) and fewer mtDNA (1.27 instead of 3.0 Mb); the totals were not re-estimated.

Placed protein is lower than the inventory for four reasons:
- residues not resolved in PDB entries;
- AlphaFold residues trimmed for pLDDT < 50 (3.4% of mass);
- copies that did not fit (mostly mitochondria, nucleolus and Golgi);
- dropped and structure-less copies (including surplus histone H2B and nucleoporins beyond the pores' need, sections 6.1 and 11).

With full-length sequences, cytosolic protein would be about 150 mg/mL.

### 10.3 Comparison with dry-mass measurements

Measured cytoplasmic dry-mass densities are around 200 mg/mL. These are rough estimates of all components, averaged over the cytoplasm (cytosol plus organelles, excluding the vacuole):

| Component | Estimate (mg/mL) | Basis |
|---|---|---|
| Macromolecules as placed | about 160 | This model (build S: cytosol 170, nucleoplasm 111, nucleolus 167, matrix 178; build R 168, 110, 164, 150; volume-weighted over the cytoplasm without the vacuole, 158 in build S and 155 in build R) |
| Full-length protein correction | about +20 | Section 10.2 |
| Lipids of internal membranes and lipid droplets | about 30–35 (build R; not re-estimated) | Build R: about 160 µm² of internal membrane plus 0.12 µm³ of lipid droplets, about 0.8 pg of lipid in 22.7 µm³ of cytoplasm (20–25 mg/mL in the larger cell, which had nearly the same membrane area). Build S has 115 µm² of internal membrane (section 3.2), so this term is smaller |
| Small molecules and ions | probably 30–40 | Free amino acids, nucleotides, glutathione, sugar phosphates, phosphate and K⁺, a few hundred mM in total |
| Glycogen and trehalose | low in log phase on glucose | Can exceed 10% of dry weight after the diauxic shift |

Together these give about 240–255 mg/mL in build R, above the measurements (225–240 mg/mL in the larger cell of builds J to Q); in build S the lipid term is smaller, which was not re-estimated. A newer measurement by refractive-index tomography of wild-type cells in glucose gives about 205–225 mg/mL with the vacuole excluded (Kunzi *et al.* 2026, a preprint; read from their Fig. 1B; check added 2026-11-12). The same study's whole-cell protein concentration, 71 mg/mL by Bradford assay normalised to the Coulter volume, is closer than before but still lower than the model's 103 mg/mL (3.5 pg in 34 µm³; 119 mg/mL at 5.0 pg in 42 µm³). **An unresolved conflict:** the dry-mass measurements would be matched by the same 3.5 pg in a cell of about 36–38 µm³, while the ribosome concentration measured by cryo-ET (23 µM; the model has 17 µM in the cytosol) asks for a smaller cytosol, about 26 µm³ of cell at this composition. Builds R and S sit between the two (build S: 17.1 µM of cytosolic ribosomes); the cell volume was set by keeping the cytosolic density, not by either measurement. The small-molecule estimate is the least certain term.

## 11. Known limitations

| Area | Limitation |
|---|---|
| Physics | Static packing (random sequential addition plus an overlap-removal squeeze): no physical interactions or dynamics and random orientations. The optional motion is illustrative (section 7). |
| Replicates | Build S is one build at master seed 0. The build is deterministic (section 13), but replicate builds at other seeds have not been made, so no number in this document has a build-to-build spread yet. Builds R and S share the copy numbers and the mRNA draw and still differ, for example, by 0.35 against 0.45 in the mitochondrial mRNA test and 1.9× against 3.3× in balanced telomere clustering (sections 6.3, 6.5). |
| Organisation | Partner proximity and network co-enrichment (section 8.4) follow assumed strengths, not measurements. Partners are near each other, not bound in specific orientations. |
| Abundances | About 2-fold per-protein uncertainty. The proteome total is set to 3.5 pg for one unbudded G1 cell, between the stoichiometric anchors (section 2.1; bulk measurements of cultures give 4–6 pg); tRNA and mRNA counts are culture means times 0.7, an assumption. 73 proteins of other genotypes, mating type or conditions are set to zero by a curated list, which is incomplete. tRNA species are split by gene copies, and their structures are predictions without modified nucleosides. Relative abundances do not match counted structures exactly: histone H3 and H2A are 0.70 and 0.77 of what the nucleosomes need (the nucleosome count is fixed, so the missing copies are simply not booked), H2B is 1.3 times (since build S the surplus is packed as free H2B; in build R it had no library entry and was left out), and nucleoporins beyond 1.25 times the pores' need are not placed. |
| Ribosomal proteins | Ribosomal proteins are 15.7% of the model's protein mass, against 29–32% in three proteomics datasets of fast-growing yeast (Metzl-Raz *et al.* 2017, Björkeroth *et al.* 2020, Xia *et al.* 2022; none is BY4741 in YPD), 15–19% in the mass-spectrometry studies collected by Ho *et al.* 2018 and 2–7% in its GFP and TAP studies. At 3.5 pg of protein, 30% would be about 460,000 ribosomes, which no count supports. The counts themselves ask for somewhat more than the model's 211,500 mature ribosomes (11,500 per µm³ of cytosol): about 14,000 per µm³ of cytoplasm by cryo-ET (Delarue *et al.* 2018) and 14,517 ± 3,193 per µm³ of cytosol (Goetz *et al.* 2025, preprint), and about 220,000 in a G1 cell by serial-section EM (Yamaguchi *et al.* 2011, abstract only). The copy-number round (IMPROVEMENTS plan U) sets the ribosomes from the cryo-ET density. |
| Geometry | A composite of two segmented cells, scaled to 34 µm³. The vacuole is one ellipsoid with the measured shape statistics of G1 cells, not a segmented shape. Organelle volume fractions are from different studies and strains (section 3.1). The cell volume follows from the protein total and an assumed cytosolic density; dry-mass measurements favour a somewhat larger cell and the cryo-ET ribosome concentration a smaller one (section 10.3). Procedural mitochondria (two tubular networks in build S) built to measured dimensions of glucose-grown cells but not to any segmented cell; nucleoid positions relative to cristae are by analogy with mammalian cells. No bud, so no cell-cycle-dependent structures. |
| Mitochondria and small compartments | Mitochondria are 3% of the cell, against 4.77% measured in glucose (Di Bartolomeo *et al.* 2020, another strain). Matrix 98.5% filled (178 mg/mL; build R 75%), intermembrane space 99.0% (build R 87%). The inner membrane, 13.1 µm² with short cristae, holds 87.4% of its proteins (carriers 37–42%; section 6.7), and only 1,940 ATP synthase dimers fit on the crista rims. The nucleolus is thinned to 0.907 of its inventory before packing and holds 84.1% of it (build R 80%). The Golgi cytosolic face keeps 0.61 of its requested copies, the rest is packed in the cytosol (section 4.4). Mitochondrially encoded subunits are inferred from their complexes. |
| Packing artefacts | Containment: 6.7 × 10⁻⁵ of sampled beads lie outside their compartment in build S, against 4.1 × 10⁻⁵ in build R (round S candidate build 1: 9.1 × 10⁻⁵); the rise has not been analysed. Cytosolic proteins are depleted to 0.16 of their density in the 4–20 nm next to Golgi and endosome membranes, where the membrane-face shell is (section 8.3). |
| Disorder | AlphaFold residues with pLDDT < 50 are deleted rather than represented, except the FG domains of 4 Nups. The simulation restores their effect on the hydrodynamic radius only (section 14.2). |
| Condensates | None are placed. A test of whether P-bodies emerge from generic interaction data was negative for specificity (section 9.1). |
| Localization | Agrees with independent N-terminal SWAT calls for 67–74% of proteins and with spatial proteomics for 85% (section 9.1); the references agree with each other for about 74%. Bulk steady state. Measured shares (CYCLoPs) only between compartment groups and only where imaging and annotation agree; within groups, and for proteins without imaging data, dual localizations are weighted by rank. 625 proteins default to the cytosol. |
| Nucleolus | Active rRNA genes are self-avoiding walks, not the real tandem array; nascent transcripts are free threads rather than compact RNPs and fit only partly (57%). The number of active genes is per genome, so the smaller nucleolus of builds R and S (0.47 µm³) is more crowded. One structure per pre-ribosome stage. No 5S rDNA (Pol III) or snoRNA genes. tRNA genes are not gathered at the nucleolus (section 6.3; the FISH test fails, 23% against 52% in build S). |
| Cytoskeleton | No septins, fimbrin or formins. The spindle-pole-body core has its measured stoichiometry but no measured copy number (section 6.2). Actin partition and filament lengths are assumptions. 0.25–0.5% of the beads of the built actin structures overlap a neighbouring built component (section 9). Only 2 of 4 intended non-kinetochore nuclear microtubules and 15 of 16 kinetochore microtubules are built (section 6.2). |
| Plasma membrane | Eisosome coats hold only 7.6% of Pil1/Lsp1; only Sur7 and Nce102 are concentrated in the furrows. The endocytic coat is partial (no clathrin) and randomly oriented, and invagination depths are assumed uniform. |
| Endocytic sites | The sites lack the ribosome exclusion zone around the actin network (181 ± 41 nm across, reaching about 100 nm beyond the tip; Picco *et al.* 2018, Buser & Drubin 2013). Invaginations are shallow: depths are drawn up to 100 nm, and the deepest in round S candidate build 1 was 92 nm, against about 100 nm on average and up to 140 nm at scission by tomography (Kukulski *et al.* 2012). The tube tips are too wide: 40 nm across, against a measured tip radius of 11.6 ± 2.6 nm, about 23 nm across (Kukulski *et al.* 2012 data in Picco *et al.* 2018). Actin per site is half the measured peak (section 6.4). |
| Nuclear pore basket | Not built. The only deposited basket model (8ZZC, integrative) resolves about a quarter of each Mlp1/Mlp2 as coarse beads, and their AlphaFold models are compact rather than the long coiled coils of the basket filaments. Mlp1, Mlp2 and Pml39 are packed free next to the nuclear envelope. |
| Contact sites | Only the nucleus–vacuole junction. The ER here is cortical and perinuclear, so ER–mitochondria (ERMES) and other cytoplasmic-ER contacts have no place. |
| Lipids and wall | Lipid compositions from isolated organelle membranes (1991), partly from non-glucose cultures; only the plasma membrane has leaflet asymmetry. Procedural cell wall. One bud scar (rim height and chitin share assumed), no birth scar; lateral-wall chitin not drawn. |
| mRNA | The mRNA between ribosomes is drawn shorter than its nucleotides, and the tails carry the rest (section 6.5). Ribosomes carry Stm1 (the idle-ribosome structure 4V88) and no tRNAs or elongation factors. Too few mRNAs are free of ribosomes: 570 of 18,136 (3%), while polysome profiling finds a median of 72% of a gene's mRNA on ribosomes (Arava *et al.* 2003, Supplementary Table 3; compared on build R, 2026-11-13). With 18,100 mRNAs a typical gene has few copies, so per-gene localization tests are noisy (mitochondrial mRNA test 0.45 / 0.66 in build S, 0.35 / 0.61 in build R with the same mRNAs; section 6.5). In a simulated box, releasing polysome ribosomes speeds 40 nm particles by 2% (measured 40–50%; section 14). |
| Missing classes | Metabolites, ions, water, glycogen, trehalose; polyphosphate only procedurally (amount assumed). |

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
| AlphaFold 3 pool predictions (in-house, unpublished) | Pair scores for partner placement | not public |
| STRING v12 | Functional and physical associations, for network-propagation profiles and loose partners | CC-BY 4.0 |
| CYCLoPs (Koh *et al.* 2015) | LOC-scores of the three wild-type GFP-collection screens (thecellvision.org/cyclops), for compartment shares | |
| Vögtle *et al.* 2017 | Submitochondrial localization of 810 proteins (Supplementary Data 5) | CC-BY 4.0 |
| Morgenstern *et al.* 2017 | High-confidence mitochondrial proteome, enrichment classes, submitochondrial classes (Table S1) | |
| Zinser *et al.* 1991 | Phospholipid classes and ergosterol of isolated organelle membranes (Tables 4, 5) | |
| EMPIAR-11416 (cryo-plasma FIB/SEM of yeast), segmented into cell, vacuole, lipid-droplet and nucleus meshes | The cell outline and lipid droplets (cell 08), the position of the vacuole (cell 08's vacuoles; their shapes up to build R) and the nucleus (cell 05) of the model cell (section 3.3) | CC0 |
| Chen *et al.* 2025, Supplementary file S8 (soft X-ray tomography of unbudded G1 BY4741 in YPD) | Number, sphericity and ellipsoid axis ratios of the vacuole (section 3.3) | |
| EMPIAR-12534 (in-cell cryo-ET of BY4741 in YPD; Chang *et al.* 2025), 50 segmented membrane meshes | Outer-membrane diameter, membrane spacing, crista area, density and depth, measured here on the YPD meshes (\`pipeline/data/validation/literature/mitochondria/\`; section 3.3) | |

### 12.2 Literature behind the parameters

| Reference | Used for |
|---|---|
| Carroll KM *et al.* (2011). *Mol. Cell. Proteomics* 10:M111.007633; Sherman F (2002). *Methods Enzymol.* 350:3–41 | Protein mass per cell (4–5 pg; 6 pg) |
| Ho B, Baryshnikova A, Brown GW (2018). *Cell Systems* 6:192–205 | Earlier protein total (molecules per cell) |
| Jorgensen P *et al.* (2002). *Science* 297:395–400 | Cell volume |
| Uchida M *et al.* (2011). *Yeast* 28:227–236 | Organelle volume fractions (soft X-ray tomography) |
| Dumoux M *et al.* (2023). *eLife* 12:e83623; EMPIAR-11416 | Cryo-plasma FIB/SEM volume of yeast; segmented cells used for nucleus, vacuole and lipid-droplet characteristics |
| Jorgensen P *et al.* (2007). *Mol. Biol. Cell* 18:3523–3532 | Nuclear-to-cell volume ratio |
| Yamaguchi M *et al.* (2011). *Microscopy* 60:321–335, doi:10.1093/jmicro/dfr052 | Serial-section EM comparison |
| Delarue M *et al.* (2018). *Cell* 174:338–349.e20, PMID 29937223 | Ribosome concentration by cryo-ET (about 14,000 per µm³ of cytoplasm) |
| Goetz SK *et al.* (2025). bioRxiv preprint, doi:10.64898/2025.12.22.696049 (Europe PMC PPR1230379; not peer reviewed) | Ribosome concentration by cryo-ET (14,517 ± 3,193 per µm³ of cytosol) |
| Orlean P (2012). *Genetics* 192:775–818 | Cell-wall thickness and composition |
| Wolf W *et al.* (2012). *PLoS ONE* 7:e39703 | Cortical ER coverage |
| Winey M *et al.* (1997). *Mol. Biol. Cell* 8:2119–2132 | NPC number |
| Ejsing CS *et al.* (2009). *PNAS* 106:2136–2141 | Lipid composition |
| Huh WK *et al.* (2003). *Nature* 425:686–691 | GFP localization |
| Kefauver JM *et al.* (2024). *Nature* 632:664–671 | Pil1/Lsp1 eisosome lattice on native membranes (8QB7, 8QB8) |
| Lynch-Day MA, Klionsky DJ (2010). *FEBS Lett.* 584:1359–1366; Yuga M *et al.* (2011). *J. Biol. Chem.* 286:13704–13713 | Cvt pathway; Cvt cargo (Ape1, Ams1, Ape4) reaches the vacuole lumen without a signal peptide |
| Laporte D *et al.* (2008). *J. Cell Biol.* 181:737–745 | 26S proteasome nuclear in growing cells |
| Bystricky K *et al.* (2004). *PNAS* 101:16495–16500 | Chromatin compaction and persistence in living yeast (fibre) |
| Arbona JM *et al.* (2017). *Genome Biol.* 18:81 | Chromatin compaction and persistence from imaging and Hi-C (fibre radius choice) |
| Pelechano V, Chávez S, Pérez-Ortín JE (2010). *PLoS One* 5:e15442 | mRNA molecules per gene (Table S1); engaged Pol II |
| Arava Y *et al.* (2003). *PNAS* 100:3889–3894 | Ribosome association of yeast mRNAs (test, not used; ribosome-free mRNAs, section 11) |
| Chereji RV, Ramachandran S, Bryson TD, Henikoff S (2018). *Genome Biol.* 19:19 | H3 chemical-cleavage dyad maps: measured nucleosome positions of the fibre |
| Rossi MJ *et al.* (2021). *Nature* 592:309–314 | Promoter nucleosome-depleted regions (NDR filter); ChIP-exo sites of chromatin proteins |
| Oberbeckmann E *et al.* (2019). *Genome Res.* 29:1996–2009 | Absolute nucleosome occupancy and count (held-out test) |
| Weiner A *et al.* (2015). *Mol. Cell* 58:371–386 | MNase nucleosome positions, unstressed (held-out test) |
| Schep AN *et al.* (2015). *Genome Res.* 25:1757–1770 | ATAC-accessible regions (held-out test) |
| Brogaard K, Xi L, Wang JP, Widom J (2012). *Nature* 486:496–501 | Chemical nucleosome map (count cross-check; sacCer2, not used) |
| Eaton ML *et al.* (2010). *Genes Dev.* 24:748–753 | Nucleosome-free regions at origins (check of the dyad calls) |
| Siow CC *et al.* (2012). *Nucleic Acids Res.* 40:D682–D686 | OriDB confirmed replication origins (MCM double hexamers, ORC) |
| Joglekar AP *et al.* (2006). *Nat. Cell Biol.* 8:581–585 | Kinetochore protein copy numbers |
| Marcand S, Gilson E, Shore D (1997). *Science* 275:986–990 | Rap1 sites every 18 bp on ~300 bp telomeric tracts (17 Rap1 per end) |
| Swygert SG *et al.* (2014). *Nat. Commun.* 5:4751; Swygert SG *et al.* (2018). *PNAS* 115:12447–12452 | Two Sir3 and one to two Sir2–Sir4 per nucleosome on reconstituted arrays (upper bounds in cells) |
| Muller EG *et al.* (2005). *Mol. Biol. Cell* 16:3341–3352 | SPB core stoichiometry Spc42:Spc29:Spc110:Cmd1:Cnm67 = 2:1:1:1:1 (as used by Viswanath *et al.* 2017) |
| Haeusler RA *et al.* (2008). *Genes Dev.* 22:2204–2214 | tRNA-gene clusters at the nucleolus in wild type (test) |
| Ikeuchi K *et al.* (2019). *EMBO J.* 38:e100276 (PDB 6I7O); Tesina P *et al.* (2020). *EMBO J.* 39:e103365 (PDB 6T83) | Collided-disome pose (6I7O used; 6T83 within 15° and 2.9 nm) |
| Meydan S, Guydosh NR (2020). *Mol. Cell* 79:588–602 | Disome footprints of 57–63 nt: at most about 5 nt between collided ribosomes |
| Zhao T *et al.* (2021). *Genome Biol.* 22:16 | 5.8% of ribosomes in collided disomes in exponential yeast (test) |
| Chan PP, Lowe TM (2016). *Nucleic Acids Res.* 44:D184–D189 | GtRNAdb mature tRNA sequences (42 isoacceptors) |
| Percudani R, Pavesi A, Ottonello S (1997). *J. Mol. Biol.* 268:322–330 | tRNA levels follow tRNA gene copy number in growing yeast |
| Passaro S *et al.* (2025). bioRxiv 10.1101/2025.06.14.659707 | Boltz-2, used to predict the tRNA structures |
| Kim SJ *et al.* (2018). *Nature* 555:475–482 | Copies of each nucleoporin per pore (Supplementary Table 2A), for the nucleoporin bookkeeping and the pore anchor |
| von der Haar T (2008). *BMC Syst. Biol.* 2:87 | 187,000 ± 56,000 ribosomes per cell (check of the protein scale) |
| Brothers M, Rine J (2022). *eLife* 11:e75653 | Sir3 ChIP-seq of cycling cells: which nucleosomes are Sir-bound |
| Goodnight D, Rine J (2020). *eLife* 9:e58910 | Sir4 ChIP-seq with untagged control (held-out test of the Sir-bound nucleosomes) |
| Teytelman L, Thurtle DM, Rine J, van Oudenaarden A (2013). *PNAS* 110:18602–18607 | Highly expressed genes give false ChIP signal (why Sir-bound nucleosomes are limited to chromosome ends) |
| Zhong J *et al.* (2016). *Genome Res.* 26:351–364 | Per-base DNase I cuts in log-phase yeast (held-out test of nucleosome positions) |
| Zelin E *et al.* (2012). *Mol. Cell* 48:459–470 | DNase-hypersensitive sites of wild-type BY4741 (held-out test) |
| Dubreuil B *et al.* (2019). *Nucleic Acids Res.* 47:D1245–D1249; Meurer M *et al.* (2018). *Nat. Methods* 15:598–600 | YeastRGB images of the C-SWAT mNeonGreen library with a nuclear-envelope marker (nucleus/cytoplasm test on raw images) |
| Das H *et al.* (2026). *Histochem. Cell Biol.* 164:8 | Absolute proteome of glucose-grown yeast (copy-number check) |
| Gade VR *et al.* (2025). *Cell Rep.* 44:116204 (preprint 10.1101/2024.11.14.623620); Kunzi M *et al.* (2026), bioRxiv 10.64898/2026.08.29.747980 | 40 nm GEM diffusion in budding yeast (test); cytoplasmic dry mass and protein concentration (check) |
| Ohsawa S *et al.* (2026). bioRxiv 10.64898/2026.09.03.749057 | Nucleus 5.3–6.3% of the cell in control strains (check) |
| Pelechano V, Wei W, Steinmetz LM (2013). *Nature* 497:127–131 | Per-gene UTR lengths (TIF-seq, YPD) |
| Subtelny AO *et al.* (2014). *Nature* 508:66–71 | Per-gene poly(A) tail lengths (PAL-seq) and translational efficiency |
| Williams CC, Jan CH, Weissman JS (2014). *Science* 346:748–751 | mRNAs translated at the mitochondrial surface (−CHX calls) |
| Jan CH, Williams CC, Weissman JS (2014). *Science* 346:1257521 | ER-enriched secretome genes (Sec63 proximity ribosome profiling) |
| Chartron JW, Hunt KCL, Frydman J (2016). *Nature* 536:224–228 | Membrane-polysome enrichment per gene (held-out test of ER mRNAs) |
| Bullitt E, Rout MP, Kilmartin JV, Akey CW (1997). *Cell* 89:1077–1086; Viswanath S *et al.* (2017). *Mol. Biol. Cell* 28:3298–3314 | Spindle pole body layers (cryo-EM; integrative model of the core) |
| Thompson M, Haeusler RA, Good PD, Engelke DR (2003). *Science* 302:1399–1401 | tRNA genes cluster at the nucleolus (check; Hi-C does not show it) |
| Saint-Georges Y *et al.* (2008). *PLoS One* 3:e2293 | Mitochondrion-localized mRNAs (held-out test) |
| Chen J *et al.* (2025). *Mol. Biol. Cell* 36:ar132; Di Bartolomeo F *et al.* (2020). *PNAS* 117:7524–7535 | Organelle volumes by soft X-ray tomography (vacuole fraction; since build S also the single near-spherical vacuole, Supplementary file S8); mitochondrial volume in glucose (checks) |
| Barbosa AD *et al.* (2015). *Mol. Biol. Cell* 26:3641–3657; Santana-Sosa S, Matos-Perdomo E, Ayra-Plasencia J, Machín F (2023). *Int. J. Mol. Sci.* 24:9829; Yang CH, Lambie EJ, Hardin J, Craft J, Snyder M (1989). *Chromosoma* 98:123–128 | Nucleolus at the envelope next to the vacuole; nucleolus opposite the spindle pole body |
| Severs NJ, Jordan EG, Williamson DH (1976). *J. Ultrastruct. Res.* 54:374–387; Kohler V, Büttner S (2021). *Contact* 4:25152564211016608; Otto FB, Thumm M (2020). *Int. J. Mol. Sci.* 21:4506 | No nuclear pores in the nucleus–vacuole junction (Severs 1976 and Pan 2000, as cited by the two reviews) |
| Egner A, Jakobs S, Hell SW (2002). *PNAS* 99:3370–3375 | Mitochondrial tubule diameter in glucose (live 4Pi microscopy) |
| Chang YT *et al.* (2025). *J. Cell Biol.* 224:e202407110 (EMPIAR-12534) | Mitochondrial membrane spacing and cristae in BY4741 in YPD (cryo-ET; meshes measured here) |
| Stoldt S *et al.* (2019). *PNAS* 116:9853–9858; Klecker T, Westermann B (2021). *Open Biol.* 11:210238 | Crista shape and placement (short, sheet-like, on two opposing sides); few cristae in glucose; ATP synthase per mitochondrial mass in glucose |
| Osman C *et al.* (2015). *PNAS* 112:E947–E956; Göke A *et al.* (2020). *Mol. Biol. Cell* 31:527–545; Seel A *et al.* (2023). *Nat. Struct. Mol. Biol.* 30:1549–1560; Kukat C *et al.* (2015). *PNAS* 112:11288–11293 | Nucleoids per µm of tubule and their spacing; mtDNA copies per cell and per fL; separate networks per cell; nucleoids in crista-free matrix (mammalian) |
| Metzl-Raz E *et al.* (2017). *eLife* 6:e28034; Björkeroth J *et al.* (2020). *PNAS* 117:21804–21812; Xia J *et al.* (2022). *Nat. Commun.* 13:2819 | Ribosomal-protein share of the proteome in fast-growing yeast (limitation, section 11) |
| Picco A *et al.* (2018). *Mol. Biol. Cell* 29:1346–1358; Buser C, Drubin DG (2013). *Microsc. Microanal.* 19:381–392 | Endocytic tip radius and ribosome exclusion zone (limitation, section 11) |
| Strádalová V *et al.* (2009). *J. Cell Sci.* 122:2887–2894; Karotki L *et al.* (2011). *J. Cell Biol.* 195:889–902 | Eisosome furrow dimensions and number |
| Bharat TAM, Hoffmann PC, Kukulski W (2018). *Structure* 26:879–886 | In-situ eisosome coat diameter |
| Kukulski W *et al.* (2012). *Cell* 150:508–520; Kukulski W *et al.* (2016). *eLife* 5:e16036 | Endocytic invagination depth at scission, vesicle size, tip radius |
| Pillet F *et al.* (2014). *BMC Biol.* 12:6 | Bud-scar diameter and profile (AFM) |
| Gotta M *et al.* (1996). *J. Cell Biol.* 134:1349; Duan Z *et al.* (2010). *Nature* 465:363; Therizols P *et al.* (2010). *PNAS* 107:2025 | Telomere foci; Rabl layout and arm-length dependence |
| Pelechano V, Chávez S, Pérez-Ortín JE (2010). *PLoS One* 5:e15442 | Engaged Pol II per gene |
| West M *et al.* (2011). *J. Cell Biol.* 193:333 | Ribosome density on the cortical ER |
| Picco A *et al.* (2015). *eLife* 4:e04535 | Molecules per endocytic site; Rvs161/167 assembly on the invagination |
| Pan X *et al.* (2000). *Mol. Biol. Cell* 11:2445; Jeong H *et al.* (2017). *PNAS* 114:E4539 | Nucleus–vacuole junction (no pores there, section 6.1); Vac8–Nvj1 structure (5XJG) |
| Kisonaite M *et al.* (2023). *Nat. Struct. Mol. Biol.* 30:670 | RAC on the 80S ribosome (7Z3N) |
| Baba M *et al.* (1997). *J. Cell Biol.* 139:1687 | Cvt vesicle size |
| Bru S *et al.* (2025). *eLife* 14:RP108181; Christ JJ *et al.* (2020). *Biotechnol. Bioeng.* 117:2089 | Vacuolar polyphosphate store; polyP chain length |
| Klis FM *et al.* (2002). *FEMS Microbiol. Rev.* 26:239–256 | Cell-wall composition (chitin 1–2%) |
| Sickmann A *et al.* (2003), *PNAS*; Reinders J *et al.* (2006), *J. Proteome Res.* | Excluded high-throughput organelle calls |
| Dolgalev G *et al.* (2023). *Int. J. Mol. Sci.*, doi:10.3390/ijms24032081 | Protein density comparison |
| Kostinski S, Reuveni S (2021). *Phys. Rev. Research*, doi:10.1103/physrevresearch.3.013020 | Ribosome count comparison |
| Warner JR (1999). *Trends Biochem. Sci.* 24:437–440 | Ribosome production rate; about 200,000 ribosomes per cell |
| French SL, Osheim YN, Cioci F, Nomura M, Beyer AL (2003). *Mol. Cell. Biol.* 23:1558–1568 | Active rRNA genes, Pol I per gene, elongation rate (Miller spreads) |
| Dammann R *et al.* (1995). *Mol. Cell. Biol.* 15:5294–5303; Merz K *et al.* (2008). *Genes Dev.* 22:1190–1204 | Active rDNA is nucleosome-depleted (psoralen accessibility; Hmo1 chromatin) |
| Kos M, Tollervey D (2010). *Mol. Cell* 37:809–820 | Co-transcriptional pre-rRNA cleavage |
| Barandun J *et al.* (2017), Cruz VE *et al.* (2022), Kater L *et al.* (2020), Malyutin AG *et al.* (2017), Scaiola A *et al.* (2018); Neyer S *et al.* (2016), *Nature* 540:607 | Pre-ribosome structures; Pol I elongation complex (5M3F) |
| Maritan M *et al.* (2022). *J. Mol. Biol.* 434:167351 | Whole-cell modelling approach |
| Weill U *et al.* (2018). *Nat. Methods* 15:617–622 | N-terminal SWAT GFP localization (independent validation, 9.1) |
| Nightingale DJH, Oliver SG, Lilley KS (2019). *Curr. Opin. Chem. Biol.* 48:86–95 | hyperLOPIT spatial proteome of yeast (validation, 9.1) |
| Gallego O *et al.* (2010). *Mol. Syst. Biol.* 6:430 | Protein–lipid arrays (membrane-anchor test, 9.1) |
| Beckmann BM *et al.* (2015). *Nat. Commun.* 6:10127; Matia-González AM *et al.* (2015). *Nat. Struct. Mol. Biol.* 22:1027–1033; Mitchell SF *et al.* (2013). *Nat. Struct. Mol. Biol.* 20:127–133 | RNA-interactome capture (condensate and diffusion tests, 9.1) |
| Rostam N *et al.* (2023). *Nat. Methods* 20:673–676 | CD-CODE condensate membership (condensate test, 9.1) |

### 12.3 PDB entries used for assemblies and templates

| Assembly or template | PDB entries |
|---|---|
| Ribosome and polysomes | 4V88, 7ZW0 |
| Pre-ribosomes and rDNA transcription | 5WLC, 7R7A, 6YLH, 5T62, 6FAI; Pol I elongation complex 5M3F |
| ER translocons | 2WW9, 7KAH |
| Nuclear pore and nucleosome | 7N9F, 1ID3 |
| Microtubules | 8QV3, 5W3F |
| Actin and actin-binding proteins | 9GO5, 5JLF, 8UXW, 8F8Q |
| Mitochondrial inner membrane | 8F39, 4B2Q, 6HU9; mitoribosome 9HM0 |
| Eisosomes | 8QB7, 8QB8 |
| RAC on the ribosome; Pol II elongation complex; Cvt cargo; nucleus–vacuole junction | 7Z3N; 1Y1W; 5JGF (Ape1 dodecamer); 5XJG |
| Other curated complexes | 4CR2, 2UV8, 4V94, 4C2M, 1WCM, 5FJ8, 7VH5, 6O7V |

The other complexes and homo-oligomers are listed per species in the viewer's inspect panel: click any molecule.

## 13. Reproducibility and data access

- **Pipeline.** The model is built by a Python pipeline (numba-accelerated) in six stages: the base structure library (stage 0), inventory, geometry, packing, validation and encoding (\`scripts/run_all.sh\`; \`scripts/README.md\` lists the stages). Build S took 87 minutes from stage 0 to the encoded tiles on 40 workers of a 64-core cluster node, 55 of them for the tile packing (all six rounds, the last with the overlap-removal squeeze); the comparisons with data (\`scripts/validate_all.sh\`, which needs the GPUs) took another 17 minutes. The real shapes come from \`pipeline/cell_meshes.npz\` (in the repository). The organisation layer of the packing (section 8.4) needs STRING (public) and the pool pair scores (not public; the derived partner lists are reproducible only with them). \`work/interactions/pair_summary.parquet\` (the pool pairs read by \`features.py\`) is an input kept in \`work/\` without a generator in the repository.
- **Stage 0, the base structure library (since 2026-11-13; \`scripts/s0_library.py\`).** The curated complexes from their RCSB assemblies (membrane complexes oriented), one AlphaFold DB model per protein of the inventory (residues with pLDDT < 50 left out, transmembrane proteins oriented from the UniProt topology), the curated PDB homo-oligomers and tRNA-Phe (1EHZ). Until build R these entries were made only by cells of a historical notebook, so the library on disk was state that no script could rebuild. Rebuilt and compared entry by entry, 5,561 of 5,647 shared entries were identical; the differences were the 82 transmembrane proteins oriented before the side rule (section 5.3) and the missing H2B entries (section 2.3). Build S was built with stage 0 (7,091 library entries); the amphiphysin heterodimer is now made by its builder (section 6.4). The library path can be set with \`YC_LIB\`.
- **Isolated builds (\`scripts/build_isolated.sh NAME [GIT_REF] [FIRST_STAGE]\`).** A build of one commit runs from a frozen git worktree of that commit in its own run folder (\`../runs/NAME/\`: \`tree/\` the worktree, \`lib/\` its own copy of the structure library, \`work/\` its own output, seeded only with download caches, and \`build.log\`), so code can be edited while a build runs and builds never share output or library. \`YC_RESUME=1\` continues a run from a later stage after a fix (the build is then marked as mixed); \`YC_LIB_FROM=empty\` starts from an empty library. Build S is the run \`runs/S3_seed0\` of commit 9568da3.
- **Master seed (\`YC_SEED\`, \`config.seed\`).** Every builder has its own random stream, numbered by its historical seed (tiles: 1000 plus the tile index). At master seed 0 the stream number itself is the seed, which reproduces the reference build; any other \`YC_SEED\` shifts all streams together, giving an independent replicate of the cell. Streams that make the structure library, the encoding or validation samples do not depend on it. Builds run with \`PYTHONHASHSEED=0\`.
- **Build stamp (\`scripts/build_stamp.py\`, \`work/build.json\`).** At the start of a build: build id, commit, whether code under \`pipeline/\`, \`scripts/\` or \`validation/\` was uncommitted, seed, start time, host, workers, every module-level constant of \`pipeline/*.py\`, and size and modification time of the inputs. After the validation run, \`build_stamp.py check\` lists result files older than the build and writes \`work/validation/MANIFEST.json\` (file, build id, age); for build S all 30 result files are of this build. With several builds a day, numbers of an earlier build had reached this document twice (2026-11-12).
- **Determinism.** The build is deterministic: the same commit, master seed, inputs and structure library give the same cell, independent of the number of workers and of the hash seed. Checked on 2026-11-13: \`main\` (build R) rebuilt by \`build_isolated.sh\` with 56 instead of 64 workers was identical to build R in all 5,912 tiles (24,484,185 instances at the same positions), all 128 pre-pass arrays and the expected totals; the seed plumbing at master seed 0 reproduced build R in all 5,912 tiles as well. Replicates at other seeds have not been built yet (section 11).
- **Source code.** The pipeline and the viewer source are in \`source/\`; \`source/README.md\` explains how to rebuild. The historical notebook and the eight modules only it used (the 2D column model and its figures) were removed in round S.
- **This repository** holds the built viewer and the encoded data.
  - The tiles under \`data/tiles/\` store every instance as a species index, a position quantised to about 0.0055 nm and an orientation quantised to about 0.004°.
  - The residue-bead library is in \`data/library/\`.
- **Full-atom coordinates.** The viewer fetches full-atom models live from RCSB and AlphaFold DB for display in Mol*.

## 14. Simulating a box of the cell (experimental)

The viewer can cut a cube of 30–200 nm out of the packed cell and run Brownian dynamics (BD) of its molecules on the GPU. While the simulation runs, the whole-cell rendering and tile streaming pause. Only the box and the membranes inside it are drawn. The aim is to watch the packed cell move on physically grounded diffusion. It is not to predict specific biochemistry. The model has excluded volume, a nonspecific attraction, screened electrostatics, specific binding at the interfaces of predicted complexes, polysome mRNA links and mean-field hydrodynamics, but no flexibility and no reactions.

The approach follows BD models of bacterial cytoplasm built from structures (McGuffee & Elcock 2010; Ando & Skolnick 2010) and of concentrated protein solutions (Mereghetti & Wade 2012; Balbo *et al.* 2013). All-atom MD of a cytoplasm (Yu *et al.* 2016) needs supercomputers for 20 ns. Rigid-body BD reaches microseconds of a 100 nm box in a browser.

Approaches considered:

| Approach | Resolution and time scale | Fit for a box of this cell in a browser |
|---|---|---|
| All-atom MD with explicit water (Yu *et al.* 2016) | atoms; tens of ns for ~100 nm | far beyond a browser GPU |
| Martini 3 coarse-grained MD (Souza *et al.* 2021); a whole minimal cell has been built in Martini (Stevens *et al.* 2023) | about 4 heavy atoms per bead, explicit CG water; ns–µs | the export (section 14.7) is the starting point for it; too heavy to run interactively |
| Elastic-network CG with implicit solvent (OpenMM custom forces) | one bead per residue; ns–µs | used for the exported OpenMM script |
| Rigid-body Brownian dynamics with implicit solvent and mean-field hydrodynamics (McGuffee & Elcock 2010; Mereghetti & Wade 2012) | residues or sphere trees; µs | **the model in the viewer** |
| Particle-based reaction–diffusion: Smoldyn (Andrews 2017), ReaDDy 2 (Hoffmann *et al.* 2019) | point particles; ms–s | no molecular shape, so it loses what the packed model provides |
| Lattice reaction–diffusion master equation: Lattice Microbes (Roberts *et al.* 2013), used for a whole minimal cell (Thornburg *et al.* 2022) | ~10 nm voxels; minutes of cell time | a natural next level for the whole cell; needs kinetic parameters this model does not have |

### 14.1 What is simulated

- **Molecules.** Every packed instance whose centre lies in the box becomes a rigid body with its packed position and orientation.
- **Schematic model.** Contacts use the library sphere tree of each species (about 30 residues per sphere). The default step is 0.2 ns, which reaches microseconds.
- **Mid model.** Contacts use clusters of about 4 residues (2 nucleotides): k-means of each species' beads, per kind, seeded by runs of the Morton-ordered beads. A cluster's radius is its bead radius plus 0.8 times the rms distance of its members from the centre. Measured along rays toward the centres of 13 species, the cluster surface then lies 0.13–0.22 nm outside the residue surface on average (0.3 nm for the ribosome), so molecules touch with a small gap instead of visibly overlapping as in the schematic model. Only clusters holding an exposed residue are kept, plus core spheres (below). The default step is 0.2 ns.
- **Residue model.** Contacts use one sphere per residue bead: radius 0.36 nm for protein, 0.6 nm for nucleotide. Only solvent-exposed beads are used, about 60% of them. A bead is exposed if a probe bead in contact along one of 26 directions overlaps none of the molecule's other beads; buried beads cannot meet another molecule before its surface beads do. The default step is 0.1 ns, which reaches hundreds of nanoseconds.
- **Core spheres (mid and residue models).** With surface spheres only, molecules are hollow shells. In a test, a strong attraction (mid model, 1.5 kT) pulled them through each other: grid cells filled far beyond their capacity, neighbours were missed, and molecules moved faster than in water. Each molecule therefore also carries core spheres: its sphere tree shrunk by 1 nm, about 1 nm inside the residue surface. Cores meet only other cores, with overlap repulsion and nothing else (no charge, attraction or contacts). They stop molecules passing through each other without adding to the cost of surface contacts. Keeping every buried bead instead (solid bodies) did the same job at about 1.5× the cost.
- **Drawing.** All models draw the residue beads.
- **Held fixed, acting only as obstacles:**
  - structures built explicitly before packing and anchored to filaments, membranes or DNA (section 6: filaments, lattices, assemblies; instance flag bit 1 from the encoder). Polysome ribosomes are the exception: they move;
  - nucleosomes and nuclear pores.
- **Membrane proteins** (membrane regions) diffuse laterally in their membrane (section 14.2), or stay fixed if that option is off.
- **Box faces.** Three boundary conditions:

  | Boundary | Behaviour | Use |
  |---|---|---|
  | Periodic | Every molecule moves. Centres are wrapped into the box and contacts use the nearest periodic image, as in standard MD of bulk solutions. | Default for boxes without membranes, such as cytoplasm or nucleoplasm. |
  | Embedded | The box is closed, but molecules in its outer 8 nm keep moving while tied to their starting positions by a harmonic tether (1 kT/nm², about 1 nm rms per axis). Molecules near the faces still see crowded surroundings, and the cell outside holds them in place. | Default for boxes containing membranes, which do not continue across the faces. |
  | Closed walls | Molecule centres reflect at the faces; parts of a molecule may stick out. | Isolated boxes. Beware the depletion layer next to the walls. |

  The first version froze an 8 nm shell completely, which immobilised about 60% of the molecules of a 60 nm box.
- **Membranes.** Membranes are walls. A grid over the box stores the signed distance to the nearest bilayer midplane, evaluated once from the same analytic membrane distance functions as the renderer. Soluble molecules are kept at least 2 nm (half a bilayer) plus their sphere radius from every bilayer, on their own side.
- **Polysomes.** Ribosomes of one polysome stay linked by their mRNA (section 14.3). The packer placed each next ribosome so that its mRNA entry lies 3–12 nm beyond the previous one's exit, within 0.7 rad of the exit direction, with its entry channel pointing back along the linker. Pairs of ribosomes in the box that meet all three conditions to the tile quantisation are linked. Links to ribosomes outside the box are dropped. A 200 nm cytoplasm box holds 125 ribosomes with 12 links joining 21 of them.
- **Not included as bodies.** mRNA threads, chromatin fibres beyond the nucleosomes, cell-wall polysaccharides and lipids are absent.

### 14.2 Diffusion coefficients

**Membrane proteins** diffuse in the tangent plane of their membrane. The plane comes from the gradient of the membrane distance grid at the centre. The protein is returned to its starting depth after every step, and it rotates only about the membrane normal. Its coefficients follow Saffman & Delbrück (1975) for a cylinder of radius *a* in a membrane of viscosity μ_m and thickness h between water of viscosity μ:

\`\`\`
D = kT / (4π μ_m h) · (ln(μ_m h / (μ a)) − 0.5772)
D_r = kT / (4π μ_m h a²)
\`\`\`

- **Parameters:** μ_m = 0.1 Pa·s for a fluid bilayer and h = 4 nm (assumptions). *a* is half the hydrodynamic radius, kept within 1–5 nm, as a rough size of the membrane-embedded part. D depends on *a* only logarithmically: about 4 µm²/s for *a* = 1.5 nm.
- **Calibration per membrane (2026-11-11; \`pipeline/membrane_diffusion.py\` → \`library/membrane_diffusion.json\`, viewer option "calibrated per membrane", on by default).** The fluid-bilayer value is multiplied by a factor per membrane: the median of measured D ÷ the Saffman–Delbrück D of the same protein (\`pipeline/data/validation/membrane_diffusion_measurements.tsv\`; FRAP of transmembrane proteins in unperturbed cells):
  - plasma membrane 5.9 × 10⁻⁴ (GFP-Sso1 and GFP-Snc1, 0.0025–0.0029 µm²/s; Valdez-Taubas & Pelham 2003);
  - vacuole 0.016 (an Sso1 anchor retargeted to the vacuole, 0.068 µm²/s; same paper);
  - ER 0.035 (Sec61 and Heh2-derived reporters, 0.15–0.18 µm²/s; Meinema *et al.* 2013);
  - nuclear envelope 0.079 (0.09–0.92 µm²/s; same paper);
  - mitochondrial membranes, Golgi/endosomes and peroxisomes: no yeast data; the ER factor is used (assumption).
- **Caveats of the calibration:** FRAP measures micrometre-scale mobility over seconds, including domains, obstacles and large immobile fractions (40–80% on the plasma membrane; Spira *et al.* 2012), and the plasma-membrane values are approximate in the source (D = r²/4t½). Applying them as one factor to the short-time motion of a 30–200 nm box makes membrane proteins as slow on nanosecond scales as on micrometre scales, which is the same compromise as the cytoplasm factor (section 14.5). Immobile fractions are not modelled.
- **Not modelled:** membrane domains as such, the cell wall pressing on the plasma membrane, and hydrodynamic coupling within the membrane.
- **Measured motion:** in a 100 nm cortical-ER box over 5 µs, membrane proteins moved a median 2.7 nm and free soluble molecules 18 nm (fluid bilayer, before the calibration). In an 80 nm cortical-ER box the median membrane-protein D is 4.3 µm²/s uncalibrated and 0.15 µm²/s calibrated (mostly ER proteins).

**Hydrodynamic radius.** R_h of each species comes from its residue beads. The Kirkwood double sum 1/R_K = ⟨1/r_ij⟩ runs over bead pairs, subsampled to 700 beads. It is scaled by 6/5, the value for a uniformly filled sphere, where R_K = 5R/6, and a 0.3 nm hydration shell is added. Elongated species get a larger R_h than globular species of the same mass.

**Trimmed disordered regions** (2026-11-11). The library leaves out AlphaFold residues with pLDDT < 50, which makes IDR-rich proteins too small hydrodynamically. For each AlphaFold monomer species, \`features.py\` regrows the trimmed IDRs (runs of ≥ 10 residues) as self-avoiding C-alpha coils from their anchor residues (virtual bond angle 2.15 rad, chosen so that free coils of 50–400 residues follow the IDP scaling laws R_g ≈ 0.254 *N*^0.522 nm and R_h ≈ 0.249 *N*^0.509 nm within about 10%; Bernadó & Blackledge 2009, Marsh & Forman-Kay 2010; the exponents are (unchecked) against the papers' text). The ratio of the full chain's Kirkwood R_h to the folded part's (\`library/rh_scale.json\`; median 1.1, up to 3) multiplies the species' R_h in the simulation.

**Diffusion in water** at 30 °C (η = 0.797 mPa·s) follows Stokes–Einstein–Debye:

\`\`\`
D_t = kT / (6π η R_h)
D_r = kT / (8π η R_h³)
\`\`\`

The Viscosity setting divides both coefficients. It stands for a fluid phase more viscous than water.

**Mean-field hydrodynamic interactions** (on by default) follow Mereghetti & Wade (2012) and Balbo *et al.* (2013). Every molecule's coefficients are multiplied by the short-time self-diffusion factors of hard spheres at the box's macromolecular volume fraction φ, using second-order expansions that include lubrication (Cichocki, Ekiel-Jeżewska & Wajnryb 1999):

\`\`\`
D_t / D_t0 = 1 − 1.832 φ − 0.219 φ²
D_r / D_r0 = 1 − 0.631 φ − 0.726 φ²
\`\`\`

- **Volume fraction.** φ is the hydrated volume of all molecules in the box (partial specific volume 0.73 mL/g plus 0.35 g/g bound water) divided by the box volume. Membranes and the lumen of other compartments count as free volume, so boxes that cut through organelles underestimate φ.
- **Typical values.** Cytoplasm boxes give φ ≈ 0.2, which makes D_t 0.63 times its value in water.
- **What it leaves out.** Explicit many-body hydrodynamics. Ando & Skolnick (2010) used Stokesian dynamics with lubrication on equivalent spheres at 250–350 mg/mL (sphere volume fractions 0.43–0.60) and found GFP at D_L/D₀ ≈ 0.08 at 300 mg/mL. That density is far above this model's (cytoplasm boxes 160–185 mg/mL, R_h-sphere volume fraction 0.32), so their number does not carry over.
- **Far-field many-body hydrodynamics tested** (2026-11-11, \`simgpu/\`): the Ewald-summed Rotne–Prager–Yamakawa mobility of all molecules of a 120 nm cytoplasm box (polydisperse, with the overlap correction of Zuk *et al.* 2014; Rotne & Prager 1969; periodic sum after Hasimoto 1959) slowed GFP-sized molecules to 0.56 of their speed in water, against 0.62 with no hydrodynamic interactions and 0.36 with the mean-field factors. Pairwise far-field hydrodynamics leaves each molecule's own mobility at its free value; the short-time slowdown comes from the near field (lubrication), which is ill-defined for the overlapping hydrodynamic spheres of shaped molecules (0.6 overlapping neighbours per molecule in the packed box). The mean-field factors stay the default.

### 14.3 Interactions and integration

**Overlap.** A harmonic penalty applies between spheres of different bodies: U = ½ k (r_i + r_j − d)² for d < r_i + r_j, with k = 20, 30 and 50 kT/nm² in the schematic, mid and residue models.

**Stickiness.** A short-range nonspecific attraction of depth ε acts per sphere pair:

\`\`\`
U = −ε cos²(π (d − r_i − r_j) / (2 · 0.6 nm))   for 0 ≤ d − r_i − r_j ≤ 0.6 nm
\`\`\`

McGuffee & Elcock (2010) needed such an attraction to reproduce GFP diffusion in *E. coli*. It is applied per sphere pair, so the effective well between two molecules is several ε, and its effect differs between models (section 14.5).

**Electrostatics** (on by default). Charged spheres interact by the screened Coulomb (Debye–Hückel) potential of two charged spheres of radii a_i, a_j:

\`\`\`
U = l_B q_i q_j exp(−κ (d − a_i − a_j)) / ((1 + κ a_i)(1 + κ a_j) d)      (kT units; d ≥ a_i + a_j, else its contact value)
\`\`\`

- **Bjerrum length** l_B = 0.72 nm (water at 30 °C).
- **Screening.** κ = 3.30 √I nm⁻¹ for ionic strength I in M. The default I = 0.3 M comes from the in vivo-like yeast assay medium of van Eunen *et al.* (2010): 300 mM K⁺, 245 mM glutamate, 50 mM phosphate, 20 mM Na⁺, 2 mM free Mg²⁺, about 0.33 M ionic strength. That gives a Debye length of 0.55 nm, so electrostatics act over about a nanometre. The Salt setting changes it.
- **Charges:**
  - Asp and Glu −1, Lys and Arg +1, His neutral at the cytosolic pH of 7.2 (Orij *et al.* 2009).
  - Each nucleotide −0.25 e. On a double helix counterion condensation leaves the fraction 1/ξ = b/l_B of the phosphate charge (Manning 1969; b = 0.17 nm).
  - A sphere carries the summed charge of the residues it represents. In the residue model, buried residues give their charge to the nearest exposed bead.
- **Cutoff.** Three Debye lengths beyond contact, with the force shifted to zero there.
- **Where the charges come from.** The structure library keeps only bead coordinates, so residue identities are recovered from the source structures by \`pipeline/simdata.py\`:
  - **Matching.** A library species' protein beads are its source's C-alpha atoms after a rigid motion (centring, membrane orientation) and the removal of some residues (low pLDDT, unused chains). Windows of four consecutive C-alphas are matched by their six internal distances. A window pose that brings 90% of the beads within 0.03 nm of a source C-alpha is accepted, and every bead takes the residue of its nearest one.
  - **Sources:** AlphaFold DB monomer models, the AlphaFold DB dimers and PDB assemblies, fetched from RCSB when not cached.
  - **Composites:** structures such as a ribosome with its translocon take the part found in the structure of the species they extend.
  - **Coverage:** 98.6% of the 3.25 M protein beads are matched. The 11 species that are not are custom-built scaffolds (actin patch, tropomyosin path, γ-TuRC, ATP synthase dimers) and four small proteins; they stay neutral.
  - **Checks:** net charges match the sequence (actin −10 in the library, about −12 from its sequence). The ribosome carries −4,017 e of bare charge (5,209 nucleotides, 2,191 basic and 999 acidic protein residues). That shrinks to about −110 e after condensation.
  - **Output:** \`library/charges.bin\`, one byte per bead.

**Specific binding** (on by default). Proteins that form a confident AlphaFold DB dimer (ipTM ≥ 0.5; the same predicted heterodimers and homodimers as the packing, section 5) attract each other at their predicted interface when both are free in the box:
- **Contacts.** Native contacts are the inter-chain residue pairs with C-alpha atoms within 0.8 nm in the dimer model, mapped onto the library beads of the two monomer species.
- **Coverage.** 371 dimer models qualify, out of 1,848 confident ones; the others lack a monomer species of one of the partners in the library. The median interface has 61 contacts.
- **Energy.** Every contact pair between two molecules of the partner species gets a well of the same cos² shape as the stickiness, of depth E/n. E is the energy of a complete interface, 20 kT by default, and n the number of distinct contact pairs at the model's resolution (residues, clusters or sphere-tree spheres), so the total does not depend on the model.
- **Orientation.** The attraction is orientation-specific, since only the predicted interface residues attract: a Gō-like interface (Levy, Wolynes & Onuchic 2004).
- **One energy for all.** Affinities are not known per pair, so E is the same for every complex. Complexes already packed as one rigid body (section 5) are unaffected.
- **Output:** \`library/contacts.json\` and \`contacts.bin\`.

**Polysome mRNA links** (on by default). Linked ribosomes (section 14.1) are joined between the mRNA exit of one and the entry of the next by a worm-like chain (Marko & Siggia 1995):

\`\`\`
F = (kT / l_p) · (1 / (4 (1 − x/L_c)²) − 1/4 + x/L_c)        (x capped at 0.95 L_c)
\`\`\`

- **Contour length.** In growing yeast, ribosomes sit one per 156 nt of ORF on average (0.64 per 100 nt; Arava *et al.* 2003). The ribosome covers the 29 nt of its channel, so a linker holds about 127 nt. At 0.49 nm per nucleotide (the contour length fitted to rU40 by Chen *et al.* 2012), L_c = 62 nm.
- **Stiffness.** Chen *et al.* (2012) measured a persistence length l_p of 1.0 nm at high salt and 2.5 nm without added salt for single-stranded rU40. 1.5 nm is used.
- **Result.** The rms end-to-end distance is about 14 nm, so a polysome stays together as a loose cluster. Every linker gets the mean length.
- **Drawing.** The links are drawn as dark brown lines between the anchors.

**Integration.** Each step:
1. **Binning.** The spheres go into two uniform grids. Grid A holds spheres up to the 90th-percentile radius, with cells of that diameter plus the attraction range. Grid B, with larger cells, holds the bigger spheres (schematic and mid) and the cores, so the common small spheres need not scan far for the rare big ones. Per-cell capacities come from the fullest cell of the starting configuration ×1.5; an earlier fixed capacity of 24 overflowed and gave asymmetric forces.
2. **Forces.** Each sphere scans the cells within its interaction range of both grids and skips:
   - cells whose box is farther than that range;
   - cells holding only spheres of its own molecule (each cell records the lowest and highest body id binned into it);
   - in the shell between contact range and the electrostatic cutoff, uncharged spheres (flagged in the cell items);
   - core spheres, unless the scanning sphere is a core itself.

   Cells can still fill beyond capacity when molecules cluster tightly. The panel reads the grid occupancy back every two seconds, and every eight command buffers while computing ahead. On any overflow it doubles the capacities before the next submission.
3. **Links.** Bond forces are computed per body.
4. **Reduction.** One workgroup per body sums its force and torque over its spheres (parallel reduction).
5. **Move.** An overdamped Euler–Maruyama step:
   \`\`\`
   dx = D_t F dt + √(2 D_t dt) ξ
   dθ = D_r τ dt + √(2 D_r dt) ξ        (kT = 1)
   \`\`\`
   The deterministic displacement is capped at 0.15 nm per step, both for the body centre and for the body surface under rotation.

All steps run in WebGPU compute shaders, up to 512 steps in one compute pass. Computing ahead, each command buffer holds about 0.1 s of GPU work, so that no submission trips the GPU watchdog.

**Speed.** Before the two-level grid, the cell skipping and the parallel per-body sums, a 60 nm residue box (surface beads, no electrostatics) ran at 279 steps/s. Half the time went to the per-body sums, which one thread per body had done serially. After: 2,145 steps/s (7.7×). A 120 nm schematic box went from 643 to 1,546 steps/s (2.4×). The rewrite left diffusion unchanged within noise (section 14.5). Steps per second with the cores, Apple-silicon laptop GPU, cytoplasm boxes:

| Model | Box | Collision spheres | Electrostatics on | Off |
|---|---|---|---|---|
| Residue | 60 nm (185 molecules) | 96,432 | 634 (0.063 µs/s) | 1,631 (0.16 µs/s) |
| Mid | 60 nm | 39,900 | 1,185 (0.24 µs/s) | 2,097 (0.42 µs/s) |
| Schematic | 120 nm (1,476 molecules) | 41,847 | 944 (0.19 µs/s) | 1,441 (0.29 µs/s) |

About half of the residue model's surface beads are charged, and each scans 1.65 nm further, so electrostatics cost 2.6× there; switching them off is the quickest speed-up.

### 14.4 Time step and convergence

Long steps overestimate diffusion. The table gives the diffusion of the GFP-sized class (R_h < 3 nm) in a 120 nm cytoplasm box, mean-field hydrodynamics on, time-averaged MSD at a 0.7 µs lag, in µm²/s:

| Stickiness ε | dt 0.02 ns | 0.05 ns | 0.1 ns | 0.2 ns | 0.5 ns |
|---|---|---|---|---|---|
| 0 | | 61.9 | 60.5 | 60.1 | 71.8 |
| 1 kT | 42.4 | 40.0 | 47.3 | 58.4 | |
| 2 kT | 13.1 | 16.0 | 24.8 | | |

- **Without attraction,** dt ≤ 0.2 ns is converged.
- **With attraction,** the 0.6 nm well needs about dt ≤ 0.05/ε ns. Contact forces summed over many sphere pairs are stiff, and a longer step lets molecules escape the well too easily. A broader 1.2 nm well did not help: more pairs fall inside it.
- **Mid model.**
  - **Without attraction,** the GFP-sized class gives 51.7, 51.6 and 51.0 µm²/s at dt 0.05, 0.1 and 0.2 ns in an 80 nm box: converged at 0.2 ns.
  - **With attraction,** the schematic rule is not enough. At 1 kT the class gave 17.9 µm²/s at dt 0.05 ns and 11.8 at 0.01 ns. Its clusters meet in many pairs per contact, which makes the summed contact forces stiff. The panel uses dt ≤ 0.01/ε ns; smaller steps were not tested.
  - **Unstable at 1.5 kT.** Even dt 0.01 ns failed there: large complexes moved faster than in water (28.6 against 23.6 µm²/s). The Stickiness slider is therefore capped at 1 kT in the mid model and 0.5 kT in the residue model (2 kT in the schematic model).
- **The panel enforces these limits.** The residue model is assumed to scale the same way (dt ≤ 0.0125/ε ns, at most 0.1 ns); this has not been measured.
- **Native contacts** make wells deeper than ε for the few bound pairs (E/n per pair: about 0.3 kT per residue pair, a few kT per sphere-tree pair at E = 20 kT). The time step is not lowered for them, so bound complexes in the schematic model may come apart somewhat too easily.

### 14.5 Calibration against yeast measurements

The readout reports, per size class, the time-averaged mean squared displacement MSD(τ) of the moving molecules over the recording, as D = MSD/(6τ). It uses the longest lag whose rms displacement stays below a quarter of the free interior of the box; longer lags feel the confinement by the box faces.

**Before electrostatics** (schematic model, 120 nm cytoplasm box of the current packing, φ = 0.198, periodic boundary, all 1,476 molecules moving, converged steps, mean-field hydrodynamics on unless noted):

| Stickiness ε | GFP-sized (R_h < 3 nm) | R_h 3–6 nm | R_h > 6 nm (ribosomes, large complexes) |
|---|---|---|---|
| 0, hydrodynamics off | 69 µm²/s (0.62 of water) | 42 | 7.6 |
| 0 | 43 (0.39) | 29 | 4.1 |
| 1 kT | 24 (0.21) | 12 | 2.4 |
| 1.5 kT | 10.0 (0.090) | 6.9 | 1.9 |
| 2 kT | 8.0 (0.072) | 6.8 | 2.4 |

The earlier frozen-shell boundary gave the same numbers within noise (62, 40, 23, 9.0 and 6.3 µm²/s for the GFP-sized class). Before the proteome-mass correction (L1; φ = 0.106, 1,051 molecules) the same box gave 77, 62, 40–42 and 13 µm²/s for the GFP-sized class. The large class holds only 15–25 molecules, so its values scatter. After the engine rewrite (section 14.3) the same box gave 41.3 and 10.4 µm²/s at ε = 0 and 1.5 kT.

**With electrostatics** (0.3 M; GFP-sized class, µm²/s; schematic in the 120 nm box, mid and residue in an 80 nm box of 446 molecules, 113 of them GFP-sized; 1 µs runs; mid and residue with core spheres):

| Stickiness ε | Schematic | Mid | Residue |
|---|---|---|---|
| 0 | 42.2 | 54.2 | 50.9 |
| 0.25 kT | | | 56.5 |
| 0.5 kT | | | 31.9 |
| 1 kT | 23.9 | 17.9 (11.8 at dt 0.01 ns) | |
| 1.25 kT | | 22.7 (unstable) | |
| 1.5 kT | 10.0 | 24.6 (unstable) | |
| 2 kT | 6.6 | | |

- **Electrostatics change little at cytoplasmic salt.** Without attraction they slowed the GFP-sized class by 0–4% (schematic 41.3 → 42.2, mid about 51 → 49). The Debye length is 0.55 nm, and most charged groups never come that close.
- **Without attraction the three models agree** within the noise of these runs: 42–54 µm²/s. The schematic model is slowest, since its sphere trees are bulkier than the residues.

**Stickiness slows molecules by clustering them.** Clusters were counted in 80 nm boxes after 1 µs. Two molecules count as in contact when any two of their surface spheres lie within 0.3 nm of touching.

| Model, ε | Largest clusters (of 446 molecules) | Single molecules | Contacts per molecule | GFP-sized D (µm²/s) |
|---|---|---|---|---|
| Schematic, 0 | 20, 18, 14 | 192 | 0.9 | 46 |
| Schematic, 1 kT | 399 | 16 | 3.8 | 21 |
| Schematic, 1.5 kT | 247, 123 | 7 | 4.4 | 10.4 |
| Schematic, 2 kT | 429 | 6 | 4.5 | 6.4 |
| Mid, 0 | 15, 10, 9 | 255 | 0.6 | 47 |
| Mid, 1 kT | 184, 109, 90 | 7 | 4.6 | 15 |
| Mid, 1.5 kT (unstable steps) | 116, 53, 51 | 4 | 12.0 | 27 |
| Residue, 0.5 kT | 285, 91, 14 | 12 | 9.3 | 25 |

- **Schematic model.** The stickiness that matches the GFP measurement (1.5–2 kT) joins most molecules into one network within a microsecond, with voids between its strands. The earlier calibration matched the number this way.
- **Mid model.** Its clusters of ~4 residues meet in many more pairs per contact than sphere-tree spheres, so the same ε binds much more strongly. At 1 kT most molecules aggregate. The rise in D at 1.25–1.5 kT is an integration artefact, not physics: those steps were too long for the stiff many-pair contacts (section 14.4).
- **Residue model.** Most molecules aggregate already at 0.5 kT.
- **What the gap really is (corrected 2026-11-11).** This section used to say that full many-body hydrodynamics closed the gap in *E. coli* (Ando & Skolnick 2010) and suggested a Viscosity of 5–6× in its place. That result was obtained at a much higher density (section 14.2), and adding far-field many-body hydrodynamics to this model's boxes changed little. The remaining gap is now calibrated on measurements instead (below).
- **Defaults.** Stickiness defaults to 0. Viscosity defaults to **2.75×**, the calibration below. The panel's readout says what each choice does for the model in use.

**Calibration on a literature benchmark** (2026-11-11; \`validation/diffusion.py\`, section 9.1). Diffusion of a molecule is modelled as D = D₀(R_h) × crowding(R_h) × *s*: D₀ from Stokes–Einstein with the R_h of the measured construct (the species that diffuses, with its trimmed IDRs restored and its fluorescent tag added as a bead model), crowding(R_h) from Brownian dynamics of packed boxes with mean-field hydrodynamics on 4 GPUs (\`simgpu/sweep.py\`: ln(D/D₀) = −1.04 − 0.099 R_h/nm in a 120 nm cytoplasm box, −1.55 − 0.142 R_h/nm in a 100 nm nucleoplasm box with its fixed chromatin), and *s* a remaining compartment factor fitted on the free tracers of the training papers: 0.365 for the cytoplasm (hence Viscosity 2.75×), 0.26 for the nucleoplasm (single-molecule tracking only). On the held-out papers the size dependence from the packed boxes beats size alone (within-paper ratios RMSE 0.11 vs 0.21 log10 units). Directly simulated 40 nm GEM particles move at D/D₀ = 0.04–0.15 in a 200 nm box, consistent with the fit; with *s* this predicts 0.25 µm²/s against 0.3 measured in both yeasts (Delarue *et al.* 2018; Garner *et al.* 2023; Molines *et al.* 2022). A binding term fixed in advance (RNA-interactome hits plus abundance-weighted pool partners) fitted the training ratios slightly better but made the held-out ones worse, so it is not used: what slows Srv2 or Ede1 (0.48 and 0.59 µm²/s, about 10× and 2× below the size prediction) is not captured by generic scores.

**Polysomes and large particles: a test that fails** (2026-11-12; \`simgpu/gade_test.py\`, links in \`simgpu/engine.py\`; nothing was fitted). Dissolving polysomes with hippuristanol speeds 40 nm GEM particles 1.4–1.5 times and 20 nm ones 1.3 times, while removing 20–25% of ribosomes in *stm1Δ* changes nothing (Gade *et al.* 2025). The server-side engine got the polysome links the browser has (neighbouring ribosomes of one mRNA tied by entropic springs), and six 200 nm boxes of the cytosol of build Q (about 85 ribosomes and 20–44 links each) were run for 10 µs with and without the links, and with 22% of the ribosomes removed; 72 tracers per condition, mean-field hydrodynamics.
- **Result:** without links the 40 nm particles are 1.02 times as fast (95% interval 0.98–1.06 at a lag of 300 ns), the 20 nm particles 0.98 times; with fewer ribosomes 1.05 times (measured 1.0).
- **What it means:** tying ribosomes together does not slow large particles in this model; ribosome crowding itself does, slightly. The measured effect needs something the simulation lacks: the mRNA as a physical chain between and beyond the ribosomes (it is not a body here), a mesh spanning several polysomes, or times beyond the 10 µs reached (tracking measures at 10–100 ms).
- **Not re-run on builds R and S:** the cytosol's density and polysome statistics are the same (sections 6.5 and 10.1).

**Measured values:**
- GFP in the fission-yeast cytoplasm: 8.6 ± 2.2 µm²/s, and 5.6 ± 2.8 µm²/s in the nucleus, by FRAP at room temperature (Delivani *et al.* 2013, checked in the text). This is roughly 0.1 of GFP in water.
- Free GFP in the budding-yeast cytoplasm by FCS: 8.7–21 µm²/s depending on the laboratory (Slaughter *et al.* 2007; Boeke *et al.* 2014; Pack *et al.* 2017). The benchmark table lists all values with quotes.
- 40 nm GEM particles in fission yeast: 0.4–0.6 µm²/s by particle tracking (Tan *et al.* 2025).
- 40 nm GEM particles in budding yeast (added 2026-11-12, held out: read from figures and placed in the test split whatever their paper hash, since they were seen after the fit): about 0.135 and 0.21 µm²/s in two experiments of Gade *et al.* 2025 (Fig. 2 of the preprint) and 0.20 in Kunzi *et al.* 2026 (preprint, Fig. 1D). The model predicts 0.25, so at this size it is 1.2–1.8× too fluid; Delarue *et al.* 2018's ~0.3, which is in the training split, sits on the other side.
- In budding yeast, ribosome crowding controls the diffusion of particles ≥ 20 nm but not of molecules ≤ 5 nm (Delarue *et al.* 2018).

**Interpretation:**
- Excluded volume, electrostatics and mean-field hydrodynamics slow GFP-sized molecules to about 0.3–0.5 of their speed in water. That is still about 2.7× faster than the benchmark's free tracers (and about 5× faster than the slowest GFP measurement), as found for bacterial models.
- A uniform nonspecific attraction can close the gap numerically, as in McGuffee & Elcock (2010). Here it does so by clustering the molecules into networks, so it is a poor stand-in for what the model lacks (many-body hydrodynamics, transient specific binding, metabolites, the mRNA and other polymers left out).
- Crowding itself matters: doubling the macromolecular density (L1) slowed GFP-sized molecules from 0.55 to 0.36 of their speed in water with no attraction.

### 14.6 Running, recording and playback

- **Placing the box.** "Place box" shows the outline at the centre of the view and frames it. While placing, the box follows the view centre: pan, orbit and zoom as usual, or click a point to bring it to the centre. The outline is drawn on top of everything, including the parts outside the cut slab.
- **Changing settings.** Changing the model, boundary, run mode, box size or the membrane option restarts a running simulation in place, with the same box.
- **Live.** Each drawn frame computes the chosen number of steps, and the molecules move as they are computed.
- **Compute first, then play.** The chosen length is computed with command buffers submitted back to back (up to eight in flight, no waiting for the display) behind a progress bar. The recording then plays in a loop. The Playback slider sets the rate, 1–60 recorded frames per second (24 by default).
- **Tracking.** Clicking a molecule in a running simulation tracks it; up to eight can be tracked, and clicking again drops one. The other molecules are drawn faded: one pixel in four, lightened, still clickable, without outlines. Each tracked molecule's centre path over the recording, up to the frame shown, is drawn as a coloured line on top of everything. A segment that crosses a periodic face is left out, since the molecule re-enters on the other side. The panel lists each tracked molecule's displacement since the start of the recording and its own diffusion estimate from successive frames.
- **Recording.** About 300 frames are recorded per run. In live mode recording continues indefinitely; past the frame budget (about 60 MB of poses) every other frame is dropped and the interval doubles, so the spacing stays uniform.
- **Throughput:** see the speed table in section 14.3. With stickiness, the shorter step lowers it in proportion. Live mode is slower, since it waits for the display.

### 14.7 Export of MD-ready files

"Export MD files" saves five files:

| File | Contents |
|---|---|
| \`yeast_box.cif\` | Coarse-grained mmCIF. One pseudo-atom per residue (CA) or nucleotide (P) bead, one chain per molecule. Chain ZZZZ holds membrane surface beads from the wall grid. Occupancy is 1 for moving molecules, 0 for fixed ones. |
| \`yeast_box.json\` | Manifest: per molecule the species, genes, UniProt, structure source (PDB entry or AlphaFold DB model), region, pose (position, quaternion), R_h, D_t, D_r, mobility; plus the settings, φ and the hydrodynamic factors. |
| \`yeast_box_trajectory.json\` | The recorded rigid-body trajectory, up to 100 frames. |
| \`run_openmm.py\` | Runs the box in OpenMM (Eastman *et al.* 2017, 2024) as a coarse-grained model, with the same boundary (periodic boundary conditions for periodic boxes, walls otherwise, plus the shell tethers of embedded boxes) and the same interactions where enabled. |
| \`trajectory_to_dcd.py\` | Expands the browser trajectory into a bead-level DCD that VMD, ChimeraX, MDAnalysis or mdtraj open together with the mmCIF. |

The manifest also lists the interactions of the browser run, at residue level: the charge of every charged bead, the mRNA links (between the beads nearest the channel anchors, with contour and persistence lengths), and the native contact pairs between every two molecules of partner species, with their share of the interface energy.

The OpenMM model has:
- an elastic network within each moving molecule (bead pairs closer than 1.2 nm, 500 kJ/mol/nm²);
- Weeks–Chandler–Andersen excluded volume between molecules, plus the same cos² attraction;
- Debye–Hückel electrostatics between the charged beads (cut off three Debye lengths beyond contact), worm-like-chain mRNA links and cos² native-contact wells, as exported;
- massless (fixed) beads for membranes, filaments and chromatin, with fixed–fixed pairs switched off;
- soft box walls;
- Langevin dynamics at 303 K with a 20 fs step after energy minimisation; DCD output.

Tested: a 60 nm plasma-membrane box gave 92,606 beads, 2,896 of them membrane beads, 18,037 moving and 266,635 network bonds. It runs on the OpenCL platform of OpenMM 8.6. On the cluster (2026-11-11) a 60 nm cytoplasm box (175,045 beads, 2.4 M network bonds) runs on CUDA at about 700–1,100 ns/day per RTX 4090. It first blew up within 80–7,200 steps: ribosome models carry about 20 bead pairs 0.001–0.1 nm apart (overlapping protein and RNA beads), which got network bonds of near-zero length. Pairs closer than 0.1 nm now get no bond, and the excluded volume continues linearly below 0.8 σ; the box then ran 50,000 steps (1 ns).

For atomistic or Martini systems, take each molecule's structure from the manifest's source and apply its pose to the structure's bead frame (CA/P positions centred as in the library).

### 14.8 Limitations

- Rigid bodies: no internal flexibility, and disordered tails stay in their AlphaFold pose.
- No reactions.
- Electrostatics are Debye–Hückel between point charges at sphere centres: no dielectric boundary, no ion correlation, no explicit Mg²⁺. In the schematic and mid models a sphere's summed charge sits at its centre.
- Specific binding exists only for pairs with a confident predicted dimer whose partners are both monomer species, with one interface energy for every pair. No other complexes, no measured affinities.
- mRNA links use the mean linker length. Free mRNA threads, the mRNA inside polysomes and links leaving the box are absent.
- Hydrodynamics only as a mean-field rescaling with one volume fraction for the whole box. The remaining slowdown is a calibrated, uniform viscosity (2.75×, section 14.5): it does not distinguish what causes it (fluid-phase viscosity, missing polymers such as mRNA and IDRs, transient binding), and specific slow-downs of individual proteins (Srv2 about 10× below its size prediction) are not reproduced.
- Membrane proteins: Saffman–Delbrück with a fluid bilayer gives about 1–4 µm²/s, while yeast plasma-membrane proteins measured by FRAP and tracking move at 0.0004–0.05 µm²/s (Bianchi *et al.* 2018; Slaughter *et al.* 2013). Since 2026-11-11 a factor per membrane, fitted on FRAP data, closes this gap on average (section 14.2), but nothing mechanistic (domains, the cell wall, immobile fractions) is modelled, and three of the eight membranes borrow the ER factor.
- Filaments and chromatin are fixed.
- No mRNA, lipids or cell-wall polymers as bodies.
- In embedded and closed boxes the faces confine the moving molecules, so diffusion over distances near the box size is not meaningful; the readout's lag limit avoids it. Periodic boxes let molecules meet their own images across a small box.
- Results depend on the time step beyond the limits in 14.4.

### 14.9 References for the simulation

| Reference | Used for |
|---|---|
| Saffman PG, Delbrück M (1975). *PNAS* 72:3111–3113 | Lateral and rotational diffusion of membrane proteins |
| Valdez-Taubas J, Pelham HRB (2003). *Curr. Biol.* 13:1636–1640; Meinema AC *et al.* (2013). *Traffic* 14:487–501; Spira F *et al.* (2012). *Nat. Cell Biol.* 14:640–648 | Membrane-protein diffusion in yeast (calibration per membrane) |
| van Eunen K *et al.* (2010). *FEBS J.* 277:749–760 | In vivo-like ionic composition of the yeast cytosol (ionic strength) |
| Orij R, Postmus J, Ter Beek A, Brul S, Smits GJ (2009). *Microbiology* 155:268–278 | Cytosolic pH 7.2 in growing yeast (residue charge states) |
| Manning GS (1969). *J. Chem. Phys.* 51:924–933 | Counterion condensation on polyelectrolytes (effective nucleotide charge) |
| Levy Y, Wolynes PG, Onuchic JN (2004). *PNAS* 101:511–516 | Gō-like models of protein binding through native interface contacts |
| Arava Y *et al.* (2003). *PNAS* 100:3889–3894 | Ribosome density on yeast mRNAs: 0.64 per 100 nt (one per 156 nt) |
| Chen H *et al.* (2012). *PNAS* 109:799–804 | Persistence length (1.0–2.5 nm) and contour length of single-stranded RNA |
| Marko JF, Siggia ED (1995). *Macromolecules* 28:8759–8770 | Worm-like-chain force–extension interpolation |
| McGuffee SR, Elcock AH (2010). *PLoS Comput. Biol.* 6:e1000694 | BD of a structure-based cytoplasm; attraction calibrated to GFP diffusion |
| Ando T, Skolnick J (2010). *PNAS* 107:18457–18462 | Crowding plus hydrodynamic interactions account for in-vivo diffusion; GFP D/D0 0.31–0.46 without and 0.08 with hydrodynamics |
| Mereghetti P, Wade RC (2012). *J. Phys. Chem. B* 116:8523–8533; Balbo J *et al.* (2013). *Biophys. J.* 104:1576–1584 | Mean-field treatment of hydrodynamic interactions in BD |
| Cichocki B, Ekiel-Jeżewska ML, Wajnryb E (1999). *J. Chem. Phys.* 111:3265–3273 | Short-time translational and rotational self-diffusion vs volume fraction |
| Yu I *et al.* (2016). *eLife* 5:e19274 | All-atom MD of a bacterial cytoplasm (scale comparison) |
| Delivani P, Chacón JM, Schroth-Diez B, Tolić-Nørrelykke IM (2013). *Bio-protocol* 3:e941 | GFP diffusion in fission-yeast cytoplasm and nucleus (FRAP) |
| Slaughter BD *et al.* (2007). *PNAS* 104:20320–20325; Boeke D *et al.* (2014). *Mol. Syst. Biol.* 10:756; Pack CG *et al.* (2017). *BMB Rep.* 50:478–483 | Free GFP and tagged proteins in budding-yeast cytoplasm (FCS; benchmark) |
| Garner RM *et al.* (2023). *Biophys. J.* 122:767–783; Molines AT *et al.* (2022). *Dev. Cell* 57:466–479 | 40 nm GEM diffusion in fission yeast (benchmark) |
| Slaughter BD *et al.* (2013). *Nat. Commun.* 4:1380; Bianchi F *et al.* (2018). *Nat. Commun.* 9:501 | Plasma-membrane protein diffusion in yeast |
| Rotne J, Prager S (1969). *J. Chem. Phys.* 50:4831–4837; Zuk PJ *et al.* (2014). *J. Fluid Mech.* 741:R5; Hasimoto H (1959). *J. Fluid Mech.* 5:317–328 | Rotne–Prager–Yamakawa mobility, unequal and overlapping spheres, periodic (Ewald) sum (\`simgpu/hydro.py\`) |
| Bernadó P, Blackledge M (2009). *Biophys. J.* 97:2839–2845; Marsh JA, Forman-Kay JD (2010). *Biophys. J.* 98:2383–2390 | Size scaling of disordered chains (IDR hydrodynamic radii) |
| Tan C, Lanz MC, Swaffer M, Skotheim J, Chang F (2025). *Mol. Biol. Cell* 36:ar51 | 40 nm GEM diffusion in fission yeast |
| Delarue M *et al.* (2018). *Cell* 174:338–349 | Ribosome crowding sets mesoscale diffusion in budding yeast |
| Eastman P *et al.* (2017). *PLoS Comput. Biol.* 13:e1005659; Eastman P *et al.* (2024). *J. Phys. Chem. B* 128:109–116 | OpenMM |
| Souza PCT *et al.* (2021). *Nat. Methods* 18:382–388; Stevens JA *et al.* (2023). *Front. Chem.* 11:1106495 | Martini 3; Martini model of an entire minimal cell |
| Andrews SS (2017). *Bioinformatics* 33:710–717; Hoffmann M, Fröhner C, Noé F (2019). *PLoS Comput. Biol.* 15:e1006830 | Particle-based reaction–diffusion (Smoldyn, ReaDDy 2) |
| Roberts E, Stone JE, Luthey-Schulten Z (2013). *J. Comput. Chem.* 34:245–255; Thornburg ZR *et al.* (2022). *Cell* 185:345–360 | Lattice Microbes; whole minimal-cell simulation |
`,v=o=>document.getElementById(o),Be=o=>o>=1e6?(o/1e6).toFixed(o>=1e7?1:2)+" M":o>=1e4?(o/1e3).toFixed(0)+" k":o.toLocaleString("en-US"),Bo=o=>o>=1e3?(o/1e3).toFixed(o>=1e4?0:2).replace(/\.?0+$/,"")+" µm":(o>=10?o.toFixed(0):o.toFixed(1))+" nm";let K,$e,Oe,No,Ue,Et,W,H,C;const B=new js,q=v("gpu"),R={cutOn:!0,plane:[0,0,1,0],slab:60,colorMode:0,regionMask:4294967295,freeHide:0,highlight:!1,selId:0,shaded:!0,organelle:!1,motion:0,l0px:1.6,l1px:2.5,minpx:.3,outline:1,showMeshes:!0,showLipids:!0,renderScale:1,memColors:!1};let ze=!1,D=!0,Tn=0,Zt=0,Eo="auto",ue;const ot={actin:["F_actin","Arp23","CP_cap","TPM1_coil"],microtubules:["MT_tubulin","gTuRC_nuc","gTuRC_cyto"],ribosomes:["ribosome_80S","ribosome_80S_RAC","ribosome_80S_Ssh1","ribosome_80S_Sec61","mitoribosome_74S_IM","CPX_1602_54S_mitochondrial_large_ribosoma","CPX_1603_37S_mitochondrial_small_ribosoma","preribo_90S","preribo_60S_nucleolar","preribo_60S_nucleoplasm","preribo_60S_cytoplasm","preribo_40S"],tRNA:["tRNA"],nucleosomes:["nucleosome"],pores:["NPC"],polymerases:["PolI","PolII","PolIII"],eisosomes:["Pil1_dimer","Lsp1_dimer"]},si=[...ot.actin,...ot.microtubules,...ot.eisosomes,"NPC","nucleosome"];let Se=-1,st=null,te=null,Pe=null;const ke=[];function St(o,e){v("loading").querySelector(".m").textContent=o,v("loading").querySelector(".pb div").style.width=`${Math.round(e*100)}%`}function Qn(o){const e=v("loading").querySelector(".err");e.style.display="block",e.innerHTML=o,v("loading").querySelector(".t").textContent="Could not start the viewer"}function Xe(o,e=2200){const t=v("toast");t.textContent=o,t.style.display="block",clearTimeout(Xe._t),Xe._t=setTimeout(()=>t.style.display="none",e)}async function ii(){if(!("gpu"in navigator)){Qn(`This browser does not expose <b>WebGPU</b>. Use a current Chrome or Edge (version 113 or later), or Safari 26 or later.
      On Linux Chrome, enable <code>chrome://flags/#enable-unsafe-webgpu</code> and <code>#enable-vulkan</code>.
      Firefox 141+ supports WebGPU on Windows.<br><br>Also check that the page is served over <code>http://localhost</code>
      (for example <code>python -m http.server</code>). Opening the file directly with <code>file://</code> does not work.`);return}try{St("Loading cell geometry",.01),[Oe,Ue]=await Promise.all([Ae("geometry.json"),Ae("bookmarks.json")]),Oe.F=await Ks(Oe,He),No=new Zs(Oe),$e=await ps(),K=await ms(St),ot.tRNA=K.types.filter(a=>a.id==="tRNA"||a.id.startsWith("tRNA_")).map(a=>a.id);const o=await He("types_tiles.bin"),e=K.types.length;Et={off:new Uint32Array(o,0,e+1),tiles:new Uint16Array(o,(e+1)*4)},St("Starting WebGPU",.7),W=new bs(q),W.onLost=a=>{Xe(a,1e4)},await W.init(K,Oe),ue=new Uint32Array(e);for(const a of Object.values(ot))for(const i of a){const h=K.byId.get(i);h!==void 0&&(ue[h]|=4)}for(const a of si){const i=K.byId.get(a);i!==void 0&&(ue[i]|=8)}W.device.queue.writeBuffer(W.typeFlags,0,ue),St("Compiling shaders and starting tile streaming",.9);let t;if(window.__YC_EMBED__){t={};for(const[a,i]of Object.entries(window.__YC_EMBED__))(a.startsWith("tiles/")||a.startsWith("library/"))&&(t[a]=lo(i))}const n=Math.max(2,Math.min(4,(navigator.hardwareConcurrency||4)-1));H=new Qs($e,W,K,Oe,n,t),H.onChange=()=>{D=!0},C=new Ws(W,H,K),bi(),xi();const s=new URLSearchParams(location.search).get("view"),r=s===null?0:isNaN(Number(s))?Math.max(0,Ue.findIndex(a=>a.name.toLowerCase().startsWith(s.toLowerCase()))):Math.min(Ue.length-1,Math.max(0,Number(s)));it(r,!1),v("loading").style.display="none",requestAnimationFrame(Wo),window.yc={R:W,T:H,cam:B,fs:R,lib:K,idx:$e,sim:C,applyBookmark:it,setType:at,setNetwork:Ut,pickAt:vt,bookmarks:Ue,setQuality:Go,describe:En,ready:!0,state:()=>({pending:H.pendingForView(),detailPending:H.detInflight,loaded:H.loaded.size,totals:H.totals(),stats:W.stats,cam:{t:B.target,d:B.dist},plane:R.plane,gpuMs:$o,gpuBusy:Lt,dirty:D}),render:()=>{D=!0}}}catch(o){console.error(o),Qn(`${String(o?.message||o)}<br><br>See the browser console for details.`)}}function Ce(){return[R.plane[0],R.plane[1],R.plane[2]]}function Me(o,e){o=he(o),R.plane=[o[0],o[1],o[2],e],Ie(),D=!0}function Ye(){const o=B.fwd;Me(o,oe(o,B.target))}let Ne=.5;try{const o=Number(localStorage.getItem("yc.nudgeStep"));o>0&&(Ne=o)}catch{}function nt(o){const e=Ce();B.target=J(B.target,O(e,o)),Me(e,R.plane[3]+o)}function ai(o){const e=q.clientWidth/q.clientHeight,t=he(o.normal);let n=he(Fe(o.up,O(t,oe(o.up,t))));return{target:[...o.target],dist:B.distForWidth(o.width,e),fwd:t,up:n}}function it(o,e=!0,t){const n=Ue[o];if(!n)return;(R.organelle||R.memColors)&&o!==0&&Ke("molecules"),n.show?(pn(n.show),tn=!0):tn&&(pn(null),tn=!1);const s=ai(n),r=he(n.normal),a=oe(r,n.target),i=Ce(),h=R.plane[3],l=R.slab;if(R.cutOn=!0,v("cutOn").checked=!0,ze=!1,v("followView").classList.remove("on"),!e){B.set(s),Me(r,a),R.slab=n.slab,Ie(),D=!0,t?.();return}B.flyTo(s,2600,()=>{Me(r,a),R.slab=n.slab,Ie(),t?.()},c=>{let d=he(J(O(i,1-c),O(r,c)));Je(J(O(i,1-c),O(r,c)))<.001&&(d=r),R.plane=[d[0],d[1],d[2],h+(a-h)*c],R.slab=Math.exp(Math.log(l)+(Math.log(n.slab)-Math.log(l))*c),Ie()}),wt(n),D=!0}function wt(o){const e=v("caption");if(!o){e.style.display="none";return}e.innerHTML=`<b>${o.name}</b>${o.caption}`,e.style.display="block",clearTimeout(wt._t),te||(wt._t=setTimeout(()=>e.style.display="none",9e3))}function _o(){te={i:0,phase:"fly",t:performance.now()},v("tourBtn").textContent="■ Stop",v("tourBtn").classList.add("on"),it(0,!0,()=>{te&&(te.phase="wait",te.t=performance.now())})}function je(){te=null,v("tourBtn").textContent="▶ Tour",v("tourBtn").classList.remove("on"),wt(null)}function ri(o){if(te){if(te.phase==="wait"&&(H.pendingForView()===0||o-te.t>14e3))te.phase="hold",te.t=o;else if(te.phase==="hold"&&o-te.t>6500){te.i=(te.i+1)%Ue.length,te.phase="fly";const e=te.i;it(e,!0,()=>{te&&te.i===e&&(te.phase="wait",te.t=performance.now())})}}}let qe=-1;function Ut(o){Se>=0&&(ue[Se]&=-2,Se=-1,Bn()),qe=o;const e=new Set(o>=0?[o,...(K.netNb[o]||[]).map(t=>t[0])]:[]);for(let t=0;t<ue.length;t++)e.has(t)?ue[t]|=1:ue[t]&=-2;R.highlight=o>=0,W.device.queue.writeBuffer(W.typeFlags,0,ue),fn(),st&&qt(),D=!0}function at(o){qe>=0&&Ut(-1),Se>=0&&(ue[Se]&=-2),Se=o,o>=0&&(ue[o]|=1),R.highlight=o>=0,W.device.queue.writeBuffer(W.typeFlags,0,ue),Bn(),D=!0}function Bn(){const o=v("selType");if(Se<0){o.style.display="none";return}const e=K.types[Se],t=H.countLoaded(Se);o.style.display="block",o.innerHTML=`<div><b>${e.gene}</b> ${e.name}</div>
    <div class="sub" style="margin:3px 0 6px">${Be(e.copies)} copies in the cell · ${Be(t)} in loaded tiles · ${Io(e)}</div>
    <div class="grid2"><button class="small" id="flyNearest">Fly to nearest</button><button class="small" id="clearType">Clear highlight</button></div>`,v("flyNearest").onclick=()=>Fo(Se),v("clearType").onclick=()=>at(-1)}function li(o){if(o=o.trim().toLowerCase(),!o)return[];const e=[];return K.types.forEach((t,n)=>{const s=t.gene.toLowerCase(),r=t.id.toLowerCase(),a=(t.uniprot||"").toLowerCase(),i=(t.name||"").toLowerCase(),h=(t.genes||"").toLowerCase();let l=0;s===o||r===o||a===o?l=100:s.startsWith(o)||r.startsWith(o)?l=60:h.split(/[;,\s]+/).includes(o)?l=55:i.includes(o)?l=30+(i.startsWith(o)?10:0):(h.includes(o)||a.startsWith(o)||(t.pdb||"").toLowerCase()===o)&&(l=20),l&&e.push({i:n,s:l+Math.log10(t.copies+1)})}),e.sort((t,n)=>n.s-t.s),e.slice(0,14).map(t=>t.i)}function Fo(o){const e=Et.off[o],t=Et.off[o+1];if(t<=e){Xe("No copies of this molecule were placed.");return}let n=-1,s=1e18;for(let d=e;d<t;d++){const f=Et.tiles[d],w=H.tileCenter(f),b=Je(Fe(w,B.target));b<s&&(s=b,n=f)}const r=H.tileCenter(n),a=K.types[o],i=Math.max(90,a.br*9),h=d=>{const f=H.nearestInTile(d,o,r);if(f<0)return;const w=H.decode(d,f),b=Ce(),m=B.state();m.target=w.p,m.dist=B.distForWidth(i,q.clientWidth/q.clientHeight),B.flyTo(m,900,()=>{Me(b,oe(b,w.p)-a.br*.4)}),R.plane[3]=oe(b,w.p)-a.br*.4;const u=H.poolIndex(d,f);Dt({kind:"inst",id:(fs|u)>>>0,info:H.instanceAt(u)})};H.loaded.get(n);const l=Ce(),c=B.state();c.target=r,c.dist=B.distForWidth(i*2.5,q.clientWidth/q.clientHeight),R.cutOn=!0,v("cutOn").checked=!0,(R.organelle||R.memColors)&&Ke("molecules"),B.flyTo(c,2200,()=>{Me(l,oe(l,r)-20);const d=H.loaded.get(n);if(d)h(d);else{const f=w=>{w.idx===n&&(H.onTileLoaded.splice(H.onTileLoaded.indexOf(f),1),h(w))};H.onTileLoaded.push(f)}}),Me(l,oe(l,r)-20)}async function vt(o,e){const t=W.width/q.clientWidth;return D=!0,W.pick(o*t,e*t)}async function Zn(o,e){const t=await vt(o,e);if(t.depth<=0)return null;const n=B.near()/t.depth,s=B.ray(o/q.clientWidth,e/q.clientHeight,q.clientWidth/q.clientHeight);return J(B.eye,O(s,n/oe(s,B.fwd)))}function Io(o){return o.src==="PDB"?"PDB "+o.pdb:o.src==="AF-multimer"?"AlphaFold complex":o.src==="pred"?"Boltz-2 prediction":"AlphaFold"}let tt={};Ae("library/diffusion.json").then(o=>{tt=o.types||{}}).catch(()=>{});let $t={};Ae("library/gene_loci.json").then(o=>{$t=o.loci||{}}).catch(()=>{});const Nn=["","I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI"];function ci(o){if(o=o.trim().toLowerCase(),o.length<2)return[];const e=[];for(const[t,n]of Object.entries($t)){const s=t.toLowerCase(),r=n[6].toLowerCase();if(s===o||r===o?e.unshift(t):(s.startsWith(o)||r.startsWith(o))&&e.push(t),e.length>30)break}return e.slice(0,4)}function di(o){const e=$t[o];if(!e)return;const t=[e[0],e[1],e[2]],n=Ce(),s=B.state();s.target=t,s.dist=B.distForWidth(260,q.clientWidth/q.clientHeight),R.cutOn=!0,v("cutOn").checked=!0,B.flyTo(s,900,()=>{Me(n,oe(n,t))}),Xe(`${o} (${e[6]}) gene locus: chromosome ${Nn[e[3]]}, ${e[4].toLocaleString()}–${e[5].toLocaleString()} bp (centre of its nucleosomes)`,4200)}let re=null,eo=!1;const lt=new Map,en=new Set;let un=()=>{};const hi=["H2A.Z","H3K4me3","H3K36me3","H3K79me3","H3K56ac","H3K9ac","H4K16ac"];function ui(o){if(!re)return eo||(eo=!0,Ae("library/genome.json").then(t=>{re=t,re.tileSet=new Set(t.tiles),re.iTypes=new Set(t.i_types),un()}).catch(()=>{})),null;if(!re.tileSet.has(o))return null;const e=lt.get(o);return e||en.has(o)?e??null:(en.add(o),Ae(`annot/${$e.ix[o]}_${$e.iy[o]}_${$e.iz[o]}.json`).then(t=>{if(t.im=new Map(t.i.map(n=>[n[0],n])),lt.set(o,t),lt.size>96){const n=lt.keys().next().value;lt.delete(n),en.delete(n)}un()}).catch(()=>{}),null)}const Q=o=>o.toLocaleString("en-US");function to(o,e,t){const n=o.filter(A=>A>0);if(n.length<5||!(e>0))return"";const s=0,r=Math.max(1,Math.ceil(Math.log10(Math.max(...n,e))*2)/2),a=Math.round((r-s)*4),i=new Array(a).fill(0),h=A=>Math.min(a-1,Math.max(0,Math.floor((Math.log10(A)-s)/(r-s)*a)));for(const A of n)i[h(A)]++;const l=300,c=46,d=2,f=(l-2*d)/a,w=Math.max(...i),b=h(e),m=n.filter(A=>A<e).length/n.length;let u="";i.forEach((A,M)=>{const _=A?Math.max(1.5,A/w*(c-4)):0,S=Math.round(10**(s+M/4)),g=Math.round(10**(s+(M+1)/4));u+=`<rect x="${(d+M*f+1).toFixed(1)}" y="${(c-_).toFixed(1)}" width="${(f-2).toFixed(1)}" height="${_.toFixed(1)}" rx="1.5" fill="${M===b?"var(--acc)":"rgba(255,255,255,0.26)"}"><title>${Q(S)}–${Q(g)} copies: ${Q(A)} ${t}</title></rect>`});const p=d+Math.min(a,Math.max(0,(Math.log10(e)-s)/(r-s)*a))*f;let x="";for(let A=0;A<=r;A++){const M=d+A*4*f,_=A<3?String(10**A):A<6?`${10**(A-3)}k`:`${10**(A-6)}M`;x+=`<line x1="${M.toFixed(1)}" x2="${M.toFixed(1)}" y1="${c}" y2="${c+3}" stroke="rgba(255,255,255,0.35)"/><text x="${M.toFixed(1)}" y="${c+13}" text-anchor="${A===0?"start":A===r?"end":"middle"}" font-size="9.5" fill="#9aa">${_}</text>`}const y=Math.floor(m*100);return`<div class="sub" style="margin-top:8px">Copies per cell (log scale): ${Q(e)}, ${y>=99?"among the highest 1%":y<1?"among the lowest":`more than ${y}%`} of the ${Q(n.length)} ${t}</div>
    <svg viewBox="0 0 ${l} ${c+16}" width="100%" style="display:block;margin-top:3px" role="img" aria-label="Histogram of copies per cell; this one has ${Q(e)}">
    ${u}<line x1="${d}" x2="${l-d}" y1="${c}" y2="${c}" stroke="rgba(255,255,255,0.35)"/>
    <line x1="${p.toFixed(1)}" x2="${p.toFixed(1)}" y1="0" y2="${c}" stroke="var(--acc)" stroke-width="2"/>${x}</svg>`}function Ct(o,e){const t=re.feat[o-1]||[],n=l=>l[3]<=6,s=l=>l[4]===l[5]?l[4]:`${l[4]} (${l[5]})`,r=t.filter(l=>l[0]<=e&&e<=l[1]).sort((l,c)=>Number(n(c))-Number(n(l)));if(r.length){const l=r.find(n);return{text:r.slice(0,3).map(d=>{const f=d[2]>0?e-d[0]+1:d[1]-e+1,w=d[1]-d[0]+1;return`in ${s(d)}, ${re.ftypes[d[3]]}${n(d)?`, bp ${Q(f)} of ${Q(w)} from its start`:""}`}).join("; "),gene:l||null}}let a=null,i=null;for(const l of t)n(l)&&(l[1]<e&&(!a||l[1]>a[1])&&(a=l),l[0]>e&&(!i||l[0]<i[0])&&(i=l));const h=[a,i].filter(Boolean).map(l=>{const c=l[2]>0?l[0]:l[1];return(l[2]>0?e<l[0]:e>l[1])?`${Q(Math.abs(c-e))} bp upstream of ${s(l)}`:null}).filter(Boolean);return{text:`intergenic${a&&i?` between ${a[4]} and ${i[4]}`:""}${h.length?` (${h.join("; ")})`:""}`,gene:null}}function Do(o){let e=-1,t,n="";if(o.kind==="free"&&(o.fk===0||o.fk===1)){const b=H.tileBySlot(o.slot);if(!b)return;e=b.idx}else if(o.kind==="inst"&&o.info?.t){if(re&&!re.iTypes.has(o.info.type))return;e=o.info.t.idx}else return;const s=ui(e);if(!s)return re&&!re.tileSet.has(e)?void 0:null;const r=s.m.length,a=s.l.length;if(o.kind==="inst"?(t=s.im.get(o.info.local),n="i"):o.fk===1?(t=s.l[o.chain],n="l"):o.chain<r?(t=s.m[o.chain],n="m"):s.x&&o.chain>=r+a&&(t=s.x[o.chain-r-a],n="x"),!t)return;const i=re.tx,h=b=>{const m=re.genes[i.g[b]];return{name:m[0],orf:m[1],u5:m[2],cds:m[3],u3:m[4],pa:m[5],L:i.n[b],nr:i.r[b],drawn:i.d[b]}},l=b=>["Transcript",`${Q(b.L)} nt: 5′ UTR ${b.u5}, CDS ${Q(b.cds)}, 3′ UTR ${b.u3}, poly(A) ${b.pa}`];if(n==="m"){const[b,m,u,p,x]=t,y=h(b),A=P=>P<y.u5?"5′ UTR":"CDS",M=P=>P<y.pa?"poly(A) tail":P<y.pa+y.u3?"3′ UTR":"CDS",_=(P,T)=>P===T?P:`${P} / ${T}`;let S="";if(m===0)S=`CDS, through ribosome ${x} of ${y.nr} (from the 5′ end)`;else if(m===1)S=`${_(M(p),M(u))}, nt ${Q(u+1)}–${Q(p+1)} from the 3′ end`;else if(m===2)S=`${_(A(u),A(p))}, nt ${Q(u+1)}–${Q(p+1)} from the 5′ cap`;else{const P=T=>{const E=y.drawn-1-T;return E<y.pa+y.u3?M(E):A(T)};S=`${_(P(u),P(p))}, nt ${Q(u+1)}–${Q(p+1)} of ${Q(y.drawn)}`}if(!re.txCount){re.txCount=new Array(re.genes.length).fill(0);for(const P of i.g)re.txCount[P]++}const g=re.txCount[i.g[b]],k=[["Gene",`${y.name}${y.name!==y.orf?` (${y.orf})`:""}`],["Copies / cell",`${Q(g)} mRNA${g===1?"":"s"} of this gene`],l(y),["Ribosomes",y.nr?`${y.nr} on this mRNA (polysome)`:"none (an untranslated mRNA)"],["This piece",S]];return m===3&&y.drawn<y.L&&k.push(["Drawn",`${Q(y.drawn)} of ${Q(y.L)} nt (the thread ran out of room; its drawn end is its 3′ end)`]),{title:`${y.name} mRNA`,line:`${S} · ${y.nr?`${y.nr} ribosomes`:"not translated"}`,rows:k,gene:y.name,orf:y.orf,hist:{vals:re.txCount,v:g,what:"genes with mRNA"},note:"Positions count from the end of the tail they are on (cap or 3′ end, where eIF4E and Pab1 sit). The mRNA between ribosomes is drawn shorter than the nucleotides it stands for, and the tails carry the rest (Model details 6.5)."}}const c=b=>`chr ${Nn[b]}`;if(n==="l"){const[b,m,u,p]=t;if(b===0||b===3){const y=Ct(m,Math.round((u+p)/2)),A=b===0?`${c(m)}:${Q(u)}${p>u?`–${Q(p)}`:""}`:`${c(m)}:~${Q(u)}`;return{line:`${A} · ${b===0?"linker between nucleosomes":"through engaged RNA polymerase II"} · ${y.text}`,rows:[["Position",`${A} (R64)`],["Genome",y.text]],gene:y.gene?.[4],orf:y.gene?.[5],note:b===0?"Linker DNA between two nucleosomes at their measured positions (Chereji et al. 2018); 147 bp are wrapped around each.":"Position of the polymerase: midpoint of its two flanking nucleosomes."}}const x=["","rDNA repeat (chr XII), an inactive copy on nucleosomes","rDNA repeat (chr XII), a 35S rRNA gene transcribed by Pol I","","mitochondrial DNA (in a nucleoid)"][b];return{line:x,rows:[["DNA",x]]}}if(n==="x"){const[,b,m]=t,u=Ct(b,m),p=u.gene?u.gene[4]:"an unannotated region";return{title:"Nascent mRNA",line:`being transcribed from ${p} by Pol II at ${c(b)}:~${Q(m)}`,rows:[["Gene",p],["Pol II at",`${c(b)}:~${Q(m)} · ${u.text}`]],gene:u.gene?.[4],orf:u.gene?.[5]}}const[,d,f,w]=t;if(d===0){const b=Ct(f,w),m=hi.filter((u,p)=>o.info.flags>>2+p&1);return{line:`${c(f)}:${Q(w)} (dyad) · ${b.text}${m.length?` · ${m.join(", ")}`:""}`,rows:[["Dyad",`${c(f)}:${Q(w)} (R64; Chereji et al. 2018)`],["Genome",b.text],["Histone marks",m.length?m.join(", "):"none enriched (Weiner et al. 2015, log2 ≥ 0.5)"]],gene:b.gene?.[4],orf:b.gene?.[5]}}if(d===1){const b=h(f),m=`ribosome ${w} of ${b.nr} from the 5′ end`;return{line:`translating ${b.name} mRNA · ${m}`,rows:[["Translating",`${b.name} mRNA${b.name!==b.orf?` (${b.orf})`:""}`],["On the mRNA",m],l(b)],gene:b.name,orf:b.orf}}if(d===2){const b=Ct(f,w),m=b.gene?b.gene[4]:null;return{line:`${m?`transcribing ${m}`:"engaged"} · ${c(f)}:~${Q(w)}`,rows:[["Position",`${c(f)}:~${Q(w)}`],["Genome",b.text]],gene:m,orf:b.gene?.[5],note:"Position: midpoint of the two nucleosomes flanking the polymerase; near gene ends the gene named can be the neighbour of the one transcribed."}}return{line:"rDNA repeat (chr XII), an inactive copy",rows:[["DNA","rDNA repeat (chr XII), inactive copy on nucleosomes"]]}}function En(o){const e=o>>>30;if(o===0)return{kind:"none",title:"Outside the cell"};if(e===1&&R.simMode&&C?.eng){const t=o&536870911,n=C.eng.sys;if(t>=n.n)return{kind:"none",title:""};const s=K.types[n.type[t]],r=C.tracked.includes(t);return{kind:"simbody",id:o,b:t,title:s.gene,sub:`${s.name} · ${n.mobile[t]?n.mem2d[t]?"in its membrane":n.tether[t]?"tethered (boundary shell)":"free":"fixed"} · click to ${r?"stop tracking":"track"}`}}if(e===1){const t=H.instanceAt(o&1073741823);return t?{kind:"inst",id:o,info:t,title:K.types[t.type].gene,sub:bt[t.region]}:{kind:"none",title:""}}if(e===2){const t=o>>>1&15,n=o&1,s=Dn[t]||"membrane";return{kind:"lipid",id:o,mem:t,leaf:n,title:"Lipid bilayer",sub:`${s} · ${t===8?"phospholipid monolayer":t===5?n?"intermembrane-space leaflet":"matrix leaflet":t===4?n?"cytoplasmic leaflet":"intermembrane-space leaflet":n?"cytoplasmic leaflet":t===0?"exoplasmic (outer) leaflet":"lumenal leaflet"}`}}if(e===3){const t=o>>>15&7,n=o&32767,s=["mRNA","DNA","β-1,3-glucan fibril","Mannan chain","pre-rRNA","FG repeats","Chitin microfibril","Polyphosphate chain"][t]??"polymer",r=["nucleotide-resolution thread (schematic path)","between nucleosomes, or an active rRNA gene; 2 bp per bead","procedural, 36 glucose per bead","procedural α-mannan on a cell-wall mannoprotein","nascent 35S pre-rRNA on an active rRNA gene, 2 nt per bead (schematic path)","disordered domain of an FG nucleoporin (NSP1, NUP49, NUP57, NUP159), one bead per residue","procedural, in the rim of a bud scar; 36 GlcNAc per bead (~18-chain microfibril)","procedural, in the vacuole lumen (300 mM phosphate units, assumption); 4 phosphates per bead"][t]??"";return{kind:"free",id:o,fk:t,chain:n,slot:o>>>18&4095,title:s,sub:r}}return o>=512&&o<528?{kind:"band",id:o,title:"Membrane (cut face)",sub:Dn[o-512]}:o>=256&&o<272?{kind:"label",id:o,title:hs[o-256],sub:"compartment (not loaded at this distance, or beyond the slab)"}:{kind:"none",title:""}}function Dt(o){st=o,o&&bn(),R.selId=o&&o.kind==="inst"?o.id>>>0:0,qt(),D=!0}function qt(){const o=v("inspect"),e=st;if(!e||e.kind==="none"){o.style.display="none";return}o.style.display="block";let t='<button class="close" id="inspClose">×</button>';const n=Do(e),s=n?n.rows.map(([i,h])=>`<div>${i}</div><div>${h}</div>`).join(""):"",r=n?.gene?K.types.findIndex(i=>i.gene===n.gene&&!i.id.startsWith("CPX")):-1,a=n?`${n.note?`<p class="sub">${n.note}</p>`:""}
      ${r>=0?`<div class="grid2"><button class="small" id="anProt" style="grid-column: 1 / 3">Highlight the ${n.gene} protein</button></div>`:""}
      ${n.orf?`<div style="margin-top:9px"><a href="https://www.yeastgenome.org/locus/${n.orf}" target="_blank">SGD: ${n.gene||n.orf}</a></div>`:""}`:"";if(e.kind==="inst"){const i=K.types[e.info.type],h=i.src==="PDB",l=i.nb,c=i.nnt;t+=`<h2>${i.gene}</h2><div>${i.name}</div><div style="margin-top:6px"><span class="badge ${h?"pdb":"af"}">${Io(i)}</span>
      ${i.oligomer>1?`<span class="badge" style="background:#666">${i.oligomer}-mer</span>`:""}</div>
      <div class="kv">
      ${i.genes&&i.genes!==i.gene?`<div>Subunits</div><div style="word-break:break-word">${i.genes.split(/;\s*/).slice(0,40).join(", ")}${i.genes.split(";").length>40?" …":""}</div>`:""}
      ${i.uniprot?`<div>UniProt</div><div>${i.uniprot}</div>`:""}
      <div>Location</div><div>${bt[e.info.region]}</div>
      ${s}
      <div>Copies / cell</div><div>${i.copies.toLocaleString("en-US")}</div>
      <div>Mass</div><div>${(i.mass/1e3).toFixed(1)} kDa</div>
      <div>Residues drawn</div><div>${(l-c).toLocaleString("en-US")} aa${c?` + ${c.toLocaleString("en-US")} nt`:""}</div>
      <div>Radius</div><div>${i.br.toFixed(1)} nm (Rg ${i.rg.toFixed(1)} nm)</div>
      ${tt[e.info.type]?`<div>Diffusion</div><div>${tt[e.info.type][0]>=10?tt[e.info.type][0].toFixed(0):tt[e.info.type][0].toFixed(1)} µm²/s predicted in the ${tt[e.info.type][1]} (calibrated on measurements; Model details 14.5)</div>`:""}
      ${!h&&i.plddt_mean!=null?`<div>Mean pLDDT</div><div>${i.plddt_mean.toFixed(1)} (${Math.round((i.plddt_frac_kept??1)*100)}% of residues kept, pLDDT ≥ 50)</div>`:""}
      <div>Position</div><div>${e.info.p.map(d=>d.toFixed(1)).join(", ")} nm</div>
      </div>
      ${to(K.types.map(d=>d.copies),i.copies,"molecular species")}
      <div class="grid2"><button class="small" id="loadAtoms">Load atomic model</button><button class="small" id="hlType">Highlight all</button>
      <button class="small" id="centerOn">Center view</button><button class="small" id="hideType">${ue[e.info.type]&2?"Show":"Hide"} this type</button>
      ${K.netNb[e.info.type]?`<button class="small${qe===e.info.type?" on":""}" id="netNb" style="grid-column: 1 / 3" title="Types whose STRING network neighbourhoods (personalized PageRank profiles) correlate most with this one">Network neighbours (${K.netNb[e.info.type].length})</button>`:""}</div>
      ${qe===e.info.type?`<div class="sub" style="margin-top:6px">${K.netNb[e.info.type].map(([d,f])=>`${K.types[d].gene.replace(/^CPX_\d+_/,"").replace(/_/g," ")} <span style="opacity:.7">${f.toFixed(2)}</span>`).join(" · ")}</div>`:""}
      <div style="margin-top:9px">${i.uniprot?`<a href="https://www.uniprot.org/uniprotkb/${i.uniprot}" target="_blank">UniProt</a>`:""}
      ${h?`<a href="https://www.rcsb.org/structure/${i.pdb}" target="_blank">RCSB ${i.pdb}</a>`:""}
      ${(i.model||i.uniprot)&&!h?`<a href="https://alphafold.ebi.ac.uk/entry/${i.model||i.uniprot}" target="_blank">AlphaFold DB</a>`:""}
      ${/^[A-Z0-9-]+$/.test(i.gene)&&!sn[i.id]&&!i.id.startsWith("CPX")?`<a href="https://www.yeastgenome.org/locus/${i.gene}" target="_blank">SGD</a>`:""}</div>${a}`}else if(t+=`<h2>${n?.title??e.title}</h2><div class="sub" style="font-size:12.5px">${n?.title?e.title+" · ":""}${e.sub||""}</div>`,n&&(t+=`<div class="kv">${s}</div>${n.hist?to(n.hist.vals,n.hist.v,n.hist.what):""}${a}`),e.kind==="lipid"){const i=Ro(e.mem,e.leaf);t+=`<div class="kv">${i.map((h,l)=>[h,l]).filter(([h])=>h>=.05).sort((h,l)=>l[0]-h[0]).map(([h,l])=>`<div>${ei[l]}</div><div>${h.toFixed(1)} mol%</div>`).join("")}</div>
        <p class="sub">Lipids are generated in your browser on the mesh of this membrane, at 0.65 nm² per lipid, and removed where they would overlap protein. Composition of this membrane${e.mem===0?" and leaflet":""} from isolated organelle membranes (Zinser et al. 1991)${e.mem===0?"; PS and most PE face the cytoplasm, sphingolipids the outside; ergosterol 40 mol% (assumption)":""}. Ergosterol is drawn as two beads.</p>`}o.innerHTML=t,v("inspClose").onclick=()=>Dt(null),r>=0&&(v("anProt").onclick=()=>at(r)),e.kind==="inst"&&(v("loadAtoms").onclick=()=>pi(e.info.type),v("hlType").onclick=()=>at(e.info.type),K.netNb[e.info.type]&&(v("netNb").onclick=()=>Ut(qe===e.info.type?-1:e.info.type)),v("centerOn").onclick=()=>{const i=B.state();i.target=e.info.p,B.flyTo(i,700)},v("hideType").onclick=()=>{ue[e.info.type]^=2,W.device.queue.writeBuffer(W.typeFlags,0,ue),qt(),D=!0})}let Rt=null,ct=null;function mi(){if(ct)return ct;const o="https://cdn.jsdelivr.net/npm/molstar@4.18.0/build/viewer/",e=t=>new Promise((n,s)=>{const r=document.createElement("link");r.rel="stylesheet",r.href=t+"molstar.css",document.head.appendChild(r);const a=document.createElement("script");a.src=t+"molstar.js",a.onload=()=>n(),a.onerror=()=>{r.remove(),a.remove(),s(new Error("Mol* failed to load"))},document.head.appendChild(a)});return ct=(window.__YC_EMBED__?e(o):e("molstar/").catch(()=>e(o))).catch(t=>{throw ct=null,t}),ct}async function pi(o){const e=K.types[o];v("molModal").style.display="flex",v("molTitle").textContent=`${e.gene} · ${e.name} · loading…`;try{await mi();const t=window.molstar;Rt?await Rt.plugin.clear():Rt=await t.Viewer.create("molHost",{layoutIsExpanded:!1,layoutShowControls:!1,layoutShowRemoteState:!1,layoutShowSequence:!0,layoutShowLog:!1,layoutShowLeftPanel:!1,viewportShowExpand:!1,viewportShowSelectionMode:!1,viewportShowAnimation:!1,collapseLeftPanel:!0,pdbProvider:"rcsb",emdbProvider:"rcsb"});let n="",s="";if(e.src==="pred")n=s=vn+`models/${e.id}.cif`;else if(e.src==="PDB"&&e.pdb){n=`https://files.rcsb.org/download/${e.pdb}-assembly1.cif`,s=`https://www.rcsb.org/structure/${e.pdb}`;try{(await fetch(n,{method:"HEAD"})).ok||(n=`https://files.rcsb.org/download/${e.pdb}.cif`)}catch{n=`https://files.rcsb.org/download/${e.pdb}.cif`}}else{const a=e.model||e.uniprot;s=`https://alphafold.ebi.ac.uk/entry/${a}`;try{const h=await(await fetch(`https://alphafold.ebi.ac.uk/api/prediction/${a}`)).json();n=(Array.isArray(h)?h[0]:h)?.cifUrl}catch{}n||(n=e.model?`https://alphafold.ebi.ac.uk/files/${e.model}-model_v1.cif`:`https://alphafold.ebi.ac.uk/files/AF-${e.uniprot}-F1-model_v6.cif`)}await Rt.loadStructureFromUrl(n,"mmcif",!1),v("molTitle").textContent=`${e.gene} · ${e.name} · ${e.src==="pred"?"predicted with Boltz-2 for this model (no modified nucleosides)":e.src==="PDB"?"PDB "+e.pdb+" (biological assembly)":e.model?`AlphaFold DB predicted ${e.oligomer}-mer (full length, incl. low-pLDDT regions)`:"AlphaFold DB model (full length, incl. low-pLDDT regions)"}`;const r=v("molLink");r.href=s,r.textContent=e.src==="pred"?"Download the model (mmCIF)":"Open source page"}catch(t){v("molTitle").textContent=`Could not load the atomic model: ${t?.message||t}. Internet access to RCSB/AlphaFold DB is needed.`}}let no=!1;function Lo(){no||(no=!0,Xe("Loading organelle surfaces (6.7 MB)…"),He("organelle_meshes.bin.gz").then(o=>{W.uploadOrganelleMeshes(o),D=!0}))}let mn=()=>{},pn=()=>{},tn=!1;function zo(){const o=R.memColors&&!R.organelle;v("modeSection").classList.toggle("on",!R.organelle&&!o),v("modeOrg").classList.toggle("on",R.organelle),v("modeMem").classList.toggle("on",o)}function fi(o){R.organelle=o,o&&Lo(),zo(),D=!0}function Ke(o){o!=="membranes"&&R.memColors&&mn(!1),fi(o==="organelles"),o==="membranes"&&mn(!0)}function Go(o){Eo=o;const e={low:[2.6,.75],med:[1.6,1],high:[1,1],ultra:[.6,1],auto:[1.6,1]},[t,n]=e[o];R.l0px=t,R.renderScale=n,W.resize(n),D=!0}function fn(){const o=v("legend"),e=(n,s)=>`<div><span class="sw" style="background:${n}"></span>${s}</div>`,t=R.colorMode;if(t===0)o.innerHTML=[e(ee[0],"cytosolic protein"),e(ee[1],"nuclear protein"),e("#EE9BC0","ribosome protein"),e("#EF9A7E","pre-ribosome"),e("#C9558F","rRNA"),e("#E27FA8","tRNA"),e("#B0437A","mRNA"),e("#F4B26A","histone"),e("#E07636","DNA"),e(ee[23],"nuclear pore"),e("#7CC6C2","FG repeats"),e("#E0893A","actin filament"),e("#D9534F","Arp2/3 complex"),e("#8E6CB0","capping protein"),e("#2FA4A0","tropomyosin"),e("#74C476","microtubule"),e("#2E8B57","γ-tubulin ring"),e(ee[14],"membrane protein"),e(ee[5],"mitochondrial"),e(ee[3],"ER lumen"),e("#F3C643","lipid head"),e("#F8E3A0","lipid tail"),e("#CDA866","β-glucan"),e("#C3C27A","mannan"),e(ee[9],"wall protein"),e("#8B88C9","large complex")].join("");else if(t===1)o.innerHTML=e("#4C8CD9","experimental (PDB)")+e("#F5A84C","AlphaFold model");else if(t>=10&&t<=16){const n=["H2A.Z (Htz1)","H3K4me3","H3K36me3","H3K79me3","H3K56ac","H3K9ac","H4K16ac"][t-10];o.innerHTML=e("#D9335A",`nucleosome enriched for ${n}`)+e("#CCC7B8","other nucleosome")+e("#E6E6E3","everything else")+'<div class="sub" style="grid-column: 1 / 3"><span>Per measured nucleosome: the mark of the nearest unstressed MNase nucleosome within 40 bp (Weiner et al. 2015), enriched at ≥ 0.5 log2. Annotation only: the marks do not change the structure.</span></div>'}else if(t===2)o.innerHTML=e("#0053D6","pLDDT > 90")+e("#65CBF3","70–90")+e("#FFDB13","50–70")+e("#B8B8B8","experimental (PDB)");else if(t===3)o.innerHTML='<div style="grid-column: 1 / 3"><div style="height:8px;border-radius:3px;background:linear-gradient(90deg,#440154,#3b528b,#21918c,#5ec962,#fde725)"></div><div style="display:flex;justify-content:space-between"><span>1</span><span>10³</span><span>10⁶ copies</span></div></div>';else if(t===6){o.innerHTML=K.netOn?'<div class="sub" style="grid-column: 1 / 3"><span>Similar colours: similar STRING network neighbourhoods (personalized PageRank profiles). Click a molecule, then <b>Network neighbours</b> to highlight its closest ones.</span></div>'+(qe>=0?`<div style="grid-column: 1 / 3"><span>Highlighted: ${K.types[qe].gene} and ${(K.netNb[qe]||[]).length} neighbours <a href="#" id="netClear" style="color:#8fc3ff">clear</a></span></div>`:"")+e("#D1D1CC","no network data"):'<div style="grid-column: 1 / 3">No network profiles in this build</div>';const n=o.querySelector("#netClear");n&&(n.onclick=s=>{s.preventDefault(),Ut(-1)})}else t===5?o.innerHTML=[0,1,2,3,4,5,6,7,8,9,14,15,18,19].map(n=>e(ee[n],bt[n])).join(""):o.innerHTML=""}function Ie(){const o=R.plane[3];v("planeD").value=String(Math.round(o)),v("planeDv").textContent=`${Number(o.toFixed(Ne<.5?2:1))} nm`;const e=R.slab;v("slab").value=String(Math.round(Math.log(e/5)/Math.log(4400/5)*1e3)),v("slabv").textContent=Bo(e),["axX","axY","axZ"].forEach((t,n)=>v(t).classList.toggle("on",Math.abs(R.plane[n])>.999)),gn=!0}function bi(){const o=K.types.length;let e=0;for(const g of K.types)e+=g.copies;if(v("subCount").textContent=`${(e/1e6).toFixed(1)} M macromolecules of ${o.toLocaleString("en-US")} species, each a real PDB or AlphaFold structure`,window.__YC_PREVIEW__){const g=$e.nInst.reduce((k,P)=>k+P,0);v("subCount").innerHTML=`<b>Preview:</b> ${$e.n} of 6,416 tiles (${(g/1e3).toFixed(0)} k molecules) from the full ${(e/1e6).toFixed(1)} M-molecule cell. The full cell is in yeast_cell_viewer.zip.`,v("modeOrg").disabled=!0,v("modeOrg").title="Organelle surfaces are only in the full app",v("modeMem").disabled=!0,v("modeMem").title="Membrane surfaces are only in the full app"}const t=v("views");Ue.forEach((g,k)=>{const P=document.createElement("button");P.textContent=`${(k+1)%10===0&&k<10?0:k<9?k+1:""}${k<10?" · ":""}${g.name}`,P.title=g.caption,P.onclick=()=>{je(),it(k),bn()},t.appendChild(P)}),v("tourBtn").onclick=()=>{te?je():(_o(),bn())};const n=v("search"),s=v("results");let r=0,a=[],i=[];const h=g=>{g<a.length?(at(a[g]),Fo(a[g])):di(i[g-a.length]),a=[],i=[],l(),n.value="",n.blur()},l=()=>{s.innerHTML=a.map((g,k)=>{const P=K.types[g];return`<div data-k="${k}" class="${k===r?"act":""}"><span class="nm"><b>${P.gene}</b> ${P.name}</span><span class="n">${Be(P.copies)}</span></div>`}).join("")+i.map((g,k)=>{const P=$t[g],T=a.length+k;return`<div data-k="${T}" class="${T===r?"act":""}"><span class="nm"><b>${g}</b> gene locus · chr ${Nn[P[3]]}:${Math.round(P[4]/1e3)} kb</span><span class="n">DNA</span></div>`}).join(""),s.querySelectorAll("div").forEach(g=>g.addEventListener("click",()=>h(Number(g.dataset.k))))};n.addEventListener("input",()=>{a=li(n.value).slice(0,10),i=ci(n.value),r=0,l()}),n.addEventListener("keydown",g=>{g.key==="ArrowDown"?(r=Math.min(a.length+i.length-1,r+1),l(),g.preventDefault()):g.key==="ArrowUp"?(r=Math.max(0,r-1),l(),g.preventDefault()):g.key==="Enter"&&a.length+i.length?h(r):g.key==="Escape"&&(a=[],i=[],l(),n.blur()),g.stopPropagation()}),v("cutOn").onchange=g=>{R.cutOn=g.target.checked,D=!0};const c=g=>()=>{const k=[0,0,0];k[g]=oe(B.fwd,[g===0?1:0,g===1?1:0,g===2?1:0])>=0?1:-1,R.cutOn=!0,v("cutOn").checked=!0,Me(k,oe(k,B.target))};v("axX").onclick=c(0),v("axY").onclick=c(1),v("axZ").onclick=c(2),v("alignView").onclick=()=>{R.cutOn=!0,v("cutOn").checked=!0,Ye()},v("followView").onclick=()=>{ze=!ze,v("followView").classList.toggle("on",ze),ze&&Ye()},v("planeD").oninput=g=>{const k=Number(g.target.value),P=Ce(),T=k-R.plane[3];B.target=J(B.target,O(P,T)),Me(P,k)};const d=v("nudgeStep");[...d.options].some(g=>Number(g.value)===Ne)&&(d.value=String(Ne)),d.onchange=()=>{Ne=Number(d.value),Ie(),d.blur();try{localStorage.setItem("yc.nudgeStep",String(Ne))}catch{}};const f=(g,k)=>{const P=v(g);let T=0,E=0;const I=()=>{clearTimeout(T),clearInterval(E)};P.addEventListener("pointerdown",F=>{F.preventDefault(),I(),Pe=null,v("sweepBtn").classList.remove("on"),k(),T=window.setTimeout(()=>{E=window.setInterval(k,60)},400)});for(const F of["pointerup","pointerleave","pointercancel"])P.addEventListener(F,I)};f("nudgeBack",()=>nt(-Ne)),f("nudgeFwd",()=>nt(Ne));const w=()=>Math.max(2,R.slab*.25),b=()=>{R.cutOn||(R.cutOn=!0,v("cutOn").checked=!0)};f("pillBack",()=>{b(),nt(-w())}),f("pillFwd",()=>{b(),nt(w())}),v("pillCut").onclick=()=>{R.cutOn=!R.cutOn,v("cutOn").checked=R.cutOn,D=!0},v("slab").oninput=g=>{const k=Number(g.target.value)/1e3;R.slab=5*Math.pow(4400/5,k),Ie(),D=!0},v("sweepBtn").onclick=()=>{if(Pe){Pe=null,v("sweepBtn").classList.remove("on");return}const g=Ce();R.cutOn=!0,v("cutOn").checked=!0;const k=Oe.G[0];Pe={t0:performance.now(),d0:-k,d1:k,dur:3e4},B.target=J(B.target,O(g,-k-R.plane[3])),v("sweepBtn").classList.add("on")},v("flipBtn").onclick=()=>{const g=Ce(),k=B.state(),P=B.target,T=O(g,-1);Me(T,-R.plane[3]),k.fwd=O(k.fwd,-1),k.target=P,B.flyTo(k,1200)},v("modeSection").onclick=()=>Ke("molecules"),v("modeOrg").onclick=()=>Ke("organelles"),v("modeMem").onclick=()=>Ke("membranes"),v("colorMode").onchange=g=>{R.colorMode=Number(g.target.value),fn(),D=!0},v("shaded").onchange=g=>{R.shaded=g.target.checked,D=!0},v("motion").onchange=g=>{R.motion=Number(g.target.value),g.target.blur(),D=!0},v("outlines").onchange=g=>{R.outline=g.target.checked?1:0,D=!0},fn();const m=v("showList"),u=g=>k=>gi(ot[g],!k),p=[["Actin","#E0893A",u("actin")],["Microtubules","#74C476",u("microtubules")],["Eisosomes","#3F8FCB",u("eisosomes")],["Ribosomes","#EE9BC0",u("ribosomes")],["tRNA","#E27FA8",u("tRNA")],["Nucleosomes","#F4B26A",u("nucleosomes")],["Nuclear pores",ee[23],u("pores")],["RNA polymerases","#9A72B8",u("polymerases")],["mRNA","#B0437A",g=>Ge(0,g)],["DNA","#E07636",g=>Ge(1,g)],["pre-rRNA","#C9558F",g=>Ge(4,g)],["FG repeats","#7CC6C2",g=>Ge(5,g)],["β-Glucan","#CDA866",g=>Ge(2,g)],["Mannan","#C3C27A",g=>Ge(3,g)],["Chitin",ee[54],g=>Ge(6,g)],["Polyphosphate",ee[55],g=>Ge(7,g)],["Lipids","#F3C643",g=>{R.showLipids=g,D=!0}],["Membranes","#F8E3A0",g=>{R.showMeshes=g,D=!0}],["Soluble proteins",ee[0],g=>Tt([0,1,2,3,4,5,6,7,8,10,11,12,13,24,25,26],g)],["Membrane proteins",ee[14],g=>Tt([14,15,16,17,18,19,20,21],g)]],x=[],y=bt.map((g,k)=>k),A=p.findIndex(([g])=>g==="Membranes"),M=()=>{R.memColors=x.every((g,k)=>g.checked===(k===A))&&!window.__YC_PREVIEW__,R.memColors&&Lo(),zo(),D=!0},_=g=>{p.forEach(([,,P],T)=>{x[T].checked=g(T),P(g(T))});const k=p.every((P,T)=>g(T));Tt(y,k),M()};p.forEach(([g,k,P],T)=>{const E=document.createElement("label");E.className="chk",E.innerHTML=`<input type="checkbox" checked><span class="sw" style="background:${k}"></span><span class="nm">${g}</span><a class="only" title="Show only ${g.toLowerCase()}">only</a>`;const I=E.querySelector("input");x.push(I),I.onchange=()=>{P(I.checked),M()},E.querySelector(".only").onclick=F=>{F.preventDefault(),F.stopPropagation(),_(N=>N===T)},m.appendChild(E)}),mn=g=>_(g?(k=>k===A):(()=>!0)),pn=g=>_(g?(k=>g.includes(p[k][0])):(()=>!0)),v("showAll").onclick=()=>_(()=>!0),v("hideAll").onclick=()=>_(()=>!1);const S=v("regionList");bt.forEach((g,k)=>{const P=document.createElement("label");P.className="chk",P.innerHTML=`<input type="checkbox" checked data-r="${k}"><span class="sw" style="background:${ee[k]}"></span>${g}`,P.querySelector("input").onchange=T=>Tt([k],T.target.checked),S.appendChild(P)}),v("quality").onchange=g=>Go(g.target.value),v("detailOn").onchange=g=>{if(H.detailEnabled=g.target.checked,!H.detailEnabled)for(const k of H.loaded.values())k.det&&H.dropDetail(k);D=!0},v("shotBtn").onclick=Uo,v("helpBtn").onclick=()=>v("helpModal").style.display="flex",v("aboutBtn").onclick=()=>v("aboutModal").style.display="flex",v("aboutBody").innerHTML=ti,v("modelBtn").onclick=()=>{v("modelBody").innerHTML||(v("modelBody").innerHTML=ni(oi)),v("modelModal").style.display="flex"},v("aboutModal").addEventListener("click",g=>{g.target.id==="aboutToModel"&&(g.preventDefault(),v("aboutModal").style.display="none",v("modelBtn").click())}),document.querySelectorAll(".modal").forEach(g=>g.addEventListener("click",k=>{(k.target===g||k.target.classList.contains("close"))&&(g.style.display="none")})),v("molClose").onclick=()=>v("molModal").style.display="none",v("panelToggle").onclick=()=>mt(v("panel").classList.contains("collapsed")),v("sheetClose").onclick=()=>mt(!1),ft.addEventListener("change",()=>mt(!ft.matches)),mt(!ft.matches&&!matchMedia("(max-height: 500px)").matches),v("caption").onclick=()=>wt(null),Ie(),_n(),wi()}const ft=matchMedia("(max-width: 760px)");function mt(o){v("panel").classList.toggle("collapsed",!o),v("panelToggle").textContent=ft.matches?"☰":o?"‹":"›"}function bn(){ft.matches&&mt(!1)}function gi(o,e){for(const t of o){const n=K.byId.get(t);n!==void 0&&(e?ue[n]|=2:ue[n]&=-3)}W.device.queue.writeBuffer(W.typeFlags,0,ue),D=!0}function Ge(o,e){e?R.freeHide&=~(1<<o):R.freeHide|=1<<o,D=!0}function Tt(o,e){for(const t of o)e?R.regionMask|=1<<t:R.regionMask&=~(1<<t);R.regionMask>>>=0,document.querySelectorAll("#regionList input").forEach(t=>{t.checked=(R.regionMask>>>Number(t.dataset.r)&1)===1}),D=!0}function _n(){let o=null;const e=new Set,t=new Map;let n=null,s=!1,r={t:0,x:0,y:0};const a=()=>{const[c,d]=[...t.values()],f=q.getBoundingClientRect();return{d:Math.hypot(d.x-c.x,d.y-c.y),mx:(c.x+d.x)/2-f.left,my:(c.y+d.y)/2-f.top,a:Math.atan2(d.y-c.y,d.x-c.x)}};q.addEventListener("contextmenu",c=>c.preventDefault()),q.addEventListener("pointerdown",c=>{if(q.setPointerCapture(c.pointerId),q.focus(),te&&je(),c.pointerType==="touch"){if(t.set(c.pointerId,{x:c.clientX,y:c.clientY}),t.size===2){o=null,s=!0,n=a();return}if(t.size>2)return;s=!1}o={x:c.clientX,y:c.clientY,btn:c.button,moved:0,shift:c.shiftKey},q.classList.add("dragging")}),q.addEventListener("pointermove",c=>{if(c.pointerType==="touch"){if(!t.has(c.pointerId))return;if(t.set(c.pointerId,{x:c.clientX,y:c.clientY}),t.size>=2&&n){const w=a();w.d>0&&n.d>0&&l(w.mx,w.my,n.d/w.d),B.pan(w.mx-n.mx,w.my-n.my,B.pxPerNmAtUnit(q.clientHeight)/B.dist);let b=w.a-n.a;b>Math.PI?b-=2*Math.PI:b<-Math.PI&&(b+=2*Math.PI),B.roll(-b),n=w,ze&&Ye(),B.anim=null,D=!0;return}}if(!o){(c.pointerType==="mouse"||c.pointerType==="pen")&&yi(c.clientX,c.clientY);return}const d=c.clientX-o.x,f=c.clientY-o.y;if(o.x=c.clientX,o.y=c.clientY,o.moved+=Math.abs(d)+Math.abs(f),o.btn===2||o.btn===1||o.shift){const w=B.pxPerNmAtUnit(q.clientHeight)/B.dist;B.pan(d,f,w)}else B.orbit(d*.005,f*.005),ze&&Ye();B.anim=null,D=!0});const i=c=>{if(c.pointerType==="touch"){if(!t.delete(c.pointerId))return;if(t.size===1){const[w]=[...t.values()];n=null,o={x:w.x,y:w.y,btn:0,moved:99,shift:!1};return}if(t.size>1){n=a();return}n=null}if(!o)return;const d=c.pointerType==="touch",f=o.moved<(d?10:4)&&o.btn===0&&!(d&&s);if(o=null,q.classList.remove("dragging"),!(c.type==="pointercancel"||!f)){if(d){const w=performance.now();if(w-r.t<320&&Math.hypot(c.clientX-r.x,c.clientY-r.y)<30){r.t=0,h(c.offsetX,c.offsetY);return}r={t:w,x:c.clientX,y:c.clientY}}if(zt&&!C.active){Zn(c.offsetX,c.offsetY).then(w=>{if(w){const b=B.state();b.target=w,B.flyTo(b,500)}});return}if(R.simMode&&C?.active){vt(c.offsetX,c.offsetY).then(w=>{w.id>>>30===1&&(C.toggleTrack(w.id&536870911),D=!0)});return}vt(c.offsetX,c.offsetY).then(w=>{const b=En(w.id);Dt(b.kind==="inst"||b.kind==="lipid"||b.kind==="free"||b.kind==="band"?b:null)})}};q.addEventListener("pointerup",i),q.addEventListener("pointercancel",i);const h=(c,d)=>{Zn(c,d).then(f=>{if(!f)return;const w=B.state();w.target=f,w.dist=Math.max(B.minDist*4,B.dist*.55),B.flyTo(w,900)})};q.addEventListener("dblclick",c=>{t.size||h(c.offsetX,c.offsetY)});const l=(c,d,f)=>{const w=q.clientWidth/q.clientHeight,b=B.ray(c/q.clientWidth,d/q.clientHeight,w),m=J(B.eye,O(b,B.dist/oe(b,B.fwd))),u=B.dist;B.zoom(f);const p=1-B.dist/u;B.target=J(B.target,O(Fe(m,B.target),p))};q.addEventListener("wheel",c=>{c.preventDefault();const d=c.deltaY||c.deltaX,f=c.deltaMode===1?d*30:d;if(c.altKey)nt(Math.sign(f)*Ne*Math.min(4,Math.abs(f)/100));else if(c.shiftKey){const w=Ce(),b=Math.max(1,Math.min(R.slab,B.dist*.05))*Math.sign(f)*Math.min(4,Math.abs(f)/60+.5);B.target=J(B.target,O(w,b)),Me(w,R.plane[3]+b)}else l(c.offsetX,c.offsetY,Math.exp(f*.0015*(c.metaKey?.125:1)));B.anim=null,te&&je(),D=!0},{passive:!1}),window.addEventListener("keydown",c=>{if(c.target.tagName==="INPUT"||c.target.tagName==="SELECT")return;const d=c.key.toLowerCase();if(e.add(d),d==="c")R.cutOn=!R.cutOn,v("cutOn").checked=R.cutOn,D=!0;else if(d==="v")R.cutOn=!0,v("cutOn").checked=!0,Ye();else if(d==="[")R.slab=Math.max(5,R.slab/1.25),Ie(),D=!0;else if(d==="]")R.slab=Math.min(4400,R.slab*1.25),Ie(),D=!0;else if(d===","||d==="."||c.code==="Comma"||c.code==="Period"){const f=d==="."||c.code==="Period";nt((f?1:-1)*(c.altKey?Ne:Math.max(2,R.slab*.25)))}else if(d==="o")Ke(R.organelle?"molecules":"organelles");else if(d==="m")Ke(R.memColors&&!R.organelle?"molecules":"membranes");else if(d==="t")te?je():_o();else if(d==="p")Uo();else if(d==="escape")Dt(null),at(-1);else if(/^[0-9]$/.test(d)){const f=d==="0"?9:Number(d)-1;f<Ue.length&&(je(),it(f))}}),window.addEventListener("keyup",c=>e.delete(c.key.toLowerCase())),window.addEventListener("blur",()=>e.clear()),_n.keys=e,window.addEventListener("resize",()=>{W.resize(R.renderScale),D=!0})}let oo=0,nn=!1,Ve=null;function yi(o,e){const t=performance.now();nn||t-oo<90||(oo=t,nn=!0,D=!0,vt(o,e).then(n=>{nn=!1,Oo(En(n.id),o,e)}))}function Oo(o,e,t){const n=v("hover");if(!o.title||o.kind==="none"){n.style.display="none",Ve=null;return}Ve={d:o,x:e,y:t};const s=Do(o);let r=`<b>${s?.title??o.title}</b>`;if(o.kind==="inst"){const a=K.types[o.info.type];r+=` · ${a.name} <span style="color:#aaa">· ${o.sub}</span>`}else s||(r+=` <span style="color:#aaa">· ${o.sub}</span>`);s&&(r+=`<br><span style="color:#ccc">${s.line}</span>`),n.innerHTML=r,n.style.display="block",n.style.left=`${e+14}px`,n.style.top=`${t+12}px`}un=()=>{Ve&&v("hover").style.display==="block"&&Oo(Ve.d,Ve.x,Ve.y),st&&(st.kind==="free"||st.kind==="inst")&&qt()};q.addEventListener("pointerleave",()=>{v("hover").style.display="none",Ve=null});let gn=!0,so=0,De=null;const dt={data:null};function wi(){const o=v("minimap");o.addEventListener("click",e=>{if(!De)return;const t=o.getBoundingClientRect(),n=((e.clientX-t.left)/t.width*2-1)*De.ext,s=(1-(e.clientY-t.top)/t.height*2)*De.ext,r=J(De.c,J(O(De.u,n),O(De.v,s))),a=B.state();a.target=r,B.flyTo(a,1200),te&&je()})}function yn(o){const e=v("minimap"),t=e.getContext("2d"),n=e.width,s=e.height,r=R.cutOn?Ce():B.fwd,a=R.cutOn?R.plane[3]:oe(B.fwd,B.target);if(gn&&o-so>250){gn=!1,so=o;let u=Fe(B.right,O(r,oe(B.right,r)));Je(u)<.001&&(u=he(Le(r,[0,0,1]))),u=he(u);let p=Le(r,u);oe(p,B.up)<0&&(p=O(p,-1));const x=2300,y=O(r,a);De={u,v:p,c:y,ext:x};const A=120,M=t.createImageData(A,A),_=[ee[80],ee[9],"#D2BE8E",ee[83],ee[3],ee[3],ee[1],ee[2],ee[4],ee[6],ee[5],"#E9C46A",ee[7],ee[8],"#C9C2DC"].map(S=>{const g=parseInt(S.slice(1),16);return[g>>16&255,g>>8&255,g&255]});_[3]=[196,214,235];for(let S=0;S<A;S++)for(let g=0;g<A;g++){const k=((g+.5)/A*2-1)*x,P=(1-(S+.5)/A*2)*x,T=J(y,J(O(u,k),O(p,P))),E=No.label(T[0],T[1],T[2]),I=_[E],F=4*(S*A+g);M.data[F]=I[0],M.data[F+1]=I[1],M.data[F+2]=I[2],M.data[F+3]=E===0?0:255}dt.data=M}if(!dt.data||!De)return;const i=yn.tmp||(yn.tmp=document.createElement("canvas"));i.width=dt.data.width,i.height=dt.data.height,i.getContext("2d").putImageData(dt.data,0,0),t.clearRect(0,0,n,s),t.imageSmoothingEnabled=!0,t.drawImage(i,0,0,n,s);const h=De,l=u=>{const p=Fe(u,h.c);return[(oe(p,h.u)/h.ext+1)/2*n,(1-oe(p,h.v)/h.ext)/2*s]},c=q.clientWidth/q.clientHeight,d=B.dist*Math.tan(B.fovy/2),f=d*c,w=[[-1,-1],[1,-1],[1,1],[-1,1]].map(([u,p])=>l(J(B.target,J(O(B.right,u*f),O(B.up,p*d)))));t.strokeStyle="#FF9400",t.lineWidth=3,t.beginPath(),w.forEach(([u,p],x)=>x?t.lineTo(u,p):t.moveTo(u,p)),t.closePath(),t.stroke();const[b,m]=l(B.target);t.fillStyle="#FF9400",t.beginPath(),t.arc(b,m,5,0,7),t.fill(),v("mmLabel").textContent=R.cutOn?`plane ${Math.round(a)} nm`:"view plane"}function Uo(){const o=R.renderScale;W.resize(2*o),R.renderScale=2*o,qo(!0),q.toBlob(e=>{if(e){const t=document.createElement("a");t.href=URL.createObjectURL(e),t.download=`yeast_cell_${new Date().toISOString().replace(/[:.]/g,"-").slice(0,19)}.png`,t.click(),setTimeout(()=>URL.revokeObjectURL(t.href),5e3)}R.renderScale=o,W.resize(o),D=!0},"image/png"),Xe("Screenshot saved (2× resolution)")}function vi(){const o=H.totals(),e=W.stats,t=ke.length>3?1e3/(ke.reduce((l,c)=>l+c)/ke.length):0,n=H.pendingForView();v("hud").innerHTML=`<b>${Be(o.inst)}</b> molecules in ${o.tiles} tiles${n?` · <span style="color:#FF9400">streaming ${n}</span>`:""}${H.starved?' · <span style="color:#FF9400">memory budget reached (thinner slab loads more)</span>':""}<br>
    drawn: <b>${Be(e.ninst)}</b> in view · ${Be(e.l0*64)} residue beads (≤) · ${Be(e.l1*16)} LOD1 · ${Be(e.l2)} LOD2 · ${Be(e.lip)} lipids · ${Be(e.fre)} polymer beads<br>
    ${t?`${t.toFixed(0)} fps · `:""}LOD bead ≥ ${R.l0px.toFixed(1)} px · GPU pools ${W.gpuMemoryMB().toFixed(0)} MB · ${(H.bytesLoaded/1048576).toFixed(0)} MB streamed${W.adapterInfo?" · "+W.adapterInfo:""}`;const s=W.width/q.clientWidth,r=B.pxPerNmAtUnit(W.height)/B.dist/s,a=120/r,i=Math.pow(10,Math.floor(Math.log10(a))),h=[1,2,5,10].map(l=>l*i).reduce((l,c)=>Math.abs(c-a)<Math.abs(l-a)?c:l);v("scalebar").querySelector(".bar").style.width=`${h*r}px`,v("scalebar").querySelector("span").textContent=Bo(h)}let $o=0,Lt=!1;function qo(o=!1){const e=performance.now();W.frame(B,R,Oe.G,H.meshDraws(),o||e-W.lastStatsT>400),W.lastStatsT=performance.now(),Lt=!0,W.device.queue.onSubmittedWorkDone().then(()=>{$o=performance.now()-e,Lt=!1})}let io=0;function Ho(){C?.active||(H.update(B,q.clientWidth/q.clientHeight,R.cutOn,R.plane,R.slab,R.organelle||R.memColors,W.height),Tn=performance.now(),(H.inflight||H.detInflight)&&(D=!0))}setInterval(()=>{H&&performance.now()-Tn>150&&Ho()},150);function Wo(o){if(requestAnimationFrame(Wo),zt&&C&&!C.active&&C.box){const n=C.box.center,s=B.target;Math.abs(n[0]-s[0])+Math.abs(n[1]-s[1])+Math.abs(n[2]-s[2])>.01&&(C.showBox([...s],C.box.size),D=!0)}const e=Zt?o-Zt:16;Zt=o;const t=_n.keys;if(t&&t.size){const n=B.dist*.9*e/1e3;let s=[0,0,0];t.has("w")&&(s=J(s,O(B.fwd,n))),t.has("s")&&(s=J(s,O(B.fwd,-n))),t.has("a")&&(s=J(s,O(B.right,-n))),t.has("d")&&(s=J(s,O(B.right,n))),t.has("q")&&(s=J(s,O(B.up,-n))),t.has("e")&&(s=J(s,O(B.up,n))),t.has("z")&&B.roll(-e*.0015),t.has("x")&&B.roll(e*.0015),(Je(s)>0||t.has("z")||t.has("x"))&&(B.fly(s),D=!0,ze&&Ye())}if(B.step(o)&&(D=!0,ze&&!te&&Ye()),Pe){const n=Math.min(1,(o-Pe.t0)/Pe.dur),s=Ce(),r=Pe.d0+(Pe.d1-Pe.d0)*n;B.target=J(B.target,O(s,r-R.plane[3])),R.plane[3]=r,Ie(),n>=1&&(Pe=null,v("sweepBtn").classList.remove("on")),D=!0}if(ri(o),R.motion>0&&(D=!0),C?.active&&(C.running&&C.replayIdx<0||C.playing||C.computing)&&(D=!0),o-Tn>90&&Ho(),D&&!Lt){if(D=!1,qo(),ke.push(e),ke.length>30&&ke.shift(),Eo==="auto"&&ke.length>=20){const n=ke.reduce((s,r)=>s+r)/ke.length;n>45&&R.l0px<3.5?(R.l0px=Math.min(3.5,R.l0px*1.15),ke.length=0):n<22&&R.l0px>.9&&(R.l0px=Math.max(.9,R.l0px/1.1),ke.length=0)}setTimeout(()=>{D=D||W.statsBusy},0)}else D||(ke.length=0);o-io>300&&(io=o,vi(),yn(o),v("pillCut").classList.toggle("on",R.cutOn),Se>=0&&H.inflight&&Bn())}ii();let ht=!0,zt=!1;function xi(){const o=()=>Number(v("simSize").value),e=()=>v("simMode").value,t=m=>{v("simStats").innerHTML=m};v("simSize").oninput=()=>{v("simSizev").textContent=`${o()} nm`,C.box&&!C.active&&C.showBox(C.box.center,o()).then(()=>{D=!0})},v("simSize").onchange=()=>{C.active&&d(!0)};const n=m=>{zt=m,v("simPlace").textContent=m?"Done placing":"Place box",v("simPlace").classList.toggle("on",m),m?t("The box sits at the centre of the view: pan, orbit and zoom as usual, or click a point to bring it to the centre. Then press Start."):C.active||t("")};v("simPlace").onclick=async()=>{if(zt){n(!1);return}if(C.active)return;await C.showBox([...B.target],o());const m=B.state();m.dist=Math.min(Math.max(m.dist,1.6*o()),4*o()),B.flyTo(m,500),n(!0),D=!0};const s=()=>{const m=C.set;v("simDt").value=String(Math.log10(m.dt)),v("simDtv").textContent=`${m.dt<.1?m.dt.toFixed(3):m.dt.toFixed(2)} ns`;const u=Cn(m.mode,m.eps);v("simDt").max=String(Math.log10(u)),v("simDt").title=`Brownian time step; at most ${u<.1?u.toFixed(3):u.toFixed(2)} ns for this model and stickiness (larger steps overestimate diffusion)`,v("simSpf").value=String(m.stepsPerFrame),v("simSpfv").textContent=`${m.stepsPerFrame}`,v("simEpsv").textContent=`${m.eps.toFixed(2)} kT`,v("simViscv").textContent=`${m.viscosity}×`,v("simPlay").textContent=C.running?"Pause":"Run live",v("simPlayRec").textContent=C.playing?"Stop playback":"Play recording",v("simAhead").textContent=C.computing?"Cancel":"Compute ahead"},r=()=>{const m=v("simEps"),u=wo[e()];m.max=String(u),Number(m.value)>u&&(m.value=String(u)),m.title=`Depth of a short-range attraction between molecular surfaces (0: excluded volume only); at most ${u} kT in this model`};r(),v("simMode").onchange=()=>{r(),v("simLen").value=e()==="residue"?"200":e()==="mid"?"1000":"5000",C.active&&d(!0)};const a=()=>({elec:v("simElec").checked,ionic:Number(v("simIonic").value),specific:v("simSpec").checked,epsGo:Number(v("simEpsGo").value),mrna:v("simMrna").checked}),i=()=>{const m=a();v("simIonicv").textContent=`${m.ionic.toFixed(2)} M`,v("simEpsGov").textContent=`${m.epsGo} kT`};for(const m of["simElec","simIonic","simSpec","simEpsGo"])v(m).oninput=()=>{i(),C.applySettings(a())};v("simMrna").onchange=()=>{C.active&&d(!0)},v("simBoundary").onchange=()=>{C.active&&d(!0)},v("simMemFluid").onchange=()=>{C.active&&d(!0)},v("simMemCal").onchange=()=>{C.active&&d(!0)},v("simHow").onchange=()=>{C.active&&d(!0)},v("simLen").value="5000";const h=m=>{const u=v("simProg");u.style.display=m===null?"none":"block",m!==null&&(u.firstElementChild.style.width=`${(100*m).toFixed(1)}%`)},l=async()=>{const m=Number(v("simLen").value),u=C.eng;h(0);const p=C.computeAhead(m,300,x=>{h(x),D=!0});s(),await p,C.eng===u&&(h(null),C.active&&C.frames.length>1&&C.play(!0),s(),D=!0)};v("simAhead").onclick=()=>{C.computing?C.abort=!0:l()},v("simPlayRec").onclick=()=>{C.play(!C.playing),s(),D=!0};let c=!1;const d=async(m=!1)=>{if(c||C.active&&!m)return;c=!0,n(!1);const u=!C.active;C.active&&(C.abort=!0,C.stop(),h(null)),await C.showBox(C.box?C.box.center:[...B.target],o()),C.applySettings((await Ko(async()=>{const{defaultSettings:p}=await Promise.resolve().then(()=>gs);return{defaultSettings:p}},void 0,import.meta.url)).defaultSettings(e())),C.applySettings({eps:Number(v("simEps").value),viscosity:Number(v("simVisc").value),hydro:v("simHI").checked,boundary:v("simBoundary").value,memFluid:v("simMemFluid").checked,memCalibrated:v("simMemCal").checked,...a()}),v("simStart").disabled=!0;try{if(await C.start(t),u){ht=R.cutOn,R.cutOn=!1,v("cutOn").checked=!1;const p=B.state();p.target=[...C.box.center],p.dist=2.2*C.box.size,B.flyTo(p,800)}R.simMode=!0,v("simRun").style.display="",v("simStart").textContent="Restart",s(),v("simHow").value==="ahead"&&(C.running=!1,l())}catch(p){t(`<span style="color:#c0392b">${p.message??p}</span>`),u||(R.simMode=!1,v("simRun").style.display="none",R.cutOn=ht,v("cutOn").checked=ht)}v("simStart").disabled=!1,c=!1,D=!0};v("simStart").onclick=()=>d(C.active),v("simDt").oninput=()=>{C.applySettings({dt:Math.pow(10,Number(v("simDt").value))}),s()},v("simSpf").oninput=()=>{C.applySettings({stepsPerFrame:Number(v("simSpf").value)}),s()},v("simEps").oninput=()=>{C.applySettings({eps:Number(v("simEps").value)}),s()},v("simHI").onchange=()=>{C.applySettings({hydro:v("simHI").checked})},v("simVisc").oninput=()=>{C.applySettings({viscosity:Number(v("simVisc").value)}),s()},v("simPlay").onclick=()=>{C.computing||(C.running=!C.running,C.running&&C.resumeLive(),s(),D=!0)},v("simLive").onclick=()=>{C.computing||(C.resumeLive(),C.running=!0,s(),D=!0)},v("simScrub").oninput=()=>{C.computing||(C.running=!1,C.replay(Number(v("simScrub").value)),s(),D=!0)},v("simExport").onclick=()=>C.export(),v("simFps").oninput=()=>{const m=Number(v("simFps").value);C.setPlayFps(m),v("simFpsv").textContent=`${m} fps`},v("simTrackClear").onclick=()=>{C.setTracked([]),D=!0},v("simExit").onclick=()=>{C.stop(),C.hideBox(),R.simMode=!1,v("simRun").style.display="none",t(""),h(null),s(),v("simStart").textContent="Start",n(!1),R.cutOn=ht,v("cutOn").checked=ht,D=!0};const f=m=>m<1e3?`${m.toFixed(1)} ns`:m<1e6?`${(m/1e3).toFixed(2)} µs`:`${(m/1e6).toFixed(3)} ms`,w=()=>{const m=v("simTracks"),u=C.active&&C.tracked.length>0;if(v("simTrackClear").style.display=u?"":"none",!u||!C.eng){m.innerHTML="";return}const p=C.eng.sys;m.innerHTML=C.trackStats().map((x,y)=>`<div><span style="display:inline-block;width:10px;height:10px;border-radius:5px;background:${So[y]};margin-right:5px"></span><b>${K.types[p.type[x.b]].gene}</b>
      moved ${x.disp.toFixed(1)} nm${Number.isFinite(x.d)?` · ${x.d.toFixed(1)} µm²/s from successive frames`:""}</div>`).join("")};C.onUpdate=()=>{w();const m=v("simScrub");m.max=String(Math.max(0,C.frames.length-1)),C.replayIdx<0?m.value=m.max:m.value=String(C.replayIdx)};let b=0;setInterval(()=>{if(!C.active||!C.eng||!C.meta)return;const m=C.meta,u=C.eng;++b%4===0&&C.pollGrid();const p=C.grid,x=p&&p.dropped>0?` · <span style="color:#c0392b">neighbour grid overflow: ${(p.dropped*100).toFixed(2)}% of spheres missed (fullest cell ${p.max} > ${p.cap})</span>`:"",y=C.replayIdx>=0?C.frames[C.replayIdx].tNs:u.timeNs;v("simScrubv").textContent=f(y);const A=C.diffusion().map(g=>`${g.label}: <b>${g.d.toFixed(2)}</b> µm²/s (in solvent ${g.d0.toFixed(1)}; n=${g.n}; lag ${f(g.lagNs)})`).join("<br>"),M=C.computing||C.playing&&C.aheadRate?C.aheadRate:C.stepsPerSec*C.set.dt/1e3,_=Number(v("simLen").value)/1e3,S=C.computing?`Computing ahead (no drawing waits): ${M.toFixed(2)} µs per second, about ${Math.max(0,(_-(u.timeNs-C.frames[0].tNs)/1e3)/Math.max(M,1e-6)).toFixed(0)} s left<br>`:C.playing?`Playing the recording (${C.frames.length} frames, ${f(C.frames[C.frames.length-1].tNs-C.frames[0].tNs)})<br>`:"";t(`${S}<b>${f(y)}</b> simulated · ${M.toFixed(2)} µs per second · ${u.sys.n} molecules
      (${m.nMobile} free${m.nMem2d?` (${m.nMem2d} of them in membranes)`:""}${m.nShell?`, ${m.nShell} tethered in the boundary shell`:""}, ${m.nFixed} fixed: ${C.set.memFluid?"":"membrane proteins, "}filaments, chromatin, nuclear pores)
      · boundary <b>${m.boundary==="periodic"?"periodic":m.boundary==="walls"?"closed walls":"embedded"}</b>${m.boundary==="periodic"&&m.hasMembrane?' <span style="color:#c0392b">(membranes do not continue across the faces)</span>':""}<br>
      Interactions: ${C.set.elec?m.charges?`electrostatics on ${m.nCharged.toLocaleString()} charged spheres (box net charge ${m.netCharge.toFixed(0)} e; Debye length ${(.303/Math.sqrt(C.set.ionic)).toFixed(2)} nm)`:"no charge data in this build":"electrostatics off"}
      · ${C.set.specific?m.contacts?`specific binding: ${m.nNativePairs} predicted complex${m.nNativePairs===1?"":"es"} can form between ${m.nNativeBodies} molecules (${C.set.epsGo} kT per interface)`:"no contact data in this build":"specific binding off"}
      · ${C.set.mrna?`${m.nLinks} mRNA link${m.nLinks===1?"":"s"} joining ${m.nPolyRibo} polysome ribosomes`:"mRNA links off"}<br>
      <span style="color:#8a8f98">${u.nSph.toLocaleString()} collision spheres · grid ${u.gridInfo.h.toFixed(1)} nm cells (capacity ${u.gridInfo.cap})${u.gridInfo.large?` + ${u.gridInfo.hB.toFixed(1)} nm cells for the ${u.gridInfo.large.toLocaleString()} spheres over ${u.gridInfo.split.toFixed(2)} nm`:""}</span><br>
      Volume fraction ${(u.sys.phi*100).toFixed(1)}% (hydrated)${x}${C.set.hydro?` · mean-field hydrodynamics: D × ${Gt(u.sys.phi).toFixed(2)}, D<sub>r</sub> × ${Ot(u.sys.phi).toFixed(2)}`:""}<br>
      Diffusion (time-averaged MSD): ${A||"…"}<br><span style="color:#8a8f98">Measured in yeast: GFP 8.6 µm²/s (fission-yeast cytoplasm, FRAP);
      budding yeast 8.7–21 µm²/s (FCS, lab-dependent); 40 nm GEM particles about 0.3 µm²/s (particle tracking). Excluded volume,
      electrostatics and mean-field hydrodynamics leave GFP-sized molecules about 2.7× too fast; the default Viscosity of
      2.75× is that calibration (fitted on training papers, checked on held-out ones; Model details, section 14.5).
      ${vo[C.set.mode]}</span>`)},500)}
