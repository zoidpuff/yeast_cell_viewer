import { Renderer, FrameState, ID_INST, ID_LIPID, ID_FREE } from './renderer';
import { SimController, TRACK_COLORS } from './sim/controller';
import { hiTrans, hiRot, dtMax, STICKY_NOTE, STICKY_MAX } from './sim/params';
import { Camera, CamState } from './camera';
import { TileManager, LoadedTile } from './tiles';
import { Geometry } from './geometry';
import { loadLibrary, loadTileIndex, fetchJson, fetchBin, Library, TileIndex, REGION_NAMES, MEM_NAMES, LABEL_NAMES, PAL, SPECIAL, b64ToBytes, DATA_BASE } from './data';
import { lipidComposition, LIPID_CLASSES } from './lipids';
import { loadFields } from './fields';
import { V3, dot, norm, cross, sub, add, scale, len } from './math';
import { ABOUT_HTML } from './about';
import { markdownToHtml } from './markdown';
import MODEL_MD from '../../../MODEL.md?raw';

const $ = <T extends HTMLElement = HTMLElement>(id: string) => document.getElementById(id) as T;
const fmt = (n: number) => n >= 1e6 ? (n / 1e6).toFixed(n >= 1e7 ? 1 : 2) + ' M' : n >= 1e4 ? (n / 1e3).toFixed(0) + ' k' : n.toLocaleString('en-US');
const fmtNm = (n: number) => n >= 1000 ? (n / 1000).toFixed(n >= 10000 ? 0 : 2).replace(/\.?0+$/, '') + ' µm' : (n >= 10 ? n.toFixed(0) : n.toFixed(1)) + ' nm';

declare global { interface Window { yc: any; } }

interface Bookmark { name: string; caption: string; target: V3; normal: V3; up: V3; width: number; slab: number; show?: string[] }

let lib: Library, idx: TileIndex, geo: any, geom: Geometry, bookmarks: Bookmark[], typesTiles: { off: Uint32Array; tiles: Uint16Array };
let R: Renderer, T: TileManager;
let sim: SimController;
const cam = new Camera();
const canvas = $<HTMLCanvasElement>('gpu');

const fs: FrameState = {
  cutOn: true, plane: [0, 0, 1, 0], slab: 60, colorMode: 0, regionMask: 0xFFFFFFFF, freeHide: 0, highlight: false,
  selId: 0, shaded: true, organelle: false, motion: 0, l0px: 1.6, l1px: 2.5, minpx: 0.3, outline: 1, showMeshes: true, showLipids: true, renderScale: 1, memColors: false,
};
let planeFollow = false;
let dirty = true, lastUpdate = 0, lastFrameT = 0, quality = 'auto';
let typeFlags: Uint32Array;
// molecule classes with their own Show switch (library type ids)
const CLASS_TYPES: Record<string, string[]> = {
  actin: ['F_actin', 'Arp23', 'CP_cap', 'TPM1_coil'],
  microtubules: ['MT_tubulin', 'gTuRC_nuc', 'gTuRC_cyto'],
  ribosomes: ['ribosome_80S', 'ribosome_80S_RAC', 'ribosome_80S_Ssh1', 'ribosome_80S_Sec61', 'mitoribosome_74S_IM', 'CPX_1602_54S_mitochondrial_large_ribosoma',
    'CPX_1603_37S_mitochondrial_small_ribosoma', 'preribo_90S', 'preribo_60S_nucleolar',
    'preribo_60S_nucleoplasm', 'preribo_60S_cytoplasm', 'preribo_40S'],
  tRNA: ['tRNA'],
  nucleosomes: ['nucleosome'],
  pores: ['NPC'],
  polymerases: ['PolI', 'PolII', 'PolIII'],
  eisosomes: ['Pil1_dimer', 'Lsp1_dimer'],
};
// connected assemblies (filament subunits, pores, chromatin): illustrative motion would pull them apart
const NO_MOTION = [...CLASS_TYPES.actin, ...CLASS_TYPES.microtubules, ...CLASS_TYPES.eisosomes, 'NPC', 'nucleosome'];
let selType = -1;
let selection: any = null;
let tour: { i: number; phase: 'fly' | 'wait' | 'hold'; t: number } | null = null;
let sweep: { t0: number; d0: number; d1: number; dur: number } | null = null;
const frameTimes: number[] = [];

function setLoading(msg: string, f: number) {
  ($('loading').querySelector('.m') as HTMLElement).textContent = msg;
  ($('loading').querySelector('.pb div') as HTMLElement).style.width = `${Math.round(f * 100)}%`;
}
function fail(msg: string) {
  const e = $('loading').querySelector('.err') as HTMLElement;
  e.style.display = 'block'; e.innerHTML = msg;
  ($('loading').querySelector('.t') as HTMLElement).textContent = 'Could not start the viewer';
}
function toast(msg: string, ms = 2200) {
  const t = $('toast'); t.textContent = msg; t.style.display = 'block';
  clearTimeout((toast as any)._t); (toast as any)._t = setTimeout(() => (t.style.display = 'none'), ms);
}

// ---------------------------------------------------------------- boot
async function boot() {
  if (!('gpu' in navigator)) {
    fail(`This browser does not expose <b>WebGPU</b>. Use a current Chrome or Edge (version 113 or later), or Safari 26 or later.
      On Linux Chrome, enable <code>chrome://flags/#enable-unsafe-webgpu</code> and <code>#enable-vulkan</code>.
      Firefox 141+ supports WebGPU on Windows.<br><br>Also check that the page is served over <code>http://localhost</code>
      (for example <code>python -m http.server</code>). Opening the file directly with <code>file://</code> does not work.`);
    return;
  }
  try {
    setLoading('Loading cell geometry', 0.01);
    [geo, bookmarks] = await Promise.all([fetchJson('geometry.json'), fetchJson('bookmarks.json')]);
    geo.F = await loadFields(geo, fetchBin);      // real cell shapes (null for the analytic cell)
    geom = new Geometry(geo);
    idx = await loadTileIndex();
    lib = await loadLibrary(setLoading);
    CLASS_TYPES.tRNA = lib.types.filter(t => t.id === 'tRNA' || t.id.startsWith('tRNA_')).map(t => t.id);   // one type per isoacceptor
    const tt = await fetchBin('types_tiles.bin');
    const NT = lib.types.length;
    typesTiles = { off: new Uint32Array(tt, 0, NT + 1), tiles: new Uint16Array(tt, (NT + 1) * 4) };
    setLoading('Starting WebGPU', 0.7);
    R = new Renderer(canvas);
    R.onLost = (m) => { toast(m, 10000); };
    await R.init(lib, geo);
    typeFlags = new Uint32Array(NT);
    for (const ids of Object.values(CLASS_TYPES)) for (const id of ids) { const i = lib.byId.get(id); if (i !== undefined) typeFlags[i] |= 4; }
    for (const id of NO_MOTION) { const i = lib.byId.get(id); if (i !== undefined) typeFlags[i] |= 8; }
    R.device.queue.writeBuffer(R.typeFlags, 0, typeFlags);
    setLoading('Compiling shaders and starting tile streaming', 0.9);
    let embedded: Record<string, Uint8Array> | undefined;
    if (window.__YC_EMBED__) {
      embedded = {};
      for (const [k, v] of Object.entries(window.__YC_EMBED__)) if (k.startsWith('tiles/') || k.startsWith('library/')) embedded[k] = b64ToBytes(v);
    }
    const nW = Math.max(2, Math.min(4, (navigator.hardwareConcurrency || 4) - 1));
    T = new TileManager(idx, R, lib, geo, nW, embedded);
    T.onChange = () => { dirty = true; };
    sim = new SimController(R, T, lib);
    buildUI();
    setupSim();
    const qv = new URLSearchParams(location.search).get('view');
    const v0 = qv === null ? 0 : (isNaN(Number(qv)) ? Math.max(0, bookmarks.findIndex(b => b.name.toLowerCase().startsWith(qv.toLowerCase()))) : Math.min(bookmarks.length - 1, Math.max(0, Number(qv))));
    applyBookmark(v0, false);
    $('loading').style.display = 'none';
    requestAnimationFrame(loop);
    window.yc = { R, T, cam, fs, lib, idx, sim, applyBookmark, setType, setNetwork, pickAt, bookmarks, setQuality, describe, ready: true, state: () => ({ pending: T.pendingForView(), detailPending: T.detInflight, loaded: T.loaded.size, totals: T.totals(), stats: R.stats, cam: { t: cam.target, d: cam.dist }, plane: fs.plane, gpuMs, gpuBusy, dirty }), render: () => { dirty = true; } };
  } catch (e: any) {
    console.error(e);
    fail(`${String(e?.message || e)}<br><br>See the browser console for details.`);
  }
}

// ---------------------------------------------------------------- plane helpers
function planeN(): V3 { return [fs.plane[0], fs.plane[1], fs.plane[2]]; }
function setPlane(n: V3, d: number) { n = norm(n); fs.plane = [n[0], n[1], n[2], d]; syncPlaneUI(); dirty = true; }
function alignToView() { const n = cam.fwd; setPlane(n, dot(n, cam.target)); }
// move the plane (and the camera with it) along its normal; nudgeStep is the user's fine step in nm
let nudgeStep = 0.5;
try { const s = Number(localStorage.getItem('yc.nudgeStep')); if (s > 0) nudgeStep = s; } catch { /* storage unavailable */ }
function movePlane(s: number) { const n = planeN(); cam.target = add(cam.target, scale(n, s)); setPlane(n, fs.plane[3] + s); }

// ---------------------------------------------------------------- bookmarks / tour
function bookmarkState(b: Bookmark): CamState {
  const aspect = canvas.clientWidth / canvas.clientHeight;
  const n = norm(b.normal);
  let up = norm(sub(b.up, scale(n, dot(b.up, n))));
  return { target: [...b.target] as V3, dist: cam.distForWidth(b.width, aspect), fwd: n, up };
}
function applyBookmark(i: number, animate = true, onDone?: () => void) {
  const b = bookmarks[i];
  if (!b) return;
  if ((fs.organelle || fs.memColors) && i !== 0) setDisplayMode('molecules');
  if (b.show) { setShowOnly(b.show); showRestricted = true; } else if (showRestricted) { setShowOnly(null); showRestricted = false; }
  const st = bookmarkState(b);
  const n1 = norm(b.normal), d1 = dot(n1, b.target);
  const n0 = planeN(), d0 = fs.plane[3], s0 = fs.slab;
  fs.cutOn = true; ($('cutOn') as HTMLInputElement).checked = true;
  planeFollow = false; $('followView').classList.remove('on');
  if (!animate) { cam.set(st); setPlane(n1, d1); fs.slab = b.slab; syncPlaneUI(); dirty = true; onDone?.(); return; }
  cam.flyTo(st, 2600, () => { setPlane(n1, d1); fs.slab = b.slab; syncPlaneUI(); onDone?.(); }, (e) => {
    let n = norm(add(scale(n0, 1 - e), scale(n1, e)));
    if (len(add(scale(n0, 1 - e), scale(n1, e))) < 1e-3) n = n1;
    fs.plane = [n[0], n[1], n[2], d0 + (d1 - d0) * e];
    fs.slab = Math.exp(Math.log(s0) + (Math.log(b.slab) - Math.log(s0)) * e);
    syncPlaneUI();
  });
  showCaption(b);
  dirty = true;
}
function showCaption(b: Bookmark | null) {
  const c = $('caption');
  if (!b) { c.style.display = 'none'; return; }
  c.innerHTML = `<b>${b.name}</b>${b.caption}`; c.style.display = 'block';
  clearTimeout((showCaption as any)._t);
  if (!tour) (showCaption as any)._t = setTimeout(() => (c.style.display = 'none'), 9000);
}
function startTour() {
  tour = { i: 0, phase: 'fly', t: performance.now() };
  $('tourBtn').textContent = '■ Stop'; $('tourBtn').classList.add('on');
  applyBookmark(0, true, () => { if (tour) { tour.phase = 'wait'; tour.t = performance.now(); } });
}
function stopTour() {
  tour = null; $('tourBtn').textContent = '▶ Tour'; $('tourBtn').classList.remove('on');
  showCaption(null);
}
function tourStep(now: number) {
  if (!tour) return;
  if (tour.phase === 'wait' && (T.pendingForView() === 0 || now - tour.t > 14000)) { tour.phase = 'hold'; tour.t = now; }
  else if (tour.phase === 'hold' && now - tour.t > 6500) {
    tour.i = (tour.i + 1) % bookmarks.length; tour.phase = 'fly';
    const i = tour.i;
    applyBookmark(i, true, () => { if (tour && tour.i === i) { tour.phase = 'wait'; tour.t = performance.now(); } });
  }
}

// ---------------------------------------------------------------- selection / search
let selNet = -1;
/** highlight a type and its strongest network neighbours (personalized PageRank profiles); -1 clears */
function setNetwork(t: number) {
  if (selType >= 0) { typeFlags[selType] &= ~1; selType = -1; renderSelType(); }
  selNet = t;
  const on = new Set<number>(t >= 0 ? [t, ...(lib.netNb[t] || []).map(x => x[0])] : []);
  for (let i = 0; i < typeFlags.length; i++) { if (on.has(i)) typeFlags[i] |= 1; else typeFlags[i] &= ~1; }
  fs.highlight = t >= 0;
  R.device.queue.writeBuffer(R.typeFlags, 0, typeFlags);
  renderLegend(); if (selection) renderInspect();
  dirty = true;
}
function setType(t: number) {
  if (selNet >= 0) setNetwork(-1);
  if (selType >= 0) typeFlags[selType] &= ~1;
  selType = t;
  if (t >= 0) typeFlags[t] |= 1;
  fs.highlight = t >= 0;
  R.device.queue.writeBuffer(R.typeFlags, 0, typeFlags);
  renderSelType();
  dirty = true;
}
function renderSelType() {
  const el = $('selType');
  if (selType < 0) { el.style.display = 'none'; return; }
  const t = lib.types[selType];
  const loaded = T.countLoaded(selType);
  el.style.display = 'block';
  el.innerHTML = `<div><b>${t.gene}</b> ${t.name}</div>
    <div class="sub" style="margin:3px 0 6px">${fmt(t.copies)} copies in the cell · ${fmt(loaded)} in loaded tiles · ${srcLabel(t)}</div>
    <div class="grid2"><button class="small" id="flyNearest">Fly to nearest</button><button class="small" id="clearType">Clear highlight</button></div>`;
  $('flyNearest').onclick = () => flyToNearest(selType);
  $('clearType').onclick = () => setType(-1);
}
function searchTypes(q: string) {
  q = q.trim().toLowerCase();
  if (!q) return [];
  const out: { i: number; s: number }[] = [];
  lib.types.forEach((t, i) => {
    const g = t.gene.toLowerCase(), id = t.id.toLowerCase(), u = (t.uniprot || '').toLowerCase(), nm = (t.name || '').toLowerCase(), gs = (t.genes || '').toLowerCase();
    let s = 0;
    if (g === q || id === q || u === q) s = 100;
    else if (g.startsWith(q) || id.startsWith(q)) s = 60;
    else if (gs.split(/[;,\s]+/).includes(q)) s = 55;
    else if (nm.includes(q)) s = 30 + (nm.startsWith(q) ? 10 : 0);
    else if (gs.includes(q) || u.startsWith(q) || (t.pdb || '').toLowerCase() === q) s = 20;
    if (s) out.push({ i, s: s + Math.log10(t.copies + 1) });
  });
  out.sort((a, b) => b.s - a.s);
  return out.slice(0, 14).map(x => x.i);
}
function flyToNearest(type: number) {
  const o0 = typesTiles.off[type], o1 = typesTiles.off[type + 1];
  if (o1 <= o0) { toast('No copies of this molecule were placed.'); return; }
  let best = -1, bd = 1e18;
  for (let k = o0; k < o1; k++) {
    const ti = typesTiles.tiles[k], c = T.tileCenter(ti);
    const d = len(sub(c, cam.target));
    if (d < bd) { bd = d; best = ti; }
  }
  const c = T.tileCenter(best);
  const t = lib.types[type];
  const width = Math.max(90, t.br * 9);
  const goTo = (tile: LoadedTile) => {
    const k = T.nearestInTile(tile, type, c);
    if (k < 0) return;
    const d = T.decode(tile, k);
    const n = planeN();
    const st = cam.state(); st.target = d.p; st.dist = cam.distForWidth(width, canvas.clientWidth / canvas.clientHeight);
    cam.flyTo(st, 900, () => { setPlane(n, dot(n, d.p) - t.br * 0.4); });
    fs.plane[3] = dot(n, d.p) - t.br * 0.4;
    const pid = T.poolIndex(tile, k);
    select({ kind: 'inst', id: (ID_INST | pid) >>> 0, info: T.instanceAt(pid) });
  };
  const loaded = T.loaded.get(best);
  const n = planeN();
  const st = cam.state(); st.target = c; st.dist = cam.distForWidth(width * 2.5, canvas.clientWidth / canvas.clientHeight);
  fs.cutOn = true; ($('cutOn') as HTMLInputElement).checked = true;
  if (fs.organelle || fs.memColors) setDisplayMode('molecules');
  cam.flyTo(st, 2200, () => {
    setPlane(n, dot(n, c) - 20);
    const L = T.loaded.get(best);
    if (L) goTo(L);
    else {
      const h = (tile: LoadedTile) => { if (tile.idx === best) { T.onTileLoaded.splice(T.onTileLoaded.indexOf(h), 1); goTo(tile); } };
      T.onTileLoaded.push(h);
    }
  });
  setPlane(n, dot(n, c) - 20);
  if (loaded) { /* handled on arrival of flight */ }
}

// ---------------------------------------------------------------- picking & inspection
async function pickAt(cx: number, cy: number) {
  const dpr = R.width / canvas.clientWidth;
  dirty = true;   // the pick is read back from the next frame (touch has no hover to trigger one)
  return R.pick(cx * dpr, cy * dpr);
}
/** world position of the surface drawn at a canvas pixel (null: background) */
async function worldAtPixel(x: number, y: number): Promise<V3 | null> {
  const r = await pickAt(x, y);
  if (r.depth <= 0) return null;
  const zv = cam.near() / r.depth;   // view-space distance along -z
  const ray = cam.ray(x / canvas.clientWidth, y / canvas.clientHeight, canvas.clientWidth / canvas.clientHeight);
  return add(cam.eye, scale(ray, zv / dot(ray, cam.fwd)));
}
function srcLabel(t: Library['types'][number]) {
  return t.src === 'PDB' ? 'PDB ' + t.pdb : t.src === 'AF-multimer' ? 'AlphaFold complex' : t.src === 'pred' ? 'Boltz-2 prediction' : 'AlphaFold';
}
// predicted in-vivo diffusion per species (library/diffusion.json, MODEL.md 14.5); absent in older builds
let diffPred: Record<string, [number, string]> = {};
fetchJson('library/diffusion.json').then(j => { diffPred = j.types || {}; }).catch(() => { /* not in this build */ });
// gene loci of the chromatin (genome_map.py): gene -> [x, y, z, chromosome, start, end, ORF]
let geneLoci: Record<string, [number, number, number, number, number, number, string]> = {};
fetchJson('library/gene_loci.json').then(j => { geneLoci = j.loci || {}; }).catch(() => { /* not in this build */ });
const ROMAN = ['', 'I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI'];
function searchLoci(q: string): string[] {
  q = q.trim().toLowerCase();
  if (q.length < 2) return [];
  const out: string[] = [];
  for (const [g, v] of Object.entries(geneLoci)) {
    const gl = g.toLowerCase(), o = v[6].toLowerCase();
    if (gl === q || o === q) out.unshift(g); else if (gl.startsWith(q) || o.startsWith(q)) out.push(g);
    if (out.length > 30) break;
  }
  return out.slice(0, 4);
}
function flyToLocus(g: string) {
  const v = geneLoci[g]; if (!v) return;
  const p: V3 = [v[0], v[1], v[2]], n = planeN();
  const st = cam.state(); st.target = p; st.dist = cam.distForWidth(260, canvas.clientWidth / canvas.clientHeight);
  fs.cutOn = true; ($('cutOn') as HTMLInputElement).checked = true;
  cam.flyTo(st, 900, () => { setPlane(n, dot(n, p)); });
  toast(`${g} (${v[6]}) gene locus: chromosome ${ROMAN[v[3]]}, ${v[4].toLocaleString()}–${v[5].toLocaleString()} bp (centre of its nucleosomes)`, 4200);
}

// ---------------------------------------------------------------- what an mRNA, DNA or chromatin instance is
// (polymer_annot.py): library/genome.json (SGD features per chromosome, gene layouts, transcripts, annotated tiles) and
// annot/<tile>.json: m per mRNA chain [transcript, part, first, last, ribosome], l per DNA chain [code, chrom, bp0, bp1],
// x per extra chain [code, chrom, bp] (chain ids nM + nL + i, worker.ts), i [instance, code, a, b]
let genome: any = null, genomeReq = false;
const annotCache = new Map<number, any>(), annotReq = new Set<number>();
let onAnnot = () => { /* set by the hover and inspect code */ };
const MARK_NAMES = ['H2A.Z', 'H3K4me3', 'H3K36me3', 'H3K79me3', 'H3K56ac', 'H3K9ac', 'H4K16ac'];
function tileAnnot(ti: number): any {      // the tile's table; null while loading or where there is none
  if (!genome) {
    if (!genomeReq) { genomeReq = true; fetchJson('library/genome.json').then(j => { genome = j; genome.tileSet = new Set(j.tiles); genome.iTypes = new Set(j.i_types); onAnnot(); }).catch(() => { /* not in this build */ }); }
    return null;
  }
  if (!genome.tileSet.has(ti)) return null;
  const a = annotCache.get(ti);
  if (a || annotReq.has(ti)) return a ?? null;
  annotReq.add(ti);
  fetchJson(`annot/${idx.ix[ti]}_${idx.iy[ti]}_${idx.iz[ti]}.json`).then(j => {
    j.im = new Map(j.i.map((r: number[]) => [r[0], r]));
    annotCache.set(ti, j);
    if (annotCache.size > 96) { const k = annotCache.keys().next().value!; annotCache.delete(k); annotReq.delete(k); }
    onAnnot();
  }).catch(() => { /* stays unannotated */ });
  return null;
}
const fmtN = (n: number) => n.toLocaleString('en-US');
/** small histogram of copies per cell over a population (log scale), with this molecule's bin and position marked */
function copyHist(vals: number[], v: number, what: string): string {
  const xs = vals.filter(x => x > 0); if (xs.length < 5 || !(v > 0)) return '';
  const lo = 0, hi = Math.max(1, Math.ceil(Math.log10(Math.max(...xs, v)) * 2) / 2), nb = Math.round((hi - lo) * 4);     // 4 bins per decade
  const bins = new Array(nb).fill(0);
  const bi = (x: number) => Math.min(nb - 1, Math.max(0, Math.floor((Math.log10(x) - lo) / (hi - lo) * nb)));
  for (const x of xs) bins[bi(x)]++;
  const W = 300, H = 46, L = 2, bw = (W - 2 * L) / nb, mx = Math.max(...bins), mine = bi(v);
  const below = xs.filter(x => x < v).length / xs.length;
  let g = '';
  bins.forEach((c, i) => {
    const h = c ? Math.max(1.5, c / mx * (H - 4)) : 0, lo_ = Math.round(10 ** (lo + i / 4)), hi_ = Math.round(10 ** (lo + (i + 1) / 4));
    g += `<rect x="${(L + i * bw + 1).toFixed(1)}" y="${(H - h).toFixed(1)}" width="${(bw - 2).toFixed(1)}" height="${h.toFixed(1)}" rx="1.5" fill="${i === mine ? 'var(--acc)' : 'rgba(255,255,255,0.26)'}"><title>${fmtN(lo_)}–${fmtN(hi_)} copies: ${fmtN(c)} ${what}</title></rect>`;
  });
  const xm = L + Math.min(nb, Math.max(0, (Math.log10(v) - lo) / (hi - lo) * nb)) * bw;
  let ticks = '';
  for (let d = 0; d <= hi; d++) {
    const x = L + d * 4 * bw, lab = d < 3 ? String(10 ** d) : d < 6 ? `${10 ** (d - 3)}k` : `${10 ** (d - 6)}M`;
    ticks += `<line x1="${x.toFixed(1)}" x2="${x.toFixed(1)}" y1="${H}" y2="${H + 3}" stroke="rgba(255,255,255,0.35)"/><text x="${x.toFixed(1)}" y="${H + 13}" text-anchor="${d === 0 ? 'start' : d === hi ? 'end' : 'middle'}" font-size="9.5" fill="#9aa">${lab}</text>`;
  }
  const pct = Math.floor(below * 100);
  return `<div class="sub" style="margin-top:8px">Copies per cell (log scale): ${fmtN(v)}, ${pct >= 99 ? 'among the highest 1%' : pct < 1 ? 'among the lowest' : `more than ${pct}%`} of the ${fmtN(xs.length)} ${what}</div>
    <svg viewBox="0 0 ${W} ${H + 16}" width="100%" style="display:block;margin-top:3px" role="img" aria-label="Histogram of copies per cell; this one has ${fmtN(v)}">
    ${g}<line x1="${L}" x2="${W - L}" y1="${H}" y2="${H}" stroke="rgba(255,255,255,0.35)"/>
    <line x1="${xm.toFixed(1)}" x2="${xm.toFixed(1)}" y1="0" y2="${H}" stroke="var(--acc)" stroke-width="2"/>${ticks}</svg>`;
}
/** features of chromosome c (1-16) at bp, or the nearest genes around it: { text, gene (ORF or RNA gene containing bp) } */
function genomeContext(c: number, bp: number): { text: string; gene: any } {
  const F: any[] = genome.feat[c - 1] || [];
  const isGene = (f: any) => f[3] <= 6;      // ORF, tRNA, snoRNA, snRNA, ncRNA, rRNA genes, pseudogene
  const lab = (f: any) => f[4] === f[5] ? f[4] : `${f[4]} (${f[5]})`;
  const inside = F.filter(f => f[0] <= bp && bp <= f[1]).sort((a, b) => Number(isGene(b)) - Number(isGene(a)));
  if (inside.length) {
    const g = inside.find(isGene);
    const parts = inside.slice(0, 3).map(f => {
      const into = f[2] > 0 ? bp - f[0] + 1 : f[1] - bp + 1, len = f[1] - f[0] + 1;
      return `in ${lab(f)}, ${genome.ftypes[f[3]]}${isGene(f) ? `, bp ${fmtN(into)} of ${fmtN(len)} from its start` : ''}`;
    });
    return { text: parts.join('; '), gene: g || null };
  }
  let up: any = null, dn: any = null;        // nearest genes left and right of bp
  for (const f of F) { if (!isGene(f)) continue; if (f[1] < bp && (!up || f[1] > up[1])) up = f; if (f[0] > bp && (!dn || f[0] < dn[0])) dn = f; }
  const near = [up, dn].filter(Boolean).map((f: any) => {     // distance to the start of a gene whose promoter side faces bp
    const start = f[2] > 0 ? f[0] : f[1], upstream = f[2] > 0 ? bp < f[0] : bp > f[1];
    return upstream ? `${fmtN(Math.abs(start - bp))} bp upstream of ${lab(f)}` : null;
  }).filter(Boolean);
  return { text: `intergenic${up && dn ? ` between ${up[4]} and ${dn[4]}` : ''}${near.length ? ` (${near.join('; ')})` : ''}`, gene: null };
}
/** identity of a picked mRNA/DNA chain or chromatin instance: { title?, line, rows, gene? }; null while loading, undefined if none */
function annotOf(d: any): any {
  let ti = -1, rec: number[] | undefined, what = '';
  if (d.kind === 'free' && (d.fk === 0 || d.fk === 1)) {
    const t = T.tileBySlot(d.slot); if (!t) return undefined;
    ti = t.idx;
  } else if (d.kind === 'inst' && d.info?.t) {
    if (genome && !genome.iTypes.has(d.info.type)) return undefined;
    ti = d.info.t.idx;
  } else return undefined;
  const a = tileAnnot(ti);
  if (!a) return genome && !genome.tileSet.has(ti) ? undefined : null;
  const nM = a.m.length, nL = a.l.length;
  if (d.kind === 'inst') { rec = a.im.get(d.info.local); what = 'i'; }
  else if (d.fk === 1) { rec = a.l[d.chain]; what = 'l'; }
  else if (d.chain < nM) { rec = a.m[d.chain]; what = 'm'; }
  else if (a.x && d.chain >= nM + nL) { rec = a.x[d.chain - nM - nL]; what = 'x'; }
  if (!rec) return undefined;
  const tx = genome.tx;
  const txInfo = (t: number) => { const g = genome.genes[tx.g[t]]; return { name: g[0], orf: g[1], u5: g[2], cds: g[3], u3: g[4], pa: g[5], L: tx.n[t], nr: tx.r[t], drawn: tx.d[t] }; };
  const layoutRow = (x: any): [string, string] => ['Transcript', `${fmtN(x.L)} nt: 5′ UTR ${x.u5}, CDS ${fmtN(x.cds)}, 3′ UTR ${x.u3}, poly(A) ${x.pa}`];
  if (what === 'm') {
    const [t, part, p0, p1, rib] = rec, x = txInfo(t);
    const r5 = (n: number) => n < x.u5 ? '5′ UTR' : 'CDS';
    const r3 = (n: number) => n < x.pa ? 'poly(A) tail' : n < x.pa + x.u3 ? '3′ UTR' : 'CDS';
    const both = (a_: string, b_: string) => a_ === b_ ? a_ : `${a_} / ${b_}`;
    let where = '';
    if (part === 0) where = `CDS, through ribosome ${rib} of ${x.nr} (from the 5′ end)`;
    else if (part === 1) where = `${both(r3(p1), r3(p0))}, nt ${fmtN(p0 + 1)}–${fmtN(p1 + 1)} from the 3′ end`;
    else if (part === 2) where = `${both(r5(p0), r5(p1))}, nt ${fmtN(p0 + 1)}–${fmtN(p1 + 1)} from the 5′ cap`;
    else {
      const reg = (n: number) => { const from3 = x.drawn - 1 - n; return from3 < x.pa + x.u3 ? r3(from3) : r5(n); };
      where = `${both(reg(p0), reg(p1))}, nt ${fmtN(p0 + 1)}–${fmtN(p1 + 1)} of ${fmtN(x.drawn)}`;
    }
    if (!genome.txCount) { genome.txCount = new Array(genome.genes.length).fill(0); for (const gi of tx.g) genome.txCount[gi]++; }
    const nCopies = genome.txCount[tx.g[t]];
    const rows: [string, string][] = [['Gene', `${x.name}${x.name !== x.orf ? ` (${x.orf})` : ''}`], ['Copies / cell', `${fmtN(nCopies)} mRNA${nCopies === 1 ? '' : 's'} of this gene`], layoutRow(x),
      ['Ribosomes', x.nr ? `${x.nr} on this mRNA (polysome)` : 'none (an untranslated mRNA)'], ['This piece', where]];
    if (part === 3 && x.drawn < x.L) rows.push(['Drawn', `${fmtN(x.drawn)} of ${fmtN(x.L)} nt (the thread ran out of room; its drawn end is its 3′ end)`]);
    return { title: `${x.name} mRNA`, line: `${where} · ${x.nr ? `${x.nr} ribosomes` : 'not translated'}`, rows, gene: x.name, orf: x.orf,
      hist: { vals: genome.txCount, v: nCopies, what: 'genes with mRNA' },
      note: 'Positions count from the end of the tail they are on (cap or 3′ end, where eIF4E and Pab1 sit). The mRNA between ribosomes is drawn shorter than the nucleotides it stands for, and the tails carry the rest (Model details 6.5).' };
  }
  const R_ = (c: number) => `chr ${ROMAN[c]}`;
  if (what === 'l') {
    const [code, c, b0, b1] = rec;
    if (code === 0 || code === 3) {
      const ctx = genomeContext(c, Math.round((b0 + b1) / 2));
      const pos = code === 0 ? `${R_(c)}:${fmtN(b0)}${b1 > b0 ? `–${fmtN(b1)}` : ''}` : `${R_(c)}:~${fmtN(b0)}`;
      return { line: `${pos} · ${code === 0 ? 'linker between nucleosomes' : 'through engaged RNA polymerase II'} · ${ctx.text}`,
        rows: [['Position', `${pos} (R64)`], ['Genome', ctx.text]], gene: ctx.gene?.[4], orf: ctx.gene?.[5],
        note: code === 0 ? 'Linker DNA between two nucleosomes at their measured positions (Chereji et al. 2018); 147 bp are wrapped around each.' : 'Position of the polymerase: midpoint of its two flanking nucleosomes.' };
    }
    const t = ['', 'rDNA repeat (chr XII), an inactive copy on nucleosomes', 'rDNA repeat (chr XII), a 35S rRNA gene transcribed by Pol I', '', 'mitochondrial DNA (in a nucleoid)'][code];
    return { line: t, rows: [['DNA', t]] };
  }
  if (what === 'x') {
    const [, c, bp] = rec, ctx = genomeContext(c, bp);
    const g = ctx.gene ? ctx.gene[4] : 'an unannotated region';
    return { title: 'Nascent mRNA', line: `being transcribed from ${g} by Pol II at ${R_(c)}:~${fmtN(bp)}`, rows: [['Gene', g], ['Pol II at', `${R_(c)}:~${fmtN(bp)} · ${ctx.text}`]], gene: ctx.gene?.[4], orf: ctx.gene?.[5] };
  }
  const [, code, ra, rb] = rec;      // instance
  if (code === 0) {
    const ctx = genomeContext(ra, rb), marks = MARK_NAMES.filter((_, i) => (d.info.flags >> (2 + i)) & 1);
    return { line: `${R_(ra)}:${fmtN(rb)} (dyad) · ${ctx.text}${marks.length ? ` · ${marks.join(', ')}` : ''}`,
      rows: [['Dyad', `${R_(ra)}:${fmtN(rb)} (R64; Chereji et al. 2018)`], ['Genome', ctx.text], ['Histone marks', marks.length ? marks.join(', ') : 'none enriched (Weiner et al. 2015, log2 ≥ 0.5)']],
      gene: ctx.gene?.[4], orf: ctx.gene?.[5] };
  }
  if (code === 1) {
    const x = txInfo(ra), where = `ribosome ${rb} of ${x.nr} from the 5′ end`;
    return { line: `translating ${x.name} mRNA · ${where}`, rows: [['Translating', `${x.name} mRNA${x.name !== x.orf ? ` (${x.orf})` : ''}`], ['On the mRNA', where], layoutRow(x)], gene: x.name, orf: x.orf };
  }
  if (code === 2) {
    const ctx = genomeContext(ra, rb), g = ctx.gene ? ctx.gene[4] : null;
    return { line: `${g ? `transcribing ${g}` : 'engaged'} · ${R_(ra)}:~${fmtN(rb)}`, rows: [['Position', `${R_(ra)}:~${fmtN(rb)}`], ['Genome', ctx.text]], gene: g, orf: ctx.gene?.[5],
      note: 'Position: midpoint of the two nucleosomes flanking the polymerase; near gene ends the gene named can be the neighbour of the one transcribed.' };
  }
  return { line: 'rDNA repeat (chr XII), an inactive copy', rows: [['DNA', 'rDNA repeat (chr XII), inactive copy on nucleosomes']] };
}

function describe(id: number): any {
  const cat = id >>> 30;
  if (id === 0) return { kind: 'none', title: 'Outside the cell' };
  if (cat === 1 && fs.simMode && sim?.eng) {    // simulation bodies carry their body index
    const b = id & 0x1FFFFFFF, sys = sim.eng.sys;     // bit 29: drawn faded
    if (b >= sys.n) return { kind: 'none', title: '' };
    const t = lib.types[sys.type[b]], tr = sim.tracked.includes(b);
    return { kind: 'simbody', id, b, title: t.gene, sub: `${t.name} · ${sys.mobile[b] ? (sys.mem2d[b] ? 'in its membrane' : sys.tether[b] ? 'tethered (boundary shell)' : 'free') : 'fixed'} · click to ${tr ? 'stop tracking' : 'track'}` };
  }
  if (cat === 1) {
    const info = T.instanceAt(id & 0x3FFFFFFF);
    if (!info) return { kind: 'none', title: '' };
    return { kind: 'inst', id, info, title: lib.types[info.type].gene, sub: REGION_NAMES[info.region] };
  }
  if (cat === 2) {
    const mem = (id >>> 1) & 15, leaf = id & 1;
    const name = MEM_NAMES[mem] || 'membrane';
    const side = mem === 8 ? 'phospholipid monolayer' : mem === 5 ? (leaf ? 'intermembrane-space leaflet' : 'matrix leaflet')
      : mem === 4 ? (leaf ? 'cytoplasmic leaflet' : 'intermembrane-space leaflet')
      : leaf ? 'cytoplasmic leaflet' : (mem === 0 ? 'exoplasmic (outer) leaflet' : 'lumenal leaflet');
    return { kind: 'lipid', id, mem, leaf, title: 'Lipid bilayer', sub: `${name} · ${side}` };
  }
  if (cat === 3) {
    const kind = (id >>> 15) & 7, chain = id & 0x7FFF;
    const nm = ['mRNA', 'DNA', 'β-1,3-glucan fibril', 'Mannan chain', 'pre-rRNA', 'FG repeats', 'Chitin microfibril', 'Polyphosphate chain'][kind] ?? 'polymer';
    const sub = ['nucleotide-resolution thread (schematic path)', 'between nucleosomes, or an active rRNA gene; 2 bp per bead', 'procedural, 36 glucose per bead',
      'procedural α-mannan on a cell-wall mannoprotein', 'nascent 35S pre-rRNA on an active rRNA gene, 2 nt per bead (schematic path)', 'disordered domain of an FG nucleoporin (NSP1, NUP49, NUP57, NUP159), one bead per residue',
      'procedural, in the rim of a bud scar; 36 GlcNAc per bead (~18-chain microfibril)',
      'procedural, in the vacuole lumen (300 mM phosphate units, assumption); 4 phosphates per bead'][kind] ?? '';
    return { kind: 'free', id, fk: kind, chain, slot: (id >>> 18) & 0xFFF, title: nm, sub };
  }
  if (id >= 0x200 && id < 0x210) return { kind: 'band', id, title: 'Membrane (cut face)', sub: MEM_NAMES[id - 0x200] };
  if (id >= 0x100 && id < 0x110) return { kind: 'label', id, title: LABEL_NAMES[id - 0x100], sub: 'compartment (not loaded at this distance, or beyond the slab)' };
  return { kind: 'none', title: '' };
}
function select(d: any) {
  selection = d;
  if (d) closeSheet();
  fs.selId = d && d.kind === 'inst' ? d.id >>> 0 : 0;
  renderInspect();
  dirty = true;
}
function renderInspect() {
  const el = $('inspect');
  const d = selection;
  if (!d || d.kind === 'none') { el.style.display = 'none'; return; }
  el.style.display = 'block';
  let h = `<button class="close" id="inspClose">×</button>`;
  const an = annotOf(d);      // mRNA, DNA and chromatin identities (polymer_annot.py)
  const anRows = an ? an.rows.map(([k, v]: [string, string]) => `<div>${k}</div><div>${v}</div>`).join('') : '';
  const protType = an?.gene ? lib.types.findIndex(t => t.gene === an.gene && !t.id.startsWith('CPX')) : -1;
  const anFoot = an ? `${an.note ? `<p class="sub">${an.note}</p>` : ''}
      ${protType >= 0 ? `<div class="grid2"><button class="small" id="anProt" style="grid-column: 1 / 3">Highlight the ${an.gene} protein</button></div>` : ''}
      ${an.orf ? `<div style="margin-top:9px"><a href="https://www.yeastgenome.org/locus/${an.orf}" target="_blank">SGD: ${an.gene || an.orf}</a></div>` : ''}` : '';
  if (d.kind === 'inst') {
    const t = lib.types[d.info.type];
    const isPDB = t.src === 'PDB';
    const nb = t.nb, nnt = t.nnt;
    h += `<h2>${t.gene}</h2><div>${t.name}</div><div style="margin-top:6px"><span class="badge ${isPDB ? 'pdb' : 'af'}">${srcLabel(t)}</span>
      ${t.oligomer > 1 ? `<span class="badge" style="background:#666">${t.oligomer}-mer</span>` : ''}</div>
      <div class="kv">
      ${t.genes && t.genes !== t.gene ? `<div>Subunits</div><div style="word-break:break-word">${t.genes.split(/;\s*/).slice(0, 40).join(', ')}${t.genes.split(';').length > 40 ? ' …' : ''}</div>` : ''}
      ${t.uniprot ? `<div>UniProt</div><div>${t.uniprot}</div>` : ''}
      <div>Location</div><div>${REGION_NAMES[d.info.region]}</div>
      ${anRows}
      <div>Copies / cell</div><div>${t.copies.toLocaleString('en-US')}</div>
      <div>Mass</div><div>${(t.mass / 1000).toFixed(1)} kDa</div>
      <div>Residues drawn</div><div>${(nb - nnt).toLocaleString('en-US')} aa${nnt ? ` + ${nnt.toLocaleString('en-US')} nt` : ''}</div>
      <div>Radius</div><div>${t.br.toFixed(1)} nm (Rg ${t.rg.toFixed(1)} nm)</div>
      ${diffPred[d.info.type] ? `<div>Diffusion</div><div>${diffPred[d.info.type][0] >= 10 ? diffPred[d.info.type][0].toFixed(0) : diffPred[d.info.type][0].toFixed(1)} µm²/s predicted in the ${diffPred[d.info.type][1]} (calibrated on measurements; Model details 14.5)</div>` : ''}
      ${!isPDB && t.plddt_mean != null ? `<div>Mean pLDDT</div><div>${t.plddt_mean.toFixed(1)} (${Math.round((t.plddt_frac_kept ?? 1) * 100)}% of residues kept, pLDDT ≥ 50)</div>` : ''}
      <div>Position</div><div>${d.info.p.map((x: number) => x.toFixed(1)).join(', ')} nm</div>
      </div>
      ${copyHist(lib.types.map(x => x.copies), t.copies, 'molecular species')}
      <div class="grid2"><button class="small" id="loadAtoms">Load atomic model</button><button class="small" id="hlType">Highlight all</button>
      <button class="small" id="centerOn">Center view</button><button class="small" id="hideType">${typeFlags[d.info.type] & 2 ? 'Show' : 'Hide'} this type</button>
      ${lib.netNb[d.info.type] ? `<button class="small${selNet === d.info.type ? ' on' : ''}" id="netNb" style="grid-column: 1 / 3" title="Types whose STRING network neighbourhoods (personalized PageRank profiles) correlate most with this one">Network neighbours (${lib.netNb[d.info.type].length})</button>` : ''}</div>
      ${selNet === d.info.type ? `<div class="sub" style="margin-top:6px">${lib.netNb[d.info.type].map(([j, r]) => `${lib.types[j].gene.replace(/^CPX_\d+_/, '').replace(/_/g, ' ')} <span style="opacity:.7">${r.toFixed(2)}</span>`).join(' · ')}</div>` : ''}
      <div style="margin-top:9px">${t.uniprot ? `<a href="https://www.uniprot.org/uniprotkb/${t.uniprot}" target="_blank">UniProt</a>` : ''}
      ${isPDB ? `<a href="https://www.rcsb.org/structure/${t.pdb}" target="_blank">RCSB ${t.pdb}</a>` : ''}
      ${(t.model || t.uniprot) && !isPDB ? `<a href="https://alphafold.ebi.ac.uk/entry/${t.model || t.uniprot}" target="_blank">AlphaFold DB</a>` : ''}
      ${/^[A-Z0-9-]+$/.test(t.gene) && !SPECIAL[t.id] && !t.id.startsWith('CPX') ? `<a href="https://www.yeastgenome.org/locus/${t.gene}" target="_blank">SGD</a>` : ''}</div>${anFoot}`;
  } else {
    h += `<h2>${an?.title ?? d.title}</h2><div class="sub" style="font-size:12.5px">${an?.title ? d.title + ' · ' : ''}${d.sub || ''}</div>`;
    if (an) h += `<div class="kv">${anRows}</div>${an.hist ? copyHist(an.hist.vals, an.hist.v, an.hist.what) : ''}${anFoot}`;
    if (d.kind === 'lipid') {
      const comp = lipidComposition(d.mem, d.leaf);
      h += `<div class="kv">${comp.map((v, c) => [v, c] as const).filter(([v]) => v >= 0.05).sort((a, b) => b[0] - a[0])
        .map(([v, c]) => `<div>${LIPID_CLASSES[c]}</div><div>${v.toFixed(1)} mol%</div>`).join('')}</div>
        <p class="sub">Lipids are generated in your browser on the mesh of this membrane, at 0.65 nm² per lipid, and removed where they would overlap protein. Composition of this membrane${d.mem === 0 ? ' and leaflet' : ''} from isolated organelle membranes (Zinser et al. 1991)${d.mem === 0 ? '; PS and most PE face the cytoplasm, sphingolipids the outside; ergosterol 40 mol% (assumption)' : ''}. Ergosterol is drawn as two beads.</p>`;
    }
  }
  el.innerHTML = h;
  $('inspClose').onclick = () => select(null);
  if (protType >= 0) $('anProt').onclick = () => setType(protType);
  if (d.kind === 'inst') {
    $('loadAtoms').onclick = () => openMolstar(d.info.type);
    $('hlType').onclick = () => setType(d.info.type);
    if (lib.netNb[d.info.type]) $('netNb').onclick = () => setNetwork(selNet === d.info.type ? -1 : d.info.type);
    $('centerOn').onclick = () => { const st = cam.state(); st.target = d.info.p; cam.flyTo(st, 700); };
    $('hideType').onclick = () => { typeFlags[d.info.type] ^= 2; R.device.queue.writeBuffer(R.typeFlags, 0, typeFlags); renderInspect(); dirty = true; };
  }
}

// ---------------------------------------------------------------- Mol*
let molViewer: any = null, molLoading: Promise<void> | null = null;
function loadMolstarLib() {
  if (molLoading) return molLoading;
  const CDN = 'https://cdn.jsdelivr.net/npm/molstar@4.18.0/build/viewer/';
  const load = (base: string) => new Promise<void>((res, rej) => {
    const l = document.createElement('link'); l.rel = 'stylesheet'; l.href = base + 'molstar.css'; document.head.appendChild(l);
    const s = document.createElement('script'); s.src = base + 'molstar.js'; s.onload = () => res(); s.onerror = () => { l.remove(); s.remove(); rej(new Error('Mol* failed to load')); };
    document.head.appendChild(s);
  });
  // local copy first (offline use); CDN fallback for the single-file preview
  molLoading = (window.__YC_EMBED__ ? load(CDN) : load('molstar/').catch(() => load(CDN))).catch((e) => { molLoading = null; throw e; });
  return molLoading;
}
async function openMolstar(type: number) {
  const t = lib.types[type];
  $('molModal').style.display = 'flex';
  $('molTitle').textContent = `${t.gene} · ${t.name} · loading…`;
  try {
    await loadMolstarLib();
    const molstar = (window as any).molstar;
    if (!molViewer) {
      molViewer = await molstar.Viewer.create('molHost', {
        layoutIsExpanded: false, layoutShowControls: false, layoutShowRemoteState: false, layoutShowSequence: true,
        layoutShowLog: false, layoutShowLeftPanel: false, viewportShowExpand: false, viewportShowSelectionMode: false, viewportShowAnimation: false,
        collapseLeftPanel: true, pdbProvider: 'rcsb', emdbProvider: 'rcsb',
      });
    } else await molViewer.plugin.clear();
    let url = '', link = '';
    if (t.src === 'pred') {                // predicted in this project (tRNAs: Boltz-2), shipped with the data
      url = link = DATA_BASE + `models/${t.id}.cif`;
    } else if (t.src === 'PDB' && t.pdb) {
      url = `https://files.rcsb.org/download/${t.pdb}-assembly1.cif`; link = `https://www.rcsb.org/structure/${t.pdb}`;
      try { const h = await fetch(url, { method: 'HEAD' }); if (!h.ok) url = `https://files.rcsb.org/download/${t.pdb}.cif`; } catch { url = `https://files.rcsb.org/download/${t.pdb}.cif`; }
    } else {
      const entry = t.model || t.uniprot;   // predicted complexes are AFDB entries of their own
      link = `https://alphafold.ebi.ac.uk/entry/${entry}`;
      try {
        const r = await fetch(`https://alphafold.ebi.ac.uk/api/prediction/${entry}`);
        const j = await r.json(); url = (Array.isArray(j) ? j[0] : j)?.cifUrl;
      } catch { /* fall through */ }
      if (!url) url = t.model ? `https://alphafold.ebi.ac.uk/files/${t.model}-model_v1.cif` : `https://alphafold.ebi.ac.uk/files/AF-${t.uniprot}-F1-model_v6.cif`;
    }
    await molViewer.loadStructureFromUrl(url, 'mmcif', false);
    $('molTitle').textContent = `${t.gene} · ${t.name} · ${t.src === 'pred' ? 'predicted with Boltz-2 for this model (no modified nucleosides)' : t.src === 'PDB' ? 'PDB ' + t.pdb + ' (biological assembly)' : t.model ? `AlphaFold DB predicted ${t.oligomer}-mer (full length, incl. low-pLDDT regions)` : 'AlphaFold DB model (full length, incl. low-pLDDT regions)'}`;
    const a = $<HTMLAnchorElement>('molLink'); a.href = link; a.textContent = t.src === 'pred' ? 'Download the model (mmCIF)' : 'Open source page';
  } catch (e: any) {
    $('molTitle').textContent = `Could not load the atomic model: ${e?.message || e}. Internet access to RCSB/AlphaFold DB is needed.`;
  }
}

// ---------------------------------------------------------------- display modes
let orgMeshesRequested = false;
function loadOrganelleMeshes() {
  if (orgMeshesRequested) return;
  orgMeshesRequested = true;
  toast('Loading organelle surfaces (6.7 MB)…');
  fetchBin('organelle_meshes.bin.gz').then(b => { R.uploadOrganelleMeshes(b); dirty = true; });
}
// membranes view = the Show list with only Membranes ticked (fs.memColors); set up with the Show list
let setMembraneView: (on: boolean) => void = () => {};
// guided views may restrict the Show list to some entries (Bookmark.show); set up with the Show list
let setShowOnly: (names: string[] | null) => void = () => {};
let showRestricted = false;
function syncModeButtons() {
  const mem = fs.memColors && !fs.organelle;
  $('modeSection').classList.toggle('on', !fs.organelle && !mem); $('modeOrg').classList.toggle('on', fs.organelle);
  $('modeMem').classList.toggle('on', mem);
}
function setOrganelle(on: boolean) {
  fs.organelle = on;
  if (on) loadOrganelleMeshes();
  syncModeButtons();
  dirty = true;
}
function setDisplayMode(m: 'molecules' | 'organelles' | 'membranes') {
  if (m !== 'membranes' && fs.memColors) setMembraneView(false);
  setOrganelle(m === 'organelles');
  if (m === 'membranes') setMembraneView(true);
}
function setQuality(q: string) {
  quality = q;
  const P: Record<string, [number, number]> = { low: [2.6, 0.75], med: [1.6, 1], high: [1.0, 1], ultra: [0.6, 1], auto: [1.6, 1] };
  const [l0, sc] = P[q];
  fs.l0px = l0; fs.renderScale = sc; R.resize(sc); dirty = true;
}
function renderLegend() {
  const el = $('legend');
  const sw = (c: string, l: string) => `<div><span class="sw" style="background:${c}"></span>${l}</div>`;
  const m = fs.colorMode;
  if (m === 0) el.innerHTML = [sw(PAL[0], 'cytosolic protein'), sw(PAL[1], 'nuclear protein'), sw('#EE9BC0', 'ribosome protein'), sw('#EF9A7E', 'pre-ribosome'), sw('#C9558F', 'rRNA'),
    sw('#E27FA8', 'tRNA'), sw('#B0437A', 'mRNA'), sw('#F4B26A', 'histone'), sw('#E07636', 'DNA'), sw(PAL[23], 'nuclear pore'), sw('#7CC6C2', 'FG repeats'), sw('#E0893A', 'actin filament'), sw('#D9534F', 'Arp2/3 complex'), sw('#8E6CB0', 'capping protein'), sw('#2FA4A0', 'tropomyosin'), sw('#74C476', 'microtubule'), sw('#2E8B57', 'γ-tubulin ring'), sw(PAL[14], 'membrane protein'),
    sw(PAL[5], 'mitochondrial'), sw(PAL[3], 'ER lumen'), sw('#F3C643', 'lipid head'), sw('#F8E3A0', 'lipid tail'), sw('#CDA866', 'β-glucan'), sw('#C3C27A', 'mannan'),
    sw(PAL[9], 'wall protein'), sw('#8B88C9', 'large complex')].join('');
  else if (m === 1) el.innerHTML = sw('#4C8CD9', 'experimental (PDB)') + sw('#F5A84C', 'AlphaFold model');
  else if (m >= 10 && m <= 16) {
    const mk = ['H2A.Z (Htz1)', 'H3K4me3', 'H3K36me3', 'H3K79me3', 'H3K56ac', 'H3K9ac', 'H4K16ac'][m - 10];
    el.innerHTML = sw('#D9335A', `nucleosome enriched for ${mk}`) + sw('#CCC7B8', 'other nucleosome') + sw('#E6E6E3', 'everything else') +
      '<div class="sub" style="grid-column: 1 / 3"><span>Per measured nucleosome: the mark of the nearest unstressed MNase nucleosome within 40 bp (Weiner et al. 2015), enriched at ≥ 0.5 log2. Annotation only: the marks do not change the structure.</span></div>';
  }
  else if (m === 2) el.innerHTML = sw('#0053D6', 'pLDDT > 90') + sw('#65CBF3', '70–90') + sw('#FFDB13', '50–70') + sw('#B8B8B8', 'experimental (PDB)');
  else if (m === 3) el.innerHTML = `<div style="grid-column: 1 / 3"><div style="height:8px;border-radius:3px;background:linear-gradient(90deg,#440154,#3b528b,#21918c,#5ec962,#fde725)"></div><div style="display:flex;justify-content:space-between"><span>1</span><span>10³</span><span>10⁶ copies</span></div></div>`;
  else if (m === 6) {
    el.innerHTML = lib.netOn
      ? '<div class="sub" style="grid-column: 1 / 3"><span>Similar colours: similar STRING network neighbourhoods (personalized PageRank profiles). Click a molecule, then <b>Network neighbours</b> to highlight its closest ones.</span></div>' +
        (selNet >= 0 ? `<div style="grid-column: 1 / 3"><span>Highlighted: ${lib.types[selNet].gene} and ${(lib.netNb[selNet] || []).length} neighbours <a href="#" id="netClear" style="color:#8fc3ff">clear</a></span></div>` : '') +
        sw('#D1D1CC', 'no network data')
      : '<div style="grid-column: 1 / 3">No network profiles in this build</div>';
    const nc = el.querySelector('#netClear'); if (nc) (nc as HTMLElement).onclick = (e) => { e.preventDefault(); setNetwork(-1); };
  }
  else if (m === 5) el.innerHTML = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 15, 18, 19].map(i => sw(PAL[i], REGION_NAMES[i])).join('');
  else el.innerHTML = '';
}

// ---------------------------------------------------------------- UI construction
function syncPlaneUI() {
  const d = fs.plane[3];
  ($('planeD') as HTMLInputElement).value = String(Math.round(d));
  $('planeDv').textContent = `${Number(d.toFixed(nudgeStep < 0.5 ? 2 : 1))} nm`;
  const s = fs.slab;
  ($('slab') as HTMLInputElement).value = String(Math.round(Math.log(s / 5) / Math.log(4400 / 5) * 1000));
  $('slabv').textContent = fmtNm(s);
  ['axX', 'axY', 'axZ'].forEach((id, k) => $(id).classList.toggle('on', Math.abs(fs.plane[k]) > 0.999));
  mmDirty = true;
}
function buildUI() {
  const NT = lib.types.length;
  let total = 0; for (const t of lib.types) total += t.copies;
  $('subCount').textContent = `${(total / 1e6).toFixed(1)} M macromolecules of ${NT.toLocaleString('en-US')} species, each a real PDB or AlphaFold structure`;
  if ((window as any).__YC_PREVIEW__) {
    const nI = idx.nInst.reduce((a, b) => a + b, 0);
    $('subCount').innerHTML = `<b>Preview:</b> ${idx.n} of 6,416 tiles (${(nI / 1e3).toFixed(0)} k molecules) from the full ${(total / 1e6).toFixed(1)} M-molecule cell. The full cell is in yeast_cell_viewer.zip.`;
    ($('modeOrg') as HTMLButtonElement).disabled = true; $('modeOrg').title = 'Organelle surfaces are only in the full app';
    ($('modeMem') as HTMLButtonElement).disabled = true; $('modeMem').title = 'Membrane surfaces are only in the full app';
  }
  // views
  const v = $('views');
  bookmarks.forEach((b, i) => {
    const btn = document.createElement('button'); btn.textContent = `${(i + 1) % 10 === 0 && i < 10 ? 0 : i < 9 ? i + 1 : ''}${i < 10 ? ' · ' : ''}${b.name}`;
    btn.title = b.caption; btn.onclick = () => { stopTour(); applyBookmark(i); closeSheet(); };
    v.appendChild(btn);
  });
  $('tourBtn').onclick = () => { if (tour) stopTour(); else { startTour(); closeSheet(); } };
  // search
  const inp = $<HTMLInputElement>('search'), res = $('results');
  let act = 0, cur: number[] = [], loci: string[] = [];
  const pick = (k: number) => {          // k < cur.length: a molecule type; otherwise a gene locus
    if (k < cur.length) { setType(cur[k]); flyToNearest(cur[k]); } else flyToLocus(loci[k - cur.length]);
    cur = []; loci = []; draw(); inp.value = ''; inp.blur();
  };
  const draw = () => {
    res.innerHTML = cur.map((i, k) => { const t = lib.types[i]; return `<div data-k="${k}" class="${k === act ? 'act' : ''}"><span class="nm"><b>${t.gene}</b> ${t.name}</span><span class="n">${fmt(t.copies)}</span></div>`; }).join('')
      + loci.map((g, j) => { const v = geneLoci[g], k = cur.length + j; return `<div data-k="${k}" class="${k === act ? 'act' : ''}"><span class="nm"><b>${g}</b> gene locus · chr ${ROMAN[v[3]]}:${Math.round(v[4] / 1000)} kb</span><span class="n">DNA</span></div>`; }).join('');
    res.querySelectorAll('div').forEach(d => d.addEventListener('click', () => pick(Number((d as HTMLElement).dataset.k))));
  };
  inp.addEventListener('input', () => { cur = searchTypes(inp.value).slice(0, 10); loci = searchLoci(inp.value); act = 0; draw(); });
  inp.addEventListener('keydown', (e) => {
    if (e.key === 'ArrowDown') { act = Math.min(cur.length + loci.length - 1, act + 1); draw(); e.preventDefault(); }
    else if (e.key === 'ArrowUp') { act = Math.max(0, act - 1); draw(); e.preventDefault(); }
    else if (e.key === 'Enter' && cur.length + loci.length) { pick(act); }
    else if (e.key === 'Escape') { cur = []; loci = []; draw(); inp.blur(); }
    e.stopPropagation();
  });
  // cut
  $<HTMLInputElement>('cutOn').onchange = (e) => { fs.cutOn = (e.target as HTMLInputElement).checked; dirty = true; };
  const ax = (k: number) => () => { const n: V3 = [0, 0, 0]; n[k] = dot(cam.fwd, [k === 0 ? 1 : 0, k === 1 ? 1 : 0, k === 2 ? 1 : 0]) >= 0 ? 1 : -1; fs.cutOn = true; ($('cutOn') as HTMLInputElement).checked = true; setPlane(n, dot(n, cam.target)); };
  $('axX').onclick = ax(0); $('axY').onclick = ax(1); $('axZ').onclick = ax(2);
  $('alignView').onclick = () => { fs.cutOn = true; ($('cutOn') as HTMLInputElement).checked = true; alignToView(); };
  $('followView').onclick = () => { planeFollow = !planeFollow; $('followView').classList.toggle('on', planeFollow); if (planeFollow) alignToView(); };
  $<HTMLInputElement>('planeD').oninput = (e) => { const d = Number((e.target as HTMLInputElement).value); const n = planeN(); const dd = d - fs.plane[3]; cam.target = add(cam.target, scale(n, dd)); setPlane(n, d); };
  const stepSel = $<HTMLSelectElement>('nudgeStep');
  if ([...stepSel.options].some(o => Number(o.value) === nudgeStep)) stepSel.value = String(nudgeStep);
  stepSel.onchange = () => {
    nudgeStep = Number(stepSel.value); syncPlaneUI(); stepSel.blur();
    try { localStorage.setItem('yc.nudgeStep', String(nudgeStep)); } catch { /* storage unavailable */ }
  };
  // -/+ buttons: one step per click, repeating while held
  const holdRepeat = (id: string, f: () => void) => {
    const b = $(id); let t0 = 0, t1 = 0;
    const stop = () => { clearTimeout(t0); clearInterval(t1); };
    b.addEventListener('pointerdown', (e) => {
      e.preventDefault(); stop(); sweep = null; $('sweepBtn').classList.remove('on');
      f();
      t0 = window.setTimeout(() => { t1 = window.setInterval(f, 60); }, 400);
    });
    for (const ev of ['pointerup', 'pointerleave', 'pointercancel']) b.addEventListener(ev, stop);
  };
  holdRepeat('nudgeBack', () => movePlane(-nudgeStep));
  holdRepeat('nudgeFwd', () => movePlane(nudgeStep));
  // touch-screen pill: coarse steps (a quarter of the slab), like the , and . keys
  const pillStep = () => Math.max(2, fs.slab * 0.25);
  const cutOnIfOff = () => { if (!fs.cutOn) { fs.cutOn = true; ($('cutOn') as HTMLInputElement).checked = true; } };
  holdRepeat('pillBack', () => { cutOnIfOff(); movePlane(-pillStep()); });
  holdRepeat('pillFwd', () => { cutOnIfOff(); movePlane(pillStep()); });
  $('pillCut').onclick = () => { fs.cutOn = !fs.cutOn; ($('cutOn') as HTMLInputElement).checked = fs.cutOn; dirty = true; };
  $<HTMLInputElement>('slab').oninput =(e) => { const x = Number((e.target as HTMLInputElement).value) / 1000; fs.slab = 5 * Math.pow(4400 / 5, x); syncPlaneUI(); dirty = true; };
  $('sweepBtn').onclick = () => {
    if (sweep) { sweep = null; $('sweepBtn').classList.remove('on'); return; }
    const n = planeN(); fs.cutOn = true; ($('cutOn') as HTMLInputElement).checked = true;
    const R0 = geo.G[0];
    sweep = { t0: performance.now(), d0: -R0, d1: R0, dur: 30000 };
    cam.target = add(cam.target, scale(n, -R0 - fs.plane[3]));
    $('sweepBtn').classList.add('on');
  };
  $('flipBtn').onclick = () => { const n = planeN(); const st = cam.state(); const tgt = cam.target; const nn = scale(n, -1) as V3; setPlane(nn, -fs.plane[3]);
    st.fwd = scale(st.fwd, -1) as V3; st.target = tgt; cam.flyTo(st, 1200); };
  // display
  $('modeSection').onclick = () => setDisplayMode('molecules');
  $('modeOrg').onclick = () => setDisplayMode('organelles');
  $('modeMem').onclick = () => setDisplayMode('membranes');
  $<HTMLSelectElement>('colorMode').onchange = (e) => {
    fs.colorMode = Number((e.target as HTMLSelectElement).value);
    renderLegend(); dirty = true;
  };
  $<HTMLInputElement>('shaded').onchange = (e) => { fs.shaded = (e.target as HTMLInputElement).checked; dirty = true; };
  $<HTMLSelectElement>('motion').onchange = (e) => { fs.motion = Number((e.target as HTMLSelectElement).value); (e.target as HTMLElement).blur(); dirty = true; };
  $<HTMLInputElement>('outlines').onchange = (e) => { fs.outline = (e.target as HTMLInputElement).checked ? 1 : 0; dirty = true; };
  renderLegend();
  // show list: named classes have their own switch (and ignore the compartment filter); "only" hides everything else
  const show = $('showList');
  const cls = (k: string) => (on: boolean) => setTypeHidden(CLASS_TYPES[k], !on);
  const classes: [string, string, (on: boolean) => void][] = [
    ['Actin', '#E0893A', cls('actin')],
    ['Microtubules', '#74C476', cls('microtubules')],
    ['Eisosomes', '#3F8FCB', cls('eisosomes')],
    ['Ribosomes', '#EE9BC0', cls('ribosomes')],
    ['tRNA', '#E27FA8', cls('tRNA')],
    ['Nucleosomes', '#F4B26A', cls('nucleosomes')],
    ['Nuclear pores', PAL[23], cls('pores')],
    ['RNA polymerases', '#9A72B8', cls('polymerases')],
    ['mRNA', '#B0437A', on => setFree(0, on)],
    ['DNA', '#E07636', on => setFree(1, on)],
    ['pre-rRNA', '#C9558F', on => setFree(4, on)],
    ['FG repeats', '#7CC6C2', on => setFree(5, on)],
    ['β-Glucan', '#CDA866', on => setFree(2, on)],
    ['Mannan', '#C3C27A', on => setFree(3, on)],
    ['Chitin', PAL[54], on => setFree(6, on)],
    ['Polyphosphate', PAL[55], on => setFree(7, on)],
    ['Lipids', '#F3C643', on => { fs.showLipids = on; dirty = true; }],
    ['Membranes', '#F8E3A0', on => { fs.showMeshes = on; dirty = true; }],
    ['Soluble proteins', PAL[0], on => setRegions([0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 24, 25, 26], on)],
    ['Membrane proteins', PAL[14], on => setRegions([14, 15, 16, 17, 18, 19, 20, 21], on)],
  ];
  const boxes: HTMLInputElement[] = [];
  const allRegions = REGION_NAMES.map((_, i) => i);
  // only membranes shown: colour them by the compartment they enclose (the section map's colours)
  const mem = classes.findIndex(([name]) => name === 'Membranes');
  // (the whole cell's membrane surfaces, as in organelle mode, without the cut: see renderer.frame)
  const syncMemColors = () => {
    fs.memColors = boxes.every((b, k) => b.checked === (k === mem)) && !(window as any).__YC_PREVIEW__;   // preview: no surfaces
    if (fs.memColors) loadOrganelleMeshes();
    syncModeButtons();
    dirty = true;
  };
  const setAll = (on: (k: number) => boolean) => {
    classes.forEach(([, , f], k) => { boxes[k].checked = on(k); f(on(k)); });
    const regionsOn = classes.every((_, k) => on(k));       // "only": every compartment off (e.g. cell-wall proteins too)
    setRegions(allRegions, regionsOn);
    syncMemColors();
  };
  classes.forEach(([name, col, f], k) => {
    const l = document.createElement('label'); l.className = 'chk';
    l.innerHTML = `<input type="checkbox" checked><span class="sw" style="background:${col}"></span><span class="nm">${name}</span><a class="only" title="Show only ${name.toLowerCase()}">only</a>`;
    const box = l.querySelector('input') as HTMLInputElement; boxes.push(box);
    box.onchange = () => { f(box.checked); syncMemColors(); };
    (l.querySelector('.only') as HTMLElement).onclick = (e) => { e.preventDefault(); e.stopPropagation(); setAll(j => j === k); };
    show.appendChild(l);
  });
  setMembraneView = (on) => setAll(on ? (j => j === mem) : (() => true));
  setShowOnly = (names) => setAll(names ? (j => names.includes(classes[j][0])) : (() => true));
  $('showAll').onclick = () => setAll(() => true);
  $('hideAll').onclick = () => setAll(() => false);
  const rl = $('regionList');
  REGION_NAMES.forEach((n, i) => {
    const l = document.createElement('label'); l.className = 'chk';
    l.innerHTML = `<input type="checkbox" checked data-r="${i}"><span class="sw" style="background:${PAL[i]}"></span>${n}`;
    (l.querySelector('input') as HTMLInputElement).onchange = (e) => setRegions([i], (e.target as HTMLInputElement).checked);
    rl.appendChild(l);
  });
  $<HTMLSelectElement>('quality').onchange = (e) => setQuality((e.target as HTMLSelectElement).value);
  $<HTMLInputElement>('detailOn').onchange = (e) => {
    T.detailEnabled = (e.target as HTMLInputElement).checked;
    if (!T.detailEnabled) for (const t of T.loaded.values()) if (t.det) T.dropDetail(t);
    dirty = true;
  };
  $('shotBtn').onclick = screenshot;
  $('helpBtn').onclick = () => ($('helpModal').style.display = 'flex');
  $('aboutBtn').onclick = () => ($('aboutModal').style.display = 'flex');
  $('aboutBody').innerHTML = ABOUT_HTML;
  $('modelBtn').onclick = () => {
    if (!$('modelBody').innerHTML) $('modelBody').innerHTML = markdownToHtml(MODEL_MD);
    $('modelModal').style.display = 'flex';
  };
  $('aboutModal').addEventListener('click', (e) => {   // "Model details" link inside the About text
    if ((e.target as HTMLElement).id === 'aboutToModel') { e.preventDefault(); $('aboutModal').style.display = 'none'; $('modelBtn').click(); }
  });
  document.querySelectorAll('.modal').forEach(m => m.addEventListener('click', (e) => {
    if (e.target === m || (e.target as HTMLElement).classList.contains('close')) (m as HTMLElement).style.display = 'none';
  }));
  $('molClose').onclick = () => ($('molModal').style.display = 'none');
  $('panelToggle').onclick = () => setPanel($('panel').classList.contains('collapsed'));
  $('sheetClose').onclick = () => setPanel(false);
  phoneMQ.addEventListener('change', () => setPanel(!phoneMQ.matches));
  setPanel(!phoneMQ.matches && !matchMedia('(max-height: 500px)').matches);   // phones start with the view, not the menu
  $('caption').onclick = () => showCaption(null);
  syncPlaneUI();
  setupInput();
  setupMinimap();
}
// phones in portrait: the panel is a bottom sheet, closed by default and after picking a view or a molecule
const phoneMQ = matchMedia('(max-width: 760px)');
function setPanel(open: boolean) {
  $('panel').classList.toggle('collapsed', !open);
  $('panelToggle').textContent = phoneMQ.matches ? '☰' : open ? '‹' : '›';
}
function closeSheet() { if (phoneMQ.matches) setPanel(false); }
function setTypeHidden(ids: string[], hidden: boolean) {
  for (const id of ids) { const i = lib.byId.get(id); if (i === undefined) continue; if (hidden) typeFlags[i] |= 2; else typeFlags[i] &= ~2; }
  R.device.queue.writeBuffer(R.typeFlags, 0, typeFlags); dirty = true;
}
function setFree(kind: number, on: boolean) { if (on) fs.freeHide &= ~(1 << kind); else fs.freeHide |= 1 << kind; dirty = true; }
function setRegions(rs: number[], on: boolean) {
  for (const r of rs) { if (on) fs.regionMask |= 1 << r; else fs.regionMask &= ~(1 << r); }
  fs.regionMask >>>= 0;
  document.querySelectorAll<HTMLInputElement>('#regionList input').forEach(inp => { inp.checked = ((fs.regionMask >>> Number(inp.dataset.r)) & 1) === 1; });
  dirty = true;
}

// ---------------------------------------------------------------- input
function setupInput() {
  let drag: null | { x: number; y: number; btn: number; moved: number; shift: boolean } = null;
  const keys = new Set<string>();
  // touch: one finger orbits (like a left drag), two fingers pinch-zoom, pan and twist-roll
  const touches = new Map<number, { x: number; y: number }>();
  let pinch: null | { d: number; mx: number; my: number; a: number } = null;
  let multi = false, lastTap = { t: 0, x: 0, y: 0 };
  const pinchState = () => {
    const [p, q] = [...touches.values()];
    const r = canvas.getBoundingClientRect();
    return { d: Math.hypot(q.x - p.x, q.y - p.y), mx: (p.x + q.x) / 2 - r.left, my: (p.y + q.y) / 2 - r.top, a: Math.atan2(q.y - p.y, q.x - p.x) };
  };
  canvas.addEventListener('contextmenu', e => e.preventDefault());
  canvas.addEventListener('pointerdown', (e) => {
    canvas.setPointerCapture(e.pointerId); canvas.focus();
    if (tour) stopTour();
    if (e.pointerType === 'touch') {
      touches.set(e.pointerId, { x: e.clientX, y: e.clientY });
      if (touches.size === 2) { drag = null; multi = true; pinch = pinchState(); return; }
      if (touches.size > 2) return;
      multi = false;
    }
    drag = { x: e.clientX, y: e.clientY, btn: e.button, moved: 0, shift: e.shiftKey };
    canvas.classList.add('dragging');
  });
  canvas.addEventListener('pointermove', (e) => {
    if (e.pointerType === 'touch') {
      if (!touches.has(e.pointerId)) return;
      touches.set(e.pointerId, { x: e.clientX, y: e.clientY });
      if (touches.size >= 2 && pinch) {
        const p = pinchState();
        if (p.d > 0 && pinch.d > 0) zoomAt(p.mx, p.my, pinch.d / p.d);
        cam.pan(p.mx - pinch.mx, p.my - pinch.my, cam.pxPerNmAtUnit(canvas.clientHeight) / cam.dist);
        let da = p.a - pinch.a;
        if (da > Math.PI) da -= 2 * Math.PI; else if (da < -Math.PI) da += 2 * Math.PI;
        cam.roll(-da);
        pinch = p;
        if (planeFollow) alignToView();
        cam.anim = null; dirty = true;
        return;
      }
    }
    if (!drag) { if (e.pointerType === 'mouse' || e.pointerType === 'pen') hover(e.clientX, e.clientY); return; }
    const dx = e.clientX - drag.x, dy = e.clientY - drag.y;
    drag.x = e.clientX; drag.y = e.clientY; drag.moved += Math.abs(dx) + Math.abs(dy);
    if (drag.btn === 2 || drag.btn === 1 || drag.shift) {
      const pxPerNm = cam.pxPerNmAtUnit(canvas.clientHeight) / cam.dist;
      cam.pan(dx, dy, pxPerNm);
    } else {
      cam.orbit(dx * 0.005, dy * 0.005);
      if (planeFollow) alignToView();
    }
    cam.anim = null; dirty = true;
  });
  const end = (e: PointerEvent) => {
    if (e.pointerType === 'touch') {
      if (!touches.delete(e.pointerId)) return;
      if (touches.size === 1) {   // second finger lifted: keep orbiting with the remaining one, no tap
        const [p] = [...touches.values()];
        pinch = null; drag = { x: p.x, y: p.y, btn: 0, moved: 99, shift: false };
        return;
      }
      if (touches.size > 1) { pinch = pinchState(); return; }
      pinch = null;
    }
    if (!drag) return;
    const touch = e.pointerType === 'touch';
    const wasClick = drag.moved < (touch ? 10 : 4) && drag.btn === 0 && !(touch && multi);
    drag = null; canvas.classList.remove('dragging');
    if (e.type === 'pointercancel' || !wasClick) return;
    if (touch) {   // double-tap flies there (no dblclick for touch)
      const now = performance.now();
      if (now - lastTap.t < 320 && Math.hypot(e.clientX - lastTap.x, e.clientY - lastTap.y) < 30) {
        lastTap.t = 0; flyToPixel(e.offsetX, e.offsetY); return;
      }
      lastTap = { t: now, x: e.clientX, y: e.clientY };
    }
    if (simPlacing && !sim.active) {   // simulation box placement: the box follows the view centre; a click moves it there
      worldAtPixel(e.offsetX, e.offsetY).then(p => { if (p) { const st = cam.state(); st.target = p; cam.flyTo(st, 500); } });
      return;
    }
    if (fs.simMode && sim?.active) {   // simulation: a click on a molecule toggles its tracking
      pickAt(e.offsetX, e.offsetY).then(r => { if (r.id >>> 30 === 1) { sim.toggleTrack(r.id & 0x1FFFFFFF); dirty = true; } });   // bit 29: faded
      return;
    }
    pickAt(e.offsetX, e.offsetY).then(r => { const d = describe(r.id); select(d.kind === 'inst' || d.kind === 'lipid' || d.kind === 'free' || d.kind === 'band' ? d : null); });
  };
  canvas.addEventListener('pointerup', end); canvas.addEventListener('pointercancel', end);
  const flyToPixel = (x: number, y: number) => {
    worldAtPixel(x, y).then(p => {
      if (!p) return;
      const st = cam.state(); st.target = p; st.dist = Math.max(cam.minDist * 4, cam.dist * 0.55);
      cam.flyTo(st, 900);
    });
  };
  canvas.addEventListener('dblclick', (e) => { if (!touches.size) flyToPixel(e.offsetX, e.offsetY); });
  // zoom by factor f (< 1 zooms in) towards a canvas pixel
  const zoomAt = (x: number, y: number, f: number) => {
    const aspect = canvas.clientWidth / canvas.clientHeight;
    const ray = cam.ray(x / canvas.clientWidth, y / canvas.clientHeight, aspect);
    const pc = add(cam.eye, scale(ray, cam.dist / dot(ray, cam.fwd)));
    const old = cam.dist; cam.zoom(f);
    const k = 1 - cam.dist / old;
    cam.target = add(cam.target, scale(sub(pc, cam.target), k));
  };
  canvas.addEventListener('wheel', (e) => {
    e.preventDefault();
    // macOS turns Shift+wheel on a mouse into horizontal scrolling: use whichever axis moved
    const raw = e.deltaY || e.deltaX;
    const dy = e.deltaMode === 1 ? raw * 30 : raw;
    if (e.altKey) {
      // fine cut-plane step: one fine step per mouse-wheel notch (~100 px), proportional for trackpads, at most 4 per event
      movePlane(Math.sign(dy) * nudgeStep * Math.min(4, Math.abs(dy) / 100));
    } else if (e.shiftKey) {
      const n = planeN(), step = Math.max(1, Math.min(fs.slab, cam.dist * 0.05)) * Math.sign(dy) * Math.min(4, Math.abs(dy) / 60 + 0.5);
      cam.target = add(cam.target, scale(n, step)); setPlane(n, fs.plane[3] + step);
    } else {
      // zoom towards the cursor (Cmd: 8x slower; Ctrl is left alone because trackpad pinch arrives as Ctrl+wheel)
      zoomAt(e.offsetX, e.offsetY, Math.exp(dy * 0.0015 * (e.metaKey ? 0.125 : 1)));
    }
    cam.anim = null; if (tour) stopTour(); dirty = true;
  }, { passive: false });
  window.addEventListener('keydown', (e) => {
    if ((e.target as HTMLElement).tagName === 'INPUT' || (e.target as HTMLElement).tagName === 'SELECT') return;
    const k = e.key.toLowerCase();
    keys.add(k);
    if (k === 'c') { fs.cutOn = !fs.cutOn; ($('cutOn') as HTMLInputElement).checked = fs.cutOn; dirty = true; }
    else if (k === 'v') { fs.cutOn = true; ($('cutOn') as HTMLInputElement).checked = true; alignToView(); }
    else if (k === '[') { fs.slab = Math.max(5, fs.slab / 1.25); syncPlaneUI(); dirty = true; }
    else if (k === ']') { fs.slab = Math.min(4400, fs.slab * 1.25); syncPlaneUI(); dirty = true; }
    else if (k === ',' || k === '.' || e.code === 'Comma' || e.code === 'Period') {   // Alt: one fine step (Option changes e.key on macOS)
      const fwd = k === '.' || e.code === 'Period';
      movePlane((fwd ? 1 : -1) * (e.altKey ? nudgeStep : Math.max(2, fs.slab * 0.25))); }
    else if (k === 'o') setDisplayMode(fs.organelle ? 'molecules' : 'organelles');
    else if (k === 'm') setDisplayMode(fs.memColors && !fs.organelle ? 'molecules' : 'membranes');
    else if (k === 't') (tour ? stopTour() : startTour());
    else if (k === 'p') screenshot();
    else if (k === 'escape') { select(null); setType(-1); }
    else if (/^[0-9]$/.test(k)) { const i = k === '0' ? 9 : Number(k) - 1; if (i < bookmarks.length) { stopTour(); applyBookmark(i); } }
  });
  window.addEventListener('keyup', (e) => keys.delete(e.key.toLowerCase()));
  window.addEventListener('blur', () => keys.clear());
  (setupInput as any).keys = keys;
  window.addEventListener('resize', () => { R.resize(fs.renderScale); dirty = true; });
}
let hoverT = 0, hoverBusy = false, hoverAt: { d: any; x: number; y: number } | null = null;
function hover(x: number, y: number) {
  const now = performance.now();
  if (hoverBusy || now - hoverT < 90) return;
  hoverT = now; hoverBusy = true;
  dirty = true;
  pickAt(x, y).then(r => { hoverBusy = false; showHover(describe(r.id), x, y); });
}
function showHover(d: any, x: number, y: number) {
  const el = $('hover');
  if (!d.title || d.kind === 'none') { el.style.display = 'none'; hoverAt = null; return; }
  hoverAt = { d, x, y };
  const an = annotOf(d);      // gene, transcript or genome position of mRNA, DNA and chromatin (loaded on demand)
  let txt = `<b>${an?.title ?? d.title}</b>`;
  if (d.kind === 'inst') { const t = lib.types[d.info.type]; txt += ` · ${t.name} <span style="color:#aaa">· ${d.sub}</span>`; }
  else if (!an) txt += ` <span style="color:#aaa">· ${d.sub}</span>`;
  if (an) txt += `<br><span style="color:#ccc">${an.line}</span>`;
  el.innerHTML = txt; el.style.display = 'block';
  el.style.left = `${x + 14}px`; el.style.top = `${y + 12}px`;
}
onAnnot = () => {
  if (hoverAt && $('hover').style.display === 'block') showHover(hoverAt.d, hoverAt.x, hoverAt.y);
  if (selection && (selection.kind === 'free' || selection.kind === 'inst')) renderInspect();
};
canvas.addEventListener('pointerleave', () => { $('hover').style.display = 'none'; hoverAt = null; });

// ---------------------------------------------------------------- minimap (analytic section through the cut plane)
let mmDirty = true, mmLastT = 0;
let mmBasis: { u: V3; v: V3; c: V3; ext: number } | null = null;
const mmImg: { data: ImageData | null } = { data: null };
function setupMinimap() {
  const c = $<HTMLCanvasElement>('minimap');
  c.addEventListener('click', (e) => {
    if (!mmBasis) return;
    const r = c.getBoundingClientRect();
    const x = ((e.clientX - r.left) / r.width * 2 - 1) * mmBasis.ext, y = (1 - (e.clientY - r.top) / r.height * 2) * mmBasis.ext;
    const p = add(mmBasis.c, add(scale(mmBasis.u, x), scale(mmBasis.v, y)));
    const st = cam.state(); st.target = p; cam.flyTo(st, 1200);
    if (tour) stopTour();
  });
}
function updateMinimap(now: number) {
  const c = $<HTMLCanvasElement>('minimap'), g = c.getContext('2d')!;
  const W = c.width, H = c.height;
  const n: V3 = fs.cutOn ? planeN() : cam.fwd;
  const d = fs.cutOn ? fs.plane[3] : dot(cam.fwd, cam.target);
  if (mmDirty && now - mmLastT > 250) {
    mmDirty = false; mmLastT = now;
    let u = sub(cam.right, scale(n, dot(cam.right, n)));
    if (len(u) < 1e-3) u = norm(cross(n, [0, 0, 1])); u = norm(u);
    let v = cross(n, u); if (dot(v, cam.up) < 0) v = scale(v, -1) as V3;
    const ext = 2300, c0 = scale(n, d) as V3;
    mmBasis = { u, v, c: c0, ext };
    const N = 120, img = g.createImageData(N, N);
    const cols = [PAL[80], PAL[9], '#D2BE8E', PAL[83], PAL[3], PAL[3], PAL[1], PAL[2], PAL[4], PAL[6], PAL[5], '#E9C46A', PAL[7], PAL[8], '#C9C2DC']
      .map(h => { const x = parseInt(h.slice(1), 16); return [(x >> 16) & 255, (x >> 8) & 255, x & 255]; });
    cols[3] = [196, 214, 235];
    for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) {
      const x = ((i + 0.5) / N * 2 - 1) * ext, y = (1 - (j + 0.5) / N * 2) * ext;
      const p = add(c0, add(scale(u, x), scale(v, y)));
      const lab = geom.label(p[0], p[1], p[2]);
      const col = cols[lab], o = 4 * (j * N + i);
      img.data[o] = col[0]; img.data[o + 1] = col[1]; img.data[o + 2] = col[2]; img.data[o + 3] = lab === 0 ? 0 : 255;
    }
    mmImg.data = img;
  }
  if (!mmImg.data || !mmBasis) return;
  const tmp = (updateMinimap as any).tmp || ((updateMinimap as any).tmp = document.createElement('canvas'));
  tmp.width = mmImg.data.width; tmp.height = mmImg.data.height; tmp.getContext('2d')!.putImageData(mmImg.data, 0, 0);
  g.clearRect(0, 0, W, H); g.imageSmoothingEnabled = true; g.drawImage(tmp, 0, 0, W, H);
  // view footprint at the target
  const B = mmBasis, toPx = (p: V3) => { const q = sub(p, B.c); return [(dot(q, B.u) / B.ext + 1) / 2 * W, (1 - dot(q, B.v) / B.ext) / 2 * H]; };
  const aspect = canvas.clientWidth / canvas.clientHeight;
  const hh = cam.dist * Math.tan(cam.fovy / 2), hw = hh * aspect;
  const corners = [[-1, -1], [1, -1], [1, 1], [-1, 1]].map(([a, b]) => toPx(add(cam.target, add(scale(cam.right, a * hw), scale(cam.up, b * hh)))));
  g.strokeStyle = '#FF9400'; g.lineWidth = 3; g.beginPath();
  corners.forEach(([x, y], i) => (i ? g.lineTo(x, y) : g.moveTo(x, y))); g.closePath(); g.stroke();
  const [tx, ty] = toPx(cam.target);
  g.fillStyle = '#FF9400'; g.beginPath(); g.arc(tx, ty, 5, 0, 7); g.fill();
  $('mmLabel').textContent = fs.cutOn ? `plane ${Math.round(d)} nm` : 'view plane';
}

// ---------------------------------------------------------------- screenshot
function screenshot() {
  const sc = fs.renderScale;
  R.resize(2 * sc);
  fs.renderScale = 2 * sc;
  renderNow(true);
  canvas.toBlob((b) => {
    if (b) {
      const a = document.createElement('a'); a.href = URL.createObjectURL(b);
      a.download = `yeast_cell_${new Date().toISOString().replace(/[:.]/g, '-').slice(0, 19)}.png`; a.click();
      setTimeout(() => URL.revokeObjectURL(a.href), 5000);
    }
    fs.renderScale = sc; R.resize(sc); dirty = true;
  }, 'image/png');
  toast('Screenshot saved (2× resolution)');
}

// ---------------------------------------------------------------- HUD / scale bar
function updateHUD() {
  const tt = T.totals(), s = R.stats;
  const fps = frameTimes.length > 3 ? 1000 / (frameTimes.reduce((a, b) => a + b) / frameTimes.length) : 0;
  const pend = T.pendingForView();
  $('hud').innerHTML = `<b>${fmt(tt.inst)}</b> molecules in ${tt.tiles} tiles${pend ? ` · <span style="color:#FF9400">streaming ${pend}</span>` : ''}${T.starved ? ' · <span style="color:#FF9400">memory budget reached (thinner slab loads more)</span>' : ''}<br>
    drawn: <b>${fmt(s.ninst)}</b> in view · ${fmt(s.l0 * 64)} residue beads (≤) · ${fmt(s.l1 * 16)} LOD1 · ${fmt(s.l2)} LOD2 · ${fmt(s.lip)} lipids · ${fmt(s.fre)} polymer beads<br>
    ${fps ? `${fps.toFixed(0)} fps · ` : ''}LOD bead ≥ ${fs.l0px.toFixed(1)} px · GPU pools ${R.gpuMemoryMB().toFixed(0)} MB · ${(T.bytesLoaded / 1048576).toFixed(0)} MB streamed${R.adapterInfo ? ' · ' + R.adapterInfo : ''}`;
  // scale bar
  const dpr = R.width / canvas.clientWidth;
  const pxPerNm = cam.pxPerNmAtUnit(R.height) / cam.dist / dpr;
  const target = 120 / pxPerNm;
  const p10 = Math.pow(10, Math.floor(Math.log10(target)));
  const nice = [1, 2, 5, 10].map(m => m * p10).reduce((a, b) => Math.abs(b - target) < Math.abs(a - target) ? b : a);
  ($('scalebar').querySelector('.bar') as HTMLElement).style.width = `${nice * pxPerNm}px`;
  ($('scalebar').querySelector('span') as HTMLElement).textContent = fmtNm(nice);
}

// ---------------------------------------------------------------- main loop
let gpuMs = 0, gpuBusy = false;
function renderNow(force = false) {
  const t0 = performance.now();
  R.frame(cam, fs, geo.G, T.meshDraws(), force || t0 - R.lastStatsT > 400);
  R.lastStatsT = performance.now();
  gpuBusy = true;
  R.device.queue.onSubmittedWorkDone().then(() => { gpuMs = performance.now() - t0; gpuBusy = false; });
}
let hudT = 0;
function tileTick() {
  if (sim?.active) return;                  // simulation mode: streaming paused, the box's tiles stay loaded
  T.update(cam, canvas.clientWidth / canvas.clientHeight, fs.cutOn, fs.plane, fs.slab, fs.organelle || fs.memColors, R.height);
  lastUpdate = performance.now();
  if (T.inflight || T.detInflight) dirty = true;
}
// streaming must not be throttled by requestAnimationFrame when frames are slow
setInterval(() => { if (T && performance.now() - lastUpdate > 150) tileTick(); }, 150);
function loop(now: number) {
  requestAnimationFrame(loop);
  if (simPlacing && sim && !sim.active && sim.box) {
    const c = sim.box.center, t = cam.target;
    if (Math.abs(c[0] - t[0]) + Math.abs(c[1] - t[1]) + Math.abs(c[2] - t[2]) > 0.01) { sim.showBox([...t] as V3, sim.box.size); dirty = true; }
  }
  const dt = lastFrameT ? now - lastFrameT : 16; lastFrameT = now;
  // keyboard flight
  const keys: Set<string> = (setupInput as any).keys;
  if (keys && keys.size) {
    const sp = cam.dist * 0.9 * dt / 1000;
    let mv: V3 = [0, 0, 0];
    if (keys.has('w')) mv = add(mv, scale(cam.fwd, sp)); if (keys.has('s')) mv = add(mv, scale(cam.fwd, -sp));
    if (keys.has('a')) mv = add(mv, scale(cam.right, -sp)); if (keys.has('d')) mv = add(mv, scale(cam.right, sp));
    if (keys.has('q')) mv = add(mv, scale(cam.up, -sp)); if (keys.has('e')) mv = add(mv, scale(cam.up, sp));
    if (keys.has('z')) cam.roll(-dt * 0.0015); if (keys.has('x')) cam.roll(dt * 0.0015);
    if (len(mv) > 0 || keys.has('z') || keys.has('x')) { cam.fly(mv); dirty = true; if (planeFollow) alignToView(); }
  }
  if (cam.step(now)) { dirty = true; if (planeFollow && !tour) alignToView(); }
  if (sweep) {
    const e = Math.min(1, (now - sweep.t0) / sweep.dur);
    const n = planeN(), d = sweep.d0 + (sweep.d1 - sweep.d0) * e;
    cam.target = add(cam.target, scale(n, d - fs.plane[3])); fs.plane[3] = d; syncPlaneUI();
    if (e >= 1) { sweep = null; $('sweepBtn').classList.remove('on'); }
    dirty = true;
  }
  tourStep(now);
  if (fs.motion > 0) dirty = true;          // illustrative motion: redraw every frame
  if (sim?.active && ((sim.running && sim.replayIdx < 0) || sim.playing || sim.computing)) dirty = true;
  if (now - lastUpdate > 90) tileTick();
  if (dirty && !gpuBusy) {
    dirty = false;
    renderNow();
    frameTimes.push(dt); if (frameTimes.length > 30) frameTimes.shift();
    // adaptive quality while interacting
    if (quality === 'auto' && frameTimes.length >= 20) {
      const avg = frameTimes.reduce((a, b) => a + b) / frameTimes.length;
      if (avg > 45 && fs.l0px < 3.5) { fs.l0px = Math.min(3.5, fs.l0px * 1.15); frameTimes.length = 0; }
      else if (avg < 22 && fs.l0px > 0.9) { fs.l0px = Math.max(0.9, fs.l0px / 1.1); frameTimes.length = 0; }
    }
    setTimeout(() => { dirty = dirty || R.statsBusy; }, 0);
  } else if (!dirty) frameTimes.length = 0;
  if (now - hudT > 300) { hudT = now; updateHUD(); updateMinimap(now); $('pillCut').classList.toggle('on', fs.cutOn); if (selType >= 0 && T.inflight) renderSelType(); }
}

boot();

// ---------------------------------------------------------------- simulation mode (sim/)
let simCut = true;
let simPlacing = false;      // the simulation box follows the view centre; clicks move the view centre
function setupSim() {
  const size = () => Number(($('simSize') as HTMLInputElement).value);
  const mode = () => ($('simMode') as HTMLSelectElement).value as 'schematic' | 'mid' | 'residue';
  const status = (msg: string) => { $('simStats').innerHTML = msg; };
  $('simSize').oninput = () => {
    $('simSizev').textContent = `${size()} nm`;
    if (sim.box && !sim.active) sim.showBox(sim.box.center, size()).then(() => { dirty = true; });
  };
  $('simSize').onchange = () => { if (sim.active) startSim(true); };
  const setPlacing = (on: boolean) => {
    simPlacing = on;
    $('simPlace').textContent = on ? 'Done placing' : 'Place box';
    $('simPlace').classList.toggle('on', on);
    if (on) status('The box sits at the centre of the view: pan, orbit and zoom as usual, or click a point to bring it to the centre. Then press Start.');
    else if (!sim.active) status('');
  };
  $('simPlace').onclick = async () => {
    if (simPlacing) { setPlacing(false); return; }
    if (sim.active) return;
    await sim.showBox([...cam.target] as V3, size());
    // frame the box (whole-cell views would show it as a dot)
    const st = cam.state(); st.dist = Math.min(Math.max(st.dist, 1.6 * size()), 4 * size()); cam.flyTo(st, 500);
    setPlacing(true); dirty = true;
  };
  const syncRun = () => {
    const s = sim.set;
    ($('simDt') as HTMLInputElement).value = String(Math.log10(s.dt)); $('simDtv').textContent = `${s.dt < 0.1 ? s.dt.toFixed(3) : s.dt.toFixed(2)} ns`;
    const dm = dtMax(s.mode, s.eps);
    ($('simDt') as HTMLInputElement).max = String(Math.log10(dm));
    $('simDt').title = `Brownian time step; at most ${dm < 0.1 ? dm.toFixed(3) : dm.toFixed(2)} ns for this model and stickiness (larger steps overestimate diffusion)`;
    ($('simSpf') as HTMLInputElement).value = String(s.stepsPerFrame); $('simSpfv').textContent = `${s.stepsPerFrame}`;
    $('simEpsv').textContent = `${s.eps.toFixed(2)} kT`; $('simViscv').textContent = `${s.viscosity}×`;
    $('simPlay').textContent = sim.running ? 'Pause' : 'Run live';
    $('simPlayRec').textContent = sim.playing ? 'Stop playback' : 'Play recording';
    $('simAhead').textContent = sim.computing ? 'Cancel' : 'Compute ahead';
  };
  const epsCeiling = () => {    // stronger stickiness is unstable in the finer models (params.STICKY_MAX)
    const e = $('simEps') as HTMLInputElement, mx = STICKY_MAX[mode()];
    e.max = String(mx); if (Number(e.value) > mx) e.value = String(mx);
    e.title = `Depth of a short-range attraction between molecular surfaces (0: excluded volume only); at most ${mx} kT in this model`;
  };
  epsCeiling();
  $('simMode').onchange = () => { epsCeiling(); ($('simLen') as HTMLSelectElement).value = mode() === 'residue' ? '200' : mode() === 'mid' ? '1000' : '5000'; if (sim.active) startSim(true); };
  const inter = () => ({ elec: ($('simElec') as HTMLInputElement).checked, ionic: Number(($('simIonic') as HTMLInputElement).value),
    specific: ($('simSpec') as HTMLInputElement).checked, epsGo: Number(($('simEpsGo') as HTMLInputElement).value), mrna: ($('simMrna') as HTMLInputElement).checked });
  const interLabels = () => { const v = inter(); $('simIonicv').textContent = `${v.ionic.toFixed(2)} M`; $('simEpsGov').textContent = `${v.epsGo} kT`; };
  for (const id of ['simElec', 'simIonic', 'simSpec', 'simEpsGo']) $(id).oninput = () => { interLabels(); sim.applySettings(inter()); };
  $('simMrna').onchange = () => { if (sim.active) startSim(true); };
  $('simBoundary').onchange = () => { if (sim.active) startSim(true); };
  $('simMemFluid').onchange = () => { if (sim.active) startSim(true); };
  $('simMemCal').onchange = () => { if (sim.active) startSim(true); };
  $('simHow').onchange = () => { if (sim.active) startSim(true); };
  ($('simLen') as HTMLSelectElement).value = '5000';
  const prog = (f: number | null) => {
    const el = $('simProg');
    el.style.display = f === null ? 'none' : 'block';
    if (f !== null) (el.firstElementChild as HTMLElement).style.width = `${(100 * f).toFixed(1)}%`;
  };
  // compute-ahead: progress bar, then play the recording in a loop
  const runAhead = async () => {
    const len = Number(($('simLen') as HTMLSelectElement).value);
    const eng = sim.eng;
    prog(0); const p = sim.computeAhead(len, 300, f => { prog(f); dirty = true; }); syncRun();
    await p;
    if (sim.eng !== eng) return;      // restarted meanwhile
    prog(null);
    if (sim.active && sim.frames.length > 1) sim.play(true);
    syncRun(); dirty = true;
  };
  $('simAhead').onclick = () => { if (sim.computing) sim.abort = true; else runAhead(); };
  $('simPlayRec').onclick = () => { sim.play(!sim.playing); syncRun(); dirty = true; };
  /** start, or restart in place after a change of model, boundary, run mode or size (same box centre) */
  let starting = false;
  const startSim = async (restart = false) => {
    if (starting || (sim.active && !restart)) return;
    starting = true;
    setPlacing(false);
    const first = !sim.active;
    if (sim.active) { sim.abort = true; sim.stop(); prog(null); }
    await sim.showBox(sim.box ? sim.box.center : [...cam.target] as V3, size());
    sim.applySettings((await import('./sim/params')).defaultSettings(mode()));
    sim.applySettings({ eps: Number(($('simEps') as HTMLInputElement).value), viscosity: Number(($('simVisc') as HTMLInputElement).value),
      hydro: ($('simHI') as HTMLInputElement).checked, boundary: ($('simBoundary') as HTMLSelectElement).value as any,
      memFluid: ($('simMemFluid') as HTMLInputElement).checked, memCalibrated: ($('simMemCal') as HTMLInputElement).checked, ...inter() });
    ($('simStart') as HTMLButtonElement).disabled = true;
    try {
      await sim.start(status);
      if (first) {
        simCut = fs.cutOn; fs.cutOn = false; ($('cutOn') as HTMLInputElement).checked = false;   // the whole box in view
        const st = cam.state(); st.target = [...sim.box!.center] as V3; st.dist = 2.2 * sim.box!.size; cam.flyTo(st, 800);   // frame the box
      }
      fs.simMode = true; $('simRun').style.display = ''; $('simStart').textContent = 'Restart'; syncRun();
      if (($('simHow') as HTMLSelectElement).value === 'ahead') { sim.running = false; runAhead(); }
    } catch (e: any) {
      status(`<span style="color:#c0392b">${e.message ?? e}</span>`);
      if (!first) { fs.simMode = false; $('simRun').style.display = 'none'; fs.cutOn = simCut; ($('cutOn') as HTMLInputElement).checked = simCut; }
    }
    ($('simStart') as HTMLButtonElement).disabled = false;
    starting = false; dirty = true;
  };
  $('simStart').onclick = () => startSim(sim.active);
  $('simDt').oninput = () => { sim.applySettings({ dt: Math.pow(10, Number(($('simDt') as HTMLInputElement).value)) }); syncRun(); };
  $('simSpf').oninput = () => { sim.applySettings({ stepsPerFrame: Number(($('simSpf') as HTMLInputElement).value) }); syncRun(); };
  $('simEps').oninput = () => { sim.applySettings({ eps: Number(($('simEps') as HTMLInputElement).value) }); syncRun(); };
  $('simHI').onchange = () => { sim.applySettings({ hydro: ($('simHI') as HTMLInputElement).checked }); };
  $('simVisc').oninput = () => { sim.applySettings({ viscosity: Number(($('simVisc') as HTMLInputElement).value) }); syncRun(); };
  $('simPlay').onclick = () => { if (sim.computing) return; sim.running = !sim.running; if (sim.running) sim.resumeLive(); syncRun(); dirty = true; };
  $('simLive').onclick = () => { if (sim.computing) return; sim.resumeLive(); sim.running = true; syncRun(); dirty = true; };
  $('simScrub').oninput = () => { if (sim.computing) return; sim.running = false; sim.replay(Number(($('simScrub') as HTMLInputElement).value)); syncRun(); dirty = true; };
  $('simExport').onclick = () => sim.export();
  $('simFps').oninput = () => { const f = Number(($('simFps') as HTMLInputElement).value); sim.setPlayFps(f); $('simFpsv').textContent = `${f} fps`; };
  $('simTrackClear').onclick = () => { sim.setTracked([]); dirty = true; };
  $('simExit').onclick = () => {
    sim.stop(); sim.hideBox(); fs.simMode = false; $('simRun').style.display = 'none'; status(''); prog(null); syncRun();
    $('simStart').textContent = 'Start'; setPlacing(false);
    fs.cutOn = simCut; ($('cutOn') as HTMLInputElement).checked = simCut; dirty = true;
  };
  const fmtT = (ns: number) => ns < 1e3 ? `${ns.toFixed(1)} ns` : ns < 1e6 ? `${(ns / 1e3).toFixed(2)} µs` : `${(ns / 1e6).toFixed(3)} ms`;
  const trackList = () => {
    const el = $('simTracks'), on = sim.active && sim.tracked.length > 0;
    $('simTrackClear').style.display = on ? '' : 'none';
    if (!on || !sim.eng) { el.innerHTML = ''; return; }
    const sys = sim.eng.sys;
    el.innerHTML = sim.trackStats().map((t, k) => `<div><span style="display:inline-block;width:10px;height:10px;border-radius:5px;background:${TRACK_COLORS[k]};margin-right:5px"></span><b>${lib.types[sys.type[t.b]].gene}</b>
      moved ${t.disp.toFixed(1)} nm${Number.isFinite(t.d) ? ` · ${t.d.toFixed(1)} µm²/s from successive frames` : ''}</div>`).join('');
  };
  sim.onUpdate = () => {
    trackList();
    const sc = $('simScrub') as HTMLInputElement;
    sc.max = String(Math.max(0, sim.frames.length - 1));
    if (sim.replayIdx < 0) sc.value = sc.max; else sc.value = String(sim.replayIdx);
  };
  let statTick = 0;
  setInterval(() => {
    if (!sim.active || !sim.eng || !sim.meta) return;
    const m = sim.meta, e = sim.eng;
    if (++statTick % 4 === 0) sim.pollGrid();
    const g = sim.grid;
    const gridNote = g && g.dropped > 0 ? ` · <span style="color:#c0392b">neighbour grid overflow: ${(g.dropped * 100).toFixed(2)}% of spheres missed (fullest cell ${g.max} > ${g.cap})</span>` : '';
    const t = sim.replayIdx >= 0 ? sim.frames[sim.replayIdx].tNs : e.timeNs;
    $('simScrubv').textContent = fmtT(t);
    const D = sim.diffusion().map(c => `${c.label}: <b>${c.d.toFixed(2)}</b> µm²/s (in solvent ${c.d0.toFixed(1)}; n=${c.n}; lag ${fmtT(c.lagNs)})`).join('<br>');
    const rate = sim.computing || (sim.playing && sim.aheadRate) ? sim.aheadRate : sim.stepsPerSec * sim.set.dt / 1e3;
    const len = Number(($('simLen') as HTMLSelectElement).value) / 1e3;
    const doing = sim.computing ? `Computing ahead (no drawing waits): ${rate.toFixed(2)} µs per second, about ${Math.max(0, (len - (e.timeNs - sim.frames[0].tNs) / 1e3) / Math.max(rate, 1e-6)).toFixed(0)} s left<br>`
      : sim.playing ? `Playing the recording (${sim.frames.length} frames, ${fmtT(sim.frames[sim.frames.length - 1].tNs - sim.frames[0].tNs)})<br>` : '';
    status(`${doing}<b>${fmtT(t)}</b> simulated · ${rate.toFixed(2)} µs per second · ${e.sys.n} molecules
      (${m.nMobile} free${m.nMem2d ? ` (${m.nMem2d} of them in membranes)` : ''}${m.nShell ? `, ${m.nShell} tethered in the boundary shell` : ''}, ${m.nFixed} fixed: ${sim.set.memFluid ? '' : 'membrane proteins, '}filaments, chromatin, nuclear pores)
      · boundary <b>${m.boundary === 'periodic' ? 'periodic' : m.boundary === 'walls' ? 'closed walls' : 'embedded'}</b>${m.boundary === 'periodic' && m.hasMembrane ? ' <span style="color:#c0392b">(membranes do not continue across the faces)</span>' : ''}<br>
      Interactions: ${sim.set.elec ? (m.charges ? `electrostatics on ${m.nCharged.toLocaleString()} charged spheres (box net charge ${m.netCharge.toFixed(0)} e; Debye length ${(0.303 / Math.sqrt(sim.set.ionic)).toFixed(2)} nm)` : 'no charge data in this build') : 'electrostatics off'}
      · ${sim.set.specific ? (m.contacts ? `specific binding: ${m.nNativePairs} predicted complex${m.nNativePairs === 1 ? '' : 'es'} can form between ${m.nNativeBodies} molecules (${sim.set.epsGo} kT per interface)` : 'no contact data in this build') : 'specific binding off'}
      · ${sim.set.mrna ? `${m.nLinks} mRNA link${m.nLinks === 1 ? '' : 's'} joining ${m.nPolyRibo} polysome ribosomes` : 'mRNA links off'}<br>
      <span style="color:#8a8f98">${e.nSph.toLocaleString()} collision spheres · grid ${e.gridInfo.h.toFixed(1)} nm cells (capacity ${e.gridInfo.cap})${e.gridInfo.large ? ` + ${e.gridInfo.hB.toFixed(1)} nm cells for the ${e.gridInfo.large.toLocaleString()} spheres over ${e.gridInfo.split.toFixed(2)} nm` : ''}</span><br>
      Volume fraction ${(e.sys.phi * 100).toFixed(1)}% (hydrated)${gridNote}${sim.set.hydro ? ` · mean-field hydrodynamics: D × ${hiTrans(e.sys.phi).toFixed(2)}, D<sub>r</sub> × ${hiRot(e.sys.phi).toFixed(2)}` : ''}<br>
      Diffusion (time-averaged MSD): ${D || '…'}<br><span style="color:#8a8f98">Measured in yeast: GFP 8.6 µm²/s (fission-yeast cytoplasm, FRAP);
      budding yeast 8.7–21 µm²/s (FCS, lab-dependent); 40 nm GEM particles about 0.3 µm²/s (particle tracking). Excluded volume,
      electrostatics and mean-field hydrodynamics leave GFP-sized molecules about 2.7× too fast; the default Viscosity of
      2.75× is that calibration (fitted on training papers, checked on held-out ones; Model details, section 14.5).
      ${STICKY_NOTE[sim.set.mode]}</span>`);
  }, 500);
}
