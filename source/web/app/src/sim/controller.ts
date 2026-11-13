// Simulation mode: choose a box, extract its molecules, run the rigid-body Brownian dynamics on the GPU while drawing
// it, record frames for replay, export MD-ready files. The whole-cell rendering and tile streaming pause meanwhile.
import type { Renderer, SimHook } from '../renderer';
import type { TileManager } from '../tiles';
import type { Library } from '../data';
import type { V3 } from '../math';
import { FILL, DRAW_SPHERES } from '../shaders';
import { SimEngine, slimG0 } from './engine';
import { extract, boxTiles, loadSimData, type SimBox, type SimMeta } from './extract';
import { defaultSettings, dtMax, type SimSettings, type TypeHydro } from './params';
import { exportFiles } from './export';

const BOX_DRAW = FILL.slice(0, FILL.indexOf('@group(0) @binding(5)')) + DRAW_SPHERES.slice(DRAW_SPHERES.indexOf('struct VOut'), DRAW_SPHERES.indexOf('@group(1) @binding(0)')) + /* wgsl */ `
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
`;

// paths of tracked molecules: beads along the segments between recorded centre positions, pixel-wide like the box
// outline, drawn on top of the molecules; a segment that crosses a periodic face is left out (the molecule re-enters
// on the other side)
const TRACK_STRIDE = 640;    // frames per molecule in the path buffer (recordings keep at most 600)
const TRACK_SUB = 12;        // beads per segment
const TRACK_MAX = 8;
export const TRACK_COLORS = ['#e6194b', '#3cb44b', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#9a6324'];
const TRACK_DRAW = FILL.slice(0, FILL.indexOf('@group(0) @binding(5)')) + DRAW_SPHERES.slice(DRAW_SPHERES.indexOf('struct VOut'), DRAW_SPHERES.indexOf('@group(1) @binding(0)')) + /* wgsl */ `
@group(1) @binding(0) var<uniform> pu: vec4u;                 // segments drawn per molecule, stride, beads per segment, molecules
@group(1) @binding(1) var<storage, read> pts: array<vec4f>;   // per molecule, per frame: centre, w = periodic box length (0: none)
const COLS = array<vec3f, 8>(${TRACK_COLORS.map(c => `vec3f(${[1, 3, 5].map(i => (parseInt(c.slice(i, i + 2), 16) / 255).toFixed(3)).join(', ')})`).join(', ')});
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
`;

export interface Frame { tNs: number; pos: Float32Array; quat: Float32Array; upos?: Float32Array }   // upos: unwrapped (periodic)

export class SimController implements SimHook {
  R: Renderer; T: TileManager; lib: Library;
  box: SimBox | null = null; set: SimSettings = defaultSettings();
  eng: SimEngine | null = null; meta: SimMeta | null = null;
  running = false; active = false; replayIdx = -1;
  frames: Frame[] = []; recording = true; maxFrames = 600; stepsSinceRec = 0; recEvery = 200;
  // compute-ahead (batches back to back, no waiting for the display) and playback of the recording
  computing = false; abort = false; aheadRate = 0;
  playing = false; playFps = 24; private playStart = 0; private playFrom = 0;
  hydroCache = new Map<number, TypeHydro>();
  boxPipe: GPURenderPipeline | null = null; boxBG: GPUBindGroup | null = null; boxG0: GPUBindGroup | null = null; boxBuf: GPUBuffer;
  onUpdate: () => void = () => {};
  stepsPerSec = 0; private lastT = 0; private lastSteps = 0; totalSteps = 0;
  // tracked molecules (body indices) and their paths
  tracked: number[] = [];
  trackPipe: GPURenderPipeline | null = null; trackBG: GPUBindGroup | null = null; trackUni: GPUBuffer; trackPts: GPUBuffer;
  constructor(R: Renderer, T: TileManager, lib: Library) {
    this.R = R; this.T = T; this.lib = lib;
    this.boxBuf = R.device.createBuffer({ size: 16, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
    this.trackUni = R.device.createBuffer({ size: 16, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
    this.trackPts = R.device.createBuffer({ size: 16 * TRACK_STRIDE * TRACK_MAX, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST });
  }

  /** track molecules (body indices; at most 8): the others are faded and the tracked ones' paths drawn */
  setTracked(ids: number[]) {
    this.tracked = [...new Set(ids)].slice(0, TRACK_MAX);
    this.eng?.setTracked(this.tracked);
    this.uploadTracks();
    this.onUpdate();
  }
  toggleTrack(b: number) {
    if (!this.eng || b < 0 || b >= this.eng.sys.n) return;
    this.setTracked(this.tracked.includes(b) ? this.tracked.filter(x => x !== b) : [...this.tracked, b]);
  }
  private uploadTracks() {
    const eng = this.eng;
    if (!eng || !this.tracked.length) return;
    const F = Math.min(this.frames.length, TRACK_STRIDE), sys = eng.sys;
    const L = sys.boundary === 'periodic' ? sys.hi[0] - sys.lo[0] : 0;
    const a = new Float32Array(4 * TRACK_STRIDE * this.tracked.length);
    this.tracked.forEach((b, k) => {
      for (let j = 0; j < F; j++) {
        const p = this.frames[j].pos, o = 4 * (k * TRACK_STRIDE + j);
        a[o] = p[3 * b]; a[o + 1] = p[3 * b + 1]; a[o + 2] = p[3 * b + 2]; a[o + 3] = L;
      }
    });
    this.R.device.queue.writeBuffer(this.trackPts, 0, a);
  }
  private async ensureTrackPipe() {
    if (this.trackPipe) return;
    const d = this.R.device;
    const m = d.createShaderModule({ code: TRACK_DRAW, label: 'simTrack' });
    const info = await m.getCompilationInfo();
    for (const e of info.messages.filter(x => x.type === 'error')) console.error(`[simTrack] ${e.lineNum}:${e.linePos} ${e.message}`);
    const l1 = d.createBindGroupLayout({ entries: [
      { binding: 0, visibility: GPUShaderStage.VERTEX, buffer: { type: 'uniform' } },
      { binding: 1, visibility: GPUShaderStage.VERTEX, buffer: { type: 'read-only-storage' } }] });
    const g0 = slimG0(this.R);
    this.trackPipe = await d.createRenderPipelineAsync({
      layout: d.createPipelineLayout({ bindGroupLayouts: [g0.layout, l1] }),
      vertex: { module: m, entryPoint: 'vsPath' },
      fragment: { module: m, entryPoint: 'fsPath', targets: [{ format: 'rgba8unorm' }, { format: 'r32uint' }] },     // id 0: no outlines over the path
      primitive: { topology: 'triangle-list', cullMode: 'none' },
      depthStencil: { format: 'depth32float', depthWriteEnabled: false, depthCompare: 'always' },     // on top: the paths stay visible through the crowd
    });
    this.trackBG = d.createBindGroup({ layout: l1, entries: [{ binding: 0, resource: { buffer: this.trackUni } }, { binding: 1, resource: { buffer: this.trackPts } }] });
    this.trackG0 = g0.bg;
  }
  private trackG0: GPUBindGroup | null = null;
  /** per tracked molecule: displacement from the start of the recording (nm, unwrapped) and its own diffusion
   *  estimate from successive recorded frames (um^2/s) */
  trackStats(): { b: number; disp: number; d: number }[] {
    const F = this.frames, last = this.replayIdx >= 0 ? this.replayIdx : F.length - 1;
    return this.tracked.map(b => {
      if (last < 1) return { b, disp: 0, d: NaN };
      const P = (f: Frame) => f.upos ?? f.pos;
      let s = 0, t = 0;
      for (let j = 1; j <= last; j++) {
        const p0 = P(F[j - 1]), p1 = P(F[j]);
        s += (p1[3 * b] - p0[3 * b]) ** 2 + (p1[3 * b + 1] - p0[3 * b + 1]) ** 2 + (p1[3 * b + 2] - p0[3 * b + 2]) ** 2; t += F[j].tNs - F[j - 1].tNs;
      }
      const a = P(F[0]), z = P(F[last]);
      return { b, disp: Math.hypot(z[3 * b] - a[3 * b], z[3 * b + 1] - a[3 * b + 1], z[3 * b + 2] - a[3 * b + 2]), d: t > 0 ? s / (6 * t) * 1e3 : NaN };
    });
  }
  /** playback rate (recorded frames per second); continues from the frame shown */
  setPlayFps(fps: number) {
    if (this.playing && this.replayIdx >= 0) { this.playFrom = this.replayIdx; this.playStart = performance.now(); }
    this.playFps = Math.max(0.5, fps);
  }

  async showBox(center: V3, size: number) {
    this.box = { center: [...center] as V3, size };
    this.R.device.queue.writeBuffer(this.boxBuf, 0, new Float32Array([...center, size / 2]));
    if (!this.boxPipe) {
      const d = this.R.device;
      const m = d.createShaderModule({ code: BOX_DRAW, label: 'simBox' });
      const l1 = d.createBindGroupLayout({ entries: [{ binding: 0, visibility: GPUShaderStage.VERTEX, buffer: { type: 'uniform' } }] });
      const g0 = slimG0(this.R); this.boxG0 = g0.bg;
      this.boxPipe = await d.createRenderPipelineAsync({
        layout: d.createPipelineLayout({ bindGroupLayouts: [g0.layout, l1] }),
        vertex: { module: m, entryPoint: 'vsEdge' },
        // outline on top of everything, not pickable
        fragment: { module: m, entryPoint: 'fsMain', targets: [{ format: 'rgba8unorm' }, { format: 'r32uint', writeMask: 0 }] },
        primitive: { topology: 'triangle-list', cullMode: 'none' },
        depthStencil: { format: 'depth32float', depthWriteEnabled: false, depthCompare: 'always' },
      });
      this.boxBG = d.createBindGroup({ layout: l1, entries: [{ binding: 0, resource: { buffer: this.boxBuf } }] });
    }
    this.R.sim = this;
  }
  hideBox() { this.box = null; if (!this.active) this.R.sim = null; }

  async start(progress: (msg: string) => void) {
    if (!this.box) throw new Error('no box');
    progress('Loading the tiles of the box');
    try { await this.T.ensureLoaded(boxTiles(this.T, this.box)); } catch (e) { this.T.unpin(); throw e; }
    progress('Loading charges and native contacts');
    const data = await loadSimData();
    progress('Building rigid bodies and diffusion coefficients');
    const { sys, meta } = extract(this.T, this.lib, this.box, this.set, this.hydroCache, data);
    if (!sys.n) { this.T.unpin(); throw new Error('the box is empty'); }
    progress(`Compiling the simulation (${sys.n} molecules, ${sys.sphLocal.length / 4} collision spheres)`);
    this.eng?.destroy();
    this.eng = new SimEngine(this.R, sys, this.set);
    await this.eng.init();
    this.meta = meta; this.frames = []; this.replayIdx = -1; this.totalSteps = 0; this.stepsSinceRec = 0; this.grid = null; this.tracked = [];
    await this.ensureTrackPipe();
    this.recEvery = this.set.recordEvery; this.playing = false;
    this.maxFrames = Math.max(100, Math.min(600, Math.floor(60e6 / (28 * sys.n))));   // ~60 MB of recorded poses
    this.active = true; this.running = true; this.R.sim = this;
    this.frames.push({ tNs: 0, pos: sys.pos.slice(), quat: sys.quat.slice() });
    this.lastT = performance.now(); this.lastSteps = 0;
  }

  stop() {
    this.running = false; this.active = false; this.playing = false; this.abort = true; this.T.unpin();
    this.eng?.destroy(); this.eng = null; this.meta = null; this.frames = [];
    this.R.sim = this.box ? this : null;
  }

  /** settings; the time step is kept at or below dtMax for the model and stickiness */
  applySettings(s: Partial<SimSettings>) {
    this.set = { ...this.set, ...s };
    this.set.dt = Math.min(this.set.dt, dtMax(this.set.mode, this.set.eps));
    this.eng?.setSettings(this.set);
  }

  /** neighbour-grid load (diagnostic), refreshed by pollGrid() */
  grid: { max: number; cap: number; dropped: number } | null = null;
  private polling = false;
  pollGrid() {
    if (!this.eng || this.polling) return;
    this.polling = true; const eng = this.eng;
    eng.gridLoad().then(g => {
      if (eng !== this.eng) return;
      this.grid = g;
      if (g.dropped > 0) eng.growGrid();         // between frames: no encoder holds the old buffer
    }).catch(() => {}).finally(() => { this.polling = false; });
  }

  // ---- SimHook
  clipBox() { return this.active && this.box ? [...this.box.center, this.box.size / 2] : null; }
  compute(enc: GPUCommandEncoder) {
    if (!this.eng) return;
    if (this.playing) { this.advancePlayback(); return; }
    if (!this.running || this.replayIdx >= 0 || this.computing) return;
    const n = Math.max(1, Math.min(512, Math.round(this.set.stepsPerFrame)));
    this.eng.step(enc, n);
    this.totalSteps += n; this.stepsSinceRec += n;
    const now = performance.now();
    if (now - this.lastT > 1000) { this.stepsPerSec = (this.totalSteps - this.lastSteps) * 1000 / (now - this.lastT); this.lastT = now; this.lastSteps = this.totalSteps; }
    if (this.recording && this.stepsSinceRec >= this.recEvery && !this.reading) {
      this.stepsSinceRec = 0; this.reading = true;
      const t = this.eng.timeNs, eng = this.eng;
      queueMicrotask(async () => {      // after this frame's command buffer is submitted
        try { const p = await eng.readPoses(); if (eng === this.eng) this.record({ tNs: t, ...p }); } catch { /* engine replaced */ }
        this.reading = false;
      });
    }
  }
  private reading = false;
  /** add a frame; past maxFrames drop every other one (uniform spacing for playback) and record half as often */
  private record(f: Frame) {
    const prev = this.frames[this.frames.length - 1], sys = this.eng?.sys;
    if (prev && sys && sys.boundary === 'periodic') {     // unwrap across the periodic faces (frames are far closer than L/2)
      const L = sys.hi[0] - sys.lo[0], pu = prev.upos ?? prev.pos, u = new Float32Array(f.pos.length);
      for (let i = 0; i < u.length; i++) { let d = f.pos[i] - prev.pos[i]; d -= L * Math.round(d / L); u[i] = pu[i] + d; }
      f.upos = u;
    }
    this.frames.push(f);
    if (this.frames.length > this.maxFrames) { this.frames = this.frames.filter((_, k) => k % 2 === 0); this.recEvery *= 2; }
    if (this.tracked.length) this.uploadTracks();
    this.onUpdate();
  }

  /** Simulate durNs as fast as the GPU allows, recording about nFrames frames, then leave the recording ready to play.
   *  Steps go out in command buffers of about 0.1 s of GPU work (at most 256 steps); at each recording point a pose copy into a staging buffer is
   *  queued in the same command buffer and read back asynchronously, with several submissions in flight, so the GPU
   *  never waits for the page. The recording restarts at the current state. */
  async computeAhead(durNs: number, nFrames: number, progress: (frac: number) => void) {
    const eng = this.eng;
    if (!eng || this.computing) return;
    const d = this.R.device, size = 32 * eng.sys.n;
    this.playing = false; if (this.replayIdx >= 0) this.resumeLive();
    this.running = false; this.computing = true; this.abort = false;
    while (this.reading) await new Promise(r => setTimeout(r, 5));
    const total = Math.max(1, Math.round(durNs / this.set.dt));
    this.recEvery = Math.max(1, Math.round(total / nFrames));
    this.frames = [{ tNs: eng.timeNs, ...(await eng.readPoses()) }]; this.uploadTracks(); this.onUpdate();
    const pool: GPUBuffer[] = [], all: GPUBuffer[] = [];
    const staging = () => { const b = pool.pop() ?? d.createBuffer({ size, usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST }); if (!all.includes(b)) all.push(b); return b; };
    const pending: { buf: GPUBuffer | null; t: number; p: Promise<void> }[] = [];
    const drain = async (keep: number) => {
      while (pending.length > keep) {
        const x = pending.shift()!; await x.p;
        if (!x.buf) continue;
        const a = new Float32Array(x.buf.getMappedRange().slice(0)); x.buf.unmap(); pool.push(x.buf);
        if (eng === this.eng) this.record({ tNs: x.t, ...eng.unpackPoses(a) });
      }
    };
    let done = 0, since = 0, batches = 0, perBatch = 16;
    const t0 = performance.now();
    try {
      while (done < total && !this.abort && eng === this.eng) {
        // command buffers of about 0.1 s of GPU work (long ones can trip the GPU watchdog); every 8th, check the grid
        if (batches > 0 && batches % 8 === 0) {
          await drain(0);
          if (this.abort || eng !== this.eng) break;
          const g = await eng.gridLoad(); this.grid = g;
          if (g.dropped > 0) eng.growGrid();
          const rate = done / Math.max(1e-3, (performance.now() - t0) / 1000);
          perBatch = Math.max(4, Math.min(256, Math.round(rate * 0.1)));
          if (this.abort || eng !== this.eng) break;      // stopped or restarted while waiting
        }
        batches++;
        const enc = d.createCommandEncoder(), recs: { buf: GPUBuffer; t: number }[] = [];
        let inBatch = 0;
        while (inBatch < perBatch && done < total) {
          const n = Math.min(this.recEvery - since, total - done, perBatch - inBatch);
          eng.step(enc, n); done += n; since += n; inBatch += n; this.totalSteps += n;
          if (since >= this.recEvery || done >= total) { since = 0; const b = staging(); eng.copyPoses(enc, b); recs.push({ buf: b, t: eng.timeNs }); }
        }
        d.queue.submit([enc.finish()]);
        for (const r of recs) pending.push({ ...r, p: r.buf.mapAsync(GPUMapMode.READ) });
        if (!recs.length) pending.push({ buf: null, t: 0, p: d.queue.onSubmittedWorkDone() });
        await drain(8);
        if (this.abort || eng !== this.eng) break;
        this.aheadRate = done * this.set.dt / 1e3 / Math.max(1e-3, (performance.now() - t0) / 1000);
        progress(done / total);
      }
      await drain(0).catch(() => {});
    } finally {
      this.computing = false;
      for (const b of all) b.destroy();
    }
  }

  // ---- playback of the recording
  play(on: boolean) {
    if (on && this.frames.length > 1) {
      this.running = false; this.playing = true;
      this.playFrom = this.replayIdx >= 0 && this.replayIdx < this.frames.length - 1 ? this.replayIdx : 0;
      this.playStart = performance.now();
    } else this.playing = false;
  }
  private advancePlayback() {
    const n = this.frames.length;
    const i = (this.playFrom + Math.floor((performance.now() - this.playStart) * this.playFps / 1000)) % n;     // loops
    if (i !== this.replayIdx) {
      this.replayIdx = i; const f = this.frames[i];
      this.eng!.writePoses(f.pos, f.quat); this.onUpdate();
    }
  }
  draw(rp: GPURenderPassEncoder) {
    if (!this.eng || !this.active) return;
    this.eng.draw(rp);
    const F = Math.min(this.frames.length, TRACK_STRIDE);
    const nSeg = Math.min(F - 1, this.replayIdx >= 0 ? this.replayIdx : F - 1);
    if (this.tracked.length && nSeg >= 1 && this.trackPipe && this.trackBG && this.trackG0) {
      this.R.device.queue.writeBuffer(this.trackUni, 0, new Uint32Array([nSeg, TRACK_STRIDE, TRACK_SUB, this.tracked.length]));
      rp.setPipeline(this.trackPipe); rp.setBindGroup(0, this.trackG0); rp.setBindGroup(1, this.trackBG);
      rp.draw(6, this.tracked.length * nSeg * TRACK_SUB);
    }
  }
  drawOverlay(rp: GPURenderPassEncoder) {
    if (this.box && this.boxPipe && this.boxBG) {
      rp.setPipeline(this.boxPipe); rp.setBindGroup(0, this.boxG0!); rp.setBindGroup(1, this.boxBG); rp.draw(6, 12 * 480);
    }
  }

  // ---- replay
  replay(i: number) {
    if (!this.eng || !this.frames.length) return;
    this.playing = false;
    if (i < 0 || i >= this.frames.length - 1 && this.running) { this.replayIdx = -1; return; }
    this.replayIdx = Math.min(i, this.frames.length - 1);
    const f = this.frames[this.replayIdx];
    this.eng.writePoses(f.pos, f.quat);
  }
  resumeLive() {
    if (!this.eng) return;
    this.playing = false;
    const f = this.frames[this.frames.length - 1];
    if (this.replayIdx >= 0 && f) this.eng.writePoses(f.pos, f.quat);
    this.replayIdx = -1;
  }

  /** translational diffusion of the mobile molecules by size class, from the time-averaged mean squared displacement
   *  over the recording: D = MSD(tau)/(6 tau) at the longest lag tau whose rms displacement stays below a quarter of the
   *  free interior of the box (longer lags feel the confinement by the box faces). Replay: up to that frame. */
  diffusion(): { label: string; d: number; n: number; d0: number; lagNs: number }[] {
    if (!this.eng || !this.meta) return [];
    const last = this.replayIdx >= 0 && !this.playing ? this.replayIdx : this.frames.length - 1;
    if (last < 1) return [];
    const F = this.frames, sys = this.eng.sys, h = this.meta.hydro;
    const free = Math.max(5, this.meta.box.size - (sys.boundary === 'embedded' ? 2 * this.set.shell : 0)), lim2 = (free / 4) ** 2;
    const P = (f: Frame) => f.upos ?? f.pos;
    const classes = [{ label: 'R_h < 3 nm (GFP-sized)', lo: 0, hi: 3 }, { label: 'R_h 3-6 nm', lo: 3, hi: 6 }, { label: 'R_h > 6 nm (ribosomes, large complexes)', lo: 6, hi: 1e9 }];
    const lags: number[] = [];
    for (let k = 1; k <= last; k = Math.max(k + 1, Math.round(k * 1.5))) lags.push(k);
    return classes.map(c => {
      const ids: number[] = []; let d0 = 0;
      for (let i = 0; i < sys.n; i++) if (sys.mobile[i] && !sys.tether[i] && h[i].rh >= c.lo && h[i].rh < c.hi) { ids.push(i); d0 += h[i].dt; }
      if (!ids.length) return null;
      let best = { d: NaN, lag: 0 };
      for (const k of lags) {
        const starts = last - k + 1, stride = Math.max(1, Math.floor(starts / 40));
        let s = 0, m = 0, tau = 0;
        for (let a = 0; a + k <= last; a += stride) {
          const f0 = F[a], f1 = F[a + k], p0 = P(f0), p1 = P(f1); tau += f1.tNs - f0.tNs;
          for (const i of ids) {
            const dx = p1[3 * i] - p0[3 * i], dy = p1[3 * i + 1] - p0[3 * i + 1], dz = p1[3 * i + 2] - p0[3 * i + 2];
            s += dx * dx + dy * dy + dz * dz;
          }
          m++;
        }
        const msd = s / (m * ids.length); tau /= m;
        if (tau <= 0) continue;
        if (best.lag && msd > lim2) break;
        best = { d: msd / (6 * tau) * 1e3, lag: tau };       // nm^2/ns -> um^2/s
      }
      return { label: c.label, d: best.d, n: ids.length, d0: d0 / ids.length * 1e3 / this.set.viscosity, lagNs: best.lag };
    }).filter((c): c is NonNullable<typeof c> => !!c);
  }

  async export() {
    if (!this.eng || !this.meta) return;
    const f = this.replayIdx >= 0 ? this.frames[this.replayIdx] : await this.eng.readPoses();
    await exportFiles(this.eng, this.meta, this.lib, this.set, f.pos, f.quat, this.frames);
  }
}
