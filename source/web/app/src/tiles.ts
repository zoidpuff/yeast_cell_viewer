// Streaming tile manager: decides which 200 nm tiles are needed for the current view and cut, loads them
// through web workers, and manages GPU pool pages (instances, lipids, free spheres) with LRU eviction.
import { TileIndex, Library, tkey } from './data';
import { Renderer, MeshDraw } from './renderer';
import { PAGE, SPAGE } from './shaders';
import { Camera } from './camera';
import { V3, dot, sub, len } from './math';
import TileWorker from './worker.ts?worker&inline';

const TILE = 200, HALF_DIAG = 173.3, MOL_MARGIN = 70;

class PageAlloc {
  used: Uint8Array; n: number; count = 0; high = 0;
  constructor(n: number) { this.n = n; this.used = new Uint8Array(n); }
  alloc(k: number): number[] | null {
    if (k <= 0) return [];
    if (this.n - this.count < k) return null;
    const out: number[] = [];
    for (let i = 0; i < this.n && out.length < k; i++) if (!this.used[i]) { this.used[i] = 1; out.push(i); }
    this.count += k;
    this.high = Math.max(this.high, out[out.length - 1] + 1);
    return out;
  }
  free(p: number[]) {
    for (const i of p) { if (this.used[i]) { this.used[i] = 0; this.count--; } }
    while (this.high > 0 && !this.used[this.high - 1]) this.high--;
  }
  get freeCount() { return this.n - this.count; }
}

export interface LoadedTile {
  idx: number; slot: number; lo: V3; inst: Uint32Array; nInst: number; pages: number[];
  basePages: number[]; nBase: number; lipPages: number[]; detPages: number[]; nLip: number; nDet: number;
  vb: GPUBuffer | null; ib: GPUBuffer | null; nIdx: number; lastUsed: number;
  det: null | { n: V3; a: number; b: number; shell: boolean }; detPending: boolean; underMesh: boolean;
  baseFreeMissing: boolean;
}

export class TileManager {
  idx: TileIndex; r: Renderer; lib: Library; G: number[];
  workers: Worker[] = [];
  state: Uint8Array;
  loaded = new Map<number, LoadedTile>();
  instA: PageAlloc; lipA: PageAlloc; freeA: PageAlloc;
  slots: number[] = [];
  pageTile: Int32Array; pageLocal: Uint32Array;
  tileState = new Uint32Array(22 * 22 * 22);
  desired: number[] = []; desiredSet = new Set<number>();
  inflight = 0; maxInflight = 12; detInflight = 0; maxDetInflight = 3;
  frameNo = 0;
  onChange: () => void = () => {};
  onTileLoaded: ((t: LoadedTile) => void)[] = [];
  bytesLoaded = 0; tilesLoadedTotal = 0; starved = false; detailStats = { lipids: 0, glucan: 0, mannan: 0 };
  detailEnabled = true;
  lastLoadErr = '';

  constructor(idx: TileIndex, r: Renderer, lib: Library, geo: any, nWorkers: number, embedded?: Record<string, Uint8Array>) {
    this.idx = idx; this.r = r; this.lib = lib; this.G = geo.G;
    this.state = new Uint8Array(idx.n);
    this.instA = new PageAlloc(r.instPages); this.lipA = new PageAlloc(r.lipPages); this.freeA = new PageAlloc(r.freePages);
    for (let s = r.maxSlots - 1; s >= 0; s--) this.slots.push(s);
    this.pageTile = new Int32Array(r.instPages).fill(-1); this.pageLocal = new Uint32Array(r.instPages);
    const volG = Array.from(idx.volGlucan);
    const base = new URL('data/', document.baseURI).href;
    for (let w = 0; w < nWorkers; w++) {
      const wk: Worker = new TileWorker();
      wk.onmessage = (e) => this.onMessage(e.data);
      wk.postMessage({ cmd: 'init', base, G: geo.G, volGlucan: volG, embedded,
        wall: { glucanPerNm3: idx.wall.glucan_beads_per_nm3, mannanPerProtein: idx.wall.mannan_beads_per_protein } });
      this.workers.push(wk);
    }
  }

  tileCenter(i: number): V3 { return [this.idx.ix[i] * TILE + 100, this.idx.iy[i] * TILE + 100, this.idx.iz[i] * TILE + 100]; }

  // ---------------------------------------------------------------- desired set
  update(cam: Camera, aspect: number, cutOn: boolean, plane: number[], slab: number, organelle: boolean, heightPx: number) {
    this.frameNo++;
    if (organelle) return;
    const eye = cam.eye, fr = cam.frustum(aspect), n = this.idx.n;
    const R = this.G[0];
    const nn: V3 = [plane[0], plane[1], plane[2]];
    const ext = 100 * (Math.abs(nn[0]) + Math.abs(nn[1]) + Math.abs(nn[2]));
    const cand: { i: number; d: number }[] = [];
    for (let i = 0; i < n; i++) {
      const c = this.tileCenter(i);
      const rel = sub(c, eye);
      const rs = HALF_DIAG + MOL_MARGIN;
      let ok = true;
      for (let k = 0; k < 5; k++) if (fr[4 * k] * rel[0] + fr[4 * k + 1] * rel[1] + fr[4 * k + 2] * rel[2] + fr[4 * k + 3] < -rs) { ok = false; break; }
      if (!ok) continue;
      if (cutOn) {
        const s = dot(nn, c) - plane[3];
        if (s + ext + MOL_MARGIN < 0 || s - ext - MOL_MARGIN > slab) continue;
      } else {
        const rc = len(c);
        if (rc < R - 150 - HALF_DIAG || rc > R + HALF_DIAG) continue;
        // front hemisphere only (the rest is hidden behind the wall)
        const toEye = sub(eye, c);
        if (dot(c, toEye) < -rc * 260) continue;
      }
      cand.push({ i, d: len(rel) });
    }
    cand.sort((a, b) => a.d - b.d);
    let budget = this.r.instPages * PAGE * 0.94, slots = this.r.maxSlots - 8;
    const des: number[] = [];
    for (const c of cand) {
      const ni = this.idx.nInst[c.i];
      if (ni > budget || slots <= 0) { this.starved = true; break; }
      budget -= Math.ceil(ni / PAGE) * PAGE; slots--;
      des.push(c.i);
    }
    if (des.length === cand.length) this.starved = false;
    this.desired = des; this.desiredSet = new Set(des);
    for (const i of des) { const t = this.loaded.get(i); if (t) t.lastUsed = this.frameNo; }
    // request loads
    for (const i of des) {
      if (this.inflight >= this.maxInflight) break;
      if (this.state[i] !== 0) continue;
      this.state[i] = 1; this.inflight++;
      this.workers[i % this.workers.length].postMessage({ cmd: 'load', idx: i, ix: this.idx.ix[i], iy: this.idx.iy[i], iz: this.idx.iz[i] });
    }
    this.updateDetail(cam, cutOn, plane, slab, heightPx);
  }

  // ---------------------------------------------------------------- procedural detail near the camera
  updateDetail(cam: Camera, cutOn: boolean, plane: number[], slab: number, heightPx: number) {
    if (!this.detailEnabled) return;
    const eye = cam.eye, pxu = cam.pxPerNmAtUnit(heightPx);
    const nn: V3 = [plane[0], plane[1], plane[2]];
    const need = Math.min(slab, 250);
    const elig: { t: LoadedTile; d: number }[] = [];
    for (const i of this.desired) {
      const t = this.loaded.get(i);
      if (!t) continue;
      const hasMem = this.idx.nTri[i] > 0, hasWall = this.idx.volGlucan[i] > 0 || (this.idx.labelMask[i] & 2) !== 0;
      if (!hasMem && !hasWall) continue;
      const c = this.tileCenter(i);
      const dist = Math.max(1, len(sub(c, eye)) - HALF_DIAG);
      const pxPerNm = pxu / dist;
      if (pxPerNm * 0.55 < 0.8) continue;
      if (!cutOn && !hasWall) continue;
      elig.push({ t, d: dist });
    }
    elig.sort((a, b) => a.d - b.d);
    const keep = new Set(elig.slice(0, 40).map(e => e.t.idx));
    // drop detail from tiles that no longer qualify (frees pool space)
    for (const t of this.loaded.values()) if (t.det && !keep.has(t.idx) && !t.detPending) this.dropDetail(t);
    for (const { t } of elig.slice(0, 40)) {
      if (this.detInflight >= this.maxDetInflight) break;
      if (t.detPending) continue;
      const D = t.det;
      let ok = false;
      if (D) {
        if (!cutOn) ok = D.shell;
        else ok = !D.shell && dot(D.n, nn) > 0.99985 && D.a <= plane[3] + 0.5 && D.b >= plane[3] + need - 0.5;
      }
      if (ok) continue;
      const m = Math.max(10, 0.35 * need);
      const range = cutOn ? { n: nn, d: plane[3], lo: -m, hi: need + m, mode: 'slab' } : { n: nn, d: plane[3], lo: 0, hi: 0, mode: 'shell' };
      t.detPending = true; this.detInflight++;
      const key = cutOn ? `${plane[3] - m}|${plane[3] + need + m}` : 'shell';
      (t as any)._req = { n: nn, a: plane[3] - m, b: plane[3] + need + m, shell: !cutOn };
      const hasMem = this.idx.nTri[t.idx] > 0, hasWall = this.idx.volGlucan[t.idx] > 0 || (this.idx.labelMask[t.idx] & 2) !== 0;
      this.workers[t.idx % this.workers.length].postMessage({
        cmd: 'detail', idx: t.idx, key, ix: this.idx.ix[t.idx], iy: this.idx.iy[t.idx], iz: this.idx.iz[t.idx],
        range, want: { lipids: hasMem && cutOn, wall: hasWall },
      });
    }
  }

  dropDetail(t: LoadedTile) {
    this.lipA.free(t.lipPages); this.freeA.free(t.detPages);
    this.clearSPages(this.r.lipPagesBuf, t.lipPages); this.clearSPages(this.r.freePagesBuf, t.detPages);
    t.lipPages = []; t.detPages = []; t.nLip = 0; t.nDet = 0; t.det = null;
    if (t.underMesh) { t.underMesh = false; this.r.setMeshSlot(t.slot, false); }
    this.syncHigh();
  }

  clearSPages(buf: GPUBuffer, pages: number[]) {
    const z = new Uint32Array(4);
    for (const p of pages) this.r.device.queue.writeBuffer(buf, p * 16, z);
  }
  syncHigh() {
    this.r.usedInstPages = this.instA.high; this.r.usedLipPages = this.lipA.high; this.r.usedFreePages = this.freeA.high;
  }

  // ---------------------------------------------------------------- eviction
  evictOne(protect: number): boolean {
    let best: LoadedTile | null = null;
    for (const t of this.loaded.values()) {
      if (t.idx === protect || this.desiredSet.has(t.idx)) continue;
      if (!best || t.lastUsed < best.lastUsed) best = t;
    }
    if (!best) return false;
    this.evict(best); return true;
  }
  evict(t: LoadedTile) {
    const q = this.r.device.queue;
    const z = new Float32Array(4);
    for (const p of t.pages) { q.writeBuffer(this.r.pageOrg, p * 16, z); this.pageTile[p] = -1; }
    this.instA.free(t.pages);
    this.freeA.free(t.basePages); this.clearSPages(this.r.freePagesBuf, t.basePages);
    this.lipA.free(t.lipPages); this.clearSPages(this.r.lipPagesBuf, t.lipPages);
    this.freeA.free(t.detPages); this.clearSPages(this.r.freePagesBuf, t.detPages);
    t.vb?.destroy(); t.ib?.destroy();
    this.slots.push(t.slot);
    this.loaded.delete(t.idx); this.state[t.idx] = 0;
    this.tileState[tkey(this.idx.ix[t.idx], this.idx.iy[t.idx], this.idx.iz[t.idx])] = 0;
    q.writeBuffer(this.r.tileState, 0, this.tileState);
    this.syncHigh();
  }
  allocWithEvict(a: PageAlloc, k: number, protect: number): number[] | null {
    let p = a.alloc(k);
    while (!p && this.evictOne(protect)) p = a.alloc(k);
    return p;
  }

  // ---------------------------------------------------------------- worker results
  onMessage(m: any) {
    if (m.cmd === 'ready') return;
    if (m.cmd === 'error') {
      console.warn('worker error', m.req, m.idx, m.error);
      this.lastLoadErr = String(m.error).split('\n')[0];
      if (m.req === 'load') { this.state[m.idx] = 0; this.inflight--; }
      if (m.req === 'detail') { this.detInflight--; const t = this.loaded.get(m.idx); if (t) t.detPending = false; }
      return;
    }
    if (m.cmd === 'loaded') { this.inflight--; this.onLoaded(m); }
    if (m.cmd === 'detail') { this.detInflight--; this.onDetail(m); }
  }

  onLoaded(m: any) {
    const i = m.idx as number;
    if (!this.desiredSet.has(i)) { this.state[i] = 0; this.onChange(); return; }   // no longer needed
    const inst: Uint32Array = m.inst, nInst = inst.length / 4;
    const q = this.r.device.queue;
    if (!this.slots.length) this.evictOne(i);
    const pages = this.allocWithEvict(this.instA, Math.ceil(nInst / PAGE), i);
    if (!pages || !this.slots.length) { this.state[i] = 0; this.starved = true; if (pages) this.instA.free(pages); return; }
    const slot = this.slots.pop()!;
    const lo: V3 = m.lo;
    const org = new ArrayBuffer(16), of = new Float32Array(org), ou = new Uint32Array(org);
    of[0] = lo[0] - 80; of[1] = lo[1] - 80; of[2] = lo[2] - 80;
    pages.forEach((p, j) => {
      const n = Math.min(PAGE, nInst - j * PAGE);
      q.writeBuffer(this.r.inst, p * PAGE * 16, inst.buffer, j * PAGE * 16, n * 16);
      ou[3] = n; q.writeBuffer(this.r.pageOrg, p * 16, org);
      this.pageTile[p] = i; this.pageLocal[p] = j * PAGE;
    });
    // mRNA + linker DNA free spheres
    const free: Float32Array = m.free, nFree = free.length / 4;
    let basePages: number[] = [], baseFreeMissing = false;
    if (nFree) {
      const bp = this.allocWithEvict(this.freeA, Math.ceil(nFree / SPAGE), i);
      if (bp) { basePages = bp; this.writeSpheres(this.r.freePool, this.r.freePagesBuf, bp, free, slot); }
      else baseFreeMissing = true;
    }
    let vb: GPUBuffer | null = null, ib: GPUBuffer | null = null;
    const nIdx = (m.ib as Uint32Array).length;
    if (nIdx) {
      vb = this.r.device.createBuffer({ size: m.vb.byteLength, usage: GPUBufferUsage.VERTEX | GPUBufferUsage.COPY_DST });
      q.writeBuffer(vb, 0, m.vb);
      ib = this.r.device.createBuffer({ size: Math.ceil(m.ib.byteLength / 4) * 4, usage: GPUBufferUsage.INDEX | GPUBufferUsage.COPY_DST });
      q.writeBuffer(ib, 0, m.ib);
      this.r.setMeshSlot(slot, false);
    }
    const t: LoadedTile = { idx: i, slot, lo, inst, nInst, pages, basePages, nBase: nFree, lipPages: [], detPages: [], nLip: 0, nDet: 0,
      vb, ib, nIdx, lastUsed: this.frameNo, det: null, detPending: false, underMesh: false, baseFreeMissing };
    this.loaded.set(i, t); this.state[i] = 2;
    this.tileState[tkey(this.idx.ix[i], this.idx.iy[i], this.idx.iz[i])] = 1;
    q.writeBuffer(this.r.tileState, 0, this.tileState);
    this.bytesLoaded += this.idx.gz[i]; this.tilesLoadedTotal++;
    this.syncHigh();
    for (const f of this.onTileLoaded.slice()) f(t);
    this.onChange();
  }

  writeSpheres(pool: GPUBuffer, pagesBuf: GPUBuffer, pages: number[], data: Float32Array, slot: number) {
    const q = this.r.device.queue, n = data.length / 4;
    pages.forEach((p, j) => {
      const c = Math.min(SPAGE, n - j * SPAGE);
      q.writeBuffer(pool, p * SPAGE * 16, data.buffer, data.byteOffset + j * SPAGE * 16, c * 16);
      q.writeBuffer(pagesBuf, p * 16, new Uint32Array([c, slot, 1, 0]));
    });
  }

  onDetail(m: any) {
    const t = this.loaded.get(m.idx);
    if (!t) return;
    t.detPending = false;
    const req = (t as any)._req;
    // replace previous detail
    this.lipA.free(t.lipPages); this.clearSPages(this.r.lipPagesBuf, t.lipPages);
    this.freeA.free(t.detPages); this.clearSPages(this.r.freePagesBuf, t.detPages);
    t.lipPages = []; t.detPages = [];
    const L: Float32Array = m.lipids, F: Float32Array = m.free;
    const nL = L.length / 4, nF = F.length / 4;
    if (nL) {
      const p = this.lipA.alloc(Math.ceil(nL / SPAGE)) ?? this.freeDetailSpace(this.lipA, Math.ceil(nL / SPAGE), t.idx);
      if (p) { t.lipPages = p; this.writeSpheres(this.r.lipPool, this.r.lipPagesBuf, p, L, t.slot); }
    }
    if (nF) {
      const p = this.freeA.alloc(Math.ceil(nF / SPAGE)) ?? this.freeDetailSpace(this.freeA, Math.ceil(nF / SPAGE), t.idx);
      if (p) { t.detPages = p; this.writeSpheres(this.r.freePool, this.r.freePagesBuf, p, F, t.slot); }
    }
    t.nLip = t.lipPages.length ? nL : 0; t.nDet = t.detPages.length ? nF : 0;
    t.det = { n: req.n, a: req.a, b: req.b, shell: req.shell };
    const under = t.nLip > 0;
    if (under !== t.underMesh && t.nIdx) { t.underMesh = under; this.r.setMeshSlot(t.slot, under); }
    this.detailStats.lipids += m.stats.lipCand;
    this.syncHigh();
    this.onChange();
  }
  // free detail of the farthest-from-use tiles until k pages are available
  freeDetailSpace(a: PageAlloc, k: number, protect: number): number[] | null {
    const cands = [...this.loaded.values()].filter(x => x.idx !== protect && (x.lipPages.length || x.detPages.length) && !this.desiredSet.has(x.idx));
    for (const c of cands) { this.dropDetail(c); const p = a.alloc(k); if (p) return p; }
    return null;
  }

  meshDraws(): MeshDraw[] {
    const out: MeshDraw[] = [];
    for (const t of this.loaded.values()) if (t.nIdx && t.vb && t.ib) out.push({ slot: t.slot, vb: t.vb, ib: t.ib, nIdx: t.nIdx, under: t.underMesh });
    return out;
  }

  // ---------------------------------------------------------------- lookups
  instanceAt(poolIdx: number) {
    const page = Math.floor(poolIdx / PAGE);
    const ti = this.pageTile[page];
    if (ti < 0) return null;
    const t = this.loaded.get(ti);
    if (!t) return null;
    const local = this.pageLocal[page] + (poolIdx % PAGE);
    return { t, local, ...this.decode(t, local) };
  }
  decode(t: LoadedTile, local: number) {
    const w = t.inst, w0 = w[4 * local], w1 = w[4 * local + 1], w2 = w[4 * local + 2];
    const QS = 360 / 65535;
    const p: V3 = [t.lo[0] - 80 + (w1 & 0xFFFF) * QS, t.lo[1] - 80 + (w1 >>> 16) * QS, t.lo[2] - 80 + (w2 & 0xFFFF) * QS];
    return { type: w0 & 0x3FFF, region: (w0 >>> 16) & 31, flags: w0 >>> 21, p };
  }
  poolIndex(t: LoadedTile, local: number) { return t.pages[Math.floor(local / PAGE)] * PAGE + (local % PAGE); }
  tileBySlot(slot: number) { for (const t of this.loaded.values()) if (t.slot === slot) return t; return null; }
  countLoaded(type: number) {
    let n = 0;
    for (const t of this.loaded.values()) { const w = t.inst; for (let k = 0; k < t.nInst; k++) if ((w[4 * k] & 0x3FFF) === type) n++; }
    return n;
  }
  nearestInTile(t: LoadedTile, type: number, to: V3) {
    let best = -1, bd = 1e18;
    for (let k = 0; k < t.nInst; k++) {
      if ((t.inst[4 * k] & 0x3FFF) !== type) continue;
      const d = this.decode(t, k), dd = (d.p[0] - to[0]) ** 2 + (d.p[1] - to[1]) ** 2 + (d.p[2] - to[2]) ** 2;
      if (dd < bd) { bd = dd; best = k; }
    }
    return best;
  }
  totals() {
    let inst = 0, lip = 0, free = 0, mesh = 0;
    for (const t of this.loaded.values()) { inst += t.nInst; lip += t.nLip; free += t.nBase + t.nDet; mesh += t.nIdx / 3; }
    return { tiles: this.loaded.size, inst, lip, free, mesh, inflight: this.inflight + this.detInflight, desired: this.desired.length };
  }
  pendingForView() {
    let n = 0; for (const i of this.desired) if (this.state[i] !== 2) n++;
    let d = 0; for (const t of this.loaded.values()) if (t.detPending) d++;
    return n + d;
  }
}
