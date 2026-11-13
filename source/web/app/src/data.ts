// Data access (fetch, or an embedded map for the single-file preview) and library preparation.
import { packColor } from './math';

declare global { interface Window { __YC_EMBED__?: Record<string, string>; } }

export const DATA_BASE = new URL('data/', document.baseURI).href;
// data/version.json changes on every deploy; it is read uncached and appended to every data URL so a browser never mixes
// cached files from two builds (type indices shift between builds)
export const DATA_VER: Promise<string> = fetch(DATA_BASE + 'version.json', { cache: 'no-store' })
  .then(r => (r.ok ? r.json() : { v: '' })).then(j => String(j.v ?? '')).catch(() => '');
export const verQuery = (v: string) => (v ? `?v=${encodeURIComponent(v)}` : '');

export function b64ToBytes(b64: string): Uint8Array {
  const bin = atob(b64); const out = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
  return out;
}

export async function gunzipIfNeeded(buf: ArrayBuffer): Promise<ArrayBuffer> {
  const u8 = new Uint8Array(buf, 0, Math.min(2, buf.byteLength));
  if (u8[0] === 0x1f && u8[1] === 0x8b) {
    const ds = new DecompressionStream('gzip');
    const stream = new Blob([buf]).stream().pipeThrough(ds);
    return await new Response(stream).arrayBuffer();
  }
  return buf;
}

export async function fetchBin(path: string, onProgress?: (n: number, total: number) => void): Promise<ArrayBuffer> {
  const emb = window.__YC_EMBED__;
  if (emb && emb[path] !== undefined) return gunzipIfNeeded(b64ToBytes(emb[path]).buffer as ArrayBuffer);
  const r = await fetch(DATA_BASE + path + verQuery(await DATA_VER));
  if (!r.ok) throw new Error(`HTTP ${r.status} for ${path}`);
  if (!onProgress || !r.body) return gunzipIfNeeded(await r.arrayBuffer());
  const total = Number(r.headers.get('content-length') || 0);
  const reader = r.body.getReader(); const chunks: Uint8Array[] = []; let n = 0;
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    chunks.push(value); n += value.length; onProgress(n, total);
  }
  const out = new Uint8Array(n); let o = 0;
  for (const c of chunks) { out.set(c, o); o += c.length; }
  return gunzipIfNeeded(out.buffer);
}
export async function fetchJson(path: string): Promise<any> {
  const b = await fetchBin(path);
  return JSON.parse(new TextDecoder().decode(b));
}

export interface TypeRec {
  id: string; gene: string; genes: string; uniprot: string; name: string; src: string; pdb: string;
  oligomer: number; nb: number; nnt: number; bo: number; so: number; sn: number; br: number; rg: number;
  mass: number; plddt_mean: number | null; plddt_frac_kept: number | null; copies: number; regions: Record<string, number>;
  model?: string;   // AlphaFold DB entry of a predicted complex (src 'AF-multimer')
}

export const REGION_NAMES = ['cytosol', 'nucleoplasm', 'nucleolus', 'ER lumen', 'vacuole lumen', 'mitochondrial matrix',
  'mitochondrial IMS', 'Golgi / endosome', 'peroxisome', 'cell wall', 'PM peripheral', 'NE peripheral', 'ER peripheral',
  'lipid droplet surface', 'plasma membrane', 'ER membrane', 'nuclear envelope', 'vacuolar membrane',
  'mito outer membrane', 'mito inner membrane', 'Golgi membrane', 'peroxisomal membrane', 'chromatin (nucleosome)', 'nuclear pore',
  'vacuole surface (cytosolic face)', 'Golgi surface (cytosolic face)', 'peroxisome surface (cytosolic face)'];
export const MEM_NAMES = ['plasma membrane', 'ER membrane', 'nuclear envelope', 'vacuolar membrane', 'mito outer membrane',
  'mito inner membrane', 'Golgi/endosome membrane', 'peroxisomal membrane', 'lipid droplet monolayer', 'cell surface'];
export const LABEL_NAMES = ['exterior', 'cell wall (mannoprotein layer)', 'cell wall (glucan/chitin layer)', 'cytosol', 'cortical ER lumen',
  'perinuclear space', 'nucleoplasm', 'nucleolus', 'vacuole', 'mito intermembrane space', 'mito matrix', 'lipid droplet',
  'Golgi / endosome', 'peroxisome', 'nuclear pore channel'];

// Goodsell-inspired palette (shared by shaders through a storage buffer)
export const PAL: Record<number, string> = {
  0: '#7EA7D3', 1: '#B993C9', 2: '#D39BC0', 3: '#A3C48C', 4: '#9FC7C0', 5: '#D9946B', 6: '#E8B48F', 7: '#B7C77E',
  8: '#C5B0E0', 9: '#9A8352', 10: '#7EA7D3', 11: '#8FA3D0', 12: '#7EA7D3', 13: '#7EA7D3', 14: '#5FAE75', 15: '#7DBE86',
  16: '#6BB596', 17: '#5BA8A0', 18: '#C27C55', 19: '#B5653F', 20: '#93B060', 21: '#A68BC9', 22: '#F4B26A', 23: '#3E9C9A',
  24: '#7EA7D3', 25: '#7EA7D3', 26: '#7EA7D3',
  32: '#F3C643', 33: '#F0BC3C', 34: '#F6D15E', 35: '#EAB542', 36: '#EFA83A', 37: '#E7A93C', 38: '#F2CA55', 39: '#E8873A',
  40: '#F3C643', 41: '#F8E3A0', 42: '#E7A93C',
  48: '#B0437A', 49: '#E07636', 50: '#CDA866', 51: '#C3C27A', 52: '#C9558F', 53: '#7CC6C2', 54: '#E2D6A8', 55: '#9B8BC4',
  56: '#5FAE75', 57: '#7DBE86', 58: '#6BB596', 59: '#5BA8A0', 60: '#C27C55', 61: '#B5653F', 62: '#93B060', 63: '#A68BC9', 64: '#E9C46A', 65: '#9A8352',
  // membranes coloured by the compartment they enclose, as in the section map (PM, ER, NE, vacuole, mito OM, IM,
  // Golgi, peroxisome, lipid droplet, cell surface); the PM is a deeper tone of the cytosol's
  66: '#A9C3E3', 67: '#A3C48C', 68: '#B993C9', 69: '#9FC7C0', 70: '#E8B48F', 71: '#D9946B', 72: '#B7C77E', 73: '#C5B0E0', 74: '#E9C46A', 75: '#9A8352',
  80: '#F6F1E7', 81: '#F1E9D3', 82: '#EEE3CC', 83: '#E4ECF5', 84: '#E8F0DF', 85: '#E8F0DF', 86: '#F1E7F1', 87: '#EBD9EA',
  88: '#E3F0EE', 89: '#F6E6DA', 90: '#F2DCCB', 91: '#F7EDC4', 92: '#EEF1DC', 93: '#EEE8F6', 94: '#ECE8F2',
};
/** bud scars (cell3d T_TORUS prims, class 9): per scar centre xyz, major radius, axis xyz, tube radius; at most 4 */
export function scarList(geo: any): number[] {
  const out: number[] = [];
  for (let i = 0; i < geo.PT.length && out.length < 32; i++) {
    if (geo.PC[i] !== 9) continue;
    const q = geo.PQ.slice(8 * i, 8 * i + 8);
    out.push(q[0], q[1], q[2], q[6], q[3], q[4], q[5], q[7]);
  }
  return out;
}

export function paletteArray(): Float32Array {
  const a = new Float32Array(96 * 4);
  for (const [k, h] of Object.entries(PAL)) {
    const v = parseInt(h.slice(1), 16), i = Number(k) * 4;
    a[i] = ((v >> 16) & 255) / 255; a[i + 1] = ((v >> 8) & 255) / 255; a[i + 2] = (v & 255) / 255; a[i + 3] = 1;
  }
  return a;
}

export const SPECIAL: Record<string, { c0: string; c1?: string; cls: string }> = {
  ribosome_80S: { c0: '#EE9BC0', c1: '#C9558F', cls: 'ribosome' },
  tRNA: { c0: '#E27FA8', c1: '#E27FA8', cls: 'tRNA' },
  nucleosome: { c0: '#F4B26A', c1: '#E07636', cls: 'nucleosome' },
  NPC: { c0: '#3E9C9A', cls: 'NPC' },
  PolI: { c0: '#9A72B8', c1: '#9A72B8', cls: 'polymerase' },
  PolII: { c0: '#9A72B8', c1: '#9A72B8', cls: 'polymerase' },
  PolIII: { c0: '#9A72B8', c1: '#9A72B8', cls: 'polymerase' },
  FAS: { c0: '#8B88C9', cls: 'complex' },
  proteasome_26S: { c0: '#8B88C9', cls: 'complex' },
  CCT: { c0: '#8B88C9', cls: 'complex' },
  F_actin: { c0: '#E0893A', c1: '#E0893A', cls: 'actin' },
  Arp23: { c0: '#D9534F', c1: '#D9534F', cls: 'actin' },
  TPM1_coil: { c0: '#2FA4A0', c1: '#2FA4A0', cls: 'actin' },
  CP_cap: { c0: '#8E6CB0', c1: '#8E6CB0', cls: 'actin' },
  MT_tubulin: { c0: '#74C476', c1: '#74C476', cls: 'microtubule' },
  Pil1_dimer: { c0: '#3F8FCB', c1: '#3F8FCB', cls: 'eisosome' },
  // Cvt cargo: aminopeptidase I dodecamers, alpha-mannosidase, receptor
  APE1_oligo: { c0: '#E0A43A', cls: 'complex' },
  AMS1_oligo: { c0: '#B9792E', cls: 'complex' },
  ATG19: { c0: '#C0473F', cls: 'complex' },
  Vac8_Nvj1: { c0: '#2E9C8F', cls: 'complex' },
  RVS161_RVS167_het: { c0: '#9B4DCA', cls: 'complex' },
  // endocytic coat (also free elsewhere)
  SLA2: { c0: '#E3B23C', cls: 'complex' },
  SLA1: { c0: '#CFA040', cls: 'complex' },
  END3: { c0: '#BB8A35', cls: 'complex' },
  Lsp1_dimer: { c0: '#7FB6E0', c1: '#7FB6E0', cls: 'eisosome' },
  gTuRC_nuc: { c0: '#2E8B57', c1: '#2E8B57', cls: 'microtubule' },
  gTuRC_cyto: { c0: '#2E8B57', c1: '#2E8B57', cls: 'microtubule' },
  ribosome_80S_Ssh1: { c0: '#EE9BC0', c1: '#C9558F', cls: 'ribosome' },
  ribosome_80S_RAC: { c0: '#EE9BC0', c1: '#C9558F', cls: 'ribosome' },
  ribosome_80S_Sec61: { c0: '#EE9BC0', c1: '#C9558F', cls: 'ribosome' },
  // ribosome assembly intermediates: assembly factors and r-proteins in coral, pre-rRNA as rRNA
  preribo_90S: { c0: '#EF9A7E', c1: '#C9558F', cls: 'pre-ribosome' },
  preribo_60S_nucleolar: { c0: '#EF9A7E', c1: '#C9558F', cls: 'pre-ribosome' },
  preribo_60S_nucleoplasm: { c0: '#EF9A7E', c1: '#C9558F', cls: 'pre-ribosome' },
  preribo_60S_cytoplasm: { c0: '#EF9A7E', c1: '#C9558F', cls: 'pre-ribosome' },
  preribo_40S: { c0: '#EF9A7E', c1: '#C9558F', cls: 'pre-ribosome' },
};

export interface Library {
  types: TypeRec[];
  beads: Float32Array;         // vec4 per bead: xyz, w = bits(kind | plddt << 8)
  st: Float32Array;            // sphere tree xyzr
  clusters: ArrayBuffer;       // per cluster: vec4f centre/radius + vec4u (offset, count, 0, 0)
  typeInfo: ArrayBuffer;       // 16 x u32 per type
  netOn: boolean;              // network profiles in this build (library/organisation.json)
  netNb: Record<string, [number, number][]>;   // per type index: strongest network neighbours [type, r]
  nClusters: number;
  byId: Map<string, number>;
}

export async function loadLibrary(progress: (msg: string, f: number) => void): Promise<Library> {
  progress('Loading structure index', 0.02);
  const tj = await fetchJson('library/types.json');
  const types: TypeRec[] = tj.types;
  const q = tj.bead_q as number;
  progress('Loading residue coordinates (2.5 M beads)', 0.06);
  const raw = await fetchBin('library/beads.bin', (n, t) => progress('Loading residue coordinates (2.5 M beads)', 0.06 + 0.5 * (t ? n / t : 0)));
  const st = new Float32Array(await fetchBin('library/spheretree.bin'));
  // optional: network-neighbourhood colour and neighbours of each type (organisation layer of the packing)
  let netCol: Record<string, string> = {}, netNb: Record<string, [number, number][]> = {}, netOn = false;
  try {
    const oj = await fetchJson('library/organisation.json');
    netCol = oj.colors || {}; netNb = oj.neighbours || {}; netOn = !!oj.colors;
  } catch { /* no organisation layer in this build */ }
  progress('Building level-of-detail clusters', 0.6);
  const nB = raw.byteLength / 8;
  const i16 = new Int16Array(raw), u8 = new Uint8Array(raw);
  const beads = new Float32Array(nB * 4), bu = new Uint32Array(beads.buffer);
  for (let i = 0; i < nB; i++) {
    beads[4 * i] = i16[4 * i] * q; beads[4 * i + 1] = i16[4 * i + 1] * q; beads[4 * i + 2] = i16[4 * i + 2] * q;
    bu[4 * i + 3] = u8[8 * i + 6] | (u8[8 * i + 7] << 8);
  }
  // 64-bead clusters (beads are Morton-ordered within each type)
  let nC = 0;
  for (const t of types) nC += Math.ceil(t.nb / 64);
  const clusters = new ArrayBuffer(nC * 32), cf = new Float32Array(clusters), cu = new Uint32Array(clusters);
  const typeInfo = new ArrayBuffer(types.length * 64), tu = new Uint32Array(typeInfo), tf = new Float32Array(typeInfo);
  let c = 0;
  types.forEach((t, ti) => {
    const c0 = c;
    for (let o = 0; o < t.nb; o += 64) {
      const n = Math.min(64, t.nb - o), b0 = t.bo + o;
      let mnx = 1e9, mny = 1e9, mnz = 1e9, mxx = -1e9, mxy = -1e9, mxz = -1e9;
      for (let k = 0; k < n; k++) {
        const j = 4 * (b0 + k), x = beads[j], y = beads[j + 1], z = beads[j + 2];
        if (x < mnx) mnx = x; if (x > mxx) mxx = x; if (y < mny) mny = y; if (y > mxy) mxy = y; if (z < mnz) mnz = z; if (z > mxz) mxz = z;
      }
      const cx = (mnx + mxx) / 2, cy = (mny + mxy) / 2, cz = (mnz + mxz) / 2;
      let r = 0;
      for (let k = 0; k < n; k++) {
        const j = 4 * (b0 + k);
        r = Math.max(r, Math.hypot(beads[j] - cx, beads[j + 1] - cy, beads[j + 2] - cz));
      }
      cf[8 * c] = cx; cf[8 * c + 1] = cy; cf[8 * c + 2] = cz; cf[8 * c + 3] = r + 0.85;
      cu[8 * c + 4] = b0; cu[8 * c + 5] = n;
      c++;
    }
    const sp = SPECIAL[t.id] ?? (t.id.startsWith('tRNA_') ? SPECIAL.tRNA : undefined);
    const isPDB = t.src === 'PDB';
    tu[16 * ti + 0] = t.bo; tu[16 * ti + 1] = t.nb; tu[16 * ti + 2] = c0; tu[16 * ti + 3] = c - c0;
    tu[16 * ti + 4] = t.so; tu[16 * ti + 5] = t.sn;
    tf[16 * ti + 6] = t.br; tf[16 * ti + 7] = Math.sqrt(5 / 3) * t.rg;
    tu[16 * ti + 8] = sp ? packColor(sp.c0) : 0; tu[16 * ti + 9] = sp && sp.c1 ? packColor(sp.c1) : 0;
    tu[16 * ti + 10] = ((isPDB ? 1 : 0) | (Math.min(t.copies, 0xFFFFFF) << 8)) >>> 0;
    tf[16 * ti + 11] = isPDB ? 255 : (t.plddt_mean ?? 255);
    tu[16 * ti + 15] = netCol[ti] ? packColor(netCol[ti]) : 0;
  });
  const byId = new Map<string, number>();
  types.forEach((t, i) => byId.set(t.id, i));
  return { types, beads, st, clusters, typeInfo, netOn, netNb, nClusters: nC, byId };
}

export interface TileIndex {
  n: number; ix: Int32Array; iy: Int32Array; iz: Int32Array; nInst: Uint32Array; nGhost: Uint32Array;
  nMr: Uint32Array; nLk: Uint32Array; nVert: Uint32Array; nTri: Uint32Array; gz: Uint32Array;
  labelMask: Uint32Array; volGlucan: Float32Array; volMannan: Float32Array; key: Map<number, number>;
  wall: any; regionNames: string[];
}
export const TGRID = 26;     // tiles per axis of the tile-state grid (-13 .. 12: +-2.6 um, encloses the real cell outline)
export const tkey = (ix: number, iy: number, iz: number) => ((ix + 13) * TGRID + (iy + 13)) * TGRID + (iz + 13);

export async function loadTileIndex(): Promise<TileIndex> {
  const j = await fetchJson('tiles/index.json');
  const cols: string[] = j.tile_columns; const T: number[][] = j.tiles; const n = T.length;
  const ci = (k: string) => cols.indexOf(k);
  const mk = <A extends Int32Array | Uint32Array | Float32Array>(C: new (n: number) => A, k: string) => {
    const a = new C(n), c = ci(k); for (let i = 0; i < n; i++) a[i] = T[i][c]; return a;
  };
  const idx: TileIndex = {
    n, ix: mk(Int32Array, 'ix'), iy: mk(Int32Array, 'iy'), iz: mk(Int32Array, 'iz'),
    nInst: mk(Uint32Array, 'n_inst'), nGhost: mk(Uint32Array, 'n_ghost'), nMr: mk(Uint32Array, 'n_mrna_beads'),
    nLk: mk(Uint32Array, 'n_link_beads'), nVert: mk(Uint32Array, 'n_vert'), nTri: mk(Uint32Array, 'n_tri'),
    gz: mk(Uint32Array, 'gz_bytes'), labelMask: mk(Uint32Array, 'label_mask'),
    volGlucan: mk(Float32Array, 'vol_glucan_nm3'), volMannan: mk(Float32Array, 'vol_mannan_nm3'),
    key: new Map(), wall: j.wall, regionNames: j.region_names,
  };
  for (let i = 0; i < n; i++) idx.key.set(tkey(idx.ix[i], idx.iy[i], idx.iz[i]), i);
  return idx;
}
