// MD-ready export of a simulation box: coarse-grained mmCIF (one pseudo-atom per residue bead, one chain per
// molecule, membrane beads as a fixed chain), a JSON manifest (molecule sources, poses, diffusion coefficients,
// settings), the recorded rigid-body trajectory, and an OpenMM script (sim/run_openmm.py) that runs the system as an
// elastic-network coarse-grained model.
import type { Library } from '../data';
import type { SimEngine } from './engine';
import type { SimMeta } from './extract';
import type { SimSettings } from './params';
import { T_K, ETA_WATER, HYDRATION, hiTrans, hiRot, NT_CHARGE, BJERRUM_NM, debyeKappa } from './params';
import { loadSimData } from './extract';
import RUN_OPENMM from './run_openmm.py?raw';
import TRAJ_TO_DCD from './trajectory_to_dcd.py?raw';

function chainId(i: number): string {
  // base-36 chain names, 1-4 characters (mmCIF allows longer; 4 keeps PDB-like tools happy)
  return i.toString(36).toUpperCase().padStart(2, '0');
}

function save(name: string, data: string | Blob) {
  const blob = typeof data === 'string' ? new Blob([data], { type: 'text/plain' }) : data;
  const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = name;
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 4000);
}

function qrot(q: ArrayLike<number>, o: number, v: [number, number, number]): [number, number, number] {
  const x = q[o], y = q[o + 1], z = q[o + 2], w = q[o + 3];
  const tx = 2 * (y * v[2] - z * v[1]), ty = 2 * (z * v[0] - x * v[2]), tz = 2 * (x * v[1] - y * v[0]);
  return [v[0] + w * tx + (y * tz - z * ty), v[1] + w * ty + (z * tx - x * tz), v[2] + w * tz + (x * ty - y * tx)];
}

/** membrane surface points of the box from the engine's wall grid (fixed beads in the export) */
async function membranePoints(eng: SimEngine): Promise<Float32Array> {
  const d = eng.d, st = eng.uniState, n = st.wd[0] * st.wd[1] * st.wd[2];
  const rb = d.createBuffer({ size: 16 * n, usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST });
  const enc = d.createCommandEncoder(); enc.copyBufferToBuffer(eng.bufs.wallGrid, 0, rb, 0, 16 * n); d.queue.submit([enc.finish()]);
  await rb.mapAsync(GPUMapMode.READ); const w = new Float32Array(rb.getMappedRange().slice(0)); rb.unmap(); rb.destroy();
  const out: number[] = [], v = st.vox, lo = eng.sys.lo;
  for (let ix = 0; ix < st.wd[0]; ix++) for (let iy = 0; iy < st.wd[1]; iy++) for (let iz = 0; iz < st.wd[2]; iz++) {
    const i = (ix * st.wd[1] + iy) * st.wd[2] + iz, dd = w[4 * i];
    if (Math.abs(dd) > 0.5 * v) continue;
    const g = [w[4 * i + 1], w[4 * i + 2], w[4 * i + 3]], gl = Math.hypot(g[0], g[1], g[2]) || 1;
    const p = [lo[0] + (ix + 0.5) * v, lo[1] + (iy + 0.5) * v, lo[2] + (iz + 0.5) * v];
    out.push(p[0] - dd * g[0] / gl, p[1] - dd * g[1] / gl, p[2] - dd * g[2] / gl);
  }
  return new Float32Array(out);
}

export async function exportFiles(eng: SimEngine, meta: SimMeta, lib: Library, set: SimSettings,
                                  pos: Float32Array, quat: Float32Array, frames: { tNs: number; pos: Float32Array; quat: Float32Array }[]) {
  const sys = eng.sys, n = sys.n;
  const lines: string[] = [
    'data_yeast_cell_box', '#',
    `_cell.length_a ${(10 * (sys.hi[0] - sys.lo[0])).toFixed(3)}`, `_cell.length_b ${(10 * (sys.hi[1] - sys.lo[1])).toFixed(3)}`,
    `_cell.length_c ${(10 * (sys.hi[2] - sys.lo[2])).toFixed(3)}`, '_cell.angle_alpha 90', '_cell.angle_beta 90', '_cell.angle_gamma 90', '#',
    'loop_', '_atom_site.group_PDB', '_atom_site.id', '_atom_site.type_symbol', '_atom_site.label_atom_id', '_atom_site.label_comp_id',
    '_atom_site.label_asym_id', '_atom_site.label_entity_id', '_atom_site.label_seq_id', '_atom_site.Cartn_x', '_atom_site.Cartn_y',
    '_atom_site.Cartn_z', '_atom_site.occupancy', '_atom_site.B_iso_or_equiv', '_atom_site.auth_seq_id', '_atom_site.auth_asym_id',
    '_atom_site.pdbx_PDB_model_num',
  ];
  let id = 1;
  const B = sys.beadLocal, K = sys.beadKind;
  for (let i = 0; i < n; i++) {
    const ch = chainId(i), mob = sys.mobile[i] ? '1.00' : '0.00';
    for (let k = 0; k < sys.beadCnt[i]; k++) {
      const b = sys.beadOff[i] + k, nuc = (K[b] & 255) === 1;
      const v = qrot(quat, 4 * i, [B[4 * b], B[4 * b + 1], B[4 * b + 2]]);
      const x = 10 * (pos[3 * i] + v[0]), y = 10 * (pos[3 * i + 1] + v[1]), z = 10 * (pos[3 * i + 2] + v[2]);
      lines.push(`ATOM ${id++} ${nuc ? 'P' : 'C'} ${nuc ? 'P' : 'CA'} ${nuc ? 'N' : 'UNK'} ${ch} 1 ${k + 1} ${x.toFixed(2)} ${y.toFixed(2)} ${z.toFixed(2)} ${mob} 0 ${k + 1} ${ch} 1`);
    }
  }
  // per-bead interaction data: charges (nucleotides condensed), mRNA links between the beads nearest the channel
  // anchors, native contacts of predicted complexes between every pair of partner molecules (residue level)
  const data = await loadSimData();
  const atomBase = new Uint32Array(n); { let a = 0; for (let i = 0; i < n; i++) { atomBase[i] = a; a += sys.beadCnt[i]; } }
  const charges: [number, number][] = [];
  if (data.charges) for (let i = 0; i < n; i++) {
    const t = lib.types[sys.type[i]];
    for (let k = 0; k < sys.beadCnt[i]; k++) {
      const z = data.charges[t.bo + k]; if (!z) continue;
      charges.push([atomBase[i] + k, (K[sys.beadOff[i] + k] & 255) === 1 ? -NT_CHARGE * z : z]);
    }
  }
  const nearestBead = (i: number, x: number, y: number, z: number) => {
    let best = 0, bd = Infinity;
    for (let k = 0; k < sys.beadCnt[i]; k++) {
      const b = 4 * (sys.beadOff[i] + k), d = (B[b] - x) ** 2 + (B[b + 1] - y) ** 2 + (B[b + 2] - z) ** 2;
      if (d < bd) { bd = d; best = k; }
    }
    return atomBase[i] + best;
  };
  const links: any[] = [];
  for (let i = 0; i < n; i++) for (let k = sys.bondOff[i]; k < sys.bondOff[i] + sys.bondCnt[i]; k++) {
    const r = sys.bonds.subarray(12 * k, 12 * k + 12);
    if (r[11] !== 1) continue;        // each link once
    const j = r[8];
    links.push({ chains: [chainId(i), chainId(j)], atoms: [nearestBead(i, r[0], r[1], r[2]), nearestBead(j, r[4], r[5], r[6])],
      anchors_local_nm: [[r[0], r[1], r[2]].map(v => +v.toFixed(3)), [r[4], r[5], r[6]].map(v => +v.toFixed(3))],
      contour_nm: +r[3].toFixed(2), persistence_nm: r[7] });
  }
  const natives: number[][] = [];
  if (data.contacts.size) {
    const byType = new Map<number, number[]>();
    for (let i = 0; i < n; i++) if (sys.mobile[i]) (byType.get(sys.type[i]) ?? byType.set(sys.type[i], []).get(sys.type[i])!).push(i);
    for (const [ta, bodiesA] of byType) for (const g of data.contacts.get(ta) ?? []) {
      if (g.other < ta) continue;                       // each type pair once (a homodimer's list holds both directions)
      const bodiesB = byType.get(g.other); if (!bodiesB) continue;
      const m = g.other === ta ? g.mine.length / 2 : g.mine.length, share = 1 / m;
      for (const i of bodiesA) for (const j of bodiesB) {
        if (i === j || (g.other === ta && j < i)) continue;
        for (let k = 0; k < m && natives.length < 300000; k++) natives.push([atomBase[i] + g.mine[k], atomBase[j] + g.theirs[k], +share.toFixed(5)]);
      }
    }
  }
  const mem = await membranePoints(eng);
  for (let k = 0; k < mem.length / 3; k++)
    lines.push(`ATOM ${id++} C M MEM ZZZZ 2 ${k + 1} ${(10 * mem[3 * k]).toFixed(2)} ${(10 * mem[3 * k + 1]).toFixed(2)} ${(10 * mem[3 * k + 2]).toFixed(2)} 0.00 0 ${k + 1} ZZZZ 1`);
  lines.push('#');
  const T = lib.types;
  const manifest = {
    format: 'yeast-cell box export v1',
    units: { length: 'nm', time: 'ns', diffusion: 'um^2/s (in water at the given temperature; divide by the viscosity factor)' },
    temperature_K: T_K, water_viscosity_Pa_s: ETA_WATER, hydration_shell_nm: HYDRATION,
    box: { lo: sys.lo, hi: sys.hi, boundary: sys.boundary, shell_nm: sys.boundary === 'embedded' ? set.shell : 0,
      tether_kT_per_nm2: sys.boundary === 'embedded' ? set.tether : 0 },
    settings: set,
    volume_fraction_hydrated: +sys.phi.toFixed(4),
    mean_field_hydrodynamics: set.hydro ? { D_t_factor: +hiTrans(sys.phi).toFixed(4), D_r_factor: +hiRot(sys.phi).toFixed(4),
      model: 'short-time self-diffusion of hard spheres at the box volume fraction, Cichocki et al. 1999 J Chem Phys 111:3265' } : null,
    note: 'Coarse-grained export of the yeast whole-cell model. Each molecule is one chain of pseudo-atoms, one per residue (CA) ' +
      'or nucleotide (P) bead of its structure, posed as in the cell (or the selected simulation frame). Chain ZZZZ holds fixed ' +
      'membrane surface beads (bilayer midplane). Mobile molecules have occupancy 1, fixed ones (membrane proteins, filaments, ' +
      'chromatin) 0. box.boundary: periodic (molecule centres wrapped into the box), walls, or embedded (molecules flagged ' +
      'tethered are held near their position by a harmonic tether of box.tether_kT_per_nm2 on their centre). To build an atomistic or Martini system, take each molecule\'s structure from "source" ' +
      'and apply its rotation (quaternion x,y,z,w) and translation (nm) to the structure\'s bead frame (CA/P centred as in the library).',
    bodies: Array.from({ length: n }, (_, i) => {
      const t = T[sys.type[i]], h = meta.hydro[i];
      return {
        chain: chainId(i), type: t.id, name: t.name, genes: t.genes, uniprot: t.uniprot || undefined,
        source: t.src === 'PDB' ? `PDB ${t.pdb}` : t.model ? `AlphaFold DB ${t.model}` : t.uniprot ? `AlphaFold DB AF-${t.uniprot}-F1` : t.src,
        region: meta.region[i], mobile: !!sys.mobile[i], tethered: !!sys.tether[i], beads: sys.beadCnt[i],
        R_h_nm: +h.rh.toFixed(3), D_t_um2_s: +(h.dt * 1e3).toFixed(3), D_r_per_us: +(h.dr * 1e3).toFixed(4),
        pos_nm: Array.from(pos.subarray(3 * i, 3 * i + 3), v => +v.toFixed(3)),
        quat_xyzw: Array.from(quat.subarray(4 * i, 4 * i + 4), v => +v.toFixed(6)),
      };
    }),
    membrane_beads: mem.length / 3,
    interactions: {
      electrostatics: set.elec && !!data.charges ? { model: 'Debye-Hueckel between charged beads (kT units)', bjerrum_nm: BJERRUM_NM,
        ionic_strength_M: set.ionic, kappa_per_nm: +debyeKappa(set.ionic).toFixed(4), nucleotide_charge: NT_CHARGE,
        note: 'charges: [atom index (0-based, CIF order), charge in e]; Lys/Arg +1, Asp/Glu -1 (pH 7.2), nucleotides scaled for counterion condensation' } : null,
      charges: set.elec ? charges : [],
      mrna_links: set.mrna ? links : [],
      mrna_link_model: 'worm-like chain (Marko-Siggia) between the beads nearest the mRNA exit of one ribosome and the entry of the next',
      native_contacts: set.specific ? natives : [],
      native_contact_model: { interface_energy_kT: set.epsGo, note: '[atom i, atom j, share]: residue contacts of a confident AlphaFold DB dimer ' +
        'between two molecules of the partner types; each pair is a well of depth share x interface_energy_kT (cos^2, 0.6 nm wide, beyond contact)' },
    },
  };
  const step = Math.max(1, Math.ceil(frames.length / 100));
  const traj = {
    note: 'Rigid-body trajectory of the browser simulation: per frame, time (ns) and for every molecule (manifest order) ' +
      'position (nm) and quaternion (x,y,z,w).',
    frames: frames.filter((_, k) => k % step === 0).map(f => ({
      t_ns: +f.tNs.toFixed(3),
      pos: Array.from(f.pos, v => +v.toFixed(2)), quat: Array.from(f.quat, v => +v.toFixed(4)),
    })),
  };
  save('yeast_box.cif', lines.join('\n') + '\n');
  save('yeast_box.json', JSON.stringify(manifest));
  save('yeast_box_trajectory.json', JSON.stringify(traj));
  save('run_openmm.py', RUN_OPENMM);
  save('trajectory_to_dcd.py', TRAJ_TO_DCD);
}
