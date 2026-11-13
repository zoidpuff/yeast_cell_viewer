// Physical parameters of the simulation models (MODEL.md section 14).
// Units: nm, ns, kT (T = 303 K).
import type { Library } from '../data';

export const T_K = 303.15;
export const KT = 1.380649e-23 * T_K;          // J
export const ETA_WATER = 0.797e-3;              // Pa s, water at 30 C
export const HYDRATION = 0.3;                   // nm added to the bead-model radius (hydration shell)
export const HALF = 2.0;                        // membrane half thickness (cell3d.HALF)

export interface TypeHydro { rh: number; dt: number; dr: number }   // nm, nm^2/ns, rad^2/ns (in water)

/** Hydrodynamic radius from the residue beads: Kirkwood double sum over bead pairs, 1/R_K = <1/r_ij>, scaled by
 *  6/5 (a filled sphere of beads gives R_K = 5R/6) plus a hydration shell. Shape-aware (elongated molecules get a
 *  larger R_h than spheres of the same mass); for globular proteins it reproduces R_h = 1.29 R_g + hydration.
 *  Translational and rotational diffusion in water from Stokes-Einstein(-Debye) with that radius. */
export function hydro(lib: Library, t: number, cache: Map<number, TypeHydro>): TypeHydro {
  const c = cache.get(t);
  if (c) return c;
  const T = lib.types[t], B = lib.beads, n = T.nb, b0 = T.bo;
  let rk = 0.5;
  if (n >= 2) {
    const m = Math.min(n, 700), step = n / m;
    let s = 0, k = 0;
    for (let i = 0; i < m; i++) {
      const a = 4 * (b0 + Math.floor(i * step));
      for (let j = i + 1; j < m; j++) {
        const b = 4 * (b0 + Math.floor(j * step));
        const r = Math.hypot(B[a] - B[b], B[a + 1] - B[b + 1], B[a + 2] - B[b + 2]);
        s += 1 / Math.max(r, 0.38); k++;
      }
    }
    rk = k ? k / s : 0.5;
  }
  const rh = Math.max(0.5, 1.2 * rk + HYDRATION);
  const dt = KT / (6 * Math.PI * ETA_WATER * rh * 1e-9) * 1e9;            // m^2/s -> nm^2/ns
  const dr = KT / (8 * Math.PI * ETA_WATER * (rh * 1e-9) ** 3) * 1e-9;     // 1/s -> 1/ns
  const out = { rh, dt, dr };
  cache.set(t, out);
  return out;
}

/** Mean-field hydrodynamic interactions, as in the BD of concentrated protein solutions of Mereghetti & Wade (2012)
 *  and Balbo et al. (2013): every molecule's short-time diffusion is scaled by that of hard spheres at the volume
 *  fraction phi of its surroundings. Second-order expansions including lubrication (Cichocki, Ekiel-Jezewska &
 *  Wajnryb 1999, J Chem Phys 111:3265); the first-order translational term is Batchelor's -1.83 phi. */
export const hiTrans = (phi: number) => Math.max(0.05, 1 - 1.832 * phi - 0.219 * phi * phi);
export const hiRot = (phi: number) => Math.max(0.05, 1 - 0.631 * phi - 0.726 * phi * phi);
/** hydrated macromolecular volume, nm^3 per Da: partial specific volume 0.73 mL/g + 0.35 g/g bound water */
export const HYD_VOL_PER_DA = (0.73 + 0.35) / 602.214;

/** Lateral diffusion of a membrane protein (Saffman & Delbrueck 1975, PNAS 72:3111): a cylinder of radius a in a
 *  membrane of viscosity MEM_VISC and thickness MEM_H between water: D = kT/(4 pi mu_m h) (ln(mu_m h/(mu a)) - gamma),
 *  D_r (about the membrane normal) = kT/(4 pi mu_m h a^2). The radius of the membrane-embedded part is taken as half
 *  the hydrodynamic radius, 1-5 nm (rough; D depends on it only logarithmically). Returns nm^2/ns and 1/ns. */
export const MEM_VISC = 0.1;      // Pa s, fluid lipid bilayer (assumption)
export const MEM_H = 4e-9;        // m, bilayer thickness
export function saffmanDelbruck(rh: number, memVisc = MEM_VISC): { dt: number; dr: number } {
  const a = Math.min(5, Math.max(1, 0.5 * rh)) * 1e-9, mh = memVisc * MEM_H;
  const dt = KT / (4 * Math.PI * mh) * (Math.log(mh / (ETA_WATER * a)) - 0.5772) * 1e9;
  const dr = KT / (4 * Math.PI * mh * a * a) * 1e-9;
  return { dt, dr };
}

/** Electrostatics: screened Coulomb (Debye-Hueckel) between charged collision spheres. Bjerrum length of water at
 *  303 K (relative permittivity 76.6): e^2 / (4 pi eps0 eps_r kT) = 0.72 nm. Debye screening kappa = sqrt(8 pi lB N_A I)
 *  = 3.30 sqrt(I / M) per nm. Default ionic strength 0.3 M: the in vivo-like yeast assay medium of van Eunen et al.
 *  (2010, FEBS J 277:749) has 300 mM K+, 245 mM glutamate, 50 mM phosphate, 20 mM Na+, 2 mM free Mg2+ (I about 0.33 M),
 *  so the Debye length is about 0.55 nm. Residue charges at the cytosolic pH of 7.2 (Orij et al. 2009, Microbiology
 *  155:268): Asp, Glu -1; Lys, Arg +1; His neutral. Nucleotides: -0.25 e (counterion condensation on double helices,
 *  Manning 1969 J Chem Phys 51:924: a charge fraction of 1/xi = b / lB remains, b = 0.17 nm per phosphate). */
export const BJERRUM_NM = 0.72;
export const NT_CHARGE = -0.25;
export const debyeKappa = (ionicM: number) => 3.30 * Math.sqrt(Math.max(1e-4, ionicM));

/** Polysome mRNA linkers (worm-like chain between the mRNA exit of one ribosome and the entry of the next). Ribosomes sit
 *  one per 156 nt of ORF on average in growing yeast (Arava et al. 2003, PNAS 100:3889); the ribosome covers the 29 nt
 *  of its channel (library model, from PDB 7ZW0), so a linker holds about 127 nt. Single-stranded RNA: 0.49 nm contour
 *  per nucleotide and persistence length between 1.0 nm (high salt) and 2.5 nm (no added salt) for rU40 (Chen et al.
 *  2012, PNAS 109:799); 1.5 nm is used. Contour length 127 x 0.49 = 62 nm, rms end-to-end distance about 14 nm. */
export const MRNA_LINK_NT = 156 - 29;
export const MRNA_NM_PER_NT = 0.49;
export const MRNA_LP = 1.5;
/** mRNA channel of the 80S library structure (ribosome_80S frame, nm; polysome.py, from PDB 7ZW0): exit (5' end) and
 *  entry (3' end) points and their outward directions; the packer places the next ribosome's entry a 3-12 nm linker
 *  beyond the previous one's exit, within 0.7 rad of the exit direction */
export const RIBO_CHANNEL = {
  exit: [-6.3208, -3.1478, 5.8007], entry: [-11.2279, -1.9972, -2.6486],
  exitDir: [-0.1301, 0.3672, 0.921], entryDir: [-0.6272, -0.7241, -0.2869], linkMin: 3.0, linkMax: 12.0, cone: 0.7,
};

export type SimMode = 'schematic' | 'mid' | 'residue';

export interface SimSettings {
  mode: SimMode;
  dt: number;            // ns per step
  stepsPerFrame: number;
  eps: number;           // kT, depth of the nonspecific attraction between surfaces (0: excluded volume only)
  viscosity: number;     // solvent viscosity relative to water (scales all diffusion coefficients)
  hydro: boolean;        // mean-field hydrodynamic interactions (scale D by the box's volume fraction)
  boundary: 'auto' | 'periodic' | 'walls' | 'embedded';   // auto: periodic without membranes in the box, else embedded
  tether: number;        // kT/nm^2, harmonic tether of the embedded shell (rms excursion sqrt(1/tether) nm per axis)
  memFluid: boolean;     // membrane proteins diffuse laterally (Saffman-Delbrueck) instead of staying fixed
  memVisc: number;       // Pa s, membrane viscosity for that
  memCalibrated: boolean; // scale it per membrane to yeast measurements (library/membrane_diffusion.json)
  kRep: number;          // kT/nm^2, overlap penalty
  shell: number;         // nm, outer layer of the box tethered in the embedded boundary
  recordEvery: number;   // steps between recorded frames
  elec: boolean;         // screened electrostatics between charged spheres
  ionic: number;         // M, ionic strength (Debye length 0.30 / sqrt(I) nm)
  specific: boolean;     // native contacts of predicted complexes (AlphaFold dimers) attract specifically
  epsGo: number;         // kT, energy of one complete native interface, shared over its contact pairs
  mrna: boolean;         // polysome ribosomes linked by their mRNA (worm-like chain)
}

export const MODES = {
  schematic: { label: 'Schematic (rigid sphere trees)', dt: 0.2, steps: 200, kRep: 20 },
  mid: { label: 'Mid (clusters of ~4 residues)', dt: 0.2, steps: 60, kRep: 30 },
  residue: { label: 'Residue detail (rigid, residue-bead contacts)', dt: 0.1, steps: 20, kRep: 50 },
};

/** Largest time step (ns) at which diffusion is converged (schematic mode, 120 nm cytoplasm box, measured: excluded
 *  volume converges by 0.2 ns; the 0.6 nm attraction well needs about 0.05/eps ns, since contact forces summed over
 *  many sphere pairs are stiff and a longer step lets molecules escape the well too easily). Residue mode (80 nm
 *  periodic box): excluded volume gives the same diffusion at 0.05, 0.1 and 0.2 ns within noise; 0.1 ns is used. With
 *  attraction the same 1/eps scaling is assumed, untested. */
export function dtMax(mode: SimMode, eps: number): number {
  const base = MODES[mode].dt;
  return eps > 0 ? Math.min(base, STICKY_DT[mode] / eps) : base;
}
/** step limit with attraction, dt <= STICKY_DT / eps ns. Mid model (80 nm box, 1 kT): GFP-sized D 17.9 um^2/s at
 *  0.05 ns but 11.8 at 0.01 ns, so 0.01 / eps; at 1.5 kT even 0.01 ns was unstable (molecules faster than in water),
 *  hence the per-model ceiling on the stickiness slider (STICKY_MAX) */
export const STICKY_DT: Record<SimMode, number> = { schematic: 0.05, mid: 0.01, residue: 0.0125 };
export const STICKY_MAX: Record<SimMode, number> = { schematic: 2, mid: 1, residue: 0.5 };

/** What stickiness does in each model (cytoplasm boxes, electrostatics and mean-field hydrodynamics on; MODEL.md 14.5).
 *  It slows GFP-sized molecules toward the 8.6 um^2/s measured in yeast mainly by clustering them. */
export const STICKY_NOTE: Record<SimMode, string> = {
  schematic: 'Stickiness reaches the measured value at 1.5–2 kT, but by joining most molecules into one network within a microsecond.',
  mid: 'Stickiness of 1 kT makes most molecules aggregate (about 12 µm²/s); stronger values are unstable and are not offered.',
  residue: 'Stickiness of 0.5 kT makes most molecules aggregate (25–32 µm²/s).',
};

/** Default viscosity: the cytoplasmic slowdown left after excluded volume and mean-field hydrodynamics, calibrated on
 *  the training papers of the diffusion benchmark (validation/diffusion.py, MODEL.md 14.5): s = 0.365, i.e. 2.74x;
 *  held-out papers: within-paper ratios RMSE 0.11 log10 units, free GFP 8.6-21 um^2/s measured vs 12.6 predicted. */
export const VISC_CALIBRATED = 2.75;

export function defaultSettings(mode: SimMode = 'schematic'): SimSettings {
  const m = MODES[mode];
  return { mode, dt: m.dt, stepsPerFrame: m.steps, eps: 0, viscosity: VISC_CALIBRATED, hydro: true, boundary: 'auto', tether: 1, memFluid: true,
    memVisc: MEM_VISC, memCalibrated: true, kRep: m.kRep, shell: 8, recordEvery: 200, elec: true, ionic: 0.3, specific: true, epsGo: 20, mrna: true };
}
