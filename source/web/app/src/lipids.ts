// Lipid composition of each membrane, for the procedural lipids (worker) and the inspector (main).
//
// Phospholipid classes (% of phospholipid) and ergosterol : phospholipid (mol/mol) of isolated organelle membranes from
// Zinser et al. 1991 J Bacteriol 173:2026, Tables 4 and 5. ER: light microsomes; Golgi/endosomes: secretory vesicles;
// nuclear envelope: nuclei; lipid-droplet monolayer: ER phospholipids without sterol (Zinser give none for it).
// Mitochondria, microsomes and peroxisomes were from lactate- or oleate-grown cells, the others from YPD.
// Plasma membrane ergosterol: the measured ratio (3.31) would make 77 mol% of the membrane sterol, more than a bilayer
// holds; 40 mol% is used instead (assumption). Plasma membrane leaflets (consensus): PS in the cytoplasmic leaflet,
// PE mostly (80%) there, sphingolipids (in "other") in the exoplasmic leaflet; other membranes are symmetric.

export const LIPID_CLASSES = ['PI', 'PE', 'PC', 'PA', 'PS', 'ergosterol', 'other (sphingolipids, PG, DAG, ...)', 'cardiolipin'];

// per membrane id (data.MEM_NAMES order): [PC, PE, PI, PS, CL, PA, other] in % of phospholipid, ergosterol : phospholipid
const ZINSER: [number[], number][] = [
  [[16.8, 20.3, 17.7, 33.6, 0.2, 3.9, 6.9], 0.40 / 0.60],   // 0 plasma membrane (sterol: assumption, see above)
  [[51.3, 33.4, 7.5, 6.6, 0.4, 0.2, 0.5], 0.26],              // 1 ER (light microsomes)
  [[44.6, 26.9, 15.1, 5.9, 0.5, 2.2, 4.3], 0.64],             // 2 nuclear envelope (nuclei; CL < 1)
  [[46.5, 19.4, 18.3, 4.4, 1.6, 2.1, 7.7], 0.18],             // 3 vacuole
  [[45.6, 32.6, 10.2, 1.2, 5.9, 4.4, 0.0], 0.01],             // 4 mitochondrial outer membrane
  [[38.4, 24.0, 16.2, 3.8, 16.1, 1.5, 0.0], 0.32],            // 5 mitochondrial inner membrane
  [[35.0, 22.3, 19.1, 12.9, 0.7, 1.2, 8.8], 0.75],            // 6 Golgi / endosomes (secretory vesicles)
  [[48.2, 22.9, 15.8, 4.5, 7.0, 1.6, 0.0], 0.39],             // 7 peroxisomes
  [[51.3, 33.4, 7.5, 6.6, 0.4, 0.2, 0.5], 0.0],               // 8 lipid-droplet monolayer
];
// share of each class in the cytoplasmic leaflet of the plasma membrane (0.5 = symmetric)
const PM_INNER: Record<number, number> = { 4: 1.0, 1: 0.8, 6: 0.0 };

/** mol % per class (LIPID_CLASSES order) of membrane mem; leaf 1 = cytoplasmic leaflet, 0 = the other one */
export function lipidComposition(mem: number, leaf: number): number[] {
  const [pl, r] = ZINSER[Math.min(mem, ZINSER.length - 1)];
  const [pc, pe, pi, ps, cl, pa, other] = pl;
  const s = pl.reduce((a, b) => a + b);
  const phos = 1 / (1 + r), erg = r / (1 + r);
  let w = [pi, pe, pc, pa, ps, 0, other, cl].map(x => x / s * phos * 100);
  w[5] = erg * 100;
  if (mem === 0) {             // move phospholipids between leaflets; the sterol share stays the same in both
    w = w.map((x, c) => c === 5 ? x : x * 2 * (leaf ? (PM_INNER[c] ?? 0.5) : 1 - (PM_INNER[c] ?? 0.5)));
    const z = w.reduce((a, b, c) => a + (c === 5 ? 0 : b), 0);
    w = w.map((x, c) => c === 5 ? x : x / z * (100 - w[5]));
  }
  return w;
}

/** cumulative class probabilities per membrane and leaflet: LCUM[mem][leaf] */
export const LCUM: number[][][] = ZINSER.map((_, m) => [0, 1].map(leaf => {
  const w = lipidComposition(m, leaf); const s = w.reduce((a, b) => a + b); let c = 0;
  return w.map(x => (c += x / s));
}));
