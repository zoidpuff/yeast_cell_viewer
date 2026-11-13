"""Calibration of membrane-protein lateral diffusion per membrane (MODEL.md 14.2, IMPROVEMENTS I14).

The simulation moves membrane proteins with the Saffman-Delbrueck coefficient of a fluid bilayer (viscosity 0.1 Pa s,
thickness 4 nm; web/app/src/sim/params.ts), about 4 um^2/s, on every membrane. Yeast measurements are far slower and
differ between membranes. Data: pipeline/data/validation/membrane_diffusion_measurements.tsv (FRAP; every value quoted
from the paper):
  plasma membrane  GFP-Sso1, GFP-Snc1 0.0025-0.0029 um^2/s (Valdez-Taubas & Pelham 2003 Curr Biol 13:1636)
  vacuole          Sso1 tail anchor retargeted to the vacuole (PNTS) 0.068 (same paper)
  ER, NE           Sec61-GFP and Heh2-derived reporters, 0.15-0.92 (Meinema et al. 2013 Traffic 14:487)
Factor per membrane = median over transmembrane proteins in unperturbed cells of D_measured / D_SD(R_h), with R_h of
the protein's library type (3 nm for engineered reporters). Membranes without yeast data (mitochondrial outer and inner,
Golgi/endosomes, peroxisomes) take the ER factor: an assumption (internal membranes without a cell wall or
eisosome/Pma1 network).
Caveats: FRAP measures micrometre-scale, long-time mobility (including domains, obstacles and immobile fractions of
40-80% on the plasma membrane), applied here as one factor on the short-time diffusion of a 30-200 nm box. The plasma
membrane value is approximate in the source (D = r^2 / 4 t_half for a 1.6 um spot).
Writes library/membrane_diffusion.json: {"scale": [PM, ER, NE, VAC, OM, IM, GOLGI, PEROX], ...}
"""
import os, sys, json, math
import numpy as np
import pandas as pd

from config import DATA

MEMBRANES = ["PM", "ER", "NE", "VAC", "MITO_OM", "MITO_IM", "GOLGI", "PEROX"]
KT = 1.380649e-23 * 303.15
ETA_WATER, MEM_VISC, MEM_H = 0.797e-3, 0.1, 4e-9
PERTURBED = ("mutant", "spheroplast", "latrunculin", "trypsin", "treated", "erg6")
TABLE = f"{DATA}/validation/membrane_diffusion_measurements.tsv"


def d_sd(rh_nm):
    """Saffman-Delbrueck D (um^2/s) as in params.ts saffmanDelbruck (radius = half R_h, within 1-5 nm)"""
    a = min(5.0, max(1.0, 0.5 * rh_nm)) * 1e-9
    mh = MEM_VISC * MEM_H
    return KT / (4 * math.pi * mh) * (math.log(mh / (ETA_WATER * a)) - 0.5772) * 1e12


def fit(lib=None, log=print):
    d = pd.read_csv(TABLE, sep="\t")
    d = d[d.D_um2_s.notna() & d.anchor_type.astype(str).str.startswith("TM")]
    d = d[~d.condition.astype(str).str.lower().str.contains("|".join(PERTURBED))]
    rh = []
    for mol in d.molecule.astype(str):
        g = mol.split(" ")[0].upper()
        t = lib.by_id.get(g) if lib is not None else None
        rh.append(lib.rh(t) if t is not None else 3.0)
    d = d.assign(rh=rh, d_sd=[d_sd(r) for r in rh])
    d["ratio"] = d.D_um2_s / d.d_sd
    memmap = {"PM": "PM", "ER": "ER", "NE": "NE", "vacuole": "VAC"}
    fac, n = {}, {}
    for m, g in d.groupby(d.membrane.map(memmap)):
        fac[m] = float(np.median(g.ratio)); n[m] = int(len(g))
    scale = [fac.get(m, fac["ER"]) for m in MEMBRANES]
    log("membrane diffusion factors: " + ", ".join(f"{m} {s:.2g}" + ("" if m in fac else " (= ER, assumed)")
                                                  for m, s in zip(MEMBRANES, scale)))
    return scale, fac, n, d


def write(data_dir, log=print):
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "simgpu"))
    import boxes
    lib = boxes.Library(data_dir)
    scale, fac, n, d = fit(lib, log)
    json.dump(dict(note="factor on the Saffman-Delbrueck D of membrane proteins per membrane; MODEL.md 14.2",
                   membranes=MEMBRANES, scale=[round(s, 6) for s in scale],
                   measured={m: dict(factor=fac[m], n=n[m]) for m in fac}, assumed_equal_to_ER=[m for m in MEMBRANES if m not in fac]),
              open(os.path.join(data_dir, "library", "membrane_diffusion.json"), "w"), indent=1)
    return scale


if __name__ == "__main__":
    from config import WEB_OUT
    s = write(WEB_OUT)
    print(dict(zip(MEMBRANES, s)))
