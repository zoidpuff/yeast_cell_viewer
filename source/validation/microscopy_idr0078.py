"""Simulated micrographs against an image screen's per-cell object counts (IMPROVEMENTS I2).

Data: Mattiazzi Usaj et al. 2020 Mol Syst Biol 16:e9243 (IDR idr0078), wild-type (his3 deletion) control cells of the
genome-wide screens: CellProfiler objects per cell for Sac6-tdTomato (actin patches), Sla1-yEGFP (endocytic sites),
Snf7-yEGFP (late endosomes) and Vph1-yEGFP (vacuoles); 33-67 k cells each. Log phase, synthetic medium, budding
population. Imaging: Opera spinning disk, 60x water (NA 1.2, n 1.333), 5 z-sections 0.8 um apart, maximum projection,
0.1076 um pixels.
Model: the same optics (confocal approximation of the spinning disk), 5 planes 0.8 um apart through the cell centre
offset by -0.4, 0, +0.4 um, along each of the three axes (9 views of the one model cell). Objects are counted two ways
inside the cell (TDH3 mask): contrast spots (microscopy.count_spots at this pixel size) and connected components above
an Otsu threshold (closer to CellProfiler's segmentation). Not the same software as the screen: compare orders of
magnitude and ranking, not exact numbers.
Output: work/validation/microscopy_idr0078.json, docs/figures/microscopy_idr0078.png
"""
import os, sys, json
import numpy as np
import pandas as pd
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
from config import WORK, DATA
import microscopy as mic

OPT = dict(na=1.2, n_imm=1.333, pixel=107.6, mode="confocal", zrange=1200.0)
MARKERS = {"SAC6": "Sac6", "SLA1": "Sla1", "SNF7": "Snf7", "VPH1": "Vph1"}


def otsu(v):
    h, e = np.histogram(v, 256); c = 0.5 * (e[1:] + e[:-1])
    w0 = np.cumsum(h); w1 = w0[-1] - w0
    m0 = np.cumsum(h * c) / np.maximum(w0, 1); m1 = ((h * c).sum() - np.cumsum(h * c)) / np.maximum(w1, 1)
    return c[np.argmax(w0 * w1 * (m0 - m1) ** 2)]


def views(cell):
    ctr = cell.pos.mean(0)
    for ax in range(3):
        for off in (-400.0, 0.0, 400.0):
            yield ax, [ctr[ax] + off + dz for dz in (-1600, -800, 0, 800, 1600)]


def main():
    cell = mic.Cell()
    tdh = cell.positions("TDH3")
    real = pd.read_csv(f"{DATA}/validation/microscopy/idr0078_mattiazziusaj2020/wt_his3_control_cells_features.csv")
    res = {}
    for gene, lab in MARKERS.items():
        pts = cell.positions(gene)
        spots, objs = [], []
        for ax, zs in views(cell):
            mip = np.max([mic.render(cell, pts, axis=ax, z=z, photons=30, seed=int(z) % 997, **OPT) for z in zs], 0)
            ref = np.max([mic.render(cell, tdh, axis=ax, z=z, **OPT) for z in zs], 0)
            mask = ref > 0.25 * np.percentile(ref, 99)
            n, _ = mic.count_spots(mip[None], mask, sigma_px=0.6, contrast=0.3, min_sep_px=2, ring=(3, 5))
            sm = ndimage.gaussian_filter(mip.astype(float), 0.7)
            lab_, k = ndimage.label((sm > otsu(sm[mask])) & mask)
            sizes = ndimage.sum(np.ones_like(sm), lab_, range(1, k + 1))
            spots.append(n); objs.append(int((sizes >= 3).sum()))
        col = [c for c in real.columns if "children" in c and "count" in c][0]
        r = real.loc[real.marker.str.contains(lab, case=False), col]
        res[gene] = dict(model_copies=int(len(pts)), model_spots=spots, model_objects=objs,
                         real_cells=int(len(r)), real_median=float(r.median()),
                         real_p10_p90=[float(r.quantile(0.1)), float(r.quantile(0.9))])
        print(gene, res[gene])
    json.dump(res, open(f"{WORK}/validation/microscopy_idr0078.json", "w"), indent=1)
    import make_figures as mf
    mf.idr0078()
    return res


if __name__ == "__main__":
    main()
