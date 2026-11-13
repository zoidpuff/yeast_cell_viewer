"""Nucleus / cytoplasm intensity ratios from raw images against the model (IMPROVEMENTS Q1; the first test on pixels).

Data: the C-SWAT library re-imaged for YeastRGB (Dubreuil et al. 2019 Nucleic Acids Res 47:D1245; library: Meurer et
al. 2018 Nat Methods 15:598): 6,145 ORFs tagged C-terminally with mNeonGreen under their native promoter, diploid
hemizygous, every strain with the nuclear-envelope marker Nup49-mCherry; spinning-disk confocal, 60x; per ORF one
16-bit grid of up to 15 x 15 single-cell crops (65 x 65 px, about 0.108 um per pixel) per channel
(RAW-C_SWAT_GFP-ORF.zip, 6 GB). Not an input of the model. The nucleus/cytosol rule of the inventory was tuned with
CYCLoPs scores on the dev half of the proteins (validation/splits.py), so the test half is the held-out number.

Per cell (all choices fixed on four sample strains before the comparison with the model):
  nucleus      marker channel smoothed (sigma 1 px), Otsu threshold, largest object, closed and filled; kept if its
               area is 60-600 px and it does not touch the crop border; interior = eroded by 2 px (inside the envelope)
  cytoplasm    ring 2-7 px outside the filled nucleus; mean of its brighter half (leaves out vacuole and outside-cell
               pixels); N = mean of the interior
Per strain: medians over cells (>= MIN_CELLS). Background: the darkest 10% of the grid. Autofluorescence: medians of
N and C over the dimmest 5% of strains, subtracted; strains whose brighter side is below 3x the cytoplasmic
autofluorescence are dropped (not detected). Measured log2(N/C) is clipped to +-CLIP.
Model: placed copies per region (types.json) -> concentration ratio (nucleoplasm + nucleolus over 2.93 um^3) /
(cytosol over 26.62 um^3), for proteins with >= 90% of their placed copies in nucleus, nucleolus or cytosol.
Blur and out-of-focus light compress measured ratios, so ranks are compared: Spearman, and the AUROC of
model-nuclear (ratio > 2) against model-cytosolic (ratio < 0.5) proteins by their measured ratio.

Abundance (check): the mean cellular mNeonGreen intensity per strain from the same screen (YeastRGB's per-ORF table)
against the model's copies per cell (PaxDb-based input), Spearman over all strains and over the brightest 70%
(the dimmest are at autofluorescence).

  python validation/image_nc.py            # needs the zip (pipeline/data/validation/microscopy/yeastrgb_dubreuil2019/raw/)
Output: work/validation/image_nc.json, image_nc.tsv (per ORF)
"""
import os, sys, io, re, json, zipfile
for _v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import multiprocessing as mp
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline")); sys.path.insert(0, HERE)
from config import WORK, DATA, NPROC  # noqa: E402

ZIP = f"{DATA}/validation/microscopy/yeastrgb_dubreuil2019/raw/RAW-C_SWAT_GFP-ORF.zip"
CROP, MIN_CELLS, CLIP = 65, 10, 5.0
AREA = (60, 600)


def measure_cell(g, r):
    """(N, C) raw means of one cell crop, or None"""
    from scipy import ndimage as ndi
    from skimage.filters import threshold_otsu
    if r.max() <= 0:
        return None
    s = ndi.gaussian_filter(r.astype(np.float32), 1.0)
    if s.max() < 3 * max(np.median(s), 1.0):           # no marker signal above the crop's level
        return None
    m = s > threshold_otsu(s)
    lab, n = ndi.label(m)
    if n == 0:
        return None
    k = 1 + int(np.argmax(ndi.sum(m, lab, range(1, n + 1))))
    m = ndi.binary_fill_holes(ndi.binary_closing(lab == k, iterations=2))
    a = int(m.sum())
    if not (AREA[0] <= a <= AREA[1]) or m[0].any() or m[-1].any() or m[:, 0].any() or m[:, -1].any():
        return None
    inner = ndi.binary_erosion(m, iterations=2)
    ring = ndi.binary_dilation(m, iterations=7) & ~ndi.binary_dilation(m, iterations=2)
    if inner.sum() < 15 or ring.sum() < 50:
        return None
    c = g[ring]
    return float(g[inner].mean()), float(c[c >= np.median(c)].mean())


def measure_grid(G, R):
    """per strain: n cells, medians of N and C above the grid background"""
    bg = float(np.median(G[G <= np.percentile(G, 10)]))
    out = []
    for i in range(0, G.shape[0] - CROP + 1, CROP):
        for j in range(0, G.shape[1] - CROP + 1, CROP):
            x = measure_cell(G[i:i + CROP, j:j + CROP].astype(np.float32), R[i:i + CROP, j:j + CROP])
            if x:
                out.append(x)
    if len(out) < MIN_CELLS:
        return len(out), np.nan, np.nan
    a = np.array(out)
    return len(out), float(np.median(a[:, 0]) - bg), float(np.median(a[:, 1]) - bg)


def _job(names):
    import tifffile
    orf, gname, rname = names
    with zipfile.ZipFile(ZIP) as z:
        G = tifffile.imread(io.BytesIO(z.read(gname))); R = tifffile.imread(io.BytesIO(z.read(rname)))
    return (orf, *measure_grid(G, R))


def measure_all():
    with zipfile.ZipFile(ZIP) as z:
        names = z.namelist()
    g = {re.search(r"\.([^.]+)\.GFP\.tif$", n).group(1): n for n in names if n.endswith(".GFP.tif")}
    r = {re.search(r"\.([^.]+)\.RFP\.tif$", n).group(1): n for n in names if n.endswith(".RFP.tif")}
    jobs = [(o, g[o], r[o]) for o in sorted(g) if o in r]
    with mp.get_context("fork").Pool(NPROC) as pool:
        rows = pool.map(_job, jobs, chunksize=8)
    return pd.DataFrame(rows, columns=["orf", "cells", "N", "C"]).set_index("orf")


def ratios(m):
    """autofluorescence-corrected log2(N/C) per detected strain"""
    m = m.dropna()
    dim = m[(m.N + m.C) <= (m.N + m.C).quantile(0.05)]
    aN, aC = float(dim.N.median()), float(dim.C.median())
    n, c = m.N - aN, m.C - aC
    det = np.maximum(n, c) > 3 * aC
    eps = 0.5 * aC
    lr = np.log2(np.maximum(n, eps) / np.maximum(c, eps)).clip(-CLIP, CLIP)
    return lr[det], dict(strains_measured=int(len(m)), detected=int(det.sum()), autofluorescence_N=aN, autofluorescence_C=aC)


def model_ratio():
    import localization as L
    df = L.model_shares()
    p = lambda c: df.get(f"placed_{c}", pd.Series(0.0, index=df.index))
    nuc, cyt = p("nucleus") + p("nucleolus"), p("cytosol")
    ok = (df.get("placed_copies", 0) > 0) & (nuc + cyt >= 0.9)
    lo = 2.0 ** -8
    lr = np.log2(np.maximum(nuc / (L.CLASS_VOL["nucleus"] + L.CLASS_VOL["nucleolus"]), lo * 1e-3)
                 / np.maximum(cyt / L.CLASS_VOL["cytosol"], lo * 1e-3)).clip(-8, 8)
    return lr[ok]


def abundance():
    from scipy.stats import spearmanr
    import localization as L
    t = json.load(open(f"{DATA}/validation/microscopy/yeastrgb_dubreuil2019/C_SWAT_GFP-by-ORF.JSON"))["data"]

    def fl(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return np.nan
    g = pd.Series({o: fl(v.get("avgGFP")) for o, v in t.items()}).dropna()
    x = pd.DataFrame(dict(gfp=g, copies=L.model_shares().copies.reindex(g.index))).dropna()
    x = x[(x.copies > 0) & (x.gfp > 0)]
    b = x[x.gfp > x.gfp.quantile(0.3)]
    return dict(strains=int(len(x)), spearman=float(spearmanr(x.gfp, x.copies).correlation),
                spearman_brightest_70pct=float(spearmanr(b.gfp, b.copies).correlation),
                log10_slope_brightest_70pct=float(np.polyfit(np.log10(b.copies), np.log10(b.gfp), 1)[0]))


def main():
    from scipy.stats import spearmanr
    from sklearn.metrics import roc_auc_score
    from splits import protein_split
    cache = f"{WORK}/validation/image_nc_raw.tsv"
    if os.path.exists(cache):
        m = pd.read_csv(cache, sep="\t", index_col=0)
    else:
        m = measure_all(); m.to_csv(cache, sep="\t")
    meas, info = ratios(m)
    mod = model_ratio()
    orfs = meas.index.intersection(mod.index)
    d = pd.DataFrame(dict(measured_log2=meas[orfs], model_log2=mod[orfs]))
    d["split"] = [protein_split(o) for o in d.index]
    res = dict(info, compared=int(len(d)), clip_log2=CLIP)
    for sp in ("dev", "test"):
        x = d[d.split == sp]
        nucl, cyto = x.model_log2 > 1, x.model_log2 < -1
        y = np.r_[np.ones(int(nucl.sum())), np.zeros(int(cyto.sum()))]
        res[sp] = dict(n=int(len(x)), spearman=float(spearmanr(x.measured_log2, x.model_log2).correlation),
                       model_nuclear=int(nucl.sum()), model_cytosolic=int(cyto.sum()),
                       auroc_nuclear_vs_cytosolic=float(roc_auc_score(y, np.r_[x.measured_log2[nucl], x.measured_log2[cyto]])),
                       measured_median_log2=dict(model_nuclear=float(x.measured_log2[nucl].median()),
                                                 model_mixed=float(x.measured_log2[~nucl & ~cyto].median()) if (~nucl & ~cyto).any() else None,
                                                 model_cytosolic=float(x.measured_log2[cyto].median())))
    res["abundance"] = abundance()
    d.to_csv(f"{WORK}/validation/image_nc.tsv", sep="\t")
    json.dump(res, open(f"{WORK}/validation/image_nc.json", "w"), indent=1)
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    main()
