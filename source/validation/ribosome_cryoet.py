"""Cytosolic ribosomes of the packed cell against cryo-electron tomography (IMPROVEMENTS G, I1).

Data (independent of the model):
  - concentration: about 23 uM cytosolic ribosomes in log-phase S. cerevisiae by in-situ cryo-ET template matching
    (Delarue et al. 2018 Cell 174:338, doi 10.1016/j.cell.2018.05.042)
  - nearest-neighbour distances: 91,155 crYOLO ribosome picks in 249 lamella tomograms (CZ cryoET Data Portal
    dataset 10007; Rangan et al. 2024 Nat Methods, doi 10.1038/s41592-024-02340-4). Picks are incomplete and not
    restricted to cytosol (their density, ~5,200 per um^3 of convex hull, is a lower bound).
Model: centres of all placed 80S ribosomes (every library type whose id starts with ribosome_80S, incl. polysomes,
RAC-bound and ER-docked) in the cytosol region. Nearest-neighbour distances are computed for the full set and for
random subsets thinned to the picks' density (incomplete detection lengthens NN distances).
Output: work/validation/ribosome_cryoet.json and docs/figures/ribosome_cryoet.png
"""
import os, sys, json, glob
import numpy as np
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
from config import WORK, DATA, WEB_OUT

VOX_A = 7.84
NA = 6.02214e23


def model_ribosomes():
    Z = np.load(f"{WORK}/contacts/instances.npz")
    types = json.load(open(f"{WEB_OUT}/library/types.json"))["types"]
    ribo = np.array([t["id"].startswith("ribosome_80S") for t in types])
    m = ribo[Z["type"].astype(np.int64)]
    return Z["pos"][m].astype(np.float64), Z["region"][m]


def cryoet_nn():
    out = []
    for f in sorted(glob.glob(f"{DATA}/validation/microscopy/czii_cryoet_10007_yeast_ribosomes/ribosome_points/*.ndjson")):
        pts = []
        for line in open(f):
            loc = json.loads(line).get("location")
            if loc:
                pts.append([loc["x"], loc["y"], loc["z"]])
        if len(pts) >= 20:
            p = np.array(pts, float) * VOX_A / 10.0
            d, _ = cKDTree(p).query(p, k=2)
            out.append(d[:, 1])
    return out


def polysome_tight(n_cyto80):
    """ribosomes in tightly packed polysomes: a chain neighbour whose mRNA entry lies within 8 nm of this ribosome's
    exit (or vice versa), the criterion of Goetz et al. 2025 (bioRxiv, doi 10.64898/2025.12.22.696049: 14.5% of
    ribosomes); Xie et al. 2024 Mol Cell report 18.6% in tightly packed polysomes. Model: pre-placed polysome chains
    (prepass.npz), mRNA channel ends of the 80S library frame (polysome.channel_path: 5' end = exit, 3' end = entry).
    Collided ribosomes (IMPROVEMENTS P5/P11, test): both members of the pairs the walk placed in the collided-disome
    pose (prepass poly_coll; polysome.COLLIDE_NT), against disome profiling in exponential yeast (Zhao et al. 2021
    Genome Biol 22:16, doi 10.1186/s13059-020-02256-0: monosome:disome peak areas 16.3:1 after RNase I, so 5.8% of
    ribosomes are in collided disomes). Build O, before the collided pose, had none: its closest neighbours were 7.7 nm
    apart at the 5th percentile, since random relative orientations clash before the ribosomes come closer."""
    from config import W3D, LIB
    z = np.load(f"{W3D}/prepass.npz")
    pos, q, ch = z["poly_pos"].astype(float), z["poly_q"].astype(float), z["poly_chain"]
    path = np.load(f"{LIB}/ribosome_80S_mrna_channel.npy")
    def rot(v):
        w, u = q[:, :1], q[:, 1:]
        t = 2 * np.cross(u, v)
        return v + w * t + np.cross(u, t)
    ex = pos + rot(np.repeat(path[:1], len(q), 0)); en = pos + rot(np.repeat(path[-1:], len(q), 0))
    d = np.linalg.norm(en[1:] - ex[:-1], axis=1); same = ch[1:] == ch[:-1]
    t = np.zeros(len(pos), bool); k = np.where(same & (d <= 8.0))[0]; t[k] = True; t[k + 1] = True
    cf = z["poly_coll"].astype(bool) if "poly_coll" in z.files else np.zeros(len(pos), bool)
    c = cf | np.r_[cf[1:], False]                      # the trailing ribosome and the one ahead of it
    # collided pairs are 4.8 nm exit-to-entry in the disome structure but 8.2 nm as placed (moved apart so the coarse
    # sphere trees clear; polysome.disome_pose), so they count as tight whatever their placed distance
    t |= c
    return dict(polysome_ribosomes=int(len(pos)), linker_nm_percentiles_5_50_95=np.percentile(d[same], [5, 50, 95]).round(2).tolist(),
                tight_fraction_of_cytosolic_80S=float(t.sum() / n_cyto80), in_chains_fraction=float(len(pos) / n_cyto80),
                goetz2025=0.145, xie2024=0.186,
                collided_fraction_of_cytosolic_80S=float(c.sum() / n_cyto80), zhao2021_disome_fraction=0.058)


def main():
    pos, reg = model_ribosomes()
    cyto = reg == 0
    import volumes
    V = volumes.measured()["cytosol"]                     # um^3 cytosol of the built geometry
    n = int(cyto.sum())
    conc_um3 = n / V
    conc_uM = conc_um3 * 1e15 / NA * 1e6
    p = pos[cyto]
    d_full = cKDTree(p).query(p, k=2)[0][:, 1]
    rng = np.random.default_rng(0)
    keep = 5200.0 / conc_um3
    q = p[rng.random(len(p)) < keep]
    d_thin = cKDTree(q).query(q, k=2)[0][:, 1]
    tomo = cryoet_nn()
    d_et = np.concatenate(tomo)
    res = dict(model_cytosolic_80S=n, model_per_um3=conc_um3, model_uM=conc_uM, delarue2018_uM=23.0,
               ratio_to_delarue=conc_uM / 23.0, model_nn_median_nm=float(np.median(d_full)),
               model_thinned_to_5200_per_um3_nn_median_nm=float(np.median(d_thin)),
               cryoet_tomograms=len(tomo), cryoet_picks=int(len(d_et)), cryoet_nn_median_nm=float(np.median(d_et)),
               cryoet_per_tomogram_median_iqr=[float(x) for x in np.percentile([np.median(t) for t in tomo], [25, 50, 75])])
    res["polysomes"] = polysome_tight(n)
    json.dump(res, open(f"{WORK}/validation/ribosome_cryoet.json", "w"), indent=1)
    print(json.dumps(res, indent=1))
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    import make_figures as mf
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    bins = np.arange(18, 120, 2)
    for d, lab, col in ((d_full, f"model, all cytosolic 80S ({conc_uM:.0f} µM)", mf.BLUE),
                        (d_thin, "model thinned to the picks' density", mf.AQUA),
                        (d_et, f"cryo-ET picks ({len(tomo)} tomograms)", mf.ORANGE)):
        h, e = np.histogram(d, bins, density=True)
        ax.plot(0.5 * (e[1:] + e[:-1]), h, color=col, lw=2, label=f"{lab}; median {np.median(d):.1f} nm")
    ax.set_xlabel("nearest-neighbour distance between ribosome centres (nm)"); ax.set_ylabel("density")
    ax.legend(fontsize=8); ax.set_title("Ribosome spacing: packed cell vs in-situ cryo-ET")
    mf.save(fig, "ribosome_cryoet")
    return res


if __name__ == "__main__":
    main()
