"""Does the model's proximity predict measured proximity (IMPROVEMENTS I4)?

Reference sets (validation/proximity_sets.py, BioGRID 5.0.263): DHFR PCA (Tarassov 2008), proximity labelling
(BioID/TurboID/APEX...), cross-linking MS, in-vivo kinase -> substrate (curated) and in-vitro array phosphorylation
(Ptacek 2005). For each bait A (kinase, PCA/BioID bait) every protein in the model is a candidate B; known partners
are positives. Scores, from simple to full model:
  abundance   copies of B per cell (the main bias of PCA and proximity labelling)
  compartment expected concentration of B in A's compartments: sum_c w_A(c) w_B(c) N_B / V_c (inventory weights;
              volumes and areas of the built cell, membranes as area x 10 nm plus their peripheral shell)
  snapshot    P(a copy of A has a copy of B within gap g) in the packed cell, separate instances only (pipeline/contacts.py)
  +assembly   snapshot or B in the same complex instance (curated complexes, Complex Portal / PDB: partly the same
              kind of evidence as the reference, so reported separately)
Metric: per-bait AUROC (ties averaged), median over baits with >= 5 positives in the model; pooled AUROC over all
bait-candidate pairs. Circularity: localization inputs (CYCLoPs, Huh 2003) are independent of PCA/BioID/XL-MS.
Output: work/validation/proximity_test.json
"""
import os, sys, json
import numpy as np
import pandas as pd
from scipy.stats import rankdata

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
from config import WORK, W3D
from contacts import Proximity

ZONES = {   # inventory compartment -> (volume regions, membrane regions) of tile_measure
    "cytosol": (["cytosol"], []), "nucleoplasm": (["nucleoplasm"], []), "nucleolus": (["nucleolus"], []),
    "ER": (["ER_lumen", "ER_peripheral"], ["mem_ER"]), "plasma_membrane": (["PM_peripheral"], ["mem_PM"]),
    "nuclear_envelope": (["NE_peripheral"], ["mem_NE"]), "mito_IM": ([], ["mem_MITO_IM"]), "mito_OM": ([], ["mem_MITO_OM"]),
    "mito_matrix": (["mito_matrix"], []), "mitochondrion": (["mito_matrix"], []), "mito_IMS": (["mito_IMS"], []),
    "golgi_endosome": (["golgi_endosome", "GOLGI_peripheral"], ["mem_GOLGI"]), "vacuole": (["vacuole", "VAC_peripheral"], ["mem_VAC"]),
    "peroxisome": (["peroxisome", "PEROX_peripheral"], ["mem_PEROX"]), "lipid_droplet": (["LD_peripheral"], []), "cell_wall": (["cell_wall"], [])}
MEM_H = 10.0


def zone_volumes():
    tm = pd.read_parquet(f"{W3D}/tile_measure.parquet")
    s = tm.sum()
    return {z: sum(s[v] for v in vol) + MEM_H * sum(s[m] for m in mem) for z, (vol, mem) in ZONES.items()}


def auroc(score, pos):
    pos = np.asarray(pos, bool)
    n1, n0 = pos.sum(), (~pos).sum()
    if n1 == 0 or n0 == 0:
        return np.nan
    r = rankdata(score)
    return float((r[pos].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def main(g=10):
    F = pd.read_parquet(f"{WORK}/features/protein_features.parquet")
    P = Proximity()
    names = P.proteins
    key_by_orf = {o: (gname if isinstance(gname, str) and gname else o) for o, gname in zip(F.orf, F.gene)}
    # model protein table aligned to contacts' protein list
    Fk = F.assign(key=[key_by_orf[o] for o in F.orf]).drop_duplicates("key").set_index("key")
    have = [n for n in names if n in Fk.index]
    idx = np.array([P.pidx[n] for n in have])
    copies = Fk.loc[have, "copies_per_cell"].values.astype(float)
    V = zone_volumes()
    zones = list(ZONES)
    W = np.array([[json.loads(w).get(z, 0.0) for z in zones] for w in Fk.loc[have, "comp_weights"]])
    W = W / np.maximum(W.sum(1, keepdims=True), 1e-12)
    vol = np.array([V[z] for z in zones])
    conc = W * copies[:, None] / vol[None, :]                 # molecules / nm^3 of each protein per zone
    Pn = P.matrix(g, "C")[np.ix_(idx, idx)]
    # same complex instance: fraction of A's instances whose type also contains B
    from scipy.sparse import diags
    PT = P.PT[idx]
    same = np.asarray((PT @ diags(P.n_inst) @ PT.T).todense()) / np.maximum(P.copies[idx][:, None], 1)
    np.fill_diagonal(same, 0)
    Pasm = 1 - (1 - np.nan_to_num(Pn)) * (1 - np.clip(same, 0, 1))
    loc = {n: k for k, n in enumerate(have)}
    sets = pd.read_csv(f"{WORK}/validation/proximity_sets.tsv", sep="\t", dtype=str)
    sets["a"] = [key_by_orf.get(o, s) for o, s in zip(sets.A_orf, sets.A)]
    sets["b"] = [key_by_orf.get(o, s) for o, s in zip(sets.B_orf, sets.B)]
    share = same > 0
    res = {}; per_bait = []
    for name, d in sets.groupby("set"):
        d = d[d.a.isin(loc) & d.b.isin(loc)]
        sym = name in ("pca", "xlms")
        pairs = set(zip(d.a, d.b)) | (set(zip(d.b, d.a)) if sym else set())
        baits = sorted({a for a, _ in pairs})
        per = {k: [] for k in ("abundance", "compartment", "snapshot", "snapshot+assembly", "snapshot, no shared complex")}
        pooled = {k: [] for k in per}; lab = []; rows = []
        for a in baits:
            ia = loc[a]
            posb = np.array([(a, b) in pairs for b in have]); posb[ia] = False
            m = np.ones(len(have), bool); m[ia] = False
            sc = {"abundance": copies, "compartment": (W[ia][None, :] * conc).sum(1),
                  "snapshot": np.nan_to_num(Pn[ia]), "snapshot+assembly": Pasm[ia]}
            nc = m & ~share[ia]                     # candidates sharing no complex type with the bait
            for k, s in sc.items():
                pooled[k].append(s[m])
                if posb[m].sum() >= 5:
                    per[k].append(auroc(s[m], posb[m]))
            if posb[nc].sum() >= 5:
                per["snapshot, no shared complex"].append(auroc(sc["snapshot"][nc], posb[nc]))
                rows.append(dict(set=name, bait=a, positives=int(posb[m].sum()), copies=float(copies[ia]),
                                 **{f"auroc_{k}": auroc(sc[k][m], posb[m]) for k in sc},
                                 auroc_nocomplex=per["snapshot, no shared complex"][-1]))
            lab.append(posb[m])
        lab = np.concatenate(lab) if lab else np.zeros(0, bool)
        rng = np.random.default_rng(0)
        ci = {k: [float(x) for x in np.percentile([np.nanmedian(rng.choice(v, len(v))) for _ in range(2000)], [2.5, 97.5])]
              for k, v in per.items() if len(v) > 1}
        res[name] = dict(pairs_in_model=len(d), baits=len(baits), baits_ge5=len(per["abundance"]),
                         median_auroc={k: float(np.nanmedian(v)) if v else None for k, v in per.items()},
                         median_auroc_ci95=ci,
                         pooled_auroc={k: auroc(np.concatenate(v), lab) if v else None for k, v in pooled.items() if pooled[k]})
        per_bait.extend(rows)
        print(name, json.dumps(res[name]))
    json.dump(dict(gap_nm=g, sets=res), open(f"{WORK}/validation/proximity_test_g{g}.json", "w"), indent=1)
    pd.DataFrame(per_bait).to_csv(f"{WORK}/validation/proximity_per_bait_g{g}.tsv", sep="\t", index=False)
    return res


if __name__ == "__main__":
    for g in (2, 10):
        print("== gap", g)
        main(g)
