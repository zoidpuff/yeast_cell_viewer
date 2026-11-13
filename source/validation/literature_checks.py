"""Checks of the model against measurements extracted from papers (2026-11-11 data upload; tables in
pipeline/data/validation/literature/*/ with verbatim quotes, see literature_numbers_README.md there).

  boeke    FCCS of endocytic proteins in the cytoplasm (Boeke et al. 2014 Mol Syst Biol 10:756): measured cytoplasmic
           concentrations vs the model's; detected vs undetected binding pairs vs the model's contact statistics
  organelles  soft X-ray tomography (Chen et al. 2025 Mol Biol Cell 36:ar132; BY4741, YPD, unbudded G1) and
           mitochondrial volume fraction (Di Bartolomeo et al. 2020 PNAS 117:7524; CEN.PK, minimal glucose)
  counts   translocon-engaged ER ribosomes (Wu et al. 2026 PLoS Comput Biol) and TORC1 per cell (Prouteau et al. 2017)
Output: work/validation/literature_checks.json
"""
import os, sys, json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
from config import WORK, DATA, WEB_OUT
from contacts import Proximity, proteins_of
from proximity_test import auroc

LIT = f"{DATA}/validation/literature"
NA = 6.02214e23
import volumes  # noqa: E402
V_CYTO_L = volumes.measured()["cytosol"] * 1e-15
CYTO_REGIONS = {"cytosol", "PM_peripheral", "NE_peripheral", "ER_peripheral", "LD_peripheral", "VAC_peripheral",
                "GOLGI_peripheral", "PEROX_peripheral"}


def instance_regions():
    import tiles3d as t3
    Z = np.load(f"{WORK}/contacts/instances.npz")
    types = json.load(open(f"{WEB_OUT}/library/types.json"))["types"]
    cnt = np.zeros((len(types), len(t3.REGION_NAMES)), np.int64)
    np.add.at(cnt, (Z["type"].astype(np.int64), Z["region"].astype(np.int64)), 1)
    names, off, mem = proteins_of(types)
    return types, cnt, names, off, mem, t3.REGION_NAMES


def boeke(types, cnt, names, off, mem, rn):
    d = pd.read_csv(f"{LIT}/boeke2014_fccs/boeke2014_pairs.tsv", sep="\t")
    pidx = {n.upper(): k for k, n in enumerate(names)}
    cyto = np.array([r in CYTO_REGIONS for r in rn])
    prot_types = {}
    for t in range(len(types)):
        for k in mem[off[t]:off[t + 1]]:
            prot_types.setdefault(k, []).append(t)

    def model_nM(g):
        k = pidx.get(str(g).upper())
        if k is None:
            return np.nan
        n = cnt[prot_types.get(k, [])][:, cyto].sum()
        return n / (V_CYTO_L * NA) * 1e9
    # measured cytoplasmic concentration per protein (mean over the rows it appears in)
    meas = {}
    for col_p, col_c in (("protein_A_green_3myeGFP", "conc_A_nM"), ("protein_B_red_3mCherry", "conc_B_nM")):
        for p, c in zip(d[col_p], d[col_c]):
            if pd.notna(c):
                meas.setdefault(p, []).append(float(c))
    conc = pd.DataFrame([dict(protein=p, measured_nM=float(np.mean(v)), model_nM=model_nM(p)) for p, v in meas.items()])
    ok = conc.dropna(); ok = ok[ok.model_nM > 0]
    conc_res = dict(proteins=len(ok), spearman=float(spearmanr(ok.measured_nM, ok.model_nM).correlation),
                    median_ratio_model_over_measured=float(np.median(ok.model_nM / ok.measured_nM)),
                    within_2x=float((abs(np.log10(ok.model_nM / ok.measured_nM)) < np.log10(2)).mean()),
                    table=ok.round(2).to_dict("records"))
    # binding: detected vs not (hetero pairs among endocytic proteins; controls excluded)
    P = Proximity()
    C2 = P.matrix(2, "C"); C10 = P.matrix(10, "C")
    from scipy.sparse import diags
    same = np.asarray((P.PT @ diags(P.n_inst) @ P.PT.T).todense()) / np.maximum(P.copies[:, None], 1)
    h = d[(d.pair_type == "hetero") & ~d.reported_interaction_TableS2.astype(str).str.startswith("control")].copy()
    h["detected"] = h.cross_correlation_detected.astype(str).str.lower().str.startswith("yes")
    rows = []
    for r in h.itertuples(index=False):
        a, b = str(r[0]).upper(), str(r[1]).upper()
        ia, ib = P.pidx.get(r[0]) if r[0] in P.pidx else P.pidx.get(a), P.pidx.get(r[1]) if r[1] in P.pidx else P.pidx.get(b)
        if ia is None or ib is None:
            continue
        rows.append(dict(A=r[0], B=r[1], detected=bool(r.detected), p2=float(np.nan_to_num(max(C2[ia, ib], C2[ib, ia]))),
                         p10=float(np.nan_to_num(max(C10[ia, ib], C10[ib, ia]))),
                         same_complex=float(max(same[ia, ib], same[ib, ia]))))
    b_ = pd.DataFrame(rows)
    bind = dict(pairs=len(b_), detected=int(b_.detected.sum()),
                auroc_p2=auroc(b_.p2, b_.detected), auroc_p10=auroc(b_.p10, b_.detected),
                auroc_same_complex=auroc(b_.same_complex, b_.detected), table=b_.round(4).to_dict("records"))
    return dict(concentrations=conc_res, binding=bind)


def organelles():
    c = pd.read_csv(f"{LIT}/chen2025_sxt/chen2025_wt_summary.tsv", sep="\t")
    v = json.load(open(f"{WORK}/results/wholecell3d_validation.json"))["density_mg_per_ml"]
    import inventory as inv, volumes
    m = volumes.measured()                               # nucleus: the target fraction (it includes the envelope lumen)
    model_frac = {"nucleus": inv.FRAC["nucleus"], "vacuole": m["vacuole"] / m["cell"],
                  "lipid droplets": m["lipid_droplet"] / m["cell"], "mitochondria": m["mitochondrion"] / m["cell"]}
    out = dict(model_volume_fraction=model_frac,
               model_density_ratio_to_cytosol={"nucleoplasm": v["nucleoplasm"] / v["cytosol"],
                                               "mito_matrix": v["mito_matrix"] / v["cytosol"]},
               chen2025_rows=c[["strain", "organelle", "metric", "mean", "sd", "n", "unit", "value_status"]].to_dict("records"))
    m = pd.read_csv(f"{LIT}/dibartolomeo2020/dibartolomeo2020_mito.tsv", sep="\t")
    out["dibartolomeo2020_rows"] = m[["strain_condition", "phase", "quantity", "value", "error", "unit", "value_status"]].to_dict("records")
    return out


def counts(types, cnt, names, off, mem, rn):
    w = pd.read_csv(f"{LIT}/wu2026_er/wu2026_er_ribosomes.tsv", sep="\t")
    p = pd.read_csv(f"{LIT}/prouteau2017/prouteau2017.tsv", sep="\t")
    ids = [t["id"] for t in types]
    docked = int(sum(cnt[i].sum() for i, t in enumerate(ids) if t in ("ribosome_80S_Sec61", "ribosome_80S_Ssh1")))
    kog = names.index("KOG1") if "KOG1" in names else None
    torc1 = int(sum(cnt[t].sum() for t in range(len(types)) if kog in set(mem[off[t]:off[t + 1]])))
    return dict(model_translocon_ribosomes=docked, wu2026=w.to_dict("records"), model_TORC1_instances=torc1,
                prouteau2017=p.to_dict("records"))


def main():
    types, cnt, names, off, mem, rn = instance_regions()
    res = dict(boeke2014=boeke(types, cnt, names, off, mem, rn), organelles=organelles(),
               counts=counts(types, cnt, names, off, mem, rn))
    json.dump(res, open(f"{WORK}/validation/literature_checks.json", "w"), indent=1, default=str)
    b = res["boeke2014"]
    print("Boeke concentrations:", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in b["concentrations"].items() if k != "table"})
    print("Boeke binding:", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in b["binding"].items() if k != "table"})
    print("model density ratios:", {k: round(x, 2) for k, x in res["organelles"]["model_density_ratio_to_cytosol"].items()})
    print("translocon ribosomes model", res["counts"]["model_translocon_ribosomes"], "| TORC1 instances model", res["counts"]["model_TORC1_instances"])
    return res


if __name__ == "__main__":
    main()
