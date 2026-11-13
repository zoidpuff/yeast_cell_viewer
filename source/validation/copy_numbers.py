"""Per-protein copy numbers of the model against two absolute proteomes (consistency check; the model's input is PaxDb).

  Ho et al. 2018 Cell Syst 6:192 (unified molecules per cell, mean over 21 studies, after GFP autofluorescence filtering;
      Table S4). Overlaps PaxDb's sources, so not independent.
  Das et al. 2026 Histochem Cell Biol 164:8 (label-free absolute proteome of BY4741 on glucose, Table S1; 2026, not in
      PaxDb or Ho).
  Chong et al. 2015 Cell 161:1413 (CYCLoPs; molecules per cell from GFP intensity of the three wild-type screens,
      calibrated by a linear regression on another abundance scale; Table 1D).
Reports rank correlation, the share within 2-fold, and the dynamic range (10th-90th percentile span) of each source on
the shared proteins. Output: work/validation/copy_numbers.json, docs/figures/copy_numbers.png
"""
import os, sys, json, warnings
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
from config import WORK, DATA

warnings.filterwarnings("ignore")
LIT = f"{DATA}/validation/literature"


def load():
    d = pd.read_excel(f"{LIT}/das2026_peroxisome/418_2026_2458_MOESM2_ESM.xlsx", sheet_name="Table S1")
    idc = [c for c in d.columns if "Majority" in str(c)][0]
    cc = [c for c in d.columns if "Mean copy" in str(c) and "Glc" in str(c)][0]
    d = pd.DataFrame(dict(orf=d[idc].astype(str).str.split(";").str[0].str.strip().str.upper(),
                          das=pd.to_numeric(d[cc], errors="coerce")))
    h = pd.read_excel(f"{LIT}/ho2018/mmc5.xlsx", header=2)
    h = pd.DataFrame(dict(orf=h["Systematic Name"].astype(str).str.upper(),
                          ho=pd.to_numeric(h["Mean molecules per cell"], errors="coerce")))
    F = pd.read_parquet(f"{WORK}/features/protein_features.parquet")
    F = pd.DataFrame(dict(orf=F.orf.str.upper(), model=F.copies_per_cell.astype(float)))
    c = pd.read_excel(f"{DATA}/validation/microscopy/cyclops_images/Table1_GFPIntensityMatrix.xls",
                      sheet_name="D. Molecules(I(g))", engine="xlrd", header=4)
    c = pd.DataFrame(dict(orf=c["ORF"].astype(str).str.upper(),
                          chong=c[["WT1", "WT2", "WT3"]].apply(pd.to_numeric, errors="coerce").mean(axis=1)))
    m = F.merge(d, on="orf").merge(h, on="orf").merge(c, on="orf", how="left")
    return m[(m.model > 0) & (m.das > 0) & (m.ho > 0)].reset_index(drop=True)


def main():
    m = load()
    L = np.log10(m[["model", "das", "ho"]])
    mc = m[m.chong > 0]; Lc = np.log10(mc[["model", "das", "chong"]])
    res = dict(proteins=len(m), totals_on_shared={c: float(m[c].sum()) for c in ("model", "das", "ho")},
               span_10_90={c: float(10 ** (L[c].quantile(0.9) - L[c].quantile(0.1))) for c in L},
               pairs={})
    for a, b in (("model", "das"), ("ho", "das"), ("model", "ho")):
        res["pairs"][f"{a}_vs_{b}"] = dict(spearman=float(spearmanr(L[a], L[b]).correlation),
                                           median_ratio=float(10 ** np.median(L[a] - L[b])),
                                           within_2x=float((abs(L[a] - L[b]) < np.log10(2)).mean()))
    res["chong2015"] = dict(proteins=len(mc), span_10_90=float(10 ** (Lc.chong.quantile(0.9) - Lc.chong.quantile(0.1))),
                            spearman_vs_das=float(spearmanr(Lc.chong, Lc.das).correlation),
                            spearman_vs_model=float(spearmanr(Lc.chong, Lc.model).correlation),
                            model_span_same_proteins=float(10 ** (Lc.model.quantile(0.9) - Lc.model.quantile(0.1))),
                            das_span_same_proteins=float(10 ** (Lc.das.quantile(0.9) - Lc.das.quantile(0.1))))
    print(json.dumps(res, indent=1))
    json.dump(res, open(f"{WORK}/validation/copy_numbers.json", "w"), indent=1)
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    import make_figures as mf
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.6))
    ax = axs[0]
    for c, lab, col in (("ho", "Ho 2018 (unified, GFP + MS)", mf.GRAY), ("model", "model (PaxDb, scaled to 5 pg)", mf.BLUE),
                        ("das", "Das 2026 (absolute MS, BY4741)", mf.ORANGE)):
        h_, e = np.histogram(L[c], np.arange(0, 7.01, 0.15), density=True)
        ax.plot(0.5 * (e[1:] + e[:-1]), h_, color=col, lw=2, label=f"{lab}; 10–90%: {res['span_10_90'][c]:.0f}×")
    ax.set_xlabel("log10 copies per cell"); ax.set_ylabel("density"); ax.legend(fontsize=7.5, loc="upper left")
    ax.set_title("Dynamic range of three proteomes")
    ax = axs[1]
    ax.hexbin(L["das"], L["model"], gridsize=50, bins="log", cmap="Blues", mincnt=1)
    ax.plot([0, 7], [0, 7], color=mf.INK2, lw=1)
    p = res["pairs"]["model_vs_das"]
    ax.text(0.03, 0.93, f"Spearman {p['spearman']:.2f}; within 2×: {p['within_2x']:.0%}", transform=ax.transAxes, fontsize=8, color=mf.INK2)
    ax.set_xlabel("log10 copies, Das 2026"); ax.set_ylabel("log10 copies, model"); ax.set_title("Model vs an independent proteome")
    mf.save(fig, "copy_numbers")
    return res


if __name__ == "__main__":
    main()
