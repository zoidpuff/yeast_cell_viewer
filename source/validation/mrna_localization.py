"""mRNAs at the mitochondrial surface (IMPROVEMENTS M4) against fractionation data (held out).

Model: per gene, the share of its mRNA molecules with a ribosome centre within NEAR_NM of the mitochondrial outer
membrane (polysome chains, prepass poly_chain / poly_chain_gene; free mRNAs count as not localized). Input of the
model: Williams et al. 2014 proximity ribosome profiling (-CHX calls, mrna_localization.py).
Data (held out): Saint-Georges et al. 2008 PLoS One 3:e2293, Table S1: % of each mRNA in the mitochondrial fraction
(MLR) for 794 genes encoding mitochondrial proteins. Their cells grew in galactose (more respiratory than the model's
glucose), so absolute shares are not comparable; the test is the ranking.
  spearman_mlr   Spearman rho between the model share and MLR % over the genes in both
  auroc_class1   AUROC of the model share for MLR class I genes (I-1, I-2: strongly localized) against classes III-IV
ER (IMPROVEMENTS M5, 2026-11-12): the same share for the cortical ER membrane, against the membrane-polysome
enrichment of Chartron, Hunt & Frydman 2016 (Nature 536:224, Supplementary Table 1; YPD, fractionation; held out, the
model's ER genes come from Jan et al. 2014 proximity profiling):
  er_spearman    Spearman rho between the model's ER share and log2 membrane-polysome enrichment
  er_auroc_top   AUROC of the model share for the top decile of enrichment against the rest
Output: work/validation/mrna_localization.json
"""
import os, sys, json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
from config import WORK, W3D, DATA
import cell3d as c3, mrna_genes as mg, mrna_localization as ml

NEAR_NM = 50.0
SG = f"{DATA}/validation/mrna/saintgeorges2008_mito_mlr/pone.0002293.s008.xls"


CHARTRON = f"{DATA}/validation/mrna/chartron2016_srp/41586_2016_BFnature19309_MOESM417_ESM.xlsx"


def model_shares(mem=c3.M_OM):
    z = np.load(f"{W3D}/prepass.npz")
    cell = c3.Cell.load(f"{W3D}/geometry.json")
    G = cell.G; glo = np.full(3, -G[0] - 10.0)
    d, _ = ml.mem_dist(z["poly_pos"].astype(np.float64), mem, *cell.args(glo, -glo, margin=50.0))
    chains, inv = np.unique(z["poly_chain"], return_inverse=True)
    near = np.zeros(len(chains), bool); np.logical_or.at(near, inv, np.abs(d) <= NEAR_NM)
    gene = z["poly_chain_gene"]
    free = pd.Series(z["mrna_gene"][z["mrna_kind"] == 0]).value_counts()
    t = pd.DataFrame(dict(gene=gene, near=near)).groupby("gene").near.agg(["sum", "count"])
    t["free"] = free.reindex(t.index).fillna(0)
    t["share"] = t["sum"] / (t["count"] + t["free"])
    gt = mg.gene_table()[["orf", "gene"]].drop_duplicates("gene").set_index("gene")
    t["orf"] = gt.orf.reindex(t.index).str.upper()
    return t.dropna(subset=["orf"]).set_index("orf"), float(np.mean(np.abs(d) <= NEAR_NM))


def main():
    t, all_near = model_shares()
    s = pd.read_excel(SG, sheet_name="Total")
    s = s.rename(columns={"SYSTEMATIC NAME": "orf", "MLR (%localisation)": "mlr", "MLR Classification": "cls"})
    s["orf"] = s.orf.str.upper()
    x = s.merge(t[["share"]], left_on="orf", right_index=True, how="inner")
    rho = spearmanr(x.mlr, x.share).correlation
    c1 = x.cls.astype(str).str.contains("I-")
    c34 = x.cls.astype(str).str.contains("III|IV")
    au = mannwhitneyu(x.share[c1], x.share[c34]).statistic / (c1.sum() * c34.sum()) if c1.any() and c34.any() else None
    tethered = t.index.isin(ml.mito_genes())
    res = dict(near_nm=NEAR_NM, ribosomes_near_om_fraction=all_near, genes_compared=int(len(x)),
               spearman_mlr=float(rho), auroc_class1_vs_3_4=None if au is None else float(au), n_class1=int(c1.sum()), n_class3_4=int(c34.sum()),
               model_share_median_tethered_genes=float(t.share[tethered].median()), model_share_median_other=float(t.share[~tethered].median()),
               mlr_median_class1=float(x.mlr[c1].median()), mlr_median_class3_4=float(x.mlr[c34].median()),
               note="Saint-Georges 2008: galactose, fractionation; model: glucose, Williams 2014 -CHX calls as input")
    # ER
    te, er_near = model_shares(c3.M_ER)
    ch = pd.read_excel(CHARTRON, sheet_name="Annotation and Enrichment")
    ch["orf"] = ch.ORF.astype(str).str.upper()
    col = [c for c in ch.columns if c.startswith("Membrane polysome enrichment")][0]
    y = ch[["orf", col, "Class"]].dropna(subset=[col]).merge(te[["share"]], left_on="orf", right_index=True, how="inner")
    top = y[col] >= y[col].quantile(0.9)
    au_er = mannwhitneyu(y.share[top], y.share[~top]).statistic / (top.sum() * (~top).sum())
    er_t = te.index.isin(ml.er_genes())
    y.rename(columns={col: "chartron_log2_membrane"}).to_csv(f"{WORK}/validation/mrna_localization_er.tsv", sep="\t", index=False)
    res.update(er_ribosomes_near_fraction=er_near, er_genes_compared=int(len(y)),
               er_spearman=float(spearmanr(y[col], y.share).correlation), er_auroc_top_decile=float(au_er),
               er_model_share_median_tethered_genes=float(te.share[er_t].median()) if er_t.any() else None,
               er_model_share_median_other=float(te.share[~er_t].median()),
               er_note="Chartron 2016: YPD, membrane polysome fractionation (log2 enrichment); model input: Jan 2014 Sec63 proximity")
    json.dump(res, open(f"{WORK}/validation/mrna_localization.json", "w"), indent=1)
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    main()
