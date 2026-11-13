"""Protein-pair reference sets for the proximity tests (IMPROVEMENTS I3/I4), from BioGRID 5.0.263 (S. cerevisiae):

  kinase_invivo   kinase -> substrate relationships curated from in-vivo evidence (BioGRID PTM relationships, mostly
                  from PhosphoGRID; Sadowski et al. 2013 Database bat026)
  kinase_invitro  in-vitro protein-array phosphorylation (Ptacek et al. 2005 Nature 438:679), BioGRID 'Biochemical
                  Activity' records of that paper -- no localization constraint in the assay
  pca             DHFR protein-fragment complementation in vivo (Tarassov et al. 2008 Science 320:1465); reports pairs
                  whose fusion termini come within ~8 nm
  prox_label      'Proximity Label-MS' records (BioID, TurboID, APEX, ...), bait -> prey
  xlms            'Cross-Linking-MS (XL-MS)' records (in-cell / lysate cross-linking)
  kid_*           Yeast KID kinase -> substrate pairs by evidence type (kid_sets)

Output: work/validation/proximity_sets.tsv (set, A, B, pmid, throughput); A is the bait / kinase where directional.
"""
import os, sys, csv
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pipeline"))

from config import DATA, WORK

BG = f"{DATA}/validation/biogrid"
ORG = f"{BG}/BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-5.0.263.tab3.txt"


def kid_sets():
    """Yeast KID (Sharifpoor et al. 2011 Genome Biol 12:R39, Additional file 3), split by evidence type. KID also scores
    'Localized to same subcellular compartment' and 'LTP Co-localization': those columns are not used, so that the sets
    do not depend on localization.
      kid_invivo     low-throughput in-vivo phosphorylation evidence (site mutants, phospho-peptide loss, phospho-shift,
                     kinase-dependent localization change excluded, in-vivo site mapping)
      kid_invitro    in-vitro evidence only (protein chips, HTP or LTP in-vitro kinase assays, in-vitro site mapping)
      kid_genetic    genetic evidence only (synthetic lethality/suppression, dosage, SGA correlation)"""
    import warnings
    warnings.filterwarnings("ignore")
    d = pd.read_excel(f"{DATA}/validation/literature/sharifpoor2011_kid/gb-2011-12-4-r39-S3.XLS", engine="xlrd", header=8)
    c = list(d.columns)
    has = lambda idx: d[[c[i] for i in idx]].notna().any(axis=1)
    invivo = has([29, 30, 31, 32, 33])
    invitro = has([4, 5, 20, 28])
    genetic = has([8, 9, 10, 11, 12, 21, 22, 23, 24])
    physical = has([6, 13, 14, 15, 17, 18, 19]) | invivo | invitro
    out = []
    for name, m in (("kid_invivo", invivo), ("kid_invitro", invitro & ~invivo & ~has([6])), ("kid_genetic", genetic & ~physical)):
        for r in d[m].itertuples(index=False):
            out.append((name, r[1], r[0], r[3], r[2], "21492431", "KID"))
    return out


def build():
    rows = []
    rel = pd.read_csv(f"{BG}/BIOGRID-PTM-RELATIONSHIPS-5.0.263.ptmrel.txt", sep="\t", dtype=str)
    rel = rel[(rel.Relationship == "kinase") & (rel["Organism ID"] == "559292")]
    ptm = pd.read_csv(f"{BG}/yeast_ptm.txt", sep="\t", dtype=str, usecols=["#PTM ID", "Official Symbol", "Systematic Name",
                                                                             "Post Translational Modification"])
    ptm = ptm.rename(columns={"#PTM ID": "ptm", "Official Symbol": "sub", "Systematic Name": "sub_orf"})
    k = rel.rename(columns={"#PTM ID": "ptm", "Official Symbol": "kin", "Systematic Name": "kin_orf"}).merge(ptm, on="ptm")
    for r in k.drop_duplicates(["kin", "sub"]).to_dict("records"):
        rows.append(("kinase_invivo", r["kin"], r["sub"], r["kin_orf"], r["sub_orf"], r["Pubmed ID"], "curated"))
    g = pd.read_csv(ORG, sep="\t", dtype=str, quoting=csv.QUOTE_NONE)
    g = g[(g["Organism ID Interactor A"] == "559292") & (g["Organism ID Interactor B"] == "559292")]
    sel = {"kinase_invitro": (g["Experimental System"] == "Biochemical Activity") & (g.Modification == "Phosphorylation")
           & (g["Publication Source"] == "PUBMED:16319894"),
           "pca": (g["Experimental System"] == "PCA") & (g["Publication Source"] == "PUBMED:18467557"),
           "prox_label": g["Experimental System"] == "Proximity Label-MS",
           "xlms": g["Experimental System"] == "Cross-Linking-MS (XL-MS)"}
    for name, m in sel.items():
        for r in g[m].itertuples(index=False):
            rows.append((name, r[7], r[8], r[5], r[6], r[14].replace("PUBMED:", ""), r[17]))
    rows += kid_sets()
    df = pd.DataFrame(rows, columns=["set", "A", "B", "A_orf", "B_orf", "pmid", "throughput"])
    df = df[df.A != df.B].drop_duplicates(["set", "A", "B", "pmid"])
    os.makedirs(f"{WORK}/validation", exist_ok=True)
    df.to_csv(f"{WORK}/validation/proximity_sets.tsv", sep="\t", index=False)
    return df


if __name__ == "__main__":
    import sys
    df = build()
    print(df.groupby("set").agg(pairs=("A", "size"), A=("A", "nunique"), B=("B", "nunique"), papers=("pmid", "nunique")))
