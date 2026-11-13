"""Stoichiometric anchors of the protein scale (external review of 2026-11-12, comment 9; a check).

The proteome is PaxDb abundances scaled to inventory.TOTAL_PROTEIN_PG of protein per cell (3.5 pg since build R). Structures whose copy number per cell is known
independently give a ratio "copies in the inventory / copies the structures need"; if the scale and the relative
abundances were right, proteins that exist only in such a structure would give ratios near 1.
  nucleoporins   per nuclear pore (inventory.N_NPC pores x measured copies per pore, Kim 2018), scaffold and FG nucleoporins
  core histones  per nucleosome (2 of each per nucleosome; measured nucleosome map plus rDNA)
  kinetochore    Ndc80 complex (8 per kinetochore) and DASH/Dam1 subunits (16 per kinetochore), 16 kinetochores
                 (Joglekar et al. 2006); note that these also have free pools in cells
  ribosomes      the inventory's 80S count against ~2 x 10^5 per cell (Warner 1999 Trends Biochem Sci 24:437;
                 von der Haar 2008 BMC Syst Biol 2:87: 187,000 +- 56,000)
Each ratio r implies a protein total of TOTAL_PROTEIN_PG / r if that anchor alone set the scale.
Output: work/validation/anchors.json
"""
import os, sys, json, re
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
from config import WORK, DATA  # noqa: E402


def main():
    import inventory as inv, ingredients as ingm, kinetochores as kin
    df, _ = inv.build()
    pool = dict(zip(df.gene, df.copies_per_cell))
    cp = pd.read_csv(f"{DATA}/complexportal_559292.tsv", sep="\t")
    gene_of_acc = dict(zip(df.uniprot, df.gene))
    ing = pd.read_csv(f"{WORK}/results/inventory/ingredients.csv").set_index("ing_id")
    res = dict(total_protein_pg=inv.TOTAL_PROTEIN_PG)

    def fam_ratios(cid, n, keep=None):
        spec = next(c[2] for c in ingm.COMPLEXES if c[0] == cid)
        fams = ingm.complex_families(cp, spec, gene_of_acc, list(pool))
        out = {}
        for f, (st, gs) in fams.items():
            if keep is None or any(keep.match(str(g)) for g in gs):
                out[f] = sum(pool.get(g, 0) for g in gs) / (st * n)
        return out
    nup = fam_ratios("NPC", inv.N_NPC, ingm.NUP_CAPPED)
    res["nucleoporins_per_pore_need"] = dict(n=len(nup), median=float(np.median(list(nup.values()))),
                                             p10_p90=[float(np.percentile(list(nup.values()), q)) for q in (10, 90)])
    his = fam_ratios("nucleosome", float(ing.loc["nucleosome", "copies_per_cell"]))
    res["core_histones_per_nucleosome_need"] = {k: float(v) for k, v in his.items()}
    res["core_histones_median"] = float(np.median(list(his.values())))
    kt = {g: pool.get(g, 0) / (kin.N_NDC80 * 16) for g in kin.NDC80C}
    kt.update({g: pool.get(g, 0) / (kin.N_DASH * 16) for g in kin.DASH})
    res["kinetochore_outer_per_need"] = dict(median=float(np.median(list(kt.values()))), ndc80_complex={g: float(kt[g]) for g in kin.NDC80C})
    ribo = float(ing[ing.index.str.startswith("ribosome_80S")].copies_per_cell.sum())
    res["ribosomes"] = dict(model=ribo, literature=[187000, 200000], ratio=[ribo / 200000, ribo / 187000])
    res["implied_protein_pg"] = dict(nucleoporins=inv.TOTAL_PROTEIN_PG / res["nucleoporins_per_pore_need"]["median"],
                                     core_histones=inv.TOTAL_PROTEIN_PG / res["core_histones_median"],
                                     ribosomes=[inv.TOTAL_PROTEIN_PG / r for r in res["ribosomes"]["ratio"]])
    json.dump(res, open(f"{WORK}/validation/anchors.json", "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
