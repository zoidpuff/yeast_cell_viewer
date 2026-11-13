"""Proximity inside mitochondria: the packed cell against quantitative BioID (IMPROVEMENTS I4).

Data: Singh et al. 2020 Mol Cell 79:1051 (doi 10.1016/j.molcel.2020.07.024), Table S1 (mmc3.xlsx): 40 BioID baits of the
mitochondrial gene-expression machinery (mitoribosome, translational activators, RNA processing), each against the
matrix control bait Kgd4 (3 replicates each); per prey a moderated log fold change and q-value. All preys are
mitochondrial, so the test asks whether the model's arrangement inside mitochondria (matrix vs inner membrane vs IMS,
membrane-tethered mitoribosomes) predicts which preys a bait labels, beyond the compartment.
Per bait: enriched preys (q.mod < 0.05 and logFC > 1) vs the other detected preys; AUROC and Spearman(logFC) of
  abundance    copies of the prey
  compartment  expected concentration of the prey in the bait's compartments (inventory weights; proximity_test.ZONES)
  snapshot     P(prey within 10 nm of a bait copy) in the packed cell (contacts.py)
Output: work/validation/singh2020_test.json
"""
import os, sys, json, warnings
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
from config import WORK, DATA
from contacts import Proximity
from proximity_test import ZONES, zone_volumes, auroc

warnings.filterwarnings("ignore")


def main(g=10):
    F = pd.read_parquet(f"{WORK}/features/protein_features.parquet")
    key = {u: (gn if isinstance(gn, str) and gn else o) for u, gn, o in zip(F.index, F.gene, F.orf)}
    gene_key = {str(gn).upper(): k for gn, k in zip(F.gene, (key[u] for u in F.index)) if isinstance(gn, str)}
    P = Proximity()
    Fk = F.assign(key=[key[u] for u in F.index]).drop_duplicates("key").set_index("key")
    zones = list(ZONES); V = zone_volumes(); vol = np.array([V[z] for z in zones])
    W = {k: np.array([json.loads(w).get(z, 0.0) for z in zones]) for k, w in zip(Fk.index, Fk.comp_weights)}
    W = {k: w / max(w.sum(), 1e-12) for k, w in W.items()}
    copies = Fk.copies_per_cell.astype(float).to_dict()
    C = P.matrix(g, "C")
    x = pd.ExcelFile(f"{DATA}/validation/proximity/singh2020_mito_bioid/mmc3.xlsx")
    rows = []
    for bait in x.sheet_names[1:]:
        d = x.parse(bait)
        b = gene_key.get(bait.upper()) or key.get(str(bait).upper())
        if b is None or b not in P.pidx or b not in W:
            continue
        d["prey"] = [key.get(u) for u in d.Uniprot]
        d = d[d.prey.isin(P.pidx) & d.prey.isin(W) & (d.prey != b)].drop_duplicates("prey")
        if len(d) < 20:
            continue
        pos = ((d["q.mod"] < 0.05) & (d.logFC > 1)).values
        ib = P.pidx[b]
        sc = {"abundance": np.array([copies[p] for p in d.prey]),
              "compartment": np.array([(W[b] * W[p] * copies[p] / vol).sum() for p in d.prey]),
              "snapshot": np.array([np.nan_to_num(C[ib, P.pidx[p]]) for p in d.prey])}
        r = dict(bait=bait, preys=len(d), enriched=int(pos.sum()))
        for k, s in sc.items():
            r[f"auroc_{k}"] = auroc(s, pos) if 3 <= pos.sum() < len(pos) else np.nan
            r[f"spearman_{k}"] = float(spearmanr(s, d.logFC).correlation)
        rows.append(r)
    t = pd.DataFrame(rows)
    summ = {c: float(t[c].median()) for c in t.columns if c.startswith(("auroc", "spearman"))}
    res = dict(source="Singh et al. 2020 Mol Cell, Table S1", baits=len(t), median=summ, per_bait=rows)
    json.dump(res, open(f"{WORK}/validation/singh2020_test.json", "w"), indent=1, default=float)
    print(len(t), "baits;", {k: round(v, 3) for k, v in summ.items()})
    return res


if __name__ == "__main__":
    main()
