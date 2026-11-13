"""Model proximity against a DHFR PCA screen with tested negatives (IMPROVEMENTS I4): PPiSeq, Liu et al. 2020
eLife 9:e62365 (1.66 M ordered pairs tested in SD + methotrexate; native promoters; fragments within ~8 nm).

Per unordered pair (A, B), scores of increasing model content (all symmetric, as expected numbers per cell):
  abundance    N_A N_B (placed particles); the PCA signal itself scales with abundance (Liu et al. 2020)
  compartment  N_A N_B sum_c w_A(c) w_B(c) / V_c        (well-mixed compartments of the inventory)
  snapshot     N_A * mean number of B within g of an A  (packed cell, separate instances; pipeline/contacts.py),
               averaged over both directions
Tests: pooled AUROC of the SD call; AUROC within deciles of N_A N_B (removes the abundance bias: does localization
or packing add information at equal abundance?); Spearman with SD fitness. 'no shared complex' drops pairs that are
subunits of one curated complex type. Nothing is tuned.
Output: work/validation/ppiseq_test.json
"""
import os, sys, json, zipfile
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
from config import WORK, DATA
from contacts import Proximity
from proximity_test import ZONES, zone_volumes, auroc

RAW = f"{DATA}/validation/proximity/liu2020_ppiseq/raw/TableS4_PPI_called_mean_fitness_pvalue_in_each_environment.csv.zip"


def load():
    z = zipfile.ZipFile(RAW)
    name = [n for n in z.namelist() if n.endswith(".csv")][0]
    d = pd.read_csv(z.open(name), usecols=["PPI", "SD_mean_fitness", "SD_PPI_called"])
    d = d[d.SD_PPI_called.notna()]
    ab = d.PPI.str.split("_", n=1, expand=True)
    d["o1"], d["o2"] = ab[0].str.upper(), ab[1].str.upper()
    d = d[d.o1 != d.o2]
    # unordered: a pair is positive if called in either orientation; fitness = max
    d["k1"] = np.where(d.o1 < d.o2, d.o1, d.o2); d["k2"] = np.where(d.o1 < d.o2, d.o2, d.o1)
    return d.groupby(["k1", "k2"], as_index=False).agg(called=("SD_PPI_called", "max"), fitness=("SD_mean_fitness", "max"))


def main(g=10):
    d = load()
    F = pd.read_parquet(f"{WORK}/features/protein_features.parquet")
    key = {o: (gn if isinstance(gn, str) and gn else o) for o, gn in zip(F.orf, F.gene)}
    P = Proximity()
    Fk = F.assign(key=[key[o] for o in F.orf]).drop_duplicates("key").set_index("key")
    d["a"], d["b"] = d.k1.map(key), d.k2.map(key)
    d = d[d.a.isin(P.pidx) & d.b.isin(P.pidx) & d.a.isin(Fk.index) & d.b.isin(Fk.index)].copy()
    ia, ib = d.a.map(P.pidx).values, d.b.map(P.pidx).values
    N = P.copies.astype(float)
    M = P.matrix(g, "M")
    zones = list(ZONES); V = zone_volumes(); vol = np.array([V[z] for z in zones])
    W = np.array([[json.loads(w).get(z, 0.0) for z in zones] for w in Fk["comp_weights"]])
    W = W / np.maximum(W.sum(1, keepdims=True), 1e-12)
    wrow = {k: i for i, k in enumerate(Fk.index)}
    wa, wb = W[d.a.map(wrow).values], W[d.b.map(wrow).values]
    from scipy.sparse import diags
    same = (P.PT @ diags(P.n_inst) @ P.PT.T).tocsr()
    d["shared_complex"] = np.asarray(same[ia, ib]).ravel() > 0
    d["abundance"] = N[ia] * N[ib]
    d["compartment"] = d.abundance * (wa * wb / vol[None, :]).sum(1)
    d["snapshot"] = 0.5 * (N[ia] * np.nan_to_num(M[ia, ib]) + N[ib] * np.nan_to_num(M[ib, ia]))
    d = d[d.abundance > 0]
    scores = ["abundance", "compartment", "snapshot"]
    out = dict(gap_nm=g, pairs=int(len(d)), positives=int(d.called.sum()), shared_complex_pairs=int(d.shared_complex.sum()))
    for label, dd in (("all", d), ("no shared complex", d[~d.shared_complex])):
        y = dd.called.values.astype(bool)
        dec = pd.qcut(np.log(dd.abundance), 10, labels=False, duplicates="drop")
        res = {}
        for s in scores:
            strat = [auroc(dd[s].values[dec == q], y[dec == q]) for q in range(dec.max() + 1)]
            res[s] = dict(pooled_auroc=auroc(dd[s].values, y), auroc_within_abundance_deciles=float(np.nanmean(strat)),
                          spearman_fitness=float(spearmanr(dd[s], dd.fitness, nan_policy="omit").correlation))
        out[label] = dict(pairs=int(len(dd)), positives=int(y.sum()), scores=res)
        print(label, len(dd), int(y.sum()), json.dumps({k: {a: round(b, 3) for a, b in v.items()} for k, v in res.items()}))
    json.dump(out, open(f"{WORK}/validation/ppiseq_test_g{g}.json", "w"), indent=1)
    return out


if __name__ == "__main__":
    for g in (5, 10):
        print("== gap", g); main(g)
