"""H5 calibration on the nucleus ONLY (nucleolus labels); P-body and stress-granule labels are not loaded here.
Score: AUROC of the proteins' enrichment in the pre-rRNA-rich phase (linear response at 0.95 x onset)."""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pipeline"))
import numpy as np
import condensates as C
e, l, ps = C.pool_lr(); eps = C.build_eps(e, l, ps); F = C.species_table()
nuc = C.Mixture(F, "nucleus", eps, rna_copies=C.PRE_RIBO, rna_nt=C.RRNA_NT)
lab = np.array([u in C.cdcode_members("Nucleolus") for u in nuc.uniprot])
print("nuclear species", len(nuc.c), "nucleolar members among them", lab.sum(), flush=True)
res = []
grid = [(0, 0, 1), (0.25, 0, 1), (0.5, 0, 1), (1, 0, 1), (2, 0, 1), (0, 0.25, 1), (0, 0.5, 1), (0, 1, 1),
        (0.25, 0.25, 1), (0.5, 0.5, 1), (1, 1, 1), (0.5, 0.25, 1), (0.25, 0.5, 1)]
for r in grid:
    r = np.array(r, float)
    g = nuc.onset(r)
    if g is None:
        print(r, "no onset", flush=True); continue
    en = C.rna_response(nuc, 0.95 * g, r)
    a = C.auroc(en, lab)
    top = [nuc.F.gene.iloc[i] for i in np.argsort(-en)[:10]]
    res.append(dict(ratios=r.tolist(), onset=g, auroc_nucleolus=a, top=top))
    print(np.round(r, 2), f"onset {g:.2f} AUROC nucleolus {a:.3f}", top, flush=True)
json.dump(res, open(os.path.join(C.OUT, "calibration_nucleus.json"), "w"), indent=1)
