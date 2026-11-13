"""H5 test (run once, parameters frozen from cond_calibrate.py): cytosol, log phase and in-silico stress.
Frozen: ratios pool:stick:RNA = 1:0:1, overall scale = the nucleus onset (same physics in both compartments)."""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pipeline"))
import numpy as np
import condensates as C
R = np.array([1.0, 0.0, 1.0])
e, l, ps = C.pool_lr(); eps = C.build_eps(e, l, ps); F = C.species_table()
nuc = C.Mixture(F, "nucleus", eps, rna_copies=C.PRE_RIBO, rna_nt=C.RRNA_NT)
g = nuc.onset(R)
pb, sg, gb = C.cdcode_members("P-body"), C.cdcode_members("Stress granule"), C.cdcode_members("G-body")
out = dict(ratios=R.tolist(), scale_nucleus_onset=g)
rng = np.random.default_rng(0)
for cond, free in (("log phase", C.FREE_MRNA_LOG), ("stress (polysome run-off)", 0.9)):
    cyt = C.Mixture(F, "cytosol", eps, rna_copies=C.MRNA_PER_CELL * free, rna_nt=C.MRNA_NT)
    g_on = cyt.onset(R)
    g_use = min(g, 0.95 * g_on) if g_on else g
    en = C.rna_response(cyt, g_use, R)
    u = cyt.uniprot
    y_pb = np.array([x in pb for x in u]); y_sg = np.array([x in sg and x not in pb for x in u]); y_gb = np.array([x in gb for x in u])
    res = dict(free_mRNA=free, onset=g_on, unstable_at_nucleus_scale=bool(g_on is not None and g_on < g), n_species=len(u),
               n_pb=int(y_pb.sum()), n_sg_only=int(y_sg.sum()), n_gbody=int(y_gb.sum()),
               auroc_pb=C.auroc(en, y_pb), auroc_sg_only=C.auroc(en, y_sg), auroc_gbody=C.auroc(en, y_gb),
               top=[cyt.F.gene.iloc[i] for i in np.argsort(-en)[:25]])
    # nulls: permuted RNA-binding evidence; permuted pool network (protein labels of E shuffled)
    null_rbp, null_pool = [], []
    rbp0, E0 = cyt.rbp.copy(), cyt.E.copy()
    for k in range(50):
        cyt.rbp = rng.permutation(rbp0); null_rbp.append(C.auroc(C.rna_response(cyt, g_use, R), y_pb)); cyt.rbp = rbp0
        p = rng.permutation(len(rbp0)); cyt.E = E0[p][:, p]; null_pool.append(C.auroc(C.rna_response(cyt, g_use, R), y_pb)); cyt.E = E0
    res["null_pb_rbp_permuted"] = [float(np.mean(null_rbp)), float(np.percentile(null_rbp, 95))]
    res["null_pb_pool_permuted"] = [float(np.mean(null_pool)), float(np.percentile(null_pool, 95))]
    # ablations
    for name, rr in (("RNA only", [0, 0, 1]), ("pool only", [1, 0, 0])):
        try:
            res[f"auroc_pb_{name}"] = C.auroc(C.rna_response(cyt, g_use, np.array(rr, float)), y_pb)
        except Exception as ex:
            res[f"auroc_pb_{name}"] = str(ex)
    out[cond] = res
    np.save(os.path.join(C.OUT, f"enrich_{cond.split()[0]}.npy"), np.c_[en])
    json.dump(dict(uniprot=u, gene=list(cyt.F.gene), enrich=en.tolist()), open(os.path.join(C.OUT, f"enrich_{cond.split()[0]}.json"), "w"))
    print(cond, json.dumps({k: v for k, v in res.items() if k != "top"}), "\n top:", res["top"], flush=True)
json.dump(out, open(os.path.join(C.OUT, "test_cytosol.json"), "w"), indent=1)
