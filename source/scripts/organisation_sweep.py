"""Choose the network-propagation parameters of the organisation layer (pipeline/organisation.py).

For STRING thresholds 400 and 700 and PPR damping 0.3-0.85, and four profile variants (raw; seed entry dropped;
seed entry dropped and degree-normalised; log), reports how well the Pearson correlation of two proteins'
personalized PageRank profiles predicts
  - shared localization: same coarse compartment (cytosol, nucleus, mitochondria, ER, ...) for proteins with a
    dominant (>= 60%) location in the inventory; also for pairs that are not directly linked;
  - co-complex membership: Complex Portal co-members against random pairs;
and how local the profiles are (median number of proteins with r > 0.3 and r > 0.5).
Usage (from the bundle root, after s1): .venv/bin/python scripts/organisation_sweep.py     (~20 min)
"""
import sys, os, json, re, numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'pipeline')); import organisation as og, config as C
from sklearn.metrics import roc_auc_score
rng = np.random.default_rng(0)
ing = pd.read_csv('work/results/inventory/ingredients.csv'); p = ing[ing.kind == 'protein']
GROUP = lambda c: ('nucleus' if c in ('nucleoplasm', 'nucleolus') else 'mito' if c.startswith('mito') else
                   'ER' if c.startswith('ER') or c == 'nuclear_envelope' else c)
loc = {}
for u, cw in zip(p.uniprot, p.comp_weights):
    w = json.loads(cw); g = {}
    for k, v in w.items(): g[GROUP(k)] = g.get(GROUP(k), 0) + v
    k, v = max(g.items(), key=lambda kv: kv[1])
    if isinstance(u, str) and v / sum(g.values()) >= 0.6: loc[u] = k
cp = pd.read_csv(f'{C.DATA}/complexportal_559292.tsv', sep='\t')
cplx = []
for s in cp['Identifiers (and stoichiometry) of molecules in complex']:
    m = set(re.findall(r'([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9])', s))
    if 2 <= len(m) <= 60: cplx.append(m)
for thr in (400, 700):
    st = og.string_edges(min_score=thr)
    nodes = sorted(set(st.a) | set(st.b)); ix = {u: i for i, u in enumerate(nodes)}; N = len(nodes)
    A = np.zeros((N, N)); a = st.a.map(ix).values; b = st.b.map(ix).values; w = st.score.values / 1000.0
    A[a, b] = w; A[b, a] = w
    deg = A.sum(0); W = A / deg[None, :]                 # column-stochastic
    adj = A > 0
    # evaluation pairs
    L = [u for u in nodes if u in loc]; Li = np.array([ix[u] for u in L]); Lg = np.array([loc[u] for u in L])
    pi = rng.integers(0, len(L), 400000); pj = rng.integers(0, len(L), 400000); k = pi != pj; pi, pj = Li[pi[k]], Li[pj[k]]
    lab_of = dict(zip(Li, Lg)); yl = np.array([lab_of[x] == lab_of[y] for x, y in zip(pi, pj)])
    na = ~adj[pi, pj]
    pos = [(ix[x], ix[y]) for c in cplx for x in c for y in c if x < y and x in ix and y in ix]
    pos = np.array(pos); ci = rng.integers(0, N, (len(pos) * 5, 2)); ci = ci[ci[:, 0] != ci[:, 1]]
    print(flush=True); print(f'STRING >= {thr}: {N} proteins, {len(st)} links; loc pairs {len(pi)} (same {yl.mean():.3f}); complex pairs {len(pos)}')
    for d in (0.3, 0.5, 0.7, 0.85):
        alpha = 1 - d
        P = alpha * np.linalg.inv(np.eye(N) - d * W)        # column i: PPR from seed i
        for var in ('raw', 'nodiag', 'deg', 'log'):
            X = P.T.copy()                                   # row i: vector of seed i
            if var in ('nodiag', 'deg', 'log'): np.fill_diagonal(X, 0); X /= X.sum(1, keepdims=True)
            if var == 'deg': X = X / np.sqrt(deg)[None, :]
            if var == 'log': X = np.log(X + 1e-6)
            Z = X - X.mean(1, keepdims=True); Z /= np.linalg.norm(Z, axis=1, keepdims=True)
            r = lambda i, j: np.einsum('ij,ij->i', Z[i], Z[j])
            rl = r(pi, pj)
            auc_l = roc_auc_score(yl, rl); auc_ln = roc_auc_score(yl[na], rl[na])
            auc_c = roc_auc_score(np.r_[np.ones(len(pos)), np.zeros(len(ci))], np.r_[r(pos[:, 0], pos[:, 1]), r(ci[:, 0], ci[:, 1])])
            R = Z @ Z.T; np.fill_diagonal(R, 0)
            print(f'  damping {d:.2f} {var:6s}: AUROC localization {auc_l:.3f} (non-adjacent {auc_ln:.3f}), co-complex {auc_c:.3f}; '
                  f'median #partners r>0.3: {np.median((R > 0.3).sum(1)):.0f}, r>0.5: {np.median((R > 0.5).sum(1)):.0f}', flush=True)
