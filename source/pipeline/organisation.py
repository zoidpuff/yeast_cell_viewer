"""Protein organisation from association data: a tunable hypothesis layer on top of random packing.

Random placement is itself an assumption (no organisation at all). This module replaces it with a weak, data-driven
one, at two scales. Neither scale is calibrated against measured protein positions; the parameters below are set so
the packing stays random-looking but measurably organised. Compartments are honoured throughout: the layer only
changes where in its compartments a protein goes and next to whom.

Network neighbourhoods (no clustering):
  - Every protein of the STRING network (v12, combined score >= STRING_MIN, edges weighted by score) gets a network
    propagation profile: its personalized PageRank (random walk with restart) vector, with damping DAMPING (the walk
    continues with probability DAMPING and restarts at the seed otherwise; a low damping keeps the profile local).
    The seed's own entry is dropped and the vector renormalised (PROFILE). A library type (monomer or complex) takes
    the mean profile of its subunits.
  - The network similarity of two types is the Pearson correlation r of their profiles: high for proteins that
    share network neighbourhoods, whether or not they are linked directly. STRING_MIN, DAMPING and the profile
    variant were chosen by how well r predicts shared localization and co-complex membership
    (scripts/organisation_sweep.py; MODEL.md section 8.4).

Large scale (across the 200 nm tiles; which proteins go where):
  - The type correlation matrix R is used as the covariance of the types' log-density fields over the cell:
    R ~ E E^T from its FIELD_K leading eigenvectors, with FIELD_K independent smooth Gaussian random fields F_k
    (random Fourier features, correlation length FIELD_L nm, unit variance), and
        log weight of type i in a tile = BETA * sum_k E_ik F_k(tile centre)
    so the log-densities of two types co-vary across the cell by BETA^2 x r_ij: network neighbours are co-enriched
    in the same places, unrelated types vary independently, and types without a network profile are spread
    uniformly. A type's copies are split over tiles in proportion to region volume (or membrane area) x
    exp(log weight). Cell totals and compartments are unchanged, and a crowding correction (wc3d) keeps the total
    macromolecular volume of each tile proportional to its compartment volume.

Meso and small scale (within a tile; who sits next to whom):
  - Contact partners: AlphaFold 3 pool pairs (best chain-pair PAE min-min <= POOL_PAE_MAX Å, paralog artefacts
    removed), proposed in contact. Strength = PAE term (1 at <= 2 Å to 0 at POOL_PAE_MAX) x (1 + r+) / 2, so pairs
    that also share a network neighbourhood pull fully and network-unsupported pairs at half strength; hubs of the
    pool graph are down-weighted.
  - Loose partners: the types most correlated with it (r >= R_LOOSE, not already contact partners), proposed
    2-15 nm away, with strength rising from 0.2 at R_LOOSE to 1 at r = 1.
  - During packing, a copy is placed next to an already placed partner (in its tile or the neighbouring halo) with a
    probability that saturates with the partner's local copy number, like mass action:
        a_p = strength_p x A_kind x n_p / (K_HALF + n_p),   P(attract) = min(P_MAX, sum_p a_p)
    The partner is chosen in proportion to a_p. The copy gets ATTRACT_TRIES proposals near it; if none passes the
    collision and compartment checks it falls back to a uniformly random position, as before.
"""
import os, json
import numpy as np
import pandas as pd
from config import DATA, WORK
from config import seed as _seed

ODIR = f"{WORK}/organisation"
STRING_DIR = f"{DATA}/string"
STRING_F = f"{STRING_DIR}/4932.protein.links.v12.0.txt.gz"
from config import POOL_SCORES                                 # summary_models.parquet of the pool predictions
POOL_PAIRS_F = f"{ODIR}/pool_pairs.parquet"                     # per protein pair: best PAE min-min (cached)

# network neighbourhoods and the large scale
STRING_MIN = 400       # STRING medium confidence; edges weighted by combined score
DAMPING = 0.5          # continuation probability of the walk (restart 1 - DAMPING)
PROFILE = "nodiag"     # seed entry dropped, renormalised
FIELD_K = 128          # leading components of the type correlation matrix that drive the fields
FIELD_L = 500.0        # nm
FIELD_NF = 128         # Fourier features per field
BETA = 0.8             # log-density standard deviation of a type fully described by the leading components
# meso and small scale
POOL_PAE_MAX = 8.0     # Å
MAX_PARTNERS = 16
R_LOOSE = 0.4          # minimum profile correlation of a loose partner
A_CONTACT = 0.6
A_LOOSE = 0.25
K_HALF = 8.0           # copies in the tile (+ halo) at which a partner reaches half its pull
P_MAX = 0.6
ATTRACT_TRIES = 40
LOOSE_GAP = (2.0, 15.0)    # nm, surface gap of loose proposals
SEED = 11
CONTACT, LOOSE = 0, 1


def _uniprot_maps():
    up = pd.read_csv(f"{DATA}/uniprot_yeast.tsv", sep="\t")
    loc2u = {}
    for e, l in zip(up.Entry, up["Gene Names (ordered locus)"].fillna("")):
        for x in str(l).replace(";", " ").split():
            loc2u.setdefault(x, e)
    return loc2u


def string_edges(min_score=STRING_MIN):
    """STRING links as UniProt pairs (a < b), combined score >= min_score"""
    loc2u = _uniprot_maps()
    s = pd.read_csv(STRING_F, sep=" ")
    s = s[s.combined_score >= min_score]
    a = s.protein1.str.split(".", n=1).str[1].map(loc2u); b = s.protein2.str.split(".", n=1).str[1].map(loc2u)
    d = pd.DataFrame(dict(a=a, b=b, score=s.combined_score.values)).dropna()
    d = d[d.a < d.b]
    return d.groupby(["a", "b"], as_index=False).score.max()


def pool_pairs():
    """best (lowest) chain-pair PAE min-min per UniProt pair over all pool and pair predictions; cached"""
    if os.path.exists(POOL_PAIRS_F):
        return pd.read_parquet(POOL_PAIRS_F)
    import pyarrow.parquet as pq, pyarrow as pa, pyarrow.compute as pc
    if not os.path.exists(POOL_SCORES):
        raise SystemExit(f"no pool scores at {POOL_SCORES}: set YC_POOL_SCORES to the pool predictions' summary_models.parquet")
    t = pq.read_table(POOL_SCORES, columns=["af3_id1", "af3_id2", "chain_pair_pae_min_min"])
    a, b = pc.utf8_upper(t["af3_id1"]), pc.utf8_upper(t["af3_id2"])
    lt = pc.less(a, b)
    t = pa.table(dict(a=pc.if_else(lt, a, b), b=pc.if_else(lt, b, a), pae=t["chain_pair_pae_min_min"]))
    g = t.group_by(["a", "b"]).aggregate([("pae", "min")]).to_pandas().rename(columns={"pae_min": "pae"})
    g = g[g.a != g.b]
    os.makedirs(ODIR, exist_ok=True); g.to_parquet(POOL_PAIRS_F)
    return g


def pool_edges(pae_max=POOL_PAE_MAX):
    """confident pool pairs without paralog artefacts: drop (A, B) if a paralog A' of A scores better with B
    (homologous interfaces, e.g. Pol II subunits against Pol III paralogs)"""
    g = pool_pairs()
    g = g[g.pae <= pae_max].copy()
    par = pd.read_csv(f"{DATA}/alliance_yeast_paralog_pairs.csv", usecols=["gene1_uniprot", "gene2_uniprot"]).dropna()
    para = {}
    for x, y in zip(par.gene1_uniprot, par.gene2_uniprot):
        para.setdefault(x, set()).add(y); para.setdefault(y, set()).add(x)
    best = {(a, b): p for a, b, p in zip(g.a, g.b, g.pae)}
    best.update({(b, a): p for (a, b), p in list(best.items())})
    keep = []
    for a, b, p in zip(g.a, g.b, g.pae):
        bad = any(best.get((a2, b), 99) < p for a2 in para.get(a, ())) or any(best.get((a, b2), 99) < p for b2 in para.get(b, ()))
        keep.append(not bad)
    return g[np.array(keep, bool)]


def lib_subunits(exp3, ing):
    """library type -> set of UniProt accessions of its subunits (complexes and heterodimers have no UniProt field
    in the inventory: their gene lists are mapped through UniProt gene names)"""
    up = pd.read_csv(f"{DATA}/uniprot_yeast.tsv", sep="\t")
    g2u = {}
    for e, prim, loc in zip(up.Entry, up["Gene Names (primary)"].fillna(""), up["Gene Names (ordered locus)"].fillna("")):
        for g in (str(prim) + ";" + str(loc)).replace(" ", ";").split(";"):
            if g:
                g2u.setdefault(g, e)
    rows = ing.set_index("ing_id")
    out = {}
    for ing_id, lib in exp3[["ing_id", "lib_id"]].drop_duplicates().itertuples(index=False):
        if ing_id not in rows.index:
            continue
        acc, genes = rows.at[ing_id, "uniprot"], rows.at[ing_id, "genes"]
        if isinstance(acc, str) and acc:
            us = [x for x in acc.split(";") if x]
        elif isinstance(genes, str):
            us = [g2u[g] for g in genes.split(";") if g in g2u]
        else:
            us = []
        if us:
            out.setdefault(lib, set()).update(us)
    return out


def ppr_profiles(edges, damping=DAMPING, profile=PROFILE):
    """personalized PageRank profile of every protein of the network: (nodes, X) with X[i] the profile of seed i
    (rows sum to 1; PROFILE "nodiag": the seed's own entry dropped)"""
    nodes = sorted(set(edges.a) | set(edges.b)); ix = {u: i for i, u in enumerate(nodes)}; N = len(nodes)
    A = np.zeros((N, N))
    a = edges.a.map(ix).values; b = edges.b.map(ix).values; w = edges.score.values / 1000.0
    A[a, b] = w; A[b, a] = w
    W = A / A.sum(0)[None, :]                                  # column-stochastic transition matrix
    X = ((1 - damping) * np.linalg.inv(np.eye(N) - damping * W)).T   # row i: stationary visits of walks restarting at i
    if profile == "nodiag":
        np.fill_diagonal(X, 0.0); X /= X.sum(1, keepdims=True)
    return nodes, X


def type_profiles(lib_ids, subunits, nodes, X):
    """mean profile over each library type's subunits, standardised for correlation: (type indices, Z) with
    Z[k] . Z[l] the Pearson correlation of types tix[k] and tix[l]"""
    ix = {u: i for i, u in enumerate(nodes)}
    tix, rows = [], []
    for i, lib in enumerate(lib_ids):
        us = [ix[u] for u in subunits.get(lib, ()) if u in ix]
        if us:
            tix.append(i); rows.append(X[us].mean(0))
    Z = np.array(rows); Z -= Z.mean(1, keepdims=True); Z /= np.linalg.norm(Z, axis=1, keepdims=True)
    return np.array(tix, np.int64), Z


def field_values(points, n_fields, seed=SEED):
    """unit-variance smooth random fields (correlation length FIELD_L) at points (nm); n x n_fields"""
    seed = _seed(seed); rng = np.random.default_rng(seed)
    K = rng.normal(0, 1.0 / FIELD_L, (n_fields, FIELD_NF, 3)); ph = rng.uniform(0, 2 * np.pi, (n_fields, FIELD_NF))
    return np.sqrt(2.0 / FIELD_NF) * np.cos(np.einsum("pc,mfc->pmf", points, K) + ph[None]).sum(2)


def field_loadings(R, k=FIELD_K):
    """E (types x k) with E E^T the best rank-k approximation of the correlation matrix R (positive part)"""
    lam, U = np.linalg.eigh(R)
    top = np.argsort(lam)[::-1][:k]
    lam = np.clip(lam[top], 0, None)
    return U[:, top] * np.sqrt(lam)[None, :], lam


def tile_logmod(lib_ids, tix, E, centres):
    """per library type with a network profile, the log weight over tiles (centres, nm): BETA * E_i . F"""
    F = field_values(centres, E.shape[1])
    return {lib_ids[i]: BETA * (F @ E[k]) for k, i in enumerate(tix)}


def partners(lib_ids, subunits, pool, tix, R):
    """partner lists per library type for the packer, as CSR arrays over lp['ids'] indices:
    off (n+1), idx, strength (0..1), kind (CONTACT / LOOSE)"""
    lib_ids = [str(x) for x in lib_ids]
    by_u = {}
    for i, lib in enumerate(lib_ids):
        for u in subunits.get(lib, ()):
            by_u.setdefault(u, []).append(i)
    row = np.full(len(lib_ids), -1, np.int64); row[tix] = np.arange(len(tix))
    r_of = lambda i, j: float(R[row[i], row[j]]) if row[i] >= 0 and row[j] >= 0 else 0.0
    shares = lambda i, j: i == j or bool(subunits[lib_ids[i]] & subunits[lib_ids[j]])
    deg = {}
    for a, b in zip(pool.a, pool.b):
        deg[a] = deg.get(a, 0) + 1; deg[b] = deg.get(b, 0) + 1
    P = [dict() for _ in lib_ids]
    # contact: pool pairs, PAE term x network support, hubs of the pool graph pull less
    for a, b, p in zip(pool.a, pool.b, pool.pae):
        if a not in by_u or b not in by_u:
            continue
        s0 = float(np.clip((POOL_PAE_MAX - p) / (POOL_PAE_MAX - 2.0), 0.05, 1.0)) / np.sqrt(max(1.0, max(deg[a], deg[b]) / 20.0))
        for i in by_u[a]:
            for j in by_u[b]:
                if shares(i, j):
                    continue
                s = s0 * (1 + max(0.0, r_of(i, j))) / 2
                for x, y in ((i, j), (j, i)):
                    if P[x].get(y, (CONTACT, -1))[1] < s:
                        P[x][y] = (CONTACT, s)
    # loose: the most correlated network neighbours
    Rm = R.copy(); np.fill_diagonal(Rm, -1)
    for k, i in enumerate(tix):
        js = np.nonzero(Rm[k] >= R_LOOSE)[0]
        js = js[np.argsort(-Rm[k, js])][:2 * MAX_PARTNERS]
        for kk in js:
            j = int(tix[kk])
            if j in P[i] or shares(i, j):
                continue
            P[i][j] = (LOOSE, float(0.2 + 0.8 * (Rm[k, kk] - R_LOOSE) / (1 - R_LOOSE)))
    off = [0]; idx = []; st = []; kd = []
    for i in range(len(lib_ids)):
        top = sorted(P[i].items(), key=lambda kv: (kv[1][0], -kv[1][1]))[:MAX_PARTNERS]   # contact first, then stronger
        for j, (kind, s) in top:
            idx.append(j); st.append(s); kd.append(kind)
        off.append(len(idx))
    return dict(off=np.array(off, np.int64), idx=np.array(idx, np.int64), strength=np.array(st, np.float64),
                kind=np.array(kd, np.int64))


def params():
    return np.array([A_CONTACT, A_LOOSE, K_HALF, P_MAX, ATTRACT_TRIES, LOOSE_GAP[0], LOOSE_GAP[1]], np.float64)


def build(exp3, ing, lib_ids, centres, log=print):
    """everything the packing needs: tile log weights per library type and partner CSR arrays; also saved to ODIR"""
    os.makedirs(ODIR, exist_ok=True)
    lib_ids = [str(x) for x in lib_ids]
    subunits = lib_subunits(exp3, ing)
    st = string_edges()
    nodes, X = ppr_profiles(st)
    tix, Z = type_profiles(lib_ids, subunits, nodes, X)
    R = Z @ Z.T
    E, lam = field_loadings(R)
    logmod = tile_logmod(lib_ids, tix, E, centres)
    pool = pool_edges()
    pr = partners(lib_ids, subunits, pool, tix, R)
    np.savez(f"{ODIR}/partners.npz", ids=np.array(lib_ids), **pr)
    np.savez_compressed(f"{ODIR}/network.npz", ids=np.array(lib_ids), tix=tix, Z=Z.astype(np.float32), E=E.astype(np.float32))
    json.dump(dict(params=dict(STRING_MIN=STRING_MIN, DAMPING=DAMPING, PROFILE=PROFILE, FIELD_K=FIELD_K, FIELD_L=FIELD_L,
                               BETA=BETA, POOL_PAE_MAX=POOL_PAE_MAX, R_LOOSE=R_LOOSE, A_CONTACT=A_CONTACT, A_LOOSE=A_LOOSE,
                               K_HALF=K_HALF, P_MAX=P_MAX)), open(f"{ODIR}/network.json", "w"))
    npart = np.diff(pr["off"])
    off_diag = R[np.triu_indices(len(R), 1)]
    log(f"organisation: {len(st)} STRING links >= {STRING_MIN} over {len(nodes)} proteins; PPR damping {DAMPING}; "
        f"{len(tix)} of {len(lib_ids)} types with a profile (pairs with r >= {R_LOOSE}: {(off_diag >= R_LOOSE).mean():.2%}); "
        f"{FIELD_K} field components hold {lam.sum() / len(R):.0%} of the correlation; {len(pool)} pool pairs (PAE <= {POOL_PAE_MAX:g}); "
        f"{(npart > 0).sum()} types with partners (median {np.median(npart[npart > 0]) if (npart > 0).any() else 0:.0f}; "
        f"{(pr['kind'] == CONTACT).sum()} contact, {(pr['kind'] == LOOSE).sum()} loose links)")
    return logmod, pr


NB_MAX, NB_MIN = 30, 0.2     # viewer: network neighbours listed per type (strongest first, r >= NB_MIN)


def encode_viewer(out_dir, log=print):
    """library/organisation.json for the viewer, keyed by viewer type index: a colour per type from a 2D map of the
    network profiles (t-SNE on 1 - r; nearby colours = similar neighbourhoods) and each type's strongest network
    neighbours (r >= NB_MIN, at most NB_MAX) for highlighting"""
    from sklearn.manifold import TSNE
    import colorsys
    tj = json.load(open(f"{out_dir}/library/types.json"))["types"]
    z = np.load(f"{ODIR}/network.npz")
    ids = [str(x) for x in z["ids"]]; tix = z["tix"]; Z = z["Z"].astype(np.float64)
    vix = {t["id"]: i for i, t in enumerate(tj)}
    keep = [k for k, i in enumerate(tix) if ids[i] in vix]
    Z = Z[keep]; vt = np.array([vix[ids[tix[k]]] for k in keep])
    R = np.clip(Z @ Z.T, -1, 1)
    Y = TSNE(2, metric="precomputed", init="random", perplexity=30, random_state=SEED).fit_transform(np.clip(1 - R, 0, 2))
    Y -= np.median(Y, 0)
    ang = np.arctan2(Y[:, 1], Y[:, 0]); rad = np.linalg.norm(Y, axis=1)
    rq = np.argsort(np.argsort(rad)) / max(1, len(rad) - 1)          # radius as a quantile: even use of saturation
    # hue: the angle's rank weighted by sqrt(copies), so the abundant neighbourhoods (glycolysis, translation) spread
    # over more of the colour wheel; the order around the map, and so similar colours for neighbours, is kept
    wt = np.sqrt(np.array([max(tj[i]["copies"], 1) for i in vt], float))
    o = np.argsort(ang); hue = np.empty(len(ang)); hue[o] = (np.cumsum(wt[o]) - 0.5 * wt[o]) / wt.sum()
    cols = {}
    for k in range(len(vt)):
        r, g, b = colorsys.hls_to_rgb(hue[k], 0.60 - 0.12 * rq[k], 0.35 + 0.5 * rq[k])
        cols[str(int(vt[k]))] = "#%02X%02X%02X" % (round(r * 255), round(g * 255), round(b * 255))
    np.fill_diagonal(R, -1)
    nb = {}
    for k in range(len(vt)):
        js = np.nonzero(R[k] >= NB_MIN)[0]
        js = js[np.argsort(-R[k, js])][:NB_MAX]
        if len(js):
            nb[str(int(vt[k]))] = [[int(vt[j]), round(float(R[k, j]), 2)] for j in js]
    json.dump(dict(colors=cols, neighbours=nb, params=json.load(open(f"{ODIR}/network.json"))["params"]),
              open(f"{out_dir}/library/organisation.json", "w"), separators=(",", ":"))
    log(f"organisation.json: {len(cols)} of {len(tj)} types with a network colour, {len(nb)} with neighbours "
        f"(median {np.median([len(v) for v in nb.values()]) if nb else 0:.0f})")
