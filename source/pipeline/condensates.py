"""Emergent condensates from generic interaction data (IMPROVEMENTS H5).

No condensate label enters the energies. Per pair of protein species i, j and per protein and RNA:
  pool    eps_ij = lam_pp * max(0, ln LR(PAE_ij)): LR is the likelihood ratio of a pair being a known physical
          interaction (IntAct) in its AlphaFold 3 pool PAE bin against all scored pairs -- includes the weak pairs
          (PAE 8-16 A) that the organisation layer leaves out
  stick   eps_ij += lam_st * s_i s_j: IDR stickers (aromatic + Arg residues in IDRs, features.py), s = n / n_95%
          capped at 1 (unnormalised counts let a few very long IDRs form single-species modes)
  (pairs of subunits of one Complex Portal complex are excluded: they move as one particle)
  RNA     eps_iR = lam_rna * (number of RNA-interactome capture screens that found i, 0-3) per binding site;
          an mRNA carries one site per 30 nt
Physics: multicomponent second-virial free energy at the compartment's real concentrations c_i (molecules/nm^3):
  f/kT = sum c_i (ln c_i - 1) + 1/2 sum B_ij c_i c_j,   B_ij = (2 pi / 3)(R_i + R_j)^3 - v_b (exp(eps_ij) - 1) n_ij
The mixture is unstable when A = I + diag(sqrt c) B diag(sqrt c) has a negative eigenvalue; the eigenvector u of the
lowest eigenvalue is the composition of the emerging condensate: enrichment of species i is u_i / sqrt(c_i).
The overall scale at which the lowest eigenvalue first reaches zero (onset) is found by bisection.
Protocol (pre-registered): the ratios lam_pp : lam_st : lam_rna are chosen ONLY on the nucleus (nucleoplasm and
nucleolus pooled, so the model does not know which proteins are nucleolar), scored by CD-CODE nucleolus membership;
then frozen. The cytosol is scored once against CD-CODE P-bodies (held out) and stress granules (negative control),
in log phase and after an in-silico stress (all mRNA released from polysomes).
"""
import os, re, json, math
import numpy as np
import pandas as pd

from config import DATA, WORK

VD = f"{DATA}/validation"
OUT = f"{WORK}/condensates"
V_B = 1.0                  # nm^3, bond volume of a pair contact (assumption; only the product with exp(eps) matters)
NT_PER_SITE = 30           # nt per RNA binding site
import volumes as _vol
V = {"cytosol": _vol.measured()["cytosol"] * 1e9,
     "nucleus": (_vol.measured()["nucleoplasm"] + _vol.measured()["nucleolus"]) * 1e9}      # nm^3, built geometry
MRNA_PER_CELL = 30000      # MODEL.md 2.2 (assumption)
MRNA_NT = 1500
FREE_MRNA_LOG = 0.2        # non-polysomal fraction of mRNA in log phase (assumption, unchecked)
RRNA_NT = 6600             # 35S pre-rRNA ~6.6 knt
PRE_RIBO = 27800           # pre-ribosomes in the nucleus (MODEL.md 2.2)


def pool_lr():
    """likelihood ratio of IntAct physical interaction per PAE bin; returns (bin edges, ln LR per bin) and the pair table"""
    ps = pd.read_parquet(f"{WORK}/interactions/pair_summary.parquet", columns=["p1", "p2", "pae_min"])
    ps["p1"] = ps.p1.str.upper(); ps["p2"] = ps.p2.str.upper()
    ps = ps[ps.p1 != ps.p2]
    ia = pd.read_csv(f"{DATA}/intact_pairs.tsv", sep="\t")
    known = set(zip(ia.acc_a, ia.acc_b)) | set(zip(ia.acc_b, ia.acc_a))
    ps["known"] = [(a, b) in known for a, b in zip(ps.p1, ps.p2)]
    edges = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 25, 40])
    ps["bin"] = np.digitize(ps.pae_min, edges) - 1
    base = ps.known.mean()
    lr = ps.groupby("bin").known.mean() / base
    lnlr = np.log(np.clip(lr.reindex(range(len(edges) - 1)).fillna(1.0).values, 1e-3, None))
    return edges, lnlr, ps


def species_table():
    """per protein (UniProt): copies per compartment group, radius, stickers, RNA-binding evidence"""
    F = pd.read_parquet(f"{WORK}/features/protein_features.parquet")
    rbp = {}
    orf_rx = r"\bY[A-P][LR]\d{3}[WC](?:-[A-Z])?\b"
    for f in ("beckmann2015_rbp.tsv", "mitchell2013_rbp.tsv"):
        for m in set(re.findall(orf_rx, open(f"{VD}/{f}").read())):
            rbp[m.upper()] = rbp.get(m.upper(), 0) + 1
    mt = pd.read_csv(f"{VD}/matia2015_rbp.tsv", sep="\t", skiprows=1, dtype=str)
    for m in set(mt.loc[mt["mRBPome"].str.strip() == "Y", "ORF (SGD)"].str.strip().str.upper()):
        rbp[m] = rbp.get(m, 0) + 1
    w = F.comp_weights.map(json.loads)
    F["n_cyto"] = F.copies_per_cell * w.map(lambda d: d.get("cytosol", 0.0))
    F["n_nuc"] = F.copies_per_cell * w.map(lambda d: d.get("nucleoplasm", 0.0) + d.get("nucleolus", 0.0))
    F["R"] = F.rh_full.fillna(F.rh_folded).fillna(3.0)
    st = F.idr_sticker.fillna(0)
    F["stick"] = np.clip(st / max(1.0, st.quantile(0.95)), 0, 1)       # 0-1, saturating at the 95th percentile
    F["rbp"] = F.orf.str.upper().map(rbp).fillna(0).astype(float)
    return F


class Mixture:
    def __init__(self, F, comp, eps_pp, rna_copies, rna_nt, rna_R=10.0, min_copies=50):
        n = F["n_cyto" if comp == "cytosol" else "n_nuc"].values
        keep = n >= min_copies
        self.F = F[keep]; self.c = n[keep] / V[comp]
        self.R = self.F.R.values; self.stick = self.F.stick.values; self.rbp = self.F.rbp.values
        idx = {u: k for k, u in enumerate(self.F.index)}
        self.E = np.zeros((len(self.c), len(self.c)))
        for (a, b), e in eps_pp.items():
            ia, ib = idx.get(a), idx.get(b)
            if ia is not None and ib is not None:
                self.E[ia, ib] = self.E[ib, ia] = max(self.E[ia, ib], e)
        self.c_rna = rna_copies / V[comp]; self.sites = rna_nt / NT_PER_SITE; self.R_rna = rna_R
        self.uniprot = list(self.F.index)

    def matrix(self, lam_pp, lam_st, lam_rna):
        R = self.R; Rr = self.R_rna
        n = len(R)
        eps = lam_pp * self.E + lam_st * np.outer(self.stick, self.stick)
        np.clip(eps, 0, 30, out=eps)
        B = (2 * np.pi / 3) * (R[:, None] + R[None, :]) ** 3 - V_B * np.expm1(eps)
        # RNA as one extra component: sites x per-site attraction
        e_r = np.clip(lam_rna * self.rbp, 0, 30)
        b_r = (2 * np.pi / 3) * (R + Rr) ** 3 - V_B * np.expm1(e_r) * self.sites
        b_rr = (2 * np.pi / 3) * (2 * Rr) ** 3
        Bf = np.zeros((n + 1, n + 1)); Bf[:n, :n] = B; Bf[:n, n] = Bf[n, :n] = b_r; Bf[n, n] = b_rr
        c = np.concatenate([self.c, [self.c_rna]])
        s = np.sqrt(c)
        A = np.eye(n + 1) + s[:, None] * Bf * s[None, :]
        return A, c

    def lowest(self, lam_pp, lam_st, lam_rna):
        A, c = self.matrix(lam_pp, lam_st, lam_rna)
        try:
            import cupy as cp
            w, U = cp.linalg.eigh(cp.asarray(A)); w = cp.asnumpy(w[:3]); U = cp.asnumpy(U[:, :3])
        except Exception:
            w, U = np.linalg.eigh(A); w = w[:3]; U = U[:, :3]
        u = U[:, 0]
        if u[np.argmax(np.abs(u))] < 0:
            u = -u
        enrich = u / np.sqrt(c)               # fold-change direction per species (last entry: RNA)
        return w, u, enrich

    def onset(self, ratios, lo=0.0, hi=40.0, tol=1e-3):
        """smallest overall scale g at which the lowest eigenvalue of A(g * ratios) reaches 0 (None if > hi)"""
        f = lambda g: self.lowest(*(g * np.asarray(ratios)))[0][0]
        if f(hi) > 0:
            return None
        while hi - lo > tol:
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
        return hi


def auroc(score, labels):
    s = np.asarray(score, float); y = np.asarray(labels, bool)
    if y.sum() == 0 or (~y).sum() == 0:
        return np.nan
    order = np.argsort(s); ranks = np.empty(len(s)); ranks[order] = np.arange(1, len(s) + 1)
    return float((ranks[y].sum() - y.sum() * (y.sum() + 1) / 2) / (y.sum() * (~y).sum()))


def cdcode_members(name):
    d = pd.read_csv(f"{VD}/cdcode_yeast_condensate_members.tsv", sep="\t", dtype=str)
    return set(d.loc[d.condensate_name == name, "uniprot_id"])


def complex_pairs():
    """protein pairs that are subunits of one curated stable complex (Complex Portal, taxon 559292): they diffuse as
    one particle, so their interface does not drive condensation"""
    cp = pd.read_csv(f"{DATA}/complexportal_559292.tsv", sep="\t", dtype=str)
    out = set()
    for mem in cp.iloc[:, 4].dropna():
        accs = sorted(set(re.findall(r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9])", mem)))
        for i in range(len(accs)):
            for j in range(i + 1, len(accs)):
                out.add((accs[i], accs[j])); out.add((accs[j], accs[i]))
    return out


def build_eps(edges, lnlr, ps):
    """pair -> weight ln LR (>0 only), pairs inside one stable complex excluded"""
    w = np.clip(lnlr, 0, None)
    b = np.clip(np.digitize(ps.pae_min.values, edges) - 1, 0, len(w) - 1)
    m = w[b] > 0
    cx = complex_pairs()
    return {(a, c): float(x) for a, c, x in zip(ps.p1.values[m], ps.p2.values[m], w[b][m]) if (a, c) not in cx}


def rna_response(mix, g, ratios):
    """enrichment of each protein in an RNA-rich phase: linear response delta c_i / c_i to a rise in the RNA chemical
    potential, ~ (A^-1)_iR sqrt(c_R / c_i), at overall scale g (below onset, where A is positive definite)"""
    A, c = mix.matrix(*(g * np.asarray(ratios)))
    try:
        import cupy as cp
        x = cp.asnumpy(cp.linalg.solve(cp.asarray(A), cp.asarray(np.eye(len(c))[:, -1])))
    except Exception:
        x = np.linalg.solve(A, np.eye(len(c))[:, -1])
    return x[:-1] * np.sqrt(c[-1] / c[:-1])
