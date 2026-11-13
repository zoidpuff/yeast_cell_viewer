"""Diffusion benchmark and per-protein diffusion model (IMPROVEMENTS H2 + H4).

Data: pipeline/data/validation/diffusion_measurements.tsv (literature, every value quoted from the paper; compiled
2026-11-11). Split by paper (validation/splits.paper_split); a compartment with no training paper keeps its test paper
in train (pre-registered rule).

Two parts, because free GFP alone spans 8.7-21 um^2/s between labs and methods:
  A  absolute: inert tracers (free GFP / NLS tracers) per compartment;
  B  relative: within-paper log ratios D_x / D_tracer of the same paper (size and interaction effects, free of lab
     and method scale).
Model: D = D0(R_h) * crowd_c(R_h) * s_c / (1 + beta * b)
  D0       Stokes-Einstein at 30 C with the construct's hydrodynamic radius: the library species that diffuses
           (oligomer or complex), trimmed IDRs restored (features.py rh_full / rh_folded), the fluorescent tag added
           as a bead model (Kirkwood)
  crowd_c  excluded volume + mean-field hydrodynamics from GPU BD of packed boxes (simgpu/sweep.py):
           ln(D/D0) = a_c + k_c R_h, per compartment
  s_c      remaining compartment slowdown (fluid-phase viscosity, missing crowders), fitted on train tracers (A)
  b        binding propensity, fixed before fitting: 1[RNA-interactome capture hit] +
           0.5 log10(1 + abundance-weighted AlphaFold 3 pool partners at PAE <= 12 A / 1e4)
  beta     fitted on train ratios (B)
Nulls: Stokes-Einstein x one constant (size only, no crowding size dependence) and the browser model (mean-field
constant factor, which in ratios equals Stokes-Einstein).
usage: python validation/diffusion.py
"""
import os, sys, re, json, math
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "pipeline")); sys.path.insert(0, os.path.join(HERE, "..", "simgpu"))
from config import DATA, WORK  # noqa: E402
from splits import paper_split  # noqa: E402

VD = f"{DATA}/validation"
KT = 1.380649e-23 * 303.15
ETA = 0.797e-3

# Curation of the constructs (from the papers' construct descriptions, not from the D values):
# id -> (library species ids or 'tracer:<n>xGFP' / 'tracer:Halo', tag, role) ; role: 'tracer' (A and ratio reference),
# 'x' (ratio to its paper's tracer), 'abs' (absolute only, no tracer in the paper). Rows not listed are excluded
# (stress, drugs, ex vivo, bound/chromatin fractions, loci, foci, heterogeneous ranges, membranes: see EXCLUDED).
AFTER_FIT = {"D095", "D096", "D097"}          # rows added after the fit (2026-11-12): test whatever their paper hash
TAG = {"GFP": 238, "Halo": 297, "mCherry": 236, "YFP3": 3 * 238, "none": 0}
CURATED = {
    # cytoplasm, S. cerevisiae
    "D001": ("tracer:1", "none", "tracer"), "D003": ("tracer:1", "none", "tracer"), "D004": ("tracer:3", "none", "x"),
    "D005": ("SRV2", "GFP", "x"), "D006": ("ARC18", "GFP", "x"), "D007": ("EDE1", "GFP", "x"), "D011": ("tracer:2", "none", "tracer"),
    "D010": ("GEM40", "none", "x"),
    "D023": ("CHS6", "GFP", "abs"), "D024": ("BUD7", "GFP", "abs"), "D025": ("CHS5_BCH1_het", "GFP", "abs"),
    "D026": ("BCH2_CHS5_het", "GFP", "abs"), "D027": ("CHS5", "GFP", "abs"),
    "D029": ("tracer:1", "none", "tracer"), "D030": ("SUP35_NM", "GFP", "x"), "D031": ("SUP35", "GFP", "x"),
    "D094": ("tracer:1", "none", "tracer"),        # added 2026-11-11 (Meinema 2013; test paper, after the fit)
    "D032": ("RNQ1", "GFP", "x"), "D034": ("HSP104_oligo", "GFP", "x"), "D036": ("SIS1", "mCherry", "x"),
    # S. pombe cytoplasm and nucleus
    "D092": ("tracer:1", "none", "tracer"), "D093": ("tracer:1", "none", "tracer"),
    "D086": ("GEM40", "none", "abs"), "D085": ("GEM40", "none", "abs"),
    # budding-yeast 40 nm GEMs, values read from figures (added 2026-11-12, IMPROVEMENTS P7; after the fit: test)
    "D095": ("GEM40", "none", "abs"), "D096": ("GEM40", "none", "abs"), "D097": ("GEM40", "none", "abs"),
    # nucleoplasm, S. cerevisiae (free fractions)
    "D062": ("tracer:Halo", "none", "tracer"), "D059": ("RAD52", "Halo", "x"), "D061": ("RAD52", "Halo", "x"),
    "D051": ("HTB1_HTZ1_het", "Halo", "abs"), "D053": ("HTB1_HTZ1_het", "Halo", "abs"), "D056": ("HTB1_HTZ1_het", "Halo", "abs"),   # free H2B travels as an H2A-H2B dimer; the library has no free H2A-H2B species (build R), H2B-H2A.Z has the same fold and size
    "D057": ("TAF1_TAF7_het", "Halo", "abs"), "D066": ("GBP2", "GFP", "abs"), "D045": ("MIG1", "GFP", "abs"),
    "D002": ("tracer:1", "none", "tracer"),
    # ER lumen
    "D072": ("tracer:1", "none", "tracer"), "D073": ("tracer:2", "none", "x"),
}
EXCLUDED_NOTE = ("stress / drugs / genetic perturbations, ex vivo mitochondria, chromatin-bound or slow components, "
                 "loci and repair foci, ranges over heterogeneous factor sets, mRNPs (no mRNA bodies in the model), "
                 "membranes (separate: Saffman-Delbrueck vs measured, reported only)")
COMP = {"cytoplasm": "cyto", "nucleoplasm": "nuc", "ER lumen": "er"}


def d0(rh_nm):
    return KT / (6 * math.pi * ETA * rh_nm * 1e-9) * 1e12            # um^2/s


def _rk(x, nsub=700):
    n = len(x); m = min(n, nsub)
    y = x[np.floor(np.arange(m) * (n / m)).astype(int)]
    d = np.sqrt(((y[:, None] - y[None]) ** 2).sum(-1)); iu = np.triu_indices(m, 1)
    return len(iu[0]) / (1.0 / np.maximum(d[iu], 0.38)).sum()


def _ball(nres, rng):
    """compact globule of nres beads (protein density ~ 1 residue per 0.13 nm^3), centred at 0"""
    R = (nres * 0.13 * 3 / (4 * math.pi)) ** (1 / 3)
    p = rng.uniform(-R, R, (nres * 3, 3)); p = p[(p ** 2).sum(1) <= R * R][:nres]
    return p, R


def construct_rh(beads, tag_res, rng=np.random.default_rng(0), n_dir=6):
    """R_h (nm) of a bead model with a fluorescent-protein tag of tag_res residues attached at its surface"""
    if tag_res == 0:
        return 1.2 * _rk(beads) + 0.3
    t, Rt = _ball(tag_res, rng)
    c = beads.mean(0); out = []
    for _ in range(n_dir):
        u = rng.normal(size=3); u /= np.linalg.norm(u)
        ext = (beads - c) @ u
        pos = c + u * (ext.max() + Rt + 0.5)
        out.append(1.2 * _rk(np.vstack([beads, t + pos])) + 0.3)
    return float(np.mean(out))


def species_beads(lib, sid, rng):
    if sid.startswith("tracer:"):
        k = sid.split(":")[1]
        if k == "Halo":
            return _ball(297, rng)[0]
        n = int(k); parts = []
        for i in range(n):                                   # tandem GFPs, ~1 nm linkers
            b, R = _ball(238, rng); parts.append(b + np.array([i * (2 * R + 1.0), 0, 0]))
        return np.vstack(parts)
    if sid == "GEM40":
        return None
    if sid == "SUP35_NM":                                    # Sup35 N+M domains (~250 aa, disordered): random coil
        import features as FT
        return FT.grow_coil(np.zeros(3), np.array([1.0, 0, 0]), 250, np.zeros((0, 3)), rng)
    t = lib.by_id.get(sid)
    if t is None:                                            # e.g. Sup35 only exists in the eRF1-eRF3 complex species
        t = next(i for i, T in enumerate(lib.types) if sid in [g.strip().upper() for g in re.split(r"[;,]", T["genes"])]
                 and T["id"].startswith("CPX"))
    return lib.bead_xyz(t).astype(np.float64)


def crowd_fits():
    """per compartment (a, k): ln(D/D0) = a + k R_h from the BD sweeps (mean-field hydrodynamics)"""
    out = {}
    for c, job in (("cyto", "cyto120_mf"), ("nuc", "nuc100_mf")):
        d = json.load(open(f"{WORK}/sweeps/{job}.json"))
        rh = np.array(d["rh"]); r = np.array(d["D"]) / np.array(d["D0"])
        A = np.vstack([np.ones_like(rh), rh]).T
        out[c] = tuple(np.linalg.lstsq(A, np.log(np.clip(r, 0.02, None)), rcond=None)[0])
    # ER lumen: 44 mg/mL (MODEL.md 10.1) -- no sweep; mean-field factor at phi ~ 0.05 and no size dependence
    out["er"] = (math.log(1 - 1.832 * 0.05), 0.0)
    return out


def binding_scores():
    F = pd.read_parquet(f"{WORK}/features/protein_features.parquet")
    rbp = set()
    # RNA-interactome capture: Beckmann et al. 2015 (Nat Commun 6:10127), Matia-Gonzalez et al. 2015 (NSMB 22:1027),
    # Mitchell et al. 2013 (NSMB 20:127); systematic ORF names anywhere in their tables
    orf_rx = r"\bY[A-P][LR]\d{3}[WC](?:-[A-Z])?\b"
    for f in ("beckmann2015_rbp.tsv", "mitchell2013_rbp.tsv"):
        rbp |= set(m.upper() for m in re.findall(orf_rx, open(f"{VD}/{f}").read()))
    mt = pd.read_csv(f"{VD}/matia2015_rbp.tsv", sep="\t", skiprows=1, dtype=str)       # all quantified; flag column
    rbp |= set(mt.loc[mt["mRBPome"].str.strip() == "Y", "ORF (SGD)"].str.strip().str.upper())
    F["rbp"] = F.orf.str.upper().isin(rbp)
    F["partner_copies"] = F.get("pool_deg_le8_copies", 0) + F.get("pool_deg_8_12_copies", 0)
    F["b"] = F.rbp.astype(float) + 0.5 * np.log10(1 + F.partner_copies / 1e4)
    return F.set_index("gene"), rbp


def table():
    import boxes as bx
    lib = bx.Library()
    F, _ = binding_scores()
    d = pd.read_csv(f"{VD}/diffusion_measurements.tsv", sep="\t", dtype=str)
    d = d[d.id.isin(CURATED)].copy()
    d["D"] = d.D_um2_s.str.extract(r"([\d.]+)").astype(float)
    d["comp"] = d.compartment.map(lambda s: COMP.get(s.split(" (")[0]))
    d["split"] = d.doi.map(paper_split)
    d.loc[d.id.isin(AFTER_FIT), "split"] = "test"      # first seen after the model was fitted: never training data
    # compartment rule: a compartment with no training paper keeps its test papers in train
    for c in d.comp.unique():
        if not ((d.comp == c) & (d.split == "train")).any():
            d.loc[d.comp == c, "split"] = "train"
    rng = np.random.default_rng(1)
    rh, b = [], []
    for i, r in d.iterrows():
        sid, tag, role = CURATED[r.id]
        if sid == "GEM40":
            rh.append(20.0); b.append(0.0); continue
        beads = species_beads(lib, sid, rng)
        x = construct_rh(beads, TAG[tag], rng)
        genes = [] if sid.startswith("tracer") else [g for g in re.split(r"[_;]", sid) if g in F.index]
        if genes and len(genes) == 1 and sid in F.index and F.loc[sid, "rh_folded"] > 0:
            x *= float(F.loc[sid, "rh_full"] / F.loc[sid, "rh_folded"])    # trimmed IDRs restored (monomers)
        rh.append(x)
        bb = [float(F.loc[g, "b"]) if not isinstance(F.loc[g, "b"], pd.Series) else float(F.loc[g, "b"].max()) for g in genes]
        b.append(max(bb) if bb else 0.0)
    d["rh"] = rh; d["b"] = b
    d["role"] = d.id.map(lambda k: CURATED[k][2])
    return d


def predict(d, cf, s, beta, crowd=True):
    out = []
    for _, r in d.iterrows():
        a, k = cf[r.comp]
        f = math.exp(a + k * r.rh) if crowd else math.exp(a + k * 3.0)
        out.append(d0(r.rh) * f * s.get(r.comp, 1.0) / (1 + beta * r.b))
    return np.array(out)


def ratios(d, pred):
    """within-paper log10 ratios to the paper's tracer in the same compartment"""
    rows = []
    for (doi, comp), g in d.groupby(["doi", "comp"]):
        tr = g[g.role == "tracer"]
        if tr.empty:
            continue
        t0 = tr.index[0]
        for i in g.index:
            if g.loc[i, "role"] == "x":
                rows.append(dict(id=g.loc[i, "id"], split=g.loc[i, "split"], meas=math.log10(g.loc[i, "D"] / g.loc[t0, "D"]),
                                 pred=math.log10(pred[d.index.get_loc(i)] / pred[d.index.get_loc(t0)])))
    return pd.DataFrame(rows)


def fit_and_score():
    d = table()
    cf = crowd_fits()
    out = {}
    for name, crowd, use_beta in (("Stokes-Einstein x constant (null)", False, False), ("+ crowding(R_h) from BD", True, False),
                                  ("+ crowding + binding", True, True)):
        # s_c: geometric-mean match of the train tracers per compartment
        base = predict(d, cf, {}, 0.0, crowd)
        s = {}
        for c in d.comp.unique():
            m = ((d.comp == c) & (d.split == "train") & (d.role == "tracer")).values
            if m.any():
                s[c] = float(np.exp(np.mean(np.log(d.D.values[m]) - np.log(base[m]))))
        beta = 0.0
        if use_beta:
            best = None
            for bt in np.concatenate([[0], np.logspace(-2, 2, 81)]):
                R = ratios(d, predict(d, cf, s, bt, crowd)); R = R[R.split == "train"]
                e = float(np.sqrt(np.mean((R.pred - R.meas) ** 2)))
                if best is None or e < best[0]:
                    best = (e, bt)
            beta = best[1]
        pred = predict(d, cf, s, beta, crowd)
        d[f"pred::{name}"] = pred
        R = ratios(d, pred)
        res = dict(s=s, beta=beta)
        for sp in ("train", "test"):
            m = (d.split == sp).values
            e_abs = np.log10(pred[m]) - np.log10(d.D.values[m])
            r = R[R.split == sp]
            res[sp] = dict(n_abs=int(m.sum()), rmse_log10_abs=float(np.sqrt(np.mean(e_abs ** 2))),
                           n_ratio=len(r), rmse_log10_ratio=float(np.sqrt(np.mean((r.pred - r.meas) ** 2))) if len(r) else None,
                           spearman_ratio=float(pd.Series(r.pred.values).corr(pd.Series(r.meas.values), method="spearman")) if len(r) > 2 else None)
        out[name] = res
    return d, out


if __name__ == "__main__":
    d, out = fit_and_score()
    os.makedirs(f"{WORK}/validation", exist_ok=True)
    d.to_csv(f"{WORK}/validation/diffusion_benchmark.tsv", sep="\t", index=False)
    json.dump(out, open(f"{WORK}/validation/diffusion_benchmark.json", "w"), indent=1)
    cols = ["id", "split", "comp", "molecule", "D", "rh", "b"] + [c for c in d.columns if c.startswith("pred::")]
    print(d[cols].to_string(float_format=lambda v: f"{v:.3g}"))
    print(json.dumps(out, indent=1))
