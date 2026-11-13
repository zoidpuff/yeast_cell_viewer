"""Localization benchmark (IMPROVEMENTS H1): the model's compartment shares per protein against independent
microscopy and spatial proteomics, on the held-out split of validation/splits.py.

Model side, per protein (systematic ORF):
  assigned  the inventory's compartment weights (inventory.build), before complexes and packing;
  placed    where its copies ended up: placed copies per region of every library type that contains the gene
            (types.json 'regions'), so packing shortfalls (nucleolus 77%, Golgi 81%, matrix 40%) show up here only.
Both are reduced to coarse classes (CLASSES). Reference side: see load_* functions; each returns ORF -> set of classes.
Scores: top-1 (the model's main class is a reference class), per-class F1, and (2026-11-12, after the Jurkat model) a
share agreement per protein, averaged plainly and weighted by protein mass (copies x mass).
usage: python validation/localization.py [--split dev|test|all] [--out DIR]
"""
import os, sys, json, re, argparse
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
from config import WEB_OUT, DATA, WORK  # noqa: E402
from splits import protein_split  # noqa: E402

CLASSES = ["cytosol", "nucleus", "nucleolus", "ER", "mitochondrion", "golgi_endosome", "vacuole", "periphery",
           "peroxisome", "lipid_droplet", "secreted"]
# inventory compartments -> classes
COMP_CLASS = {"cytosol": "cytosol", "nucleoplasm": "nucleus", "nucleolus": "nucleolus", "ER": "ER",
              "nuclear_envelope": "ER", "mito_IM": "mitochondrion", "mito_OM": "mitochondrion", "mito_IMS": "mitochondrion",
              "mito_matrix": "mitochondrion", "mitochondrion": "mitochondrion", "golgi_endosome": "golgi_endosome",
              "vacuole": "vacuole", "plasma_membrane": "periphery", "peroxisome": "peroxisome",
              "lipid_droplet": "lipid_droplet", "cell_wall": "secreted"}
# packing regions (tiles3d.REGION_NAMES + nucleosome 22, NPC 23) -> classes
REGION_CLASS = ["cytosol", "nucleus", "nucleolus", "ER", "vacuole", "mitochondrion", "mitochondrion", "golgi_endosome",
                "peroxisome", "secreted", "periphery", "ER", "ER", "lipid_droplet",
                "periphery", "ER", "ER", "vacuole", "mitochondrion", "mitochondrion", "golgi_endosome", "peroxisome",
                "nucleus", "ER", "vacuole", "golgi_endosome", "peroxisome"]


def gene_to_orf():
    sgd = pd.read_csv(f"{DATA}/sgd/SGD_features.tab", sep="\t", header=None)
    m = {}
    for orf, std, typ in zip(sgd[3], sgd[4], sgd[1]):
        if typ == "ORF" and isinstance(orf, str):
            m[orf.upper()] = orf
            if isinstance(std, str):
                m[std.upper()] = orf
    return m


def model_shares():
    """DataFrame index ORF, columns assigned_<class>, placed_<class> (shares summing to 1 where defined), copies"""
    import inventory as I
    inv, _ = I.build()
    g2o = gene_to_orf()
    rows = {}
    for orf, cw, n, md in zip(inv.orf, inv.comp_weights, inv.copies_per_cell, inv.mass_da):
        w = json.loads(cw); s = {}
        for c, v in w.items():
            k = COMP_CLASS.get(c)
            if k:
                s[k] = s.get(k, 0.0) + v
        z = sum(s.values()) or 1.0
        rows[orf] = {f"assigned_{k}": v / z for k, v in s.items()}
        rows[orf]["copies"] = n; rows[orf]["mass_da"] = md
    types = json.load(open(f"{WEB_OUT}/library/types.json"))["types"]
    placed = {}
    for t in types:
        genes = [g.strip().upper() for g in re.split(r"[;,]", t["genes"]) if g.strip()]
        for g in genes:
            orf = g2o.get(g)
            if not orf:
                continue
            d = placed.setdefault(orf, {})
            for r, c in t["regions"].items():
                k = REGION_CLASS[int(r)] if int(r) < len(REGION_CLASS) else None
                if k:
                    d[k] = d.get(k, 0.0) + c
    for orf, d in placed.items():
        z = sum(d.values())
        if z > 0:
            rows.setdefault(orf, {}).update({f"placed_{k}": v / z for k, v in d.items()})
            rows[orf]["placed_copies"] = z
    df = pd.DataFrame.from_dict(rows, orient="index").fillna(0.0)
    df.index.name = "orf"
    return df


# Volumes (um^3) that turn shares into concentrations, as microscopy sees them (MODEL.md 3.1 table; the plasma
# membrane as a 0.25 um diffraction-limited rim of its 55 um^2; ER = cortical + perinuclear lumen). Added after the
# first look at dev (2026-11-11): microscopy calls where a protein is concentrated, not where most copies are; a
# protein at equal concentration in nucleus and cytosol has 92% of its copies in the cytosol.
import volumes as _vol  # noqa: E402
_m = _vol.measured()          # the built geometry (until 2026-11-12 a copy of the MODEL.md table, with stale mitochondria)
CLASS_VOL = {"cytosol": _m["cytosol"], "nucleus": _m["nucleoplasm"], "nucleolus": _m["nucleolus"], "ER": _m["ER"],
             "mitochondrion": _m["mitochondrion"], "golgi_endosome": _m["golgi_endosome"], "vacuole": _m["vacuole"],
             "periphery": _m["area"]["PM"] * 0.25, "peroxisome": _m["peroxisome"],
             "lipid_droplet": _m["lipid_droplet"], "secreted": _m["cell_wall"]}
CONC_CALL = 0.5           # classes within this factor of the highest concentration count as model calls


def conc_view(df, prefix):
    """shares -> concentrations normalised to the maximum, per protein (same column names)"""
    cols = [f"{prefix}_{c}" for c in CLASSES if f"{prefix}_{c}" in df.columns]
    c = df[cols].copy()
    for col in cols:
        c[col] = c[col] / CLASS_VOL[col[len(prefix) + 1:]]
    mx = c.max(axis=1).replace(0, np.nan)
    return c.div(mx, axis=0).fillna(0.0)


def top_class(df, prefix):
    cols = [f"{prefix}_{c}" for c in CLASSES if f"{prefix}_{c}" in df.columns]
    sub = df[cols]
    top = sub.idxmax(axis=1).str.replace(prefix + "_", "", regex=False)
    top[sub.sum(axis=1) == 0] = None
    return top


def score(model_df, ref, prefix, thresh=0.2, name="", mode="share"):
    """compare model shares with reference class sets (ORF -> set). mode 'share': top class = most copies, calls =
    share >= thresh; mode 'conc': top = highest concentration, calls = within CONC_CALL of it.
    Returns dict of metrics and per-protein table"""
    raw = model_df
    if mode == "conc":
        model_df = conc_view(model_df, prefix); thresh = CONC_CALL
    orfs = [o for o in ref if o in model_df.index and model_df.loc[o, [c for c in model_df.columns if c.startswith(prefix + "_")]].sum() > 0]
    top = top_class(model_df.loc[orfs], prefix)
    hit = np.array([top[o] in ref[o] for o in orfs])
    # multi-label: model classes with share >= thresh
    tp = {c: 0 for c in CLASSES}; fp = dict(tp); fn = dict(tp)
    per = []
    for o in orfs:
        m = {c for c in CLASSES if model_df.loc[o].get(f"{prefix}_{c}", 0.0) >= thresh}
        r = ref[o]
        for c in CLASSES:
            if c in m and c in r: tp[c] += 1
            elif c in m: fp[c] += 1
            elif c in r: fn[c] += 1
        # share agreement (Jurkat model's Tier-1 score, 2026-11-12): 1 - 0.5 x model share outside the reference
        # classes - 0.5 x fraction of reference classes holding < 5% of the model copies; on shares, in both modes
        sh = np.array([raw.loc[o].get(f"{prefix}_{c}", 0.0) for c in CLASSES]); sh = sh / max(sh.sum(), 1e-12)
        inref = np.array([c in r for c in CLASSES])
        agree = 1 - 0.5 * sh[~inref].sum() - 0.5 * (sh[inref] < 0.05).mean() if inref.any() else np.nan
        per.append(dict(orf=o, top=top[o], ref="|".join(sorted(r)), model="|".join(sorted(m)), hit=top[o] in r,
                        jaccard=len(m & r) / max(1, len(m | r)), agree=agree,
                        mass=float(raw.loc[o].get("copies", 0.0) * raw.loc[o].get("mass_da", 0.0))))
    f1 = {}
    for c in CLASSES:
        p = tp[c] / (tp[c] + fp[c]) if tp[c] + fp[c] else np.nan
        rc = tp[c] / (tp[c] + fn[c]) if tp[c] + fn[c] else np.nan
        f1[c] = dict(n_ref=tp[c] + fn[c], precision=p, recall=rc, f1=(2 * p * rc / (p + rc)) if p == p and rc == rc and p + rc > 0 else np.nan)
    rng = np.random.default_rng(0)
    boot = [hit[rng.integers(0, len(hit), len(hit))].mean() for _ in range(1000)] if len(hit) else [np.nan]
    macro = np.nanmean([v["f1"] for v in f1.values() if v["n_ref"] >= 10])
    P = pd.DataFrame(per)
    ag = P.agree.values if len(P) else np.array([np.nan]); w = P.mass.values if len(P) else np.array([1.0])
    ok = np.isfinite(ag)
    return dict(name=name, prefix=prefix, mode=mode, n=len(orfs), top1=float(hit.mean()) if len(hit) else np.nan,
                top1_ci=(float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))),
                agree=float(ag[ok].mean()) if ok.any() else np.nan,
                agree_mass_weighted=float((ag[ok] * w[ok]).sum() / max(w[ok].sum(), 1e-12)) if ok.any() else np.nan,
                macro_f1=float(macro), per_class=f1), P


# ------------------------------------------------------------------ references (independent of CYCLoPs and Huh 2003)
VD = f"{DATA}/validation"
SWAT_MAP = {"cytosol": {"cytosol"}, "nucleus": {"nucleus"}, "nucleolus": {"nucleolus"}, "er": {"ER"},
            "nuclear periphery": {"ER"}, "mitochondria": {"mitochondrion"}, "vacuole": {"vacuole"},
            "vacuole membrane": {"vacuole"}, "cell periphery": {"periphery"}, "bud": {"periphery"},
            "bud neck": {"periphery"}, "punctate": {"golgi_endosome", "peroxisome", "lipid_droplet"}}
LOPIT_MAP = {"Cytosol": {"cytosol"}, "Ribosome": {"cytosol"}, "Nucleus": {"nucleus", "nucleolus"}, "Histone": {"nucleus"},
             "ER": {"ER"}, "Golgi": {"golgi_endosome"}, "Mitochondrion": {"mitochondrion"}, "Vacuole": {"vacuole"},
             "V ATPase": {"vacuole"}, "Plasma membrane": {"periphery"},
             "Proteasome catalytic subunit": {"cytosol", "nucleus"}, "Proteasome regulatory subunit": {"cytosol", "nucleus"}}


def load_swat(column="NATIVEpr-GFP localization"):
    """Weill et al. 2018 (Nat Methods 15:617) SWAT N-terminal GFP calls, Table S1. NATIVEpr: native promoter;
    'NOP1pr-GFP localization ' (with the trailing space): constitutive NOP1 promoter. The 'Known C' Localization'
    column (Huh + LoQAtE) is not used."""
    d = pd.read_csv(f"{VD}/weill2018_S1.tsv", sep="\t", dtype=str)
    out = {}
    for orf, v in zip(d["ORF"], d[column]):
        if not isinstance(v, str):
            continue
        s = set()
        for t in v.split(","):
            s |= SWAT_MAP.get(t.strip().lower(), set())
        if s:
            out[orf.strip().upper()] = s
    return out


def load_hyperlopit():
    """hyperLOPIT of yeast (Nightingale et al. 2019, Curr Opin Chem Biol 48:86; pRolocdata yeast2018): SVM-predicted
    location after the authors' score threshold ('unknown' dropped), plus the marker proteins"""
    d = pd.read_csv(f"{VD}/hyperlopit2018_fData.tsv", sep="\t", dtype=str)
    out = {}
    for desc, pred, mk in zip(d["Description"], d["predicted.location"], d["markers"]):
        m = re.search(r"GN=(\S+)", desc or "")
        if not m:
            continue
        lab = mk if mk != "unknown" else pred
        s = LOPIT_MAP.get(lab)
        if s:
            out[m.group(1).upper()] = s
    g2o = gene_to_orf()
    return {g2o[k]: v for k, v in out.items() if k in g2o}


def load_cyclops_fit():
    """CYCLoPs WT LOC-scores (an INPUT of the model): classes with mean score >= 0.3 -- agreement is a training fit"""
    import inventory as I
    mp = {"cytosol": "cytosol", "nucleoplasm": "nucleus", "nucleolus": "nucleolus", "nuclear_envelope": "ER", "ER": "ER",
          "golgi_endosome": "golgi_endosome", "mitochondrion": "mitochondrion", "peroxisome": "peroxisome",
          "vacuole": "vacuole", "plasma_membrane": "periphery", "periphery": "periphery"}
    out = {}
    for orf, sc in I.cyclops_scores().items():
        s = {mp[I.CYCLOPS_MAP[k]] for k, v in sc.items() if v >= 0.3 and I.CYCLOPS_MAP[k] in mp}
        if s:
            out[orf.upper()] = s
    return out


def run(split="dev", out=None):
    M = model_shares()
    refs = {"SWAT_native": load_swat(), "SWAT_NOP1": load_swat("NOP1pr-GFP localization "),
            "hyperLOPIT": load_hyperlopit(), "CYCLoPs (input: fit, not a test)": load_cyclops_fit()}
    res, tables = [], {}
    for name, ref in refs.items():
        if split != "all":
            ref = {o: v for o, v in ref.items() if protein_split(o) == split}
        for prefix in ("assigned", "placed"):
            for mode in ("share", "conc"):
                r, per = score(M, ref, prefix, name=name, mode=mode)
                res.append(r); tables[(name, prefix, mode)] = per
    if out:
        os.makedirs(out, exist_ok=True)
        json.dump(res, open(f"{out}/localization_{split}.json", "w"), indent=1, default=float)
        for (n, p, m), t in tables.items():
            t.to_csv(f"{out}/per_protein_{split}_{re.sub(r'[^A-Za-z0-9]+', '_', n)}_{p}_{m}.tsv", sep="\t", index=False)
    return res, tables, M


def confusion(per):
    """reference class -> model top class counts (single-class references only)"""
    s = per[~per.ref.str.contains(r"\|")]
    return pd.crosstab(s.ref, s.top)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--split", default="dev"); ap.add_argument("--out", default=f"{WORK}/validation")
    a = ap.parse_args()
    res, tables, _ = run(a.split, a.out)
    for r in res:
        print(f"{r['name']:34s} {r['prefix']:8s} {r['mode']:5s} n={r['n']:5d} top1={r['top1']:.3f} [{r['top1_ci'][0]:.3f},{r['top1_ci'][1]:.3f}] "
              f"agree={r['agree']:.3f} (mass-weighted {r['agree_mass_weighted']:.3f}) macroF1={r['macro_f1']:.3f}")
