"""Membrane association of soluble proteins from lipid anchors and lipid-binding domains (IMPROVEMENTS H6).

Evidence per ORF (no localization data enters):
  lipid    UniProt lipidation features: N-myristoyl glycine, S-palmitoyl cysteine, S-farnesyl / S-geranylgeranyl
           cysteine (GPI anchors are left to the secretory rules)
  domain   membrane-binding domains (InterPro/Pfam: PH, PX, FYVE, C2, BAR/F-BAR, ENTH/ANTH, GRAM, CRAL-TRIO,
           oxysterol-binding, SMP) or a lipid-array hit (Gallego et al. 2010, Mol Syst Biol 6:430)
Target membrane when the annotation names none: by anchor chemistry and lipid preference --
  geranylgeranyl (Rab/Rho GTPases on vesicles), PI3P binders (FYVE, PX) -> Golgi/endosome; myristoyl, palmitoyl,
  farnesyl, PI(4,5)P2 binders, everything else -> plasma membrane.
inventory.build() moves the fraction BOUND_F[evidence] of a soluble protein's copies to the peripheral shell of that
membrane (MODEL.md 4.4); BOUND_F was fitted on the dev half of the SWAT comparison (validation/tune_anchor.py).
"""
import re
import pandas as pd

from config import DATA

VD = f"{DATA}/validation"
MEM_DOMAINS = re.compile(r"PH|PX|FYVE|C2|BAR|ENTH|ANTH|GRAM|CRAL|Sec14|Oxysterol|SMP")
ENDO_DOMAINS = re.compile(r"FYVE|PX")
MEMBRANE_COMPS = ("plasma_membrane", "golgi_endosome", "vacuole", "ER", "nuclear_envelope", "peroxisome",
                  "mito_OM", "lipid_droplet")


def anchors():
    """ORF -> (evidence 'lipid'|'domain', default target compartment)"""
    out = {}
    u = pd.read_csv(f"{VD}/uniprot_features.tsv", sep="\t", dtype=str)
    for orf, lip in zip(u["Gene Names (ordered locus)"], u["Lipidation"]):
        if not isinstance(orf, str) or not isinstance(lip, str):
            continue
        orf = orf.split()[0].upper()
        notes = re.findall(r'/note="([^"]+)"', lip)
        if any("GPI-anchor" in n for n in notes):
            continue
        kinds = {k for n in notes for k in ("myristoyl", "palmitoyl", "farnesyl", "geranylgeranyl") if k in n}
        if kinds:
            out[orf] = ("lipid", "golgi_endosome" if kinds == {"geranylgeranyl"} else "plasma_membrane")
    m = pd.read_csv(f"{VD}/membrane_binding_domains_yeast.tsv", sep="\t", dtype=str)
    for orf, dom in zip(m["locus"], m["domain_hits"]):
        if isinstance(orf, str) and isinstance(dom, str) and MEM_DOMAINS.search(dom) and orf.upper() not in out:
            out[orf.upper()] = ("domain", "golgi_endosome" if ENDO_DOMAINS.search(dom) else "plasma_membrane")
    g = pd.read_csv(f"{VD}/gallego2010_lipid.tsv", sep="\t", skiprows=1, dtype=str)
    lipids = [c for c in g.columns[3:] if c.strip()]
    for _, r in g.iterrows():
        orf = str(r["ORF"]).upper()
        if orf in out:
            continue
        hits = {c for c in lipids if str(r[c]).strip() not in ("0", "", "nan")}
        if hits:
            out[orf] = ("domain", "golgi_endosome" if any(h.startswith("PtdIns(3)P") for h in hits) and
                        not any(h.startswith("PtdIns(4,5)P2") for h in hits) else "plasma_membrane")
    return out


def apply(w, ev, f):
    """move fraction f of a soluble protein's weights to a membrane's peripheral shell: an annotated membrane
    compartment if w has one (largest), else the anchor's default. w: dict compartment -> weight (sums to 1)"""
    kind, default = ev
    mem = {c: v for c, v in w.items() if c in MEMBRANE_COMPS}
    if sum(mem.values()) >= f:
        return w
    tgt = max(mem, key=mem.get) if mem else default
    out = {c: v * (1 - f) for c, v in w.items()}
    out[tgt] = out.get(tgt, 0.0) + f
    return out
