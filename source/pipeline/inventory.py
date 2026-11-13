"""Build a composition inventory for a haploid S. cerevisiae cell (S288C, log phase, rich glucose).

Layer 1 of the Maritan et al. 2021-style recipe: copies per cell, compartment, size.
"""
import re, json
import numpy as np
import pandas as pd

from config import DATA as D
# Proteome total: protein MASS in this cell; the relative abundances stay PaxDb's.
# Bulk measurements give 4-6 pg per haploid cell (Carroll et al. 2011 Mol Cell Proteomics 10:M111.007633, BNID 110550;
# Sherman 2002 Methods Enzymol 350:3, BNID 105082): averages over a growing culture, used here until 2026-11-12
# (5.0 pg in the culture's mean volume of 42 fL). The model cell is one unbudded G1 cell with one genome, and
# structures counted per cell say it holds less (validation/anchors.py, at 5.0 pg): core histones 1.3x what one
# genome's nucleosomes need (so 3.9 pg), nucleoporins 1.5x what the pores need (3.2 pg). 3.5 pg lies between the two;
# the ribosome count it gives (about 211,000, against ~187,000-200,000 measured) was not used and stays a check.
# The cell volume shrinks with it (CELL_VOL below) so that concentrations stay where they were.
POP_MEAN_PROTEIN_PG = 5.0
TOTAL_PROTEIN_PG = 3.5
CELL_SHARE = TOTAL_PROTEIN_PG / POP_MEAN_PROTEIN_PG     # this cell as a share of the population-average cell: also
                                                        # scales the tRNA and mRNA counts, which were per average cell

# ---------------- geometry (literature-anchored; see report) ----------------
CELL_VOL = 34.0          # um^3 (2026-11-12; was 42, the mean log-phase volume of Jorgensen 2002): the volume at which
                         # the cytosol keeps its macromolecular density (about 170 mg/mL) with 3.5 pg of protein and the
                         # fractions below, so the diffusion calibration stays valid. Between the culture mean (42) and
                         # unbudded G1 cells by soft X-ray tomography (24.9 um3, Chen 2025). Dry-mass density alone would
                         # put it at 36-38, the cryo-ET ribosome concentration at about 26 (MODEL.md 3.1).
WALL_THICK = 0.11        # um, lower bound of 110-200 nm (Orlean 2012)
FRAC = dict(nucleus=0.07, vacuole=0.172, mitochondrion=0.03, lipid_droplet=0.005)
# nucleus 7% (2026-11-12; was 8%): the nucleus-to-cell ratio is constant across cell sizes at about 7% (fluorescence,
#   Jorgensen 2007); soft X-ray tomography of unbudded G1 BY4741 in YPD gives 5.9 +- 1.3% (Chen 2025) and 5.3-6.3%
#   (Ohsawa 2026, preprint), the two segmented cryo-FIB/SEM nuclei 8.0% and 9.7% (cut cells, biased up), serial-section
#   EM 10.5% (Yamaguchi 2011)
# vacuole 17.2% (2026-11-12; was 9%): Chen 2025, the measurement matched in strain, medium and stage (17.2 +- 6.0%, 287
#   cells). The five segmented cells gave 6-13% (a lower bound: the segmentation flattens vacuoles), SXT 7% (Uchida 2011),
#   fluorescence 25% (Joyner 2016).
# lipid droplets 0.5%: median of five segmented cells (0.3-0.8%, depositor labels); SXT 1% (Uchida 2011)
# mitochondria 3% (2026-11-11; was 2%): measurements span 1.7% (serial-section EM, Yamaguchi 2011) to 4.77% (3D
#   fluorescence in glucose, Di Bartolomeo 2020 PNAS; diffraction inflates 300 nm tubules; CEN.PK). 3% is where the
#   requested matrix protein reaches 1.1x the cytosol's density, the soft X-ray absorption ratio of mitochondria to
#   cytosol (Chen 2025 MBoC): a fit to that ratio. At 2% the matrix asked for ~270 mg/mL and only 40% fitted.
NUCLEOLUS_OF_NUC = 0.20  # Uchida 2011
ENDOMEMBRANE_FRAC = 0.01 # assumption: Golgi+endosomes+peroxisomes (not segmented by SXT)
CER_COVER = 0.48         # fraction of PM covered by cortical ER (Wolf 2012)
NPC_PER_UM2 = 10.8       # pores per um2 of nuclear surface: the 117 pores (Winey 1997, mean over the cell cycle; 65-182,
                         # rising from G1) on the earlier 3.36 um3 nucleus. The count now follows the nucleus surface
                         # (an assumption; the G1 density of that study could not be read here: unchecked)
N_NPC = int(round(NPC_PER_UM2 * 4 * np.pi * (3 * FRAC["nucleus"] * CELL_VOL / (4 * np.pi)) ** (2 / 3)))


def geometry():
    R = (3 * CELL_VOL / (4 * np.pi)) ** (1 / 3)
    Rp = R - WALL_THICK
    g = {"cell_radius_um": R, "protoplast_radius_um": Rp, "cell_volume_um3": CELL_VOL}
    v = {}
    v["cell_wall"] = CELL_VOL - 4 / 3 * np.pi * Rp ** 3
    nuc = FRAC["nucleus"] * CELL_VOL
    v["nucleolus"] = NUCLEOLUS_OF_NUC * nuc
    v["nucleoplasm"] = nuc - v["nucleolus"]
    v["vacuole"] = FRAC["vacuole"] * CELL_VOL
    v["mitochondrion"] = FRAC["mitochondrion"] * CELL_VOL
    v["lipid_droplet"] = FRAC["lipid_droplet"] * CELL_VOL
    v["golgi_endosome"] = ENDOMEMBRANE_FRAC * CELL_VOL * 0.8
    v["peroxisome"] = ENDOMEMBRANE_FRAC * CELL_VOL * 0.2
    pm_area = 4 * np.pi * Rp ** 2
    cer_area = CER_COVER * pm_area
    v["ER"] = cer_area * 0.024 + 4 * np.pi * ((nuc * 3 / 4 / np.pi) ** (1 / 3)) ** 2 * 0.03  # cER sheet 24 nm + NE lumen ~30 nm
    v["nuclear_envelope"] = 0.0
    v["plasma_membrane"] = 0.0
    v["cytosol"] = 4 / 3 * np.pi * Rp ** 3 - sum(v[k] for k in v if k not in ("cell_wall",))
    rn = (nuc * 3 / 4 / np.pi) ** (1 / 3)
    areas = {"plasma_membrane": pm_area, "nuclear_envelope": 4 * np.pi * rn ** 2 * 2,
             "ER": cer_area * 2, "vacuole": 4 * np.pi * (v["vacuole"] * 3 / 4 / np.pi) ** (2 / 3)}
    g.update(volumes_um3=v, areas_um2=areas, nucleus_radius_um=rn, n_npc=N_NPC,
             cer_cover=CER_COVER, wall_thickness_um=WALL_THICK)
    return g


# ---------------- proteins this cell does not have (external review of 2026-11-12, comments 2 and 3) ----------------
# The abundances are PaxDb's integrated dataset, which mixes strains, media and stages. The model cell is BY4741
# (MATa his3-delta1 leu2-delta0 met15-delta0 ura3-delta0), haploid, unbudded G1-like, in glucose. Curated lists, not
# a measurement: genes the genotype deletes, alpha-specific genes, meiosis/sporulation-specific genes (a haploid
# cannot enter meiosis), glucose-repressed genes, and proteins degraded before G1. Their copies are set to zero before
# the proteome is scaled to its mass. A condition-matched abundance set would replace the last three lists.
NOT_IN_THIS_CELL = {
    "genotype": ["HIS3", "LEU2", "MET17", "URA3"],
    "alpha-specific": ["MF(ALPHA)1", "MF(ALPHA)2", "STE3", "SAG1", "MATALPHA1", "MATALPHA2", "HMLALPHA1", "HMLALPHA2"],
    "meiosis and sporulation": ["SPO11", "DMC1", "ZIP1", "ZIP2", "ZIP3", "HOP1", "HOP2", "MND1", "RED1", "REC8", "MEK1",
                                "MSH4", "MSH5", "NDT80", "IME1", "IME2", "SPO13", "SPO16", "SPS100", "SSP1", "SSP2",
                                "MPC54", "SPO21", "SPO74", "ADY3", "ADY4", "MEI4", "MEI5", "SAE3", "MER1", "REC102",
                                "REC104", "REC114", "SPR3", "SPR28", "DON1", "SMK1", "SPS1", "GIP1", "SPS2", "SPS22",
                                "SPO19", "SPR2"],
    "glucose-repressed": ["GAL1", "GAL2", "GAL7", "GAL10", "PCK1", "FBP1", "MLS1", "ICL1", "CYC7", "ADH2", "JEN1"],
    "absent in G1": ["MCD1", "CDC5", "SGO1", "PDS1", "CLB1", "CLB2", "CLB3", "CLB4"],
}
# "spindle pole" in UniProt or GO: only the components of the spindle pole body itself go to the nuclear envelope
# (spb.py then draws the plaque from them); for every other protein the term is ignored and its other annotations
# decide (calmodulin had 43,000 copies spread over the whole envelope; review comment 5)
SPB_CORE = {"SPC42", "SPC29", "SPC110", "CNM67", "NUD1", "SPC72", "SFI1", "KAR1", "CDC31", "MPS2", "MPS3", "BBP1",
            "NBP1", "TUB4", "SPC97", "SPC98"}


def _spindle_pole(a, gene):
    out = [("nuclear_envelope" if gene in SPB_CORE else None) if c == "spindle_pole" else c for c in a]
    return list(dict.fromkeys(c for c in out if c))


# ---------------- localization parsing ----------------
RULES = [  # (regex on lowercase text, compartment) - first match per location token
    (r"cell wall|secreted|extracellular|periplasm|cell surface", "cell_wall"),
    (r"nucleolus", "nucleolus"),
    (r"nuclear pore|nucleus membrane|nucleus inner membrane|nucleus outer membrane|nuclear envelope", "nuclear_envelope"),
    (r"spindle pole", "spindle_pole"),          # resolved per protein in build(): see SPB_CORE
    # mitochondrial sub-compartments (inventory v3); the generic term is the fallback
    (r"mitochondri\w* inner membrane|inner mitochondrial membrane|mitochondrial respiratory chain", "mito_IM"),
    (r"mitochondri\w* outer membrane|outer mitochondrial membrane", "mito_OM"),
    (r"mitochondri\w* intermembrane space", "mito_IMS"),
    (r"mitochondri\w* matrix|mitochondrial nucleoid", "mito_matrix"),
    (r"mitochondri", "mitochondrion"),
    (r"lipid droplet|lipid particle", "lipid_droplet"),
    (r"peroxisom", "peroxisome"),
    (r"vacuol", "vacuole"),
    (r"golgi|endosome|vesicle|prevacuolar|multivesicular|clathrin|trans-golgi|autophag", "golgi_endosome"),
    (r"endoplasmic reticulum|microsome", "ER"),
    (r"cell membrane|plasma membrane|eisosome|cell cortex", "plasma_membrane"),
    (r"nucleus|nucleoplasm|chromosome|chromatin|telomere|kinetochore", "nucleoplasm"),
    (r"cytoplasm|cytosol|cytoskeleton|p-body|stress granule|bud|ribosome|polysome|actin|microtubule", "cytosol"),
]


# Large-scale organelle proteomics papers that over-assign abundant cytosolic proteins
# (Sickmann 2003, Reinders 2006 mitochondrial proteomes). Location tokens supported ONLY by these are dropped.
HT_ORGANELLE_PROTEOMICS = {"14576278", "16823961"}


HUH2003 = "14562095"  # Huh et al. 2003 genome-wide GFP localization (steady-state bulk localization)
HT_MIN_ENTRIES = 40   # a paper cited as location evidence in this many UniProt entries counts as high-throughput
_PMID_ENTRIES = None


def _pmid_entries():
    """PubMed id -> number of UniProt entries citing it as subcellular-location evidence"""
    global _PMID_ENTRIES
    if _PMID_ENTRIES is None:
        import collections
        u = pd.read_csv(f"{D}/uniprot_yeast.tsv", sep="\t", usecols=["Subcellular location [CC]"])
        c = collections.Counter()
        for s in u.iloc[:, 0].dropna():
            c.update(set(re.findall(r"PubMed:(\d+)", s.split("Note=")[0])))
        _PMID_ENTRIES = c
    return _PMID_ENTRIES


def parse_uniprot_cc(s):
    """Return ordered list of compartments. If any token is supported by Huh et al. 2003 GFP data, only tokens with
    Huh support or experimental evidence from a dedicated (not high-throughput) study are kept, in UniProt order: the
    bulk steady-state localization without proteomics-only calls (GPD1 keeps its cytoplasm, PNC1 its cytoplasm and
    nucleus)."""
    if not isinstance(s, str):
        return []
    s = s.split("Note=")[0]
    s = s.replace("SUBCELLULAR LOCATION:", "")
    out, huh, has_huh = [], [], False
    ht = _pmid_entries()
    for m in re.finditer(r"([^.;{}]+?)\s*(\{[^}]*\})?\s*[.;]", s):
        t, ev = m.group(1).strip().lower(), (m.group(2) or "")
        if not t:
            continue
        pm = set(re.findall(r"PubMed:(\d+)", ev))
        if pm and pm <= HT_ORGANELLE_PROTEOMICS and "ECO:0000250" not in ev:
            continue
        for rx, c in RULES:
            if re.search(rx, t):
                out.append(c)
                has_huh |= HUH2003 in pm
                if HUH2003 in pm or ("ECO:0000269" in ev and any(ht[p] < HT_MIN_ENTRIES for p in pm)):
                    huh.append(c)
                break
    allc = out
    out = list(dict.fromkeys(huh if has_huh else out))
    return refine_mito(out, allc)


MITO_SUB = ("mito_IM", "mito_OM", "mito_IMS", "mito_matrix")
MITO_SECONDARY = 0.1   # weight factor of a secondary mitochondrial location behind a cytosolic/nuclear primary one


def refine_mito(comps, extra):
    """replace the generic 'mitochondrion' by the sub-compartments named in extra (e.g. a Huh 2003 GFP call of
    'mitochondrion' while UniProt also says 'mitochondrion inner membrane')"""
    if "mitochondrion" not in comps:
        return comps
    sub = [c for c in dict.fromkeys(extra) if c in MITO_SUB]
    if not sub:
        return comps
    out = []
    for c in comps:
        out += sub if c == "mitochondrion" else [c]
    return list(dict.fromkeys(out))


SOLUBLE = ("cytosol", "nucleoplasm", "nucleolus")


def tm_go_fallback(s):
    """membrane-bearing compartments named in GO CC (with 'cell periphery' read as plasma membrane), in GO order"""
    if not isinstance(s, str):
        return []
    out = []
    for t in s.lower().split(";"):
        t = re.sub(r"\[go:\d+\]", "", t).strip()
        if re.search(r"cell periphery", t):
            out.append("plasma_membrane"); continue
        for rx, c in RULES:
            if re.search(rx, t):
                if c not in SOLUBLE and c != "cell_wall":
                    out.append(c)
                break
    return list(dict.fromkeys(out))


def parse_go_cc(s):
    if not isinstance(s, str):
        return []
    out = []
    for t in s.lower().split(";"):
        t = re.sub(r"\[go:\d+\]", "", t).strip()
        if re.search(r"complex$", t) and not re.search(r"pore|ribosom|proteasome", t):
            continue
        for rx, c in RULES:
            if re.search(rx, t):
                out.append(c)
                break
    out = list(dict.fromkeys(out))
    # GO CC (no evidence filter available here) inherits high-throughput organelle calls for abundant
    # cytosolic enzymes; if cytosol is annotated, keep only cyto/nuclear terms.
    if "cytosol" in out:
        out = [c for c in out if c in ("cytosol", "nucleoplasm", "nucleolus")]
    return out


def weights(comps, vol, mito_hc=False):
    """Rank weights (1, 1/2, 1/4, ...) by annotation order (first-listed = primary); the cytosol and
    nucleoplasm shares are pooled and re-split by volume (equal concentration across the nuclear pore).
    mito_hc: in the high-confidence mitochondrial proteome (Morgenstern 2017)"""
    if not comps:
        return {"cytosol": 1.0}
    raw = {c: 0.5 ** i for i, c in enumerate(comps)}
    # dual-localised cytosolic/nuclear proteins with a secondary mitochondrial location (e.g. SOD1, ADK1 "Cytoplasm.
    # Mitochondrion intermembrane space", typically noted as a small percent): the mitochondrial pool is minor
    if comps[0] in ("cytosol", "nucleoplasm"):
        for c in comps[1:]:
            if c in MITO_SUB or c == "mitochondrion":
                raw[c] *= MITO_SECONDARY
    # and the reverse: a high-confidence mitochondrial protein with a secondary location elsewhere (PET9, the major
    # ADP/ATP carrier, "Peroxisome membrane" from one study): the other pool is minor
    if mito_hc and (comps[0] in MITO_SUB or comps[0] == "mitochondrion"):
        for c in comps[1:]:
            if not (c in MITO_SUB or c == "mitochondrion"):
                raw[c] *= MITO_SECONDARY
    cn = [c for c in comps if c in ("cytosol", "nucleoplasm")]
    if cn:
        pooled = sum(raw[c] for c in cn); tot = sum(vol[c] for c in cn)
        for c in cn:
            raw[c] = pooled * vol[c] / tot
    z = sum(raw.values())
    return {c: v / z for c, v in raw.items()}


# ---------------- measured localization (replaces annotation ranks where available) ----------------
VOGTLE_F = f"{D}/vogtle2017_submito.xlsx"     # Vogtle et al. 2017 Nat Commun 8:290, Supplementary Data 5
VOGTLE_SUB = {"Outer membrane": "mito_OM", "Inner membrane": "mito_IM", "Inter-membrane space": "mito_IMS", "Matrix": "mito_matrix"}


def vogtle_submito():
    """ORF -> mitochondrial sub-compartment, for proteins with exactly one SILAC-based assignment (810 proteins)"""
    v = pd.read_excel(VOGTLE_F, header=1)
    out = {}
    for _, r in v.iterrows():
        subs = [c for k, c in VOGTLE_SUB.items() if isinstance(r[k], str)]
        if len(subs) == 1:
            out[r["ORF"]] = subs[0]
    return out


def apply_vogtle(comps, sub):
    """replace every mitochondrial token by the measured sub-compartment"""
    if not any(c == "mitochondrion" or c in MITO_SUB for c in comps):
        return comps
    return list(dict.fromkeys(sub if (c == "mitochondrion" or c in MITO_SUB) else c for c in comps))


MORGENSTERN_F = f"{D}/morgenstern2017/TableS1_mmc2.xlsx"   # Morgenstern et al. 2017 Cell Rep 19:2836, Table S1


def morgenstern():
    """ORF -> (in the high-confidence mitochondrial proteome, enrichment class of the pure/crude comparison,
    submitochondrial class 'Matrix/IM' | 'IMS/IM' | 'OM' | None)"""
    m = pd.read_excel(MORGENSTERN_F, sheet_name="Results - Summary of Study", header=1)
    return {r["Gene Systematic Name"]: (r["High-confidence mito. proteome"] == 1, r["Classes"] if isinstance(r["Classes"], str) else None,
                                        r["Submitochondrial localization5"] if isinstance(r["Submitochondrial localization5"], str) else None)
            for _, r in m.iterrows()}


def apply_morgenstern(comps, rec, ntm, go):
    """(new comps, tag): drop the mitochondrial tokens of proteins that are outside the high-confidence proteome and
    depleted from pure mitochondria (Class 3; e.g. COPI coat, Sec7, Glt1); give generic 'mitochondrion' tokens the
    measured sub-compartment (Matrix/IM and IMS/IM: soluble proteins to the matrix or IMS, membrane proteins to the IM)"""
    hc, cls, sub = rec
    has_mito = any(c == "mitochondrion" or c in MITO_SUB for c in comps)
    if not has_mito:
        return comps, None
    if not hc and cls == "Class 3":
        rest = [c for c in comps if not (c == "mitochondrion" or c in MITO_SUB)]
        if not rest:
            rest = [c for c in parse_go_cc(go) if not (c == "mitochondrion" or c in MITO_SUB)] or ["cytosol"]
        return rest, "+Morgenstern_not_mito"
    if "mitochondrion" in comps and sub:
        tgt = {"OM": "mito_OM", "Matrix/IM": "mito_IM" if ntm > 0 else "mito_matrix",
               "IMS/IM": "mito_IM" if ntm > 0 else "mito_IMS"}.get(sub)
        if tgt:
            return list(dict.fromkeys(tgt if c == "mitochondrion" else c for c in comps)), "+Morgenstern_sub"
    return comps, None


CYCLOPS_FILES = [f"{D}/cyclops/WT{k}_LOCscore.xls" for k in (1, 2, 3)]   # Koh et al. 2015 G3, three wild-type screens
CYCLOPS_MIN = 0.05          # LOC-scores below this are ignored
CYCLOPS_MAP = {"Cytoplasm": "cytosol", "Nucleus": "nucleoplasm", "Spindle Pole": "nucleoplasm", "Nucleolus": "nucleolus",
               "Nuclear Periphery": "nuclear_envelope", "ER": "ER", "Golgi": "golgi_endosome", "Endosome": "golgi_endosome",
               "Mitochondria": "mitochondrion", "Peroxisome": "peroxisome", "Vacuole/Vacuolar Membrane": "vacuole",
               "Cortical Patches": "plasma_membrane", "Cell Periphery": "periphery", "Bud": "periphery",
               "Budneck": "periphery", "Budsite": "periphery"}
MEMBRANE_COMPS = ("plasma_membrane", "nuclear_envelope", "ER", "golgi_endosome", "vacuole", "peroxisome", "mitochondrion",
                  "lipid_droplet") + MITO_SUB
# compartment groups: CYCLoPs sets the share of each group; the split within a group (cytosol vs nucleus, nucleoplasm vs
# nucleolus, ER vs nuclear envelope, mitochondrial sub-compartments) keeps the annotation-based rule. Imaging is reliable
# at the level of groups, less so within them: the classifier calls many abundant cytosolic enzymes "Nucleus" (ERG10,
# TKL1), and LOC-scores are fractions of cells, not concentrations.
_GROUPS = {"cytosol": "nucleocytoplasm", "nucleoplasm": "nucleocytoplasm", "nucleolus": "nucleocytoplasm",
           "nuclear_envelope": "ER", "ER": "ER"}
GROUP = lambda c: "mitochondrion" if c in MITO_SUB else _GROUPS.get(c, c)
CYCLOPS_NEW_MIN = 0.5       # a group the annotation does not list is added only from this LOC-score on, and never the
                            # nucleocytoplasm (where GFP fusions that fail to import end up, e.g. Cyc1) or mitochondria
                            # for proteins outside the high-confidence mitochondrial proteome (punctate-class confusion)
RPROT = re.compile(r"^(RP[LS]\d+|RPP[012])[AB]?$")    # GFP-tagged r-proteins accumulate in the nucleus (not assembled)


def cyclops_scores():
    """ORF -> {CYCLoPs class: mean LOC-score over the wild-type screens that imaged the protein}"""
    d = pd.concat([pd.read_excel(f, header=3).rename(columns={"ORF": "orf", "Gene name": "gene"}) for f in CYCLOPS_FILES])
    cls = [c for c in d.columns if c in CYCLOPS_MAP]
    m = d.groupby("orf")[cls].mean()
    return {orf: {c: float(v) for c, v in row.items() if v >= CYCLOPS_MIN} for orf, row in m.iterrows()}


# Nucleus vs cytosol within the nucleocytoplasm group (IMPROVEMENTS H7a, 2026-11-11). The annotation rule pools the
# two and splits them by volume (equal concentration), which put proteins that imaging shows enriched in the nucleus
# 92% in the cytosol. The CYCLoPs Nucleus-minus-Cytoplasm LOC-score difference d predicts the independent N-terminal
# SWAT calls (Weill et al. 2018; dev half: d > 0.5 -> 91% nuclear-only), so the concentration ratio nucleus/cytosol
# is set to exp(NUC_K * d) for proteins with both scores. NUC_K was fitted on the dev half of the SWAT comparison.
NUC_K = 5.0


def nucleocyto_split(scores, inner, vol, k=None):
    """re-split the cytosol/nucleoplasm weights of `inner` (a group's annotation weights) by imaging enrichment"""
    k = NUC_K if k is None else k
    if not scores or "cytosol" not in inner or "nucleoplasm" not in inner or k <= 0:
        return inner
    d = scores.get("Nucleus", 0.0) - scores.get("Cytoplasm", 0.0)
    tot = inner["cytosol"] + inner["nucleoplasm"]
    en = float(np.exp(k * d))
    cn = en * vol["nucleoplasm"]; cc = vol["cytosol"]
    out = dict(inner)
    out["nucleoplasm"] = tot * cn / (cn + cc); out["cytosol"] = tot * cc / (cn + cc)
    return out


def cyclops_weights(scores, comps, ntm, vol, mito_ok=True):
    """compartment weights from CYCLoPs LOC-scores (fraction of cells showing the protein in a compartment), or None
    when they cannot be used. Each compartment group gets a share proportional to its summed scores; within a group
    the annotation-based weights (weights()) decide. 'Cell periphery' and bud calls count as plasma membrane for
    membrane proteins and proteins annotated there, and as cytosol otherwise (abundant cytosolic proteins look
    peripheral because the vacuole presses the cytoplasm against the cortex). Groups the annotation does not list are
    added only under CYCLOPS_NEW_MIN's conditions, through the same rules. Not used when the scores and the annotation share no
    group (C-terminal GFP mislocalises some secreted and lumenal proteins), for membrane proteins without a
    membrane-bearing call, or for secreted cell-wall proteins."""
    if not scores or "cell_wall" in comps:
        return None
    pm_ok = ntm > 0 or "plasma_membrane" in comps
    cls_c = {}
    for k, v in scores.items():
        c = CYCLOPS_MAP[k]
        if c == "periphery":
            c = "plasma_membrane" if pm_ok else "cytosol"
        cls_c[c] = cls_c.get(c, 0.0) + v
    gs = {}
    for c, v in cls_c.items():
        gs[GROUP(c)] = gs.get(GROUP(c), 0.0) + v
    ann_g = {}
    for c, f in (weights(comps, vol) if comps else {}).items():
        ann_g.setdefault(GROUP(c), {})[c] = f
    if comps and not set(ann_g) & set(gs):
        return None
    new_ok = lambda g, v: v >= CYCLOPS_NEW_MIN and g != "nucleocytoplasm" and (g != "mitochondrion" or mito_ok)
    gs = {g: v for g, v in gs.items() if g in ann_g or new_ok(g, v)}
    if not gs or (ntm > 0 and not any(g != "nucleocytoplasm" for g in gs)):
        return None
    z = sum(gs.values()); w = {}
    for g, v in gs.items():
        inner = ann_g[g] if g in ann_g else weights([c for c in cls_c if GROUP(c) == g], vol)
        if g == "nucleocytoplasm":
            inner = nucleocyto_split(scores, inner, vol)
        zi = sum(inner.values())
        for c, f in inner.items():
            w[c] = w.get(c, 0.0) + v / z * f / zi
    return w


def subunit_oligomer(s):
    if not isinstance(s, str):
        return 1
    s = s.lower()
    m = {"homodimer": 2, "homotrimer": 3, "homotetramer": 4, "homopentamer": 5, "homohexamer": 6,
         "homoheptamer": 7, "homooctamer": 8, "homodecamer": 10, "homododecamer": 12}
    for k, n in m.items():
        if re.search(r"\b" + k + r"\b", s.split(".")[0] + "." + s[:200]):
            return n
    return 1


# curated locations where the UniProt annotation is plainly incomplete
LOC_FIX = {
    # eisosome core BAR protein (GO eisosome; Pil1/Lsp1 lattice on the PM, Kefauver et al. 2024); UniProt lists only
    # "lipid droplet", from an organelle proteomics screen
    "PIL1": ["plasma_membrane"],
}


# membrane-bound fraction of soluble proteins with a lipid anchor / a membrane-binding domain (H6; fitted on dev)
BOUND_F = {"lipid": 0.0, "domain": 0.0}   # off: lipid anchors 0.7 helped dev by 0.1 point but not test (2026-11-11); domains hurt dev


def build():
    g = geometry()
    vol = dict(g["volumes_um3"])
    # Volumes for the concentration rules below (equal concentration of dual cytosol/nucleus proteins, the fitted nucleus
    # rule) as the cell is built: the lumen of the nuclear envelope lies inside the nucleus of geometry(), which
    # books it under "ER" and takes it off the cytosol. Until 2026-11-13 the nucleoplasm was 1.90 um3 here and 1.56
    # as built, so "equal concentration" came out 1.24 times enriched in the nucleus (stack review). Now 1.65 and
    # 18.44 against 1.56 and 18.42 built (the rest is pore channels and the smoothing of the nuclear surface).
    ne_lumen = 4 * np.pi * g["nucleus_radius_um"] ** 2 * 0.03
    vol["nucleoplasm"] -= ne_lumen; vol["cytosol"] += ne_lumen
    pax = pd.read_csv(f"{D}/paxdb_integrated.txt", sep="\t", comment="#",
                      names=["gene", "string_id", "ppm"])
    pax["orf"] = pax.string_id.str.replace("4932.", "", regex=False)
    pax = pax.groupby("orf", as_index=False).ppm.sum()
    u = pd.read_csv(f"{D}/uniprot_yeast.tsv", sep="\t")
    u.columns = ["uniprot", "reviewed", "orf", "gene", "length", "mass_da", "cc", "go_cc", "tm", "subunit", "name"]
    u["orf"] = u.orf.astype(str).str.split(r"[ ;]").str[0]
    u = u.drop_duplicates("orf")
    df = u.merge(pax, on="orf", how="left")
    matched_ppm = df.ppm.sum()
    df["ppm"] = df.ppm.fillna(0.0)
    absent = {g: k for k, gs in NOT_IN_THIS_CELL.items() for g in gs}
    m_abs = df.gene.isin(list(absent))
    absent_ppm = df.loc[m_abs].groupby(df.gene.map(absent)).ppm.agg(["count", "sum"])
    df.loc[m_abs, "ppm"] = 0.0
    # rescale so that total = Ho et al. 2018 total proteome size
    frac = df.ppm / df.ppm.sum()
    n_total = TOTAL_PROTEIN_PG * 1e-12 * 6.02214e23 / float((frac * df.mass_da.fillna(df.mass_da.median())).sum())
    df["copies_per_cell"] = frac * n_total
    df["n_tm"] = df.tm.fillna("").str.count("TRANSMEM")
    # integral (multi-pass) proteins without TRANSMEM ranges are outer-membrane beta-barrels (TOM40, SAM50, MDM10,
    # MDM34): flag them as membrane proteins; membrane_orient uses their beta-strands as membrane segments
    barrel = df.cc.fillna("").str.contains("Multi-pass membrane protein") & (df.n_tm == 0)
    df.loc[barrel, "n_tm"] = 1
    locs, src = [], []
    for cc, go, ntm, gene_ in zip(df.cc, df.go_cc, df.n_tm, df.gene):
        a = _spindle_pole(parse_uniprot_cc(cc), gene_)
        if "mitochondrion" in a:          # UniProt generic: take the sub-compartment from GO when it names one
            gc = parse_go_cc(go)
            if sum(c in MITO_SUB for c in gc) == 1:
                a = refine_mito(a, gc)
        s_ = "UniProt" if a else None
        if not a:
            a = _spindle_pole(parse_go_cc(go), gene_); s_ = "GO" if a else "default_cytosol"
        if ntm > 0 and all(c in SOLUBLE for c in a):
            # TM protein with only soluble annotations (UniProt often just says "Membrane"): membrane locations from GO
            gm = tm_go_fallback(go)
            if gm:
                a, s_ = gm, s_ + "+GO_membrane"
        locs.append(a); src.append(s_)
    for i, gene in enumerate(df.gene):
        if gene in LOC_FIX:
            locs[i] = list(LOC_FIX[gene]); src[i] = "curated"
    # measured sub-mitochondrial locations (Vogtle 2017) replace annotated or generic ones
    vg = vogtle_submito()
    for i, orf in enumerate(df.orf):
        if orf in vg:
            new = apply_vogtle(locs[i], vg[orf])
            if new != locs[i]:
                locs[i] = new; src[i] += "+Vogtle"
    # high-confidence mitochondrial proteome and submitochondrial classes (Morgenstern 2017); Vogtle takes precedence
    mg = morgenstern(); not_mito = set()
    for i, (orf, ntm, go) in enumerate(zip(df.orf, df.n_tm, df.go_cc)):
        if orf in mg and orf not in vg:
            new, tag = apply_morgenstern(locs[i], mg[orf], ntm, go)
            if tag:
                locs[i] = new; src[i] += tag
                if tag == "+Morgenstern_not_mito":
                    not_mito.add(orf)
    df["compartments"] = ["|".join(x) if x else "cytosol" for x in locs]
    df["oligomer_n"] = df.subunit.map(subunit_oligomer)
    # measured shares (CYCLoPs) replace the annotation ranks where the two sources agree
    cy = cyclops_scores()
    w = []
    for i, (orf, ntm, gene) in enumerate(zip(df.orf, df.n_tm, df.gene.fillna(""))):
        cw = None if RPROT.match(gene) else cyclops_weights(cy.get(orf), locs[i], ntm, vol,
                                                            mito_ok=orf in mg and mg[orf][0] and orf not in not_mito)
        if cw is None:
            w.append(weights(locs[i], vol, mito_hc=orf in mg and mg[orf][0] and orf not in not_mito))
            if orf in cy:
                src[i] += "+CYCLoPs_conflict"
        else:
            w.append(cw); src[i] += "+CYCLoPs"
    # membrane association from lipid anchors and lipid-binding domains (membrane_anchor.py, IMPROVEMENTS H6)
    import membrane_anchor as MA
    anc = MA.anchors()
    for i, (orf, ntm) in enumerate(zip(df.orf, df.n_tm)):
        ev = anc.get(str(orf).upper())
        if ev is None or ntm > 0 or "cell_wall" in locs[i] or BOUND_F.get(ev[0], 0) <= 0:
            continue
        new = MA.apply(w[i], ev, BOUND_F[ev[0]])
        if new is not w[i]:
            w[i] = new; src[i] += f"+anchor_{ev[0]}"
    df["loc_source"] = src
    df["comp_weights"] = [json.dumps({k: round(v, 4) for k, v in x.items()}) for x in w]
    df["gene"] = df.gene.fillna(df.orf)
    keep = ["orf", "gene", "uniprot", "name", "length", "mass_da", "copies_per_cell", "ppm",
            "compartments", "comp_weights", "loc_source", "n_tm", "oligomer_n", "subunit"]
    df = df[keep].sort_values("copies_per_cell", ascending=False).reset_index(drop=True)
    df.attrs["paxdb_ppm_matched"] = matched_ppm
    df.attrs["not_in_this_cell"] = {k: dict(proteins=int(r["count"]), share_of_ppm=float(r["sum"] / matched_ppm))
                                    for k, r in absent_ppm.iterrows()}
    return df, g


def per_compartment(df, g):
    rows = {}
    for cw, n, m in zip(df.comp_weights, df.copies_per_cell, df.mass_da):
        for c, w in json.loads(cw).items():
            r = rows.setdefault(c, [0.0, 0.0])
            r[0] += n * w; r[1] += n * w * m
    out = pd.DataFrame([(c, a, b) for c, (a, b) in rows.items()], columns=["compartment", "protein_copies", "protein_mass_da"])
    out["volume_um3"] = out.compartment.map(g["volumes_um3"])
    out["area_um2"] = out.compartment.map(g["areas_um2"])
    out["protein_mg_per_ml"] = out.protein_mass_da / 6.022e23 / (out.volume_um3 * 1e-15)  # g/L == mg/mL
    return out.sort_values("protein_copies", ascending=False)
