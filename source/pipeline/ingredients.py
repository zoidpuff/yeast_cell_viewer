"""Assemble packing ingredients: curated macromolecular complexes (PDB assemblies) + residual
proteins (AlphaFold monomers, or PDB homo-oligomers for abundant proteins)."""
import json, re
import numpy as np
import pandas as pd

from config import DATA as D
GENOME_BP = 12.07e6      # S288C nuclear genome (chromosomal, excl. rDNA repeat expansion)
NRL_BP = 165             # assumed nucleosome repeat length

# curated complexes: (id, label, complex-portal ids or explicit gene->stoich, PDB assembly, compartment weights, membrane)
COMPLEXES = [
    ("ribosome_80S", "80S ribosome", "RIBO", "4V88", {"cytosol": 0.97, "ER": 0.03}, None),
    # nucleus/cytosol split from the subunits' measured weights (CYCLoPs, MODEL.md 4.2): the 26S proteasome is
    # assembled in and mostly localized to the nucleus in growing cells (Laporte et al. 2008 J Cell Biol 181:737)
    ("proteasome_26S", "26S proteasome", ["CPX-2262"], "4CR2", ("from_subunits", ("nucleoplasm", "cytosol")), None),
    ("FAS", "Fatty acid synthase (a6b6)", {"FAS1": 6, "FAS2": 6}, "2UV8", {"cytosol": 1.0}, None),
    ("CCT", "TRiC/CCT chaperonin", ["CPX-2156"], "4V94", {"cytosol": 1.0}, None),
    ("PolII", "RNA polymerase II", ["CPX-2662"], "1WCM", {"nucleoplasm": 1.0}, None),
    ("PolI", "RNA polymerase I", ["CPX-1664"], "4C2M", {"nucleolus": 1.0}, None),
    ("PolIII", "RNA polymerase III", ["CPX-2660"], "5FJ8", {"nucleoplasm": 1.0}, None),
    ("NPC", "Nuclear pore complex", "NPC_KIM2018", "7N9F", {"nuclear_envelope": 1.0}, "NE"),
    ("Pma1_hexamer", "Pma1 H+-ATPase hexamer", {"PMA": (6, ["PMA1", "PMA2"])}, "7VH5", {"plasma_membrane": 1.0}, "PM"),
    ("V-ATPase", "Vacuolar V-ATPase V1Vo", {"VMA1": 3, "VMA2": 3, "VMA4": 3, "VMA10": 3, "VMA8": 1, "VMA7": 1,
                                             "VMA5": 1, "VMA13": 1, "VPH1": 1, "VMA3": 8, "VMA11": 1, "VMA16": 1,
                                             "VMA6": 1, "VMA9": 1}, "6O7V", {"vacuole": 1.0}, "VAC"),
    ("nucleosome", "Nucleosome core particle", {"H3": (2, ["HHT1; HHT2"]), "H4": (2, ["HHF1; HHF2"]),
                                                "H2A": (2, ["HTA1", "HTA2"]), "H2B": (2, ["HTB1", "HTB2"])},
     "1ID3", {"nucleoplasm": 1.0}, None),
]
# chains of the 80S structure 4V88 beyond the RPL/RPS/RPP proteins (found by the PDB audit of 2026-11-12: both were
# also placed as free proteins, so counted twice). Asc1 (RACK1) is a 40S component, one per ribosome. Stm1 is in the
# crystal because it clamps the idle ribosomes of the starved cells it came from; its 188,000 copies are fewer than the
# ribosomes, so all of them are booked to ribosomes and none is placed free (the library structure still has one Stm1
# in every ribosome; a translating 80S without it is a library change not made yet).
RIBO_EXTRA = ("ASC1", "STM1")
# nucleoporins beyond the assembled pores (review comment 8): PaxDb gives about three pores' worth per pore. Scaffold
# and FG nucleoporins exist almost only in pores, so their free pool is capped at NUP_FREE of the copies in pores (an
# assumption) and the surplus is not placed (reported by s1). Sec13, Seh1 (also COPII and SEA), Dyn2 and the basket
# proteins Mlp1/2 are left alone.
NUP_CAPPED = re.compile(r"^(NUP\d+|NSP1|NIC96|NDC1|POM\d+|ASM4|GLE1)$")
NUP_FREE = 0.25
import inventory as _inv
FIXED_COPIES = {"NPC": _inv.N_NPC, "nucleosome": GENOME_BP / NRL_BP}        # pores: inventory.N_NPC (by areal density)
NUC_FROM_MAP = True      # nucleosomes on the arms = measured dyads (nucleosome_map.py; 66,245, 10% below GENOME_BP / NRL_BP)


# copies of each nucleoporin per pore, measured by quantitative mass spectrometry and calibrated imaging (Kim et al.
# 2018 Nature 555:475, Supplementary Table 2A, column "Copies per NPC"; Nup145 is one gene for Nup145N and Nup145C).
# Until 2026-11-12 the pore took its subunits from Complex Portal CPX-824, which gives no stoichiometry, so one copy
# of each nucleoporin per pore was booked (117 per cell instead of 936-5,616) and the rest was placed as free protein.
NPC_COPIES = {"GLE1": 8, "NUP42": 8, "NUP116": 16, "NUP100": 16, "DYN2": 16, "NUP82": 16, "NUP159": 16, "NSP1": 48,
              "NIC96": 32, "NUP49": 32, "NUP57": 32, "NUP84": 16, "NUP85": 16, "NUP120": 16, "NUP133": 16, "NUP145": 16,
              "SEH1": 16, "SEC13": 16, "NUP157": 16, "NUP170": 16, "NUP188": 16, "NUP192": 16, "NUP53": 16, "ASM4": 16,
              "NDC1": 16, "POM34": 16, "POM152": 16, "NUP1": 8, "NUP60": 16, "MLP1": 8, "MLP2": 8}


def complex_families(cp, spec, gene_of_acc, genes):
    """Return {family: (stoich, [paralog genes])}."""
    if spec == "NPC_KIM2018":
        return {g: (n, [g]) for g, n in NPC_COPIES.items() if g in genes}
    if spec == "RIBO":
        fam = {}
        for g in genes:
            if re.match(r"^(RP[LS]\d+|RPP[012])[AB]?$", g):
                fam.setdefault(re.sub(r"[AB]$", "", g), []).append(g)
            elif g in RIBO_EXTRA:
                fam[g] = [g]
        return {f: (1, gs) for f, gs in fam.items()}
    if isinstance(spec, dict):
        return {k: (v if isinstance(v, tuple) else (v, [k])) for k, v in spec.items()}
    out = {}
    for cid in spec:
        row = cp[cp["#Complex ac"] == cid].iloc[0]
        for tok in row["Identifiers (and stoichiometry) of molecules in complex"].split("|"):
            m = re.match(r"([A-Z0-9]+)(?:-\w+)?\((\d+)\)", tok)
            if m and gene_of_acc.get(m.group(1)):
                g = gene_of_acc[m.group(1)]
                out[g] = (max(1, int(m.group(2))), [g])
    return out


def subunit_weights(df, genes, keep):
    """copy-weighted mean compartment weights of a complex's subunits, restricted to `keep` and renormalised"""
    acc = {k: 0.0 for k in keep}
    for g in genes:
        r = df[df.gene == g]
        if len(r) == 0:
            continue
        w = json.loads(r.comp_weights.iloc[0]); n = float(r.copies_per_cell.iloc[0])
        for k in keep:
            acc[k] += n * w.get(k, 0.0)
    z = sum(acc.values())
    return {k: round(v / z, 4) for k, v in acc.items()} if z > 0 else {k: 1.0 / len(keep) for k in keep}


def build_ingredients(df):
    cp = pd.read_csv(f"{D}/complexportal_559292.tsv", sep="\t")
    gene_of_acc = dict(zip(df.uniprot, df.gene))
    pool = dict(zip(df.gene, df.copies_per_cell))
    mass = dict(zip(df.gene, df.mass_da))
    rows, consumed = [], {}
    for cid, label, spec, pdb, cw, memb in COMPLEXES:
        fams = complex_families(cp, spec, gene_of_acc, list(pool))
        fam_copies = {f: sum(pool.get(g, 0) for g in gs) / st for f, (st, gs) in fams.items()}
        if cid in FIXED_COPIES:
            copies = FIXED_COPIES[cid]
            if cid == "nucleosome" and NUC_FROM_MAP:
                import nucleosome_map
                copies = nucleosome_map.arm_count()
            if cid == "nucleosome":         # plus the nucleosomal rDNA (inactive repeats and spacers), not in GENOME_BP
                import rdna
                copies += rdna.nucleosomes_in_rdna()
        elif cid == "Pma1_hexamer":
            copies = fam_copies["PMA"]
        else:
            vals = [v for v in fam_copies.values() if v > 0]
            copies = float(np.median(vals)) if vals else 0.0
        m_tot = 0.0
        for f, (st, gs) in fams.items():
            tot = sum(pool.get(g, 0) for g in gs)
            m_tot += st * np.mean([mass.get(g, 0) for g in gs])
            for g in gs:
                share = pool.get(g, 0) / tot if tot > 0 else 1 / len(gs)
                consumed[g] = consumed.get(g, 0) + copies * st * share
        genes_all = sorted(g for _, gs in fams.values() for g in gs)
        if isinstance(cw, tuple) and cw[0] == "from_subunits":
            cw = subunit_weights(df, genes_all, cw[1])
        rows.append(dict(ing_id=cid, label=label, kind="complex", genes=";".join(genes_all), n_subunit_types=len(fams),
                         copies_per_cell=copies, comp_weights=json.dumps(cw), mass_da=m_tot, pdb=pdb,
                         membrane=memb or "", oligomer_n=1, uniprot=""))
    # residual monomer pool
    npc_genes = set(next(x["genes"] for x in rows if x["ing_id"] == "NPC").split(";"))
    surplus = [0, 0.0, 0.0]
    for r in df.itertuples():
        left = max(0.0, r.copies_per_cell - consumed.get(r.gene, 0.0))
        if r.gene in npc_genes and NUP_CAPPED.match(str(r.gene)) and left > NUP_FREE * consumed.get(r.gene, 0.0):
            cap = NUP_FREE * consumed.get(r.gene, 0.0)
            surplus[0] += 1; surplus[1] += left - cap; surplus[2] += (left - cap) * r.mass_da; left = cap
        if left < 1:
            continue
        rows.append(dict(ing_id=r.gene, label=r.name, kind="protein", genes=r.gene, n_subunit_types=1,
                         copies_per_cell=left, comp_weights=r.comp_weights, mass_da=r.mass_da, pdb="",
                         membrane="", oligomer_n=r.oligomer_n, uniprot=r.uniprot, n_tm=r.n_tm))
    ing = pd.DataFrame(rows)
    ing["n_tm"] = ing.n_tm.fillna(0).astype(int)
    print(f"nucleoporin surplus not placed: {surplus[1]:.0f} copies of {surplus[0]} nucleoporins, {surplus[2] / 1e9:.2f} GDa "
          f"(free pool capped at {NUP_FREE:.0%} of the copies in pores)", flush=True)
    return ing, consumed
