"""Pre-ribosomal particles: ribosome assembly intermediates from the nucleolus to the cytoplasm.

A growing yeast cell makes about 2,000-4,000 ribosomes per minute (Warner 1999). The small subunit assembles on the
nascent 35S pre-rRNA in the nucleolus as the 90S pre-ribosome (SSU processome), and leaves the nucleus early as a
pre-40S particle that matures in the cytoplasm. The large subunit passes through nucleolar, nucleoplasmic and
cytoplasmic pre-60S states. One cryo-EM structure per stage, with all of its RNA (pre-rRNA, U3 snoRNA):

  preribo_90S              5WLC  complete SSU processome (Barandun et al. 2017)              nucleolus
  preribo_60S_nucleolar    7R7A  nucleolar pre-60S, state E1 (Cruz et al. 2022)              nucleolus
  preribo_60S_nucleoplasm  6YLH  Rix1-Rea1 nucleoplasmic pre-60S (Kater et al. 2020)        nucleoplasm
  preribo_60S_cytoplasm    5T62  Nmd3-Lsg1-Tif6 cytoplasmic pre-60S (Malyutin et al. 2017)  cytosol
  preribo_40S              6FAI  cytoplasmic pre-40S (Scaiola et al. 2018)                  cytosol, partly nucleoplasm

Copies: the number of particles in a stage is set by its assembly factors, the same way as for other complexes: the
robust median over factors found only in that stage's structure (its markers), capped at twice the 25th percentile
over all of its factors. The particles consume their factors from the protein pools (complexes built later get what
is left). Ribosomal proteins are not consumed from the pool; instead the mature 80S count drops by the r-protein
content of the particles (as a fraction of an 80S), so ribosomal-protein mass is conserved.
"""
import os, re, json
import numpy as np
import pandas as pd
import gemmi
import extend_complexes as xc
import assemblies as A
import membrane_orient as mo
import build_structs as bs
from config import LIB

RIBO = "ribosome_80S"
STAGES = [   # (ing_id, label, PDB, compartment weights)
    ("preribo_90S", "90S pre-ribosome (SSU processome)", "5WLC", {"nucleolus": 1.0}),
    ("preribo_60S_nucleolar", "Nucleolar pre-60S ribosome (state E)", "7R7A", {"nucleolus": 1.0}),
    ("preribo_60S_nucleoplasm", "Nucleoplasmic pre-60S ribosome (Rix1-Rea1)", "6YLH", {"nucleoplasm": 1.0}),
    ("preribo_60S_cytoplasm", "Cytoplasmic pre-60S ribosome (Nmd3, Lsg1)", "5T62", {"cytosol": 1.0}),
    # pre-40S particles are exported soon after the 90S is released; assumption: 15% of them are still nuclear
    ("preribo_40S", "Pre-40S ribosome", "6FAI", {"cytosol": 0.85, "nucleoplasm": 0.15}),
]
RPROT = re.compile(r"^(RP[LS]\d+|RPP[012])[AB]?$|^ASC1$")
INFO_F = f"{LIB}/preribosomes.json"


def _content(pid, gene_of):
    """per structure: beads (all chains, RNA included), residue counts per gene, nucleotides"""
    st = gemmi.read_structure(xc.fetch_assembly(pid)); st.remove_ligands_and_waters(); st.setup_entities()
    st.remove_empty_chains()
    xyz, kind, acc, pos = mo.beads_with_positions(st, align=A.entity_align(pid))
    ea = A.entity_uniprot(pid)
    chains = {}
    for ch in st[0]:
        ent = st.get_entity_of(ch.get_polymer())
        if ent is None or ent.polymer_type == gemmi.PolymerType.Rna:
            continue
        for a in ea.get(ent.name, ()):
            g = gene_of.get(a)
            if g:
                chains[g] = chains.get(g, 0) + 1
    res = {}
    for a in acc[kind == 0]:
        g = gene_of.get(a)
        if g:
            res[g] = res.get(g, 0) + 1
    return xyz, kind, chains, res


def build(df, log=print):
    """library entries; returns {ing_id: dict(sub={gene: chains}, n_nt, rprot_res, n_beads)}"""
    if os.path.exists(INFO_F) and all(os.path.exists(f"{LIB}/{s[0]}.npz") for s in STAGES):
        return json.load(open(INFO_F))
    gene_of = dict(zip(df.uniprot, df.gene))
    info = {}
    for lid, label, pid, _ in STAGES:
        xyz, kind, chains, res = _content(pid, gene_of)
        bs.save(lid, xyz - xyz.mean(0), kind, extra={"pdb": pid}, per=30)
        info[lid] = dict(sub=chains, n_nt=int((kind == 1).sum()), n_aa=int((kind == 0).sum()),
                         rprot_res=int(sum(n for g, n in res.items() if RPROT.match(g))))
        log(f"  {lid} ({pid}): {len(chains)} protein species, {info[lid]['n_aa']} residues, {info[lid]['n_nt']} nt, "
            f"radius {np.linalg.norm(xyz - xyz.mean(0), axis=1).max():.1f} nm")
    json.dump(info, open(INFO_F, "w"))
    return info


def assign(ing, df, log=print):
    """add the pre-ribosome ingredients, consuming their assembly factors and the matching share of mature ribosomes"""
    info = build(df, log=log)
    mass = dict(zip(df.gene, df.mass_da))
    pool = {r.ing_id: r.copies_per_cell for r in ing[ing.kind == "protein"].itertuples()}
    factors = {lid: {g: s for g, s in info[lid]["sub"].items() if not RPROT.match(g)} for lid in info}
    ribo_res = int((np.load(f"{LIB}/{RIBO}.npz")["kind"] == 0).sum())
    rows, used, ribo_equiv = [], {}, 0.0
    for lid, label, pid, cw in STAGES:
        fac = factors[lid]
        others = set().union(*(set(factors[o]) for o in factors if o != lid))
        markers = {g: s for g, s in fac.items() if g not in others and pool.get(g, 0) > 0} or fac
        n = min(A.complex_copies(pool, markers), A.complex_copies(pool, fac))
        for g, s in fac.items():
            take = min(pool.get(g, 0), n * s)
            pool[g] = pool.get(g, 0) - take; used[g] = used.get(g, 0) + take
        f_ribo = info[lid]["rprot_res"] / ribo_res
        ribo_equiv += n * f_ribo
        genes = sorted(info[lid]["sub"])
        rows.append(dict(ing_id=lid, label=label, kind="complex", genes=";".join(genes), n_subunit_types=len(genes),
                         copies_per_cell=n, comp_weights=json.dumps(cw), pdb=pid, membrane="", oligomer_n=1, uniprot="", n_tm=0,
                         mass_da=sum(mass.get(g, 0.0) * s for g, s in info[lid]["sub"].items()) + 330.0 * info[lid]["n_nt"]))
        log(f"  {lid}: {n:.0f} particles from {len(markers)} stage markers "
            f"({', '.join(sorted(markers)[:6])}{'...' if len(markers) > 6 else ''}); r-proteins = {f_ribo:.2f} of an 80S")
    ing = ing.copy()
    for g, k in used.items():
        ing.loc[(ing.ing_id == g) & (ing.kind == "protein"), "copies_per_cell"] -= k
    r = ing.ing_id == RIBO
    n80 = float(ing.loc[r, "copies_per_cell"].iloc[0])
    ing.loc[r, "copies_per_cell"] = n80 - ribo_equiv
    ing = ing[ing.copies_per_cell > 0.5]
    n_pre = sum(x["copies_per_cell"] for x in rows)
    log(f"pre-ribosomes: {n_pre:.0f} particles ({n_pre / n80:.1%} of {n80:.0f} ribosomes); mature 80S reduced by "
        f"{ribo_equiv:.0f} (their r-protein content); {sum(used.values()):.0f} factor copies consumed")
    return pd.concat([ing, pd.DataFrame(rows)], ignore_index=True), dict(n=n_pre, ribo_equiv=ribo_equiv, used=used)
