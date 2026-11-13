"""Mitoribosomes on the inner membrane (literature item L6).

In yeast, the vast majority of mitoribosomes are tethered to the inner membrane, with the polypeptide exit tunnel
facing it, through Mba1 and the 21S rRNA expansion segment 96-ES1 (cryo-ET of isolated mitochondria: Pfeffer et al.
2015 Nat Commun 6:6019). They make the hydrophobic core subunits of the respiratory chain and ATP synthase, which are
inserted co-translationally by Oxa1.

Structure: the complete 74S monosome of 9HM0 (yeast mitoribosome with mRNA; the translational activator Aep3 is
left out). Membrane frame: the exit-tunnel axis (from the large subunit's centre to the centroid of the exit-site
proteins uL23m MRP20, uL24m MRPL40 and uL29m MRPL4) is the membrane normal, pointing from the matrix into the membrane
(the matrix is the positive, N side of the inner membrane in the library frame), and the monosome is lowered until its
nearest bead is GAP nm above the matrix-side leaflet (2 nm half-thickness). No residue enters the bilayer, so the
packer treats it as a membrane assembly that sits on the membrane.

Copies: MEMBRANE_FRACTION of the smaller of the large- and small-subunit pools (assumption from "vast majority"); the
remaining subunits stay free in the matrix, as before.
"""
import os, json
import numpy as np
import pandas as pd
import gemmi
import extend_complexes as xc
import assemblies as A
import membrane_orient as mo
import build_structs as bs
from config import LIB

PDB = "9HM0"
LIB_ID = "mitoribosome_74S_IM"
LSU, SSU = "CPX_1602_54S_mitochondrial_large_ribosoma", "CPX_1603_37S_mitochondrial_small_ribosoma"
EXIT_GENES = ("MRP20", "MRPL40", "MRPL4")      # uL23m, uL24m, uL29m around the tunnel exit
DROP_GENES = ("AEP3",)
MEMBRANE_FRACTION = 0.9
HALF = 2.0                                      # membrane half-thickness (cell3d.HALF)
GAP = 1.0                                       # nm between the nearest bead and the leaflet surface
INFO_F = f"{LIB}/{LIB_ID}.json"


def build(acc_of, log=print):
    """library entry in the inner-membrane frame; returns info dict"""
    if os.path.exists(INFO_F) and os.path.exists(f"{LIB}/{LIB_ID}.npz"):
        return json.load(open(INFO_F))
    st = gemmi.read_structure(xc.fetch_assembly(PDB)); st.remove_ligands_and_waters(); st.setup_entities()
    st.remove_empty_chains()
    xyz, kind, acc, pos = mo.beads_with_positions(st, align=A.entity_align(PDB))
    drop = np.isin(acc, [acc_of[g] for g in DROP_GENES if g in acc_of])
    xyz, kind, acc = xyz[~drop], kind[~drop], acc[~drop]
    exit_acc = [acc_of[g] for g in EXIT_GENES]
    ex = np.isin(acc, exit_acc)
    c = xyz.mean(0); e = xyz[ex].mean(0)
    n = (c - e) / np.linalg.norm(c - e)             # from the membrane into the matrix (+z, N side)
    # orthonormal frame with z = n
    a = np.array([1.0, 0, 0]) if abs(n[0]) < 0.9 else np.array([0, 1.0, 0])
    u = np.cross(n, a); u /= np.linalg.norm(u); v = np.cross(n, u)
    X = np.c_[(xyz - e) @ u, (xyz - e) @ v, (xyz - e) @ n]
    shift = (HALF + GAP) - X[:, 2].min()            # nearest bead GAP above the leaflet; everything at z > 0 (matrix)
    X[:, 2] += shift
    bs.save(LIB_ID, X, kind, extra=dict(pdb=PDB, membrane="TM", tm=np.zeros(len(X), bool)), per=30)
    info = dict(n_beads=int(len(X)), n_nt=int((kind == 1).sum()), exit_z=float(shift), centre_z=float(X[:, 2].mean()),
                exit_residues=int(ex.sum()))
    near = X[:, 2] < X[:, 2].min() + 3.0
    info["nearest_kind_nt"] = float((kind[near] == 1).mean())
    json.dump(info, open(INFO_F, "w"))
    log(f"  {LIB_ID} ({PDB}): {len(X)} beads ({info['n_nt']} nt); exit site {shift - HALF:.1f} nm above the leaflet, "
        f"centre {info['centre_z']:.1f} nm from the midplane; beads within 3 nm of the membrane-facing side: "
        f"{info['nearest_kind_nt']:.0%} RNA")
    return info


def assign(ing, df, log=print):
    """replace MEMBRANE_FRACTION of the mitoribosomal subunit pools by membrane-bound 74S monosomes"""
    acc_of = dict(zip(df.gene, df.uniprot))
    info = build(acc_of, log=log)
    i = ing.set_index("ing_id")
    if LSU not in i.index or SSU not in i.index:
        log("  mitoribosomes: subunit complexes not in the inventory; skipped")
        return ing, None
    nl, ns = float(i.loc[LSU, "copies_per_cell"]), float(i.loc[SSU, "copies_per_cell"])
    n = MEMBRANE_FRACTION * min(nl, ns)
    genes = sorted(set(i.loc[LSU, "genes"].split(";")) | set(i.loc[SSU, "genes"].split(";")))
    row = dict(ing_id=LIB_ID, label="74S mitoribosome on the inner membrane", kind="complex", genes=";".join(g for g in genes if g),
               n_subunit_types=len(genes), copies_per_cell=n, comp_weights=json.dumps({"mito_IM": 1.0}),
               mass_da=float(i.loc[LSU, "mass_da"]) + float(i.loc[SSU, "mass_da"]), pdb=PDB, membrane="TM",
               oligomer_n=1, uniprot="", n_tm=1)
    ing = ing.copy()
    ing.loc[ing.ing_id == LSU, "copies_per_cell"] = nl - n
    ing.loc[ing.ing_id == SSU, "copies_per_cell"] = ns - n
    ing = ing[ing.copies_per_cell > 0.5]
    log(f"mitoribosomes: {n:.0f} membrane-bound 74S ({MEMBRANE_FRACTION:.0%} of min(LSU {nl:.0f}, SSU {ns:.0f})); "
        f"free subunits left: LSU {nl - n:.0f}, SSU {ns - n:.0f}")
    return pd.concat([ing, pd.DataFrame([row])], ignore_index=True), dict(n=n, info=info)
